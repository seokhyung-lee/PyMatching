// Copyright 2026 PyMatching fork contributors.
// Modified PyMatching file: Stim-free reusable whole-graph reweight/decode API.
// Licensed under the Apache License, Version 2.0. See LICENSE.

#include "pymatching/core/matching_region.h"
#include "pymatching/core/matching_reuse.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

#include "pymatching/sparse_blossom/flooder/graph_flooder.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "pymatching/sparse_blossom/search/search_flooder.h"

namespace pm::core {
namespace {

using EdgeKey = std::pair<uint32_t, int64_t>;

EdgeKey canonical_key(const ClassEdge& edge, size_t num_detectors) {
    if (edge.node1 >= num_detectors) {
        throw std::invalid_argument("class edge node1 is out of range");
    }
    if (edge.node2 == BOUNDARY_NODE) {
        return {edge.node1, BOUNDARY_NODE};
    }
    if (edge.node2 < 0 || static_cast<uint64_t>(edge.node2) >= num_detectors) {
        throw std::invalid_argument("class edge node2 is out of range");
    }
    if (edge.node1 == static_cast<uint32_t>(edge.node2)) {
        throw std::invalid_argument("self-loop class edges are not supported");
    }
    return {std::min(edge.node1, static_cast<uint32_t>(edge.node2)),
            std::max(edge.node1, static_cast<uint32_t>(edge.node2))};
}

void process_timeline(Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    if (!mwpm.flooder.queue.empty()) {
        throw std::logic_error("matcher queue was not cleared before decode");
    }
    mwpm.flooder.queue.cur_time = 0;
    for (uint64_t detection : detection_events) {
        auto* node = &mwpm.flooder.graph.nodes[detection];
        if (!node->region_that_arrived) {
            mwpm.create_detection_event(node);
        } else if (node->reached_from_source != node) {
            throw std::logic_error("retained detection has the wrong source");
        }
    }
    while (true) {
        auto event = mwpm.flooder.run_until_next_mwpm_notification();
        if (event.event_type == NO_EVENT) {
            break;
        }
        mwpm.process_event(event);
    }
    if (mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size()) {
        mwpm.reset();
        throw std::invalid_argument("no perfect matching exists for the supplied syndrome");
    }
}

void extract_match_edges(Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    mwpm.flooder.match_edges.clear();
    for (uint64_t detection : detection_events) {
        auto& node = mwpm.flooder.graph.nodes[detection];
        if (node.region_that_arrived != nullptr) {
            mwpm.shatter_blossom_and_extract_match_edges(node.region_that_arrived_top, mwpm.flooder.match_edges);
        }
    }
}

}  // namespace

struct MatchingRegion::Impl {
    struct TopologyEdge {
        EdgeKey key;
        std::vector<size_t> class_indices;
        weight_int* matching_weight = nullptr;
        weight_int* matching_reverse_weight = nullptr;
        weight_int* search_weight = nullptr;
        weight_int* search_reverse_weight = nullptr;
        size_t selected_class = 0;
        uint32_t external_weight = 0;
    };

    size_t detector_count;
    std::vector<ClassEdge> class_edges;
    std::vector<TopologyEdge> topology_edges;
    std::vector<std::vector<size_t>> search_topology_index;
    Mwpm matcher;
    const MatchingOptions options;
    bool mature = false;
    std::vector<uint8_t> previous_syndrome;
    std::vector<double> previous_float_weights;
    std::vector<uint32_t> retained_quantized;
    double previous_scale = 0;
    double previous_maximum = 0;

    Impl(size_t num_detectors, std::vector<ClassEdge> edges, MatchingOptions selected_options)
        : detector_count(num_detectors), class_edges(std::move(edges)), options(selected_options) {
        std::map<EdgeKey, std::vector<size_t>> groups;
        for (size_t i = 0; i < class_edges.size(); ++i) {
            groups[canonical_key(class_edges[i], detector_count)].push_back(i);
        }
        MatchingGraph matching_graph(detector_count, 0);
        SearchGraph search_graph(detector_count);
        for (const auto& [key, indices] : groups) {
            topology_edges.push_back({key, indices});
            if (key.second == BOUNDARY_NODE) {
                matching_graph.add_boundary_edge(key.first, 0, {});
                search_graph.add_boundary_edge(key.first, 0, {});
            } else {
                matching_graph.add_edge(key.first, static_cast<size_t>(key.second), 0, {});
                search_graph.add_edge(key.first, static_cast<size_t>(key.second), 0, {});
            }
        }
        matcher = Mwpm(GraphFlooder(std::move(matching_graph)), SearchFlooder(std::move(search_graph)));
        for (const auto& node : matcher.search_flooder.graph.nodes) {
            search_topology_index.emplace_back(node.neighbors.size(), SIZE_MAX);
        }
        bind_weight_handles();
        for (const auto& indices : search_topology_index) {
            if (std::any_of(indices.begin(), indices.end(), [&](size_t i) { return i >= topology_edges.size(); })) {
                throw std::logic_error("unmapped search topology edge");
            }
        }
    }

    void bind_weight_handles() {
        for (auto& edge : topology_edges) {
            const auto [u, v] = edge.key;
            auto& matching_u = matcher.flooder.graph.nodes[u];
            auto* matching_v = v == BOUNDARY_NODE ? nullptr : &matcher.flooder.graph.nodes[v];
            size_t matching_index = matching_u.index_of_neighbor(matching_v);
            edge.matching_weight = &matching_u.neighbor_weights[matching_index];
            if (matching_v != nullptr) {
                edge.matching_reverse_weight = &matching_v->neighbor_weights[matching_v->index_of_neighbor(&matching_u)];
            }
            auto& search_u = matcher.search_flooder.graph.nodes[u];
            auto* search_v = v == BOUNDARY_NODE ? nullptr : &matcher.search_flooder.graph.nodes[v];
            size_t search_index = search_u.index_of_neighbor(search_v);
            edge.search_weight = &search_u.neighbor_weights[search_index];
            const size_t topology_id = static_cast<size_t>(&edge - topology_edges.data());
            search_topology_index[u][search_index] = topology_id;
            if (search_v != nullptr) {
                const size_t reverse = search_v->index_of_neighbor(&search_u);
                edge.search_reverse_weight = &search_v->neighbor_weights[reverse];
                search_topology_index[static_cast<size_t>(v)][reverse] = topology_id;
            }
        }
    }

    void validate_syndrome(const std::vector<uint8_t>& syndrome) const {
        if (syndrome.size() != detector_count) {
            throw std::invalid_argument("adjusted_syndrome length does not match num_detectors");
        }
        if (std::any_of(syndrome.begin(), syndrome.end(), [](uint8_t bit) { return bit > 1; })) {
            throw std::invalid_argument("adjusted_syndrome entries must be zero or one");
        }
    }

    bool update_weights(const std::vector<uint32_t>& weights, bool reusable, std::vector<GraphFillRegion*>& dirty) {
        size_t changed_edges = 0;
        for (auto& edge : topology_edges) {
            size_t selected = edge.class_indices.front();
            uint32_t selected_weight = weights[selected];
            for (size_t class_index : edge.class_indices) {
                uint32_t candidate = weights[class_index];
                if (candidate < selected_weight || (candidate == selected_weight && class_index < selected)) {
                    selected = class_index;
                    selected_weight = candidate;
                }
            }
            if (selected_weight > MAX_EXACT_EDGE_WEIGHT) {
                throw std::invalid_argument("quantized edge weight exceeds MAX_EXACT_EDGE_WEIGHT");
            }
            if (reusable && edge.external_weight != selected_weight) {
                // A zero-cost plateau can carry a retained contact whose
                // endpoints belong to other matched regions. Endpoint ownership
                // alone does not invalidate that witness when the plateau opens.
                if ((edge.external_weight == 0 && selected_weight > 0) || ++changed_edges > 64) {
                    reusable = false;
                } else {
                    auto& u = matcher.flooder.graph.nodes[edge.key.first];
                    auto* v = edge.key.second == BOUNDARY_NODE ? nullptr : &matcher.flooder.graph.nodes[edge.key.second];
                    auto* a = u.region_that_arrived_top;
                    auto* b = v ? v->region_that_arrived_top : nullptr;
                    const int64_t ra = a ? a->radius.y_intercept() + u.wrapped_radius_cached : 0;
                    const int64_t rb = b ? b->radius.y_intercept() + v->wrapped_radius_cached : 0;
                    // The original v4 strict-gap rule applies only across
                    // distinct regions, and to both old and new costs.
                    const bool separated = a != b && int64_t(edge.external_weight) * 2 > ra + rb &&
                                           int64_t(selected_weight) * 2 > ra + rb;
                    if (!separated) {
                        detail::invalidate_matched_node(matcher, edge.key.first, dirty);
                        if (v) detail::invalidate_matched_node(matcher, edge.key.second, dirty);
                    }
                }
            }
            // Equal-cost parallel classes may exchange identity. Update it even
            // when all four graph-weight stores can be skipped.
            edge.selected_class = selected;
            if (options.sparse_updates && edge.external_weight == selected_weight) continue;
            const weight_int internal_weight = selected_weight * 2U;
            *edge.matching_weight = internal_weight;
            *edge.search_weight = internal_weight;
            if (edge.matching_reverse_weight != nullptr) {
                *edge.matching_reverse_weight = internal_weight;
                *edge.search_reverse_weight = internal_weight;
            }
            edge.external_weight = selected_weight;
        }
        return reusable;
    }

    void prepare_matching(bool reusable, const std::vector<uint8_t>& syndrome, std::vector<GraphFillRegion*>& dirty) {
        if (!reusable) {
            matcher.reset_for_reuse();
            return;
        }
        for (size_t i = 0; i < syndrome.size(); ++i) {
            if (syndrome[i] != previous_syndrome[i]) detail::invalidate_matched_node(matcher, i, dirty);
        }
        detail::erase_dirty_regions(matcher, dirty);
        matcher.flooder.queue.clear();
        matcher.flooder.queue.cur_time = 0;
        for (auto& node : matcher.flooder.graph.nodes) node.node_event_tracker.clear();
    }

    DecodeResult decode_quantized(
        const std::vector<uint32_t>& weights,
        const std::vector<uint8_t>& syndrome,
        double normalising_constant) {
        if (weights.size() != class_edges.size()) {
            throw std::invalid_argument("weights length does not match num_class_edges");
        }
        validate_syndrome(syndrome);
        bool reusable = options.reuse_state && mature;
        // A failed solve never leaves a reusable state. The next valid call
        // resets it before entering the matcher, even after a partial update.
        mature = false;
        if (!options.reuse_state) matcher.reset_for_reuse();
        std::vector<GraphFillRegion*> dirty;
        if (reusable) dirty.reserve(128);
        reusable = update_weights(weights, reusable, dirty);
        if (options.reuse_state) prepare_matching(reusable, syndrome, dirty);
        std::vector<uint64_t> detections;
        for (size_t i = 0; i < syndrome.size(); ++i) {
            if (syndrome[i]) {
                detections.push_back(i);
            }
        }
        process_timeline(matcher, detections);
        if (options.reuse_state) detail::extract_retained_matches(matcher, detections);
        else extract_match_edges(matcher, detections);

        DecodeResult result{std::vector<uint8_t>(class_edges.size(), 0), 0, normalising_constant};
        for (const auto& match_edge : matcher.flooder.match_edges) {
            const size_t from = match_edge.loc_from - &matcher.flooder.graph.nodes[0];
            const size_t to = match_edge.loc_to == nullptr
                                  ? SIZE_MAX
                                  : static_cast<size_t>(match_edge.loc_to - &matcher.flooder.graph.nodes[0]);
            matcher.search_flooder.iter_edges_on_shortest_path_from_middle(
                from, to, [&](const SearchGraphEdge& path_edge) {
                    uint32_t u = static_cast<uint32_t>(path_edge.detector_node - &matcher.search_flooder.graph.nodes[0]);
                    auto& topology_edge = topology_edges[search_topology_index[u][path_edge.neighbor_index]];
                    result.selected_edges[topology_edge.selected_class] ^= 1;
                    if (result.objective >
                        std::numeric_limits<int64_t>::max() - topology_edge.external_weight) {
                        throw std::overflow_error("matching objective exceeds i64");
                    }
                    result.objective += topology_edge.external_weight;
                });
        }
        if (options.reuse_state) {
            previous_syndrome = syndrome;
            mature = true;
        }
        return result;
    }
};

MatchingRegion::MatchingRegion(size_t num_detectors, std::vector<ClassEdge> class_edges, MatchingOptions options)
    : impl_(std::make_unique<Impl>(num_detectors, std::move(class_edges), options)) {
}
MatchingRegion::~MatchingRegion() = default;
MatchingRegion::MatchingRegion(MatchingRegion&&) noexcept = default;
MatchingRegion& MatchingRegion::operator=(MatchingRegion&&) noexcept = default;

size_t MatchingRegion::num_detectors() const {
    return impl_->detector_count;
}

size_t MatchingRegion::num_class_edges() const {
    return impl_->class_edges.size();
}

void MatchingRegion::reset_matching_state() {
    impl_->mature = false;
    impl_->previous_syndrome.clear();
    impl_->matcher.reset_for_reuse();
}

std::vector<uint64_t> MatchingRegion::matching_state_snapshot() const {
    if (!impl_->options.reuse_state) return {};
    if (!impl_->mature) throw std::logic_error("matching snapshot requires a successful retained solve");
    auto words = detail::frozen_matching_snapshot(impl_->matcher);
    words.push_back(impl_->previous_syndrome.size());
    for (auto bit : impl_->previous_syndrome) words.push_back(bit);
    words.push_back(impl_->topology_edges.size());
    for (const auto& edge : impl_->topology_edges) {
        words.push_back(edge.external_weight);
        words.push_back(edge.selected_class);
    }
    // Exact quantization memoization changes computation cost only. Keeping it
    // in this key would let an earlier shot affect recurrence in mixed f64/u32
    // sequences, even though every current quantized input is identical.
    return words;
}

DecodeResult MatchingRegion::reweight_f64_and_decode(
    const std::vector<double>& weights, const std::vector<uint8_t>& adjusted_syndrome) try {
    if (weights.size() != impl_->class_edges.size()) {
        throw std::invalid_argument("weights length does not match num_class_edges");
    }
    double maximum = 0;
    bool all_integral = true;
    for (double weight : weights) {
        if (!std::isfinite(weight) || weight < 0) {
            throw std::invalid_argument("f64 weights must be finite and nonnegative");
        }
        maximum = std::max(maximum, weight);
        all_integral &= std::round(weight) == weight;
    }
    const double scale = (all_integral && maximum <= MAX_EXACT_EDGE_WEIGHT) || maximum == 0
                             ? 1.0
                             : static_cast<double>(MAX_EXACT_EDGE_WEIGHT) / maximum;
    const bool scale_overflowed = !std::isfinite(scale);
    if (!impl_->options.sparse_updates) {
        std::vector<uint32_t> quantized;
        quantized.reserve(weights.size());
        for (double weight : weights) {
            const double scaled = scale_overflowed
                                      ? (weight / maximum) * static_cast<double>(MAX_EXACT_EDGE_WEIGHT)
                                      : weight * scale;
            quantized.push_back(static_cast<uint32_t>(std::round(scaled)));
        }
        return impl_->decode_quantized(quantized, adjusted_syndrome, scale);
    }
    auto& quantized = impl_->retained_quantized;
    const bool reusable_scale = impl_->previous_float_weights.size() == weights.size() &&
                                impl_->previous_scale == scale &&
                                (!scale_overflowed || impl_->previous_maximum == maximum);
    quantized.resize(weights.size());
    for (size_t i = 0; i < weights.size(); ++i) {
        const double weight = weights[i];
        // Exact computation reuse, never an epsilon-based weight approximation.
        if (reusable_scale && weight == impl_->previous_float_weights[i]) continue;
        const double scaled = scale_overflowed
                                  ? (weight / maximum) * static_cast<double>(MAX_EXACT_EDGE_WEIGHT)
                                  : weight * scale;
        quantized[i] = static_cast<uint32_t>(std::round(scaled));
    }
    impl_->previous_float_weights = weights;
    impl_->previous_scale = scale;
    impl_->previous_maximum = maximum;
    return impl_->decode_quantized(quantized, adjusted_syndrome, scale);
} catch (...) {
    impl_->mature = false;
    throw;
}

DecodeResult MatchingRegion::reweight_u32_and_decode(
    const std::vector<uint32_t>& weights, const std::vector<uint8_t>& adjusted_syndrome) try {
    for (uint32_t weight : weights) {
        if (weight > MAX_EXACT_EDGE_WEIGHT) {
            throw std::invalid_argument("u32 weight exceeds MAX_EXACT_EDGE_WEIGHT");
        }
    }
    return impl_->decode_quantized(weights, adjusted_syndrome, 1.0);
} catch (...) {
    impl_->mature = false;
    throw;
}

}  // namespace pm::core
