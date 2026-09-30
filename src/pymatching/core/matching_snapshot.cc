// Copyright 2026 PyMatching fork contributors.
// Added file: exact comparison of frozen retained matching state.
// Licensed under the Apache License, Version 2.0. See LICENSE.
#include "pymatching/core/matching_reuse.h"

#include <stdexcept>
#include <unordered_map>

namespace pm::core::detail {
namespace {

// Index zero denotes null. Live-region IDs follow detector-ordered roots and
// ordered-child preorder, independent of retained physical arena capacity.
class SnapshotWriter {
    const Mwpm& matcher;
    std::unordered_map<const GraphFillRegion*, uint64_t> region_ids;
    std::vector<const GraphFillRegion*> ordered_regions;

    void identify(const GraphFillRegion* region) {
        if (region_ids.contains(region)) return;
        region_ids.emplace(region, ordered_regions.size() + 1);
        ordered_regions.push_back(region);
        for (const auto& child : region->blossom_children) identify(child.region);
    }

    uint64_t node_id(const DetectorNode* node) const {
        return node ? static_cast<uint64_t>(node - matcher.flooder.graph.nodes.data()) + 1 : 0;
    }

    uint64_t region_id(const GraphFillRegion* region) const {
        return region ? region_ids.at(region) : 0;
    }

    void edge(const CompressedEdge& value) {
        words.insert(words.end(), {node_id(value.loc_from), node_id(value.loc_to), value.obs_mask});
    }

    void region(const GraphFillRegion& value) {
        if (!value.radius.is_frozen() || value.alt_tree_node ||
            value.shrink_event_tracker.has_desired_time || value.shrink_event_tracker.has_queued_time) {
            throw std::logic_error("matching snapshot requires frozen idle regions");
        }
        words.insert(words.end(), {region_id(value.blossom_parent), region_id(value.blossom_parent_top),
                                 static_cast<uint64_t>(value.radius.y_intercept())});
        // A contracted child's inactive match may be stale or uninitialized.
        // Its internal pairing is represented by the ordered blossom cycle.
        if (!value.blossom_parent) {
            words.push_back(region_id(value.match.region));
            edge(value.match.edge);
        }
        words.push_back(value.blossom_children.size());
        for (const auto& child : value.blossom_children) {
            words.push_back(region_id(child.region));
            edge(child.edge);
        }
        words.push_back(value.shell_area.size());
        for (const auto* node : value.shell_area) words.push_back(node_id(node));
    }

    void regions() {
        const auto& arena = matcher.flooder.region_arena;
        if (ordered_regions.size() != arena.allocated.size() - arena.available.size()) {
            throw std::logic_error("matching snapshot has unreachable live regions");
        }
        words.push_back(ordered_regions.size());
        for (const auto* value : ordered_regions) region(*value);
        const auto& trees = matcher.node_arena;
        if (trees.allocated.size() != trees.available.size()) {
            throw std::logic_error("matching snapshot contains an unfinished alternating tree");
        }
        // Source-ordered invalidation removes address-order tie breaking. Arena
        // slot/free-list identities and spare capacity then affect allocation
        // cost only, not future matching decisions. Including them would make
        // recurrence depend on how much memory an earlier physical shot used.
    }

   public:
    std::vector<uint64_t> words;

    explicit SnapshotWriter(const Mwpm& value) : matcher(value) {
        for (const auto& node : matcher.flooder.graph.nodes) {
            if (node.reached_from_source == &node && node.region_that_arrived_top) identify(node.region_that_arrived_top);
        }
    }

    std::vector<uint64_t> write() {
        if (!matcher.flooder.queue.empty() || !matcher.search_flooder.queue.empty() ||
            !matcher.search_flooder.reached_nodes.empty()) {
            throw std::logic_error("matching snapshot requires drained matching and search queues");
        }
        regions();
        words.push_back(matcher.flooder.graph.nodes.size());
        for (const auto& node : matcher.flooder.graph.nodes) {
            words.insert(words.end(), {region_id(node.region_that_arrived), region_id(node.region_that_arrived_top),
                node_id(node.reached_from_source), static_cast<uint64_t>(node.radius_of_arrival),
                static_cast<uint64_t>(node.wrapped_radius_cached), node.observables_crossed_from_source});
        }
        // Immutable adjacency is fixed by this MatchingRegion. Matching queue
        // time/node trackers are reset before reuse. Inactive shrink timestamps,
        // extraction scratch and reset search distances/targets are never read
        // before being overwritten, so their stale storage is not semantic state.
        return std::move(words);
    }
};

}  // namespace

std::vector<uint64_t> frozen_matching_snapshot(const Mwpm& matcher) {
    return SnapshotWriter(matcher).write();
}

}  // namespace pm::core::detail
