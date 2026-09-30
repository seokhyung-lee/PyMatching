// Copyright 2026 PyMatching fork contributors.
// Added file: whole matched-pair invalidation and read-only blossom extraction.
// Licensed under the Apache License, Version 2.0. See LICENSE.
#include "pymatching/core/matching_reuse.h"

#include <algorithm>
#include <stdexcept>

namespace pm::core::detail {
namespace {

// Read-only counterpart to the native pair-and-shatter extraction order.
void expand_internal(GraphFillRegion* region, DetectorNode* endpoint, std::vector<CompressedEdge>& out) {
    if (region->blossom_children.empty()) return;
    auto* child = endpoint->region_that_arrived;
    while (child && child->blossom_parent != region) child = child->blossom_parent;
    if (!child) throw std::logic_error("external endpoint is not inside blossom");
    const auto& children = region->blossom_children;
    size_t base = 0;
    while (base < children.size() && children[base].region != child) ++base;
    if (base == children.size()) throw std::logic_error("missing blossom base");
    for (size_t offset = 1; offset < children.size(); offset += 2) {
        const auto& a = children[(base + offset) % children.size()];
        const auto& b = children[(base + offset + 1) % children.size()];
        expand_internal(a.region, a.edge.loc_from, out);
        expand_internal(b.region, a.edge.loc_to, out);
        out.push_back(a.edge);
    }
    expand_internal(child, endpoint, out);
}

DetectorNode* first_source(GraphFillRegion* region) {
    if (region->blossom_children.empty()) return region->shell_area.front();
    auto* result = first_source(region->blossom_children.front().region);
    for (const auto& child : region->blossom_children) result = std::min(result, first_source(child.region));
    return result;
}

void delete_hierarchy(Mwpm& matcher, GraphFillRegion* region) {
    for (const auto& child : region->blossom_children) delete_hierarchy(matcher, child.region);
    region->cleanup_shell_area();
    matcher.flooder.region_arena.del(region);
}

}  // namespace

void extract_retained_matches(Mwpm& matcher, const std::vector<uint64_t>& detections) {
    auto& out = matcher.flooder.match_edges;
    out.clear();
    std::vector<GraphFillRegion*> seen;
    seen.reserve(detections.size());
    for (auto detection : detections) {
        auto* region = matcher.flooder.graph.nodes[detection].region_that_arrived_top;
        if (!region || std::find(seen.begin(), seen.end(), region) != seen.end()) continue;
        seen.push_back(region);
        expand_internal(region, region->match.edge.loc_from, out);
        if (region->match.region) {
            seen.push_back(region->match.region);
            expand_internal(region->match.region, region->match.edge.loc_to, out);
        }
        out.push_back(region->match.edge);
    }
}

void invalidate_matched_node(Mwpm& matcher, size_t node, std::vector<GraphFillRegion*>& dirty) {
    auto* region = matcher.flooder.graph.nodes[node].region_that_arrived_top;
    if (region) {
        dirty.push_back(region);
        if (region->match.region) dirty.push_back(region->match.region);
    }
}

void erase_dirty_regions(Mwpm& matcher, std::vector<GraphFillRegion*>& dirty) {
    // All sources belong to one fixed detector vector. Its order is a logical
    // detector order, unlike sorting separately allocated region addresses.
    std::sort(dirty.begin(), dirty.end(), [](auto* a, auto* b) { return first_source(a) < first_source(b); });
    dirty.erase(std::unique(dirty.begin(), dirty.end()), dirty.end());
    for (auto* region : dirty) delete_hierarchy(matcher, region);
}

}  // namespace pm::core::detail
