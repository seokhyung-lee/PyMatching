// Copyright 2026 PyMatching fork contributors.
// Added file: snapshot field, inactive-storage and order sensitivity tests.
// Licensed under the Apache License, Version 2.0. See LICENSE.
#include "pymatching/core/matching_reuse.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using namespace pm;
using pm::core::detail::frozen_matching_snapshot;

namespace {
Mwpm fixture() {
    MatchingGraph graph(5, 0);
    SearchGraph search(5);
    Mwpm matcher{GraphFlooder(std::move(graph)), SearchFlooder(std::move(search))};
    std::vector<GraphFillRegion*> leaves;
    for (size_t i = 0; i < 4; ++i) {
        auto* region = matcher.flooder.region_arena.alloc_default_constructed();
        auto* node = &matcher.flooder.graph.nodes[i];
        region->radius = VaryingCT::frozen(2);
        region->shell_area.push_back(node);
        node->region_that_arrived = node->region_that_arrived_top = region;
        node->reached_from_source = node;
        leaves.push_back(region);
    }
    auto& nodes = matcher.flooder.graph.nodes;
    auto* blossom = matcher.flooder.region_arena.alloc_default_constructed();
    blossom->radius = VaryingCT::frozen(1);
    for (size_t i = 0; i < 3; ++i) {
        blossom->blossom_children.push_back({leaves[i], {&nodes[i], &nodes[(i + 1) % 3], 0}});
        leaves[i]->wrap_into_blossom(blossom);
    }
    blossom->add_match(leaves[3], {&nodes[0], &nodes[3], 0});
    leaves[0]->shell_area.push_back(&nodes[4]);
    nodes[4].region_that_arrived = leaves[0];
    nodes[4].region_that_arrived_top = blossom;
    nodes[4].reached_from_source = &nodes[0];
    nodes[4].wrapped_radius_cached = 1;
    for (size_t i = 0; i < 2; ++i) {
        auto* unused = matcher.flooder.region_arena.alloc_default_constructed();
        matcher.flooder.region_arena.del(unused);
    }
    return matcher;
}

template <typename Apply, typename Restore>
void changes_key(Mwpm& matcher, Apply apply, Restore restore) {
    const auto before = frozen_matching_snapshot(matcher);
    apply();
    assert(frozen_matching_snapshot(matcher) != before);
    restore();
    assert(frozen_matching_snapshot(matcher) == before);
}
}  // namespace

int main() {
    auto matcher = fixture();
    auto other = fixture();
    assert(frozen_matching_snapshot(matcher) == frozen_matching_snapshot(other));
    auto& nodes = matcher.flooder.graph.nodes;
    auto* leaf = nodes[0].region_that_arrived;
    auto* blossom = nodes[0].region_that_arrived_top;
    changes_key(matcher, [&] { leaf->radius = VaryingCT::frozen(3); }, [&] { leaf->radius = VaryingCT::frozen(2); });
    changes_key(matcher, [&] { ++nodes[0].radius_of_arrival; }, [&] { --nodes[0].radius_of_arrival; });
    changes_key(matcher, [&] { ++nodes[0].wrapped_radius_cached; }, [&] { --nodes[0].wrapped_radius_cached; });
    changes_key(matcher, [&] { nodes[4].reached_from_source = &nodes[1]; }, [&] { nodes[4].reached_from_source = &nodes[0]; });
    changes_key(matcher, [&] { ++nodes[4].observables_crossed_from_source; }, [&] { --nodes[4].observables_crossed_from_source; });
    changes_key(matcher, [&] { std::reverse(leaf->shell_area.begin(), leaf->shell_area.end()); },
                         [&] { std::reverse(leaf->shell_area.begin(), leaf->shell_area.end()); });
    changes_key(matcher, [&] { std::swap(blossom->blossom_children[0], blossom->blossom_children[1]); },
                         [&] { std::swap(blossom->blossom_children[0], blossom->blossom_children[1]); });
    changes_key(matcher, [&] { ++blossom->match.edge.obs_mask; }, [&] { --blossom->match.edge.obs_mask; });
    auto& arena = matcher.flooder.region_arena;
    auto* unused_a = arena.alloc_default_constructed();
    auto* unused_b = arena.alloc_default_constructed();
    arena.del(unused_a);
    arena.del(unused_b);
    const auto before = frozen_matching_snapshot(matcher);
    std::reverse(arena.available.begin(), arena.available.end());
    std::reverse(arena.allocated.begin(), arena.allocated.end());
    assert(frozen_matching_snapshot(matcher) == before);
    // These fields are inactive or overwritten before the next semantic read.
    leaf->match = Match{nullptr, {nullptr, nullptr, 123}};
    leaf->shrink_event_tracker.desired_time = cyclic_time_int{123};
    leaf->shrink_event_tracker.queued_time = cyclic_time_int{456};
    nodes[0].node_event_tracker.has_queued_time = true;
    matcher.flooder.queue.cur_time = 100;
    matcher.search_flooder.graph.nodes[0].distance_from_source = 321;
    assert(frozen_matching_snapshot(matcher) == before);
    leaf->shrink_event_tracker.has_desired_time = true;
    bool rejected = false;
    try { frozen_matching_snapshot(matcher); } catch (const std::logic_error&) { rejected = true; }
    assert(rejected);
    std::cout << "snapshot field/order/lifecycle checks passed\n";
}
