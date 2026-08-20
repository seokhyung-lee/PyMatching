// Copyright 2026 PyMatching fork contributors.
// Modified PyMatching file: dependency-free native contract tests.
// Licensed under the Apache License, Version 2.0. See LICENSE.

#include "pymatching/core/matching_region.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using pm::core::BOUNDARY_NODE;
using pm::core::ClassEdge;
using pm::core::DecodeResult;
using pm::core::MatchingRegion;
using pm::core::MAX_EXACT_EDGE_WEIGHT;

namespace {

const std::vector<ClassEdge> EDGES = {
    {0, 1}, {1, 2}, {2, 3}, {0, 3}, {0, BOUNDARY_NODE}, {3, BOUNDARY_NODE}, {1, 0}};

void assert_same(const DecodeResult& actual, const DecodeResult& expected) {
    assert(actual.selected_edges == expected.selected_edges);
    assert(actual.objective == expected.objective);
    assert(actual.normalising_constant == expected.normalising_constant);
}

DecodeResult fresh_u32(const std::vector<uint32_t>& weights, const std::vector<uint8_t>& syndrome) {
    MatchingRegion fresh(4, EDGES);
    return fresh.reweight_u32_and_decode(weights, syndrome);
}

DecodeResult fresh_f64(const std::vector<double>& weights, const std::vector<uint8_t>& syndrome) {
    MatchingRegion fresh(4, EDGES);
    return fresh.reweight_f64_and_decode(weights, syndrome);
}

void test_fresh_rebuild_equivalence() {
    MatchingRegion reused(4, EDGES);
    const std::vector<std::vector<uint32_t>> weights = {
        {9, 1, 1, 9, 2, 8, 3},
        {1, 9, 9, 1, 8, 2, 4},
        {MAX_EXACT_EDGE_WEIGHT, 0, MAX_EXACT_EDGE_WEIGHT, 0, 0, MAX_EXACT_EDGE_WEIGHT, 0},
        {7, 4, 6, 2, 9, 5, 1},
    };
    const std::vector<std::vector<uint8_t>> syndromes = {
        {1, 0, 0, 0}, {1, 1, 0, 0}, {1, 0, 1, 0}, {0, 0, 0, 1}};
    for (size_t i = 0; i < weights.size(); ++i) {
        assert_same(reused.reweight_u32_and_decode(weights[i], syndromes[i]), fresh_u32(weights[i], syndromes[i]));
    }

    const std::vector<double> f64_a = {0.3, 1.1, 1.2, 4.0, 2.0, 7.0, 0.2};
    const std::vector<double> f64_b = {8.0, 0.25, 0.25, 1.0, 9.0, 3.0, 2.0};
    assert_same(reused.reweight_f64_and_decode(f64_a, {1, 1, 0, 0}), fresh_f64(f64_a, {1, 1, 0, 0}));
    assert_same(reused.reweight_f64_and_decode(f64_b, {1, 0, 1, 0}), fresh_f64(f64_b, {1, 0, 1, 0}));
}

void test_parallel_boundary_and_objective() {
    MatchingRegion region(4, EDGES);
    auto parallel = region.reweight_u32_and_decode({9, 50, 50, 50, 50, 50, 2}, {1, 1, 0, 0});
    assert(parallel.objective == 2);
    assert(parallel.selected_edges[6] == 1);
    assert(parallel.selected_edges[0] == 0);

    auto boundary = region.reweight_u32_and_decode({20, 20, 20, 20, 3, 20, 20}, {1, 0, 0, 0});
    assert(boundary.objective == 3);
    assert(boundary.selected_edges[4] == 1);

    auto integral_f64 = region.reweight_f64_and_decode({9, 50, 50, 50, 50, 50, 2}, {1, 1, 0, 0});
    assert(integral_f64.objective == 2);
    assert(integral_f64.normalising_constant == 1.0);

    MatchingRegion maximum_region(1, {{0, BOUNDARY_NODE}});
    auto maximum = maximum_region.reweight_u32_and_decode({MAX_EXACT_EDGE_WEIGHT}, {1});
    assert(maximum.objective == MAX_EXACT_EDGE_WEIGHT);
    assert(maximum.selected_edges == std::vector<uint8_t>{1});

    auto nonintegral = maximum_region.reweight_f64_and_decode({0.5}, {1});
    assert(nonintegral.objective == MAX_EXACT_EDGE_WEIGHT);
    assert(nonintegral.normalising_constant == 2.0 * MAX_EXACT_EDGE_WEIGHT);

    auto large_integral = maximum_region.reweight_f64_and_decode(
        {2.0 * static_cast<double>(MAX_EXACT_EDGE_WEIGHT)}, {1});
    assert(large_integral.objective == MAX_EXACT_EDGE_WEIGHT);
    assert(large_integral.normalising_constant == 0.5);

    auto huge_finite = maximum_region.reweight_f64_and_decode({1e100}, {1});
    assert(huge_finite.objective == MAX_EXACT_EDGE_WEIGHT);
    assert(std::isfinite(huge_finite.normalising_constant));
    assert(huge_finite.normalising_constant > 0);
    assert(
        std::abs(
            huge_finite.normalising_constant -
            static_cast<double>(MAX_EXACT_EDGE_WEIGHT) / 1e100) <=
        std::numeric_limits<double>::epsilon() * huge_finite.normalising_constant);

    auto tiny_finite = maximum_region.reweight_f64_and_decode({1e-307}, {1});
    assert(tiny_finite.objective == MAX_EXACT_EDGE_WEIGHT);
    assert(tiny_finite.normalising_constant == std::numeric_limits<double>::infinity());
}

void test_changing_shortest_paths_and_repetition() {
    MatchingRegion region(4, EDGES);
    const std::vector<uint8_t> syndrome = {1, 0, 1, 0};
    auto via_one = region.reweight_u32_and_decode({1, 1, 20, 20, 50, 50, 30}, syndrome);
    assert(via_one.objective == 2);
    assert(via_one.selected_edges[0] == 1 && via_one.selected_edges[1] == 1);

    auto via_three = region.reweight_u32_and_decode({20, 20, 1, 1, 50, 50, 30}, syndrome);
    assert(via_three.objective == 2);
    assert(via_three.selected_edges[2] == 1 && via_three.selected_edges[3] == 1);

    auto repeated = region.reweight_u32_and_decode({20, 20, 1, 1, 50, 50, 30}, syndrome);
    assert_same(repeated, via_three);
}

void test_exception_recovery() {
    MatchingRegion no_boundary(3, {{0, 1}, {1, 2}, {0, 2}});
    bool failed = false;
    try {
        (void)no_boundary.reweight_u32_and_decode({1, 1, 1}, {1, 0, 0});
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    assert(failed);
    auto recovered = no_boundary.reweight_u32_and_decode({2, 1, 9}, {1, 0, 1});
    assert(recovered.objective == 3);

    failed = false;
    try {
        (void)no_boundary.reweight_f64_and_decode(
            {1.0, std::numeric_limits<double>::quiet_NaN(), 2.0}, {1, 0, 1});
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    assert(failed);
    assert_same(
        no_boundary.reweight_u32_and_decode({2, 1, 9}, {1, 0, 1}),
        MatchingRegion(3, {{0, 1}, {1, 2}, {0, 2}}).reweight_u32_and_decode({2, 1, 9}, {1, 0, 1}));
}

}  // namespace

int main() {
    test_fresh_rebuild_equivalence();
    test_parallel_boundary_and_objective();
    test_changing_shortest_paths_and_repetition();
    test_exception_recovery();
    std::cout << "pymatching_core_tests: all assertions passed\n";
}
