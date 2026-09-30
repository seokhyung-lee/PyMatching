// Copyright 2026 PyMatching fork contributors.
// Added file: exact reuse, warm sequence and error-lifecycle contract tests.
// Licensed under the Apache License, Version 2.0. See LICENSE.
#include "pymatching/core/matching_region.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

using namespace pm::core;

namespace {
const std::vector<ClassEdge> EDGES{{0, 1}, {1, 2}, {2, 3}, {0, 3}, {0, -1}, {3, -1}, {1, 0}};

void verify(const DecodeResult& result, const DecodeResult& reference, size_t n,
            const std::vector<ClassEdge>& edges, const std::vector<uint8_t>& syndrome, bool strict) {
    if (result.objective != reference.objective)
        std::cerr << "objective actual=" << result.objective << " expected=" << reference.objective << '\n';
    assert(result.objective == reference.objective);
    assert(result.normalising_constant == reference.normalising_constant ||
           std::abs(result.normalising_constant - reference.normalising_constant) <= 1e-10);
    if (strict) assert(result.selected_edges == reference.selected_edges);
    assert(result.selected_edges.size() == edges.size());
    std::vector<uint8_t> actual(n, 0);
    for (size_t i = 0; i < edges.size(); ++i) {
        assert(result.selected_edges[i] <= 1);
        actual[edges[i].node1] ^= result.selected_edges[i];
        if (edges[i].node2 != BOUNDARY_NODE) actual[edges[i].node2] ^= result.selected_edges[i];
    }
    assert(actual == syndrome);
}

void verify_cost(const DecodeResult& result, const std::vector<uint32_t>& weights) {
    int64_t cost = 0;
    for (size_t i = 0; i < weights.size(); ++i) cost += int64_t(weights[i]) * result.selected_edges[i];
    assert(result.objective == cost);
}

void exact_tiny(MatchingOptions options) {
    const std::vector<ClassEdge> edges{{0, 1}, {1, 2}, {2, 3}, {0, 3}, {0, -1}, {3, -1}, {0, 2}};
    MatchingRegion candidate(4, edges, options);
    for (uint32_t iteration = 0; iteration < 96; ++iteration) {
        std::vector<uint32_t> weights(edges.size());
        for (size_t e = 0; e < edges.size(); ++e)
            weights[e] = iteration % 8 == 0 ? (e % 2 ? 0 : MAX_EXACT_EDGE_WEIGHT) : (e * 7 + iteration * 13) % 11;
        for (uint32_t code = 0; code < 16; ++code) {
            std::vector<uint8_t> syndrome(4);
            for (size_t v = 0; v < 4; ++v) syndrome[v] = (code >> v) & 1;
            int64_t optimum = INT64_MAX;
            for (uint32_t mask = 0; mask < (1U << edges.size()); ++mask) {
                int64_t cost = 0;
                std::vector<uint8_t> parity(4, 0);
                for (size_t e = 0; e < edges.size(); ++e) if ((mask >> e) & 1) {
                    cost += weights[e];
                    parity[edges[e].node1] ^= 1;
                    if (edges[e].node2 != -1) parity[edges[e].node2] ^= 1;
                }
                if (parity == syndrome) optimum = std::min(optimum, cost);
            }
            const auto actual = candidate.reweight_u32_and_decode(weights, syndrome);
            verify(actual, {actual.selected_edges, optimum, 1}, 4, edges, syndrome, false);
            verify_cost(actual, weights);
            if (options.reuse_state) assert(!candidate.matching_state_snapshot().empty());
        }
    }
}

void mixed_updates(MatchingOptions options, uint32_t seed, size_t graphs, bool stress = false) {
    std::mt19937 random(seed);
    for (size_t graph = 0; graph < graphs; ++graph) {
        size_t n = 8 + graph % (stress ? 89 : 25);
        std::vector<ClassEdge> edges;
        for (size_t i = 0; i < n; ++i) edges.push_back({uint32_t(i), int64_t((i + 1) % n)});
        for (size_t i = 0; i < n; ++i) {
            uint32_t u = random() % n, v = random() % n;
            if (u != v) edges.push_back({u, v});
        }
        if (graph % 2) edges.push_back({0, -1});
        MatchingRegion cold(n, edges), candidate(n, edges, options);
        std::vector<uint32_t> weights(edges.size());
        std::vector<uint8_t> syndrome(n);
        for (auto& w : weights) w = random() % 101;
        for (size_t iteration = 0; iteration < 500; ++iteration) {
            const auto previous_weights = weights;
            const auto previous_syndrome = syndrome;
            for (size_t j = 0; j < 1 + iteration % 4; ++j)
                weights[random() % weights.size()] = iteration % 97 == 0 ? MAX_EXACT_EDGE_WEIGHT : random() % 101;
            if (iteration % 3) {
                const auto& edge = edges[random() % edges.size()];
                syndrome[edge.node1] ^= 1;
                if (edge.node2 != -1) syndrome[edge.node2] ^= 1;
            }
            const auto actual = candidate.reweight_u32_and_decode(weights, syndrome);
            const auto expected = cold.reweight_u32_and_decode(weights, syndrome);
            if (actual.objective != expected.objective) {
                std::cerr << "seed=" << seed << " graph=" << graph << " iteration=" << iteration << '\n';
                for (size_t e = 0; e < edges.size(); ++e)
                    std::cerr << "edge " << edges[e].node1 << ' ' << edges[e].node2 << ' '
                              << previous_weights[e] << " -> " << weights[e] << '\n';
                std::cerr << "syndrome ";
                for (size_t v = 0; v < n; ++v) std::cerr << int(previous_syndrome[v]) << int(syndrome[v]) << ' ';
                std::cerr << '\n';
            }
            verify(actual, expected, n, edges, syndrome, !options.reuse_state);
            verify_cost(actual, weights);
            if (options.reuse_state) assert(!candidate.matching_state_snapshot().empty());
        }
    }
}

void floating_interleaving(MatchingOptions options) {
    MatchingRegion reference(4, EDGES), candidate(4, EDGES, options);
    std::vector<double> base{0.3, 1.1, 1.2, 4.0, 2.0, 7.0, 0.2};
    for (size_t step = 0; step < 256; ++step) {
        base[step % base.size()] = double((step * 13) % 97) / 17.0;
        const double scales[]{1.0, 1.0, 1.0, 1e100, 1e-307, 2e-307, 0.0};
        auto weights = base;
        for (auto& weight : weights) weight *= scales[(step / 8) % 7];
        std::vector<uint8_t> syndrome(4);
        for (size_t i = 0; i < 4; ++i) syndrome[i] = (step >> i) & 1;
        auto check = [&] {
            verify(candidate.reweight_f64_and_decode(weights, syndrome),
                   reference.reweight_f64_and_decode(weights, syndrome), 4, EDGES, syndrome, !options.reuse_state);
        };
        check();
        if (options.reuse_state) {
            const auto snapshot = candidate.matching_state_snapshot();
            check();
            assert(candidate.matching_state_snapshot() == snapshot);
        }
        const std::vector<uint32_t> integers{1, 1, 0, 0, 3, 2, 1};
        verify(candidate.reweight_u32_and_decode(integers, syndrome),
               reference.reweight_u32_and_decode(integers, syndrome), 4, EDGES, syndrome, !options.reuse_state);
        check();
    }
    // Equal infinite scales do not imply equal normalizations.
    for (double maximum : {1e-307, 2e-307, 1e-307}) {
        std::vector<double> weights{1e-308, maximum, 1e-308, maximum, maximum, maximum, maximum};
        verify(candidate.reweight_f64_and_decode(weights, {1, 1, 0, 0}),
               reference.reweight_f64_and_decode(weights, {1, 1, 0, 0}), 4, EDGES, {1, 1, 0, 0}, !options.reuse_state);
    }
}

void independent_sequences_and_recovery(MatchingOptions options) {
    MatchingRegion reused(4, EDGES, options);
    const std::vector<uint32_t> a{2, 1, 3, 1, 8, 4, 2}, b{0, 4, 0, 3, 1, 2, 0};
    for (size_t sequence = 0; sequence < 12; ++sequence) {
        MatchingRegion fresh(4, EDGES, options);
        reused.reset_matching_state();
        for (const auto& weights : {a, b, a}) {
            const std::vector<uint8_t> syndrome{1, 0, uint8_t(sequence % 2), 0};
            verify(reused.reweight_u32_and_decode(weights, syndrome), fresh.reweight_u32_and_decode(weights, syndrome),
                   4, EDGES, syndrome, true);
        }
        bool rejected = false;
        try { reused.reweight_f64_and_decode({1, 2, 3, 4, 5, 6, 7}, {2, 0, 0, 0}); }
        catch (const std::invalid_argument&) { rejected = true; }
        assert(rejected);
        if (options.reuse_state) {
            bool snapshot_rejected = false;
            try { reused.matching_state_snapshot(); }
            catch (const std::logic_error&) { snapshot_rejected = true; }
            assert(snapshot_rejected);
        }
        MatchingRegion recovery(4, EDGES, options);
        verify(reused.reweight_u32_and_decode(a, {1, 0, 1, 0}), recovery.reweight_u32_and_decode(a, {1, 0, 1, 0}),
               4, EDGES, {1, 0, 1, 0}, true);
    }
    MatchingRegion impossible(2, {}, options);
    bool rejected = false;
    try { impossible.reweight_u32_and_decode({}, {1, 0}); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    assert(impossible.reweight_u32_and_decode({}, {0, 0}).objective == 0);

    // A widespread update takes the original >64-topology-edge cold path.
    std::vector<ClassEdge> many;
    for (uint32_t i = 0; i < 80; ++i) many.push_back({i, -1});
    MatchingRegion many_reused(80, many, options), many_cold(80, many);
    for (uint32_t cost : {1, 2, 0, 3}) {
        std::vector<uint32_t> weights(80, cost);
        std::vector<uint8_t> syndrome(80, 1);
        verify(many_reused.reweight_u32_and_decode(weights, syndrome), many_cold.reweight_u32_and_decode(weights, syndrome),
               80, many, syndrome, true);
    }
}

void snapshot_equality_is_shot_independent() {
  for (bool sparse : {false, true}) {
    MatchingRegion fresh(1, {{0, -1}}, {sparse, true}), prewarmed(1, {{0, -1}}, {sparse, true});
    prewarmed.reweight_f64_and_decode({1}, {1});
    prewarmed.reset_matching_state();
    std::vector<std::vector<uint64_t>> a, b;
    for (uint8_t syndrome : {0, 1, 0, 1, 0}) {
        if (syndrome) {
            fresh.reweight_f64_and_decode({1}, {syndrome});
            prewarmed.reweight_f64_and_decode({1}, {syndrome});
        } else {
            fresh.reweight_u32_and_decode({1}, {syndrome});
            prewarmed.reweight_u32_and_decode({1}, {syndrome});
        }
        a.push_back(fresh.matching_state_snapshot());
        b.push_back(prewarmed.matching_state_snapshot());
    }
    for (size_t i = 0; i < a.size(); ++i) for (size_t j = 0; j < a.size(); ++j)
        assert((a[i] == a[j]) == (b[i] == b[j]));
  }
}
}  // namespace

int main() {
    snapshot_equality_is_shot_independent();
    for (bool sparse : {false, true}) for (bool retained : {false, true}) {
        MatchingOptions options{sparse, retained};
        std::cout << "testing sparse=" << sparse << " reuse=" << retained << " tiny\n" << std::flush;
        exact_tiny(options);
        std::cout << "floating\n" << std::flush;
        floating_interleaving(options);
        std::cout << "lifecycle\n" << std::flush;
        independent_sequences_and_recovery(options);
        std::cout << "mixed\n" << std::flush;
        // Includes all earlier failed sequence coordinates and long histories.
        mixed_updates(options, 20260929, 40);
        mixed_updates(options, 20260930, 200);
        mixed_updates(options, 20261002, 100, true);
        std::cout << "sparse=" << sparse << " reuse=" << retained << " tiny=1536 mixed=170000 all_verified=1\n";
    }
}
