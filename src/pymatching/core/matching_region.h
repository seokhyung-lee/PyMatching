// Copyright 2026 PyMatching fork contributors.
// Modified PyMatching file: Stim-free reusable whole-graph reweight/decode API.
// Licensed under the Apache License, Version 2.0. See LICENSE.

#ifndef PYMATCHING_CORE_MATCHING_REGION_H
#define PYMATCHING_CORE_MATCHING_REGION_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace pm::core {

inline constexpr uint32_t MAX_EXACT_EDGE_WEIGHT = (1U << 24U) - 1U;
inline constexpr int64_t BOUNDARY_NODE = -1;

struct ClassEdge {
    uint32_t node1;
    int64_t node2;
};

struct DecodeResult {
    std::vector<uint8_t> selected_edges;
    int64_t objective;
    // May be positive infinity when the exact f64-to-integer scale is larger
    // than the f64 range. The selected edges and integer objective remain valid.
    double normalising_constant;
};

class MatchingRegion {
   public:
    MatchingRegion(size_t num_detectors, std::vector<ClassEdge> class_edges);
    ~MatchingRegion();
    MatchingRegion(MatchingRegion&&) noexcept;
    MatchingRegion& operator=(MatchingRegion&&) noexcept;
    MatchingRegion(const MatchingRegion&) = delete;
    MatchingRegion& operator=(const MatchingRegion&) = delete;

    size_t num_detectors() const;
    size_t num_class_edges() const;
    DecodeResult reweight_f64_and_decode(
        const std::vector<double>& weights, const std::vector<uint8_t>& adjusted_syndrome);
    DecodeResult reweight_u32_and_decode(
        const std::vector<uint32_t>& weights, const std::vector<uint8_t>& adjusted_syndrome);

   private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace pm::core

#endif
