// Copyright 2026 PyMatching fork contributors.
// Added file: conservative retained matching lifecycle and exact-state helpers.
// Licensed under the Apache License, Version 2.0. See LICENSE.
#ifndef PYMATCHING_CORE_MATCHING_REUSE_H
#define PYMATCHING_CORE_MATCHING_REUSE_H

#include <cstdint>
#include <vector>
#include "pymatching/sparse_blossom/matcher/mwpm.h"

namespace pm::core::detail {

void extract_retained_matches(Mwpm& matcher, const std::vector<uint64_t>& detections);
void invalidate_matched_node(Mwpm& matcher, size_t node, std::vector<GraphFillRegion*>& dirty);
void erase_dirty_regions(Mwpm& matcher, std::vector<GraphFillRegion*>& dirty);
std::vector<uint64_t> frozen_matching_snapshot(const Mwpm& matcher);

}  // namespace pm::core::detail
#endif
