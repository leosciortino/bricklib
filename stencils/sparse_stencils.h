#pragma once

#include <functional>

using SparseIndicator = std::function<bool(long, long, long)>;
using SparseStencilLauncher =
    std::function<void(unsigned,
                       const Brick<Dim<BDIM>, Dim<VFOLD>> &,
                       const Brick<Dim<BDIM>, Dim<VFOLD>> &)>;
using SparseReferenceStencil =
    std::function<bElem(Brick<Dim<BDIM>, Dim<VFOLD>> &,
                        unsigned, long, long, long)>;

void run_sparsecu(const SparseIndicator &indicator,
                  const SparseStencilLauncher &launch_stencil,
                  const SparseReferenceStencil &reference_stencil);

void launch_star_sparsecu(unsigned stencil_radius,
                          unsigned active_brick_count,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bIn,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bOut);

void launch_cube_sparsecu(unsigned stencil_radius,
                          unsigned active_brick_count,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bIn,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bOut);
