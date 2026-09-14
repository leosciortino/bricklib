//
// Created by Leo Sciortino on 9/8/26 Adapted from Tuowen Zhao's code.
//

#include "brick-hip.h"
#include "brick.h"
#include "macro_coeffs.h"
#include "stencils/stencils.h"
#include "stencils.hip.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "bricksetup.h"
#include "multiarray.h"
#include "brickcompare.h"
#include "hipvfold.h"
#include "vecscatter.h"

__global__ void
d3star7_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                    Brick <Dim<BDIM>, Dim<VFOLD>> bOut) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("python/star_7.py", VSVEC, (BDIM), (VFOLD), b);
}

__global__ void
d3star13_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                     Brick <Dim<BDIM>, Dim<VFOLD>> bOut) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("python/star_13.py", VSVEC, (BDIM), (VFOLD), b);
}

__global__ void
d3star19_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                     Brick <Dim<BDIM>, Dim<VFOLD>> bOut) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("python/star_19.py", VSVEC, (BDIM), (VFOLD), b);
}

__global__ void
d3star25_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                     Brick <Dim<BDIM>, Dim<VFOLD>> bOut) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("python/star_25.py", VSVEC, (BDIM), (VFOLD), b);
}

__global__ void
d3cube27_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                     Brick <Dim<BDIM>, Dim<VFOLD>> bOut) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("python/cube_27.py", VSVEC, (BDIM), (VFOLD), b);
}

__global__ void
d3cube125_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                      Brick <Dim<BDIM>, Dim<VFOLD>> bOut) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("python/cube_125.py", VSVEC, (BDIM), (VFOLD), b);
}

#undef bIn
#undef bOut

void launch_star_sparsecu(unsigned stencil_radius,
                          unsigned active_brick_count,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bIn,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bOut) {
  dim3 block(active_brick_count), thread(64);
  switch (stencil_radius) {
    case 1:
      d3star7_brick_trans<<<block, thread>>>(active_brick_count, bIn, bOut);
      break;
    case 2:
      d3star13_brick_trans<<<block, thread>>>(active_brick_count, bIn, bOut);
      break;
    case 3:
      d3star19_brick_trans<<<block, thread>>>(active_brick_count, bIn, bOut);
      break;
    case 4:
      d3star25_brick_trans<<<block, thread>>>(active_brick_count, bIn, bOut);
      break;
    default:
      throw std::invalid_argument("star stencil radius must be between 1 and 4");
  }
}

void launch_cube_sparsecu(unsigned stencil_radius,
                          unsigned active_brick_count,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bIn,
                          const Brick<Dim<BDIM>, Dim<VFOLD>> &bOut) {
  dim3 block(active_brick_count), thread(64);
  switch (stencil_radius) {
    case 1:
      d3cube27_brick_trans<<<block, thread>>>(active_brick_count, bIn, bOut);
      break;
    case 2:
      d3cube125_brick_trans<<<block, thread>>>(active_brick_count, bIn, bOut);
      break;
    default:
      throw std::invalid_argument("cube stencil radius must be 1 or 2");
  }
}

void run_sparsecu(const SparseIndicator &indicator,
                  const SparseStencilLauncher &launch_stencil,
                  const SparseReferenceStencil &reference_stencil) {
  if (!indicator)
    throw std::invalid_argument("sparse indicator must be callable");
  if (!launch_stencil)
    throw std::invalid_argument("sparse stencil launcher must be callable");
  if (!reference_stencil)
    throw std::invalid_argument("sparse reference stencil must be callable");
  // Maps x, y, z brick indices to a flattened brick index.
  auto grid_index = [](int x, int y, int z) {
    return (z * STRIDEB + y) * STRIDEB + x;
  };
  std::vector<char> active(STRIDEB * STRIDEB * STRIDEB, 0);
  for (int i = 0; i < N; i++)
    for (int j = 0 ; j < N; j++)
      for (int k = 0; k < N; k++)
        // Mark a brick active when the indicator selects any point in it.
        if (indicator(i, j, k))
          active[grid_index(GB + i / TILE, GB + j / TILE, GB + k / TILE)] = 1;
   
  // sparse_id: coordinate -> brickID
  // coordinate: brickid -> coordinate

  struct index {
    int x;
    int y;
    int z;
  };

  // maps grid_index(x, y, z) ie flattened logical coord ->  sparse brick coordinate \in 0..active_bricks
  std::vector<int> logical_coord_to_sparse_id(STRIDEB * STRIDEB * STRIDEB, -1);

  // maps brick coordinate-> x, y, z
  std::vector<index> sparse_id_to_logical_coord;

  for (int i = GB; i < GB + NB; i++)
    for (int j = GB; j < GB + NB; j++)
      for (int k = GB; k < GB + NB; k++) {
        if (active[grid_index(i, j, k)]){
          // maps from brick coordinate x, y, z to the brick coordinate \in 0..active bricks
          logical_coord_to_sparse_id[grid_index(i, j, k)] = sparse_id_to_logical_coord.size();
          sparse_id_to_logical_coord.push_back({i, j, k});
        }
    }

  unsigned active_brick_count = sparse_id_to_logical_coord.size();
  
  // Create BrickInfo for active_brick_count number of bricks
  BrickInfo<3> bInfo(active_brick_count);

  // Now we iterate over all of the sparse ids and add surrounding bricks
  for (unsigned i = 0; i < active_brick_count; i++) {
    unsigned n = 0;
    for (int dz = -1; dz <= 1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          int x = sparse_id_to_logical_coord[i].x;
          int y = sparse_id_to_logical_coord[i].y;
          int z = sparse_id_to_logical_coord[i].z;
          int b_neighbor_logical_coords = grid_index(x + dx, y + dy, z + dz);
          if (logical_coord_to_sparse_id[b_neighbor_logical_coords] == -1) {
            logical_coord_to_sparse_id[b_neighbor_logical_coords] = sparse_id_to_logical_coord.size();
            sparse_id_to_logical_coord.push_back({x + dx, y + dy, z + dz});
          }
          unsigned neighbor_id = logical_coord_to_sparse_id[b_neighbor_logical_coords];
          if (neighbor_id == -1) {
            throw std::runtime_error("missing neighboring brick id");
          }
          // add neighbors to brick info
          bInfo.adj[i][n++] = neighbor_id;
        }
  }
  
  // Total number of bricks which includes "active bricks" and boundry "zero bricks"
  long num_total_bricks = sparse_id_to_logical_coord.size();
  

  bElem *in_ptr = randomArray({active_brick_count, TILE, TILE, TILE});
  // moveBrickInfo
  BrickInfo<3> *bInfo_dev;
  BrickInfo<3> _bInfo_dev = movBrickInfo(bInfo, hipMemcpyHostToDevice);
  {
    unsigned size = sizeof(BrickInfo < 3 > );
    hipMalloc(&bInfo_dev, size);
    hipMemcpy(bInfo_dev, &_bInfo_dev, size, hipMemcpyHostToDevice);
  }

  // Create BrickStorage
  auto bSize  = cal_size<BDIM>::value;
  // allocate enought for bsize by two so we have an in and out array
  auto bStorage = BrickStorage::allocate(num_total_bricks, bSize * 2);

  Brick <Dim<BDIM>, Dim<VFOLD>> bIn(&bInfo, bStorage, 0);
  Brick <Dim<BDIM>, Dim<VFOLD>> bOut(&bInfo, bStorage, bSize);

  std::fill_n(
    bStorage.dat.get(),
    bStorage.chunks * bStorage.step,
    0.0);

  for (unsigned b = 0; b < active_brick_count; ++b)
    for (long k = 0; k < TILE; ++k)
      for (long j = 0; j < TILE; ++j)
        for (long i = 0; i < TILE; ++i) {
          long pos =
              ((b * TILE + k) * TILE + j) * TILE + i;

          bIn[b][k][j][i] = in_ptr[pos];
        }
  
  
  // Expected output
  std::vector<bElem> expected(active_brick_count * bSize);
  for (unsigned b = 0; b < active_brick_count; ++b)
    for (long k = 0; k < TILE; ++k)
      for (long j = 0; j < TILE; ++j)
        for (long i = 0; i < TILE; ++i) {
          unsigned pos = ((b * TILE + k) * TILE + j) * TILE + i;
          expected[pos] = reference_stencil(bIn, b, k, j, i);
        }

  // movBrickStorage
  BrickStorage bStorage_dev = movBrickStorage(bStorage, hipMemcpyHostToDevice);
  Brick<Dim<BDIM>, Dim<VFOLD>> bIn_dev(bInfo_dev, bStorage_dev, 0);
  Brick<Dim<BDIM>, Dim<VFOLD>> bOut_dev(bInfo_dev, bStorage_dev, bSize);

  auto brick_func_trans = [&]() -> void {
    launch_stencil(active_brick_count, bIn_dev, bOut_dev);
  };

  std::cout << "Trans: " << hiptime_func(brick_func_trans) << std::endl;

  hipMemcpy(bStorage.dat.get(), bStorage_dev.dat.get(),
            bStorage.chunks * bStorage.step * sizeof(bElem),
            hipMemcpyDeviceToHost);
  hipDeviceSynchronize();

  bool matches = true;
  for (unsigned b = 0; b < active_brick_count; ++b)
    for (long k = 0; k < TILE; ++k)
      for (long j = 0; j < TILE; ++j)
        for (long i = 0; i < TILE; ++i) {
          unsigned pos = ((b * TILE + k) * TILE + j) * TILE + i;
          bElem actual = bOut[b][k][j][i];
          bElem reference = expected[pos];
          bElem diff = std::abs(actual - reference);
          bool close = diff < BRICK_TOLERANCE ||
                       diff < (std::abs(actual) + std::abs(reference)) * BRICK_TOLERANCE;
          matches = matches && close;
        }

  free(in_ptr);
  free(bInfo.adj);
  hipFree(_bInfo_dev.adj);
  hipFree(bInfo_dev);
  if (!matches)
    throw std::runtime_error("sparse brick result mismatch");

 
}
