//
// Created by Tuowen Zhao on 12/5/18.
//

#include "brick-cuda.h"
#include "brick.h"
#include "stencils/stencils.h"
#include "stencils_cu.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "bricksetup.h"
#include "multiarray.h"
#include "brickcompare.h"
#include "cudavfold.h"
#include "vecscatter.h"

__global__ void
d3pt7_brick_trans(unsigned active_count, Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                  Brick <Dim<BDIM>, Dim<VFOLD>> bOut,
                  bElem *coeff) {
  unsigned b = blockIdx.x;
  if (b >= active_count)
    return;
  brick("7pt.py", VSVEC, (BDIM), (VFOLD), b);
}

#undef bIn
#undef bOut

// This function comparse the sparse implemntation to the cond implementaiton
void d3pt7_sparsecu(double line_width) {

  if (line_width < 0.0)
    throw std::invalid_argument("line width must be non-negative");

    // Maps x, y, z, brick indices to the flattened brick indices
  auto grid_index = [](int x, int y, int z) {
    return (z * STRIDEB + y) * STRIDEB + x;
  };
  // returns true if there the (x, y, z) is less than line_width away from the line x = y = z
  // We use the equation (x−y)^2+(y−z)^2+(z−x)^2 ≤ 3d^2
  auto in_line = [line_width](long x, long y, long z) {
    return (x-y)*(x-y) + (y-z)*(y-z) + (z-x)*(z-x) <= 3 * line_width * line_width;
  };

  std::vector<char> active(STRIDEB * STRIDEB * STRIDEB, 0);
  for (int i = 0; i < N; i++)
    for (int j = 0 ; j < N; j++)
      for (int k = 0; k < N; k++)
        // If our x, y, z are within the specificed line disance, set the corresponding brick to active
        if (in_line(i, j, k))
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

  
  // Need a way given a f: brickid -> neighbor
  // so we need a way to go from brickid -> coor
  // then coord -> neighbors
  // then get neighbors brick ids 
 

  
  

  // Total number of bricks which includes "active bricks" and boundry "zero bricks"
  long num_total_bricks = sparse_id_to_logical_coord.size();
  // long num_zero_boundary_bricks = num_total_bricks-active_brick_count;
  
  // unsigned size = TILE * TILE * TILE * active_brick_count * sizeof(bElem);
  bElem *in_ptr = randomArray({active_brick_count, TILE, TILE, TILE});
  // bElem *out_ptr = zeroArray({active_brick_count, TILE, TILE, TILE});
  // bElem *zero_boundary_vals = zeroArray({num_zero_boundary_bricks, TILE, TILE, TILE});
  bElem *coeff_dev;
  {
    unsigned size = 129 * sizeof(bElem);
    cudaMalloc(&coeff_dev, size);
    cudaMemcpy(coeff_dev, coeff, size, cudaMemcpyHostToDevice);
  }


  // // Allocate in and out device arrays
  // bElem *in_dev, *out_dev;
  // {
  //   cudaMalloc(&in_dev, size);
  //   cudaMemcpy(in_dev, in_ptr, size, cudaMemcpyHostToDevice);
  // }
  // {
  //   cudaMalloc(&out_dev, size);
  //   cudaMemcpy(out_dev, out_ptr, size, cudaMemcpyHostToDevice);
  // }

  // moveBrickInfo
  BrickInfo<3> *bInfo_dev;
  BrickInfo<3> _bInfo_dev = movBrickInfo(bInfo, cudaMemcpyHostToDevice);
  {
    unsigned size = sizeof(BrickInfo < 3 > );
    cudaMalloc(&bInfo_dev, size);
    cudaMemcpy(bInfo_dev, &_bInfo_dev, size, cudaMemcpyHostToDevice);
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
          expected[pos] = coeff[5] * bIn[b][k + 1][j][i] + coeff[6] * bIn[b][k - 1][j][i] +
                          coeff[3] * bIn[b][k][j + 1][i] + coeff[4] * bIn[b][k][j - 1][i] +
                          coeff[1] * bIn[b][k][j][i + 1] + coeff[2] * bIn[b][k][j][i - 1] +
                          coeff[0] * bIn[b][k][j][i];
        }

  // movBrickStorage
  BrickStorage bStorage_dev = movBrickStorage(bStorage, cudaMemcpyHostToDevice);
  
  auto brick_func_trans = [active_brick_count, &bInfo_dev, &bStorage_dev, &coeff_dev]() -> void {
    auto bSize = cal_size<BDIM>::value;
    Brick <Dim<BDIM>, Dim<VFOLD>> bIn(bInfo_dev, bStorage_dev, 0);
    Brick <Dim<BDIM>, Dim<VFOLD>> bOut(bInfo_dev, bStorage_dev, bSize);
    dim3 block(active_brick_count), thread(32);
    d3pt7_brick_trans << < block, thread >> > (active_brick_count, bIn, bOut, coeff_dev);
  };

  std::cout << "d3pt7" << std::endl;
  std::cout << "Trans: " << cutime_func(brick_func_trans) << std::endl;

  cudaMemcpy(bStorage.dat.get(), bStorage_dev.dat.get(),
             bStorage.chunks * bStorage.step * sizeof(bElem),
             cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

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
  cudaFree(_bInfo_dev.adj);
  cudaFree(bInfo_dev);
  cudaFree(coeff_dev);

  if (!matches)
    throw std::runtime_error("sparse brick result mismatch");

 
}

__global__ void
d3pt7_brick_trans(unsigned (*grid)[STRIDEB][STRIDEB], Brick <Dim<BDIM>, Dim<VFOLD>> bIn,
                  Brick <Dim<BDIM>, Dim<VFOLD>> bOut,
                  bElem *coeff) {
  long tk = GB + blockIdx.z;
  long tj = GB + blockIdx.y;
  long ti = GB + blockIdx.x;
  unsigned b = grid[tk][tj][ti];
  brick("7pt.py", VSVEC, (BDIM), (VFOLD), b);
}
