//
// Created by Tuowen Zhao on 12/5/18.
//

#include <iostream>
#include <string>
#include <cmath>

#include "stencils/macro_coeffs.h"
#include "stencils/stencils.hip.h"

int main(int argc, char **argv) {
  if (argc > 5) {
    std::cerr << "usage: " << argv[0]
              << " [star|cube] [radius] [width] [line|circle]" << std::endl;
    return 1;
  }

  const std::string stencil_name = argc > 1 ? argv[1] : "star";
  const std::string indicator_name = argc > 4 ? argv[4] : "line";
  unsigned long parsed_stencil_radius;
  double line_width;
  try {
    parsed_stencil_radius = argc > 2 ? std::stoul(argv[2]) : 1;
    line_width = argc > 3 ? std::stod(argv[3]) : 1.0;
  } catch (const std::exception &) {
    std::cerr << "radius and line width must be numeric" << std::endl;
    return 1;
  }

  if (parsed_stencil_radius < 1 || parsed_stencil_radius > 4) {
    std::cerr << "stencil radius must be between 1 and 4" << std::endl;
    return 1;
  }
  const unsigned stencil_radius =
      static_cast<unsigned>(parsed_stencil_radius);

  if (line_width < 0.0) {
    std::cerr << "width must be non-negative" << std::endl;
    return 1;
  }
  if (indicator_name != "line" && indicator_name != "circle") {
    std::cerr << "indicator must be either 'line' or 'circle'" << std::endl;
    return 1;
  }

  // Returns true when (x, y, z) is within line_width of x = y = z.
  const auto in_line = [line_width](long x, long y, long z) {
    return (x - y) * (x - y) + (y - z) * (y - z) +
               (z - x) * (z - x) <
           3 * line_width * line_width;
  };

  // Returns true when (x, y, z) is within circle_width of a sphere
  // centered at (cx, cy, cz) with radius circle_radius.
  const double circle_width = line_width;
  const double circle_radius = N / 4.0;
  const double cx = (N - 1) / 2.0;
  const double cy = cx;
  const double cz = cx;
  const auto in_circle =
      [circle_width, circle_radius, cx, cy, cz](long x, long y, long z) {
        const double dx = x - cx;
        const double dy = y - cy;
        const double dz = z - cz;

        const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

        return std::abs(dist - circle_radius) < circle_width;
      };
  const SparseIndicator indicator =
      indicator_name == "circle" ? SparseIndicator(in_circle)
                                 : SparseIndicator(in_line);
  SparseStencilLauncher launch_stencil;
  SparseReferenceStencil reference_stencil;
  using BrickType = Brick<Dim<BDIM>, Dim<VFOLD>>;

  if (stencil_name == "star") {
    launch_stencil = [stencil_radius](
        unsigned active_brick_count,
        const Brick<Dim<BDIM>, Dim<VFOLD>> &bIn,
        const Brick<Dim<BDIM>, Dim<VFOLD>> &bOut) {
      launch_star_sparsecu(stencil_radius, active_brick_count, bIn, bOut);
    };
    switch (stencil_radius) {
      case 1:
        reference_stencil = &stencil_reference::star_7<BrickType>;
        break;
      case 2:
        reference_stencil = &stencil_reference::star_13<BrickType>;
        break;
      case 3:
        reference_stencil = &stencil_reference::star_19<BrickType>;
        break;
      case 4:
        reference_stencil = &stencil_reference::star_25<BrickType>;
        break;
    }
  } else if (stencil_name == "cube") {
    if (stencil_radius > 2) {
      std::cerr << "cube stencil radius must be 1 or 2" << std::endl;
      return 1;
    }
    launch_stencil = [stencil_radius](
        unsigned active_brick_count,
        const Brick<Dim<BDIM>, Dim<VFOLD>> &bIn,
        const Brick<Dim<BDIM>, Dim<VFOLD>> &bOut) {
      launch_cube_sparsecu(stencil_radius, active_brick_count, bIn, bOut);
    };
    reference_stencil = stencil_radius == 1
                            ? &stencil_reference::cube_27<BrickType>
                            : &stencil_reference::cube_125<BrickType>;
  } else {
    std::cerr << "stencil name must be either 'star' or 'cube'" << std::endl;
    return 1;
  }

  std::cout << "3D " << stencil_name << " radius " << stencil_radius
            << ", " << indicator_name << " width " << line_width
            << std::endl;
  run_sparsecu(indicator, launch_stencil, reference_stencil);
  return 0;
}
