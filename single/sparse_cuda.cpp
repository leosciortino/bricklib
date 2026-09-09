//
// Created by Tuowen Zhao on 12/5/18.
//

#include <iostream>
#include "stencils/stencils_cu.h"
#include <random>
#include "brick.h"
#include "stencils/stencils.h"

bElem *coeff;

int main(int argc, char **argv) {
  double line_width = argc > 1 ? std::stod(argv[1]) : 1.0;
  if (line_width < 0.0) {
    std::cerr << "line width must be non-negative" << std::endl;
    return 1;
  }

  coeff = (bElem *) malloc(129 * sizeof(bElem));
  std::random_device r;
  std::mt19937_64 mt(r());
  std::uniform_real_distribution<bElem> u(0, 1);

  for (int i = 0; i < 129; ++i)
    coeff[i] = u(mt);

  d3pt7_sparsecu(line_width);
  free(coeff);
  return 0;
}
