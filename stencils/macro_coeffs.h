#ifndef BRICK_MACROS_COEFFS_H
#define BRICK_MACROS_COEFFS_H

#include "brick.h"

//STAR_STENCIL_RADIUS-valid values are 1,2,3,4
//CUBE_STENCIL_RADIUS-valid values are 1,2
#ifdef STENCIL_RADIUS_1
#define STAR_STENCIL_RADIUS 1
#define CUBE_STENCIL_RADIUS 1
#elif STENCIL_RADIUS_2
#define STAR_STENCIL_RADIUS 2
#define CUBE_STENCIL_RADIUS 2
#elif STENCIL_RADIUS_3
#define STAR_STENCIL_RADIUS 3
#define CUBE_STENCIL_RADIUS 1
#elif STENCIL_RADIUS_4
#define STAR_STENCIL_RADIUS 4
#define CUBE_STENCIL_RADIUS 1
#endif

#define MPI_BETA 0.1
#define MPI_ALPHA 0.4

#define MPI_A0 0.1
#define MPI_A1 0.06
#define MPI_A2 0.045
#define MPI_A3 0.03
#define MPI_A4 0.015

#define MPI_B0 0.4
#define MPI_B1 0.07
#define MPI_B2 0.03

#define MPI_C0 0.1
#define MPI_C1 0.04
#define MPI_C2 0.03
#define MPI_C3 0.01
#define MPI_C4 0.006
#define MPI_C5 0.004
#define MPI_C6 0.005
#define MPI_C7 0.002
#define MPI_C8 0.003
#define MPI_C9 0.001

namespace stencil_reference {

template <typename BrickType>
inline bElem star_7(BrickType &bIn, unsigned b,
                    long k, long j, long i) {
    return (bIn[b][k + 1][j][i] + bIn[b][k - 1][j][i] + 
            bIn[b][k][j + 1][i] + bIn[b][k][j - 1][i] + 
            bIn[b][k][j][i + 1] + bIn[b][k][j][i - 1]) * MPI_BETA + 
            bIn[b][k][j][i] * MPI_ALPHA;

}

template <typename BrickType>
inline bElem star_13(BrickType &bIn, unsigned b,
                    long k, long j, long i) {
    return (bIn[b][k + 2][j][i] + bIn[b][k - 2][j][i] + 
            bIn[b][k][j + 2][i] + bIn[b][k][j - 2][i] + 
            bIn[b][k][j][i + 2] + bIn[b][k][j][i - 2]) * MPI_B2 + 
            (bIn[b][k + 1][j][i] + bIn[b][k - 1][j][i] + 
            bIn[b][k][j + 1][i] + bIn[b][k][j - 1][i] + 
            bIn[b][k][j][i + 1] + bIn[b][k][j][i - 1]) * MPI_B1 + 
            bIn[b][k][j][i] * MPI_B0;
}

template <typename BrickType>
inline bElem star_19(BrickType &bIn, unsigned b,
                    long k, long j, long i) {
    return (bIn[b][k + 3][j][i] + bIn[b][k - 3][j][i] + 
            bIn[b][k][j + 3][i] + bIn[b][k][j - 3][i] + 
            bIn[b][k][j][i + 3] + bIn[b][k][j][i - 3]) * MPI_A3 + 
            (bIn[b][k + 2][j][i] + bIn[b][k - 2][j][i] + 
            bIn[b][k][j + 2][i] + bIn[b][k][j - 2][i] + 
            bIn[b][k][j][i + 2] + bIn[b][k][j][i - 2]) * MPI_A2 + 
            (bIn[b][k + 1][j][i] + bIn[b][k - 1][j][i] + 
            bIn[b][k][j + 1][i] + bIn[b][k][j - 1][i] + 
            bIn[b][k][j][i + 1] + bIn[b][k][j][i - 1]) * MPI_A1 + 
            bIn[b][k][j][i] * MPI_A0;
}
template <typename BrickType>
inline bElem star_25(BrickType &bIn, unsigned b,
                    long k, long j, long i) {
    return (bIn[b][k + 4][j][i] + bIn[b][k - 4][j][i] + 
            bIn[b][k][j + 4][i] + bIn[b][k][j - 4][i] + 
            bIn[b][k][j][i + 4] + bIn[b][k][j][i - 4]) * MPI_A4 + 
            (bIn[b][k + 3][j][i] + bIn[b][k - 3][j][i] + 
            bIn[b][k][j + 3][i] + bIn[b][k][j - 3][i] + 
            bIn[b][k][j][i + 3] + bIn[b][k][j][i - 3]) * MPI_A3 + 
            (bIn[b][k + 2][j][i] + bIn[b][k - 2][j][i] + 
            bIn[b][k][j + 2][i] + bIn[b][k][j - 2][i] + 
            bIn[b][k][j][i + 2] + bIn[b][k][j][i - 2]) * MPI_A2 + 
            (bIn[b][k + 1][j][i] + bIn[b][k - 1][j][i] + 
            bIn[b][k][j + 1][i] + bIn[b][k][j - 1][i] + 
            bIn[b][k][j][i + 1] + bIn[b][k][j][i - 1]) * MPI_A1 + 
            bIn[b][k][j][i] * MPI_A0;
}


template <typename BrickType>
inline bElem cube_27(BrickType &bIn, unsigned b,
                    long k, long j, long i) {
    return (MPI_C0 * bIn[b][k][j][i] + 
       MPI_C1 * (bIn[b][k + 1][j][i] + 
                 bIn[b][k - 1][j][i] + 
                 bIn[b][k][j + 1][i] + 
                 bIn[b][k][j - 1][i] + 
                 bIn[b][k][j][i + 1] + 
                 bIn[b][k][j][i - 1]) +
       MPI_C3 * (bIn[b][k + 1][j + 1][i] + 
                 bIn[b][k - 1][j + 1][i] + 
                 bIn[b][k + 1][j - 1][i] + 
                 bIn[b][k - 1][j - 1][i] + 
                 bIn[b][k + 1][j][i + 1] + 
                 bIn[b][k - 1][j][i + 1] + 
                 bIn[b][k + 1][j][i - 1] + 
                 bIn[b][k - 1][j][i - 1] + 
                 bIn[b][k][j + 1][i + 1] + 
                 bIn[b][k][j - 1][i + 1] + 
                 bIn[b][k][j + 1][i - 1] + 
                 bIn[b][k][j - 1][i - 1]) +
       MPI_C6 * (bIn[b][k + 1][j + 1][i + 1] + 
                 bIn[b][k - 1][j + 1][i + 1] + 
                 bIn[b][k + 1][j - 1][i + 1] + 
                 bIn[b][k - 1][j - 1][i + 1] + 
                 bIn[b][k + 1][j + 1][i - 1] + 
                 bIn[b][k - 1][j + 1][i - 1] + 
                 bIn[b][k + 1][j - 1][i - 1] + 
                 bIn[b][k - 1][j - 1][i - 1]) );
}

template <typename BrickType>
inline bElem cube_125(BrickType &bIn, unsigned b,
                    long k, long j, long i) {
    return (MPI_C0 * bIn[b][k][j][i] +
       MPI_C1 * (bIn[b][k + 1][j][i] + 
                 bIn[b][k - 1][j][i] + 
                 bIn[b][k][j + 1][i] + 
                 bIn[b][k][j - 1][i] + 
                 bIn[b][k][j][i + 1] + 
                 bIn[b][k][j][i - 1]) + 
       MPI_C2 * (bIn[b][k + 2][j][i] + 
                 bIn[b][k - 2][j][i] + 
                 bIn[b][k][j + 2][i] + 
                 bIn[b][k][j - 2][i] + 
                 bIn[b][k][j][i + 2] + 
                 bIn[b][k][j][i - 2]) + 
       MPI_C3 * (bIn[b][k + 1][j + 1][i] + 
                 bIn[b][k - 1][j + 1][i] + 
                 bIn[b][k + 1][j - 1][i] + 
                 bIn[b][k - 1][j - 1][i] + 
                 bIn[b][k + 1][j][i + 1] + 
                 bIn[b][k - 1][j][i + 1] + 
                 bIn[b][k + 1][j][i - 1] + 
                 bIn[b][k - 1][j][i - 1] + 
                 bIn[b][k][j + 1][i + 1] + 
                 bIn[b][k][j - 1][i + 1] + 
                 bIn[b][k][j + 1][i - 1] + 
                 bIn[b][k][j - 1][i - 1]) + 
       MPI_C4 * (bIn[b][k + 1][j + 2][i] + 
                 bIn[b][k - 1][j + 2][i] + 
                 bIn[b][k + 1][j - 2][i] + 
                 bIn[b][k - 1][j - 2][i] + 
                 bIn[b][k + 1][j][i + 2] + 
                 bIn[b][k - 1][j][i + 2] + 
                 bIn[b][k + 1][j][i - 2] + 
                 bIn[b][k - 1][j][i - 2] + 
                 bIn[b][k][j + 1][i + 2] + 
                 bIn[b][k][j - 1][i + 2] + 
                 bIn[b][k][j + 1][i - 2] + 
                 bIn[b][k][j - 1][i - 2] + 
                 bIn[b][k + 2][j + 1][i] + 
                 bIn[b][k - 2][j + 1][i] + 
                 bIn[b][k + 2][j - 1][i] + 
                 bIn[b][k - 2][j - 1][i] + 
                 bIn[b][k + 2][j][i + 1] + 
                 bIn[b][k - 2][j][i + 1] + 
                 bIn[b][k + 2][j][i - 1] + 
                 bIn[b][k - 2][j][i - 1] + 
                 bIn[b][k][j + 2][i + 1] + 
                 bIn[b][k][j - 2][i + 1] + 
                 bIn[b][k][j + 2][i - 1] + 
                 bIn[b][k][j - 2][i - 1]) + 
       MPI_C5 * (bIn[b][k + 2][j + 2][i] + 
                 bIn[b][k - 2][j + 2][i] + 
                 bIn[b][k + 2][j - 2][i] + 
                 bIn[b][k - 2][j - 2][i] + 
                 bIn[b][k + 2][j][i + 2] + 
                 bIn[b][k - 2][j][i + 2] + 
                 bIn[b][k + 2][j][i - 2] + 
                 bIn[b][k - 2][j][i - 2] + 
                 bIn[b][k][j + 2][i + 2] + 
                 bIn[b][k][j - 2][i + 2] +
                 bIn[b][k][j + 2][i - 2] + 
                 bIn[b][k][j - 2][i - 2]) + 
       MPI_C6 * (bIn[b][k + 1][j + 1][i + 1] + 
                 bIn[b][k - 1][j + 1][i + 1] + 
                 bIn[b][k + 1][j - 1][i + 1] + 
                 bIn[b][k - 1][j - 1][i + 1] + 
                 bIn[b][k + 1][j + 1][i - 1] + 
                 bIn[b][k - 1][j + 1][i - 1] + 
                 bIn[b][k + 1][j - 1][i - 1] + 
                 bIn[b][k - 1][j - 1][i - 1]) +
       MPI_C7 * (bIn[b][k + 1][j + 1][i + 2] + 
                 bIn[b][k - 1][j + 1][i + 2] + 
                 bIn[b][k + 1][j - 1][i + 2] + 
                 bIn[b][k - 1][j - 1][i + 2] + 
                 bIn[b][k + 1][j + 1][i - 2] + 
                 bIn[b][k - 1][j + 1][i - 2] + 
                 bIn[b][k + 1][j - 1][i - 2] + 
                 bIn[b][k - 1][j - 1][i - 2] + 
                 bIn[b][k + 1][j + 2][i + 1] + 
                 bIn[b][k - 1][j + 2][i + 1] + 
                 bIn[b][k + 1][j - 2][i + 1] + 
                 bIn[b][k - 1][j - 2][i + 1] + 
                 bIn[b][k + 1][j + 2][i - 1] + 
                 bIn[b][k - 1][j + 2][i - 1] + 
                 bIn[b][k + 1][j - 2][i - 1] + 
                 bIn[b][k - 1][j - 2][i - 1] + 
                 bIn[b][k + 2][j + 1][i + 1] + 
                 bIn[b][k - 2][j + 1][i + 1] + 
                 bIn[b][k + 2][j - 1][i + 1] + 
                 bIn[b][k - 2][j - 1][i + 1] + 
                 bIn[b][k + 2][j + 1][i - 1] + 
                 bIn[b][k - 2][j + 1][i - 1] + 
                 bIn[b][k + 2][j - 1][i - 1] + 
                 bIn[b][k - 2][j - 1][i - 1]) + 
       MPI_C8 * (bIn[b][k + 2][j + 2][i + 1] + 
                 bIn[b][k - 2][j + 2][i + 1] + 
                 bIn[b][k + 2][j - 2][i + 1] + 
                 bIn[b][k - 2][j - 2][i + 1] + 
                 bIn[b][k + 2][j + 2][i - 1] + 
                 bIn[b][k - 2][j + 2][i - 1] + 
                 bIn[b][k + 2][j - 2][i - 1] + 
                 bIn[b][k - 2][j - 2][i - 1] + 
                 bIn[b][k + 2][j + 1][i + 2] + 
                 bIn[b][k - 2][j + 1][i + 2] + 
                 bIn[b][k + 2][j - 1][i + 2] + 
                 bIn[b][k - 2][j - 1][i + 2] + 
                 bIn[b][k + 2][j + 1][i - 2] + 
                 bIn[b][k - 2][j + 1][i - 2] + 
                 bIn[b][k + 2][j - 1][i - 2] + 
                 bIn[b][k - 2][j - 1][i - 2] + 
                 bIn[b][k + 1][j + 2][i + 2] + 
                 bIn[b][k - 1][j + 2][i + 2] + 
                 bIn[b][k + 1][j - 2][i + 2] + 
                 bIn[b][k - 1][j - 2][i + 2] + 
                 bIn[b][k + 1][j + 2][i - 2] + 
                 bIn[b][k - 1][j + 2][i - 2] + 
                 bIn[b][k + 1][j - 2][i - 2] + 
                 bIn[b][k - 1][j - 2][i - 2]) + 
       MPI_C9 * (bIn[b][k + 2][j + 2][i + 2] + 
                 bIn[b][k - 2][j + 2][i + 2] + 
                 bIn[b][k + 2][j - 2][i + 2] + 
                 bIn[b][k - 2][j - 2][i + 2] + 
                 bIn[b][k + 2][j + 2][i - 2] + 
                 bIn[b][k - 2][j + 2][i - 2] + 
                 bIn[b][k + 2][j - 2][i - 2] +
                 bIn[b][k - 2][j - 2][i - 2]) );
}

} // namespace stencil_reference

#endif // BRICK_MACROS_COEFFS_H
