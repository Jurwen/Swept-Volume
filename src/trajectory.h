//
//  trajectory.h
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/5/24.
//

#ifndef trajectory_h
#define trajectory_h

#include "adaptive_column_grid.h"

double sweptFunc();

Eigen::RowVector3d trajLine(double t, std::span<const Scalar, 3> coords){
    Eigen::RowVector3d translation = (1 - t) * Eigen::RowVector3d(1./3., 0., 0.) + t * Eigen::RowVector3d(0, 1./3., 0.) + Eigen::RowVector3d(1./3., 1./3. , 1./2.);
    return Eigen::Map<const Eigen::RowVector3d>(coords.data()) - translation;
}


///hard-coded example 1: sphere traversing through a straight line:
std::pair<Scalar, Eigen::RowVector4d> sphereLine(Eigen::RowVector4d input) {
    // Extract input values {xx, yy, zz, tt}
    Scalar xx = input[0];
    Scalar yy = input[1];
    Scalar zz = input[2];
    Scalar tt = input[3];

    Scalar value = -0.09
        + std::pow((-1.0 / 3.0 + 1.0 / 3.0 * (-1 + tt) + xx), 2)
        + std::pow((-0.333333 - tt / 3.0 + yy), 2)
        + std::pow((-0.5 + zz), 2);

    Eigen::RowVector4d gradient;
    gradient[0] = 2 * (-1.0 / 3.0 + 1.0 / 3.0 * (-1 + tt) + xx); // Partial derivative w.r.t. xx
    gradient[1] = 2 * (-0.333333 - tt / 3.0 + yy);               // Partial derivative w.r.t. yy
    gradient[2] = 2 * (-0.5 + zz);                              // Partial derivative w.r.t. zz
    gradient[3] = 2.0 / 3.0 * (-1.0 / 3.0 + 1.0 / 3.0 * (-1 + tt) + xx)
                  - 2.0 / 3.0 * (-0.333333 - tt / 3.0 + yy);     // Partial derivative w.r.t. tt

    return {value, gradient};
}

std::pair<Scalar, Eigen::RowVector4d> sphereLoopDLoop(Eigen::RowVector4d input) {
    // Extract inputs {xx, yy, zz, tt}
    Scalar xx = input[0];
    Scalar yy = input[1];
    Scalar zz = input[2];
    Scalar tt = input[3];

    // Compute intermediate terms
    Scalar t1 = -0.175 - 4 * std::pow(1 - tt, 2) * tt + 2 * (1 - tt) * std::pow(tt, 2) - (2 * std::pow(tt, 3)) / 3 + xx;
    Scalar t2 = -0.2 - 2 * std::pow(1 - tt, 2) * tt - 2 * (1 - tt) * std::pow(tt, 2) + yy;
    Scalar t3 = -0.5 + zz;

    // Compute the scalar value (first row)
    Scalar value = -0.024025 + std::pow(t1, 2) + std::pow(t2, 2) + std::pow(t3, 2);

    // Compute the gradient (second row)
    Eigen::RowVector4d gradient;

    // Gradient with respect to xx
    gradient[0] = 2 * t1;

    // Gradient with respect to yy
    gradient[1] = 2 * t2;

    // Gradient with respect to zz
    gradient[2] = 2 * t3;

    // Gradient with respect to tt
    gradient[3] = 2 * (
        t1 * (-4 * std::pow(1 - tt, 2) + 12 * (1 - tt) * tt - 4 * std::pow(tt, 2)) +
        t2 * (-2 * std::pow(1 - tt, 2) + 2 * std::pow(tt, 2))
     );
     

    return {value, gradient};
}

std::pair<Scalar, Eigen::RowVector4d> flippingDonut(Eigen::RowVector4d inputs) {
    Scalar value;
    Eigen::RowVector4d gradient;
    // Unpack the inputs
    Scalar xx = inputs[0];
    Scalar yy = inputs[1];
    Scalar zz = inputs[2];
    Scalar tt = inputs[3];
    
    // Constants
    const Scalar pi = 3.14159;
    
    // Precomputed terms
    Scalar cos_pi_tt = std::cos(pi * tt);
    Scalar sin_pi_tt = std::sin(pi * tt);
    Scalar term_zz = -0.51 - 0.01 * tt + zz;
    Scalar term_tt = -0.5 - 0.01 * (1 - tt);
    Scalar term_tt2 = -0.25 - 0.01 * (1 - tt) - 0.51 * tt;
    
    Scalar sqrt_inner = std::sqrt(0.0 + term_zz * term_zz +
                                  std::pow(-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt +
                                           term_tt2 * sin_pi_tt + xx * sin_pi_tt, 2));
    Scalar sqrt_outer = -0.2 + sqrt_inner;
    
    Scalar term1 = -0.01 + term_tt2 * cos_pi_tt + xx * cos_pi_tt -
    term_tt * sin_pi_tt - yy * sin_pi_tt;
    
    // Compute scalar value
    value = -0.0025 + std::pow(term1, 2) + std::pow(sqrt_outer, 2);
    
    // Compute gradient
    // Gradient w.r.t. xx
    gradient[0] = 2 * cos_pi_tt * term1 +
    (2 * sin_pi_tt *
     (-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt +
      term_tt2 * sin_pi_tt + xx * sin_pi_tt) *
     sqrt_outer) /
    sqrt_inner;
    
    // Gradient w.r.t. yy
    gradient[1] = -2 * sin_pi_tt * term1 +
    (2 * cos_pi_tt *
     (-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt +
      term_tt2 * sin_pi_tt + xx * sin_pi_tt) *
     sqrt_outer) /
    sqrt_inner;
    
    // Gradient w.r.t. zz
    gradient[2] = (2 * term_zz * sqrt_outer) / sqrt_inner;
    
    // Gradient w.r.t. tt
    gradient[3] = 2 * (-pi * term_tt * cos_pi_tt - pi * yy * cos_pi_tt -
                       0.01 * sin_pi_tt - pi * term_tt2 * sin_pi_tt -
                       pi * xx * sin_pi_tt) *
    term1 +
    ((-0.02 * term_zz +
      2 * (-0.01 + term_tt * cos_pi_tt + yy * cos_pi_tt +
           term_tt2 * sin_pi_tt + xx * sin_pi_tt) *
      (0.01 * cos_pi_tt + pi * term_tt2 * cos_pi_tt +
       pi * xx * cos_pi_tt - 0.5 * sin_pi_tt -
       pi * term_tt * sin_pi_tt - pi * yy * sin_pi_tt)) *
     sqrt_outer) /
    sqrt_inner;
    return {value, gradient};
}

#endif /* trajectory_h */
