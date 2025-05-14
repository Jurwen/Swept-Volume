//
//  trajectory.h
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/5/24.
//

#ifndef trajectory_h
#define trajectory_h

#include <numbers>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/signed_distance.h>

void trajLine3D(double t, Eigen::RowVector3d& xt, Eigen::RowVector3d& vt) {
    // Define the fixed vectors
    Eigen::RowVector3d start(0.01, 0.01, 0.0);
    Eigen::RowVector3d end(0.51, 0.0, 0.01);
    Eigen::RowVector3d offset(0.25, 0.5, 0.5);
    
    // Compute the linear interpolation and add the offset
    xt = (1.0 - t) * start + t * end + offset;
    vt = end - start;
}

void trajLine3D2(double t, Eigen::RowVector3d& xt, Eigen::RowVector3d& vt) {
    // Define the fixed vectors
    Eigen::RowVector3d start(0.01, 0.01, 0.0);
    Eigen::RowVector3d end(0.0, 0.01, 0.51);
    Eigen::RowVector3d offset(0.5, 0.5, 0.25);
    
    // Compute the linear interpolation and add the offset
    xt = (1.0 - t) * start + t * end + offset;
    vt = end - start;
}

void trajBezier(double t, Eigen::RowVector3d& xt, Eigen::RowVector3d& vt) {
    // Compute position values
    xt(0) = 0.175 + 4 * (1 - t) * (1 - t) * t - 2 * (1 - t) * t * t + (2.0 / 3.0) * t * t * t;
    xt(1) = 0.2 + 2 * (1 - t) * (1 - t) * t + 2 * (1 - t) * t * t;
    xt(2) = 0.5;
    
    // Compute velocity values (derivatives)
    vt(0) = 4 * (1 - t) * (1 - t) - 12 * (1 - t) * t + 4 * t * t;
    vt(1) = 2 * (1 - t) * (1 - t) - 2 * t * t;
    vt(2) = 0.0; // No velocity change in the z-direction
}

void trajLineRot3D(double t, Eigen::Matrix3d& Rt, Eigen::Matrix3d& VRt, const int rotNum) {
    const Scalar pi = 3.14159;
    // Compute sine and cosine of theta
    double cosTheta = std::cos(t * rotNum * pi);
    double sinTheta = std::sin(t * rotNum * pi);
    double dTheta_dt = rotNum * pi;
    // Rotation matrix Rz(theta)
    Rt << cosTheta, -sinTheta, 0,
    sinTheta,  cosTheta, 0,
    0,        0,        1;
    
    VRt << -sinTheta * dTheta_dt, -cosTheta * dTheta_dt, 0,
    cosTheta * dTheta_dt, -sinTheta * dTheta_dt, 0,
    0,        0,        0;
}

void trajLineRot3Dx(double t, Eigen::Matrix3d& Rt, Eigen::Matrix3d& VRt, const int rotNum) {
    const Scalar pi = 3.14159;
    // Compute sine and cosine of theta
    double cosTheta = std::cos(t * rotNum * pi);
    double sinTheta = std::sin(t * rotNum * pi);
    double dTheta_dt = rotNum * pi;
    // Rotation matrix Rz(theta)
    Rt <<
    1, 0, 0,
    0, cosTheta, -sinTheta,
    0, sinTheta,  cosTheta;
    
    VRt <<
    0,        0,        0,
    0, -sinTheta * dTheta_dt, -cosTheta * dTheta_dt,
    0, cosTheta * dTheta_dt, -sinTheta * dTheta_dt;
}

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

std::pair<Scalar, Eigen::RowVector4d> sphereLoopDLoop(Eigen::RowVector4d inputs) {
    Scalar value;
    Eigen::RowVector4d gradient;
    // Unpack the inputs
    Scalar xx = inputs[0];
    Scalar yy = inputs[1];
    Scalar zz = inputs[2];
    Scalar tt = inputs[3];
    
    // Scalar computation
    Scalar term1 = -0.175 - 4 * std::pow(1 - tt, 2) * tt + 2 * (1 - tt) * std::pow(tt, 2) - (2 * std::pow(tt, 3)) / 3 + xx;
    Scalar term2 = -0.2 - 2 * std::pow(1 - tt, 2) * tt - 2 * (1 - tt) * std::pow(tt, 2) + yy;
    Scalar term3 = -0.5 + zz;
    
    value = -0.024025 + std::pow(term1, 2) + std::pow(term2, 2) + std::pow(term3, 2);
    
    // Gradient computation
    gradient[0] = 2 * term1; // Partial derivative w.r.t. xx
    gradient[1] = 2 * term2; // Partial derivative w.r.t. yy
    gradient[2] = 2 * term3; // Partial derivative w.r.t. zz
    
    Scalar d_term1_d_tt = -4 * std::pow(1 - tt, 2) + 12 * (1 - tt) * tt - 4 * std::pow(tt, 2);
    Scalar d_term2_d_tt = -2 * std::pow(1 - tt, 2) + 2 * std::pow(tt, 2);
    
    gradient[3] = 2 * d_term1_d_tt * term1 + 2 * d_term2_d_tt * term2;
    
    
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
    gradient[3] = 2 * (-0.5 * cos_pi_tt -pi * term_tt * cos_pi_tt - pi * yy * cos_pi_tt -
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

// Function to compute the scalar value and gradient
std::pair<Scalar, Eigen::RowVector4d> flippingDonutFullTurn(Eigen::RowVector4d inputs) {
    Scalar value;
    Eigen::RowVector4d gradient;
    // Unpack the inputs
    Scalar xx = inputs[0];
    Scalar yy = inputs[1];
    Scalar zz = inputs[2];
    Scalar tt = inputs[3];
    
    // Constants
    const Scalar pi2 = 6.28319; // 2 * π
    
    // Precomputed terms
    Scalar cos_pi2_tt = std::cos(pi2 * tt);
    Scalar sin_pi2_tt = std::sin(pi2 * tt);
    Scalar term_zz = -0.51 - 0.01 * tt + zz;
    Scalar term_tt = -0.5 - 0.01 * (1 - tt);
    Scalar term_tt2 = -0.25 - 0.01 * (1 - tt) - 0.51 * tt;
    
    Scalar sqrt_inner = std::sqrt(0.0 + term_zz * term_zz +
                                  std::pow(-0.01 + term_tt * cos_pi2_tt + yy * cos_pi2_tt +
                                           term_tt2 * sin_pi2_tt + xx * sin_pi2_tt, 2));
    Scalar sqrt_outer = -0.2 + sqrt_inner;
    
    Scalar term1 = -0.01 + term_tt2 * cos_pi2_tt + xx * cos_pi2_tt -
    term_tt * sin_pi2_tt - yy * sin_pi2_tt;
    
    // Compute scalar value
    value = -0.0025 + std::pow(term1, 2) + std::pow(sqrt_outer, 2);
    
    // Compute gradient
    // Gradient w.r.t. xx
    gradient[0] = 2 * cos_pi2_tt * term1 +
    (2 * sin_pi2_tt *
     (-0.01 + term_tt * cos_pi2_tt + yy * cos_pi2_tt +
      term_tt2 * sin_pi2_tt + xx * sin_pi2_tt) *
     sqrt_outer) /
    sqrt_inner;
    
    // Gradient w.r.t. yy
    gradient[1] = -2 * sin_pi2_tt * term1 +
    (2 * cos_pi2_tt *
     (-0.01 + term_tt * cos_pi2_tt + yy * cos_pi2_tt +
      term_tt2 * sin_pi2_tt + xx * sin_pi2_tt) *
     sqrt_outer) /
    sqrt_inner;
    
    // Gradient w.r.t. zz
    gradient[2] = (2 * term_zz * sqrt_outer) / sqrt_inner;
    
    // Gradient w.r.t. tt
    gradient[3] = 2 * (-0.5 * cos_pi2_tt - pi2 * term_tt * cos_pi2_tt - pi2 * yy * cos_pi2_tt -
                       0.01 * sin_pi2_tt - pi2 * term_tt2 * sin_pi2_tt -
                       pi2 * xx * sin_pi2_tt) *
    term1 +
    ((-0.02 * term_zz +
      2 * (-0.01 + term_tt * cos_pi2_tt + yy * cos_pi2_tt +
           term_tt2 * sin_pi2_tt + xx * sin_pi2_tt) *
      (0.01 * cos_pi2_tt + pi2 * term_tt2 * cos_pi2_tt +
       pi2 * xx * cos_pi2_tt - 0.5 * sin_pi2_tt -
       pi2 * term_tt * sin_pi2_tt - pi2 * yy * sin_pi2_tt)) *
     sqrt_outer) /
    sqrt_inner;
    return {value, gradient};
}

std::pair<Scalar, Eigen::RowVector4d> flippingDonut2(Eigen::RowVector4d input) {
    Scalar x = input(0);
    Scalar y = input(1);
    Scalar z = input(2);
    Scalar t = input(3);
    
    constexpr Scalar pi = 3.14159;
    Scalar cos_pt = std::cos(pi * t);
    Scalar sin_pt = std::sin(pi * t);
    
    // Base terms
    Scalar A = -0.3 - 0.3 * t;
    Scalar B = -0.2 - 0.4 * t;
    Scalar zTerm = -0.45 + z;
    
    // Construct nested expressions
    Scalar nested1 = A * cos_pt + y * cos_pt + B * sin_pt + x * sin_pt;
    Scalar nested2 = B * cos_pt + x * cos_pt - A * sin_pt - y * sin_pt;
    
    Scalar inner1 = 0.6 * nested1 + 0.8 * nested2;
    Scalar firstTerm = inner1 * inner1;
    
    Scalar termY = nested1 - 0.6 * inner1;
    
    Scalar sqrtInner = std::sqrt(zTerm * zTerm + (nested2 - 0.8 * inner1) * (nested2 - 0.8 * inner1) + termY * termY);
    Scalar secondTerm = std::pow(-0.2 + sqrtInner, 2);
    
    Scalar value = -0.0025 + firstTerm + secondTerm;
    
    // Gradient
    Eigen::RowVector4d grad;
    
    // Reuse intermediate terms for chain rule
    Scalar d_inner1_dx = 0.8 * cos_pt + 0.6 * sin_pt;
    Scalar d_inner1_dy = 0.6 * cos_pt - 0.8 * sin_pt;
    
    Scalar d_nested2_dx = cos_pt;
    Scalar d_nested2_dy = -sin_pt;
    Scalar d_termY_dx = sin_pt - 0.6 * (0.6 * sin_pt + 0.8 * d_nested2_dx);
    Scalar d_termY_dy = cos_pt - 0.6 * (0.6 * cos_pt + 0.8 * d_nested2_dy);
    
    Scalar d_sqrt_dx = ((nested2 - 0.8 * inner1) * (d_nested2_dx - 0.8 * d_inner1_dx) + termY * d_termY_dx) / sqrtInner;
    Scalar d_sqrt_dy = ((nested2 - 0.8 * inner1) * (d_nested2_dy - 0.8 * d_inner1_dy) + termY * d_termY_dy) / sqrtInner;
    
    grad(0) = 2 * d_inner1_dx * inner1 + 2 * (-0.2 + sqrtInner) * d_sqrt_dx;
    grad(1) = 2 * d_inner1_dy * inner1 + 2 * (-0.2 + sqrtInner) * d_sqrt_dy;
    grad(2) = 2 * zTerm * (-0.2 + sqrtInner) / sqrtInner;
    
    // Derivatives with respect to t
    Scalar dA_dt = -0.3;
    Scalar dB_dt = -0.4;
    Scalar d_cos_pt = -pi * sin_pt;
    Scalar d_sin_pt = pi * cos_pt;
    
    Scalar d_nested1_dt = dA_dt * cos_pt + A * d_cos_pt + y * d_cos_pt +
    dB_dt * sin_pt + B * d_sin_pt + x * d_sin_pt;
    
    Scalar d_nested2_dt = dB_dt * cos_pt + B * d_cos_pt + x * d_cos_pt -
    dA_dt * sin_pt - A * d_sin_pt - y * d_sin_pt;
    
    Scalar d_inner1_dt = 0.6 * d_nested1_dt + 0.8 * d_nested2_dt;
    
    Scalar d_termY_dt = d_nested1_dt - 0.6 * (0.6 * d_nested1_dt + 0.8 * d_nested2_dt);
    Scalar d_sqrt_dt = ((nested2 - 0.8 * inner1) * (d_nested2_dt - 0.8 * d_inner1_dt) + termY * d_termY_dt) / sqrtInner;
    
    grad(3) = 2 * d_inner1_dt * inner1 + 2 * (-0.2 + sqrtInner) * d_sqrt_dt;
    
    return {value, grad};
}

std::pair<Scalar, Eigen::RowVector4d> rotatingSphere(Eigen::RowVector4d inputs) {
    // Unpack the inputs
    Scalar xx = inputs(0);
    Scalar yy = inputs(1);
    Scalar zz = inputs(2);
    Scalar tt = inputs(3);
    
    // Constants
    const Scalar pi2 = 6.28319; // 2 * π
    const Scalar coeff = 2.51327; // π * 0.8
    
    // Precomputed terms
    Scalar cos_pi2_tt = std::cos(pi2 * tt);
    Scalar sin_pi2_tt = std::sin(pi2 * tt);
    
    // Compute scalar value
    Scalar value = -0.0625 + std::pow(-0.5 + zz, 2) +
    std::pow(-0.5 + yy + 0.2 * cos_pi2_tt, 2) +
    std::pow(-0.5 + xx - 0.2 * sin_pi2_tt, 2);
    
    // Compute gradient
    Eigen::RowVector4d gradient;
    gradient(0) = 2 * (-0.5 + xx - 0.2 * sin_pi2_tt);               // Gradient w.r.t. xx
    gradient(1) = 2 * (-0.5 + yy + 0.2 * cos_pi2_tt);               // Gradient w.r.t. yy
    gradient(2) = 2 * (-0.5 + zz);                                  // Gradient w.r.t. zz
    gradient(3) = -coeff * cos_pi2_tt * (-0.5 + xx - 0.2 * sin_pi2_tt) // Gradient w.r.t. tt
    - coeff * (-0.5 + yy + 0.2 * cos_pi2_tt) * sin_pi2_tt;
    
    // Return the value and gradient as a pair
    return {value, gradient};
}

std::pair<Scalar, Eigen::RowVector4d> rotatingSphere2(Eigen::RowVector4d input) {
    Scalar x = input(0);
    Scalar y = input(1);
    Scalar z = input(2);
    Scalar t = input(3);
    
    constexpr Scalar pi = std::numbers::pi;
    Scalar cos_pt = std::cos(2 * pi * t);
    Scalar sin_pt = std::sin(2 * pi * t);
    
    // Value calculation
    Scalar value = -0.1225 + std::pow(-0.5 + z, 2) + std::pow(-0.5 + y + 0.25 * cos_pt, 2) + std::pow(-0.5 + x - 0.25 * sin_pt, 2);
    
    // Gradient calculation
    Eigen::RowVector4d grad;
    grad(0) = 2 * (-0.5 + x - 0.25 * sin_pt);
    grad(1) = 2 * (-0.5 + y + 0.25 * cos_pt);
    grad(2) = 2 * (-0.5 + z);
    grad(3) = -pi * cos_pt * (-0.5 + x - 0.25 * sin_pt) - pi * sin_pt * (-0.5 + y + 0.25 * cos_pt);
    
    return {value, grad};
}


std::pair<Scalar, Eigen::RowVector4d> rotatingSpherewLift(Eigen::RowVector4d inputs) {
    // Unpack the inputs
    Scalar xx = inputs(0);
    Scalar yy = inputs(1);
    Scalar zz = inputs(2);
    Scalar tt = inputs(3);
    
    // Constants
    const Scalar pi3 = 9.42478;  // 3 * π
    const Scalar coeff = 4.71239; // 1.5 * π
    const Scalar drift = -0.15;   // Coefficient for tt in zz term
    const Scalar zz_offset = -0.4; // Updated offset for zz
    
    // Precomputed terms
    Scalar cos_pi3_tt = std::cos(pi3 * tt);
    Scalar sin_pi3_tt = std::sin(pi3 * tt);
    Scalar term_zz = zz_offset + zz + drift * tt;
    
    // Compute scalar value
    Scalar value = -0.0225 + std::pow(term_zz, 2) +
    std::pow(-0.5 + yy + 0.25 * cos_pi3_tt, 2) +
    std::pow(-0.5 + xx - 0.25 * sin_pi3_tt, 2);
    
    // Compute gradient
    Eigen::RowVector4d gradient;
    gradient(0) = 2 * (-0.5 + xx - 0.25 * sin_pi3_tt);               // Gradient w.r.t. xx
    gradient(1) = 2 * (-0.5 + yy + 0.25 * cos_pi3_tt);               // Gradient w.r.t. yy
    gradient(2) = 2 * term_zz;                                       // Gradient w.r.t. zz
    gradient(3) = 2 * drift * term_zz                                    // Gradient w.r.t. tt
    - coeff * cos_pi3_tt * (-0.5 + xx - 0.25 * sin_pi3_tt)
    - coeff * (-0.5 + yy + 0.25 * cos_pi3_tt) * sin_pi3_tt;
    
    // Return the value and gradient as a pair
    return {value, gradient};
}

std::pair<Scalar, Eigen::RowVector4d> tearDropLine(Eigen::RowVector4d inputs) {
    // Unpack the inputs
    Scalar xx = inputs(0);
    Scalar yy = inputs(1);
    Scalar zz = inputs(2);
    Scalar tt = inputs(3);
    
    // Precomputed terms
    Scalar term_xx = -0.5 - 0.01 * (1.0 - tt) - 0.01 * tt + xx;  // Term involving xx
    Scalar term_yy = -0.25 - 0.01 * (1.0 - tt) - 0.51 * tt + yy; // Term involving yy
    Scalar term_zz = -0.5 - 0.01 * (1.0 - tt) - 0.01 * tt + zz;  // Term involving zz
    
    // Compute scalar value
    Scalar value = -128.0 * std::pow(term_xx, 4) -
    512.0 * std::pow(term_xx, 5) +
    16.0 * std::pow(term_yy, 2) +
    16.0 * std::pow(term_zz, 2);
    
    // Compute gradient
    Eigen::RowVector4d gradient;
    gradient(0) = -512.0 * std::pow(term_xx, 3) - 2560.0 * std::pow(term_xx, 4); // Gradient w.r.t. xx
    gradient(1) = 32.0 * term_yy;                                               // Gradient w.r.t. yy
    gradient(2) = 32.0 * term_zz;                                               // Gradient w.r.t. zz
    gradient(3) = -16.0 * term_yy;                                              // Gradient w.r.t. tt
    
    // Return the value and gradient as a pair
    return {value, gradient};
    
}

std::pair<Scalar, Eigen::RowVector4d> ellipsoidSine(Eigen::RowVector4d inputs) {
    // Unpack the inputs
    Scalar xx = inputs(0);
    Scalar yy = inputs(1);
    Scalar zz = inputs(2);
    Scalar tt = inputs(3);
    
    // Constants
    const Scalar pi2 = 6.28319;  // 2 * π
    const Scalar sin_coeff = 0.3; // Sine coefficient for yy
    const Scalar cos_coeff = 113.097; // Coefficient from gradient
    const Scalar scalar_coeff_xx = 30.0; // Coefficient for xx term
    const Scalar scalar_coeff_yy = 30.0; // Coefficient for yy term
    
    // Precomputed terms
    Scalar term_xx = -0.2 - 0.6 * tt + xx;                  // Term involving xx and tt
    Scalar term_yy = -0.5 + yy - sin_coeff * std::sin(pi2 * tt); // Term involving yy and sine
    Scalar term_zz = -0.5 + zz;                             // Term involving zz
    
    // Compute scalar value
    Scalar value = -0.2025 +
    scalar_coeff_xx * std::pow(term_xx, 2) +
    std::pow(term_zz, 2) +
    scalar_coeff_yy * std::pow(term_yy, 2);
    
    // Compute gradient
    Eigen::RowVector4d gradient;
    gradient(0) = 2 * scalar_coeff_xx * term_xx;                        // Gradient w.r.t. xx
    gradient(1) = 2 * scalar_coeff_yy * term_yy;                        // Gradient w.r.t. yy
    gradient(2) = 2 * term_zz;                                          // Gradient w.r.t. zz
    gradient(3) = -36.0 * term_xx -
    cos_coeff * std::cos(pi2 * tt) * term_yy;             // Gradient w.r.t. tt
    
    // Return the value and gradient as a pair
    return {value, gradient};
}

std::pair<Scalar, Eigen::RowVector4d> ellipsoidLine(Eigen::RowVector4d inputs) {
    // Unpack inputs
    Scalar xx = inputs(0);
    Scalar yy = inputs(1);
    Scalar zz = inputs(2);
    Scalar tt = inputs(3);
    
    // Precomputed terms
    Scalar term_xx = -1.0 / 3.0 + (1.0 / 3.0) * (-1.0 + tt) + xx;
    Scalar term_yy = -0.333333 - tt / 3.0 + yy;
    Scalar term_zz = -0.5 + zz;
    
    // Compute scalar value
    Scalar value = -0.2025 +
    10.0 * std::pow(term_xx, 2) +
    10.0 * std::pow(term_yy, 2) +
    std::pow(term_zz, 2);
    
    // Compute gradient
    Eigen::RowVector4d gradient;
    gradient(0) = 20.0 * term_xx;  // ∂/∂xx
    gradient(1) = 20.0 * term_yy;  // ∂/∂yy
    gradient(2) = 2.0 * term_zz;   // ∂/∂zz
    gradient(3) = (20.0 / 3.0) * term_xx - (20.0 / 3.0) * term_yy;  // ∂/∂tt
    
    // Return the value and gradient as a pair
    return {value, gradient};
}
#endif /* trajectory_h */
