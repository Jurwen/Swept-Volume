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
#include <stf/stf.h>

#include <filesystem>

#include "MeshSDF.h"

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
    Eigen::RowVector3d start(-0.11, 0.01, 0.0);
    Eigen::RowVector3d end(0.61, 0.01, 0.01);
    Eigen::RowVector3d offset(0.25, 0.5, 0.5);
    
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

std::pair<Scalar, Eigen::RowVector4d> bezier(Eigen::RowVector4d inputs) {
    static stf::ImplicitTorus base_shape(0.07, 0.03, {0.0, 0.0, 0.0});
    // stf::ImplicitSphere base_shape(0.07, {0.0, 0.0, 0.0});
    static stf::PolyBezier<3> bezier({{0.2, 0.2, 0.3}, {1.4, 0.8, 0.3}, {-0.4, 0.8, 0.3}, {0.8, 0.2, 0.3}});
    static stf::SweepFunction<3> sweep_function(base_shape, bezier);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> elbow(Eigen::RowVector4d inputs) {
    //stf::ImplicitSphere base_shape(0.2, {0.0, 0.0, 0.0});
    static stf::ImplicitTorus base_shape(0.2, 0.05, {0.0, 0.0, 0.0});
    static stf::Polyline<3> polyline({{0.3, 0.3, 0.3}, {0.7, 0.3, 0.3}, {0.7, 0.7, 0.3}});
    static stf::SweepFunction<3> sweep_function(base_shape, polyline);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> blend_sphere_torus(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere(0.07, {0.0, 0.0, 0.0});
    static stf::ImplicitTorus torus(0.1, 0.04, {0.0, 0.0, 0.0});
    static stf::Polyline<3> polyline({{0.2, 0.5, 0.5}, {0.8, 0.5, 0.5}});
    // stf::ImplicitSphere base_shape(0.05, {0.0, 0.0, 0.0});
    static stf::SweepFunction<3> sphere_sweep(sphere, polyline);
    static stf::SweepFunction<3> torus_sweep(torus, polyline);
    static stf::InterpolateFunction<3> sweep_function(sphere_sweep, torus_sweep);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> blend_spheres(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere(0.2, {0.0, 0.0, 0.0}, 2);
    static stf::ImplicitSphere sphere2(0.2, {0.0, 0.3, 0.0}, 2);
    static stf::ImplicitSphere sphere3(0.2, {0.0, -0.3, 0.0}, 2);
    static stf::ImplicitUnion two_spheres(sphere2, sphere3, 0.03);

    static stf::Polyline<3> polyline({{0.2, 0.5, 0.5}, {0.8, 0.5, 0.5}});
    static stf::SweepFunction<3> sphere_sweep(sphere, polyline);
    static stf::SweepFunction<3> two_sphere_sweep(two_spheres, polyline);
    static stf::InterpolateFunction<3> sweep_function(sphere_sweep, two_sphere_sweep);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}


std::pair<Scalar, Eigen::RowVector4d> sphere_spiral(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere(0.05, {0.0, 0.0, 0.0});
    // clang-format off
    static auto polybezier = stf::PolyBezier<3>::from_samples({
        { 0.500000, 0.050000, 0.500000},
        { 0.522888, 0.056137, 0.570442},
        { 0.381791, 0.074382, 0.585884},
        { 0.326728, 0.104237, 0.374110},
        { 0.585411, 0.144887, 0.237132},
        { 0.831076, 0.195223, 0.500000},
        { 0.616414, 0.253873, 0.858287},
        { 0.166606, 0.319237, 0.742225},
        { 0.147082, 0.389532, 0.243590},
        { 0.638583, 0.462839, 0.073486},
        { 0.948463, 0.537161, 0.500000},
        { 0.634803, 0.610468, 0.914880},
        { 0.166606, 0.680763, 0.742225},
        { 0.195223, 0.746127, 0.278567},
        { 0.602308, 0.804777, 0.185128},
        { 0.776396, 0.855113, 0.500000},
        { 0.566184, 0.895763, 0.703694},
        { 0.381791, 0.925618, 0.585884},
        { 0.440078, 0.943863, 0.456464},
        { 0.500000, 0.950000, 0.500000},
        { 0.425932, 0.943863, 0.500000},
    });
    // clang-format on
    static stf::SweepFunction<3> sweep_function(sphere, polybezier);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.finite_difference_gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> brush_stroke(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere base_shape(0.045, {0.0, 0.0, 0.0});

    // clang-format off
    static std::vector<std::array<stf::Scalar, 3>> samples{
        {0.3090, 0.6504, 0.5},
        {0.3074, 0.6294, 0.5},
        {0.3000, 0.5009, 0.5},
        {0.3934, 0.4126, 0.5},
        {0.4897, 0.3216, 0.5},
        {0.6325, 0.3306, 0.5},
        {0.6360, 0.3432, 0.5},
        {0.6389, 0.3537, 0.5},
        {0.5428, 0.3618, 0.5},
        {0.4755, 0.4415, 0.5},
        {0.3973, 0.5340, 0.5},
        {0.4045, 0.6679, 0.5},
        {0.4223, 0.6732, 0.5},
        {0.4402, 0.6784, 0.5},
        {0.4594, 0.5506, 0.5},
        {0.5693, 0.4865, 0.5},
        {0.6152, 0.4597, 0.5},
        {0.6628, 0.4525, 0.5},
        {0.7000, 0.4514, 0.5},
    };
    // clang-format on
    static stf::PolyBezier<3> transform(samples, false);
    static stf::SweepFunction<3> sweep_function(base_shape, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> brush_stroke_blending(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere0(0.045, {0.0, 0.0, 0.0});
    static stf::ImplicitSphere sphere1(0.02, {0.0, 0.0, 0.03});
    static stf::ImplicitSphere sphere2(0.02, {0.0, 0.0, -0.03});
    static stf::ImplicitUnion<3> base_shape(sphere1, sphere2, 0.002);

    // clang-format off
    static std::vector<std::array<stf::Scalar, 3>> samples{
        {0.3090, 0.6504, 0.32},
        {0.3074, 0.6294, 0.34},
        {0.3000, 0.5009, 0.36},
        {0.3934, 0.4126, 0.38},
        {0.4897, 0.3216, 0.4},
        {0.6325, 0.3306, 0.42},
        {0.6360, 0.3432, 0.44},
        {0.6389, 0.3537, 0.46},
        {0.5428, 0.3618, 0.48},
        {0.4755, 0.4415, 0.5},
        {0.3973, 0.5340, 0.52},
        {0.4045, 0.6679, 0.54},
        {0.4223, 0.6732, 0.56},
        {0.4402, 0.6784, 0.58},
        {0.4594, 0.5506, 0.6},
        {0.5693, 0.4865, 0.62},
        {0.6152, 0.4597, 0.64},
        {0.6628, 0.4525, 0.66},
        {0.7000, 0.4514, 0.68},
    };
    // clang-format on
    static stf::PolyBezier<3> stroke(samples, false);
    static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 360 * 3);
    static stf::Compose<3> rotating_stroke(stroke, rotation);
    //static stf::SweepFunction<3> sweep_function(base_shape, transform);
    static stf::SweepFunction<3> sweep1(sphere0, stroke);
    static stf::SweepFunction<3> sweep2(base_shape, stroke);
    static stf::InterpolateFunction<3> blend(
        sweep1,
        sweep2,
        [](stf::Scalar t) { return (std::sin(t * 3 * 2 * M_PI - M_PI / 2) + 1) / 2; },
        [](stf::Scalar t) { return 3 * M_PI * std::cos(t * 3 * 2 * M_PI - M_PI / 2); });
    auto& sweep_function = blend;

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}


std::pair<Scalar, Eigen::RowVector4d> knot(Eigen::RowVector4d inputs) {
    // static stf::ImplicitSphere base_shape(0.025, {0.0, 0.0, 0.0});
    static stf::ImplicitCapsule<3> base_shape(0.01, {-0.1, 0.0, 0.0}, {0.1, 0, 0});
    static std::vector<std::array<stf::Scalar, 3>> samples{
        { 0.9000, 0.5000, 0.5000 },
        { 0.9000, 0.5832, 0.6464 },
        { 0.6598, 0.5920, 0.7375 },
        { 0.5174, 0.5985, 0.6732 },
        { 0.3750, 0.6050, 0.6089 },
        { 0.4716, 0.6607, 0.3911 },
        { 0.4060, 0.5342, 0.3268 },
        { 0.3405, 0.4077, 0.2625 },
        { 0.2279, 0.1952, 0.3536 },
        { 0.3000, 0.1536, 0.5000 },
        { 0.3721, 0.1120, 0.6464 },
        { 0.4998, 0.3156, 0.7375 },
        { 0.5766, 0.4357, 0.6732 },
        { 0.6535, 0.5557, 0.6089 },
        { 0.6535, 0.4443, 0.3911 },
        { 0.5766, 0.5643, 0.3268 },
        { 0.4998, 0.6844, 0.2625 },
        { 0.3721, 0.8880, 0.3536 },
        { 0.3000, 0.8464, 0.5000 },
        { 0.2279, 0.8048, 0.6464 },
        { 0.3405, 0.5923, 0.7375 },
        { 0.4060, 0.4658, 0.6732 },
        { 0.4716, 0.3392, 0.6089 },
        { 0.3750, 0.3950, 0.3911 },
        { 0.5174, 0.4015, 0.3268 },
        { 0.6598, 0.4080, 0.2625 },
        { 0.9000, 0.4168, 0.3536 },
        { 0.9000, 0.5000, 0.5000 },
    };
    static auto curve = stf::PolyBezier<3>(samples, false);
    //auto curve = stf::Polyline<3>(samples);
    static stf::SweepFunction<3> sweep_function(base_shape, curve);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> concentric_rings(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere1(0.05, {0.6, 0.5, 0.5});
    static stf::ImplicitSphere sphere2(0.05, {0.71, 0.5, 0.5});
    static stf::ImplicitUnion base_shape(sphere1, sphere2);

    static stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    static stf::SweepFunction<3> sweep_function(base_shape, rotation);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> spinning_rod(Eigen::RowVector4d inputs) {
    static stf::ImplicitCapsule<3> rod(0.02, {0.5, 0.5, 0.5}, {0.7, 0.5, 0.5});

    static stf::Rotation<3> spin({0.6, 0.5, 0.5}, {0, 1, 0}, 360 * 3);
    static stf::Rotation<3> rotation({0.5, 0.5, 0.5}, {0, 0, 1});
    static stf::Compose<3> transform(rotation, spin);
    static stf::SweepFunction<3> sweep_function(rod, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> letter_L(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere(0.04, {0.0, 0.0, 0.0});
    static stf::PolyBezier<3> curve(
        {
            {0.6941,0.4189,0.45},
            {0.6457,0.3864,0.4532520646014161},
            {0.5952,0.3504,0.4567053138602828},
            {0.5448,0.3504,0.45951681021477053},
            {0.5076,0.3504,0.461587292026215},
            {0.4752,0.3696,0.46368795156949405},
            {0.4447,0.3899,0.46573106669305236},
            {0.4110,0.4126,0.4679918240358986},
            {0.3755,0.4392,0.4704676882285434},
            {0.3422,0.4392,0.47232022458615164},
            {0.3180,0.4392,0.47367148639993645},
            {0.3000,0.4204,0.4751204549163221},
            {0.3000,0.3993,0.4762973603670379},
            {0.3000,0.3782,0.4774742658177537},
            {0.3192,0.3555,0.4791290732784213},
            {0.3567,0.3555,0.4812213496352494},
            {0.3893,0.3555,0.48304119417478214},
            {0.4193,0.3729,0.48497004222424356},
            {0.4439,0.3966,0.48687369938850344},
            {0.4670,0.4196,0.48869220574687294},
            {0.5081,0.4804,0.49278055529388515},
            {0.5323,0.5170,0.49522667680735677},
            {0.5946,0.6112,0.5015202760795618},
            {0.6257,0.6496,0.5042756679944657},
            {0.6683,0.6496,0.5066512734412809},
            {0.6855,0.6496,0.5076102334381605},
            {0.7000,0.6363,0.508705395789782},
            {0.7000,0.6159,0.5098387121497305},
            {0.7000,0.5741,0.5121707285057785},
            {0.6157,0.4908,0.5187744847481601},
            {0.5025,0.4873,0.5250870538782416},
            {0.4029,0.4842,0.5306437666433996},
            {0.3395,0.5440,0.5355028764354673},
            {0.3395,0.5933,0.5382489891538043},
            {0.3395,0.6245,0.5399925527844943},
            {0.3649,0.6488,0.541950304765777},
            {0.4017,0.6488,0.5439989920318379},
            {0.4498,0.6488,0.5466797211140239},
            {0.4799,0.6022,0.5497688624782774},
            {0.4834,0.6002,0.55},
        },
        false);

    static stf::SweepFunction<3> sweep_function(sphere, curve);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> letter_L_blend(Eigen::RowVector4d inputs) {
    static stf::ImplicitSphere sphere(0.04, {0.0, 0.0, 0.0});
    static stf::ImplicitTorus torus(0.03, 0.015, {0.0, 0.0, 0.0});
    static stf::PolyBezier<3> curve(
        {
            {0.6941,0.4189,0.45},
            {0.6457,0.3864,0.4532520646014161},
            {0.5952,0.3504,0.4567053138602828},
            {0.5448,0.3504,0.45951681021477053},
            {0.5076,0.3504,0.461587292026215},
            {0.4752,0.3696,0.46368795156949405},
            {0.4447,0.3899,0.46573106669305236},
            {0.4110,0.4126,0.4679918240358986},
            {0.3755,0.4392,0.4704676882285434},
            {0.3422,0.4392,0.47232022458615164},
            {0.3180,0.4392,0.47367148639993645},
            {0.3000,0.4204,0.4751204549163221},
            {0.3000,0.3993,0.4762973603670379},
            {0.3000,0.3782,0.4774742658177537},
            {0.3192,0.3555,0.4791290732784213},
            {0.3567,0.3555,0.4812213496352494},
            {0.3893,0.3555,0.48304119417478214},
            {0.4193,0.3729,0.48497004222424356},
            {0.4439,0.3966,0.48687369938850344},
            {0.4670,0.4196,0.48869220574687294},
            {0.5081,0.4804,0.49278055529388515},
            {0.5323,0.5170,0.49522667680735677},
            {0.5946,0.6112,0.5015202760795618},
            {0.6257,0.6496,0.5042756679944657},
            {0.6683,0.6496,0.5066512734412809},
            {0.6855,0.6496,0.5076102334381605},
            {0.7000,0.6363,0.508705395789782},
            {0.7000,0.6159,0.5098387121497305},
            {0.7000,0.5741,0.5121707285057785},
            {0.6157,0.4908,0.5187744847481601},
            {0.5025,0.4873,0.5250870538782416},
            {0.4029,0.4842,0.5306437666433996},
            {0.3395,0.5440,0.5355028764354673},
            {0.3395,0.5933,0.5382489891538043},
            {0.3395,0.6245,0.5399925527844943},
            {0.3649,0.6488,0.541950304765777},
            {0.4017,0.6488,0.5439989920318379},
            {0.4498,0.6488,0.5466797211140239},
            {0.4799,0.6022,0.5497688624782774},
            {0.4834,0.6002,0.55},
        },
        false);

    static stf::SweepFunction<3> sphere_sweep(sphere, curve);
    static stf::Rotation<3> torus_rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 360 * 3);
    static stf::Compose<3> torus_curve(curve, torus_rotation);
    static stf::SweepFunction<3> torus_sweep(torus, torus_curve);
    static stf::OffsetFunction<3> offset_function(
        torus_sweep,
        [](stf::Scalar t) { return - 0.01 * std::sin(t * 6 * M_PI) - 0.01; },
        [](stf::Scalar t) { return - 0.01 * std::cos(t * 6 * M_PI) * 6 * M_PI; });
    auto& sweep_function = offset_function;
    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> torus_rotation(Eigen::RowVector4d inputs) {
    static stf::ImplicitTorus base_shape(0.2, 0.04, {0.25, 0.5, 0.5});
    static stf::Rotation<3> rotation({0.25, 0.5, 0.5}, {1, 0, 0});
    static stf::Translation<3> translation({-0.5, 0, 0});
    static stf::Compose<3> flip(translation, rotation);
    static stf::SweepFunction<3> sweep_function(base_shape, flip);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}


std::pair<Scalar, Eigen::RowVector4d> loopDloop_with_offset(Eigen::RowVector4d inputs) {
    static stf::ImplicitTorus base_shape(0.07, 0.03, {0.0, 0.0, 0.0});
    static stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.3}, {-0.4, 0.8, 0.3}, {0.8, 0.2, 0.3}});
    static stf::SweepFunction<3> torus_sweep(base_shape, bezier);
    static stf::OffsetFunction<3> offset_function(
        torus_sweep,
        [](stf::Scalar t) { return -0.02 * std::cos(t * 2 * M_PI) - 0.02; },
        [](stf::Scalar t) { return 0.02 * std::sin(t * 2 * M_PI) * 2 * M_PI; });
    auto& sweep_function = offset_function;

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> loopDloop_with_offset_v2(Eigen::RowVector4d inputs) {
    static stf::ImplicitTorus base_shape(0.07, 0.03, {0.0, 0.0, 0.0});
    static stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.3}, {-0.4, 0.8, 0.3}, {0.8, 0.2, 0.3}});
    static stf::SweepFunction<3> torus_sweep(base_shape, bezier);
    static stf::OffsetFunction<3> offset_function(
        torus_sweep,
        [](stf::Scalar t) { return -0.04 * t; },
        [](stf::Scalar t) { return -0.04 ; });
    auto& sweep_function = offset_function;

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> doghead(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "doghead_800_shifted.xyz",
            data_dir / "vipss_data" / "doghead_800_shifted_coeff",
            {0.0, 0.0, 0.0},
            0.2, true);
    static stf::Translation<3> translation({-0.5, 0.0, 0.0});
    static stf::Rotation<3> rotation_Y({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}, 180);
    static stf::Rotation<3> rotation_X({0.25, 0.5, 0.5}, {1.0, 0.0, 0.0}, 360);
    static stf::Rotation<3> rotation_Z({0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}, 180);
    static stf::Compose<3> transform(translation, rotation_Y);
    static stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.4}, {-0.4, 0.8, 0.5}, {0.8, 0.2, 0.6}});
    static stf::SweepFunction<3> sweep_function(base_shape, bezier);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> star_S(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    std::string filename = (data_dir / "meshes" / "star2.obj").string();

    static MeshSDF sdf(filename, {0.0, 0.0, 0.0}, 0.1, -0.02);
    //static stf::ImplicitCapsule<3> sdf(0.025, {0.0, -0.06, 0.0}, {0, 0.06, 0});
    static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, 360 * 5);
    //static stf::Translation<3> translation({0.0, -0.2, 0.0});
    static stf::PolyBezier<3> bezier(
        {{0.25, 0.5, 0.25}, {1.5, 0.5, 0.4}, {-0.5, 0.5, 0.6}, {0.75, 0.5, 0.75}});
    static stf::Compose<3> transform(bezier, rotation);
    static stf::SweepFunction<3> sweep_function(sdf, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> star_D(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    std::string filename = (data_dir / "meshes" / "star2.obj").string();

    static MeshSDF sdf(filename, {0.0, 0.0, 0.0}, 0.1, -0.02);
    //static stf::ImplicitCapsule<3> sdf(0.025, {0.0, -0.06, 0.0}, {0, 0.06, 0});
    static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, 360 * 5);
    static stf::PolyBezier<3> bezier(
        {{0.25, 0.5, 0.25}, {0.75, 0.5, 0.4}, {0.75, 0.5, 0.6}, {0.25, 0.5, 0.75}});
    static stf::Compose<3> transform(bezier, rotation);
    static stf::SweepFunction<3> sdf_sweep(sdf, transform);

    static stf::ImplicitSphere sphere(0.05, {0.3, 0.5, 0.8});
    static stf::Translation<3> translation({0.0, 0.0, 0.6});
    static stf::SweepFunction<3> sphere_sweep(sphere, translation);

    static stf::UnionFunction<3> sweep_function(sdf_sweep, sphere_sweep, 0.02);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> star_F(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    std::string filename = (data_dir / "meshes" / "star2.obj").string();

    static MeshSDF sdf(filename, {0.0, 0.0, 0.0}, 0.1, -0.02);
    //static stf::ImplicitCapsule<3> sdf(0.025, {0.0, -0.06, 0.0}, {0, 0.06, 0});
    static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, 360 * 5);
    static stf::PolyBezier<3> bezier(
        {{0.4, 0.5, 0.25}, {0.4, 0.5, 0.6}, {0.3, 0.5, 0.75}, {0.75, 0.5, 0.75}});
    static stf::Compose<3> transform(bezier, rotation);
    static stf::SweepFunction<3> sdf_sweep(sdf, transform);

    static stf::ImplicitSphere sphere(0.05, {0.6, 0.5, 0.5});
    static stf::Translation<3> translation({0.3, 0.0, 0.0});
    static stf::SweepFunction<3> sphere_sweep(sphere, translation);

    static stf::UnionFunction<3> sweep_function(sdf_sweep, sphere_sweep, 0.01);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> star_I(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    std::string filename = (data_dir / "meshes" / "star2.obj").string();

    static MeshSDF sdf(filename, {0.5, 0.5, 0.25}, 0.1, -0.01);
    static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 300);
    static stf::Translation<3> translation({0.0, 0.0, -0.5});
    static stf::Compose<3> transform(translation, rotation);
    static stf::SweepFunction<3> sweep_function(sdf, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.finite_difference_gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> fertility(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "fertility.xyz",
            data_dir / "vipss_data" / "fertility_coeff.txt",
            {0.5, 0.5, 0.5},
            0.3, false);
    static stf::Translation<3> translation({-0.2, 0.0, 0.0});
    static stf::Rotation<3> rotation_Y({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}, 180);
    static stf::Rotation<3> rotation_X({0.25, 0.5, 0.5}, {1.0, 0.0, 0.0}, 360);
    static stf::Rotation<3> rotation_Z({0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}, 180);
    static stf::Compose<3> transform(translation, rotation_Y);
    static stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.4}, {-0.4, 0.8, 0.5}, {0.8, 0.2, 0.6}});
    static stf::SweepFunction<3> sweep_function(base_shape, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> fertility_v2(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "fertility.xyz",
            data_dir / "vipss_data" / "fertility_coeff.txt",
            {0.5, 0.25, 0.5},
            0.3, false);
    static stf::Rotation<3> rotation_Y({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}, 180 * 3);
    static stf::Translation<3> translation({0.0, -0.5, 0.0});
    static stf::Compose<3> transform(translation, rotation_Y);
    static stf::SweepFunction<3> sdf_sweep(base_shape, transform);
    static stf::OffsetFunction<3> sweep_function(
        sdf_sweep,
        [](stf::Scalar t) { return -0.1 * (1-t); },
        [](stf::Scalar t) { return 0.1 ; });

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> fertility_v3(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    // clang-format off
    static std::vector<std::array<stf::Scalar, 3>> samples{
        {0.3090, 0.5, 0.6504},
        {0.3074, 0.5, 0.6294},
        {0.3000, 0.5, 0.5009},
        {0.3934, 0.5, 0.4126},
        {0.4897, 0.5, 0.3216},
        {0.6325, 0.5, 0.3306},
        {0.6360, 0.5, 0.3432},
        {0.6389, 0.5, 0.3537},
        {0.5428, 0.5, 0.3618},
        {0.4755, 0.5, 0.4415},
        {0.3973, 0.5, 0.5340},
        {0.4045, 0.5, 0.6679},
        {0.4223, 0.5, 0.6732},
        {0.4402, 0.5, 0.6784},
        {0.4594, 0.5, 0.5506},
        {0.5693, 0.5, 0.4865},
        {0.6152, 0.5, 0.4597},
        {0.6628, 0.5, 0.4525},
        {0.7000, 0.5, 0.4514},
    };
    // clang-format on

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "fertility.xyz",
            data_dir / "vipss_data" / "fertility_coeff.txt",
            {0.0, 0.0, 0.0},
            0.1, false);
    static stf::PolyBezier<3> bezier(samples, false);
    static stf::SweepFunction<3> sdf_sweep(base_shape, bezier);
    static stf::OffsetFunction<3> sweep_function(
        sdf_sweep,
        [](stf::Scalar t) { return -0.1 * t; },
        [](stf::Scalar t) { return -0.1; });

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> fertility_v4(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "fertility.xyz",
            data_dir / "vipss_data" / "fertility_coeff.txt",
            {0.5, 0.25, 0.5},
            0.3, false);
    static stf::Rotation<3> rotation_Y({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}, 180 * 3);
    static stf::Translation<3> translation({0.0, -0.5, 0.0});
    static stf::Compose<3> transform(translation, rotation_Y);
    static stf::SweepFunction<3> sdf_sweep(base_shape, transform);
    static stf::OffsetFunction<3> sweep_function(
        sdf_sweep,
        [](stf::Scalar t) { return 0.1 * (1-t); },
        [](stf::Scalar t) { return -0.1 ; });

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> fertility_v5(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "fertility.xyz",
            data_dir / "vipss_data" / "fertility_coeff.txt",
            {0.5, 0.25, 0.5},
            0.15, false);
    static stf::Rotation<3> rotation_Y({0.5, 0.5, 0.5}, {0.0, 1.0, 0.0}, 180 * 3);
    static stf::Translation<3> translation({0.0, -0.5, 0.0});
    static stf::Scale<3> scale({0.5, 0.5, 0.5}, {0.5, 0.25, 0.5});
    static stf::Compose<3> scale_rotate(rotation_Y, scale);
    static stf::Compose<3> transform(translation, scale_rotate);
    static stf::SweepFunction<3> sweep_function(base_shape, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> fertility_v6(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);

    static stf::Duchon base_shape(
            data_dir / "vipss_data" / "fertility.xyz",
            data_dir / "vipss_data" / "fertility_coeff.txt",
            {0.0, 0.0, 0.0},
            0.2, false);
    static stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.4}, {-0.4, 0.8, 0.5}, {0.8, 0.2, 0.6}});
    static stf::SweepFunction<3> sdf_sweep(base_shape, bezier);
    static stf::OffsetFunction<3> sweep_function(
        sdf_sweep,
        [](stf::Scalar t) { return 0.1 * (1-t); },
        [](stf::Scalar t) { return -0.1 ; });

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> bunny_blend(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    std::string filename = (data_dir / "meshes" / "bunny.obj").string();

    static MeshSDF sdf(filename, {0.25, 0.5, 0.5}, 0.2, -0.02);
    static stf::Rotation<3> rotation_Z({0.5, 0.5, 0.5}, {0.0, 0.0, 1.0}, 360);
    static stf::Rotation<3> rotation_X({0.25, 0.5, 0.5}, {1.0, 0.0, 0.0}, 360);
    static stf::Compose<3> transform(rotation_Z, rotation_X);
    static stf::SweepFunction<3> sweep_function(sdf, transform);

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}


std::pair<Scalar, Eigen::RowVector4d> loopDloop_with_offset_v3(Eigen::RowVector4d inputs) {
    static stf::ImplicitTorus base_shape(0.07, 0.04, {0.0, 0.0, 0.0});
    static stf::PolyBezier<3> bezier(
        {{0.2, 0.2, 0.3}, {1.4, 0.8, 0.3}, {-0.4, 0.8, 0.3}, {0.8, 0.2, 0.3}});
    static stf::SweepFunction<3> torus_sweep(base_shape, bezier);
    static stf::OffsetFunction<3> offset_function(
        torus_sweep,
        [](stf::Scalar t) { return -0.04 * t; },
        [](stf::Scalar t) { return -0.04 ; });
    auto& sweep_function = offset_function;

    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.finite_difference_gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> VIPSS_blend(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    static stf::Duchon doghead(
                               data_dir / "test" / "dog_head" / "doghead_800.xyz",
                               data_dir / "test" / "dog_head" / "doghead_800_coeff",
                               {0, 0, 0}, 1, true);
    static stf::Duchon kitten(
                               data_dir / "test" / "kitten" / "kitten_893.xyz",
                               data_dir / "test" / "kitten" / "kitten_893.coeff",
                               {0, 0, 0}, 1, true);
    static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 90);
    static stf::PolyBezier<3> bezier({{5, 4., 3.}, {5, 7., 4.5}, {5, 2.5, 4.5}, {5, 5.5, 3.}});
    stf::Compose<3> bezierRot(bezier, rotation);
    static stf::SweepFunction<3> dog_sweep(doghead, bezierRot);
    static stf::SweepFunction<3> kitten_sweep(kitten, bezier);
    static stf::InterpolateFunction<3> sweep_function(kitten_sweep, dog_sweep);
    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.finite_difference_gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

std::pair<Scalar, Eigen::RowVector4d> VIPSS_blend2(Eigen::RowVector4d inputs) {
    std::filesystem::path data_dir(DATA_DIR);
    static stf::Duchon doghead(
                               data_dir / "vipss_data" / "doghead_800.xyz",
                               data_dir / "vipss_data" / "doghead_800_coeff",
                               {5, 5, 3}, 1, true);
    static stf::Duchon kitten(
                               data_dir / "vipss_data" / "kitten_893.xyz",
                               data_dir / "vipss_data" / "kitten_893_coeff",
                               {5, 5, 3}, 1, true);
    //static stf::Rotation<3> rotation({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 90);
    //static stf::PolyBezier<3> bezier({{5, 4., 3.}, {5, 7., 4.5}, {5, 2.5, 4.5}, {5, 5.5, 3.}});
    
    //New Trajectory:
    static stf::Rotation<3> rotation({4.1, 5, 3}, {0.0, 1.0, 0.0}, 180+90);
    static stf::Rotation<3> rotation_dog({5, 5, 3}, {0.0, 1.0, 0.0}, 45);
//    static stf::Translation<3> translation({0, 0.0, 0});
    static stf::Translation<3> translation({-2, 0.0, -1.5});
//    stf::Compose<3> bezierRot(bezier, rotation);
//    static stf::SweepFunction<3> dog_sweep(doghead, bezierRot);
//    static stf::SweepFunction<3> kitten_sweep(kitten, bezier);
//    static stf::InterpolateFunction<3> sweep_function(kitten_sweep, dog_sweep);
    stf::Compose<3> bezierRot(translation, rotation);
    stf::Compose<3> rot_dog(bezierRot, rotation_dog);
    static stf::SweepFunction<3> dog_sweep(doghead, rot_dog);
    static stf::SweepFunction<3> kitten_sweep(kitten, bezierRot);
    static stf::InterpolateFunction<3> sweep_function(kitten_sweep, dog_sweep);
    Scalar value = sweep_function.value({inputs(0), inputs(1), inputs(2)}, inputs(3));
    auto gradient = sweep_function.finite_difference_gradient({inputs(0), inputs(1), inputs(2)}, inputs(3));
    return {value, Eigen::RowVector4d(gradient[0], gradient[1], gradient[2], gradient[3])};
}

#endif /* trajectory_h */
