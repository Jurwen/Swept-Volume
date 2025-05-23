//
//  ref_crit.cpp
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/10/24.
//
//#define no_inside
#include <ref_crit.h>

///Stores the 2d origin for the convex hull check happened in zero-crossing criteria.
std::array<double, 2> query_2d = {0.0, 0.0}; // X, Y

/// bezier ordinates of a 4D cubic simplex
const Eigen::Matrix<double, 35, 5> ls {
    /// Vertices
    {3, 0, 0, 0, 0}, {0, 3, 0, 0, 0}, {0, 0, 3, 0, 0}, {0, 0, 0, 3, 0}, {0, 0, 0, 0, 3},
    
    /// Edges
    {2, 1, 0, 0, 0}, {2, 0, 1, 0, 0}, {2, 0, 0, 1, 0}, {2, 0, 0, 0, 1},
    {1, 2, 0, 0, 0}, {0, 2, 1, 0, 0}, {0, 2, 0, 1, 0}, {0, 2, 0, 0, 1},
    {1, 0, 2, 0, 0}, {0, 1, 2, 0, 0}, {0, 0, 2, 1, 0}, {0, 0, 2, 0, 1},
    {1, 0, 0, 2, 0}, {0, 1, 0, 2, 0}, {0, 0, 1, 2, 0}, {0, 0, 0, 2, 1},
    {1, 0, 0, 0, 2}, {0, 1, 0, 0, 2}, {0, 0, 1, 0, 2}, {0, 0, 0, 1, 2},
    
    /// Faces
    {1, 1, 1, 0, 0}, {1, 1, 0, 1, 0}, {1, 1, 0, 0, 1}, {1, 0, 1, 1, 0},
    {1, 0, 1, 0, 1}, {1, 0, 0, 1, 1}, {0, 1, 1, 1, 0}, {0, 1, 1, 0, 1},
    {0, 1, 0, 1, 1}, {0, 0, 1, 1, 1}
};

/// Hard coded parameters for taking 
const std::array<std::array<size_t, 2>, 15> derMatrix {{{0, 8}, {5, 27}, {6, 29}, {7, 30}, {8, 21}, {9, 12}, {25, 32}, {26, 33}, {27, 22}, {13, 16}, {28, 34}, {29, 23}, {17, 20}, {30, 24}, {21, 4}}};

const std::array<std::array<size_t, 5>, 35> elevMatrix {{{0, 15, 15, 15, 15}, {15, 5, 15, 15, 15}, {15, 15, 9, 15, 15}, {15,
    15, 15, 12, 15}, {15, 15, 15, 15, 14}, {1, 0, 15, 15, 15}, {2, 15,
        0, 15, 15}, {3, 15, 15, 0, 15}, {4, 15, 15, 15, 0}, {5, 1, 15, 15,
            15}, {15, 6, 5, 15, 15}, {15, 7, 15, 5, 15}, {15, 8, 15, 15, 5}, {9,
                15, 2, 15, 15}, {15, 9, 6, 15, 15}, {15, 15, 10, 9, 15}, {15, 15,
                    11, 15, 9}, {12, 15, 15, 3, 15}, {15, 12, 15, 7, 15}, {15, 15, 12,
                        10, 15}, {15, 15, 15, 13, 12}, {14, 15, 15, 15, 4}, {15, 14, 15, 15,
                            8}, {15, 15, 14, 15, 11}, {15, 15, 15, 14, 13}, {6, 2, 1, 15,
                                15}, {7, 3, 15, 1, 15}, {8, 4, 15, 15, 1}, {10, 15, 3, 2, 15}, {11,
                                    15, 4, 15, 2}, {13, 15, 15, 4, 3}, {15, 10, 7, 6, 15}, {15, 11, 8,
                                        15, 6}, {15, 13, 15, 8, 7}, {15, 15, 13, 11, 10}}};

/// returns a `bool` value that `true` represents positive and `false` represents negative of the input value `x`.
bool get_sign(double x) {
    return x > 0;
}

std::array<int, 16> topFIndices = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 20, 21, 23, 26};
std::array<int, 16> botFIndices = {5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 26, 27, 28, 29};

bool bezier4D(
              const std::array<Eigen::RowVector4d, 5>& pts,
              const std::array<double, 5>& vals,
              const std::array<Eigen::RowVector4d, 5>& grads,
              Eigen::RowVector<double, 35>& bezier,
              bool& inside)
{
    // Decompose inputs
    auto& p1 = pts[0];
    auto& p2 = pts[1];
    auto& p3 = pts[2];
    auto& p4 = pts[3];
    auto& p5 = pts[4];
    
    double v1 = vals[0];
    double v2 = vals[1];
    double v3 = vals[2];
    double v4 = vals[3];
    double v5 = vals[4];
    
    const auto& g1 = grads[0];
    const auto& g2 = grads[1];
    const auto& g3 = grads[2];
    const auto& g4 = grads[3];
    const auto& g5 = grads[4];
    
    auto p12 = p1 - p2;
    auto p13 = p1 - p3;
    auto p14 = p1 - p4;
    auto p15 = p1 - p5;
    auto p23 = p2 - p3;
    auto p24 = p2 - p4;
    auto p25 = p2 - p5;
    auto p34 = p3 - p4;
    auto p35 = p3 - p5;
    auto p45 = p4 - p5;
    // Compute edge values
    std::array<double, 4> v1s = {
        v1 - g1.dot(p12) / 3,
        v1 - g1.dot(p13) / 3,
        v1 - g1.dot(p14) / 3,
        v1 - g1.dot(p15) / 3
    };
    
    std::array<double, 4> v2s = {
        v2 + g2.dot(p12) / 3,
        v2 - g2.dot(p23) / 3,
        v2 - g2.dot(p24) / 3,
        v2 - g2.dot(p25) / 3
    };
    
    std::array<double, 4> v3s = {
        v3 + g3.dot(p13) / 3,
        v3 + g3.dot(p23) / 3,
        v3 - g3.dot(p34) / 3,
        v3 - g3.dot(p35) / 3
    };
    
    std::array<double, 4> v4s = {
        v4 + g4.dot(p14) / 3,
        v4 + g4.dot(p24) / 3,
        v4 + g4.dot(p34) / 3,
        v4 - g4.dot(p45) / 3
    };
    
    std::array<double, 4> v5s = {
        v5 + g5.dot(p15) / 3,
        v5 + g5.dot(p25) / 3,
        v5 + g5.dot(p35) / 3,
        v5 + g5.dot(p45) / 3
    };
    
    // Compute face values
    double e1 = (v1s[0] + v1s[1] + v2s[0] + v2s[1] + v3s[0] + v3s[1]) / 6;
    double face1 = e1 + (e1 - (v1 + v2 + v3) / 3) / 2;
    double e2 = (v1s[0] + v1s[2] + v2s[0] + v2s[2] + v4s[0] + v4s[1]) / 6;
    double face2 = e2 + (e2 - (v1 + v2 + v4) / 3) / 2;
    
    double e3 = (v1s[0] + v1s[3] + v2s[0] + v2s[3] + v5s[0] + v5s[1]) / 6;
    double face3 = e3 + (e3 - (v1 + v2 + v5) / 3) / 2;
    
    double e4 = (v1s[1] + v1s[2] + v3s[0] + v3s[2] + v4s[0] + v4s[2]) / 6;
    double face4 = e4 + (e4 - (v1 + v3 + v4) / 3) / 2;
    
    double e5 = (v1s[1] + v1s[3] + v3s[0] + v3s[3] + v5s[0] + v5s[2]) / 6;
    double face5 = e5 + (e5 - (v1 + v3 + v5) / 3) / 2;
    
    double e6 = (v1s[2] + v1s[3] + v4s[0] + v4s[3] + v5s[0] + v5s[3]) / 6;
    double face6 = e6 + (e6 - (v1 + v4 + v5) / 3) / 2;
    
    double e7 = (v2s[1] + v2s[2] + v3s[1] + v3s[2] + v4s[1] + v4s[2]) / 6;
    double face7 = e7 + (e7 - (v2 + v3 + v4) / 3) / 2;
    
    double e8 = (v2s[1] + v2s[3] + v3s[1] + v3s[3] + v5s[1] + v5s[2]) / 6;
    double face8 = e8 + (e8 - (v2 + v3 + v5) / 3) / 2;
    
    double e9 = (v2s[2] + v2s[3] + v4s[1] + v4s[3] + v5s[1] + v5s[3]) / 6;
    double face9 = e9 + (e9 - (v2 + v4 + v5) / 3) / 2;
    
    double e10 = (v3s[2] + v3s[3] + v4s[2] + v4s[3] + v5s[2] + v5s[3]) / 6;
    double face10 = e10 + (e10 - (v3 + v4 + v5) / 3) / 2;
    // Combine results into a single row vector
    bezier << v1, v2, v3, v4, v5,
    v1s[0], v1s[1], v1s[2], v1s[3],
    v2s[0], v2s[1], v2s[2], v2s[3],
    v3s[0], v3s[1], v3s[2], v3s[3],
    v4s[0], v4s[1], v4s[2], v4s[3],
    v5s[0], v5s[1], v5s[2], v5s[3],
    face1, face2, face3, face4, face5, face6, face7, face8, face9, face10;
    //    bezier.segment(25, 10) << face1, face2, face3, face4, face5,
    //    face6, face7, face8, face9, face10;
#ifndef no_inside
    inside = inside || bezier({0,1,2,3,5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 25, 26, 28, 31}).maxCoeff() <= 0;
    inside = inside || bezier({1,2,3,4,10, 11, 12, 14, 15, 16, 18, 19, 20, 22, 23, 24, 31, 32, 33, 34}).maxCoeff() <= 0;
#endif
    if (get_sign(bezier.maxCoeff()) == get_sign(bezier.minCoeff())){
        return false;
    }
    return true;
}

void  adjugate(const Eigen::Matrix<double, 4, 4>& mat,
               Eigen::Matrix4d& adjugate) {
    //Eigen::Matrix4d adjugate;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Eigen::Matrix3d minor;
            
            int rowOffset = 0;
            for (int row = 0; row < 4; ++row) {
                if (row == i) {
                    rowOffset = 1;
                    continue;
                }
                int colOffset = 0;
                for (int col = 0; col < 4; ++col) {
                    if (col == j) {
                        colOffset = 1;
                        continue;
                    }
                    minor(row - rowOffset, col - colOffset) = mat(row, col);
                }
            }
            // Compute the cofactor
            double cofactor = std::pow(-1, i + j) * minor.determinant();
            adjugate(j, i) = cofactor; // Note the transpose here
        }
    }
    
    //return adjugate;
}

void bezierElev(Eigen::RowVector<double, 16> ords,
                Eigen::RowVector<double, 35>& bezier){
    for (size_t i = 0; i < 35; i++){
        Eigen::Vector<double, 5> elevRow;
        elevRow << ords[elevMatrix[i][0]], ords[elevMatrix[i][1]], ords[elevMatrix[i][2]], ords[elevMatrix[i][3]], ords[elevMatrix[i][4]];
        bezier[i] = (ls.row(i) * elevRow);
    }
    bezier = bezier / 3.0;
    //return bezierGrad / 3.0;
}


bool bezierDerOrds(const Eigen::RowVector<double, 35>& ords,
                    const std::array<Eigen::RowVector4d, 5> verts,
                    Eigen::RowVector<double, 35>& bezierGrad)
{
    Eigen::RowVector<double, 16> vals;
    double norm = (verts[0] - verts[4]).norm();
    // Loop through rows of derMatrix
    for (int i = 0; i < 15; i++) {
        vals[i] = (ords[derMatrix[i][0]] - ords[derMatrix[i][1]]) * 3.0 / norm; // Normalize by the norm of verts
    }
    vals[15] = 0;
    bezierElev(vals, bezierGrad);
    return get_sign(bezierGrad.maxCoeff()) == get_sign(bezierGrad.minCoeff());
}

std::array<double, 70> parse_convex_points2d(const Eigen::Matrix<double, 2, 35> valList) {
    std::array<double, 70> transposed;
    Eigen::MatrixXd::Map(transposed.data(), 2, 35) = valList;
    return transposed;
}

bool outHullClip2D(Eigen::Matrix<double, 2, 35> pts){
    const double eps = 0.0000001;
    bool r1, r2;
    double t;
    auto perp = [](Eigen::Vector2d data){
        return Eigen::Vector2d{-data[1], data[0]};
    };
    std::array<Eigen::Vector2d, 2> range = {-perp(pts.col(0)), perp(pts.col(0))};
    for (size_t i = 1; i < 35; i++){
        t = range[0].dot(pts.col(i));
        if (t > eps){
            r1 = true;
        }else if (t < -eps){
            r1 = false;
        }else if (perp(range[0]).dot(pts.col(i)) > 0){
            r1 = true;
        }else{
            return false;
        }
        t = range[1].dot(pts.col(i));
        if (t > eps){
            r2 = true;
        }else if (t < -eps){
            r2 = false;
        }else if (perp(range[0]).dot(pts.col(i)) < 0){
            r2 = true;
        }else{
            return false;
        }
        
        if (!r1 && !r2){
            return false;
        }else if (!r1){
            range[0] = -perp(pts.col(i));
        }else if (!r2){
            range[1] = perp(pts.col(i));
        }
    }
    return true;
}

/// See header
bool refineFt(
              const std::array<vertex4d, 5>& verts,
              const double threshold,
              bool& inside,
              bool& choice,
              bool& zeroX,
              std::array<double, timer_amount>& profileTimer){
    auto& p1 = verts[0].coord;
    auto& p2 = verts[1].coord;
    auto& p3 = verts[2].coord;
    auto& p4 = verts[3].coord;
    auto& p5 = verts[4].coord;
    
    auto& v1 = verts[0].valGradList.first;
    auto& v2 = verts[1].valGradList.first;
    auto& v3 = verts[2].valGradList.first;
    auto& v4 = verts[3].valGradList.first;
    auto& v5 = verts[4].valGradList.first;
    
    auto& g1 = verts[0].valGradList.second;
    auto& g2 = verts[1].valGradList.second;
    auto& g3 = verts[2].valGradList.second;
    auto& g4 = verts[3].valGradList.second;
    auto& g5 = verts[4].valGradList.second;
    Eigen::RowVector<double, 35> bezierVals;
    if (!bezier4D({p1, p2, p3, p4, p5}, {v1, v2, v3, v4, v5}, {g1, g2, g3, g4, g5}, bezierVals, inside)){
        return false;
    };
    Eigen::RowVector<double, 35> bezierGrad;
    if (bezierDerOrds(bezierVals, {p1, p2, p3, p4, p5}, bezierGrad)){
        return false;
    }
    Eigen::Matrix<double, 2, 35> nPoints_eigen;
    nPoints_eigen << bezierVals, bezierGrad;
    zeroX = !outHullClip2D(nPoints_eigen);
    if (zeroX){
        auto vec1 = p2 - p1, vec2 = p3 - p1, vec3 = p4 - p1, vec4 = p5 - p1;
        Eigen::Matrix4d vec;
        vec << vec1, vec2, vec3, vec4;
        Eigen::Matrix4d adj;
        adjugate(vec, adj);
        auto v1 = bezierGrad[0];
        auto v2 = bezierGrad[1];
        auto v3 = bezierGrad[2];
        auto v4 = bezierGrad[3];
        auto v5 = bezierGrad[4];
        Eigen::RowVector4d gradList = Eigen::RowVector4d(v2-v1, v3-v1, v4-v1, v5-v1) * adj;
        Eigen::RowVector<double, 30> error = ((bezierGrad.tail(30) - (bezierGrad.head(5) * ls.bottomRows(30).transpose()) / 3) * vec.determinant()).array().abs();
        if (error.maxCoeff() > threshold * gradList.norm()){
            Eigen::RowVector<double, 16> topFError = error(topFIndices);
            Eigen::RowVector<double, 16> botFError = error(botFIndices);
            choice = std::max(error[3], error[16]) > std::min(topFError.maxCoeff(), botFError.maxCoeff());
            return true;
        }
    }
    return false;
}

/// Construct the values of one function at the bezier control points within a tet.
///
/// @param[in] pts          The vertex cooridantes of four tet vertices.
/// @param[in] vals         The function values at four tet vertices.
/// @param[in] grads            The total derivative of the functions in x, y, z, t direction at four tet vertices.
/// @param[out] bezier          The eigen vector of 20 bezier values.
///
/// @return         If 20 bezier values contain zero-crossing.
bool bezier3D(
              const std::array<Eigen::RowVector4d, 4>& pts,
              const Eigen::RowVector4d& vals,
              const std::array<Eigen::RowVector4d, 4>& grads,
              Eigen::RowVector<double, 20>& bezier)
{
    // Decompose inputs
    const auto& p1 = pts[0];
    const auto& p2 = pts[1];
    const auto& p3 = pts[2];
    const auto& p4 = pts[3];
    
    const double& v1 = vals[0];
    const double& v2 = vals[1];
    const double& v3 = vals[2];
    const double& v4 = vals[3];
    
    const auto& g1 = grads[0];
    const auto& g2 = grads[1];
    const auto& g3 = grads[2];
    const auto& g4 = grads[3];
    
    // Compute edge values
    std::array<double, 4> v1s = {
        v1 + g1.dot(p2 - p1) / 3,
        v1 + g1.dot(p3 - p1) / 3,
        v1 + g1.dot(p4 - p1) / 3
    };
    
    std::array<double, 4> v2s = {
        v2 + g2.dot(p3 - p2) / 3,
        v2 + g2.dot(p4 - p2) / 3,
        v2 + g2.dot(p1 - p2) / 3
    };
    
    std::array<double, 4> v3s = {
        v3 + g3.dot(p4 - p3) / 3,
        v3 + g3.dot(p1 - p3) / 3,
        v3 + g3.dot(p2 - p3) / 3
    };
    
    std::array<double, 4> v4s = {
        v4 + g4.dot(p1 - p4) / 3,
        v4 + g4.dot(p2 - p4) / 3,
        v4 + g4.dot(p3 - p4) / 3
    };
    // Compute face values
    double face1 = (9 * (v2s[0] + v2s[1] + v3s[0] + v3s[2] + v4s[1] + v4s[2]) / 6 - v2 - v3 - v4)/ 6;
    double face2 =(9 * (v1s[1] + v1s[2] + v3s[0] + v3s[1] + v4s[0] + v4s[2]) / 6 - v1 - v3 - v4)/ 6;
    double face3 =(9 * (v1s[0] + v1s[2] + v2s[1] + v2s[2] + v4s[0] + v4s[1]) / 6 - v1 - v2 - v4)/ 6;
    double face4 =(9 * (v1s[0] + v1s[1] + v2s[0] + v2s[2] + v3s[1] + v3s[2]) / 6 - v1 - v2 - v3)/ 6;
    
    // Combine results into a single row vector
    bezier << v1, v2, v3, v4,
    v1s[0], v1s[1], v1s[2],
    v2s[0], v2s[1], v2s[2],
    v3s[0], v3s[1], v3s[2],
    v4s[0], v4s[1], v4s[2],
    face1, face2, face3, face4;
    
    if (get_sign(bezier.maxCoeff()) == get_sign(bezier.minCoeff())){
        return false;
    }
    return true;
}

/// Construct the value differences between linear interpolations and bezier approximations at 16 bezier control points (excluding control points at tet vertices)
/// @param[in] valList          The eigen vector of 20 bezier values.
///
/// @return         The value differences at 16 control points.
Eigen::Vector<double, 16> bezierDiff(const Eigen::Vector<double,20> valList)
{
    /// Constant coefficient to obtain linear interpolated values at each bezier control points
    const Eigen::Matrix<double, 16, 4> linear_coeff {{2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}, {0, 2, 1, 0},{0, 2, 0, 1}, {1, 2, 0, 0}, {0, 0, 2, 1}, {1, 0, 2, 0},{0, 1, 2, 0}, {1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2},{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};
    Eigen::Vector<double, 16> linear_val = (linear_coeff * valList.head(4)) / 3;
    return valList.tail(16) - linear_val;
}

/// See header
bool refine3D(std::array<Eigen::RowVector4d, 4> pts,
              Eigen::RowVector<double, 4> vals,
              std::array<Eigen::RowVector4d, 4> grads,
              const double threshold){
    Eigen::RowVector<double, 20> bezierVals;
    if (!bezier3D(pts, vals, grads, bezierVals)){
        return false;
    };
    Eigen::Vector3d vec1 = pts[1].head(3) - pts[0].head(3),
    vec2 = pts[2].head(3) - pts[0].head(3),
    vec3 = pts[3].head(3) - pts[0].head(3);
    Eigen::Matrix3d vec;
    vec << vec1, vec2, vec3;
    double D = vec.determinant();
    Eigen::Matrix3d crossMatrix;
    crossMatrix << vec2.cross(vec3), vec3.cross(vec1), vec1.cross(vec2);
    Eigen::Vector3d unNormF = Eigen::RowVector3d(vals(1)-vals(0), vals(2)-vals(0), vals(3)-vals(0)) * crossMatrix.transpose();
    double lhs = (bezierDiff(bezierVals) * D).array().abs().maxCoeff();
    double rhs = threshold * unNormF.norm();
    if (lhs > rhs) {
        return true;
    }
    return false;
}
