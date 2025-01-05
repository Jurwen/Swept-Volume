//
//  ref_crit.cpp
//  adaptive_column_grid
//
//  Created by Yiwen Ju on 12/10/24.
//

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

Eigen::RowVector<double, 35> bezier4D(
                                      const std::array<Eigen::RowVector4d, 5>& pts,
                                      const std::array<double, 5>& vals,
                                      const std::array<Eigen::RowVector4d, 5>& grads)
{
    // Decompose inputs
    const auto& p1 = pts[0];
    const auto& p2 = pts[1];
    const auto& p3 = pts[2];
    const auto& p4 = pts[3];
    const auto& p5 = pts[4];
    
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
    
    // Compute edge values
    std::array<double, 4> v1s = {
        v1 + g1.dot(p2 - p1) / 3,
        v1 + g1.dot(p3 - p1) / 3,
        v1 + g1.dot(p4 - p1) / 3,
        v1 + g1.dot(p5 - p1) / 3
    };
    
    std::array<double, 4> v2s = {
        v2 + g2.dot(p1 - p2) / 3,
        v2 + g2.dot(p3 - p2) / 3,
        v2 + g2.dot(p4 - p2) / 3,
        v2 + g2.dot(p5 - p2) / 3
    };
    
    std::array<double, 4> v3s = {
        v3 + g3.dot(p1 - p3) / 3,
        v3 + g3.dot(p2 - p3) / 3,
        v3 + g3.dot(p4 - p3) / 3,
        v3 + g3.dot(p5 - p3) / 3
    };
    
    std::array<double, 4> v4s = {
        v4 + g4.dot(p1 - p4) / 3,
        v4 + g4.dot(p2 - p4) / 3,
        v4 + g4.dot(p3 - p4) / 3,
        v4 + g4.dot(p5 - p4) / 3
    };
    
    std::array<double, 4> v5s = {
        v5 + g5.dot(p1 - p5) / 3,
        v5 + g5.dot(p2 - p5) / 3,
        v5 + g5.dot(p3 - p5) / 3,
        v5 + g5.dot(p4 - p5) / 3
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
    Eigen::RowVectorXd result(35);
    result << v1, v2, v3, v4, v5,
    v1s[0], v1s[1], v1s[2], v1s[3],
    v2s[0], v2s[1], v2s[2], v2s[3],
    v3s[0], v3s[1], v3s[2], v3s[3],
    v4s[0], v4s[1], v4s[2], v4s[3],
    v5s[0], v5s[1], v5s[2], v5s[3],
    face1, face2, face3, face4, face5, face6, face7, face8, face9, face10;
    
    return result;
}

Eigen::Matrix4d adjugate(const Eigen::Matrix<double, 4, 4>& mat) {
    Eigen::Matrix4d adjugate;
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
    
    return adjugate;
}

const Eigen::RowVector<double, 35> bezierElev(Eigen::RowVector<double, 15> ords){
    Eigen::RowVector<double, 16> zeroOrds;
    zeroOrds << ords, 0;
    Eigen::RowVector<double, 35> bezierGrad(35);
    for (size_t i = 0; i < 35; i++){
        Eigen::Vector<double, 5> elevRow;
        elevRow << zeroOrds[elevMatrix[i][0]], zeroOrds[elevMatrix[i][1]], zeroOrds[elevMatrix[i][2]], zeroOrds[elevMatrix[i][3]], zeroOrds[elevMatrix[i][4]];
        bezierGrad[i] = ls.row(i) * elevRow;
    }
    return bezierGrad / 3.0;
}


Eigen::RowVector<double, 35> bezierDerOrds(const Eigen::VectorXd& ords, const std::array<Eigen::RowVector4d, 5> verts) {
    Eigen::RowVector<double, 15> vals;
    double norm = (verts[0] - verts[4]).norm();
    // Loop through rows of derMatrix
    for (int i = 0; i < 15; i++) {
        vals[i] = (ords[derMatrix[i][0]] - ords[derMatrix[i][1]]) * 3.0 / norm; // Normalize by the norm of verts
    }
    return bezierElev(vals);
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

std::array<int, 16> topFIndices = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14, 20, 21, 23, 26};
std::array<int, 16> botFIndices = {5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 26, 27, 28, 29};

bool refineContour(
                   const std::array<vertex4d, 5> verts,
                   const double threshold,
                   bool& inside,
                   bool& choice,
                   std::array<double, timer_amount>& profileTimer){
//    Timer ref_crit_timer(ref_crit, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
//    Timer zeroX_crit_one_timer(zeroX_crit_one, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    std::array<Eigen::RowVector4d, 5> pts;
    std::array<double, 5> vals;
    std::array<Eigen::RowVector4d, 5> grads;
    for (size_t i = 0; i < 5; i++){
        pts[i] = verts[i].coord;
        vals[i] = verts[i].valGradList.first;
        grads[i] = verts[i].valGradList.second;
    }
    Eigen::RowVector<double, 35> bezierVals = bezier4D(pts, vals, grads);
    inside = inside || (std::max(bezierVals.head(4).maxCoeff(), (bezierVals.tail(30))(topFIndices).maxCoeff()) <= 0);
    inside = inside || (std::max(bezierVals({1,2,3,4}).maxCoeff(), (bezierVals.tail(30))(botFIndices).maxCoeff()) <= 0);
    if (get_sign(bezierVals.maxCoeff()) == get_sign(bezierVals.minCoeff())){
        return false;
    }
    Eigen::RowVector<double, 35> bezierGrad = bezierDerOrds(bezierVals, pts);
    if (get_sign(bezierGrad.maxCoeff()) == get_sign(bezierGrad.minCoeff())){
        return false;
    }
//    zeroX_crit_one_timer.Stop();
//    Timer zeroX_crit_timer(zeroX_crit_two, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    Eigen::Matrix<double, 2, 35> nPoints_eigen;
    nPoints_eigen << bezierVals, bezierGrad;
    std::array<double, 70> nPoints = parse_convex_points2d(nPoints_eigen);
//    bool zeroX = convex_hull_membership::contains<2, double>(nPoints, query_2d);
    bool zeroX = !outHullClip2D(nPoints_eigen);
//    zeroX_crit_timer.Stop();
//    Timer distance_crit_timer(distance_crit, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    if (zeroX){
        //return true;
        auto vec1 = pts[1] - pts[0], vec2 = pts[2] - pts[0], vec3 = pts[3] - pts[0], vec4 = pts[4] - pts[0];
        Eigen::Matrix4d vec;
        vec << vec1, vec2, vec3, vec4;
        auto adj = adjugate(vec);
        std::array<Eigen::RowVector4d, 2> gradList;
        double v1 = vals[0];
        double v2 = vals[1];
        double v3 = vals[2];
        double v4 = vals[3];
        double v5 = vals[4];
        gradList[0] = Eigen::RowVector4d(v2-v1, v3-v1, v4-v1, v5-v1) * adj;
        v1 = bezierGrad[0];
        v2 = bezierGrad[1];
        v3 = bezierGrad[2];
        v4 = bezierGrad[3];
        v5 = bezierGrad[4];
        gradList[1] = Eigen::RowVector4d(v2-v1, v3-v1, v4-v1, v5-v1) * adj;
        Eigen::Matrix<double, 2, 30> diffList;
        Eigen::RowVector<double, 30>diff1 = bezierVals.tail(30) - (bezierVals.head(5) * ls.bottomRows(30).transpose()) / 3;
        Eigen::RowVector<double, 30>diff2 = bezierGrad.tail(30) - (bezierGrad.head(5) * ls.bottomRows(30).transpose()) / 3;
        diffList << diff1, diff2;
        double D = vec.determinant();
        double gradNorm = (gradList[1][3] * gradList[0] - gradList[0][3] * gradList[1]).norm();
        Eigen::RowVector<double, 30> error = (Eigen::RowVector2d{gradList[1][3], -1 * gradList[0][3]} * diffList).array().abs();
        if (std::abs(error.maxCoeff() * D) > std::abs(threshold * gradNorm)){
            Eigen::RowVector<double, 16> topFError = error(topFIndices);
            Eigen::RowVector<double, 16> botFError = error(botFIndices);
            choice = std::max(error[3], error[16]) > std::min(topFError.maxCoeff(), botFError.maxCoeff());
            return true;
        }
    }
//    distance_crit_timer.Stop();
//    ref_crit_timer.Stop();
    return false;
}

/// Here begins the subdivision criteria for 3D caps.
Eigen::Vector<double, 20> bezierConstruct(const Eigen::RowVector4d vals,
                                          const Eigen::Matrix<double, 4, 3> grads,
                                          const Eigen::Matrix<double, 3, 6> vec)
{
    Eigen::RowVector3d v0s, v1s, v2s, v3s;
    v0s = grads.row(0) * vec(Eigen::all, {0, 1, 2}) / 3;
    v0s.array() += vals(0);
    v1s =  grads.row(1) * vec(Eigen::all, {3, 4, 0}) / 3;
    v1s = v1s.asDiagonal() * Eigen::Vector3d({1, 1, -1});
    v1s.array() += vals(1);
    v2s = grads.row(2) * vec(Eigen::all, {5, 1, 3}) / 3;
    v2s = v2s.asDiagonal() * Eigen::Vector3d({1, -1, -1});
    v2s.array() += vals(2);
    v3s = grads.row(3) * vec(Eigen::all, {2, 4, 5}) / 3;
    v3s *= -1;
    v3s.array() += vals(3);
    
    
    double vMid0 = (9 * (v1s(0) + v1s(1) + v2s(0) + v2s(2) + v3s(1) + v3s(2)) / 6 - vals(1) - vals(2) - vals(3))/ 6;
    double vMid1 =(9 * (v0s[1] + v0s[2] + v2s[0] + v2s[1] + v3s[0] + v3s[2]) / 6 - vals(0) - vals(2) - vals(3))/ 6;
    double vMid2 =(9 * (v0s[0] + v0s[2] + v1s[1] + v1s[2] + v3s[0] + v3s[1]) / 6 - vals(0) - vals(1) - vals(3))/ 6;
    double vMid3 =(9 * (v0s[0] + v0s[1] + v1s[0] + v1s[2] + v2s[1] + v2s[2]) / 6 - vals(0) - vals(1) - vals(2))/ 6;
    Eigen::RowVector<double, 20> valList;
    valList << vals, v0s, v1s, v2s, v3s, vMid0, vMid1, vMid2, vMid3;
    return valList;
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

bool refineCap(
               const std::array<vertex4d, 4> verts,
               const double threshold,
               bool& zeroX)
{
    Eigen::Matrix<double, 4, 3> pts;
    std::array<Eigen::RowVector4d,4> tet_info;
    for (size_t i = 0; i < 4; i++){
        pts.row(i) = verts[i].coord.head(3);
        tet_info[i] = Eigen::RowVector4d{verts[i].valGradList.first, verts[i].valGradList.second(0), verts[i].valGradList.second(1), verts[i].valGradList.second(2)};
    }
    Eigen::RowVector<double, 20> valList;
    Eigen::RowVector<double, 16> diffList;
    Eigen::Vector3d eigenVec1 = pts.row(1) - pts.row(0), eigenVec2 = pts.row(2) - pts.row(0), eigenVec3 = pts.row(3) - pts.row(0), eigenVec4 = pts.row(2) - pts.row(1), eigenVec5 = pts.row(3) - pts.row(1), eigenVec6 = pts.row(3) - pts.row(2);
    Eigen::Matrix<double, 3, 6> vec;
    vec << eigenVec1, eigenVec2, eigenVec3, eigenVec4, eigenVec5, eigenVec6;
    double D = vec.leftCols(3).determinant();
    double sqD = D*D;
    Eigen::Matrix3d crossMatrix;
    crossMatrix << eigenVec2.cross(eigenVec3), eigenVec3.cross(eigenVec1), eigenVec1.cross(eigenVec2);
    //single function linearity check:
    Eigen::Matrix4d func_info;
    func_info << tet_info[0], tet_info[1], tet_info[2], tet_info[3];
    Eigen::RowVector4d vals = func_info.col(0);
    Eigen::Matrix<double, 4, 3> grads_eigen = func_info.rightCols(3);
    valList = bezierConstruct(vals, grads_eigen, vec);
    if (get_sign(valList.maxCoeff()) != get_sign(valList.minCoeff())){
        zeroX = true;
        Eigen::Vector3d unNormF = Eigen::RowVector3d(vals(1)-vals(0), vals(2)-vals(0), vals(3)-vals(0)) * crossMatrix.transpose();
        diffList = bezierDiff(valList);
        double error = std::max(diffList.maxCoeff(), -diffList.minCoeff());
        double lhs = error * error * sqD;
        double rhs = threshold * threshold * unNormF.squaredNorm();
        if (lhs > rhs) {
            return true;
        }
    }
    return false;
}
