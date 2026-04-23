#include <span>
#include <queue>
#include <optional>
#include <mtet/io.h>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>
#include <igl/parallel_for.h>
#include <igl/readOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <igl/random_points_on_mesh.h>
#include <igl/writePLY.h>
#include <igl/signed_distance.h>
#include <igl/sparse_voxel_grid.h>
#include <igl/upsample.h>
#include <igl/get_seconds.h>
#include <igl/facet_adjacency_matrix.h>
#include <igl/barycentric_coordinates.h>
#include <igl/grid.h>
#include <igl/connected_components.h>
#include <igl/polygon_corners.h>
#include <igl/per_face_normals.h>
#include <igl/slice.h>
#include <igl/per_corner_normals.h>
#include <igl/swept_volume_signed_distance.h>

#include "init_grid.h"
#include "io.h"
#include "col_gridgen.h"
#include "trajectory.h"
#include "timer.h"

#include <sweep/generalized_sweep.h>
#include <lagrange/views.h>
#include <lagrange/utils/invalid.h>

#include <catch2/catch.hpp>
#include <random>

TEST_CASE("graident check through finite differences regarding analytical functions", "[Analytical][examples]") {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> polynomialFunc;

    const double delta = 0.000001;
    std::function<Eigen::RowVector4d(Eigen::RowVector4d)> finite_diff_grad = [&](Eigen::RowVector4d data)->Eigen::RowVector4d
    {
        std::pair<Scalar,Eigen::RowVector4d> valGrad;
        Eigen::RowVector4d grad;
        valGrad = polynomialFunc(data);
        grad << (polynomialFunc(data + Eigen::RowVector4d{delta, 0, 0, 0}).first - polynomialFunc(data - Eigen::RowVector4d{delta, 0, 0, 0}).first) / (2 * delta),
        (polynomialFunc(data + Eigen::RowVector4d{0, delta, 0, 0}).first - polynomialFunc(data - Eigen::RowVector4d{0, delta, 0, 0}).first) / (2 * delta),
        (polynomialFunc(data + Eigen::RowVector4d{0, 0, delta, 0}).first - polynomialFunc(data - Eigen::RowVector4d{0, 0, delta, 0}).first) / (2 * delta),
        (polynomialFunc(data + Eigen::RowVector4d{0, 0, 0, delta}).first - polynomialFunc(data - Eigen::RowVector4d{0, 0, 0, delta}).first) / (2 * delta);
        return grad;
    };
    const double threshold = 0.02;
    std::function<std::array<int, 4>(int)> check_grad = [&](int num)->std::array<int, 4>
    {
        std::array<int, 4> ret{0, 0, 0, 0};
        Eigen::RowVector4d randP, poly_grad, finite_diff;
        int diff_amount_x = 0, diff_amount_y = 0, diff_amount_z = 0, diff_amount_t = 0;
        for (size_t i = 0; i < num; i ++){
            randP << dis(gen), dis(gen), dis(gen), dis(gen);
            poly_grad = polynomialFunc(randP).second;
            finite_diff = finite_diff_grad(randP);
            double diff = (poly_grad - finite_diff).cwiseAbs().maxCoeff();
            if (diff > threshold){
                Eigen::RowVector4d diffList = (poly_grad - finite_diff).cwiseAbs();
                for (size_t i = 0; i < 4; i++){
                    if (diff == diffList(i)){
                        ret[i]++;
                    }
                }
                std::cout << "Random Points Coordinates: " << randP << std::endl;
                std::cout << "Polynomial Computed Gradient: " << poly_grad << std::endl;
                std::cout << "Finite Difference Gradient: " << finite_diff << std::endl;
                std::cout << "---------------------------------------------------------------------------" << std::endl;
            }
        }
        return ret;
    };

    SECTION("sphere through line:") {
        polynomialFunc = sphereLine;
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    SECTION("flipping donut:") {
        polynomialFunc = flippingDonut;
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    SECTION("flipping donut full turn:") {
        polynomialFunc = flippingDonutFullTurn;
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    SECTION("sphere Loop D Loop:") {
        polynomialFunc = sphereLoopDLoop;
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    SECTION("rotating sphere with lift:") {
        polynomialFunc = rotatingSpherewLift;
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    SECTION("ellipsoid through sinusoidal curve:") {
        polynomialFunc = ellipsoidSine;
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
}


TEST_CASE("graident check through finite differences regarding libigl functions", "[LIBIGL][examples]") {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::string function_file = std::string(DATA_DIR) + "/test/sphere_lvl6.obj";
    igl::read_triangle_mesh(function_file,V,F);
    igl::AABB<Eigen::MatrixXd,3> tree;
    tree.init(V,F);
    igl::FastWindingNumberBVH fwn_bvh;
    int order = 2;
    igl::fast_winding_number(V,F,order,fwn_bvh);
    int rotation = 0;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> libiglFunc = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
    {
        Scalar value;
        Eigen::RowVector4d gradient;
        const double iso = 0.001;
        Eigen::RowVector3d P = data.head(3);
        double t = data[3];
        Eigen::RowVector3d running_closest_point = V.row(0);
        double running_sign = 1.0;
        int i;
        double s,sqrd,sqrd2,s2;
        Eigen::Matrix3d VRt,Rt;
        Eigen::RowVector3d xt,vt,pos,c,c2,xyz_grad,point_velocity;
        trajLine3D(t, xt, vt);
        trajLineRot3D(t, Rt, VRt, rotation);

        pos = ((Rt.inverse())*((P - xt).transpose())).transpose();
        // fast winding number
        Eigen::VectorXd w;
        igl::fast_winding_number(fwn_bvh,2.0,pos,w);
        s = 1.-2.*w(0);
        sqrd = tree.squared_distance(V,F,pos,i,c);
//        std::cout << "c: " << c << std::endl;
//        std::cout << "xt: " << xt << std::endl;
//        std::cout << "vt: " << vt << std::endl;
//        std::cout << "Rt: " << Rt << std::endl;
//        std::cout << "VRt: " << VRt << std::endl;
        value = s*sqrt(sqrd);
        Eigen::RowVector3d cp = c - pos;
//        std::cout << "cp: " << cp << std::endl;
        cp.normalize();
//        std::cout << "normalized cp: " << cp << std::endl;
        //std::cout << cp << std::endl;
        xyz_grad  = (-s) * cp * Rt.inverse();
        gradient << xyz_grad;
        //std::cout << xyz_grad << std::endl;
        point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
//        std::cout << VRt * (Rt.inverse()) << std::endl;
//        std::cout << "point velocity: " << point_velocity << std::endl;
        gradient(3) = (-s) * cp.dot(point_velocity);
        //std::cout << s * cp.dot(point_velocity) << std::endl;
        return {value, gradient};
    };
    
    const double delta = 0.000001;
    std::function<Eigen::RowVector4d(Eigen::RowVector4d)> finite_diff_grad = [&](Eigen::RowVector4d data)->Eigen::RowVector4d
    {
        std::pair<Scalar,Eigen::RowVector4d> valGrad;
        Eigen::RowVector4d grad;
        valGrad = libiglFunc(data);
        grad << (libiglFunc(data + Eigen::RowVector4d{delta, 0, 0, 0}).first - libiglFunc(data - Eigen::RowVector4d{delta, 0, 0, 0}).first) / (2 * delta),
        (libiglFunc(data + Eigen::RowVector4d{0, delta, 0, 0}).first - libiglFunc(data - Eigen::RowVector4d{0, delta, 0, 0}).first) / (2 * delta),
        (libiglFunc(data + Eigen::RowVector4d{0, 0, delta, 0}).first - libiglFunc(data - Eigen::RowVector4d{0, 0, delta, 0}).first) / (2 * delta),
        (libiglFunc(data + Eigen::RowVector4d{0, 0, 0, delta}).first - libiglFunc(data - Eigen::RowVector4d{0, 0, 0, delta}).first) / (2 * delta);
        return grad;
    };
    const double threshold = 0.02;
    std::function<std::array<int, 4>(int)> check_grad = [&](int num)->std::array<int, 4>
    {
        std::array<int, 4> ret{0, 0, 0, 0};
        Eigen::RowVector4d randP, libigl_grad, finite_diff;
        int diff_amount_x = 0, diff_amount_y = 0, diff_amount_z = 0, diff_amount_t = 0;
        for (size_t i = 0; i < num; i ++){
            //0.847431 0.598953 0.516438 0.865228
            //randP << dis(gen), dis(gen), dis(gen), dis(gen);
            randP << 0.847431, 0.598953, 0.516438, (double)i / num;
            libigl_grad = libiglFunc(randP).second;
            finite_diff = finite_diff_grad(randP);
            double diff = (libigl_grad - finite_diff).cwiseAbs().maxCoeff();
            if (diff > threshold){
                Eigen::RowVector4d diffList = (libigl_grad - finite_diff).cwiseAbs();
                for (size_t i = 0; i < 4; i++){
                    if (diff == diffList(i)){
                        ret[i]++;
                    }
                }
                std::cout << "Random Points Coordinates: " << randP << std::endl;
                std::cout << "Libigl Computed Gradient: " << libigl_grad << std::endl;
                std::cout << "Finite Difference Gradient: " << finite_diff << std::endl;
                std::cout << "---------------------------------------------------------------------------" << std::endl;
            }
        }
        return ret;
    };
    
    SECTION("sphere through line:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    function_file = std::string(DATA_DIR) + "/test/torus_lvl500.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("torus through line:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    
    function_file = std::string(DATA_DIR) + "/test/bunny.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("translating bunny:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }

    function_file = std::string(DATA_DIR) + "/test/kingkong.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("translating kingkong:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }

    function_file = std::string(DATA_DIR) + "/test/fandisk.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("translating fandisk:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    
    rotation = 1;
    function_file = std::string(DATA_DIR) + "/test/sphere_lvl6.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("sphere through line with rotation:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    function_file = std::string(DATA_DIR) + "/test/torus_lvl500.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("flipping donut:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }

    function_file = std::string(DATA_DIR) + "/test/torus_lvl400_noTurn.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    libiglFunc = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
    {
        Scalar value;
        Eigen::RowVector4d gradient;
        const double iso = 0.001;
        Eigen::RowVector3d P = data.head(3);
        double t = data[3];
        Eigen::RowVector3d running_closest_point = V.row(0);
        double running_sign = 1.0;
        int i;
        double s,sqrd,sqrd2,s2;
        Eigen::Matrix3d VRt,Rt;
        Eigen::RowVector3d xt,vt,pos,c,c2,xyz_grad,point_velocity;
        trajLine3D2(t, xt, vt);
        trajLineRot3Dx(t, Rt, VRt, rotation);
        pos = ((Rt.inverse())*((P - xt).transpose())).transpose();
        // fast winding number
        Eigen::VectorXd w;
        igl::fast_winding_number(fwn_bvh,2.0,pos,w);
        s = 1.-2.*w(0);
        sqrd = tree.squared_distance(V,F,pos,i,c);
        value = s*sqrt(sqrd);
        Eigen::RowVector3d cp = c - pos;
        cp.normalize();
        //std::cout << cp << std::endl;
        xyz_grad  = (-s) * cp * Rt.inverse();
        gradient << xyz_grad;
        //std::cout << xyz_grad << std::endl;
        point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
        gradient(3) = (-s) * cp.dot(point_velocity);
        //std::cout << s * cp.dot(point_velocity) << std::endl;
        return {value, gradient};
    };
    SECTION("flipping donut 2:") {
        igl::WindingNumberAABB<double, int> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
}

// ---------------------------------------------------------------------------
// envelope_edge_id post-conditions
//
// Run the full pipeline on a small example and check two invariants of the
// edge-label inheritance pass:
//   (1) Every valid `envelope_edge_id[eid]` places the arrangement edge
//       within `geom_tol = 1e-6 * bbox_diag` of the claimed input edge.
//   (2) Every invalid `envelope_edge_id[eid]` with |F_U| >= 1 does NOT
//       lie on any input edge of any parent envelope facet.
// ---------------------------------------------------------------------------

namespace {

double dist_point_to_segment_test(
    const Eigen::Matrix<double, 3, 1>& P,
    const Eigen::Matrix<double, 3, 1>& A,
    const Eigen::Matrix<double, 3, 1>& B)
{
    Eigen::Matrix<double, 3, 1> AB = B - A;
    double ab2 = AB.squaredNorm();
    if (ab2 < 1e-30) return (P - A).norm();
    double t = (P - A).dot(AB) / ab2;
    t = std::max(0.0, std::min(1.0, t));
    return (P - (A + t * AB)).norm();
}

void check_envelope_edge_id_post_conditions(const sweep::SweepResult& r)
{
    using Index = sweep::Index;
    const auto& env = r.envelope;
    const auto& arr = r.arrangement;
    const Index invalid_id = lagrange::invalid<Index>();

    REQUIRE(arr.has_attribute("envelope_facet_id"));
    REQUIRE(arr.has_attribute("envelope_edge_id"));
    auto envelope_facet_id =
        lagrange::attribute_vector_view<Index>(arr, "envelope_facet_id");
    auto envelope_edge_id =
        lagrange::attribute_vector_view<Index>(arr, "envelope_edge_id");

    auto V_env = lagrange::vertex_view(env).template cast<double>();
    auto V_arr = lagrange::vertex_view(arr).template cast<double>();

    Eigen::Matrix<double, 3, 1> bb_min = V_env.colwise().minCoeff().transpose();
    Eigen::Matrix<double, 3, 1> bb_max = V_env.colwise().maxCoeff().transpose();
    const double bbox_diag = (bb_max - bb_min).norm();
    // Must match `envelope_edge_id_rel_tol` in src/post_processing.h.
    const double geom_tol = 1e-6 * bbox_diag;

    auto max_endpoint_dist = [&](Index u0, Index u1, Index e_in) -> double {
        auto [iv0, iv1] = env.get_edge_vertices(e_in);
        Eigen::Matrix<double, 3, 1> A = V_env.row(iv0).transpose();
        Eigen::Matrix<double, 3, 1> B = V_env.row(iv1).transpose();
        Eigen::Matrix<double, 3, 1> P0 = V_arr.row(u0).transpose();
        Eigen::Matrix<double, 3, 1> P1 = V_arr.row(u1).transpose();
        return std::max(
            dist_point_to_segment_test(P0, A, B),
            dist_point_to_segment_test(P1, A, B));
    };

    const Index num_edges = arr.get_num_edges();
    std::vector<Index> parents;
    std::vector<Index> f_u;
    parents.reserve(8);
    f_u.reserve(8);

    size_t checked_valid = 0;

    for (Index eid = 0; eid < num_edges; eid++) {
        parents.clear();
        arr.foreach_facet_around_edge(eid, [&](Index fid) {
            parents.push_back(envelope_facet_id[fid]);
        });
        f_u.clear();
        for (size_t i = 0; i < parents.size(); i++) {
            size_t occ = 0;
            for (size_t j = 0; j < parents.size(); j++) {
                if (parents[j] == parents[i]) occ++;
            }
            if (occ == 1) f_u.push_back(parents[i]);
        }

        auto [u0, u1] = arr.get_edge_vertices(eid);
        const Index label = envelope_edge_id[eid];

        if (label != invalid_id) {
            // Post-condition (1): claimed parent edge holds geometrically.
            double d = max_endpoint_dist(u0, u1, label);
            INFO("eid=" << eid << " label=" << label << " dist=" << d
                        << " tol=" << geom_tol);
            REQUIRE(d <= geom_tol);
            checked_valid++;
        } else if (!f_u.empty()) {
            // Post-condition (2): no parent edge of any F_U facet is within tol.
            for (Index p : f_u) {
                Index cb = env.get_facet_corner_begin(p);
                for (Index k = 0; k < 3; k++) {
                    Index v0 = env.get_corner_vertex(cb + k);
                    Index v1 = env.get_corner_vertex(cb + (k + 1) % 3);
                    Index e_in = env.find_edge_from_vertices(v0, v1);
                    if (e_in == invalid_id) continue;
                    double d = max_endpoint_dist(u0, u1, e_in);
                    INFO("eid=" << eid << " (invalid) parent=" << p
                                << " candidate e_in=" << e_in << " dist=" << d
                                << " tol=" << geom_tol);
                    REQUIRE(d > geom_tol);
                }
            }
        }
    }

    REQUIRE(checked_valid > 0);
}

} // namespace

TEST_CASE("envelope_edge_id post-conditions", "[EdgeLabel][examples]")
{
    struct Case
    {
        std::string name;
        std::string sweep_yaml;
        std::string config_yaml;
    };

    std::vector<Case> cases = {
        {"simple",
         std::string(EXAMPLE_DIR) + "/simple/sweep.yaml",
         std::string(EXAMPLE_DIR) + "/simple/config.yaml"},
        {"brush_stroke_blend",
         std::string(EXAMPLE_DIR) + "/brush_stroke_blend/sweep.yaml",
         std::string(EXAMPLE_DIR) + "/brush_stroke_blend/config.yaml"},
    };

    for (const auto& c : cases) {
        SECTION(c.name)
        {
            auto r = sweep::generalized_sweep_from_config(c.sweep_yaml, c.config_yaml);
            check_envelope_edge_id_post_conditions(r);
        }
    }
}
