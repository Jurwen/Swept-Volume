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
    std::string function_file = std::string(TEST_FILE) + "/test/sphere_lvl6.obj";
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
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    function_file = std::string(TEST_FILE) + "/test/torus_lvl500.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("torus through line:") {
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    
    function_file = std::string(TEST_FILE) + "/test/bunny.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("translating bunny:") {
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }

    function_file = std::string(TEST_FILE) + "/test/kingkong.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("translating kingkong:") {
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }

    function_file = std::string(TEST_FILE) + "/test/fandisk.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("translating fandisk:") {
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
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
    function_file = std::string(TEST_FILE) + "/test/sphere_lvl6.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("sphere through line with rotation:") {
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }
    function_file = std::string(TEST_FILE) + "/test/torus_lvl500.obj";
    igl::read_triangle_mesh(function_file,V,F);
    tree.init(V,F);
    igl::fast_winding_number(V,F,order,fwn_bvh);
    SECTION("flipping donut:") {
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
        hier.set_mesh(V,F);
        hier.grow();
        std::array<int, 4> error_num = check_grad(100);
        std::cout << "error in x: " << error_num[0]  << " error in y: " << error_num[1] << " error in z: " << error_num[2] << " error in t: " << error_num[3] << std::endl;
        REQUIRE(error_num[0] == 0);
        REQUIRE(error_num[1] == 0);
        REQUIRE(error_num[2] == 0);
        REQUIRE(error_num[3] == 0);
    }

    function_file = std::string(TEST_FILE) + "/test/torus_lvl400_noTurn.obj";
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
        igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
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
