#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <Marching4D.h>
#include <mtet/io.h>
#include <implicit_functions.h>
#include <chrono>
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

int main(int argc, const char *argv[])
{
    struct
    {
        std::string grid_file;
        std::string function_file;
        double threshold;
        int max_splits = std::numeric_limits<int>::max();
        int rot = 0;
    } args;
    CLI::App app{"Generalized Swept Volume"};
    app.add_option("grid", args.grid_file, "Initial grid file")->required();
    app.add_option("function", args.function_file, "Implicit function file")->required();
    app.add_option("-t,--threshold", args.threshold, "Threshold value");
    app.add_option("-m,--max-splits", args.max_splits, "Maximum number of splits");
    app.add_option("-r,--rotation-number", args.rot, "Number of rotations");
    
    CLI11_PARSE(app, argc, argv);
    // Read initial grid
    mtet::MTetMesh grid;
    if (args.grid_file.find(".json") != std::string::npos){
        grid = init_grid::load_tet_mesh(args.grid_file);
        mtet::save_mesh("init.msh", grid);
        grid = mtet::load_mesh("init.msh");
    } else {
        grid = mtet::load_mesh(args.grid_file);
    }
    
    int max_splits = args.max_splits;
    std::string function_file = args.function_file;
    double threshold = args.threshold;
    int rotation = args.rot;
    
    /// main function:
    /// the lambda function for function evaluations
    ///  @param[in] data            The 4D coordinate
    ///  @return    A std::pari<Scalar, Eigen::RowVector4d> of the value and the gradients at this 4D point
    ///
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(args.function_file,V,F);
    igl::AABB<Eigen::MatrixXd,3> tree;
    tree.init(V,F);
    igl::FastWindingNumberBVH fwn_bvh;
    int order = 2;
    igl::fast_winding_number(V,F,order,fwn_bvh);
    igl::WindingNumberAABB<Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi> hier;
    hier.set_mesh(V,F);
    hier.grow();
    /// Definition of the implicit function using explicit meshes
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
        value = s*sqrt(sqrd) - iso;
        Eigen::RowVector3d cp = pos - c;
        cp.normalize();
        //std::cout << cp << std::endl;
        xyz_grad  = s * cp * Rt.inverse();
        gradient << xyz_grad;
        //std::cout << xyz_grad << std::endl;
        point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
        gradient(3) =  s * cp.dot(point_velocity);
        //std::cout << s * cp.dot(point_velocity) << std::endl;
        return {value, gradient};
    };
    
//    std::cout << libiglFunc(Eigen::RowVector4d({0.5, 0.5, 0.5, 0})).first << std::endl;
//    std::cout << libiglFunc(Eigen::RowVector4d({0.5, 0.5, 0.5, 0})).second << std::endl;
//    const int level = 21;
//    std::vector<double> value_list;
//    std::vector<Eigen::RowVector4d> grad_list;
//    value_list.reserve(pow(level, 4));
//    grad_list.reserve(pow(level, 4));
//    Eigen::RowVector4d input;
//    for (size_t l = 0; l < level; l++){
//        input(3) = (double)l / (double)(level - 1);
//        for (size_t i = 0; i < level; i++){
//            input(0) = (double)i / (double)(level - 1);
//            for (size_t j = 0; j < level; j++){
//                input(1) = (double)j / (double)(level - 1);
//                for (size_t k = 0; k < level; k++){
//                    input(2) = (double)k / (double)(level - 1);
//                    auto ret = libiglFunc(input);
//                    value_list.emplace_back(ret.first);
//                    grad_list.emplace_back(ret.second);
//                }
//            }
//        }
//    }
//    std::string test_file = "test_grad_field.json";
//    if (std::filesystem::exists(test_file.c_str())){
//        std::filesystem::remove(test_file.c_str());
//    }
//    using json = nlohmann::json;
//    std::ofstream fout(test_file.c_str(),std::ios::app);
//    nlohmann::json jOut;
//    jOut.push_back(nlohmann::json(value_list));
//    jOut.push_back(nlohmann::json(grad_list));
//    fout << jOut.dump(4, ' ', true, nlohmann::json::error_handler_t::replace) << std::endl;
//    fout.close();
//
//
//
    
    
    auto implicit_sweep = [&](Eigen::RowVector4d data){
        return sphereLine(data);
    };
//    vertexCol newVert;
//    constexpr mtet::Scalar arr[3] = { 0.25, 0.5, 0.5 };
//    initNewVert(newVert, arr, implicit_sweep, threshold);
//    std::string test_file = "vert_time_sampling.json";
//    if (std::filesystem::exists(test_file.c_str())){
//        std::filesystem::remove(test_file.c_str());
//    }
//    using json = nlohmann::json;
//    std::ofstream fout(test_file.c_str(),std::ios::app);
//    nlohmann::json jOut;
//    for (auto vert : newVert.vert4dList){
//        jOut.push_back(nlohmann::json(vert.coord));
//    }
//    fout << jOut.dump(4, ' ', true, nlohmann::json::error_handler_t::replace) << std::endl;
//    fout.close();
    
    
    

    ///
    ///
    ///Grid generation:
    vertExtrude vertexMap;
    tetExtrude_test insideMap;
    tetExtrude cell5Map;
    std::array<double, timer_amount> profileTimer = {0, 0, 0, 0, 0, 0, 0, 0};
    auto starterTime = std::chrono::high_resolution_clock::now();
    if (!gridRefine(grid, vertexMap, insideMap, cell5Map, implicit_sweep, threshold, max_splits, profileTimer)){
        throw std::runtime_error("ERROR: grid generation failed");
        return 0;
    };
    auto stopperTime = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    auto duration = end - start;
    double ms = duration * 0.001;

    /// save the grid output for discretization tool
    save_mesh_json("grid.json", grid);
    /// write grid and active tets
    mtet::save_mesh("tet_grid.msh", grid);
    /// Start of meshing the swept volume
    std::vector<std::array<double, 4>> verts;
    std::vector<std::array<size_t, 5>> simps;
    std::vector<std::array<size_t, 4>> ulsimp;
    std::vector<std::array<size_t, 4>> llsimp;
    std::vector<double> values;
    //Timer meshing_timer(meshing, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
    /// Convert column-wise data structure to lists of 4D verts and tets
    convert_4d_grid(grid, vertexMap, cell5Map, verts, simps, ulsimp, llsimp, values);
    //meshing_timer.Stop();

    std::vector<std::array<double, 4>> output_vertices_4d;
    std::vector<std::array<size_t, 4>> output_tets_3d;
    std::vector<std::array<size_t, 6>> output_prisms_3d;

    std::vector<double> values_3d;
    std::vector<std::array<double, 4> > gradients_3d;
    std::vector<std::array<double, 3>> output_vertices;
    std::vector<std::array<size_t, 3>> output_triangles;
    /// First marching using the time-derivative function
    marching4D::Marching4SimplexUL(verts, simps, ulsimp, llsimp, values, output_vertices_4d, output_tets_3d
                                   , output_prisms_3d
                                   );
    
//    // build mesh
//    mtet::MTetMesh mesh;
//    std::vector<mtet::VertexId> vertex_ids;
//    vertex_ids.reserve(output_vertices_4d.size());
//    for (auto &v: output_vertices_4d) {
//        vertex_ids.push_back(mesh.add_vertex(v[0], v[1], v[2]));
//    }
//    for (auto &t: output_tets_3d) {
//        mesh.add_tet(vertex_ids[t[0]], vertex_ids[t[1]], vertex_ids[t[2]], vertex_ids[t[3]]);
//    }
//    mesh.initialize_connectivity();
//    mtet::save_mesh("1-tier-grid.msh", mesh);
    /// Re-evaluate the implicit function
    values_3d.reserve(output_vertices_4d.size());
    gradients_3d.reserve(output_vertices_4d.size());
    for (size_t i = 0; i < output_vertices_4d.size(); i++){
        Eigen::RowVector4d vert_4d;
        vert_4d << output_vertices_4d[i][0], output_vertices_4d[i][1], output_vertices_4d[i][2], output_vertices_4d[i][3];
        std::pair<Scalar, Eigen::RowVector4d> valGrad = implicit_sweep(vert_4d);
        values_3d.emplace_back(valGrad.first);
        gradients_3d.emplace_back(std::array<double, 4>{valGrad.second(0), valGrad.second(1), valGrad.second(2), valGrad.second(3)});
    }
    
    
    
    
    /// Second marching
    //std::cout << output_prisms_3d.size() << std::endl;
    marching4D::MarchingTetPrism(output_vertices_4d, output_tets_3d, output_prisms_3d, values_3d, gradients_3d, output_vertices, output_triangles, false);
    std::cout << output_vertices.size() <<" " << output_triangles.size() << " ";
//    std::string debug_file = "debug.json";
//    if (std::filesystem::exists(debug_file.c_str())){
//        std::filesystem::remove(debug_file.c_str());
//    }
//    using json = nlohmann::json;
//    std::ofstream fout(debug_file.c_str(),std::ios::app);
//    json jOut;
//    jOut.push_back(json(verts));
//    jOut.push_back(json(simps));
//    jOut.push_back(json(ulsimp));
//    jOut.push_back(json(llsimp));
//    jOut.push_back(json(values));
//    jOut.push_back(json(output_vertices_4d));
//    jOut.push_back(json(output_tets_3d));
//    jOut.push_back(json(output_prisms_3d));
//    jOut.push_back(json(values_3d));
//    jOut.push_back(json(gradients_3d));
//    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
//    fout.close();
//    std::string debug_output_file = "debug_output.json";
//    if (std::filesystem::exists(debug_output_file.c_str())){
//        std::filesystem::remove(debug_output_file.c_str());
//    }
//    using json = nlohmann::json;
//    std::ofstream fout_output(debug_output_file.c_str(),std::ios::app);
//    json jOut_output;
//    jOut_output.push_back(json(output_vertices));
//    jOut_output.push_back(json(output_triangles));
//    fout_output << jOut_output.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
//    fout_output.close();
//    for (size_t i = 0; i < timer_amount; i++){
//        std::cout
//        << time_label[i] << " "
//        << profileTimer[i] << " ";
//    }
//    std::cout << ms << " ";
//
////  saving 4D grid
//    if (!save_4d_grid("grid4D.json", verts, simps, ulsimp, llsimp)){
//        throw std::runtime_error ("Error: save 4D grid failed.");
//    }
//  saving surface meshes
    if (!save_surface_mesh("mesh.json", output_vertices, output_triangles)){
        throw std::runtime_error ("Error: save surface mesh failed.");
    }

    return 0;
}
