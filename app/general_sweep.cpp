#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <Marching4D.h>
#include <mtet/io.h>
#include <implicit_functions.h>
#include <chrono>
#include <igl/read_triangle_mesh.h>
#include <igl/signed_distance.h>
#include <mtetcol/contour.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/io.h>
//#include <arrangement/Arrangement.h>
//#include <orient_triangles.h>
//#include <igl/write_triangle_mesh.h>¯

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
    ///libigl input using mesh files (unstable gradients, need high resolution mesh input):
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
        value = s*sqrt(sqrd);
        Eigen::RowVector3d cp = c - pos;
        cp.normalize();
        //std::cout << cp << std::endl;
        xyz_grad  = (-s) * cp * Rt.inverse();
        gradient << xyz_grad;
        //std::cout << xyz_grad << std::endl;
        point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
        gradient(3) =  (-s) * cp.dot(point_velocity);
        //std::cout << s * cp.dot(point_velocity) << std::endl;
        return {value, gradient};
    };
    
    std::function<Scalar(Eigen::RowVector4d)> libiglVal = [&](Eigen::RowVector4d data)->Scalar
    {
        Scalar value;
        Eigen::RowVector4d gradient;
        Eigen::RowVector3d P = data.head(3);
        double t = data[3];
        Eigen::RowVector3d running_closest_point = V.row(0);
        double running_sign = 1.0;
        int i;
        double s,sqrd,sqrd2,s2;
        Eigen::Matrix3d VRt,Rt;
        Eigen::RowVector3d xt,vt;
        trajLine3D(t, xt, vt);
        trajLineRot3D(t, Rt, VRt, rotation);
        
        Eigen::RowVector3d pos, c;
        pos = ((Rt.inverse())*((P - xt).transpose())).transpose();
        // fast winding number
        Eigen::VectorXd w;
        igl::fast_winding_number(fwn_bvh,2.0,pos,w);
        s = 1.-2.*w(0);
        sqrd = tree.squared_distance(V,F,pos,i,c);
        value = s*sqrt(sqrd);
        return value;
    };
    
    const double delta = 0.0000001;
    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> libiglNumFunc = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
    {
        Scalar value = libiglVal(data);
        Eigen::RowVector4d grad;
        grad << (libiglVal(data + Eigen::RowVector4d{delta, 0, 0, 0}) - libiglVal(data - Eigen::RowVector4d{delta, 0, 0, 0})) / (2 * delta),
        (libiglVal(data + Eigen::RowVector4d{0, delta, 0, 0}) - libiglVal(data - Eigen::RowVector4d{0, delta, 0, 0})) / (2 * delta),
        (libiglVal(data + Eigen::RowVector4d{0, 0, delta, 0}) - libiglVal(data - Eigen::RowVector4d{0, 0, delta, 0})) / (2 * delta),
        (libiglVal(data + Eigen::RowVector4d{0, 0, 0, delta}) - libiglVal(data - Eigen::RowVector4d{0, 0, 0, delta})) / (2 * delta);
        return {value, grad};
    };
    
    auto implicit_sweep = [&](Eigen:: RowVector4d data){
        return flippingDonutFullTurn(data);
    };
    

    ///
    ///
    ///Grid generation:
    vertExtrude vertexMap;
    insidenessMap insideMap;
    std::array<double, timer_amount> profileTimer = {0, 0, 0, 0, 0, 0, 0, 0};
    auto starterTime = std::chrono::high_resolution_clock::now();
    if (!gridRefine(grid, vertexMap, insideMap, implicit_sweep, threshold, max_splits, profileTimer)){
        throw std::runtime_error("ERROR: grid generation failed");
        return 0;
    };

    /// save the grid output for discretization tool
    save_mesh_json("grid.json", grid);
    /// write grid and active tets
    mtet::save_mesh("tet_grid.msh", grid);
    
    
    Scalar iso_value = 0.0;
    bool cyclic = false;
    std::vector<mtetcol::Scalar> verts;
    std::vector<mtetcol::Index> simps;
    std::vector<std::vector<double>> time;
    std::vector<std::vector<double>> values;
    convert_4d_grid_mtetcol(grid, vertexMap,
                            verts,
                            simps,
                            time,
                            values,
                            cyclic);
    auto stopperTime = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    auto duration = end - start;
    double ms = duration * 0.001;
    std::cout << "Grid-building time: " << ms * 0.001 << std::endl;
    std::function<std::span<double>(size_t)> time_func = [&](size_t index)->std::span<double>{
        return time[index];
    };
    std::function<std::span<double>(size_t)> values_func = [&](size_t index)->std::span<double>{
        return values[index];
    };
    mtetcol::SimplicialColumn<4> columns;
    columns.set_vertices(verts);
    columns.set_simplices(simps);
    columns.set_time_samples(time_func, values_func);
    
    auto contour = columns.extract_contour(iso_value, cyclic);
    
    // Triangulate cycles if needed
    contour.triangulate_cycles();
    
    size_t num_contour_vertices = contour.get_num_vertices();
    std::vector<double> function_values(num_contour_vertices);
    for (size_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        function_values[i] = implicit_sweep(Eigen::RowVector4d{pos[0], pos[1], pos[2], pos[3]}).first;
    }
    // Extract isocontour
    auto isocontour = contour.isocontour(function_values);
    isocontour.triangulate_cycles();
    stopperTime = std::chrono::high_resolution_clock::now();
    end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    duration = end - start;
    ms = duration * 0.001;
    std::cout << "Total time: " << ms * 0.001 << std::endl;
    mtetcol::save_contour("contour.msh", isocontour);
    
    
//    /// Mathematica isosurfacing output:
//    std::vector<std::array<double, 3>> verts_math;
//    std::vector<std::array<size_t, 4>> simps_math;
//    std::vector<std::vector<double>> time_math;
//    std::vector<std::vector<double>> values_math;
//    convert_4d_grid_col(grid, vertexMap,
//                        verts_math,
//                        simps_math,
//                        time_math,
//                        values_math);
//    std::string column_iso_file = "column_iso.json";
//    if (std::filesystem::exists(column_iso_file.c_str())){
//        std::filesystem::remove(column_iso_file.c_str());
//    }
//    using json = nlohmann::json;
//    std::ofstream fout(column_iso_file.c_str(),std::ios::app);
//    json jOut;
//    jOut.push_back(json(verts_math));
//    jOut.push_back(json(simps_math));
//    jOut.push_back(json(time_math));
//    jOut.push_back(json(values_math));
//    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
//    fout.close();
//    /// End of Mathematica output
    
/*
//    arrangement::MatrixFr vertices;    // nx3 Vertex matrix
//    arrangement::MatrixIr faces;       // mx3 Face matrix
//    arrangement::VectorI face_labels;  // mx1 Face labels
//    for (size_t i = 0; i < output_vertices.size(); i++){
//        Eigen::Vector3d coord;
//        coord << output_vertices[i][0], output_vertices[i][1], output_vertices[i][2];
//        vertices << coord;
//    }
//    for (size_t i = 0; i < output_triangles.size(); i++){
//        Eigen::Vector3i faceid;
//        faceid << output_triangles[i][0], output_triangles[i][1], output_triangles[i][2];
//        faces << faceid;
//    }
//    //igl::read_triangle_mesh("../data/three_cubes_randOriented.obj", vertices, faces);
//    face_labels = Eigen::VectorXi::LinSpaced(faces.rows(), 0, faces.rows() - 1);
//
//    // convert faces to vector of arrays
//    std::vector<std::array<size_t, 3>> triangles(faces.rows());
//    for (size_t i = 0; i < faces.rows(); i++) {
//        triangles[i][0] = faces(i, 0);
//        triangles[i][1] = faces(i, 1);
//        triangles[i][2] = faces(i, 2);
//    }
//
//    // orient the triangles
//    std::vector<bool> kept;
//    sweepLabeling::OrientManifoldTriangles(triangles, kept);
//
//    // convert the triangles back to faces
//    for (size_t i = 0; i < triangles.size(); i++) {
//        faces(i, 0) = triangles[i][0];
//        faces(i, 1) = triangles[i][1];
//        faces(i, 2) = triangles[i][2];
//    }
//
//    // create a mesh arrangement engine.
//    auto engine = arrangement::Arrangement::create_mesh_arrangement(vertices, faces, face_labels);
//
//    engine->run();
//
//    const auto& out_vertices = engine->get_vertices();
//    auto out_faces = engine->get_faces();
//
//    // parent face index for each output face
//    const auto& parent_ids = engine->get_out_face_labels();
//
//    // compute cells for each output faces and recover face orientation
//    const auto& patches = engine->get_patches();
//    const auto& patch_cells = engine->get_cells();
//    std::vector<std::array<size_t,2>> cells(out_faces.rows());
//    for (size_t i=0; i < out_faces.rows(); i++) {
//        auto parent_id = parent_ids(i);
//        auto patch_id = patches(i);
//        if (kept[parent_id]) {
//            cells[i][0] = patch_cells(patch_id,0);
//            cells[i][1] = patch_cells(patch_id,1);
//        } else {
//            cells[i][0] = patch_cells(patch_id,1);
//            cells[i][1] = patch_cells(patch_id,0);
//            std::swap(out_faces(i,0), out_faces(i,1));
//        }
//    }
//
//    // write the mesh to a file
//    igl::write_triangle_mesh("output.obj", out_vertices, out_faces);
    
    */
    return 0;
}
