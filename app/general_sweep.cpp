#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <mtet/io.h>
#include <chrono>
#include <igl/read_triangle_mesh.h>
#include <igl/signed_distance.h>
#include <mtetcol/contour.h>
#include <mtetcol/simplicial_column.h>
#include <mtetcol/io.h>
#include <igl/write_triangle_mesh.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include "init_grid.h"
#include "io.h"
#include "col_gridgen.h"
#include "trajectory.h"
#include "post_processing.h"
#include "timer.h"

#define SAVE_CONTOUR
//#define batch_stats
int main(int argc, const char *argv[])
{
    struct
    {
        std::string grid_file;
        std::string output_path;
        std::string function_file = "";
        double threshold = 0.0005;
        double traj_threshold = 0.01;
        int max_splits = std::numeric_limits<int>::max();
        int rot = 0;
        bool without_snapping = false;
        bool without_opt_triangulation = false;
    } args;
    CLI::App app{"Generalized Swept Volume"};
    app.add_option("grid", args.grid_file, "Initial grid file")->required();
    app.add_option("output", args.output_path, "Output path")->required();
    app.add_option("-f,--function", args.function_file, "Implicit function file");
    app.add_option("-t,--threshold", args.threshold, "Threshold value");
    app.add_option("--tt, --traj_threshold", args.traj_threshold, "Threshold value for trajectory");
    app.add_option("-m,--max-splits", args.max_splits, "Maximum number of splits");
    app.add_option("-r,--rotation-number", args.rot, "Number of rotations");
    app.add_flag("--without-snapping", args.without_snapping, "Disable vertex snapping in iso-surfacing step");
    app.add_flag("--without-optimal-triangulation", args.without_opt_triangulation,
                 "Disable optimal triangulation in iso-surfacing triangulation step");
    
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
    std::string output_path = args.output_path;
    int max_splits = args.max_splits;
    std::string function_file = args.function_file;
    double threshold = args.threshold;
    double traj_threshold = args.traj_threshold;
    int rotation = args.rot;
    const int dim = 4;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::AABB<Eigen::MatrixXd,3> tree;
    igl::FastWindingNumberBVH fwn_bvh;
    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> implicit_sweep;
    if (args.function_file != ""){
        if (args.function_file.find(".obj") != std::string::npos){
            igl::read_triangle_mesh(args.function_file,V,F);
            tree.init(V,F);
            int order = 2;
            igl::fast_winding_number(V, F, order, fwn_bvh);\
            igl::WindingNumberAABB<double, int> hier;
            hier.set_mesh(V,F);
            hier.grow();
            /// main function:
            /// the lambda function for function evaluations
            ///  @param[in] data            The 4D coordinate
            ///  @return    A std::pari<Scalar, Eigen::RowVector4d> of the value and the gradients at this 4D point
            ///
            ///libigl input using mesh files (unstable gradients, need high resolution mesh input):
            implicit_sweep = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
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
                xyz_grad  = (-s) * cp * Rt.inverse();
                gradient.template head<3>() << xyz_grad;
                point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
                gradient(3) =  (-s) * cp.dot(point_velocity);
                std::cout << s * cp.dot(point_velocity) << std::endl;
                std::cout << "value: " << value << std::endl;
                return {value, gradient};
            };
        } else if (args.function_file == "elbow") {
            implicit_sweep = elbow;
        } else if (args.function_file == "bezier") {
            implicit_sweep = bezier;
        } else if (args.function_file == "blend_spheres") {
            implicit_sweep = blend_spheres;
        } else if (args.function_file == "blend_sphere_torus") {
            implicit_sweep = blend_sphere_torus;
        } else if (args.function_file == "sphere_spiral") {
            implicit_sweep = sphere_spiral;
        } else if (args.function_file == "knot") {
            implicit_sweep = knot;
        } else if (args.function_file == "brush_stroke") {
            implicit_sweep = brush_stroke;
        } else if (args.function_file == "brush_stroke_blending") {
            implicit_sweep = brush_stroke_blending;
        } else if (args.function_file == "concentric_rings") {
            implicit_sweep = concentric_rings;
        } else if (args.function_file == "spinning_rod") {
            implicit_sweep = spinning_rod;
        } else if (args.function_file == "letter_L") {
            implicit_sweep = letter_L;
        } else if (args.function_file == "letter_L_blend") {
            implicit_sweep = letter_L_blend;
        } else if (args.function_file == "torus_rotation") {
            implicit_sweep = torus_rotation;
        } else if (args.function_file == "loopDloop_with_offset") {
            implicit_sweep = loopDloop_with_offset;
        } else if (args.function_file == "loopDloop_with_offset_v2") {
            implicit_sweep = loopDloop_with_offset_v2;
        } else if (args.function_file == "doghead") {
            implicit_sweep = doghead;
        } else if (args.function_file == "star_S") {
            implicit_sweep = star_S;
        } else if (args.function_file == "star_D") {
            implicit_sweep = star_D;
        } else if (args.function_file == "star_F") {
            implicit_sweep = star_F;
        } else if (args.function_file == "star_I") {
            implicit_sweep = star_I;
        } else if (args.function_file == "fertility") {
            implicit_sweep = fertility;
        } else if (args.function_file == "fertility_v2") {
            implicit_sweep = fertility_v2;
        } else if (args.function_file == "fertility_v3") {
            implicit_sweep = fertility_v3;
        } else if (args.function_file == "fertility_v4") {
            implicit_sweep = fertility_v4;
        } else if (args.function_file == "fertility_v5") {
            implicit_sweep = fertility_v5;
        } else if (args.function_file == "fertility_v6") {
            implicit_sweep = fertility_v6;
        } else if (args.function_file == "bunny_blend") {
            implicit_sweep = bunny_blend;
        } else if (args.function_file == "loopDloop_with_offset_v3") {
            implicit_sweep = loopDloop_with_offset_v3;
        } else if (args.function_file == "VIPSS_blend") {
            implicit_sweep = VIPSS_blend;
        } else {
            throw std::runtime_error("ERROR: file format not supported");
        }
    }else{
        /// use hard coded models as default/testing purpose.
        implicit_sweep = [&](Eigen:: RowVector4d data){
            return flippingDonutFullTurn(data);
        };
    }
    ///
    ///
    ///Grid generation:
    vertExtrude vertexMap;
    insidenessMap insideMap;
    std::cout << "Start to generate the background grid..." << std::endl;
    std::array<double, timer_amount> profileTimer = {0, 0, 0, 0, 0, 0, 0, 0};
    auto starterTime = std::chrono::high_resolution_clock::now();
    spdlog::set_level(spdlog::level::off);
    if (!gridRefine(grid, vertexMap, insideMap, implicit_sweep, threshold, traj_threshold, max_splits, profileTimer)){
        throw std::runtime_error("ERROR: grid generation failed");
        return 0;
    };
    spdlog::set_level(spdlog::level::info);
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
    auto grid_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Grid-building time: " << (grid_end - start) * 0.000001 << std::endl;
    std::function<std::span<double>(size_t)> time_func = [&](size_t index)->std::span<double>{
        return time[index];
    };
    std::function<std::span<double>(size_t)> values_func = [&](size_t index)->std::span<double>{
        return values[index];
    };
    mtetcol::SimplicialColumn<dim> columns;
    columns.set_vertices(verts);
    columns.set_simplices(simps);
    columns.set_time_samples(time_func, values_func);
    
    auto contour = columns.extract_contour(iso_value, cyclic);
    stopperTime = std::chrono::high_resolution_clock::now();
    auto surface_1_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Surfacing time: " << (surface_1_end - grid_end) * 0.000001 << std::endl;
    size_t num_contour_vertices = contour.get_num_vertices();
    std::vector<double> function_values(num_contour_vertices);
    std::vector<double> gradient_values(num_contour_vertices * dim);
    for (uint32_t i = 0; i < num_contour_vertices; ++i) {
        auto pos = contour.get_vertex(i);
        auto pos_eval = implicit_sweep(Eigen::RowVector4d{pos[0], pos[1], pos[2], pos[3]});
        function_values[i] = pos_eval.first;
        gradient_values[dim * i] = pos_eval.second[0];
        gradient_values[dim * i + 1] = pos_eval.second[1];
        gradient_values[dim * i + 2] = pos_eval.second[2];
        gradient_values[dim * i + 3] = pos_eval.second[3];
    }
    // Extract isocontour
    auto isocontour = contour.isocontour(function_values, gradient_values, !args.without_snapping);
    isocontour.triangulate_cycles(!args.without_opt_triangulation);
    if (!std::filesystem::exists(output_path)) {
        // Attempt to create the directory
        if (std::filesystem::create_directory(output_path)) {
            std::cout << "Directory created successfully." << std::endl;
        } else {
            std::cerr << "Failed to create directory." << std::endl;
        }
    } else {
        std::cout << "Directory already exists. Removing all contents..." << std::endl;
        for (const auto& entry : std::filesystem::directory_iterator(output_path)) {
            std::error_code ec;
            std::filesystem::remove_all(entry.path(), ec);
            if (ec) {
                std::cerr << "Error removing " << entry.path() << ": " << ec.message() << std::endl;
            }
        }
    }
    stopperTime = std::chrono::high_resolution_clock::now();
    auto surface_2_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Surfacing time: " << (surface_2_end - surface_1_end) * 0.000001 << std::endl;
    
#ifdef SAVE_CONTOUR
    mtetcol::save_contour(output_path + "/temporal_grid.obj", contour);
    mtetcol::save_contour(output_path + "/contour.msh", isocontour);
    
    /// Mathematica isosurfacing output:
    std::vector<std::array<double, 3>> verts_math;
    std::vector<std::array<size_t, 4>> simps_math;
    std::vector<std::vector<double>> time_math;
    std::vector<std::vector<double>> values_math;
    convert_4d_grid_col(grid, vertexMap,
                        verts_math,
                        simps_math,
                        time_math,
                        values_math);
    std::string column_iso_file = "column_iso.json";
    {
        if (std::filesystem::exists(column_iso_file.c_str())){
            std::filesystem::remove(column_iso_file.c_str());
        }
        using json = nlohmann::json;
        std::ofstream fout(column_iso_file.c_str(),std::ios::app);
        json jOut;
        jOut.push_back(json(verts_math));
        jOut.push_back(json(simps_math));
        jOut.push_back(json(time_math));
        jOut.push_back(json(values_math));
        fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
        fout.close();
    }
    /// End of Mathematica output
#endif
    arrangement::MatrixFr vertices;    // nx3 Vertex matrix
    arrangement::MatrixIr faces;       // mx3 Face matrix
    size_t num_vertices = isocontour.get_num_vertices();
    size_t num_triangles = isocontour.get_num_cycles();
    vertices.resize(num_vertices, 3);
    faces.resize(num_triangles, 3);
    for (size_t i = 0; i < num_vertices; i++){
        Eigen::RowVector3d coord;
        auto pos = isocontour.get_vertex(i);
        vertices(i,0) = pos[0];
        vertices(i,1) = pos[1];
        vertices(i,2) = pos[2];
    }
    for (size_t i = 0; i < num_triangles; i++){
        size_t ind = 0;
        auto tris = isocontour.get_cycle(i);
        assert(tris.size() == 3);
        for (auto si : tris) {
            mtetcol::Index seg_id = index(si);
            bool seg_ori = mtetcol::orientation(si);
            auto seg = isocontour.get_segment(seg_id);
            faces(i, ind) = (seg_ori ? seg[0] : seg[1]);
            ind ++;
        }
    }
    arrangement::MatrixFr out_vertices;  // Output vertex matrix
    std::vector<arrangement::MatrixIr> out_faces;     // Output face matrix
    compute_sweep_volume(vertices, faces, out_vertices, out_faces, output_path);
    stopperTime = std::chrono::high_resolution_clock::now();
    auto arrangement_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Arrangement time: " << (arrangement_end - surface_2_end) * 0.000001 << std::endl;
    igl::write_triangle_mesh(output_path + "/mesh" + ".obj", vertices, faces);
        for (size_t i = 0; i < out_faces.size(); i++){
            igl::write_triangle_mesh(output_path + "/" + std::to_string(i) + ".obj", out_vertices, out_faces[i]);
        }
    for (size_t i = 0; i < out_faces.size(); i++){
        igl::write_triangle_mesh(output_path + "/" + std::to_string(i) + ".obj", out_vertices, out_faces[i]);
    }
#ifdef batch_stats
    std::string stats_file = "stats.json";
    {
        std::vector<std::array<double, 3>> verts_math;
        std::vector<std::array<size_t, 4>> simps_math;
        std::vector<std::vector<double>> time_math;
        std::vector<std::vector<double>> values_math;
        convert_4d_grid_col(grid, vertexMap,
                            verts_math,
                            simps_math,
                            time_math,
                            values_math);
        int num_non_tet_poly = 0;
        for (uint32_t i = 0; i < contour.get_num_polyhedra(); ++i) {
            if (!contour.is_polyhedron_regular(i)) num_non_tet_poly++;
        }
        int after_tris = 0;
        for (size_t i = 0; i < out_faces.size(); i++) after_tris += out_faces[i].rows();
        using json = nlohmann::json;
        std::ofstream fout(stats_file.c_str(),std::ios::app);
        json jOut;
        jOut["grid tets"] = grid.get_num_tets();
        jOut["contour cells"] = contour.get_num_polyhedra();
        jOut["non_tet polys"] = num_non_tet_poly;
        jOut["tris before arrangement"] = num_triangles;
        jOut["tris after arrangement"] = after_tris;
        fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
        fout.close();
    }
#endif
    return 0;
}
