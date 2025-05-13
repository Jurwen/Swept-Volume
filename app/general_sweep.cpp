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
#include <arrangement/Arrangement.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <numeric>   // for std::iota
#include <utility>   // for std::pair
#include <limits>    // for std::numeric_limits

#include "init_grid.h"
#include "io.h"
#include "col_gridgen.h"
#include "trajectory.h"
#include "timer.h"




struct pairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        auto a = static_cast<uint32_t>(p.first);
        auto b = static_cast<uint32_t>(p.second);
        uint64_t high = static_cast<uint64_t>(a) << 32;
        uint64_t low  = static_cast<uint64_t>(b);
        uint64_t packed = high | low;
        return static_cast<std::size_t>(packed);
    }
};


// Compute the “validPatch” list as in your Mathematica code.
// - cellNum: number of cells
// - volInfo0: length‑cellNum vector of volumes
// - cellData: list of {u,v} zero‑based adjacency pairs
//
// Returns: a vector of indices i into cellData for which
//          hash[min(u,v)][max(u,v)] was never set to true.
std::pair<std::vector<bool>, std::vector<std::vector<bool>>>
computeValidPatch(size_t cellNum,
                  const std::vector<double>& volInfo0,
                  const arrangement::MatrixIr& cellData)
{
    std::vector<std::vector<bool>> adj(cellNum, std::vector<bool>(cellNum, false));
    std::vector<bool> valid(cellNum, false), merged(cellNum, true);
    std::vector<std::vector<bool>> hash(cellNum, std::vector<bool>(cellNum, false));

    //    We want, for each index i, a “rank” so that largest volInfo0 → rank 0,
    //    next largest → rank 1, etc.
    std::vector<int> order(cellNum);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b){
                  return volInfo0[a] > volInfo0[b];   // descending
              });
    std::vector<int> rankDesc(cellNum);
    for(int pos = 0; pos < cellNum; ++pos) {
        rankDesc[ order[pos] ] = pos;
    }
    const double threshold = 1e-5;
    for(int i = 0; i < cellNum; ++i) {
        if (volInfo0[i] > threshold) valid[i] = true;
    }
    // build adjacency from cellData
    for( size_t e = 0; e < cellData.rows(); e++) {
        int u = cellData(e, 0), v = cellData(e, 1);
        if (u>=0 && u<cellNum && v>=0 && v<cellNum) {
            adj[u][v] = adj[v][u] = true;
        }
    }

    // handle = Xnor/@Transpose[{valid,merged}]
    auto allHandleTrue = [&](){
        for(int i=0;i<cellNum;++i)
            if (valid[i] != merged[i]) return false;
        return true;
    };

    // While[Not@AllTrue@handle, … do merge steps …]
    while(!allHandleTrue()) {
        int target = 0;
        while(target < cellNum && valid[target] == merged[target])
            ++target;
        // find its neighbors: mergeCell = Flatten@Position[adj[[target]],True]
        std::vector<int> neigh;
        neigh.reserve(cellNum);
        for(int j=0; j<cellNum; ++j) {
            if (adj[target][j]) neigh.push_back(j);
        }
        if (neigh.empty()) break;  // nothing to merge with?

        auto it = std::min_element(
            neigh.begin(), neigh.end(),
            [&](int a,int b){ return rankDesc[a] < rankDesc[b]; }
        );
        int mergeCell = *it;

        // union adjacency rows & columns:
        for(int k=0;k<cellNum;++k) {
            adj[mergeCell][k] = adj[mergeCell][k] || adj[target][k];
            adj[k][mergeCell] = adj[k][mergeCell] || adj[k][target];
        }

        hash[target][mergeCell] = true;
        hash[mergeCell][target] = true;

        merged[target] = false;
    }

    // 9) Build validPatch: those edges {u,v} for which hash[min,max]==false
    std::vector<bool> validPatch(cellData.rows(), false);
    for(int i = 0; i < cellData.rows(); ++i) {
        int u = cellData(i, 0),
            v = cellData(i, 1);
        if (! hash[u][v]) {
            validPatch[i] = true;
        }
    }
    return {validPatch, hash};
}


/**
 * Compute the surface of the sweep volume.
 *
 * @param vertices Input vertices of the sweep generator.
 * @param faces Input faces of the sweep generator.
 * @param out_vertices Output vertices of the sweep volume surface.
 * @param out_faces Output faces of the sweep volume surface.
 */
void compute_sweep_volume(const arrangement::MatrixFr& vertices, const arrangement::MatrixIr& faces,
                          arrangement::MatrixFr& out_vertices, std::vector<arrangement::MatrixIr>& out_faces, std::string output_path) {
    // assume the input faces have already been oriented based on second order time derivatives
    // then the sweep volume consists of arrangement cells with positive winding number
    
    // Initialize face labels
    arrangement::VectorI face_labels = Eigen::VectorXi::LinSpaced(faces.rows(), 0, faces.rows() - 1);
    // create a mesh arrangement engine.
    auto engine = arrangement::Arrangement::create_mesh_arrangement(vertices, faces, face_labels);
    // Alternatively, create a fast mesh arrangement engine.
    // face_labels.setZero();  // fast arrangement cannot use too large face labels
    // auto engine = arrangement::Arrangement::create_fast_arrangement(vertices, faces, face_labels);
    // or, create a Geogram arrangement engine
    // the cell ids associated with patches can often be MAX_INT64.
    // auto engine = arrangement::Arrangement::create_geogram_arrangement(vertices, faces, face_labels);
    
// The following is based on James' arrangement code in Python: https://github.com/qnzhou/arrangement-benchmark/blob/main/python/arrangement/__main__.py
    engine->run();
    const auto& cell_data       = engine->get_cells();           // (#facets x 2) array
    const auto& patches         = engine->get_patches();         // list of facet indices
    const auto& parent_facets   = engine->get_out_face_labels();   // size = #facets
    const auto& winding_number  = engine->get_winding_number();  // (#facets x 2) array
    out_vertices = engine->get_vertices();
    const auto& arrangement_faces = engine->get_faces();
    size_t num_cells = engine->get_num_cells();
    size_t num_patches = engine->get_num_patches();
    size_t num_facets = arrangement_faces.rows();
    std::cout << "Printing arrangement stats: " << std::endl;
    std::cout << "num_cells: " << num_cells << std::endl;
    std::cout << "num_patches: " << num_patches << std::endl;
    std::cout << "num_vertices: " << engine->get_vertices().rows() << std::endl;
    std::cout << "num_faces: " << arrangement_faces.rows() << std::endl;
    std::vector<int> wind_list(num_cells);
    std::vector<double> volInfo(num_cells);
    std::vector<size_t> cellIt;
    for (size_t i = 0; i < num_cells; ++i) {
        // collect all facets incident on cell i
        std::vector<int> active_facets;
        active_facets.reserve(num_patches);
        int active_facets_count = 0;
        bool active_orientations;
        double meshVol = 0.0;
        for (int facet_idx = 0; facet_idx < patches.size(); ++facet_idx) {
            int c0 = cell_data(patches[facet_idx], 0);
            int c1 = cell_data(patches[facet_idx], 1);
            if (c0 == i || c1 == i) {
                active_facets.emplace_back(facet_idx);
                active_facets_count++;
                int j = arrangement_faces(facet_idx, 0), k = arrangement_faces(facet_idx, 1), l = arrangement_faces(facet_idx, 2);
                Eigen::Vector3d vj = out_vertices.row(j);
                Eigen::Vector3d vk = out_vertices.row(k);
                Eigen::Vector3d vl = out_vertices.row(l);
                meshVol += vj.dot( vk.cross(vl) );
                if (active_facets_count == 1){
                    active_orientations = (c0 == i);
                }
            }
        }
        meshVol = std::abs(meshVol) / 6.0;
        volInfo[i] = meshVol;
        active_facets.resize(active_facets_count);
        if (active_facets.empty()) {
            continue;
        }
        
        // pick winding from the first active facet, based on its orientation
        int active_wind;
        if (!active_orientations) {
            active_wind = winding_number(active_facets[0], 1);
        } else {
            active_wind = winding_number(active_facets[0], 0);
        }
        wind_list[i] = active_wind;
    }
    std::vector<bool> valid(num_cells, false);
    const double threshold = 1e-5;
    for(int i = 0; i < num_cells; ++i) {
        if (volInfo[i] > threshold) valid[i] = true;
        if (valid[i] && wind_list[i] == 0){
            cellIt.push_back(i);
            //std::cout << i << " ";
        }
    }
    //std::cout << std::endl;
    
    // Pruning of error cells
    const auto& prunedInfo = computeValidPatch(num_cells, volInfo, cell_data);
    const std::vector<std::vector<bool>>& hash = prunedInfo.second;
    const std::vector<bool>& valid_patchInd = prunedInfo.first;

    // Check edge valences for finding the feature lines and corners
    std::vector<std::array<std::array<double, 3>, 2>> feature_lines;
    std::vector<std::array<double, 3>> corners;
    std::unordered_map<std::pair<int, int>, int, pairHash> edge_valence;
    std::unordered_map<int, int> vert_valence;
    auto push_valance = [&](int v0, int v1){
        if (v0 < v1){
            edge_valence[{v0, v1}] ++;
        } else{
            edge_valence[{v1, v0}] ++;
        }
    };
    auto push_vert_valance = [&](int v){
        vert_valence[v] ++;
    };
    for (size_t i = 0; i < arrangement_faces.rows(); i++){
        if (valid_patchInd[patches(i, 0)]){
            const auto& face = arrangement_faces.row(i);
            push_valance(face[0], face[1]);
            push_valance(face[1], face[2]);
            push_valance(face[0], face[2]);
        }
    }
    auto check_valence = [&](int v0, int v1, int v2){
        std::vector<int> vert_list = {v0, v1, v2};
        std::sort(vert_list.begin(), vert_list.end());
        if (edge_valence[{vert_list[0], vert_list[1]}] > 2){
            edge_valence[{vert_list[0], vert_list[1]}] = 2;
            push_vert_valance(vert_list[0]);
            push_vert_valance(vert_list[1]);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(vert_list[0], 0), out_vertices(vert_list[0], 1), out_vertices(vert_list[0], 2)}, std::array<double, 3>{out_vertices(vert_list[1], 0), out_vertices(vert_list[1], 1), out_vertices(vert_list[1], 2)}});
        }
        if (edge_valence[{vert_list[0], vert_list[2]}] > 2){
            edge_valence[{vert_list[0], vert_list[2]}] = 2;
            push_vert_valance(vert_list[0]);
            push_vert_valance(vert_list[2]);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(vert_list[0], 0), out_vertices(vert_list[0], 1), out_vertices(vert_list[0], 2)}, std::array<double, 3>{out_vertices(vert_list[2], 0), out_vertices(vert_list[2], 1), out_vertices(vert_list[2], 2)}});
        }
        if (edge_valence[{vert_list[1], vert_list[2]}] > 2){
            edge_valence[{vert_list[1], vert_list[2]}] = 2;
            push_vert_valance(vert_list[1]);
            push_vert_valance(vert_list[2]);
            feature_lines.push_back(std::array<std::array<double, 3>, 2>{std::array<double, 3>{out_vertices(vert_list[1], 0), out_vertices(vert_list[1], 1), out_vertices(vert_list[1], 2)}, std::array<double, 3>{out_vertices(vert_list[2], 0), out_vertices(vert_list[2], 1), out_vertices(vert_list[2], 2)}});
        }
    };
    std::vector<std::pair<size_t, bool>> patch_list;
    for (size_t i = 0; i < cellIt.size(); ++i){
        auto cluster = hash[cellIt[i]];
        cluster[cellIt[i]] = true;
        arrangement::MatrixIr patch_faces;
        patch_faces.resize(arrangement_faces.rows(), 3);
        size_t face_count = 0;
        for (size_t f = 0; f < num_facets; ++f){
            auto patch_id = patches(f);
            if (cluster[cell_data(patch_id, 0)] != cluster[cell_data(patch_id, 1)]){
                if (cluster[cell_data(patch_id, 0)]){
                    patch_faces.row(face_count) = (arrangement_faces.row(f));
                }else{
                    patch_faces.row(face_count) = (arrangement_faces.row(f).reverse());
                }
                check_valence(patch_faces(face_count, 0), patch_faces(face_count, 1), patch_faces(face_count, 2));
                face_count++;
            }
        }
        patch_faces.conservativeResize(face_count, 3);
        out_faces.emplace_back(patch_faces);
    }
    for (int i = 0; i < out_vertices.rows(); i++){
        if (vert_valence[i] != 2 && vert_valence[i] == 1){
            corners.push_back({out_vertices(i, 0), out_vertices(i, 1), out_vertices(i, 2)});
        }
    }
    std::string feature_lines_file = output_path + "/features.json";
    if (std::filesystem::exists(feature_lines_file.c_str())){
        std::filesystem::remove(feature_lines_file.c_str());
    }
    using json = nlohmann::json;
    std::ofstream fout(feature_lines_file.c_str(),std::ios::app);
    json jOut;
    jOut.push_back(json(feature_lines));
    jOut.push_back(json(corners));
    fout << jOut.dump(4, ' ', true, json::error_handler_t::replace) << std::endl;
    fout.close();
}



int main(int argc, const char *argv[])
{
    struct
    {
        std::string grid_file;
        std::string output_path;
        std::string function_file;
        double threshold;
        int max_splits = std::numeric_limits<int>::max();
        int rot = 0;
    } args;
    CLI::App app{"Generalized Swept Volume"};
    app.add_option("grid", args.grid_file, "Initial grid file")->required();
    app.add_option("output", args.output_path, "Output path")->required();
    app.add_option("-f,--function", args.function_file, "Implicit function file");
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
    std::string output_path = args.output_path;
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
//    Eigen::MatrixXd V;
//    Eigen::MatrixXi F;
//    igl::read_triangle_mesh(args.function_file,V,F);
//    igl::AABB<Eigen::MatrixXd,3> tree;
//    tree.init(V,F);
//    igl::FastWindingNumberBVH fwn_bvh;
//    int order = 2;
//    igl::fast_winding_number(V,F,order,fwn_bvh);
//    igl::WindingNumberAABB<double, int> hier;
//    hier.set_mesh(V,F);
//    hier.grow();
    /// Definition of the implicit function using explicit meshes
//    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> libiglFunc = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
//    {
//        Scalar value;
//        Eigen::RowVector4d gradient;
//        const double iso = 0.001;
//        Eigen::RowVector3d P = data.head(3);
//        double t = data[3];
//        Eigen::RowVector3d running_closest_point = V.row(0);
//        double running_sign = 1.0;
//        int i;
//        double s,sqrd,sqrd2,s2;
//        Eigen::Matrix3d VRt,Rt;
//        Eigen::RowVector3d xt,vt,pos,c,c2,xyz_grad,point_velocity;
//        trajLine3D(t, xt, vt);
//        trajLineRot3D(t, Rt, VRt, rotation);
//
//        pos = ((Rt.inverse())*((P - xt).transpose())).transpose();
//        // fast winding number
//        Eigen::VectorXd w;
//        igl::fast_winding_number(fwn_bvh,2.0,pos,w);
//        s = 1.-2.*w(0);
//        sqrd = tree.squared_distance(V,F,pos,i,c);
//        value = s*sqrt(sqrd);
//        Eigen::RowVector3d cp = c - pos;
//        cp.normalize();
//        //std::cout << cp << std::endl;
//        xyz_grad  = (-s) * cp * Rt.inverse();
//        gradient << xyz_grad;
//        //std::cout << xyz_grad << std::endl;
//        point_velocity = (-Rt.inverse()*VRt*Rt.inverse()*(P.transpose() - xt.transpose()) - Rt.inverse()*vt.transpose()).transpose();
//        gradient(3) =  (-s) * cp.dot(point_velocity);
//        //std::cout << s * cp.dot(point_velocity) << std::endl;
//        return {value, gradient};
//    };
//    
//    std::function<Scalar(Eigen::RowVector4d)> libiglVal = [&](Eigen::RowVector4d data)->Scalar
//    {
//        Scalar value;
//        Eigen::RowVector4d gradient;
//        Eigen::RowVector3d P = data.head(3);
//        double t = data[3];
//        Eigen::RowVector3d running_closest_point = V.row(0);
//        double running_sign = 1.0;
//        int i;
//        double s,sqrd,sqrd2,s2;
//        Eigen::Matrix3d VRt,Rt;
//        Eigen::RowVector3d xt,vt;
//        trajLine3D(t, xt, vt);
//        trajLineRot3D(t, Rt, VRt, rotation);
//        
//        Eigen::RowVector3d pos, c;
//        pos = ((Rt.inverse())*((P - xt).transpose())).transpose();
//        // fast winding number
//        Eigen::VectorXd w;
//        igl::fast_winding_number(fwn_bvh,2.0,pos,w);
//        s = 1.-2.*w(0);
//        sqrd = tree.squared_distance(V,F,pos,i,c);
//        value = s*sqrt(sqrd);
//        return value;
//    };
//    
//    const double delta = 0.0000001;
//    std::function<std::pair<Scalar, Eigen::RowVector4d>(Eigen::RowVector4d)> libiglNumFunc = [&](Eigen::RowVector4d data)->std::pair<Scalar, Eigen::RowVector4d>
//    {
//        Scalar value = libiglVal(data);
//        Eigen::RowVector4d grad;
//        grad << (libiglVal(data + Eigen::RowVector4d{delta, 0, 0, 0}) - libiglVal(data - Eigen::RowVector4d{delta, 0, 0, 0})) / (2 * delta),
//        (libiglVal(data + Eigen::RowVector4d{0, delta, 0, 0}) - libiglVal(data - Eigen::RowVector4d{0, delta, 0, 0})) / (2 * delta),
//        (libiglVal(data + Eigen::RowVector4d{0, 0, delta, 0}) - libiglVal(data - Eigen::RowVector4d{0, 0, delta, 0})) / (2 * delta),
//        (libiglVal(data + Eigen::RowVector4d{0, 0, 0, delta}) - libiglVal(data - Eigen::RowVector4d{0, 0, 0, delta})) / (2 * delta);
//        return {value, grad};
//    };
    
    auto implicit_sweep = [&](Eigen:: RowVector4d data){
        return flippingDonut(data);
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
    auto grid_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Grid-building time: " << (grid_end - start) * 0.000001 << std::endl;
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
//    int tet_num = 0;
//    for (size_t i = 0; i < contour.get_num_polyhedra(); ++i){
//        if (contour.is_polyhedron_regular(i)) tet_num++;
//    }
//    std::cout << contour.get_num_vertices() << " " << contour.get_num_polyhedra() << " " << tet_num << std::endl;
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
    auto surface_end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    std::cout << "Surfacing time: " << (surface_end - grid_end) * 0.000001 << std::endl;
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
    /// End of Mathematica output
    
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
    std::cout << "Arrangement time: " << (arrangement_end - surface_end) * 0.000001 << std::endl;
    igl::write_triangle_mesh(output_path + "/mesh" + ".obj", vertices, faces);
    for (size_t i = 0; i < out_faces.size(); i++){
        igl::write_triangle_mesh(output_path + "/" + std::to_string(i) + ".obj", out_vertices, out_faces[i]);
    }
    
    
    return 0;
}
