#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <mtet/io.h>
#include <implicit_functions.h>
#include <chrono>

#include "init_grid.h"
#include "io.h"
#include "col_gridgen.h"
#include "trajectory.h"
#include "ref_crit.h"
#include "timer.h"

int main(int argc, const char *argv[])
{
    struct
    {
        std::string grid_file;
        std::string function_file;
        double threshold;
        int max_splits = std::numeric_limits<int>::max();;
    } args;
    CLI::App app{"Generalized Swept Volume"};
    app.add_option("grid", args.grid_file, "Initial grid file")->required();
    app.add_option("function", args.function_file, "Implicit function file")->required();
    app.add_option("-t,--threshold", args.threshold, "Threshold value");
    app.add_option("-m,--max-splits", args.max_splits, "Maximum number of splits");
    
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
    
    /// main function:
    
    /// Read implicit function
//    std::vector<std::unique_ptr<ImplicitFunction<double>>> functions;
//    load_functions(function_file, functions);
//    std::unique_ptr<ImplicitFunction<double>> &object = functions[0];
//    auto trajFunc = trajLine;
    
    ///
    /// the lambda function for function evaluations
    ///  @param[in] data            The 4D coordinate
    ///  @return    A std::pari<Scalar, Eigen::RowVector4d> of the value and the gradients at this 4D point
    
    auto implicit_sweep = [&](Eigen::RowVector4d data){
        return flippingDonut(data);
    };
    ///
    ///
    ///Grid generation:
    vertExtrude vertexMap;
    tetExtrude cell5Map;
    std::array<double, timer_amount> profileTimer = {0,0};
    auto starterTime = std::chrono::high_resolution_clock::now();
    if (!gridRefine(grid, vertexMap, cell5Map, implicit_sweep, threshold, max_splits, profileTimer)){
        throw std::runtime_error("ERROR: grid generation failed");
        return 0;
    };
    std::cout << profileTimer[0] << " " << profileTimer[1] << " ";
    auto stopperTime = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(starterTime).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(stopperTime).time_since_epoch().count();
    auto duration = end - start;
    double ms = duration * 0.001;
    std::cout << ms << " ";

    /// save the grid output for discretization tool
    save_mesh_json("grid.json", grid);
    /// write grid and active tets
    mtet::save_mesh("tet_grid.msh", grid);

//    std::cout << "saving 4D grid..." << std::endl;
    if (!save_4d_grid("grid4D.json",grid,vertexMap,cell5Map)){
        throw std::runtime_error ("Error: save 4D grid failed.");
    }

    return 0;
}
