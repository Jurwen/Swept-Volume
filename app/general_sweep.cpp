#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <mtet/io.h>
#include <implicit_functions.h>

#include "init_grid.h"
#include "io.h"
#include "col_gridgen.h"
#include "trajectory.h"

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
    } else {
        grid = mtet::load_mesh(args.grid_file);
    }
    
    int max_splits = args.max_splits;
    std::string function_file = args.function_file;
    double threshold = args.threshold;
    
    /// main function:
    
    /// Read implicit function
    std::vector<std::unique_ptr<ImplicitFunction<double>>> functions;
    load_functions(function_file, functions);
    std::unique_ptr<ImplicitFunction<double>> &object = functions[0];
    auto trajFunc = trajLine;
    
    ///
    /// the lambda function for function evaluations
    ///  @param[in] data            The 3D coordinate
    ///  @param[in] funcNum         The number of functions
    ///
    ///  @return `Eigen::RowVector4d`. It represents the value at 0th index and gradients at {1, 2, 3} index.
    auto implicit_sweep = [&](std::span<const Scalar, 4> data){
        Eigen::RowVector3d traversed = trajFunc(data[3], data.first<3>());
        Eigen::RowVector4d eval;
        eval[0] = object->evaluate_gradient(traversed[0], traversed[1], traversed[2], eval[1], eval[2], eval[3]);
        return Eigen::RowVector4d{-1 * eval[0], -1 * eval[1], -1 * eval[2], -1 * eval[3]};
    };
    
    ///
    ///
    ///Grid generation:
    vertexCol timeMap;
    tetCol cell5Map;
    if (!gridRefine(grid, timeMap, cell5Map, implicit_sweep, threshold)){
        throw std::runtime_error("ERROR: grid generation failed");
        return 0;
    };
    
    /// save the grid output for discretization tool
    save_mesh_json("grid.json", grid);
    /// write grid and active tets
    mtet::save_mesh("tet_grid.msh", grid);
    
    std::cout << "saving 4D grid..." << std::endl;
    if (!save_4d_grid("grid4D.json",grid,timeMap,cell5Map)){
        throw std::runtime_error ("Error: save 4D grid failed.");
    }

    return 0;
}
