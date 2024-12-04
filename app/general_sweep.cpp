#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>
#include <mtet/io.h>

#include "init_grid.h"
#include "io.h"

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
    ///
    vertexCol timeMap;
    tetCol cell5Map;
    init_grid::init5CGrid(3, grid, 1024, timeMap, cell5Map);
    
    /// save the grid output for discretization tool
    save_mesh_json("grid.json", grid);
    /// write grid and active tets
    mtet::save_mesh("tet_grid.msh", grid);
    if (!save_4d_grid("grid4D.json",grid,timeMap,cell5Map)){
        std::cout << "Error: save 4D grid failed." << std::endl;
    }

    return 0;
}
