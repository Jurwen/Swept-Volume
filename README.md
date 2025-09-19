# Swept-Volume

This code implements the ACM SIGGRAPH ASIA 2025 paper: Lifted Surfacing of Generalized Sweep Volumes

Given any sweep represented as a smooth time-varying implicit function satisfying a genericity assumption, this algorithm produces a watertight and intersection-free surface that faithfully captures the geometric and topological features.

## Build

Use the following command to build: 

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
The program `general_sweep` will be generated in the build file. 

### Dependency

Currently, there are some third-party libraries that are not yet open sourced. We are waiting for some internal authorizations. Please stay tuned for the full release coming soon.

## Usage

To use the `general_sweep` tool, you must provide an initial grid file and output path as required arguments, along with any desired options. The trajectory functions are defined in `trajectory.h`. 

```bash
./general_sweep <grid> <output> [OPTIONS]
```


## 4D Implicit Function Framework: 

The input of this program is any generalized sweep that is represented by a smooth 4D implicit function. Currently, we provide a series of pre-defined functions which include all paper examples. You can specify an implicit function file or use one of the predefined function names. Unfortunately, we do not provide a GUI or a user-friendly tool for defining such functions at this moment. If you want to specify your own sweep, please refer to this [repository](https://github.com/qnzhou/space-time-functions) and specific use cases in `trajectory.h` for details.

### Positional Arguments

- `grid` : The path to the initial grid file that will be used for grid generation. This file can be either a `.msh` or `.json` file. When a `.json` file is provided, it will be converted to a `.msh` file internally. The format to specify an initial grid can be found [here](https://github.com/Jurwen/Swept-Volume/blob/main/data/test/grid_1.json);
- `output` : The output directory path where all generated files will be saved. The tool will create this directory if it doesn't exist. Output files include:
  - `contour.msh` : The mesh before arrangement (if `SAVE_CONTOUR` is enabled)
  - `0.obj`, `1.obj`, ... : Separated cell components with 0 winding number.
  - `features.json` : Feature lines (first slot) and points (second slot) 

### Options

- `-h, --help` : Show the help message and exit the program.
- `-f, --function <file>` : Optional. Specify an implicit function file or predefined function name. Can be:
  - A predefined function name (e.g., `fertility_v4`, `kitten_dog`, `letter_L_blend`, `ball_genus_roll`, `tangle_chair_S`, `star_S`, etc.)
  - See `trajectory.h` for all available predefined functions
- `-t, --threshold <value>` : Set the threshold value for grid generation (default: 0.0005). Lower values produce coarser grids. This is a DOUBLE value that controls the precision level.
- `--tt, --traj-threshold <value>` : Set the threshold value for trajectory processing (default: 0.005). This is a DOUBLE value that controls trajectory precision.
- `-m, --max-splits <number>` : Set the maximum number of splits for grid generation to avoid infinite subdivision (default: unlimited). This is a sanity parameter to prevent degeneracies.
- `--without-snapping` : Disable vertex snapping in the iso-surfacing step.
- `--without-optimal-triangulation` : Disable optimal triangulation in the iso-surfacing triangulation step. 
