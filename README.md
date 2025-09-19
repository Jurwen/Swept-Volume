# Lifted Surfacing of Generalized Sweep Volumes

This code implements the ACM SIGGRAPH ASIA 2025 paper: Lifted Surfacing of Generalized Sweep Volumes

<img width="1600" alt="ball-rolling" src="https://github.com/user-attachments/assets/4751b437-091b-4626-976e-3c0a74132838" />
>Top: A wire-like ball rolls forward while offsetting, changing its genus from 41 to 29. Bottom: The time-colored sweep boundary and sharp creases (2nd row), in transparency (3rd row, creases hidden for clarity), and a cut-away view (bottom row).

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


### 4D Implicit Function Framework: 

The input of this program is any generalized sweep that is represented by a smooth 4D implicit function. Currently, we provide a series of pre-defined functions which include all paper examples. You can specify an implicit function file or use one of the predefined function names. Unfortunately, we do not provide a GUI or a user-friendly tool for defining such functions at this moment. If you want to specify your own sweep, please refer to this [repository](https://github.com/qnzhou/space-time-functions) and specific use cases in `trajectory.h` for details.

### Positional Arguments

- `grid` : The path to the initial grid file that will be used for grid generation. This file can be either a `.msh` or `.json` file. When a `.json` file is provided, it will be converted to a `.msh` file internally. The format to specify an initial grid can be found [here](https://github.com/Jurwen/Swept-Volume/blob/main/data/test/grid_1.json);
- `output` : The output directory path where all generated files will be saved. The tool will create this directory if it doesn't exist. Output files include:
  - `contour.msh` : The mesh before arrangement (if `SAVE_CONTOUR` is enabled)
  - `0.obj`, `1.obj`, ... : Separated cell components with 0 winding number.
  - `features.json` : Feature lines (first slot) and points (second slot) 

### Options

- `-h, --help` : Show the help message and exit the program.
- `-f, --function <file>` : Specify an implicit function file or predefined function name. Can be:
  - A predefined function name (e.g., `fertility_v4`, `kitten_dog`, `letter_L_blend`, `ball_genus_roll`, `tangle_chair_S`, `star_S`, etc.)
  - See `trajectory.h` for all available predefined functions
- `-t, --threshold <value>` : Set the threshold value for grid generation (default: 0.0005). Lower values produce coarser grids. This is a DOUBLE value that controls the precision level.
- `--tt, --traj-threshold <value>` : Set the threshold value for trajectory processing (default: 0.005). This is a DOUBLE value that controls trajectory precision.
- `-m, --max-splits <number>` : Set the maximum number of splits for grid generation to avoid infinite subdivision (default: unlimited). This is a sanity parameter to prevent degeneracies.
- `--without-snapping` : Disable vertex snapping in the iso-surfacing step.
- `--without-optimal-triangulation` : Disable optimal triangulation in the iso-surfacing triangulation step. 

## Example:

The following is an example of how to use the `general_sweep` tool with common options:

```bash
./general_sweep ../data/test/grid_1.json ../output/brush_stroke_example -t 0.0005 --tt 0.005 -f brush_stroke_blending
```

### Parameter Breakdown:

This example command demonstrates how to generate a swept volume using the `brush_stroke_blending` function. Let's examine each parameter:

#### Required Positional Arguments:
- **`../data/test/grid_1.json`** : The initial grid file that defines the starting tetrahedral mesh structure. This JSON file contains the initial spatial discretization that will be refined during the swept volume computation.

- **`../output/brush_stroke_example`** : The output directory where all generated files will be saved. The tool will create this directory if it doesn't exist.

#### Optional Parameters:
- **`-t 0.0005`** : Sets the **grid refinement threshold** to 0.0005. This parameter controls how finely the algorithm subdivides the initial grid based on the implicit function's gradient magnitude. A smaller value (like 0.0005) means:
  - Higher precision in capturing surface details
  - More computational time and memory usage
  - Better preservation of sharp features and fine geometric details
  - The algorithm will subdivide grid cells more aggressively where the function changes rapidly

- **`--tt 0.005`** : Sets the **trajectory threshold** to 0.005. This parameter controls the precision of trajectory processing, which is crucial for:
  - Temporal discretization of the 4D sweep
  - Determining how finely to sample the time dimension of the sweep
  - Capturing temporal variations in the implicit function
  - A smaller value provides better temporal resolution but increases computation time

- **`-f brush_stroke_blending`** : Specifies the **implicit function** to use. The `brush_stroke_blending` function represents a particular sweep pattern that:
  - Creates a brush stroke-like swept volume with blending effects
  - Demonstrates complex topological changes during the sweep
  - Shows how multiple geometric elements can be blended together in the swept volume, as the base brush morphs between one sphere and two spheres.
  - Is one of the predefined functions available in `trajectory.h`

### Expected Output:
When this command runs successfully, you'll find in the `../output/brush_stroke_example/` directory:
- **`0.obj`, `1.obj`, ...** : Individual mesh components representing different parts of the swept envelope with 0 winding number
- **`features.json`** : Feature lines and points that capture the topological structure of the swept volume
- **`contour.msh`** : The intermediate envelope mesh before pruning, containing per-vertex info of `time` and `is_regular` (if enabled)

This example showcases how different threshold values balance between computational efficiency and geometric fidelity, making it suitable for both quick prototyping and high-quality final results.
