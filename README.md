# Swept-Volume

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

Currently, all the packages dependencies are available.

## Usage

To use the `general_sweep` tool, you must provide an initial grid file as arguments, along with any desired options. The functions are now hard-coded in `trajectory.h`. Feel free to change it to any lambda function at `implicit_sweep`.  build

```bash
./general_sweep <grid> <function> [OPTIONS]
```

### Positional Arguments

- `grid` : The path to the initial grid file that will be used for gridgen. This file can either be a `.msh` or `.json` file. 
Examples of grid files can be found in the `data/grid` directory.

### Options

- `-h, --help` : Show the help message and exit the program.
- `function` : Optional. A mesh file that uses SDF as implicit functions. This serves the purpose to compare with Silvia's code.
- `-t, --threshold` : Set the threshold value for the grid generation of the swept volume. The lower the value, the coarser the grid. This is a `DOUBLE` value that defines the precision level of the gridgen.
- `-m, --max-splits` : Set the maximum number of splits for the grid generation in order to avoid degeneracies, which lead to inifinite subdivision. This is a sanity parameter. 
