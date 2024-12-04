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

To use the `gridgen` tool, you must provide an initial grid file and implicit function file as arguments, along with any desired options.

```bash
./general_sweep <grid> <function> [OPTIONS]
```

### Positional Arguments

- `grid` : The path to the initial grid file that will be used for gridgen. This file can either be a `.msh` or `.json` file. 
Examples of grid files can be found in the `data/grid` directory.
- `function` : The path to the implicit function file that is used as the obejct to sweep.

### Options

- `-h, --help` : Show the help message and exit the program.
- `-t, --threshold` : Set the threshold value for the grid generation of the swept volume. The lower the value, the coarser the grid. This is a `DOUBLE` value that defines the precision level of the gridgen.
