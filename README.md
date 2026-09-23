# Tensor Mapper Documentation

## Overview

Tensor Mapper is a C++ application for transferring data fields between different mesh formats. It supports various mapping algorithms and can handle VTK unstructured grids (.vtu), Abaqus input files (.inp), and NumPy arrays (.npy).

## Supported File Formats

### Input Formats
- **VTU**: VTK unstructured grid format
- **INP**: Abaqus input files with nodes, elements, and element sets
- **NPY**: NumPy arrays (for material/voxel data) in metres

### Output Formats
- **VTU**: VTK unstructured grid format
- **INP**: Abaqus input files with field data

## Command Line Usage

```bash
tensor_mapper <source_file> <target_file> <output_file> [options]
```

### Command Line Options

#### Mapping Methods
- `--direct`: Direct (nearest point) mapping (default)
- `--npoint`: N-point averaging
- `--proximity`: Proximity-weighted mapping
- `--gaussian`: Gaussian weighting 

#### Mapping Parameters
- `--radius <value>`: Search radius for proximity mapping and Gaussian mapping
- `--neighbors <n>`: Number of neighbors for N-point mapping (default: 5)
- `--sigma <value>`: Gaussian sigma parameter

#### Smoothing Options
- `--smooth`: Enable Laplacian smoothing
- `--smooth-iterations <n>`: Number of smoothing iterations (default: 1)
- `--relax <factor>`: Relaxation factor for smoothing (default: 0.5)

#### Data Selection
- `--data-name <name>`: Specify which data array to map
- `--target-location <nodal|elemental>`: Target data location
  - `nodal`: Maps to mesh points/vertices
  - `elemental`: Maps to mesh cells/elements
  - NPY material mapping always outputs elemental data

#### NPY File Options
- `--origin <x,y,z>`: Origin coordinates for NPY files (default: 0,0,0)
- `--spacing <dx,dy,dz>`: Spacing between voxels (default: 1,1,1)
- `--material <id> <name>`: Define material ID and name mappings

## Mapping Methods

### Direct Mapping
Maps each target point/element to the nearest source point/element.

**Use case**: Simple transfer when meshes are similar
**Parameters**: None

### N-Point Mapping
Averages values from N nearest neighbors with equal weights.

**Parameters**: `--neighbors <n>`

### Proximity Mapping
Uses inverse distance weighting for all points within a specified radius.

**Use case**: Interpolation with distance based weighting
**Parameters**: `--radius <value>`

### Gaussian Mapping
Applies Gaussian-weighted averaging where weights are determined by the function: weight = exp(-distance²/(2*sigma²))

Use case: Natural smooth transitions with gradual influence decay

Parameters:

--radius <value>: Maximum search distance (auto-calculates as 3*sigma if not specified)

--sigma <value>: Controls the spread of the Gaussian bell curve

## Smoothing

Laplacian smoothing reduces noise and creates smoother field distributions:

```
newValue = originalValue + relaxationFactor * (neighborAverage - originalValue)
```

- **High relaxation factor**: Favors neighbor values (more smoothing)
- **Low relaxation factor**: Preserves original data better

## File Format Details

### VTU Files
Standard VTK unstructured grid format supporting:
- Point and cell data arrays
- Multiple data components
- Various element types

### INP Files
Abaqus input format with support for:
- Node definitions (`*NODE`)
- Element definitions (`*ELEMENT, TYPE=...`)
- Element sets (`*ELSET`)
- Field data output as distribution tables

### NPY Files
NumPy arrays representing 3D voxel data:
- Must be 3D arrays (depth × height × width)
- Supports uint8 material IDs
- Requires origin and spacing parameters for coordinate mapping

## Material Mapping (NPY to Mesh)

When using NPY files as source, the application:

1. **Loads voxel data**: Reads 3D material ID array
2. **Maps coordinates**: Converts voxel indices to world coordinates
3. **Cell intersection**: Determines which voxels fall within each mesh element
4. **Counts materials**: Tallies material IDs per element
5. **Creates output arrays**:
   - `MaterialName_Count`: Count of each material per element
   - `DominantMaterial`: ID of most prevalent material per element

### Coordinate System
- NPY arrays use (x,y,z) indexing
- Mesh coordinates are mapped with x↔z flipping for proper orientation
- Origin and spacing convert voxel indices to world coordinates

## Dependencies

- **VTK**: Visualization Toolkit for mesh handling
- **OpenMP**: Parallel processing support
- **ZLIB**: Compression support for NPY files
- **cnpy**: NumPy file reading library

## Run instructions
Just run the exe from the folder, make sure dlls are in same folder.

## Examples

```bash
# NPY material mapping with smoothing
tensor_mapper tubec_points.npy midtow_out_tubec.vtu output_smoothed.vtu --smooth --smooth-iterations 1 --origin -0.029369011121814493,-0.033022958968174136,0.7303319676666113 --spacing 0.00003734269999999991,0.000037342699999999916,0.000037342699999994034 --material 1 voids


# Direct mapping
tensor_mapper source.inp target.vtu output.vtu --direct --data-name CellVoidCount --target-location nodal

# N-point averaging
tensor_mapper data.vtu mesh.inp result.inp --npoint --neighbors 5 --target-location elemental

# Proximity mapping
tensor_mapper mesh1.vtu mesh2.inp result.vtu --proximity --radius 2.0 --target-location nodal --data-name CellVoidCount

# Gaussian mapping
tensor_mapper source.vtu target.inp result.vtu --gaussian --radius 5.0 --sigma 1.5 --target-location nodal

# Direct mapping with smoothing
tensor_mapper mesh1.vtu mesh2.vtu result.vtu --direct --smooth --smooth-iterations 3 --relax 0.4 --target-location elemental --data-name CellVoidCount
```

## Build Instructions (only needed if you are making from the source code)

```bash
# Configure with CMake

mkdir build
cd build

cmake ..

# Build
cmake --build . --config Release
```
Then copy the contents (exe and dlls) of the release folder to where you want to run it from. 

