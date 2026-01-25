# FlatLand
FlatLand is a modular C++ utility designed for the geometric analysis and visualization of 3D models via orthographic projection and rasterization.

While primarily focused on computing the projected silhouette and visible surface area of 3D meshes (OBJ format), the latest version introduces scalar field mapping. This allows users to project arbitrary data (e.g., temperature, stress, pressure) defined on nodes or faces onto the 2D viewing plane, generating average value statistics and heatmap visualizations.

The tool features a custom Z-buffer implementation for depth testing and edge-function rasterization for high-fidelity coverage calculations. It is optimized for batch processing, supporting runtime precision switching (Float vs Double) and per-view resolution settings.

## Key Features
 - **Projected Area Analysis**: Accurate calculation of visible surface area from arbitrary viewing angles.
 - **Scalar Data Mapping**: Load external data files (node-based or face-based) to compute visible average values and generate heatmap images.
 - **Flexible Batch Mode**: Process thousands of views with a single load of the mesh geometry. Supports overriding resolution and data files per-view.
 - **Runtime Precision**: Switch between float (speed) and double (accuracy) arithmetic via CLI flags.
 - **Visual Output**: export .ppm heatmaps of the projected mesh.
 - **JSON Integration**: Structured output for easy parsing in automated pipelines.

## Build Instructions
The tool relies on standard C++17 and CMake (3.14+). No third-party dependencies are required.
```shell
# Clone the repository
git clone https://github.com/Interfluo/flatland.git
cd flatland

# Build using the provided helper script
./run_build.sh
```
Alternatively, you can build manually via CMake:
```shell
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
```

To verify the build:
```shell
./build/flatland --help
```

## Usage
```shell
./flatland <mesh.obj> [options]
```
### Core options
| Flag | Arguments | Description |
| :--- | :-------- | :---------- |
|  -v  | x y z     | Add a single view direction. Can be repeated multiple times. | 
|  -b  | file      | Load views from a batch file (see format below). | 
|  -o  | prefix    | "Save heatmap images with this filename prefix (e.g., out_0000.ppm)." | 
|  -j  | none      | Output results as structured JSON to stdout. | 
|  -h  | none      | Show help message. |

### Simulation Parameters
| Flag | Arguments | Description |
| :--- | :-------- | :---------- |
| -r        | val | Set default pixel resolution (default: 0.001). Used if a batch entry doesn't specify one.| 
| -d        | file |  "Path to a default scalar data file (text file, one value per line)."| 
| -p        | mode |  "Set math precision: float (default, faster) or double (higher accuracy)."| 
| --no-cull |none |  Disable backface culling (renders all matching faces).| 

### Batch File Format
The batch file (-b) allows for defining mixed-resolution or mixed-data jobs. Each line represents a view.
```text
<nx> <ny> <nz> [resolution] [data_file_path]
```
 - resolution and data_file_path are optional.
 - If omitted, the global defaults (set via -r and -d) are used.

**Example** views.txt:
```text
# Standard view using global defaults
1.0 0.0 0.0

# High-resolution view (0.005 pixel size)
0.0 1.0 0.0 0.005

# View with specific scalar data file overlay
0.0 0.0 1.0 0.002 ./pressure_load_case_2.txt
```

### Examples
**1. Basic Area Calculation**: Calculate the projected area of bunny.obj from the X-axis using default settings.
```shell
./flatland bunny.obj -v 1 0 0
```
**2. Heatmap Generation with Double Precision**: Map values from temps.txt onto the mesh, render a heatmap to heat_0000.ppm, using high-precision math.
```shell
./flatland engine.obj -v 0 1 1 -d temps.txt -p double -o heat
```
**3. JSON Pipeline**: Run a batch of views defined in views.txt, outputting strict JSON for a Python script to ingest.
```shell
./flatland part.obj -b views.txt -j > results.json
```

## Benchmarks

### Benchmark 1: Bunny 
In this benchmark, we test the impact of varying resolution (pixel area) against computed projected area and note execution times. We use the [Stanford Bunny](https://docs.pyvista.org/api/examples/_autosummary/pyvista.examples.downloads.download_bunny.html#pyvista.examples.downloads.download_bunny)  geometry.

<div align="center" style="background-color: white; padding: 10px;">
  <img src="https://github.com/user-attachments/assets/eaefb062-e3fc-403a-8bed-a19a79b1a1dd" width="80%" />
  <p>Pixel size vs image quality</p>
</div>

The projected area calculation converges as expected, and execution times remain fast even for the highest resolution cases.

<table align="center" style="background-color: white; border-collapse: collapse; border: none;">
  <tr>
    <!-- Left Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/a425288f-d986-4660-914c-7a2a5031622b" width="100%" />
      <p style="color: black; margin-top: 10px;">Projected Area error for Bunnyn</p>
    </td>
    <!-- Right Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/2c9f6436-a825-4480-8ccc-06277a2f9989" width="100%" />
      <p style="color: black; margin-top: 10px;">Execution time for Bunny</p>
    </td>
  </tr>
</table>

*Note, that the lowest resolution image is 4x4 pixels while highest resolution image is 12070x15402 pixels.*

### Benchmark 2: Sphere
This benchmark exercises the tool on an [Icosphere](https://docs.pyvista.org/api/utilities/_autosummary/pyvista.icosphere
) of varying subdivision levels. This is a vital validation case since the area of a unit sphere is known to be $\pi$ from any view angle.

<div align="center" style="background-color: white; padding: 10px;">
  <img src="https://github.com/user-attachments/assets/de64a054-27bb-41b9-83e1-f89fdeecfbec" width="80%" />
  <p>Varying Icosphere subdivision levels</p>
</div>

The results demonstrate exponential convergence towards the correct analytical answer as mesh and pixel resolution increase.

<table align="center" style="background-color: white; border-collapse: collapse; border: none;">
  <tr>
    <!-- Left Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/8c09d227-fee3-42be-9b14-3b63958bc537" width="100%" />
      <p style="color: black; margin-top: 10px;">Projected Area error for Icosphere</p>
    </td>
    <!-- Right Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/17a51518-9d43-4bd4-857f-b6e737c47705" width="100%" />
      <p style="color: black; margin-top: 10px;">Execution time for Icosphere</p>
    </td>
  </tr>
</table>


## Notes / Tips

1. **Batch Processing**: For large datasets, always use the batch file (-b) option. Repeatedly calling the tool from the CLI incurs overhead from re-loading and parsing the OBJ file every time.

2. **Resolution**: Use the image output (-o) to perform a sweep of pixel sizes when setting up a new case. Ensure you are not under-resolving thin geometry features.

3. **Culling**: Running without back-face culling (--no-cull) is supported but rarely necessary unless you are dealing with non-manifold or open surface meshes where "inside" faces contribute to the silhouette.

