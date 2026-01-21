# Flat-land
FlatLand is a modular C++ utility designed for the geometric analysis and visualization of 3D models via orthographic projection and rasterization. The tool primarily focuses on computing the projected silhouette and visible surface area of 3D meshes (loaded via the OBJ file format) from arbitrary viewing angles. It utilizes a custom Z-buffer implementation to handle occlusions and depth testing, ensuring that only the forward-most surfaces contribute to the final calculations.

Technically, the tool functions by transforming 3D triangles into a 2D coordinate system defined by a user-specified "view normal." Once projected, the software performs high-fidelity rasterization using barycentric coordinates to determine pixel coverage. It calculates the total "covered area" by mapping triangles onto a discrete image grid, where each pixel’s status (occupied or empty) is tracked. This makes the tool particularly useful for engineering or architectural applications where understanding the footprint or visible profile of a complex object is necessary.

The software is built for efficiency and batch processing, supporting features like automated output naming and modular math structures for 3D vector operations. Beyond numerical data—such as visible triangle counts and total projected area—FlatLand can generate visual output in the form of PPM (Portable Pixmap) images. These images represent the silhouette of the object, providing a clear, high-contrast visual representation of the mesh from the perspective of the calculated view frame.

## Build Instructions
Tool doesn't rely on any third party libraries so it is extremely straight forward to build:
```shell
git clone https://github.com/Interfluo/flatland.git
cd src/
./run_build.sh
```

This should create a build folder, to check that everything is working:
```shell
cd build
./flatland --help
```
which should output:
```shell
Usage: ./flatland <input.obj> [options]

Computes visible projected area of an OBJ mesh from one or more view directions.

Options:
  -h, --help                Show this help
  -v, --verbose             Verbose output
  --view <x> <y> <z>        Single view normal (default: 1 0 0)
  --views-file <file>       Text file with multiple view normals (x y z per line)
  -p, --pixel-size <val>    Pixel size (default: 0.001)
  -o, --output <prefix>     Output PPM prefix (batch: prefix_0000.ppm etc.)
  --no-image                Disable PPM output
  --json                    JSON output
  --no-bfc                  Disable back-face culling
```

## A few notes / tips
1. For running large numbers of cases I recommend using the views-file option (see examples folder for an example of one), repeatedly calling the tool from the CLI adds a lot of overhead from reading in the mesh file.
2. I recommend using the image output to do a sweep of pixel-size when initially setting up a case, make sure that you aren't under resolving the geometry. (Worst case scenario running quite high resolution should be pretty fast anyways)
3. I generally recommend running without back-face culling, this is a bit faster and shouldn't make a difference form most situations.

## Benchmark 1: Bunny 
In this first benchmark we test the impact of varying resolution (pixel area) against computed projected area and note execution times. 
We use the stanford bunny geometry for fun.

<div align="center" style="background-color: white; padding: 10px;">
  <img src="https://github.com/user-attachments/assets/eaefb062-e3fc-403a-8bed-a19a79b1a1dd" width="80%" />
  <p>Pixel size vs image quality</p>
</div>

Primary take-away is that the projected area calculation converges as expected and that execution times are quite fast even for the highest resolution case. 
<table align="center" style="background-color: white; border-collapse: collapse; border: none;">
  <tr>
    <!-- Left Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/a425288f-d986-4660-914c-7a2a5031622b" width="100%" />
      <p style="color: black; margin-top: 10px;">Projected Area error vs pixel area for Bunnyn</p>
    </td>
    <!-- Right Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/2c9f6436-a825-4480-8ccc-06277a2f9989" width="100%" />
      <p style="color: black; margin-top: 10px;">Execution time vs pixel area for Bunny</p>
    </td>
  </tr>
</table>

Note, the lowest resolution image is 4x4 pixels while highest resolution image is 12070x15402 pixels.

## Benchmark 2: Sphere
For this benchamrk we excercise the tool on Icosphere of varying resolutions as shown below. This is a great test since we know that the area of a unit sphere is $\pi$ from any view angle. 
<div align="center" style="background-color: white; padding: 10px;">
  <img src="https://github.com/user-attachments/assets/de64a054-27bb-41b9-83e1-f89fdeecfbec" width="80%" />
  <p>Varying IcoSphere subdivision levels</p>
</div>

In this benchamrk we see that the cases exponentially converge towards the correct answer a mesh and pixel resolution increase.
<table align="center" style="background-color: white; border-collapse: collapse; border: none;">
  <tr>
    <!-- Left Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/8c09d227-fee3-42be-9b14-3b63958bc537" width="100%" />
      <p style="color: black; margin-top: 10px;">Projected Area error vs pixel area for Bunnyn</p>
    </td>
    <!-- Right Image -->
    <td align="center" style="background-color: white; border: none; padding: 10px; width: 50%;">
      <img src="https://github.com/user-attachments/assets/17a51518-9d43-4bd4-857f-b6e737c47705" width="100%" />
      <p style="color: black; margin-top: 10px;">Execution time vs pixel area for Bunny</p>
    </td>
  </tr>
</table>

I have additionally done some sanity checks with various other geometries and all cases seem to work well. 
