# Flat-land
Computes the projected area of an input obj mesh from any view angle using a rasterization based technique with dynamic windowing. Tool has been validated and benchmarked on a few datasets shown below.

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

<div align="center" style="background-color: white; padding: 10px;">
  <img src="https://github.com/user-attachments/assets/de64a054-27bb-41b9-83e1-f89fdeecfbec" width="80%" />
  <p>Varying IcoSphere subdivision levels</p>
</div>

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
