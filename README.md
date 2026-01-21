# Flat-land
Computes the projected area of an input obj mesh from any view angle using a rasterization based technique with dynamic windowing. Tool has been validated and benchmarked on a few datasets, whitepaper is to come.

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
