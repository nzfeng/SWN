# Winding Numbers on Discrete Surfaces

C++ demo for "[Winding Numbers on Discrete Surfaces](https://nzfeng.github.io/research/WNoDS/index.html)" by [Nicole Feng](https://nzfeng.github.io/index.html), [Mark Gillespie](https://markjgillespie.com/), and [Keenan Crane](https://www.cs.cmu.edu/~kmcrane/), presented at SIGGRAPH 2023.

Paper PDF (4.4mb): [link](https://nzfeng.github.io/research/WNoDS/WNoDS.pdf)
Talk (10 minutes): [link]()

![teaser image](media/teaser.png)

If this code contributes to academic work, please cite as:
```bibtex
@article{Feng:2023:WND,
    author = {Feng, Nicole and Gillespie, Mark and Crane, Keenan},
    title = {Winding Numbers on Discrete Surfaces,
    journal = {ACM Trans. Graph.},
    volume = {42},
    number = {4},
    month = aug,
    year = {2023},
    articleno = {},
    doi = {},
    publisher = {ACM},
    address = {New York, NY, USA}
}
```
TODO: Fill in article numbers.

# Getting started

```
git clone --recursive https://github.com/nzfeng/SWN.git
cd SWN
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8 # or however many cores you have or want to use
bin/main /path/to/mesh /path/to/curve --viz
```

A Polyscope GUI will open:

![Screenshot of Polyscope GUI](media/GUI.png)

# Usage

Once the executable `bin/main` has been built, running `bin/main /path/to/mesh /path/to/curve` will output the input mesh with texture coordinates containing the winding number solution. 

You can pass several solve options to the command line, options which are also shown in the GUI:

|flag | purpose|
| ------------- |-------------|
|`--outputFilename=output.obj`| File to save output mesh to, along with homogeneous texture coordinates |
|`--doHomologyCorrection`=true| Correct for nonbounding curve components |
|`--allowRemeshing`=true| Allow the input surface mesh to be re-meshed |
|`--viz`| Show the GUI |
|`--version`, `-v`| Version info |
|`--help`, `-h`| Display help |

## File formats
The input mesh may be an `obj`, `ply`, `off`, or `stl`. See [the geometry-central website](https://geometry-central.net/surface/utilities/io/#reading-meshes) for up-to-date information on supported file types.

The input curve can be specified either as a discrete _1-chain_ or a _dual 1-chain_. Specifying the curve as a 1-chain encodes the curve with _tangential orientation_, while a dual 1-chain encodes _normal orientation_. Using a dual 1-chain is only necessary if the surface is non-orientable.

<!-- Ordinarily, we assume that jumps increase in the direction obtained by rotating the tangent 90 degrees counter-clockwise. On
a nonorientable surface, however, there is no consistent notion of counter-clockwise—even though curves can still meaningfully bound regions -->

## Homogenous texture coordinates

For better visualization of the solution around curve endpoints, we perform projective interpolation (see Figure 4 from the paper.) 

By default, we output texture coordinates in _homogeneous coordinates_, where rather than the standard uv coordinates, we output 3-dimensional uvw texture coordinates. To visualize these textures you can interpolate the 3d coordinates linearly across each triangle. Then, for each pixel you perform a homogeneous divide, dividing the first two coordinates by the last coordinate to obtain the final uv texture coordinates. This can be done e.g. in a shader or via shader nodes in Blender (see `Example.blend` for an example).

## Remeshing

Depending on the curve input, the surface may need to be re-meshed for optimal output.

SWN is formulated for curves that conform to mesh edges. However, you can still specify curves generically as sequences of barycentric points along the surface, with the condition that the curve is continuous and linear within each triangle face. If `--allowRemeshing=true`, the surface will be re-meshed so that the curve lies entirely along mesh edges (with no change to the curve geometry.) If `--allowRemeshing=false`, the program will use a Poisson formulation of SWN to give valid output, but homology correction is no longer possible.

Second, since curve endpoints are omitted from the solve (Section 2.3.2 in the paper), curves that span only one mesh edge will be ignored. If `--allowRemeshing=true`, such edges will be subdivided so that these parts of the curve will not be ignored.

# Visualization

The `render/` directory contains an example Blender file (`Example.blend`) that can load and visualize meshes and curves, with the SWN solution. 

<!-- The blender file should open to a Python script in the `Scripting` workspace. You can load your own uniformized mesh by changing the mesh name in the script and clicking on `Run Script`. This will load your model and apply a correctly-interpolated checkerboard texture. -->

![Screenshot of the provided Blender file](media/BlenderFile.png)
TODO: example Blender file