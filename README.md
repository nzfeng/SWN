# Winding Numbers on Discrete Surfaces

**WARNING:** This repo is still under construction! The full code release will be soon, probably before SIGGRAPH 2023.

C++ demo for "[Winding Numbers on Discrete Surfaces](https://nzfeng.github.io/research/WNoDS/index.html)" by [Nicole Feng](https://nzfeng.github.io/index.html), [Mark Gillespie](https://markjgillespie.com/), and [Keenan Crane](https://www.cs.cmu.edu/~kmcrane/), presented at SIGGRAPH 2023.

Paper PDF (4.4mb): [link](https://nzfeng.github.io/research/WNoDS/WNoDS.pdf)

Talk (10 minutes): [link]()

![teaser image](media/teaser.png)

If this code contributes to academic work, please cite as:
```bibtex
@article{Feng:2023:WND,
    author = {Feng, Nicole and Gillespie, Mark and Crane, Keenan},
    title = {Winding Numbers on Discrete Surfaces},
    year = {2023},
    issue_date = {August 2023},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    volume = {42},
    number = {4},
    issn = {0730-0301},
    url = {https://doi.org/10.1145/3592401},
    doi = {10.1145/3592401},
    journal = {ACM Trans. Graph.},
    month = {jul},
    articleno = {36}
}
```

# Getting started

## Gurobi
The program relies on Gurobi to solve linear programs, so you must first install [Gurobi](https://www.gurobi.com/). Those affiliated with an university can get the academic version for free. Otherwise, a free trial is available.

TODO: editing custom path in `cmake/modules/FindGUROBI.cmake`
<!-- [How do I resolve "undefined reference" errors while linking Gurobi in C++?](https://support.gurobi.com/hc/en-us/articles/360039093112) -->

## Running the program
```
git clone --recursive https://github.com/nzfeng/SWN.git
cd SWN
mkdir build
cd build
cmake -DCOMISO_BUILD_WITHOUT_BOOST=On -DCMAKE_BUILD_TYPE=Release ..
make -j8 # or however many cores you have or want to use
bin/main /path/to/mesh --c=/path/to/curve --viz
```

A Polyscope GUI will open:

![Screenshot of Polyscope GUI](media/GUI.png)

# Usage

Once the executable `bin/main` has been built, running `bin/main /path/to/mesh /path/to/curve` will output the input mesh with texture coordinates containing the winding number solution. 

You can pass several solve options to the command line, options which are also shown in the GUI:

|flag | purpose|
| ------------- |-------------|
|`--curveFilename=input.[obj,txt]`| Filepath to input curve. |
|`--allowRemeshing`=true| Allow the input surface mesh to be re-meshed |
|`--outputFilename=output.obj`| File to save output mesh to, along with homogeneous texture coordinates |
|`--correctNonbounding`=true| Correct for nonbounding curve components |
|`--approximateResidual`=false| Use reduced-size linear program to approximate residual function, instead of solving a more expensive LP. |
<!--|`--viz`| Show the GUI |-->
|`--verbose`, `-V`| Verbose output |
|`--version`, `-v`| Version info |
|`--help`, `-h`| Display help |

TODO: CL flags for all output?

## Curve input
TODO

## File formats
The input mesh may be an `obj`, `ply`, `off`, or `stl`. See [the geometry-central website](https://geometry-central.net/surface/utilities/io/#reading-meshes) for up-to-date information on supported file types.

The input curve can be specified either as a discrete _1-chain_ or a _dual 1-chain_. Specifying the curve as a 1-chain encodes the curve with _tangential orientation_, while a dual 1-chain encodes _normal orientation_. Using a dual 1-chain is only necessary if the surface is non-orientable.

<!-- Ordinarily, we assume that jumps increase in the direction obtained by rotating the tangent 90 degrees counter-clockwise. On
a nonorientable surface, however, there is no consistent notion of counter-clockwise—even though curves can still meaningfully bound regions -->

## Re-meshing

Depending on the curve input, the surface may need to be re-meshed for optimal output.

First, SWN is formulated for curves that conform to mesh edges. However, you can still specify curves generically as sequences of barycentric points along the surface, with the condition that the curve is continuous and linear within each triangle face. If `--allowRemeshing=true`, the surface will be re-meshed so that the curve lies entirely along mesh edges (with no change to the curve geometry.) If `--allowRemeshing=false`, the program will use a Poisson formulation of SWN to give valid output, but homology correction is no longer possible.

Second, since curve endpoints are omitted from the solve (Section 2.3.2 in the paper), curves that span only one mesh edge will be ignored. If `--allowRemeshing=true`, such edges will be subdivided so that these parts of the curve will not be ignored.

## Intrinsic re-meshing

The GUI can perform intrinsic re-meshing, for example generating an intrinsic Delaunay triangulation for better numerical behavior. If you choose to solve on an intrinsic mesh, exporting the solution will export TODO

Notes: 
* Only manifold meshes can be intrinsically re-triangulated.
* Delaunay refinement may not terminate with a minimum angle value >30 degrees.

# Output

From the GUI menu, you can export the solution, as well as other intermediate functions, as an OBJ file containing both the mesh and texture coordinates. 

From the GUI menu, you can also export curves as [OBJ line elements](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Line_elements). In most 3D graphics software, how smooth a curve is rendered depends on its connectivity (i.e., there can be noticeable gaps between curves specified as separate segments.) Therefore, this program will automatically greedily compute the connected components of the curve before exporting, so that the curve can be written with maximum connectivity and hence yield maximum smoothness during rendering. 

Whenever one calls an export function, the data that gets exported is from most recently computed solution.

## Homogenous texture coordinates

For better visualization of the solution around curve endpoints, the Polyscope shader performs projective interpolation (see Figure 4 from the paper.) 

By default, we output texture coordinates in _homogeneous coordinates_, where rather than the standard uv coordinates, we output 3-dimensional uvw texture coordinates. To visualize these textures you can interpolate the 3d coordinates linearly across each triangle. Then, for each pixel you perform a homogeneous divide, dividing the first two coordinates by the last coordinate to obtain the final uv texture coordinates. This can be done e.g. in a shader or via shader nodes in Blender (see `Example.blend` for an example).

## Visualization

The `render/` directory contains an example Blender file (`Example.blend`) that can load and visualize meshes and curves, with the SWN solution. 

<!-- The blender file should open to a Python script in the `Scripting` workspace. You can load your own uniformized mesh by changing the mesh name in the script and clicking on `Run Script`. This will load your model and apply a correctly-interpolated checkerboard texture. -->

![Screenshot of the provided Blender file](media/BlenderFile.png)
TODO: example Blender file