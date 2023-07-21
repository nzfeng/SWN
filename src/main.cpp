#include "geometrycentral/surface/common_subdivision.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/mutation_manager.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <chrono>
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

#include "surface_winding_numbers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// == intrinsic triangulation stuff
std::unique_ptr<ManifoldSurfaceMesh> manifoldMesh;
std::unique_ptr<VertexPositionGeometry> manifoldGeom;
std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri;

// == Polyscope stuff
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceMesh* psCsMesh; // common subdivision
polyscope::SurfaceGraphQuantity* intEdgeQ;

// == the SWN solver
std::unique_ptr<SurfaceWindingNumbersSolver> SWNSolver;
std::unique_ptr<SurfaceWindingNumbersSolver> intrinsicSolver;

// == solve parameters
bool DO_HOMOLOGY_CORRECTION = true;
float REFINE_AREA_THRESH = std::numeric_limits<float>::infinity();
float REFINE_ANGLE_THRESH = 30.;
int MAX_INSERTIONS = -1;
bool USING_MANIFOLD_MESH = false; // if using the re-meshed (manifold) mesh
bool USE_SPECIAL_BASES = true;
int MAX_EDGE_SPLITS = 3;

// == program parameters
std::string MESHNAME = "input mesh";
std::string MESH_FILEPATH, CURVE_FILEPATH;
enum SolverMode { OriginalMesh = 0, IntrinsicMesh };
int SOLVER_MODE = SolverMode::OriginalMesh;
bool VIS_INTRINSIC_MESH = false;

// == curve data
std::vector<SurfacePoint> CURVE_NODES;
std::vector<std::array<size_t, 2>> CURVE_EDGES;
std::vector<Halfedge> curveHalfedges;
std::vector<Halfedge> curveHalfedgesOnManifold;
std::vector<Halfedge> curveHalfedgesOnIntrinsic;


SurfaceMesh& getMesh() {
    return USING_MANIFOLD_MESH ? *manifoldMesh : *mesh;
}

VertexPositionGeometry& getGeom() {
    return USING_MANIFOLD_MESH ? *manifoldGeom : *geometry;
}

void setManifoldMesh() {

    USING_MANIFOLD_MESH = true;
    SWNSolver.reset(new SurfaceWindingNumbersSolver(*manifoldGeom));
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, manifoldGeom->vertexPositions, manifoldMesh->getFaceVertexList(),
                                            polyscopePermutations(*manifoldMesh));
}

void ensureHaveManifoldMesh() {

    if (!mesh->isManifold() || !mesh->isOriented())
        throw std::logic_error(
            "SWN: Mesh must be manifold and orientable to convert to geometrycentral::ManifoldSurfaceMesh.");

    if (manifoldMesh == nullptr) {
        manifoldMesh = mesh->toManifoldMesh();
        manifoldMesh->compress();
    }
    if (manifoldGeom == nullptr) {
        manifoldGeom = geometry->reinterpretTo(*manifoldMesh);
        manifoldGeom->refreshQuantities();
    }

    psMesh = polyscope::registerSurfaceMesh(MESHNAME, manifoldGeom->inputVertexPositions,
                                            manifoldMesh->getFaceVertexList(), polyscopePermutations(*manifoldMesh));
}

/*
 * If an intrinsic triangulation has not yet been created, create one. Otherwise do nothing.
 * The intrinsic triangulation is always based on the original input mesh.
 */
void ensureHaveIntrinsicTriangulation() {

    if (intTri != nullptr) return;

    ensureHaveManifoldMesh();

    intTri = std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>(
        new IntegerCoordinatesIntrinsicTriangulation(*manifoldMesh, *manifoldGeom));

    // If curve is given as a set of mesh edges, set these edges as "marked" in the intrinsic triangulation so they
    // never get flipped.
    if (curveHalfedges.size() + curveHalfedgesOnManifold.size() > 0) {

        EdgeData<bool> markedEdges(*intTri->intrinsicMesh, false);

        // The intrinsic triangulation should at first just be a copy of the input mesh.
        if (curveHalfedgesOnManifold.size() == 0) {
            std::vector<SurfacePoint> curveNodesOnManifold;
            for (const auto& pt : CURVE_NODES) {
                curveNodesOnManifold.push_back(reinterpretTo(pt, *manifoldMesh));
            }
            curveHalfedgesOnManifold = setCurveHalfedges(curveNodesOnManifold, CURVE_EDGES);
        }

        setManifoldMesh();
        getGeom().requireEdgeIndices();
        for (Halfedge he : curveHalfedgesOnManifold) {
            size_t eIdx = getGeom().edgeIndices[he.edge()];
            markedEdges[intTri->intrinsicMesh->edge(eIdx)] = true;
        }
        getGeom().unrequireEdgeIndices();

        intTri->setMarkedEdges(markedEdges);
    }
}

void ensureHaveIntrinsicSolver() {

    ensureHaveIntrinsicTriangulation();
    intrinsicSolver.reset(new SurfaceWindingNumbersSolver(*intTri));
}

void visualizeIntrinsicEdges() {

    if (intTri == nullptr) return;

    std::vector<std::vector<Vector3>> result;
    EdgeData<std::vector<SurfacePoint>> traces = intTri->traceAllIntrinsicEdgesAlongInput();
    for (Edge e : intTri->intrinsicMesh->edges()) {
        result.emplace_back();
        std::vector<Vector3>& thisResult = result.back();
        // Convert to 3D positions.
        std::vector<Vector3> positions;
        for (const auto& pt : traces[e]) {
            positions.push_back(pt.interpolate(manifoldGeom->vertexPositions));
        }
        // Add the points to the list
        thisResult.insert(std::end(thisResult), std::begin(positions), std::end(positions));
    }

    intEdgeQ = psMesh->addSurfaceGraphQuantity("intrinsic edges", result);
    intEdgeQ->setEnabled(true);
    intEdgeQ->setColor(polyscope::render::RGB_ORANGE);
    intEdgeQ->setRadius(0.0005);
}

void functionCallback() {

    // Solve on original mesh;
    ImGui::RadioButton("Original mesh", &SOLVER_MODE, SolverMode::OriginalMesh);
    // Solve on an intrinsic mesh
    ImGui::RadioButton("Intrinsic mesh", &SOLVER_MODE, SolverMode::IntrinsicMesh);

    if (ImGui::TreeNode("Advanced solving options")) {
        ImGui::Checkbox("Do homology correction", &DO_HOMOLOGY_CORRECTION);
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("Intrinsic mesh improvement")) {
        if (ImGui::Checkbox("Show intrinsic edges", &VIS_INTRINSIC_MESH)) {
            ensureHaveIntrinsicTriangulation();
            visualizeIntrinsicEdges();
        }

        if (ImGui::Button("Flip to Delaunay")) {
            ensureHaveIntrinsicTriangulation();
            intTri->flipToDelaunay();
            visualizeIntrinsicEdges();
        }

        ImGui::InputFloat("Angle thresh", &REFINE_ANGLE_THRESH);
        ImGui::InputFloat("Area thresh", &REFINE_AREA_THRESH);
        ImGui::InputInt("Max insert", &MAX_INSERTIONS);
        if (ImGui::Button("Delaunay refine")) {
            ensureHaveIntrinsicTriangulation();
            intTri->delaunayRefine(REFINE_ANGLE_THRESH, REFINE_AREA_THRESH, MAX_INSERTIONS);
            visualizeIntrinsicEdges();
        }

        ImGui::TreePop();
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Solve for surface winding number on a triangle mesh.");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
    args::Positional<std::string> curveFilename(parser, "curve", "A curve file.");
    args::Group group(parser);

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    if (!inputFilename) {
        std::cerr << "Please specify a mesh file as argument" << std::endl;
        return EXIT_FAILURE;
    }

    // Load mesh
    MESH_FILEPATH = args::get(inputFilename);
    std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);
    // Read line objects, if they exist in mesh file.
    readLines(*mesh, MESH_FILEPATH, CURVE_NODES, CURVE_EDGES);

    // Initialize polyscope
    polyscope::init();

    polyscope::state::userCallback = functionCallback;

    // Initialize solver.
    SWNSolver = std::unique_ptr<SurfaceWindingNumbersSolver>(new SurfaceWindingNumbersSolver(*geometry));

    // Register the mesh with polyscope
    MESHNAME = polyscope::guessNiceNameFromPath(MESH_FILEPATH);
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Read curve.
    if (curveFilename) {
        std::string CURVE_FILEPATH = args::get(curveFilename);
        readCurves(*mesh, CURVE_FILEPATH, CURVE_NODES, CURVE_EDGES);
    }

    // Display curve.
    displayCurves(*geometry, CURVE_NODES, CURVE_EDGES);

    // // TODO: convert curves to halfedges, if applicable
    // if (curveHalfedges.size() == 0) curveHalfedges = setCurveHalfedges(curveNodes, curveEdges, true);
    // if (curveNodes.size() + inputCurveEdges.size() + curveEdges.size() + curveHalfedges.size() != 0) {
    //     if (mesh->isManifold()) {
    //         HERE();
    //         // std::cerr << "Remeshing if necessary..." << std::endl;
    //         // bool wasRemeshed = remeshIfNecessary();
    //         // std::cerr << "Splitting edges if necessary..." << std::endl;
    //         // if (wasRemeshed) {
    //         //     std::cerr << "mesh split along curve" << std::endl;
    //         //     splitEdgesIfNecessary(curveHalfedgesOnManifold);
    //         //     std::cerr << manifoldMesh->nVertices() << " " << manifoldMesh->nEdges() << " "
    //         //               << manifoldMesh->nFaces() << std::endl;
    //         //     std::cerr << "edges split" << std::endl;
    //         // } else {
    //         ensureHaveManifoldMesh();

    //         // std::vector<Halfedge> halfedgesOnManifold;
    //         // geometry->requireVertexIndices();
    //         // for (Halfedge he : curveHalfedges) {
    //         //     Vertex vA = manifoldMesh->vertex(geometry->vertexIndices[he.tailVertex()]);
    //         //     Vertex vB = manifoldMesh->vertex(geometry->vertexIndices[he.tipVertex()]);
    //         //     std::cerr << he.tailVertex() << " " << he.tipVertex() << " " << vA << " " << vB << std::endl;
    //         //     halfedgesOnManifold.push_back(determineHalfedgeFromVertices(vA, vB));
    //         // }
    //         // geometry->requireVertexIndices();

    //         // for (Halfedge he : curveHalfedges) {
    //         //     halfedgesOnManifold.push_back(manifoldMesh->halfedge(he.getIndex()));
    //         // }

    //         std::vector<SurfacePoint> curveNodesOnManifold;
    //         for (const auto& pt : curveNodes) {
    //             curveNodesOnManifold.push_back(reinterpretTo(pt, *manifoldMesh));
    //         }
    //         curveHalfedgesOnManifold = setCurveHalfedges(curveNodesOnManifold, curveEdges, true);
    //         std::cout << "curveHalfedgesOnManifold.size(): " << curveHalfedgesOnManifold.size() << std::endl;

    //         // splitEdgesIfNecessary(curveHalfedgesOnManifold);
    //         // std::cerr << "edges split" << std::endl;
    //         // }
    //     }
    // }

    // if (mesh->isManifold()) {
    //     ensureHaveIntrinsicSolver();
    // }

    polyscope::show();

    return EXIT_SUCCESS;
}