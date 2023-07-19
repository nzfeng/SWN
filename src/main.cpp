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

// intrinsic triangulation stuff
std::unique_ptr<ManifoldSurfaceMesh> manifoldMesh;
std::unique_ptr<VertexPositionGeometry> manifoldGeom;
std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri;

// Polyscope stuff
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceMesh* psCsMesh; // common subdivision
polyscope::SurfaceGraphQuantity* intEdgeQ;

// The SWN solver
std::unique_ptr<SurfaceWindingNumbersSolver> SWNSolver;

// Parameters
std::string MESHNAME = "input mesh";
std::string MESHNAME = "input mesh";
enum SolverMode { OriginalMesh = 0, IntrinsicMesh };
int SOLVER_MODE = SolverMode::OriginalMesh;
bool DO_HOMOLOGY_CORRECTION = true;
float REFINE_AREA_THRESH = std::numeric_limits<float>::infinity();
float REFINE_ANGLE_THRESH = 30.;
int MAX_INSERTIONS = -1;
bool USING_MANIFOLD_MESH = false; // if using the re-meshed (manifold) mesh
bool USE_SPECIAL_BASES = true;
int MAX_EDGE_SPLITS = 3;

std::vector<Halfedge> curveHalfedgesOnManifold;
std::vector<Halfedge> curveHalfedgesOnIntrinsic;

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

    switchPolyscopeMeshToDisplay(false, true);
    intEdgeQ = psMesh->addSurfaceGraphQuantity("intrinsic edges", result);
    intEdgeQ->setEnabled(true);
    intEdgeQ->setColor(polyscope::render::RGB_ORANGE);
    intEdgeQ->setRadius(0.0005);
}

void functionCallback {

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
    std::string loadedFilename = args::get(inputFilename);
    MESH_FILEPATH = loadedFilename;
    std::tie(mesh, geometry) = readSurfaceMesh(loadedFilename);

    // Initialize polyscope
    polyscope::init();

    polyscope::state::userCallback = functionCallback;

    // Initialize solver.
    SWNSolver = std::unique_ptr<SurfaceWindingNumbersSolver>(new SurfaceWindingNumbersSolver(*geometry));

    // Register the mesh with polyscope
    MESHNAME = polyscope::guessNiceNameFromPath(loadedFilename);
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // TODO: Read curve.
    if (curveFilename) {
        std::string curveFilepath = args::get(curveFilename);
        CURVE_FILEPATH = curveFilepath;
        std::cerr << "Reading input curves..." << std::endl;
        readInputCurves(*mesh, curveFilepath, curveNodes, inputCurveEdges, curveEdges, curveHalfedges, MARKED_EDGES);
        std::cerr << "Setting input halfedges..." << std::endl;
        if (curveHalfedges.size() == 0) curveHalfedges = setCurveHalfedges(curveNodes, curveEdges, true);
        std::cerr << curveHalfedges.size() << std::endl;
        if (curveNodes.size() + inputCurveEdges.size() + curveEdges.size() + curveHalfedges.size() != 0) {
            if (mesh->isManifold()) {
                HERE();
                // std::cerr << "Remeshing if necessary..." << std::endl;
                // bool wasRemeshed = remeshIfNecessary();
                // std::cerr << "Splitting edges if necessary..." << std::endl;
                // if (wasRemeshed) {
                //     std::cerr << "mesh split along curve" << std::endl;
                //     splitEdgesIfNecessary(curveHalfedgesOnManifold);
                //     std::cerr << manifoldMesh->nVertices() << " " << manifoldMesh->nEdges() << " "
                //               << manifoldMesh->nFaces() << std::endl;
                //     std::cerr << "edges split" << std::endl;
                // } else {
                ensureHaveManifoldMesh();

                // std::vector<Halfedge> halfedgesOnManifold;
                // geometry->requireVertexIndices();
                // for (Halfedge he : curveHalfedges) {
                //     Vertex vA = manifoldMesh->vertex(geometry->vertexIndices[he.tailVertex()]);
                //     Vertex vB = manifoldMesh->vertex(geometry->vertexIndices[he.tipVertex()]);
                //     std::cerr << he.tailVertex() << " " << he.tipVertex() << " " << vA << " " << vB << std::endl;
                //     halfedgesOnManifold.push_back(determineHalfedgeFromVertices(vA, vB));
                // }
                // geometry->requireVertexIndices();

                // for (Halfedge he : curveHalfedges) {
                //     halfedgesOnManifold.push_back(manifoldMesh->halfedge(he.getIndex()));
                // }

                std::vector<SurfacePoint> curveNodesOnManifold;
                for (const auto& pt : curveNodes) {
                    curveNodesOnManifold.push_back(reinterpretTo(pt, *manifoldMesh));
                }
                curveHalfedgesOnManifold = setCurveHalfedges(curveNodesOnManifold, curveEdges, true);
                std::cout << "curveHalfedgesOnManifold.size(): " << curveHalfedgesOnManifold.size() << std::endl;

                // splitEdgesIfNecessary(curveHalfedgesOnManifold);
                // std::cerr << "edges split" << std::endl;
                // }
            }
        }
    }

    if (mesh->isManifold()) {
        ensureHaveIntrinsicSolver();
    }

    std::cerr << "Displaying curves..." << std::endl;
    displayCurves(getGeom(), getCurveHalfedges(), curveNodes, curveEdges, psMesh, MARKED_EDGES);

    polyscope::show();

    return EXIT_SUCCESS;
}