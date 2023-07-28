#include "geometrycentral/surface/transfer_functions.h"

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

// == This mesh & geometry are used to store the initial input. They should never get changed throughout the program.
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// == Intrinsic triangulation stuff. All remeshing is performed on the manifold mesh.
std::unique_ptr<ManifoldSurfaceMesh> manifoldMesh;
std::unique_ptr<VertexPositionGeometry> manifoldGeom;
std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri;

// == Polyscope stuff
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceMesh* psCsMesh; // common subdivision
polyscope::CurveNetwork* intEdgeQ;

// == SWN solver
std::unique_ptr<SurfaceWindingNumbersSolver> SWNSolver;
std::unique_ptr<SurfaceWindingNumbersSolver> intrinsicSolver;

// == solve parameters
bool ALLOW_REMESHING = true;
bool DO_HOMOLOGY_CORRECTION = true;
bool APPROXIMATE_RESIDUAL = false;
std::string OUTPUT_FILENAME;
float EPSILON = 1e-2;
float REFINE_AREA_THRESH = std::numeric_limits<float>::infinity();
float REFINE_ANGLE_THRESH = 25.;
int MAX_INSERTIONS = -1;
bool USING_MANIFOLD_MESH = false; // if using the re-meshed (manifold) mesh
bool USE_SPECIAL_BASES = true;
int MAX_EDGE_SPLITS = 3;

// == program parameters
std::string MESHNAME = "input mesh";
std::string DATA_DIR, MESHROOT, MESH_FILEPATH, CURVE_FILEPATH;
enum SolverMode { ExtrinsicMesh = 0, IntrinsicMesh };
int SOLVER_MODE = SolverMode::ExtrinsicMesh;
bool VIS_INTRINSIC_MESH = false;
bool USING_GUI = true;
enum CurveExportMode { CLOSED_BOUNDING = 0, CLOSED_NONBOUNDING, CLOSED_ALL, BOUNDING, NONBOUNDING };
enum CurveMode { DUAL = 0, PRIMAL, PRIMAL_MANIFOLD, BARYCENTRIC };
struct SolverState {
    int solverMode;
    int curveMode;
};
SolverState LAST_SOLVE;

// == curve data
std::vector<SurfacePoint> CURVE_NODES;
std::vector<std::array<size_t, 2>> CURVE_EDGES;
std::vector<std::array<Face, 2>> DUAL_CHAIN;
std::vector<Halfedge> curveHalfedges;
std::vector<Halfedge> curveHalfedgesOnManifold;
std::vector<Halfedge> curveHalfedgesOnIntrinsic;

// == intermediate data
std::vector<Halfedge> BOUNDING_PARTS, NONBOUNDING_PARTS, BOUNDING_LOOPS, NONBOUNDING_LOOPS;
std::vector<SurfacePoint> CONTOUR_NODES;
std::vector<std::array<size_t, 2>> CONTOUR_EDGES;


SurfaceMesh& getMesh() {
    return USING_MANIFOLD_MESH ? *manifoldMesh : *mesh;
}

VertexPositionGeometry& getGeom() {
    return USING_MANIFOLD_MESH ? *manifoldGeom : *geometry;
}

std::vector<Halfedge>& getCurveHalfedges() {
    return USING_MANIFOLD_MESH ? curveHalfedgesOnManifold : curveHalfedges;
}

void setManifoldMesh() {

    USING_MANIFOLD_MESH = true;
    SWNSolver.reset(new SurfaceWindingNumbersSolver(*manifoldGeom));
    if (USING_GUI) {
        psMesh =
            polyscope::registerSurfaceMesh(MESHNAME, manifoldGeom->vertexPositions, manifoldMesh->getFaceVertexList());
        psMesh->setAllPermutations(polyscopePermutations(*manifoldMesh));
    }
    // TODO: reinterpret dual chain
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
}

void reinterpretCurveToManifoldMesh() {

    ensureHaveManifoldMesh();
    if (curveHalfedgesOnManifold.size() == 0) {
        std::vector<SurfacePoint> curveNodesOnManifold;
        for (const auto& pt : CURVE_NODES) {
            curveNodesOnManifold.push_back(reinterpretTo(pt, *manifoldMesh));
        }
        curveHalfedgesOnManifold = convertToHalfedges(curveNodesOnManifold, CURVE_EDGES);
    }
}

/*
 * Re-mesh so that the curve is constrained to mesh edges in the new mesh. Return curve as a series of halfedges
 * in the new mesh.
 *
 * Any intrinsic triangulation built on top of the new mesh should mark edges that lie on the curve (using
 * setMarkedEdges()) so that such edges won't be flipped.
 */
std::vector<Halfedge> remeshToConforming() {

    ensureHaveManifoldMesh();

    MutationManager mm(*manifoldMesh, *manifoldGeom);

    // Need to re-interpret curveNodes in terms of manifold mesh.
    std::vector<SurfacePoint> newCurveNodes;
    for (const auto& pt : CURVE_NODES) {
        newCurveNodes.push_back(reinterpretTo(pt, *manifoldMesh));
    }

    std::vector<Halfedge> heOnCurve = mm.cutAlongPath(newCurveNodes, CURVE_EDGES);
    std::cerr << "Mesh cut along " << heOnCurve.size() << " new edges." << std::endl;

    // Keep track of the halfedges in the new mesh that segments correspond to, before compression.
    EdgeData<int> chain(*manifoldMesh, convertToChain(manifoldGeom->mesh, heOnCurve));
    manifoldMesh->compress();

    // Now need to triangulate the newly polygonal faces. To be safe, just call it on every face.
    for (Face f : manifoldMesh->faces()) {
        manifoldMesh->triangulate(f);
    }
    std::cerr << "Mesh re-triangulated." << std::endl;

    manifoldMesh->compress();

    std::vector<Halfedge> newCurve;
    manifoldGeom->requireEdgeIndices();
    for (Edge e : manifoldMesh->edges()) {
        int chainVal = chain[manifoldGeom->edgeIndices[e]];
        if (chainVal != 0) {
            Halfedge he = (chainVal > 0) ? e.halfedge() : e.halfedge().twin();
            for (int i = 0; i < abs(chainVal); i++) newCurve.push_back(he);
        }
    }
    manifoldGeom->unrequireEdgeIndices();

    return newCurve; // on manifold mesh
}

/*
 * If curve is not constrained to mesh edges, re-mesh so that it is.
 */
bool remeshIfNecessary() {

    if (curveHalfedges.size() > 0) return false;

    curveHalfedgesOnManifold = remeshToConforming();

    manifoldMesh->validateConnectivity();

    setManifoldMesh();

    return true;
}

bool splitEdgesIfNecessary(const std::vector<Halfedge>& halfedges) {

    if (halfedges.size() == 0) throw std::logic_error("splitEdgesIfNecessary() called at the wrong time!");

    ensureHaveManifoldMesh();

    EdgeData<int> chain(*manifoldMesh, 0);
    for (Halfedge he : halfedges) {
        Edge e = he.edge();
        int sign = he.orientation() ? 1 : -1;
        chain[e] += sign;
    }

    // Determine curve components.
    std::vector<std::vector<Halfedge>> components = getCurveComponents(*manifoldGeom, halfedges);

    MutationManager mm(*manifoldMesh, *manifoldGeom);

    bool didWeSplit = false;
    int maxEdges = 1 << MAX_EDGE_SPLITS;
    for (const auto& component : components) {
        if (component.size() >= maxEdges) continue;
        int nFold = std::ceil(maxEdges / component.size());
        // get highest bit
        int nSplits = 0;
        while (nFold >>= 1) ++nSplits;

        didWeSplit = true;
        for (const Halfedge& he : component) {
            Edge e = he.edge();
            int chainVal = chain[e];
            Halfedge newHe = mm.splitEdge(e, 0.5);
            Edge otherNewEdge =
                !e.isBoundary() ? newHe.twin().next().twin().next().twin().edge() : newHe.twin().next().twin().edge();
            chain[newHe.edge()] = chainVal;
            chain[otherNewEdge] = chainVal;
            nSplits -= 1;

            std::vector<Edge> edgesToSplit = {otherNewEdge, newHe.edge()};
            std::vector<Edge> nextEdgesToSplit;
            while (nSplits > 0) {
                while (edgesToSplit.size() > 0) {
                    Edge e = edgesToSplit.back();
                    edgesToSplit.pop_back();
                    Halfedge newHe = mm.splitEdge(e, 0.5);
                    Edge otherNewEdge = !e.isBoundary() ? newHe.twin().next().twin().next().twin().edge()
                                                        : newHe.twin().next().twin().edge();
                    chain[newHe.edge()] = chainVal;
                    chain[otherNewEdge] = chainVal;
                    nextEdgesToSplit.push_back(otherNewEdge);
                    nextEdgesToSplit.push_back(newHe.edge());
                }
                nSplits -= 1;
                edgesToSplit = nextEdgesToSplit;
                nextEdgesToSplit.clear();
            }
        }
    }

    curveHalfedgesOnManifold.clear();
    for (Edge e : manifoldMesh->edges()) {
        int val = chain[e];
        Halfedge he = e.halfedge();
        if (val > 0) {
            for (int i = 0; i < val; i++) curveHalfedgesOnManifold.push_back(he);
        } else {
            for (int i = 0; i < abs(val); i++) curveHalfedgesOnManifold.push_back(he.twin());
        }
    }

    manifoldMesh->compress();
    manifoldGeom->refreshQuantities();
    manifoldMesh->validateConnectivity();

    setManifoldMesh();

    return didWeSplit;
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
        reinterpretCurveToManifoldMesh();

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

void resetCurveOnIntrinsicTriangulation() {
    curveHalfedgesOnIntrinsic = determineHalfedgesInIntrinsicTriangulation(*intTri, curveHalfedgesOnManifold);
}

void ensureHaveIntrinsicSolver() {

    ensureHaveIntrinsicTriangulation();
    intrinsicSolver.reset(new SurfaceWindingNumbersSolver(*intTri));
}

void visualizeIntrinsicEdges() {

    if (intTri == nullptr) return;

    std::vector<Vector3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    EdgeData<std::vector<SurfacePoint>> traces = intTri->traceAllIntrinsicEdgesAlongInput();
    for (Edge e : intTri->intrinsicMesh->edges()) {

        // Convert to 3D positions.
        size_t N = nodes.size();
        size_t M = traces[e].size();
        std::vector<Vector3> positions;
        for (const auto& pt : traces[e]) {
            nodes.push_back(pt.interpolate(manifoldGeom->vertexPositions));
        }
        for (size_t i = 0; i < M - 1; i++) edges.push_back({N + i, N + i + 1});
    }

    intEdgeQ = polyscope::registerCurveNetwork("intrinsic edges", nodes, edges);
    intEdgeQ->setEnabled(VIS_INTRINSIC_MESH);
    intEdgeQ->setColor(polyscope::render::RGB_ORANGE);
    intEdgeQ->setRadius(0.0005);
}

/*
 * Display CornerData quantity on the common subdivision of an intrinsic triangulation.
 */
void visualizeIntrinsicSolution(const CornerData<double>& w, const std::string& name) {

    CommonSubdivision& cs = intTri->getCommonSubdivision();
    cs.constructMesh();
    // Linearly interpolate data from intrinsic mesh to the common subdivision.
    CornerData<double> cs_w = interpolateAcrossB(cs, w);

    ManifoldSurfaceMesh& csMesh = *cs.mesh;
    VertexData<Vector3> csPositions = cs.interpolateAcrossA(manifoldGeom->vertexPositions);

    // Register and display with Polyscope
    psCsMesh = polyscope::registerSurfaceMesh("common subdivision", csPositions, csMesh.getFaceVertexList(),
                                              polyscopePermutations(csMesh));
    psCsMesh->setEnabled(true);
    psCsMesh->addCornerScalarQuantity(name, cs_w)->setEnabled(true);
}

void exportCurves(int curveExportMode) {

    switch (LAST_SOLVE.solverMode) {
        case (SolverMode::ExtrinsicMesh): {
            switch (LAST_SOLVE.curveMode) {
                case (CurveMode::DUAL): {
                    // DUAL_CHAIN
                    break;
                }
                case (CurveMode::PRIMAL): {
                    // curveHalfedges
                    switch (curveExportMode) {
                        case (CurveExportMode::CLOSED_BOUNDING): {
                            std::vector<Halfedge> loops = getCompletedLoops(curveHalfedges, SWNSolver->wFunction);
                            std::vector<std::vector<Halfedge>> components = getCurveComponents(*geometry, loops);
                            exportCurvesAsOBJ(geometry->vertexPositions, components, DATA_DIR + "BoundingLoops.obj");
                            break;
                        }
                        case (CurveExportMode::CLOSED_NONBOUNDING): {
                            std::vector<Halfedge> loops = getCompletedLoops(curveHalfedges, SWNSolver->vFunction);
                            std::vector<std::vector<Halfedge>> components = getCurveComponents(*geometry, loops);
                            exportCurvesAsOBJ(geometry->vertexPositions, components, DATA_DIR + "NonBoundingLoops.obj");
                            break;
                        }
                        case (CurveExportMode::CLOSED_ALL): {
                            // TODO
                            break;
                        }
                        case (CurveExportMode::BOUNDING): {
                            // TODO
                            break;
                        }
                        case (CurveExportMode::NONBOUNDING): {
                            // TODO
                            break;
                        }
                    }
                    break;
                }
                case (CurveMode::PRIMAL_MANIFOLD): {
                    // curveHalfedgesOnManifold
                    break;
                }
                case (CurveMode::BARYCENTRIC): {
                    // CURVE_NODES, CURVE_EDGES
                    break;
                }
            }
            break;
        }
        case (SolverMode::IntrinsicMesh): {
            // curveHalfedgesOnIntrinsic;
            break;
        }
    }

    // switch (curveExportMode) {
    //     case (CurveExportMode::CLOSED_BOUNDING): {
    //         filename = DATA_DIR + "BoundingLoops.obj";
    //         break;
    //     }
    //     case (CurveExportMode::CLOSED_NONBOUNDING): {
    //         filename = DATA_DIR + "NonBoundingLoops.obj";
    //         break;
    //     }
    //     case (CurveExportMode::CLOSED_ALL): {
    //         filename = DATA_DIR + "CompletedCurve.obj";
    //         break;
    //     }
    //     case (CurveExportMode::BOUNDING): {
    //         filename = DATA_DIR + "BoundingParts.obj";
    //         break;
    //     }
    //     case (CurveExportMode::NONBOUNDING): {
    //         filename = DATA_DIR + "NonBoundingLoops.obj";
    //         break;
    //     }
    // }
}

void solve() {

    switch (SOLVER_MODE) {
        case (SolverMode::ExtrinsicMesh): {
            if (USING_GUI) displayCurves(getGeom(), getCurveHalfedges(), CURVE_NODES, CURVE_EDGES, DUAL_CHAIN);
            SWNSolver->doHomologyCorrection = DO_HOMOLOGY_CORRECTION;
            SWNSolver->approximateResidual = APPROXIMATE_RESIDUAL;
            SWNSolver->epsilon = EPSILON;
            LAST_SOLVE.solverMode = SolverMode::ExtrinsicMesh;
            CornerData<double> w;
            if (DUAL_CHAIN.size() > 0) {
                std::cerr << "m2" << std::endl;
                w = SWNSolver->solve(DUAL_CHAIN);
                LAST_SOLVE.curveMode = CurveMode::DUAL;
            } else if (!USING_MANIFOLD_MESH && curveHalfedges.size() > 0) {
                std::cerr << "m3" << std::endl;
                w = SWNSolver->solve(curveHalfedges);
                LAST_SOLVE.curveMode = CurveMode::PRIMAL;
            } else if (USING_MANIFOLD_MESH && curveHalfedgesOnManifold.size() > 0) {
                std::cerr << "m4" << std::endl;
                w = SWNSolver->solve(curveHalfedgesOnManifold);
                LAST_SOLVE.curveMode = CurveMode::PRIMAL_MANIFOLD;
            } else {
                std::cerr << "m5" << std::endl;
                w = SWNSolver->solve(CURVE_NODES, CURVE_EDGES);
                LAST_SOLVE.curveMode = CurveMode::BARYCENTRIC;
            }
            if (USING_GUI) psMesh->addCornerScalarQuantity("w", w)->setEnabled(true);
            Vector<double> wVector = w.toVector();
            std::cerr << "w min: " << wVector.minCoeff() << "\tw max: " << wVector.maxCoeff() << std::endl;
            break;
        }
        case (SolverMode::IntrinsicMesh): {
            ensureHaveIntrinsicSolver();
            resetCurveOnIntrinsicTriangulation();
            intrinsicSolver->doHomologyCorrection = DO_HOMOLOGY_CORRECTION;
            intrinsicSolver->approximateResidual = APPROXIMATE_RESIDUAL;
            intrinsicSolver->epsilon = EPSILON;
            LAST_SOLVE.solverMode = SolverMode::IntrinsicMesh;
            CornerData<double> w = intrinsicSolver->solve(curveHalfedgesOnIntrinsic);
            Vector<double> wVector = w.toVector();
            std::cerr << "w min: " << wVector.minCoeff() << "\tw max: " << wVector.maxCoeff() << std::endl;
            if (USING_GUI) visualizeIntrinsicSolution(w, "w");
            // transferBtoA(intTri, f_intrinsic, TransferMethod::L2);
            break;
        }
    }
}

void functionCallback() {

    if (ImGui::Button("Solve!")) {

        // remove everything except for input curves
        psMesh->removeAllQuantities();

        solve();
    }

    // Solve on original mesh;
    ImGui::RadioButton("Extrinsic mesh", &SOLVER_MODE, SolverMode::ExtrinsicMesh);
    // Solve on an intrinsic mesh
    ImGui::RadioButton("Intrinsic mesh", &SOLVER_MODE, SolverMode::IntrinsicMesh);

    if (ImGui::TreeNode("Advanced solving options")) {
        ImGui::Checkbox("Solve for nonbounding loops", &DO_HOMOLOGY_CORRECTION);
        ImGui::Checkbox("Use reduced-size linear program", &APPROXIMATE_RESIDUAL);
        ImGui::InputFloat("epsilon", &EPSILON);
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

    if (ImGui::TreeNode("Export")) {
        if (ImGui::TreeNode("Functions")) {
            if (ImGui::Button("Export winding number function")) {
                // TODO
            }
            if (ImGui::Button("Export initial jump harmonic function u")) {
                // TODO
            }
            if (ImGui::Button("Export residual function v")) {
                // if (SOLVER_MODE == SolverMode::ExtrinsicMesh && SWNSolver->vFunction.getMesh() != nullptr)
            }
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Curves")) {
            if (ImGui::Button("Export input curves")) {
                std::vector<std::vector<std::array<size_t, 2>>> connectedComponents =
                    getCurveComponents(*mesh, CURVE_NODES, CURVE_EDGES);
                exportCurvesAsOBJ(geometry->vertexPositions, CURVE_NODES, connectedComponents,
                                  DATA_DIR + "InputCurves.obj");
            }
            if (ImGui::Button("Export completed bounding loops")) {
                exportCurves(CurveExportMode::CLOSED_BOUNDING);
            }
            if (ImGui::Button("Export completed nonbounding loops")) {
                exportCurves(CurveExportMode::CLOSED_NONBOUNDING);
            }
            if (ImGui::Button("Export all completed curves")) {
                exportCurves(CurveExportMode::CLOSED_ALL);
            }
            if (ImGui::Button("Export bounding parts")) {
                exportCurves(CurveExportMode::BOUNDING);
            }
            if (ImGui::Button("Export nonbounding parts")) {
                exportCurves(CurveExportMode::NONBOUNDING);
            }
            ImGui::TreePop();
        }
        ImGui::TreePop();
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Solve for surface winding number on a triangle mesh.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
    args::ValueFlag<std::string> curveFilename(parser, "curve", "An input curve file", {"curve", "c"});
    args::ValueFlag<std::string> allowRemeshing(parser, "allowRemeshing", "Allow re-meshing of the input mesh.",
                                                {"allowRemeshing", "m"});
    args::ValueFlag<std::string> outputFilename(parser, "outputFilename", "File to save output mesh to.",
                                                {"outputFilename", "o"});
    args::ValueFlag<std::string> doHomologyCorrection(parser, "identifyNonbounding", "", {"identifyNonbounding", "h"});
    args::ValueFlag<std::string> approximateResidual(parser, "approximateResidual", "", {"approximateResidual", "r"});

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
    DATA_DIR = getHomeDirectory(MESH_FILEPATH); // extract home directory
    OUTPUT_FILENAME = DATA_DIR + "output.obj";
    std::tie(mesh, geometry) = readSurfaceMesh(MESH_FILEPATH);
    // Read line objects, if they exist in mesh file.
    readLines(*mesh, MESH_FILEPATH, CURVE_NODES, CURVE_EDGES, DUAL_CHAIN);

    // Read curve & other arguments.
    if (curveFilename) {
        std::string CURVE_FILEPATH = args::get(curveFilename);
        readCurves(*mesh, CURVE_FILEPATH, CURVE_NODES, CURVE_EDGES, DUAL_CHAIN);
    }
    if (allowRemeshing) {
        ALLOW_REMESHING = isStringTrue(args::get(allowRemeshing));
    }
    if (outputFilename) {
        OUTPUT_FILENAME = args::get(outputFilename);
    }
    if (doHomologyCorrection) {
        DO_HOMOLOGY_CORRECTION = isStringTrue(args::get(doHomologyCorrection));
    }
    if (approximateResidual) {
        APPROXIMATE_RESIDUAL = isStringTrue(args::get(approximateResidual));
    }

    // Initialize solver.
    SWNSolver = std::unique_ptr<SurfaceWindingNumbersSolver>(new SurfaceWindingNumbersSolver(*geometry));

    // Convert input curve to other forms, if applicable.
    if (CURVE_EDGES.size() > 0) {
        curveHalfedges = convertToHalfedges(CURVE_NODES, CURVE_EDGES);
        if (ALLOW_REMESHING && mesh->isManifold() && curveHalfedges.size() == 0) {
            if (remeshIfNecessary()) std::cerr << "Mesh re-meshed so curve is conforming." << std::endl;
        }
    }

    if (mesh->isManifold()) {
        reinterpretCurveToManifoldMesh();
        if (ALLOW_REMESHING && curveHalfedgesOnManifold.size() > 0 && splitEdgesIfNecessary(curveHalfedgesOnManifold)) {
            std::cerr << "Mesh re-meshed so curve contains no single-edge components." << std::endl;
        }
        ensureHaveIntrinsicSolver();
    }

    if (USING_GUI) {
        // Initialize polyscope
        polyscope::init();

        polyscope::state::userCallback = functionCallback;

        // Register the mesh with polyscope
        MESHROOT = polyscope::guessNiceNameFromPath(MESH_FILEPATH);
        psMesh =
            polyscope::registerSurfaceMesh(MESHNAME, getGeom().inputVertexPositions, getMesh().getFaceVertexList());
        psMesh->setAllPermutations(polyscopePermutations(getMesh()));

        // Display curve.
        displayCurves(getGeom(), getCurveHalfedges(), CURVE_NODES, CURVE_EDGES, DUAL_CHAIN);

        polyscope::show();
    }

    return EXIT_SUCCESS;
}