#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "ortools/linear_solver/linear_solver.h"
#include "utils.h"

#include <complex>
#include <queue>
#include <set>

// for debugging
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

// Implements the Surface Winding Numbers (SWN) algorithm from
//      Winding Number on Discrete Surfaces
//      Nicole Feng, Mark Gillespie, Keenan Crane
//      SIGGRAPH 2023

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace operations_research;

class SurfaceWindingNumbersSolver {

  public:
    // === Constructor
    SurfaceWindingNumbersSolver(IntrinsicGeometryInterface& geom);

    // === Solve

    // The only reason the solve() functions aren't const is because they update the intermediate quantities stored in
    // member variables.

    /* Curve(s) is specified as a primal 1-chain. Edges of the curve that pass through non-manifold vertices are
     * omitted. (If this is undesirable, consider using a dual 1-chain; see below.)
     */
    CornerData<double> solve(const Vector<double>& primalChain);

    /* Curve(s) is specified as an unordered collection of oriented mesh edges. */
    CornerData<double> solve(const std::vector<Halfedge>& curve);

    /* Curve is specified as an ordered sequence of mesh vertices. */
    CornerData<double> solve(const std::vector<Vertex>& curve);

    /* Multiple curves are specified as ordered sequences of vertices. */
    CornerData<double> solve(const std::vector<std::vector<Vertex>>& curves);

    /* Curve(s) is specified as a *dual* 1-chain, where each edge of the curve is specified as a pair of dual vertices
     * (mesh faces) indicating an oriented dual edge from the first face to the second.
     *
     * Specifying the curve as a dual 1-chain is the only method that achieves full generality, in that it can be used
     * to un-ambiguously specify the curve's sides and orientations on non-manifold or non-orientable meshes. */
    CornerData<double> solve(const std::vector<std::array<Face, 2>>& curve);

    void setLPSolver(const std::string& solverID);

    // === Parameters
    double epsilon = 1e-2;
    bool doHomologyCorrection = true;
    bool approximateResidual = false;
    bool verbose = false;

    // === Intermediate solve quantities
    CornerData<double> uFunction, vFunction, wFunction;
    EdgeData<double> gFunction; // g := Jv

  private:
    // === Members
    SurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

    bool simplyConnected;
    std::string LP_SOLVER = "SCIP";
    HalfedgeData<double> halfedgeCotanWeights;

    // === Solvers
    // the 2-form Laplacian used to get co-exact component of a primal 1-form
    std::unique_ptr<PositiveDefiniteSolver<double>> coexactSolver;

    // === DEC operators used throughout
    SparseMatrix<double> d0, d0T, hodge1, hodge1Inv, d1, d1T;

    // === Algorithm steps
    CornerData<double> computeReducedCoordinates(const Vector<double>& chain,
                                                 const std::vector<Vertex>& interiorVertices,
                                                 const std::map<Vertex, Halfedge>& outgoingHalfedgeOnCurve) const;

    CornerData<double> solveJumpLaplaceEquation(const std::vector<Vertex>& interiorVertices,
                                                const VertexData<bool>& isInteriorEndpoint,
                                                const CornerData<double>& reducedCoordinates) const;

    Vector<double> DarbouxDerivative(const VertexData<bool>& isInteriorEndpoint, const CornerData<double>& u) const;

    SparseMatrix<double> buildLaplacian(const VertexData<bool>& isInteriorEndpoint, const VertexData<size_t>& DOFindex,
                                        const size_t& nDOFs) const;

    Vector<double> buildJumpLaplaceRHS(const std::vector<Vertex>& interiorVertices,
                                       const VertexData<bool>& isInteriorEndpoint,
                                       const CornerData<double>& reducedCoordinates, const VertexData<size_t>& DOFindex,
                                       const size_t& nDOFs) const;

    Vector<double> harmonicComponent(const Vector<double>& omega) const;

    CornerData<double> residualFunction(const Vector<double>& chain, const Vector<double>& gamma) const;

    CornerData<double> integrateLocally(const Vector<double>& gamma) const;

    CornerData<double> solveLinearProgram(const Vector<double>& chain, const CornerData<double>& vInit) const;

    CornerData<double> approximateResidualFunction(const Vector<double>& chain,
                                                   const std::vector<std::pair<Vertex, bool>>& endpoints,
                                                   const Vector<double>& gamma) const;

    Vector<double> dijkstraCompleteCurve(const Vector<double>& chain,
                                         const std::vector<std::pair<Vertex, bool>>& curveEndpoints) const;

    std::vector<Halfedge> dijkstraPath(IntrinsicGeometryInterface& geom, const Vertex& startVert,
                                       const std::set<Vertex>& endVerts) const;

    CornerData<double> subtractJumpDerivative(const Vector<double>& chain, const std::vector<Vertex>& interiorVertices,
                                              const VertexData<bool>& isInteriorEndpoint,
                                              const std::map<Vertex, Halfedge>& outgoingHalfedgeOnCurve,
                                              const CornerData<double>& resid) const;

    // === Auxiliary functions

    void ensureHaveCoexactSolver();

    Vector<double> computeCoExactComponent(const Vector<double>& omega) const;

    // === Utility functions

    Vector<double> convertToChain(const std::vector<Halfedge>& curve) const;

    Vector<double> convertToChain(const std::vector<Vertex>& curve) const;

    Vector<double> jumpDerivative(const CornerData<double>& f) const;
};
