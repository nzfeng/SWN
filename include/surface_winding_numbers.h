#pragma once

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/barycentric_vector.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "utils.h"

#include <complex>
#include <set>

#include <CoMISo/Config/config.hh>
#include <CoMISo/NSolver/CBCSolver.hh>
#include <CoMISo/NSolver/COMISOSolver.hh>
#include <CoMISo/NSolver/CPLEXSolver.hh>
#include <CoMISo/NSolver/GUROBISolver.hh>
// #include <CoMISo/NSolver/LPSolveSolver.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/LinearProblem.hh>
#include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/NSolver/NProblemInterface.hh>
#include <CoMISo/NSolver/VariableType.hh>
#include <CoMISo/Utils/StopWatch.hh>

// Implements the Surface Winding Numbers (SWN) algorithm from
//      Winding Number on Discrete Surfaces
//      Nicole Feng, Mark Gillespie, Keenan Crane
//      SIGGRAPH 2023

using namespace geometrycentral;
using namespace geometrycentral::surface;

class SurfaceWindingNumbersSolver {

  public:
    // === Constructor
    SurfaceWindingNumbersSolver(IntrinsicGeometryInterface& geom);

    // === Solve
    CornerData<double> solve(const std::vector<Halfedge>& curve, bool doHomologyCorrection = true,
                             bool useSpecialBases = true);

    VertexData<double> solve(const std::vector<BarycentricVector>& curve, bool doHomologyCorrection = true,
                             bool useWhitneyElements = false);

    VertexData<double> solve(const std::vector<SurfacePoint>& curveNodes,
                             const std::vector<std::array<size_t, 2>>& curveEdges, bool doHomologyCorrection = true,
                             bool useWhitneyElements = false);

    // === Edit parameters
    void setEpsilon(double epsilon);

  private:
    // === Members
    SurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

    double EPSILON = 1e-2;
    double SPECIAL_VAL = std::numeric_limits<double>::quiet_NaN();
    bool USE_SPECIAL_BASES = true;

    // === Curve data
    std::vector<Halfedge> curveHalfedges;        // the curve, represented as a set of halfedges
    std::vector<BarycentricVector> curveVectors; // the curve, represented as a set of BarycentricVectors

    // === Solvers
    // factors the standard cotan Laplacian
    std::unique_ptr<PositiveDefiniteSolver<double>> cotanLaplaceSolver;
    // the 2-form Laplacian used to get co-exact component of a primal 1-form
    std::unique_ptr<PositiveDefiniteSolver<double>> coexactSolver;
    std::unique_ptr<PositiveDefiniteSolver<double>> crouzeixRaviartPoissonSolver; // use Crouzeix-Raviart elements
    std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;                // use Whitney elements

    // === DEC operators used throughout
    SparseMatrix<double> d0, d0T, hodge1, hodge1Inv, d1, d1T;
    HalfedgeData<double> halfedgeCotanWeights;

    // === Helpers, utility / auxiliary functions
    Vector<double> shiftVector(const Vector<double>& v) const;
    Eigen::SparseMatrix<double, Eigen::RowMajor> convertToRowMajor(const SparseMatrix<double>& M) const;
    void shiftWedgeVectorByAverageValueAlongCurves(Vector<double>& wedgeVals,
                                                   const std::vector<Eigen::Triplet<double>>& jumpPairs) const;
    void shiftToGetIntegers(Vector<double>& u_wedges, const std::vector<Eigen::Triplet<double>>& jumpPairs) const;
    Vector<double> shiftToGetIntegers(const size_t& N, const CornerData<double>& u_corners,
                                      const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                      const CornerData<size_t>& wedgeIndices) const;
    double averageJumpValueAlongAllCurves(Vector<double>& wedgeVals,
                                          const std::vector<Eigen::Triplet<double>>& jumpPairs) const;
    double interpolateCornerData(const CornerData<double>& cornerVals, Face f) const;
    bool isValidWedgeIndex(size_t wIdx) const;
    bool isCurveClosed(const std::vector<Halfedge>& curve);
    Vector<double> convertToChain(const std::vector<Halfedge>& curve);
    Vector<double> cornerVector(const Vector<double>& v, const CornerData<size_t>& wedgeIndices);
    std::vector<BarycentricVector> gradient(const Vector<double>& omega);

    std::vector<double> computeJumpLengths(const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                           const std::vector<Eigen::Triplet<double>>& jumpPairsDedup,
                                           const CornerData<size_t>& wedgeIndices,
                                           const std::vector<Vertex>& bVertices);

    std::vector<Halfedge> dijkstraCancelBoundary(const std::vector<Halfedge>& curve);
    std::vector<Halfedge> dijkstraCompleteCurve(const std::vector<Halfedge>& curve);

    std::vector<Halfedge> completeCurve(const std::vector<Halfedge>& curve, const Vector<double>& gamma);

    void ensureHaveCotanLaplaceSolver();
    void ensureHaveHodgeDecompositionSolvers();
    Vector<double> computeExactComponent(const Vector<double>& omega);
    Vector<double> computeCoExactComponent(const Vector<double>& omega);

    VertexData<double> edgeMidpointDataToVertexData(const Vector<double>& u);
    std::vector<BarycentricVector>
    convertCurveToBarycentricVectors(const std::vector<SurfacePoint>& curveNodes,
                                     const std::vector<std::array<size_t, 2>>& curveEdges) const;

    // void ensureHaveLaplaceSolver();
    void ensureHaveCrouzeixRaviartPoissonSolver();
    void ensureHavePoissonSolver();

    // === Functions for solving the Poisson equation
    void discretizeSegmentInWhitneyBasis(const BarycentricVector& seg, Vector<double>& N_Gamma);
    void discretizeSegmentInCrouzeixRaviartBasis(const BarycentricVector& seg, Vector<double>& N_Gamma);

    Vector<double> discretizeRHS(bool useWhitneyElements);
    // Solve on input mesh using Whitney elements.
    VertexData<double> solvePoissonEquation();
    // Solve on input using Crouzei-Raviart elements.
    VertexData<double> solveCrouzeixRaviartPoissonEquation();

    // === Functions for solving the jump equation
    CornerData<size_t> indexWedges(const std::vector<Halfedge>& curve, size_t& iN, size_t& N,
                                   std::vector<Eigen::Triplet<double>>& wedgePairs, std::vector<Vertex>& bVertices,
                                   std::vector<Vertex>& curveEndpoints, bool useSpecialBases = true);
    std::vector<Eigen::Triplet<double>> computeCumulativeJumps(const std::vector<Halfedge>& curve,
                                                               const std::vector<Vertex>& bVertices,
                                                               const CornerData<size_t>& wedgeIndices);

    CornerData<double> solveJumpEquation(const size_t& N, const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                         const std::vector<Eigen::Triplet<double>>& qJumpPairs,
                                         const CornerData<size_t>& wedgeIndices,
                                         const std::vector<Vertex>& curveEndpoints);

    Vector<double> computeOmega(const CornerData<double>& u, const std::vector<Vertex>& curveEndpoints);

    void computeOmegaAndHodgeDecomposition(const CornerData<double>& u, const std::vector<Vertex>& curveEndpoints,
                                           Vector<double>& omega, Vector<double>& dAlpha, Vector<double>& deltaBeta,
                                           Vector<double>& gamma);
    std::vector<Eigen::Triplet<double>> getFinalJumps(const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                                      const CornerData<double>& v_corners,
                                                      const CornerData<size_t>& wedgeIndices, size_t N) const;

    Vector<double> computeIntegratedDivergenceOverBoundaryCells(const Vector<double>& gamma,
                                                                const std::vector<Vertex>& bVertices,
                                                                const CornerData<size_t>& wedgeIndices, size_t N) const;
    Vector<double> computeIntegratedDivergence(const Vector<double>& omega, const CornerData<size_t>& wedgeIndices,
                                               size_t N) const;
    Vector<double> computeIntegratedDivergence(const HalfedgeData<double>& omega,
                                               const CornerData<size_t>& wedgeIndices, size_t N) const;

    CornerData<double> integrateExactly(const Vector<double>& gamma, std::vector<Halfedge>& branchCut);

    std::vector<std::pair<size_t, double>> matchUpJumpPairs(const std::vector<Eigen::Triplet<double>>& jumpPairs_u,
                                                            const CornerData<size_t>& wedgeIndices_u, const size_t& N_u,
                                                            const std::vector<Eigen::Triplet<double>>& jumpPairs_v,
                                                            const CornerData<size_t>& wedgeIndices_v,
                                                            const size_t& N_v);

    Vector<double> integrateGammaAndMinimizeJumpsL1(const std::vector<Halfedge>& completedCurve,
                                                    const Vector<double>& gamma, CornerData<size_t>& wedgeIndices_v);

    CornerData<double> integrateExactlyBranchCutAnywhere(const Vector<double>& gamma,
                                                         const std::vector<Halfedge>& curve);

    CornerData<double> reducedLinearProgram(const Vector<double>& gamma, const std::vector<Halfedge>& curve);

    CornerData<double> integrateExactlyPerFace(const Vector<double>& gamma, const FaceData<double>& shiftPerFace);
};
