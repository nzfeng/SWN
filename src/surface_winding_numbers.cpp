#include "surface_winding_numbers.h"

// ==== SOLVE

/*
 * Input: A discrete dual 1-chain ∈ Z^|E| encoding the input curve Γ.
 *
 * Output: The winding number function on corners. If useSpecialBases = true, corners adjacent to interior endpoints
 * will have placeholder values SPECIAL_VAL; their values are to be interpolated using the scheme described in
 * Section 2.3.2.
 */
CornerData<double> SurfaceWindingNumbersSolver::solve(const Vector<int>& dualChain, bool doHomologyCorrection,
                                                      bool useSpecialBases) const {
    USE_SPECIAL_BASES = useSpecialBases;
    // TODO: compute interior endpoints, boundary = VertexData<int>, VertexData<bool> isInteriorEndpoint

    // DEC operators

    // Determine whether the mesh is simply-connected.
    bool simplyConnected = false;
    if (mesh.isManifold() && mesh.nConnectedComponents() == 1) {
        size_t chi = mesh.nVertices() + mesh.nFaces() - mesh.nEdges();
        simplyConnected = (chi == 1 || chi == 2);
    }

    computeReducedCoordinates();
    CornerData<double> w = solveJumpEquation();

    if (!simplyConnected && doHomologyCorrection) {
        Vector<double> gamma = DarbouxDerivative(isInteriorEndpoint);
        if (isCurveClosed(dualChain)) gamma = harmonicComponent(gamma);
        residualFunction(); // IntegrateLocally(), ComputeRelativeJumps(), SolveLinearProgram(), RecoverSolution();
        subtractJumpDerivative();
        CornerData<double> w = solveJumpEquation();
    }
    return w;
}

CornerData<double> SurfaceWindingNumbersSolver::solve(const std::vector<Halfedge>& curve, bool doHomologyCorrection,
                                                      bool useSpecialBases) const {

    // TODO: Convert the input curve to a dual 1-chain, then call generic solver.
}


// ==== ALGORITHM STEPS


/*
 * Input: The curve Γ, represented by relevant pre-computed quantities (its interior vertices).
 *
 * Output: A function c on corners, that expresses values at corners relative to a reference value at a corner adjacent
 * to the same vertex.
 */
CornerData<double>
SurfaceWindingNumbersSolver::computeReducedCoordinates(const std::vector<Vertex>& interiorVertices) const {

    CornerData<double> reducedCoordinates(mesh, 0);
    for (const auto& vi : interiorVertices) {
        if (!vi.isManifold()) continue;
        // TODO
    }
    return reducedCoordinates;
}

/*
 * Input: The curve Γ, represented by relevant pre-computed quantities: its interior vertices, interior endpoints, and
 * reduced coordinates.
 *
 * Output: A function u. If doing projective interpolation (Section 2.3.2), values at corners adjacent to interior
 * endpoints are left undefined, to be interpolated later.
 */
CornerData<double> SurfaceWindingNumbersSolver::solveJumpEquation(const std::vector<Vertex>& interiorVertices,
                                                                  const VertexData<bool>& isInteriorEndpoint,
                                                                  const CornerData<double>& reducedCoordinates) const {

    SparseMatrix<double> L = buildLaplacian(isInteriorEndpoint);
    Vector<double> b = buildJumpLaplaceRHS(interiorVertices, isInteriorEndpoint, reducedCoordinates);
    shiftDiagonal(L, 1e-8); // hack to ensure L is PD and not just PSD
    Vector<double> u0 = solvePositiveDefinite(L, b);
    // Apply shifts to recover u.
    CornerData<double> u(mesh);
    geom.requireVertexIndices();
    for (Vertex v : mesh.vertices()) {
        size_t vIdx = geom.vertexIndices[v];
        for (Corner c : v.adjacentCorners()) u[c] = u0[vIdx] + reducedCoordinates[c];
    }
    geom.unrequireVertexIndices();
    return u;
}

/*
 * Input: The interior endpoints of a curve Γ.
 *
 * Output: The Darboux derivative Du = ω ∈ R^|E|, represented as a discrete primal 1-form. Values of ω at edges incident
 * on interior endpoints are 0 by definition (Section 2.4.2.)
 */
Vector<double> SurfaceWindingNumbersSolver::DarbouxDerivative(const VertexData<bool>& isInteriorEndpoint) const {

    geom.requireEdgeIndices();
    Vector<double> omega = Vector<double>::Zero(mesh.nEdges());

    for (Edge e : mesh.edges()) {

        // For every edge connecting to a vertex at a curve endpoint, set the value to zero.
        if (isInteriorEndpoint[e.firstVertex()] || isInteriorEndpoint[e.secondVertex()]) continue;

        size_t i = geom.edgeIndices[e];
        Halfedge he = e.halfedge();
        Corner c1 = he.corner();
        Corner c2 = he.next().corner();
        Corner d2 = c2;
        Corner d1 = c1;
        if (!e.isBoundary()) {
            d2 = he.twin().corner();
            d1 = he.twin().next().corner();
        }
        omega[i] = 0.5 * ((u[c2] - u[c1]) + (u[d2] - u[d1]));
    }
    geom.unrequireEdgeIndices();
    return omega;
}

/*
 * Input: The interior endpoints of a curve Γ.
 *
 * Output: The standard cotan Laplacian on V^*, the set of vertices minus interior endpoints.
 */
SparseMatrix<double> SurfaceWindingNumbersSolver::buildLaplacian(const VertexData<bool>& isInteriorEndpoint) const {

    // Omit curve endpoints from the system. First map existing vertices to their new DOF indices.
    VertexData<size_t> DOFindex(mesh);
    VertexData<bool> isCurveEndpoint(mesh, false);
    for (const Vertex& v : curveEndpoints) {
        isCurveEndpoint[v] = true;
    }
    size_t nDOF = 0;
    for (Vertex v : mesh.vertices()) {
        if (!isCurveEndpoint[v]) {
            DOFindex[v] = nDOF;
            nDOF++;
        }
    }

    // Build Laplace matrix.
    FaceData<bool> skippedFace(mesh, false);
    HalfedgeData<bool> skipHalfedge(mesh, false);
    SparseMatrix<double> L(nDOF, nDOF);
    std::vector<Eigen::Triplet<double>> tripletList;
    geom.requireHalfedgeCotanWeights();
    for (Face f : mesh.faces()) {

        // bool skipFace = false;
        // for (Vertex v : f.adjacentVertices()) {
        //     if (isCurveEndpoint[v]) {
        //         skipFace = true;
        //         skippedFace[f] = true;
        //         for (Halfedge he : f.adjacentHalfedges()) skipHalfedge[he] = true;
        //         break;
        //     }
        // }
        // if (skipFace) continue;

        for (Halfedge he : f.adjacentHalfedges()) {
            if (isCurveEndpoint[he.tailVertex()] || isCurveEndpoint[he.tipVertex()]) {
                skipHalfedge[he] = true;
                continue;
            }

            size_t i = DOFindex[he.tailVertex()];
            size_t j = DOFindex[he.tipVertex()];
            double w = halfedgeCotanWeights[he];
            if (std::isnan(w) || !std::isfinite(w)) w = 1;
            tripletList.emplace_back(i, i, w);
            tripletList.emplace_back(j, j, w);
            tripletList.emplace_back(i, j, -w);
            tripletList.emplace_back(j, i, -w);
        }
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
    return L;
}

/*
 * Input: The interior vertices and interior endpoints of a curve Γ.
 *
 * Output: The r.h.s. |V^*|-vector encoding the jump conditions for the jump Laplace equation.
 */
SurfaceWindingNumbersSolver::buildJumpLaplaceRHS(const std::vector<Vertex>& interiorVertices,
                                                 const VertexData<bool>& isInteriorEndpoint) const {

    // Build RHS (with DOFs at curve endpoints omitted)
    Vector<double> RHS = Vector<double>::Zero(nDOF);
    size_t m = qJumpPairs.size();
    for (size_t i = 0; i < m; i++) {
        size_t w_i = qJumpPairs[i].row();
        size_t w_j = qJumpPairs[i].col();
        double jumpVal = qJumpPairs[i].value();
        size_t wIdx = w_j;
        Vertex v = wedgeToVertex[wIdx];
        size_t vIdx = DOFindex[v];
        for (Corner c : v.adjacentCorners()) {
            if (wedgeIndices[c] == wIdx) {
                Halfedge heA = c.halfedge();
                Halfedge heB = heA.next().next();

                if (!skipHalfedge[heA]) {
                    RHS[vIdx] -= halfedgeCotanWeights[heA] * jumpVal;
                    RHS[DOFindex[heA.tipVertex()]] += halfedgeCotanWeights[heA] * jumpVal;
                }
                if (!skipHalfedge[heB]) {
                    RHS[vIdx] -= halfedgeCotanWeights[heB] * jumpVal;
                    RHS[DOFindex[heB.tailVertex()]] += halfedgeCotanWeights[heB] * jumpVal;
                }
            }
        }
    }
}

/*
 * Input: A discrete primal 1-form ω ∈ R^|E| with no exact component.
 *
 * Output: The harmonic component γ ∈ R^|E| of ω, also a primal 1-form.
 */
Vector<double> SurfaceWindingNumbersSolver::harmonicComponent(const Vector<double>& omega) const {

    ensureHaveCoexactSolver();
    Vector<double> deltaBeta = computeCoExactComponent(omega);
    Vector<double> gamma = omega - deltaBeta;
    return gamma;
}

SurfaceWindingNumbersSolver::residualFunction() const {
    integrateLocally();
    computeRelativeJumps();
    recoverSolution();
}

SurfaceWindingNumbersSolver::integrateLocally() const {}

SurfaceWindingNumbersSolver::computeRelativeJumps() const {}

SurfaceWindingNumbersSolver::recoverSolution() const {}

/*
 * Input: The interior endpoints of curve Γ, residual function v on corners, and reduced coordinates associated with Γ.
 *
 * Output: Updated reduced coordinates encoding new jump constraints for the jump Laplace equation.
 */
SurfaceWindingNumbersSolver::subtractJumpDerivative() const {}


// ==== AUXILIARY FUNCTIONS


/*
 * Build the solvers for performing Hodge decomposition on a primal 1-form.
 *
 * The output of this function only remains valid if the mesh remains un-mutated! (which is ensured since the mesh, etc.
 * variables are private; if these change, the user must re-construct a new SurfaceWindingNumbersSolver.)
 */
void SurfaceWindingNumbersSolver::ensureHaveCoexactSolver() {

    if (coexactSolver == nullptr) {
        geom.requireDECOperators();
        d0 = geom.d0;
        d0T = d0.transpose();
        d1 = geom.d1;
        d1T = d1.transpose();
        hodge1 = geom.hodge1;
        hodge1Inv = geom.hodge1Inverse;
        geom.requireEdgeIndices();
        const EdgeData<size_t>& eIdx = geom.edgeIndices;
        for (Edge e : mesh.edges()) {
            double wInv = halfedgeCotanWeights[e.halfedge()];
            if (!e.isBoundary()) wInv += halfedgeCotanWeights[e.halfedge().twin()];
            double w = 1. / wInv;
            if (!std::isfinite(w) || std::isnan(w) || w < 0) {
                std::cout << "bad w: " << w << vendl;
                w = 1;
            } else if (w > 1000) {
                std::cout << "really big w: " << w << vendl;
                w = 1000;
            }
            hodge1Inv.coeffRef(eIdx[e], eIdx[e]) = w;
        }
        geom.unrequireEdgeIndices();
        SparseMatrix<double> B = d1 * hodge1Inv * d1T;
        shiftDiagonal(B, 1e-8);
        coexactSolver.reset(new PositiveDefiniteSolver<double>(B));
        geom.unrequireDECOperators();
    }
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Returns the coexact component 𝛿β of ω, as a primal 1-form.
 */
Vector<double> SurfaceWindingNumbersSolver::computeCoExactComponent(const Vector<double>& omega) {

    Vector<double> rhs = d1 * omega;
    Vector<double> betaTilde = coexactSolver->solve(rhs);
    return hodge1Inv * d1T * betaTilde;
}


// ==== UTILITIES

// ==== CONSTRUCTOR


SurfaceWindingNumbersSolver::SurfaceWindingNumbersSolver(IntrinsicGeometryInterface& geom_)
    : mesh(geom_.mesh), geom(geom_) {

    if (!mesh.isTriangular()) throw std::logic_error("Mesh must be triangular to run SWN.");
}


/*
 * Shift vector so that its minimum value is 0.
 */
Vector<double> SurfaceWindingNumbersSolver::shiftVector(const Vector<double>& v) const {

    Vector<double> u = v;
    double shift = v.minCoeff();
    for (int i = 0; i < v.size(); i++) {
        u[i] -= shift;
    }
    return u;
}

/*
 * Convert an Eigen SparseMatrix in column-major storage order to row-major.
 */
Eigen::SparseMatrix<double, Eigen::RowMajor>
SurfaceWindingNumbersSolver::convertToRowMajor(const SparseMatrix<double>& M) const {

    size_t nRows = M.rows();
    size_t nCols = M.cols();
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(nRows, nCols);
    std::vector<Eigen::Triplet<double>> tripletList;
    // For each pair of adjacent corners, add their cotan weights to the matrix.
    for (size_t i = 0; i < nCols; i++) {
        for (SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
            tripletList.emplace_back(it.row(), it.col(), it.value());
        }
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

/*
 * Compute a shift vector that will result in the input curve being contoured at a levelset that partitions the domain
 * into roughly integer values. (Don't necessarily want to shift so that the input curve gets contoured at
 * half-integers; consider two duplicate contractible curves.)
 */
void SurfaceWindingNumbersSolver::shiftToGetIntegers(Vector<double>& u_wedges,
                                                     const std::vector<Eigen::Triplet<double>>& jumpPairs) const {
    double avgValOnOneSide = 0.;
    size_t nValidJumps = 0;
    size_t m = jumpPairs.size();
    for (size_t i = 0; i < m; i++) {
        size_t w_i = jumpPairs[i].row();
        size_t w_j = jumpPairs[i].col();
        double jumpVal = jumpPairs[i].value();
        if (abs(jumpVal) > 0.5) {
            avgValOnOneSide += (jumpVal > 0) ? u_wedges[w_i] : u_wedges[w_j]; // take outside
            nValidJumps += 1;
        }
    }
    if (nValidJumps == 0) return;
    avgValOnOneSide /= nValidJumps;

    // Doesn't matter than we shift to the nearest int. Just need to shift to some int.
    double someInt = roundToNearestInt(avgValOnOneSide);
    double shift = someInt - avgValOnOneSide;
    u_wedges += shift * Vector<double>::Ones(u_wedges.size());
}

Vector<double> SurfaceWindingNumbersSolver::shiftToGetIntegers(const size_t& N, const CornerData<double>& u_corners,
                                                               const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                                               const CornerData<size_t>& wedgeIndices) const {

    Vector<double> u_wedges(N);
    for (Corner c : mesh.corners()) {
        size_t wIdx = wedgeIndices[c];
        if (isValidWedgeIndex(wIdx)) u_wedges[wIdx] = u_corners[c];
    }
    shiftToGetIntegers(u_wedges, jumpPairs);
    return u_wedges;
}

/*
 * Shift corner data by the average value along the curves. That way, curves will be at approximately the 0.5-levelset.
 */
void SurfaceWindingNumbersSolver::shiftWedgeVectorByAverageValueAlongCurves(
    Vector<double>& wedgeVals, const std::vector<Eigen::Triplet<double>>& jumpPairs) const {

    double shift = 0.;
    size_t m = jumpPairs.size();
    for (size_t i = 0; i < m; i++) {
        size_t w_i = jumpPairs[i].row();
        size_t w_j = jumpPairs[i].col();
        shift += 0.5 * (wedgeVals[w_i] + wedgeVals[w_j]);
    }
    shift /= m;
    wedgeVals -= Vector<double>::Ones(wedgeVals.size()) * shift;
}

/*
 * Compute average jump value along *all* curves.
 * Warning: This does not (yet) compute the average jump value for each connected component!
 */
double SurfaceWindingNumbersSolver::averageJumpValueAlongAllCurves(
    Vector<double>& wedgeVals, const std::vector<Eigen::Triplet<double>>& jumpPairs) const {

    double avg = 0.;
    size_t m = jumpPairs.size();
    for (size_t i = 0; i < m; i++) {
        size_t w_i = jumpPairs[i].row();
        size_t w_j = jumpPairs[i].col();
        avg += abs(wedgeVals[w_i] - wedgeVals[w_j]);
    }
    avg /= m;
    return avg;
}

/*
 * Interpolate the given scalar corner data at the barycenter of the given face.
 */
double SurfaceWindingNumbersSolver::interpolateCornerData(const CornerData<double>& cornerVals, Face f) const {

    double val = 0.;
    for (Corner c : f.adjacentCorners()) {
        val += cornerVals[c];
    }
    return val / f.degree();
}

bool SurfaceWindingNumbersSolver::isValidWedgeIndex(size_t wIdx) const {
    return wIdx < mesh.nCorners();
}

bool SurfaceWindingNumbersSolver::isCurveClosed(const std::vector<Halfedge>& curve) {

    geom.requireDECOperators();
    SparseMatrix<double> B = geom.d0.transpose();
    geom.unrequireDECOperators();

    geom.requireEdgeIndices();
    Vector<double> chain = Vector<double>::Zero(mesh.nEdges());
    for (const Halfedge& he : curve) {
        chain[geom.edgeIndices[he.edge()]] += (he.orientation()) ? 1 : -1;
    }
    geom.unrequireEdgeIndices();

    double eps = 1e-5;
    Vector<double> boundary = B * chain;
    if (abs(boundary.norm()) < eps) return true;

    // Check if the curve is closed relative to the boundary.
    geom.requireVertexIndices();
    for (Vertex v : mesh.vertices()) {
        size_t vIdx = geom.vertexIndices[v];
        if (abs(boundary[vIdx]) > eps && !v.isBoundary()) return false;
    }
    geom.unrequireVertexIndices();
    return true;
}

Vector<double> SurfaceWindingNumbersSolver::convertToChain(const std::vector<Halfedge>& curve) {

    Vector<double> chain = Vector<double>::Zero(mesh.nEdges());
    geom.requireEdgeIndices();
    for (const Halfedge& he : curve) {
        size_t eIdx = geom.edgeIndices[he.edge()];
        chain[eIdx] += he.orientation() ? 1 : -1;
    }
    geom.unrequireEdgeIndices();
    return chain;
}

/*
 * Convert a vector of values on wedges to a vector of values on corners.
 * Use a special value
 */
Vector<double> SurfaceWindingNumbersSolver::cornerVector(const Vector<double>& v,
                                                         const CornerData<size_t>& wedgeIndices) {

    geom.requireCornerIndices();
    Vector<double> u(mesh.nCorners());
    for (Corner c : mesh.corners()) {
        size_t wIdx = wedgeIndices[c];
        size_t cIdx = geom.cornerIndices[c];
        if (isValidWedgeIndex(wIdx)) {
            u[cIdx] = v[wIdx];
        } else {
            u[cIdx] = SPECIAL_VAL;
        }
    }
    geom.unrequireCornerIndices();
    return u;
}

/*
 * Given a 1-form ω, return the gradient vector (a single vector, since the gradient of a linear function in 2D is
 * constant per face), corresponding to ω restricted to each face. If ω is harmonic, we can in fact represent it
 * via a single vector per face on the entire mesh.
 */
std::vector<BarycentricVector> SurfaceWindingNumbersSolver::gradient(const Vector<double>& omega) {

    std::vector<BarycentricVector> grad(mesh.nFaces());

    // Whitney-interpolate ω at the barycenter of each face.
    geom.requireFaceAreas();
    // Store barycentric vectors representing the halfedges of each face
    std::vector<Vector3> heVecs = {Vector3{-1., 1., 0.}, Vector3{0., -1., 1.}, Vector3{1., -1., 0.}};
    for (Face f : mesh.faces()) {
        size_t fIdx = geom.faceIndices[f];
        double A = geom.faceAreas[f];
        // gradient vectors of the piecewise hat functions (Whitney 0-forms), except the 1/2A factor
        std::vector<BarycentricVector> vecs;
        // 1-form values around the current face
        std::vector<double> oneFormVals;
        int k = 0;
        for (Halfedge he : f.adjacentHalfedges()) {
            // vi = he.tailVertex(); vj = he.tipVertex();
            size_t eIdx = geom.edgeIndices[he.edge()];
            double val = he.orientation() ? omega[eIdx] : -omega[eIdx];
            oneFormVals.push_back(val);
            BarycentricVector heVec = BarycentricVector(f, heVecs[k]);
            heVec = heVec.rotated90(geom);
            vecs.push_back(heVec);
            k++;
        }
        BarycentricVector g = BarycentricVector(f);
        for (int i = 0; i < 3; i++) {
            g += (oneFormVals[(i + 1) % 3] + oneFormVals[(i + 2) % 3]) * vecs[i];
        }
        g /= 6. * A;

        grad[fIdx] = g;
    }
    geom.unrequireFaceAreas();

    return grad;
}

/*
 * Determine the "length" of the jump between each jump (wedge) pair.
 */
std::vector<double> SurfaceWindingNumbersSolver::computeJumpLengths(
    const std::vector<Eigen::Triplet<double>>& jumpPairs, const std::vector<Eigen::Triplet<double>>& jumpPairsDedup,
    const CornerData<size_t>& wedgeIndices, const std::vector<Vertex>& bVertices) {

    std::map<std::pair<size_t, size_t>, double> lengthMap;
    geom.requireEdgeLengths();
    for (const Vertex& v : bVertices) {
        for (Halfedge he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) continue;
            Corner cB = he.corner();               // inside corner
            Corner cA = he.twin().next().corner(); // outside corner
            size_t wA = wedgeIndices[cA];
            size_t wB = wedgeIndices[cB];
            if (wA != wB) {
                double length = geom.edgeLengths[he.edge()];
                std::pair<size_t, size_t> keyA = std::make_pair(wA, wB);
                std::pair<size_t, size_t> keyB = std::make_pair(wB, wA);
                auto searchA = lengthMap.find(keyA);
                auto searchB = lengthMap.find(keyB);
                if (searchA == lengthMap.end()) {
                    lengthMap[keyA] = 0.5 * length;
                } else {
                    searchA->second += 0.5 * length;
                }
                if (searchB == lengthMap.end()) {
                    lengthMap[keyB] = 0.5 * length;
                } else {
                    searchB->second += 0.5 * length;
                }
            }
        }
    }
    geom.unrequireEdgeLengths();

    // Match up lengths to jumpPairsDedup
    size_t m = jumpPairsDedup.size();
    std::vector<double> jumpLengths(m, 0.);
    for (size_t i = 0; i < m; i++) {
        const Eigen::Triplet<double>& jumpPair = jumpPairsDedup[i];
        std::pair<size_t, size_t> key = std::make_pair(jumpPair.row(), jumpPair.col());
        jumpLengths[i] = lengthMap[key];
    }
    return jumpLengths;
}

void SurfaceWindingNumbersSolver::ensureHaveCotanLaplaceSolver() {

    if (cotanLaplaceSolver != nullptr) return;

    geom.requireCotanLaplacian();
    SparseMatrix<double>& C = geom.cotanLaplacian;
    shiftDiagonal(C, 1e-8);
    cotanLaplaceSolver.reset(new PositiveDefiniteSolver<double>(C));
    geom.unrequireCotanLaplacian();
}


// ==== JUMP EQUATION

/*
 * Index wedges, where a "wedge" is a union of corners. The first <iN> indices will belong to interior wedges, and
 * the next <bN> will belong to boundary wedges. Interior wedges are those that are adjacent to an interior vertex,
 * where "interior" refers to vertices that are not on the curve. Curve endpoints, however, still count as
 * "interior".
 *
 * This function returns a map from corners to wedge indices.
 *
 * <iN> and <N> keep track of the # of interior and total wedges, respectively.
 *
 * <wedgePairs> keeps track of which wedges are "paired". This information is used to construct the constraints in
 * the linear system later. The value of the winding number function should jump by the corresponding value from the
 * first wedge to the second wedge in the pair.
 *
 * <bVertices> keeps track of which vertices are "boundary", i.e. those that lie on the curve but are not curve
 * endpoints. This information is used to iterate over boundary wedges/corners during the homology correction step
 * later.
 *
 * <curveEndpoints> returns a vector of which vertices are at curve endpoints.
 *
 * NOTE: This function always orders wedge pairs such that the jump from w_i to w_j is positive.
 *
 * Warning: This function only works if the curve goes through manifold elements.
 */
CornerData<size_t> SurfaceWindingNumbersSolver::indexWedges(const std::vector<Halfedge>& curve, size_t& iN, size_t& N,
                                                            std::vector<Eigen::Triplet<double>>& wedgePairs,
                                                            std::vector<Vertex>& bVertices,
                                                            std::vector<Vertex>& curveEndpoints, bool useSpecialBases) {
    Vector<double> chain = convertToChain(curve);
    geom.requireDECOperators();
    SparseMatrix<double> B1 = geom.d0.transpose();
    Vector<double> chainBoundary = B1 * chain;
    geom.unrequireDECOperators();

    if (psMesh) psMesh->addVertexScalarQuantity("partialGamma", chainBoundary);

    bVertices.clear();
    N = 0;

    // In the end, only corners around curve endpoints should have invalid (i.e. > mesh.nCorners()) wedge indices.
    CornerData<size_t> wedgeIndices(mesh, mesh.nCorners() + 100);

    // Determine interior/boundary vertices (where "boundary" => boundary on the "cut" mesh).
    // Also store an arbitrary outgoing cut halfedge per boundary vertex (for use while indexing boundary wedges, below)
    std::map<Vertex, Halfedge> outgoingCutHalfedge;
    geom.requireEdgeIndices();
    for (Vertex v : mesh.vertices()) {

        size_t nIncidentCuts = 0;
        for (Halfedge he : v.outgoingHalfedges()) {
            if (chain[geom.edgeIndices[he.edge()]] != 0) {
                nIncidentCuts += 1;
                outgoingCutHalfedge.insert(std::make_pair(v, he));
            }
        }

        // Do nothing if curve runs through a nonmanifold vertex (basically just omit that part of the curve.)
        // This gets rid of parts of the curve that are on nonmanifold edges.
        if (abs(chainBoundary[v.getIndex()]) > 1e-5 && v.isManifold() && !v.isBoundary() && useSpecialBases) {
            curveEndpoints.push_back(v);
        } else if (v.isManifold() &&
                   ((nIncidentCuts > 1 && !v.isBoundary()) || (nIncidentCuts > 0 && v.isBoundary()))) {
            bVertices.push_back(v); // vertex is on curve
        } else {
            for (Corner c : v.adjacentCorners()) wedgeIndices[c] = N;
            N++;
        }
    }
    iN = N;

    // Index boundary wedges.
    wedgePairs.clear();
    for (const Vertex& v : bVertices) {

        // 1. Index all wedges on this vertex first.
        Halfedge start = outgoingCutHalfedge[v]; // get an initial cut edge on this vertex

        // Go CW around the vertex starting from <start>.
        Halfedge curr = start;
        do {
            do {
                curr = curr.twin().next();       // go CW around vertex
                wedgeIndices[curr.corner()] = N; // each corner in this wedge gets the same wedge index
                // until we reach the next cut edge, or hit the boundary
            } while (chain[geom.edgeIndices[curr.edge()]] == 0 && !curr.edge().isBoundary());
            N++;
        } while (curr != start);

        // 2. Then visit all wedges incident on this vertex, and assign wedge pairs & jump values.
        curr = start;
        size_t startWIdx = wedgeIndices[curr.corner()];
        std::vector<Eigen::Triplet<double>> jumpPairs;
        double totalJump = 0.;
        do {
            int jumpVal = chain[geom.edgeIndices[curr.edge()]];
            if (!curr.orientation()) jumpVal = -jumpVal;
            if (jumpVal != 0) {
                size_t w_i = wedgeIndices[curr.corner()];
                size_t w_j = wedgeIndices[curr.twin().next().corner()];
                if (jumpVal > 0) {
                    jumpPairs.emplace_back(w_j, w_i, jumpVal);
                } else {
                    jumpPairs.emplace_back(w_i, w_j, -jumpVal);
                }
                totalJump += jumpVal;
            }
            curr = curr.twin().next();
        } while (curr != start);
        // If the jump values around a vertex are compatible (i.e. the jumps sum to zero), then one of the jump pairs
        // around the vertex is redundant. TODO: For now, also do this even if the jumps are incompatible.
        // jumpPairs.pop_back();

        // // If the jump values around a vertex are incompatible (i.e. the jumps don't sum to zero), then distribute
        // the
        // // "excess" evenly across all jump pairs around this vertex. TODO: This may change
        // else {
        //     size_t nJumpPairs = jumpPairs.size();
        //     for (size_t i = 0; i < nJumpPairs; i++) {
        //         double jumpVal = jumpPairs[i].value();
        //         double newJumpVal = jumpVal - signs[i] * totalJump / nJumpPairs;
        //         jumpPairs[i] = {jumpPairs[i].row(), jumpPairs[i].col(), newJumpVal};
        //     }
        // }
        wedgePairs.insert(wedgePairs.end(), jumpPairs.begin(), jumpPairs.end());
    }
    geom.unrequireEdgeIndices();

    return wedgeIndices;
}

/*
 * Given a 1-form γ, compute the integrated divergence (Neumann BCs) on the cut mesh.
 */
Vector<double> SurfaceWindingNumbersSolver::computeIntegratedDivergenceOverBoundaryCells(
    const Vector<double>& gamma, const std::vector<Vertex>& bVertices, const CornerData<size_t>& wedgeIndices,
    size_t N) const {

    Vector<double> rhs = Vector<double>::Zero(N);
    geom.requireEdgeIndices();
    for (const auto& v : bVertices) {
        for (Corner c : v.adjacentCorners()) {
            // Get the two halfedges bounding this corner, on this one side of the "cut" mesh.
            Halfedge heA = c.halfedge();
            Halfedge heB = heA.next().next();
            size_t eIdxA = geom.edgeIndices[heA.edge()];
            size_t eIdxB = geom.edgeIndices[heB.edge()];
            double wA = halfedgeCotanWeights[heA] * gamma[eIdxA];
            double wB = halfedgeCotanWeights[heB] * gamma[eIdxB];
            size_t wIdx = wedgeIndices[c];
            rhs[wIdx] += (v == heA.edge().secondVertex()) ? wA : -wA; // if flow is "into" v
            rhs[wIdx] += (v == heB.edge().secondVertex()) ? wB : -wB;
        }
    }
    geom.unrequireEdgeIndices();

    return rhs;
}

/*
 * Given a 1-form ω, integrate its divergence over the whole mesh.
 * Note: possibly temp function
 */
Vector<double> SurfaceWindingNumbersSolver::computeIntegratedDivergence(const Vector<double>& omega,
                                                                        const CornerData<size_t>& wedgeIndices,
                                                                        size_t N) const {

    geom.requireEdgeIndices();
    geom.requireHalfedgeCotanWeights();
    // geom.requireEdgeCotanWeights();
    Vector<double> rhs = Vector<double>::Zero(N);
    for (const auto& v : mesh.vertices()) {
        for (Corner c : v.adjacentCorners()) {
            // Get the two halfedges bounding this corner, on this one side of the "cut" mesh.
            Halfedge heA = c.halfedge();
            Halfedge heB = heA.next().next();
            size_t eIdxA = geom.edgeIndices[heA.edge()];
            size_t eIdxB = geom.edgeIndices[heB.edge()];
            double wA = halfedgeCotanWeights[heA] * omega[eIdxA];
            double wB = halfedgeCotanWeights[heB] * omega[eIdxB];
            // double wA = 0.5 * geom.edgeCotanWeights[heA.edge()] * omega[eIdxA];
            // double wB = 0.5 * geom.edgeCotanWeights[heB.edge()] * omega[eIdxB];
            size_t wIdx = wedgeIndices[c];
            rhs[wIdx] += (v == heA.edge().secondVertex()) ? wA : -wA; // if flow is "into" v
            rhs[wIdx] += (v == heB.edge().secondVertex()) ? wB : -wB;
        }
    }
    // geom.unrequireEdgeCotanWeights();
    geom.unrequireHalfedgeCotanWeights();
    geom.unrequireEdgeIndices();
    return rhs;
}

/*
 * Integrate gamma while minimizing the L1 norm of the jumps, via linear program.
 */
Vector<double> SurfaceWindingNumbersSolver::integrateGammaAndMinimizeJumpsL1(const std::vector<Halfedge>& closedCurve,
                                                                             const Vector<double>& gamma,
                                                                             CornerData<size_t>& wedgeIndices_v) {

    size_t iN, N;
    std::vector<Eigen::Triplet<double>> jumpPairs_v;
    std::vector<Vertex> bVertices;
    std::vector<Vertex> curveEndpoints;
    wedgeIndices_v = indexWedges(closedCurve, iN, N, jumpPairs_v, bVertices, curveEndpoints);

    // Get rid of redundant jump pairs.
    std::vector<Eigen::Triplet<double>> newJumpPairs;
    size_t idx = 0;
    while (idx < jumpPairs_v.size()) {
        newJumpPairs.push_back(jumpPairs_v[idx]);
        size_t w_i = jumpPairs_v[idx].row();
        size_t w_j = jumpPairs_v[idx].col();
        if (idx + 1 > jumpPairs_v.size()) break;
        // All duplicate jump pairs should be next to each other.
        size_t w_i_next = jumpPairs_v[idx + 1].row();
        size_t w_j_next = jumpPairs_v[idx + 1].col();
        // All jump pairs should be oriented s.t. jump is positive.
        if (w_i == w_i_next && w_j == w_j_next) {
            idx += 2; // skip the next pair!
        } else {
            idx += 1;
        }
    }

    std::vector<std::vector<Corner>> wedgeToCorners(N);
    for (Corner c : mesh.corners()) {
        size_t wIdx = wedgeIndices_v[c];
        if (isValidWedgeIndex(wIdx)) wedgeToCorners[wIdx].push_back(c);
    }
    // CornerData<double> newWedgePairs(mesh, 0);
    // CornerData<double> wedgePairs(mesh, 0);
    // for (size_t i = 0; i < newJumpPairs.size(); i++) {
    //     size_t w_i = newJumpPairs[i].row();
    //     size_t w_j = newJumpPairs[i].col();
    //     for (Corner c : wedgeToCorners[w_i]) newWedgePairs[c] = -1;
    //     for (Corner c : wedgeToCorners[w_j]) newWedgePairs[c] = 1;
    // }
    // for (size_t i = 0; i < jumpPairs_v.size(); i++) {
    //     size_t w_i = jumpPairs_v[i].row();
    //     size_t w_j = jumpPairs_v[i].col();
    //     for (Corner c : wedgeToCorners[w_i]) wedgePairs[c] = -1;
    //     for (Corner c : wedgeToCorners[w_j]) wedgePairs[c] = 1;
    // }
    // if (psMesh) psMesh->addCornerScalarQuantity("wedge pairs", wedgePairs);
    // if (psMesh) psMesh->addCornerScalarQuantity("new wedge pairs", newWedgePairs);

    size_t m = newJumpPairs.size();
    size_t numVars = N + m; // N DOFs + m slack variables
    std::cerr << "Get an instance of a LinearProblem for integrateGammaAndMinimizeJumpsL1()..." << std::endl;
    COMISO::LinearProblem lp(numVars);

    // Set up objective using jump pairs.
    std::cerr << "Setting up objective..." << std::endl;
    std::vector<double> jumpLengths = computeJumpLengths(jumpPairs_v, newJumpPairs, wedgeIndices_v, bVertices);
    Vector<double> checkLengths = Vector<double>::Zero(N);
    for (size_t i = 0; i < m; i++) {
        size_t w_i = newJumpPairs[i].row();
        size_t w_j = newJumpPairs[i].col();
        if (!isValidWedgeIndex(w_i) || !isValidWedgeIndex(w_j)) throw std::logic_error("bad wedge index");
        double jumpVal = newJumpPairs[i].value();
        double coeff = jumpLengths[i] / abs(jumpVal);
        lp.coeffs()[N + i] = coeff;
        checkLengths[w_i] -= coeff;
        checkLengths[w_j] += coeff;
    }
    CornerData<double> checkWedgeLengths(mesh, cornerVector(checkLengths, wedgeIndices_v));
    if (psMesh) psMesh->addCornerScalarQuantity("check lengths", checkWedgeLengths);

    std::cerr << "Setting up variables..." << std::endl;
    std::vector<COMISO::PairIndexVtype> X(numVars);
    std::vector<COMISO::NConstraintInterface*> constraints;
    // Set up solution vector
    for (size_t j = 0; j < numVars; j++) {
        X[j] = COMISO::PairIndexVtype(j, COMISO::Real);
    }
    std::cerr << "Setting up constraints..." << std::endl;
    // Add constraint to enforce that the slack variables are positive.
    for (size_t i = 0; i < m; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(N + i) = 1.;
        COMISO::LinearConstraint* lc =
            new COMISO::LinearConstraint(coeffs, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lc);
    }
    // Add constraints to define the slack variables.
    for (size_t i = 0; i < m; i++) {
        size_t w_i = newJumpPairs[i].row();
        size_t w_j = newJumpPairs[i].col();

        COMISO::LinearConstraint::SVectorNC coeffsA(numVars);
        coeffsA.coeffRef(w_i) = -1.;
        coeffsA.coeffRef(w_j) = 1.;
        coeffsA.coeffRef(N + i) = -1.;
        COMISO::LinearConstraint* lcA =
            new COMISO::LinearConstraint(coeffsA, 0., COMISO::LinearConstraint::NC_LESS_EQUAL);
        constraints.push_back(lcA);

        COMISO::LinearConstraint::SVectorNC coeffsB(numVars);
        coeffsB.coeffRef(w_i) = -1.;
        coeffsB.coeffRef(w_j) = 1.;
        coeffsB.coeffRef(N + i) = 1.;
        COMISO::LinearConstraint* lcB =
            new COMISO::LinearConstraint(coeffsB, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcB);
    }
    // Add constraint that each jump cannot exceed what it was originally in u.
    // Don't need to match up jumps, since this function is only used the "closed curve" algorithm, where the wedge
    // pairs in the original/completed curve are the same.
    for (size_t i = 0; i < m; i++) {
        double jumpVal = newJumpPairs[i].value();
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(N + i) = 1.;
        COMISO::LinearConstraint* lc =
            new COMISO::LinearConstraint(coeffs, -abs(jumpVal), COMISO::LinearConstraint::NC_LESS_EQUAL);
        constraints.push_back(lc);
    }

    // Add constraint to enforce integration. Use first-order constraint that dv = h.
    geom.requireEdgeIndices();
    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        double gamma_i = gamma[geom.edgeIndices[e]];
        Corner cA_i = he.corner();
        Corner cA_j = he.next().corner();
        size_t wA_i = wedgeIndices_v[cA_i];
        size_t wA_j = wedgeIndices_v[cA_j];
        if (!isValidWedgeIndex(wA_i) || !isValidWedgeIndex(wA_j)) throw std::logic_error("bad wedge index");
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(wA_j) = 1.;
        coeffs.coeffRef(wA_i) = -1.;
        COMISO::LinearConstraint* lc =
            new COMISO::LinearConstraint(coeffs, -gamma_i, COMISO::LinearConstraint::NC_EQUAL);
        constraints.push_back(lc);

        if (!e.isBoundary()) {
            Corner cB_i = he.twin().next().corner();
            Corner cB_j = he.twin().corner();
            size_t wB_i = wedgeIndices_v[cB_i];
            size_t wB_j = wedgeIndices_v[cB_j];
            if (!isValidWedgeIndex(wB_i) || !isValidWedgeIndex(wB_j)) throw std::logic_error("bad wedge index");
            if (wB_i != wA_i || wB_j != wA_j) {
                COMISO::LinearConstraint::SVectorNC coeffsB(numVars);
                coeffsB.coeffRef(wB_j) = 1.;
                coeffsB.coeffRef(wB_i) = -1.;
                COMISO::LinearConstraint* lcB =
                    new COMISO::LinearConstraint(coeffsB, -gamma_i, COMISO::LinearConstraint::NC_EQUAL);
                constraints.push_back(lcB);
            }
        }
    }
    geom.unrequireEdgeIndices();

#if (COMISO_GUROBI_AVAILABLE)
    std::cout << "Getting GUROBI solver... " << std::endl;
    COMISO::GUROBISolver gsol;

    std::cout << "Solve..." << std::endl;
    double time_limit = 1e10; // time limit in seconds
    gsol.solve(&lp, constraints, X, time_limit);
    std::cerr << "Solved" << std::endl;
#endif

    Vector<double> v_vector(N);
    for (size_t i = 0; i < N; i++) {
        v_vector[i] = lp.x()[i];
    }
    std::cerr << "integrateGammaAndMinimizeJumpsL1() solved" << std::endl;
    std::cerr << "Objective is " << lp.eval_f(lp.x().data()) << std::endl;

    // Delete pointers
    for (size_t i = 0; i < constraints.size(); i++) {
        delete constraints[i];
    }

    // // Check jump values...
    // for (size_t i = 0; i < m; i++) {
    //     size_t w_i = newJumpPairs[i].row();
    //     size_t w_j = newJumpPairs[i].col();
    //     std::cerr << abs(v_vector[w_i] - v_vector[w_j]) << std::endl;
    // }

    return v_vector;
}

/*
 * Try reduced-size LP, by separating application of shortest path heuristic +
 * minimizing L1 norm of the jumps.
 */
CornerData<double> SurfaceWindingNumbersSolver::reducedLinearProgram(const Vector<double>& gamma,
                                                                     const std::vector<Halfedge>& curve) {

    geom.requireEdgeIndices();

    // Complete curve using a shortest-path heuristic.
    // std::vector<Halfedge> completedCurve =
    //     dijkstraCompleteCurve(curve); // disallow curves from connecting to
    //     themselves
    std::vector<Halfedge> completedCurve = dijkstraCancelBoundary(curve);
    EdgeData<bool> dijkCurve(mesh, false);
    for (const Halfedge& he : completedCurve) dijkCurve[he.edge()] = true;
    if (psMesh) psMesh->addEdgeScalarQuantity("Dijkstra-completed curve", dijkCurve);

    // Detect connected components.
    // In the same loop, integrate gamma within each connected component.
    Vector<double> chain = convertToChain(completedCurve);
    double eps = 1e-5;
    int regionLabel = 0;
    FaceData<int> visitedFace(mesh, 0);
    CornerData<double> vPreShift(mesh, 0);
    for (Face seedFace : mesh.faces()) {

        if (visitedFace[seedFace] > 0) continue;

        // Locally integrate within face.
        Halfedge he0 = seedFace.halfedge();
        Halfedge he1 = he0.next();
        Halfedge he2 = he1.next();
        vPreShift[he1.corner()] =
            he0.orientation() ? gamma[geom.edgeIndices[he0.edge()]] : -gamma[geom.edgeIndices[he0.edge()]];
        vPreShift[he2.corner()] = vPreShift[he1.corner()];
        vPreShift[he2.corner()] +=
            he1.orientation() ? gamma[geom.edgeIndices[he1.edge()]] : -gamma[geom.edgeIndices[he1.edge()]];

        // BFS
        regionLabel += 1;
        std::vector<Face> bag = {seedFace};
        while (bag.size() > 0) {
            Face currFace = bag.back();
            bag.pop_back();
            visitedFace[currFace] = regionLabel;
            for (Halfedge he : currFace.adjacentHalfedges()) {
                Face neighbor = he.twin().face();
                if (abs(chain[geom.edgeIndices[he.edge()]]) < eps && visitedFace[neighbor] == 0) {
                    bag.push_back(neighbor);
                    // Locally integrate within face, and match across edge.
                    he0 = he.twin().next();
                    he1 = he0.next();
                    he2 = he1.next();
                    vPreShift[he0.corner()] = vPreShift[he.corner()];
                    vPreShift[he1.corner()] = vPreShift[he0.corner()];
                    vPreShift[he1.corner()] += (he0.orientation() ? gamma[geom.edgeIndices[he0.edge()]]
                                                                  : -gamma[geom.edgeIndices[he0.edge()]]);
                    vPreShift[he2.corner()] = vPreShift[he1.corner()];
                    vPreShift[he2.corner()] += (he1.orientation() ? gamma[geom.edgeIndices[he1.edge()]]
                                                                  : -gamma[geom.edgeIndices[he1.edge()]]);
                }
            }
        }
    }
    geom.unrequireEdgeIndices();
    int nComponents = regionLabel;
    if (psMesh) psMesh->addFaceScalarQuantity("connected components", visitedFace);
    if (psMesh) psMesh->addCornerScalarQuantity("v pre-shift", vPreShift);

    EdgeData<size_t> bEdgeIdx(mesh, 0); // dense re-indexing of edges on component boundaries
    std::vector<Edge> bEdges;
    size_t reIdx = 0;
    for (Edge e : mesh.edges()) {
        if (e.isBoundary()) continue;
        Halfedge he = e.halfedge();
        if (visitedFace[he.face()] != visitedFace[he.twin().face()]) {
            bEdgeIdx[e] = reIdx;
            bEdges.push_back(e);
            reIdx++;
        }
    }

    // Run reduced-size LP, where [# of DOFs] = [# of connected components].
    size_t numVars = nComponents + reIdx; // DOFs + slack variables
    std::cerr << "Get an instance of a LinearProblem for reducedLinearProgram()..." << std::endl;
    COMISO::LinearProblem lp(numVars);

    std::cerr << "Setting up objective..." << std::endl;
    geom.requireEdgeLengths();
    for (size_t i = 0; i < reIdx; i++) {
        // |total shift| = |existing diff. across bEdge i + component shift}
        lp.coeffs()[nComponents + i] = geom.edgeLengths[bEdges[i]];
    }
    geom.unrequireEdgeLengths();

    std::cerr << "Setting up variables..." << std::endl;
    std::vector<COMISO::PairIndexVtype> X(numVars);
    for (size_t j = 0; j < numVars; j++) {
        X[j] = COMISO::PairIndexVtype(j, COMISO::Real);
    }

    std::cerr << "Setting up constraints..." << std::endl;
    std::vector<COMISO::NConstraintInterface*> constraints;

    // Add constraint to enforce that the slack variables are positive.
    for (size_t i = 0; i < reIdx; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(nComponents + i) = 1.;
        COMISO::LinearConstraint* lc =
            new COMISO::LinearConstraint(coeffs, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lc);
    }

    // Add constraints on slack variables so that for each edge that was
    // originally in the curve, the jump in v is at most jump in u.
    Vector<double> origChain = convertToChain(curve);
    for (size_t i = 0; i < reIdx; i++) {
        double chainCoeff = origChain[geom.edgeIndices[bEdges[i]]];
        if (abs(chainCoeff) > eps) {
            COMISO::LinearConstraint::SVectorNC coeffs(numVars);
            coeffs.coeffRef(nComponents + i) = 1.;
            COMISO::LinearConstraint* lc =
                new COMISO::LinearConstraint(coeffs, -abs(chainCoeff), COMISO::LinearConstraint::NC_LESS_EQUAL);
            constraints.push_back(lc);
        }
    }

    // Add constraints to define the slack variables
    for (size_t i = 0; i < reIdx; i++) {
        Edge e = bEdges[i];
        Halfedge he = e.halfedge();
        Corner c0 = he.corner();
        Corner c1 = he.next().corner();
        Corner d0 = he.twin().next().corner();
        Corner d1 = he.twin().next().next().next().corner();
        int r0 = visitedFace[he.face()] - 1;
        int r1 = visitedFace[he.twin().face()] - 1;
        double initShift = 0.5 * ((vPreShift[c0] - vPreShift[d0]) + (vPreShift[c1] - vPreShift[d1]));

        COMISO::LinearConstraint::SVectorNC coeffsA(numVars);
        coeffsA.coeffRef(r0) = -1.;
        coeffsA.coeffRef(r1) = 1.;
        coeffsA.coeffRef(nComponents + i) = -1.;
        COMISO::LinearConstraint* lcA =
            new COMISO::LinearConstraint(coeffsA, -initShift, COMISO::LinearConstraint::NC_LESS_EQUAL);
        constraints.push_back(lcA);

        COMISO::LinearConstraint::SVectorNC coeffsB(numVars);
        coeffsB.coeffRef(r0) = -1.;
        coeffsB.coeffRef(r1) = 1.;
        coeffsB.coeffRef(nComponents + i) = 1.;
        COMISO::LinearConstraint* lcB =
            new COMISO::LinearConstraint(coeffsB, -initShift, COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcB);
    }

#if (COMISO_GUROBI_AVAILABLE)
    std::cout << "Getting GUROBI solver... " << std::endl;
    COMISO::GUROBISolver gsol;

    std::cout << "Solve..." << std::endl;
    double time_limit = 1e10; // time limit in seconds
    gsol.solve(&lp, constraints, X, time_limit);
    std::cerr << "Solved" << std::endl;
#endif

    // Reconstruct solution using shifts solved for from LP.
    CornerData<double> v_corners = vPreShift; // v post-shift
    for (size_t i = 0; i < nComponents; i++) {
        double shift = lp.x()[i];
        std::cerr << shift << std::endl;
        for (Face f : mesh.faces()) {
            if (visitedFace[f] != i + 1) continue;
            for (Corner c : f.adjacentCorners()) v_corners[c] += shift;
        }
    }
    std::cerr << "reducedLinearProgram() solved" << std::endl;
    std::cerr << "Objective is " << lp.eval_f(lp.x().data()) << std::endl;

    // Delete pointers
    for (size_t i = 0; i < constraints.size(); i++) {
        delete constraints[i];
    }

    return v_corners;
}

/*
 * Try new algo for "curve completion". Integrate up h exactly, but allow v to jump across any edge. We will try to
 * *encourage* v to jump across the original curve while minimizing jumps elsewhere.
 */
CornerData<double> SurfaceWindingNumbersSolver::integrateExactlyBranchCutAnywhere(const Vector<double>& gamma,
                                                                                  const std::vector<Halfedge>& curve) {

    return reducedLinearProgram(gamma, curve);

    Vector<double> chain = convertToChain(curve);
    CornerData<double> vPreShift = integrateExactlyPerFace(gamma, FaceData<double>(mesh, 0));
    // if (psMesh) psMesh->addEdgeScalarQuantity("CHAIN", chain);

    size_t F = mesh.nFaces();
    size_t E = mesh.nEdges();
    size_t numVars = F + E; // DOFs + slack variables
    std::cerr << "Get an instance of a LinearProblem for integrateExactlyBranchCutAnywhere()..." << std::endl;
    COMISO::LinearProblem lp(numVars);

    // TODO: Try weighting edges inversely proportional to their distance from the input curve

    // Set up objective.
    std::cerr << "Setting up objective..." << std::endl;
    geom.requireEdgeLengths();
    double eps = 1e-2;
    for (size_t i = 0; i < E; i++) {
        double coeff = (abs(chain[i]) < 1e-5) ? 1. : eps;
        lp.coeffs()[F + i] = coeff * geom.edgeLengths[i];
    }
    geom.unrequireEdgeLengths();

    std::cerr << "Setting up variables..." << std::endl;
    std::vector<COMISO::PairIndexVtype> X(numVars);
    std::vector<COMISO::NConstraintInterface*> constraints;
    // Set up solution vector
    for (size_t j = 0; j < numVars; j++) {
        X[j] = COMISO::PairIndexVtype(j, COMISO::Real);
    }
    std::cerr << "Setting up constraints..." << std::endl;
    // Add constraint to enforce that the slack variables are positive.
    for (size_t i = 0; i < E; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(F + i) = 1.;
        COMISO::LinearConstraint* lc =
            new COMISO::LinearConstraint(coeffs, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lc);

        // // Add constraints on slack variables so that for each edge that was originally in the curve, the jump in v
        // is
        // // at most jump in u.
        // if (abs(chain[i]) > 1e-5) {
        //     COMISO::LinearConstraint::SVectorNC coeffsB(numVars);
        //     coeffsB.coeffRef(F + i) = 1.;
        //     COMISO::LinearConstraint* lcB =
        //         new COMISO::LinearConstraint(coeffsB, -abs(chain[i]), COMISO::LinearConstraint::NC_LESS_EQUAL);
        //     constraints.push_back(lcB);
        // }
    }
    // Add constraints to define the slack variables.
    geom.requireFaceIndices();
    for (size_t i = 0; i < E; i++) {
        // Compute the difference in face values across edge i
        Edge ei = mesh.edge(i);
        if (ei.isBoundary()) continue;
        Halfedge he = ei.halfedge();
        // Warning: the below 2 lines assume a manifold mesh
        size_t fA = geom.faceIndices[he.face()];
        size_t fB = geom.faceIndices[he.twin().face()];
        double diff = vPreShift[he.corner()] - vPreShift[he.twin().next().corner()];

        // -z -Cijmu <= 0
        COMISO::LinearConstraint::SVectorNC coeffsA(numVars);
        coeffsA.coeffRef(fA) = -1.;
        coeffsA.coeffRef(fB) = 1.;
        coeffsA.coeffRef(F + i) = -1.;
        COMISO::LinearConstraint* lcA =
            new COMISO::LinearConstraint(coeffsA, -diff, COMISO::LinearConstraint::NC_LESS_EQUAL);
        constraints.push_back(lcA);

        // z - Cijmu >= 0
        COMISO::LinearConstraint::SVectorNC coeffsB(numVars);
        coeffsB.coeffRef(fA) = -1.;
        coeffsB.coeffRef(fB) = 1.;
        coeffsB.coeffRef(F + i) = 1.;
        COMISO::LinearConstraint* lcB =
            new COMISO::LinearConstraint(coeffsB, -diff, COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcB);

        // Gij Cijmu >= 0
        double gij = chain(i);
        COMISO::LinearConstraint::SVectorNC coeffsC(numVars);
        coeffsC.coeffRef(fA) = 1. * gij;
        coeffsC.coeffRef(fB) = -1. * gij;
        COMISO::LinearConstraint* lcC =
            new COMISO::LinearConstraint(coeffsC, diff * gij, COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcC);

        // Gij Cijmu <= Gij * Gij
        COMISO::LinearConstraint::SVectorNC coeffsD(numVars);
        coeffsD.coeffRef(fA) = 1. * gij;
        coeffsD.coeffRef(fB) = -1. * gij;
        COMISO::LinearConstraint* lcD =
            new COMISO::LinearConstraint(coeffsD, diff * gij - gij * gij, COMISO::LinearConstraint::NC_LESS_EQUAL);
        constraints.push_back(lcD);
    }
    geom.unrequireFaceIndices();

#if (COMISO_GUROBI_AVAILABLE)
    std::cout << "Getting GUROBI solver... " << std::endl;
    COMISO::GUROBISolver gsol;

    std::cout << "Solve..." << std::endl;
    double time_limit = 1e10; // time limit in seconds
    gsol.solve(&lp, constraints, X, time_limit);
    std::cerr << "Solved" << std::endl;
#endif

    // Reconstruct solution using shifts solved for from LP.
    CornerData<double> v_corners = vPreShift; // v post-shift
    for (size_t i = 0; i < F; i++) {
        Face f = mesh.face(i);
        double shift = lp.x()[i];
        for (Corner c : f.adjacentCorners()) v_corners[c] += shift;
    }
    std::cerr << "integrateExactlyBranchCutAnywhere() solved" << std::endl;
    std::cerr << "Objective is " << lp.eval_f(lp.x().data()) << std::endl;

    // Delete pointers
    for (size_t i = 0; i < constraints.size(); i++) {
        delete constraints[i];
    }

    return v_corners;
}

CornerData<double> SurfaceWindingNumbersSolver::integrateExactlyPerFace(const Vector<double>& gamma,
                                                                        const FaceData<double>& shiftPerFace) {
    CornerData<double> v_corners(mesh, 0);
    geom.requireEdgeIndices();
    for (Face f : mesh.faces()) {
        double shift = shiftPerFace[f];
        for (Corner c : f.adjacentCorners()) v_corners[c] += shift;
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t eIdx = geom.edgeIndices[he.edge()];
            Corner cA = he.corner();
            Corner cB = he.next().corner();
            double diff = he.orientation() ? gamma[eIdx] : -gamma[eIdx];
            v_corners[cB] = v_corners[cA] + diff;
        }
    }
    geom.unrequireEdgeIndices();
    return v_corners;
}

/*
 * <bVertices> = interior vertices of the curve that are also manifold
 */
std::vector<Eigen::Triplet<double>> SurfaceWindingNumbersSolver::computeCumulativeJumps(
    const std::vector<Halfedge>& curve, const std::vector<Vertex>& bVertices, const CornerData<size_t>& wedgeIndices) {


    std::vector<Eigen::Triplet<double>> jumpPairs;
    Vector<double> chain = convertToChain(curve);
    double eps = 1e-5;
    geom.requireEdgeIndices();

    for (const Vertex& v : bVertices) {

        // start at wedge of v.halfedge() [important for boundary vertices, since this is the first wedge]
        size_t initWedge = wedgeIndices[v.halfedge().corner()];
        double totalJump = 0;
        Halfedge start = v.halfedge().next().next().twin();

        // TODO: bugs happen if boundary edges are ever part of the curve --- filter those out

        // Go CW around the vertex starting from <start>.
        Halfedge curr = start;
        do {
            size_t eIdx = geom.edgeIndices[curr.edge()];
            if (abs(chain[eIdx]) > eps) {
                double jump = (curr.orientation()) ? chain[eIdx] : -chain[eIdx];
                size_t w_in = wedgeIndices[curr.corner()];
                size_t w_out = wedgeIndices[curr.twin().next().corner()];
                totalJump += jump;
                jumpPairs.emplace_back(initWedge, w_in, totalJump);
            }
            curr = curr.next().next().twin();
        } while (curr != start && !curr.edge().isBoundary());
    }
    geom.unrequireEdgeIndices();


    if (psMesh) {
        CornerData<double> cumulativeJumps(mesh, 0);
        Vector<double> wedgeCumulativeJump(mesh.nCorners(), 5000);
        for (Eigen::Triplet<double> jp : jumpPairs) {
            wedgeCumulativeJump(jp.col()) = jp.value();
        }
        for (Corner c : mesh.corners()) {
            if (isValidWedgeIndex(wedgeIndices[c])) {
                cumulativeJumps[c] = wedgeCumulativeJump(wedgeIndices[c]);
            }
        }
        psMesh->addCornerScalarQuantity("cumulative jump", cumulativeJumps);
        // psMesh->addEdgeScalarQuantity("chain", chain);
    }

    return jumpPairs;
}

/*
 * Solve the jump equation with the given jumps encoded in <jumpPairs>.
 * Use COV trick to get linear system with fewer variables than if we solved saddle point system.
 */
CornerData<double> SurfaceWindingNumbersSolver::solveJumpEquation(const size_t& N,
                                                                  const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                                                  const std::vector<Eigen::Triplet<double>>& qJumpPairs,
                                                                  const CornerData<size_t>& wedgeIndices,
                                                                  const std::vector<Vertex>& curveEndpoints) {
    // Vector<double> jumpVals;
    // SparseMatrix<double> L = buildWedgeLaplace(N, wedgeIndices);
    // shiftDiagonal(L, 1e-8);
    // SparseMatrix<double> B = buildJumpMatrix(N, jumpPairs, jumpVals);
    // size_t m = B.rows();
    // SparseMatrix<double> Z(m, m);
    // SparseMatrix<double> LHS1 = horizontalStack<double>({L, B.transpose()});
    // SparseMatrix<double> LHS2 = horizontalStack<double>({B, Z});
    // SparseMatrix<double> LHS = verticalStack<double>({LHS1, LHS2});
    // Vector<double> RHS = Vector<double>::Zero(N + m);
    // RHS.tail(m) = jumpVals;
    // Vector<double> result = solveSquare(LHS, RHS);
    // Vector<double> u_vector = result.head(N);
    // std::cerr << "Jump eqn. residual: " << (LHS * result - RHS).squaredNorm() << std::endl;
    // // shiftToGetIntegers(u_vector, jumpPairs);
    // CornerData<double> u(mesh, cornerVector(u_vector, wedgeIndices));
    // return u;

    // Omit curve endpoints from the system. First map existing vertices to their new DOF indices.
    VertexData<size_t> DOFindex(mesh);
    VertexData<bool> isCurveEndpoint(mesh, false);
    for (const Vertex& v : curveEndpoints) {
        isCurveEndpoint[v] = true;
    }
    size_t nDOF = 0;
    for (Vertex v : mesh.vertices()) {
        if (!isCurveEndpoint[v]) {
            DOFindex[v] = nDOF;
            nDOF++;
        }
    }

    // Build Laplace matrix.
    FaceData<bool> skippedFace(mesh, false);
    HalfedgeData<bool> skipHalfedge(mesh, false);
    SparseMatrix<double> L(nDOF, nDOF);
    std::vector<Eigen::Triplet<double>> tripletList;
    geom.requireHalfedgeCotanWeights();
    for (Face f : mesh.faces()) {

        // bool skipFace = false;
        // for (Vertex v : f.adjacentVertices()) {
        //     if (isCurveEndpoint[v]) {
        //         skipFace = true;
        //         skippedFace[f] = true;
        //         for (Halfedge he : f.adjacentHalfedges()) skipHalfedge[he] = true;
        //         break;
        //     }
        // }
        // if (skipFace) continue;

        for (Halfedge he : f.adjacentHalfedges()) {
            if (isCurveEndpoint[he.tailVertex()] || isCurveEndpoint[he.tipVertex()]) {
                skipHalfedge[he] = true;
                continue;
            }

            size_t i = DOFindex[he.tailVertex()];
            size_t j = DOFindex[he.tipVertex()];
            double w = halfedgeCotanWeights[he];
            if (std::isnan(w) || !std::isfinite(w)) w = 1;
            tripletList.emplace_back(i, i, w);
            tripletList.emplace_back(j, j, w);
            tripletList.emplace_back(i, j, -w);
            tripletList.emplace_back(j, i, -w);
        }
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
    checkSymmetric(L);
    checkHermitian(L);
    shiftDiagonal(L, 1e-8);
    std::cerr << "built cotan Laplace" << std::endl;

    // Use the boxed equation in Section 2.3 of TACD06. We have to generalize their formula to an arbitrary number of
    // wedges incident on each vertex, by keeping track of the jumps relative to only a single wedge.

    // Keep track of what vertices each wedge is a part of.
    std::vector<Vertex> wedgeToVertex(N);
    for (Vertex v : mesh.vertices()) {
        for (Corner c : v.adjacentCorners()) {
            size_t wIdx = wedgeIndices[c];
            if (isValidWedgeIndex(wIdx)) wedgeToVertex[wIdx] = v;
        }
    }

    // Form "cumulative" jump pairs. All wedge jump pairs around a
    // vertex should be adjacent to each other. Note: This assumes that the input jumpPairs include jumps across all
    // halfedges of the curve, even redundant jumps.
    // std::vector<Eigen::Triplet<double>> qJumpPairs;
    // size_t idx = 0;
    // while (idx < jumpPairs.size()) {
    //     size_t wi = jumpPairs[idx].row();
    //     size_t wj = jumpPairs[idx].col();
    //     double totalJump = jumpPairs[idx].value();
    //     double lastJump = totalJump; // last jump from wi to wj
    //     qJumpPairs.emplace_back(wi, wj, totalJump);
    //     Vertex currVertex = wedgeToVertex[wi];
    //     assert(wedgeToVertex[wi] == wedgeToVertex[wj]);
    //     size_t offset = 1;
    //     size_t wi_last = wi;
    //     size_t wj_last = wj;
    //     while (idx + offset < jumpPairs.size()) {
    //         size_t wi_next = jumpPairs[idx + offset].row();
    //         size_t wj_next = jumpPairs[idx + offset].col();
    //         size_t jumpVal = jumpPairs[idx + offset].value();
    //         Vertex nextVertex = wedgeToVertex[wi_next];
    //         assert(wedgeToVertex[wi_next] == wedgeToVertex[wj_next]);
    //         if (nextVertex != currVertex) {
    //             currVertex = nextVertex;
    //             break;
    //         }
    //         // When we initially indexed the wedges, we should have gone in order around a vertex, s.t. wedge pairs
    //         are
    //         // consecutive.
    //         assert(wi_next == wi_last || wi_next == wj_last || wj_next == wi_last || wj_next == wj_last);
    //         // If there are only 2 wedges in this vertex, no need to add the other redundant wedge.
    //         if (wi_next == wi_last && wj_next == wj_last) {
    //             offset++;
    //             break;
    //         }
    //         // Else, we should form the rest of the "cumulative jumps".
    //         size_t jumpDst;
    //         if (wi_next == wj_last) {
    //             totalJump += jumpVal;
    //             jumpDst = wj_next;
    //         } else if (wi_next == wi_last) {
    //             totalJump = totalJump - lastJump + jumpVal;
    //             jumpDst = wj_next;
    //         } else if (wj_next == wi_last) {
    //             totalJump = totalJump - lastJump - jumpVal;
    //             jumpDst = wi_next;
    //         } else if (wj_next == wj_last) {
    //             totalJump -= jumpVal;
    //             jumpDst = wi_next;
    //         } else {
    //             throw std::logic_error("Something went wrong when computing cumulative jumps.");
    //         }
    //         // Break if we've already circled around to our original jump pair.
    //         if (jumpDst == wi || jumpDst == wj) {
    //             offset++;
    //             break;
    //         } else {
    //             qJumpPairs.emplace_back(wi, jumpDst, totalJump);
    //         }

    //         wi_last = wi_next;
    //         wj_last = wj_next;
    //         lastJump = jumpVal;
    //         offset++;
    //     }
    //     idx += offset;
    // }

    // Build RHS (with DOFs at curve endpoints omitted)
    Vector<double> RHS = Vector<double>::Zero(nDOF);
    size_t m = qJumpPairs.size();
    for (size_t i = 0; i < m; i++) {
        size_t w_i = qJumpPairs[i].row();
        size_t w_j = qJumpPairs[i].col();
        double jumpVal = qJumpPairs[i].value();
        size_t wIdx = w_j;
        Vertex v = wedgeToVertex[wIdx];
        size_t vIdx = DOFindex[v];
        for (Corner c : v.adjacentCorners()) {
            if (wedgeIndices[c] == wIdx) {
                Halfedge heA = c.halfedge();
                Halfedge heB = heA.next().next();

                if (!skipHalfedge[heA]) {
                    RHS[vIdx] -= halfedgeCotanWeights[heA] * jumpVal;
                    RHS[DOFindex[heA.tipVertex()]] += halfedgeCotanWeights[heA] * jumpVal;
                }
                if (!skipHalfedge[heB]) {
                    RHS[vIdx] -= halfedgeCotanWeights[heB] * jumpVal;
                    RHS[DOFindex[heB.tailVertex()]] += halfedgeCotanWeights[heB] * jumpVal;
                }
            }
        }
    }
    geom.unrequireHalfedgeCotanWeights();
    // std::cerr << "rhs built" << std::endl;

    // std::cerr << "Solving cotan Laplace..." << std::endl;
    Vector<double> u_vector = solvePositiveDefinite(L, RHS);

    // Go back and assign the right values to the omitted wedge DOFs.
    // std::cerr << "Assigning right values to corners..." << std::endl;
    CornerData<double> u_corners(mesh);
    for (Vertex v : mesh.vertices()) {
        size_t vIdx = DOFindex[v];
        for (Corner c : v.adjacentCorners()) {
            u_corners[c] = u_vector[vIdx];
        }
    }
    // if (psMesh) psMesh->addCornerScalarQuantity("u test", u_corners);
    for (size_t i = 0; i < m; i++) {
        size_t w_i = qJumpPairs[i].row();
        size_t w_j = qJumpPairs[i].col();
        double jumpVal = qJumpPairs[i].value();
        size_t wIdx = w_j;
        Vertex v = wedgeToVertex[wIdx];
        size_t vIdx = DOFindex[v];
        for (Corner c : v.adjacentCorners()) {
            if (wedgeIndices[c] == wIdx) {
                u_corners[c] = u_vector[vIdx] + jumpVal;
            }
        }
    }
    // if (psMesh) psMesh->addCornerScalarQuantity("u (corrected)", u_corners);

    Vector<double> u_wedges = shiftToGetIntegers(N, u_corners, jumpPairs, wedgeIndices);
    CornerData<double> u(mesh, cornerVector(u_wedges, wedgeIndices));

    // CornerData<double> u(mesh, cornerVector(shiftVector(u_corners.toVector()), wedgeIndices));

    return u;
}

Vector<double> SurfaceWindingNumbersSolver::computeOmega(const CornerData<double>& u,
                                                         const std::vector<Vertex>& curveEndpoints) {

    // if (psMesh) {
    //     std::vector<std::pair<size_t, int>> vCount;
    //     for (const Vertex& v : curveEndpoints) {
    //         vCount.emplace_back(v.getIndex(), 1);
    //     }
    //     psMesh->addVertexCountQuantity("curve endpoints", vCount);
    // }

    geom.requireEdgeIndices();
    Vector<double> omega = Vector<double>::Zero(mesh.nEdges()); // omega := Du
    VertexData<bool> isCurveEndpoint(mesh, false);
    for (const Vertex& v : curveEndpoints) isCurveEndpoint[v] = true;
    for (Edge e : mesh.edges()) {
        // For every edge connecting to a vertex at a curve endpoint, set the value to zero.
        if (isCurveEndpoint[e.firstVertex()] || isCurveEndpoint[e.secondVertex()]) continue;
        size_t i = geom.edgeIndices[e];
        Halfedge he = e.halfedge();
        Corner c1 = he.corner();
        Corner c2 = he.next().corner();
        Corner d2 = c2;
        Corner d1 = c1;
        if (!e.isBoundary()) {
            d2 = he.twin().corner();
            d1 = he.twin().next().corner();
        }
        omega[i] = 0.5 * ((u[c2] - u[c1]) + (u[d2] - u[d1]));
    }
    geom.unrequireEdgeIndices();
    return omega;
}

std::vector<std::vector<Halfedge>> SurfaceWindingNumbersSolver::computeHomologyBasis() {

    std::vector<std::vector<Halfedge>> generators;
    if (mesh.hasBoundary()) {
        generators = buildPrimalRelativeGenerators();
    } else {
        TreeCotree treeCotree(mesh, geom);
        treeCotree.buildPrimalGenerators();
        generators = treeCotree.primalGenerators;
    }

    return generators;
}

/*
 * Compute a homology basis (where generators are given as a set of primal edges), and a cohomology basis (where 1-forms
 * are given as values on primal edges.)
 */
std::tuple<std::vector<std::vector<Halfedge>>, std::vector<Vector<double>>>
SurfaceWindingNumbersSolver::computeHomologyCohomologyBasis() {

    // Compute the homology generators on *primal edges* of the mesh.
    // (More specifically, tree-cotree computes a homotopy basis. But every homotopy basis is also a homology basis.)
    std::vector<std::vector<Halfedge>> generators = computeHomologyBasis();

    // For each generator, solve for the corresponding harmonic 1-form by solving for a harmonic (scalar) potential that
    // has a unit jump across the curve. (Such a potential cannot be defined as being constant 1 on one side and
    // 0 on the other; we have to make sure the potential is indeed harmonic everywhere, up to a unit-jump across the
    // curve.)
    size_t nG = generators.size();
    size_t E = mesh.nEdges();
    std::vector<Vector<double>> oneForms(nG);
    geom.requireEdgeIndices();
    geom.requireHalfedgeCotanWeights();
    for (size_t i = 0; i < nG; i++) {
        size_t iN, N;
        std::vector<Eigen::Triplet<double>> jumpPairs;
        std::vector<Vertex> bVertices;
        std::vector<Vertex> curveEndpoints; // should be none
        CornerData<size_t> wedgeIndices = indexWedges(generators[i], iN, N, jumpPairs, bVertices, curveEndpoints, true);
        assert(curveEndpoints.size() == 0);
        std::vector<Eigen::Triplet<double>> qJumpPairs = computeCumulativeJumps(generators[i], bVertices, wedgeIndices);
        CornerData<double> u = solveJumpEquation(N, jumpPairs, qJumpPairs, wedgeIndices, curveEndpoints);
        Vector<double> omega(E);
        for (Edge e : mesh.edges()) {
            size_t eIdx = geom.edgeIndices[e];
            Halfedge he = e.halfedge();
            Corner cA = he.corner();
            Corner cB = he.next().corner();
            omega[eIdx] = u[cB] - u[cA];
        }
        oneForms[i] = omega;
    }
    geom.unrequireHalfedgeCotanWeights();

    // for (size_t i = 0; i < nG; i++) {
    //     if (psMesh) {
    //         psMesh->addOneFormIntrinsicVectorQuantity("orig omega " + std::to_string(i), tempOneForms[i],
    //                                                   polyscopeEdgeOrientations(mesh));
    //     }
    // }
    geom.unrequireEdgeIndices();

    // Visualize generators & corresponding 1-forms
    for (size_t i = 0; i < nG; i++) {
        EdgeData<bool> eta(mesh, false);
        for (const Halfedge& he : generators[i]) eta[he.edge().getIndex()] = true;
        if (psMesh) {
            psMesh->addEdgeScalarQuantity("eta " + std::to_string(i), eta);
            psMesh->addOneFormIntrinsicVectorQuantity("omega " + std::to_string(i), oneForms[i],
                                                      polyscopeEdgeOrientations(mesh));
        }
    }
    return std::tuple<std::vector<std::vector<Halfedge>>, std::vector<Vector<double>>>{generators, oneForms};
}

/*
 * Compute a relative homology basis (where generators are given as a set of primal edges), and a cohomology basis
 * (where 1-forms are given as values on primal edges.)
 */
std::tuple<std::vector<std::vector<Halfedge>>, std::vector<Vector<double>>>
SurfaceWindingNumbersSolver::computeRelativeHomologyCohomologyBasis() {

    std::vector<std::vector<Halfedge>> generators = buildPrimalRelativeGenerators();
    size_t nG = generators.size();
    std::vector<Vector<double>> oneForms(nG);

    // *Warning*: We are assuming that the elements in the homology basis are generated by connecting the first boundary
    // loop to every other boundary loop. We solve for the corresponding cohomology basis by solving linear systems: a
    // Dirichlet problem where we assign 1 along all boundary loops, and 0 along one.

    // Set interior vertices.
    size_t V = mesh.nVertices();
    Vector<bool> setAMembership(V);
    geom.requireVertexIndices();
    geom.requireCotanLaplacian();
    for (Vertex v : mesh.vertices()) setAMembership[geom.vertexIndices[v]] = !v.isBoundary();
    SparseMatrix<double> C = geom.cotanLaplacian;
    shiftDiagonal(C, 1e-8);

    // Construct the decomposition
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(C, setAMembership, true);

    // Determine bcVals
    Vector<double> iVals, bcVals;
    Vector<double> uVals = Vector<double>::Zero(V);
    Vector<double> rhsValsA = Vector<double>::Zero(mesh.nInteriorVertices());
    BoundaryLoop bl_source;
    // don't know an easier way to get the first loop
    for (BoundaryLoop bl : mesh.boundaryLoops()) {
        bl_source = bl;
        break;
    }
    for (Vertex v : bl_source.adjacentVertices()) uVals[geom.vertexIndices[v]] = 0;
    decomposeVector(decomp, uVals, iVals, bcVals);

    // Solve problem
    Vector<double> combinedRHS = rhsValsA - decomp.AB * bcVals;
    Vector<double> Aresult = solvePositiveDefinite(decomp.AA, combinedRHS);

    // Combine the two boundary conditions and interior solution to a full vector
    Vector<double> omega_vector = reassembleVector(decomp, Aresult, bcVals);
    if (psMesh) psMesh->addOneFormIntrinsicVectorQuantity("omega", omega_vector, polyscopeEdgeOrientations(mesh));

    geom.unrequireCotanLaplacian();
    geom.unrequireVertexIndices();

    // Visualize generators & corresponding 1-forms.
    for (size_t i = 0; i < nG; i++) {
        EdgeData<bool> eta(mesh, false);
        for (const Halfedge& he : generators[i]) eta[he.edge().getIndex()] = true;
        if (psMesh) {
            psMesh->addEdgeScalarQuantity("eta " + std::to_string(i), eta);
        }
    }
    return std::tuple<std::vector<std::vector<Halfedge>>, std::vector<Vector<double>>>{generators, oneForms};
}


/*
 * TODO: Get rid of "junk loops" from the input curve, which is represented as a 1-chain. "Junk loops" are components
 * that do not contain any edge of the original curve.
 */


/*
 * Compute the Darboux derivative omega := du, and perform Hodge decomposition.
 */
void SurfaceWindingNumbersSolver::computeOmegaAndHodgeDecomposition(const CornerData<double>& u,
                                                                    const std::vector<Vertex>& curveEndpoints,
                                                                    Vector<double>& omega, Vector<double>& dAlpha,
                                                                    Vector<double>& deltaBeta, Vector<double>& gamma) {

    VertexData<bool> isCurveEndpoint(mesh, false);
    for (const Vertex& v : curveEndpoints) isCurveEndpoint[v] = true;

    geom.requireEdgeIndices();
    omega = Vector<double>::Zero(mesh.nEdges());
    for (Edge e : mesh.edges()) {
        // For every edge connecting to a vertex at a curve endpoint, set the value to zero.
        if (isCurveEndpoint[e.firstVertex()] || isCurveEndpoint[e.secondVertex()]) continue;

        size_t i = geom.edgeIndices[e];
        Halfedge he = e.halfedge();
        Corner c1 = he.corner();
        Corner c2 = he.next().corner();
        Corner d2 = c2;
        Corner d1 = c1;
        if (!he.edge().isBoundary()) {
            d2 = he.twin().corner();
            d1 = he.twin().next().corner();
        }
        omega[i] = 0.5 * ((u[c2] - u[c1]) + (u[d2] - u[d1]));
    }
    std::cerr << "omega computed" << std::endl;

    // Hodge-decompose.
    ensureHaveHodgeDecompositionSolvers();
    std::cerr << "got hodge solvers" << std::endl;
    dAlpha = computeExactComponent(omega);
    std::cerr << "dAlpha computed" << std::endl;
    deltaBeta = computeCoExactComponent(omega);
    std::cerr << "deltaBeta computed" << std::endl;
    gamma = omega - dAlpha - deltaBeta;
    geom.unrequireEdgeIndices();
}

/*
 * Get final jumps by reading off jumps from v across the original curve.
 */
std::vector<Eigen::Triplet<double>>
SurfaceWindingNumbersSolver::getFinalJumps(const std::vector<Eigen::Triplet<double>>& jumpPairs,
                                           const CornerData<double>& v_corners, const CornerData<size_t>& wedgeIndices,
                                           size_t N) const {

    // Map each wedge in u to one of its constituent corners.
    std::vector<Corner> wedgeToCorner(N);
    for (Corner c : mesh.corners()) {
        if (isValidWedgeIndex(wedgeIndices[c])) wedgeToCorner[wedgeIndices[c]] = c;
    }

    // Read off the jumps from (u - v) across original parts of the curve, to prescribe them as jumps in one last
    // jump equation solve.
    std::cerr << "New jump vals: ";
    std::vector<Eigen::Triplet<double>> newJumpPairs;
    for (const auto& jumpPair : jumpPairs) {
        size_t w_i = jumpPair.row();
        size_t w_j = jumpPair.col();
        double jumpVal = jumpPair.value(); // orig. jump val in u

        Corner c_i = wedgeToCorner[w_i]; // get any corner of the wedge
        Corner c_j = wedgeToCorner[w_j];
        // I am assuming that the wedges used to solve for v are always a superset of the wedges used to solve for u.
        double newJumpVal = jumpVal - (v_corners[c_j] - v_corners[c_i]);
        // if (abs(v_corners[c_j] - v_corners[c_i]) > 0.1)  newJumpVal = 0;
        // std::cerr << newJumpVal << " ";
        newJumpPairs.emplace_back(w_i, w_j, newJumpVal);
    }
    std::cerr << std::endl;
    return newJumpPairs;
}


// ==== POISSON EQUATION

/*
 * Value at each vertex is the average of the edge values.
 */
VertexData<double> SurfaceWindingNumbersSolver::edgeMidpointDataToVertexData(const Vector<double>& u) {

    VertexData<double> w(mesh, 0);
    for (Vertex v : mesh.vertices()) {
        for (Edge e : v.adjacentEdges()) {
            w[v] += u[geom.edgeIndices[e]];
        }
        w[v] /= v.degree();
    }
    return w;
}

/*
 * Build the solver for when we solve the Poisson equation on the input mesh using Crouzeix-Raviart basis elements.
 *
 * The output of this function only remains valid if the mesh remains un-mutated! (which is ensured since the mesh, etc.
 * variables are private; if these change, the user must re-construct a new SurfaceWindingNumbersSolver.)
 */
void SurfaceWindingNumbersSolver::ensureHaveCrouzeixRaviartPoissonSolver() {

    if (crouzeixRaviartPoissonSolver != nullptr) return;

    geom.requireCrouzeixRaviartLaplacian();
    SparseMatrix<double>& L = geom.crouzeixRaviartLaplacian;
    crouzeixRaviartPoissonSolver.reset(new PositiveDefiniteSolver<double>(L));

    geom.unrequireCrouzeixRaviartLaplacian();
}

/*
 * Build the solver for when we solve the Poisson equation on the input mesh using the standard Whitney basis
 * elements (piecewise linear hat functions.)
 *
 * The output of this function only remains valid if the mesh remains un-mutated! (which is ensured since the mesh, etc.
 * variables are private; if these change, the user must re-construct a new SurfaceWindingNumbersSolver.)
 */
void SurfaceWindingNumbersSolver::ensureHavePoissonSolver() {

    if (poissonSolver != nullptr) return;

    geom.requireCotanLaplacian();
    SparseMatrix<double>& L = geom.cotanLaplacian;
    poissonSolver.reset(new PositiveDefiniteSolver<double>(L));

    geom.unrequireCotanLaplacian();
}

/*
 * Convert a curve in <curveNodes>, <curveEdges> form to a series of BarycentricVector.
 */
std::vector<BarycentricVector> SurfaceWindingNumbersSolver::convertCurveToBarycentricVectors(
    const std::vector<SurfacePoint>& curveNodes, const std::vector<std::array<size_t, 2>>& curveEdges) const {
    std::vector<BarycentricVector> curve;
    for (const auto& seg : curveEdges) {
        curve.emplace_back(BarycentricVector(curveNodes[seg[0]], curveNodes[seg[1]]));
    }
    return curve;
}

/*
 * Given a linear segment within a face of the mesh (represented by its endpoints, and the face(s) it lies within), sum
 * its contribution onto each *vertex element* (this is when we're solving the Poisson equation using Whitney elements.)
 *
 * Even though the computations for Crouzeix-Raviart vs. Whitney elements are basically the same, let's have separate
 * functions to avoid evaluating an extra if-statement for each segment.
 */
void SurfaceWindingNumbersSolver::discretizeSegmentInWhitneyBasis(const BarycentricVector& seg,
                                                                  Vector<double>& N_Gamma) {

    std::vector<Face> faces;
    switch (seg.type) {
        case BarycentricVectorType::Face:
            faces.push_back(seg.face);
            break;
        case BarycentricVectorType::Edge:
            for (Face f : seg.edge.adjacentFaces()) {
                faces.push_back(f);
            }
            break;
        default:
            return;
            break;
    }
    if (faces.size() == 0) return; // shouldn't get here except for numerical error reasons maybe
    double split = 1. / faces.size();

    for (const Face& f : faces) {
        double A = geom.faceAreas[f];
        for (Halfedge he : f.adjacentHalfedges()) {
            // get vertex across from this halfedge (assuming this is a triangle mesh)
            Vertex v = he.next().tipVertex();
            size_t vIdx = geom.vertexIndices[v];
            // Evaluate the inner product between this segment and the halfedge vector.
            Vector2 edgeCoords = he.orientation() ? Vector2{-1., 1.} : Vector2{1., -1.};
            BarycentricVector heVec = BarycentricVector(he.edge(), edgeCoords);
            N_Gamma[vIdx] += split * dot(geom, seg, heVec) * 0.5 / A;
        }
    }
}

/*
 * Given a linear segment within a face of the mesh (represented by its endpoints, and the face(s) it lies within), sum
 * its contribution onto each *edge element* (this is when we're solving the Poisson equation using Crouzeix-Ravaiart
 * elements.)
 *
 * Even though the computations for Crouzeix-Raviart vs. Whitney elements are basically the same, let's have separate
 * functions to avoid evaluating an extra if-statement for each segment.
 */
void SurfaceWindingNumbersSolver::discretizeSegmentInCrouzeixRaviartBasis(const BarycentricVector& seg,
                                                                          Vector<double>& N_Gamma) {

    std::vector<Face> faces;
    switch (seg.type) {
        case BarycentricVectorType::Face:
            faces.push_back(seg.face);
            break;
        case BarycentricVectorType::Edge:
            for (Face f : seg.edge.adjacentFaces()) {
                faces.push_back(f);
            }
            break;
        default:
            return;
            break;
    }
    if (faces.size() == 0) return; // shouldn't get here except for numerical error reasons maybe
    double split = 1. / faces.size();

    for (const Face& f : faces) {
        double A = geom.faceAreas[f];
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t eIdx = geom.edgeIndices[he.edge()];
            // Evaluate the inner product between this segment and the halfedge vector.
            Vector2 edgeCoords = he.orientation() ? Vector2{-1., 1.} : Vector2{1., -1.};
            BarycentricVector heVec = BarycentricVector(he.edge(), edgeCoords);
            N_Gamma[eIdx] += split * dot(geom, seg, heVec) / A;
        }
    }
}


/*
 * Discretize the RHS of the Poisson equation.
 */
Vector<double> SurfaceWindingNumbersSolver::discretizeRHS(bool useWhitneyElements) {

    Vector<double> N_Gamma;

    if (useWhitneyElements) {
        N_Gamma = Vector<double>::Zero(mesh.nVertices());
        for (const auto& seg : curveVectors) {
            discretizeSegmentInWhitneyBasis(seg, N_Gamma);
        }
    } else {
        N_Gamma = Vector<double>::Zero(mesh.nEdges());
        for (const auto& seg : curveVectors) {
            discretizeSegmentInCrouzeixRaviartBasis(seg, N_Gamma);
        }
    }
    return N_Gamma;
}

/*
 * Solve SWN on a (possibly non-manifold) surface mesh via a Poisson equation, discretized using the usual piecewise
 * linear hat functions.
 */
VertexData<double> SurfaceWindingNumbersSolver::solvePoissonEquation() {

    ensureHavePoissonSolver();

    bool useWhitneyElements = true;
    Vector<double> N_Gamma = discretizeRHS(useWhitneyElements);

    Vector<double> u_vec = poissonSolver->solve(N_Gamma);
    u_vec = shiftVector(u_vec);

    VertexData<double> u(mesh);
    u.fromVector(u_vec);
    return u;
}


/*
 * Solve SWN on a (possibly non-manifold) surface mesh via a Poisson equation, discretized using Crouzeix-Raviart
 * basis elements.
 *
 * The Crouzeix-Raviart Laplacian, like the cotan Laplacian, by default is built by summing up
 * contributions per-face; this should work well enough for non-manifold meshes (except in the case that the curve goes
 * through a non-manifold edge, in which case the output of SWN is ambiguous.)
 */
VertexData<double> SurfaceWindingNumbersSolver::solveCrouzeixRaviartPoissonEquation() {

    ensureHaveCrouzeixRaviartPoissonSolver();

    bool useWhitneyElements = false;
    Vector<double> N_Gamma = discretizeRHS(useWhitneyElements);

    Vector<double> u_vec = crouzeixRaviartPoissonSolver->solve(N_Gamma);
    u_vec = shiftVector(u_vec);

    VertexData<double> u = edgeMidpointDataToVertexData(u_vec);
    return u;
}

// ==== SOLVE

/*
 * Complete the curve using Dijkstra, but it is less restrictive than the dijkstraCompleteCurve() function, which
 * prevents curve components from connecting to themselves.
 */
std::vector<Halfedge> SurfaceWindingNumbersSolver::dijkstraCancelBoundary(const std::vector<Halfedge>& curve) {

    double eps = 1e-5;
    geom.requireDECOperators();
    SparseMatrix<double> B = geom.d0.transpose();
    geom.unrequireDECOperators();
    Vector<double> chain = Vector<double>::Zero(mesh.nEdges());
    for (const Halfedge& he : curve) {
        size_t i = he.edge().getIndex();
        double val = he.orientation() ? 1 : -1;
        chain[i] = val;
    }
    Vector<double> boundary = B * chain;
    if (boundary.norm() < eps) return curve; // the input curve already has no boundary

    std::vector<std::pair<Vertex, bool>> endpoints;
    geom.requireVertexIndices();
    for (Vertex v : mesh.vertices()) {
        size_t vIdx = geom.vertexIndices[v];
        if (abs(boundary[vIdx]) > eps) {
            endpoints.emplace_back(v, boundary[vIdx] > 0);
        }
    }
    geom.unrequireVertexIndices();

    // Create custom EdgeLengthGeometry
    EdgeData<double> dijkstraWeights(mesh);
    geom.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
        double length = geom.edgeLengths[e];
        Halfedge he = e.halfedge();
        dijkstraWeights[e] = length;
    }
    geom.unrequireEdgeLengths();
    // Remove original curve edges from the graph.
    const double infinity = std::numeric_limits<double>::infinity();
    for (const Halfedge& he : curve) {
        dijkstraWeights[he.edge()] = infinity;
    }
    EdgeLengthGeometry dijkstraMesh(mesh, dijkstraWeights);

    // Only connect endpoints of opposite sign.
    Vertex startVert;
    bool sgn;
    std::vector<Halfedge> completedCurve = curve;
    while (endpoints.size() > 0) {
        std::tie(startVert, sgn) = endpoints.back();
        endpoints.pop_back();
        std::set<Vertex> endVerts;
        for (auto& tup : endpoints) {
            if (tup.second != sgn) endVerts.insert(tup.first);
        }
        std::vector<Halfedge> path = dijkstraPath(dijkstraMesh, startVert, endVerts);
        Vertex endVert = path.back().tipVertex();
        if (!sgn) {
            std::reverse(path.begin(), path.end());
            for (Halfedge& he : path) {
                he = he.twin();
            }
        }
        completedCurve.insert(completedCurve.end(), path.begin(), path.end());
        // Determine which endpoint we ended at, so we can delete it.
        for (size_t i = 0; i < endpoints.size(); i++) {
            if (endpoints[i].first == endVert) {
                endpoints.erase(endpoints.begin() + i);
                break;
            }
        }
    }

    // Make sure that the complete curve has no boundary.
    Vector<double> completedChain = Vector<double>::Zero(mesh.nEdges());
    for (const Halfedge& he : completedCurve) {
        size_t i = he.edge().getIndex();
        double val = he.orientation() ? 1 : -1;
        completedChain[i] = val;
    }
    std::cerr << "boundary of completed chain: " << (B * completedChain).norm() << std::endl;
    assert((B * completedChain).norm() < eps);

    return completedCurve;
}

/*
 * Complete the curve using a version of Dijkstra. Prevents curve components from connecting to themselves.
 */
std::vector<Halfedge> SurfaceWindingNumbersSolver::dijkstraCompleteCurve(const std::vector<Halfedge>& curve) {

    // Determine the curve endpoints by determining all connected curve components. We need the latter information to
    // prevent curves from being completed by a reversed copy of themselves.
    std::vector<std::vector<Halfedge>> curveComponents = getCurveComponents(geom, curve);
    std::vector<std::pair<Vertex, bool>> curveEndpoints; // store index of the vertex and its sign
    for (const auto& cpt : curveComponents) {
        // It is possible that a curve component is closed, i.e. doesn't have endpoints; check.
        Vertex vertexA = cpt.front().tailVertex();
        Vertex vertexB = cpt.back().tipVertex();
        if (vertexA != vertexB) {
            curveEndpoints.emplace_back(vertexA, false);
            curveEndpoints.emplace_back(vertexB, true);
        }
    }
    std::cerr << "# curve components: " << curveComponents.size() << std::endl;
    std::cerr << "# pairs of curve endpoints: " << curveEndpoints.size() << std::endl;
    if (curveEndpoints.size() == 0) return curve;

    // Visualize curve endpoints to make sure we got them right.
    std::vector<std::pair<size_t, int>> vertexCountValues;
    for (auto& tup : curveEndpoints) {
        vertexCountValues.emplace_back(tup.first.getIndex(), !tup.second ? -1 : 1);
    }
    if (psMesh != NULL) psMesh->addVertexCountQuantity("endpoints", vertexCountValues);

    std::vector<Halfedge> completedCurve = curve;
    // Create custom edgelengthgeometry
    EdgeData<double> dijkstraWeights(mesh);
    geom.requireEdgeLengths();
    for (Edge e : mesh.edges()) {
        double length = geom.edgeLengths[e];
        Halfedge he = e.halfedge();
        dijkstraWeights[e] = length;
    }
    geom.unrequireEdgeLengths();
    // Remove original curve edges from the graph.
    const double infinity = std::numeric_limits<double>::infinity();
    for (const auto& cpt : curveComponents) {
        for (const Halfedge& he : cpt) {
            dijkstraWeights[he.edge()] = infinity;
        }
    }
    EdgeLengthGeometry dijkstraMesh(mesh, dijkstraWeights);

    // Only connect endpoints of opposite sign.
    std::vector<std::pair<Vertex, bool>> tempEndpoints = curveEndpoints;
    while (tempEndpoints.size() > 0) {
        Vertex startVert;
        bool sgn;
        std::tie(startVert, sgn) = tempEndpoints.back();
        tempEndpoints.pop_back();
        std::set<Vertex> endVerts;
        assert(tempEndpoints.size() > 0);
        for (auto& tup : tempEndpoints) {
            if (tup.second != sgn) endVerts.insert(tup.first);
        }
        std::vector<Halfedge> path = dijkstraPath(dijkstraMesh, startVert, endVerts);
        assert(path.size() > 0);
        Vertex endVert = path.back().tipVertex();
        if (!sgn) {
            std::reverse(path.begin(), path.end());
            for (Halfedge& he : path) {
                he = he.twin();
            }
        }
        completedCurve.insert(completedCurve.end(), path.begin(), path.end());
        // Determine which endpoint we ended at, so we can delete it.
        for (size_t i = 0; i < tempEndpoints.size(); i++) {
            if (tempEndpoints[i].first == endVert) {
                tempEndpoints.erase(tempEndpoints.begin() + i);
                break;
            }
        }
    }
    return completedCurve;
}

std::vector<Halfedge> SurfaceWindingNumbersSolver::ordinaryOHCP(const std::vector<Halfedge>& curve) {

    std::vector<Halfedge> homologousCurve = curve;

    size_t E = mesh.nEdges();
    size_t F = mesh.nFaces();
    size_t numVars = 2 * E + 2 * F;

    std::cerr << "Get an instance of a LinearProblem for ordinary OHCP..." << std::endl;
    COMISO::LinearProblem lp(numVars);

    // Set up objective.
    geom.requireEdgeLengths();
    for (size_t i = 0; i < E; i++) {
        double length = geom.edgeLengths[i];
        lp.coeffs()[i] = length;
        lp.coeffs()[E + i] = length;
    }
    geom.unrequireEdgeLengths();

    std::cerr << "Setting up variables..." << std::endl;

    std::vector<COMISO::PairIndexVtype> X(numVars);
    std::vector<COMISO::NConstraintInterface*> constraints;
    // x+ and x-, y+ and y-
    for (size_t j = 0; j < numVars; j++) {
        X[j] = COMISO::PairIndexVtype(j, COMISO::Real);
    }
    // Add constraint to enforce that the positive component of each chain coefficient is positive
    for (size_t i = 0; i < numVars; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(i) = 1.;
        COMISO::LinearConstraint* lcA =
            new COMISO::LinearConstraint(coeffs, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcA);
    }


    geom.requireEdgeIndices();
    // Convert input curve into discrete 1-chain.
    Vector<double> c = Vector<double>::Zero(E);
    for (const Halfedge& he : homologousCurve) {
        c[geom.edgeIndices[he.edge()]] += (he.orientation()) ? 1 : -1;
    }
    geom.unrequireEdgeIndices();

    // Set up constraint matrix [I -I -B B]
    std::cerr << "Setting up constraints..." << std::endl;
    geom.requireDECOperators();
    SparseMatrix<double>& BT = geom.d1; // == transpose of boundary matrix in col-major order
    geom.unrequireDECOperators();
    for (size_t i = 0; i < E; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        double rhs = c[i];
        coeffs.coeffRef(i) = 1;
        coeffs.coeffRef(E + i) = -1;
        // Going down columns of B^T is equivalent to going through rows of B.
        for (Eigen::SparseMatrix<double>::InnerIterator it(BT, i); it; ++it) {
            coeffs.coeffRef(2 * E + it.row()) = -it.value();
            coeffs.coeffRef(2 * E + F + it.row()) = it.value();
        }
        COMISO::LinearConstraint* lc = new COMISO::LinearConstraint(coeffs, -rhs, COMISO::LinearConstraint::NC_EQUAL);
        constraints.push_back(lc);
    }

    // Check that the input curve has no boundary.
    SparseMatrix<double> B = BT.transpose();
    std::cerr << "ordinary OHCP input 1-chain boundary norm (want this to be 0): " << (B * c).norm() << std::endl;


#if (COMISO_GUROBI_AVAILABLE)
    std::cout << "Getting GUROBI solver... " << std::endl;
    COMISO::GUROBISolver gsol;

    std::cout << "Solve..." << std::endl;
    double time_limit = 1e10; // time limit in seconds
    gsol.solve(&lp, constraints, X, time_limit);
    std::cerr << "Solved" << std::endl;
#endif

    std::vector<Halfedge> completedCurve;
    Vector<double> x_vector(E);
    for (size_t i = 0; i < E; i++) {
        double chainCoeff = lp.x()[i] - lp.x()[E + i];
        x_vector[i] = chainCoeff;
        if (abs(chainCoeff) > 1e-5) {
            Halfedge he = mesh.edge(i).halfedge();
            if (chainCoeff < 0) he = he.twin();
            for (size_t j = 0; j < (size_t)(abs(chainCoeff)); j++) completedCurve.push_back(he);
        }
    }
    std::cerr << "Ordinary OHCP solved" << std::endl;
    std::cerr << "Objective is " << lp.eval_f(lp.x().data()) << std::endl;

    // // Visualize the 1-chain x
    // if (psMesh) psMesh->addEdgeScalarQuantity("1-chain x (ordinary OHCP)", x_vector);

    // Delete pointers
    for (size_t i = 0; i < constraints.size(); i++) {
        delete constraints[i];
    }

    return completedCurve;
}

/*
 * Find the shortest curve homologous to the input curve, such that it is a subset of the given curve.
 */
std::vector<Halfedge> SurfaceWindingNumbersSolver::OHCP_subset(const std::vector<Halfedge>& curve) {

    size_t E = mesh.nEdges();
    size_t F = mesh.nFaces();
    size_t numVars = 2 * E + 2 * F;

    std::cerr << "Get an instance of a LinearProblem for OHBP..." << std::endl;
    COMISO::LinearProblem lp(numVars);

    // Set up objective.
    geom.requireEdgeLengths();
    for (size_t i = 0; i < E; i++) {
        double length = geom.edgeLengths[i];
        lp.coeffs()[i] = length;
        lp.coeffs()[E + i] = length;
    }
    geom.unrequireEdgeLengths();

    std::cerr << "Setting up variables..." << std::endl;

    std::vector<COMISO::PairIndexVtype> X(numVars);
    std::vector<COMISO::NConstraintInterface*> constraints;
    // x+ and x-, y+ and y-
    for (size_t j = 0; j < numVars; j++) {
        X[j] = COMISO::PairIndexVtype(j, COMISO::Real);
    }
    // Add constraint to enforce that the positive component of each chain coefficient is positive
    for (size_t i = 0; i < numVars; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(i) = 1.;
        COMISO::LinearConstraint* lcA =
            new COMISO::LinearConstraint(coeffs, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcA);
    }


    geom.requireEdgeIndices();
    // Convert input curve into discrete 1-chain.
    Vector<double> c = Vector<double>::Zero(E);
    for (const Halfedge& he : curve) {
        c[geom.edgeIndices[he.edge()]] += (he.orientation()) ? 1 : -1;
    }
    geom.unrequireEdgeIndices();

    // Set up constraint matrix [I -I -B B]
    // Also add constraints to enforce that the output chain is a subset of the input chain.
    std::cerr << "Setting up constraints..." << std::endl;
    geom.requireDECOperators();
    SparseMatrix<double>& BT = geom.d1; // == transpose of boundary matrix in col-major order
    geom.unrequireDECOperators();
    for (size_t i = 0; i < E; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        double rhs = c[i];
        coeffs.coeffRef(i) = 1;
        coeffs.coeffRef(E + i) = -1;
        // Going down columns of B^T is equivalent to going through rows of B.
        for (Eigen::SparseMatrix<double>::InnerIterator it(BT, i); it; ++it) {
            coeffs.coeffRef(2 * E + it.row()) = -it.value();
            coeffs.coeffRef(2 * E + F + it.row()) = it.value();
        }
        COMISO::LinearConstraint* lc = new COMISO::LinearConstraint(coeffs, -rhs, COMISO::LinearConstraint::NC_EQUAL);
        constraints.push_back(lc);

        // Constraint that the output chain is a subset of the input chain.
        if (abs(rhs) < 1e-5) {
            COMISO::LinearConstraint::SVectorNC coeffsA(numVars);
            coeffsA.coeffRef(i) = 1;
            coeffsA.coeffRef(E + i) = -1;
            COMISO::LinearConstraint* lcA =
                new COMISO::LinearConstraint(coeffsA, 0., COMISO::LinearConstraint::NC_EQUAL);
            constraints.push_back(lcA);
        } else if (rhs > 0.) {
            COMISO::LinearConstraint::SVectorNC coeffsA(numVars);
            coeffsA.coeffRef(i) = 1;
            coeffsA.coeffRef(E + i) = -1;
            COMISO::LinearConstraint* lcA =
                new COMISO::LinearConstraint(coeffsA, -rhs, COMISO::LinearConstraint::NC_LESS_EQUAL);
            COMISO::LinearConstraint* lcB =
                new COMISO::LinearConstraint(coeffsA, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
            constraints.push_back(lcA);
            constraints.push_back(lcB);
        } else {
            COMISO::LinearConstraint::SVectorNC coeffsA(numVars);
            coeffsA.coeffRef(i) = 1;
            coeffsA.coeffRef(E + i) = -1;
            COMISO::LinearConstraint* lcA =
                new COMISO::LinearConstraint(coeffsA, -rhs, COMISO::LinearConstraint::NC_GREATER_EQUAL);
            COMISO::LinearConstraint* lcB =
                new COMISO::LinearConstraint(coeffsA, 0., COMISO::LinearConstraint::NC_LESS_EQUAL);
            constraints.push_back(lcA);
            constraints.push_back(lcB);
        }
    }

    // // Check that the input curve has no boundary.
    // SparseMatrix<double> B = BT.transpose();
    // std::cerr << "OHCP_subset input 1-chain boundary norm (want this to be 0): " << (B * c).norm() << std::endl;


#if (COMISO_GUROBI_AVAILABLE)
    std::cout << "Getting GUROBI solver... " << std::endl;
    COMISO::GUROBISolver gsol;

    std::cout << "Solve..." << std::endl;
    double time_limit = 1e10; // time limit in seconds
    gsol.solve(&lp, constraints, X, time_limit);
    std::cerr << "Solved" << std::endl;
#endif

    std::vector<Halfedge> completedCurve;
    Vector<double> x_vector(E);
    for (size_t i = 0; i < E; i++) {
        double chainCoeff = lp.x()[i] - lp.x()[E + i];
        x_vector[i] = chainCoeff;
        if (abs(chainCoeff) > 1e-5) {
            Halfedge he = mesh.edge(i).halfedge();
            if (chainCoeff < 0) he = he.twin();
            for (size_t j = 0; j < (size_t)(abs(chainCoeff)); j++) completedCurve.push_back(he);
        }
    }
    std::cerr << "OHCP_subset solved" << std::endl;
    std::cerr << "Objective is " << lp.eval_f(lp.x().data()) << std::endl;

    // Delete pointers
    for (size_t i = 0; i < constraints.size(); i++) {
        delete constraints[i];
    }

    return completedCurve;
}

/*
 * One possibility we talked about (on 11/10/22) is to integrate gamma exactly via search, resulting in a jump in the
 * "right" homology class. Then run OHCP (while fixing some values of the solution 1-chain vector) to get the shortest
 * homologous chain that contains every edge of the original input curve. This may not handle "splitting", i.e. we will
 * always start with homology coefficients at most one, but we may want to complete the curve s.t. we obtain two
 * generator loops. We could also just form little contractible loops from a dashed line forming a generator loop, while
 * the actual generator loop is elsewhere (like over the skinniest part of the torus)... I sent a Slack to
 * surface-winding-numbers on 11/11/22.
 *
 * Edit: I just realized a nice thing about OHCP; it's guaranteed to return closed loops, because it returns a
 * 1-chain in the form x = c + By, i.e. the closed input 1-chain plus the boundary of a 2-chain, which is always closed.
 * This prevents the possibility of the completed curve only ever completing 1 dashed generator loop, even if there are
 * two. (Weird things happen if the input curve is open -- but then again, all open curves should be null-homologous
 * anyway. We also should not get open curves from the exact integration step.)
 *
 * Edit: I wonder what happens if we minimize the L2 norm instead. This is slightly different than finding the shortest
 * cycle.
 */
std::vector<Halfedge> SurfaceWindingNumbersSolver::completeCurveOHCP(const std::vector<Halfedge>& curve,
                                                                     const Vector<double>& gamma) {

    std::vector<Halfedge> homologousCurve;
    // // First get the desired homology class of curve by integrating gamma via search.
    // CornerData<double> integratedGamma = integrateExactly(gamma, homologousCurve);
    // if (psMesh != NULL) {
    //     EdgeData<bool> edgeOnCurve(mesh, false);
    //     for (const Halfedge& he : homologousCurve) {
    //         edgeOnCurve[he.edge()] = true;
    //     }
    //     psMesh->addEdgeScalarQuantity("branch cut", edgeOnCurve);
    //     psMesh->addCornerScalarQuantity("integrated gamma on faces", integratedGamma);
    // }

    // The branch cut might not be closed (i.e. it might have boundary.) To correct this issue, all we have to do is
    // simply connect each pair of opposite-sign boundary points using Dijkstra, omitting any already-marked edges.
    // homologousCurve = dijkstraCancelBoundary(homologousCurve);

    // Decompose gamma into a cohomology basis. Delete components that are small. Take the ceiling of
    // nonzero coefficients. Use the corresponding curve as the homology class.
    std::vector<std::vector<Halfedge>> homologyBasis;
    std::vector<Vector<double>> cohomologyBasis;
    std::tie(homologyBasis, cohomologyBasis) = computeHomologyCohomologyBasis();
    size_t nG = cohomologyBasis.size();
    Eigen::VectorXd gammaInBasis(nG);
    Vector<double> chain = Vector<double>::Zero(mesh.nEdges());
    std::vector<Halfedge> tempCurve;
    // We can't simply project gamma onto the 1-forms in the cohomology basis, because the basis isn't necessarily
    // orthogonal. Instead, we have to solve a (small) 2g x 2g linear system.
    Eigen::MatrixXd cMat(nG, nG);
    Eigen::VectorXd cRHS(nG);
    for (size_t i = 0; i < nG; i++) {
        for (size_t j = 0; j < nG; j++) {
            cMat(i, j) = cohomologyBasis[i].dot(cohomologyBasis[j]);
        }
        cRHS(i) = cohomologyBasis[i].dot(gamma);
    }
    gammaInBasis = cMat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(cRHS);

    geom.requireEdgeIndices();
    for (size_t i = 0; i < nG; i++) {
        std::cerr << "gammaInBasis[i]: " << gammaInBasis[i] << std::endl;
        if (abs(gammaInBasis[i]) > 1e-3) {
            size_t nCopies = std::ceil(abs(gammaInBasis[i]));
            std::vector<Halfedge> basisCurve = homologyBasis[i];
            if (gammaInBasis[i] < 0) {
                for (Halfedge& he : basisCurve) he = he.twin();
            }
            for (const Halfedge& he : basisCurve) {
                double val = he.orientation() ? 1 : -1;
                val *= nCopies; // because nCopies is a size_t, can't directly use -nCopies above
                chain[geom.edgeIndices[he.edge()]] += val;
            }
            for (size_t j = 0; j < nCopies; j++) {
                tempCurve.insert(tempCurve.end(), basisCurve.begin(), basisCurve.end());
            }
        }
    }
    geom.unrequireEdgeIndices();
    homologousCurve = tempCurve;
    // Visualize 1-forms and generators to debug
    if (psMesh) {
        for (size_t i = 0; i < cohomologyBasis.size(); i++) {
            EdgeData<bool> edgeOnGen(mesh, false);
            for (const Halfedge& he : homologyBasis[i]) edgeOnGen[he.edge()] = true;
            psMesh->addEdgeScalarQuantity("generator " + std::to_string(i), edgeOnGen);
            psMesh->addOneFormIntrinsicVectorQuantity("omega " + std::to_string(i), cohomologyBasis[i],
                                                      polyscopeEdgeOrientations(mesh));
        }
    }

    // Run ordinary OHCP just to make sure everything is behaving as expected.
    std::vector<Halfedge> shortenedCurve = ordinaryOHCP(homologousCurve);
    EdgeData<int> edgeOnShortCurve(mesh, 0);
    for (const Halfedge& he : shortenedCurve) {
        edgeOnShortCurve[he.edge()] += 1;
    }
    if (psMesh) psMesh->addEdgeScalarQuantity("OHCP shortened curve", edgeOnShortCurve);

    return OHCP_supset(chain, curve);
}

/*
 * Run a variant of OHCP. Find the shortest chain homologous to the input chain, that contains the input curve as a
 * *subset*.
 */
std::vector<Halfedge> SurfaceWindingNumbersSolver::OHCP_supset(const Vector<double>& chain,
                                                               const std::vector<Halfedge>& curve) {

    size_t E = mesh.nEdges();
    size_t F = mesh.nFaces();
    size_t numVars = 2 * E + 2 * F;

    std::cerr << "Get an instance of a LinearProblem for OHPP..." << std::endl;
    COMISO::LinearProblem lp(numVars);

    // Set up objective.
    geom.requireEdgeLengths();
    for (size_t i = 0; i < E; i++) {
        double length = geom.edgeLengths[i];
        lp.coeffs()[i] = length;
        lp.coeffs()[E + i] = length;
    }
    geom.unrequireEdgeLengths();

    std::cerr << "Setting up variables..." << std::endl;

    std::vector<COMISO::PairIndexVtype> X(numVars);
    std::vector<COMISO::NConstraintInterface*> constraints;
    // x+ and x-, y+ and y-
    for (size_t j = 0; j < numVars; j++) {
        X[j] = COMISO::PairIndexVtype(j, COMISO::Real);
    }
    // Add constraint to enforce that the positive component of each chain coefficient is positive
    for (size_t i = 0; i < numVars; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        coeffs.coeffRef(i) = 1.;
        COMISO::LinearConstraint* lcA =
            new COMISO::LinearConstraint(coeffs, 0., COMISO::LinearConstraint::NC_GREATER_EQUAL);
        constraints.push_back(lcA);
    }


    geom.requireEdgeIndices();

    // Convert the original curve into a discrete 1-chain.
    Vector<double> c_orig = Vector<double>::Zero(E);
    for (const Halfedge& he : curve) {
        c_orig[geom.edgeIndices[he.edge()]] += (he.orientation()) ? 1 : -1;
    }
    // if (psMesh) psMesh->addEdgeScalarQuantity("input 1-chain", c);
    geom.unrequireEdgeIndices();

    // Set up constraint matrix [I -I -B B]
    std::cerr << "Setting up constraints..." << std::endl;
    geom.requireDECOperators();
    SparseMatrix<double>& BT = geom.d1; // == transpose of boundary matrix in col-major order
    SparseMatrix<double> B1 = geom.d0.transpose();
    geom.unrequireDECOperators();
    for (size_t i = 0; i < E; i++) {
        COMISO::LinearConstraint::SVectorNC coeffs(numVars);
        // Enforcing certain values of the 1-chain should preserve total unimodularity (we're adding zero-columns)
        // Also now there's no constraint on these values of x, so they will be zero (not itself a problem, just have to
        // remember to set these values before returning from the solution.)
        double rhs = chain[i] - c_orig[i];
        // If edge i is not contained in the original curve, enforce the usual constraint. Otherwise the constraint just
        // says that the boundary of the solution 2-chain has the desired i-th coefficient.
        if (abs(c_orig[i]) < 1e-5) {
            coeffs.coeffRef(i) = 1;
            coeffs.coeffRef(E + i) = -1;
        }
        // Going down columns of B^T is equivalent to going through rows of B.
        for (Eigen::SparseMatrix<double>::InnerIterator it(BT, i); it; ++it) {
            coeffs.coeffRef(2 * E + it.row()) = -it.value();
            coeffs.coeffRef(2 * E + F + it.row()) = it.value();
        }
        COMISO::LinearConstraint* lc = new COMISO::LinearConstraint(coeffs, -rhs, COMISO::LinearConstraint::NC_EQUAL);
        constraints.push_back(lc);
    }

#if (COMISO_GUROBI_AVAILABLE)
    std::cout << "Getting GUROBI solver... " << std::endl;
    COMISO::GUROBISolver gsol;

    std::cout << "Solve..." << std::endl;
    double time_limit = 1e10; // time limit in seconds
    gsol.solve(&lp, constraints, X, time_limit);
    std::cerr << "Solved" << std::endl;
#endif

    Vector<double> x_vector = c_orig;
    std::vector<Halfedge> completedCurve = curve;
    for (size_t i = 0; i < E; i++) {
        double chainCoeff = lp.x()[i] - lp.x()[E + i];
        x_vector[i] += chainCoeff;
        if (abs(chainCoeff) > 1e-5) {
            Halfedge he = mesh.edge(i).halfedge();
            if (chainCoeff < 0) he = he.twin();
            completedCurve.push_back(he);
        }
    }
    std::cerr << "OHPP solved" << std::endl;
    std::cerr << "Objective is " << lp.eval_f(lp.x().data()) << std::endl;

    // Delete pointers
    for (size_t i = 0; i < constraints.size(); i++) {
        delete constraints[i];
    }

    // return x_vector;
    return completedCurve;
}

/*
 * Maybe try DFS instead of BFS. During the search, forbid crossing edges of Gamma. During the search (i.e.
 * integration), if we ever arrive at a point where we are going to *not* jump over Gamma, backtrack and try a
 * different neighboring triangle. Keep backtracking and doing the search until we've run out of options for queuing
 * new triangle faces. Of course we could still integrate gamma such that segments that are supposed to be part of a
 * generator loop get completed into small contractible loops. This is less likely to happen the more "filled in"
 * the generator loop is, but is also in theory kind of random (up to the order in which we iterate around faces and
 * do our search...)
 */

/*
 * Generic curve completion code. Swap out the final algorithm as needed.
 */
std::vector<Halfedge> SurfaceWindingNumbersSolver::completeCurve(const std::vector<Halfedge>& curve,
                                                                 const Vector<double>& gamma) {

    // return dijkstraCompleteCurve(curve);
    return completeCurveOHCP(curve, gamma);
}

CornerData<double> SurfaceWindingNumbersSolver::integrateExactly(const Vector<double>& gamma,
                                                                 std::vector<Halfedge>& branchCut) {

    double infinity = std::numeric_limits<double>::infinity();
    VertexData<bool> vertexVisited(mesh, false);
    FaceData<bool> faceVisited(mesh, false);
    CornerData<bool> cornerVisited(mesh, true);
    VertexData<double> vertexVals(mesh, 0.);
    CornerData<double> cornerVals(mesh, 0.);

    geom.requireEdgeIndices();
    double eps = 1e-5;
    // only need this outer loop if the mesh is comprised of multiple disconnected components
    for (Face start : mesh.faces()) {
        if (faceVisited[start]) continue;

        faceVisited[start] = true;
        std::vector<Face> queue = {start};

        // Fill in corner values for the root face. gamma should be locally integrable within each face.
        for (Halfedge he : start.adjacentHalfedges()) {
            double gamma_i = gamma[geom.edgeIndices[he.edge()]];
            if (!he.orientation()) gamma_i *= -1;
            Corner currC = he.corner();
            Corner otherC = he.next().corner();
            double otherSol = cornerVals[currC] + gamma_i;
            cornerVals[otherC] = otherSol;
        }

        while (queue.size() > 0) {
            Face curr = queue.back();
            queue.pop_back();

            // Visit adjacent faces.
            for (Halfedge he : curr.adjacentHalfedges()) {
                if (he.edge().isBoundary()) continue;
                Face otherF = he.twin().face();

                if (faceVisited[otherF]) {
                    // orient halfedge so that greater values are on the left side.
                    double diffA = cornerVals[he.corner()] - cornerVals[he.twin().next().corner()];
                    double diffB = cornerVals[he.next().corner()] - cornerVals[he.twin().corner()];
                    double diff = 0.5 * (diffA + diffB);
                    if (abs(diff) < eps) continue;
                    Halfedge cutHalfedge = (diff > 0) ? he : he.twin();
                    branchCut.push_back(cutHalfedge);
                    continue;
                }

                // Fill in corner values for this face.
                for (Halfedge otherHe : otherF.adjacentHalfedges()) {
                    double gamma_i = gamma[geom.edgeIndices[otherHe.edge()]];
                    if (!otherHe.orientation()) gamma_i *= -1;
                    Corner currC = otherHe.corner();
                    Corner otherC = otherHe.next().corner();
                    double otherSol = cornerVals[currC] + gamma_i;
                    cornerVals[otherC] = otherSol;
                }
                // If the corner values along the shared halfedge can be matched via a constant shift, then do that.
                double shiftA = cornerVals[he.corner()] - cornerVals[he.twin().next().corner()];
                double shiftB = cornerVals[he.next().corner()] - cornerVals[he.twin().corner()];
                assert(abs(shiftA - shiftB) < eps);
                for (Corner c : otherF.adjacentCorners()) cornerVals[c] += shiftA;
                faceVisited[otherF] = true;
                queue.push_back(otherF);
            }
        }
    }
    geom.unrequireEdgeIndices();

    // Make sure that this curve has no boundary.
    geom.requireDECOperators();
    SparseMatrix<double> B = geom.d0.transpose();
    Vector<double> chain = Vector<double>::Zero(mesh.nEdges());
    for (const Halfedge& he : branchCut) {
        size_t i = he.edge().getIndex();
        double val = he.orientation() ? 1 : -1;
        chain[i] = val;
    }
    Vector<double> boundary = B * chain;
    std::cerr << "branch cut boundary norm (want this to be 0): " << boundary.norm() << std::endl;
    geom.unrequireDECOperators();

    // Construct 1-chain using all jumps (coefficients may be fractional.)
    Vector<double> allJumpChain = Vector<double>::Zero(mesh.nEdges());
    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        double diffA = cornerVals[he.corner()] - cornerVals[he.twin().next().corner()];
        double diffB = cornerVals[he.next().corner()] - cornerVals[he.twin().corner()];
        double diff = 0.5 * (diffA + diffB);
        allJumpChain[e.getIndex()] = diff;
    }
    std::cerr << "all-jumps chain boundary norm (want this to be 0): " << (B * allJumpChain).norm() << std::endl;

    return cornerVals;
}

/*
 * The jumpPairs formed from the input curve will be a subset of the jumpPairs formed from the completed curve.
 *
 * Returns a vector that maps each jumpPair in the completed curve (v) to the corresponding jumpPair in the input curve
 * (u), and the (absolute value of the) jumpVal in u. *Warning*: This function assumes that wedge pairs are always
 * ordered so that the jump is positive going from w_i to w_j.
 *
 * Not sure of a more efficient way to do this matching...
 */
std::vector<std::pair<size_t, double>>
SurfaceWindingNumbersSolver::matchUpJumpPairs(const std::vector<Eigen::Triplet<double>>& jumpPairs_u,
                                              const CornerData<size_t>& wedgeIndices_u, const size_t& N_u,
                                              const std::vector<Eigen::Triplet<double>>& jumpPairs_v,
                                              const CornerData<size_t>& wedgeIndices_v, const size_t& N_v) {

    // Map each wedge to its set of corners
    std::vector<std::set<size_t>> wedgeToCorners_u(N_u);
    std::vector<std::set<size_t>> wedgeToCorners_v(N_v);
    geom.requireCornerIndices();
    for (Corner c : mesh.corners()) {
        size_t cIdx = geom.cornerIndices[c];
        size_t wIdx_u = wedgeIndices_u[c];
        size_t wIdx_v = wedgeIndices_v[c];

        if (isValidWedgeIndex(wIdx_u)) wedgeToCorners_u[wIdx_u].insert(cIdx);
        if (isValidWedgeIndex(wIdx_v)) wedgeToCorners_v[wIdx_v].insert(cIdx);
    }
    geom.unrequireCornerIndices();

    // Now match up jump pairs in v to jump pairs in u; the following vector keeps track of the i-th jumpPair_u that
    // corresponds to a jumpPair in v, and the (absolute value of the) jumpVal in u. Not sure of a more efficient way to
    // do this...
    std::vector<std::pair<size_t, double>> jumpCaps;
    for (size_t j = 0; j < jumpPairs_v.size(); j++) {
        size_t wi_v = jumpPairs_v[j].row();
        size_t wj_v = jumpPairs_v[j].col();
        bool foundMatch = false;
        for (size_t i = 0; i < jumpPairs_u.size(); i++) {
            size_t wi_u = jumpPairs_u[i].row();
            size_t wj_u = jumpPairs_u[i].col();
            double jumpVal = abs(jumpPairs_u[i].value());
            assert(jumpPairs_u[i].value() > 0);

            bool matchA = (wedgeToCorners_u[wi_u] == wedgeToCorners_v[wi_v]) &&
                          (wedgeToCorners_u[wj_u] == wedgeToCorners_v[wj_v]);
            // bool matchB = (wedgeToCorners_u[wi_u] == wedgeToCorners_v[wj_v]) &&
            //               (wedgeToCorners_u[wj_u] == wedgeToCorners_v[wi_v]);

            if (matchA) {
                jumpCaps.emplace_back(j, jumpVal);
                foundMatch = true;
                break;
            }
        }
        // can't enforce this, otherwise we are restricting v to only jump over wedges in u
        // if (!foundMatch) jumpCaps.emplace_back(j, 0.);
    }
    assert(jumpCaps.size() == jumpPairs_v.size());
    return jumpCaps;
}

CornerData<double> SurfaceWindingNumbersSolver::solve(CornerData<double>& u, const std::vector<Halfedge>& curve,
                                                      bool doHomologyCorrection, bool useSpecialBases) {

    // Index wedges.
    size_t iN, N;
    std::vector<Eigen::Triplet<double>> jumpPairs;
    std::vector<Vertex> bVertices;
    std::vector<Vertex> curveEndpoints;
    CornerData<size_t> wedgeIndices = indexWedges(curve, iN, N, jumpPairs, bVertices, curveEndpoints, useSpecialBases);

    if (psMesh) {
        psMesh->addCornerScalarQuantity("wedge index", wedgeIndices);
    }

    bool simplyConnected = false;
    if (mesh.isManifold() && mesh.nConnectedComponents() == 1) {
        int chi = ((int)mesh.nVertices() + (int)mesh.nFaces()) - (int)mesh.nEdges();
        simplyConnected = (chi == 1 || chi == 2);
    }
    if (simplyConnected) doHomologyCorrection = false;
    std::cout << "simplyConnected: " << (simplyConnected ? "true" : "false") << std::endl;

    // TODO: Subsequent steps be more efficient if I first converted the input curve into a chain (don't have to iterate
    // over possibly many copies of edges, multiple times; can just use the pre-computed chain coefficient.) Downside is
    // that I would have to change the code for indexWedges(), and call convertToChain()


    // Solve for u.
    std::cerr << "Solving jump equation..." << std::endl;
    std::vector<Eigen::Triplet<double>> qJumpPairs = computeCumulativeJumps(curve, bVertices, wedgeIndices);
    std::cerr << "cumulative jump pairs computed" << std::endl;
    u = solveJumpEquation(N, jumpPairs, qJumpPairs, wedgeIndices, curveEndpoints);
    shiftSoMinimumValueisZero(u);
    std::cerr << "jump equation solved" << std::endl;
    std::pair<double, double> uMapRange = getMinMax(u);
    double uRange = uMapRange.second - uMapRange.first;
    std::cerr << "u min: " << uMapRange.first << "\tu max: " << uMapRange.second << std::endl;
    if (psMesh) psMesh->addCornerScalarQuantity("u", u)->setMapRange(uMapRange);
    if (psMesh) plotIntegerRoundedFunction(u, curve, psMesh, "u");

    if (!doHomologyCorrection) {
        std::pair<double, double> wMapRange = getMinMax(u);
        std::cerr << "w min: " << wMapRange.first << "\tw max: " << wMapRange.second << std::endl;

        if (psMesh) plotIntegerRoundedFunction(u, curve, psMesh, "w");

        return u;
    }

    CornerData<double> w;
    if (isCurveClosed(curve)) {
        Vector<double> omega = computeOmega(u, curveEndpoints);
        std::cerr << "|Du|: " << omega.norm() << std::endl;
        if (psMesh) psMesh->addOneFormIntrinsicVectorQuantity("Du", omega, polyscopeEdgeOrientations(mesh));
        Vector<double> v_vector = integrateGammaAndMinimizeJumpsL1(curve, omega, wedgeIndices);
        CornerData<double> v_corners(mesh, cornerVector(v_vector, wedgeIndices));
        shiftSoMinimumValueisZero(v_corners);
        if (psMesh) psMesh->addCornerScalarQuantity("v", v_corners);
        w = u - v_corners;
        if (psMesh) psMesh->addCornerScalarQuantity("w", w);
    } else {
        Vector<double> omega, dAlpha, deltaBeta, gamma;
        computeOmegaAndHodgeDecomposition(u, curveEndpoints, omega, dAlpha, deltaBeta, gamma);
        std::cerr << "Hodge decomposition done" << std::endl;
        if (psMesh) psMesh->addOneFormIntrinsicVectorQuantity("gamma", gamma, polyscopeEdgeOrientations(mesh));

        CornerData<double> v_corners = integrateExactlyBranchCutAnywhere(gamma, curve);
        shiftSoMinimumValueisZero(v_corners);
        if (psMesh) psMesh->addCornerScalarQuantity("v", v_corners);
        if (psMesh) psMesh->addCornerScalarQuantity("u-v", u - v_corners);

        std::vector<Eigen::Triplet<double>> newJumpPairs = getFinalJumps(qJumpPairs, v_corners, wedgeIndices, N);
        std::cerr << "Solving final jump equation..." << std::endl;
        w = solveJumpEquation(N, jumpPairs, newJumpPairs, wedgeIndices, curveEndpoints);
        shiftSoMinimumValueisZero(w);
        std::pair<double, double> wMapRange = getMinMax(w);
        wMapRange.second = wMapRange.first + uRange;
        if (psMesh) psMesh->addCornerScalarQuantity("w", w)->setMapRange(wMapRange);
        std::cerr << "w min: " << wMapRange.first << "\tw max: " << wMapRange.second << std::endl;
    }

    std::pair<double, double> wMapRange = getMinMax(w);
    std::cerr << "w min: " << wMapRange.first << "\tw max: " << wMapRange.second << std::endl;

    if (psMesh) plotIntegerRoundedFunction(w, curve, psMesh, "w");


    return w;
}

CornerData<double> SurfaceWindingNumbersSolver::solveFinalAndReturnAllData(
    const std::vector<Halfedge>& curve, CornerData<double>& u, CornerData<double>& v_corners,
    CornerData<double>& w_rounded, Vector<double>& gamma, bool doHomologyCorrection, bool useSpecialBases) {

    // Index wedges.
    size_t iN, N;
    std::vector<Eigen::Triplet<double>> jumpPairs;
    std::vector<Vertex> bVertices;
    std::vector<Vertex> curveEndpoints;
    CornerData<size_t> wedgeIndices = indexWedges(curve, iN, N, jumpPairs, bVertices, curveEndpoints, useSpecialBases);
    CornerData<double> w;
    gamma = Vector<double>::Zero(mesh.nEdges());

    if (isCurveClosed(curve) || !doHomologyCorrection) {

        // Solve for u.
        std::vector<Eigen::Triplet<double>> qJumpPairs = computeCumulativeJumps(curve, bVertices, wedgeIndices);
        std::cerr << "cumulative jump pairs computed" << std::endl;
        u = solveJumpEquation(N, jumpPairs, qJumpPairs, wedgeIndices, curveEndpoints);
        if (psMesh) psMesh->addCornerScalarQuantity("u", u);
        w = u;
        v_corners = u;

        if (doHomologyCorrection) {
            // Compute h := Du, and integrate h s.t. the L1 norm of the length of the jumps are minimized.
            Vector<double> omega = computeOmega(u, curveEndpoints);
            if (psMesh) psMesh->addOneFormIntrinsicVectorQuantity("Du", omega, polyscopeEdgeOrientations(mesh));
            gamma = omega;
            Vector<double> v_vector = integrateGammaAndMinimizeJumpsL1(curve, omega, wedgeIndices);
            v_corners.fromVector(cornerVector(v_vector, wedgeIndices));
            if (psMesh) psMesh->addCornerScalarQuantity("v", v_corners);
            w = u - v_corners;
        }
    } else if (doHomologyCorrection) {

        // Compute u just for visual debugging purposes (we don't actually need to compute it)
        std::vector<Eigen::Triplet<double>> qJumpPairs = computeCumulativeJumps(curve, bVertices, wedgeIndices);
        std::cerr << "cumulative jump pairs computed" << std::endl;
        u = solveJumpEquation(N, jumpPairs, qJumpPairs, wedgeIndices, curveEndpoints);
        std::pair<double, double> uMapRange = getMinMax(u);
        double uRange = uMapRange.second - uMapRange.first;
        std::cerr << "u min: " << uMapRange.first << "\tu max: " << uMapRange.second << std::endl;
        // if (psMesh) psMesh->addCornerScalarQuantity("u", u)->setMapRange(uMapRange);


        Vector<double> omega, dAlpha, deltaBeta;
        computeOmegaAndHodgeDecomposition(u, curveEndpoints, omega, dAlpha, deltaBeta, gamma);

        // if (psMesh) psMesh->addOneFormIntrinsicVectorQuantity("gamma", gamma, polyscopeEdgeOrientations(mesh));

        v_corners = integrateExactlyBranchCutAnywhere(gamma, curve);
        if (psMesh) psMesh->addCornerScalarQuantity("v", v_corners);
        if (psMesh) psMesh->addCornerScalarQuantity("u-v", u - v_corners);

        std::cerr << "Getting final jumps..." << std::endl;
        std::vector<Eigen::Triplet<double>> newJumpPairs = getFinalJumps(qJumpPairs, v_corners, wedgeIndices, N);
        std::cerr << "Solving final jump equation..." << std::endl;
        CornerData<double> w_corners = solveJumpEquation(N, jumpPairs, newJumpPairs, wedgeIndices, curveEndpoints);
        w = w_corners;
    }

    std::pair<double, double> wMapRange = getMinMax(w);
    std::cerr << "w min: " << wMapRange.first << "\tw max: " << wMapRange.second << std::endl;

    w_rounded = integerRounded(w, curve);
    // w_rounded = integerRounded(u, curve);

    return w;
}

/*
 * Same as the other solve() function that solves the Poisson equation, but when the curve is already given in the
 * form of BarycentricVectors.
 */
VertexData<double> SurfaceWindingNumbersSolver::solve(const std::vector<BarycentricVector>& curve,
                                                      bool doHomologyCorrection, bool useWhitneyElements) {
    curveVectors = curve;

    VertexData<double> w(mesh);

    geom.requireEdgeIndices();
    geom.requireVertexIndices();
    geom.requireFaceAreas();

    if (!useWhitneyElements) {
        w = solveCrouzeixRaviartPoissonEquation();
    } else {
        w = solvePoissonEquation();
    }

    if (doHomologyCorrection) {

        // Compute ω.
        geom.requireFaceIndices();
        Vector<double> omega = Vector<double>::Zero(mesh.nEdges());
        for (Edge e : mesh.edges()) {
            size_t i = geom.edgeIndices[e];
            double theta_i = 2. * PI * w[e.firstVertex()];
            double theta_j = 2. * PI * w[e.secondVertex()];

            std::complex<double> z_j, z_i;
            z_j = Vector2::fromAngle(theta_j);
            z_i = Vector2::fromAngle(theta_i);
            omega[i] = std::arg(z_j / z_i);
        }

        // Solve for γ.
        ensureHaveHodgeDecompositionSolvers();
        Vector<double> dAlpha = computeExactComponent(omega);
        Vector<double> deltaBeta = computeCoExactComponent(omega);
        Vector<double> gamma = omega - dAlpha - deltaBeta;
        std::vector<BarycentricVector> gammaVecs = gradient(gamma);

        if (!useWhitneyElements) {
            // TODO
        } else {
            // TODO
        }
        geom.unrequireFaceIndices();
    }

    geom.unrequireEdgeIndices();
    geom.unrequireVertexIndices();
    geom.unrequireFaceAreas();

    return w;
}

/*
 * Solve SWN on a (possibly non-manifold) surface mesh via a Poisson equation, with or without homology correction
 * (true by default), or using Whitney or Crouzeix-Raviart elements (C-R by default.)
 *
 * <curveEdges> indexes into <curveNodes>; each element describes an edge represented by two nodes, which *must* share
 * a common face.
 *
 * TODO: Possibly provide the option to allow mesh mutation, and re-mesh s.t. the curve lies on edges.
 */
VertexData<double> SurfaceWindingNumbersSolver::solve(const std::vector<SurfacePoint>& curveNodes,
                                                      const std::vector<std::array<size_t, 2>>& curveEdges,
                                                      bool doHomologyCorrection, bool useWhitneyElements) {

    std::vector<BarycentricVector> curve = convertCurveToBarycentricVectors(curveNodes, curveEdges);

    return solve(curve, doHomologyCorrection, useWhitneyElements);
}