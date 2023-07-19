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