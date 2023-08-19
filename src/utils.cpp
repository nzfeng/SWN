#include "utils.h"

// ===================== NUMERICAL

int mod(int a, int b) {
    return (b + (a % b)) % b;
}

/* Round to nearest int. Add +- 0.5 because compiler always truncates. */
int roundToNearestInteger(double x) {

    double y = x + 0.5 - (x < 0);
    return (int)y;
}

bool isSpecial(double x) {
    return std::isnan(x);
}

std::tuple<double, double> minMax(const Vector<double>& vec) {

    const double inf = std::numeric_limits<double>::infinity();
    double maxVal = -inf;
    double minVal = inf;
    for (size_t i = 0; i < vec.size(); i++) {
        double val = vec[i];
        if (isSpecial(val)) continue;
        maxVal = std::max(maxVal, val);
        minVal = std::min(minVal, val);
    }
    return std::make_tuple(minVal, maxVal);
}

/* Compute the global shift (part of the rounding procedure) described in Section 3.5.1. */
double shift(const CornerData<double>& func, const std::vector<Halfedge>& curve) {

    SurfaceMesh* mesh = func.getMesh();
    Vector<int> chain = convertToChain(*mesh, curve);
    double avgShiftOnOneSide = 0.;
    size_t nValidJumps = 0;
    for (Edge e : mesh->edges()) {
        size_t eIdx = e.getIndex();
        if (!e.isBoundary() && chain[eIdx] != 0) {
            Halfedge he = (chain[eIdx] > 0) ? e.halfedge() : e.halfedge().twin();
            Corner cA = he.corner();
            Corner cB = he.next().corner();
            if (isSpecial(func[cA]) || isSpecial(func[cB])) continue;
            double avg = 0.5 * (func[cA] + func[cB]);
            double avgInt = roundToNearestInteger(avg);
            double shift = abs(avgInt - avg);
            avgShiftOnOneSide += shift;
            nValidJumps += 1;
        }
    }
    avgShiftOnOneSide = nValidJumps == 0 ? 0 : avgShiftOnOneSide / nValidJumps;
    return avgShiftOnOneSide;
}

/* Implement the rounding procedure described in Section 3.5.1. */
FaceData<double> round(const CornerData<double>& func, const std::vector<Halfedge>& curve) {

    // Compute global shift.
    double avgShiftOnOneSide = shift(func, curve);

    // Threshold on a per-face basis
    SurfaceMesh* mesh = func.getMesh();
    FaceData<double> F(*mesh);
    for (Face f : mesh->faces()) {
        double avgVal = 0.;
        double nValidCorners = 0.;
        for (Corner c : f.adjacentCorners()) {
            if (isSpecial(func[c])) continue;
            avgVal += func[c] + avgShiftOnOneSide;
            nValidCorners += 1;
        }
        avgVal /= nValidCorners;
        F[f] = roundToNearestInteger(avgVal);
    }
    return F;
}


// ===================== OPERATORS


/* Apparently geometry-central can only build D as part of geom.DECoperators(), but strictly speaking D doesn't need
 * geometry. So just build the matrix myself here. */
template <typename T>
SparseMatrix<T> b1(SurfaceMesh& mesh) {
    SparseMatrix<T> B(mesh.nVertices(), mesh.nEdges());
    std::vector<Eigen::Triplet<T>> tripletList;
    for (Edge e : mesh.edges()) {
        size_t eIdx = e.getIndex();
        tripletList.emplace_back(e.firstVertex().getIndex(), eIdx, -1);
        tripletList.emplace_back(e.secondVertex().getIndex(), eIdx, 1);
    }
    B.setFromTriplets(tripletList.begin(), tripletList.end());
    return B;
}

template <typename T>
SparseMatrix<T> b2(SurfaceMesh& mesh) {
    SparseMatrix<T> B(mesh.nEdges(), mesh.nFaces());
    std::vector<Eigen::Triplet<T>> tripletList;
    for (Face f : mesh.faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            tripletList.emplace_back(he.edge().getIndex(), f.getIndex(), (he.orientation() ? 1 : -1));
        }
    }
    B.setFromTriplets(tripletList.begin(), tripletList.end());
    return B;
}

CornerData<Vector2> toProjectiveCoordinates(const CornerData<double>& u) {

    SurfaceMesh* mesh = u.getMesh();
    CornerData<Vector2> u_hom(*mesh);
    for (const Corner& c : mesh->corners()) {
        if (isSpecial(u[c])) {
            u_hom[c] = {0., 0.};
        } else {
            u_hom[c] = {u[c], 1.0};
        }
    }
    return u_hom;
}

CornerData<double> fromProjectiveCoordinates(const CornerData<Vector2>& u_hom) {

    SurfaceMesh* mesh = u_hom.getMesh();
    CornerData<double> u(*mesh);
    double eps = 1e-5;
    for (const Corner& c : mesh->corners()) {
        u[c] = (abs(u_hom[c][1]) < eps) ? SPECIAL_VAL : u_hom[c][0] / u_hom[c][1];
    }
    return u;
}

// ===================== I/O

std::string getHomeDirectory(const std::string& filepath) {

    std::string dir(filepath.begin(), filepath.begin() + filepath.find_last_of("/") + 1);
    return dir;
}

bool isStringTrue(const std::string& input) {

    std::string data = input;
    std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c) { return std::tolower(c); });
    return (data == "true" || data == "on");
}

/* Read in curve, encoded as OBJ line objects. */
void readLines(SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
               std::vector<std::array<size_t, 2>>& curveEdges, std::vector<std::array<Face, 2>>& dualChain,
               int offset) {

    std::ifstream curr_file(filepath.c_str());
    std::string line;
    std::string X; // re-written variable for holding the leading (indicator) char
    size_t idx;    // re-written variable for holding element index

    if (curr_file.is_open()) {

        while (!curr_file.eof()) {

            getline(curr_file, line);

            // Ignore any newlines
            if (line == "") {
                continue;
            }

            std::istringstream iss(line);
            iss >> X; // get the first char in the line

            if (X == "l") {
                std::vector<SurfacePoint> nodes;
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx += offset; // OBJ elements are 1-indexed, whereas geometry-central uses 0-indexing
                    nodes.emplace_back(mesh.vertex(idx));
                }
                size_t N = nodes.size();
                size_t M = curveNodes.size();
                for (size_t i = 0; i < N - 1; i++) {
                    curveEdges.push_back({M + i, M + (i + 1)});
                }
                curveNodes.insert(curveNodes.end(), nodes.begin(), nodes.end());
            } else if (X == "c") {
                std::vector<size_t> indices;
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx += offset;
                    indices.push_back(idx);
                }
                size_t N = indices.size();
                for (size_t i = 0; i < N - 1; i++) {
                    dualChain.push_back({mesh.face(indices[i]), mesh.face(indices[i + 1])});
                }
            } else if (X == "#offset") {
                iss >> offset;
            }
        }
        curr_file.close();
    } else {
        std::cerr << "Could not open file <" << filepath << ">" << std::endl;
    }
}

/* Read in curve data from a separate file, encoded as segments between barycentric points. */
void readCurves(SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
                std::vector<std::array<size_t, 2>>& curveEdges, std::vector<std::array<Face, 2>>& dualChain,
                int offset) {

    std::ifstream curr_file(filepath.c_str());
    std::string line;
    std::string X; // re-written variable for holding the leading (indicator) char
    size_t idx;    // re-written variable for holding element index
    double a, b, c;

    if (curr_file.is_open()) {

        while (!curr_file.eof()) {

            getline(curr_file, line);

            // Ignore any newlines
            if (line == "") {
                continue;
            }

            std::istringstream iss(line);
            iss >> X; // get the first char in the line

            if (X == "v") {
                iss >> idx;
                idx += offset;
                curveNodes.emplace_back(mesh.vertex(idx));
            } else if (X == "e") {
                // WARNING: The edge index is assumed to agree with geometry-central's edge-indexing scheme. Depending
                // on which software you used to export curve data on the mesh, your edge indices may disagree with the
                // indices geometry-central computes from the vertices.
                iss >> idx >> a;
                idx += offset;
                curveNodes.emplace_back(mesh.edge(idx), a);
            } else if (X == "f") {
                iss >> idx >> a >> b >> c;
                idx += offset;
                curveNodes.emplace_back(mesh.face(idx), Vector3{a, b, c});
            } else if (X == "l") {
                std::vector<size_t> indices;
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx += offset;
                    indices.push_back(idx);
                }
                size_t N = indices.size();
                for (size_t i = 0; i < N - 1; i++) {
                    curveEdges.push_back({indices[i], indices[i + 1]});
                }
            } else if (X == "c") {
                std::vector<size_t> indices;
                while (true) {
                    if (iss.eof()) break;
                    iss >> idx;
                    idx += offset;
                    indices.push_back(idx);
                }
                size_t N = indices.size();
                for (size_t i = 0; i < N - 1; i++) {
                    dualChain.push_back({mesh.face(indices[i]), mesh.face(indices[i + 1])});
                }
            } else if (X == "#offset") {
                iss >> offset;
            }
        }
        curr_file.close();
    } else {
        std::cerr << "Could not open file <" << filepath << ">" << std::endl;
    }
}

void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<SurfacePoint>& curveNodes,
                  const std::vector<std::vector<std::array<size_t, 2>>>& curveEdges, const std::string& filename) {

    std::fstream f;
    f.open(filename, std::ios::out | std::ios::trunc);

    Vector3 pos;
    size_t a, b;
    if (f.is_open()) {
        for (const SurfacePoint& pt : curveNodes) {
            pos = pt.interpolate(vertexPositions);
            f << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }
        for (auto curve : curveEdges) {
            // Assume that each curve is connected
            f << "l ";
            for (auto seg : curve) {
                a = seg[0];
                b = seg[1];
                f << seg[0] + 1 << " "; // OBJ is 1-indexed
            }
            f << b + 1; // OBJ is 1-indexed
            f << "\n";
        }

        f.close();
        std::cerr << "File " << filename << " written succesfully." << std::endl;
    } else {
        std::cerr << "Could not save curves '" << filename << "'!" << std::endl;
    }
}

void exportCurves(const VertexData<Vector3>& vertexPositions, const std::vector<std::vector<Halfedge>>& curveHalfedges,
                  const std::string& filename) {

    std::vector<SurfacePoint> curveNodes;
    std::vector<std::vector<std::array<size_t, 2>>> curveEdges;
    for (const auto& curve : curveHalfedges) {
        curveEdges.emplace_back();
        auto& currCurve = curveEdges.back();
        for (const Halfedge& he : curve) {
            curveNodes.emplace_back(he.tailVertex());
            curveNodes.emplace_back(he.tipVertex());
            size_t N = curveNodes.size();
            currCurve.push_back({N - 2, N - 1});
        }
    }
    exportCurves(vertexPositions, curveNodes, curveEdges, filename);
}

void exportFunction(EmbeddedGeometryInterface& geom, const CornerData<double>& u, const std::string& filename) {

    CornerData<Vector2> texCoords = toProjectiveCoordinates(u);
    WavefrontOBJ::write(filename, geom, texCoords);
    std::cerr << "File " << filename << " written." << std::endl;
}

/* Export on common subdivision. */
void exportFunction(IntegerCoordinatesIntrinsicTriangulation& intTri, VertexPositionGeometry& manifoldGeom,
                    const CornerData<double>& u, const std::string& filename) {

    CornerData<Vector2> texCoords = toProjectiveCoordinates(u);

    CommonSubdivision& cs = intTri.getCommonSubdivision();
    cs.constructMesh(true, true);

    // Interpolate homogenous coordinates from intrinsic mesh to the common subdivision
    CornerData<Vector2> cs_u = interpolateVector2AcrossB(cs, texCoords);
    cs_u = toProjectiveCoordinates(fromProjectiveCoordinates(cs_u));

    // Export mesh
    ManifoldSurfaceMesh& csMesh = *cs.mesh;
    VertexData<Vector3> csPositions = cs.interpolateAcrossA(manifoldGeom.vertexPositions);
    VertexPositionGeometry geom(csMesh, csPositions);
    writeSurfaceMesh(csMesh, geom, cs_u, filename);
    std::cerr << "File " << filename << " written." << std::endl;
}


// ===================== CURVE MANIPULATION

Halfedge determineHalfedgeFromVertices(const Vertex& vA, const Vertex& vB) {

    for (Halfedge he : vA.outgoingHalfedges()) {
        if (he.tipVertex() == vB) return he;
    }
    return Halfedge();
}

/*
 * A SurfacePoint on the original input mesh, re-interpret to a new mesh with the same topology.
 */
SurfacePoint reinterpretTo(const SurfacePoint& p, SurfaceMesh& otherMesh) {

    switch (p.type) {
        case (SurfacePointType::Vertex):
            return SurfacePoint(otherMesh.vertex(p.vertex.getIndex()));
            break;
        case (SurfacePointType::Edge):
            return SurfacePoint(otherMesh.edge(p.edge.getIndex()), p.tEdge);
            break;
        case (SurfacePointType::Face):
            return SurfacePoint(otherMesh.face(p.face.getIndex()), p.faceCoords);
            break;
    }
    throw std::logic_error("Bad switch");
}

SurfaceMesh* getMesh(const SurfacePoint& p) {
    switch (p.type) {
        case (SurfacePointType::Vertex): {
            return p.vertex.getMesh();
            break;
        }
        case (SurfacePointType::Edge): {
            return p.edge.getMesh();
            break;
        }
        case (SurfacePointType::Face): {
            return p.face.getMesh();
            break;
        }
    }
    return nullptr;
}

/*
 * Convert curve represented using barycentric points to a collection of halfedges, if all curve edges are constrained
 * to mesh edges.
 *
 * Return an empty vector if curve is not constrained to mesh edges.
 */
std::vector<Halfedge> convertToHalfedges(const std::vector<SurfacePoint>& curveNodes,
                                         const std::vector<std::array<size_t, 2>>& curveEdges) {

    std::vector<Halfedge> curveHalfedges;
    if (curveNodes.size() == 0 || curveEdges.size() == 0) return curveHalfedges;

    SurfaceMesh* mesh = getMesh(curveNodes[0]);

    if (mesh->isOriented()) {
        for (const auto& seg : curveEdges) {
            const SurfacePoint& pt1 = curveNodes[seg[0]];
            const SurfacePoint& pt2 = curveNodes[seg[1]];

            // For some reason, very occasionally this code fails. I checked that pt2 is indeed among the adjacent
            // vertices of pt1; the == operator among vertices just doesn't work sometimes, I guess.
            Edge e = sharedEdge(pt1, pt2);
            if (e == Edge()) {
                std::cerr << "Input curve segment between <" << pt1 << "> and <" << pt2
                          << "> is not constrained to mesh edges." << std::endl;
                return std::vector<Halfedge>();
            }
            if (pt1.type != SurfacePointType::Vertex || pt2.type != SurfacePointType::Vertex)
                return std::vector<Halfedge>();

            Halfedge he = (pt1.vertex == e.firstVertex()) ? e.halfedge() : e.halfedge().twin();
            curveHalfedges.push_back(he);

            // if (pt1.type != SurfacePointType::Vertex || pt2.type != SurfacePointType::Vertex) {
            //     return std::vector<Halfedge>();
            // }

            // Halfedge he = Halfedge();
            // for (Halfedge heA : pt1.vertex.outgoingHalfedges()) {
            //     if (heA.tipVertex() == pt2.vertex) {
            //         he = heA;
            //         break;
            //     }
            // }

            // if (he != Halfedge()) {
            //     curveHalfedges.push_back(he);
            // } else {
            //     std::cerr << "Input curve segment between <" << pt1 << "> and <" << pt2
            //               << "> is not constrained to mesh edges." << std::endl;
            //     return std::vector<Halfedge>();
            // }
        }
    } else {
        for (const auto& seg : curveEdges) {
            const SurfacePoint& pt1 = curveNodes[seg[0]];
            const SurfacePoint& pt2 = curveNodes[seg[1]];

            if (pt1.type != SurfacePointType::Vertex || pt2.type != SurfacePointType::Vertex) {
                return std::vector<Halfedge>();
            }

            bool found = false;
            for (Edge e : mesh->edges()) {
                if (e.firstVertex() == pt1.vertex && e.secondVertex() == pt2.vertex) {
                    curveHalfedges.push_back(e.halfedge());
                    found = true;
                    break;
                } else if (e.firstVertex() == pt2.vertex && e.secondVertex() == pt1.vertex) {
                    curveHalfedges.push_back(e.halfedge().twin());
                    found = true;
                    break;
                }
            }

            if (!found) {
                std::cerr << "Input curve is not constrained to mesh edges." << std::endl;
                return std::vector<Halfedge>();
            }
        }
    }

    return curveHalfedges;
}

Vector<int> convertToChain(const SurfaceMesh& mesh, const std::vector<Halfedge>& curve) {

    Vector<int> chain = Vector<int>::Zero(mesh.nEdges());
    for (const Halfedge& he : curve) {
        chain[he.edge().getIndex()] += he.orientation() ? 1 : -1;
    }
    return chain;
}

/*
 * Given a collection of halfedges, determine all the connected curve components it represents.
 *
 * This function is useful for extracting smooth curves when we go to export curves as OBJ and render them. If the
 * curve is exported as a collection of (disconnected) edges, then the resulting curve when thickened and rendered
 * in Blender will have unsightly gaps.
 */
std::vector<std::vector<Halfedge>> getCurveComponents(IntrinsicGeometryInterface& geom,
                                                      const std::vector<Halfedge>& curveHalfedges, bool useEndpoints) {

    SurfaceMesh& mesh = geom.mesh;
    Vector<int> chain = convertToChain(mesh, curveHalfedges);
    SparseMatrix<int> B = b1<int>(mesh);
    Vector<int> boundary = B * chain;
    double eps = 1e-5;
    size_t E = mesh.nEdges();
    std::vector<std::vector<Halfedge>> curves;
    geom.requireEdgeIndices();
    while (chain.norm() > eps) {
        for (size_t i = 0; i < E; i++) {
            if (chain[i] == 0) continue;

            // Get a starting halfedge.
            Halfedge startHe = mesh.edge(i).halfedge();
            int sign = 1;
            if (chain[i] < 0) {
                startHe = startHe.twin();
                sign = -1;
            }
            chain[i] -= sign;

            curves.emplace_back();
            std::vector<Halfedge>& currCurve = curves.back();
            currCurve.push_back(startHe);

            // Attach edges to endpoints until we can't.
            Halfedge currHe = startHe;
            while (true) {
                Vertex v = currHe.tipVertex();
                if (useEndpoints && boundary[v.getIndex()] != 0) break;
                bool didWeFindOne = false;
                for (Halfedge he : v.outgoingHalfedges()) {
                    size_t eIdx = geom.edgeIndices[he.edge()];
                    bool isValid = (chain[eIdx] > 0 && he.orientation()) || (chain[eIdx] < 0 && !he.orientation());
                    if (isValid) {
                        currHe = he;
                        currCurve.push_back(he);
                        chain[eIdx] -= (chain[eIdx] > 0) ? 1 : -1;
                        didWeFindOne = true;
                        break;
                    }
                }
                if (!didWeFindOne) break;
            }

            currHe = startHe;
            while (true) {
                Vertex v = currHe.tailVertex();
                if (useEndpoints && boundary[v.getIndex()] != 0) break;
                bool didWeFindOne = false;
                for (Halfedge he : v.incomingHalfedges()) {
                    size_t eIdx = geom.edgeIndices[he.edge()];
                    bool isValid = (chain[eIdx] > 0 && he.orientation()) || (chain[eIdx] < 0 && !he.orientation());
                    if (isValid) {
                        currHe = he;
                        currCurve.insert(currCurve.begin(), he);
                        chain[eIdx] -= (chain[eIdx] > 0) ? 1 : -1;
                        didWeFindOne = true;
                        break;
                    }
                }
                if (!didWeFindOne) break;
            }
        }
    }
    geom.unrequireEdgeIndices();

    return curves;
}

/*
 * Same as getCurveComponents(), but for when the curve is not constrained to mesh edges and is given as a series of
 * SurfacePoints.
 *
 * Warning: This function assumes that curve nodes have the same index iff they have the same position.
 */
std::vector<std::vector<std::array<size_t, 2>>>
getCurveComponents(SurfaceMesh& mesh, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges) {

    std::vector<std::array<size_t, 2>> edgesToAdd = curveEdges;
    std::vector<std::vector<std::array<size_t, 2>>> curves;
    size_t nSegs = curveEdges.size();
    while (edgesToAdd.size() > 0) {
        std::array<size_t, 2> startSeg = edgesToAdd.back();
        edgesToAdd.pop_back();
        curves.emplace_back();
        std::vector<std::array<size_t, 2>>& currCurve = curves.back();
        currCurve.push_back(startSeg);

        // Add segs to the front end until we can't.
        std::array<size_t, 2> currSeg = startSeg;
        while (true) {
            const SurfacePoint& front = curveNodes[currSeg[1]];
            bool didWeFindOne = false;
            for (size_t i = 0; i < edgesToAdd.size(); i++) {
                std::array<size_t, 2> otherSeg = edgesToAdd[i];
                if (curveNodes[otherSeg[0]] == front) {
                    currSeg = otherSeg;
                    currCurve.push_back(otherSeg);
                    edgesToAdd.erase(edgesToAdd.begin() + i);
                    didWeFindOne = true;
                    break;
                }
            }
            if (!didWeFindOne) break;
        }

        // Add segs to the back end until we can't.
        currSeg = startSeg;
        while (true) {
            const SurfacePoint& back = curveNodes[currSeg[0]];
            bool didWeFindOne = false;
            for (size_t i = 0; i < edgesToAdd.size(); i++) {
                std::array<size_t, 2> otherSeg = edgesToAdd[i];
                if (curveNodes[otherSeg[1]] == back) {
                    currSeg = otherSeg;
                    currCurve.insert(currCurve.begin(), otherSeg);
                    edgesToAdd.erase(edgesToAdd.begin() + i);
                    didWeFindOne = true;
                    break;
                }
            }
            if (!didWeFindOne) break;
        }
    }
    return curves;
}

/*
 * Given the input curve and computed residual function, return:
 *     - the bounding parts of the input curve
 *     - the nonbounding parts of the input curve
 * as edges in the mesh.
 */
std::tuple<std::vector<Halfedge>, std::vector<Halfedge>>
getCurveDecomposition(const std::vector<Halfedge>& curveHalfedges, const CornerData<double>& vFunc, double epsilon) {

    SurfaceMesh* mesh = vFunc.getMesh();
    if (mesh == nullptr) return std::make_tuple(std::vector<Halfedge>(), std::vector<Halfedge>());

    std::vector<Halfedge> nonbound, bound;
    std::vector<Halfedge> nbLoops = getJumpLocus(curveHalfedges, vFunc, epsilon);
    Vector<int> chain = convertToChain(*mesh, nbLoops);
    for (const Halfedge& he : curveHalfedges) {
        int coeff = chain[he.edge().getIndex()];
        if (coeff == 0 || (sgn(coeff) > 0) != he.orientation()) {
            bound.push_back(he);
        } else {
            nonbound.push_back(he);
        }
    }
    return std::make_tuple(bound, nonbound);
}

std::vector<Halfedge> getJumpLocus(const std::vector<Halfedge>& curveHalfedges, const FaceData<double>& func,
                                   double epsilon) {

    SurfaceMesh* mesh = func.getMesh();
    if (mesh == nullptr) return std::vector<Halfedge>();

    // Get (oriented) boundary.
    SparseMatrix<double> B = b2<double>(*mesh);
    Vector<double> w = func.toVector();
    Vector<double> Bw = B * w;

    std::vector<Halfedge> curve;
    for (Edge e : mesh->edges()) {
        double c = Bw[e.getIndex()];
        Halfedge he = e.halfedge();
        if (abs(c) > epsilon) curve.push_back(c > 0. ? he : he.twin());
    }
    return curve;
}

std::vector<Halfedge> getJumpLocus(const std::vector<Halfedge>& curveHalfedges, const CornerData<double>& func,
                                   double epsilon) {

    SurfaceMesh* mesh = func.getMesh();
    if (mesh == nullptr) return std::vector<Halfedge>();

    // Get (oriented) boundary.
    std::vector<Halfedge> curve;
    for (Edge e : mesh->edges()) {
        if (e.isBoundary()) continue;
        Halfedge he = e.halfedge();
        double funcA = 0.5 * (func[he.corner()] + func[he.next().corner()]);
        double funcB = 0.5 * (func[he.twin().corner()] + func[he.twin().next().corner()]);
        double c = funcA - funcB;
        if (abs(c) > epsilon) curve.push_back(c > 0. ? he : he.twin());
    }
    return curve;
}

/*
 * Given the input curve and a function, extract the function's jump locus (restricted to mesh edges) as a curve.
 */
std::vector<Halfedge> getCompletedBoundingLoops(const std::vector<Halfedge>& curveHalfedges,
                                                const CornerData<double>& func, double epsilon) {
    SurfaceMesh* mesh = func.getMesh();
    if (mesh == nullptr) return std::vector<Halfedge>();

    // Apply rounding procedure.
    FaceData<double> W = round(func, curveHalfedges);
    return getJumpLocus(curveHalfedges, W, epsilon);
}

// ===================== MESH MUTATION

/*
 * Build a new intrinsic triangulatiion, and set the given curve halfedges as marked edges.
 */
void resetIntrinsicTriangulationAndMarkEdges(std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>& intTri,
                                             ManifoldSurfaceMesh& manifoldMesh, VertexPositionGeometry& manifoldGeom,
                                             const std::vector<Halfedge>& curveHalfedgesOnManifold) {

    intTri = std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation>(
        new IntegerCoordinatesIntrinsicTriangulation(manifoldMesh, manifoldGeom));

    // If curve is given as a set of mesh edges, set these edges as "marked" in the intrinsic triangulation so they
    // never get flipped.
    manifoldGeom.requireEdgeIndices();
    EdgeData<bool> markedEdges(*intTri->intrinsicMesh, false);
    // The intrinsic triangulation should at first just be a copy of the input mesh.
    for (Halfedge he : curveHalfedgesOnManifold) {
        size_t eIdx = manifoldGeom.edgeIndices[he.edge()];
        markedEdges[intTri->intrinsicMesh->edge(eIdx)] = true;
    }
    manifoldGeom.unrequireEdgeIndices();
    intTri->setMarkedEdges(markedEdges);
}

/*
 * Given a curve represented as a series of halfedges in the original input mesh, return the set of corresponding
 * halfedges in the (possibly refined) intrinsic triangulation.
 */
std::vector<Halfedge>
determineHalfedgesInIntrinsicTriangulation(IntegerCoordinatesIntrinsicTriangulation& intTri,
                                           const std::vector<Halfedge>& curveHalfedgesOnManifold) {

    // For some reason, need to call this first or else everything fails. Intrinsic mesh isn't built/indexed right.
    intTri.traceAllIntrinsicEdgesAlongInput();

    // Use traceInputHalfedgeAlongIntrinsic(), and verify that all SurfacePoints are vertex-type.
    std::vector<Halfedge> intCurve;
    for (size_t heIdx = 0; heIdx < curveHalfedgesOnManifold.size(); heIdx++) {
        Halfedge he = curveHalfedgesOnManifold[heIdx];
        std::vector<SurfacePoint> pts = intTri.traceInputHalfedgeAlongIntrinsic(he);
        for (SurfacePoint pt : pts) {
            assert(pt.type == SurfacePointType::Vertex);
        }
        size_t n = pts.size();
        for (size_t i = 0; i < n - 1; i++) {
            Vertex vA = pts[i].vertex;
            Vertex vB = pts[i + 1].vertex;
            for (Halfedge outHe : vA.outgoingHalfedges()) {
                if (outHe.tipVertex() == vB) {
                    intCurve.push_back(outHe);
                    break;
                }
            }
        }
    }
    return intCurve;
}


// ===================== VISUALIZATION

void displayCurves(const VertexPositionGeometry& geometry, const std::vector<Halfedge>& curveHalfedges,
                   const std::vector<SurfacePoint>& curveNodes, const std::vector<std::array<size_t, 2>>& curveEdges,
                   const std::vector<std::array<Face, 2>>& dualChain) {

    if (curveEdges.size() > 0 && &geometry.mesh == getMesh(curveNodes[0])) {
        std::vector<Vector3> nodes;
        std::vector<std::array<size_t, 2>> edges;
        for (const auto& pair : curveEdges) {
            size_t N = nodes.size();
            nodes.push_back(curveNodes[pair[0]].interpolate(geometry.vertexPositions));
            nodes.push_back(curveNodes[pair[1]].interpolate(geometry.vertexPositions));
            edges.push_back({N, N + 1});
        }
        polyscope::registerCurveNetwork("input curve", nodes, edges)->setColor({0, 0, 0})->setEnabled(true);
    }
    if (dualChain.size() > 0) {
        std::vector<Vector3> nodes;
        std::vector<std::array<size_t, 2>> edges;
        SurfaceMesh& mesh = geometry.mesh;
        polyscope::SurfaceMesh* psMesh = polyscope::getSurfaceMesh("input mesh");
        glm::vec3 baseColor = psMesh->getSurfaceColor();
        FaceData<glm::vec3> fColors(mesh, baseColor);
        for (const auto& pair : dualChain) {
            fColors[pair[0]] = {1., 0.388235, 0.278431}; // orange for "outside"
            fColors[pair[1]] = {0.294118, 0., 0.509804}; // indigo for "inside"
        }
        psMesh->addFaceColorQuantity("input dual chain", fColors)->setEnabled(true);
        polyscope::registerCurveNetwork("input curve edges", nodes, edges)->setColor({0, 0, 0})->setEnabled(true);
    }
    if (curveHalfedges.size() > 0) {
        std::vector<Vector3> nodes;
        std::vector<std::array<size_t, 2>> edges;
        for (const Halfedge& he : curveHalfedges) {
            size_t N = nodes.size();
            nodes.push_back(geometry.vertexPositions[he.tailVertex()]);
            nodes.push_back(geometry.vertexPositions[he.tipVertex()]);
            edges.push_back({N, N + 1});
        }
        polyscope::registerCurveNetwork("input curve", nodes, edges)->setColor({0, 0, 0})->setEnabled(true);
    }
}


CornerData<Vector2> interpolateVector2AcrossB(CommonSubdivision& cs, const CornerData<Vector2>& dataB) {

    CornerData<Vector2> interp(*cs.mesh);
    for (Vertex v : cs.mesh->vertices()) {
        SurfacePoint posB = cs.sourcePoints[v]->posB;
        if (posB.type == SurfacePointType::Face) {
            // this vertex in the common subdivison lies within a face in the intrinsic mesh
            SurfacePoint pB_face = posB.inSomeFace();
            Face fB = pB_face.face;
            Vector2 val = {0., 0.};
            size_t idx = 0;
            for (Halfedge he : fB.adjacentHalfedges()) {
                Corner cB = he.corner();
                val += dataB[cB] * pB_face.faceCoords[idx];
                idx++;
            }
            for (Corner c : v.adjacentCorners()) {
                interp[c] = val;
            }
        } else {
            for (Corner c : v.adjacentCorners()) {
                Face f = c.face();
                Face fB = cs.sourceFaceB[f];
                SurfacePoint pB_face = posB.inFace(fB);
                Vector2 val = {0., 0.};
                size_t idx = 0;
                for (Halfedge he : fB.adjacentHalfedges()) {
                    Corner cB = he.corner();
                    val += dataB[cB] * pB_face.faceCoords[idx];
                    idx++;
                }
                interp[c] = val;
            }
        }
    }
    return interp;
}