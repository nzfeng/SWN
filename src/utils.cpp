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


// ===================== I/O

std::string getHomeDirectory(const std::string& filepath) {

    std::string dir(filepath.begin(), filepath.begin() + filepath.find_last_of("/"));
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

// TODO: texture exporting

// TODO: nice curve exporting


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

    SurfacePoint pt0 = curveNodes[0];
    SurfaceMesh* mesh;
    if (pt0.type == SurfacePointType::Vertex) {
        mesh = pt0.vertex.getMesh();
    } else if (pt0.type == SurfacePointType::Edge) {
        mesh = pt0.edge.getMesh();
    } else if (pt0.type == SurfacePointType::Face) {
        mesh = pt0.face.getMesh();
    }

    if (mesh->isOriented()) {
        for (const auto& seg : curveEdges) {
            const SurfacePoint& pt1 = curveNodes[seg[0]];
            const SurfacePoint& pt2 = curveNodes[seg[1]];

            // // For some reason, very occasionally this code fails. I checked that pt2 is indeed among the adjacent
            // vertices of pt1; the == operator among vertices just doesn't work sometimes, I guess.

            // Edge e = sharedEdge(pt1, pt2);
            // if (e == Edge()) std::cerr << "no common edge found" << std::endl;
            // if (e == Edge()) return std::vector<Halfedge>();
            // if (pt1.type != SurfacePointType::Vertex || pt2.type != SurfacePointType::Vertex)
            //     return std::vector<Halfedge>();

            // Halfedge he = (pt1.vertex == e.firstVertex()) ? e.halfedge() : e.halfedge().twin();
            // curveHalfedges.push_back(he);

            if (pt1.type != SurfacePointType::Vertex || pt2.type != SurfacePointType::Vertex) {
                return std::vector<Halfedge>();
            }

            Halfedge he = Halfedge();
            for (Halfedge heA : pt1.vertex.outgoingHalfedges()) {
                if (heA.tipVertex() == pt2.vertex) {
                    he = heA;
                    break;
                }
            }

            if (he != Halfedge()) {
                curveHalfedges.push_back(he);
            } else {
                std::cerr << "Input curve segment between <" << pt1 << "> and <" << pt2
                          << "> is not constrained to mesh edges." << std::endl;
                return std::vector<Halfedge>();
            }
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

Vector<int> convertToChain(IntrinsicGeometryInterface& geom, const std::vector<Halfedge>& curve) {

    SurfaceMesh& mesh = geom.mesh;
    Vector<int> chain = Vector<int>::Zero(mesh.nEdges());
    geom.requireEdgeIndices();
    for (const Halfedge& he : curve) {
        size_t eIdx = geom.edgeIndices[he.edge()];
        chain[eIdx] += he.orientation() ? 1 : -1;
    }
    geom.unrequireEdgeIndices();
    return chain;
}

/*
 * Given a collection of halfedges, determine all the connected curve components it represents.
 * This function is useful for extracting smooth curves when we go to export curves as OBJ and render them. If the
 * curve is exported as a collection of (disconnected) edges, then the resulting curve when thickened and rendered
 * in Blender will have unsightly gaps.
 */
std::vector<std::vector<Halfedge>> getCurveComponents(IntrinsicGeometryInterface& geom,
                                                      const std::vector<Halfedge>& curveHalfedges) {

    SurfaceMesh& mesh = geom.mesh;
    Vector<int> chain = convertToChain(geom, curveHalfedges);
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
        for (const SurfacePoint& p : curveNodes) nodes.push_back(p.interpolate(geometry.vertexPositions));
        polyscope::registerCurveNetwork("input curve", nodes, curveEdges)->setColor({0, 0, 0})->setEnabled(true);
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