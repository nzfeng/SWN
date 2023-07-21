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

/*
 * Convert curve represented using barycentric points to a collection of halfedges, if all curve edges are constrained
 * to mesh edges.
 *
 * Return an empty vector if curve is not constrained to mesh edges.
 */
std::vector<Halfedge> setCurveHalfedges(const std::vector<SurfacePoint>& curveNodes,
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


/* Read in curve, encoded as OBJ line objects. */
void readLines(SurfaceMesh& mesh, const std::string& filepath, std::vector<SurfacePoint>& curveNodes,
               std::vector<std::array<size_t, 2>>& curveEdges, int offset) {

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
                std::vector<std::array<size_t, 2>>& curveEdges, int offset) {

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
            } else if (X == "#offset") {
                iss >> offset;
            }
        }
        curr_file.close();
    } else {
        std::cerr << "Could not open file <" << filepath << ">" << std::endl;
    }
}


// ===================== MESH MUTATION


// ===================== VISUALIZATION

void displayCurves(const VertexPositionGeometry& geometry, const std::vector<SurfacePoint>& curveNodes,
                   const std::vector<std::array<size_t, 2>>& curveEdges) {

    std::vector<Vector3> nodes;
    for (const SurfacePoint& p : curveNodes) nodes.push_back(p.interpolate(geometry.vertexPositions));
    polyscope::registerCurveNetwork("input curve", nodes, curveEdges)->setColor({0, 0, 0})->setEnabled(true);
}