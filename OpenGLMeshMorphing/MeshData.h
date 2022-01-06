#include "VBO.h"
#include "Mesh.h"

struct VertexData {
	Vertex vertex;
	int index;
	int eqClass;
};

struct EdgeData {
	VertexData v1;
	VertexData v2;

	bool equals(EdgeData &e) const {
		if (v1.eqClass == e.v1.eqClass) {
			return v2.eqClass == e.v2.eqClass;
		} else if (v1.eqClass == e.v2.eqClass) {
			return v2.eqClass == e.v1.eqClass;
		}
		return false;
	}
};

class MeshData {
private:
	const float EPSILON = 0.001f;

	Mesh mesh;

	int vertexCount;
	int edgesCount;

	bool initialized = false;

	void initEdges();
	void initVertices();
	void initBorder();

public:
	int getVertexCount() { return vertexCount; }
	int getEdgesCount() { return edgesCount; }

	VertexData* vertices;
	EdgeData* edges;

	std::vector<EdgeData> border;

	MeshData(Mesh& mesh);
	void init();
};