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

	bool isBorder;

	bool equals(EdgeData &e) const {
		if (v1.eqClass == e.v1.eqClass) {
			return v2.eqClass == e.v2.eqClass;
		} else if (v1.eqClass == e.v2.eqClass) {
			return v2.eqClass == e.v1.eqClass;
		}
		return false;
	}
};

struct UniqueEdgeData {
	int v1;
	int v2;

	bool isBorder;

	bool equals(UniqueEdgeData& e) {
		return (v1 == e.v1 && v2 == e.v2) 
			|| (v1 == e.v2 && v2 == e.v1);
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
	void initUniqueEdges();

public:
	int getVertexCount() { return vertexCount; }
	int getEdgesCount() { return edgesCount; }

	VertexData* vertices;
	EdgeData* edges;

	std::vector<EdgeData> border;
	std::vector<UniqueEdgeData> uniqueEdges;

	MeshData(Mesh& mesh);
	void init();
};