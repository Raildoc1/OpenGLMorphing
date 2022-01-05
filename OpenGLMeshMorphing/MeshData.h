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

public:
	int getVertexCount() { return vertexCount; }
	int getEdgesCount() { return edgesCount; }

	VertexData* vertices;
	EdgeData* edges;

	MeshData(Mesh& mesh);
	void init();
};