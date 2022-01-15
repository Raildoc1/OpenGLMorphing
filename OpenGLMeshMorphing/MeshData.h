#include "VBO.h"
#include "Mesh.h"
#include <map>

struct VertexData {
	Vertex vertex;
	int index;
	int eqClass;
	bool isBorder;
};

struct MapEntity {
	glm::vec2 image;
	bool locked;
};

struct EdgeData {
	VertexData v1;
	VertexData v2;

	float length;

	bool isBorder;

	bool equals(EdgeData &e) const {
		if (v1.eqClass == e.v1.eqClass) {
			return v2.eqClass == e.v2.eqClass;
		} else if (v1.eqClass == e.v2.eqClass) {
			return v2.eqClass == e.v1.eqClass;
		}
		return false;
	}

	bool adjacent(EdgeData& e) {
		if (this->equals(e)) {
			return false;
		}

		if (v1.eqClass == e.v1.eqClass || v1.eqClass == e.v2.eqClass) {
			return true;
		}

		if (v2.eqClass == e.v1.eqClass || v2.eqClass == e.v2.eqClass) {
			return true;
		}

		return false;
	}
};

struct UniqueEdgeData {
	int v1;
	int v2;

	glm::vec3 pos1;
	glm::vec3 pos2;

	bool isBorder;

	bool equals(UniqueEdgeData& e) {
		return (v1 == e.v1 && v2 == e.v2) 
			|| (v1 == e.v2 && v2 == e.v1);
	}
};

class MeshData {
private:
	const float EPSILON = 0.001f;
	const int MAX_ITER = 10000;

	Mesh mesh;

	int vertexCount;
	int edgesCount;

	float borderLength = 0.0f;
	float lastEnergy = 0.0f;

	std::vector<int> fixedIndices;
	std::map<int, glm::vec2> derivatives;

	bool initialized = false;

	void initEdges();
	void initVertices();
	void initFixedIndices();
	void initBorder();
	void sortBorder();
	void initUniqueEdges();
	void initHarmonicK();
	void initMap();

	int edgeInTriangle(int e1, int e2, int t1, int t2, int t3);

	float calculateMapEnergy();

public:
	int getVertexCount() { return vertexCount; }
	int getEdgesCount() { return edgesCount; }
	float getBorderLength() { return borderLength; }

	VertexData* vertices;
	EdgeData* edges;

	std::vector<EdgeData> border;
	std::vector<UniqueEdgeData> uniqueEdges;
	std::map<int, MapEntity> map;

	float** k;

	MeshData(Mesh& mesh);
	~MeshData();
	void init();

	float tickMap();
	void harmonizeMap();
};