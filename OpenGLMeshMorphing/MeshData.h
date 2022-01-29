#include "VBO.h"
#include "Mesh.h"
#include <map>

enum class VertexType { Source, Target, Extra, Merged };

struct VertexData {
	Vertex vertex;
	int index;
	int eqClass;
	bool isBorder;
};

struct MapEntity {
	glm::vec2 image;
	bool locked;
	bool border;
	float phi;
};

struct BorderVertex {
	int eqClass;
	float phi;
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

	VertexType type;

	int baseEqClass1;
	int baseEqClass2;

	bool isBorder;

	bool equals(UniqueEdgeData& e) {
		return (v1 == e.v1 && v2 == e.v2) 
			|| (v1 == e.v2 && v2 == e.v1);
	}

	bool adjacent(UniqueEdgeData& e) {
		if (this->equals(e)) {
			return false;
		}

		if (v1 == e.v1 || v1 == e.v2) {
			return true;
		}

		if (v2 == e.v1 || v2 == e.v2) {
			return true;
		}

		return false;
	}
};

class MeshData {
private:
	const float EPSILON = 0.001f;
	const int MAX_ITER = 10000;

	Mesh mesh;

	int vertexCount;
	int edgesCount;
	int indicesCount;

	float borderLength = 0.0f;
	float lastEnergy = 0.0f;

	std::vector<int> fixedIndices;
	std::map<int, glm::vec2> derivatives;
	std::vector<BorderVertex> borderVertices;

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
	int getIndicesCount() { return indicesCount; }

	VertexData* vertices;
	EdgeData* edges;
	int* triangles;

	std::vector<EdgeData> border;
	std::vector<UniqueEdgeData> uniqueEdges;
	std::map<int, MapEntity> map;

	float** k;

	MeshData(Mesh& mesh);
	~MeshData();
	void init();

	float tickMap();
	void harmonizeMap();
	glm::vec3 findVertexPos(glm::vec2 mapPos);
	glm::vec3 findBorderPos(float phi);
};