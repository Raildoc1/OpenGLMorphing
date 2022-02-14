#include "VBO.h"
#include "Mesh.h"
#include <map>
#include <set>

enum class VertexType { Source, Target, Extra, Merged };
enum class CoeffType { One, DHC, MVC, Kanai };

struct VertexData {
	Vertex vertex;
	int index;
	int eqClass;
	int setIndex;
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

	bool isTriangle(UniqueEdgeData& e1, UniqueEdgeData& e2) {
		if (v1 == e1.v1) {
			if (e2.v1 == e1.v2 && e2.v2 == v2) {
				return true;
			} else if (e2.v1 == v2 && e2.v2 == e1.v2) {
				return true;
			}
		} else if (v1 == e1.v2) {
			if (e2.v1 == e1.v1 && e2.v2 == v2) {
				return true;
			} else if (e2.v1 == v2 && e2.v2 == e1.v1) {
				return true;
			}
		}

		if (v1 == e2.v1) {
			if (e1.v1 == e2.v2 && e1.v2 == v2) {
				return true;
			} else if (e1.v1 == v2 && e1.v2 == e2.v2) {
				return true;
			}
		} else if (v1 == e2.v2) {
			if (e1.v1 == e2.v1 && e1.v2 == v2) {
				return true;
			} else if (e1.v1 == v2 && e1.v2 == e2.v1) {
				return true;
			}
		}

		return false;
	}
};

class MeshData {
private:
	const float EPSILON = 0.00001f;
	const int MAX_ITER = 25000;

	Mesh mesh;

	CoeffType coeffType = CoeffType::Kanai;

	int vertexCount;
	int edgesCount;
	int indicesCount;

	float borderLength = 0.0f;
	float lastEnergy = 0.0f;
	float rotation = 0.0f;

	std::vector<int> fixedIndices;
	std::map<int, glm::vec2> derivatives;
	std::vector<BorderVertex> borderVertices;
	std::set<int> eqClassesSet;

	float** A;
	float* U;
	float* V;
	int looseVerteicesAmount;

	bool initialized = false;
	bool invertBorder = false;

	void initEdges();
	void initVertices();
	void initFixedIndices();
	void initBorder();
	void sortBorder();
	void initUniqueEdges();
	void initHarmonicK();
	void initLambda();
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
	std::vector<int> vertexSetList;

	float** k;
	float** lambda;

	MeshData(Mesh& mesh, float rotation, bool invertBorder);
	~MeshData();
	void init();

	float tickMap();
	void harmonizeMap();
	glm::vec3 findVertexPos(glm::vec2 mapPos);
	glm::vec3 findBorderPos(float phi);
};