#include "VBO.h"
#include "Mesh.h"
#include <map>
#include <set>

enum class VertexType { Unknown, Source, Target, Merged, Removed };
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

struct UniqueVertexData {
	VertexType type;
	int eqClass;
	bool isBorder;

	UniqueVertexData(VertexType type, int eqClass, bool isBorder) :
		type(type),
		eqClass(eqClass),
		isBorder(isBorder)
	{}

	UniqueVertexData(int eqClass, bool isBorder) :
		type(VertexType::Unknown),
		eqClass(eqClass),
		isBorder(isBorder)
	{}

	bool operator==(const UniqueVertexData& rhs) {
		return eqClass == rhs.eqClass && type == rhs.type;
	}

	operator int() const { return eqClass; }
};

struct UniqueEdgeData {
	UniqueVertexData v1;
	UniqueVertexData v2;

	VertexType type;

	UniqueEdgeData(UniqueVertexData v1, UniqueVertexData v2) :
		v1(v1),
		v2(v2)
	{}

	UniqueEdgeData(VertexType type, UniqueVertexData v1, UniqueVertexData v2) :
		type(type),
		v1(v1),
		v2(v2)
	{}
	
	bool isBorder() const {
		return v1.isBorder && v2.isBorder;
	}

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
		if (v1.eqClass == e1.v1.eqClass) {
			if (e2.v1.eqClass == e1.v2.eqClass && e2.v2.eqClass == v2.eqClass) {
				return true;
			} else if (e2.v1.eqClass == v2.eqClass && e2.v2.eqClass == e1.v2.eqClass) {
				return true;
			}
		} else if (v1.eqClass == e1.v2.eqClass) {
			if (e2.v1.eqClass == e1.v1.eqClass && e2.v2.eqClass == v2.eqClass) {
				return true;
			} else if (e2.v1.eqClass == v2.eqClass && e2.v2.eqClass == e1.v1.eqClass) {
				return true;
			}
		}

		if (v1.eqClass == e2.v1.eqClass) {
			if (e1.v1.eqClass == e2.v2.eqClass && e1.v2.eqClass == v2.eqClass) {
				return true;
			} else if (e1.v1.eqClass == v2.eqClass && e1.v2.eqClass == e2.v2.eqClass) {
				return true;
			}
		} else if (v1.eqClass == e2.v2.eqClass) {
			if (e1.v1.eqClass == e2.v1.eqClass && e1.v2.eqClass == v2.eqClass) {
				return true;
			} else if (e1.v1.eqClass == v2.eqClass && e1.v2.eqClass == e2.v1.eqClass) {
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

	CoeffType coeffType = CoeffType::MVC;

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

	bool isBorder(const UniqueEdgeData& e) { return vertices[e.v1.eqClass].isBorder && vertices[e.v2.eqClass].isBorder; };

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