#include "VBO.h"
#include "Mesh.h"
#include <map>
#include <set>

using std::vector, std::map, std::set;

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

struct TriangleVertex {
	int index;
	glm::vec2 image;
	bool isBorder;

	TriangleVertex(int index, glm::vec2 image, bool isBorder) :
		index(index),
		image(image),
		isBorder(isBorder)
	{ }

	TriangleVertex(VertexData vertexData, glm::vec2 image) :
		index(vertexData.eqClass),
		image(image),
		isBorder(vertexData.isBorder)
	{ }
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

	int triangle1 = -1;
	int triangle2 = -1;

	UniqueEdgeData() :
		v1(UniqueVertexData(-1, false)),
		v2(UniqueVertexData(-1, false)),
		type(VertexType::Unknown)
	{}

	UniqueEdgeData(UniqueVertexData v1, UniqueVertexData v2) :
		v1(v1),
		v2(v2),
		type(VertexType::Unknown)
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

struct Triangle {
	TriangleVertex a;
	TriangleVertex b;
	TriangleVertex c;

	vector<UniqueEdgeData> edges;

	Triangle(TriangleVertex a, TriangleVertex b, TriangleVertex c) : a(a), b(b), c(c)
	{
		invAr = 1.0f / (-b.image.y * c.image.x + a.image.y * (-b.image.x + c.image.x) + a.image.x * (b.image.y - c.image.y) + b.image.x * c.image.y);

		s1 = a.image.y * c.image.x - a.image.x * c.image.y;
		s2 = c.image.y - a.image.y;
		s3 = a.image.x - c.image.x;

		t1 = a.image.x * b.image.y - a.image.y * b.image.x;
		t2 = a.image.y - b.image.y;
		t3 = b.image.x - a.image.x;

		UniqueVertexData v1(a.index, a.isBorder);
		UniqueVertexData v2(b.index, b.isBorder);
		UniqueVertexData v3(c.index, c.isBorder);

		edges = vector<UniqueEdgeData>();
		edges.push_back(UniqueEdgeData(v1, v2));
		edges.push_back(UniqueEdgeData(v2, v3));
		edges.push_back(UniqueEdgeData(v3, v1));
	}

	void getBarycentricCoordinates(glm::vec2 point, float* s, float* t) {
		*s = invAr * (s1 + s2 * point.x + s3 * point.y);
		*t = invAr * (t1 + t2 * point.x + t3 * point.y);
	}

	bool equals(Triangle& e) const {
		if (a.index == e.a.index) {
			if (b.index == e.b.index) {
				return c.index == e.c.index;
			}
			if (b.index == e.c.index) {
				return c.index == e.b.index;
			}
		}
		if (a.index == e.b.index) {
			if (b.index == e.a.index) {
				return c.index == e.c.index;
			}
			if (b.index == e.c.index) {
				return c.index == e.a.index;
			}
		}
		if (a.index == e.c.index) {
			if (b.index == e.a.index) {
				return c.index == e.b.index;
			}
			if (b.index == e.b.index) {
				return c.index == e.a.index;
			}
		}
		return false;
	}

private:
	float invAr;
	float s1, s2, s3;
	float t1, t2, t3;
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

	vector<int> fixedIndices;
	map<int, glm::vec2> derivatives;
	vector<BorderVertex> borderVertices;
	set<int> eqClassesSet;

	int looseVerteicesAmount;

	bool initialized = false;
	bool invertBorder = false;

	void initEdges();
	void initVertices();
	void initBorder();
	void sortBorder();
	void initUniqueEdges();
	void initHarmonicK();
	void initLambda();
	void initMap();
	void initTriangles();
	void initFixedIndices();

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

	vector<EdgeData> border;
	vector<UniqueEdgeData> uniqueEdges;
	map<int, MapEntity> unitCircleMap;
	vector<int> vertexSetList;
	vector<Triangle> triangles;
	vector<vector<int>> meshMatrix;

	vector<vector<vector<int>>> edgesToTriangles;

	float** k;
	float** lambda;

	MeshData(Mesh& mesh, float rotation, bool invertBorder);
	~MeshData();
	void init();

	float tickMap();
	void harmonizeMap();
	void addTriangleToEdgeMap(int v1, int v2, int triangle);
	int getOppositeTriangle(int v1, int v2, int triangle);
	glm::vec3 findVertexPos(glm::vec2 mapPos, int* triangle);
	glm::vec3 findBorderPos(float phi);
	int getBorderTriangle(float phi);
};