#pragma once
#include "MeshData.h"
#include "SuperMesh.h"
#include <map>

using std::vector, std::map;

struct BorderEntity
{
	int eqClass;
	float phi;
	VertexType type;
};

struct MorphEntity
{
	VertexType vertexType;
	glm::vec3 srcPos;
	glm::vec3 tarPos;
	glm::vec3 srcNorm;
	glm::vec3 tarNorm;

	MorphEntity() :
		vertexType(VertexType::Unknown),
		srcPos(glm::vec3()),
		tarPos(glm::vec3()),
		srcNorm(glm::vec3()),
		tarNorm(glm::vec3())
	{ }

	MorphEntity(VertexType vertexType) :
		vertexType(vertexType),
		srcPos(glm::vec3()),
		tarPos(glm::vec3()),
		srcNorm(glm::vec3()),
		tarNorm(glm::vec3())
	{ }

	MorphEntity(VertexType vertexType, glm::vec3 srcPos, glm::vec3 tarPos, glm::vec3 srcNorm, glm::vec3 tarNorm) :
		vertexType(vertexType),
		srcPos(srcPos),
		tarPos(tarPos),
		srcNorm(srcNorm),
		tarNorm(tarNorm)
	{ }
};

struct IntersectionEntity
{
	int vertexIndex;
	glm::vec2 image;

	IntersectionEntity(int vertexIndex, glm::vec2 image) :
		vertexIndex(vertexIndex),
		image(image)
	{ }
};

struct ClockwiseComparator {

	std::map<int, MapEntity>* map;

	ClockwiseComparator(std::map<int, MapEntity>* map) : map(map) { }

	bool operator () (const UniqueEdgeData &e1, const UniqueEdgeData& e2)
	{
		if (e1.v1.eqClass != e2.v1.eqClass) {
			std::cout << "edges have different origins!" << std::endl;
			std::cout << "--- " << std::string(e1) << " " << std::string(e2) << std::endl;
		}

		glm::vec2 a = (*map)[e1.v2].image - (*map)[e1.v1].image;
		glm::vec2 b = (*map)[e2.v2].image - (*map)[e1.v1].image;

		float angleA = atan2(a.x, a.y);
		float angleB = atan2(b.x, b.y);

		return angleA > angleB;
	}
};

struct Feature {
	int srcEqClass;
	int tarEqClass;

	Feature(int srcEqClass, int tarEqClass) :
		srcEqClass(srcEqClass),
		tarEqClass(tarEqClass)
	{ }
};

class HarmonicMapper
{
private:
	MeshData* source;
	MeshData* target;

	bool initialized = false;

	int lastVertexIndex = 0;
	int firstExtraIndex = 0;
	int firstTargetIndex = 0;
	bool extraIndexesReached = false;

	float mergeDistance = 0.001f;
	int nextVertexIndex;

	std::vector <SuperVertex> superVertices = std::vector <SuperVertex>();
	std::vector <GLuint> superIndices = std::vector <GLuint>();
public:
	std::map<int, MapEntity> map;
	std::map<int, MorphEntity> finalMorphMap;
	vector<UniqueEdgeData> uniqueEdges;
	vector<BorderEntity> border;

	vector<vector<IntersectionEntity>> sourceIntersections;
	vector<vector<IntersectionEntity>> targetIntersections;

	HarmonicMapper(MeshData& source, MeshData& target);

	void init(vector<Feature>& features);
	void init();
	void initMap();
	void adjustFeatures(vector<Feature>& features);
	void mergeMaps();
	void initEdges();
	void fixMapBound();
	void fixed_fixIntersections();
	void mergeCloseVertices();
	bool fixIntersection(int i0, int i1, int j0, int j1, bool moveBound, int* bound1, int* bound2);
	void fixUniqueEdges();
	void retriangulate();
	void fast_retriangulate();
	void clearMap();

	void Equalize(int v1, int v2);
	bool hasEdge(int v1, int v2);
	bool isClockWise(glm::vec2 a, glm::vec2 b, glm::vec2 c);
	void markEdgeRemoved(int v);

	SuperMesh* generateSuperMesh();
};

