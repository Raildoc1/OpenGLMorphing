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

	MorphEntity() :
		vertexType(VertexType::Unknown),
		srcPos(glm::vec3()),
		tarPos(glm::vec3())
	{ }

	MorphEntity(VertexType vertexType) :
		vertexType(vertexType),
		srcPos(glm::vec3()),
		tarPos(glm::vec3())
	{ }

	MorphEntity(VertexType vertexType, glm::vec3 srcPos, glm::vec3 tarPos) :
		vertexType(vertexType),
		srcPos(srcPos),
		tarPos(tarPos)
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

public:
	std::map<int, MapEntity> map;
	std::map<int, MorphEntity> finalMorphMap;
	vector<UniqueEdgeData> uniqueEdges;
	vector<BorderEntity> border;

	vector<vector<IntersectionEntity>> sourceIntersections;
	vector<vector<IntersectionEntity>> targetIntersections;

	HarmonicMapper(MeshData& source, MeshData& target);

	void init();
	void initMap();
	void mergeMaps();
	void initEdges();
	void fixMapBound();
	void fixed_fixIntersections();
	void mergeCloseVertices();
	bool fixIntersection(int i0, int i1, int j0, int j1, bool moveBound, int* bound1, int* bound2);
	void fixUniqueEdges();
	void fast_retriangulate();
	void clearMap();

	void Equalize(int v1, int v2);
	bool hasEdge(int v1, int v2);
	bool isClockWise(glm::vec2 a, glm::vec2 b, glm::vec2 c);
	void markEdgeRemoved(int v);

	SuperMesh* generateSuperMesh();
};

