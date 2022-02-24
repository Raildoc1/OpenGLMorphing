#pragma once
#include "MeshData.h"
#include "SuperMesh.h"
#include <map>

struct BorderEntity
{
	int eqClass;
	float phi;
	VertexType type;
};

struct MorphEntity
{
	VertexType vertexType;
	int baseEqClass;
	glm::vec3 srcPos;
	glm::vec3 tarPos;
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

public:
	std::map<int, MapEntity> map;
	std::map<int, MorphEntity> finalMorphMap;
	std::vector<UniqueEdgeData> uniqueEdges;
	std::vector<BorderEntity> border;

	HarmonicMapper(MeshData &source, MeshData &target);

	static bool TryFindIntersection(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d, glm::vec2* intersection, bool exclusively);

	void init();
	void initMap();
	void initEdges();
	void fixMapBound();
	void fixIntersections();
	void mergeCloseVertices();
	bool fixIntersection(int i0, int i1, int j0, int j1, bool moveBound, int* bound1, int* bound2);
	void fixUniqueEdges();
	void retriangulate();
	void fast_retriangulate();
	void clearMap();

	void Equalize(int v1, int v2);
	bool hasEdge(int v1, int v2);

	SuperMesh* generateSuperMesh();
};

