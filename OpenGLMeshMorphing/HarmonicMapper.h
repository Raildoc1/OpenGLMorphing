#pragma once
#include "MeshData.h"
#include <map>

struct BorderEntity
{
	int eqClass;
	float phi;
};

class HarmonicMapper
{
private:
	MeshData* source;
	MeshData* target;

	bool initialized = false;

	int lastVertexIndex = 0;
	int firstExtraIndex = 0;
	bool extraIndexesReached = false;

public:
	std::map<int, MapEntity> map;
	std::vector<UniqueEdgeData> uniqueEdges;
	std::vector<BorderEntity> border;

	HarmonicMapper(MeshData &source, MeshData &target);

	static bool TryFindIntersection(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d, glm::vec2* intersection);

	void init();
	void initMap();
	void initEdges();
	void fixMapBound();
	void fixIntersections();
	bool fixIntersection();
	void fixUniqueEdges();

	void Equalize(int v1, int v2);
};

