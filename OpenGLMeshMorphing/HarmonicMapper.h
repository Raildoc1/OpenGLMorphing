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

public:
	std::map<int, MapEntity> map;
	std::vector<UniqueEdgeData> uniqueEdges;
	std::vector<BorderEntity> border;

	HarmonicMapper(MeshData &source, MeshData &target);

	void init();
	void initMap();
	void initEdges();
	void fixMapBound();
	void fixIntersections();

	void Equalize(int v1, int v2);
};

