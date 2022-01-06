#pragma once
#include "MeshData.h"

struct MapEntity {
	int eqClass;
	glm::vec2 image;
	bool isBorder;
};

class HarmonicMapper
{
private:
	MeshData* source;
	MeshData* target;

	bool initialized = false;

public:
	HarmonicMapper(MeshData &source, MeshData &target);

	std::vector<MapEntity> sourceMap;
	std::vector<MapEntity> targetMap;

	void init();

	void initMap(MeshData* mesh, std::vector<MapEntity> &map);
};

