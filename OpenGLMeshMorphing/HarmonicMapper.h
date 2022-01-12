#pragma once
#include "MeshData.h"
#include <map>

struct MapEntity {
	glm::vec2 image;
	bool locked;
};

class HarmonicMapper
{
private:
	MeshData* source;
	MeshData* target;

	bool initialized = false;

	float calculateEnergy();

public:
	HarmonicMapper(MeshData &source, MeshData &target);

	std::map<int, MapEntity> sourceMap;
	std::map<int, MapEntity> targetMap;

	void init();

	void initMap(MeshData* mesh, std::map<int, MapEntity> &map);
};

