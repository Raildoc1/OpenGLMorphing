#pragma once
#include "MeshData.h"
#include <map>

class HarmonicMapper
{
private:
	MeshData* source;
	MeshData* target;

	bool initialized = false;

	float calculateEnergy();

public:
	HarmonicMapper(MeshData &source, MeshData &target);

	void init();
};

