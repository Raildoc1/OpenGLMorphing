#include "HarmonicMapper.h"

HarmonicMapper::HarmonicMapper(MeshData &source, MeshData &target)
{
	this->source = &source;
	this->target = &target;

	sourceMap = std::vector<MapEntity>();
	targetMap = std::vector<MapEntity>();
}

void HarmonicMapper::init()
{
	initMap(source, sourceMap);
	initMap(source, sourceMap);

	initialized = true;
}

void HarmonicMapper::initMap(MeshData* mesh, std::vector<MapEntity>& map)
{
	for (size_t i = 0; i < mesh->getVertexCount(); i++)
	{
		if (mesh->vertices[i].isBorder) {
			std::cout << i << " is border vertex" << std::endl;
		}
	}
}
