#include "HarmonicMapper.h"

HarmonicMapper::HarmonicMapper(MeshData &source, MeshData &target)
{
	this->source = &source;
	this->target = &target;

	sourceMap = std::map<int, MapEntity>();
	targetMap = std::map<int, MapEntity>();
}

void HarmonicMapper::init()
{
	initMap(source, sourceMap);
	initMap(source, sourceMap);

	initialized = true;
}

void HarmonicMapper::initMap(MeshData* mesh, std::map<int, MapEntity> &map)
{
	map = std::map<int, MapEntity>();

	std::cout << "Border vertices:" << std::endl;

	float currentLength = 0.0f;
	float dPi = glm::pi<float>() * 2.0f;
	float borderLength = mesh->getBorderLength();

	for (size_t i = 0; i < mesh->border.size(); i++)
	{
		float t = dPi * currentLength / borderLength;

		MapEntity e;
		e.locked = true;
		e.image = glm::vec2(glm::cos(t), glm::sin(t));

		currentLength += mesh->border[i].length;

		map[mesh->border[i].v1.eqClass] = e;
	}

	std::cout << "Border length: " << mesh->getBorderLength() << std::endl;

	for (size_t i = 0; i < mesh->getVertexCount(); i++)
	{
		int eqClass = mesh->vertices[i].eqClass;
		bool duplicate = false;

		for (size_t i = 0; i < map.size(); i++)
		{
			if (map.count(eqClass) > 0) {
				duplicate = true;
				break;
			}
		}

		if (duplicate) {
			continue;
		}

		MapEntity e;
		e.image = glm::vec2(0.0f, 0.0f);
		e.locked = false;

		map[eqClass] = e;
	}

}
