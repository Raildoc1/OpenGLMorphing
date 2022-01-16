#include "HarmonicMapper.h"

HarmonicMapper::HarmonicMapper(MeshData &source, MeshData &target)
{
	this->source = &source;
	this->target = &target;

	map = std::map<int, MapEntity>();
	uniqueEdges = std::vector<UniqueEdgeData>();
	border = std::vector<BorderEntity>();
}

void HarmonicMapper::init()
{
	initMap();
	initEdges();
	fixMapBound();
	fixIntersections();
	initialized = true;
}

void HarmonicMapper::initMap()
{
	for (auto const& x : source->map)
	{
		map[x.first] = x.second;
	}

	int vertexCount = source->getVertexCount();

	for (auto const& x : target->map)
	{
		map[x.first + vertexCount] = x.second;
	}
}

void HarmonicMapper::initEdges()
{
	for (auto const& x : source->uniqueEdges)
	{
		UniqueEdgeData temp = x;
		uniqueEdges.push_back(temp);
	}

	int vertexCount = source->getVertexCount();

	for (auto const& x : target->uniqueEdges)
	{
		UniqueEdgeData temp = x;
		temp.v1 += vertexCount;
		temp.v2 += vertexCount;
		uniqueEdges.push_back(temp);
	}
}

void HarmonicMapper::fixMapBound()
{

	for (auto const& x : uniqueEdges)
	{
		if (!x.isBorder)
		{
			continue;
		}

		BorderEntity e;
		e.eqClass = x.v1;
		e.phi = map[x.v1].phi;

		if (border.size() == 0) {
			border.push_back(e);
			continue;
		}

		if (e.phi <= border[0].phi) {
			auto it = border.begin();
			border.insert(it, e);
			continue;
		}

		if (e.phi >= border[border.size() - 1].phi) {
			border.push_back(e);
			continue;
		}

		for (int i = 0; i < border.size() - 1; i++)
		{
			if (border[i].phi < e.phi && e.phi < border[i + 1].phi) {
				border.insert(border.begin() + i + 1, e);
				break;
			}
		}
	}

	for (auto it = uniqueEdges.begin(); it != uniqueEdges.end();)
	{
		if (!(*it).isBorder)
		{
			it++;
			continue;
		}

		it = uniqueEdges.erase(it);
	}

	for (size_t i = 0; i < border.size(); i++)
	{
		std::cout << border[i].eqClass << " - " << border[i].phi << std::endl;

		UniqueEdgeData e;

		int j = (i + 1) % (border.size());

		e.v1 = border[i].eqClass;
		e.v2 = border[j].eqClass;

		e.isBorder = true;

		uniqueEdges.push_back(e);
	}
}

void HarmonicMapper::fixIntersections()
{

}

void HarmonicMapper::Equalize(int v1, int v2)
{
	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		if (uniqueEdges[i].v1 == v2) {
			uniqueEdges[i].v1 = v1;
		}
		if (uniqueEdges[i].v2 == v2) {
			uniqueEdges[i].v2 = v1;
		}
	}


	// TODO: uniq
	/*for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		for (size_t j = 0; j < uniqueEdges.size(); j++)
		{
			if (i == j) {
				continue;
			}


		}
	}*/
}