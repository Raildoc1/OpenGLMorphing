#include "HarmonicMapper.h"
#include <glm/gtx/string_cast.hpp>

HarmonicMapper::HarmonicMapper(MeshData &source, MeshData &target)
{
	this->source = &source;
	this->target = &target;

	map = std::map<int, MapEntity>();
	uniqueEdges = std::vector<UniqueEdgeData>();
	border = std::vector<BorderEntity>();

	lastVertexIndex = this->source->getVertexCount() + this->target->getVertexCount() - 1;
	firstExtraIndex = lastVertexIndex + 1;
}

bool HarmonicMapper::TryFindIntersection(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d, glm::vec2* intersection)
{
	float v1 = (d.x - c.x) * (a.y - c.y) - (d.y - c.y) * (a.x - c.x);
	float v2 = (d.x - c.x) * (b.y - c.y) - (d.y - c.y) * (b.x - c.x);
	float v3 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	float v4 = (b.x - a.x) * (d.y - a.y) - (b.y - a.y) * (d.x - a.x);

	float z1 = (b.y - a.y) / (b.x - a.x);
	float z2 = (d.y - c.y) / (d.x - c.x);

	*intersection = glm::vec2();

	intersection->x = (c.y - a.y - z2 * c.x + z1 * a.x) / (z1 - z2);
	intersection->y = (b.y - a.y) * (intersection->x - a.x) / (b.x - a.x) + a.y;

	return (v1 * v2 < 0.0f) && (v3 * v4 < 0.0f);
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
	fixUniqueEdges();
	while (fixIntersection()) { }
}

bool HarmonicMapper::fixIntersection()
{
	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		if (uniqueEdges[i].v1 == -1) {
			continue;
		}

		bool intersectionFound = false;

		if (i > firstExtraIndex) {
			extraIndexesReached = true;
		}

		size_t j = i + 1;

		if (extraIndexesReached && j < firstExtraIndex) {
			j = firstExtraIndex;
		}

		for (; j < uniqueEdges.size(); j++)
		{
			if (i == j) {
				continue;
			}

			if (uniqueEdges[j].v1 == -1) {
				continue;
			}

			int a = uniqueEdges[i].v1;
			int b = uniqueEdges[i].v2;
			int c = uniqueEdges[j].v1;
			int d = uniqueEdges[j].v2;

			glm::vec2 intersection;

			if (intersectionFound = TryFindIntersection(map[a].image, map[b].image, map[c].image, map[d].image, &intersection)) {

				//std::cout << "intersect: (" << to_string(map[a].image) << ", " << to_string(map[b].image) << ") - (" << to_string(map[c].image) << ", " << to_string(map[d].image) << ")" << std::endl;
				//std::cout << "intersect: (" << a << ", " << b << ") - (" << c << ", " << d << ")" << std::endl;

				MapEntity e;
				e.border = false;
				e.locked = false;
				e.phi = 0.0f;
				e.image = intersection;

				int centerIndex = ++lastVertexIndex;

				map[centerIndex] = e;

				UniqueEdgeData e1;
				UniqueEdgeData e2;
				UniqueEdgeData e3;
				UniqueEdgeData e4;

				glm::vec3 center1 = (uniqueEdges[i].pos1 + uniqueEdges[i].pos2) / 2.0f;
				glm::vec3 center2 = (uniqueEdges[j].pos1 + uniqueEdges[j].pos2) / 2.0f;

				e1.v1 = a;
				e1.v2 = centerIndex;
				e1.isBorder = uniqueEdges[i].isBorder;
				e1.pos1 = uniqueEdges[i].pos1;
				e1.pos2 = center1;

				e2.v1 = centerIndex;
				e2.v2 = b;
				e2.isBorder = uniqueEdges[i].isBorder;
				e2.pos1 = center1;
				e2.pos2 = uniqueEdges[i].pos2;

				e3.v1 = c;
				e3.v2 = centerIndex;
				e3.isBorder = uniqueEdges[j].isBorder;
				e3.pos1 = uniqueEdges[j].pos1;
				e3.pos2 = center2;

				e4.v1 = centerIndex;
				e4.v2 = d;
				e4.isBorder = uniqueEdges[j].isBorder;
				e4.pos1 = center2;
				e1.pos2 = uniqueEdges[j].pos2;

				uniqueEdges[i].v1 = -1;
				uniqueEdges[j].v1 = -1;

				uniqueEdges.push_back(e1);
				uniqueEdges.push_back(e2);
				uniqueEdges.push_back(e3);
				uniqueEdges.push_back(e4);

				fixUniqueEdges();
				return true;
			}
		}
	}

	return false;
}

void HarmonicMapper::fixUniqueEdges()
{
	for (size_t i = 0; i < uniqueEdges.size() - 1; i++)
	{
		for (size_t j = i + 1; j < uniqueEdges.size(); j++)
		{
			if (uniqueEdges[i].equals(uniqueEdges[j])) {
				uniqueEdges[j].v1 = -1;
			}
		}
	}

	for (auto it = uniqueEdges.begin(); it != uniqueEdges.end();)
	{
		if (it->v1 == -1)
		{
			it = uniqueEdges.erase(it);
			continue;
		}

		it++;
	}
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

	fixUniqueEdges();
}