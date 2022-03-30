#include "HarmonicMapper.h"
#include <glm/gtx/string_cast.hpp>
#include <filesystem>

void print_time(const clock_t& prev, const std::string msg) {
	std::cout << std::setw(30) << msg << std::setw(30) << float(clock() - prev) / CLOCKS_PER_SEC << " seconds." << std::endl;
}

HarmonicMapper::HarmonicMapper(MeshData& source, MeshData& target)
{
	this->source = &source;
	this->target = &target;

	map = std::map<int, MapEntity>();
	finalMorphMap = std::map<int, MorphEntity>();
	uniqueEdges = std::vector<UniqueEdgeData>();
	border = std::vector<BorderEntity>();

	int source_vertex_count = source.getVertexCount();
	int unique_edges_count = source.uniqueEdges.size();

	source_edges = (std::vector<UniqueEdgeData>*)calloc(source_vertex_count, sizeof(std::vector<UniqueEdgeData>));

	for (size_t i = 0; i < source_vertex_count; i++)
	{
		source_edges[i] = std::vector<UniqueEdgeData>();
	}

	for (size_t i = 0; i < unique_edges_count; i++)
	{
		int v1 = source.uniqueEdges[i].v1;
		int v2 = source.uniqueEdges[i].v2;

		source_edges[v1].push_back(source.uniqueEdges[i]);
		source_edges[v2].push_back(source.uniqueEdges[i]);
	}

	used_edges = (bool**)calloc(source_vertex_count, sizeof(bool*));
	for (size_t i = 0; i < source_vertex_count; i++)
	{
		used_edges[i] = (bool*)calloc(source_vertex_count, sizeof(bool));
	}

	lastVertexIndex = this->source->getVertexCount() + this->target->getVertexCount() - 1;
	nextVertexIndex = lastVertexIndex + 1;

	sourceIntersections = std::vector<std::vector<IntersectionEntity>>(source.uniqueEdges.size(), std::vector<IntersectionEntity>());
	targetIntersections = std::vector<std::vector<IntersectionEntity>>(target.uniqueEdges.size(), std::vector<IntersectionEntity>());
}

bool HarmonicMapper::TryFindIntersection(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d, glm::vec2* intersection, bool exclusively)
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

	if ((v1 * v2 >= 0.0f) || (v3 * v4 >= 0.0f)) {
		//std::cout << v1 * v2 << ", " << v3 * v4 << std::endl;
	}

	if (exclusively) {
		return (v1 * v2 < 0.0f) && (v3 * v4 < 0.0f);
	}

	return ((v1 * v2 < 0.0f) && (v3 * v4 <= 0.0f)) || ((v1 * v2 <= 0.0f) && (v3 * v4 < 0.0f));
}

void HarmonicMapper::init()
{
	clock_t time_stamp = clock();

	initMap();
	print_time(time_stamp, "initMap finished");
	time_stamp = clock();

	initEdges();
	print_time(time_stamp, "initEdges finished");
	time_stamp = clock();

	mergeMaps();
	print_time(time_stamp, "mergeMaps finished");
	time_stamp = clock();

	//fixIntersections();
	fixed_fixIntersections();
	print_time(time_stamp, "fixIntersections finished");
	time_stamp = clock();

	fixMapBound();
	print_time(time_stamp, "fixMapBound finished");
	time_stamp = clock();

	clearMap();
	print_time(time_stamp, "clearMap finished");
	time_stamp = clock();

	//fast_retriangulate();
	print_time(time_stamp, "retriangulate finished");
	time_stamp = clock();

	initialized = true;
}

void HarmonicMapper::initMap()
{
	//float rotation = glm::pi<float>();
	float rotation = 0.0f;

	for (auto const& x : source->unitCircleMap)
	{
		int i = x.first;

		map[i] = x.second;

		glm::vec3 rotatedPosition = glm::rotateZ(glm::vec3(x.second.image.x, x.second.image.y, 0.0f), -rotation);
		int triangle;

		glm::vec3 tarPos;

		if (!x.second.border) {
			tarPos = target->findVertexPos(x.second.image, &triangle);
		}
		else {
			tarPos = target->findBorderPos(x.second.phi);
		}

		finalMorphMap[i] = MorphEntity(
			VertexType::Source,
			source->vertices[i].vertex.position,
			tarPos
		);
	}

	int vertexCount = source->getVertexCount();

	for (auto const& x : target->unitCircleMap)
	{
		int i = x.first;

		map[i + vertexCount] = x.second;

		glm::vec3 rotatedPosition = glm::rotateZ(glm::vec3(x.second.image.x, x.second.image.y, 0.0f), rotation);
		int triangle;

		glm::vec3 srcPos;

		if (!x.second.border) {
			srcPos = source->findVertexPos(x.second.image, &triangle);
		}
		else {
			srcPos = source->findBorderPos(x.second.phi);
		}

		finalMorphMap[i + vertexCount] = MorphEntity(
			VertexType::Target,
			srcPos,
			target->vertices[i].vertex.position
		);
	}
}

void HarmonicMapper::initEdges()
{
	for (auto const& x : source->uniqueEdges)
	{
		UniqueEdgeData temp = x;
		temp.v1.type = temp.v2.type = temp.type = VertexType::Source;
		temp.type = VertexType::Source;
		uniqueEdges.push_back(temp);
	}

	int vertexCount = source->getVertexCount();
	firstTargetIndex = source->uniqueEdges.size();

	for (auto const& x : target->uniqueEdges)
	{
		UniqueEdgeData temp = x;
		temp.v1.eqClass += vertexCount;
		temp.v2.eqClass += vertexCount;
		temp.v1.type = temp.v2.type = temp.type = VertexType::Target;
		temp.type = VertexType::Target;
		uniqueEdges.push_back(temp);
	}

	firstExtraIndex = source->uniqueEdges.size() + target->uniqueEdges.size();
}

void HarmonicMapper::fixMapBound()
{
	for (auto const& x : uniqueEdges)
	{
		if (!x.isBorder())
		{
			continue;
		}

		BorderEntity e;
		e.eqClass = x.v1.eqClass;
		e.phi = map[x.v1.eqClass].phi;
		e.type = finalMorphMap[x.v1.eqClass].vertexType;

		if (border.size() == 0) {
			border.push_back(e);
			continue;
		}

		if (e.phi < border[0].phi) {
			auto it = border.begin();
			border.insert(it, e);
			continue;
		}
		else if (e.phi == border[0].phi) {
			Equalize(border[0].eqClass, e.eqClass);
			continue;
		}

		if (e.phi > border[border.size() - 1].phi) {
			border.push_back(e);
			continue;
		}
		else if (e.phi == border[border.size() - 1].phi) {
			border[border.size() - 1].type = VertexType::Merged;
			Equalize(border[border.size() - 1].eqClass, e.eqClass);
			continue;
		}

		for (int i = 0; i < border.size() - 1; i++)
		{
			if (border[i].phi < e.phi && e.phi < border[i + 1].phi) {
				border.insert(border.begin() + i + 1, e);
				break;
			}
			else if (border[i].phi == e.phi) {
				border[i].type = VertexType::Merged;
				Equalize(border[i].eqClass, e.eqClass);
			}
		}
	}

	fixUniqueEdges();

	for (auto it = uniqueEdges.begin(); it != uniqueEdges.end();)
	{
		if (!(*it).isBorder())
		{
			it++;
			continue;
		}

		it = uniqueEdges.erase(it);
	}

	for (size_t i = 0; i < border.size(); i++)
	{
		//std::cout << border[i].eqClass << " - " << border[i].phi << std::endl;


		int j = (i + 1) % (border.size());

		UniqueVertexData v1 = UniqueVertexData(finalMorphMap[border[i].eqClass].vertexType, border[i].eqClass, true);
		UniqueVertexData v2 = UniqueVertexData(finalMorphMap[border[j].eqClass].vertexType, border[j].eqClass, true);
		UniqueEdgeData e(v1, v2);

		uniqueEdges.push_back(e);
	}

	//std::cout << "Border length: " << border.size() << std::endl;
}

void HarmonicMapper::fixIntersections()
{
	fixUniqueEdges();
	std::cout << "looking for intersections..." << std::endl;
	std::cout << "Edges amount = " << uniqueEdges.size() << std::endl;

	std::cout << "Checking (0, " << firstTargetIndex << ") -> (" << firstTargetIndex << ", " << firstExtraIndex << ")" << std::endl;

	int n = 0;

	while (fixIntersection(0, firstTargetIndex, firstTargetIndex, firstExtraIndex - 1, false, nullptr, nullptr)) {
		n++;
	}

	std::cout << "first step finished with " << uniqueEdges.size() << " edges!" << std::endl;
	fixUniqueEdges();

	int fixAmount = 0;
	int bound1 = 0;
	int bound2 = firstExtraIndex;
	while (fixIntersection(bound1, uniqueEdges.size(), bound2, uniqueEdges.size(), true, &bound1, &bound2)) {
		n++;
		if (++fixAmount % 100 == 0) {
			//std::cout << fixAmount << std::endl;
		}
	}

	fixUniqueEdges();
	std::cout << n << " intersections fixed!" << std::endl;
	std::cout << "Edges amount = " << uniqueEdges.size() << std::endl;
}

void HarmonicMapper::mergeCloseVertices() {
	for (size_t i = 0; i < uniqueEdges.size() - 1; i++)
	{
		if (uniqueEdges[i].isBorder()) {
			continue;
		}

		for (size_t j = i + 1; j < uniqueEdges.size(); j++)
		{
			if (uniqueEdges[j].isBorder()) {
				continue;
			}

			int v1 = uniqueEdges[i].v1.eqClass;
			int u1 = uniqueEdges[i].v2.eqClass;
			int v2 = uniqueEdges[j].v1.eqClass;
			int u2 = uniqueEdges[j].v2.eqClass;

			bool iIsBorder = uniqueEdges[i].isBorder();
			bool jIsBorder = uniqueEdges[j].isBorder();

			glm::vec2 v1_map = map[v1].image;
			glm::vec2 u1_map = map[u1].image;
			glm::vec2 v2_map = map[v2].image;
			glm::vec2 u2_map = map[u2].image;

			if (glm::distance(v1_map, v2_map) < mergeDistance) {
				if (iIsBorder) {
					map[v2].image = map[v1].image;
				}
				else if (jIsBorder) {
					map[v1].image = map[v2].image;
				}
				else {
					glm::vec2 avgPos = (map[v1].image + map[v2].image) / 2.0f;
					map[v1].image = map[v2].image = avgPos;
				}
				Equalize(v1, v2);
			}
			else if (glm::distance(v1_map, u2_map) < mergeDistance) {
				if (iIsBorder) {
					map[u2].image = map[v1].image;
				}
				else if (jIsBorder) {
					map[v1].image = map[u2].image;
				}
				else {
					glm::vec2 avgPos = (map[v1].image + map[u2].image) / 2.0f;
					map[v1].image = map[u2].image = avgPos;
				}
				Equalize(v1, u2);
			}

			if (glm::distance(u1_map, v2_map) < mergeDistance) {
				if (iIsBorder) {
					map[u1].image = map[v2].image;
				}
				else if (jIsBorder) {
					map[v2].image = map[u1].image;
				}
				else {
					glm::vec2 avgPos = (map[u1].image + map[v2].image) / 2.0f;
					map[u1].image = map[v2].image = avgPos;
				}
				Equalize(u1, v2);
			}
			else if (glm::distance(u1_map, u2_map) < mergeDistance) {
				if (iIsBorder) {
					map[u1].image = map[u2].image;
				}
				else if (jIsBorder) {
					map[u2].image = map[u1].image;
				}
				else {
					glm::vec2 avgPos = (map[u1].image + map[u2].image) / 2.0f;
					map[u1].image = map[u2].image = avgPos;
				}
				Equalize(u1, u2);
			}
		}
	}

	fixUniqueEdges();
}

bool HarmonicMapper::fixIntersection(int i0, int i1, int j0, int j1, bool moveBound = false, int* bound1 = nullptr, int* bound2 = nullptr)
{
	for (size_t i = i0; i < i1; i++)
	{
		if (uniqueEdges[i].v1.type == VertexType::Removed) {
			//std::cout << "this edge must be deleted" << std::endl;
			continue;
		}


		for (size_t j = j0; j < j1; j++)
		{
			if (i == j) {
				continue;
			}

			if (uniqueEdges[j].v1.type == VertexType::Removed) {
				//std::cout << "this edge must be deleted" << std::endl;
				continue;
			}

			auto a = uniqueEdges[i].v1;
			auto b = uniqueEdges[i].v2;
			auto c = uniqueEdges[j].v1;
			auto d = uniqueEdges[j].v2;

			glm::vec2 intersection;

			if (TryFindIntersection(map[a].image, map[b].image, map[c].image, map[d].image, &intersection, true)) {

				//std::cout << "intersect: (" << to_string(map[a].image) << ", " << to_string(map[b].image) << ") - (" << to_string(map[c].image) << ", " << to_string(map[d].image) << ")" << std::endl;
				//std::cout << "intersect: (" << a << ", " << b << ") - (" << c << ", " << d << ")" << std::endl;

				MapEntity e;
				e.border = false;
				e.locked = false;
				e.phi = 0.0f;
				e.image = intersection;

				int centerIndex = ++lastVertexIndex;

				map[centerIndex] = e;

				glm::vec3 aPos;
				glm::vec3 bPos;
				glm::vec3 cPos;
				glm::vec3 dPos;

				float t1 = glm::distance(intersection, map[a].image) / glm::distance(map[b].image, map[a].image);
				float t2 = glm::distance(intersection, map[c].image) / glm::distance(map[d].image, map[c].image);

				MorphEntity me;
				me.vertexType = VertexType::Merged;

				VertexType abType = VertexType::Unknown;
				VertexType cdType = VertexType::Unknown;

				if (uniqueEdges[i].v1.type == VertexType::Source
					|| uniqueEdges[i].v2.type == VertexType::Source
					|| uniqueEdges[i].type == VertexType::Source
					|| uniqueEdges[j].v1.type == VertexType::Target
					|| uniqueEdges[j].v2.type == VertexType::Target
					|| uniqueEdges[j].type == VertexType::Target) {
					aPos = finalMorphMap[a].srcPos;
					bPos = finalMorphMap[b].srcPos;
					cPos = finalMorphMap[c].tarPos;
					dPos = finalMorphMap[d].tarPos;

					me.srcPos = aPos * (1.0f - t1) + bPos * t1;
					me.tarPos = cPos * (1.0f - t2) + dPos * t2;

					abType = VertexType::Source;
					cdType = VertexType::Target;

				}
				else if (uniqueEdges[i].v1.type == VertexType::Target
					|| uniqueEdges[i].v2.type == VertexType::Target
					|| uniqueEdges[i].type == VertexType::Target
					|| uniqueEdges[j].v1.type == VertexType::Source
					|| uniqueEdges[j].v2.type == VertexType::Source
					|| uniqueEdges[j].type == VertexType::Source) {
					aPos = finalMorphMap[a].tarPos;
					bPos = finalMorphMap[b].tarPos;
					cPos = finalMorphMap[c].srcPos;
					dPos = finalMorphMap[d].srcPos;

					me.tarPos = aPos * (1.0f - t1) + bPos * t1;
					me.srcPos = cPos * (1.0f - t2) + dPos * t2;

					abType = VertexType::Target;
					cdType = VertexType::Source;
				}
				else {
					std::cout << "********************* " << i << ": type1 = " << (int)uniqueEdges[i].v1.type << "; type2 = " << (int)uniqueEdges[i].v2.type << std::endl;
				}

				finalMorphMap[centerIndex] = me;

				auto center = UniqueVertexData(VertexType::Merged, centerIndex, false);

				UniqueEdgeData e1 = UniqueEdgeData(abType, a, center);
				UniqueEdgeData e2 = UniqueEdgeData(abType, center, b);
				UniqueEdgeData e3 = UniqueEdgeData(cdType, c, center);
				UniqueEdgeData e4 = UniqueEdgeData(cdType, center, d);

				markEdgeRemoved(i);
				markEdgeRemoved(j);

				int delta = 0;

				if (i < firstExtraIndex) {
					delta++;
				}

				if (j < firstExtraIndex) {
					delta++;
				}

				firstExtraIndex -= delta;

				uniqueEdges.push_back(e1);
				uniqueEdges.push_back(e2);
				uniqueEdges.push_back(e3);
				uniqueEdges.push_back(e4);
				return true;
			}
		}

		if (moveBound) {
			(*bound1)++;
		}
	}

	if (moveBound) {
		(*bound2)++;
	}

	return false;
}

void HarmonicMapper::fixUniqueEdges()
{
	std::cout << "Looking for doubles: " << uniqueEdges.size() << " -> ";

	for (size_t i = 0; i < uniqueEdges.size() - 1; i++)
	{
		if (uniqueEdges[i].v1 == uniqueEdges[i].v2) {
			markEdgeRemoved(i);
			continue;
		}

		for (size_t j = i + 1; j < uniqueEdges.size(); j++)
		{
			if (uniqueEdges[i].equals(uniqueEdges[j])) {
				markEdgeRemoved(j);
				//std::cout << "remove " << j << "; firstTargetIndex = " << firstTargetIndex << "; firstExtraIndex = " << firstExtraIndex << std::endl;
			}
		}
	}

	int i = 0;

	for (auto it = uniqueEdges.begin(); it != uniqueEdges.end();)
	{
		if (i < firstTargetIndex) {
			firstTargetIndex--;
		}

		if (i < firstExtraIndex) {
			firstExtraIndex--;
		}

		if (it->v1.type == VertexType::Removed)
		{
			it = uniqueEdges.erase(it);
			continue;
		}

		it++;
		i++;
	}

	std::cout << uniqueEdges.size() << std::endl;
}

void HarmonicMapper::mergeMaps()
{
	int n = 0;
	UniqueVertexData v1 = source->uniqueEdges[0].v1;
	std::vector<UniqueEdgeData> work_list = std::vector<UniqueEdgeData>();
	std::vector<UniqueEdgeData> candidate_list = std::vector<UniqueEdgeData>();

	for (auto& e : source_edges[v1])
	{
		work_list.push_back(e);
		used_edges[e.v1][e.v2] = used_edges[e.v2][e.v1] = true;
	}

	while (!work_list.empty()) {
		UniqueEdgeData e_a = work_list.back();
		work_list.pop_back();

		UniqueVertexData v1 = e_a.v1;
		UniqueVertexData v2 = e_a.v2;

		if (e_a.isBorder()) {
			continue;
		}

		int triangleIndex = -1;

		if (v2.isBorder) {
			UniqueVertexData temp = v1;
			v1 = v2;
			v2 = temp;
		}

		if (v1.isBorder) {
			triangleIndex = target->getBorderTriangle(source->unitCircleMap[v1.eqClass].phi);
		}
		else {
			target->findVertexPos(source->unitCircleMap[v1.eqClass].image, &triangleIndex);
		}

		if (triangleIndex < 0) {
			std::cout << "no triangle found for vertex { index = " << v1.eqClass << ", isBorder = " << v1.isBorder << "}" << std::endl;
		}

		if (triangleIndex >= 0) {
			candidate_list.insert(
				std::end(candidate_list),
				std::begin(target->triangles[triangleIndex].edges),
				std::end(target->triangles[triangleIndex].edges)
			);
		}

		while (!candidate_list.empty()) {
			UniqueEdgeData e_b = candidate_list.back();
			candidate_list.pop_back();

			if (e_b.isBorder()) {
				continue;
			}

			glm::vec2 intersection;

			if (TryFindIntersection(
				source->unitCircleMap[v1].image, source->unitCircleMap[v2].image,
				target->unitCircleMap[e_b.v1].image, target->unitCircleMap[e_b.v2].image,
				&intersection, true
			)) {
				int intersectionIndex = nextVertexIndex++;
				IntersectionEntity in = IntersectionEntity(intersectionIndex, intersection);

				sourceIntersections[source->meshMatrix[e_a.v1][e_a.v2]].push_back(in);
				targetIntersections[target->meshMatrix[e_b.v1][e_b.v2]].push_back(in);

				finalMorphMap[intersectionIndex] = MorphEntity(VertexType::Merged); // TODO: final morph map!

				map[intersectionIndex] = MapEntity(intersection);

				n++;
				int oppositeTriangle = target->getOppositeTriangle(e_b.v1, e_b.v2, triangleIndex);
				triangleIndex = oppositeTriangle;

				if (oppositeTriangle == -1) {
					continue;
				}

				for (UniqueEdgeData& e : target->triangles[oppositeTriangle].edges)
				{
					if (e.equals(e_b)) {
						continue;
					}
					candidate_list.push_back(e);
				}
			}
		}

		for (auto& e : source_edges[v2])
		{
			if (used_edges[e.v1][e.v2]) {
				continue;
			}

			work_list.push_back(e);
			used_edges[e.v1][e.v2] = used_edges[e.v2][e.v1] = true;
		}
	}
	std::cout << "intersections amount = " << n << std::endl;
}

void HarmonicMapper::retriangulate()
{
	std::cout << "retriangulating..." << std::endl;

	int i = 0;

	for (auto const& x : map)
	{
		for (auto const& y : map)
		{
			if (x.first == y.first) {
				continue;
			}

			if (hasEdge(x.first, y.first)) {
				continue;
			}

			bool intersectionFound = false;

			for (auto const& e : uniqueEdges)
			{
				glm::vec2 intersection;
				if (TryFindIntersection(x.second.image, y.second.image, map[e.v1].image, map[e.v2].image, &intersection, false)) {
					intersectionFound = true;
					break;
				}
			}

			if (!intersectionFound) {
				UniqueVertexData v1data = UniqueVertexData(finalMorphMap[x.first].vertexType, x.first, false);
				UniqueVertexData v2data = UniqueVertexData(finalMorphMap[y.first].vertexType, y.first, false);
				UniqueEdgeData e = UniqueEdgeData(v1data, v2data);

				uniqueEdges.push_back(e);

				//std::cout << i++ << std::endl;
			}
		}
	}

	std::cout << "retriangulating finished successfully!" << std::endl;
}

void HarmonicMapper::fast_retriangulate()
{
	std::cout << "retriangulating..." << std::endl;

	int n = 0;

	for (size_t i = 0; i < uniqueEdges.size() - 2; i++)
	{
		int triangles_desired_amount = uniqueEdges[i].isBorder() ? 1 : 2;

		for (size_t j = i + 1; j < uniqueEdges.size() - 1; j++)
		{
			if (!uniqueEdges[i].adjacent(uniqueEdges[j])) {
				continue;
			}

			int v1 = -1;
			int v2 = -1;

			if (uniqueEdges[i].v1 == uniqueEdges[j].v1) {
				v1 = uniqueEdges[i].v2;
				v2 = uniqueEdges[j].v2;
			}
			else if (uniqueEdges[i].v1 == uniqueEdges[j].v2) {
				v1 = uniqueEdges[i].v2;
				v2 = uniqueEdges[j].v1;
			}
			else if (uniqueEdges[i].v2 == uniqueEdges[j].v1) {
				v1 = uniqueEdges[i].v1;
				v2 = uniqueEdges[j].v2;
			}
			else {
				v1 = uniqueEdges[i].v1;
				v2 = uniqueEdges[j].v1;
			}

			if (v1 == -1 || v2 == -1) {
				std::cout << "[ERROR]: edge not found: v1 = " << v1 << ", v2 = " << v2 << std::endl;
				continue;
			}

			if (hasEdge(v1, v2)) {
				triangles_desired_amount--;
				if (triangles_desired_amount <= 0) {
					break;
				}
				continue;
			}

			bool intersectionFound = false;

			for (auto const& e : uniqueEdges)
			{
				auto temp = e;
				if (uniqueEdges[i].equals(temp) || uniqueEdges[j].equals(temp)) {
					continue;
				}

				glm::vec2 intersection;
				if (TryFindIntersection(map[v1].image, map[v2].image, map[e.v1].image, map[e.v2].image, &intersection, false)) {
					//std::cout << "***out***" << std::endl;
					intersectionFound = true;
					break;
				}
			}

			if (!intersectionFound) {
				UniqueVertexData v1data = UniqueVertexData(finalMorphMap[v1].vertexType, v1, false);
				UniqueVertexData v2data = UniqueVertexData(finalMorphMap[v2].vertexType, v2, false);
				UniqueEdgeData e = UniqueEdgeData(v1data, v2data);

				uniqueEdges.push_back(e);

				//std::cout << n++ << std::endl;

				triangles_desired_amount--;
				if (triangles_desired_amount <= 0) {
					break;
				}
			}

		}
	}

	std::cout << "retriangulating finished successfully!" << std::endl;
}

void HarmonicMapper::clearMap()
{
	//std::cout << "map size before = " << map.size() << std::endl;
	for (auto it = map.begin(); it != map.end();)
	{
		bool found = false;
		for (auto const& e : uniqueEdges)
		{
			if (it->first == e.v1 || it->first == e.v2) {
				found = true;
				break;
			}
		}

		if (!found)
		{
			it = map.erase(it);
			continue;
		}

		//std::cout << it->first << " -> " << to_string(it->second.image) << std::endl;

		it++;
	}
	//std::cout << "map size after = " << map.size() << std::endl;
}

void HarmonicMapper::fixed_fixIntersections()
{
	struct LineSortComparator {
		LineSortComparator(glm::vec2 base) : base(base) { }
		bool operator () (IntersectionEntity a, IntersectionEntity b) {
			return glm::distance2(base, a.image) < glm::distance2(base, b.image);
		}
		glm::vec2 base;
	};

	for (size_t i = 0; i < sourceIntersections.size(); i++)
	{
		if (sourceIntersections[i].empty()) {
			continue;
		}

		vector<IntersectionEntity> intersections = sourceIntersections[i];
		UniqueEdgeData e = source->uniqueEdges[i];
		std::sort(intersections.begin(), intersections.end(), LineSortComparator(source->unitCircleMap[e.v1].image));

		uniqueEdges[i].type = VertexType::Removed;

		//std::cout << "[" << i << "]" << (std::string)(source->uniqueEdges[i]) << " " << (std::string)(uniqueEdges[i]) << std::endl;
		//std::cout << (std::string)(source->uniqueEdges[0]) << " " << (std::string)(uniqueEdges[0]) << std::endl;

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Source,
				UniqueVertexData(VertexType::Source, e.v1.eqClass, e.v1.isBorder),
				UniqueVertexData(VertexType::Merged, intersections[0].vertexIndex, false)
			)
		);

		for (size_t i = 0; i < intersections.size() - 1; i++)
		{
			uniqueEdges.push_back(
				UniqueEdgeData(
					VertexType::Source,
					UniqueVertexData(VertexType::Merged, intersections[i].vertexIndex, false),
					UniqueVertexData(VertexType::Merged, intersections[i + 1].vertexIndex, false)
				)
			);
		}

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Source,
				UniqueVertexData(VertexType::Merged, intersections[intersections.size() - 1].vertexIndex, false),
				UniqueVertexData(VertexType::Source, e.v2.eqClass, e.v2.isBorder)
			)
		);
	}

	// TODO: target intersections!

	int vertexCount = source->getVertexCount();

	for (size_t i = 0; i < targetIntersections.size(); i++)
	{
		if (targetIntersections[i].empty()) {
			continue;
		}

		vector<IntersectionEntity> intersections = targetIntersections[i];
		UniqueEdgeData e = target->uniqueEdges[i];
		std::sort(intersections.begin(), intersections.end(), LineSortComparator(target->unitCircleMap[e.v1].image));

		uniqueEdges[i + source->uniqueEdges.size()].type = VertexType::Removed;

		//std::cout << (std::string)(target->uniqueEdges[i]) << " " << (std::string)(uniqueEdges[i + vertexCount]) << std::endl;

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Target,
				UniqueVertexData(VertexType::Target, e.v1.eqClass + vertexCount, e.v1.isBorder),
				UniqueVertexData(VertexType::Merged, intersections[0].vertexIndex, false)
			)
		);

		for (size_t i = 0; i < intersections.size() - 1; i++)
		{
			uniqueEdges.push_back(
				UniqueEdgeData(
					VertexType::Target,
					UniqueVertexData(VertexType::Merged, intersections[i].vertexIndex, false),
					UniqueVertexData(VertexType::Merged, intersections[i + 1].vertexIndex, false)
				)
			);
		}

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Target,
				UniqueVertexData(VertexType::Merged, intersections[intersections.size() - 1].vertexIndex, false),
				UniqueVertexData(VertexType::Target, e.v2.eqClass + vertexCount, e.v2.isBorder)
			)
		);
	}

	for (auto it = uniqueEdges.begin(); it != uniqueEdges.end();)
	{
		if (it->type == VertexType::Removed)
		{
			it = uniqueEdges.erase(it);
			continue;
		}

		it++;
	}
}

void HarmonicMapper::Equalize(int v1, int v2)
{
	auto type1 = finalMorphMap[v1].vertexType;
	auto type2 = finalMorphMap[v2].vertexType;

	VertexType mergedType = VertexType::Merged;

	if (type1 == type2) {
		mergedType = type1;
	}

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		if (uniqueEdges[i].v1 == v2) {
			uniqueEdges[i].v1.eqClass = v1;
			uniqueEdges[i].v1.type = mergedType;
		}
		if (uniqueEdges[i].v2 == v2) {
			uniqueEdges[i].v2.eqClass = v1;
			uniqueEdges[i].v2.type = mergedType;
		}
		//uniqueEdges[i].type = VertexType::Merged;
	}

	//fixUniqueEdges();
}

bool HarmonicMapper::hasEdge(int v1, int v2) {

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		if (uniqueEdges[i].v1 == v1 && uniqueEdges[i].v2 == v2) {
			return true;
		}
		if (uniqueEdges[i].v1 == v2 && uniqueEdges[i].v2 == v1) {
			return true;
		}
	}
	return false;
}

SuperMesh* HarmonicMapper::generateSuperMesh() {

	std::cout << "generating super mesh..." << std::endl;
	std::cout << "uniqueEdges.size() = " << uniqueEdges.size() << std::endl;
	fixUniqueEdges();
	std::cout << "uniqueEdges.size() = " << uniqueEdges.size() << std::endl;
	std::vector <SuperVertex>* vertices = new std::vector <SuperVertex>();
	std::vector <GLuint>* indices = new std::vector <GLuint>();
	std::cout << "vectors created successfully!" << std::endl;

	int n = 0;

	for (size_t i = 0; i < uniqueEdges.size() - 2; i++)
	{
		auto e1 = uniqueEdges[i];

		int triangles_desired_amount = e1.isBorder() ? 1 : 2;

		for (size_t j = i + 1; j < uniqueEdges.size() - 1; j++)
		{
			if (i == j) {
				continue;
			}

			auto e2 = uniqueEdges[j];
			for (size_t k = j + 1; k < uniqueEdges.size(); k++)
			{
				if (k == i || k == j) {
					continue;
				}

				auto e3 = uniqueEdges[k];

				if (e1.isTriangle(e2, e3)) {
					int v3_ind = 0;
					if (e3.v1.eqClass != e1.v1.eqClass && e3.v1.eqClass != e1.v2.eqClass) {
						v3_ind = e3.v1;
					}
					else {
						v3_ind = e3.v2;
					}

					SuperVertex v1;
					v1.position1 = finalMorphMap[e1.v1].srcPos;
					v1.position2 = finalMorphMap[e1.v1].tarPos;

					SuperVertex v2;
					v2.position1 = finalMorphMap[e1.v2].srcPos;
					v2.position2 = finalMorphMap[e1.v2].tarPos;

					SuperVertex v3;
					v3.position1 = finalMorphMap[v3_ind].srcPos;
					v3.position2 = finalMorphMap[v3_ind].tarPos;

					vertices->push_back(v1);
					indices->push_back(vertices->size() - 1);

					if (isClockWise(map[e1.v1].image, map[e1.v2].image, map[v3_ind].image)) {
						vertices->push_back(v2);
						indices->push_back(vertices->size() - 1);
						vertices->push_back(v3);
						indices->push_back(vertices->size() - 1);
					}
					else {
						vertices->push_back(v3);
						indices->push_back(vertices->size() - 1);
						vertices->push_back(v2);
						indices->push_back(vertices->size() - 1);
					}

					if (n++ % 100 == 0) {
						//std::cout << n << std::endl;
					}

					if ((--triangles_desired_amount) <= 0) {
						goto out;
					}

					break;
				}
			}
		}
	out:;
	}

	std::cout << "Super mesh vertices amount = " << vertices->size() << std::endl;
	/*for (auto const& v : *vertices)
	{
		std::cout << to_string(v.position1) << std::endl;
	}*/

	std::cout << "Super mesh indices amount = " << indices->size() << std::endl;

	SuperMesh* superMesh = new SuperMesh(*vertices, *indices);

	std::cout << "Super mesh generated successfully!" << std::endl;

	return superMesh;
}

bool HarmonicMapper::isClockWise(glm::vec2 a, glm::vec2 b, glm::vec2 c) {
	return ((b.x * c.y) - (b.y * c.x)) - ((a.x * c.y) - (a.y * c.x)) + ((a.x * b.y) - (a.y * b.x)) < 0.0f;
}

void HarmonicMapper::markEdgeRemoved(int v)
{
	uniqueEdges[v].v1.type = VertexType::Removed;
	uniqueEdges[v].v2.type = VertexType::Removed;
	uniqueEdges[v].type = VertexType::Removed;
}