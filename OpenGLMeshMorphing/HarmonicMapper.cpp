#include "HarmonicMapper.h"
#include <glm/gtx/string_cast.hpp>
#include <filesystem>
#include "Utils.h"

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
	
	lastVertexIndex = this->source->getVertexCount() + this->target->getVertexCount() - 1;
	nextVertexIndex = lastVertexIndex + 1;

	sourceIntersections = std::vector<std::vector<IntersectionEntity>>(source.uniqueEdges.size(), std::vector<IntersectionEntity>());
	targetIntersections = std::vector<std::vector<IntersectionEntity>>(target.uniqueEdges.size(), std::vector<IntersectionEntity>());
}

void HarmonicMapper::init(vector<Feature>& features)
{
	adjustFeatures(features);
	init();
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

	//clearMap();
	print_time(time_stamp, "clearMap finished");
	time_stamp = clock();

	//fast_retriangulate();
	retriangulate();
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
		glm::vec3 tarNorm;

		if (!x.second.border) {
			tarPos = target->findVertexPos(x.second.image, triangle, tarNorm);
		}
		else {
			tarPos = target->findBorderPos(x.second.phi, tarNorm);
		}

		finalMorphMap[i] = MorphEntity(
			VertexType::Source,
			source->vertices[i].vertex.position,
			tarPos,
			source->vertices[i].vertex.normal,
			tarNorm
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
		glm::vec3 srcNorm;

		if (!x.second.border) {
			srcPos = source->findVertexPos(x.second.image, triangle, srcNorm);
		}
		else {
			srcPos = source->findBorderPos(x.second.phi, srcNorm);
		}

		finalMorphMap[i + vertexCount] = MorphEntity(
			VertexType::Target,
			srcPos,
			target->vertices[i].vertex.position,
			srcNorm,
			target->vertices[i].vertex.normal
		);
	}
}

void HarmonicMapper::adjustFeatures(vector<Feature>& features)
{
	vector<glm::vec2> w = vector<glm::vec2>(features.size());
	vector<float> sourceD = vector<float>(features.size(), 0.5f);
	vector<float> targetD = vector<float>(features.size(), 0.5f);
	vector<glm::vec2> sourceT = vector<glm::vec2>(features.size());
	vector<glm::vec2> targetT = vector<glm::vec2>(features.size());
	vector<glm::vec2> sourceFeaturePositions = vector<glm::vec2>(features.size());
	vector<glm::vec2> targetFeaturePositions = vector<glm::vec2>(features.size());

	for (size_t i = 0; i < features.size(); i++)
	{
		sourceFeaturePositions[i] = source->unitCircleMap[features[i].srcEqClass].image;
		targetFeaturePositions[i] = target->unitCircleMap[features[i].tarEqClass].image;
		w[i] = 0.5f * (sourceFeaturePositions[i] + targetFeaturePositions[i]);
		sourceT[i] = w[i] - sourceFeaturePositions[i];
		targetT[i] = w[i] - targetFeaturePositions[i];
	}

	std::map<int, MapEntity> srcMap = source->unitCircleMap;
	std::map<int, MapEntity> tarMap = target->unitCircleMap;

	vector<int> changedVertices = vector<int>();
	bool changed = false;
	float d = 1.0f;

	while (d > 0.0f) {
		do {
			changed = false;
			for (size_t i = 0; i < features.size(); i++)
			{
				sourceFeaturePositions[i] = source->unitCircleMap[features[i].srcEqClass].image;
				sourceT[i] = w[i] - sourceFeaturePositions[i];
			}

			for (size_t i = 0; i < features.size(); i++)
			{
				if (glm::length(sourceT[i]) < 0.001f) {
					continue;
				}

				float t = 1.0f;

				do {
					changedVertices.clear();
					srcMap = source->unitCircleMap;

					std::cout << to_string(sourceFeaturePositions[i]) << " ->  " << to_string(w[i]) << " ";

					for (auto& m : srcMap) {

						if (m.second.border) {
							continue;
						}

						float featureDistance = glm::distance(m.second.image, sourceFeaturePositions[i]);

						if (featureDistance > d) {
							continue;
						}

						m.second.image += t * sourceT[i] * (d - featureDistance);
						changedVertices.push_back(m.first);
					}
					t *= 0.9f;
					std::cout << "map valid = " << source->checkMapValidity(srcMap, changedVertices) << std::endl;
				} while ((changedVertices.size() > 0 && !source->checkMapValidity(srcMap, changedVertices)) || t > 0.01f);

				source->unitCircleMap = srcMap;
				sourceFeaturePositions[i] = source->unitCircleMap[features[i].srcEqClass].image;
				glm::vec2 newT = w[i] - sourceFeaturePositions[i];

				if (glm::distance(newT, sourceT[i]) > 0.0001f) {
					changed = true;
				}
			}
		} while (changed);
		d -= 0.11f;
	}

	d = 1.0f;

	while (d > 0.0f) {
		do {
			changed = false;
			for (size_t i = 0; i < features.size(); i++)
			{
				targetFeaturePositions[i] = target->unitCircleMap[features[i].tarEqClass].image;
				targetT[i] = w[i] - targetFeaturePositions[i];
			}

			for (size_t i = 0; i < features.size(); i++)
			{
				if (glm::length(targetT[i]) < 0.001f) {
					continue;
				}

				std::cout << glm::length(targetT[i]) << std::endl;

				float t = 1.0f;

				do {
					changedVertices.clear();
					tarMap = target->unitCircleMap;

					std::cout << to_string(targetFeaturePositions[i]) << " ->  " << to_string(w[i]) << " ";

					for (auto& m : tarMap) {
						if (m.second.border) {
							continue;
						}

						float featureDistance = glm::distance(m.second.image, targetFeaturePositions[i]);

						if (featureDistance > d) {
							continue;
						}

						m.second.image += t * targetT[i] * (d - featureDistance);
						changedVertices.push_back(m.first);
					}
					t *= 0.9f;
					std::cout << "map valid = " << target->checkMapValidity(tarMap, changedVertices) << std::endl;
				} while (changedVertices.size() > 0 && !target->checkMapValidity(tarMap, changedVertices) && t > 0.01f);

				target->unitCircleMap = tarMap;
				targetFeaturePositions[i] = target->unitCircleMap[features[i].tarEqClass].image;
				glm::vec2 newT = w[i] - targetFeaturePositions[i];

				if (glm::distance(newT, targetT[i]) > 0.0001f) {
					changed = true;
				}
			}
		} while (changed);
		d -= 0.11f;
	}

	source->initTriangles();
	target->initTriangles();

	std::cout << "features adjusted successfully!" << std::endl;
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
		if (!x.isBorder)
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
		if (!(*it).isBorder)
		{
			it++;
			continue;
		}

		it = uniqueEdges.erase(it);
	}

	for (size_t i = 0; i < border.size(); i++)
	{
		int j = (i + 1) % (border.size());

		UniqueVertexData v1 = UniqueVertexData(finalMorphMap[border[i].eqClass].vertexType, border[i].eqClass, true);
		UniqueVertexData v2 = UniqueVertexData(finalMorphMap[border[j].eqClass].vertexType, border[j].eqClass, true);
		UniqueEdgeData e(v1, v2, true);

		uniqueEdges.push_back(e);
	}
}

void HarmonicMapper::mergeCloseVertices() {
	for (size_t i = 0; i < uniqueEdges.size() - 1; i++)
	{
		if (uniqueEdges[i].isBorder) {
			continue;
		}

		for (size_t j = i + 1; j < uniqueEdges.size(); j++)
		{
			if (uniqueEdges[j].isBorder) {
				continue;
			}

			int v1 = uniqueEdges[i].v1.eqClass;
			int u1 = uniqueEdges[i].v2.eqClass;
			int v2 = uniqueEdges[j].v1.eqClass;
			int u2 = uniqueEdges[j].v2.eqClass;

			bool iIsBorder = uniqueEdges[i].isBorder;
			bool jIsBorder = uniqueEdges[j].isBorder;

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
	int sourceVertexCount = source->getVertexCount();
	UniqueVertexData v1 = source->uniqueEdges[0].v1;

	vector<UniqueEdgeData> work_list = vector<UniqueEdgeData>();
	vector<UniqueEdgeData> candidate_list = vector<UniqueEdgeData>();
	vector<vector<bool>> used_edges = vector<vector<bool>>(sourceVertexCount, vector<bool>(sourceVertexCount, false));

	for (auto& e : source->adjacencyList[v1])
	{
		work_list.push_back(e);
		used_edges[e.v1][e.v2] = used_edges[e.v2][e.v1] = true;
	}

	while (!work_list.empty()) {
		UniqueEdgeData e_a = work_list.back();
		work_list.pop_back();

		UniqueVertexData v1 = e_a.v1;
		UniqueVertexData v2 = e_a.v2;

		if (e_a.isBorder) {
			continue;
		}

		int triangleIndex = -1;
		bool swapped = false;

		if (v2.isBorder) {
			UniqueVertexData temp = v1;
			v1 = v2;
			v2 = temp;
			swapped = true;
		}

		glm::vec3 normal;

		if (v1.isBorder) {
			triangleIndex = target->getBorderTriangle(source->unitCircleMap[v1.eqClass].phi);
		}
		else {
			target->findVertexPos(source->unitCircleMap[v1.eqClass].image, triangleIndex, normal);
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

			if (e_b.isBorder) {
				continue;
			}

			glm::vec2 intersection;

			if (Utils::tryFindIntersection(
				source->unitCircleMap[v1].image, source->unitCircleMap[v2].image,
				target->unitCircleMap[e_b.v1].image, target->unitCircleMap[e_b.v2].image,
				&intersection, true
			)) {
				int intersectionIndex = nextVertexIndex++;
				IntersectionEntity in = IntersectionEntity(intersectionIndex, intersection);

				sourceIntersections[source->meshMatrix[e_a.v1][e_a.v2]].push_back(in);
				targetIntersections[target->meshMatrix[e_b.v1][e_b.v2]].push_back(in);

				float t_src = glm::distance(intersection, map[e_a.v1].image) / glm::distance(map[e_a.v1].image, map[e_a.v2].image);
				float t_tar = glm::distance(intersection, map[e_b.v1 + sourceVertexCount].image) / glm::distance(map[e_b.v1 + sourceVertexCount].image, map[e_b.v2 + sourceVertexCount].image);

				finalMorphMap[intersectionIndex] = MorphEntity(
					VertexType::Merged,
					finalMorphMap[e_a.v1].srcPos * (1 - t_src) + finalMorphMap[e_a.v2].srcPos * t_src,
					finalMorphMap[e_b.v1 + sourceVertexCount].tarPos * (1 - t_tar) + finalMorphMap[e_b.v2 + sourceVertexCount].tarPos * t_tar,
					glm::normalize(finalMorphMap[e_a.v1].srcNorm * (1 - t_src) + finalMorphMap[e_a.v2].srcNorm * t_src),
					glm::normalize(finalMorphMap[e_b.v1 + sourceVertexCount].tarNorm * (1 - t_tar) + finalMorphMap[e_b.v2 + sourceVertexCount].tarNorm * t_tar)
				);
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

		int nextOrigin = swapped ? v1.eqClass : v2.eqClass;
		for (auto& e : source->adjacencyList[nextOrigin])
		{
			if (used_edges[e.v1][e.v2]) {
				continue;
			}

			work_list.push_back(e);
			used_edges[e.v1][e.v2] = used_edges[e.v2][e.v1] = true;
		}
	}
	std::cout << "intersections amount = " << n << std::endl;

	for (size_t i = 0; i < source->uniqueEdges.size(); i++)
	{
		if (!used_edges[source->uniqueEdges[i].v1][source->uniqueEdges[i].v2]) {
			std::cout << "edge " << std::string(source->uniqueEdges[i])
				<< " never used ( isBorder = " << source->uniqueEdges[i].isBorder
				<< " v1.isBorder = " << source->uniqueEdges[i].v1.isBorder
				<< " v2.isBorder = " << source->uniqueEdges[i].v2.isBorder << " )" << std::endl;
		}
	}
}

void HarmonicMapper::retriangulate()
{
	int superVertexCount = nextVertexIndex;
	vector<vector<int>> edgesUsed = vector<vector<int>>(superVertexCount, vector<int>(superVertexCount, 0));
	vector<bool> verticesUsed = vector<bool>(superVertexCount, false);
	vector<vector<UniqueEdgeData>> edges = vector<vector<UniqueEdgeData>>(superVertexCount, vector<UniqueEdgeData>());
	ClockwiseComparator comparator = ClockwiseComparator(&map);
	superVertices.clear();
	superIndices.clear();

	for (auto& e : uniqueEdges)
	{
		edges[e.v1].push_back(e);
		edges[e.v2].push_back(e.turn());
	}

	for (auto& e : edges)
	{
		std::sort(e.begin(), e.end(), comparator);
	}

	int v0 = uniqueEdges[0].v1.eqClass;
	vector<int> workList = vector<int>();
	vector<UniqueEdgeData> polygonList = vector<UniqueEdgeData>();
	workList.push_back(v0);
	verticesUsed[v0] = true;

	while (!workList.empty())
	{
		int origin = workList.back();
		workList.pop_back();

		for (auto& e : edges[origin])
		{
			int currentVertex = e.v2;
			int prevVertex = e.v1;
			polygonList.push_back(e);
			//std::cout << origin << " : ";
			while (currentVertex != origin) {
				//std::cout << currentVertex << " ";

				if (edgesUsed[polygonList.back().v1][polygonList.back().v2] > 1
					|| edgesUsed[polygonList.back().v2][polygonList.back().v1] > 1
					|| polygonList.size() > 6) {
					//std::cout << currentVertex << " abort ";
					polygonList.clear();
					break;
				}

				UniqueEdgeData currentEdge;

				for (size_t i = 0; i < edges[currentVertex].size(); i++)
				{
					if (edges[currentVertex][i].v2 == prevVertex) {
						prevVertex = currentVertex;
						currentEdge = edges[currentVertex][(i + 1) % edges[currentVertex].size()];
						polygonList.push_back(currentEdge);
						break;
					}
				}

				currentVertex = currentEdge.v2;
			}

			/*for (auto& p : polygonList)
			{
				std::cout << (std::string)(p) << " ";
			}*/

			for (auto& p : polygonList)
			{
				if (edgesUsed[p.v1][p.v2] > 1 || edgesUsed[p.v2][p.v1] > 1) {
					polygonList.clear();
					break;
				}
			}

			for (auto& p : polygonList)
			{
				edgesUsed[p.v1][p.v2]++;
				edgesUsed[p.v2][p.v1]++;
			}

			SuperVertex originSuperVertex(finalMorphMap[origin].srcPos, finalMorphMap[origin].tarPos, finalMorphMap[origin].srcNorm, finalMorphMap[origin].tarNorm);

			int polygonEdgesAmount = polygonList.size();
			for (int i = 1; i < polygonEdgesAmount - 1; i++)
			{
				SuperVertex v1(
					finalMorphMap[polygonList[i].v1].srcPos,
					finalMorphMap[polygonList[i].v1].tarPos,
					finalMorphMap[polygonList[i].v1].srcNorm,
					finalMorphMap[polygonList[i].v1].tarNorm
				);
				SuperVertex v2(
					finalMorphMap[polygonList[i + 1].v1].srcPos,
					finalMorphMap[polygonList[i + 1].v1].tarPos,
					finalMorphMap[polygonList[i + 1].v1].srcNorm,
					finalMorphMap[polygonList[i + 1].v1].tarNorm
				);

				superVertices.push_back(originSuperVertex);
				superIndices.push_back(superVertices.size() - 1);
				superVertices.push_back(v1);
				superIndices.push_back(superVertices.size() - 1);
				superVertices.push_back(v2);
				superIndices.push_back(superVertices.size() - 1);
			}

			//std::cout << polygonList.size() << std::endl;
			polygonList.clear();

			if (!verticesUsed[e.v2]) {
				workList.push_back(e.v2);
				verticesUsed[e.v2] = true;
			}
		}
	}

	for (auto& e : edges)
	{
		if (e.empty()) {
			continue;
		}
		if (!verticesUsed[e[0].v1]) {
			std::cout << "vertex" << e[0].v1 << " was never used!" << std::endl;
		}
	}
}

void HarmonicMapper::fast_retriangulate()
{
	std::cout << "retriangulating..." << std::endl;

	int n = 0;

	for (size_t i = 0; i < uniqueEdges.size() - 2; i++)
	{
		int triangles_desired_amount = uniqueEdges[i].isBorder ? 1 : 2;

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
				if (Utils::tryFindIntersection(map[v1].image, map[v2].image, map[e.v1].image, map[e.v2].image, &intersection, false)) {
					intersectionFound = true;
					break;
				}
			}

			if (!intersectionFound) {
				UniqueVertexData v1data = UniqueVertexData(finalMorphMap[v1].vertexType, v1, false);
				UniqueVertexData v2data = UniqueVertexData(finalMorphMap[v2].vertexType, v2, false);
				UniqueEdgeData e = UniqueEdgeData(v1data, v2data, false);

				uniqueEdges.push_back(e);

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
		it++;
	}
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
				UniqueVertexData(VertexType::Merged, intersections[0].vertexIndex, false),
				false
			)
		);

		for (size_t i = 0; i < intersections.size() - 1; i++)
		{
			uniqueEdges.push_back(
				UniqueEdgeData(
					VertexType::Source,
					UniqueVertexData(VertexType::Merged, intersections[i].vertexIndex, false),
					UniqueVertexData(VertexType::Merged, intersections[i + 1].vertexIndex, false),
					false
				)
			);
		}

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Source,
				UniqueVertexData(VertexType::Merged, intersections[intersections.size() - 1].vertexIndex, false),
				UniqueVertexData(VertexType::Source, e.v2.eqClass, e.v2.isBorder),
				false
			)
		);
	}

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

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Target,
				UniqueVertexData(VertexType::Target, e.v1.eqClass + vertexCount, e.v1.isBorder),
				UniqueVertexData(VertexType::Merged, intersections[0].vertexIndex, false),
				false
			)
		);

		for (size_t i = 0; i < intersections.size() - 1; i++)
		{
			uniqueEdges.push_back(
				UniqueEdgeData(
					VertexType::Target,
					UniqueVertexData(VertexType::Merged, intersections[i].vertexIndex, false),
					UniqueVertexData(VertexType::Merged, intersections[i + 1].vertexIndex, false),
					false
				)
			);
		}

		uniqueEdges.push_back(
			UniqueEdgeData(
				VertexType::Target,
				UniqueVertexData(VertexType::Merged, intersections[intersections.size() - 1].vertexIndex, false),
				UniqueVertexData(VertexType::Target, e.v2.eqClass + vertexCount, e.v2.isBorder),
				false
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
	}
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

//SuperMesh* HarmonicMapper::generateSuperMesh() {
//
//	std::cout << "generating super mesh..." << std::endl;
//	std::cout << "uniqueEdges.size() = " << uniqueEdges.size() << std::endl;
//	fixUniqueEdges();
//	std::cout << "uniqueEdges.size() = " << uniqueEdges.size() << std::endl;
//	std::vector <SuperVertex>* vertices = new std::vector <SuperVertex>();
//	std::vector <GLuint>* indices = new std::vector <GLuint>();
//	std::cout << "vectors created successfully!" << std::endl;
//
//	int n = 0;
//
//	for (size_t i = 0; i < uniqueEdges.size() - 2; i++)
//	{
//		auto e1 = uniqueEdges[i];
//
//		int triangles_desired_amount = e1.isBorder ? 1 : 2;
//
//		for (size_t j = i + 1; j < uniqueEdges.size() - 1; j++)
//		{
//			if (i == j) {
//				continue;
//			}
//
//			auto e2 = uniqueEdges[j];
//			for (size_t k = j + 1; k < uniqueEdges.size(); k++)
//			{
//				if (k == i || k == j) {
//					continue;
//				}
//
//				auto e3 = uniqueEdges[k];
//
//				if (e1.isTriangle(e2, e3)) {
//					int v3_ind = 0;
//					if (e3.v1.eqClass != e1.v1.eqClass && e3.v1.eqClass != e1.v2.eqClass) {
//						v3_ind = e3.v1;
//					}
//					else {
//						v3_ind = e3.v2;
//					}
//
//					SuperVertex v1(finalMorphMap[e1.v1].srcPos, finalMorphMap[e1.v1].tarPos);
//					SuperVertex v2(finalMorphMap[e1.v2].srcPos, finalMorphMap[e1.v2].tarPos);
//					SuperVertex v3(finalMorphMap[v3_ind].srcPos, finalMorphMap[v3_ind].tarPos);
//
//					vertices->push_back(v1);
//					indices->push_back(vertices->size() - 1);
//
//					if (isClockWise(map[e1.v1].image, map[e1.v2].image, map[v3_ind].image)) {
//						vertices->push_back(v2);
//						indices->push_back(vertices->size() - 1);
//						vertices->push_back(v3);
//						indices->push_back(vertices->size() - 1);
//					}
//					else {
//						vertices->push_back(v3);
//						indices->push_back(vertices->size() - 1);
//						vertices->push_back(v2);
//						indices->push_back(vertices->size() - 1);
//					}
//
//					if (n++ % 100 == 0) {
//						//std::cout << n << std::endl;
//					}
//
//					if ((--triangles_desired_amount) <= 0) {
//						goto out;
//					}
//
//					break;
//				}
//			}
//		}
//	out:;
//	}
//
//	std::cout << "Super mesh vertices amount = " << vertices->size() << std::endl;
//	std::cout << "Super mesh indices amount = " << indices->size() << std::endl;
//
//	SuperMesh* superMesh = new SuperMesh(*vertices, *indices);
//
//	std::cout << "Super mesh generated successfully!" << std::endl;
//
//	return superMesh;
//}

SuperMesh* HarmonicMapper::generateSuperMesh()
{
	std::cout << "Super mesh vertices amount = " << superVertices.size() << std::endl;
	std::cout << "Super mesh indices amount = " << superIndices.size() << std::endl;

	SuperMesh* superMesh = new SuperMesh(superVertices, superIndices);

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