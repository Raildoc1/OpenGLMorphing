#include "MeshData.h"

MeshData::MeshData(Mesh& mesh)
{
	this->mesh = mesh;
	this->vertexCount = mesh.vertices.size();
	this->edgesCount = mesh.indices.size();

	vertices = (VertexData*)calloc(vertexCount, sizeof(VertexData));
	edges = (EdgeData*)calloc(edgesCount, sizeof(EdgeData));
	k = (float**)calloc(vertexCount, sizeof(float*));
	triangles = (int*)calloc(mesh.indices.size(), sizeof(int));

	for (size_t i = 0; i < vertexCount; i++)
	{
		k[i] = (float*)calloc(vertexCount, sizeof(float));
	}

	border = std::vector<EdgeData>();
	uniqueEdges = std::vector<UniqueEdgeData>();
	fixedIndices = std::vector<int>();
	map = std::map<int, MapEntity>();
	derivatives = std::map<int, glm::vec2>();
	borderVertices = std::vector<BorderVertex>();
}

MeshData::~MeshData()
{
	for (size_t i = 0; i < vertexCount; i++)
	{
		delete k[i];
	}

	delete k;
	delete vertices;
	delete edges;
}

void MeshData::init()
{
	initVertices();
	initFixedIndices();
	initEdges();
	initBorder();
	initUniqueEdges();
	initHarmonicK();
	initMap();
	lastEnergy = calculateMapEnergy();

	initialized = true;
}

void MeshData::initEdges()
{
	for (size_t i = 0; i < edgesCount; i += 3)
	{
		EdgeData e1;
		EdgeData e2;
		EdgeData e3;

		e1.v1 = vertices[mesh.indices[i + 0]];
		e1.v2 = vertices[mesh.indices[i + 1]];
		e1.length = glm::distance(e1.v1.vertex.position, e1.v2.vertex.position);

		e2.v1 = vertices[mesh.indices[i + 1]];
		e2.v2 = vertices[mesh.indices[i + 2]];
		e2.length = glm::distance(e2.v1.vertex.position, e2.v2.vertex.position);

		e3.v1 = vertices[mesh.indices[i + 2]];
		e3.v2 = vertices[mesh.indices[i + 0]];
		e3.length = glm::distance(e3.v1.vertex.position, e3.v2.vertex.position);

		edges[i + 0] = e1;
		edges[i + 1] = e2;
		edges[i + 2] = e3;
	}
}

void MeshData::initVertices()
{
	for (size_t i = 0; i < vertexCount; i++)
	{
		VertexData data;
		data.vertex = mesh.vertices[i];
		data.index = i;
		data.eqClass = -1;
		data.isBorder = false;

		for (size_t j = 0; j < vertexCount; j++)
		{
			if (glm::distance2(mesh.vertices[i].position, mesh.vertices[j].position) < EPSILON) {
				data.eqClass = j;
				break;
			}
		}

		vertices[i] = data;
	}
}

void MeshData::initFixedIndices()
{
	indicesCount = mesh.indices.size();
	for (size_t i = 0; i < indicesCount; i++)
	{
		triangles[i] = vertices[mesh.indices[i]].eqClass;
		fixedIndices.push_back(vertices[mesh.indices[i]].eqClass);
	}
}

void MeshData::initBorder()
{
	for (size_t i = 0; i < edgesCount; i++)
	{
		bool isDuplicate = false;

		for (size_t j = 0; j < edgesCount; j++)
		{
			if (i == j) {
				continue;
			}

			if (edges[i].equals(edges[j])) {
				/*std::cout << i << ": " << "(" << edges[i].v1.eqClass << ", " << edges[i].v2.eqClass << ")" << " == " 
					<< i << ": " << "(" << edges[j].v1.eqClass << ", " << edges[j].v2.eqClass << ")" << std::endl;*/
				isDuplicate = true;
				break;
			}
		}

		bool isBorder = !isDuplicate;
		edges[i].isBorder = isBorder;

		if (isBorder) {
			border.push_back(edges[i]);
			vertices[edges[i].v1.eqClass].isBorder = true;
			vertices[edges[i].v2.eqClass].isBorder = true;
		}
	}
	sortBorder();

	borderLength = 0.0f;

	for (size_t i = 0; i < border.size(); i++)
	{
		borderLength += border[i].length;
	}
}

void MeshData::sortBorder()
{
	for (size_t i = 0; i < border.size() - 1; i++)
	{
		for (size_t j = i + 1; j < border.size(); j++)
		{
			if (border[i].adjacent(border[j]) && (i + 1 != j)) {
				EdgeData t = border[i + 1];
				border[i + 1] = border[j];
				border[j] = t;
				break;
			}
		}
	}
}

void MeshData::initUniqueEdges()
{
	for (size_t i = 0; i < edgesCount; i++)
	{
		UniqueEdgeData e;

		e.v1 = edges[i].v1.eqClass;
		e.v2 = edges[i].v2.eqClass;
		e.pos1 = edges[i].v1.vertex.position;
		e.pos2 = edges[i].v2.vertex.position;
		e.isBorder = edges[i].isBorder;

		bool isDuplicate = false;
		for (size_t i = 0; i < uniqueEdges.size(); i++)
		{
			if (uniqueEdges[i].equals(e)) {
				isDuplicate = true;
				break;
			}
		}

		if (!isDuplicate) {
			uniqueEdges.push_back(e);
		}
	}
}

void MeshData::initHarmonicK()
{
	for (size_t i = 0; i < vertexCount; i++)
	{
		for (size_t j = 0; j < vertexCount; j++)
		{
			k[i][j] = 0.0f;
		}
	}

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		if (uniqueEdges[i].isBorder) {
			continue;
		}

		bool firstFound = false;
		int k1 = -1;
		int k2 = -1;

		for (size_t j = 0; j < mesh.indices.size(); j += 3)
		{
			int v1 = vertices[mesh.indices[j + 0]].eqClass;
			int v2 = vertices[mesh.indices[j + 1]].eqClass;
			int v3 = vertices[mesh.indices[j + 2]].eqClass;

			int e1 = uniqueEdges[i].v1;
			int e2 = uniqueEdges[i].v2;

			int k0 = edgeInTriangle(e1, e2, v1, v2, v3);

			if (k0 != -1) {
				if (!firstFound) {
					firstFound = true;
					k1 = k0;
				} else {
					k2 = k0;
					break;
				}
			}
		}
		
		if (k1 == -1 || k2 == -1) {
			std::cerr << "Triangle not found!" << std::endl;
		}

		glm::vec3 vk1 = vertices[k1].vertex.position;
		glm::vec3 vk2 = vertices[k2].vertex.position;
		glm::vec3 vi = uniqueEdges[i].pos1;
		glm::vec3 vj = uniqueEdges[i].pos2;

		std::cout << "vertices[" << k1 << "].vertex.position = " << "(" << vk1.x << ", " << vk1.y << ")" << std::endl;

		float lik1 = glm::distance2(vi, vk1);
		float ljk1 = glm::distance2(vj, vk1);
		float lik2 = glm::distance2(vi, vk2);
		float ljk2 = glm::distance2(vj, vk2);
		float lji = glm::distance2(vi, vj);

		glm::vec3 ik1 = vk1 - vi;
		glm::vec3 jk1 = vk1 - vj;
		glm::vec3 ik2 = vk2 - vi;
		glm::vec3 jk2 = vk2 - vj;

		float Aijk1 = glm::cross(ik1, jk1).length() / 2.0f;
		float Aijk2 = glm::cross(ik2, jk2).length() / 2.0f;

		std::cout << "lik1 = (" << vi.x <<  ", " << vi.y << ") to (" << vk1.x << ", " << vk1.y << ") = " << lik1 << "; ljk1 = " << ljk1 << "; lik2 = " << lik2 << "; ljk2 = " << ljk2 << "; lji = " << lji << std::endl;
		std::cout << "Aijk1 = " << Aijk1 << "; Aijk2 = " << Aijk2 << std::endl;

		k[uniqueEdges[i].v1][uniqueEdges[i].v2] = k[uniqueEdges[i].v2][uniqueEdges[i].v1] = ((lik1 + ljk1 - lji) / Aijk1) + ((lik2 + ljk2 - lji) / Aijk2);
	}

	for (size_t i = 0; i < vertexCount; i++)
	{
		for (size_t j = 0; j < vertexCount; j++)
		{
			if (i == j) {
				k[i][j] = 0;
			}

			k[i][j] = k[vertices[i].eqClass][vertices[j].eqClass];
		}
	}
}

int MeshData::edgeInTriangle(int e1, int e2, int t1, int t2, int t3)
{
	if ((e1 == t1 && e2 == t2) ||
		(e2 == t1 && e1 == t2)) {
		return t3;
	}

	if ((e1 == t1 && e2 == t3) ||
		(e2 == t1 && e1 == t3)) {
		return t2;
	}

	if ((e1 == t3 && e2 == t2) ||
		(e2 == t3 && e1 == t2)) {
		return t1;
	}

	return -1;
}

void MeshData::initMap()
{
	std::cout << "Border vertices:" << std::endl;

	float currentLength = 0.0f;
	float dPi = glm::pi<float>() * 2.0f;
	float borderLength = getBorderLength();

	for (size_t i = 0; i < border.size(); i++)
	{
		float t = dPi * currentLength / borderLength;

		MapEntity e;
		e.locked = true;
		e.image = glm::vec2(glm::cos(t), glm::sin(t));
		e.border = true;
		e.phi = t;

		currentLength += border[i].length;

		map[border[i].v1.eqClass] = e;

		BorderVertex v;
		v.eqClass = border[i].v1.eqClass;
		v.phi = t;

		borderVertices.push_back(v);

		std::cout << "phi = " << t << std::endl;
	}

	//std::cout << "Border length: " << getBorderLength() << std::endl;

	for (size_t i = 0; i < getVertexCount(); i++)
	{
		int eqClass = vertices[i].eqClass;
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

		float t = dPi * eqClass / (float)vertexCount;
		//std::cout << "t[" << eqClass << "] = " << "(" << t << std::endl;

		MapEntity e;
		//e.image = 0.5f * glm::vec2(glm::cos(t), glm::sin(t));
		e.image = glm::vec2(0.0f, 0.0f);
		e.locked = false;
		e.border = false;
		e.phi = 0.0f;

		//std::cout << "image = " << e.image.x << ", " << e.image.y << std::endl;

		map[eqClass] = e;
	}

}

float MeshData::calculateMapEnergy()
{
	if (!initialized) {
		return 0.0f;
	}

	float result = 0.0f;

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		derivatives[uniqueEdges[i].v1] = glm::vec2(0.0f, 0.0f);
		derivatives[uniqueEdges[i].v2] = glm::vec2(0.0f, 0.0f);
	}

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		int v1 = uniqueEdges[i].v1;
		int v2 = uniqueEdges[i].v2;
		glm::vec2 delta = k[v1][v2] * (map[v1].image - map[v2].image);

		derivatives[v1] += delta;
		derivatives[v2] -= delta;

		//std::cout << "derivatives[" << v1 << "] = " << "(" << derivatives[v1].x << ", " << derivatives[v1].y << ") - k[" << v1 << "][" << v2 << "] = " << k[v1][v2] << std::endl;

		result += k[v1][v2] * glm::distance2(map[v1].image, map[v2].image);
	}

	return result / 2.0f;
}

float MeshData::tickMap()
{
	if (!initialized) {
		return 0.0f;
	}

	float energy = calculateMapEnergy();
	float energyDelta = energy - lastEnergy;

	lastEnergy = energy;

	for (auto const& x : map)
	{
		if (x.second.locked) {
			continue;
		}

		map[x.first].image -= 0.001f * derivatives[x.first];
	}

	return energyDelta;
}

void MeshData::harmonizeMap()
{
	int iterations = 0;
	float delta;
	do
	{
		iterations++;
		delta = tickMap();

		//std::cout << "iterations = " << iterations << "; delta = " << delta << std::endl;

		if (delta < 0.0f) {
			delta = -delta;
		}

	} while (iterations < MAX_ITER && delta > 0.01f);

	//std::cout << "iterations = " << iterations << "; delta = " << delta << std::endl;
}

glm::vec3 MeshData::findVertexPos(glm::vec2 mapPos) {
	for (size_t i = 0; i < indicesCount; i += 3)
	{
		glm::vec2 p1 = map[triangles[i + 0]].image;
		glm::vec2 p2 = map[triangles[i + 1]].image;
		glm::vec2 p3 = map[triangles[i + 2]].image;

		float area = 0.5 * (-p2.y * p3.x + p1.y * (-p2.x + p3.x) + p1.x * (p2.y - p3.y) + p2.x * p3.y);

		float s = 1 / (2 * area) * (p1.y * p3.x - p1.x * p3.y + (p3.y - p1.y) * mapPos.x + (p1.x - p3.x) * mapPos.y);
		float t = 1 / (2 * area) * (p1.x * p2.y - p1.y * p2.x + (p1.y - p2.y) * mapPos.x + (p2.x - p1.x) * mapPos.y);

		glm::vec3 v1 = vertices[triangles[i + 0]].vertex.position;
		glm::vec3 v2 = vertices[triangles[i + 1]].vertex.position;
		glm::vec3 v3 = vertices[triangles[i + 2]].vertex.position;

		if (s > 0 && t > 0 && 1 - s - t > 0) {
			return v1 + (v2 - v1) * s + (v3 - v1) * t;
		}
	}
}

glm::vec3 MeshData::findBorderPos(float phi) {
	for (size_t i = 0; i < borderVertices.size(); i++)
	{
		int j = (i + 1) % borderVertices.size();

		float phi1 = borderVertices[i].phi;
		float phi2 = borderVertices[j].phi;

		if (phi1 < phi && phi < phi2) {
			float t = (phi - phi1) / (phi2 - phi1);
			return vertices[borderVertices[i].eqClass].vertex.position * t +
				vertices[borderVertices[j].eqClass].vertex.position * (1.0f - t);
		} else if (borderVertices[i].phi == phi) {
			return vertices[borderVertices[i].eqClass].vertex.position;
		}
	}

	std::cout << "border position isn't found :c" << std::endl;

	return glm::vec3();
}