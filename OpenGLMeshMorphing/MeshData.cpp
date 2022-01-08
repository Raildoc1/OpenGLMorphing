#include "MeshData.h"

MeshData::MeshData(Mesh& mesh)
{
	this->mesh = mesh;
	this->vertexCount = mesh.vertices.size();
	this->edgesCount = mesh.indices.size();

	vertices = (VertexData*)calloc(vertexCount, sizeof(VertexData));
	edges = (EdgeData*)calloc(edgesCount, sizeof(EdgeData));
	k = (float**)calloc(vertexCount, sizeof(float*));

	for (size_t i = 0; i < vertexCount; i++)
	{
		k[i] = (float*)calloc(vertexCount, sizeof(float));
	}

	border = std::vector<EdgeData>();
	uniqueEdges = std::vector<UniqueEdgeData>();
	fixedIndices = std::vector<int>();
}

void MeshData::init()
{
	initVertices();
	initFixedIndices();
	initEdges();
	initBorder();
	initUniqueEdges();
	initHarmonicK();

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

		e2.v1 = vertices[mesh.indices[i + 1]];
		e2.v2 = vertices[mesh.indices[i + 2]];

		e3.v1 = vertices[mesh.indices[i + 2]];
		e3.v2 = vertices[mesh.indices[i + 0]];

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
	for (size_t i = 0; i < mesh.indices.size(); i++)
	{
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
			int v1 = mesh.indices[j + 0];
			int v2 = mesh.indices[j + 1];
			int v3 = mesh.indices[j + 2];

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

		glm::vec3 vk1 = vertices[k1].vertex.position;
		glm::vec3 vk2 = vertices[k2].vertex.position;
		glm::vec3 vi = uniqueEdges[i].pos1;
		glm::vec3 vj = uniqueEdges[i].pos2;

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