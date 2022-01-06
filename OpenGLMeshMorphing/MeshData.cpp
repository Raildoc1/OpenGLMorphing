#include "MeshData.h"

MeshData::MeshData(Mesh& mesh)
{
	this->mesh = mesh;
	this->vertexCount = mesh.vertices.size();
	this->edgesCount = mesh.indices.size();

	vertices = (VertexData*)calloc(vertexCount, sizeof(VertexData));
	edges = (EdgeData*)calloc(edgesCount, sizeof(EdgeData));

	border = std::vector<EdgeData>();
}

void MeshData::init()
{
	initVertices();
	initEdges();
	initBorder();

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
				std::cout << i << ": " << "(" << edges[i].v1.eqClass << ", " << edges[i].v2.eqClass << ")" << " == " 
					<< i << ": " << "(" << edges[j].v1.eqClass << ", " << edges[j].v2.eqClass << ")" << std::endl;
				isDuplicate = true;
				break;
			}
		}

		if (!isDuplicate) {
			border.push_back(edges[i]);
		}
	}
}