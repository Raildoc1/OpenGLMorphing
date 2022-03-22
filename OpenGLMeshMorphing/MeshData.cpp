#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include "MeshData.h"
#include <filesystem>
#include <Eigen/Dense>

void print_time_stamp(const clock_t& prev, const std::string msg) {
	std::cout << std::setw(30) << msg << std::setw(30) << float(clock() - prev) / CLOCKS_PER_SEC << " seconds." << std::endl;
}

float* solveSLE(float** A, float* b, int n) {
	float** AA = (float**)calloc(n, sizeof(float*));
	for (size_t i = 0; i < n; i++)
	{
		AA[i] = (float*)calloc(n + 1, sizeof(float));
		memcpy(AA[i], A[i], n * sizeof(float));
		AA[i][n] = b[i];
	}

	for (int i = 0; i < n - 1; i++) {
		for (int h = i + 1; h < n; h++)
		{
			float t = AA[h][i] / AA[i][i];
			for (int j = 0; j <= n; j++)
			{
				AA[h][j] = AA[h][j] - t * AA[i][j];
			}
		}
	}

	float* result = (float*)calloc(n, sizeof(float));

	for (int i = n - 1; i >= 0; i--)
	{
		result[i] = AA[i][n];
		for (int j = n - 1; j > i; j--)
		{
			result[i] = result[i] - AA[i][j] * result[j];
		}
		result[i] = result[i] / AA[i][i];
	}

	return result;
}

MeshData::MeshData(Mesh& mesh, float rotation = 0.0f, bool invertBorder = false)
{
	this->mesh = mesh;
	this->vertexCount = mesh.vertices.size();
	this->edgesCount = mesh.indices.size();
	this->rotation = rotation;
	this->invertBorder = invertBorder;

	vertices = (VertexData*)calloc(vertexCount, sizeof(VertexData));
	edges = (EdgeData*)calloc(edgesCount, sizeof(EdgeData));
	k = (float**)calloc(vertexCount, sizeof(float*));
	lambda = (float**)calloc(vertexCount, sizeof(float*));
	triangles = (int*)calloc(mesh.indices.size(), sizeof(int));

	for (size_t i = 0; i < vertexCount; i++)
	{
		k[i] = (float*)calloc(vertexCount, sizeof(float));
		lambda[i] = (float*)calloc(vertexCount, sizeof(float));
	}

	border = std::vector<EdgeData>();
	uniqueEdges = std::vector<UniqueEdgeData>();
	fixedIndices = std::vector<int>();
	map = std::map<int, MapEntity>();
	derivatives = std::map<int, glm::vec2>();
	borderVertices = std::vector<BorderVertex>();
	eqClassesSet = std::set<int>();
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
	clock_t time_stamp = clock();

	initVertices();
	print_time_stamp(time_stamp, "initVertices finished");
	time_stamp = clock();

	initFixedIndices();
	print_time_stamp(time_stamp, "initFixedIndices finished");
	time_stamp = clock();

	initEdges();
	print_time_stamp(time_stamp, "initEdges finished");
	time_stamp = clock();

	initBorder();
	print_time_stamp(time_stamp, "initBorder finished");
	time_stamp = clock();

	initUniqueEdges();
	print_time_stamp(time_stamp, "initUniqueEdges finished");
	time_stamp = clock();

	initHarmonicK();
	print_time_stamp(time_stamp, "initHarmonicK finished");
	time_stamp = clock();

	initLambda();
	print_time_stamp(time_stamp, "initLambda finished");
	time_stamp = clock();

	initMap();
	print_time_stamp(time_stamp, "initMap finished");
	time_stamp = clock();

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
			glm::vec3 pos1 = mesh.vertices[i].position;
			glm::vec3 pos2 = mesh.vertices[j].position;
			if (glm::distance2(pos1, pos2) < EPSILON) {
				data.eqClass = j;
				eqClassesSet.insert(j);
				break;
			}
		}

		vertices[i] = data;
	}

	std::vector<int> v;
	v.reserve(eqClassesSet.size());
	std::copy(eqClassesSet.begin(), eqClassesSet.end(), std::back_inserter(v));

	for (size_t i = 0; i < vertexCount; i++)
	{
		for (size_t j = 0; j < v.size(); j++)
		{
			if (vertices[i].eqClass == v[j]) {
				vertices[i].setIndex = j;
				break;
			}
		}
	}

	std::cout << "set size = " << eqClassesSet.size() << std::endl;
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
			eqClassesSet.erase(edges[i].v1.eqClass);
			eqClassesSet.erase(edges[i].v2.eqClass);
		}
	}
	sortBorder();

	borderLength = 0.0f;

	for (size_t i = 0; i < border.size(); i++)
	{
		borderLength += border[i].length;
	}

	vertexSetList.reserve(eqClassesSet.size());
	std::copy(eqClassesSet.begin(), eqClassesSet.end(), std::back_inserter(vertexSetList));

	for (size_t i = 0; i < vertexCount; i++)
	{
		vertices[i].setIndex = -1;
		for (size_t j = 0; j < vertexSetList.size(); j++)
		{
			if (vertices[i].eqClass == vertexSetList[j]) {
				vertices[i].setIndex = j;
				break;
			}
		}
	}

	looseVerteicesAmount = eqClassesSet.size();

	std::cout << "loose verteices amount = " << looseVerteicesAmount << std::endl;
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
		UniqueVertexData v1 = UniqueVertexData(edges[i].v1.eqClass, edges[i].v1.isBorder);
		UniqueVertexData v2 = UniqueVertexData(edges[i].v2.eqClass, edges[i].v2.isBorder);
		UniqueEdgeData e = UniqueEdgeData(v1, v2);

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
		bool firstFound = false;
		int k1 = -1;
		int k2 = -1;

		for (size_t j = 0; j < mesh.indices.size(); j += 3)
		{
			int v1 = vertices[mesh.indices[j + 0]].eqClass;
			int v2 = vertices[mesh.indices[j + 1]].eqClass;
			int v3 = vertices[mesh.indices[j + 2]].eqClass;

			int e1 = uniqueEdges[i].v1.eqClass;
			int e2 = uniqueEdges[i].v2.eqClass;

			int k0 = edgeInTriangle(e1, e2, v1, v2, v3);

			if (k0 != -1) {
				if (!firstFound) {
					firstFound = true;
					k1 = k0;
				}
				else {
					k2 = k0;
					break;
				}
			}
		}

		if (k1 == -1 && k2 == -1) {
			std::cerr << "Triangle not found!" << std::endl;
		}
		else {
			float result = 0.0f;

			glm::vec3 vi = vertices[uniqueEdges[i].v1.eqClass].vertex.position;
			glm::vec3 vj = vertices[uniqueEdges[i].v2.eqClass].vertex.position;

			float lji = glm::distance2(vi, vj);

			if (k1 != -1) {
				glm::vec3 vk1 = vertices[k1].vertex.position;

				if (coeffType == CoeffType::Kanai) {

					float lik1 = glm::distance2(vi, vk1);
					float ljk1 = glm::distance2(vj, vk1);

					glm::vec3 ik1 = vk1 - vi;
					glm::vec3 jk1 = vk1 - vj;

					float Aijk1 = glm::cross(ik1, jk1).length() / 2.0f;

					result += (lik1 + ljk1 - lji) / Aijk1;

					if (result < 0.0f) {
						//std::cout << "Alert!" << std::endl;
						//std::cout << "lji = " << lji << std::endl;
						//std::cout << "lik1 = " << lik1 << std::endl;
						//std::cout << "ljk1 = " << ljk1 << std::endl;
					}
				}
				else if (coeffType == CoeffType::MVC) {
					float dotij = glm::dot(glm::normalize(vj - vi), glm::normalize(vk1 - vi));
					float alpha = glm::acos(dotij);

					result += glm::tan(alpha / 2.0f) / glm::distance(vi, vj);

					if (result < 0.0f) {
						//std::cout << "Alert!" << std::endl;
					}
				}
				else if (coeffType == CoeffType::DHC) {
					float dotij = glm::dot(glm::normalize(vj - vk1), glm::normalize(vi - vk1));
					float gamma = glm::acos(dotij);

					result += 1.0f / glm::tan(gamma);
				}
			}

			if (k2 != -1) {
				glm::vec3 vk2 = vertices[k2].vertex.position;

				if (coeffType == CoeffType::Kanai) {
					float lik2 = glm::distance2(vi, vk2);
					float ljk2 = glm::distance2(vj, vk2);

					glm::vec3 ik2 = vk2 - vi;
					glm::vec3 jk2 = vk2 - vj;

					float Aijk2 = glm::cross(ik2, jk2).length() / 2.0f;

					result += (lik2 + ljk2 - lji) / Aijk2;

					if (result < 0.0f) {
						//std::cout << "Alert!" << std::endl;
						//std::cout << "lji = " << lji << std::endl;
						//std::cout << "lik2 = " << lik2 << std::endl;
						//std::cout << "ljk2 = " << ljk2 << std::endl;
					}
				}
				else if (coeffType == CoeffType::MVC) {
					float dotji = glm::dot(glm::normalize(vj - vi), glm::normalize(vk2 - vi));
					float beta = glm::acos(dotji);

					result += glm::tan(beta / 2.0f) / glm::distance(vi, vj);

					if (result < 0.0f) {
						//std::cout << "Alert!" << std::endl;
					}
				}
				else if (coeffType == CoeffType::DHC) {
					float dotji = glm::dot(glm::normalize(vj - vk2), glm::normalize(vi - vk2));
					float gamma = glm::acos(dotji);

					result += 1.0f / glm::tan(gamma);
				}
			}

			if (coeffType == CoeffType::One) {
				k[uniqueEdges[i].v1.eqClass][uniqueEdges[i].v2.eqClass] = k[uniqueEdges[i].v2.eqClass][uniqueEdges[i].v1.eqClass] = 1.0f;
			} else {
				k[uniqueEdges[i].v1.eqClass][uniqueEdges[i].v2.eqClass] = k[uniqueEdges[i].v2.eqClass][uniqueEdges[i].v1.eqClass] = result;
			}
		}
	}

	for (size_t i = 0; i < vertexCount; i++)
	{
		for (size_t j = 0; j < vertexCount; j++)
		{
			if (i == j) {
				k[i][j] = 0;
			}

			k[i][j] = k[vertices[i].eqClass][vertices[j].eqClass];
			//std::cout << std::setw(6) << std::setprecision(2) << k[i][j] << " ";
		}
		//std::cout << std::endl;
	}
	//std::cout << std::endl;
}

void MeshData::initLambda()
{
	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		int i0 = uniqueEdges[i].v1.eqClass;
		int i1 = uniqueEdges[i].v2.eqClass;
		float current_k = k[i0][i1];

		lambda[i0][i1] += current_k;
		lambda[i1][i0] += current_k;

		for (size_t j = 0; j < uniqueEdges.size(); j++)
		{
			if (i == j) {
				continue;
			}

			int j0 = uniqueEdges[j].v1.eqClass;
			int j1 = uniqueEdges[j].v2.eqClass;

			if (j0 == i0 || j0 == i1) {
				lambda[j0][j1] += current_k;
				/*if (j0 == 1080 && j1 == 966) {
					std::cout << lambda[j0][j1] << " " << current_k << std::endl;
				}*/
			}
			else if (j1 == i0 || j1 == i1) {
				lambda[j1][j0] += current_k;
				/*if (j1 == 1080 && j0 == 966) {
					std::cout << lambda[j0][j1] << " " << current_k << std::endl;
				}*/
			}
		}

	}

	//std::cout << "lambdas: " << std::endl;
	for (size_t i = 0; i < vertexCount; i++)
	{
		for (size_t j = 0; j < vertexCount; j++)
		{
			if (lambda[i][j] == 0.0f) {
				//std::cout << std::setw(6) << std::setprecision(2) << lambda[i][j] << " ";
				continue;
			}

			lambda[i][j] = k[i][j] / lambda[i][j];
			//std::cout << std::setw(6) << std::setprecision(2) << lambda[i][j] << " ";
			/*if (i == 1080 && j == 966) {
				std::cout << "lambda = " << lambda[i][j] << " " << k[i][j] << std::endl;
			}*/
			if (lambda[i][j] < 0.0f) {
				//std::cout << "lambda[" << i << "][" << j << "] = " << lambda[i][j] << std::endl;
			}
			//std::cout << "lambda[" << i << "][" << j << "] = " << lambda[i][j] << " ";
		}
		//std::cout << std::endl;
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
	Eigen::MatrixXf A = Eigen::MatrixXf::Identity(looseVerteicesAmount, looseVerteicesAmount);
	Eigen::VectorXf U = Eigen::VectorXf::Zero(looseVerteicesAmount);
	Eigen::VectorXf V = Eigen::VectorXf::Zero(looseVerteicesAmount);

	float currentLength = 0.0f;
	float dPi = glm::pi<float>() * 2.0f;
	float borderLength = getBorderLength();

	for (size_t i = 0; i < border.size(); i++)
	{
		size_t j = i;

		if (invertBorder) {
			j = border.size() - i - 1;
		}

		float t = dPi * currentLength / borderLength + rotation;

		while (t < 0.0f) {
			t += dPi;
		}

		while (t >= dPi) {
			t -= dPi;
		}

		MapEntity e;
		e.locked = true;
		e.image = glm::vec2(glm::cos(t), glm::sin(t));
		e.border = true;
		e.phi = t;

		currentLength += border[j].length;

		map[border[j].v1.eqClass] = e;

		BorderVertex v;
		v.eqClass = border[j].v1.eqClass;
		v.phi = t;

		borderVertices.push_back(v);
	}

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		int i0 = uniqueEdges[i].v1.eqClass;
		int i1 = uniqueEdges[i].v2.eqClass;

		VertexData v1 = vertices[i0];
		VertexData v2 = vertices[i1];

		if (v1.isBorder && v2.isBorder) {
			continue;
		}

		if (v1.isBorder && !v2.isBorder) {
			U[v2.setIndex] += lambda[v2.eqClass][v1.eqClass] * map[v1.eqClass].image.x;
			V[v2.setIndex] += lambda[v2.eqClass][v1.eqClass] * map[v1.eqClass].image.y;
		}
		else if (!v1.isBorder && v2.isBorder) {
			U[v1.setIndex] += lambda[v1.eqClass][v2.eqClass] * map[v2.eqClass].image.x;
			V[v1.setIndex] += lambda[v1.eqClass][v2.eqClass] * map[v2.eqClass].image.y;
		}
		else {
			A(v1.setIndex, v2.setIndex) = -lambda[v1.eqClass][v2.eqClass];
			A(v2.setIndex, v1.setIndex) = -lambda[v2.eqClass][v1.eqClass];
		}
	}

	std::cout << "loose vertices amount = " << looseVerteicesAmount << std::endl;

	clock_t time_stamp = clock();

	Eigen::HouseholderQR<Eigen::MatrixXf> hqr = A.householderQr();
	Eigen::VectorXf u = hqr.solve(U);
	Eigen::VectorXf v = hqr.solve(V);

	print_time_stamp(time_stamp, "partialPivLu SLE solutions found");

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

		MapEntity e;
		e.image = glm::vec2(u[vertices[eqClass].setIndex], v[vertices[eqClass].setIndex]);
		e.locked = false;
		e.border = false;
		e.phi = 0.0f;

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
		derivatives[uniqueEdges[i].v1.eqClass] = glm::vec2(0.0f, 0.0f);
		derivatives[uniqueEdges[i].v2.eqClass] = glm::vec2(0.0f, 0.0f);
	}

	for (size_t i = 0; i < uniqueEdges.size(); i++)
	{
		int v1 = uniqueEdges[i].v1.eqClass;
		int v2 = uniqueEdges[i].v2.eqClass;
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

		map[x.first].image -= 0.05f * derivatives[x.first];
	}

	std::cout << energy << std::endl;

	return energyDelta;
}

void MeshData::harmonizeMap()
{
	std::cout << "harmonizeMap" << std::endl;

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

	} while (iterations < MAX_ITER && delta > 0.001f);

	//std::cout << "iterations = " << iterations << "; delta = " << delta << std::endl;
}

glm::vec3 MeshData::findVertexPos(glm::vec2 mapPos) {

	float area = 0.0f;
	float s = 0.0f;
	float t = 0.0f;

	for (size_t i = 0; i < indicesCount; i += 3)
	{
		glm::vec2 p1 = map[fixedIndices[i + 0]].image;
		//std::cout << "read1 " << i << std::endl;
		glm::vec2 p2 = map[fixedIndices[i + 1]].image;
		//std::cout << "read1 " << i + 1 << std::endl;
		glm::vec2 p3 = map[fixedIndices[i + 2]].image;
		//std::cout << "read1 " << i + 2 << std::endl;

		area = 0.5 * (-p2.y * p3.x + p1.y * (-p2.x + p3.x) + p1.x * (p2.y - p3.y) + p2.x * p3.y);

		s = 1 / (2 * area) * (p1.y * p3.x - p1.x * p3.y + (p3.y - p1.y) * mapPos.x + (p1.x - p3.x) * mapPos.y);
		t = 1 / (2 * area) * (p1.x * p2.y - p1.y * p2.x + (p1.y - p2.y) * mapPos.x + (p2.x - p1.x) * mapPos.y);

		glm::vec3 v1 = vertices[fixedIndices[i + 0]].vertex.position;
		//std::cout << "read2 " << i << std::endl;
		glm::vec3 v2 = vertices[fixedIndices[i + 1]].vertex.position;
		//std::cout << "read2 " << i + 1 << std::endl;
		glm::vec3 v3 = vertices[fixedIndices[i + 2]].vertex.position;
		//std::cout << "read2 " << i + 2 << std::endl;

		if (s > -EPSILON && t > -EPSILON && 1 - s - t > -EPSILON) {

			if (area < 0.0f) {
				std::cout << "ERROR: findVertexPos: area is negative! {area = " << area << ", s = " << s << ", t = " << t << "}" << std::endl;
			}

			if (s < 0.0f || t < 0.0f) {
				std::cout << "ERROR: findVertexPos: coordinates are negative! {area = " << area << ", s = " << s << ", t = " << t << "}" << std::endl;
			}

			//std::cout << "{area = " << area << ", s = " << s << ", t = " << t << ", s + t = " << s + t << "}" << std::endl;

			return v1 + (v2 - v1) * s + (v3 - v1) * t;
		}
	}

	std::cout << "ERROR: findVertexPos: could not find vertex position! {area = " << area << ", s = " << s << ", t = " << t << "}" << std::endl;
	return glm::vec3(0.0f);
}

glm::vec3 MeshData::findBorderPos(float phi) {

	float pi = glm::pi<float>();

	std::cout << "looking for position..." << std::endl;

	int i = 0;
	while (i++ < 2) {
		for (size_t i = 0; i < borderVertices.size(); i++)
		{
			int j = (i + 1) % borderVertices.size();

			float phi1 = borderVertices[i].phi;
			float phi2 = borderVertices[j].phi;

			while (phi1 >= pi * 2.0f) {
				phi1 -= pi * 2.0f;
			}

			while (phi2 >= pi * 2.0f) {
				phi2 -= pi * 2.0f;
			}

			while (phi1 < 0.0f) {
				phi1 += pi * 2.0f;
			}

			while (phi2 < 0.0f) {
				phi2 += pi * 2.0f;
			}

			while (phi2 < phi1) {
				phi2 += pi * 2.0f;
			}

			//std::cout << phi1 << " < " << phi << " < " << phi2 << std::endl;
			if (phi1 <= phi && phi <= phi2) {
				float t = (phi - phi1) / (phi2 - phi1);
				return vertices[borderVertices[j].eqClass].vertex.position * t +
					vertices[borderVertices[i].eqClass].vertex.position * (1.0f - t);
			}
		}

		phi += pi * 2.0f;
	}

	std::cout << "border position isn't found :c" << std::endl;

	return glm::vec3();
}