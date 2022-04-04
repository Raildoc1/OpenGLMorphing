#pragma once
#include<glm/glm.hpp>
#include <vector>

using std::vector, glm::vec2;

class Utils
{
public:
	static bool tryFindIntersection(vec2 a, vec2 b, vec2 c, vec2 d, vec2* intersection, bool exclusively);
	static bool pointInPolygon(vec2& point, vector<vec2>& polygon);
};

