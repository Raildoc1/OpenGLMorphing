#pragma once
#include<glm/glm.hpp>

class Utils
{
public:
	static bool tryFindIntersection(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d, glm::vec2* intersection, bool exclusively);
};

