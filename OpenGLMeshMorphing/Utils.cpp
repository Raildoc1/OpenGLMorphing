#include "Utils.h"
#include <iostream>

bool Utils::tryFindIntersection(vec2 a, vec2 b, vec2 c, vec2 d, vec2* intersection, bool exclusively)
{
	float v1 = (d.x - c.x) * (a.y - c.y) - (d.y - c.y) * (a.x - c.x);
	float v2 = (d.x - c.x) * (b.y - c.y) - (d.y - c.y) * (b.x - c.x);
	float v3 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
	float v4 = (b.x - a.x) * (d.y - a.y) - (b.y - a.y) * (d.x - a.x);

	float z1 = (b.y - a.y) / (b.x - a.x);
	float z2 = (d.y - c.y) / (d.x - c.x);

	*intersection = vec2();

	intersection->x = (c.y - a.y - z2 * c.x + z1 * a.x) / (z1 - z2);
	intersection->y = (b.y - a.y) * (intersection->x - a.x) / (b.x - a.x) + a.y;

	if (exclusively) {
		return (v1 * v2 < 0.0f) && (v3 * v4 < 0.0f);
	}

	return ((v1 * v2 < 0.0f) && (v3 * v4 <= 0.0f)) || ((v1 * v2 <= 0.0f) && (v3 * v4 < 0.0f));
}

bool Utils::pointInPolygon(vec2& point, vector<vec2>& polygon)
{
	int intersections = 0;

	for (size_t i = 0; i < polygon.size(); i++)
	{
		size_t j = (i + 1) % polygon.size();
		if (polygon[i].x > point.x && polygon[j].x > point.x) {
			continue;
		}
		if (polygon[i].x < point.x && polygon[j].x < point.x) {
			continue;
		}
		glm::vec2 v1 = polygon[i] - point;
		glm::vec2 v2 = polygon[j] - point;
		float intersection = (v1.y * v2.x - v1.x * v2.y) / (v2.x - v1.x);
		if (intersection > 0.0f) {
			intersections++;
		}
	}

	return intersections % 2 == 1;
}
