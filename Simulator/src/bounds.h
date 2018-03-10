#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "globals.h"

class Bounds
{
public:
	Point3f min, max;

	Bounds()
	{
        double inf = std::numeric_limits<T>::infinity();
		min << inf, inf, inf;
		max << -inf, -inf, -inf;
	}

	Bounds(const Point3f& min, const Point3f& max) : min(min), max(max)
	{}

	Bounds(const Point3f& p) : min(p), max(p)
	{}

	// Returns a vector representing the diagonal of the box
	Vector3f Diagonal() const
    {
        return max - min;
    }

	// Returns the index of the axis with the largest length
	int MaximumExtent() const;

	// Transforms this Bounds3f by the input Transform and also returns a Bounds3f representing our new boundaries.
	// Transforming a Bounds3f is equivalent to creating a Bounds3f that encompasses the non-axis-aligned box resulting from
	// transforming this Bounds3f's eight corners
	Bounds ApplyTransform(const Matrix4f& tr);

	bool Intersect(const Ray& r , float* t) const;
	bool Contains(const Vector3f& pos);
};

bool Intersect_AABB_with_AABB(const Bounds& a, const Bounds& b);
Bounds Union(const Bounds& b1, const Bounds& b2);
Bounds Union(const Bounds& b1, const Point3f& p);
Bounds Union(const Bounds& b1, const Point4f& p);