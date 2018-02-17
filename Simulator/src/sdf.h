#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <stdlib.h>

using T = double;

namespace SDF
{
    // p must be in local space
    float inline sdPlane(Eigen::Matrix<T, 3, 1> p, Eigen::Matrix<T, 4, 1> n)
    {
        // n must be normalized
        Eigen::Matrix<T, 3, 1> _n(n[0], n[1], n[2]);
        return p.dot(_n) + n[3];
    }

    // p must be in local space
    float inline sdSphere(Eigen::Matrix<T, 3, 1> p, float s)
    {
        return p.norm() - s;
    }
}
