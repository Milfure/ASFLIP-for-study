#pragma once

#include <cmath>
#include "Utilize/Util.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"

template<class Type>
using Vector3 = Eigen::Vector3<Type>;

using Vector3f = Vector3<real>;
using Vector3i = Vector3<int32>;

template<class Type>
using Matrix3 = Eigen::Matrix3<Type>;

using Matrix3f = Matrix3<real>;
using Matrix3i = Matrix3<int32>;

template<int count>
Vector3f sqr(Vector3f vec)
{
	Vector3f result = Vector3f::Constant(1.f);
	for (size_t c = 0; c < count; c++)
	{
		for (size_t i = 0; i < 3; i++)
		{
			result[i] *= vec[i];
		}
	}
	return result;
}

template<class Type>
using Jacobian = Eigen::JacobiSVD<Type>;

using Jacobian3 = Jacobian<Matrix3f>;

template<class Type>
using BDC = Eigen::BDCSVD<Type>;

using BDC3 = BDC<Matrix3f>;