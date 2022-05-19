#include <iostream>

#include "Math/Math.h"

int main(void)
{
	Eigen::MatrixX<real> mat(3,3);

	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			mat(i, j) = (real)3 * i + j;
		}
	}

	std::cout << mat << std::endl;

	Jacobian<Eigen::MatrixX<real>> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

	std::cout << svd.computeU() << " " << svd.computeV() << std::endl;

	return 0;
}