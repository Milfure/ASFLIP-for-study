#pragma once

#include <cassert>

#include "Eigen/Core"
#include "Eigen/SVD"

#include "Utilize/Util.h"

template<class Type>
class Matrix final
{
public:
	using Jacobian = Eigen::JacobiSVD<Eigen::MatrixX<Type>>;

public:
	Matrix();
	Matrix(size_t row, size_t col);
	Matrix(Type value);
	Matrix(const Eigen::MatrixX<Type>& mat);
	Matrix(Eigen::MatrixX<Type>&& mat);
	Matrix(const Matrix<Type>& other);
	Matrix(Matrix<Type>&& other) noexcept;

	~Matrix();
	Matrix<Type>& operator=(const Matrix<Type>& other);
	Matrix<Type>& operator=(Matrix<Type>&& other) noexcept;

	Matrix<Type> operator+(const Type constant);
	Matrix<Type> operator-(const Type constant);
	Matrix<Type> operator*(const Type constant);
	Matrix<Type> operator/(const Type constant);

	Matrix<Type>& operator+=(const Type constant);
	Matrix<Type>& operator-=(const Type constant);
	Matrix<Type>& operator*=(const Type constant);
	Matrix<Type>& operator/=(const Type constant);

	Matrix<Type> operator+(const Matrix<Type>& other);
	Matrix<Type> operator-(const Matrix<Type>& other);
	Matrix<Type> operator*(const Matrix<Type>& other);
	Matrix<Type> operator/(const Matrix<Type>& other);

	Matrix<Type>& operator+=(const Matrix<Type>& other);
	Matrix<Type>& operator-=(const Matrix<Type>& other);
	Matrix<Type>& operator*=(const Matrix<Type>& other);
	Matrix<Type>& operator/=(const Matrix<Type>& other);

	Matrix<Type> Diagonal();
	Matrix<Type> Product();
	Matrix<Type> Exp();
	Matrix<Type> Identity();
	Matrix<Type> Inverse();

	bool IsDiagonal() const;
	bool IsIdentity() const;
	bool IsOrthogonal() const;
	bool IsUnitary() const;
	
	void JacobianSVD(Matrix<Type>& u, Matrix<Type>& sig, Matrix<Type>& v);
	void BdcSVD(Matrix<Type>& u, Matrix<Type>& sig, Matrix<Type>& v);

	Type Norm() const;
	Matrix<Type> Normalize();
	Matrix<Type> Normalized();

	Type Determinant();
	Type Trace();

	Matrix<Type> Pow();
	Matrix<Type> Sqrt();

	void Resize(size_t row, size_t col);
	void Zero();

private:
	Eigen::MatrixX<Type> mMat;
};

template<class Type>
Matrix<Type>::Matrix()
{

}

template<class Type>
Matrix<Type>::Matrix(size_t row, size_t col)
	: mMat(row, col)
{

}

template<class Type>
Matrix<Type>::Matrix(Type value)
{

}

template<class Type>
Matrix<Type>::Matrix(const Eigen::MatrixX<Type>& mat)
	: mMat(mat)
{

}

template<class Type>
Matrix<Type>::Matrix(Eigen::MatrixX<Type>&& mat)
	: mMat(std::move(mat))
{

}

template<class Type>
Matrix<Type>::Matrix(const Matrix<Type>& other)
	: mMat(other.mMat)
{

}

template<class Type>
Matrix<Type>::Matrix(Matrix<Type>&& other) noexcept
	: mMat(std::move(other.mMat))
{

}

template<class Type>
Matrix<Type>::~Matrix()
{
	
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator=(const Matrix<Type>& other)
{
	mMat = other.mMat;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator=(Matrix<Type>&& other) noexcept
{
	mMat = std::move(other.mMat);
}

template<class Type>
Matrix<Type> Matrix<Type>::operator+(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			tmp(row, col) = mMat(row, col) + constant;
		}
	}

	return tmp;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator-(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			tmp(row, col) = mMat(row, col) + constant;
		}
	}

	return tmp;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator*(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			tmp(row, col) = mMat(row, col) * constant;
		}
	}

	return tmp;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator/(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			tmp(row, col) = mMat(row, col) / constant;
		}
	}

	return tmp;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator+=(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			mMat(row, col) += constant;
		}
	}

	return *this;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator-=(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			mMat(row, col) -= constant;
		}
	}

	return *this;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator*=(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			mMat(row, col) *= constant;
		}
	}

	return *this;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator/=(const Type constant)
{
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; row < maxRow; row++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			mMat(row, col) /= constant;
		}
	}

	return *this;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator+(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == other.mMat.cols());
	return mMat + other.mMat;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator-(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == other.mMat.cols());
	return mMat - other.mMat;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator*(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == other.mMat.cols());
	
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; i < maxRow; i++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			tmp(row, col) = mMat(row, col) * other.mMat(row, col);
		}
	}

	return tmp;
}

template<class Type>
Matrix<Type> Matrix<Type>::operator/(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == mMat.cols());
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	Matrix<Type> tmp(maxRow, maxCol);

	for (size_t row = 0; i < maxRow; i++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			tmp(row, col) = mMat(row, col) / other.mMat(row, col);
		}
	}

	return tmp;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator+=(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == other.mMat.cols());
	mMat += other.mMat;
	return *this;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator-=(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == other.mMat.cols());
	mMat -= other.mMat;
	return *this;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator*=(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == other.mMat.cols());

	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	for (size_t row = 0; i < maxRow; i++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			mMat(row, col) *= other.mMat(row, col);
		}
	}

	return *this;
}

template<class Type>
Matrix<Type>& Matrix<Type>::operator/=(const Matrix<Type>& other)
{
	assert(mMat.rows() == other.mMat.rows() && mMat.cols() == mMat.cols());
	size_t maxRow = mMat.rows();
	size_t maxCol = mMat.cols();

	for (size_t row = 0; i < maxRow; i++)
	{
		for (size_t col = 0; col < maxCol; col++)
		{
			mMat(row, col) /= other.mMat(row, col);
		}
	}

	return *this;
}