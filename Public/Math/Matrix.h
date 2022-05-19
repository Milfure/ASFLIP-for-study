#pragma once

#include "Eigen/Dense"
#include "Eigen/SVD"

template<class Type>
class Matrix final
{
private:
	using Jacobian = Eigen::JacobiSVD<Eigen::MatrixX<Type>>;

public:
	Matrix();
	Matrix(size_t row, size_t col);
	Matrix(Type value);
	Matrix(const Matrix<Type>& other);
	Matrix(Matrix<Type>&& other) noexcept;

	~Matrix();
	Matrix<Type>& operator=(const Matrix<Type>& other);
	Matrix<Type>& operator=(Matrix<Type>&& other) noexcept;

	Matrix<Type> operator+(const Matrix<Type>& other);
	Matrix<Type> operator-(const Matrix<Type>& other);
	Matrix<Type> operator*(const Matrix<Type>& other);
	Matrix<Type> operator/(const Matrix<Type>& other);

	Matrix<Type>&& operator+(const Matrix<Type>& other);
	Matrix<Type>&& operator-(const Matrix<Type>& other);
	Matrix<Type>&& operator*(const Matrix<Type>& other);
	Matrix<Type>&& operator/(const Matrix<Type>& other);

	Matrix<Type>& operator+=(const Matrix<Type>& other);
	Matrix<Type>& operator-=(const Matrix<Type>& other);
	Matrix<Type>& operator*=(const Matrix<Type>& other);
	Matrix<Type>& operator/=(const Matrix<Type>& other);

	Matrix<Type> Transpose();
	Matrix<Type> SVD();
	Matrix<Type> Inverse();
	Type Trace();

private:
	Eigen::MatrixX<Type> mMat;
};