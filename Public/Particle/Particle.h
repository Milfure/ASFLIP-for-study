#pragma once

#include "Math/Math.h"

class Particle
{
public:
	Particle();

	virtual ~Particle();

	void SetPosition(Vector3f x);
	void AddPosition(const Vector3f x);
	Vector3f GetPosition() const;

	void SetVelocity(Vector3f v);
	void AddVelocity(const Vector3f v);
	Vector3f GetVelocity() const;

	void SetAffine(Matrix3f affine);
	Matrix3f GetAffine() const;

	void SetDeformation(Matrix3f f);
	Matrix3f GetDeformation() const;

	void SetJacobian(real j);
	real GetJacobian() const;

	void SetVolume(real volume);
	real GetVolume() const;

	void SetDensity(real density);
	real GetDensity() const;

	real GetMass() const;

protected:	
	/*volume and mass*/
	real mVolume;
	real mDensity;
	real mMass;

	Vector3f Size;
	Vector3f Center; //center of size

	Vector3f V; //velocity
	Vector3f X; //position
	Matrix3f C; //affine velocity field
	Matrix3f F; //deformation gradient
	real J; //plastic deformation or Jacobian determinant
};