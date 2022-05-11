#include "Particle/Particle.h"

Particle::Particle()
{

}

Particle::~Particle()
{

}

void Particle::SetPosition(Vector3f x)
{
	X = x;
}

void Particle::AddPosition(const Vector3f x)
{
	X += x;
}

Vector3f Particle::GetPosition() const
{
	return X;
}

void Particle::SetVelocity(Vector3f v)
{
	V = v;
}

void Particle::AddVelocity(const Vector3f v)
{
	V += v;
}

Vector3f Particle::GetVelocity() const
{
	return V;
}

void Particle::SetAffine(Matrix3f affine)
{
	C = affine;
}

Matrix3f Particle::GetAffine() const
{
	return C;
}

void Particle::SetDeformation(Matrix3f f)
{
	F = f;
}

Matrix3f Particle::GetDeformation() const
{
	return F;
}

void Particle::SetJacobian(real j)
{
	J = j;
}

real Particle::GetJacobian() const
{
	return J;
}

void Particle::SetVolume(real volume)
{
	mVolume = volume;
	mDensity = mMass / mVolume;
}

real Particle::GetVolume() const
{
	return mVolume;
}

void Particle::SetDensity(real density)
{
	mDensity = density;
	mVolume = mMass / mDensity;
}

real Particle::GetDensity() const
{
	return mDensity;
}

real Particle::GetMass() const
{
	return mMass;
}