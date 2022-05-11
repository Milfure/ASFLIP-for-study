#include "Grid/Grid.h"

Grid::Grid()
	: mV(Vector3f::Zero())
	, mV0(Vector3f::Zero())
	, mMass(0.f)
{

}

Grid::~Grid()
{

}

void Grid::SetSize(uint16 size)
{
	mSize = size;
}

uint16 Grid::GetSize() const
{
	return mSize;
}

void Grid::SetCurrentVelocity(Vector3f v)
{
	mV = v;
}

void Grid::AddCurrentVelocity(Vector3f v)
{
	mV += v;
}

Vector3f Grid::GetCurrentVelocity() const
{
	return mV;
}

void Grid::SetPrevVelocity(Vector3f v)
{
	mV0 = v;
}

void Grid::AddPrevVelocity(Vector3f v)
{
	mV0 += v;
}

Vector3f Grid::GetPrevVelocity() const
{
	return mV0;
}

void Grid::SetMass(real mass)
{
	mMass = mass;	
}

void Grid::AddMass(real mass)
{
	mMass += mass;
}

real Grid::GetMass() const
{
	return mMass;
}

void Grid::Reset()
{
	for (size_t i = 0; i < 3; i++)
	{
		mV[i] = 0.f;
		mV0[i] = 0.f;
	}
	mMass = 0.f;
}

