#pragma once

#include "Math/Math.h"

class Grid
{
public:
	Grid();
	Grid(uint16 size);

	virtual ~Grid();

	void SetSize(uint16 size);
	uint16 GetSize() const;

	void SetCurrentVelocity(Vector3f v);
	void AddCurrentVelocity(Vector3f v);
	Vector3f GetCurrentVelocity() const;

	void SetPrevVelocity(Vector3f v);
	void AddPrevVelocity(Vector3f v);
	Vector3f GetPrevVelocity() const;

	void SetMass(real mass);
	void AddMass(real mass);
	real GetMass() const;

	void Reset();

private:
	uint16 mSize = 96;

	Vector3f mV; //velocity
	Vector3f mV0; // previous velocity
	real mMass; //mass
};