#pragma once

#include <vector>

#include "Math/Math.h"
#include "Utilize/Util.h"
#include "Grid/GridManager.h"

class Particle;

class ASFLIP final
{
public:
	ASFLIP();

	~ASFLIP();

	Vector3f WorldSpaceToMaterialSpace(Vector3f& x, Vector3f& translation, Matrix3f& rotation);

	void ProjectDruckerParager(Matrix3f& sig, real& j);
	Matrix3f NeoHookeanElasticity(Matrix3f& U, Matrix3f& sig);

	void Update();

	void P2G();
	void ExternalForceAndCollision();
	void G2P();

private:
	/*For calculate ASFLIP*/
	real VelocityAdjustment = 0.99f;
	real PositionAdjustmentMin = 0.f;
	real PositionAdjustmentMax = 1.f;
	real AffineStretching = 1.f;
	real AffineRotation = 1.f;
	real ParticleCollision = 1.f;
	real AdvParams[6];

	/*simulating environment setting*/
	real Quality = 1.f;
	uint16 GridSize = 96;
	uint16 GridNum = 96;
	real Dx = 1.f / GridSize; //Delta x
	real InvDx = (real)GridSize; // Inverse delta x

	/*delta time setting*/
	real FrameDt = 4e-3f;
	real Dt = 1e-4f / Quality;

	/*Young's modulus and Poisson's ratio*/
	real E = 5e5f;
	real Nu = 0.3f;

	/*Bulk modulus and Shear modulus*/
	real Kappa0 = E / (3 * (1 - Nu * 2));
	real Mu0 = E / (2 * (1 + Nu));

	/*Plasticity parameter*/
	real FrictionAngle = 40.f;
	real SinPi = sin(FrictionAngle / 180.f * 3.141592653f);
	real MaterialFriction = 1.633f * SinPi / (3.f - SinPi);
	real VolumeRecoveryRate = 0.5f;

	Vector3f Gravity;

	std::vector<Particle> Particles;
	GridManager mGridManager;
};