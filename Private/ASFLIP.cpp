#include "ASFLIP.h"
#include "Grid/Grid.h"
#include "Particle/Particle.h"

ASFLIP::ASFLIP()
{

}

ASFLIP::~ASFLIP()
{

}

Vector3f ASFLIP::WorldSpaceToMaterialSpace(Vector3f& x, Vector3f& translation, Matrix3f& rotation)
{
	Vector3f tmp = x - translation;
	Vector3f X = rotation.transpose() * tmp;
	return x;
}

void ASFLIP::ProjectDruckerParager(Matrix3f& s, real& j)
{
	real JSe = s.determinant();
	for (size_t i = 0; i < 3; i++)
	{
		s(i, i) = std::max(1e-6f, abs(s(i, i) * j));
	}

	if (s.trace() >= 1.f)
	{
		for (size_t i = 0; i < 3; i++)
		{
			s(i, i) = 1.f;
		}
		j *= std::pow(std::max(1e-6f, JSe), VolumeRecoveryRate);
	}

	else
	{
		j = 1.f;
		real Je = std::max(1e-6f, s.trace());
		real sqrs0 = s(0, 0) * s(0, 0);
		real sqrs1 = s(1, 1) * s(1, 1);
		real sqrs2 = s(2, 2) * s(2, 2);
		real traceB2 = (sqrs0 + sqrs1 + sqrs2) / 2.f;
		real Je2 = Je * Je;
		real yieldThreshold = -MaterialFriction * Kappa0 * 0.5f * (Je2 - 1.f);
		real devB0 = sqrs0 - traceB2;
		real devB1 = sqrs1 - traceB2;
		real devB2 = sqrs2 - traceB2;
		real norm2DevB = devB0 * devB0 + devB1 * devB1 + devB2 * devB2;
		real muNormDevBBar = Mu0 * std::sqrt(norm2DevB / Je);

		if (muNormDevBBar > yieldThreshold)
		{
			real detB = sqrs0 * sqrs1 * sqrs2;
			real detDevB = devB0 * devB1 * devB2;
			real lambda2 = yieldThreshold / std::max(1e-6f, muNormDevBBar);
			real lambda1 = std::sqrt(std::max(0.f, detB - lambda2 * lambda2 * detDevB));
			s(0, 0) = std::sqrt(std::abs(lambda1 + lambda2 * devB0));
			s(1, 1) = std::sqrt(std::abs(lambda1 + lambda2 * devB1));
			s(2, 2) = std::sqrt(std::abs(lambda1 + lambda2 * devB2));
		}
	}
}

Matrix3f ASFLIP::NeoHookeanElasticity(Matrix3f& U, Matrix3f& sig)
{
	real J = sig.determinant();
	real mu_J_1_2 = Mu0 * sqrt(J);
	real JPrime = Kappa0 * 0.5f * (J * J - 1.f);
	real sqrS_1_2 = sig.trace() / 2.f;
	Matrix3f stress;
	for (size_t index = 0; index < 3; index++)
	{
		stress(index, index) = (sig(index, index) * sig(index, index) - sqrS_1_2) * mu_J_1_2;
	}
	stress = U * stress * U.transpose();
	for (size_t index = 0; index < 3; index++)
	{
		stress(index, index) += JPrime;
	}

	return stress;
}

void ASFLIP::Update()
{
	mGridManager.RestGrids();

	P2G();
	ExternalForceAndCollision();
	G2P();
}

void ASFLIP::P2G()
{
	real affineStretching = AffineStretching;
	real affineRotation = AffineRotation;
	real rc0 = (affineStretching + affineRotation) * 0.5f;
	real rc1 = (affineStretching - affineRotation) * 0.5f;

	for (Particle& particle : Particles)
	{
		Vector3f Xp = particle.GetPosition();
		Vector3f Vp = particle.GetVelocity();
		Vector3i base = (Xp * InvDx - Vector3f::Constant(0.5f)).cast<int32>();
		Vector3f Fx = Xp * InvDx - base.cast<real>();

		Vector3f w[3];
		/*w[0] = Vector3f::Constant(0.5f) * sqr<2>(Vector3f::Constant(1.5f) - Fx);
		w[1] = Vector3f::Constant(0.75f) - sqr<2>(Fx - Vector3f::Constant(1.f));
		w[2] = Vector3f::Constant(0.5f) * sqr<2>(Fx - Vector3f::Constant(0.5f));*/

		Matrix3f Fp = particle.GetDeformation();
		Matrix3f Cp = particle.GetAffine();

		Fp = (Matrix3f::Identity() + Dt * Cp) * Fp;
		real Jp = particle.GetJacobian();

		Jacobian3 svd(Fp);
		Matrix3f U, sig, V;
		U = svd.matrixU();
		sig = svd.singularValues().asDiagonal();
		V = svd.matrixV();

		ProjectDruckerParager(sig, Jp);

		Fp = U * sig * V.transpose();

		particle.SetDeformation(Fp);

		Matrix3f stress = NeoHookeanElasticity(U, sig);
		stress = (-Dt * particle.GetVolume() * 4.f * InvDx * InvDx) * stress;

		Matrix3f affineWithoutStress = particle.GetMass() * (Cp * rc0 + Cp.transpose() * rc1);
		Matrix3f affine = stress + affineWithoutStress;

		Grid*** grids = mGridManager.GetGrids();
		for (size_t x = 0; x < 3; x++)
		{
			for (size_t y = 0; y < 3; y++)
			{
				for (size_t z = 0; z < z; z++)
				{
					Grid& grid = grids[base[0] + x][base[1] + y][base[2] + z];
					Vector3i offset(x, y, z);
					Vector3f dpos = (offset.cast<real>() - Fx) * Dx;
					real weight = w[x][0] + w[y][1] + w[z][2];
					grid.AddCurrentVelocity(Vector3f(weight * (particle.GetMass() * Vp + affine * dpos)));
					grid.AddPrevVelocity(Vector3f(weight * (particle.GetMass() * Vp + affineWithoutStress * dpos)));
					grid.AddMass(weight * particle.GetMass());
				}
			}
		}
	}
}

void ASFLIP::ExternalForceAndCollision()
{
	Grid*** grids = mGridManager.GetGrids();
	uint16 girdNum = mGridManager.GetGridNum();
	for (size_t x = 0; x < GridNum; x++)
	{
		for (size_t y = 0; y < GridNum; y++)
		{
			for (size_t z = 0; z < GridNum; z++)
			{
				Grid& grid = grids[x][y][z];
				real gMass = grid.GetMass();
				if (gMass > 0)
				{
					Vector3f gV = grid.GetCurrentVelocity();
					gV *= 1.f / gMass;
					gV += Dt * Gravity;
					grid.SetPrevVelocity(grid.GetPrevVelocity() * (1.f / gMass));

					if (x < 3 && gV[0] < 0)
					{
						gV[0] = 0.f;
					}
				}
			}
		}
	}
}

void ASFLIP::G2P()
{
	/*for(Particle& particle : Particles)
	{
		Vector3f Xp = particle.GetPosition();
		Vector3i base = (Xp * InvDx - Vector3f::Constant(0.5f)).cast<int32>();
		Vector3f Fx = Xp * InvDx - base.cast<real>();

		Vector3f w[3];
		w[0] = Vector3f::Constant(0.5f) * sqr<2>(Vector3f::Constant(1.5f) - Fx);
		w[1] = Vector3f::Constant(0.75f) - sqr<2>(Fx - Vector3f::Constant(1.f));
		w[2] = Vector3f::Constant(0.5f) * sqr<2>(Fx - Vector3f::Constant(0.5f));

		Vector3f newV = Vector3f::Constant(0.f);
		Matrix3f newC = Matrix3f::Constant(0.f);
		
		Grid*** grids = mGridManager.GetGrids();
		for (size_t x = 0; x < 3; x++)
		{
			for (size_t y = 0; y < 3; y++)
			{
				for (size_t z = 0; z < 3; z++)
				{
					Grid& grid = grids[base[0] + x][base[1] + y][base[2] + z];
					Vector3f dpos = Vector3i(x, y, z).cast<real>() - Fx;
					Vector3f gV = grid.GetCurrentVelocity();
					real weight = w[x][0] + w[y][1] + w[z][2];
					newV += weight * gV;
					newC += 4.f * InvDx * weight * gV.cross(dpos);
				}
			}
		}
		if (VelocityAdjustment > 0.f)
		{
			Vector3f Vp = particle.GetVelocity();
			real posAdj = PositionAdjustmentMax;

			if (posAdj > 0.f && ParticleCollision > 0.f)
			{

			}

			if (PositionAdjustmentMin < posAdj)
			{
				real logdJ = newC.trace() * Dt;
				real J = particle.GetDeformation().determinant();
				if (log(std::max(1e-15f, J)) + logdJ < -0.001f)
				{
					posAdj = PositionAdjustmentMin;
				}
			}
			Vector3f oldV = Vector3f::Constant(0.f);

			for (size_t x = 0; x < 3; x++)
			{
				for (size_t y = 0; y < 3; y++)
				{
					for (size_t z = 0; z < 3; z++)
					{
						Grid& grid = grids[base[0] + x][base[1] + y][base[2] + z];
						Vector3f gV0 = grid.GetPrevVelocity();
						real weight = w[x][0] * w[y][1] * w[z][2];
						oldV += weight * gV0;
					}
				}
			}

			Vector3f diffV = Vp - oldV;
			particle.SetVelocity(newV + VelocityAdjustment * diffV);
			particle.SetPosition(Xp + (newV + posAdj * VelocityAdjustment * diffV) * Dt);
		}
		else
		{
			particle.SetVelocity(newV);
			particle.SetPosition(Xp + newV * Dt);
		}
		particle.SetAffine(newC);
	}*/
}

