#pragma once

#include "Grid.h"
#include "ParticleArray.h"
#include "Particle.h"

#include <cmath>
#include <vector>

using namespace std;

namespace pfc 
{
	template<GridTypes gridType>
	class CurrentDeposition
	{
	public:
		template<class T_Particle>
		void operator()(Grid<FP, gridType>* grid, T_Particle* particle) {};

		template<class T_ParticleArray>
		void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {};
	};

	template<GridTypes gridType>
	class FirstOrderCurrentDeposition : public CurrentDeposition
	{
	public:
		template<class T_ParticleArray>
		void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {
			grid->zeroizeJ();
			oneParticleDeposition(grid, &particle);
		}


		template<class T_ParticleArray>
		void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {
			typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
			
			grid->zeroizeJ();

			for (int i = 0; i < particleArray->size(); i++) {
				ParticleProxyType particle = (*particleArray)[i];
				oneParticleDeposition(grid, &particle);
			}
		}

	private:
		template<class T_Particle>
		void oneParticleDeposition(Grid<FP, gridType>* grid, T_Particle* particle)
		{
			FP3 particlePosition = particle->getPosition();
			Int3 idx;
			FP3 internalCoords;

			FP3 JBeforeDeposition = (particle->getVelocity() * particle->getCharge()) / (grid->steps.x * grid->steps.y * grid->steps.z);

			// i, j, k
			grid->getGridCoords(particlePosition, FP3(0, 0, 0), idx, internalCoords); // don't know where find J shift

			Jx(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * JBeforeDeposition.z;

			// i + 1, j, k 
			grid->getGridCoords(particlePosition + FP3((FP)1 * grid->steps.x, 0, 0), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * ((FP)1 - internalCoords.z) * JBeforeDeposition.z;

			// i, j + 1, k
			grid->getGridCoords(particlePosition + FP3(0, (FP)1 * grid->steps.y, 0), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * ((FP)1 - internalCoords.z) * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * ((FP)1 - internalCoords.z) * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * ((FP)1 - internalCoords.z) * JBeforeDeposition.z;

			// i, j, k + 1
			grid->getGridCoords(particlePosition + FP3(0, 0, (FP)1 * grid->steps.z), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * internalCoords.z * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * internalCoords.z * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * ((FP)1 - internalCoords.y) * internalCoords.z * JBeforeDeposition.z;

			// i + 1, j + 1, k
			grid->getGridCoords(particlePosition + FP3((FP)1 * grid->steps.x, (FP)1 * grid->steps.y, 0), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += internalCoords.x * internalCoords.y * ((FP)1 - internalCoords.z) * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += internalCoords.x * internalCoords.y * ((FP)1 - internalCoords.z) * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += internalCoords.x * internalCoords.y * ((FP)1 - internalCoords.z) * JBeforeDeposition.z;

			// i + 1, j, k + 1
			grid->getGridCoords(particlePosition + FP3((FP)1 * grid->steps.x, 0, (FP)1 * grid->steps.z), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * internalCoords.z * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * internalCoords.z * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += internalCoords.x * ((FP)1 - internalCoords.y) * internalCoords.z * JBeforeDeposition.z;

			// i, j + 1, k + 1
			grid->getGridCoords(particlePosition + FP3(0, (FP)1 * grid->steps.y, (FP)1 * grid->steps.z), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * internalCoords.z * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * internalCoords.z * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += ((FP)1 - internalCoords.x) * internalCoords.y * internalCoords.z * JBeforeDeposition.z;

			// i + 1, j + 1, k + 1 
			grid->getGridCoords(particlePosition + FP3((FP)1 * grid->steps.x, (FP)1 * grid->steps.y, (FP)1 * grid->steps.z), FP3(0, 0, 0), idx, internalCoords);

			Jx(idx.x, idx.y, idx.z) += internalCoords.x * internalCoords.y * internalCoords.z * JBeforeDeposition.x;
			Jy(idx.x, idx.y, idx.z) += internalCoords.x * internalCoords.y * internalCoords.z * JBeforeDeposition.y;
			Jz(idx.x, idx.y, idx.z) += internalCoords.x * internalCoords.y * internalCoords.z * JBeforeDeposition.z;
		}
	};
}