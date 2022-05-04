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
	class CurrentsWeighing
	{
	public:
		template<class T_ParticleArray>
		void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {};
	};

	template<GridTypes gridType>
	class FirstOrderCurrentsWeighing : public CurrentsWeighing
	{
	public:
		template<class T_ParticleArray>
		void operator()(Grid<FP, gridType>* grid, T_ParticleArray* particleArray) {
			typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
			FP3 nullCoords;
			for (int i = 0; i < grid->numCells.x; i++)
				for (int j = 0; j < grid->numCells.y; j++)
					for (int k = 0; k < grid->numCells.z; k++) {
						FP3 coords = grid->JxPosition(i, j, k);
						grid->Jx(i, j, k) = nullCoords;
						coords = grid->JyPosition(i, j, k);
						grid->Jy(i, j, k) = nullCoords;
						coords = grid->JzPosition(i, j, k);
						grid->Jz(i, j, k) = nullCoords;
					}

			for (int i = 0; i < particleArray->size(); i++) {
				ParticleProxyType particle = (*particleArray)[i];
				oneParticleWeighing(grid, &particle);
			}
		}

	private:
		template<class T_Particle>
		void oneParticleWeighing(Grid<FP, gridType>* grid, T_Particle* particle)
		{
			int i = trunc((particle->getPosition().x - grid->origin.x) / grid->steps.x);
			int j = trunc((particle->getPosition().y - grid->origin.y) / grid->steps.y);
			int k = trunc((particle->getPosition().z - grid->origin.z) / grid->steps.z);

			FP3 JBeforeWeighing = (particle->getVelocity() * particle->getCharge()) / (grid->steps.x * grid->steps.y * grid->steps.z);

			std::vector<FP3> newPosition;

			int cellNodeNumber = 2; // the number of cell nodes in each direction x,y,z

			for (int _i = 0; _i < cellNodeNumber; _i++)
				for (int _j = 0; _j < cellNodeNumber; _j++)
					for (int _k = 0; _k < cellNodeNumber; _k++) {
						FP3 cellNode((particle->getPosition().x + _i * grid->steps.x - (grid->origin.x + i * grid->steps.x)) / grid->steps.x,
							(particle->getPosition().y + _j * grid->steps.y - (grid->origin.y + j * grid->steps.y)) / grid->steps.y,
							(particle->getPosition().z + _k * grid->steps.z - (grid->origin.z + k * grid->steps.z)) / grid->steps.z);
						newPosition.push_back(cellNode);
					}

			for (int _i = 0; _i < cellNodeNumber; _i++)
				for (int _j = 0; _j < cellNodeNumber; _j++)
					for (int _k = 0; _k < cellNodeNumber; _k++) {
						// index of cell node in vector newPosition
						int newPositionIndex = _i * cellNodeNumber * cellNodeNumber + _j * cellNodeNumber + _k;
						FP3 currentDensityDeposit(JBeforeWeighing);

						if (_i) {
							currentDensityDeposit *= newPosition[newPositionIndex].x;
						}
						else {
							currentDensityDeposit *= (FP)1 - newPosition[newPositionIndex].x;
						}

						if (_j) {
							currentDensityDeposit *= newPosition[newPositionIndex].y;
						}
						else {
							currentDensityDeposit *= (FP)1 - newPosition[newPositionIndex].y;
						}

						if (_k) {
							currentDensityDeposit *= newPosition[newPositionIndex].z;
						}
						else {
							currentDensityDeposit *= (FP)1 - newPosition[newPositionIndex].z;
						}

						grid->Jx(i + _i, j + _j, k + _k) += currentDensityDeposit.x;
						grid->Jy(i + _i, j + _j, k + _k) += currentDensityDeposit.y;
						grid->Jz(i + _i, j + _j, k + _k) += currentDensityDeposit.z;
					}
		}
	};
}