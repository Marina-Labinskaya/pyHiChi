#pragma once

#include "Dimension.h"
#include "Grid.h"
#include "Particle.h"
#include "Vectors.h"
#include "VectorsProxy.h"
#include <random>

namespace pfc {
    class ParticleGenerator {
        typedef typename VectorTypeHelper<Three, Real>::Type PositionType;
        typedef typename VectorTypeHelper<Three, Real>::Type MomentumType;
        typedef ParticleTypes TypeIndexType;
        typedef Real WeightType;
     public:

        template<class T_ParticleArray>
        void operator()(T_ParticleArray* particleArray, double(*f)(double, double, double), double T0,
            const FP3& minCoords, const FP3& maxCoords, int counter = 1, const MomentumType& momentum = FP3(0, 0, 0), WeightType weight = 1,
            TypeIndexType typeIndex = ParticleTypes::Electron);
    private:
        FP3 GetParticleRandomPosition(double(*f)(double, double, double), const FP3& minCoord, const FP3& maxCoord, int counter) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> rd_posX(minCoord.x, maxCoord.x);
            std::uniform_real_distribution<> rd_posY(minCoord.y, maxCoord.y);
            std::uniform_real_distribution<> rd_posZ(minCoord.z, maxCoord.z);
            std::uniform_real_distribution<> u;
            bool pos_correct = false;// , posY_correct = false, posZ_correct = false;
            FP3 pos;
            while (!pos_correct) {
                pos.x = rd_posX(gen);
                pos.y = rd_posY(gen);
                pos.z = rd_posZ(gen);
                if ((FP)u(gen) < (FP)(f(pos.x, pos.y, pos.z) * (maxCoord.x - minCoord.x)*(maxCoord.y - minCoord.y)* (maxCoord.z - minCoord.z) / counter))
                    pos_correct = true;
            }

            /*while (!posY_correct) {
                pos.y = rd_posY(gen);
                if ((FP)u(gen) < (FP)(f(pos.x, pos.y, pos.z) * (maxCoord.y - minCoord.y) / counter))
                    posY_correct = true;
            }

            while (!posZ_correct) {
                pos.z = rd_posZ(gen);
                if ((FP)u(gen) < (FP)(f(pos.x, pos.y, pos.z) * (maxCoord.z - minCoord.z) / counter))
                    posZ_correct = true;
            }*/

            return pos;
        }
        MomentumType GetParticleRandomMomentum(const MomentumType& momentum, double T0) {
            double alpha = 1.5 / (std::pow(T0 / (constants::electronMass * constants::c * constants::c) + 1, 2) - 1);
            double sigma = std::sqrt(0.5 / alpha) * constants::electronMass * constants::c;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> rd_momentum(0.0, 1.0);
            MomentumType newMomentum(momentum.x + sigma*rd_momentum(gen), momentum.y + sigma*rd_momentum(gen), momentum.z + sigma*rd_momentum(gen));
            return newMomentum;
        }
    };

    template<class T_ParticleArray>
    void ParticleGenerator::operator()(T_ParticleArray* particleArray, double(*f)(double, double, double), double T0,
        const FP3& minCoords, const FP3& maxCoords, int counter, const MomentumType& momentum, WeightType weight,
        TypeIndexType typeIndex) {

        for (int i = 0; i < counter; ++i) {
            PositionType ParticlePosition(GetParticleRandomPosition(f, minCoords, maxCoords, counter));
            Particle3d NewParticle(ParticlePosition, GetParticleRandomMomentum(momentum, T0), weight, typeIndex);
            particleArray->pushBack(NewParticle);
        }
    }
}
