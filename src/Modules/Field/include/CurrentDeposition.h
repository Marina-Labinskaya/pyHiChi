#pragma once

#include "Constants.h"
#include "Grid.h"
#include "FormFactor.h"
#include "ParticleArray.h"
#include "Particle.h"
#include <iostream>
#include <vector>

namespace pfc
{
    template<class TGrid, class DerivedClass, class LocalDerivedClass>
    class CurrentDeposition
    {
    public:
        enum class ZeroizeJ {
            NOT_USE_ZEROIZEJ, USE_ZEROIZEJ
        };

        CurrentDeposition(double _dt) : halfDt(0.5 * _dt) {}

        template<class T_Particle>
        void operator()(TGrid* grid, const T_Particle& particle,
            CurrentDeposition::ZeroizeJ UsingZeroizeJ = CurrentDeposition::ZeroizeJ::NOT_USE_ZEROIZEJ) {
            if (UsingZeroizeJ == ZeroizeJ::USE_ZEROIZEJ)
                grid->zeroizeJ();
            static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle);
        }

//         template<class T_ParticleArray>
//         void operator()(TGrid* grid, T_ParticleArray* particleArray) {
//             typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
//             grid->zeroizeJ();
//             LocalDerivedClass Depositor(Int3(0,0,0),
//                 static_cast<DerivedClass*>(this), grid);
// #pragma omp declare reduction (add_currents : LocalDerivedClass : omp_out.addCurrents()) \
//     initializer(omp_priv = omp_orig)
// #pragma omp parallel for reduction(add_currents: Depositor)
//             for (int i = 0; i < particleArray->size(); i++) {
//                 ParticleProxyType particle = (*particleArray)[i];
//                 Int3 baseGridIdx = grid->getBaseIndex(
//                     particle.getPosition() - (particle.getVelocity() * halfDt));
//                 if (Depositor.baseGridIdx != baseGridIdx) {
//                     Depositor.addCurrents();
//                     Depositor.setNewGridCell(baseGridIdx);
//                 }
//                 static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle, Depositor);
//             }
//         }

        template<class T_ParticleArray>
        void operator()(TGrid* grid, T_ParticleArray* particleArray) {
            typedef typename T_ParticleArray::ParticleProxyType ParticleProxyType;
            grid->zeroizeJ();
            int num_threads = omp_get_max_threads();
            std::vector<LocalDerivedClass> Depositor(num_threads, LocalDerivedClass(Int3(0,0,0),
               static_cast<DerivedClass*>(this), grid));
#pragma omp parallel for
            for (int i = 0; i < particleArray->size(); i++) {
                int thread_num = omp_get_thread_num();
                ParticleProxyType particle = (*particleArray)[i];
                Int3 baseGridIdx = grid->getBaseIndex(
                    particle.getPosition() - (particle.getVelocity() * halfDt));
                if (Depositor[thread_num].baseGridIdx != baseGridIdx) {
                    Depositor[thread_num].addCurrents();
                    Depositor[thread_num].setNewGridCell(baseGridIdx);
                }
                static_cast<DerivedClass*>(this)->depositOneParticle(grid, &particle, Depositor[thread_num]);
            }
#pragma omp parallel for
           for (int i = 0; i < num_threads; ++i)
               Depositor[i].addCurrents();
        }


        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle) {}

        const double halfDt;
    };

    template<class TGrid, class TCurrentDeposition, int BlockSize>
    class LocalDeposition
    {
    public:

        static const int blockSize = BlockSize;

        LocalDeposition(const Int3& _baseGridIdx, TCurrentDeposition* currentDeposition, TGrid* _grid) :
            baseGridIdx(_baseGridIdx), grid(_grid), halfDt(currentDeposition->halfDt)
        {            
            blockOffset = (blockSize + 1) / 2;
            startGridIdx = baseGridIdx - Int3(blockOffset, blockOffset, blockOffset);

            for (int i = 0; i < this->blockSize; ++i)
                for (int j = 0; j < this->blockSize; ++j)
                    for (int k = 0; k < this->blockSize; ++k)
                    {
                        Jx[i][j][k] = 0;
                        Jy[i][j][k] = 0;
                        Jz[i][j][k] = 0;
                    }
        }

        void setNewGridCell(const Int3& _baseGridIdx) {
            baseGridIdx = _baseGridIdx;
            startGridIdx = baseGridIdx - Int3(blockOffset, blockOffset, blockOffset);
        }

        void addCurrents()
        {
            for (int i = 0; i < blockSize; ++i)
                for (int j = 0; j < blockSize; ++j)
                    for (int k = 0; k < blockSize; ++k)
                    {
                        Int3 temp = startGridIdx + Int3(i, j, k);
                        Int3 gridIdx = Int3(temp.x % grid->numCells.x, temp.y % grid->numCells.y, temp.z % grid->numCells.z);
                        #pragma omp atomic
                        grid->Jx(gridIdx) += Jx[i][j][k];
                        #pragma omp atomic
                        grid->Jy(gridIdx) += Jy[i][j][k];
                        #pragma omp atomic
                        grid->Jz(gridIdx) += Jz[i][j][k];
                    }

            for (int i = 0; i < this->blockSize; ++i)
                for (int j = 0; j < this->blockSize; ++j)
                    for (int k = 0; k < this->blockSize; ++k)
                    {
                        Jx[i][j][k] = 0;
                        Jy[i][j][k] = 0;
                        Jz[i][j][k] = 0;
                    }
        }

        template<class TParticle>
        void depositCurrent(const TParticle& particle) {}
        Int3 baseGridIdx;
    protected:

        FP Jx[blockSize][blockSize][blockSize];
        FP Jy[blockSize][blockSize][blockSize];
        FP Jz[blockSize][blockSize][blockSize];
        TGrid* grid;
        Int3 startGridIdx;
        int blockOffset;
        const double halfDt;
    };

    template<class TGrid>
    class CurrentDepositionCIC;

    template<class TGrid>
    class LocalDepositionCIC : public LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 4>
    {
    public:
        LocalDepositionCIC() : LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 4>() {}
        LocalDepositionCIC(const LocalDepositionCIC& deposition):
            LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 4>(
            static_cast<LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 4>>(deposition)) {} 
        LocalDepositionCIC(const Int3& blockIdx, CurrentDepositionCIC<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionCIC<TGrid>, 4>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 particlePosition = particle->getPosition() - (particle->getVelocity() * this->halfDt);
            FP3 current = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) /
                this->grid->steps.volume();

            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            idxJx = this->grid->getIndexJx(particlePosition) - this->startGridIdx;
            internalCoordsJx = this->grid->getInternalCoordsJx(particlePosition);

            idxJy = this->grid->getIndexJy(particlePosition) - this->startGridIdx;
            internalCoordsJy = this->grid->getInternalCoordsJy(particlePosition);

            idxJz = this->grid->getIndexJz(particlePosition) - this->startGridIdx;
            internalCoordsJz = this->grid->getInternalCoordsJz(particlePosition);

            FormFactorCIC formFactorJx, formFactorJy, formFactorJz;
            
            formFactorJx(internalCoordsJx);
            formFactorJy(internalCoordsJy);
            formFactorJz(internalCoordsJz);
            depositComponent(idxJx, current.x, formFactorJx, this->Jx);
            depositComponent(idxJy, current.y, formFactorJy, this->Jy);
            depositComponent(idxJz, current.z, formFactorJz, this->Jz);
        }

    private:
        void depositComponent(const Int3 & idx, const FP & value, FormFactorCIC& formFactor,
            FP current[4][4][4])
        {
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j++) {
                    for (int k = 0; k <= 1; k++) {
                        current[idx.x + i][idx.y + j][idx.z + k] += 
                            (formFactor.c[0][i] *
                             formFactor.c[1][j] *
                             formFactor.c[2][k]
                            ) * value;
                    }
                }
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionCIC : public CurrentDeposition<TGrid, CurrentDepositionCIC<TGrid>,
        LocalDepositionCIC<TGrid>>
    {
    public:
        CurrentDepositionCIC(double _dt) : CurrentDeposition<TGrid, CurrentDepositionCIC<TGrid>,
            LocalDepositionCIC<TGrid>>(_dt) {}

        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle, LocalDepositionCIC<TGrid>& Depositor)
        {
            //FP3 particlePosition = particle->getPosition() - (particle->getVelocity() * this->dt / 2.0);
            //Int3 baseGridIdx = grid->getBaseIndex(particlePosition);
            //LocalDepositionCIC<TGrid> Depositor(baseGridIdx, this, grid);
            Depositor.depositCurrent(particle);
        }
    };

    template<class TGrid>
    class CurrentDepositionTSC;

    template<class TGrid>
    class LocalDepositionTSC : public LocalDeposition<TGrid, CurrentDepositionTSC<TGrid>, 5>
    {
    public:

        LocalDepositionTSC(const Int3 & blockIdx, CurrentDepositionTSC<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionTSC<TGrid>, 5>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 particlePosition = particle->getPosition() - (particle->getVelocity() * this->halfDt);
            FP3 current = (particle->getVelocity() * particle->getCharge() * particle->getWeight()) /
                this->grid->steps.volume();

            Int3 idxJx, idxJy, idxJz;
            FP3 internalCoordsJx, internalCoordsJy, internalCoordsJz;
            idxJx = this->grid->getClosestIndexJx(particlePosition) - this->startGridIdx;
            internalCoordsJx = this->grid->getClosestInternalCoordsJx(particlePosition);

            idxJy = this->grid->getClosestIndexJy(particlePosition) - this->startGridIdx;
            internalCoordsJy = this->grid->getClosestInternalCoordsJy(particlePosition);

            idxJz = this->grid->getClosestIndexJz(particlePosition) - this->startGridIdx;
            internalCoordsJz = this->grid->getClosestInternalCoordsJz(particlePosition);

            FormFactorTSC formFactorJx, formFactorJy, formFactorJz;
            
            formFactorJx(internalCoordsJx);
            formFactorJy(internalCoordsJy);
            formFactorJz(internalCoordsJz);
            depositComponent(idxJx, current.x, formFactorJx, this->Jx);
            depositComponent(idxJy, current.y, formFactorJy, this->Jy);
            depositComponent(idxJz, current.z, formFactorJz, this->Jz);
        }

    private:

        void depositComponent(const Int3 & idx,
            const FP & value, FormFactorTSC& formFactor, FP current[5][5][5])
        {
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    for (int k = -1; k <= 1; k++) {
                        current[idx.x + i][idx.y + j][idx.z + k] += 
                            ( formFactor.c[0][i + 1] 
                            * formFactor.c[1][j + 1]
                            * formFactor.c[2][k + 1])
                            * value;
                    }
                }
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionTSC : public CurrentDeposition<TGrid, CurrentDepositionTSC<TGrid>,
        LocalDepositionTSC<TGrid>>
    {
    public:
        CurrentDepositionTSC(double _dt) : CurrentDeposition<TGrid, CurrentDepositionTSC<TGrid>,
            LocalDepositionTSC<TGrid>>(_dt) {}

        template<class T_Particle>
        void depositOneParticle(TGrid* grid, T_Particle* particle, LocalDepositionTSC<TGrid>& Depositor)
        {
            
            //FP3 particlePosition = particle->getPosition() - (particle->getVelocity() * this->dt / 2.0);
            //Int3 baseGridIdx = grid->getBaseIndex(particlePosition);
            //LocalDepositionTSC<TGrid> Depositor(baseGridIdx, this, grid);
            Depositor.depositCurrent(particle);
        }
    };
}
