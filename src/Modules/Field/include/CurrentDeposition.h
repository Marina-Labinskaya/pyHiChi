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
                    particle.getPosition());// - (particle.getVelocity() * halfDt));
                if (Depositor[thread_num].getBaseGridIdx() != baseGridIdx) {
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
        void depositOneParticle(TGrid* grid, T_Particle* particle, LocalDerivedClass& Depositor)
        {
            Depositor.depositCurrent(particle);
        }

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
            blockOffset = (blockSize) / 2;
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
                        #pragma omp atomic
                        grid->Jx(remainder(startGridIdx + Int3(i, j, k), grid->numCells)) += Jx[i][j][k];
                        #pragma omp atomic
                        grid->Jy(remainder(startGridIdx + Int3(i, j, k), grid->numCells)) += Jy[i][j][k];
                        #pragma omp atomic
                        grid->Jz(remainder(startGridIdx + Int3(i, j, k), grid->numCells)) += Jz[i][j][k];
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
        Int3 getBaseGridIdx() { return baseGridIdx; }
    protected:

        FP Jx[blockSize][blockSize][blockSize];
        FP Jy[blockSize][blockSize][blockSize];
        FP Jz[blockSize][blockSize][blockSize];
        TGrid* grid;
        Int3 startGridIdx;
        int blockOffset;
        const double halfDt;
        Int3 baseGridIdx;
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
    };

    template<class TGrid>
    class CurrentDepositionVB;

    template<class TGrid>
    class LocalDepositionVB : public LocalDeposition<TGrid, CurrentDepositionVB<TGrid>, 3>
    {
    public:

        LocalDepositionVB(
            const Int3 & blockIdx, CurrentDepositionVB<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionVB<TGrid>, 3>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            //FP3 velocity = particle->getVelocity();
            //FP3 oldCoords = particle->getPosition() - (particle->getVelocity() * (2 * this->halfDt));
            //FP3 coords = particle->getPosition();
            //Int3 oldLocalOrigin = this->grid->getClosestBaseIndex(oldCoords);
            //Int3 localOrigin = this->grid->getClosestBaseIndex(particle->getPosition());
            FP3 origin = this->grid->origin + this->startGridIdx * this->grid->steps;
            FP charge = particle->getCharge() * particle->getWeight();// / this->grid->steps.volume() / (2 * this->halfDt);
            FP3 coords = particle->getPosition();
            FP3 velocity = particle->getVelocity();
            FP3 oldCoords = coords - velocity * (2 * this->halfDt);
            FP3 middleCellOrigin = origin + 0.5 * this->grid->steps;
            Int3 oldLocalOrigin = truncate((oldCoords - middleCellOrigin) / this->grid->steps);
            Int3 localOrigin = truncate((coords - middleCellOrigin) * this->grid->steps);
            depositCurrentRecursive(coords, oldCoords, localOrigin, oldLocalOrigin, velocity, charge);
        }

    private:
        void depositCurrentRecursive(FP3 coords, FP3 oldCoords, Int3 localOrigin, Int3 oldLocalOrigin, FP3 velocity,
            FP charge)
        {
            FP3 origin = this->grid->origin + this->startGridIdx * this->grid->steps;
            /* Find the axis that is first intersected.
            Take max time as it is measured for backward movement. */
            FP maxTimeToIntersection = 0;
            int earliestIntersectionAxis = -1;
            for (int d = 0; d < 3; d++)
            {
                if (localOrigin[d] != oldLocalOrigin[d])
                {
                    FP newOldCoords = origin[d] + 0.5 * this->grid->steps[d] +
                        std::max(oldLocalOrigin[d], localOrigin[d]) * this->grid->steps[d];
                    FP timeToIntersection = (coords[d] - newOldCoords) / velocity[d];
                    if (timeToIntersection > maxTimeToIntersection)
                    {
                        maxTimeToIntersection = timeToIntersection;
                        earliestIntersectionAxis = d;
                    }
                }
            }
            if (earliestIntersectionAxis >= 0)
            {
                /* Split particle along axis dim0; dim1 and dim2 are two other axes
                (0 = x, 1 = y, 2 = z) */
                const int dim0 = earliestIntersectionAxis;
                const int dim1 = (dim0 + 1) % 3;
                const int dim2 = (dim0 + 2) % 3;
                FP3 newOldCoords;
                newOldCoords[dim0] = (origin + 0.5 * this->grid->steps)[dim0] +
                    std::max(oldLocalOrigin[dim0], localOrigin[dim0]) * this->grid->steps[dim0];
                newOldCoords[dim1] = coords[dim1] - velocity[dim1] * maxTimeToIntersection;
                newOldCoords[dim2] = coords[dim2] - velocity[dim2] * maxTimeToIntersection;
                Int3 newOldLocalOrigin = oldLocalOrigin;
                newOldLocalOrigin[dim0] = localOrigin[dim0];
                depositCurrentRecursive(coords, newOldCoords, localOrigin, newOldLocalOrigin, velocity, charge);
                coords = newOldCoords;
                localOrigin = oldLocalOrigin;
            }

            FP3 delta = (coords - oldCoords) / this->grid->steps;
            FP3 midway = (0.5 * (coords + oldCoords) - origin) / this->grid->steps -
                FP3(localOrigin) - FP3(0.5, 0.5, 0.5);

            FP3 iMidway = FP3(1.0, 1.0, 1.0) - midway;
            FP3 csteps = charge * this->grid->steps;
            FP3 delcs = delta * csteps;
            FP3 csdelmid, csdelimid;
            csdelmid[0] = delcs.x * midway.z;
            csdelmid[1] = delcs.y * midway.z;
            csdelmid[2] = delcs.z * midway.x;
            csdelimid[0] = delcs.x * iMidway.z;
            csdelimid[1] = delcs.y * iMidway.z;
            csdelimid[2] = delcs.z * iMidway.x;
            FP3 csdel = csteps * delta.x * delta.y * delta.z / FP3(12, 12, 12);
            Int3 idx = localOrigin;
            Int3 iIdx = idx + Int3(1, 1, 1);        

            this->Jx[iIdx.x][iIdx.y][iIdx.z] += csdelmid[0] * midway.y + csdel.x;
            this->Jx[iIdx.x][idx.y][iIdx.z] += csdelmid[0] * iMidway.y - csdel.x;
            this->Jx[iIdx.x][iIdx.y][idx.z] += csdelimid[0] * midway.y - csdel.x;
            this->Jx[iIdx.x][idx.y][idx.z] += csdelimid[0] * iMidway.y + csdel.x;

            this->Jy[iIdx.x][iIdx.y][iIdx.z] += csdelmid[1] * midway.x + csdel.y;
            this->Jy[idx.x][iIdx.y][iIdx.z] += csdelmid[1] * iMidway.x - csdel.y;
            this->Jy[iIdx.x][iIdx.y][idx.z] += csdelimid[1] * midway.x - csdel.y;
            this->Jy[idx.x][iIdx.y][idx.z] += csdelimid[1] * iMidway.x + csdel.y;

            this->Jz[iIdx.x][iIdx.y][iIdx.z] += csdelmid[2] * midway.y + csdel.z;
            this->Jz[iIdx.x][idx.y][iIdx.z] += csdelmid[2] * iMidway.y - csdel.z;
            this->Jz[idx.x][iIdx.y][iIdx.z] += csdelimid[2] * midway.y - csdel.z;
            this->Jz[idx.x][idx.y][iIdx.z] += csdelimid[2] * iMidway.y + csdel.z;
        }
    };

    template<class TGrid>
    class CurrentDepositionVB :
        public CurrentDeposition<TGrid, CurrentDepositionVB<TGrid>, LocalDepositionVB<TGrid>>
    {
    public:
        CurrentDepositionVB(double _dt) : 
            CurrentDeposition<TGrid, CurrentDepositionVB<TGrid>, LocalDepositionVB<TGrid>>(_dt) {}
    };

    template<class TGrid>
    class CurrentDepositionZ1;

    template<class TGrid>
    class LocalDepositionZ1 : public LocalDeposition<TGrid, CurrentDepositionZ1<TGrid>, 3>
    {
    public:

        LocalDepositionZ1(
            const Int3 & blockIdx, CurrentDepositionZ1<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionZ1<TGrid>, 3>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 origin = this->grid->origin + this->startGridIdx * this->grid->steps;
            FP3 coords = particle->getPosition();
            FP3 oldCoords = particle->getPosition() - particle->getVelocity() * (2 * this->halfDt);
            Int3 localOrigin = truncate((coords - origin) / this->grid->steps + FP3(0.5, 0.5, 0.5));
            Int3 oldLocalOrigin = truncate((oldCoords - origin) / this->grid->steps + FP3(0.5, 0.5, 0.5));
            //std::cout << truncate((coords - origin) / this->grid->steps) << std::endl;
            //std::cout << "oldLocalOrigin: " << oldLocalOrigin << std::endl;
            //std::cout << "localOrigin: " << localOrigin << std::endl;
            FP charge = particle->getCharge() * particle->getWeight() / this->grid->steps.volume() / (2 * this->halfDt);
            FP3 r;
            for(int i = 0; i < 3; i++)
            {
                if(oldLocalOrigin[i] == localOrigin[i]) r[i] = (coords[i] + oldCoords[i]) * (FP)0.5;
                else r[i] = (oldLocalOrigin[i] + localOrigin[i]) * this->grid->steps[i] * (FP)0.5 + origin[i];
            }

            FP3 F[2], W1[2], W2[2];
            W1[0] = ((oldCoords + r) * (FP)0.5 - origin) / this->grid->steps - (FP3)oldLocalOrigin + FP3(0.5, 0.5, 0.5);
            W1[1] = ((coords + r) * (FP)0.5 - origin) / this->grid->steps - (FP3)localOrigin + FP3(0.5, 0.5, 0.5);
            W2[0] = FP3(1.0, 1.0, 1.0) - W1[0];
            W2[1] = FP3(1.0, 1.0, 1.0) - W1[1];
            F[0] = charge * (r - oldCoords);
            F[1] = charge * (coords - r);

            Int3 idx[2];
            idx[0] = oldLocalOrigin;
            idx[1] = localOrigin;
            
            // std::cout << "startGridIdx: " << this->startGridIdx << std::endl;
            // std::cout << "baseGridIdx: " << this->baseGridIdx << std::endl;
            // std::cout << "idx[0]: " << idx[0] << std::endl;
            // std::cout << "idx[1]: " << idx[1] << std::endl;
            for(int i = 0; i < 2; i++)
            {
                this->Jx[idx[i].x][idx[i].y]    [idx[i].z]     += F[i].x * W1[i].y * W1[i].z;
                this->Jx[idx[i].x][idx[i].y - 1][idx[i].z]     += F[i].x * W2[i].y * W1[i].z;
                this->Jx[idx[i].x][idx[i].y]    [idx[i].z - 1] += F[i].x * W1[i].y * W2[i].z;
                this->Jx[idx[i].x][idx[i].y - 1][idx[i].z - 1] += F[i].x * W2[i].y * W2[i].z;

                this->Jy[idx[i].x]    [idx[i].y][idx[i].z]     += F[i].y * W1[i].x * W1[i].z;
                this->Jy[idx[i].x - 1][idx[i].y][idx[i].z]     += F[i].y * W2[i].x * W1[i].z;
                this->Jy[idx[i].x]    [idx[i].y][idx[i].z - 1] += F[i].y * W1[i].x * W2[i].z;
                this->Jy[idx[i].x - 1][idx[i].y][idx[i].z - 1] += F[i].y * W2[i].x * W2[i].z;

                this->Jz[idx[i].x]    [idx[i].y]    [idx[i].z] += F[i].z * W1[i].x * W1[i].y;
                this->Jz[idx[i].x - 1][idx[i].y]    [idx[i].z] += F[i].z * W2[i].x * W1[i].y;
                this->Jz[idx[i].x]    [idx[i].y - 1][idx[i].z] += F[i].z * W1[i].x * W2[i].y;
                this->Jz[idx[i].x - 1][idx[i].y - 1][idx[i].z] += F[i].z * W2[i].x * W2[i].y;
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionZ1 : public CurrentDeposition<TGrid, CurrentDepositionZ1<TGrid>,
        LocalDepositionZ1<TGrid>>
    {
    public:
        CurrentDepositionZ1(double _dt) : CurrentDeposition<TGrid, CurrentDepositionZ1<TGrid>,
            LocalDepositionZ1<TGrid>>(_dt) {}
    };

    template<class TGrid>
    class CurrentDepositionZ2;

    template<class TGrid>
    class LocalDepositionZ2 : public LocalDeposition<TGrid, CurrentDepositionZ2<TGrid>, 5>
    {
    public:
        LocalDepositionZ2(
            const Int3 & blockIdx, CurrentDepositionZ2<TGrid>* currentDeposition, TGrid* _grid)
            : LocalDeposition<TGrid, CurrentDepositionZ2<TGrid>, 5>(blockIdx, currentDeposition, _grid) {}

        template<class T_Particle>
        void depositCurrent(T_Particle* particle)
        {
            FP3 coords = particle->getPosition();
            FP3 minCoords = this->grid->origin + this->startGridIdx * this->grid->steps;
            
            FP3 oldCoords = coords - particle->getVelocity() * (2 * this->halfDt);
            Int3 localOrigin = truncate((coords - minCoords) / this->grid->steps);
            Int3 oldLocalOrigin = truncate((oldCoords - minCoords) / this->grid->steps);
            FP charge = particle->getCharge() * particle->getWeight() / this->grid->steps.volume();

            FP3 r;
            for (int i = 0; i < 3; i++)
            {
                if (oldLocalOrigin[i] == localOrigin[i]) r[i] = (coords[i] + oldCoords[i]) * (FP)0.5;
                else r[i] = std::max(oldLocalOrigin[i], localOrigin[i]) * this->grid->steps[i] + minCoords[i];
            }

            FP3 W[2];
            W[0] = ((oldCoords + r) * (FP)0.5 - minCoords) / this->grid->steps - (FP3)oldLocalOrigin -
                FP3(0.5, 0.5, 0.5);
            W[1] = ((coords + r) * (FP)0.5 - minCoords) / this->grid->steps - (FP3)localOrigin -
                FP3(0.5, 0.5, 0.5);

            FP3 F[2][2];
            FP3 cv = charge * particle->getVelocity() * (FP)0.5;
            F[0][0] = cv * (FP3(0.5, 0.5, 0.5) - W[0]);
            F[1][0] = cv * (FP3(0.5, 0.5, 0.5) - W[1]);
            F[0][1] = cv * (FP3(0.5, 0.5, 0.5) + W[0]);
            F[1][1] = cv * (FP3(0.5, 0.5, 0.5) + W[1]);

            FP3 W1[2], W2[2], W3[2];
            W1[0] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) - W[0]) * (FP3(0.5, 0.5, 0.5) - W[0]);
            W1[1] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) - W[1]) * (FP3(0.5, 0.5, 0.5) - W[1]);

            W2[0] = FP3(0.75, 0.75, 0.75) - W[0] * W[0];
            W2[1] = FP3(0.75, 0.75, 0.75) - W[1] * W[1];

            W3[0] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) + W[0]) * (FP3(0.5, 0.5, 0.5) + W[0]);
            W3[1] = (FP)0.5 * (FP3(0.5, 0.5, 0.5) + W[1]) * (FP3(0.5, 0.5, 0.5) + W[1]);

            Int3 idx[2];
            idx[0] = oldLocalOrigin;
            idx[1] = localOrigin;
            // std::cout << minCoords << std::endl;
            // std::cout << "oldLocalOrigin: " << oldLocalOrigin << std::endl;
            // std::cout << "localOrigin: " << localOrigin << std::endl;
            // std::cout << "startGridIdx: " << this->startGridIdx << std::endl;
            // std::cout << "baseGridIdx: " << this->baseGridIdx << std::endl;
            // std::cout << "idx[0]: " << idx[0] << std::endl;
            // std::cout << "idx[1]: " << idx[1] << std::endl;

            for(int i = 0; i < 2; i++)
            for(int j = 0; j < 2; j++)
            {
                this->Jx[idx[i].x + j][idx[i].y - 1][idx[i].z - 1] += F[i][j].x * W1[i].y * W1[i].z;
                this->Jx[idx[i].x + j][idx[i].y - 1][idx[i].z]     += F[i][j].x * W1[i].y * W2[i].z;
                this->Jx[idx[i].x + j][idx[i].y - 1][idx[i].z + 1] += F[i][j].x * W1[i].y * W3[i].z;
                this->Jx[idx[i].x + j][idx[i].y]    [idx[i].z - 1] += F[i][j].x * W2[i].y * W1[i].z;
                this->Jx[idx[i].x + j][idx[i].y]    [idx[i].z]     += F[i][j].x * W2[i].y * W2[i].z;
                this->Jx[idx[i].x + j][idx[i].y]    [idx[i].z + 1] += F[i][j].x * W2[i].y * W3[i].z;
                this->Jx[idx[i].x + j][idx[i].y + 1][idx[i].z - 1] += F[i][j].x * W3[i].y * W1[i].z;
                this->Jx[idx[i].x + j][idx[i].y + 1][idx[i].z]     += F[i][j].x * W3[i].y * W2[i].z;
                this->Jx[idx[i].x + j][idx[i].y + 1][idx[i].z + 1] += F[i][j].x * W3[i].y * W3[i].z;

                this->Jy[idx[i].x - 1][idx[i].y + j][idx[i].z - 1] += F[i][j].y * W1[i].x * W1[i].z;
                this->Jy[idx[i].x - 1][idx[i].y + j][idx[i].z]     += F[i][j].y * W1[i].x * W2[i].z;
                this->Jy[idx[i].x - 1][idx[i].y + j][idx[i].z + 1] += F[i][j].y * W1[i].x * W3[i].z;
                this->Jy[idx[i].x]    [idx[i].y + j][idx[i].z - 1] += F[i][j].y * W2[i].x * W1[i].z;
                this->Jy[idx[i].x]    [idx[i].y + j][idx[i].z]     += F[i][j].y * W2[i].x * W2[i].z;
                this->Jy[idx[i].x]    [idx[i].y + j][idx[i].z + 1] += F[i][j].y * W2[i].x * W3[i].z;
                this->Jy[idx[i].x + 1][idx[i].y + j][idx[i].z - 1] += F[i][j].y * W3[i].x * W1[i].z;
                this->Jy[idx[i].x + 1][idx[i].y + j][idx[i].z]     += F[i][j].y * W3[i].x * W2[i].z;
                this->Jy[idx[i].x + 1][idx[i].y + j][idx[i].z + 1] += F[i][j].y * W3[i].x * W3[i].z;

                this->Jz[idx[i].x - 1][idx[i].y - 1][idx[i].z + j] += F[i][j].z * W1[i].y * W1[i].x;
                this->Jz[idx[i].x]    [idx[i].y - 1][idx[i].z + j] += F[i][j].z * W1[i].y * W2[i].x;
                this->Jz[idx[i].x + 1][idx[i].y - 1][idx[i].z + j] += F[i][j].z * W1[i].y * W3[i].x;
                this->Jz[idx[i].x - 1][idx[i].y]    [idx[i].z + j] += F[i][j].z * W2[i].y * W1[i].x;
                this->Jz[idx[i].x]    [idx[i].y]    [idx[i].z + j] += F[i][j].z * W2[i].y * W2[i].x;
                this->Jz[idx[i].x + 1][idx[i].y]    [idx[i].z + j] += F[i][j].z * W2[i].y * W3[i].x;
                this->Jz[idx[i].x - 1][idx[i].y + 1][idx[i].z + j] += F[i][j].z * W3[i].y * W1[i].x;
                this->Jz[idx[i].x]    [idx[i].y + 1][idx[i].z + j] += F[i][j].z * W3[i].y * W2[i].x;
                this->Jz[idx[i].x + 1][idx[i].y + 1][idx[i].z + j] += F[i][j].z * W3[i].y * W3[i].x;
            }
        }
    };

    template<class TGrid>
    class CurrentDepositionZ2 : public CurrentDeposition<TGrid, CurrentDepositionZ2<TGrid>,
        LocalDepositionZ2<TGrid>>
    {
    public:
        CurrentDepositionZ2(double _dt) : CurrentDeposition<TGrid, CurrentDepositionZ2<TGrid>,
            LocalDepositionZ2<TGrid>>(_dt) {}
    };

}
