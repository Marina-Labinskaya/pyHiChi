#pragma once

#include "FieldSolver.h"
#include "Grid.h"

#include <iostream>

namespace pfc
{
    template<GridTypes gridTypes>
    class FieldSolver;
    template<GridTypes gridTypes>
    class RealFieldSolver;

    template <GridTypes gridTypes>
    class CurrentBoundaryConditions
    {
    public:
        CurrentBoundaryConditions(RealFieldSolver<gridTypes>* fieldSolver = 0);
        virtual void updateCurrentBoundaries();
    protected:
        RealFieldSolver<gridTypes>* fieldSolver;
        // major index is index of edge, minor index is index of component
        std::function<FP(FP, FP, FP, FP)> jLeft[3][3];
        std::function<FP(FP, FP, FP, FP)> jRight[3][3];
        FP3 leftCoeff;
        FP3 rightCoeff;
    };

    template<GridTypes gridTypes>
    CurrentBoundaryConditions<gridTypes>::CurrentBoundaryConditions(RealFieldSolver<gridTypes>* _fieldSolver)
    {
        fieldSolver = _fieldSolver;
        leftCoeff = FP3(0, 0, 0);
        rightCoeff = FP3(0, 0, 0);
        for (int d = 0; d < 3; ++d)
        {
            //in seq if fieldGeneration in area
            leftCoeff[d] = 1;
            rightCoeff[d] = 1;
        }
    }

    template<GridTypes gridTypes>
    void CurrentBoundaryConditions<gridTypes>::updateCurrentBoundaries()
    {
        Grid<FP, gridTypes> * grid = fieldSolver->grid;
        const FP time = fieldSolver->globalTime + fieldSolver->timeShiftB;
        const FP cdt = constants::c * fieldSolver->dt;
        const FP3 norm_coeffs = FP3(cdt, cdt, cdt) / grid->steps;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;
            int begin1 = fieldSolver->internalEAreaBegin[dim1];
            int begin2 = fieldSolver->internalEAreaBegin[dim2];
            int end1 = fieldSolver->internalEAreaEnd[dim1];
            int end2 = fieldSolver->internalEAreaEnd[dim2];
            //OMP_FOR_COLLAPSE()
            for (int j = begin1; j < end1; j++)
                for (int k = begin2; k < end2; k++)
                {
                    // Adjust indexes for symmetry of generation coordinates, use B indexes for consistency
                    /// Видимо, не нужно добавлять здесь 1, а ниже вычитать 1 вместо 2
                    Int3 index;
                    index[dim0] = fieldSolver->internalBAreaBegin[dim0] + 1;
                    index[dim1] = j;
                    index[dim2] = k;
                    FP3 jxCoords, jyCoords, jzCoords;
                    Int3 indexes[3] = { index, index, index };
                    indexes[dim0][dim0]++;
                    jxCoords = grid->JxPosition(indexes[0].x, indexes[0].y,
                        indexes[0].z);
                    jyCoords = grid->JyPosition(indexes[1].x, indexes[1].y,
                        indexes[1].z);
                    jzCoords = grid->JzPosition(indexes[2].x, indexes[2].y,
                        indexes[2].z);
                    FP coeff = leftCoeff[dim0] * norm_coeffs[dim0];
                    grid->Jx(indexes[0]) += coeff * jLeft[dim0][0](
                        jxCoords.x, jxCoords.y, jxCoords.z, time);
                    grid->Jy(indexes[1]) += coeff * jLeft[dim0][1](
                        jyCoords.x, jyCoords.y, jyCoords.z, time);
                    grid->Jz(indexes[2]) += coeff * jLeft[dim0][2](
                        jzCoords.x, jzCoords.y, jzCoords.z, time);

                    index[dim0] = fieldSolver->internalBAreaEnd[dim0] - 2;
                    jxCoords += grid->JxPosition(index.x, index.y, index.z);
                    jyCoords += grid->JyPosition(index.x, index.y, index.z);
                    jzCoords += grid->JzPosition(index.x, index.y, index.z);
                    coeff = rightCoeff[dim0] * norm_coeffs[dim0];
                    grid->Jx(index) += coeff * jRight[dim0][0](
                        jxCoords.x, jxCoords.y, jxCoords.z, time);
                    grid->Jy(index) += coeff * jRight[dim0][1](
                        jyCoords.x, jyCoords.y, jyCoords.z, time);
                    grid->Jz(index) += coeff * jRight[dim0][2](
                        jzCoords.x, jzCoords.y, jzCoords.z, time);
                }
        }
    }

    template<GridTypes gridTypes>
    class PeriodicalCurrentBoundaryConditions : public CurrentBoundaryConditions<gridTypes>
    {
    public:
        PeriodicalCurrentBoundaryConditions(RealFieldSolver<gridTypes>* fieldSolver = 0) :
            CurrentBoundaryConditions<gridTypes>(fieldSolver) {}
        virtual void updateCurrentBoundaries();
    };

    template<GridTypes gridTypes>
    void PeriodicalCurrentBoundaryConditions<gridTypes>::updateCurrentBoundaries()
    {
        Grid<FP, gridTypes>* grid = this->fieldSolver->grid;
        for (int dim0 = 0; dim0 < grid->dimensionality; dim0++)
        {
            int dim1 = (dim0 + 1) % 3;
            int dim2 = (dim0 + 2) % 3;
            int begin1 = 0;
            int begin2 = 0;
            int end1 = grid->numCells[dim1];
            int end2 = grid->numCells[dim2];

            //OMP_FOR_COLLAPSE()
            for (int j = begin1; j < end1; j++)
                for (int k = begin2; k < end2; k++)
                {
                    // Adjust indexes for symmetry of generation coordinates
                    Int3 indexL, indexR;
                    indexL[dim0] = 0;
                    indexL[dim1] = j;
                    indexL[dim2] = k;
                    indexR[dim0] = indexL[dim0] + grid->numInternalCells[dim0];
                    indexR[dim1] = j;
                    indexR[dim2] = k;

                    const int numExchangeCells = grid->getNumExternalLeftCells()[dim0];

                    for (int i = 0; i < numExchangeCells; i++) {
                        indexL[dim0] += i;
                        indexR[dim0] += i;

                        grid->Jx(indexR) += grid->Jx(indexL);
                        grid->Jy(indexR) += grid->Jy(indexL);
                        grid->Jz(indexR) += grid->Jz(indexL);

                        grid->Jx(indexL) = grid->Jx(indexR);
                        grid->Jy(indexL) = grid->Jy(indexR);
                        grid->Jz(indexL) = grid->Jz(indexR);
                    }

                    indexL[dim0] = grid->getNumExternalLeftCells()[dim0];
                    indexR[dim0] = indexL[dim0] + grid->numInternalCells[dim0];

                    for (int i = 0; i < numExchangeCells; i++) {
                        indexL[dim0] += i;
                        indexR[dim0] += i;

                        grid->Jx(indexL) += grid->Jx(indexR);
                        grid->Jy(indexL) += grid->Jy(indexR);
                        grid->Jz(indexL) += grid->Jz(indexR);

                        grid->Jx(indexR) = grid->Jx(indexL);
                        grid->Jy(indexR) = grid->Jy(indexL);
                        grid->Jz(indexR) = grid->Jz(indexL);
                    }
                }
        }
    }
}