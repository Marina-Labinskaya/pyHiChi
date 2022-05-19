#include "TestingUtility.h"

#include "CurrentDeposition.h"
#include "Constants.h"
#include "FDTD.h"
#include "FieldValue.h"
#include "Particle.h"
//#include "Pusher.h"
#include <random>

using namespace pfc;

TEST(TCurrentDeposition, CurrentDepositionForOneParticleIsRight) {
    FP stepx = 0.1, stepy = 0.2, stepz = 0.1;

    Int3 numInternalCells = Int3(10, 7, 12);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepx, stepy, stepz);
    Int3 globalGridDims = numInternalCells;

    Int3 idxX, idxY, idxZ;
    FP3 internalCoordsX, internalCoordsY, internalCoordsZ;
    Real charge = constants::electronCharge;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);

    ScalarField<FP> expectJx(grid.numCells), expectJy(grid.numCells), expectJz(grid.numCells);

    expectJx.zeroize();
    expectJy.zeroize();
    expectJz.zeroize();

    FP x = (FP)0.3 * grid.steps.x;
    FP y = (FP)0.1 * grid.steps.y;
    FP z = (FP)0.4 * grid.steps.z;
    Particle3d::PositionType position(x, y, z);
    Particle3d::MomentumType momentum(FP3(15.0, 7.0, 10.0));
    Particle3d particle(position, momentum);

    FP3 JBefore = (particle.getVelocity() * particle.getCharge()) / (steps.x * steps.y * steps.z);

    FirstOrderCurrentDeposition<YeeGridType> currentdeposition;
    currentdeposition(&grid, particle);

    grid.getIndexJxCoords(position, idxX, internalCoordsX);
    currentdeposition.CurrentDensityAfterDeposition(expectJx, idxX, internalCoordsX, JBefore.x);

    grid.getIndexJyCoords(position, idxY, internalCoordsY);
    currentdeposition.CurrentDensityAfterDeposition(expectJy, idxY, internalCoordsY, JBefore.y);

    grid.getIndexJzCoords(position, idxZ, internalCoordsZ);
    currentdeposition.CurrentDensityAfterDeposition(expectJz, idxZ, internalCoordsZ, JBefore.z);

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                EXPECT_NEAR (grid.Jx(i, j, k), expectJx(i, j, k), 0.00001);
                EXPECT_NEAR (grid.Jy(i, j, k), expectJy(i, j, k), 0.00001);
                EXPECT_NEAR (grid.Jz(i, j, k), expectJz(i, j, k), 0.00001);
            }
}

TEST(TCurrentDeposition, CurrentDepositionForManyParticlesIsRight) {
    FP stepX = 0.1, stepY = 0.2, stepZ = 0.1;
    FP minpX = 1.0, minpY = 1.0, minpZ = 1.0;
    FP maxpX = 10.0, maxpY = 10.0, maxpZ = 10.0;
    FP minX = 0.0, minY = 0.0, minZ = 0.0;
    FP maxX, maxY, maxZ;

    Int3 numInternalCells = Int3(10, 7, 12);
    FP3 minCoords = FP3(0.0, 0.0, 0.0);
    FP3 steps = FP3(stepX, stepY, stepZ);
    Int3 globalGridDims = numInternalCells;

    Int3 idxX, idxY, idxZ;
    FP3 internalCoordsX, internalCoordsY, internalCoordsZ;

    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
    maxX = minX + grid.numInternalCells.x * stepX;
    maxY = minY + grid.numInternalCells.y * stepY;
    maxZ = minZ + grid.numInternalCells.z * stepZ;

    ScalarField<FP> expectJx(grid.numCells), expectJy(grid.numCells), expectJz(grid.numCells);
    expectJx.zeroize();
    expectJy.zeroize();
    expectJz.zeroize();

    ParticleArray3d particleArray;
    int numParticles = 10;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disX(minX, maxX);
    std::uniform_real_distribution<> disY(minY, maxY);
    std::uniform_real_distribution<> disZ(minZ, maxZ);

    std::uniform_real_distribution<> dispX(minpX, maxpX);
    std::uniform_real_distribution<> dispY(minpY, maxpY);
    std::uniform_real_distribution<> dispZ(minpZ, maxpZ);

    for (int i = 0; i < numParticles; i++)
    {
        Particle3d::PositionType position(disX(gen), disY(gen), disZ(gen));
        Particle3d::MomentumType momentum(dispX(gen), dispY(gen), dispZ(gen));
        Particle3d particle(position, momentum);
        particleArray.pushBack(particle);
    }

    FirstOrderCurrentDeposition<YeeGridType> currentdeposition;
    currentdeposition(&grid, &particleArray);

    for (int i = 0; i < numParticles; i++)
    {
        FP3 JBefore = (particleArray[i].getVelocity() * particleArray[i].getCharge()) / (steps.x * steps.y * steps.z);

        grid.getIndexJxCoords(particleArray[i].getPosition(), idxX, internalCoordsX);
        currentdeposition.CurrentDensityAfterDeposition(expectJx, idxX, internalCoordsX, JBefore.x);

        grid.getIndexJyCoords(particleArray[i].getPosition(), idxY, internalCoordsY);
        currentdeposition.CurrentDensityAfterDeposition(expectJy, idxY, internalCoordsY, JBefore.y);

        grid.getIndexJzCoords(particleArray[i].getPosition(), idxZ, internalCoordsZ);
        currentdeposition.CurrentDensityAfterDeposition(expectJz, idxZ, internalCoordsZ, JBefore.z);
    }

    for (int i = 0; i < grid.numCells.x; ++i)
        for (int j = 0; j < grid.numCells.y; ++j)
            for (int k = 0; k < grid.numCells.z; ++k) {
                EXPECT_NEAR(grid.Jx(i, j, k), expectJx(i, j, k), 0.00001);
                EXPECT_NEAR(grid.Jy(i, j, k), expectJy(i, j, k), 0.00001);
                EXPECT_NEAR(grid.Jz(i, j, k), expectJz(i, j, k), 0.00001);
            }
}

//TEST(TCurrentDeposition, particleCanGoThroughACycle) {
//    // create grid
//    FP stepX = 0.01, stepY = 0.02, stepZ = 0.01;
//    Int3 numInternalCells = Int3(10, 20, 10);
//    FP3 minCoords = FP3(0.0, 0.0, 0.0);
//    FP3 steps = FP3(stepX, stepY, stepZ);
//    Int3 globalGridDims = Int3(0, 0, 0);
//    YeeGrid grid(numInternalCells, minCoords, steps, globalGridDims);
//
//    // field solver
//    const double dt = 0.001;
//    FDTD fdtd(&grid, dt);
//    for (int i = 0; i < 100; ++i) {
//        fdtd.updateFields();
//    }
//
//    // interpolate
//    grid.setInterpolationType(Interpolation_CIC);
//    //grid.getFields();
//    // here interpolation foe fields in Ex, ...
//
//    // pusher
//    ParticleArray3d particleArray;
//    FP Ex = 0.0, Ey = 0.0, Ez = 0.0;
//    FP Bx = 1.0, By = 1.0, Bz = 1.0;
//    int numParticles = 10;
//
//    FP minpX = 1.0, minpY = 1.0, minpZ = 1.0;
//    FP maxpX = 10.0, maxpY = 10.0, maxpZ = 10.0;
//    FP minX = 0.0, minY = 0.0, minZ = 0.0;
//    FP maxX, maxY, maxZ;
//
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> disX(minX, maxX);
//    std::uniform_real_distribution<> disY(minY, maxY);
//    std::uniform_real_distribution<> disZ(minZ, maxZ);
//
//    std::uniform_real_distribution<> dispX(minpX, maxpX);
//    std::uniform_real_distribution<> dispY(minpY, maxpY);
//    std::uniform_real_distribution<> dispZ(minpZ, maxpZ);
//
//    for (int i = 0; i < numParticles; i++)
//    {
//        Particle3d::PositionType position(disX(gen), disY(gen), disZ(gen));
//        Particle3d::MomentumType momentum(dispX(gen), dispY(gen), dispZ(gen));
//        Particle3d particle(position, momentum);
//        particleArray.pushBack(particle);
//    }
//
//    /*BorisPusher scalarPusher;
//    FP timeStep = 0.01;
//    auto field = ValueField(Ex, Ey, Ez, Bx, By, Bz);
//    scalarPusher(&particleArray, field, timeStep);*/
//
//    // current deposition
//    FirstOrderCurrentDeposition<YeeGridType> currentdeposition;
//    currentdeposition(&grid, &particleArray);
//}
