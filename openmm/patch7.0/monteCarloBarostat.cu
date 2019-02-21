/**
 * Scale the particle positions with each axis independent
 */

extern "C" __global__ void scalePositions(float scaleX, float scaleY, float scaleZ, int numMolecules, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, real4* __restrict__ posq,
        const int* __restrict__ moleculeAtoms, const int* __restrict__ moleculeStartIndex) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numMolecules; index += blockDim.x*gridDim.x) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        // Find the center of each molecule.

        real3 center = make_real3(0, 0, 0);
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            center.x += pos.x;
            center.y += pos.y;
            center.z += pos.z;
        }
        real invNumAtoms = RECIP(numAtoms);
        center.x *= invNumAtoms;
        center.y *= invNumAtoms;
        center.z *= invNumAtoms;

//VO 2018 changes below
        real3 oldCenter = center;
        // Move it into the first periodic box
        APPLY_PERIODIC_TO_POS(center);
        // Scale molecule center
        center.x*=scaleX;
        center.y*=scaleY;
        center.z*=scaleZ;
        // Restore to original position
        RESTORE_PERIODIC(center,oldCenter);
        real3 delta = make_real3(center.x-oldCenter.x, center.y-oldCenter.y, center.z-oldCenter.z);
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}
