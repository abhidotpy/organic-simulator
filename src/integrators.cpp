
#include <cmath>
#include "xmath.h"
#include "system.h"
#include "integrators.h"

void VerletIntegrator :: initial_step()
{
    int np = system.N;
    for(int i = 0; i < np; i++)
    {
        system.atom.vx[i] += system.atom.fx[i] * dt / 2.0;
        system.atom.vy[i] += system.atom.fy[i] * dt / 2.0;
        system.atom.vz[i] += system.atom.fz[i] * dt / 2.0;
        
        system.atom.rx[i] += ( system.atom.vx[i] * dt );
        system.atom.ry[i] += ( system.atom.vy[i] * dt );
        system.atom.rz[i] += ( system.atom.vz[i] * dt );

        if (system.box.is_periodic)
        {
            system.atom.rx[i] -= std::rint(system.atom.rx[i] / system.box.box_length) * system.box.box_length;
            system.atom.ry[i] -= std::rint(system.atom.ry[i] / system.box.box_length) * system.box.box_length;
            system.atom.rz[i] -= std::rint(system.atom.rz[i] / system.box.box_length) * system.box.box_length;
        }
    }
}

void VerletIntegrator :: final_step()
{
    int np = system.N;
    for(int i = 0; i < np; i++)
    {
        system.atom.vx[i] += system.atom.fx[i] * dt / 2.0;
        system.atom.vy[i] += system.atom.fy[i] * dt / 2.0;
        system.atom.vz[i] += system.atom.fz[i] * dt / 2.0;
    }
}
