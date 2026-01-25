#include <cmath>
#include "xmath.h"
#include "system.h"
#include "force.h"

void LennardJones::calculate_forces(int i, int j)
{
    Vector3 rij = Vector3(0.0, 0.0, 0.0);
    Vector3 fij = Vector3(0.0, 0.0, 0.0);
    
    rij.x = system.atom.rx[i] - system.atom.rx[j];
    rij.y = system.atom.ry[i] - system.atom.ry[j];
    rij.z = system.atom.rz[i] - system.atom.rz[j];

    if (system.box.is_periodic)
    {
        rij.x -= system.box.box_length * rint(rij.x / system.box.box_length);
        rij.y -= system.box.box_length * rint(rij.y / system.box.box_length);
        rij.z -= system.box.box_length * rint(rij.z / system.box.box_length);
    }
    
    double rij_sq = rij.mag2();
    double sr2_lj = pow(sigma, 2) / rij_sq;
    system.potential_energy += 4 * epsilon * ( pow(sr2_lj, 6) - pow(sr2_lj, 3) );

    double vir = 24.0 * epsilon * ( 2.0 * pow(sr2_lj, 6) - pow(sr2_lj, 3) );
    fij = vir * rij / rij_sq;
    system.virial += vir;
    
    system.atom.fx[i] += fij.x;
    system.atom.fy[i] += fij.y;
    system.atom.fz[i] += fij.z;

    system.atom.fx[j] -= fij.x;
    system.atom.fy[j] -= fij.y;
    system.atom.fz[j] -= fij.z;
}

void LennardJones12::calculate_forces(int i, int j)
{
    Vector3 rij = Vector3(0.0, 0.0, 0.0);
    Vector3 fij = Vector3(0.0, 0.0, 0.0);
    
    rij.x = system.atom.rx[i] - system.atom.rx[j];
    rij.y = system.atom.ry[i] - system.atom.ry[j];
    rij.z = system.atom.rz[i] - system.atom.rz[j];

    if (system.box.is_periodic)
    {
        rij.x -= system.box.box_length * rint(rij.x / system.box.box_length);
        rij.y -= system.box.box_length * rint(rij.y / system.box.box_length);
        rij.z -= system.box.box_length * rint(rij.z / system.box.box_length);
    }
    
    double rij_sq = rij.mag2();
    double sr2_lj = pow(sigma, 2) / rij_sq;
    system.potential_energy += epsilon * ( pow(sr2_lj, 6) - 2.0 * pow(sr2_lj, 3) );

    system.virial += 12.0 * epsilon * ( pow(sr2_lj, 6) - pow(sr2_lj, 3) );
    fij = system.virial * rij / rij_sq;
    
    system.atom.fx[i] += fij.x;
    system.atom.fy[i] += fij.y;
    system.atom.fz[i] += fij.z;

    system.atom.fx[j] -= fij.x;
    system.atom.fy[j] -= fij.y;
    system.atom.fz[j] -= fij.z;
}

void LennardJones_WCA::calculate_forces(int i, int j)
{
    Vector3 rij = Vector3(0.0, 0.0, 0.0);
    Vector3 fij = Vector3(0.0, 0.0, 0.0);
    
    rij.x = system.atom.rx[i] - system.atom.rx[j];
    rij.y = system.atom.ry[i] - system.atom.ry[j];
    rij.z = system.atom.rz[i] - system.atom.rz[j];

    if (system.box.is_periodic)
    {
        rij.x -= system.box.box_length * rint(rij.x / system.box.box_length);
        rij.y -= system.box.box_length * rint(rij.y / system.box.box_length);
        rij.z -= system.box.box_length * rint(rij.z / system.box.box_length);
    }
    
    double rij_sq = rij.mag2();
    double sigma_sq = pow(sigma, 2);
    double sr2_lj = sigma_sq / rij_sq;

    if ( rij_sq <= sigma_sq )
    {
        system.potential_energy += epsilon * ( pow(sr2_lj, 6) - 2.0 * pow(sr2_lj, 3) ) + 1.0;
        system.virial += 12.0 * epsilon * ( pow(sr2_lj, 6) - pow(sr2_lj, 3) );
    }
    else
    {
        system.potential_energy += 0.0;
        system.virial += 0.0;
    }
    
    fij = system.virial * rij / rij_sq;
    
    system.atom.fx[i] += fij.x;
    system.atom.fy[i] += fij.y;
    system.atom.fz[i] += fij.z;

    system.atom.fx[j] -= fij.x;
    system.atom.fy[j] -= fij.y;
    system.atom.fz[j] -= fij.z;
}

void MorsePotential::calculate_forces(int i, int j)
{
        Vector3 rij = Vector3(0.0, 0.0, 0.0);
        Vector3 fij = Vector3(0.0, 0.0, 0.0);
        
        rij.x = system.atom.rx[i] - system.atom.rx[j];
        rij.y = system.atom.ry[i] - system.atom.ry[j];
        rij.z = system.atom.rz[i] - system.atom.rz[j];

        if (system.box.is_periodic)
        {
            rij.x -= system.box.box_length * rint(rij.x / system.box.box_length);
            rij.y -= system.box.box_length * rint(rij.y / system.box.box_length);
            rij.z -= system.box.box_length * rint(rij.z / system.box.box_length);
        }
        
        double rp = rij.mag();

        system.potential_energy += D * pow((1.0 - exp(-alpha * (rp - r0))), 2.0) - D;
        fij = -2.0 * alpha * D * exp(-alpha * (rp - r0)) * (1.0 - exp(-alpha * (rp - r0))) * rij / rp;
        
        system.atom.fx[i] += fij.x;
        system.atom.fy[i] += fij.y;
        system.atom.fz[i] += fij.z;

        system.atom.fx[j] -= fij.x;
        system.atom.fy[j] -= fij.y;
        system.atom.fz[j] -= fij.z;
    }
