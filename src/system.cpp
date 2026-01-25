#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include "xmath.h"
#include "system.h"

void System :: set_num_atoms(int np)
{
    N = np;
    num_atoms_set = true;

    DOF = 3 * N;

    if (N > 0) 
    {
        atom.rx.resize(N); atom.ry.resize(N); atom.rz.resize(N);
        atom.vx.resize(N); atom.vy.resize(N); atom.vz.resize(N);
        atom.fx.resize(N); atom.fy.resize(N); atom.fz.resize(N);
    }
}

void System :: set_density(double dens)
{
    density_set = true;
    density = dens;
}

void System :: add_atom_position(int id, double x, double y, double z)
{
    if (num_atoms_set == false)
    {
        std::cerr << "Error: System has no atoms. Set number of atoms first.\n";
        return;
    }
    else
    {
        try
        {
            atom.rx.at(id) = x;
            atom.ry.at(id) = y;
            atom.rz.at(id) = z;
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "Error: Atom index " << id << " does not exist \
            in system of " << N << " atoms.\n";
        }
    }
}

void System :: add_atom_velocity(int id, double vx, double vy, double vz)
{
    if (num_atoms_set == false)
    {
        std::cerr << "Error: System has no atoms. Set number of atoms first.\n";
        return;
    }
    else
    {
        try
        {
            atom.vx.at(id) = vx;
            atom.vy.at(id) = vy;
            atom.vz.at(id) = vz;
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "Error: Particle index " << id << " does not exist \
            in system of " << N << " atoms.\n";
        }
    }
}

void System :: set_velocity_to_temperature(double temperature, bool remove_com_velocity)
{
    if (num_atoms_set == false)
    {
        std::cout << "Error: System has no atoms. Set number of atoms first.\n";
        return;
    }

    double target_velocity = std::sqrt(temperature);
    double vx_gen, vy_gen, vz_gen;
    double vx_sum, vy_sum, vz_sum;
    vx_sum = 0.0; 
    vy_sum = 0.0; 
    vz_sum = 0.0;

    for (int i = 0; i < N; i++)
    {
        atom.vx.at(i) = target_velocity * std::normal_distribution<double>(0.0, 1.0)(rand_gen);
        atom.vy.at(i) = target_velocity * std::normal_distribution<double>(0.0, 1.0)(rand_gen);
        atom.vz.at(i) = target_velocity * std::normal_distribution<double>(0.0, 1.0)(rand_gen);
    }

    if (remove_com_velocity)
    {
        DOF -= 3;
        for (int i = 0; i < N; i++)
        {
            vx_sum += atom.vx.at(i);
            vy_sum += atom.vy.at(i);
            vz_sum += atom.vz.at(i);
        }

        vx_sum /= double(N);
        vy_sum /= double(N);
        vz_sum /= double(N);

        for (int i = 0; i < N; i++)
        {
            atom.vx.at(i) -= vx_sum; 
            atom.vy.at(i) -= vy_sum; 
            atom.vz.at(i) -= vz_sum; 
        }
    }
}

void System :: set_rotational_dof()
{
    rotations_enabled = true;

    DOF += 3 * N;

    atom. qw.resize(N); atom. qx.resize(N); atom. qy.resize(N); atom. qz.resize(N);
    atom.qdw.resize(N); atom.qdx.resize(N); atom.qdy.resize(N); atom.qdz.resize(N);

    atom.lx.resize(N); atom.ly.resize(N); atom.lz.resize(N);
    atom.wx.resize(N); atom.wy.resize(N); atom.wz.resize(N);
    atom.tx.resize(N); atom.ty.resize(N); atom.tz.resize(N);
}


void System :: initialize_forces()
{
    std::fill(atom.fx.begin(), atom.fx.end(), 0.0);
    std::fill(atom.fy.begin(), atom.fy.end(), 0.0);
    std::fill(atom.fz.begin(), atom.fz.end(), 0.0);

    std::fill(atom.tx.begin(), atom.tx.end(), 0.0);
    std::fill(atom.ty.begin(), atom.ty.end(), 0.0);
    std::fill(atom.tz.begin(), atom.tz.end(), 0.0);

    potential_energy = 0.0;
    kinetic_energy = 0.0;
    temperature = 0.0;
    pressure = 0.0;
    virial = 0.0;
}

void System :: generate_fcc_lattice(int unit_cells)
{
    fcc_lattice_generated = true;
    int num_part = 4 * unit_cells * unit_cells * unit_cells;
    set_num_atoms(num_part);

    if(!box_lengths_set)
    {
        if (!density_set)
        {
            std::cerr << "Density not set. Set density before generating FCC lattice.\n";
            exit(1);
        }
        
        box_lengths_set = true;
        box.volume = double(N) / density;
        box.box_length = pow(box.volume, 1.0 / 3.0);

        box.box_xlen = box.box_length;
        box.box_ylen = box.box_length;
        box.box_zlen = box.box_length;
    }

    double cell = box.box_length / double(unit_cells);
    double half_cell = cell / 2.0;
    double rroot3  = 0.5773503;

    atom.rx.at(0) = 0.0;
    atom.ry.at(0) = 0.0;
    atom.rz.at(0) = 0.0;

    atom.rx.at(1) = half_cell;
    atom.ry.at(1) = half_cell;
    atom.rz.at(1) = 0.0;

    atom.rx.at(2) = 0.0;
    atom.ry.at(2) = half_cell;
    atom.rz.at(2) = half_cell;

    atom.rx.at(3) = half_cell;
    atom.ry.at(3) = 0.0;
    atom.rz.at(3) = half_cell;

    int M = 0;
    for (int i = 0; i < unit_cells; i++)
    {
        for (int j = 0; j < unit_cells; j++)
        {
            for (int k = 0; k < unit_cells; k++)
            {
                for (int iref = 0; iref < 4; iref++)
                {
                    atom.rx.at(iref + M) = atom.rx.at(iref) + cell * i;
                    atom.ry.at(iref + M) = atom.ry.at(iref) + cell * j;
                    atom.rz.at(iref + M) = atom.rz.at(iref) + cell * k;
                }
                M += 4;
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        atom.rx.at(i) = atom.rx.at(i) - box.box_length / 2.0;
        atom.ry.at(i) = atom.ry.at(i) - box.box_length / 2.0;
        atom.rz.at(i) = atom.rz.at(i) - box.box_length / 2.0;
    }
}

void System :: initialize()
{
    std::cout << "This is SUCROSE v1.0\n";

    if (num_atoms_set)
    {
        std::cout << "Number of atoms in system set to " << N << "\n";
    }
    else
    {
        std::cerr << "No atoms in system.\n";
        exit(1);
    }

    if (density_set)
    {
        std::cout << "Density of system set to " << density << " atoms per cubic unit cell\n";
    }
    else
    {
        std::cerr << "Density not set for system. Cannot determine simualtion box size.\n";
        exit(1);
    }

    if (box_lengths_set)
    {
        std::cout << "Box length of system set to " << box.box_length << "\n";
    }
    else
    {
        box_lengths_set = true;
        box.volume = double(N) / density;
        box.box_length = pow(box.volume, 1.0 / 3.0);

        box.box_xlen = box.box_length;
        box.box_ylen = box.box_length;
        box.box_zlen = box.box_length;
        std::cout << "Box length of system set to " << box.box_length << "\n";
    }

    if (box.is_periodic)
        std::cout << "Periodic boundary conditions enabled for current system\n";
    else
        std::cout << "Periodic boundary conditions disabled for current system\n";

    if (rotations_enabled)
        std::cout << "Rotational degrees of freedom enabled for current system\n";
    else
        std::cout << "Rotational degrees of freedom disabled for current system\n";

    if (default_mass)
            atom.mass = 1.0;

        if (default_charge)
            atom.charge = 0.0;

    if (fcc_lattice_generated)
        std::cout << "FCC lattice generated with " << N << " atoms\n";
    
    initialize_forces();

    std::cout << "System initialization complete...\n\n";
}


void System :: calculate_state_variables()
{
    kinetic_energy = 0.0;
    temperature = 0.0;
    pressure = 0.0;

    for (int i = 0; i < N; i++)
    {
        kinetic_energy += 0.5 * atom.mass * (atom.vx[i] * atom.vx[i] + atom.vy[i] * atom.vy[i] + atom.vz[i] * atom.vz[i]);
    }

    total_energy = potential_energy + kinetic_energy;

    temperature = 2.0 * kinetic_energy / double(DOF);
    
    if (box.volume > 0)
    {
        pressure = (density * temperature) + (virial / (3.0 * box.volume));
    }

    potential_energy /= double(N);
    kinetic_energy /= double(N);
    total_energy /= double(N);
}
