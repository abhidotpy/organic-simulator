#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include "xmath.h"

struct Atom
{
    double mass;
    double charge;

    std::array <double, 3> radius;
    std::array <double, 3> shape; 

    std::vector <double> rx, ry, rz;
    std::vector <double> ex, ey, ez;
    std::vector <double> vx, vy, vz;
    std::vector <double> fx, fy, fz;

    std::vector <double> qw, qx, qy, qz;
    std::vector <double> qdw, qdx, qdy, qdz;
    std::vector <double> lx, ly, lz;
    std::vector <double> wx, wy, wz;
    std::vector <double> tx, ty, tz;
};

struct Box
{
    double box_length;
    double volume;
    double box_xlen, box_ylen, box_zlen;
    bool is_periodic;
};

class System
{
    public:
    int N;
    int DOF;
    double potential_energy;
    double kinetic_energy;
    double total_energy;
    double temperature;
    double pressure;
    double density;
    double virial;
    
    bool default_mass = true, default_charge = true;
    bool rotations_enabled = false;

    bool box_lengths_set = false; 
    bool num_atoms_set = false;
    bool density_set = false;

    bool fcc_lattice_generated = false;

    Atom atom;
    Box box;
    std::mt19937 rand_gen;

    System() :  default_mass(true), default_charge(true), 
                box_lengths_set(false), num_atoms_set(false), density_set(false),
                rotations_enabled(false), fcc_lattice_generated(false), 
                N(0), DOF(0), potential_energy(0.0), kinetic_energy(0.0), 
                total_energy(0.0), temperature(0.0), pressure(0.0), density(0.0), virial(0.0)
    {
        std::random_device rd;
        rand_gen = std::mt19937(rd());
        atom.mass = 1.0;
        atom.charge = 0.0;

        box.box_length = 0.0;
        box.volume = 0.0;
        box.box_xlen = 0.0;
        box.box_ylen = 0.0;
        box.box_zlen = 0.0;
        box.is_periodic = false;
    }

    ~System() {}

    void set_system_mass(double mass)
    {
        default_mass = false;
        atom.mass = mass;
    }

    void set_system_charge(double charge)
    {   
        default_charge = false;
        atom.charge = charge;
    }

    void set_periodic_boundaries()
    {
        box.is_periodic = true;
    }

    void set_num_atoms(int np);

    void set_density(double dens);

    void add_atom_position(int id, double x, double y, double z);

    void add_atom_velocity(int id, double vx, double vy, double vz);

    void set_velocity_to_temperature(double temperature, bool remove_com_velocity=true);
    
    void set_rotational_dof();

    void initialize_forces();

    void generate_fcc_lattice(int unit_cells);

    void initialize();

    void calculate_state_variables();
};
