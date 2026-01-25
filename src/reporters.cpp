#include <iostream>
#include <fstream>
#include <iomanip>
#include "system.h"
#include "reporters.h"

StateDataReporter::StateDataReporter(System& sys, const std::string& filename) : Reporter(sys, filename)
{
    file << "Step\t"
         << "Potential_energy\t"
         << "Kinetic_energy\t\t"
         << "Total_energy\t\t"
         << "Temperature\t\t\t"
         << "Pressure\t\t\t"
         << "Density\n";
}

void StateDataReporter::report(int step)
{
    file << step << "\t\t" << std::setprecision(10) << std::fixed
         << system.potential_energy << "\t\t"
         << system.kinetic_energy << "\t\t"
         << system.total_energy << "\t\t"
         << system.temperature << "\t\t"
         << system.pressure << "\t\t"
         << system.density << "\n";
}

TrajectoryReporter::TrajectoryReporter(System& sys, const std::string& filename) : Reporter(sys, filename) {}

void TrajectoryReporter :: report(int step)
{
    int np = system.N;
    file << "ITEM: TIMESTEP\n"
            << step << "\n"
            << "ITEM: NUMBER OF ATOMS\n" 
            << system.N 
            << "\nITEM: BOX BOUNDS pp pp pp\n"
            << -system.box.box_length / 2 << "\t" << system.box.box_length / 2 << "\n"
            << -system.box.box_length / 2 << "\t" << system.box.box_length / 2 << "\n"
            << -system.box.box_length / 2 << "\t" << system.box.box_length / 2 << "\n"
            << "ITEM: ATOMS id radius x y z vx vy vz fx fy fz\n";

    for (int i = 0; i < system.N; i++)
    {
        file << i + 1 << "\t" << 0.5 << "\t"
                << system.atom.rx[i] << "\t" << system.atom.ry[i] << "\t" << system.atom.rz[i] << "\t" 
                << system.atom.vx[i] << "\t" << system.atom.vy[i] << "\t" << system.atom.vz[i] << "\t"
                << system.atom.fx[i] << "\t" << system.atom.fy[i] << "\t" << system.atom.fz[i] << "\n";
    }
}

ConfigurationReporter::ConfigurationReporter(System& sys, const std::string& filename) : Reporter(sys, filename) {}

void ConfigurationReporter::report(int step)
{
    int np = system.N;
        for (int i = 0; i < np; i++)
        {
            file << std::setprecision(10) << std::fixed
                 << system.atom.rx[i] << "\t" << system.atom.ry[i] << "\t" << system.atom.rz[i] << "\t" 
                 << system.atom.vx[i] << "\t" << system.atom.vy[i] << "\t" << system.atom.vz[i] << "\t"
                 << system.atom.fx[i] << "\t" << system.atom.fy[i] << "\t" << system.atom.fz[i] << "\n";
        }
        file << "\n";
}
