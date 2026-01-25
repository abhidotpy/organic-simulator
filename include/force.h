
#pragma once

#include <cmath>
#include "xmath.h"
#include "system.h"

class Force
{
    protected:
    System& system;

    Force(System& sys) : system(sys) {}
    
    public:
    virtual void calculate_forces(int i, int j) = 0;
};

class LennardJones : public Force
{
    private:
    double epsilon;
    double sigma;

    public:
    LennardJones(System& sys, double epsilon, double sigma) : Force(sys), epsilon(epsilon), sigma(sigma) {}

    void calculate_forces(int i, int j);

};

class LennardJones12 : public Force
{
    private:
        double epsilon;
        double sigma;

    public:
    LennardJones12(System& sys, double epsilon, double sigma) : Force(sys), epsilon(epsilon), sigma(sigma) {}
    
    void calculate_forces(int i, int j);
};

class LennardJones_WCA : public Force
{
    private:
        double epsilon;
        double sigma;
    public:
    LennardJones_WCA(System& sys, double epsilon, double sigma) : Force(sys), epsilon(epsilon), sigma(sigma) {}
    
    void calculate_forces(int i, int j);
};

class MorsePotential : public Force
{
    private:
        double D;
        double r0;
        double alpha;

    public:
    MorsePotential(System& sys, double dissociation_energy, double equilibrium_distance, double width)
                     : Force(sys), D(dissociation_energy), r0(equilibrium_distance), alpha(width) {}
    
    void calculate_forces(int i, int j);
};

class GayBernePotential : public Force
{
    private:

};
