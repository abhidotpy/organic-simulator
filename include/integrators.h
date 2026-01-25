
#pragma once

#include <cmath>
#include "xmath.h"
#include "system.h"


class Integrator
{
    protected:
    System& system;
    double dt;

    public:
    Integrator(System& sys, double dt) : system(sys), dt(dt) {}

    ~Integrator() {}

    double get_timestep() { return dt; }

    void set_timestep(double dt) { this->dt = dt; }

    virtual void initial_step() = 0;
    virtual void final_step() = 0;
};

class VerletIntegrator : public Integrator
{
    public:
    VerletIntegrator(System& sys, double dt) : Integrator(sys, dt) {}
    
    void initial_step();
    void final_step();

};
