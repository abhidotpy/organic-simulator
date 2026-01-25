#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include "system.h"

class Reporter
{
    protected:
    System& system;
    std::ofstream file;

    public:
    Reporter(System& sys, const std::string& filename) : system(sys), file(filename) {}

    virtual ~Reporter()
    {
        if (file.is_open())
        {
            file.close();
        }
    }

    virtual void report(int step) = 0;
};

class StateDataReporter : public Reporter
{
    public:
    StateDataReporter(System& sys, const std::string& filename);
    ~StateDataReporter()
    {
        if (file.is_open())
        {
            file.close();
        }
    }
    void report(int step);
};

class TrajectoryReporter : public Reporter
{
    public:
    TrajectoryReporter(System& sys, const std::string& filename);
    ~TrajectoryReporter()
    {
        if (file.is_open())
        {
            file.close();
        }
    }
    void report(int step);
};

class ConfigurationReporter : public Reporter
{
    public:
    ConfigurationReporter(System& sys, const std::string& filename);
    ~ConfigurationReporter()
    {
        if (file.is_open())
        {
            file.close();
        }
    }
    void report(int step);
};
