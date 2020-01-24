// Header file for the ICs and RHS functions for the Colman 2013 model

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>

class FileSys {
public:
    std::string Pacemap;
    std::string FB_map;
    std::string Fcell_map;
    std::string Geometry_file;
    std::string Fibre_theta;
    std::string Fibre_phi;
    std::string Stim_file;
    std::string SAN_type;
    std::string Stim_time_file;
    std::string Stim_amp_file;
    std::string Apicobasal_file;
    std::string RVIndex_file;
    FileSys();
    ~FileSys();
};

class Simulation_Config: public FileSys {
public:
    int BCL;
    double Total_time;
    float dt;
    float ISO;
    float Ach;
    float Diff_Scale;
    double S2;
    int S1_number;

    int Model_type;
    int model_out;
    int IKur_type;
    int region;
    int region_3D;
    int AF_model;
    int mutation;
    int FB_type;
    int FB_number;
    double Ggap;
    int tau_type;

    std::string mutation_char;
    std::string region_char;
    std::string tau_type_char;
    std::string ICs;
    std::string Stim_type;
    std::string Model_type_char;

    
    Simulation_Config();
    Simulation_Config(int argc, char *argv[]);  // another constructor
    ~Simulation_Config();
    void Initilise();
    void Config_handling(int argc, char *argv[]);
    void Report_Config();
    void Report_All();
};




#endif
