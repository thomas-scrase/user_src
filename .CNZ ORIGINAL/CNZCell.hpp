/* CNZ_HPP
*   implementing single cell model of CNZ
    coupled with Rice et al contraction model

*   Haibo NI qiangzi.ni@gmail.com
*   Date: 21 January 2015 (Wednesday)
    15:04
*/


#ifndef CNZ_HPP
#define CNZ_HPP

#include <iostream>
#include <exception>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "ISAC_Container.hpp"
#include "EnumSimulationCtrl.hpp"  // parameter input enum types;
#include "SingleCellParameter.hpp"  // handling all the parameters regarding different cell types and mutations.
#include "RiceMyofilament.hpp"
#include "ExplicitSolver.hpp"
#include <string>
/*enum AtriaCellType {RA=0, PM, CT, RAA, AS, AVR, BB, LA, LAA, PV, SAN_C, SAN_P};
enum MuationType {WT=0, D322H, E48G, A305T, Y155C, D469E, P488S};
enum AFType {NONE=0, AF1, AF2, AF3, AF4};*/


class CNZ_Cell
{
public:
    CNZ_Cell(AtriaCellType cell_type = RA, MuationType mutation_type = WT, AFType AF_model = NONE, bool usr_IsSAC = true, bool UseICFromFile = false);
    CNZ_Cell(std::vector<double> & vec, AtriaCellType cell_type = RA, MuationType mutation_type = WT, AFType AF_model = NONE, bool usr_IsSAC = true);
    // CNZ_Cell(std::string cell_type = "RA", std::string mutation_type = "WT", AFType AF_model = NONE, bool usr_IsSAC = true, bool UseICFromFile = true);
    ~CNZ_Cell() {};
    double SolveElectricalModel(const double dt);
    void SingleElectricalTimeStep(const double usr_stim, const double usr_dt);
    void RunTest(int Beats = 100, double BCL = 1000.0, double usr_dt = 0.005);
    void InitialiseElectricalStates();
    void InitialiseElectricalStatesFromFile(const char *file_name);
    void SingleTimeStep(const double usr_stim, const double usr_dt);
    void OutputInitialConditions(const char *file_name);
    void InitialiseElectricalStatesFromVector(std::vector<double> vec);

    double state[41];    // states of CNZ model;, updated in SolveElectricalModel();
    ISAC_Container SAC;  // stretch activated channel
    bool UseSAC;
    double Vm;
    double Vmo;
    double temp1, temp2;
    SingleCellPara para;
    double m_t;
    AtriaCellType m_celltype;
    RiceMyofilament force_model;

    double m_IKur    ; // IKur / Cm;
    double m_ICaL    ; // ICaL / Cm;
    double m_INCX    ; // INaCa / Cm;
    double m_Jrel_ss ; // Jrelss;
    double m_ISAC    ; // SAC.ISAC / Cm;

    //access functions via pointer to allow coupling to tissue solver
    double& membrane_potential(){return Vm;}
    double& membrane_potential_previous(){return Vmo;}
    double& sarcomere_length(){return force_model.SL;}
    double& sarcomere_rest(){return force_model.SLrest;}

    void set_membrane_potential(const double new_Vm){Vm = new_Vm; state[0] = new_Vm;}
    void set_membrane_potential_previous(const double new_Vmo){Vmo=new_Vmo;}
    // void set_sarcomere_length(){return force_model.SL;}
    // void set_sarcomere_rest(){return force_model.SLrest;}
};

/*Nxb       = y_myofilament[0];
XBpreR    = y_myofilament[1];
XBpostR   = y_myofilament[2];
SL        = y_myofilament[3];
x_XBpostR = y_myofilament[4];
x_XBpreR  = y_myofilament[5];
TropCaL   = y_myofilament[6];
TropCaH   = y_myofilament[7];
intf      = y_myofilament[8];*/
#endif