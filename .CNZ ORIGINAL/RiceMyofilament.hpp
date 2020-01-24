/*
* A class implement Rice Myofilament model
* Rice et al. 2008

* Haibo Ni
* qiangzi.ni@gmail.com
*/
#ifndef RICEMYOFILAMENT_HPP
#define RICEMYOFILAMENT_HPP

#include <cmath>
#include <algorithm>    // std::max
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "RiceMyofilamentConstants.hpp"
#include "ExplicitSolver.hpp"

#define sign(a) ((a) < (0.) ? (-1.0) : (1.0))
#define heav(a) ((a) < (0.) ? (0.0) : (1.0))

class RiceMyofilament
{
public:
    RiceMyofilament(double usr_SLset = 2.2, double usr_SLrest = 1.9);
    ~RiceMyofilament(){};
    void InitialiseMyofilamentStates();
    double MyofilamentODE( const double Cai);
    void FESolve(const double dt);
    void SetSLrest(double usr_SLrest);
    void SetSLset(double usr_SLset);
    double GetStrain();
    void OutputState();
    double SolveForceSingleTimeStep(const double Cai, const double usr_dt);
	double MyofilamentODETemperature(const double Cai);  // for future use, ie fitting experimental
    double N, dN;   //
    double XBprer, dXBprer;   //
    double XBpostr, dXBpostr ;  //
    double SL, dSL;  //
    double xXBpostr, dxXBpostr;  //
    double xXBprer, dxXBprer ;   //
    double TRPNCaL, dTRPNCaL ;  // Low-affinity Ca-bound troponin (uM)
    double TRPNCaH, dTRPNCaH ;  // High-affinity Ca-bound troponin (uM)
    double intf,   dintf;  // integrateddouble  force

    double SLset ;         // initial length
    double SLrest;       //   (um) rest SL length for 0 passive force
    double Force, NormalizedForce;
    double EngineeringStrain;
};



#endif