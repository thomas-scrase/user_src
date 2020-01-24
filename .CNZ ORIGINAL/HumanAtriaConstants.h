/* static constants */


#ifndef HUMAN_ATRIA_CONSTANTS_H
#define HUMAN_ATRIA_CONSTANTS_H

// #include <math.h>
#include <cmath>

static const double ISO = 0.0;
static const double AF = 0;
// static const double RA = 0;
//-------------------------------------------------------------------------
// Computed variables
//-------------------------------------------------------------------------
// Enviromental states

static const double Frdy = 96485.0;
static const double R = 8314.0;					//defined in CNZ_Constants.hpp
static const double Temp = 310.0;				//defined in CNZ_Constants.hpp
static const double FoRT = Frdy / (R * Temp);
static const double Cmem = 1.10e-10;  // was 1.10e-10 here by haibo
// cell geometry
static const double cellRadius = 10.25;
static const double cellLength = 100.0;
static const double Vcell = 3.14159265358979 * pow(cellRadius, 2.0) * cellLength * 1.0e-15;   // units in liter L
static const double Vjunc = 0.0539 * 0.01 * Vcell;  // from orginal grandi et al model
// static const double Vjunc = 0.02 * Vcell; // model from ORd model.
static const double Vmyo = 0.65 * Vcell;
static const double Vsr = 0.035 * Vcell; // Haibo, change 29/05/2013
static const double Vsl = 0.02 * Vcell;
// the values in the non-C code. by haibo
static const double DcaJuncSL = 1.64e-6;       // in cm
static const double DcaSLcyto = 1.22e-6;
static const double DnaJuncSL = 1.09e-5;
static const double DnaSLcyto = 1.79e-5;

// comment out by haibo, 29. 05. 2013
/*DcaJuncSL = 1.2205e-5;
DcaSLcyto = 2.8914e-6;
DnaJuncSL = 2.7121e-7;
DnaSLcyto = 1.2722e-6;
*/


static const double junctionLength = 160.0e-3;
static const double junctionRadius = 15.0e-3;
static const double distSLcyto = 0.45;
static const double distJuncSL = 0.5;
static const double SAjunc = 20150.0 * 3.14159265358979 * 2.0 * junctionLength * junctionRadius; // Haibo
static const double SAsl = 3.14159265358979 * 2.0 * cellRadius * cellLength; // Haibo

// fractional parameters:

static const double Fjunc = 0.11;
static const double Fjunc_CaL = 0.9;
static const double Fsl_CaL = 1.0 - Fjunc_CaL;
static const double Fsl = 1.0 - Fjunc;

// fixed ion concerntrations.
static const double Cli = 15;   // Intracellular Cl  [mM]
static const double Clo = 150;  // Extracellular Cl  [mM]
static const double GB_Ko  = 5.4;   // Extracellular K   [mM]
static const double GB_Nao = 140;  // Extracellular Na  [mM]
static const double GB_Cao = 1.8;  // Extracellular Ca  [mM]
static const double Mgi = 1.0;    // Intracellular Mg  [mM]


// Na transport parameters
static const double GB_GNa   = 23.0 * (1.0 - 0.1*AF);  // [mS/uF]
static const double GB_GNab  = 0.597e-3;    // [mS/uF]
static const double IbarNaK  = 1.26;     // [uA/uF]    // 1.8*0.7
static const double KmNaip   = 11.0 * (1.0 - 0.25 * ISO);         // [mM]11 by haibo, from matlab code
static const double KmKo     = 1.5;         // [mM]1.5
static const double Q10NaK   = 1.63;
static const double Q10KmNai = 1.39;

// K current parameters:
static const double pNaK = 0.01833;
static const double gkp  = 0.002;

// Cl current parameters
static const double GClCa = 0.5 * 0.109625;
static const double KdClCa = 100.0e-3;
static const double GClB = 1.0 * 9.0e-3;

// I_Ca parameters
/*pNa = 0.75e-8;     //   [cm/sec]  was    pNa = 3.375e-9; here, haibo
pCa = 2.7e-4;      //  [cm/sec]
pK = 1.35e-7;      //   [cm/sec]    was    pK = 6.075e-8; by haibo*/
static const double pNa = (1.0 + 0.5*ISO) * (1.0 - 0.5*AF) * 0.75e-8;       // [cm/sec]
static const double pCa = (1.0 + 0.5*ISO) * (1.0 - 0.5*AF) * 2.7e-4;       // [cm/sec]
static const double pK  = (1.0 + 0.5*ISO) * (1.0 - 0.5*AF) * 1.35e-7;        // [cm/sec]
static const double Q10CaL = 1.8;      //

static const double IbarNCX = (1.0 + 0.4*AF) * 3.15;    //  % [uA/uF]4.5 in ventricle
static const double KmCai = 3.59e-3;    // [mM]
static const double KmCao = 1.3;        // [mM]
static const double KmNai = 12.29;      // [mM]
static const double KmNao = 87.5;       // [mM]
static const double ksat = 0.27;        // [none]   was 0.32 here, by haibo
static const double nu = 0.35;          // [none]   // from SK code, not in atrial paper. Haibo.
static const double Kdact = 0.384e-3;   // [mM] was 0.225 e-3 haibo
static const double Q10NCX = 1.57;      // [none]
static const double IbarSLCaP = 0.0471; //  was 0.0673
static const double KmPCa = 0.5e-3;     // [mM]
static const double GCaB = 6.0643e-4;   // [uA/uF]  was    GCaB = 5.513e-4; haibo
static const double Q10SLCaP = 2.35;    // [none]


// SR SR_Fluxes
static const double Q10SRCaP = 2.6;
static const double Vmax_SRCaP = 5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
static const double Kmf = (2.5 - 1.25*ISO) * 0.246e-3;          // [mM] was    Kmf = 2.0*0.246e-3; // Haibo
static const double Kmr = 1.7;               // [mM]L cytosol
static const double hillSRCaP = 1.787;       // [mM]
static const double ks = 25.0;                 // [1/ms]
static const double koCa = 10.0 + 20.0*AF + 10.0 * ISO * (1.0 - AF);               // [mM^-2 1/ms]   %default 10   modified 20
static const double kom = 0.06;              // [1/ms]
static const double kiCa = 0.5;              // [1/mM/ms]
static const double kim = 0.005;             // [1/ms]
static const double ec50SR = 0.45;           // [mM]

// % Buffering parameters
// % koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
static const double Bmax_Naj = 7.561;
static const double Bmax_Nasl = 1.65;
static const double koff_na = 1.0e-3;
static const double kon_na = 0.1e-3;
static const double Bmax_TnClow = 70.0e-3;
static const double koff_tncl =  (1.0 + 0.5 * ISO) * 19.6e-3;  // according to the matlab code
static const double kon_tncl = 32.7;
static const double Bmax_TnChigh = 140.0e-3;

static const double kon_sll = 100.0;
static const double Bmax_SLlowj = 4.6e-3 * Vmyo / Vjunc * 0.1;
static const double koff_sll = 1300.0e-3;
static const double kon_slh = 100.0;
static const double Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1;
static const double koff_slh = 30.0e-3;

static const double Bmax_SLlowsl = 37.4e-3 * Vmyo / Vsl;
static const double Bmax_SLhighsl = 13.4e-3 * Vmyo / Vsl;
static const double kon_tnchca = 2.37;
static const double koff_tnchca = 0.032e-3;
static const double kon_tnchmg = 3.0e-3;
static const double koff_tnchmg = 3.33e-3;
static const double kon_cam = 34.0;
static const double Bmax_CaM = 24.0e-3;
static const double koff_cam = 238.0e-3;
static const double kon_myoca = 13.8;
static const double Bmax_myosin = 140.0e-3;
static const double koff_myoca = 0.46e-3;
static const double kon_myomg = 0.0157;
static const double koff_myomg = 0.057e-3;
static const double kon_sr = 100.0;
static const double Bmax_SR = 19.0 * 0.9e-3;
static const double koff_sr = 60.0e-3;
static const double fcaCaj = 0.0;
static const double fcaCaMSL = 0.0;
static const double Qpow = (Temp - 310.0) / 10.0;




static const double GB_ECl = (1.0 / FoRT * log(Cli / Clo)); //
static const double gks_junc = 0.0035;
static const double gks_sl = 0.0035;

/* epi ventricle parameters */
static const double GtoSlow = 0.0156; // 0.0156 (mS/microF)
static const double GtoFast = 0.165;   // 0.1144  (mS/microF) Haibo
static const double iksmultiply = 1.0;  // gks  0.0035  (mS/microF)
static const double ik1multiply = 1.0; //  gk1 0.35     (mS/microF)
static const double ina_c = 0.0; // ina ss
/**************************************/

static const double kon_csqn = 100.0;
static const double Bmax_Csqn = 140.0e-3 * Vmyo / Vsr;
static const double koff_csqn = 65.0;
static const double MaxSR = 15.0;
static const double MinSR = 1.0;

static const double J_ca_juncsl = 8.2413e-13;
static const double J_ca_slmyo = 3.7243e-12;
static const double J_na_juncsl = 1.8313e-14;
//  J_na_juncsl = 1.8313e-10; // CHANGE: Simulink and C
static const double J_na_slmyo = 1.6386e-12;
static const double Cm = Cmem;									//defined in CNZ_Constants.hpp

static const double Vsr_junc = Vsr * 0.08; // 0.08 of Vsr as used in ORd Model. could try more like CRN model and others see any improvements.
static const double Vsr_network = Vsr - Vsr_junc;  // the rest were assigned as network SR volume.



#endif