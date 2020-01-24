
#include <math.h>
#include "HumanAtriaConstants.h"

class ICaL_class
{
public:
    double ICaL_j, ICaL_sl;
    
};


double  GB_ICaL_vclamp(double *Y, double V, double dt, double Caj, double Casl, ICaL_class & ICaLtemp) {

    // const double Q10CaL = 1.8;      //
    // const double Frdy = 96485.0;
    // const double R = 8314.0;
    // const double Temp = 310.0;
    // const double FoRT = Frdy / (R * Temp);
    // const double GB_Ko  = 5.4;   // Extracellular K   [mM]
    // const double GB_Nao = 140;  // Extracellular Na  [mM]
    // const double GB_Cao = 1.8;  // Extracellular Ca  [mM]
    // const double fcaCaMSL = 0.0;
    // const double Qpow = (Temp - 310.0) / 10.0;
    // const double Fjunc = 0.11;
    // const double Fjunc_CaL = 0.9;
    // const double Fsl_CaL = 1.0 - Fjunc_CaL;
    // const double Fsl = 1.0 - Fjunc;
/*    const double Caj = 1.737475e-4;
    const double Casl = 1.031812e-4;*/

    double pNa = 0.75e-8;       // [cm/sec]
    double pCa = 2.7e-4;       // [cm/sec]
    double pK  = 1.35e-7;        // [cm/sec]
    float fcaCaj = 0.0;

    /* atrial ICaL */
    double dss = 1.0 / (1.0 + exp(-(V + 9.0) / 6.0));
    double taud = 1.0 * dss * (1.0 - exp(-(V + 9.0) / 6.0)) / (0.035 * (V + 9.0));
    double fss = 1.0 / (1.0 + exp((V + 30.0) / 7.0)) + 0.2 / (1.0 + exp((50.0 - V) / 20.0));
    double tauf = 1.0 / (0.0197 * exp(-pow(0.0337 * (V + 25.0), 2.0)) + 0.02);
    // dY[10] = (dss - Y[10]) / taud;   // derivative of d
    // dY[11] = (fss - Y[11]) / tauf;   // derivative of f
    // dY[12] = 1.7 * Y[1] * (1.0 - Y[12]) - 11.9e-3 * Y[12]; // derivative of fca_bj
    // dY[13] = 1.7 * Casl * (1.0 - Y[13]) - 11.9e-3 * Y[13]; // derivative of fca_bsl
    double d_dt = (dss - Y[0]) / taud;   // derivative of d
    double f_dt = (fss - Y[1]) / tauf;   // derivative of f
    double fca_bj_dt = 1.7 * Caj * (1.0 - Y[2]) - 11.9e-3 * Y[2]; // derivative of fca_bj
    double fca_bsl_dt = 1.7 * Casl * (1.0 - Y[3]) - 11.9e-3 * Y[3]; // derivative of fca_bsl


    double ibarca_j = pCa * 4.0 * V * Frdy * FoRT * (0.341 * Caj * exp(2.0 * V * FoRT) - 0.341 * GB_Cao) / (exp(2.0 * V * FoRT) - 1.0);
    double ibarca_sl = pCa * 4.0 * V * Frdy * FoRT * (0.341 * Casl * exp(2.0 * V * FoRT) - 0.341 * GB_Cao) / (exp(2.0 * V * FoRT) - 1.0);
    double ibark = pK * V * Frdy * FoRT * (0.75 * Y[27] * exp(V * FoRT) - 0.75 * GB_Ko) / (exp(V * FoRT) - 1.0);
    double ibarna_j = pNa * V * Frdy * FoRT * (0.75 * Y[31] * exp(V * FoRT) - 0.75 * GB_Nao) / (exp(V * FoRT) - 1.0);
    double ibarna_sl = pNa * V * Frdy * FoRT * (0.75 * Y[32] * exp(V * FoRT) - 0.75 * GB_Nao) / (exp(V * FoRT) - 1.0);

    double ICaL_Ca_j = Fjunc_CaL * ibarca_j * Y[0] * Y[1] * (1.0 - Y[2] + fcaCaj) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
    double ICaL_Ca_sl = Fsl_CaL * ibarca_sl * Y[0] * Y[1] * (1.0 - Y[3] + fcaCaMSL) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
    /***********************/
    Y[0] = Y[0] + d_dt * dt;
    Y[1] = Y[1] + f_dt * dt;
    Y[2] = Y[2] + fca_bj_dt * dt;
    Y[3] = Y[3] + fca_bsl_dt * dt;

    double ICaL_K = ibark * Y[0] * Y[1] * (Fjunc_CaL * (fcaCaj + (1.0 - Y[2])) + Fsl_CaL * (fcaCaMSL + (1.0 - Y[3]))) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
    double ICaL_Na_j = Fjunc_CaL * ibarna_j * Y[0] * Y[1] * (1.0 - Y[2] + fcaCaj) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
    double ICaL_Na_sl = Fsl_CaL * ibarna_sl * Y[0] * Y[1] * (1.0 - Y[3] + fcaCaMSL) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
    // I_CaNa = ICaL_Na_j + ICaL_Na_sl;
    // I_Catot = I_Ca + I_CaK + I_CaNa;
    // ICaL = ICaL_Ca_j + ICaL_Ca_sl + ICaL_K + ICaL_Na_j + ICaL_Na_sl;
    ICaLtemp.ICaL_j = ICaL_Ca_j;
    ICaLtemp.ICaL_sl = ICaL_Ca_sl;
    return ICaL_Ca_j + ICaL_Ca_sl + ICaL_K + ICaL_Na_j + ICaL_Na_sl;
}


void model_initialise(double *Y) {
    Y[0] = (7.175662e-6);  // 10: d (dimensionless) (in I_Ca)
    Y[1] = (1.00)   ; // 11: f (dimensionless) (in I_Ca)
    Y[2] = (2.421991e-2);   // 12: f_Ca_Bj (dimensionless) (in I_Ca)
    Y[3] = (1.452605e-2);   // 13: f_Ca_Bsl (dimensionless) (in I_Ca)
}