#ifndef CNZ_CPP
#define CNZ_CPP
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "CNZCell.hpp"
#include "CNZ_Constants.hpp"
#include "stimulus.h"


CNZ_Cell::CNZ_Cell( std::vector<double> & vec, AtriaCellType cell_type , MuationType mutation_type, AFType AF_model, bool usr_IsSAC)
    : UseSAC(usr_IsSAC),
      m_celltype(cell_type),
      para(cell_type, mutation_type, AF_model),
      force_model()
{

    InitialiseElectricalStates();
    m_t = 0.0;
    InitialiseElectricalStatesFromVector(vec);
    Vm = state[0];
    // Vmo = state[0];

}

CNZ_Cell::CNZ_Cell(AtriaCellType cell_type, MuationType mutation_type, AFType AF_model, bool usr_IsSAC, bool UseICFromFile)
    : UseSAC(usr_IsSAC),
      m_celltype(cell_type),
      para(cell_type, mutation_type, AF_model),
      force_model()
{
    InitialiseElectricalStates();
    m_t = 0.0;

    if (UseICFromFile) {
        auto outICs_name = std::string("ICs/") + GetCellTypeToString(cell_type) + std::string("__") +  GetMutationTypeToString(mutation_type) + "_ICs.dat";
        InitialiseElectricalStatesFromFile(outICs_name.c_str());
    }
    Vm = state[0];
    // Vmo = state[0];
}



double CNZ_Cell::SolveElectricalModel(const double dt) {
    double Ena, Ek, Eca, Ekf, Enaf;
    double alpha, beta, tau, inf, a, b;
    double INa, IKr, IKs, ICaL, IK1, Iab, IbK, IbCa;
    double IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT;
    double fnak, sigma;
    double naidot, kidot, caidot;
    double V         = state[0];
    double m         = state[1];
    double h         = state[2];
    double j         = state[3];
    double d         = state[4];
    double f         = state[5];
    double xr        = state[6];
    double xs        = state[7];
    double nai       = state[8];
    double cai       = state[9];
    double ki        = state[10];
    double fca       = state[11];
    double itr       = state[12];
    double its       = state[13];
    double isusr     = state[14];
    double isuss     = state[15];
    double Cass      = state[16];
    double CaSR1     = state[17];
    double CaSR2     = state[18];
    double SERCACa   = state[19];
    double SERCACass = state[20];
    double RyRoss    = state[21];
    double RyRcss    = state[22];
    double RyRass    = state[23];
    double RyRo3     = state[24];
    double RyRc3     = state[25];
    double RyRa3     = state[26];
    double dd        = state[27];
    double ff        = state[28];
    double rkv       = state[29];
    double skv       = state[30];
    double kif       = state[31];
    double naif      = state[32];
    double Vmf       = state[33];
    double CNZ_a     = state[34];
    double CNZ_i     = state[35];
    double If_y      = state[38];

    // incorporation of INa state dependent drug effects
    double BA = state[39];   // Blockade of INa channel,
    double BI = state[40];    // blockade of INa channel.

    double dBAdt = 0.0;//0.0*para.drug_INa_Ka * para.drug_INa_concen * m*m*m*h*j * ( 1 - BA - BI) - 0.0*para.drug_INa_La*BA ;// * exp(-z*V*F/R/T)
    double dBIdt = 0.0;//0.0*para.drug_INa_Ki * para.drug_INa_concen * (1-h) * (1- BA - BI) - 0.0*para.drug_INa_Li * BI; // * exp(-z*V*F/R/T)
    state[39] += 0.0;//dt * dBAdt;
    state[40] = 0.0;// += 0.0;//dt * dBIdt;

    // Reversal Potentials0.0;//
    Ena = 26.71 * log(CRN_nac / nai);
    Ek = 26.71 * log(CRN_kc / ki);
    Eca = 13.35 * log(CRN_cac / cai);
    Ekf = 26.71 * log(kof / kif);
    Enaf = 26.71 * log(naof / naif);

    //INa
    INa = para.GNa * Cm * CRN_gna * (1 - BA - BI) * m * m * m * h * j * (V - Ena);
    // INa = para.GNa * Cm * CRN_gna * m * m * m * h * j * (V - Ena);

    // m gate
    alpha =  0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13)));
    if (fabs(V + 47.13) < 1e-10) alpha = 3.2;
    beta = 0.08 * exp(-V / 11.0);

    tau = 1 / (alpha + beta);
    inf = alpha * tau;

    m = inf + (m - inf) * exp(-dt / tau); // steady state approx due to fast tau

    // j gate
    if (V >= -40.0)
    {
        alpha  = 0.0;
        beta = 0.3 * exp(-2.535e-7 * V) / (1 + exp(-0.1 * (V + 32)));
    }
    else
    {
        alpha = (-1.2714e5 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) * (V + 37.78) / (1 + exp(0.311 * (V + 79.23)));
        beta = 0.1212 * exp(-0.01052 * V) / (1 + exp(-0.1378 * (V + 40.14)));
    }
    tau = 1 / (alpha + beta);
    inf = alpha * tau;
    j = inf + (j - inf) * exp(-dt / tau);

    // h gate
    if (V >= -40.0)
    {
        alpha  = 0.0;
        beta = 1 / (0.13 * (1 + exp((V + 10.66) / -11.1)));
    }
    else
    {
        alpha = 0.135 * exp((V + 80) / -6.8);
        beta = 3.56 * exp(0.079 * V) + 3.1e5 * exp(0.35 * V);
    }
    tau = 1 / (alpha + beta);
    inf = alpha * tau;

    h = inf + (h - inf) * exp(-dt / tau);

    // CRN IKr

    IKr = para.GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4));
    a = 0.0003 * (V + 14.1) / (1 - exp((V + 14.1) / -5));
    b = 0.000073898 * (V - 3.3328) / (exp((V - 3.3328) / 5.1237) - 1);
    if (fabs(V + 14.1) < 1e-10) a = 0.0015;
    if (fabs(V - 3.3328) < 1e-10) b = 3.7836118e-4;  //
    tau = 1 / (a + b);
    inf = 1 / (1 + exp((V + para.IKr_ac_shift + 14.1) / (-6.5 * para.IKr_ac_grad)));
    xr = inf + (xr - inf) * exp(-dt / tau);
    // end CRN IKr
    // IKs
    IKs = para.GKs * Cm * CRN_gks * xs * xs * (V - Ek);

    // xs gate
    a = 0.00004 * (V - 19.9) / (1 - exp((V - 19.9) / -17));
    b = 0.000035 * (V - 19.9) / (exp((V - 19.9) / 9) - 1);
    if (fabs(V - 19.9) < 1e-10) /* denominator = 0 */
    {
        a = 0.00068;
        b = 0.000315;
    }
    tau = 0.5 / (a + b); // note lagrer taus may be more accurate
    inf = sqrt(1 / (1 + exp((V - 19.9 - para.IKs_shift) / (-12.7 * para.IKs_grad))));
    xs = inf + (xs - inf) * exp(-dt / tau);

    //ICaL
    ICaL = 2.125 * para.GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL);

    // fca gate
    inf = 1 / (1 + (Cass / 0.00035));
    tau = 2.0;
    fca = inf + (fca - inf) * exp(-dt / tau);

    // d gate
    a = 1 / (1 + exp((V + 10) / -6.24));
    tau = a * (1 - exp((V + 10) / -6.24)) / (0.035 * (V + 10));
    if (fabs(V + 10) < 1e-10) tau = a * 4.579;
    inf = 1 / (1 + exp((V + 10) / -8));
    d = inf + (d - inf) * exp(-dt / tau);

    // f gate
    inf = exp(-(V + 28) / 6.9) / (1 + exp(-(V + 28) / 6.9));
    tau = 1.5 * 2 * 3 / (0.0197 * exp(-0.0337 * 0.0337 * (V + 10) * (V + 10)) + 0.02);
    f = inf + (f - inf) * exp(-dt / tau);

    //Ito
    Ito = para.Gto * Cm * MAL_gto * itr * its * (V - Ek);

    // r gate
    inf = 1.0 / (1.0 + exp((V - 1.0) / -11.0));
    tau = (0.0035 * exp(-(V / 30.0) * 2) + 0.0015);
    itr = inf + (itr - inf) * exp(-(dt / 1000) / tau);

    // s gate
    inf = 1.0 / (1.0 + exp((V + 40.5) / 11.5));
    tau = (0.025635 * exp (-((V + 52.45) / 15.8827) * 2) + 0.01414);
    its = inf + (its - inf) * exp(-(dt / 1000) / tau);

    //isusr
    inf = 1 / (1 + exp(-(V + para.IKur_ac_shift + 6) / (8.6 * para.IKur_ac_grad)));
    tau = (0.009 / (1.0 + exp((V + 5) / 12)) + 0.0005);
    isusr = inf + (isusr - inf) * exp(-(dt / 1000) / tau);

    //isuss
    inf = 1 / (1 + exp((V + para.IKur_inac_shift + 7.5) / (10 * para.IKur_inac_grad)));
    tau = (0.59 / (1 + exp((V + 60.0) / 10)) + 3.05);
    isuss = inf + (isuss - inf) * exp(-(dt / 1000) / tau);

    //IKur
    if (para.m_mutation_type == D322H)
    {
        double IKur_K0 = -8.26597;
        double IKur_c = 4.5128;
        double IKur_x0 = 1.899769;
        double IKur_y0 = 20.5232;

        IKur_c = 3.6887;
        IKur_x0 = 2.84400335;
        IKur_y0 = 15.2672201;
        IKur = Cm * CNZ_gkur * (IKur_c + IKur_x0 / (1.0 + exp((V - IKur_y0) / (IKur_K0)))) * CNZ_a * CNZ_i * (V - Ek);
        inf = ((1.0) / (1 + exp((V + 10.22) / (-6.81))) );
    } else {

        IKur = /*((1 - para.IKur_type_CNZ) * ( Cm * MAL_gkur * isusr * isuss * (V - Ek))) +
           (para.IKur_type_CNZ) **/
            Cm * CNZ_gkur   * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597))))
            * CNZ_a * CNZ_i * (V - Ek);
        inf = ((para.IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + para.IKur_ac1_shift) / (-5.75418 * para.IKur_ac1_grad))) )
              * ((para.IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + para.IKur_ac2_shift) / (-11.51037561 * para.IKur_ac2_grad))))
              + para.IKur_ac_add;
    }    
    IKur = para.para_IKur_cond * para.GKur * IKur;
    double K_Q10 = 3.5308257834747638;


    inf = 1 / (1 + exp(-(V - (-6 + para.para_IKur_Vhchange)) / (8.6 * para.para_IKur_slope)));
    // para_IKur_Vhchange, para_IKur_slope, 
    //CNZ_a


    tau =  (45.6666746826 / (1 + exp((V + 11.2306497073) / 11.5254705962)) + 4.26753514993)
           * (0.262186042981 / (1 + exp((V + 35.8658312707) / (-3.87510627762))) + 0.291755017928); //
    tau = para.para_IKur_timeconstants*tau / K_Q10 ;

    CNZ_a = inf + (CNZ_a - inf) * exp(-(dt) / tau);

    //CNZ_i
    inf = (para.IKur_inac_mult * 0.52424)
          / (1.0 + exp((V + 15.1142 + para.IKur_inac_shift) / (7.567021 * para.IKur_inac_grad)))
          + 0.4580778 + para.IKur_inac_add;
    tau = 2328 / (1 + exp(((V) - 9.435) / (3.5827))) + 1739.139;
    tau = tau / K_Q10;

    CNZ_i = inf + (CNZ_i - inf) * exp(-(dt) / tau);

    // If
    If = 0.0; //para.Gf * Cm * CRN_gf * If_y * (V - (-22));

    //If y gate
    // inf = 1 / (1 + exp((V + 90.95 + para.If_vshift) / (10.1 * para.If_grad))   );
    // tau = 1 / (1.2783E-04 * exp(-V / 9.2424) + 121.6092 * exp(V / 9.2424)   );

    // If_y = inf + (If_y - inf) * exp(-dt / tau);

    //ICaT
    ICaT = 0.0;// para.GCaT * Cm * gcaT * ff * dd * (V - EcaT);

    // dd gate
    // a = 1.0680 * exp((V + 26.3) / 30.0);
    // b  = 1.0680 * exp(-(V + 26.3) / 30.0);
    // tau = 1.0 / (a + b);
    // inf = 1.0 / (1.0 + exp(-(V + 37.0) / 6.8));
    // dd = inf + (dd - inf) * exp(-dt / tau);

    // // ff gate
    // a = 0.0153 * exp(-(V + 71.0) / 83.3);
    // b  = 0.0150 * exp((V + 71.0) / 15.38);
    // tau = 1.0 / (a + b);
    // inf = 1.0 / (1.0 + exp((V + 71.0) / 9.0));
    // ff = inf + (ff - inf) * exp(-dt / tau);

    //IK1
    IK1 = para.GK1 * Cm * CRN_gk1 * (V - Ek + para.IK1_v_shift)
          / (1 + exp(0.07 * (V + 80.0 + para.IK1_v_shift)));

    //Iab
    Iab = Cm * 0.0003879 * (V + 69.6) / (1 - 0.8377 * exp((V + 49.06) / 1056));

    //Ibk
    IbK = Cm * CRN_gbk * (V - Ek);

    //IbCa
    IbCa = para.Gbca * Cm * CRN_gbca * (V - Eca);

    // IbNa
    IbNa = Cm * CRN_gbna * (V - Ena);

    //ICap
    ICap = 1.4 * para.GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);

    //INaCa
    INaCa = para.GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
            (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
            (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
             exp(V * (CRN_gammalr - 1) * F / (R * T)));

    //INaK
    sigma = (exp(CRN_nac / 67.3) - 1) / 7.0;
    fnak = 1 / (1 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0365 * sigma * exp(-V * F / (R * T)));
    INaK = Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(CRN_kmnai / nai, 1.5));
    double IGap = 0.0;

#ifdef FIBROSIS

    if (FB_number != 0) {

        double Ikv, gkv, infr, infs, taur, taus;
        double IK1f, gk1f, ak1f, bk1f;
        double INaKf, INaKf_max, Vrevf, NaKBf, kmKf, kmNaf;
        double Ibnaf, gbnaf;
        double Iftot;

        //IKv
        gkv = FB_Gkv * 0.25;

        Ikv = Cmf * gkv * rkv * skv * (Vmf - Ekf);

        infr = (1 / (1 + exp(-(Vmf + 20 - 25 + FB_Ikv_shift) / 11)));
        infs = 1 / (1 + exp((Vmf + 23 - 25 + FB_Ikv_shift) / 7));

        taur = 1000 * (0.0203 + 0.1380 * exp(pow(-((Vmf + 20) / 25.9), 2))); // Check which units this is in, likely seconds and we want ms
        taus = 1000 * (1.574 + 5.268 * exp(pow(-((Vmf + 23) / 22.7), 2))); // ^^ hence why multiplied by 1000

        rkv = infr + (rkv - infr) * exp(-dt / taur);
        skv = infs + (skv - infs) * exp(-dt / taus);

        // Ik1
        ak1f = 0.1 / (1 + exp(0.06 * (Vmf - Ekf - 200.0 ) )  );
        bk1f = (3 * exp(0.0002 * (Vmf - Ekf + 100.0)) + exp(0.1 * (Vmf - Ekf - 10.0))) / (1.0 + exp(-0.5 * (Vmf - Ekf)));
        gk1f = 0.4822;
        IK1f = gk1f * (ak1f / (ak1f + bk1f)) * (Vmf - Ekf);

        // INaK
        Vrevf = -150.0;
        NaKBf = -200.0;
        kmKf = 1.0;
        kmNaf = 11.0;

        INaKf_max = para.FB_GNaK * 1.644;

        gbnaf = 0.0095;
        Ibnaf = gbnaf * (Vmf - Enaf);

        Iftot = Cmf * (Ikv + IK1f + INaKf + Ibnaf);

        IGap = GGAP * (V - Vmf);

        Vmf = Vmf - dt * ((Iftot - IGap) / Cmf);

    } // end FB if
    else IGap = 0;
#endif


    // compute ISAC;
    double strain = (force_model.SL - force_model.SLrest) / force_model.SLrest;
    if (UseSAC)
    {
        SAC.Compute_ISAC(strain, V, Cm);
    }
    // Concentrations
    naidot = (-3 * INaK - 3 * INaCa - IbNa - INa - SAC.Isac_Na) / (F * CRN_vi);
    kidot = (2 * INaK - IK1 - Ito - IKur - IKr - IKs - IbK - SAC.Isac_K) / (F * CRN_vi);
    nai = nai + dt * naidot;
    ki = ki + dt * kidot;

    //Caclium handling
    double betass;
    betass = pow(
                 ( 1 + SLlow * KdSLlow / pow(((Cass) + KdSLlow), 2)
                   + SLhigh * KdSLhigh / pow(((Cass) + KdSLhigh), 2)
                   + BCa * KdBCa / pow(((Cass) + KdBCa), 2)  )
                 , (-1));

    double betai, gammai;
    betai = pow(( 1 + BCa * KdBCa / pow((cai + KdBCa), 2)  ), (-1));
    gammai = BCa * KdBCa / pow((cai + KdBCa), 2);

    double betaSR1, betaSR2;
    betaSR1 = pow( ( 1 + CSQN * KdCSQN / pow((CaSR1 + KdCSQN), 2) ), (-1));
    betaSR2 = pow( ( 1 + CSQN * KdCSQN / pow((CaSR2 + KdCSQN), 2) ), (-1));

    double Jj_nj;
    Jj_nj = DCa * Aj_nj / xj_nj * ((Cass) - cai) * 1e-6;

    double J_SERCASR, J_bulkSERCA;
    double J_SERCASRss, J_bulkSERCAss;

    J_SERCASR =  (-k3 * pow(CaSR1, 2) * (cpumps - SERCACa) + k4 * (SERCACa)) * Vnonjunct3 * 2;
    J_bulkSERCA = (k1 * pow(cai, 2) * (cpumps - (SERCACa)) - k2 * (SERCACa)) * Vnonjunct3 * 2;
    J_SERCASRss = (-k3 * pow(CaSR2, 2) * (cpumps - (SERCACass)) + k4 * (SERCACass)) * Vss * 2;
    J_bulkSERCAss = (k1 * pow((Cass), 2) * (cpumps - (SERCACass)) - k2 * (SERCACass)) * Vss * 2;

    double RyRtauadapt = 1;

    double RyRtauactss = 5e-3;
    double RyRtauinactss = 15e-3;
    double RyRtauact = 18.75e-3;
    double RyRtauinact = 87.5e-3;

    double nuss = 625 * Vss;
    double RyRSRCass = (1 - 1 / (1 +  exp((CaSR2 - 0.3) / 0.1)));
    double RyRainfss = 0.505 - 0.427 / (1 + exp((( Cass + (para.fRyR * Cass) ) * 1000 - 0.29) / 0.082));
    double RyRoinfss = (1 - 1 / (1 +  exp(((Cass + (para.fRyR * Cass)   ) * 1000 - ((RyRass) + 0.22)) / 0.03)));
    double RyRcinfss = (1 / (1 + exp(((Cass + (para.fRyR * Cass )) * 1000 - ((RyRass) + 0.02)) / 0.01)));
    double Jrelss = nuss * ( (RyRoss) ) * (RyRcss) * RyRSRCass * ( CaSR2 -  (Cass) );

    double nu3 = 1 * Vnonjunct3;
    double RyRSRCa3 = (1 - 1 / (1 +  exp((CaSR1 - 0.3) / 0.1)));
    double RyRainf3 =  0.505 - 0.427 / (1 + exp(( (cai + ( para.fRyR * cai)  ) * 1000 - 0.29) / 0.082));
    double RyRoinf3 = (1 - 1 / (1 +  exp(( (cai + ( para.fRyR * cai) ) * 1000 - ((RyRa3) + 0.22)) / 0.03)));
    double RyRcinf3 = (1 / (1 +  exp(( (cai + (para.fRyR * cai ) ) * 1000 - ((RyRa3) + 0.02)) / 0.01)));
    double Jrel3 = nu3 * ( (RyRo3) ) * (RyRc3) * RyRSRCa3 * ( CaSR1 -  cai );

    Jrelss = para.fIRel * Jrelss;
    Jrel3 = para.fIRel * Jrel3;

    double JSRCaleak3 = para.GSR_leak * kSRleak * ( CaSR1 - cai ) * Vnonjunct3;
    double JSRCaleakss = para.GSR_leak * kSRleak * ( CaSR2 - (Cass) ) * Vss;

    double O_JSERCASR = J_SERCASR;
    double O_JSERCASRss = J_SERCASRss;
    double O_JBULKSERCA = J_bulkSERCA;
    double O_JBULKSERCAss = J_bulkSERCAss;
    double O_JRELss = Jrelss;
    double O_JREL = Jrel3;
    double O_JSRLEAK = JSRCaleak3;
    double O_JSRLEAKss = JSRCaleakss;

    double JCa, JCass;
    JCa = -para.BULK_CONST * J_bulkSERCA + JSRCaleak3 + Jrel3 + Jj_nj;
    JCass = -Jj_nj + JSRCaleakss - para.BULK_CONST * J_bulkSERCAss + Jrelss;

    double JSRCa1, JSRCa2;
    JSRCa1 = J_SERCASR - JSRCaleak3 - Jrel3;
    JSRCa2 = J_SERCASRss - JSRCaleakss - Jrelss;

    double dy;

    dy = 0.5 * (-J_SERCASR + J_bulkSERCA) / Vnonjunct3;
    SERCACa = Foward_Euler(dt / 1000, dy, SERCACa);
    dy = 0.5 * (-J_SERCASRss + J_bulkSERCAss) / Vss;
    SERCACass = Foward_Euler(dt / 1000, dy, SERCACass);

    RyRoss = Euler_inf(dt / 1000, RyRoss, RyRoinfss, RyRtauactss);
    RyRcss = Euler_inf(dt / 1000, RyRcss, RyRcinfss, RyRtauinactss);
    RyRass = Euler_inf(dt / 1000, RyRass, RyRainfss, RyRtauadapt);
    RyRo3 = Euler_inf(dt / 1000, RyRo3, RyRoinf3, RyRtauact);
    RyRc3 = Euler_inf(dt / 1000, RyRc3, RyRcinf3, RyRtauinact);
    RyRa3 = Euler_inf(dt / 1000, RyRa3, RyRainf3, RyRtauadapt);

    dy =  betass * ( JCass / Vss + ((-( para.RyR * ICaL) - IbCa - ICap - ICaT + 2 * INaCa - SAC.Isac_Ca)) / (2 * Vss * 1000 * F) );
    Cass = Foward_Euler(dt / 1000, dy, Cass);

    dy = JCa / Vnonjunct3 * betai;
    cai = Foward_Euler(dt / 1000, dy, cai);

    dy =  betaSR1 * (DCaSR) * ( (CaSR2 - 2 * CaSR1 + CaSR1) / (dx * dx) + (CaSR1 - CaSR2) / (2 * 3 * (dx * dx)) ) + JSRCa1 / VSR3 * betaSR1;

    CaSR1 = Foward_Euler(dt, dy, CaSR1);

    dy = betaSR2 * (DCaSR) * ( (CaSR2 - 2 * CaSR2 + CaSR1) / (dx * dx) + (CaSR2 - CaSR1) / (2 * 4 * (dx * dx)) ) + JSRCa2 / VSR4 * betaSR2;

    CaSR2 = Foward_Euler(dt, dy, CaSR2);

    // return state variables
    state[1] = m;
    state[2] = h;
    state[3] = j;
    state[4] = d;
    state[5] = f;
    state[6] = xr;
    state[7] = xs;
    state[8] = nai;  // not updating nai,
    state[9] = cai;
    state[10] = ki;
    state[11] = fca;
    state[12] = itr;
    state[13] = its;
    state[14] = isusr;
    state[15] = isuss;
    state[16] = Cass;
    state[17] = CaSR1;
    state[18] = CaSR2;
    state[19] = SERCACa;
    state[20] = SERCACass;
    state[21] = RyRoss;
    state[22] = RyRcss;
    state[23] = RyRass;
    state[24] = RyRo3;
    state[25] = RyRc3;
    state[26] = RyRa3;
    state[27] = dd;
    state[28] = ff;

    state[29] = rkv;
    state[30] = skv;
    state[31] = kif;
    state[32] = naif;
    state[33] = Vmf;
    state[34] = CNZ_a;
    state[35] = CNZ_i;
    state[38] = If_y;

    // keep log of currents; Fri 15 Apr 2016 11:51:29 BST
    temp1 = INa / Cm;
    temp2 = IKur / Cm;



    m_IKur    = IKur/Cm;
    m_ICaL    = ICaL/Cm;
    m_INCX    = INaCa/Cm;
    m_Jrel_ss = Jrelss;
    m_ISAC    = SAC.ISAC/Cm;

    double Iion_tot = (INa + Ito + IKur + IKr + IKs + ICaL
                       + IK1 + IbK + IbNa + IbCa + Iab + INaCa + INaK
                       + ICap + If + ICaT + SAC.ISAC
#ifdef FIBROSIS
                       + (FB_number * IGap)
#endif
                      ) / Cm;

    return Iion_tot;


} // end Return Iion_tot

void CNZ_Cell::SingleElectricalTimeStep(const double usr_stim, const double usr_dt) {
    double Itot = SolveElectricalModel(usr_dt);
    Vmo = Vm;
    state[0] = state[0] - usr_dt * ( Itot + usr_stim);  // note double & Vm = state[0];
    Vm = state[0];
    // state[0] = Vm;   // update Vm in state
}



void CNZ_Cell::InitialiseElectricalStates() {
    // initial conditions for static variant
    state[0] = -76.079842; // V
    state[1] = 0.006676; // m
    state[2] = 0.896736; // h
    state[3] = 0.918836; // j
    state[4] = 0.000259; // d
    state[5] = 0.999059; // f
    state[6] = 0.000072; // xr
    state[7] = 0.022846; // xs
    state[8] = 11.170000; // nai
    state[9] = 0.000098;  // Cai
    state[10] = 115.632438; // ki
    state[11] = 0.770253; // fca
    state[12] = 0.000905; // Itr
    state[13] = 0.956638; // Its
    state[14] = 0.000289; // Isusr
    state[15] = 0.998950; // Isuss
    state[16] = 0.000104; // Cass
    state[17] = 0.437859; // CaSR1
    state[18] = 0.434244; // CaSR2
    state[19] = 0.002432; // SERCACa
    state[20] = 0.002443; // SERCACass
    state[21] = 0.000412; // RyRoss
    state[22] = 0.967156; // RyRcss
    state[23] = 0.118222; // RyRass
    state[24] = 0.000365; // RyRo3
    state[25] = 0.977008; // RyRc3
    state[26] = 0.115451; // RyRa3
    state[27] = 0.003182; // dd
    state[28] = 0.637475; // ff
    // FB data
    state[29] = 0.011694; // rkv
    state[30] = 0.996878; // skv
    state[31] = 129.434900; // kif
    state[32] = 8.554700; // naif
    state[33] = -43.806331; // vmf

    state[34] = 0; // CNZ_a
    state[35] = 1; // CNZ_i
    state[36] = 0; //ACh_i
    state[37] = 0; // ACh_j
    state[38] = 0; // if y
    state[39] = 0; // BA in INa state-dependent block
    state[40] = 0; // BI in INa state-dependent block
}



void CNZ_Cell::InitialiseElectricalStatesFromVector(std::vector<double> vec) {

    int i = 0;
    for ( i = 0; i < 41; ++i)
    {
        state[i] = vec[i];
    }
    i = 40;
    force_model.N        = vec[41];
    force_model.XBprer   = vec[42];
    force_model.XBpostr  = vec[43];
    force_model.SL       = vec[44];
    force_model.xXBpostr = vec[45];
    force_model.xXBprer  = vec[46];
    force_model.TRPNCaL = vec[47];
    force_model.TRPNCaH  = vec[48];
    force_model.intf     = vec[49];
    // std::cout << force_model.N  << std::endl;
    // std::cout << force_model.intf  << std::endl;
}



void CNZ_Cell::InitialiseElectricalStatesFromFile(const char *file_name) {
    std::fstream readin(file_name, std::ios_base::in);
    if (!readin.is_open())
    {
        std::cerr << "Failed to open initial data file " << file_name << std::endl;
        std::exit(0);
    }
    int i = 0;

    for (int i = 0; i < 41; ++i)
    {
        readin >> state[i];
    }

    // readin initial conditions from the same file
    readin >> force_model.N
           >> force_model.XBprer
           >> force_model.XBpostr
           >> force_model.SL
           >> force_model.xXBpostr
           >> force_model.xXBprer
           >> force_model.TRPNCaL
           >> force_model.TRPNCaH
           >> force_model.intf ;
    readin.close();
}

void CNZ_Cell::OutputInitialConditions(const char *file_name) {
    std::fstream outfile(file_name, std::ios_base::out);
    if (!outfile.is_open())
    {
        std::cerr << "Failed to open initial data file " << file_name << std::endl;
        std::exit(0);
    }

    for (int i = 0; i < 41; ++i)
    {
        outfile << state[i] << std::endl;
    }


    // outfile initial conditions from the same file
    outfile << force_model.N    << std::endl
            << force_model.XBprer  << std::endl
            << force_model.XBpostr  << std::endl
            << force_model.SL  << std::endl
            << force_model.xXBpostr  << std::endl
            << force_model.xXBprer  << std::endl
            << force_model.TRPNCaL  << std::endl
            << force_model.TRPNCaH  << std::endl
            << force_model.intf << std::endl ;
    outfile.close();

}



void CNZ_Cell::SingleTimeStep(const double usr_stim, const double usr_dt) {
    SingleElectricalTimeStep(usr_stim, usr_dt);
    double Cai = state[9];
    double dCai_feeback = force_model.SolveForceSingleTimeStep(Cai, usr_dt);  // dCai_feeback for future use
    m_t += usr_dt;   // member t
}



void CNZ_Cell::RunTest(int Beats, double BCL, double usr_dt) {
    double Total_time = Beats * (BCL);
    unsigned long long count = 0;
    double dt;
    dt = usr_dt;
    std::ofstream out("output.dat", std::ofstream::out);
    double t;
    for (t = 0.0; t <= Total_time; t += dt)
    {
        double istim = S1(0.0, -20.0, BCL, t, 2.0);

        SingleTimeStep(istim, usr_dt);
        count ++;
        if (count % 100 == 0) {
            // Output_Voltage();
            out << t << " "
                << Vm << " "
                << state[9] << " "
                << force_model.SL  << " "
                << force_model.GetStrain() << std::endl;     // strain is computed as engineering strain: (SL - SL_ini) / SL_ini
        }
    }
    out.close();
}




#endif