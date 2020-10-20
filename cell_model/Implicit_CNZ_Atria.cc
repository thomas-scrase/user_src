#include "Implicit_CNZ_Atria.h"

namespace oomph{

void ImplicitCNZAtria::return_initial_membrane_potential(double &v,
                                                        const unsigned &cell_type)
{
    v = -76.079842; // V
}

bool ImplicitCNZAtria::return_initial_state_variable(const unsigned &n,
                                                    double &v,
                                                    const unsigned &cell_type)
{
    switch(n){
        case 0 : v = 0.006676; // m
                break;
        case 1 : v = 0.896736; // h
                break;
        case 2 : v = 0.918836; // j
                break;
        case 3 : v = 0.000259; // d
                break;
        case 4 : v = 0.999059; // f
                break;
        case 5 : v = 0.000072; // xr
                break;
        case 6 : v = 0.022846; // xs
                break;
        case 7 : v = 11.170000; // nai
                break;
        case 8 : v = 0.000098;  // Cai
                break;
        case 9 : v = 115.632438; // ki
                break;
        case 10 : v = 0.770253; // fca
                break;
        case 11 : v = 0.000905; // Itr
                break;
        case 12 : v = 0.956638; // Its
                break;
        case 13 : v = 0.000289; // Isusr
                break;
        case 14 : v = 0.998950; // Isuss
                break;
        case 15 : v = 0.000104; // Cass
                break;
        case 16 : v = 0.437859; // CaSR1
                break;
        case 17 : v = 0.434244; // CaSR2
                break;
        case 18 : v = 0.002432; // SERCACa
                break;
        case 19 : v = 0.002443; // SERCACass
                break;
        case 20 : v = 0.000412; // RyRoss
                break;
        case 21 : v = 0.967156; // RyRcss
                break;
        case 22 : v = 0.118222; // RyRass
                break;
        case 23 : v = 0.000365; // RyRo3
                break;
        case 24 : v = 0.977008; // RyRc3
                break;
        case 25 : v = 0.115451; // RyRa3
                break;
        case 26 : v = 0.003182; // dd
                break;
        case 27 : v = 0.637475; // ff
                break;
        // FB data
        case 28 : v = 0.011694; // rkv
                break;
        case 29 : v = 0.996878; // skv
                break;
        case 30 : v = 129.434900; // kif
                break;
        case 31 : v = 8.554700; // naif
                break;
        case 32 : v = -43.806331; // vmf
                break;
        case 33 : v = 0; // CNZ_a
                break;
        case 34 : v = 1; // CNZ_i
                break;
        case 35 : v = 0; //ACh_i
                break;
        case 36 : v = 0; // ACh_j
                break;
        case 37 : v = 0; // if y
                break;
        case 38 : v = 0; // BA in INa state-dependent block
                break;
        case 39 : v = 0; // BI in INa state-dependent block
                break;
        case 40 : v = 0.97;  // N
                break;
        case 41 : v = 0.01;   // XBprer
                break;
        case 42 : v = 0.01;   // XBpostr
                break;
        case 43 : v = SLset;   // SL
                break;
        case 44 : v = x_0;   // xXBpostr
                break;
        case 45 : v = 0.0;   // xXBprer
                break;
        case 46 : v = 0.01447254;   // TropCaL
                break;
        case 47 : v = 0.2320947;   // TropCaH
                break;
        case 48 : v = 0.0;   // intergral of force, normallised
                break;
        default : return false;
    }
    return true;
}

double ImplicitCNZAtria::cm(CellState &state)
{
    return 1.0;
}

void ImplicitCNZAtria::fill_in_generic_residual_contribution_cell_base( CellState &state,
                                                                    Vector<double> &residuals,
                                                                    DenseMatrix<double> &jacobian,
                                                                    unsigned flag)
{
    double Ena, Ek, Eca, Ekf, Enaf;
    double alpha, beta, tau, inf, a, b;
    double INa, IKr, IKs, ICaL, IK1, Iab, IbK, IbCa;
    double IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT;
    double ISAC, Isac_Ca, Isac_Na, Isac_K; //Stretch activated channels
    double fnak, sigma;
    double naidot, kidot, caidot;
    double V         = state.get_vm();
    double m         = state.get_var(0,0);
    double h         = state.get_var(0,1);
    double j         = state.get_var(0,2);
    double d         = state.get_var(0,3);
    double f         = state.get_var(0,4);
    double xr        = state.get_var(0,5);
    double xs        = state.get_var(0,6);
    double nai       = state.get_var(0,7);
    double cai       = state.get_var(0,8);
    double ki        = state.get_var(0,9);
    double fca       = state.get_var(0,10);
    double itr       = state.get_var(0,11);
    double its       = state.get_var(0,12);
    double isusr     = state.get_var(0,13);
    double isuss     = state.get_var(0,14);
    double Cass      = state.get_var(0,15);
    double CaSR1     = state.get_var(0,16);
    double CaSR2     = state.get_var(0,17);
    double SERCACa   = state.get_var(0,18);
    double SERCACass = state.get_var(0,19);
    double RyRoss    = state.get_var(0,20);
    double RyRcss    = state.get_var(0,21);
    double RyRass    = state.get_var(0,22);
    double RyRo3     = state.get_var(0,23);
    double RyRc3     = state.get_var(0,24);
    double RyRa3     = state.get_var(0,25);
    double dd        = state.get_var(0,26);
    double ff        = state.get_var(0,27);
    double rkv       = state.get_var(0,28);
    double skv       = state.get_var(0,29);
    double kif       = state.get_var(0,30);
    double naif      = state.get_var(0,31);
    double Vmf       = state.get_var(0,32);
    double CNZ_a     = state.get_var(0,33);
    double CNZ_i     = state.get_var(0,34);
    double ACh_i     = state.get_var(0,35);
    double ACh_j     = state.get_var(0,36);
    double If_y      = state.get_var(0,37);
    // incorporation of INa state dependent drug effects
    double BA        = state.get_var(0,38);   // Blockade of INa channel,
    double BI        = state.get_var(0,39);    // blockade of INa channel.

    double N         = state.get_var(0,40);
    double XBprer    = state.get_var(0,41);
    double XBpostr   = state.get_var(0,42);
    double SL        = state.get_var(0,43);
    double xXBpostr  = state.get_var(0,44);
    double xXBprer   = state.get_var(0,45);
    double TRPNCaL   = state.get_var(0,46);
    double TRPNCaH   = state.get_var(0,47);
    double intf      = state.get_var(0,48);

    double dt = state.get_dt();

    unsigned cell_type = state.get_cell_type();

    std::string AF_model = "NONE";

    bool UseSAC = true;


    //////////////////////////////////////////////////////////////////////////////////////
    //Cell and mutation type parameters

    ////////////////////////////////
    //Declare the parameters
    ////////////////////////////////
    double GNa;
    double tau_ina_act;
    double tau_ina_inact;
    double drug_INa_Ka, drug_INa_La, drug_INa_Ki,drug_INa_Li;
    double drug_IKur_KO,     drug_IKur_ZKO,     drug_IKur_KC,     drug_IKur_ZKC,     drug_IKur_LO,     drug_IKur_ZLO,     drug_IKur_LC,     drug_IKur_ZLC;
    double drug_IKur_concen;
    double drug_INa_concen;

    double INa_Act_Shift;
    double Gto;
    double ITo_vshift;
    double Ito_Act_Shift;
    
    double GCaL;
    double ICaL_inac_shift;
    double ICaL_ac_shift;
    double ICaL_tau_rate;
    
    double GKs;
    double IKs_ac_grad;
    double IKs_shift;
    double IKs_grad;

    double GKr;
    double IKr_ac_grad;
    double IKr_ac_shift;
    double IKr_grad;

    double GK1;
    double IK1_v_shift;

    double Gf;
    double If_grad;
    double If_vshift;
    
    double GNaCa;
    double GNaK;
    double GCaT;
    double GCap;    
    double Gbca;

    double BULK_CONST;
    double GSR_leak;
    double RyR;
    double fIRel;
    double fRyR;

    double GKur;
    double IKur_cnd_add;
    double IKur_cnd_mult;
    double IKur_cnd_shift;
    double IKur_cnd_grad;
    double IKur_ac_shift;
    double IKur_ac_grad;
    double IKur_ac1_mult;
    double IKur_ac1_shift;
    double IKur_ac1_grad;
    double IKur_ac2_mult;
    double IKur_ac2_shift;
    double IKur_ac2_grad;
    double IKur_ac_add;
    double IKur_inac_mult;
    double IKur_inac_shift;
    double IKur_inac_grad;
    double IKur_inac_add;
    double Simple_Ikur_ac_shift;
    double Simple_Ikur_ac_grad;
    double Simple_Ikur_inac_shift;
    double Simple_Ikur_inac_grad;
    double Simple_GKur;

    int IKur_Model_type;
    int IKur_type_CNZ;
    int tau_type;
    float Fcell;
    float Zindex;
    int phase;
    float SAN_Gna;
    float SAN_Gto;
    float SAN_GcaL;
    float SAN_Gks;
    float SAN_Gkr;
    float SAN_Gf;
    float SAN_Gbna;
    float SAN_Gbca;
    int   SAN_periph;
    int state_len;
    int state_len2;
    // AtriaCellType m_cell_type;
    // IKurMuationType m_mutation_type;
    // AFType m_AF_model;

    float para_IKur_cond, para_IKur_Vhchange, para_IKur_slope, para_IKur_timeconstants;

    ////////////////////////////////
    // Default scaling factors
    ////////////////////////////////
    GNa = GCaL = GKur = GKr = GK1 = GKs = Gto =  GNaCa = GNaK = GCap = fIRel = RyR = GSR_leak = GCaT = Gbca = 1.0;
    IKur_ac_grad = IKur_inac_grad = IKs_grad = IKr_ac_grad = 1;
    BULK_CONST = 1.25;
    fRyR = IKur_ac_shift = IKur_inac_shift = IKs_shift =  IK1_v_shift = IKr_ac_shift  =  0.0;
    IKur_ac1_grad = IKur_ac2_grad = 1;
    IKur_ac1_mult = IKur_inac_mult = IKur_ac2_mult = 1;
    IKur_ac_add = IKur_inac_add = 0;
    IKur_ac1_shift = IKur_ac2_shift = 0;
    If_vshift = 0;
    If_grad = 1;
    Gf = 0.4;
    // Istim = 0;
    Fcell = 0;
    Zindex = 0;
    phase = 0;
    SAN_Gna = 0;
    SAN_Gto = 0;
    SAN_GcaL = 0;
    SAN_Gks = 0;
    SAN_Gkr = 0;
    SAN_Gf = 0;
    SAN_Gbna = 0;
    SAN_Gbca = 0;
    SAN_periph = 0;
    IKur_type_CNZ = 1;
    Simple_Ikur_ac_shift = 0;
    Simple_Ikur_ac_grad = 1;
    Simple_Ikur_inac_shift = 0;
    Simple_Ikur_inac_grad = 1;
    Simple_GKur = 1;
    IKr_grad = 1;

    drug_INa_concen = 0.0; // uM
    drug_INa_Ka = 0.1/1000.0; // ms-1*uM-1
    drug_INa_Ki = 0.1/1000.0; // ms-1 * uM -1
    drug_INa_La = 0.1; // ms-1
    drug_INa_Li = 0.01; // ms-1

    drug_IKur_KO = 0;
    drug_IKur_ZKO = 0;
    drug_IKur_KC = 0;
    drug_IKur_ZKC = 0;

    drug_IKur_LO = 0;
    drug_IKur_ZLO = 0;
    drug_IKur_LC = 0;
    drug_IKur_ZLC = 0;

    drug_IKur_concen = 0.0;


    INa_Act_Shift=0;
    Ito_Act_Shift=0.0;


    para_IKur_cond = 1.0;
    para_IKur_Vhchange = 0.0;
    para_IKur_slope = 1.0;
    para_IKur_timeconstants =1.0;



    ////////////////////////////////
    //Set cell type para
    ////////////////////////////////
    // m_cell_type = cell_type;
    If_vshift = 0;
    If_grad = 1;
    Gf = 0.4;

    //      0 RA, 1 PM, 2 CT, 3 RAA, 4 AS, 5 AVR, 6 BB, 7 LA, 8 LAA, 9 PV, 10 SAN_C, 11 SAN_P
    if (cell_type == 0/*RA*/) {
        BULK_CONST = 1.20;//1.25;
    }
    else  if (cell_type == 1/*PM*/) { // PM
        GCaL = 0.94 * GCaL;
    }
    else if (cell_type == 2/*CT*/) { // CT
        GCaL = 1.68;//2;  // note 1.68 = Biophysically based (Feng et al), 2 better replicates APD differences (1.68 is still within expt range)
        Gto = 1.35 * Gto;
    }
    else if (cell_type == 3/*RAA*/ || cell_type == 4/*AS*/) { // RAA / AS
        BULK_CONST = 1.5;
        Gto = 0.53;//0.4; // or 0.53
        if (cell_type == 4/*AS*/) { // AS
            GCaL = 0.4;//0.25;
            GKur = 0.667;
            Gto = 0.4 * 0.4;
            IK1_v_shift = -6;
            GNa = 1.3;
        }
        /*if (IKur_type_CNZ == 0){
            IKur_ac_shift = -14;
            IKur_inac_shift = 25;
            IKur_ac_grad = 2.1317;
            IKur_inac_grad = 1.769;
        }*/
    }
    else if (cell_type == 5/*AVR*/) { // AVR
        GCaL = 0.67;
        GKr = 1.63;
        Gto = 0.6;
    }
    else if (cell_type == 6/*BB*/) { // BB
        GCaL = 2.32;//2.2;//1.6;
        Gto = 1.17;//0.8;
        GK1 = 0.85;
        BULK_CONST = 1.6;//1.22;

    }
    else if (cell_type == 7/*LA*/) { // LA
        Gf = 1;
        GKr = 1.6;

    }
    else if (cell_type == 8/*LAA*/) { // LAA
        Gf = 1;
        GKr = 1.6;
        Gto = 0.53;

        if (IKur_type_CNZ == 0) {
            IKur_ac_shift = -14;
            IKur_inac_shift = 25;
            IKur_ac_grad = 2.1317;
            IKur_inac_grad = 1.769;
        }
        GKur = 0.8;
    }
    /*else if (cell_type == 9PV) { // PV

        IKr_ac_shift = -4.2;
        IKr_ac_grad = 0.938;
        IKs_shift = -9.6;
        IKs_grad = 1.047;
        IK1_v_shift = -15;
        GK1 = 0.65;

    }
    else if (cell_type == 10) { // PV new (rough ehrlich) -- note wrong on Gkr and Gks
        GKr = 0.75 * 1.6;
        Gto = 0.8 * 0.53;
        GKs = 0.7 * 1.8;
        IKs_shift = -13.43;
        IKs_grad = 0.5;//0.75;
        BULK_CONST = 1.3;
        GK1 = 0.7;
        GCaL = 0.65;//0.75;
        IK1_v_shift = -10;
    }*/
    else if (cell_type == 9/*PV*/) { // PV new Ehrlich + Cha (Primarily Ehrlich) - USE THIS PV MODEL
        GKr = 1.5 * 1.6; //1.6*1.6;
        Gto = 0.75;//0.72*0.53;
        GKs = 1.5;//1.75*1.8;
        GCaL = 0.7;
        GK1 = 0.62;//0.69;
    }
    else if (cell_type == 10/*SAN_C*/) { // SAN C
        GNa = 0.06;
        GKur = 0.26;
        GK1 = 0.42;
        GCaT = 4.5 * 17; //5.2*17;
        Gf = 4.4;
        GNaCa = 0.5;
        Gto = 0.4;
        GCaL = 0.68;
        GKs = 0.69;
        GKr = 0.45;
        BULK_CONST = 3;
        If_grad = 1.1;
    }
    else if (cell_type == 11/*SAN_P*/) { // SAN P

        GNa = 0.3;
        GKur = 0.26;
        GK1 = 0.3;//0.42;
        GCaT = 4.5 * 17;
        Gf = 3.93 * 4.4;
        GNaCa = 0.5;
        Gto = 0.8;//0.92;
        GCaL = 2.28;
        GKs = 2;//4.19;
        GKr = 2.75;
        BULK_CONST = 3;
        If_grad = 1.075;
    }


    ////////////////////////////////
    //Set IKur Mutation Para
    ////////////////////////////////
    // m_mutation_type = mutation_type;

    if (*Mutation_pt == "WT") { /*WT*/
        IKur_ac1_mult   = 1.0;
        IKur_ac1_shift  = 0.0;
        IKur_ac1_grad   = 1.0;
        IKur_ac2_mult   = 1.0;
        IKur_ac2_shift  = 0.0;
        IKur_ac2_grad   = 1.0;
        IKur_ac_add     = 0.0;
        IKur_inac_mult  = 1.0;
        IKur_inac_shift = 0.0;
        IKur_inac_grad  = 1.0;
        IKur_inac_add   = 0.0;
        GKur *= 1.0;
    }
    else if (*Mutation_pt == "D322H") {  /*D322H*/
        // IKur_ac1_mult   = 1.0;
        // IKur_ac1_shift  = -9.3932825747999988;
        // IKur_ac1_grad   = 1.3785501205210822;
        // IKur_ac2_mult   = 1.0;
        // IKur_ac2_shift  = 15.229362644299998;
        // IKur_ac2_grad   = 0.25822502847411422;
        // IKur_ac_add     = 0.0;

        // IKur_ac1_shift  = -9.3932825747999988;
        // IKur_ac1_grad   = 1.3785501205210822;
        // IKur_ac2_shift  = 15.229362644299998;
        // IKur_ac2_grad   = 0.25822502847411422;

        IKur_ac1_mult   = 1.0;
        IKur_ac2_mult   = 1.0;
        IKur_ac_add     = 0.0;
        IKur_ac1_shift = -7.6684;
        IKur_ac1_grad = 1.18;
        IKur_ac2_shift = 23.59;
        IKur_ac2_grad = 0.14173;

        
        // Inactivation
        IKur_inac_mult  = 0.87254304449870279;
        IKur_inac_shift = -9.6150807654200001;
        IKur_inac_grad  = 0.80075124328979652;
        IKur_inac_add   = 0.073931087206000057;

        // GKur *= 1.720232;
        // GKur *= 1.820232;
        GKur = 1.74186924; /// lastest D322H Model


        Simple_Ikur_ac_shift += -3.2641;
        Simple_Ikur_ac_grad *= 0.8;
        Simple_Ikur_inac_shift += 9.61;
        Simple_Ikur_inac_grad *= 0.8;
        Simple_GKur *= 1.7995;

    }
    else if (*Mutation_pt == "E48G") { /*E48G*/
        IKur_ac1_mult   = 1.0;
        IKur_ac1_shift  =  -8.914978286109998;
        IKur_ac1_grad   = 1.687979145508135 ;
        IKur_ac2_mult   = 1.0;
        IKur_ac2_shift  = 14.395569738899999;
        IKur_ac2_grad   = 0.312429308790384;
        IKur_ac_add     = 0.0;
        // // Inactivation
        IKur_inac_mult  = 0.935206196986113;
        IKur_inac_shift = -4.02829;
        IKur_inac_grad  = 0.940280253452448;
        IKur_inac_add   = 0.039197146604;

        GKur *= 1.286906;


        Simple_Ikur_ac_shift += -3.02;
        Simple_Ikur_ac_grad *= 0.95;
        Simple_Ikur_inac_shift += 4.02;
        Simple_Ikur_inac_grad *= 0.94;
        Simple_GKur *= 1.32;
    }
    else if (*Mutation_pt == "A305T") {  /*A305T*/
        IKur_ac1_mult   = 1.0;
        IKur_ac1_shift  = -5.135031219399998 ;
        IKur_ac1_grad   = 2.40452882662690 ;
        IKur_ac2_mult   = 1.0;
        IKur_ac2_shift  = 13.643107896399998;
        IKur_ac2_grad   = 0.387170393885174;
        IKur_ac_add     = 0.0;
        // Inactivation
        IKur_inac_mult  = 0.874602534570044;
        IKur_inac_shift = -6.032311251379999;
        IKur_inac_grad  = 0.748490912918043;
        IKur_inac_add   = 0.070656047487;

        GKur *= 1.388228;

        Simple_Ikur_ac_shift += -4.068;
        Simple_Ikur_ac_grad *= 1.102;
        Simple_Ikur_inac_shift += 6.03;
        Simple_Ikur_inac_grad *= 0.74;
        Simple_GKur *= 1.445;
    }
    else if (*Mutation_pt == "Y155C") {   /*Y155C*/
        IKur_ac1_mult   = 1.0;
        IKur_ac1_shift  = -4.430626802899999 ;
        IKur_ac1_grad   = 1.651529981648123 ;
        IKur_ac2_mult   = 1.0;
        IKur_ac2_shift  = 3.3770968185;
        IKur_ac2_grad   =  0.517356087786261;
        IKur_ac_add     = 0.0;
        // Inactivation
        IKur_inac_mult  = 0.936476668495346;
        IKur_inac_shift = -5.012779270500001;
        IKur_inac_grad  = 0.781775743434569;
        IKur_inac_add   = 0.034583585385;

        GKur *= 0.461283;
        Simple_Ikur_ac_shift += -.89;
        Simple_Ikur_ac_grad *= 0.75;
        Simple_Ikur_inac_shift += 5.01;
        Simple_Ikur_inac_grad *= 0.78;
        Simple_GKur *= 0.4745;
    }
    else if (*Mutation_pt == "D469E") {  /*D469E*/

        IKur_ac1_mult   = 1.0;
        IKur_ac1_shift  = 0.0 ;
        IKur_ac1_grad   = 1.0;
        IKur_ac2_mult   = 1.0;
        IKur_ac2_shift  = 0.0;
        IKur_ac2_grad   =  1.0;
        IKur_ac_add     = 0.0;
        // Inactivation
        IKur_inac_mult  =  1.174004519012284;
        IKur_inac_shift = -4.5454272743;
        IKur_inac_grad  = 0.842432589996777;
        IKur_inac_add   = -0.09291551303;

        GKur *= 0.563448;
        
        Simple_Ikur_inac_shift += 4.5455;
        Simple_Ikur_inac_grad *= 0.84;
        Simple_GKur *= 0.5457;
    }
    else if (*Mutation_pt == "P488S") {
        /*P488S*/   /*note, P488s and D469E used the same activation if as the WT, and
                                                        P488S used the inactivation of P488S/WT, due to the absense of
                                                        P488S inactvation data*/

        IKur_ac1_mult   = 1.0;
        IKur_ac1_shift  = 0.0 ;
        IKur_ac1_grad   = 1.0;
        IKur_ac2_mult   = 1.0;
        IKur_ac2_shift  = 0.0;
        IKur_ac2_grad   =  1.0;
        IKur_ac_add     = 0.0;
        // Inactivation
        IKur_inac_mult  =  1.212942827258507;
        IKur_inac_shift = 1.4151025432;
        IKur_inac_grad  = 0.687940393974062;
        IKur_inac_add   = -0.110508855073;

        GKur *= 0.038593;


        Simple_Ikur_ac_grad *= 0.83;
        Simple_Ikur_inac_shift += -1.4251;
        Simple_Ikur_inac_grad *= 0.68;
        Simple_GKur *= 0.0384;

    } else {
        std::cerr << "wrong IKur mutation type \n";
        std::exit(0);
    }

    ////////////////////////////////
    //Set AF type Para
    ////////////////////////////////
    // m_AF_model = AF_model;
    if (AF_model == "NONE") { // SR condition, do not change anything
        GCaL = GCaL;
        Gto = Gto;
        GK1 = GK1;

    } else if (AF_model == "AF1") { // Bosch based model  // based on Henggui's paper. in cardiovascular research 2005.
        GCaL = 0.26 * GCaL;
        Gto = 0.3 * Gto;
        GK1 = 3.35 * GK1;
        Ito_Act_Shift= 16;
        INa_Act_Shift=1.6;
    }
    else if (AF_model == "AF2") { // Workman based model
        GCaL = 0.35 * GCaL;
        Gto = 0.35 * Gto;
        GK1 = 1.75 * GK1;
    }
    else if (AF_model == "AF3") { // Wagoner based model
        GCaL = 0.37 * GCaL;
        Gto = 0.34 * Gto;
        GKur = 0.51 * GKur;
        GK1 = 2.06 * GK1;
    }
    else if (AF_model == "AF4") { // Generic
        GCaL = 0.3 * GCaL;
        Gto = 0.35 * Gto;
        GK1 = 2 * GK1;
        GKur = 0.5 * GKur;
        GNaCa = 1.55 * GNaCa;
        BULK_CONST = 1.5; // Increase according to shanmugam et al, AF_SR_updatke
        GSR_leak = 1.25;
        fIRel = 3 * fIRel;
        GKs = 2 * GKs;
        RyR = 2.5;
    }
    else {
        std::cerr << "wrong AF model type!!! \n";
        std::exit(0);
    }




    //End cell and mutation type parameters
    //////////////////////////////////////////////////////////////////////////////////////


    // Reversal Potentials0.0;//
    Ena = 26.71 * log(CRN_nac / nai);
    Ek = 26.71 * log(CRN_kc / ki);
    Eca = 13.35 * log(CRN_cac / cai);
    Ekf = 26.71 * log(kof / kif);
    Enaf = 26.71 * log(naof / naif);

    //INa
    INa = GNa * Cm * CRN_gna * (1 - BA - BI) * m * m * m * h * j * (V - Ena);
    // INa = GNa * Cm * CRN_gna * m * m * m * h * j * (V - Ena);

    // m gate
    alpha =  0.32 * (V + 47.13) / (1 - exp(-0.1 * (V + 47.13)));
    if (fabs(V + 47.13) < 1e-10) alpha = 3.2;
    beta = 0.08 * exp(-V / 11.0);

    tau = 1 / (alpha + beta);
    inf = alpha * tau;

    // m = inf + (m - inf) * exp(-dt / tau); // steady state approx due to fast tau
    residuals[0] -= state.get_var(1,0) - (inf - m)/tau;

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
    // j = inf + (j - inf) * exp(-dt / tau);
    residuals[2] -= state.get_var(1,0) - (inf - j)/tau;

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
    residuals[1] -= state.get_var(1,1) - (inf - h)/tau;

    // CRN IKr

    IKr = GKr * Cm * CRN_gkr * xr * (V - Ek) / (1 + exp((V + 15) / 22.4));
    a = 0.0003 * (V + 14.1) / (1 - exp((V + 14.1) / -5));
    b = 0.000073898 * (V - 3.3328) / (exp((V - 3.3328) / 5.1237) - 1);
    if (fabs(V + 14.1) < 1e-10) a = 0.0015;
    if (fabs(V - 3.3328) < 1e-10) b = 3.7836118e-4;  //
    tau = 1 / (a + b);
    inf = 1 / (1 + exp((V + IKr_ac_shift + 14.1) / (-6.5 * IKr_ac_grad)));
    // xr = inf + (xr - inf) * exp(-dt / tau);
    residuals[5] -= state.get_var(1,5) - (inf - xr)/tau;
    // end CRN IKr
    // IKs
    IKs = GKs * Cm * CRN_gks * xs * xs * (V - Ek);

    // xs gate
    a = 0.00004 * (V - 19.9) / (1 - exp((V - 19.9) / -17));
    b = 0.000035 * (V - 19.9) / (exp((V - 19.9) / 9) - 1);
    if (fabs(V - 19.9) < 1e-10) /* denominator = 0 */
    {
    a = 0.00068;
    b = 0.000315;
    }
    tau = 0.5 / (a + b); // note lagrer taus may be more accurate
    inf = sqrt(1 / (1 + exp((V - 19.9 - IKs_shift) / (-12.7 * IKs_grad))));
    // xs = inf + (xs - inf) * exp(-dt / tau);
    residuals[6] -= state.get_var(1,6) - (inf - xs)/tau;

    //ICaL
    ICaL = 2.125 * GCaL * Cm * CRN_gcaL * d * f * fca * (V - CRN_ErL);

    // fca gate
    inf = 1 / (1 + (Cass / 0.00035));
    tau = 2.0;
    // fca = inf + (fca - inf) * exp(-dt / tau);
    residuals[10] -= state.get_var(1,10) - (inf - fca)/tau;

    // d gate
    a = 1 / (1 + exp((V + 10) / -6.24));
    tau = a * (1 - exp((V + 10) / -6.24)) / (0.035 * (V + 10));
    if (fabs(V + 10) < 1e-10) tau = a * 4.579;
    inf = 1 / (1 + exp((V + 10) / -8));
    // d = inf + (d - inf) * exp(-dt / tau);
    residuals[3] -= state.get_var(1,3) -(inf - d)/tau;

    // f gate
    inf = exp(-(V + 28) / 6.9) / (1 + exp(-(V + 28) / 6.9));
    tau = 1.5 * 2 * 3 / (0.0197 * exp(-0.0337 * 0.0337 * (V + 10) * (V + 10)) + 0.02);
    // f = inf + (f - inf) * exp(-dt / tau);
    residuals[4] -= state.get_var(1,4) - (inf - f)/tau;

    //Ito
    Ito = Gto * Cm * MAL_gto * itr * its * (V - Ek);

    // r gate
    inf = 1.0 / (1.0 + exp((V - 1.0) / -11.0));
    tau = (0.0035 * exp(-(V / 30.0) * 2) + 0.0015);
    // itr = inf + (itr - inf) * exp(-(dt / 1000) / tau);
    residuals[11] -= state.get_var(1,11) - (inf - itr)/tau;

    // s gate
    inf = 1.0 / (1.0 + exp((V + 40.5) / 11.5));
    tau = (0.025635 * exp (-((V + 52.45) / 15.8827) * 2) + 0.01414);
    // its = inf + (its - inf) * exp(-(dt / 1000) / tau);
    residuals[12] -= state.get_var(1,12) - (inf - its)/tau;

    //isusr
    inf = 1 / (1 + exp(-(V + IKur_ac_shift + 6) / (8.6 * IKur_ac_grad)));
    tau = (0.009 / (1.0 + exp((V + 5) / 12)) + 0.0005);
    // isusr = inf + (isusr - inf) * exp(-(dt / 1000) / tau);
    residuals[13] -= state.get_var(1,13) - (inf - isusr)/tau;

    //isuss
    inf = 1 / (1 + exp((V + IKur_inac_shift + 7.5) / (10 * IKur_inac_grad)));
    tau = (0.59 / (1 + exp((V + 60.0) / 10)) + 3.05);
    // isuss = inf + (isuss - inf) * exp(-(dt / 1000) / tau);
    residuals[14] -= state.get_var(1,14) - (inf - isuss)/tau;
    //IKur



    if (*Mutation_pt == "D322H")
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

    IKur = /*((1 - IKur_type_CNZ) * ( Cm * MAL_gkur * isusr * isuss * (V - Ek))) +
    (IKur_type_CNZ) **/
    Cm * CNZ_gkur   * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597))))
    * CNZ_a * CNZ_i * (V - Ek);
    inf = ((IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + IKur_ac1_shift) / (-5.75418 * IKur_ac1_grad))) )
      * ((IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + IKur_ac2_shift) / (-11.51037561 * IKur_ac2_grad))))
      + IKur_ac_add;
    }


    IKur = para_IKur_cond * GKur * IKur;
    double K_Q10 = 3.5308257834747638;


    inf = 1 / (1 + exp(-(V - (-6 + para_IKur_Vhchange)) / (8.6 * para_IKur_slope)));
    // para_IKur_Vhchange, para_IKur_slope, 
    //CNZ_a


    tau =  (45.6666746826 / (1 + exp((V + 11.2306497073) / 11.5254705962)) + 4.26753514993)
    * (0.262186042981 / (1 + exp((V + 35.8658312707) / (-3.87510627762))) + 0.291755017928); //
    tau = para_IKur_timeconstants*tau / K_Q10 ;

    // CNZ_a = inf + (CNZ_a - inf) * exp(-(dt) / tau);
    residuals[33] -= state.get_var(1,33) - (inf - CNZ_a)/tau;

    //CNZ_i
    inf = (IKur_inac_mult * 0.52424)
    / (1.0 + exp((V + 15.1142 + IKur_inac_shift) / (7.567021 * IKur_inac_grad)))
    + 0.4580778 + IKur_inac_add;
    tau = 2328 / (1 + exp(((V) - 9.435) / (3.5827))) + 1739.139;
    tau = tau / K_Q10;

    // CNZ_i = inf + (CNZ_i - inf) * exp(-(dt) / tau);
    residuals[34] -= state.get_var(1,34) - (inf - CNZ_i)/tau;

    // If
    If = 0.0; //Gf * Cm * CRN_gf * If_y * (V - (-22));

    //If y gate
    // inf = 1 / (1 + exp((V + 90.95 + If_vshift) / (10.1 * If_grad))   );
    // tau = 1 / (1.2783E-04 * exp(-V / 9.2424) + 121.6092 * exp(V / 9.2424)   );

    // If_y = inf + (If_y - inf) * exp(-dt / tau);
    residuals[37] -= state.get_var(1,37);

    //ICaT
    ICaT = 0.0;// GCaT * Cm * gcaT * ff * dd * (V - EcaT);

    // dd gate
    // a = 1.0680 * exp((V + 26.3) / 30.0);
    // b  = 1.0680 * exp(-(V + 26.3) / 30.0);
    // tau = 1.0 / (a + b);
    // inf = 1.0 / (1.0 + exp(-(V + 37.0) / 6.8));
    // dd = inf + (dd - inf) * exp(-dt / tau);
    residuals[26] -= state.get_var(1,26);

    // // ff gate
    // a = 0.0153 * exp(-(V + 71.0) / 83.3);
    // b  = 0.0150 * exp((V + 71.0) / 15.38);
    // tau = 1.0 / (a + b);
    // inf = 1.0 / (1.0 + exp((V + 71.0) / 9.0));
    // ff = inf + (ff - inf) * exp(-dt / tau);
    residuals[27] -= state.get_var(1,27);


    //IK1
    IK1 = GK1 * Cm * CRN_gk1 * (V - Ek + IK1_v_shift)
    / (1 + exp(0.07 * (V + 80.0 + IK1_v_shift)));

    //Iab
    Iab = Cm * 0.0003879 * (V + 69.6) / (1 - 0.8377 * exp((V + 49.06) / 1056));

    //Ibk
    IbK = Cm * CRN_gbk * (V - Ek);

    //IbCa
    IbCa = Gbca * Cm * CRN_gbca * (V - Eca);

    // IbNa
    IbNa = Cm * CRN_gbna * (V - Ena);

    //ICap
    ICap = 1.4 * GCap * Cm * CRN_icapbar * Cass / (Cass + CRN_kmcap);

    //INaCa
    INaCa = GNaCa * Cm * CRN_knacalr / (pow(CRN_kmnalr, 3.0) + pow(CRN_nac, 3.0)) / (CRN_kmcalr + CRN_cac) /
    (1 + CRN_ksatlr * exp((CRN_gammalr - 1) * V * F / (R * T))) *
    (nai * nai * nai * CRN_cac * exp(V * CRN_gammalr * F / (R * T)) - CRN_nac * CRN_nac * CRN_nac * Cass *
     exp(V * (CRN_gammalr - 1) * F / (R * T)));

    //INaK
    sigma = (exp(CRN_nac / 67.3) - 1) / 7.0;
    fnak = 1 / (1 + 0.1245 * exp(-0.1 * V * F / (R * T)) + 0.0365 * sigma * exp(-V * F / (R * T)));
    INaK = Cm * CRN_inakbar * fnak * CRN_kc / (CRN_kc + CRN_kmko) / (1 + pow(CRN_kmnai / nai, 1.5));
    double IGap = 0.0;

    // #ifdef FIBROSIS
    double FB_number = 0.0;
    if (FB_number != 0) {

    double Ikv, gkv, infr, infs, taur, taus;
    double IK1f, gk1f, ak1f, bk1f;
    double INaKf, INaKf_max, Vrevf, NaKBf, kmKf, kmNaf;
    double Ibnaf, gbnaf;
    double Iftot;

    //Some things which weren't declared anywhere
    double FB_Gkv = 0.0;
    double FB_Ikv_shift = 0.0;
    double FB_GNaK = 0.0;
    double GGAP = 0.0;
    ///

    //IKv
    gkv = FB_Gkv * 0.25;

    Ikv = Cmf * gkv * rkv * skv * (Vmf - Ekf);

    infr = (1 / (1 + exp(-(Vmf + 20 - 25 + FB_Ikv_shift) / 11)));
    infs = 1 / (1 + exp((Vmf + 23 - 25 + FB_Ikv_shift) / 7));

    taur = 1000 * (0.0203 + 0.1380 * exp(pow(-((Vmf + 20) / 25.9), 2))); // Check which units this is in, likely seconds and we want ms
    taus = 1000 * (1.574 + 5.268 * exp(pow(-((Vmf + 23) / 22.7), 2))); // ^^ hence why multiplied by 1000

    // rkv = infr + (rkv - infr) * exp(-dt / taur);
    residuals[28] -= state.get_var(1,28) - (infr - rkv)/taur;
    // skv = infs + (skv - infs) * exp(-dt / taus);
    residuals[29] -= state.get_var(1,29) - (infs - skv)/taus;

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

    INaKf_max = FB_GNaK * 1.644;

    gbnaf = 0.0095;
    Ibnaf = gbnaf * (Vmf - Enaf);

    Iftot = Cmf * (Ikv + IK1f + INaKf + Ibnaf);

    IGap = GGAP * (V - Vmf);

    // Vmf = Vmf - dt * ((Iftot - IGap) / Cmf);
    residuals[32] -= state.get_var(1,32) - (-((Iftot - IGap) / Cmf));

    } // end FB if
    else{
        IGap = 0;
        residuals[28] -= state.get_var(1,28);
        residuals[29] -= state.get_var(1,29);
        residuals[32] -= state.get_var(1,32);
    }
    // #endif


    // compute ISAC;
    double strain = (SL - SLrest) / SLrest;
    if (UseSAC)
    {
        double PK = 1.0;
        double PNa = 1.0;
        double PCa = 1.0;
        double GSac = 0.00485 * 1.25; // need to be updated  ms/uF
        double Pm;
        const double Ke = 0.026590;// new fitting, was 0.02; // Dr Ismail's thesis
        const double strain_half = 0.163;// using 0.233; // Modeling of Stretch-Activated Sarcolemmal Channels in Smooth Muscle Cells
        // http://link.springer.com/chapter/10.1007%2F978-3-642-03882-2_197#page-1
        // 0.163 Role of Stretch-activated Channels in the Heart: Action Potential and Ca2+ Transients
        // in http://www.ncbi.nlm.nih.gov/books/NBK7490/
        const double Esac = -1.0; //mv see  Dr Ismail's thesis
        Pm = 1.0 / (1 + exp(-(strain - strain_half) / Ke));
        ISAC = Cm * (GSac * Pm * (V - Esac));
        int total_Sca = PNa + PCa + PK;
        Isac_Ca = PCa * ISAC / total_Sca;
        Isac_Na = PNa * ISAC / total_Sca;
        Isac_K = PK * ISAC / total_Sca;
    }
    // Concentrations
    naidot = (-3 * INaK - 3 * INaCa - IbNa - INa - Isac_Na) / (F * CRN_vi);
    kidot = (2 * INaK - IK1 - Ito - IKur - IKr - IKs - IbK - Isac_K) / (F * CRN_vi);
    // nai = nai + dt * naidot;
    residuals[7] -= state.get_var(1,7) - naidot;
    ki = ki + dt * kidot;
    residuals[9] -= state.get_var(1,9) - kidot;

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
    double RyRainfss = 0.505 - 0.427 / (1 + exp((( Cass + (fRyR * Cass) ) * 1000 - 0.29) / 0.082));
    double RyRoinfss = (1 - 1 / (1 +  exp(((Cass + (fRyR * Cass)   ) * 1000 - ((RyRass) + 0.22)) / 0.03)));
    double RyRcinfss = (1 / (1 + exp(((Cass + (fRyR * Cass )) * 1000 - ((RyRass) + 0.02)) / 0.01)));
    double Jrelss = nuss * ( (RyRoss) ) * (RyRcss) * RyRSRCass * ( CaSR2 -  (Cass) );

    double nu3 = 1 * Vnonjunct3;
    double RyRSRCa3 = (1 - 1 / (1 +  exp((CaSR1 - 0.3) / 0.1)));
    double RyRainf3 =  0.505 - 0.427 / (1 + exp(( (cai + ( fRyR * cai)  ) * 1000 - 0.29) / 0.082));
    double RyRoinf3 = (1 - 1 / (1 +  exp(( (cai + ( fRyR * cai) ) * 1000 - ((RyRa3) + 0.22)) / 0.03)));
    double RyRcinf3 = (1 / (1 +  exp(( (cai + (fRyR * cai ) ) * 1000 - ((RyRa3) + 0.02)) / 0.01)));
    double Jrel3 = nu3 * ( (RyRo3) ) * (RyRc3) * RyRSRCa3 * ( CaSR1 -  cai );

    Jrelss = fIRel * Jrelss;
    Jrel3 = fIRel * Jrel3;

    double JSRCaleak3 = GSR_leak * kSRleak * ( CaSR1 - cai ) * Vnonjunct3;
    double JSRCaleakss = GSR_leak * kSRleak * ( CaSR2 - (Cass) ) * Vss;

    double O_JSERCASR = J_SERCASR;
    double O_JSERCASRss = J_SERCASRss;
    double O_JBULKSERCA = J_bulkSERCA;
    double O_JBULKSERCAss = J_bulkSERCAss;
    double O_JRELss = Jrelss;
    double O_JREL = Jrel3;
    double O_JSRLEAK = JSRCaleak3;
    double O_JSRLEAKss = JSRCaleakss;

    double JCa, JCass;
    JCa = -BULK_CONST * J_bulkSERCA + JSRCaleak3 + Jrel3 + Jj_nj;
    JCass = -Jj_nj + JSRCaleakss - BULK_CONST * J_bulkSERCAss + Jrelss;

    double JSRCa1, JSRCa2;
    JSRCa1 = J_SERCASR - JSRCaleak3 - Jrel3;
    JSRCa2 = J_SERCASRss - JSRCaleakss - Jrelss;

    double dy;

    // dy = 0.5 * (-J_SERCASR + J_bulkSERCA) / Vnonjunct3;
    SERCACa += dt / 1000*dy;
    residuals[18] -= state.get_var(1,18) - 1000.0*dy;
    // dy = 0.5 * (-J_SERCASRss + J_bulkSERCAss) / Vss;
    SERCACass += dt / 1000*dy;
    residuals[19] -= state.get_var(1,19) - 1000.0*dy;

// Euler_inf(temp_type dt, temp_type gate, temp_type inf, temp_type tau) {
//     return (inf + (gate - inf) * exp(-(dt) / tau));
// }

    // RyRoss = RyRoinfss + (RyRoss - RyRoinfss) * exp(-(dt / 1000.0) / RyRtauactss);//Euler_inf(dt / 1000, RyRoss, RyRoinfss, RyRtauactss);
    residuals[20] -= state.get_var(1,20) - (RyRoinfss - RyRoss)/(1000.0*RyRtauactss);
    // RyRcss = RyRcinfss + (RyRcss - RyRcinfss) * exp(-(dt / 1000.0) / RyRtauinactss);//Euler_inf(dt / 1000, RyRcss, RyRcinfss, RyRtauinactss);
    residuals[21] -= state.get_var(1,21) - (RyRcinfss - RyRcss)/(1000.0*RyRtauinactss);
    // RyRass = RyRainfss + (RyRass - RyRainfss) * exp(-(dt / 1000.0) / RyRtauadapt);//Euler_inf(dt / 1000, RyRass, RyRainfss, RyRtauadapt);
    residuals[22] -= state.get_var(1,22) - (RyRainfss - RyRass)/(1000.0*RyRtauadapt);
    // RyRo3 = RyRoinf3 +(RyRo3 - RyRoinf3) * exp(-(dt / 1000.0) / RyRtauact);//Euler_inf(dt / 1000, RyRo3, RyRoinf3, RyRtauact);
    residuals[23] -= state.get_var(1,23) - (RyRoinf3 - RyRo3)/(1000.0*RyRtauact);
    // RyRc3 = RyRcinf3 +(RyRc3 - RyRcinf3) * exp(-(dt / 1000.0) / RyRtauinact);//Euler_inf(dt / 1000, RyRc3, RyRcinf3, RyRtauinact);
    residuals[24] -= state.get_var(1,24) - (RyRcinf3 - RyRc3)/(1000.0*RyRtauinact);
    // RyRa3 = RyRainf3 +(RyRa3 - RyRainf3) * exp(-(dt / 1000.0) / RyRtauadapt);//Euler_inf(dt / 1000, RyRa3, RyRainf3, RyRtauadapt);
    residuals[25] -= state.get_var(1,25) - (RyRainf3 - RyRa3)/(1000.0*RyRtauadapt);

    dy =  betass * ( JCass / Vss + ((-( RyR * ICaL) - IbCa - ICap - ICaT + 2 * INaCa - Isac_Ca)) / (2 * Vss * 1000 * F) );
    // Cass += dt / 1000*dy;
    residuals[15] -= state.get_var(1,15) - 1000.0*dy;

    dy = JCa / Vnonjunct3 * betai;
    // cai += dt / 1000*dy;
    residuals[8] -= state.get_var(1,8) - 1000.0*dy;

    dy =  betaSR1 * (DCaSR) * ( (CaSR2 - 2 * CaSR1 + CaSR1) / (dx * dx) + (CaSR1 - CaSR2) / (2 * 3 * (dx * dx)) ) + JSRCa1 / VSR3 * betaSR1;

    // CaSR1 += dt*dy;
    residuals[16] -= state.get_var(1,16) - dy;

    dy = betaSR2 * (DCaSR) * ( (CaSR2 - 2 * CaSR2 + CaSR1) / (dx * dx) + (CaSR2 - CaSR1) / (2 * 4 * (dx * dx)) ) + JSRCa2 / VSR4 * betaSR2;

    // CaSR2 += dt*dy;
    residuals[17] -= state.get_var(1,17) - dy;

    //Compute force model

    double P        = 1 - N - XBprer - XBpostr; //
    //   Compute single overlap fractions
    double sovr_ze   = std::min(len_thick / 2, SL / 2);       //   z-line end
    double sovr_cle = std::max(SL / 2 - (SL - len_thin), len_hbare / 2); //   centerline of end
    double len_sovr  = sovr_ze - sovr_cle;               //   single overlap length
    double SOVFThick = len_sovr * 2 / (len_thick - len_hbare); //   thick filament overlap frac
    double SOVFThin  = len_sovr / len_thin;              //   thin filament overlap frac

    //   Compute combined Ca binding to high- (w/XB) and low- (no XB) sites
    double Tropreg = (1 - SOVFThin) * TRPNCaL + SOVFThin * TRPNCaH;
    double permtot = sqrt(1 / (1 + pow((perm50 / Tropreg), nperm)));
    double inprmt = std::min(1.0 / permtot, 100.0);

    //   Adjustments for Ca activation, temperature, SL, stress and strain
    double konT    = kon * pow(Qkon, ((Temp - 310.0) / 10.0));
    double koffLT  = koffL * pow(Qkoff, ((Temp - 310.0) / 10.0)) * koffmod;
    double koffHT  = koffH * pow(Qkoff, ((Temp - 310.0) / 10.0)) * koffmod;
    double kn_pT   = kn_p * permtot * pow(Qkn_p, ((Temp - 310.0) / 10.0));
    double kp_nT   = kp_n * inprmt * pow(Qkp_n, ((Temp - 310.0) / 10.0));
    double fappT   = fapp * xbmodsp * pow(Qfapp, ((Temp - 310.0) / 10.0));
    double gapslmd = 1 + (1 - SOVFThick) * gslmod;
    double gappT   = gapp * gapslmd * xbmodsp * pow(Qgapp, ((Temp - 310.0) / 10.0));
    double hfmd    = exp(-sign(xXBprer) * hfmdc * ((xXBprer / x_0) * (xXBprer / x_0)));
    double hbmd    = exp(sign((xXBpostr - x_0)) * hbmdc * (((xXBpostr - x_0) / x_0) * ((xXBpostr - x_0) / x_0)));
    double hfT     = hf * hfmd * xbmodsp * pow(Qhf, ((Temp - 310.0) / 10.0));
    double hbT     = hb * hbmd * xbmodsp * pow(Qhb, ((Temp - 310.0) / 10.0));
    double gxbmd   = heav(x_0 - xXBpostr) * exp(sigmap * ((x_0 - xXBpostr) / x_0) * ((x_0 - xXBpostr) / x_0))
                     + (1 - heav(x_0 - xXBpostr)) * exp(sigman * (((xXBpostr - x_0) / x_0) * ((xXBpostr - x_0) / x_0)));
    double gxbT    = gxb * gxbmd * xbmodsp * pow(Qgxb, ((Temp - 310.0) / 10.0));

    //   Regulation and corssbridge cycling state derivatives
    double dTRPNCaL  = konT * cai * (1 - TRPNCaL) - koffLT * TRPNCaL;
    double dTRPNCaH  = konT * cai * (1 - TRPNCaH) - koffHT * TRPNCaH;
    // double dN_NoXB   = -kn_pT * N_NoXB + kp_nT * P_NoXB;
    // double dP_NoXB   = -kp_nT * P_NoXB + kn_pT * N_NoXB;
    double dN        = -kn_pT * N + kp_nT * P;
    // double //   dP      = -kp_nT*P + kn_pT*N - fappT*P + gappT*XBprer + gxbT*XBpostr;
    double dXBprer   = fappT * P - gappT * XBprer - hfT * XBprer + hbT * XBpostr;
    double dXBpostr  = hfT * XBprer - hbT * XBpostr - gxbT * XBpostr;
    // double dP        = -(dN + dXBprer + dXBpostr);

    //   steady-state fractions in XBprer and XBpostr using King-Altman rule
    double SSXBprer = (hb * fapp + gxb * fapp) / (gxb * hf + fapp * hf + gxb * gapp + hb * fapp + hb * gapp + gxb * fapp);
    double SSXBpostr = fapp * hf / (gxb * hf + fapp * hf + gxb * gapp + hb * fapp + hb * gapp + gxb * fapp);
    //   normalization for scaling active and passive force (maximal force)
    double Fnordv = kxb * x_0 * SSXBpostr;

    //   Calculate Forces (active, passive, preload, afterload)
    double force = kxb * SOVFThick * (xXBpostr * XBpostr + xXBprer * XBprer);
    double active = force / Fnordv;
    //Force = force;
    // NormalizedForce = active;
    double ppforce_c = heav(SL - SL_c) * PCon_c * (exp(PExp_c * fabs(SL - SL_c)) - 1);

    double ppforce_t = sign(SL - SLrest) * PCon_t *(exp(PExp_t *fabs((SL - SLrest))) - 1);
    // if (Rice_Usr_para.singlecell == False) ppforce_t += std::max(0, (PCon_c * exp(PExp_c * (SL - SL_c))));

    // Assume initial is generated by a preload force to counteract passive force
    // Preload force is computed here
    double PreloadF = sign((SLset - SLrest)) * PCon_t *(exp(PExp_t *fabs((SLset - SLrest))) - 1.0);
    // if (Rice_Usr_para.singlecell == False) PreloadF += std::max(0, (PCon_c * exp(PExp_c * (SLset - SL_c))));

    double ppforce = ppforce_t + ppforce_c;
    // ppforce += std::max(0, (PCon_c * exp(PExp_c * (SL - SL_c))));
    double preload = sign(SLset - SLrest) * PCon_t *(exp(PExp_t *fabs(SLset - SLrest)) - 1);
    // printf("%f\n", preload);
    preload = PreloadF;
    // preload += std::max(0, (PCon_c * exp(PExp_c * (SLset - SL_c))));
    double afterload = 0;  // either static constant or due to series elastic element
    /*if (Rice_Usr_para.SEon == True) {
        afterload = KSE * (SLset - SL);
    }
    */

    double dintf = /*(1 - Rice_Usr_para.SEon_LengthClamp) **/ (-ppforce + preload - active + afterload); //   total force
    double dSL = /*(1 - Rice_Usr_para.SEon_LengthClamp)**/((intf + (SLset - SL) * visc) / massf) * heav(SL - SLmin) * heav(SLmax - SL);
    //   Mean strain of strongly-bound states due to SL motion and XB cycling
    double dutyprer  = (hbT * fappT + gxbT * fappT)  //   duty fractions using the
                       / (fappT * hfT + gxbT * hfT + gxbT * gappT + hbT * fappT + hbT * gappT + gxbT * fappT);
    double dutypostr = fappT * hfT              //   King-Alman Rule
                       / (fappT * hfT + gxbT * hfT + gxbT * gappT + hbT * fappT + hbT * gappT + gxbT * fappT);
    double dxXBprer = dSL / 2.0 + xPsi / dutyprer * (-xXBprer * fappT + (xXBpostr - x_0 - xXBprer) * hbT);
    double dxXBpostr = dSL / 2.0 + xPsi / dutypostr * (x_0 + xXBprer - xXBpostr) * hfT;

    //   Ca buffering by low-affinity troponin C (LTRPNCa)
    double FrSBXB    = (XBpostr + XBprer) / (SSXBpostr + SSXBprer);
    double dFrSBXB   = (dXBpostr + dXBprer) / (SSXBpostr + SSXBprer);

    double dsovr_ze  = -dSL / 2.0 * heav(len_thick - SL);
    double dsovr_cle = -dSL / 2.0 * heav((2.0 * len_thin - SL) - len_hbare);
    double dlen_sovr = dsovr_ze - dsovr_cle;
    double dSOVFThin = dlen_sovr / len_thin;
    double dSOVFThick = 2.0 * dlen_sovr / (len_thick - len_hbare);

    double TropToT = Trop_conc * ((1 - SOVFThin) * TRPNCaL
                                  + SOVFThin * (FrSBXB * TRPNCaH + (1 - FrSBXB) * TRPNCaL));
    double dTropToT = Trop_conc * (-dSOVFThin * TRPNCaL + (1 - SOVFThin) * dTRPNCaL
                                   + dSOVFThin * (FrSBXB * TRPNCaH + (1 - FrSBXB) * TRPNCaL)
                                   + SOVFThin * (dFrSBXB * TRPNCaH + FrSBXB * dTRPNCaH - dFrSBXB * TRPNCaL
                                           + (1 - FrSBXB) * dTRPNCaL));
    /*double dforce = kxb * dSOVFThick * (xXBpostr * XBpostr + xXBprer * XBprer)
                    + kxb * SOVFThick * (dxXBpostr * XBpostr + xXBpostr * dXBpostr
                                         + dxXBprer * XBprer + xXBprer * dXBprer);*/
    
    // N        += dt*dN;
    residuals[40] -= state.get_var(1,40) - dN;
    // XBprer   += dt*dXBprer;
    residuals[41] -= state.get_var(1,41) - dXBprer;
    // XBpostr  += dt*dXBpostr;
    residuals[42] -= state.get_var(1,42) - dXBpostr;
    // SL       += dt*dSL;
    residuals[43] -= state.get_var(1,43) - dSL;
    // xXBpostr += dt*dxXBpostr;
    residuals[44] -= state.get_var(1,44) - dxXBpostr;
    // xXBprer  += dt*dxXBprer;
    residuals[45] -= state.get_var(1,45) - dxXBprer;
    // TRPNCaL  += dt*dTRPNCaL;
    residuals[46] -= state.get_var(1,46) - dTRPNCaL;
    // TRPNCaH  += dt*dTRPNCaH;
    residuals[47] -= state.get_var(1,47) - dTRPNCaH;
    // intf     += dt*dintf;
    residuals[48] -= state.get_var(1,48) - dintf;

    double dCai_feeback = dTropToT;


    //End solve force model


    //Handle some variables which have no equations in the original code
    residuals[30] -= state.get_var(1,30);
    residuals[31] -= state.get_var(1,31);
    residuals[35] -= state.get_var(1,35);
    residuals[36] -= state.get_var(1,36);
    residuals[38] -= state.get_var(1,38);
    residuals[39] -= state.get_var(1,39);


    // return state variables
    // new_state[0] = m;
    // new_state[1] = h;
    // new_state[2] = j;
    // new_state[3] = d;
    // new_state[4] = f;
    // new_state[5] = xr;
    // new_state[6] = xs;
    // new_state[7] = nai;
    // new_state[8] = cai;
    // new_state[9] = ki;
    // new_state[10] = fca;
    // new_state[11] = itr;
    // new_state[12] = its;
    // new_state[13] = isusr;
    // new_state[14] = isuss;
    // new_state[15] = Cass;
    // new_state[16] = CaSR1;
    // new_state[17] = CaSR2;
    // new_state[18] = SERCACa;
    // new_state[19] = SERCACass;
    // new_state[20] = RyRoss;
    // new_state[21] = RyRcss;
    // new_state[22] = RyRass;
    // new_state[23] = RyRo3;
    // new_state[24] = RyRc3;
    // new_state[25] = RyRa3;
    // new_state[26] = dd;
    // new_state[27] = ff;
    // new_state[28] = rkv;
    // new_state[29] = skv;
    // new_state[30] = kif;
    // new_state[31] = naif;
    // new_state[32] = Vmf;
    // new_state[33] = CNZ_a;
    // new_state[34] = CNZ_i;
    // new_state[35] = ACh_i;
    // new_state[36] = ACh_j;
    // new_state[37] = If_y;
    // // incorporation of INa state dependent drug effects
    // new_state[38] = BA;   // Blockade of INa channel,
    // new_state[39] = BA;    // blockade of INa channel.

    // new_state[40] = N;
    // new_state[41] = XBprer;
    // new_state[42] = XBpostr;
    // new_state[43] = SL;
    // new_state[44] = xXBpostr;
    // new_state[45] = xXBprer;
    // new_state[46] = TRPNCaL;
    // new_state[47] = TRPNCaH;
    // new_state[48] = intf;

    // keep log of currents; Fri 15 Apr 2016 11:51:29 BST
    // temp1 = INa / Cm;
    // temp2 = IKur / Cm;



    // m_IKur    = IKur/Cm;
    // m_ICaL    = ICaL/Cm;
    // m_INCX    = INaCa/Cm;
    // m_Jrel_ss = Jrelss;
    // m_ISAC    = ISAC/Cm;

    double Iion_tot = (INa + Ito + IKur + IKr + IKs + ICaL
               + IK1 + IbK + IbNa + IbCa + Iab + INaCa + INaK
               + ICap + If + ICaT + ISAC
    // #ifdef FIBROSIS
               + (FB_number * IGap)
    // #endif
              ) / Cm;

    state.set_membrane_current(Iion_tot);

    state.set_active_strain((SL - SLset) / SLset);

    // for(unsigned i=0; i<)

} // end Return Iion_tot



ImplicitCNZAtria::ImplicitCNZAtria(){
    //The largest timestep such that this model converges
    Intrinsic_dt = 0.02;
    //CNZ constants
    CRN_vcell = 20100.0; /* um3 */
    CRN_vi = 13668;
    CRN_vup = 1109.52/*0.0552*vcell*/;
    CRN_vrel = 96.48/*0.0048*vcell*/;
    T = 310; /* 37 Celcius */
    CRN_Tfac = 3;
    CRN_Csp = 1e+6; /* pF/cm2 */
    F = 96.4867; /* coul/mmol */
    R = 8.3143; /* J K-1 mol-1 */
    CRN_kb = 5.4; /* mM */
    CRN_nab = 140; /* mM */
    CRN_cab = 1.8; /* mM */
    CRN_nac = 140;
    CRN_cac = 1.8;
    CRN_kc = 5.4;
    CRN_gna = 7.8; /* nS/pF */
    CRN_gto = 0.1652; /* nS/pF */
    CRN_gkr = 0.029411765; /* nS/pF */
    CRN_gks = 0.12941176; /* nS/pF */
    CRN_gcaL = 0.1294; /* nS/pF */
    CRN_ErL = 65.0; /* mV */
    CRN_gk1 = 0.09; /* nS/pF */
    CRN_gbna = 0.0006744375; /* nS/pF */
    CRN_gbk = 0.0;
    CRN_gbca = 0.001131; /* nS/pF */
    CRN_inakbar = 0.59933874; /* pA/pF */
    CRN_kmnai = 10.0; /* mM */
    CRN_kmko = 1.5; /* mM */
    CRN_icapbar = 0.275; /* pA/pF */
    CRN_kmcap = 0.0005; /* mM */
    CRN_knacalr = 1600.0; /* pA/pF */
    CRN_kmnalr = 87.5; /* mM */
    CRN_kmcalr = 1.38; /* mM */
    CRN_ksatlr = 0.1;
    CRN_gammalr = 0.35;
    CRN_trpnbar = 0.070; /* mM */
    CRN_cmdnbar = 0.050; /* mM */
    CRN_csqnbar = 10; /* mM */
    CRN_kmcsqn = 0.8; /* mM */
    CRN_kmcmdn = 0.00238; /* mM */
    CRN_kmtrpn = 0.0005; /* mM */
    CRN_grelbar = 30.0; /* ms-1 */
    CRN_kmup = 0.00092; /* mM */
    CRN_iupbar = 0.005; /* mM/ms */
    CRN_caupmax = 15.0; /* mM */
    CRN_kupleak = 0.00033333336/*iupbar/caupmax*/; /* ms-1 */
    CRN_tautr = 180.0; /* ms */
    CRN_gf_na = 0.0944;////0.10/* 0.09*/; /* nS/pF */
    CRN_gf_k = /*0.0752*/ 0.0752;
    CRN_gf = 0.025;//0.0752;
    CRN_Ef = -22.0; /* mV */
    gcaT = (0.15 * 0.22) / 17.62;
    EcaT = 45.0;
    MAL_gto = 0.75471 * 0.1962; //nS/pF **** in manuscript it is given in nS, so has been /100 here
    MAL_gkur = 0.05874;
    CNZ_gkur = 0.006398;

    // New parameters

     GONG_gto = 0.103;
     gf = 0.0385;

    // KM model intracellular model parameters
    Vss = 2 * 4.99232e-5;
    rjunct = 6.5;
    lcell = 122.051;

    dx = 1.625;
    Aj_nj = 2492.324412;// = M_PI*rjunct*2*lcell*0.5; // atea between junct and non junct
    xj_nj = 0.822500;// = 0.02/2 + dx/2; // diffusion distance from center to center of junct to first njunct
    xj_nj_Nai = 3.260000;// = 0.02/2 + 2*dx; // diffusion distance from center of junct to center of njunct (between 2nd and 3rd njunct)

    Vnonjunct3 = 6 * 0.002531;
    VSR3 = 2 * 0.000057;
    VSR4 = 2 * 0.000080;
    Vcytosol = 0.008150;
    Vnonjunct_Nai = 0.008100;

    BCa = 24e-3;
    SLlow = 165;
    SLhigh = 13;

    KdBCa = 2.38e-3;
    KdSLlow = 1.1;
    KdSLhigh = 13e-3;

    CSQN =  6.7;
    KdCSQN = 0.8;

    BNa = 1.1319;
    KdBNa = 10;

    DCa = 780;// % µm2/s
    DCaSR = 44;//
    DCaBm = 25; //% µm2/s
    DNa = 0.12;

    SERCAKmf = 0.25e-3;
    SERCAKmr = 1.8;
    k4 = 7.5 ;// % pump rate
    k3 = 2.314815; //= k4 / SERCAKmr^2;
    k1 = 7500000.000000;// = 1000^2 * k4;
    k2 = 0.468750;// = k1 * SERCAKmf^2;
    cpumps = 40e-3;
    kSRleak = 6e-3;

    // Fibroblast stuff
    Cmf = 6.3; // pF
    vif = 0.00000137;//0.00000137;//0.00000000137; // um3  (from 0.00137 nL -> 0.00000137 nm3 -> 0.00000000137 um3
    naof = 130.0110; //mM
    kof = 5.3581; //mM
    Cm = 100.0;




    //Rice myofilament constants
    #define p0 1.754251548964904
    #define p1 0.905622641626625
    #define p2 0.499437793063966
    #define p3 0.400000001127317
    #define p4 1.000000000000000
    #define p5 1.000000000000000
    #define p6 1.981229252026256
    #define p7 0.511387864324941
    #define p8 0.023420000000000
    // rolling back to its original value.
    SLmax   = 2.4;// belus 2010. fig6, was 2.4;        //   (um) maximum sarcomere length
    SLmin   = 1.4;//1.4;        //   (um) minimum sarcomere length
    len_thin  = 1.2;//1.2;      //   (um) thin filament length
    len_thick = 1.65; // //1.65;     //   (um) thick filament length
    len_hbare = 0.1;      //   (um) length of bare portion of thick filament
    //   Temperature Dependence
    Qkon = 1.5;
    Qkoff = 1.3;
    // static const double Qkoff = 1.4;
    Qkn_p = 1.6;
    Qkp_n = 1.6;
    Qfapp = 6.25;
    // static const double Qgapp = 6.25;
    Qgapp = 2.5;
    Qhf = 6.25;
    Qhb = 6.25;
    Qgxb = 6.25;
    //   Ca binding to troponin
    kon     = p0 * 50/*e-3*/;    //   (1/[ms uM])
    koffL   = p1 * p0 * 250e-3; //   (1/ms)
    koffH   = p1 * p0 * 25e-3;  //   (1/ms)
    perm50  = p6 * 0.5;      //   perm variable that controls n to p transition
    nperm   = p7 * 15;       //     in Hill-like fashion
    kn_p    = 500e-3;     //   (1/ms)
    kp_n    = 50e-3;      //   (1/ms)
    koffmod = 1.0;        //   mod to change species
    //   Thin filament regulation and crossbridge cycling
    fapp    = p2 * 500e-3;   //   (1/ms) XB on rate
    gapp    = p3 * 70e-3;    //   (1/ms) XB off rate
    gslmod  = 6;          //   controls SL effect on gapp
    hf      = p2 * 2000e-3;  //   (1/ms) rate between pre-force and force states
    hfmdc   = 5;          //
    hb      = p3 * 400e-3;   //   (1/ms) rate between pre-force and force states
    hbmdc   = 0;          //
    gxb     = p3 * 70e-3;    //   (1/ms) ATP consuming transition rate
    sigmap  = 8;          //   distortion dependence of STP using transition gxb
    sigman  = 1;          //
    xbmodsp = 4.0 / 3.0;      //   mouse specific modification for XB cycling rates
    //   Mean strain of strongly-bound states
    x_0     = 0.007;      //   (um) strain induced by head rotation
    xPsi    = 2;          //   scaling factor balancing SL motion and XB cycling
    //   Normalized active a nd passive force
    PCon_t  = 0.002;      //   (norm Force) passive force due to titin
    PExp_t  = 10;         //     these apply to trabeculae and single cells only
    SL_c    = 2.25;       //   (um) resting length for collagen
    PCon_c  = 0.02;       //   (norm Force) passive force due to collagen
    PExp_c  = 70;         //     these apply to trabeculae and single cells only
    //   Calculation of complete muscle response
    massf   = 0.000025e6;  //   ([norm Force ms^2]/um) muscle mass
    visc    = 0.003e3;    //   ([norm Force ms]/um) muscle viscosity
    KSE     = 1;          //   (norm Force/um) series elastic element
    kxb     = 120;        //   (mN/mm^2) maximal force
    Trop_conc = 70e-3;       //   (uM) troponin concentration
    Temp = 310.0/*-17.0*/;

    SLset = 2.2;         // initial length
    SLrest = 1.9;       //   (um) rest SL length for 0 passive force
}

}//end namespace