#include "CNZ_Atria.h"

namespace oomph{

    CNZCell::CNZCell()   :   CellModelBase()
    {
        //The required storage
        //..without variables with zero residual
        this->Required_Storage = 45;
        
        //Overload the requests which are not the default value
        Requires_Strain = true;
        Requires_Fibrosis = true;

        //Constants
        CRN_vcell = 20100.0; /* um3 */
        CRN_vi = 13668;
        CRN_vup = 1109.52/*0.0552*vcell*/;
        CRN_vrel = 96.48/*0.0048*vcell*/;
        T = 310; /* 37 Celcius */
        CRN_Tfac = 3;
        CRN_Csp = 1e+6; /* pF/cm2 */
        F = 96.4867; /* coul/mmol */
        R = 8.3143; /* J K-1 mol-1 */
        FoRT = 0.0374351923;    /*  F/RT */
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
        GGAP = 0.0;
        FB_Ikv_shift = 0.0;
        FB_Gkv = 0.0;
        //Membrane capacitance
        Cm = 100.0; //pF
        //ISAC parameters
        //densities of the different ISAC types
        ISAC_pNa = 1.0;
        ISAC_pK = 1.0;
        ISAC_pCa = 1.0;
        ISAC_total_Sca = (ISAC_pNa + ISAC_pK + ISAC_pCa);
        //parameters
        ISAC_Esac = -1.0;
        ISAC_Ke = 0.026590;
        ISAC_strain_half = 0.163;
        ISAC_GSac = 0.0;//6.0625e-3;
        //Switch parameters from SingleCellParameters.cpp
        IKur_type_CNZ = 1;


        //BEGIN Assign force model parameters
        double p0 = 1.754251548964904;
        double p1 = 0.905622641626625;
        double p2 = 0.499437793063966;
        double p3 = 0.400000001127317;
        double p4 = 1.000000000000000;
        double p5 = 1.000000000000000;
        double p6 = 1.981229252026256;
        double p7 = 0.511387864324941;
        double p8 = 0.023420000000000;
        // rolling back to its original value.
        SLmax   = 2.4;// belus 2010. fig6, was 2.4;        //   (um) maximum sarcomere length
        SLmin   = 1.4;//1.4;        //   (um) minimum sarcomere length
        len_thin  = 1.2;//1.2;      //   (um) thin filament length
        len_thick = 1.65; // //1.65;     //   (um) thick filament length
        len_hbare = 0.1;      //   (um) length of bare portion of thick filament
        //   Temperature Dependence
        Qkon = 1.5;
        Qkoff = 1.3;
        // Qkoff = 1.4;
        Qkn_p = 1.6;
        Qkp_n = 1.6;
        Qfapp = 6.25;
        // Qgapp = 6.25;
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
        //END Assign force model parameters
    }






    //====================================================================
    //====================================================================
    // Define the channel currents
    //====================================================================
    //====================================================================
    void CNZCell::INa_current(CellState &state){
        state.set_ina_current(Cm*
                            get_GNa(state.cell_type())* //!!!!!
                            CRN_gna*                                //Blockade, originally from state[39] and state[40]
                            pow(state.var(0,0),3)*
                            state.var(0,1)*
                            state.var(0,2)*
                            (state.vm() - state.ena()));
    }

    void CNZCell::IKr_current(CellState &state){
        state.set_ikr_current(Cm*
                            get_GKr(state.cell_type())*
                            CRN_gkr*
                            state.var(0,5)*
                            (state.vm() - state.ek())/
                            (1.0 + exp((state.vm() + 15.0) / 22.4)));
    }

    void CNZCell::IKs_current(CellState &state){
        state.set_iks_current(Cm*
                            get_GKs(state.cell_type())*
                            CRN_gks*
                            pow(
                                state.var(0,6)
                                ,2.0
                                )*
                            (state.vm() - state.ek()));
    }

    void CNZCell::ICaL_current(CellState &state){
        state.set_ical_current(Cm*
                            2.125*
                            get_GCaL(state.cell_type())*
                            CRN_gcaL*
                            state.var(0,3)*
                            state.var(0,4)*
                            state.var(0,10)*
                            (state.vm() - state.eca()));
    }

    void CNZCell::IK1_current(CellState &state){
        state.set_ik1_current(Cm*
                            get_GK1(state.cell_type())*
                            CRN_gk1*
                            (state.vm() - state.ek() + get_IK1_v_shift(state.cell_type()))/
                            (1.0 + exp(0.07 * (state.vm() + 80.0 + get_IK1_v_shift(state.cell_type())))));
    }

    void CNZCell::Iab_current(CellState &state){
        state.set_iab_current(Cm*
                            0.0003879*
                            (state.vm() + 69.6)/
                            (1.0 - 0.8377*exp((state.vm() + 49.06)/1056.0)));
    }

    void CNZCell::IbK_current(CellState &state){
        state.set_ibk_current(Cm*
                            CRN_gbk*
                            (state.vm() - state.ek()));
    }

    void CNZCell::IbCa_current(CellState &state){
        state.set_ibca_current(Cm*
                            get_Gbca(state.cell_type())*
                            CRN_gbca*
                            (state.vm() - state.eca()));
    }

    void CNZCell::IbNa_current(CellState &state){
        state.set_ibna_current(Cm*
                            CRN_gbna*
                            (state.vm() - state.ena()));
    }

    void CNZCell::ICap_current(CellState &state){
        state.set_icap_current(Cm*
                            1.4*
                            get_GCap(state.cell_type())*
                            CRN_icapbar*
                            state.var(0,15)/
                            (state.var(0,15) + CRN_kmcap ));
    }

    void CNZCell::INaCa_current(CellState &state){
        state.set_inaca_current(Cm*
                            get_GNaCa(state.cell_type())*
                            CRN_knacalr/
                            (pow(CRN_kmnalr,3.0)+pow(state.na_o(),3.0))/
                            (CRN_kmcalr+state.ca_o())/
                            (1.0+CRN_ksatlr*exp((CRN_gammalr-1.0)*state.vm()*FoRT))*
                            (pow(state.var(0,7),3.0)*state.ca_o()*exp(state.vm()*CRN_gammalr*FoRT)-
                                pow(state.na_o(),3.0)*state.var(0,15)*exp(state.vm()*(CRN_gammalr-1.0)*FoRT)));
    }

    void CNZCell::INaK_current(CellState &state){
        state.set_inak_current(Cm*
                            CRN_inakbar/
                            ((1.0+0.1245*exp(-0.1*state.vm()*FoRT)+0.0365*((exp(state.na_o()/67.3)-1.0)/7.0)*exp(-state.vm()*FoRT))*
                            (1.0 + CRN_kmko / state.k_o())*
                            (1.0 + pow(CRN_kmnai/state.var(0,7),1.5))));
    }

    void CNZCell::Ito_current(CellState &state){
        state.set_ito_current(Cm*
                            get_Gto(state.cell_type())*
                            MAL_gto*
                            state.var(0,11)*
                            state.var(0,12)*
                            (state.vm() - state.ek()));
    }

    void CNZCell::IKur_current(CellState &state){
        state.set_ikur_current(Cm*
                            get_IKur_cond(state.cell_type(),this->mutation())*
                            get_GKur(state.cell_type(),this->mutation())*
                            CNZ_gkur*
                            (get_IKur_c(state.cell_type(),this->mutation()) + get_IKur_x0(state.cell_type(),this->mutation())/(1.0 + exp((state.vm() - get_Ikur_y0(state.cell_type(),this->mutation()))/(-8.26597))))*
                            state.var(0,28)*
                            state.var(0,29)*
                            (state.vm() - state.ek()));
    }

    void CNZCell::If_current(CellState &state){
                /*
                Cm*
                get_Gf(get_state.cell_type()_at_node_CNZ(n))*
                CRN_gf*
                state.var(0,30)*
                (state.vm() - Rev_Pot);
                */
    }

    void CNZCell::IK1f_current(CellState &state){
        double ak1f = 0.1 / (1 + exp(0.06 * (state.var(0,44) - state.ekf() - 200.0 ) )  );
        double bk1f = (3 * exp(0.0002 * (state.var(0,44) - state.ekf() + 100.0)) + exp(0.1 * (state.var(0,44) - state.ekf() - 10.0))) / (1.0 + exp(-0.5 * (state.var(0,44) - state.ekf())));
        double gk1f = 0.4822;
        // state.set_ik1f_current(gk1f * (ak1f / (ak1f + bk1f)) * (state.var(0,44) - state.ekf()));
        state.set_ik1f_current(0.0);

    }

    void CNZCell::IbNaf_current(CellState &state){
        double gbnaf = 0.0095;
        state.set_ibnaf_current(gbnaf * (state.var(0,44) - state.enaf()));
    }

    void CNZCell::ICaT_current(CellState &state){
        // return  0.0;
                /*
                Cm*
                get_GCaT(state.cell_type())*
                gcaT*
                state.var(0,27)*
                state.var(0,26)*
                (state.vm() - EcaT);
                */
    }

    void CNZCell::IGap_current(CellState &state){
        state.set_igap_current(GGAP * (state.vm() - state.var(0,44)));
    }

    inline void CNZCell::IKv_current(CellState &state){
        double gkv = FB_Gkv * 0.25;
        state.set_ikv_current(Cmf * gkv * state.var(0,40) * state.var(0,41) * (state.var(0,44) - state.ekf()));
    }

    void CNZCell::ISAC_current(CellState &state){
        double ISAC = Cm*
                        ISAC_GSac*
                        (state.vm() - ISAC_Esac)/
                        (1.0 + exp(-(state.stress() - ISAC_strain_half)/ISAC_Ke));

        state.set_isac_na_current(ISAC_pNa*ISAC/ISAC_total_Sca);
        state.set_isac_ca_current(ISAC_pCa*ISAC/ISAC_total_Sca);
        state.set_isac_k_current(ISAC_pK*ISAC/ISAC_total_Sca);
    }

    //fibrosis currents

    //====================================================================
    //====================================================================
    // Define the channel currents
    //====================================================================
    //====================================================================

    void CNZCell::ENa_reversal(CellState &state)
    {
        state.set_ena(26.71*log(state.na_o() / state.var(0,7)));
    }
    void CNZCell::EK_reversal(CellState &state)
    {
        state.set_ek(26.71*log(state.k_o() / state.var(0,9)));
    }
    void CNZCell::ECa_reversal(CellState &state)
    {
        state.set_eca(13.35*log(state.ca_o() / state.var(0,8)));
    }

    void CNZCell::EKf_reversal(CellState &state)
    {
        state.set_ekf(26.71 * log(kof / state.var(0,42)));
    }
    void CNZCell::Enaf_reversal(CellState &state)
    {
        state.set_enaf(26.71 * log(naof / state.var(0,43)));
    }






    //====================================================================
    //====================================================================
    // Define the channel residuals
    //====================================================================
    //====================================================================
    inline void CNZCell::INa_current_residual(const CellState &state, Vector<double> &residuals){
        double alpha, beta;
        //m
        //Calculate rates
        if(std::fabs(state.vm() + 47.13) > 1e-10){alpha = 0.32 * (state.vm() + 47.13) / (1.0 - exp(-0.1 * (state.vm() + 47.13)));}
        else{alpha = 3.2;}
        beta = 0.08 * exp(-state.vm() / 11.0);

        residuals[0] -= state.var(1,0) + state.var(0,0)*(alpha + beta) - alpha;

        //h
        //Calculate rates
        if(state.vm()>=-40.0){
            alpha = 0.0;
            beta = 1.0 / (0.13 * (1.0 + exp((state.vm() + 10.66) / -11.1)));
        }
        else{
            alpha = 0.135 * exp((state.vm() + 80.0) / -6.8);
            beta = 3.56 * exp(0.079 * state.vm()) + 3.1e5 * exp(0.35 * state.vm());
        }

        residuals[1] -= state.var(1,1) + state.var(0,1)*(alpha + beta) - alpha;

        //j
        //Calculate rates
        if (state.vm() >= -40.0)
        {
            alpha  = 0.0;
            beta = 0.3 * exp(-2.535e-7 * state.vm()) / (1.0 + exp(-0.1 * (state.vm() + 32.0)));
        }
        else
        {
            alpha = (-1.2714e5 * exp(0.2444 * state.vm()) - 3.474e-5 * exp(-0.04391 * state.vm())) * (state.vm() + 37.78) / (1.0 + exp(0.311 * (state.vm() + 79.23)));
            beta = 0.1212 * exp(-0.01052 * state.vm()) / (1.0 + exp(-0.1378 * (state.vm() + 40.14)));
        }
        residuals[2] -= state.var(1,2) + state.var(0,2)*(alpha + beta) - alpha;
    }

    inline void CNZCell::IKr_current_residual(const CellState &state, Vector<double> &residuals){
        double alpha, beta;
        //xr
        //Calculate rates
        if (std::fabs(state.vm() + 14.1) > 1e-10){alpha = 0.0003 * (state.vm() + 14.1) / (1.0 - exp(-(state.vm() + 14.1) / 5.0));}
        else{alpha = 0.0015;}
        if (fabs(state.vm() - 3.3328) > 1e-10){beta = 0.000073898 * (state.vm() - 3.3328) / (exp((state.vm() - 3.3328) / 5.1237) - 1.0);}
        else{beta = 3.7836118e-4;}
        residuals[5] += state.var(1,5) - ( 1.0 / ( 1.0 + exp( -(state.vm() + get_IKr_ac_shift(state.cell_type(), this->mutation()) + 14.1) / (6.5 * get_IKr_ac_grad(state.cell_type(), this->mutation()))) ) - state.var(0,5) )*(alpha + beta);
    }

    inline void CNZCell::IKs_current_residual(const CellState &state, Vector<double> &residuals){
        double alpha, beta;
        //xs
        //Calculate rates
        if (std::fabs(state.vm() - 19.9) > 1e-10){
            
            alpha = 0.00004 * (state.vm() - 19.9) / (1.0 - exp(-(state.vm() - 19.9) / 17.0));
            beta = 0.000035 * (state.vm() - 19.9) / (exp((state.vm() - 19.9) / 9.0) - 1.0);
        }
        else{
            alpha = 0.00068;
            beta = 0.000315;
        }
        residuals[6] += state.var(1,6) - ( sqrt(1.0 / (1.0 + exp(-(state.vm() - 19.9 - get_IKs_shift(state.cell_type(), this->mutation())) / (12.7 * get_IKs_grad(state.cell_type(), this->mutation()))))) - state.var(0,6) )*2.0*(alpha + beta);
    }

    inline void CNZCell::ICaL_current_residual(const CellState &state, Vector<double> &residuals){
        //d
 // and jacobian without rates         
        if (std::fabs(state.vm() + 10.0) > 1e-10)
        {
            residuals[3] += state.var(1,3) - 0.035*(state.vm() + 10.0)*( (1.0+exp(-(state.vm()+10.0)/6.24)) / (1.0-exp(-(state.vm()+10.0)/6.24)) )*(1.0/(1.0+exp(-(state.vm()+10.0)/8.0)) - state.var(0,3));
        }
        else
        {   
            residuals[3] += state.var(1,3) - (1.0/4.579) * (1.0 + exp((state.vm() + 10.0) / -6.24)) * ( 1.0 / (1.0 + exp((state.vm() + 10.0) / -8.0)) - state.var(0,3));
        }
        //f
        residuals[4] += state.var(1,4) - ( exp(-(state.vm()+28.0)/6.9)/(1.0+exp(-(state.vm()+28.0)/6.9)) - state.var(0,4) ) / (9.0 / (0.0197 * exp(-pow(0.0337*(state.vm() + 10),2.0)) + 0.02));
        //fca
        residuals[10] += state.var(1,10) - ( 1.0 / (1.0 + (state.var(0,15) / 0.00035)) - state.var(0,10) )/(2.0);
    }

    //====================================================================
    //Either no variables for these channels or variables shared with other channels
    inline void CNZCell::IK1_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::Iab_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::IbK_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::IbCa_current_residual(const CellState &state, Vector<double> &residuals){

    }
    inline void CNZCell::IbNa_current_residual(const CellState &state, Vector<double> &residuals){

    }
    inline void CNZCell::ICap_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::INaCa_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::INaK_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::IGap_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    inline void CNZCell::ISAC_current_residual(const CellState &state, Vector<double> &residuals) {

    }
    //End no variables for these channels or variables shared with other channels
    //====================================================================

    inline void CNZCell::Ito_current_residual(const CellState &state, Vector<double> &residuals){
        //itr //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[11] += state.var(1,11) - ( 1.0 / (1.0 + exp(-(state.vm() - 1.0)/11.0)) - state.var(0,11) ) / ((0.0035 * exp(-(state.vm() / 15.0)) + 0.0015)*1000.0);
        //its //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[12] += state.var(1,12) - ( 1.0 / (1.0 + exp((state.vm() + 40.5) / 11.5)) - state.var(0,12) ) / ((0.025635 * exp (-(state.vm() + 52.45) / 7.94135) + 0.01414)*1000.0);
    }

    inline void CNZCell::IKur_current_residual(const CellState &state, Vector<double> &residuals){
        //CNZ_a
        residuals[28] += state.var(1,28) - ( 1.0 / (1.0 + exp(-(state.vm() - (-6.0 + get_IKur_Vhchange(state.cell_type(), this->mutation()) )) / (8.6 * get_IKur_slope(state.cell_type(),this->mutation()) ))) - state.var(0,28) ) / ( ( get_IKur_timeconstants(state.cell_type(),this->mutation())*(45.6666746826 / (1.0 + exp((state.vm() + 11.2306497073) / 11.5254705962)) + 4.26753514993)*(0.262186042981 / (1.0 + exp((state.vm() + 35.8658312707) / (-3.87510627762))) + 0.291755017928) ) / 3.5308257834747638 );
        //CNZ_i
        residuals[29] += state.var(1,29) - ( (get_IKur_inac_mult(state.cell_type(),this->mutation()) * 0.52424) / (1.0 + exp((state.vm() + 15.1142 + get_IKur_inac_shift(state.cell_type(),this->mutation()) ) / (7.567021 * get_IKur_inac_grad(state.cell_type(),this->mutation()) ))) + 0.4580778 + get_IKur_inac_add(state.cell_type(),this->mutation())  - state.var(0,29) ) / ( (2328.0 / (1.0 + exp((state.vm() - 9.435) / (3.5827))) + 1739.139) / 3.5308257834747638 );
    }

    inline void CNZCell::If_current_residual(const CellState &state, Vector<double> &residuals){
        //If_y
        residuals[30] += state.var(1,30) - ( 1.0 / (1.0 + exp((state.vm() + 90.95 + get_If_vshift(state.cell_type(),this->mutation())) / (10.1 * get_If_grad(state.cell_type(),this->mutation()))) ) - state.var(0,30) ) * (1.2783E-04 * exp(-state.vm() / 9.2424) + 121.6092 * exp(state.vm() / 9.2424));
    }

    inline void CNZCell::ICaT_current_residual(const CellState &state, Vector<double> &residuals){
        //dd
        residuals[26] += state.var(1,26) - ( 1.0 / (1.0 + exp(-(state.vm() + 37.0) / 6.8)) - state.var(0,26) ) * ( 1.0680*exp((state.vm() + 26.3)/30.0) + 1.0680*exp(-(state.vm() + 26.3)/30.0) );
        //ff
        residuals[27] += state.var(1,27) - ( 1.0 / (1.0 + exp((state.vm() + 71.0) / 9.0)) - state.var(0,27) ) * ( 0.0153*exp(-(state.vm() + 71.0)/83.3) + 0.0150*exp((state.vm() + 71.0)/15.38) );
        
    }

    //====================================================================
    //====================================================================
    // Define the concentration residuals
    //====================================================================
    //====================================================================

    //The current residual functions
    inline void CNZCell::Ca_i_residual(const CellState &state, Vector<double> &residuals){
    }

    inline void CNZCell::Na_i_residual(const CellState &state, Vector<double> &residuals){
        //nai
        residuals[7] += state.var(1,7) - (-3.0 * state.inak() - 3.0 * state.inaca() - state.ibna() - state.ina() - state.isac_na()) / (F *CRN_vi);
    }

    inline void CNZCell::K_i_residual(const CellState &state, Vector<double> &residuals){
        //ki       
        residuals[9] += state.var(1,9) - (2.0 * state.inak() - state.ik1() - state.ito() - state.ikur()- state.ikr() - state.iks() - state.ibk() - state.isac_k()) / (F * CRN_vi);
    }

    inline void CNZCell::Calcium_dynamics_residual(const CellState &state, Vector<double> &residuals){
        //values required for changes in various calcium concentrations
        double betass;
        betass = pow(
                     ( 1 + SLlow * KdSLlow / pow((state.var(0,15) + KdSLlow), 2)
                       + SLhigh * KdSLhigh / pow((state.var(0,15) + KdSLhigh), 2)
                       + BCa * KdBCa / pow(((state.var(0,15)) + KdBCa), 2)  )
                     , (-1));
        double betai, gammai;
        betai = pow(( 1 + BCa * KdBCa / pow((state.var(0,8) + KdBCa), 2)  ), (-1));
        gammai = BCa * KdBCa / pow((state.var(0,8) + KdBCa), 2);
        double betaSR1, betaSR2;
        betaSR1 = pow( ( 1 + CSQN * KdCSQN / pow((state.var(0,16) + KdCSQN), 2) ), (-1));
        betaSR2 = pow( ( 1 + CSQN * KdCSQN / pow((state.var(0,17) + KdCSQN), 2) ), (-1));
        double Jj_nj;
        Jj_nj = DCa * Aj_nj / xj_nj * (state.var(0,15) - state.var(0,8)) * 1e-6;
        double J_SERCASR, J_bulkSERCA;
        double J_SERCASRss, J_bulkSERCAss;
        J_SERCASR =  (-k3 * pow(state.var(0,16), 2) * (cpumps - state.var(0,18)) + k4 * state.var(0,18)) * Vnonjunct3 * 2;
        J_bulkSERCA = (k1 * pow(state.var(0,8), 2) * (cpumps - state.var(0,18)) - k2 * state.var(0,18)) * Vnonjunct3 * 2;
        J_SERCASRss = (-k3 * pow(state.var(0,17), 2) * (cpumps - state.var(0,19)) + k4 * state.var(0,19)) * Vss * 2;
        J_bulkSERCAss = (k1 * pow((state.var(0,15)), 2) * (cpumps - state.var(0,19)) - k2 * state.var(0,19)) * Vss * 2;
        double RyRtauadapt = 1.0;
        double RyRtauactss = 5e-3;
        double RyRtauinactss = 15e-3;
        double RyRtauact = 18.75e-3;
        double RyRtauinact = 87.5e-3;
        double nuss = 625.0 * Vss;
        double RyRSRCass = (1 - 1 / (1 +  exp((state.var(0,17) - 0.3) / 0.1)));
        double RyRainfss = 0.505 - 0.427 / (1 + exp(( ( state.var(0,15) + (get_fRyR(state.cell_type(), this->mutation()) * state.var(0,15)) ) * 1000 - 0.29) / 0.082));
        double RyRoinfss = (1 - 1 / (1 +  exp(((state.var(0,15) + (get_fRyR(state.cell_type(), this->mutation()) * state.var(0,15))   ) * 1000 - ((state.var(0,22)) + 0.22)) / 0.03)));
        double RyRcinfss = (1 / (1 + exp(((state.var(0,15) + (get_fRyR(state.cell_type(), this->mutation()) * state.var(0,15) )) * 1000 - ((state.var(0,22)) + 0.02)) / 0.01)));
        double Jrelss = nuss * ( (state.var(0,20)) ) * (state.var(0,21)) * RyRSRCass * ( state.var(0,17) -  (state.var(0,15)) );
        double nu3 = 1 * Vnonjunct3;
        double RyRSRCa3 = (1 - 1 / (1 +  exp((state.var(0,16) - 0.3) / 0.1)));
        double RyRainf3 =  0.505 - 0.427 / (1 + exp(( (state.var(0,8) + ( get_fRyR(state.cell_type(), this->mutation()) * state.var(0,8))  ) * 1000 - 0.29) / 0.082));
        double RyRoinf3 = (1 - 1 / (1 +  exp(( (state.var(0,8) + ( get_fRyR(state.cell_type(), this->mutation()) * state.var(0,8)) ) * 1000 - ((state.var(0,25)) + 0.22)) / 0.03)));
        double RyRcinf3 = (1 / (1 +  exp(( (state.var(0,8) + (get_fRyR(state.cell_type(), this->mutation()) * state.var(0,8) ) ) * 1000 - ((state.var(0,25)) + 0.02)) / 0.01)));
        double Jrel3 = nu3 * ( (state.var(0,23)) ) * (state.var(0,24)) * RyRSRCa3 * ( state.var(0,16) -  state.var(0,8) );
        Jrelss = get_fIRel(state.cell_type(), this->mutation()) * Jrelss;
        Jrel3 = get_fIRel(state.cell_type(), this->mutation()) * Jrel3;
        double JSRCaleak3 = get_GSR_leak(state.cell_type(), this->mutation()) * kSRleak * ( state.var(0,16) - state.var(0,8) ) * Vnonjunct3;
        double JSRCaleakss = get_GSR_leak(state.cell_type(), this->mutation()) * kSRleak * ( state.var(0,17) - (state.var(0,15)) ) * Vss;
        double JCa, JCass;
        JCa = -get_BULK_CONST(state.cell_type(), this->mutation()) * J_bulkSERCA + JSRCaleak3 + Jrel3 + Jj_nj;
        JCass = -Jj_nj + JSRCaleakss - get_BULK_CONST(state.cell_type(), this->mutation()) * J_bulkSERCAss + Jrelss;
        double JSRCa1, JSRCa2;
        JSRCa1 = J_SERCASR - JSRCaleak3 - Jrel3;
        JSRCa2 = J_SERCASRss - JSRCaleakss - Jrelss;
        //cai
        residuals[8] += state.var(1,8) - JCa / (Vnonjunct3 * betai * 1000.0);
        //isusr //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[13] += state.var(1,13) - ( 1.0 / (1.0 + exp(-(state.vm() + get_IKur_ac_shift(state.cell_type(), this->mutation()) + 6.0) / (8.6 * get_IKur_ac_grad(state.cell_type(), this->mutation())))) - state.var(0,13) ) / ((0.009 / (1.0 + exp((state.vm() + 5.0) / 12.0)) + 0.0005)*1000.0);
        //isuss //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[14] +=state.var(1,14)- ( 1.0 / (1.0 + exp((state.vm() + get_IKur_inac_shift(state.cell_type(),this->mutation()) + 7.5) / (10.0 * get_IKur_inac_grad(state.cell_type(),this->mutation())))) - state.var(0,14) ) / ((0.59 / (1.0 + exp((state.vm() + 60.0) / 10.0)) + 3.05)*1000.0);
        //Cass //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[15] += state.var(1,15) - betass * ( JCass / Vss + ((-( get_RyR(state.cell_type(),this->mutation()) * state.ical()) - state.ibca() - state.icap() - state.icat() + 2 * state.inaca() - state.isac_ca())) / (2000.0 * Vss * F) )/1000.0;
        //CaSR1
        residuals[16] +=    state.var(1,16) - 
                                    betaSR1 *
                                    (
                                        DCaSR *
                                        (
                                            ( state.var(0,17) - state.var(0,16) ) / ( pow(dx,2) )
                                            +
                                            ( state.var(0,16) - state.var(0,17) ) / ( 6.0 * pow(dx,2) )
                                        )
                                        +
                                        JSRCa1 / VSR3
                                    );
        //CaSR2
        residuals[17] +=    state.var(1,17) -
                                    betaSR2 *
                                    (
                                        DCaSR *
                                        (
                                            ( state.var(0,16) - state.var(0,17) ) / ( pow(dx,2) )
                                            +
                                            ( state.var(0,17) - state.var(0,16) ) / ( 8.0 * pow(dx,2) )
                                        )
                                        +
                                        JSRCa2 / VSR4
                                    );
        //SERCACa //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[18] += state.var(1,18) - 0.5 * (-J_SERCASR + J_bulkSERCA) / (Vnonjunct3*1000.0);
        //SERCACass //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[19] += state.var(1,19) - 0.5 * (-J_SERCASRss + J_bulkSERCAss) / (Vss*1000.0);
        //RyRoss //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[20] += state.var(1,20) - ( RyRoinfss - state.var(0,20) ) / (RyRtauactss*1000.0);
        //RyRcss //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[21] += state.var(1,21) - ( RyRcinfss - state.var(0,21) ) / (RyRtauinactss*1000.0);
        //RyRass //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[22] += state.var(1,22) - ( RyRainfss - state.var(0,22) ) / (RyRtauadapt*1000.0);
        //RyRo3 //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[23] += state.var(1,23) - ( RyRoinf3 - state.var(0,23) ) / (RyRtauact*1000.0);
        //RyRc3 //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[24] += state.var(1,24) - ( RyRcinf3 - state.var(0,24) ) / (RyRtauinact*1000.0);
        //RyRa3 //!!!!! 1000.0 in tau is from CNZCell.cpp
        residuals[25] += state.var(1,25) - ( RyRainf3 - state.var(0,25) ) / (RyRtauadapt*1000.0);
    }

    //Contains two channel gating variables and potassium and sodium fibrosis concentrations
    //  and fibrosis transmembrane potential
    inline void CNZCell::fibrosis_residuals(const CellState &state, Vector<double> &residuals){    
        //preallocate the tau and inf variables
        double inf, tau;
        //rkv
        inf = (1 / (1 + exp(-(state.var(0,44) + 20 - 25 + FB_Ikv_shift) / 11)));
        tau = 1000 * (0.0203 + 0.1380 * exp(pow(-((state.var(0,44) + 20) / 25.9), 2)));
        residuals[40] += state.var(1,40) - (inf - state.var(0,40))/tau;

        //skv
        inf = 1 / (1 + exp((state.var(0,44) + 23 - 25 + FB_Ikv_shift) / 7));
        tau = 1000 * (1.574 + 5.268 * exp(pow(-((state.var(0,44) + 23) / 22.7), 2)));
        residuals[41] += state.var(1,41) - (inf - state.var(0,41))/tau;

        //kif
        //NO RESIDUAL IN ORIGINAL CNZCell model
        residuals[42] += state.var(1,42);//state.var(1,42);

        //naif
        //NO RESIDUAL IN ORIGINAL CNZCell model
        residuals[43] += state.var(1,43);//state.var(1,43);

        //Vmf
        residuals[44] += state.var(1,44) + state.ikv() + state.ik1f() /*+ state.inakf()*/ + state.ibnaf() - state.igap()/Cmf;
        //NO FUNCTON FOR INAKF IN ORIGINAL CODE
    }

    inline void CNZCell::rice_force_model_residuals(const CellState &state, Vector<double> &residuals){
        // BEGIN FORCE MODEL

        //Begin calculate variables for force model
        // Copied and only edited for variable names from original CNZCell
        //for readability preallocate the memory for the force variables
        double N =          state.var(0,31);
        double XBprer =     state.var(0,32);
        double XBpostr =    state.var(0,33);
        double SL =         state.var(0,34);
        double xXBpostr =   state.var(0,35);
        double xXBprer =    state.var(0,36);
        double TRPNCaL =    state.var(0,37);
        double TRPNCaH =    state.var(0,38);
        double intf =       state.var(0,39);
        //call the resting sarcomma length only once
        double SLset =      get_SLset(state.cell_type(),this->mutation());
        double SLrest =     get_SLrest(state.cell_type(),this->mutation());
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
        double konT    = kon /** pow(Qkon, ((Temp - 310.0) / 10.0))*/;
        double koffLT  = koffL /** pow(Qkoff, ((Temp - 310.0) / 10.0))*/ * koffmod;
        double koffHT  = koffH /** pow(Qkoff, ((Temp - 310.0) / 10.0))*/ * koffmod;
        double kn_pT   = kn_p * permtot /** pow(Qkn_p, ((Temp - 310.0) / 10.0))*/;
        double kp_nT   = kp_n * inprmt /** pow(Qkp_n, ((Temp - 310.0) / 10.0))*/;
        double fappT   = fapp * xbmodsp /** pow(Qfapp, ((Temp - 310.0) / 10.0))*/;
        double gapslmd = 1 + (1 - SOVFThick) * gslmod;
        double gappT   = gapp * gapslmd * xbmodsp /** pow(Qgapp, ((Temp - 310.0) / 10.0))*/;
        double hfmd    = exp(-rice_sign(xXBprer) * hfmdc * ((xXBprer / x_0) * (xXBprer / x_0)));
        double hbmd    = exp(rice_sign((xXBpostr - x_0)) * hbmdc * (((xXBpostr - x_0) / x_0) * ((xXBpostr - x_0) / x_0)));
        double hfT     = hf * hfmd * xbmodsp /** pow(Qhf, ((Temp - 310.0) / 10.0))*/;
        double hbT     = hb * hbmd * xbmodsp /** pow(Qhb, ((Temp - 310.0) / 10.0))*/;
        double gxbmd   = rice_heav(x_0 - xXBpostr) * exp(sigmap * ((x_0 - xXBpostr) / x_0) * ((x_0 - xXBpostr) / x_0))
                         + (1 - rice_heav(x_0 - xXBpostr)) * exp(sigman * (((xXBpostr - x_0) / x_0) * ((xXBpostr - x_0) / x_0)));
        double gxbT    = gxb * gxbmd * xbmodsp /** pow(Qgxb, ((Temp - 310.0) / 10.0))*/;
        //   Regulation and corssbridge cycling state derivatives
        double dTRPNCaL  = konT * state.var(0,8) * (1 - TRPNCaL) - koffLT * TRPNCaL;
        double dTRPNCaH  = konT * state.var(0,8) * (1 - TRPNCaH) - koffHT * TRPNCaH;
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
        double ppforce_c = rice_heav(SL - SL_c) * PCon_c * (exp(PExp_c * fabs(SL - SL_c)) - 1);
        double ppforce_t = rice_sign(SL - SLrest) * PCon_t *(exp(PExp_t *fabs((SL - SLrest))) - 1);
        // if (Rice_Usr_para.singlecell == False) ppforce_t += std::max(0, (PCon_c * exp(PExp_c * (SL - SL_c))));
        // Assume initial is generated by a preload force to counteract passive force
        // Preload force is computed here
        double PreloadF = rice_sign((SLset - SLrest)) * PCon_t *(exp(PExp_t *fabs((SLset - SLrest))) - 1.0);
        // if (Rice_Usr_para.singlecell == False) PreloadF += std::max(0, (PCon_c * exp(PExp_c * (SLset - SL_c))));
        double ppforce = ppforce_t + ppforce_c;
        // ppforce += std::max(0, (PCon_c * exp(PExp_c * (SL - SL_c))));
        double preload = rice_sign(SLset - SLrest) * PCon_t *(exp(PExp_t *fabs(SLset - SLrest)) - 1);
        // preload += std::max(0, (PCon_c * exp(PExp_c * (SLset - SL_c))));
        double afterload = 0;  // either static constant or due to series elastic element
        double dintf = /*(1 - Rice_Usr_para.SEon_LengthClamp) **/ (-ppforce + preload - active + afterload); //   total force
        double dSL = /*(1 - Rice_Usr_para.SEon_LengthClamp)**/((intf + (SLset - SL) * visc) / massf) * rice_heav(SL - SLmin) * rice_heav(SLmax - SL);
        //   Mean strain of strongly-bound states due to SL motion and XB cycling
        double dutyprer  = (hbT * fappT + gxbT * fappT)  //   duty fractions using the
                       / (fappT * hfT + gxbT * hfT + gxbT * gappT + hbT * fappT + hbT * gappT + gxbT * fappT);
        double dutypostr = fappT * hfT              //   King-Alman Rule
                       / (fappT * hfT + gxbT * hfT + gxbT * gappT + hbT * fappT + hbT * gappT + gxbT * fappT);
        double dxXBprer = dSL / 2.0 + xPsi / dutyprer * (-xXBprer * fappT + (xXBpostr - x_0 - xXBprer) * hbT);
        double dxXBpostr = dSL / 2.0 + xPsi / dutypostr * (x_0 + xXBprer - xXBpostr) * hfT;
        //FOR RETURNING dCai_feedback WHICH DOES NOTHING IN THE ORIGINAL CNZCell CODE
        //   Ca buffering by low-affinity troponin C (LTRPNCa)
        double FrSBXB    = (XBpostr + XBprer) / (SSXBpostr + SSXBprer);
        double dFrSBXB   = (dXBpostr + dXBprer) / (SSXBpostr + SSXBprer);
        double dsovr_ze  = -dSL / 2.0 * rice_heav(len_thick - SL);
        double dsovr_cle = -dSL / 2.0 * rice_heav((2.0 * len_thin - SL) - len_hbare);
        double dlen_sovr = dsovr_ze - dsovr_cle;
        double dSOVFThin = dlen_sovr / len_thin;
        double dSOVFThick = 2.0 * dlen_sovr / (len_thick - len_hbare);
        double TropToT = Trop_conc * ((1 - SOVFThin) * TRPNCaL
                                      + SOVFThin * (FrSBXB * TRPNCaH + (1 - FrSBXB) * TRPNCaL));
        double dTropToT = Trop_conc * (-dSOVFThin * TRPNCaL + (1 - SOVFThin) * dTRPNCaL
                                       + dSOVFThin * (FrSBXB * TRPNCaH + (1 - FrSBXB) * TRPNCaL)
                                       + SOVFThin * (dFrSBXB * TRPNCaH + FrSBXB * dTRPNCaH - dFrSBXB * TRPNCaL
                                               + (1 - FrSBXB) * dTRPNCaL));
        //END FOR RETURNING dCai_feedback WHICH DOES NOTHING IN THE ORIGINAL CNZCell CODE

        //rice_N
        residuals[31] += state.var(1,31) - dN;
        //rice_XBprer
        residuals[32] += state.var(1,32) - dXBprer;
        //rice_XBpostr
        residuals[33] += state.var(1,33) - dXBpostr;
        //rice_SL
        residuals[34] += state.var(1,34) - dSL;
        //rice_xXBpostr
        residuals[35] += state.var(1,35) - dxXBpostr;
        //rice_xXBprer
        residuals[36] += state.var(1,36) - dxXBprer;
        //rice_TRPNCaL
        residuals[37] += state.var(1,37) - dTRPNCaL;
        //rice_TRPNCaH
        residuals[38] += state.var(1,38) - dTRPNCaH;
        //rice_intf
        residuals[39] += state.var(1,39) - dintf;
        // END FORCE MODEL
    }



    





    //====================================================================
    //====================================================================
    // Fill in residuals
    //====================================================================
    //====================================================================
    void CNZCell::fill_in_generic_residual_contribution_cell_base(CellState &state,
                                                                    Vector<double> &residuals,
                                                                    DenseMatrix<double> &jacobian,
                                                                    unsigned flag)
    {   
        //Calculate reversal potentials
        ENa_reversal(state);
        ECa_reversal(state);
        EK_reversal(state);
        EKf_reversal(state);
        Enaf_reversal(state);

        //Calculate channel currents
        INa_current(state);
        IKr_current(state);
        IKs_current(state);
        ICaL_current(state);
        IK1_current(state);
        Iab_current(state);
        IbK_current(state);
        IbCa_current(state);
        IbNa_current(state);
        ICap_current(state);
        INaCa_current(state);
        INaK_current(state);
        Ito_current(state);
        IKur_current(state);
        If_current(state);
        ICaT_current(state);
        IGap_current(state);

        IK1f_current(state);
        IbNaf_current(state);
        IKv_current(state);

        ISAC_current(state);
        
        //Call the channel residual functions
        INa_current_residual(state, residuals);
        IKr_current_residual(state, residuals);
        IKs_current_residual(state, residuals);
        ICaL_current_residual(state, residuals);
        IK1_current_residual(state, residuals);
        Iab_current_residual(state, residuals);
        IbK_current_residual(state, residuals);
        IbCa_current_residual(state, residuals);
        IbNa_current_residual(state, residuals);
        ICap_current_residual(state, residuals);
        INaCa_current_residual(state, residuals);
        INaK_current_residual(state, residuals);
        Ito_current_residual(state, residuals);
        IKur_current_residual(state, residuals);
        If_current_residual(state, residuals);
        ICaT_current_residual(state, residuals);
        IGap_current_residual(state, residuals);
        ISAC_current_residual(state, residuals);

        //Call the fibrosis residual functions
        fibrosis_residuals(state, residuals);

        //Call the force model residual functions
        rice_force_model_residuals(state, residuals);

        //Call the current residual functions
        Ca_i_residual(state, residuals);
        Na_i_residual(state, residuals);
        K_i_residual(state, residuals);
        Calcium_dynamics_residual(state, residuals);
    }

    inline void CNZCell::return_initial_membrane_potential(double &v, const unsigned &cell_type){
        v = -76.079842;
    }

    inline bool CNZCell::return_initial_value(const unsigned &n, double &v, const unsigned &cell_type){
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
            case 28 : v = 0; // CNZ_a
                    break;
            case 29 : v = 1; // CNZ_i
                    break;
            case 30 : v = 0; // if y
                    break;
            case 31 : v = 0.97;  // N
                    break;
            case 32 : v = 0.01;   // XBprer
                    break;
            case 33 : v = 0.01;   // XBpostr
                    break;
            case 34 : v = 1.9;   // SL
                    break;
            case 35 : v = 0.007;   // xXBpostr
                    break;
            case 36 : v = 0.0;   // xXBprer
                    break;
            case 37 : v = 0.01447254;   // TropCaL
                    break;
            case 38 : v = 0.2320947;   // TropCaH
                    break;
            case 39 : v = 0.0;   // intergral of force, normallised
                    break;
            // FB data
            case 40 : v = 0.011694; // rkv
                    break;
            case 41 : v = 0.996878; // skv
                    break;
            case 42 : v = 129.434900; // kif
                    break;
            case 43 : v = 8.554700; // naif
                    break;
            case 44 : v = -43.806331; // vmf
                    break;

            default : return false;
        }
        return true;
    }

}