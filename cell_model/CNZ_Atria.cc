#include "cell_model_elements.h"

namespace oomph{

    CNZCell::CNZCell()   :   CellModelBase()
    {
        //The required storage
        //..without variables with zero residual
        this->Required_Storage = 40;


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


        //====================================================================
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
        //====================================================================
    }

    //====================================================================
    //====================================================================
    // Define the channel currents
    //====================================================================
    //====================================================================
    double CNZCell::INa_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_GNa(cell_type)* //!!!!!
                CRN_gna*                                //Blockade, originally from state[39] and state[40]
                pow(node_var(node,m_index_CNZCell(),local_ind),3)*
                node_var(node,h_index_CNZCell(),local_ind)*
                node_var(node,j_index_CNZCell(),local_ind)*
                (Vm - Rev_Pot);
    }

    double CNZCell::IKr_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_GKr(cell_type)*
                CRN_gkr*
                node_var(node,xr_index_CNZCell(),local_ind)*
                (Vm - Rev_Pot)/
                (1.0 + exp((Vm + 15.0) / 22.4));
    }

    double CNZCell::IKs_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_GKs(cell_type)*
                CRN_gks*
                pow(
                    node_var(node,xs_index_CNZCell(),local_ind)
                    ,2.0
                    )*
                (Vm - Rev_Pot);
    }

    double CNZCell::ICaL_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                2.125*
                get_GCaL(cell_type)*
                CRN_gcaL*
                node_var(node,d_index_CNZCell(),local_ind)*
                node_var(node,f_index_CNZCell(),local_ind)*
                node_var(node,fca_index_CNZCell(),local_ind)*
                (Vm - Rev_Pot);
    }

    double CNZCell::IK1_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_GK1(cell_type)*
                CRN_gk1*
                (Vm - Rev_Pot + get_IK1_v_shift(cell_type))/
                (1.0 + exp(0.07 * (Vm + 80.0 + get_IK1_v_shift(cell_type))));
    }

    double CNZCell::Iab_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                0.0003879*
                (Vm + 69.6)/
                (1.0 - 0.8377*exp((Vm + 49.06)/1056.0));
    }

    double CNZCell::IbK_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                CRN_gbk*
                (Vm - Rev_Pot);
    }

    double CNZCell::IbCa_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_Gbca(cell_type)*
                CRN_gbca*
                (Vm - Rev_Pot);
    }

    double CNZCell::IbNa_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                CRN_gbna*
                (Vm - Rev_Pot);
    }

    double CNZCell::ICap_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                1.4*
                get_GCap(cell_type)*
                CRN_icapbar*
                node_var(node, Cass_index_CNZCell(), local_ind)/
                (node_var(node, Cass_index_CNZCell(),local_ind) + CRN_kmcap );
    }

    double CNZCell::INaCa_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_GNaCa(cell_type)*
                CRN_knacalr/
                (pow(CRN_kmnalr,3.0)+pow(Ext_conc[0],3.0))/
                (CRN_kmcalr+Ext_conc[1])/
                (1.0+CRN_ksatlr*exp((CRN_gammalr-1.0)*Vm*FoRT))*
                (pow(node_var(node,nai_index_CNZCell(),local_ind),3.0)*Ext_conc[1]*exp(Vm*CRN_gammalr*FoRT)-
                    pow(Ext_conc[0],3.0)*node_var(node, Cass_index_CNZCell(),local_ind)*exp(Vm*(CRN_gammalr-1.0)*FoRT));
    }

    double CNZCell::INaK_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                CRN_inakbar/
                (
                    (1.0+0.1245*exp(-0.1*Vm*FoRT)+0.0365*((exp(Ext_conc[0]/67.3)-1.0)/7.0)*exp(-Vm*FoRT))*
                    (1.0 + CRN_kmko / Ext_conc[2])*
                    (1.0 + pow(CRN_kmnai/node_var(node,nai_index_CNZCell(),local_ind),1.5))
                );
    }

    double CNZCell::Ito_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_Gto(cell_type)*
                MAL_gto*
                node_var(node,itr_index_CNZCell(),local_ind)*
                node_var(node,its_index_CNZCell(),local_ind)*
                (Vm - Rev_Pot);
    }

    double CNZCell::IKur_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  Cm*
                get_IKur_cond(cell_type,mut_type)*
                get_GKur(cell_type,mut_type)*
                CNZ_gkur*
                (get_IKur_c(cell_type,mut_type) + get_IKur_x0(cell_type,mut_type)/(1.0 + exp((Vm - get_Ikur_y0(cell_type,mut_type))/(-8.26597))))*
                node_var(node,CNZ_a_index_CNZCell(),local_ind)*
                node_var(node,CNZ_i_index_CNZCell(),local_ind)*
                (Vm - Rev_Pot);
    }

    double CNZCell::If_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  0.0;
                /*
                Cm*
                get_Gf(get_cell_type_at_node_CNZ(n))*
                CRN_gf*
                node_var(node,If_y_index_CNZ(),local_ind)*
                (Vm - Rev_Pot);
                */
    }

    double CNZCell::ICaT_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return  0.0;
                /*
                Cm*
                get_GCaT(cell_type)*
                gcaT*
                node_var(node,ff_index_CNZCell(),local_ind)*
                node_var(node,dd_index_CNZCell(),local_ind)*
                (Vm - EcaT);
                */
    }

    //!!!!!
    //MISSING PARAMETER GGAP IN ORIGINAL CNZCell.cpp

    double CNZCell::IGap_current_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &Rev_Pot) const {
        return 0.0;
    }


    double CNZCell::ISAC_ISAC_CNZCell(Node* node,
                                        const Vector<double> &Ext_conc,
                                        const Vector<unsigned> &local_ind,
                                        const unsigned &cell_type,
                                        const unsigned &mut_type,
                                        const unsigned &fibrosis,
                                        const double &Vm,
                                        const double &strain) const {
        return  Cm*
                ISAC_GSac*
                (Vm - ISAC_Esac)/
                (1.0 + exp(-(strain - ISAC_strain_half)/ISAC_Ke));
    }

    double CNZCell::ISAC_Na_current_CNZCell(Node* node,
                                            const Vector<double> &Ext_conc,
                                            const Vector<unsigned> &local_ind,
                                            const unsigned &cell_type,
                                            const unsigned &mut_type,
                                            const unsigned &fibrosis,
                                            const double &Vm,
                                            const double &ISAC_) const {
        return  ISAC_pNa*
                ISAC_/
                ISAC_total_Sca;

    }

    double CNZCell::ISAC_K_current_CNZCell(Node* node,
                                            const Vector<double> &Ext_conc,
                                            const Vector<unsigned> &local_ind,
                                            const unsigned &cell_type,
                                            const unsigned &mut_type,
                                            const unsigned &fibrosis,
                                            const double &Vm,
                                            const double &ISAC_) const {
        return  ISAC_pK*
                ISAC_/
                ISAC_total_Sca;
    }

    double CNZCell::ISAC_Ca_current_CNZCell(Node* node,
                                            const Vector<double> &Ext_conc,
                                            const Vector<unsigned> &local_ind,
                                            const unsigned &cell_type,
                                            const unsigned &mut_type,
                                            const unsigned &fibrosis,
                                            const double &Vm,
                                            const double &ISAC_) const {
        return  ISAC_pCa*
                ISAC_/
                ISAC_total_Sca;
    }

    //====================================================================
    //====================================================================
    // Calculate the sub residual for CNZCell
    //====================================================================
    //====================================================================
    void CNZCell::fill_in_generic_residual_contribution_cell_base( Node* node,
                            const double& Vm,
                            const double& strain,
                            const Vector<double> &Ext_conc,
                            const Vector<unsigned> &local_ind,
                            const unsigned &cell_type,
                            const unsigned &mut_type,
                            const unsigned &fibrosis,
                            Vector<double> &residuals,
                            DenseMatrix<double> &jacobian,
                            unsigned flag)
    {   
        // std::cout << "entered residual in cell_model" << std::endl;
        //====================================================================
        //====================================================================
        // Calculate some variables before residuals is calculated
        //====================================================================
        //====================================================================
        //Preallocate memory for the reversal potentials
        double Ena, Eca, Ek;
        //Preallocate memory for the ISAC scale
        double ISAC_;
        //Preallocate memory for the channel currents used when calcualting concentration changes (all of them)
        double INa, IKr, IKs, ICaL, IK1, Iab, IbK, IbCa, IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT, IGap, ISAC_Na, ISAC_K, ISAC_Ca;
        //Preallocate memory for alpha and beta
        double alpha, beta;
        //Calculate reversal potentials
        get_reversal_Na(Ext_conc[0], node_var(node,nai_index_CNZCell(),local_ind), Ena);
        get_reversal_Ca(Ext_conc[1], node_var(node,cai_index_CNZCell(),local_ind), Eca);
        get_reversal_K(Ext_conc[2], node_var(node,ki_index_CNZCell(),local_ind), Ek);

        //Get the ISAC scale
        ISAC_ = ISAC_ISAC_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, strain);

        //Calculate channel currents
        INa = INa_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ena);
        IKr = IKr_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ek);
        IKs = IKs_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ek);
        ICaL = ICaL_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, CRN_ErL);
        IK1 = IK1_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ek);
        Iab = Iab_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);          //!!!!! NO NERNST POTENTIAL USED IN CNZ CODE
        IbK = IbK_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ek);
        IbCa = IbCa_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Eca);
        IbNa = IbNa_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ena);
        ICap = ICap_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);        //!!!!! NO NERNST POTENTIAL USED IN CNZ CODE
        INaCa = INaCa_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);      //!!!!! NO NERNST POTENTIAL USED IN CNZ CODE
        INaK = INaK_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);        //!!!!! NO NERNST POTENTIAL USED IN CNZ CODE
        Ito = Ito_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ek);
        IKur = IKur_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, Ek);
        If = If_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);            //!!!!! NO NERNST POTENTIAL USED IN CNZ CODE
        ICaT = ICaT_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);        //!!!!! NO NERNST POTENTIAL USED IN CNZ CODE

        IGap = IGap_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, 0.0);        //!!!!! NOT IMPLEMENTED YET

        ISAC_Na = ISAC_Na_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, ISAC_);
        ISAC_K = ISAC_K_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, ISAC_);
        ISAC_Ca = ISAC_Ca_current_CNZCell(node, Ext_conc, local_ind, cell_type, mut_type, fibrosis, Vm, ISAC_);



        //====================================================================
        //====================================================================
        // Loop through the local equations
        //====================================================================
        //====================================================================

        //====================================================================
        //m
        //====================================================================
        int var_ind = m_index_CNZCell();
        //Calculate rates
        if(std::fabs(Vm + 47.13) > 1e-10){alpha = 0.32 * (Vm + 47.13) / (1.0 - exp(-0.1 * (Vm + 47.13)));}
        else{alpha = 3.2;}
        beta = 0.08 * exp(-Vm / 11.0);
        //Populate residuals
        // residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) + node_var(node,var_ind,local_ind)*(alpha + beta) - alpha;

        // Populate residual with step from previous value to next one, to be used with identity jacobian
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) + node_var(node,var_ind,local_ind)*(alpha + beta) - alpha;
        //Populate jacobian
        if(flag){
            // jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0) + alpha + beta;

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0);
        }
        //====================================================================
        //h
        //====================================================================
        var_ind = h_index_CNZCell();
        //Calculate rates
        if(Vm>=-40.0){
            alpha = 0.0;
            beta = 1.0 / (0.13 * (1.0 + exp((Vm + 10.66) / -11.1)));
        }
        else{
            alpha = 0.135 * exp((Vm + 80.0) / -6.8);
            beta = 3.56 * exp(0.079 * Vm) + 3.1e5 * exp(0.35 * Vm);
        }
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) + node_var(node,var_ind,local_ind)*(alpha + beta) - alpha;

        // Populate residual with step from previous value to next one, to be used with identity jacobian
        // residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) + node_var(node,var_ind,local_ind)*(alpha + beta) - alpha
        //                         + ( 1 - node->time_stepper_pt()->weight(1,0) )*node->value(0,local_ind[var_ind]) - ( 1 + node->time_stepper_pt()->weight(1,1) )*node->value(1,local_ind[var_ind]);

        //Populate jacobian
        if(flag){
            // jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0) + alpha + beta;

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //j
        //====================================================================
        var_ind = j_index_CNZCell();
        //Calculate rates
        if (Vm >= -40.0)
        {
            alpha  = 0.0;
            beta = 0.3 * exp(-2.535e-7 * Vm) / (1.0 + exp(-0.1 * (Vm + 32.0)));
        }
        else
        {
            alpha = (-1.2714e5 * exp(0.2444 * Vm) - 3.474e-5 * exp(-0.04391 * Vm)) * (Vm + 37.78) / (1.0 + exp(0.311 * (Vm + 79.23)));
            beta = 0.1212 * exp(-0.01052 * Vm) / (1.0 + exp(-0.1378 * (Vm + 40.14)));
        }
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) + node_var(node,var_ind,local_ind)*(alpha + beta) - alpha;

        // Populate residual with step from previous value to next one, to be used with identity jacobian
        // residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) + node_var(node,var_ind,local_ind)*(alpha + beta) - alpha
        //                         + ( 1 - node->time_stepper_pt()->weight(1,0) )*node->value(0,local_ind[var_ind]) - ( 1 + node->time_stepper_pt()->weight(1,1) )*node->value(1,local_ind[var_ind]);

        //Populate jacobian
        if(flag){
            // jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0) + alpha + beta;

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0);
        }
        //====================================================================
        //d
        //====================================================================
        var_ind = d_index_CNZCell();
        //Populate residuals and jacobian without rates         
        if (std::fabs(Vm + 10.0) > 1e-10)
        {
            //Populate residuals
            residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - 0.035*(Vm + 10.0)*( (1.0+exp(-(Vm+10.0)/6.24)) / (1.0-exp(-(Vm+10.0)/6.24)) )*(1.0/(1.0+exp(-(Vm+10.0)/8.0)) - node_var(node,var_ind,local_ind));
            //Populate jacobian
            if(flag){
               // jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0) + 0.035*(Vm+10.0)*(1.0+exp(-(Vm+10.0)/8.0))/(1.0-exp(-(Vm+10.0)/6.24));   

               // Jacobian is identity
                jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
           }
        }
        else
        {   
            //Populate residuals
            residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - (1.0/4.579) * (1.0 + exp((Vm + 10.0) / -6.24)) * ( 1.0 / (1.0 + exp((Vm + 10.0) / -8.0)) - node_var(node,var_ind,local_ind));
            //Populate jacobian
            if(flag){
                // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) + 5.60897e-3;

                // Jacobian is identity
                jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
            }
        }

        //====================================================================
        //f
        //====================================================================
        var_ind = f_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( exp(-(Vm+28.0)/6.9)/(1.0+exp(-(Vm+28.0)/6.9)) - node_var(node,var_ind,local_ind) ) / (9.0 / (0.0197 * exp(-pow(0.0337*(Vm + 10),2.0)) + 0.02));
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) - 1.0/(1.5 * 2.0 * 3.0 / (0.0197 * exp(-0.0337 * 0.0337 * (Vm + 10.0) * (Vm + 10.0)) + 0.02));

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //xr
        //====================================================================
        var_ind = xr_index_CNZCell();
        //Calculate rates
        if (std::fabs(Vm + 14.1) > 1e-10){alpha = 0.0003 * (Vm + 14.1) / (1.0 - exp(-(Vm + 14.1) / 5.0));}
        else{alpha = 0.0015;}
        if (fabs(Vm - 3.3328) > 1e-10){beta = 0.000073898 * (Vm - 3.3328) / (exp((Vm - 3.3328) / 5.1237) - 1.0);}
        else{beta = 3.7836118e-4;}
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / ( 1.0 + exp( -(Vm + get_IKr_ac_shift(cell_type, mut_type) + 14.1) / (6.5 * get_IKr_ac_grad(cell_type, mut_type))) ) - node_var(node,var_ind,local_ind) )*(alpha + beta);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) - node_var(node,var_ind,local_ind)*(alpha + beta);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //xs
        //====================================================================
        var_ind = xs_index_CNZCell();
        //Calculate rates
        if (std::fabs(Vm - 19.9) > 1e-10){
            
            alpha = 0.00004 * (Vm - 19.9) / (1.0 - exp(-(Vm - 19.9) / 17.0));
            beta = 0.000035 * (Vm - 19.9) / (exp((Vm - 19.9) / 9.0) - 1.0);
        }
        else{
            alpha = 0.00068;
            beta = 0.000315;
        }
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( sqrt(1.0 / (1.0 + exp(-(Vm - 19.9 - get_IKs_shift(cell_type, mut_type)) / (12.7 * get_IKs_grad(cell_type, mut_type))))) - node_var(node,var_ind,local_ind) )*2.0*(alpha + beta);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) - node_var(node,var_ind,local_ind)*2.0*(alpha + beta);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //nai
        //====================================================================
        var_ind = nai_index_CNZCell();
        
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - (-3.0 * INaK - 3.0 * INaCa - IbNa - INa - ISAC_Na) / (F *CRN_vi);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //ki
        //====================================================================
        var_ind = ki_index_CNZCell();
        
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - (2.0 * INaK - IK1 - Ito - IKur - IKr - IKs - IbK - ISAC_K) / (F * CRN_vi);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0); //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //values required for changes in various calcium concentrations
        //====================================================================
        double betass;
        betass = pow(
                     ( 1 + SLlow * KdSLlow / pow((node_var(node,Cass_index_CNZCell(),local_ind) + KdSLlow), 2)
                       + SLhigh * KdSLhigh / pow((node_var(node,Cass_index_CNZCell(),local_ind) + KdSLhigh), 2)
                       + BCa * KdBCa / pow(((node_var(node,Cass_index_CNZCell(),local_ind)) + KdBCa), 2)  )
                     , (-1));

        double betai, gammai;
        betai = pow(( 1 + BCa * KdBCa / pow((node_var(node,cai_index_CNZCell(),local_ind) + KdBCa), 2)  ), (-1));
        gammai = BCa * KdBCa / pow((node_var(node,cai_index_CNZCell(),local_ind) + KdBCa), 2);

        double betaSR1, betaSR2;
        betaSR1 = pow( ( 1 + CSQN * KdCSQN / pow((node_var(node,CaSR1_index_CNZCell(),local_ind) + KdCSQN), 2) ), (-1));
        betaSR2 = pow( ( 1 + CSQN * KdCSQN / pow((node_var(node,CaSR2_index_CNZCell(),local_ind) + KdCSQN), 2) ), (-1));

        double Jj_nj;
        Jj_nj = DCa * Aj_nj / xj_nj * (node_var(node,Cass_index_CNZCell(),local_ind) - node_var(node,cai_index_CNZCell(),local_ind)) * 1e-6;

        double J_SERCASR, J_bulkSERCA;
        double J_SERCASRss, J_bulkSERCAss;

        J_SERCASR =  (-k3 * pow(node_var(node,CaSR1_index_CNZCell(),local_ind), 2) * (cpumps - node_var(node,SERCACa_index_CNZCell(),local_ind)) + k4 * node_var(node,SERCACa_index_CNZCell(),local_ind)) * Vnonjunct3 * 2;
        J_bulkSERCA = (k1 * pow(node_var(node,cai_index_CNZCell(),local_ind), 2) * (cpumps - node_var(node,SERCACa_index_CNZCell(),local_ind)) - k2 * node_var(node,SERCACa_index_CNZCell(),local_ind)) * Vnonjunct3 * 2;
        J_SERCASRss = (-k3 * pow(node_var(node,CaSR2_index_CNZCell(),local_ind), 2) * (cpumps - node_var(node,SERCACass_index_CNZCell(),local_ind)) + k4 * node_var(node,SERCACass_index_CNZCell(),local_ind)) * Vss * 2;
        J_bulkSERCAss = (k1 * pow((node_var(node,Cass_index_CNZCell(),local_ind)), 2) * (cpumps - node_var(node,SERCACass_index_CNZCell(),local_ind)) - k2 * node_var(node,SERCACass_index_CNZCell(),local_ind)) * Vss * 2;

        double RyRtauadapt = 1.0;
        double RyRtauactss = 5e-3;
        double RyRtauinactss = 15e-3;
        double RyRtauact = 18.75e-3;
        double RyRtauinact = 87.5e-3;

        double nuss = 625.0 * Vss;
        double RyRSRCass = (1 - 1 / (1 +  exp((node_var(node,CaSR2_index_CNZCell(),local_ind) - 0.3) / 0.1)));
        double RyRainfss = 0.505 - 0.427 / (1 + exp(( ( node_var(node,Cass_index_CNZCell(),local_ind) + (get_fRyR(cell_type, mut_type) * node_var(node,Cass_index_CNZCell(),local_ind)) ) * 1000 - 0.29) / 0.082));
        double RyRoinfss = (1 - 1 / (1 +  exp(((node_var(node,Cass_index_CNZCell(),local_ind) + (get_fRyR(cell_type, mut_type) * node_var(node,Cass_index_CNZCell(),local_ind))   ) * 1000 - ((node_var(node,RyRass_index_CNZCell(),local_ind)) + 0.22)) / 0.03)));
        double RyRcinfss = (1 / (1 + exp(((node_var(node,Cass_index_CNZCell(),local_ind) + (get_fRyR(cell_type, mut_type) * node_var(node,Cass_index_CNZCell(),local_ind) )) * 1000 - ((node_var(node,RyRass_index_CNZCell(),local_ind)) + 0.02)) / 0.01)));
        double Jrelss = nuss * ( (node_var(node,RyRoss_index_CNZCell(),local_ind)) ) * (node_var(node,RyRcss_index_CNZCell(),local_ind)) * RyRSRCass * ( node_var(node,CaSR2_index_CNZCell(),local_ind) -  (node_var(node,Cass_index_CNZCell(),local_ind)) );

        double nu3 = 1 * Vnonjunct3;
        double RyRSRCa3 = (1 - 1 / (1 +  exp((node_var(node,CaSR1_index_CNZCell(),local_ind) - 0.3) / 0.1)));
        double RyRainf3 =  0.505 - 0.427 / (1 + exp(( (node_var(node,cai_index_CNZCell(),local_ind) + ( get_fRyR(cell_type, mut_type) * node_var(node,cai_index_CNZCell(),local_ind))  ) * 1000 - 0.29) / 0.082));
        double RyRoinf3 = (1 - 1 / (1 +  exp(( (node_var(node,cai_index_CNZCell(),local_ind) + ( get_fRyR(cell_type, mut_type) * node_var(node,cai_index_CNZCell(),local_ind)) ) * 1000 - ((node_var(node,RyRa3_index_CNZCell(),local_ind)) + 0.22)) / 0.03)));
        double RyRcinf3 = (1 / (1 +  exp(( (node_var(node,cai_index_CNZCell(),local_ind) + (get_fRyR(cell_type, mut_type) * node_var(node,cai_index_CNZCell(),local_ind) ) ) * 1000 - ((node_var(node,RyRa3_index_CNZCell(),local_ind)) + 0.02)) / 0.01)));
        double Jrel3 = nu3 * ( (node_var(node,RyRo3_index_CNZCell(),local_ind)) ) * (node_var(node,RyRc3_index_CNZCell(),local_ind)) * RyRSRCa3 * ( node_var(node,CaSR1_index_CNZCell(),local_ind) -  node_var(node,cai_index_CNZCell(),local_ind) );

        Jrelss = get_fIRel(cell_type, mut_type) * Jrelss;
        Jrel3 = get_fIRel(cell_type, mut_type) * Jrel3;

        double JSRCaleak3 = get_GSR_leak(cell_type, mut_type) * kSRleak * ( node_var(node,CaSR1_index_CNZCell(),local_ind) - node_var(node,cai_index_CNZCell(),local_ind) ) * Vnonjunct3;
        double JSRCaleakss = get_GSR_leak(cell_type, mut_type) * kSRleak * ( node_var(node,CaSR2_index_CNZCell(),local_ind) - (node_var(node,Cass_index_CNZCell(),local_ind)) ) * Vss;

        double JCa, JCass;
        JCa = -get_BULK_CONST(cell_type, mut_type) * J_bulkSERCA + JSRCaleak3 + Jrel3 + Jj_nj;
        JCass = -Jj_nj + JSRCaleakss - get_BULK_CONST(cell_type, mut_type) * J_bulkSERCAss + Jrelss;

        double JSRCa1, JSRCa2;
        JSRCa1 = J_SERCASR - JSRCaleak3 - Jrel3;
        JSRCa2 = J_SERCASRss - JSRCaleakss - Jrelss;


        //====================================================================
        //cai
        //====================================================================
        var_ind = cai_index_CNZCell();

        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - JCa / (Vnonjunct3 * betai * 1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //fca
        //====================================================================
        var_ind = fca_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + (node_var(node,Cass_index_CNZCell(),local_ind) / 0.00035)) - node_var(node,var_ind,local_ind) )/(2.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) + 0.5;

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        // jacobian(var_ind,node_var(node,Cass_index_CNZCell(),local_ind)) += 0.00035/(1 + pow(node_var(node,Cass_index_CNZCell(),local_ind),2.0));

        //====================================================================
        //itr
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = itr_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp(-(Vm - 1.0)/11.0)) - node_var(node,var_ind,local_ind) ) / ((0.0035 * exp(-(Vm / 15.0)) + 0.0015)*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) + 1.0/(0.0035 * exp(-(Vm / 30.0) * 2.0) + 0.0015);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //its
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = its_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp((Vm + 40.5) / 11.5)) - node_var(node,var_ind,local_ind) ) / ((0.025635 * exp (-(Vm + 52.45) / 7.94135) + 0.01414)*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) + 1.0/(0.025635 * exp (-((Vm + 52.45) / 15.8827) * 2.0) + 0.01414);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //isusr
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = isusr_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp(-(Vm + get_IKur_ac_shift(cell_type, mut_type) + 6.0) / (8.6 * get_IKur_ac_grad(cell_type, mut_type)))) - node_var(node,var_ind,local_ind) ) / ((0.009 / (1.0 + exp((Vm + 5.0) / 12.0)) + 0.0005)*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) + 1.0/(0.009 / (1.0 + exp((Vm + 5.0) / 12.0)) + 0.0005);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //isuss
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = isuss_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp((Vm + get_IKur_inac_shift(cell_type,mut_type) + 7.5) / (10.0 * get_IKur_inac_grad(cell_type,mut_type)))) - node_var(node,var_ind,local_ind) ) / ((0.59 / (1.0 + exp((Vm + 60.0) / 10.0)) + 3.05)*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0) + 1.0/(0.59 / (1.0 + exp((Vm + 60.0) / 10.0)) + 3.05);

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //Cass
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = Cass_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - betass * ( JCass / Vss + ((-( get_RyR(cell_type,mut_type) * ICaL) - IbCa - ICap - ICaT + 2 * INaCa - ISAC_Ca)) / (2000.0 * Vss * F) )/1000.0;
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //CaSR1
        //====================================================================
        var_ind = CaSR1_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] +=    node_var_derivative(node, var_ind, local_ind) - 

                                    betaSR1 *
                                    (
                                        DCaSR *
                                        (
                                            ( node_var(node,CaSR2_index_CNZCell(),local_ind) - node_var(node,var_ind,local_ind) ) / ( pow(dx,2) )
                                            +
                                            ( node_var(node,var_ind,local_ind) - node_var(node,CaSR2_index_CNZCell(),local_ind) ) / ( 6.0 * pow(dx,2) )
                                        )
                                        +
                                        JSRCa1 / VSR3
                                    );
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }   
        //====================================================================
        //CaSR2
        //====================================================================
        var_ind = CaSR2_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] +=    node_var_derivative(node, var_ind, local_ind) -

                                    betaSR2 *
                                    (
                                        DCaSR *
                                        (
                                            ( node_var(node,CaSR1_index_CNZCell(),local_ind) - node_var(node,var_ind,local_ind) ) / ( pow(dx,2) )
                                            +
                                            ( node_var(node,var_ind,local_ind) - node_var(node,CaSR1_index_CNZCell(),local_ind) ) / ( 8.0 * pow(dx,2) )
                                        )
                                        +
                                        JSRCa2 / VSR4
                                    );
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //SERCACa
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = SERCACa_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - 0.5 * (-J_SERCASR + J_bulkSERCA) / (Vnonjunct3*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //SERCACass
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = SERCACass_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - 0.5 * (-J_SERCASRss + J_bulkSERCAss) / (Vss*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //RyRoss
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = RyRoss_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( RyRoinfss - node_var(node,var_ind,local_ind) ) / (RyRtauactss*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //RyRcss
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = RyRcss_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( RyRcinfss - node_var(node,var_ind,local_ind) ) / (RyRtauinactss*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //RyRass
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = RyRass_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( RyRainfss - node_var(node,var_ind,local_ind) ) / (RyRtauadapt*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //RyRo3
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = RyRo3_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( RyRoinf3 - node_var(node,var_ind,local_ind) ) / (RyRtauact*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //RyRc3
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = RyRc3_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( RyRcinf3 - node_var(node,var_ind,local_ind) ) / (RyRtauinact*1000.0);
        //Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //RyRa3
        //==================================================================== //!!!!! 1000.0 in tau is from CNZCell.cpp
        var_ind = RyRa3_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( RyRainf3 - node_var(node,var_ind,local_ind) ) / (RyRtauadapt*1000.0);
        // Populate jacobian
        if(flag){
            // jacobian(var_ind,var_ind) += node->time_stepper_pt()->weight(1,0);    //!!!!! COMPLETE THIS

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }

        //====================================================================
        //dd
        //====================================================================
        var_ind = dd_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp(-(Vm + 37.0) / 6.8)) - node_var(node,var_ind,local_ind) ) * ( 1.0680*exp((Vm + 26.3)/30.0) + 1.0680*exp(-(Vm + 26.3)/30.0) );
        // //Populate jacobian
        if(flag){
            // jacobian(var_ind, var_ind) += 1.0;

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        // //====================================================================
        // //ff
        // //====================================================================
        var_ind = ff_index_CNZCell();
        //Populate residuals
        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp((Vm + 71.0) / 9.0)) - node_var(node,var_ind,local_ind) ) * ( 0.0153*exp(-(Vm + 71.0)/83.3) + 0.0150*exp((Vm + 71.0)/15.38) );
        // //Populate jacobian
        if(flag){
            // jacobian(var_ind, var_ind) += 1.0;

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        
        // //====================================================================
        // //rkv
        // //====================================================================
        // var_ind = rkv_index_CNZCell();
        // //Populate residuals
        // residuals[local_ind[var_ind]] += 0.0;
        // //Populate jacobian
        // // jacobian(var_ind, var_ind) += 1.0;

        // //====================================================================
        // //skv
        // //====================================================================
        // var_ind = skv_index_CNZCell();
        // //Populate residuals
        // residuals[local_ind[var_ind]] += 0.0;
        // //Populate jacobian
        // // jacobian(var_ind, var_ind) += 1.0;

        // //====================================================================
        // //kif
        // //====================================================================
        // var_ind = kif_index_CNZCell();
        // //Populate residuals
        // residuals[local_ind[var_ind]] += 0.0;
        // //Populate jacobian
        // // jacobian(var_ind, var_ind) += 1.0;

        // //====================================================================
        // //naif
        // //====================================================================
        // var_ind = naif_index_CNZCell();
        // //Populate residuals
        // residuals[local_ind[var_ind]] += 0.0;
        // //Populate jacobian
        // // jacobian(var_ind, var_ind) += 1.0;

        // //====================================================================
        // //Vmf
        // //====================================================================
        // var_ind = Vmf_index_CNZCell();
        // //Populate residuals
        // residuals[local_ind[var_ind]] += 0.0;
        // //Populate jacobian
        // // jacobian(var_ind, var_ind) += 1.0;

        //====================================================================
        //CNZ_a
        //====================================================================
        var_ind = CNZ_a_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp(-(Vm - (-6.0 + get_IKur_Vhchange(cell_type, mut_type) )) / (8.6 * get_IKur_slope(cell_type,mut_type) ))) - node_var(node,var_ind,local_ind) ) / ( ( get_IKur_timeconstants(cell_type,mut_type)*(45.6666746826 / (1.0 + exp((Vm + 11.2306497073) / 11.5254705962)) + 4.26753514993)*(0.262186042981 / (1.0 + exp((Vm + 35.8658312707) / (-3.87510627762))) + 0.291755017928) ) / 3.5308257834747638 );
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //CNZ_i
        //====================================================================
        var_ind = CNZ_i_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( (get_IKur_inac_mult(cell_type,mut_type) * 0.52424) / (1.0 + exp((Vm + 15.1142 + get_IKur_inac_shift(cell_type,mut_type) ) / (7.567021 * get_IKur_inac_grad(cell_type,mut_type) ))) + 0.4580778 + get_IKur_inac_add(cell_type,mut_type)  - node_var(node,var_ind,local_ind) ) / ( (2328.0 / (1.0 + exp((Vm - 9.435) / (3.5827))) + 1739.139) / 3.5308257834747638 );
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //If_y
        //====================================================================
        var_ind = If_y_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - ( 1.0 / (1.0 + exp((Vm + 90.95 + get_If_vshift(cell_type,mut_type)) / (10.1 * get_If_grad(cell_type,mut_type))) ) - node_var(node,var_ind,local_ind) ) * (1.2783E-04 * exp(-Vm / 9.2424) + 121.6092 * exp(Vm / 9.2424));
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) +=  node->time_stepper_pt()->weight(1,0); 
        }

        
        //====================================================================
        //====================================================================
        // BEGIN FORCE MODEL
        //====================================================================
        //====================================================================

        //====================================================================
        //Begin calculate variables for force model
        // Copied and only edited for variable names from original CNZCell
        //====================================================================
        //for readability preallocate the memory for the force variables
        double N =          node_var(node,rice_N_index_CNZCell(),local_ind);
        double XBprer =     node_var(node,rice_XBprer_index_CNZCell(),local_ind);
        double XBpostr =    node_var(node,rice_XBpostr_index_CNZCell(),local_ind);
        double SL =         node_var(node,rice_SL_index_CNZCell(),local_ind);
        double xXBpostr =   node_var(node,rice_xXBpostr_index_CNZCell(),local_ind);
        double xXBprer =    node_var(node,rice_xXBprer_index_CNZCell(),local_ind);
        double TRPNCaL =    node_var(node,rice_TRPNCaL_index_CNZCell(),local_ind);
        double TRPNCaH =    node_var(node,rice_TRPNCaH_index_CNZCell(),local_ind);
        double intf =       node_var(node,intf_index_CNZCell(),local_ind);

        //call the resting sarcomma length only once
        double SLset =      get_SLset(cell_type,mut_type);
        double SLrest =     get_SLrest(cell_type,mut_type);

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
        double dTRPNCaL  = konT * node_var(node,cai_index_CNZCell(),local_ind) * (1 - TRPNCaL) - koffLT * TRPNCaL;
        double dTRPNCaH  = konT * node_var(node,cai_index_CNZCell(),local_ind) * (1 - TRPNCaH) - koffHT * TRPNCaH;
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

        //====================================================================
        //rice_N
        //====================================================================
        var_ind = rice_N_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dN;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_XBprer
        //====================================================================
        var_ind = rice_XBprer_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dXBprer;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_XBpostr
        //====================================================================
        var_ind = rice_XBpostr_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dXBpostr;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_SL
        //====================================================================
        var_ind = rice_SL_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dSL;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_xXBpostr
        //====================================================================
        var_ind = rice_xXBpostr_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dxXBpostr;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_xXBprer
        //====================================================================
        var_ind = rice_xXBprer_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dxXBprer;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_TRPNCaL
        //====================================================================
        var_ind = rice_TRPNCaL_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dTRPNCaL;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_TRPNCaH
        //====================================================================
        var_ind = rice_TRPNCaH_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dTRPNCaH;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //rice_intf
        //====================================================================
        var_ind = intf_index_CNZCell();

        residuals[local_ind[var_ind]] += node_var_derivative(node, var_ind, local_ind) - dintf;
        if(flag){

            // Jacobian is identity
            jacobian(var_ind, var_ind) += node->time_stepper_pt()->weight(1,0); 
        }
        //====================================================================
        //====================================================================
        // END FORCE MODEL
        //====================================================================
        //====================================================================
    }
}