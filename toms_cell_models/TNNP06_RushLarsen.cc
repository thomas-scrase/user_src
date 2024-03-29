#include "TNNP06_RushLarsen.h"


namespace oomph{

TNNP06RushLarsen::TNNP06RushLarsen(const unsigned& number_of_backup_values) : CellModelBaseFullyPartitioned(number_of_backup_values)
{
	cell_type = 100;


	Tent_R = 8314.472;
	Tent_F = 96485.3415;
	Tent_T = 310.0;
	Tent_Ko = 5.4;
	Tent_Cao = 2.0;
	Tent_Nao = 140.0;
	Tent_Vc = 0.016404;
	Tent_Vsr = 0.001094;
	Tent_Vss = 0.00005468;
	Tent_Bufc = 0.2;
	Tent_Kbufc = 0.001;
	Tent_Bufsr = 10.;
	Tent_Kbufsr = 0.3;
	Tent_Bufss = 0.4;
	Tent_Kbufss = 0.00025;
	Tent_Vmaxup = 0.006375;
	Tent_Kup = 0.00025;
	Tent_Vrel = 0.102;
	Tent_k1_ = 0.15;
	Tent_k2_ = 0.045;
	Tent_k3 = 0.060;
	Tent_k4 = 0.005;
	Tent_EC = 1.5;
	Tent_maxsr = 2.5;
	Tent_minsr = 1.;
	Tent_Vleak = 0.00036;
	Tent_Vxfer = 0.0038;
	Tent_CAPACITANCE = 0.185;
	Tent_pKNa  = 0.03;
	Tent_GbNa  = 0.00029;
	Tent_KmK   = 1.0;
	Tent_KmNa  = 40.0;
	Tent_knak  = 2.724;
	Tent_GCaL  = 0.00003980;
	Tent_GbCa  = 0.000592;
	Tent_knaca = 1000;
	Tent_KmNai = 87.5;
	Tent_KmCa  = 1.38;
	Tent_ksat  = 0.1;
	Tent_n     = 0.35;
	Tent_GpCa  = 0.1238;
	Tent_KpCa  = 0.0005;
	Tent_GpK   = 0.0146;
	Tent_RTONF = (Tent_R * Tent_T / Tent_F);
	Tent_inverseVcF2  =  (1. / (2 * Tent_Vc * Tent_F));
	Tent_inverseVcF   =  (1. / (Tent_Vc * Tent_F));
	Tent_inversevssF2 = (1. / (2 * Tent_Vss * Tent_F));

	//Assign the names of the variables used by this model
	Names_Of_Cell_Variables =
	{
		"Cai",
		"Nai",
		"Ki",
		"CaSRf",
		"CaSS",
		"sm",
		"sh",
		"sj",
		"sd",
		"sf",
		"sf2",
		"sfcass",
		"sr",
		"ss",
		"sxr1",
		"sxr2",
		"sxs",
		"y",
		"sRR"
	};

	Names_Of_Output_Data =
	{

	};

	FinalizeConstruction();
}

TNNP06RushLarsen::~TNNP06RushLarsen()
{

}






double TNNP06RushLarsen::get_initial_state_variable(const unsigned &v)
{
	// double STATE[20];
	if (cell_type == 108 || cell_type == 107 || cell_type == 106) {
		switch(v){
			case Vm_TT : return -74.7890522727;
			case Cai_TT: return 0.0001720623; // Cai
			case Nai_TT: return 8.5447311020;    // Nai
			case Ki_TT: return 136.9896086978;   // Ki
			case CaSRf_TT: return 3.2830723338;     // CaSRf
			case CaSS_TT: return 0.0006146554; // CaSS
			case sm_TT: return 0.0145766758;     // sm
			case sh_TT: return 0.2979720207;    // sh
			case sj_TT: return 0.0692509548;    // sj
			case sd_TT: return 0.0001356656;     // sd
			case sf_TT: return 0.5943228461;    // sf
			case sf2_TT: return 0.8265709174;     // sf2
			case sfcass_TT: return 0.9767040566;     // sfcass
			case sr_TT: return 0.0006830833;     // sr
			case ss_TT: return 0.9717098312;     // ss
			case sxr1_TT: return 0.4663168269;     // sxr1
			case sxr2_TT: return 0.3657472179;     // sxr2
			case sxs_TT: return 0.0486609588;     // sxs
			case y_TT: return 0.0184308075;     // y
			case sRR_TT: return 0.8199969443;    // sRR
			// case 19: return 0.0;     // sOO
		}
	}
	else{
		switch(v){
			case Vm_TT : return -86.2;
			case Cai_TT: return 0.00007; // Cai
			case Nai_TT: return 7.67;    // Nai
			case Ki_TT: return 138.3;   // Ki
			case CaSRf_TT: return 1.3;     // CaSRf
			case CaSS_TT: return 0.00007; // CaSS
			case sm_TT: return 0.0;     // sm
			case sh_TT: return 0.75;    // sh
			case sj_TT: return 0.75;    // sj
			case sd_TT: return 0.0;     // sd
			case sf_TT: return 1.0;     // sf
			case sf2_TT: return 1.0;     // sf2
			case sfcass_TT: return 1.0;     // sfcass
			case sr_TT: return 0.0;     // sr
			case ss_TT: return 1.0;     // ss
			case sxr1_TT: return 0.0;     // sxr1
			case sxr2_TT: return 0.0;     // sxr2
			case sxs_TT: return 0.0;     // sxs
			case y_TT: return 0.0;     // y
			case sRR_TT: return 1.0;     // sRR
			// case 19: return 0.0;     // sOO	
		}	
	}
}

void TNNP06RushLarsen::TakeTimestep(const double& dt, double* state)
{
	double Vm 	= state[Vm_TT];
	double Cai	= state[Cai_TT];
	double Nai	= state[Nai_TT];
	double Ki		= state[Ki_TT];
	double CaSRf	= state[CaSRf_TT];
	double CaSS	= state[CaSS_TT];
	double sm		= state[sm_TT];
	double sh		= state[sh_TT];
	double sj		= state[sj_TT];
	double sd		= state[sd_TT];
	double sf		= state[sf_TT];
	double sf2	= state[sf2_TT];
	double sfcass	= state[sfcass_TT];
	double sr		= state[sr_TT];
	double ss		= state[ss_TT];
	double sxr1	= state[sxr1_TT];
	double sxr2	= state[sxr2_TT];
	double sxs	= state[sxs_TT];
	double y		= state[y_TT];
	double sRR	= state[sRR_TT];


	const double ZIndex = 1.0;


	static double Tent_inverseVcF2  =  (1. / (2 * Tent_Vc * Tent_F));
	static double Tent_inverseVcF   =  (1. / (Tent_Vc * Tent_F));
	static double Tent_inversevssF2 = (1. / (2 * Tent_Vss * Tent_F));

	// double Ek, Ena, Eks, Eca;
	// double INa, GNa, AM, BM, TAU_M, M_INF, AH_1, BH_1, TAU_H, AH_2, BH_2, H_INF, AJ_1, BJ_1, TAU_J, AJ_2, BJ_2, J_INF;
	// double ICaL, D_INF, Ad, Bd, Cd, TAU_D, F_INF, Af, Bf, Cf, TAU_F, F2_INF, Af2, Bf2, Cf2, TAU_F2, FCaSS_INF, TAU_FCaSS;
	// double Ito, Gto, R_INF, S_INF, TAU_R, TAU_S;
	// double IKr, Gkr, Xr1_INF, axr1, bxr1, TAU_Xr1, Xr2_INF, axr2, bxr2, TAU_Xr2;
	// double IKs, Gks, Xs_INF, Axs, Bxs, TAU_Xs;
	// double IK1, GK1, Ak1, Bk1, rec_iK1;
	// double INaCa;
	// double INaK, rec_iNaK;
	// double Isus, Gsus, susa;
	// double If, ifk, ifna, Gfna, Gfk, y_inf, Ay, By, tau_y;
	// double IpCa;
	// double IpK, rec_ipK;
	// double IbNa;
	// double IbCa;
	// double kCaSR, k1, k2, dRR;
	// double Irel, Ileak, Iup, Ixfer;
	double CaCSQN, bjsr, cjsr;
	double CaSSBuf, bcss, ccss;
	double CaBuf, bc, cc;
	double dNai, dKi, dCai, dCaSR, dCaSS;

	double Ek = Tent_RTONF * (log((Tent_Ko / Ki)));
	double Ena = Tent_RTONF * (log((Tent_Nao / Nai)));
	double Eks = Tent_RTONF * (log((Tent_Ko + Tent_pKNa * Tent_Nao) / (Ki + Tent_pKNa * Nai)));
	// Eks = Tent_RTONF * (log((Tent_Ko + Tent_pKNa * Tent_Nao) / (Ki + Tent_pKNa * 7.67)));
	double Eca = 0.5 * Tent_RTONF * (log((Tent_Cao / Cai)));

	/*if (cell_type == 13 || cell_type == 14 || cell_type == 15)
	    GNa = 14.838;
	else
	    GNa = 130.5744;*/  // the original implementation of 108 and VW cells were in the other way round regarding the INa conductance..
	// change to the true values: here: haibo Sun 10 May 2015 09:44:30 BST


	// if (cell_type == 108 || cell_type == 107 || cell_type == 106)
	// 	GNa = 130.5744;
	// else
	// 	GNa = 14.838;
	double GNa = (cell_type == 108 || cell_type == 107 || cell_type == 106) ? 130.5744 : 14.838;

	double AM = 1. / (1. + exp((-60. - Vm) / 5.));
	double BM = 0.1 / (1. + exp((Vm + 35.) / 5.)) + 0.10 / (1. + exp((Vm - 50.) / 200.));
	double TAU_M = AM * BM;
	double M_INF = 1. / ((1. + exp((-56.86 - Vm) / 9.03)) * (1. + exp((-56.86 - Vm) / 9.03)));
	// if (Vm >= -40.) {
	// 	AH_1 = 0.;
	// 	BH_1 = (0.77 / (0.13 * (1. + exp(-(Vm + 10.66) / 11.1))));
	// 	TAU_H = 1.0 / (AH_1 + BH_1);
	// } else {
	// 	AH_2 = (0.057 * exp(-(Vm + 80.) / 6.8));
	// 	BH_2 = (2.7 * exp(0.079 * Vm) + (3.1e5) * exp(0.3485 * Vm));
	// 	TAU_H = 1.0 / (AH_2 + BH_2);
	// }
	double AH_1 = (Vm >= -40.)*0. + (!(Vm >= -40.))*(0.057 * exp(-(Vm + 80.) / 6.8));
	double BH_1 = (Vm >= -40.)*(0.77 / (0.13 * (1. + exp(-(Vm + 10.66) / 11.1)))) + (!(Vm >= -40.))*(2.7 * exp(0.079 * Vm) + (3.1e5) * exp(0.3485 * Vm));

	double TAU_H = 1.0 / (AH_1 + BH_1);


	double H_INF = 1. / ((1. + exp((Vm + 71.55) / 7.43)) * (1. + exp((Vm + 71.55) / 7.43)));
	// if (Vm >= -40.) {
	// 	AJ_1 = 0.;
	// 	BJ_1 = (0.6 * exp((0.057) * Vm) / (1. + exp(-0.1 * (Vm + 32.))));
	// 	TAU_J = 1.0 / (AJ_1 + BJ_1);
	// } else {
	// 	AJ_2 = (((-2.5428e4) * exp(0.2444 * Vm) - (6.948e-6) * exp(-0.04391 * Vm)) * (Vm + 37.78) / (1. + exp(0.311 * (Vm + 79.23))));
	// 	BJ_2 = (0.02424 * exp(-0.01052 * Vm) / (1. + exp(-0.1378 * (Vm + 40.14))));
	// 	TAU_J = 1.0 / (AJ_2 + BJ_2);
	// }
	double AJ_1 = (Vm >= -40.) ? 0. : (((-2.5428e4) * exp(0.2444 * Vm) - (6.948e-6) * exp(-0.04391 * Vm)) * (Vm + 37.78) / (1. + exp(0.311 * (Vm + 79.23))));
	double BJ_1 = (Vm >= -40.) ? (0.6 * exp((0.057) * Vm) / (1. + exp(-0.1 * (Vm + 32.)))) : (0.02424 * exp(-0.01052 * Vm) / (1. + exp(-0.1378 * (Vm + 40.14))));
	double TAU_J = 1.0 / (AJ_1 + BJ_1);


	double J_INF = H_INF;
	// sm = M_INF - (M_INF - sm) * exp(-dt / TAU_M);
	// sh = H_INF - (H_INF - sh) * exp(-dt / TAU_H);
	// sj = J_INF - (J_INF - sj) * exp(-dt / TAU_J);
	double INa = GNa * sm * sm * sm * sh * sj * (Vm - Ena);

	double D_INF = 1. / (1. + exp((-8 - Vm) / 7.5));
	double Ad = 1.4 / (1. + exp((-35 - Vm) / 13)) + 0.25;
	double Bd = 1.4 / (1. + exp((Vm + 5) / 5));
	double Cd = 1. / (1. + exp((50 - Vm) / 20));
	double TAU_D = Ad * Bd + Cd;
	double F_INF = 1. / (1. + exp((Vm + 20) / 7));
	double Af = 1102.5 * exp(-(Vm + 27) * (Vm + 27) / 225);
	double Bf = 200. / (1 + exp((13 - Vm) / 10.));
	double Cf = (180. / (1 + exp((Vm + 30) / 10))) + 20;
	double TAU_F = Af + Bf + Cf;
	double F2_INF = 0.67 / (1. + exp((Vm + 35) / 7)) + 0.33;
	double Af2 = 562 * exp(-(Vm + 27) * (Vm + 27) / 240);
	double Bf2 = 31 / (1. + exp((25 - Vm) / 10));
	double Cf2 = 80 / (1. + exp((Vm + 30) / 10));
	double TAU_F2 = Af2 + Bf2 + Cf2;
	double FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
	double TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;
	// sd = D_INF - (D_INF - sd) * exp(-dt / TAU_D);
	// sf = F_INF - (F_INF - sf) * exp(-dt / TAU_F);
	// sf2 = F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2);
	// sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS);
	double ICaL = Tent_GCaL * sd * sf * sf2 * sfcass * 4 * (Vm - 15) * (Tent_F * Tent_F / (Tent_R * Tent_T)) * (0.25 * exp(2 * (Vm - 15) * Tent_F / (Tent_R * Tent_T)) * CaSS - Tent_Cao) / (exp(2 * (Vm - 15) * Tent_F / (Tent_R * Tent_T)) - 1.);

	double Gto;
	double R_INF, S_INF, TAU_R, TAU_S;
	if (cell_type == 108 || cell_type == 107 || cell_type == 106) {
		Gto = 0.08184;
		R_INF = 1. / (1. + exp((20 - Vm) / 13.));
		S_INF = 1. / (1. + exp((Vm + 27) / 13.));
		TAU_R = 10.45 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 7.3;
		TAU_S = 85.*exp(-(Vm + 25.) * (Vm + 25.) / 320.) + 5. / (1. + exp((Vm - 40.) / 5.)) + 42.;




		/* formulation from StewardP, human Purkinje fibre cells */
		R_INF = 1. / (1. + exp((20 - Vm) / 13.));
		S_INF = 1. / (1. + exp((Vm + 27) / 13.));
		TAU_R = 1.1 * 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8 + 6.5;
		TAU_S = 85.*exp(-(Vm + 45. - 20.) * (Vm + 45. - 20.) / 320.) + 5. / (1. + exp((Vm - 20. - 20.) / 5.)) + 3. + 39.;

	} else if (cell_type == 102) {
		Gto = 0.073;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 28) / 5.));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 1000.*exp(-(Vm + 67) * (Vm + 67) / 1000.) + 8.;
	} else if (cell_type == 101 || cell_type == 100) {
		Gto = 0.294;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20) / 5.));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.0 * exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (cell_type == 105) {
		Gto = 0.073 * 4;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 28) / 5.));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 1000.*exp(-(Vm + 67) * (Vm + 67) / 1000.) + 8.;
	} else {
		Gto = 0.294 * 0.8;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20) / 5.));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.*exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	}
	// ss = S_INF - (S_INF - ss) * exp(-dt / TAU_S);
	// sr = R_INF - (R_INF - sr) * exp(-dt / TAU_R);

	double Ito = Gto * sr * ss * (Vm - Ek);

	//APEX BASE HETERO
	// first remove heterogeneity in Ito, by haibo time: Feb 27 2015
	/*if (cell_type == 102 || cell_type == 101 || cell_type == 100) {
	    Ito = Ito * ((-1.0 * (ZIndex - 10) / 318) + 2.0);
	}*/

	// if (cell_type == 108 || cell_type == 107 || cell_type == 106)
	// 	Gkr = 0.0918;
	// else
	// 	Gkr = 0.153;
	double Gkr = (cell_type == 108 || cell_type == 107 || cell_type == 106) ? 0.0918 : 0.153;

	double Xr1_INF = 1. / (1. + exp((-26. - Vm) / 7.));
	double axr1 = 450. / (1. + exp((-45. - Vm) / 10.));
	double bxr1 = 6. / (1. + exp((Vm - (-30.)) / 11.5));
	double TAU_Xr1 = axr1 * bxr1;
	double Xr2_INF = 1. / (1. + exp((Vm - (-88.)) / 24.));
	double axr2 = 3. / (1. + exp((-60. - Vm) / 20.));
	double bxr2 = 1.12 / (1. + exp((Vm - 60.) / 20.));
	double TAU_Xr2 = axr2 * bxr2;
	// sxr1 = Xr1_INF - (Xr1_INF - sxr1) * exp(-dt / TAU_Xr1);
	// sxr2 = Xr2_INF - (Xr2_INF - sxr2) * exp(-dt / TAU_Xr2);


	double IKr = Gkr * sqrt(Tent_Ko / 5.4) * sxr1 * sxr2 * (Vm - Ek);

	double Gks;
	if (cell_type == 108 || cell_type == 107 || cell_type == 106)
		Gks = 0.2352;
	else if (cell_type == 101)
		Gks = 0.098;
	else if (cell_type == 102 || cell_type == 100)
		Gks = 0.392;
	else if (cell_type == 104)
		Gks = 0.098 * 2;
	else
		Gks = 0.392 * 2;
	double Xs_INF = 1. / (1. + exp((-5. - Vm) / 14.));
	double Axs = (1400. / (sqrt(1. + exp((5. - Vm) / 6))));
	double Bxs = (1. / (1. + exp((Vm - 35.) / 15.)));
	double TAU_Xs = Axs * Bxs + 80;
	// sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);

	double IKs = Gks * sxs * sxs * (Vm - Eks);

	//APEX BASE HETERO
	// if (cell_type == 102 || cell_type == 101 || cell_type == 100) {
	// 	// IKs = IKs * ((-1.0 * (ZIndex - 10) / 318) + 2.0);  // John's version
	// 	IKs = IKs * ZIndex;   //  by haibo
	// }

	IKs  = (cell_type == 102 || cell_type == 101 || cell_type == 100) ? IKs * ZIndex : IKs;

	double GK1, Ak1, Bk1, rec_iK1, IK1;
	if (cell_type != 108 && cell_type != 107 && cell_type != 106) {
		GK1 = 5.405;
		Ak1 = 0.1 / (1. + exp(0.06 * (Vm - Ek - 200)));
		Bk1 = (3.*exp(0.0002 * (Vm - Ek + 100)) + exp(0.1 * (Vm - Ek - 10))) / (1. + exp(-0.5 * (Vm - Ek)));
		rec_iK1 = Ak1 / (Ak1 + Bk1);
		IK1 = GK1 * rec_iK1 * (Vm - Ek);
	} else {
		GK1 = 0.065;
		rec_iK1 = 1. / (1. + exp(0.1 * (Vm + 75.44)));
		IK1 = GK1 * rec_iK1 * (Vm - Ek - 8.0);
	}

	double Gsus, susa, Isus;
	if (cell_type == 108 || cell_type == 107 || cell_type == 106) {
		Gsus = 0.0227;
		susa = 1. / (1. + exp((5 - Vm) / 17.));
		Isus = Gsus * susa * (Vm - Ek);
	} else
		Isus = 0.0;

	double Gfna, Gfk, y_inf, Ay, By, tau_y, ifna, ifk, If;
	if (cell_type == 108 || cell_type == 107 || cell_type == 106) {
		Gfna = 0.0145654;
		Gfk = 0.0234346;
		y_inf = 1.0 / (1.0 + exp((Vm + 80.6) / 6.8));
		Ay = exp(-2.9 - (0.04 * Vm));
		By = exp(3.6 + (0.11 * Vm));
		tau_y = 4000. / (Ay + By);
		// y = y_inf - (y_inf - y) * exp(-dt / tau_y);
		ifna = Gfna * y * (Vm - Ena);
		ifk = Gfk * y * (Vm - Ek);
		If = ifk + ifna;
	} else {
		//this was not in the original, but to minimise changes we must set these parameters such that the derivatives
		//	of these variables are zero rather than nans
		Gfna = 0.0145654;
		Gfk = 0.0234346;
		y_inf = y;
		Ay = 1.0;
		By = 1.0;
		tau_y = 1.0;


		ifna = 0.0;
		ifk = 0.0;
		If = 0.0;
	}

	double INaCa = Tent_knaca * (1. / (Tent_KmNai * Tent_KmNai * Tent_KmNai + Tent_Nao * Tent_Nao * Tent_Nao)) * (1. / (Tent_KmCa + Tent_Cao)) * (1. / (1 + Tent_ksat * exp((Tent_n - 1) * Vm * Tent_F / (Tent_R * Tent_T)))) * (exp(Tent_n * Vm * Tent_F / (Tent_R * Tent_T)) * Nai * Nai * Nai * Tent_Cao - exp((Tent_n - 1) * Vm * Tent_F / (Tent_R * Tent_T)) * Tent_Nao * Tent_Nao * Tent_Nao * Cai * 2.5);

	double rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * Vm * Tent_F / (Tent_R * Tent_T)) + 0.0353 * exp(-Vm * Tent_F / (Tent_R * Tent_T))));
	double INaK = Tent_knak * (Tent_Ko / (Tent_Ko + Tent_KmK)) * (Nai / (Nai + Tent_KmNa)) * rec_iNaK;

	double IpCa = Tent_GpCa * Cai / (Tent_KpCa + Cai);

	double rec_ipK = 1. / (1. + exp((25 - Vm) / 5.98));
	double IpK = Tent_GpK * rec_ipK * (Vm - Ek);

	double IbNa = Tent_GbNa * (Vm - Ena);

	double IbCa = Tent_GbCa * (Vm - Eca);

	double kCaSR = Tent_maxsr - ((Tent_maxsr - Tent_minsr) / (1 + (Tent_EC / CaSRf) * (Tent_EC / CaSRf)));
	double k1 = Tent_k1_ / kCaSR;
	double k2 = Tent_k2_ * kCaSR;
	double dRR = Tent_k4 * (1 - sRR) - k2 * CaSS * sRR;
	// sOO = k1 * CaSS * CaSS * sRR / (Tent_k3 + k1 * CaSS * CaSS);
	const double sOO = k1 * CaSS * CaSS * sRR / (Tent_k3 + k1 * CaSS * CaSS);

	double Irel = Tent_Vrel * sOO * (CaSRf - CaSS);
	double Ileak = Tent_Vleak * (CaSRf - Cai);
	double Iup = Tent_Vmaxup / (1. + ((Tent_Kup * Tent_Kup) / (Cai * Cai)));
	double Ixfer = Tent_Vxfer * (CaSS - Cai);

	CaCSQN = Tent_Bufsr * CaSRf / (CaSRf + Tent_Kbufsr);
	dCaSR = dt * (Iup - Irel - Ileak);
	bjsr = Tent_Bufsr - CaCSQN - dCaSR - CaSRf + Tent_Kbufsr;
	cjsr = Tent_Kbufsr * (CaCSQN + dCaSR + CaSRf);
	CaSRf = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;
	// double dCaSR = (Iup - Irel - Ileak)/(1.0  + (Tent_Kbufsr*Tent_Bufsr)/pow(CaSRf + Tent_Kbufsr, 2.0));

	CaSSBuf = Tent_Bufss * CaSS / (CaSS + Tent_Kbufss);
	dCaSS = dt * (-Ixfer * (Tent_Vc / Tent_Vss) + Irel * (Tent_Vsr / Tent_Vss) + (-ICaL * Tent_inversevssF2 * Tent_CAPACITANCE));
	bcss = Tent_Bufss - CaSSBuf - dCaSS - CaSS + Tent_Kbufss;
	ccss = Tent_Kbufss * (CaSSBuf + dCaSS + CaSS);
	CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;
	// double dCaSS =(-Ixfer * (Tent_Vc / Tent_Vss) + Irel * (Tent_Vsr / Tent_Vss) + (-ICaL * Tent_inversevssF2 * Tent_CAPACITANCE))/(1.0 + (Tent_Bufss*Tent_Kbufss)/pow(CaSS+Tent_Kbufss, 2.0));


	CaBuf = Tent_Bufc * Cai / (Cai + Tent_Kbufc);
	dCai = dt * ((-(IbCa + IpCa - 2 * INaCa) * Tent_inverseVcF2 * Tent_CAPACITANCE) - (Iup - Ileak) * (Tent_Vsr / Tent_Vc) + Ixfer);
	bc = Tent_Bufc - CaBuf - dCai - Cai + Tent_Kbufc;
	cc = Tent_Kbufc * (CaBuf + dCai + Cai);
	Cai = (sqrt(bc * bc + 4 * cc) - bc) / 2;
	// double dCai = ((-(IbCa + IpCa - 2 * INaCa) * Tent_inverseVcF2 * Tent_CAPACITANCE) - (Iup - Ileak) * (Tent_Vsr / Tent_Vc) + Ixfer)/(1.0 + (Tent_Bufc*Tent_Kbufc)/pow(Cai + Tent_Kbufc, 2.0));
	

	dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa + ifna) * Tent_inverseVcF * Tent_CAPACITANCE;

	dKi = -(IK1 + Ito + Isus + IKr + IKs - 2 * INaK + IpK + ifk) * Tent_inverseVcF * Tent_CAPACITANCE;

	state[Vm_TT] 			+= -dt*(IKr + IKs + IK1 + Ito + Isus + (INa) + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + (If) + get_stimulus(this->time()));
	state[Cai_TT] 			+= dt*dCai;
	state[Nai_TT]			+= dt*dNai;
	state[Ki_TT]			+= dt*dKi;
	state[CaSRf_TT]			= CaSRf;
	state[CaSS_TT]			= CaSS;

	state[sm_TT]			= M_INF+(sm-M_INF)*exp(-dt/TAU_M);
	state[sh_TT]			= H_INF+(sh-H_INF)*exp(-dt/TAU_H);
	state[sj_TT]			= J_INF+(sj-J_INF)*exp(-dt/TAU_J);
	state[sd_TT]			= D_INF+(sd-D_INF)*exp(-dt/TAU_D);
	state[sf_TT]			= F_INF+(sf-F_INF)*exp(-dt/TAU_F);
	state[sf2_TT]			= F2_INF+(sf2-F2_INF)*exp(-dt/TAU_F2);
	state[sfcass_TT]		= FCaSS_INF+(sfcass-FCaSS_INF)*exp(-dt/TAU_FCaSS);
	state[sr_TT]			= S_INF+(sr-S_INF)*exp(-dt/TAU_S);
	state[ss_TT]			= S_INF+(ss-S_INF)*exp(-dt/TAU_S);
	state[sxr1_TT]			= Xr1_INF+(sxr1-Xr1_INF)*exp(-dt/TAU_Xr1);
	state[sxr2_TT]			= Xr2_INF+(sxr2-Xr2_INF)*exp(-dt/TAU_Xr2);
	state[sxs_TT]			= Xs_INF+(sxs-Xs_INF)*exp(-dt/TAU_Xs);
	state[y_TT]				= y_inf+(y-y_inf)*exp(-dt/tau_y);
	state[sRR_TT]			+= dt*dRR;

}


void TNNP06RushLarsen::get_output(double *state, double *out)
{
	//Intentionally empty
}

}; //End namespace