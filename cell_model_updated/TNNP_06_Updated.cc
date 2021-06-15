#include "TNNP_06_Updated.h"

namespace oomph{
	
ExplicitTNNP06VentUpdated::ExplicitTNNP06VentUpdated(){ //: CellModelBaseUpdated(){
	//Assign the values of the physical constants
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
		"sRR"//,
		// "sOO"
	};
	Names_Of_Other_Parameters =
	{
		"ZIndex"
	};
	Names_Of_Other_Variables =
	{

	};
	Names_Of_Output_Data =
	{
		"ActiveStrain"
		// "EKs",
		// "IKs"
	};

	FinalizeConstruction();
}



double ExplicitTNNP06VentUpdated::return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
{
	double STATE[20];
	if (cell_type == 108 || cell_type == 107 || cell_type == 106) {
		STATE[0]  = 0.0001720623; // Cai
		STATE[1]  = 8.5447311020;    // Nai
		STATE[2]  = 136.9896086978;   // Ki
		STATE[3]  = 3.2830723338;     // CaSRf
		STATE[4]  = 0.0006146554; // CaSS
		STATE[5]  = 0.0145766758;     // sm
		STATE[6]  = 0.2979720207;    // sh
		STATE[7]  = 0.0692509548;    // sj
		STATE[8]  = 0.0001356656;     // sd
		STATE[9]  = 0.5943228461;    // sf
		STATE[10] = 0.8265709174;     // sf2
		STATE[11] = 0.9767040566;     // sfcass
		STATE[12] = 0.0006830833;     // sr
		STATE[13] = 0.9717098312;     // ss
		STATE[14] = 0.4663168269;     // sxr1
		STATE[15] = 0.3657472179;     // sxr2
		STATE[16] = 0.0486609588;     // sxs
		STATE[17] = 0.0184308075;     // y
		STATE[18] = 0.8199969443;    // sRR
		STATE[19] = 0.0;     // sOO
	}
	else{
		STATE[0]  = 0.00007; // Cai
		STATE[1]  = 7.67;    // Nai
		STATE[2]  = 138.3;   // Ki
		STATE[3]  = 1.3;     // CaSRf
		STATE[4]  = 0.00007; // CaSS
		STATE[5]  = 0.0;     // sm
		STATE[6]  = 0.75;    // sh
		STATE[7]  = 0.75;    // sj
		STATE[8]  = 0.0;     // sd
		STATE[9]  = 1.0;     // sf
		STATE[10] = 1.0;     // sf2
		STATE[11] = 1.0;     // sfcass
		STATE[12] = 0.0;     // sr
		STATE[13] = 1.0;     // ss
		STATE[14] = 0.0;     // sxr1
		STATE[15] = 0.0;     // sxr2
		STATE[16] = 0.0;     // sxs
		STATE[17] = 0.0;     // y
		STATE[18] = 1.0;     // sRR
		STATE[19] = 0.0;     // sOO
	}

	if(v>=0 && v < 20){
		return STATE[v];
	}
}

double ExplicitTNNP06VentUpdated::return_initial_membrane_potential(const unsigned &cell_type)
{
	if (cell_type == 108 || cell_type == 107 || cell_type == 106){
		return -74.7890522727;
	}
	else{
		return -86.2;
	}
}

void ExplicitTNNP06VentUpdated::Calculate_Derivatives(const double &Vm,
											const Vector<double> &CellVariables,
											const double &t,
											const unsigned &cell_type,
											const double &Istim,
											const Vector<double> &Other_Parameters,
											const Vector<double> &Other_Variables,
											
											Vector<double> &Variable_Derivatives,
											double &Iion)
{
	//Get all the cell variables
	const double Cai	= CellVariables[Cai_TT];
	const double Nai	= CellVariables[Nai_TT];
	const double Ki		= CellVariables[Ki_TT];
	const double CaSRf	= CellVariables[CaSRf_TT];
	const double CaSS	= CellVariables[CaSS_TT];
	const double sm		= CellVariables[sm_TT];
	const double sh		= CellVariables[sh_TT];
	const double sj		= CellVariables[sj_TT];
	const double sd		= CellVariables[sd_TT];
	const double sf		= CellVariables[sf_TT];
	const double sf2	= CellVariables[sf2_TT];
	const double sfcass	= CellVariables[sfcass_TT];
	const double sr		= CellVariables[sr_TT];
	const double ss		= CellVariables[ss_TT];
	const double sxr1	= CellVariables[sxr1_TT];
	const double sxr2	= CellVariables[sxr2_TT];
	const double sxs	= CellVariables[sxs_TT];
	const double y		= CellVariables[y_TT];
	const double sRR	= CellVariables[sRR_TT];
	// const double sOO	= CellVariables[sOO_TT];

	// std::cout << "Cai " << Cai << std::endl;
	// std::cout << "Nai " << Nai << std::endl;
	// std::cout << "Ki " << Ki << std::endl;
	// std::cout << "CaSRf " << CaSRf << std::endl;
	// std::cout << "CaSS " << CaSS << std::endl;
	// std::cout << "sm " << sm << std::endl;
	// std::cout << "sh " << sh << std::endl;
	// std::cout << "sj " << sj << std::endl;
	// std::cout << "sd " << sd << std::endl;
	// std::cout << "sf " << sf << std::endl;
	// std::cout << "sf2 " << sf2 << std::endl;
	// std::cout << "sfcass " << sfcass << std::endl;
	// std::cout << "sr " << sr << std::endl;
	// std::cout << "ss " << ss << std::endl;
	// std::cout << "sxr1 " << sxr1 << std::endl;
	// std::cout << "sxr2 " << sxr2 << std::endl;
	// std::cout << "sxs " << sxs << std::endl;
	// std::cout << "y " << y << std::endl;
	// std::cout << "sRR " << sRR << std::endl;
	// std::cout << "sOO " << sOO << std::endl;

	const double ZIndex = Other_Parameters[ZIndex_TT];


	// static double Tent_RTONF = (Tent_R * Tent_T / Tent_F);
	static double Tent_inverseVcF2  =  (1. / (2 * Tent_Vc * Tent_F));
	static double Tent_inverseVcF   =  (1. / (Tent_Vc * Tent_F));
	static double Tent_inversevssF2 = (1. / (2 * Tent_Vss * Tent_F));

	double Ek, Ena, Eks, Eca;
	double INa, GNa, AM, BM, TAU_M, M_INF, AH_1, BH_1, TAU_H, AH_2, BH_2, H_INF, AJ_1, BJ_1, TAU_J, AJ_2, BJ_2, J_INF;
	double ICaL, D_INF, Ad, Bd, Cd, TAU_D, F_INF, Af, Bf, Cf, TAU_F, F2_INF, Af2, Bf2, Cf2, TAU_F2, FCaSS_INF, TAU_FCaSS;
	double Ito, Gto, R_INF, S_INF, TAU_R, TAU_S;
	double IKr, Gkr, Xr1_INF, axr1, bxr1, TAU_Xr1, Xr2_INF, axr2, bxr2, TAU_Xr2;
	double IKs, Gks, Xs_INF, Axs, Bxs, TAU_Xs;
	double IK1, GK1, Ak1, Bk1, rec_iK1;
	double INaCa;
	double INaK, rec_iNaK;
	double Isus, Gsus, susa;
	double If, ifk, ifna, Gfna, Gfk, y_inf, Ay, By, tau_y;
	double IpCa;
	double IpK, rec_ipK;
	double IbNa;
	double IbCa;
	double kCaSR, k1, k2, dRR;
	double Irel, Ileak, Iup, Ixfer;
	double CaCSQN, bjsr, cjsr;
	double CaSSBuf, bcss, ccss;
	double CaBuf, bc, cc;
	double dNai, dKi, dCai, dCaSR, dCaSS;

	Ek = Tent_RTONF * (log((Tent_Ko / Ki)));
	Ena = Tent_RTONF * (log((Tent_Nao / Nai)));
	Eks = Tent_RTONF * (log((Tent_Ko + Tent_pKNa * Tent_Nao) / (Ki + Tent_pKNa * Nai)));
	// Eks = Tent_RTONF * (log((Tent_Ko + Tent_pKNa * Tent_Nao) / (Ki + Tent_pKNa * 7.67)));
	Eca = 0.5 * Tent_RTONF * (log((Tent_Cao / Cai)));

	/*if (cell_type == 13 || cell_type == 14 || cell_type == 15)
	    GNa = 14.838;
	else
	    GNa = 130.5744;*/  // the original implementation of 108 and VW cells were in the other way round regarding the INa conductance..
	// change to the true values: here: haibo Sun 10 May 2015 09:44:30 BST


	if (cell_type == 108 || cell_type == 107 || cell_type == 106)
		GNa = 130.5744;
	else
		GNa = 14.838;
	AM = 1. / (1. + exp((-60. - Vm) / 5.));
	BM = 0.1 / (1. + exp((Vm + 35.) / 5.)) + 0.10 / (1. + exp((Vm - 50.) / 200.));
	TAU_M = AM * BM;
	M_INF = 1. / ((1. + exp((-56.86 - Vm) / 9.03)) * (1. + exp((-56.86 - Vm) / 9.03)));
	if (Vm >= -40.) {
		AH_1 = 0.;
		BH_1 = (0.77 / (0.13 * (1. + exp(-(Vm + 10.66) / 11.1))));
		TAU_H = 1.0 / (AH_1 + BH_1);
	} else {
		AH_2 = (0.057 * exp(-(Vm + 80.) / 6.8));
		BH_2 = (2.7 * exp(0.079 * Vm) + (3.1e5) * exp(0.3485 * Vm));
		TAU_H = 1.0 / (AH_2 + BH_2);
	}
	H_INF = 1. / ((1. + exp((Vm + 71.55) / 7.43)) * (1. + exp((Vm + 71.55) / 7.43)));
	if (Vm >= -40.) {
		AJ_1 = 0.;
		BJ_1 = (0.6 * exp((0.057) * Vm) / (1. + exp(-0.1 * (Vm + 32.))));
		TAU_J = 1.0 / (AJ_1 + BJ_1);
	} else {
		AJ_2 = (((-2.5428e4) * exp(0.2444 * Vm) - (6.948e-6) * exp(-0.04391 * Vm)) * (Vm + 37.78) / (1. + exp(0.311 * (Vm + 79.23))));
		BJ_2 = (0.02424 * exp(-0.01052 * Vm) / (1. + exp(-0.1378 * (Vm + 40.14))));
		TAU_J = 1.0 / (AJ_2 + BJ_2);
	}
	J_INF = H_INF;
	// sm = M_INF - (M_INF - sm) * exp(-dt / TAU_M);
	// sh = H_INF - (H_INF - sh) * exp(-dt / TAU_H);
	// sj = J_INF - (J_INF - sj) * exp(-dt / TAU_J);
	INa = GNa * sm * sm * sm * sh * sj * (Vm - Ena);

	D_INF = 1. / (1. + exp((-8 - Vm) / 7.5));
	Ad = 1.4 / (1. + exp((-35 - Vm) / 13)) + 0.25;
	Bd = 1.4 / (1. + exp((Vm + 5) / 5));
	Cd = 1. / (1. + exp((50 - Vm) / 20));
	TAU_D = Ad * Bd + Cd;
	F_INF = 1. / (1. + exp((Vm + 20) / 7));
	Af = 1102.5 * exp(-(Vm + 27) * (Vm + 27) / 225);
	Bf = 200. / (1 + exp((13 - Vm) / 10.));
	Cf = (180. / (1 + exp((Vm + 30) / 10))) + 20;
	TAU_F = Af + Bf + Cf;
	F2_INF = 0.67 / (1. + exp((Vm + 35) / 7)) + 0.33;
	Af2 = 562 * exp(-(Vm + 27) * (Vm + 27) / 240);
	Bf2 = 31 / (1. + exp((25 - Vm) / 10));
	Cf2 = 80 / (1. + exp((Vm + 30) / 10));
	TAU_F2 = Af2 + Bf2 + Cf2;
	FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
	TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;
	// sd = D_INF - (D_INF - sd) * exp(-dt / TAU_D);
	// sf = F_INF - (F_INF - sf) * exp(-dt / TAU_F);
	// sf2 = F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2);
	// sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS);
	ICaL = Tent_GCaL * sd * sf * sf2 * sfcass * 4 * (Vm - 15) * (Tent_F * Tent_F / (Tent_R * Tent_T)) * (0.25 * exp(2 * (Vm - 15) * Tent_F / (Tent_R * Tent_T)) * CaSS - Tent_Cao) / (exp(2 * (Vm - 15) * Tent_F / (Tent_R * Tent_T)) - 1.);

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

	Ito = Gto * sr * ss * (Vm - Ek);

	//APEX BASE HETERO
	// first remove heterogeneity in Ito, by haibo time: Feb 27 2015
	/*if (cell_type == 102 || cell_type == 101 || cell_type == 100) {
	    Ito = Ito * ((-1.0 * (ZIndex - 10) / 318) + 2.0);
	}*/

	if (cell_type == 108 || cell_type == 107 || cell_type == 106)
		Gkr = 0.0918;
	else
		Gkr = 0.153;
	Xr1_INF = 1. / (1. + exp((-26. - Vm) / 7.));
	axr1 = 450. / (1. + exp((-45. - Vm) / 10.));
	bxr1 = 6. / (1. + exp((Vm - (-30.)) / 11.5));
	TAU_Xr1 = axr1 * bxr1;
	Xr2_INF = 1. / (1. + exp((Vm - (-88.)) / 24.));
	axr2 = 3. / (1. + exp((-60. - Vm) / 20.));
	bxr2 = 1.12 / (1. + exp((Vm - 60.) / 20.));
	TAU_Xr2 = axr2 * bxr2;
	// sxr1 = Xr1_INF - (Xr1_INF - sxr1) * exp(-dt / TAU_Xr1);
	// sxr2 = Xr2_INF - (Xr2_INF - sxr2) * exp(-dt / TAU_Xr2);


	IKr = Gkr * sqrt(Tent_Ko / 5.4) * sxr1 * sxr2 * (Vm - Ek);

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
	Xs_INF = 1. / (1. + exp((-5. - Vm) / 14.));
	Axs = (1400. / (sqrt(1. + exp((5. - Vm) / 6))));
	Bxs = (1. / (1. + exp((Vm - 35.) / 15.)));
	TAU_Xs = Axs * Bxs + 80;
	// sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);

	IKs = Gks * sxs * sxs * (Vm - Eks);

	//APEX BASE HETERO
	if (cell_type == 102 || cell_type == 101 || cell_type == 100) {
		// IKs = IKs * ((-1.0 * (ZIndex - 10) / 318) + 2.0);  // John's version
		IKs = IKs * ZIndex;   //  by haibo
	}

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

	if (cell_type == 108 || cell_type == 107 || cell_type == 106) {
		Gsus = 0.0227;
		susa = 1. / (1. + exp((5 - Vm) / 17.));
		Isus = Gsus * susa * (Vm - Ek);
	} else
		Isus = 0.0;

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

	INaCa = Tent_knaca * (1. / (Tent_KmNai * Tent_KmNai * Tent_KmNai + Tent_Nao * Tent_Nao * Tent_Nao)) * (1. / (Tent_KmCa + Tent_Cao)) * (1. / (1 + Tent_ksat * exp((Tent_n - 1) * Vm * Tent_F / (Tent_R * Tent_T)))) * (exp(Tent_n * Vm * Tent_F / (Tent_R * Tent_T)) * Nai * Nai * Nai * Tent_Cao - exp((Tent_n - 1) * Vm * Tent_F / (Tent_R * Tent_T)) * Tent_Nao * Tent_Nao * Tent_Nao * Cai * 2.5);

	rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * Vm * Tent_F / (Tent_R * Tent_T)) + 0.0353 * exp(-Vm * Tent_F / (Tent_R * Tent_T))));
	INaK = Tent_knak * (Tent_Ko / (Tent_Ko + Tent_KmK)) * (Nai / (Nai + Tent_KmNa)) * rec_iNaK;

	IpCa = Tent_GpCa * Cai / (Tent_KpCa + Cai);

	rec_ipK = 1. / (1. + exp((25 - Vm) / 5.98));
	IpK = Tent_GpK * rec_ipK * (Vm - Ek);

	IbNa = Tent_GbNa * (Vm - Ena);

	IbCa = Tent_GbCa * (Vm - Eca);

	kCaSR = Tent_maxsr - ((Tent_maxsr - Tent_minsr) / (1 + (Tent_EC / CaSRf) * (Tent_EC / CaSRf)));
	k1 = Tent_k1_ / kCaSR;
	k2 = Tent_k2_ * kCaSR;
	dRR = Tent_k4 * (1 - sRR) - k2 * CaSS * sRR;
	// sRR += dt * dRR;
	// sOO = k1 * CaSS * CaSS * sRR / (Tent_k3 + k1 * CaSS * CaSS);
	const double sOO = k1 * CaSS * CaSS * sRR / (Tent_k3 + k1 * CaSS * CaSS);

	Irel = Tent_Vrel * sOO * (CaSRf - CaSS);
	Ileak = Tent_Vleak * (CaSRf - Cai);
	Iup = Tent_Vmaxup / (1. + ((Tent_Kup * Tent_Kup) / (Cai * Cai)));
	Ixfer = Tent_Vxfer * (CaSS - Cai);

	// CaCSQN = Tent_Bufsr * CaSRf / (CaSRf + Tent_Kbufsr);
	// dCaSR = dt * (Iup - Irel - Ileak);
	// bjsr = Tent_Bufsr - CaCSQN - dCaSR - CaSRf + Tent_Kbufsr;
	// cjsr = Tent_Kbufsr * (CaCSQN + dCaSR + CaSRf);
	// // CaSRf = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;
	dCaSR = (Iup - Irel - Ileak)/(1.0  + (Tent_Kbufsr*Tent_Bufsr)/pow(CaSRf + Tent_Kbufsr, 2.0));

	// CaSSBuf = Tent_Bufss * CaSS / (CaSS + Tent_Kbufss);
	// dCaSS = dt * (-Ixfer * (Tent_Vc / Tent_Vss) + Irel * (Tent_Vsr / Tent_Vss) + (-ICaL * Tent_inversevssF2 * Tent_CAPACITANCE));
	// bcss = Tent_Bufss - CaSSBuf - dCaSS - CaSS + Tent_Kbufss;
	// ccss = Tent_Kbufss * (CaSSBuf + dCaSS + CaSS);
	// // CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;
	dCaSS =(-Ixfer * (Tent_Vc / Tent_Vss) + Irel * (Tent_Vsr / Tent_Vss) + (-ICaL * Tent_inversevssF2 * Tent_CAPACITANCE))/(1.0 + (Tent_Bufss*Tent_Kbufss)/pow(CaSS+Tent_Kbufss, 2.0));


	// CaBuf = Tent_Bufc * Cai / (Cai + Tent_Kbufc);
	// dCai = dt * ((-(IbCa + IpCa - 2 * INaCa) * Tent_inverseVcF2 * Tent_CAPACITANCE) - (Iup - Ileak) * (Tent_Vsr / Tent_Vc) + Ixfer);
	// bc = Tent_Bufc - CaBuf - dCai - Cai + Tent_Kbufc;
	// cc = Tent_Kbufc * (CaBuf + dCai + Cai);
	// // Cai = (sqrt(bc * bc + 4 * cc) - bc) / 2;
	dCai = ((-(IbCa + IpCa - 2 * INaCa) * Tent_inverseVcF2 * Tent_CAPACITANCE) - (Iup - Ileak) * (Tent_Vsr / Tent_Vc) + Ixfer)/(1.0 + (Tent_Bufc*Tent_Kbufc)/pow(Cai + Tent_Kbufc, 2.0));
	

	dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa + ifna) * Tent_inverseVcF * Tent_CAPACITANCE;
	// Nai += dt * dNai;

	dKi = -(IK1 + Ito + Isus + IKr + IKs - 2 * INaK + IpK + ifk) * Tent_inverseVcF * Tent_CAPACITANCE;
	// Ki += dt * dKi;


	//Set the values of the variable derivatives
	// Variable_Derivatives[Cai_TT] 	= ((sqrt(bc * bc + 4 * cc) - bc) / 2 - Cai)/dt;
	// Variable_Derivatives[Nai_TT]		= dNai;
	// Variable_Derivatives[Ki_TT]		= dKi;
	// Variable_Derivatives[CaSRf_TT]	= ((sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2 - CaSRf)/dt;
	// Variable_Derivatives[CaSS_TT]	= ((sqrt(bcss * bcss + 4 * ccss) - bcss) / 2 - CaSS)/dt;
	// Variable_Derivatives[sm_TT]		= (M_INF - (M_INF - sm) * exp(-dt / TAU_M) - sm)/dt;
	// Variable_Derivatives[sh_TT]		= (H_INF - (H_INF - sh) * exp(-dt / TAU_H) - sh)/dt;
	// Variable_Derivatives[sj_TT]		= (J_INF - (J_INF - sj) * exp(-dt / TAU_J) - sj)/dt;
	// Variable_Derivatives[sd_TT]		= (D_INF - (D_INF - sd) * exp(-dt / TAU_D) - sd)/dt;
	// Variable_Derivatives[sf_TT]		= (F_INF - (F_INF - sf) * exp(-dt / TAU_F) - sf)/dt;
	// Variable_Derivatives[sf2_TT]		= (F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2) - sf2)/dt;
	// Variable_Derivatives[sfcass_TT]	= (FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS) - sfcass)/dt;
	// Variable_Derivatives[sr_TT]		= (S_INF - (S_INF - sr) * exp(-dt / TAU_S) - sr)/dt;
	// Variable_Derivatives[ss_TT]		= (S_INF - (S_INF - ss) * exp(-dt / TAU_S) - ss)/dt;
	// Variable_Derivatives[sxr1_TT]	= (Xr1_INF - (Xr1_INF - sxr1) * exp(-dt / TAU_Xr1) - sxr1)/dt;
	// Variable_Derivatives[sxr2_TT]	= (Xr2_INF - (Xr2_INF - sxr2) * exp(-dt / TAU_Xr2) - sxr2)/dt;

	// Variable_Derivatives[sr_TT]		= (S_INF - sr)/ TAU_S;
	// Variable_Derivatives[ss_TT]		= (S_INF - ss)/ TAU_S;
	// Variable_Derivatives[sxr1_TT]	= (Xr1_INF - sxr1)/ TAU_Xr1;
	// Variable_Derivatives[sxr2_TT]	= (Xr2_INF - sxr2)/ TAU_Xr2;
	// Variable_Derivatives[sxs_TT]		= (Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs) - sxs)/dt;
	// Variable_Derivatives[y_TT]		= (y_inf - (y_inf - y) * exp(-dt / tau_y) - y)/dt;
	// Variable_Derivatives[sRR_TT]		= dRR;
	// Variable_Derivatives[sOO_TT]		= (k1 * CaSS * CaSS * sRR / (Tent_k3 + k1 * CaSS * CaSS) -sOO)/dt;


	//Set the values of the variable derivatives
	
	


	// Variable_Derivatives[Cai_TT] 	= ((sqrt(bc * bc + 4 * cc) - bc) / 2 - Cai)/dt;
	// Variable_Derivatives[CaSRf_TT]	= ((sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2 - CaSRf)/dt;
	// Variable_Derivatives[CaSS_TT]	= ((sqrt(bcss * bcss + 4 * ccss) - bcss) / 2 - CaSS)/dt;
	Variable_Derivatives[Cai_TT] 		= dCai;			//((sqrt(bc * bc + 4 * cc) - bc) / 2 - Cai)/dt;
	Variable_Derivatives[Nai_TT]		= dNai;
	Variable_Derivatives[Ki_TT]			= dKi;
	Variable_Derivatives[CaSRf_TT]		= dCaSR;		//((sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2 - CaSRf)/dt;
	Variable_Derivatives[CaSS_TT]		= dCaSS;		//((sqrt(bcss * bcss + 4 * ccss) - bcss) / 2 - CaSS)/dt;
	Variable_Derivatives[sm_TT]			= (1.0 / TAU_M)*(M_INF - sm);
	Variable_Derivatives[sh_TT]			= (1.0 / TAU_H)*(H_INF - sh);
	Variable_Derivatives[sj_TT]			= (1.0 / TAU_J)*(J_INF - sj);
	Variable_Derivatives[sd_TT]			= (1.0 / TAU_D)*(D_INF - sd);
	Variable_Derivatives[sf_TT]			= (1.0 / TAU_F)*(F_INF - sf);
	Variable_Derivatives[sf2_TT]		= (1.0 / TAU_F2)*(F2_INF - sf2);
	Variable_Derivatives[sfcass_TT]		= (1.0 / TAU_FCaSS)*(FCaSS_INF - sfcass);
	Variable_Derivatives[sr_TT]			= (1.0 / TAU_S)*(S_INF - sr);
	Variable_Derivatives[ss_TT]			= (1.0 / TAU_S)*(S_INF - ss);
	Variable_Derivatives[sxr1_TT]		= (1.0 / TAU_Xr1)*(Xr1_INF - sxr1);
	Variable_Derivatives[sxr2_TT]		= (1.0 / TAU_Xr2)*(Xr2_INF - sxr2);
	Variable_Derivatives[sxs_TT]		= (1.0 / TAU_Xs)*(Xs_INF - sxs);
	Variable_Derivatives[y_TT]			= (1.0 / tau_y)*(y_inf - y);
	Variable_Derivatives[sRR_TT]		= dRR;

	//Set the value of the total ionic current
	Iion = -(IKr + IKs + IK1 + Ito + Isus + (INa) + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + (If) + Istim);

	for(unsigned i=0; i<Num_Cell_Vars; i++){
		if(!std::isfinite(Variable_Derivatives[i])){
			// oomph_info << "Shuffling " << i << "th derivative due to NaN" << std::endl;
			Variable_Derivatives[i] = ((double)std::rand()/(RAND_MAX));
		}
	}
	if(!std::isfinite(Iion)){
		// oomph_info << "Shuffling iion due to NaN" << std::endl;
		Iion = ((double)std::rand()/(RAND_MAX));
	}
}





void ExplicitTNNP06VentUpdated::get_optional_output(const double &Vm,
												const Vector<double> &CellVariables,
												const double &t,
												const unsigned &cell_type,
												const double &Istim,
												const Vector<double> &Other_Parameters,
												const Vector<double> &Other_Variables,

												Vector<double> &Out) const
{
	// const double sxs	= CellVariables[sxs_TT];

	// const double Nai	= CellVariables[Nai_TT];
	// const double Ki		= CellVariables[Ki_TT];

	// const double Eks = Tent_RTONF * (log((Tent_Ko + Tent_pKNa * Tent_Nao) / (Ki + Tent_pKNa * Nai)));

	// // double iks_right_ventricle_factor;
 // //    if (cell_type == 103 || cell_type == 105 || cell_type == 104)
	// // 	iks_right_ventricle_factor = 2.0;
	// // else 
	// // 	iks_right_ventricle_factor = 1.0;
	
 // //    //OURS
 // //    double Gks;
 // //    if (cell_type == 100 || cell_type == 103)		Gks = 0.392 * iks_right_ventricle_factor; 
 // //    else if (cell_type == 102 || cell_type == 105)  Gks = 0.392 * iks_right_ventricle_factor;
 // //    else if (cell_type == 101 || cell_type == 104) Gks = 0.098 * iks_right_ventricle_factor; 

	// double Gks;
	// if (cell_type == 108 || cell_type == 107 || cell_type == 106)
	// 	Gks = 0.2352;
	// else if (cell_type == 101)
	// 	Gks = 0.098;
	// else if (cell_type == 102 || cell_type == 100)
	// 	Gks = 0.392;
	// else if (cell_type == 104)
	// 	Gks = 0.098 * 2;
	// else
	// 	Gks = 0.392 * 2;

	// // Xs_INF = 1. / (1. + exp((-5. - Vm) / 14.));
	// // Axs = (1400. / (sqrt(1. + exp((5. - Vm) / 6))));
	// // Bxs = (1. / (1. + exp((Vm - 35.) / 15.)));
	// // TAU_Xs = Axs * Bxs + 80;
	
	// // sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);
	
	// Out[EKs_TT] = Eks;
 //    Out[IKs_TT] = Gks * sxs * sxs * (Vm - Eks);
}






} //end namespace
