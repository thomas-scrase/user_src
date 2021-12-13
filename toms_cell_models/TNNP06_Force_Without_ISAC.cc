#include "TNNP06_Force_Without_ISAC.h"


namespace oomph{

TNNP06_Force_Without_ISAC::TNNP06_Force_Without_ISAC()
{
	//The active strain is the first (and only) other variable the model generates
	active_strain_index = 0;


	// constants
	TTCell_Ko = 5.4;
	TTCell_Cao = 2.0;
	TTCell_Nao = 140.0;
	TTCell_Vc = 0.016404;
	TTCell_Vsr = 0.001094;
	TTCell_Vss = 0.00005468;
	TTCell_Bufc = 0.2;
	TTCell_Kbufc = 0.001;
	TTCell_Bufsr = 10.;
	TTCell_Kbufsr = 0.3;
	TTCell_Bufss = 0.4;
	TTCell_Kbufss = 0.00025;
	TTCell_Vmaxup = 0.006375;
	TTCell_Kup = 0.00025;
	TTCell_Vrel = 0.102; //40.8;
	TTCell_k1_ = 0.15;
	TTCell_k2_ = 0.045;
	TTCell_k3 = 0.060;
	TTCell_k4 = 0.005; //0.000015;
	TTCell_EC = 1.5;
	TTCell_maxsr = 2.5;
	TTCell_minsr = 1.;
	TTCell_Vleak = 0.00036;
	TTCell_Vxfer = 0.0038;
	TTCell_R = 8314.472;
	TTCell_F = 96485.3415;
	TTCell_T = 310.0;
	TTCell_CAPACITANCE = 0.185;
	TTCell_pKNa = 0.03;
	TTCell_GbNa = 0.00029;
	TTCell_KmK = 1.0;
	TTCell_KmNa = 40.0;
	TTCell_knak = 2.724;
	TTCell_GCaL = 0.00003980;
	TTCell_GbCa = 0.000592;
	TTCell_knaca = 1000;
	TTCell_KmNai = 87.5;
	TTCell_KmCa = 1.38;
	TTCell_ksat = 0.1;
	TTCell_n = 0.35;
	TTCell_GpCa = 0.1238;
	TTCell_KpCa = 0.0005;
	TTCell_GpK = 0.0146;
	TTCell_RTONF = (TTCell_R *TTCell_T) / TTCell_F;
	TTCell_inverseVcF2 = 1.0 / (2.0 * TTCell_Vc *TTCell_F);
	TTCell_inverseVcF = 1.0 / (TTCell_Vc *TTCell_F);
	TTCell_inversevssF2 = 1.0 / (2.0 * TTCell_Vss *TTCell_F);
	//Rice myofilament constants
	p0 = 1.754251548964904;
	p1 = 0.905622641626625;
	p2 = 0.499437793063966;
	p3 = 0.400000001127317;
	p4 = 1.000000000000000;
	p5 = 1.000000000000000;
	p6 = 1.981229252026256;
	p7 = 0.511387864324941;
	p8 = 0.023420000000000;
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



	//Assign the names of the variables used by this model
	Names_Of_Cell_Variables =
	{
		"Vm",
		"cai",
		"CaSR",
		"CaSS",
		"Nai",
		"Ki",
		"sm",
		"sh",
		"sj",
		"sxr1",
		"sxr2",
		"sxs",
		"ss",
		"sr",
		"sd",
		"sf",
		"sf2",
		"sfcass",
		"sRR",
		"mNaL",
		"hNaL",
		"wt_C3",
		"wt_O",
		"wt_C1",
		"wt_C2",
		"wt_I",
		"N",
		"XBprer",
		"XBpostr",
		"SL",
		"xXBpostr",
		"xXBprer",
		"TRPNCaL",
		"TRPNCaH",
		"intf"
	};
	Names_Of_Other_Parameters =
	{
		"ABIndex",
		"IS_Index",
		"RV_Index"
	};
	Names_Of_Other_Variables =
	{

	};
	Names_Of_Output_Data =
	{
		"Strain"
	};

	FinalizeConstruction();
}

TNNP06_Force_Without_ISAC::~TNNP06_Force_Without_ISAC()
{

}

double TNNP06_Force_Without_ISAC::return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
{
	const double Vars[50] = {8.37604e-05, 2.43747, 0.000187031, 8.66537, 135.875, 0.00157112, 0.756077, 0.755924, 0, 1,
			6.89887e-24, 1.74656e-18, 0.156871, 0.369442, 0.326288, 0.128122, 0.0189051, 9.24159e-05, 0.000163242, 9.61494e-05, 1.89167e-05, 2.04175e-08, 2.40518e-08, 7.09812e-09, 2.00552e-12, 1.18375e-12, 7.40385e-17,
			0.999998, 2.25571e-08, 3.19131e-05, 0.988368, 0.999517, 0.999982, 0.991112, 0.00029464, 0.427938, 0.998997, 5.46719e-05, 0.000619172, 0.000306603, 2.2339e-05, 0.999848, 2.1554e-05, 0.000114638, 2.19965, 0.00700003, 1.84508e-08, 0.0181649, 0.156304, -0.00105964};

	return Vars[v];
}

double TNNP06_Force_Without_ISAC::return_initial_membrane_potential(const unsigned &cell_type)
{
	return -85.6435;
}


void TNNP06_Force_Without_ISAC::Calculate_Derivatives(const Boost_State_Type &Variables,
									const double &t,
									const unsigned &cell_type,
									const double &Istim,
									const Vector<double> &Other_Parameters,
									const Vector<double> &Other_Variables,
									Vector<double> &Variable_Derivatives,
									double &Iion)
{
	const double cai 	  = Variables[cai_TT]; // 0.00007;
	const double CaSR   = Variables[CaSR_TT]; // 1.3;
	const double CaSS   = Variables[CaSS_TT]; // 0.00007;
	const double Nai    = Variables[Nai_TT]; // 7.67;
	const double Ki     = Variables[Ki_TT]; // 138.3;
	const double sm     = Variables[sm_TT]; // 0.0;
	const double sh     = Variables[sh_TT]; // 0.75;
	const double sj     = Variables[sj_TT]; // 0.75;
	const double sxr1   = Variables[sxr1_TT]; // 0.0;
	const double sxr2   = Variables[sxr2_TT]; // 1.0;
	const double sxs    = Variables[sxs_TT]; // 0.0;
	const double ss     = Variables[ss_TT]; // 1.0;
	const double sr     = Variables[sr_TT]; // 0.0;
	const double sd     = Variables[sd_TT]; // 0.0;
	const double sf     = Variables[sf_TT]; // 1.0;
	const double sf2    = Variables[sf2_TT]; // 1.0;
	const double sfcass = Variables[sfcass_TT]; // 1.0;
	const double sRR    = Variables[sRR_TT]; // 1.0;
	const double mNaL   = Variables[mNaL_TT];
	const double hNaL   = Variables[hNaL_TT];
	const double wt_C3  = Variables[wt_C3_TT];  // 1.0;
	const double wt_O   = Variables[wt_O_TT];  // 0.0;
	const double wt_C1  = Variables[wt_C1_TT];  // 0.0;
	const double wt_C2  = Variables[wt_C2_TT];  // 0.0;
	const double wt_I   = Variables[wt_I_TT];  // 0.0;

	///Rice myofillament variables
	const double N         = Variables[N_TT];
	const double XBprer    = Variables[XBprer_TT];
	const double XBpostr   = Variables[XBpostr_TT];
	const double SL        = Variables[SL_TT];
	const double xXBpostr  = Variables[xXBpostr_TT];
	const double xXBprer   = Variables[xXBprer_TT];
	const double TRPNCaL   = Variables[TRPNCaL_TT];
	const double TRPNCaH   = Variables[TRPNCaH_TT];
	const double intf      = Variables[intf_TT];

	const double Vm = Variables[Num_Cell_Vars];

	double ABIndex = Other_Parameters[0];
	double IS_Index = Other_Parameters[1];
	double RV_Index = Other_Parameters[2];


	double GTo_ABh = 1.0;
	double InAcTo_Vhalf_ABh = 0.0;
	double InAcTo_Vk_ABh = 1.0;
	double GKs_ABh = 1.0;
	double TauKs_ABh = 1.0;
	double AcKs_Vhalf_ABh = 0.0;
	double AB_IKr_ActVhalf = 0.0;


	if (ABIndex < 0.0 or ABIndex > 1.0) {
		std::cerr << "Wrong Apcial-Basal Ratio, at ABIndex = " << ABIndex << std:: endl;
		std::cerr << "Program Existing... " << std::endl;
		std::exit(0);
	}
	/* V1, ABIndex   == 0.0 -> Apical Cells */
	/*GTo_ABh          = (1.0 + (1 - ABIndex) * (29.6 / 16.5 - 1));      // Apical cells: 29.6/16.5; Basal Cells: 1.0
	InAcTo_Vhalf_ABh = -(1 - ABIndex) * 4.0;                          // Apical Cells: -4; Basal Cells: 0.0;
	InAcTo_Vk_ABh    = 1.0 + (1 - ABIndex) * ((4.5 - 3.4) / 3.4);   // Apical cells: 4.5/3.4; Basal Cells: 1.0
	GKs_ABh          = 1.0 + (1 - ABIndex) * (5.6 / 2.1 - 1.0);      // Apical Cells: 5.6/2.1; Basal Cells: 1.0
	// TauKs_ABh        = 1.0 + (1 - ABIndex) * ((358.0 - 516.0) / 516.0);      // Apical Cells:358/516 ; Basal Cells: 1.0
	AB_IKr_ActVhalf = (1-ABIndex) * 4.9;									// Apical cells: 4.9; Basal cells 0.0
	*/


	// longer action potential

	/*GTo_ABh          = 1.0 -  ABIndex * (1.0 - 16.5 / 29.6); // Apical cells: 1.0; Basal Cells: 16.5/29.6
	InAcTo_Vhalf_ABh = ABIndex * 4.0;                          // Apical Cells: 0.0; Basal Cells: 4.0;
	InAcTo_Vk_ABh    = 1.0 -  ABIndex * (1.0 - 3.4 / 4.5); // Apical cells: 1.0; Basal Cells: 3.4/4.5
	GKs_ABh          = 1.0 -  ABIndex * (1.0 - 2.1 / 4.2);//1.0 -  ABIndex * (1.0 - 2.1 / 5.6);  // Apical Cells: 1.0; Basal Cells: 2.1/5.6
	TauKs_ABh = 1.0 + ABIndex * (516.0 / 318.0 - 1); // Apical Cells:1.0 ; Basal Cells: 516/358
	AB_IKr_ActVhalf = -ABIndex * 4.9; // Apical cells: 0.0; Basal cells -4.9*/



	//version three , find the middel value;
	/*{double x;
		x = (29.6-16.5)/(29.6+16.5);
		GTo_ABh          = 1.0 -  (ABIndex-0.5) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
		InAcTo_Vhalf_ABh = (ABIndex-0.5) * 4.0;                          // Apical Cells: -2.0; Basal Cells: 2.0;
		x = (4.5-3.4)/(4.5+3.4);
		InAcTo_Vk_ABh    = 1.0 -  (ABIndex-0.5) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
		x = (5.6-2.1)/(5.6+2.1);
		GKs_ABh          =  1.0 -  (ABIndex-0.5) * x+0.1;     // Apical Cells: 1.0; Basal Cells: 2.1/5.6
		x = (516-318)/(516+318);
		TauKs_ABh = 1.0 +  (ABIndex-0.5) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
		// TauKs_ABh = 1.0 - ABIndex*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
		AB_IKr_ActVhalf = -(ABIndex-0.5) * 4.9; // Apical cells: 2.45; Basal cells -2.45
	} */


	// AB_IKr_ActVhalf = -ABIndex * 4.9; // Apical cells: 0.0; Basal cells -4.9

	// version four, middle 2.
	/*GTo_ABh          = 1.0 -  ABIndex * (1.0 - 16.5 / 29.6); // Apical cells: 1.0; Basal Cells: 16.5/29.6
	InAcTo_Vhalf_ABh = ABIndex * 4.0;                          // Apical Cells: 0.0; Basal Cells: 4.0;
	InAcTo_Vk_ABh    = 1.0 -  ABIndex * (1.0 - 3.4 / 4.5); // Apical cells: 1.0; Basal Cells: 3.4/4.5
	GKs_ABh          = 1.5 -  ABIndex * (1.0 - 2.1 / 5.6);//1.0 -  ABIndex * (1.0 - 2.1 / 5.6);  // Apical Cells: 1.0; Basal Cells: 2.1/5.6
	TauKs_ABh = 1.0 + ABIndex * (516.0 / 318.0 - 1); // Apical Cells:1.0 ; Basal Cells: 516/358*/
	{
		double x;
		double A = 0.4;  // was  0.5 for the LV, here change to 0.35, to increase APD, see what happens.

		if (cell_type == 105 or cell_type == 104 or cell_type == 103)
		{
			A = 0.4;
		}

		if (cell_type == 103 or cell_type == 100) {
			A = 0.2;
		}
		x = (29.6 / 16.5 - 1) / (A + 29.6 / 16.5 * (1 - A) );
		GTo_ABh          = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
		InAcTo_Vhalf_ABh = (ABIndex - 1.0) * 4.0;                        // Apical Cells: -2.0; Basal Cells: 2.0;
		x = (4.5 / 3.4 - 1) / (A + 4.5 / 3.4 * (1 - A) );
		InAcTo_Vk_ABh    = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
		x = (5.6 / 2.1 - 1) / (A + 5.6 / 2.1 * (1 - A) );
		GKs_ABh          =  1.0 -  (ABIndex - A) * x;   // Apical Cells: 1.0; Basal Cells: 2.1/5.6
		// x = (516-318)/(516+318);
		x = (358.0 / 516.0 - 1) / (A + 358.0 / 516.0 * (1 - A) );
		// TauKs_ABh = 1.0 -  (ABIndex-1.0) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
		// TauKs_ABh = 1.0 - ABIndex*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
		AB_IKr_ActVhalf = -(ABIndex - 1.0) * 4.9; // Apical cells: 2.45; Basal cells -2.45
	}





	/*{
		double x;
		double A = 0.5;

		if (cell_type == 105 or cell_type == 104 or cell_type == 103)
		{
			A = 0.4;
		}
		x = (29.6 / 16.5 - 1) / (A + 29.6 / 16.5 * (1 - A) );
		GTo_ABh          = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 16.5/29.6
		InAcTo_Vhalf_ABh = (ABIndex - 1.0) * 4.0;                        // Apical Cells: -2.0; Basal Cells: 2.0;
		x = (4.5 / 3.4 - 1) / (A + 4.5 / 3.4 * (1 - A) );
		InAcTo_Vk_ABh    = 1.0 -  (ABIndex - A) * x; // Apical cells: 1.0; Basal Cells: 3.4/4.5
		x = (5.6 / 2.1 - 1) / (A + 5.6 / 2.1 * (1 - A) );
		GKs_ABh          =  1.0 -  (ABIndex - A) * x;   // Apical Cells: 1.0; Basal Cells: 2.1/5.6
		// x = (516-318)/(516+318);
		x = (358.0 / 516.0 - 1) / (A + 358.0 / 516.0 * (1 - A) );
		// TauKs_ABh = 1.0 -  (ABIndex-1.0) * x; // Apical Cells:1.0 ; Basal Cells: 516/358
		// TauKs_ABh = 1.0 - ABIndex*(516.0/318.0-1);// Apical Cells:1.0 ; Basal Cells: 516/358
		AB_IKr_ActVhalf = -(ABIndex - 1.0) * 4.9; // Apical cells: 2.45; Basal cells -2.45
	}
	*/

	/* end of Apical-Basal heterogeneity settings */
	/*if (cell_type == 105 or cell_type == 104 or cell_type == 103)
	{
		if (ABIndex < 0.5)
		{
			ABIndex = 0.5;
		}
		GTo_ABh          = (1.0 + (1 - ABIndex) * (29.6 / 16.5 - 1));      // Apical cells: 29.6/16.5; Basal Cells: 1.0
		InAcTo_Vhalf_ABh = -(1 - ABIndex) * 4.0;                          // Apical Cells: -4; Basal Cells: 0.0;

	}*/

	/*implementation of Ischemia effects: */

	/*reference:
	Ivan V. Kazbanov, et al. 2014. Plos One. Effect of Global Cardiac Ischemia on Human Ventricular
											 Fibrillation: Insights from a Multi-scale Mechanistic
											 Model of the Human Heart */
	double Ischemia_TTCell_Ko = TTCell_Ko + IS_Index * 5.4;
	double Ek                 = TTCell_RTONF * (log((Ischemia_TTCell_Ko / Ki)));
	float f_atp               = 0.0055 * IS_Index;
	double GKATP              = 155.0;
	double IKATP              = GKATP * f_atp * pow(Ischemia_TTCell_Ko / 5.4, 0.3) * (Vm - Ek) / (40 + 3.5 * exp(0.025 * Vm));
	float Acidosis_factor     = 1.0  - 0.2 * IS_Index; // Acidosis_factor is the scaling factor to multiply INa and ICaL, to represent the effect of Acidosis following ischemia...



	/*end of ischemia implementation */



	double Ena = TTCell_RTONF * (log((TTCell_Nao / Nai)));
	double Eks = TTCell_RTONF * (log((Ischemia_TTCell_Ko + TTCell_pKNa * TTCell_Nao) / (Ki + TTCell_pKNa * Nai)));
	// std::cout << TTCell_RTONF << "\t" << TTCell_Cao << "\t" << cai << std::endl;
	double Eca = 0.5 * TTCell_RTONF * (log((TTCell_Cao / cai)));
	double Ak1 = 0.1 / (1. + exp(0.06 * (Vm - Ek - 200)));
	double Bk1 = (3.0 * exp(0.0002 * (Vm - Ek + 100)) + exp(0.1 * (Vm - Ek - 10))) / (1. + exp(-0.5 * (Vm - Ek)));
	double rec_iK1 = Ak1 / (Ak1 + Bk1);
	double rec_iNaK = (1. / (1. + 0.1245 * exp(-0.1 * Vm * TTCell_F / (TTCell_R * TTCell_T)) + 0.0353 * exp(-Vm * TTCell_F / (TTCell_R * TTCell_T))));
	double rec_ipK = 1. / (1. + exp((25 - Vm) / 5.98));

	double R_INF;
	double S_INF;
	double TAU_R;
	double TAU_S;
	double GNaL, Gks, Gto;

	double Mfactor = 2.0;

	if (cell_type == 100) {
		GNaL = 0.0065;
		Gks = 0.392 * 1.4;  // was 0.392, changed in the same way as the ORd model did, to increase the APD difference in EPI and ENDO by haibo
		Gto = 0.294;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.*exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (cell_type == 102) {
		GNaL = 0.0065;
		Gks = 0.392;
		Gto = 0.073;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 28 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 1000.*exp(-(Vm + 67) * (Vm + 67) / 1000.) + 8.;
	} else if (cell_type == 101) {
		GNaL = 0.0095;
		Gks = 0.098 * Mfactor;
		Gto = 0.294;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.0 * exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (cell_type == 105) {
		GNaL = 0.0065;
		Gks = 0.392 * (1.0 + 0.8 * RV_Index); // was 2.0
		Gto = 0.073 * (1.0 + 0.8 * RV_Index); // was 3.0
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 28 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 1000.*exp(-(Vm + 67) * (Vm + 67) / 1000.) + 8.;
	} else if (cell_type == 104) {
		GNaL = 0.0095;
		Gks = 0.098 * Mfactor *  (1.0 + 0.8 * RV_Index); // was  0.098 * 2.0
		Gto = 0.294 * (1.0 + 0.8 * RV_Index); // was 3.0;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.0 * exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (cell_type == 103) {
		GNaL = 0.0065;
		Gks = 0.392 * 1.4 * (1.0 + 0.8 * RV_Index); // was 0.392, changed in the same way as the ORd model did, by haibo // was 0.392 * 1.4 * 2.0
		Gto = 0.294 * (1.0 + 0.8 * RV_Index); // was 3.0
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.*exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else {
		std::cerr << "Cell type wrong!" << std::endl;
		// std::exit(0);
		throw OomphLibError("Cell type wrong!",
	                       OOMPH_CURRENT_FUNCTION,
	                       OOMPH_EXCEPTION_LOCATION);
	}
	// cout << Gto << endl;
	/* Ito */
	double Ito = GTo_ABh * Gto * sr * ss * (Vm - Ek);
	// ss = S_INF - (S_INF - ss) * exp(-dt / TAU_S);
	// sr = R_INF - (R_INF - sr) * exp(-dt / TAU_R);
	

	/* INa */
	double GNa = 14.838;
	if (cell_type == 104 or cell_type == 101)
	{
		// GNa *= 1.5;  // tranmural hereterogeneity includes INa
	}
	double INa = Acidosis_factor * GNa * sm * sm * sm * sh * sj * (Vm - Ena);
	// double AM = 1. / (1. + exp((-60. - Vm) / 5.));
	// double BM = 0.1 / (1. + exp((Vm + 35.) / 5.)) + 0.10 / (1. + exp((Vm - 50.) / 200.));
	double TAU_M = (1. / (1. + exp((-60. - Vm) / 5.))) * ( 0.1 / (1. + exp((Vm + 35.) / 5.)) + 0.10 / (1. + exp((Vm - 50.) / 200.))) ;
	double TAU_H;
	double TAU_J;
	double M_INF = 1. / ((1. + exp((-56.86 - Vm) / 9.03)) * (1. + exp((-56.86 - Vm) / 9.03)));
	if (Vm >= -40.) {
		double AH_1 = 0.;
		double BH_1 = (0.77 / (0.13 * (1. + exp(-(Vm + 10.66) / 11.1))));
		TAU_H = 1.0 / (AH_1 + BH_1);
	} else {
		double AH_2 = (0.057 * exp(-(Vm + 80.) / 6.8));
		double BH_2 = (2.7 * exp(0.079 * Vm) + (3.1e5) * exp(0.3485 * Vm));
		TAU_H = 1.0 / (AH_2 + BH_2);
	}
	double H_INF = 1. / ((1. + exp((Vm + 71.55) / 7.43)) * (1. + exp((Vm + 71.55) / 7.43)));
	if (Vm >= -40.) {
		double AJ_1 = 0.;
		double BJ_1 = (0.6 * exp((0.057) * Vm) / (1. + exp(-0.1 * (Vm + 32.))));
		TAU_J = 1.0 / (AJ_1 + BJ_1);
	} else {
		double AJ_2 = (((-2.5428e4) * exp(0.2444 * Vm) - (6.948e-6) *
		                exp(-0.04391 * Vm)) * (Vm + 37.78) /
		               (1. + exp(0.311 * (Vm + 79.23))));
		double BJ_2 = (0.02424 * exp(-0.01052 * Vm) / (1. + exp(-0.1378 * (Vm + 40.14))));
		TAU_J = 1.0 / (AJ_2 + BJ_2);
	}
	double J_INF = H_INF;

	// sm = M_INF - (M_INF - sm) * exp(-dt / TAU_M);
	// sh = H_INF - (H_INF - sh) * exp(-dt / TAU_H);
	// sj = J_INF - (J_INF - sj) * exp(-dt / TAU_J);
	
	

	/* IKr */  // CalculateTTIKr
	/*double Xr1_INF;
	double Xr2_INF;
	double TAU_Xr1;
	double TAU_Xr2;
	double Gkr = 0.153;
	Xr1_INF = 1. / (1. + exp((-26. - Vm) / 7.));
	double axr1 = 450. / (1. + exp((-45. - Vm) / 10.));
	double bxr1 = 6. / (1. + exp((Vm - (-30.)) / 11.5));
	TAU_Xr1 = axr1 * bxr1;
	Xr2_INF = 1. / (1. + exp((Vm - (-88.)) / 24.));
	double axr2 = 3. / (1. + exp((-60. - Vm) / 20.));
	double bxr2 = 1.12 / (1. + exp((Vm - 60.) / 20.));
	TAU_Xr2 = axr2 * bxr2;
	sxr1 = Xr1_INF - (Xr1_INF - sxr1) * exp(-dt / TAU_Xr1);
	sxr2 = Xr2_INF - (Xr2_INF - sxr2) * exp(-dt / TAU_Xr2);*/
	// double IKr =  Gkr * sqrt(Ischemia_TTCell_Ko / 5.4) * sxr1 * sxr2 * (Vm - Ek);
	// CalculateMarkovIKr();

	/* Markov IKr Model */
	double Gkr = 0.153;
	double m_epiMidRatio = 1.5;
	double epi_factor   = 1.8 * m_epiMidRatio;
	double endo_factor  = 1.8;
	double mcell_factor = 1.8;
	double wt_a1 =  2.172;
	double wt_b1 =  1.077;
	double wt_a2 =  0.00655   * exp(0.5 * 0.05547153 * (Vm - 36. - AB_IKr_ActVhalf));
	double wt_a  =  0.00555   * exp(0.05547153 * (Vm - 12. - AB_IKr_ActVhalf));
	double wt_b  =  0.002357  * exp(-0.036588 * (Vm - AB_IKr_ActVhalf));
	double wt_b2 = 0.65 * 0.0029357 * exp(0.69 * -0.02158 * (Vm - AB_IKr_ActVhalf));
	double wt_ai = 0.11 * 0.439     * exp(1.7 *  -0.02352 * (Vm + 25. - AB_IKr_ActVhalf)) * (4.5 / Ischemia_TTCell_Ko);
	double wt_bi = 0.4 *  0.656     * exp(      0.000942 * (Vm - AB_IKr_ActVhalf)) * ((pow((4.5 / Ischemia_TTCell_Ko), 0.3)));
	double wt_mu = (wt_ai * wt_b2) / wt_bi;
	if (cell_type == 102 or cell_type == 105)
		Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * endo_factor;
	else if (cell_type == 101 or cell_type == 104)
		Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * mcell_factor;
	else
		Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * epi_factor;

	double IKr = Gkr * wt_O * (Vm - Ek);

	

	// {	
		double wt_dC3 = (wt_b * wt_C2) - (wt_a * wt_C3);
		double wt_dC2 = -((wt_b + wt_a1) * wt_C2) + (wt_a * wt_C3) + (wt_b1 * wt_C1);
		double wt_dC1 = -((wt_b1 + wt_a2 + wt_a2) * wt_C1) + (wt_a1 * wt_C2) + (wt_b2 * wt_O) + (wt_mu * wt_I);
		double wt_dO  =  -((wt_b2 + wt_bi) * wt_O) + (wt_a2 * wt_C1) + (wt_ai * wt_I);
		double wt_dI  = -((wt_mu + wt_ai) * wt_I) + (wt_a2 * wt_C1) + (wt_bi * wt_O);

		// wt_O  += dt * wt_dO;
		// wt_C1 += dt * wt_dC1;
		// wt_C2 += dt * wt_dC2;
		// wt_C3 += dt * wt_dC3;
		// wt_I  += dt * wt_dI;
	// }
	
	

	/* IKs */
	double IKs = GKs_ABh * Gks * sxs * sxs * (Vm - Eks);  // with apicalbasal heterogeneity wit direct scaling factor, by haibo
	double Xs_INF = 1. / (1. + exp((-5. - Vm) / 14.));
	// double Axs = (1400. / (sqrt(1. + exp((5. - Vm) / 6))));
	// double Bxs = (1. / (1. + exp((Vm - 35.) / 15.)));
	double TAU_Xs = (1400. / (sqrt(1. + exp((5. - Vm) / 6)))) * (1. / (1. + exp((Vm - 35.) / 15.))) + 80;
	TAU_Xs = TauKs_ABh * TAU_Xs;
	// sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);

	
	

	/* ICaL */
	double ICaL = Acidosis_factor * TTCell_GCaL * sd * sf * sf2 * sfcass * 4 * (Vm - 15) * (TTCell_F * TTCell_F / (TTCell_R * TTCell_T)) *
	              (0.25 * exp(2 * (Vm - 15) * TTCell_F / (TTCell_R * TTCell_T)) * CaSS - TTCell_Cao) / (exp(2 * (Vm - 15) * TTCell_F / (TTCell_R * TTCell_T)) - 1.);

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
	double Af2 = 600 * exp(-(Vm + 25) * (Vm + 25) / 170);
	double Bf2 = 31 / (1. + exp((25 - Vm) / 10));
	double Cf2 = 16 / (1. + exp((Vm + 30) / 10));
	double TAU_F2 = Af2 + Bf2 + Cf2;

	// sd = D_INF - (D_INF - sd) * exp(-dt / TAU_D);
	// sf = F_INF - (F_INF - sf) * exp(-dt / TAU_F);
	// sf2 = F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2);

	double FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
	double TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;
	// sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS);

	
	//INaL gates

	// double mNaL_INF, hNaL_INF, TAU_mNaL, TAU_hNaL;
	// need to update INaL later...
	// mNaL = mNaL_INF - (mNaL_INF - mNaL) * exp(-dt / TAU_mNaL);
	// hNaL = hNaL_INF - (hNaL_INF - hNaL) * exp(-dt / TAU_hNaL);
	// double INaL = GNaL * mNaL * mNaL * mNaL * hNaL * (Vm - Ena);
	double INaL = 0.0;



	/*INaL was taken from ORd 2011 model. */
	INaL = (1 + 10 * IS_Index) * GNaL * mNaL * hNaL * (Vm - Ena) ;
	double TMP_INF = 1.0 / (1.0 + exp((-(Vm + 42.85)) / 5.264));

	// mNaL = TMP_INF - (TMP_INF - mNaL) * exp(-dt / TAU_M);  // double tmL = TAU_M;
	double TMP_INF_1 = 1.0 / (1.0 + exp((Vm + 87.61) / 7.488));

	// hNaL = TMP_INF - (TMP_INF - hNaL) * exp(-dt / 200.0);  // double thL = 200.0;
	/*double hLssp = 1.0 / (1.0 + exp((Vm + 93.81) / 7.488));
	double thLp = 3.0 * thL;
	hLp = hLssp - (hLssp - hLp) * exp(-dt / thLp);*/
	// double fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));

	






	double IbNa = TTCell_GbNa * (Vm - Ena);
	// std::cout << TTCell_GbCa << "\t" << Vm << "\t" << Eca << std::endl;
	double IbCa = TTCell_GbCa * (Vm - Eca);
	double IpK = TTCell_GpK * rec_ipK * (Vm - Ek);
	double IpCa = TTCell_GpCa * cai / (TTCell_KpCa + cai);
	double INaK = (1 - 0.35 * IS_Index) * TTCell_knak * (Ischemia_TTCell_Ko / (Ischemia_TTCell_Ko + TTCell_KmK)) * (Nai / (Nai + TTCell_KmNa)) * rec_iNaK;
	double INaCa = TTCell_knaca * (1. / (TTCell_KmNai * TTCell_KmNai * TTCell_KmNai + TTCell_Nao * TTCell_Nao * TTCell_Nao)) * (1. / (TTCell_KmCa + TTCell_Cao)) *
	               (1. / (1 + TTCell_ksat * exp((TTCell_n - 1) * Vm * TTCell_F / (TTCell_R * TTCell_T)))) *
	               (exp(TTCell_n * Vm * TTCell_F / (TTCell_R * TTCell_T)) * Nai * Nai * Nai * TTCell_Cao -
	                exp((TTCell_n - 1) * Vm * TTCell_F / (TTCell_R * TTCell_T)) * TTCell_Nao * TTCell_Nao * TTCell_Nao * cai * 2.5);

	double GK1 = 5.405;
	if (cell_type == 100)
		GK1 *= 1.2; // according to the ORd model supple, added by haibo
	double IK1 = GK1 * rec_iK1 * (Vm - Ek);
	// double sItot = IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL + Istim + IKATP;

	/* update concentrations */
	double kCaSR = TTCell_maxsr - ((TTCell_maxsr - TTCell_minsr) / (1 + (TTCell_EC / CaSR) * (TTCell_EC / CaSR)));
	double k1 = TTCell_k1_ / kCaSR;
	double k2 = TTCell_k2_ * kCaSR;
	double dRR = TTCell_k4 * (1 - sRR) - k2 * CaSS * sRR;
	// sRR += dt * dRR;
	double sOO = k1 * CaSS * CaSS * sRR / (TTCell_k3 + k1 * CaSS * CaSS);

	double Irel = TTCell_Vrel * sOO * (CaSR - CaSS);
	double Ileak = TTCell_Vleak * (CaSR - cai);
	double Iup = TTCell_Vmaxup / (1. + ((TTCell_Kup * TTCell_Kup) / (cai * cai)));
	double Ixfer = TTCell_Vxfer * (CaSS - cai);

	// {	double CaCSQN = TTCell_Bufsr * CaSR / (CaSR + TTCell_Kbufsr);
		// double dCaSR = dt * (Iup - Irel - Ileak);
		// double bjsr = TTCell_Bufsr - CaCSQN - dCaSR - CaSR + TTCell_Kbufsr;
		// double cjsr = TTCell_Kbufsr * (CaCSQN + dCaSR + CaSR);
		// CaSR = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;
	// }


	// {	double CaSSBuf = TTCell_Bufss * CaSS / (CaSS + TTCell_Kbufss);
	// 	double dCaSS = dt * (-Ixfer * (TTCell_Vc / TTCell_Vss) + Irel * (TTCell_Vsr / TTCell_Vss) + (-ICaL * TTCell_inversevssF2 * TTCell_CAPACITANCE));
	// 	double bcss = TTCell_Bufss - CaSSBuf - dCaSS - CaSS + TTCell_Kbufss;
	// 	double ccss = TTCell_Kbufss * (CaSSBuf + dCaSS + CaSS);
	// 	// CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;
	// }


	// double CaBuf = TTCell_Bufc * cai / (cai + TTCell_Kbufc);
	// double dCai = dt * ((-(IbCa + IpCa - 2 * INaCa) * TTCell_inverseVcF2 * TTCell_CAPACITANCE) - (Iup - Ileak) * (TTCell_Vsr / TTCell_Vc) + Ixfer);
	// double bc = TTCell_Bufc - CaBuf - dCai - cai + TTCell_Kbufc;
	// double cc = TTCell_Kbufc * (CaBuf + dCai + cai);
	// cai = (sqrt(bc * bc + 4 * cc) - bc) / 2;


	// double dNai = -(INa + INaL + IbNa + 3 * INaK + 3 * INaCa) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	// Nai += dt * dNai;

	// double dKi = -(/*Istim +*/ IK1 + Ito + IKr + IKs - 2 * INaK + IpK + IKATP) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	// Ki += dt * dKi;





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
    // double dxXBpostr = dSL / 2.0 + xPsi / dutypostr * (x_0 + xXBprer - xXBpostr) * hfT;
    double dxXBpostr = dSL / 2.0 + xPsi / (fappT / (fappT * hfT + gxbT * hfT + gxbT * gappT + hbT * fappT + hbT * gappT + gxbT * fappT)) * (x_0 + xXBprer - xXBpostr);

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
    // XBprer   += dt*dXBprer;
    // XBpostr  += dt*dXBpostr;
    // SL       += dt*dSL;
    // xXBpostr += dt*dxXBpostr;
    // xXBprer  += dt*dxXBprer;
    // TRPNCaL  += dt*dTRPNCaL;
    // TRPNCaH  += dt*dTRPNCaH;
    // intf     += dt*dintf;

    double dCai_feeback = dTropToT;

    //End solve force model

	Variable_Derivatives[cai_TT] = ((-(IbCa + IpCa - 2 * INaCa) * TTCell_inverseVcF2 * TTCell_CAPACITANCE) - (Iup - Ileak) * (TTCell_Vsr / TTCell_Vc) + Ixfer)/(1.0 + (TTCell_Bufc*TTCell_Bufc)/pow(cai + TTCell_Bufc, 2.0));// - dTropToT/1000.0;
	Variable_Derivatives[CaSR_TT] = (Iup - Irel - Ileak)/(1.0  + (TTCell_Kbufsr*TTCell_Bufsr)/pow(CaSR + TTCell_Kbufsr, 2.0));
	Variable_Derivatives[CaSS_TT] = (-Ixfer * (TTCell_Vc / TTCell_Vss) + Irel * (TTCell_Vsr / TTCell_Vss) + (-ICaL * TTCell_inversevssF2 * TTCell_CAPACITANCE))/(1.0 + (TTCell_Bufss*TTCell_Kbufss)/pow(CaSS+TTCell_Kbufss, 2.0));
	Variable_Derivatives[Nai_TT] = -(INa + INaL + IbNa + 3 * INaK + 3 * INaCa) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	Variable_Derivatives[Ki_TT] = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK + IKATP) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	Variable_Derivatives[sm_TT] = (1.0/TAU_M)*(M_INF - sm);
	Variable_Derivatives[sh_TT] = (1.0/TAU_H)*(H_INF - sh);
	Variable_Derivatives[sj_TT] = (1.0/TAU_J)*(J_INF - sj);
	Variable_Derivatives[sxr1_TT] = 0.0; //Corresponds to the H-H type model of IKr, zero since we are using MC type model
	Variable_Derivatives[sxr2_TT] = 0.0; //Corresponds to the H-H type model of IKr, zero since we are using MC type model
	Variable_Derivatives[sxs_TT] = (1.0/TAU_Xs)*(Xs_INF - sxs);
	Variable_Derivatives[ss_TT] = (1.0/TAU_S)*(S_INF - ss);
	Variable_Derivatives[sr_TT] = (1.0/TAU_R)*(R_INF - sr);
	Variable_Derivatives[sd_TT] = (1.0/TAU_D)*(D_INF - sd);
	Variable_Derivatives[sf_TT] = (1.0/TAU_F)*(F_INF - sf);
	Variable_Derivatives[sf2_TT] = (1.0/TAU_F2)*(F2_INF - sf2);
	Variable_Derivatives[sfcass_TT] = (1.0/TAU_FCaSS)*(FCaSS_INF - sfcass);
	Variable_Derivatives[sRR_TT] = dRR;
	Variable_Derivatives[mNaL_TT] = (1.0/TAU_M)*(TMP_INF - mNaL);
	Variable_Derivatives[hNaL_TT] = (1.0/200.0)*(TMP_INF_1 - hNaL);
	Variable_Derivatives[wt_O_TT] = wt_dO; 
	Variable_Derivatives[wt_C1_TT] = wt_dC1; 
	Variable_Derivatives[wt_C2_TT] = wt_dC2; 
	Variable_Derivatives[wt_C3_TT] =wt_dC3;
	Variable_Derivatives[wt_I_TT] = wt_dI;

	Variable_Derivatives[N_TT] = 0.0;
	Variable_Derivatives[XBprer_TT] = 0.0;
	Variable_Derivatives[XBpostr_TT] = 0.0;
	Variable_Derivatives[SL_TT] = 0.0;
	Variable_Derivatives[xXBpostr_TT] = 0.0;
	Variable_Derivatives[xXBprer_TT] = 0.0;
	Variable_Derivatives[TRPNCaL_TT] = 0.0;
	Variable_Derivatives[TRPNCaH_TT] = 0.0;
	Variable_Derivatives[intf_TT] = 0.0;


	// Variable_Derivatives[N_TT] = dN;
	// Variable_Derivatives[XBprer_TT] = dXBprer;
	// Variable_Derivatives[XBpostr_TT] = dXBpostr;
	// Variable_Derivatives[SL_TT] = dSL;
	// Variable_Derivatives[xXBpostr_TT] = dxXBpostr;
	// Variable_Derivatives[xXBprer_TT] = dxXBprer;
	// Variable_Derivatives[TRPNCaL_TT] = dTRPNCaL;
	// Variable_Derivatives[TRPNCaH_TT] = dTRPNCaH;
	// Variable_Derivatives[intf_TT] = dintf;


	Iion = -(IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL + Istim + IKATP);


	// //If a derivative is non-finite then set it to some random value
	// for(unsigned i=0; i<Num_Cell_Vars; i++){
	// 	if(!std::isfinite(Variable_Derivatives[i]))
	// 	{
	// 		oomph_info << "Non-finite at index " << i << std::endl;



	// 		// oomph_info << dSL << std::endl;
	// 	 //    oomph_info << xPsi << std::endl;
	// 	 //    oomph_info << dutypostr << std::endl;
	// 	 //    oomph_info << x_0 << std::endl;
	// 	 //    oomph_info << xXBprer << std::endl;
	// 	 //    oomph_info << xXBpostr << std::endl;
	// 	 //    oomph_info << hfT << std::endl;

	// 	 //    oomph_info << "==========" << std::endl;

	// 	 //    oomph_info << fappT << std::endl;
	// 	 //    oomph_info << hfT << std::endl;
	// 	 //    oomph_info << gxbT << std::endl;
	// 	 //    oomph_info << gappT << std::endl;
	// 	 //    oomph_info << hbT << std::endl;

	// 	 //    oomph_info << "==========" << std::endl;

	// 	 //    oomph_info << hf << std::endl;
	// 		// oomph_info << hfmd << std::endl;
	// 		// oomph_info << xbmodsp << std::endl;
	// 		// oomph_info << Qhf << std::endl;
	// 		// oomph_info << Temp << std::endl;

	// 		// oomph_info << "==========" << std::endl;

	// 	 //    oomph_info << xXBprer << std::endl;
	// 	 //    oomph_info << hfmdc << std::endl;
	// 	 //    oomph_info << x_0 << std::endl;

	// 		exit(0);
	// 	}
	// }
	// //Do the same for the ionic current
	// if(!std::isfinite(Iion))
	// {
	// 	oomph_info << "Non-finite in Iion" << std::endl;
	// 	exit(0);
	// }
}


void TNNP06_Force_Without_ISAC::get_optional_output(const Boost_State_Type &Variables,
								const double &t,
								const unsigned &cell_type,
								const double &Istim,
								const Vector<double> &Other_Parameters,
								const Vector<double> &Other_Variables,
								Vector<double> &Out)
{
	// Out[0] = (SL - SLrest)/SLrest;
	Out[0] = (Variables[SL_TT] - SLset)/SLset;
}

}; //End namespace