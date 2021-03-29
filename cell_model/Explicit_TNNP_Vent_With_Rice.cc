#include "Explicit_TNNP_Vent_With_Rice.h"

namespace oomph{
ExplicitTNNPVentWithRice::ExplicitTNNPVentWithRice() : CellModelBase(){

	//The largest timestep such that this model converges
	Intrinsic_dt = 0.02;

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

    //End assign rice myofilament constants
}

bool ExplicitTNNPVentWithRice::compatible_cell_types(const unsigned& cell_type){
    switch(cell_type){
        case 100 : return true; //LVEPI
        case 101 : return true; //LVMCELL
        case 102 : return true; //LVENDO
        case 103 : return true; //RVEPI
        case 104 : return true; //RVMCELL
        case 105 : return true; //RVENDO
        case 106 : return true; //PFI
        case 107 : return true; //PFMB
        case 108 : return true; //PF
    }
    return false;
}

//return the initial membrane potential for the cell_typeth cell type
inline void ExplicitTNNPVentWithRice::return_initial_membrane_potential(double &v, const unsigned &cell_type){
	// TNNP_explicit_With_Rice_ICs, from timestep 0.01ms, and S1 protocol of 50 of length 1000ms.
	// Stimulus was of strength -200, with a duration of 0.5ms
	double Potents[6] = {
		//Cell Type 100
		-85.6435,
		//Cell Type 101
		-85.3843,
		//Cell Type 102
		-85.5237,
		//Cell Type 103
		-85.5265,
		//Cell Type 104
		-85.3843,
		//Cell Type 105
		-85.5237
	};

	if(cell_type>=100 && cell_type<106){v = Potents[cell_type-100];}
	// v = -86.2;
}

//Return the initial condition for the nth variable and cell_typeth cell type
inline bool ExplicitTNNPVentWithRice::return_initial_state_variable(const unsigned &n, double &v, const unsigned &cell_type){
	// TNNP_explicit_With_Rice_ICs, from timestep 0.005ms, and S1 protocol of 50 of length 1000ms.
	// Stimulus was of strength -200, with a duration of 0.5ms
	double Vars[6][50] = 
	{
		//Cell Type 100
		{8.37604e-05, 2.43747, 0.000187031, 8.66537, 135.875, 0.00157112, 0.756077, 0.755924, 0, 1,
			6.89887e-24, 1.74656e-18, 0.156871, 0.369442, 0.326288, 0.128122, 0.0189051, 9.24159e-05, 0.000163242, 9.61494e-05, 1.89167e-05, 2.04175e-08, 2.40518e-08, 7.09812e-09, 2.00552e-12, 1.18375e-12, 7.40385e-17, //MC IKS
			0.999998, 2.25571e-08, 3.19131e-05, 0.988368, 0.999517, 0.999982, 0.991112, 0.00029464, 0.427938, 0.998997, 5.46719e-05, 0.000619172, 0.000306603, 2.2339e-05, 0.999848, 2.1554e-05, 0.000114638, 2.19965, 0.00700003, 1.84508e-08, 0.0181649, 0.156304, -0.00105964},
		//Cell Type 101
		{9.49936e-05, 3.02001, 0.000202782, 8.72602, 135.707, 0.00166011, 0.749115, 0.748807, 0, 1,
			1.41158e-23, 3.41548e-18, 0.115371, 0.330199, 0.354426, 0.169183, 0.0303846, 8.64031e-05, 0.000185485, 0.000132811, 3.18036e-05, 2.42679e-08, 3.47526e-08, 1.24833e-08, 3.03125e-12, 2.17772e-12, 1.42485e-16, //MC IKS
			0.999996, 2.35533e-08, 3.30352e-05, 0.981555, 0.999499, 0.999977, 0.989834, 0.000309505, 0.415331, 0.998728, 0.000146077, 0.000713499, 0.000352412, 6.01965e-05, 0.999771, 3.2432e-05, 0.000172512, 2.19947, 0.00700005, 3.06805e-08, 0.0205513, 0.173646, -0.00159578},
		//Cell Type 102
		{8.45909e-05, 2.48002, 0.00018953, 8.57782, 135.955, 0.00161163, 0.752879, 0.75267, 0, 1,
			8.82515e-24, 2.18799e-18, 0.144564, 0.3594, 0.335085, 0.138911, 0.0216485, 9.1776e-05, 0.000171134, 0.000106416, 2.2113e-05, 2.18503e-08, 2.71744e-08, 8.47026e-09, 2.31308e-12, 1.442e-12, 9.20701e-17, //MC IKS
			0.699156, 2.3012e-08, 3.24267e-05, 0.98727, 0.999509, 0.999981, 0.990825, 0.000301417, 0.422925, 0.998918, 8.02978e-05, 0.000648029, 0.000320643, 3.29355e-05, 0.999843, 2.22677e-05, 0.000118438, 2.19963, 0.00700003, 1.9706e-08, 0.0183418, 0.157616, -0.00109481},
		//Cell Type 103
		{8.39527e-05, 2.44599, 0.00018794, 8.64708, 135.892, 0.00161068, 0.752954, 0.752759, 0, 1,
			7.80986e-24, 1.93723e-18, 0.155313, 0.368233, 0.32741, 0.129432, 0.0192315, 9.39584e-05, 0.000167084, 9.90781e-05, 1.96288e-05, 2.13166e-08, 2.52808e-08, 7.51285e-09, 2.15022e-12, 1.27801e-12, 8.15355e-17, //MC IKS
			0.999998, 2.30013e-08, 3.24146e-05, 0.988345, 0.999509, 0.999982, 0.991101, 0.000301257, 0.423664, 0.998981, 5.71666e-05, 0.00062761, 0.000310764, 2.34694e-05, 0.999847, 2.17167e-05, 0.000115503, 2.19964, 0.00700003, 1.86343e-08, 0.0182058, 0.156607, -0.00106766},
		//Cell Type 104
		{9.49936e-05, 3.02001, 0.000202782, 8.72601, 135.707, 0.00166011, 0.749115, 0.748807, 0, 1,
			1.41159e-23, 3.4155e-18, 0.115371, 0.330199, 0.354426, 0.169183, 0.0303846, 8.64032e-05, 0.000185485, 0.000132811, 3.18036e-05, 2.4268e-08, 3.47527e-08, 1.24834e-08, 3.03126e-12, 2.17773e-12, 1.42486e-16, //MC IKS
			0.999996, 2.35534e-08, 3.30353e-05, 0.981555, 0.999499, 0.999977, 0.989834, 0.000309506, 0.415331, 0.998728, 0.000146077, 0.0007135, 0.000352412, 6.01966e-05, 0.999771, 3.2432e-05, 0.000172512, 2.19947, 0.00700005, 3.06805e-08, 0.0205513, 0.173646, -0.00159578},
		//Cell Type 105
		{8.45909e-05, 2.48002, 0.00018953, 8.57782, 135.955, 0.00161164, 0.752878, 0.75267, 0, 1,
			8.82522e-24, 2.18801e-18, 0.144564, 0.3594, 0.335085, 0.138911, 0.0216485, 9.17761e-05, 0.000171134, 0.000106417, 2.2113e-05, 2.18503e-08, 2.71745e-08, 8.47028e-09, 2.31309e-12, 1.442e-12, 9.20706e-17, //MC IKS
			0.699156, 2.30121e-08, 3.24267e-05, 0.98727, 0.999509, 0.999981, 0.990825, 0.000301417, 0.422925, 0.998918, 8.02979e-05, 0.00064803, 0.000320643, 3.29355e-05, 0.999843, 2.22677e-05, 0.000118438, 2.19963, 0.00700003, 1.9706e-08, 0.0183418, 0.157616, -0.00109481}
	};
	if(n>=0 && n<required_nodal_variables() && cell_type>=100 && cell_type<106){v = Vars[cell_type-100][n]; return true;}
	else{return false;}
}//End get initial conditions


//a virtual function to extract black box nodal parameters from the cell state
//	we implement this as virtual so that in the case this class is used in a
//	combined cell model class. Then it can be overloaded so that black box nodal
//	parameters associated with each model have a unique index
inline void ExplicitTNNPVentWithRice::extract_black_box_parameters_ExplicitTNNPVentWithRice(double &abindex,
																double &isindex,
																double &rvindex,
																CellState &Cellstate){
	abindex = Cellstate.get_black_box_nodal_parameters(0);
	isindex = Cellstate.get_black_box_nodal_parameters(1);
	rvindex = Cellstate.get_black_box_nodal_parameters(2);
}

void ExplicitTNNPVentWithRice::explicit_timestep(CellState &Cellstate, Vector<double> &new_state){
	double Vm  	  = Cellstate.get_vm(); // -86.2;
	double cai 	  = new_state[0]; // 0.00007;
	double CaSR   = new_state[1]; // 1.3;
	double CaSS   = new_state[2]; // 0.00007;
	double Nai    = new_state[3]; // 7.67;
	double Ki     = new_state[4]; // 138.3;
	double sm     = new_state[5]; // 0.0;
	double sh     = new_state[6]; // 0.75;
	double sj     = new_state[7]; // 0.75;
	double sxr1   = new_state[8]; // 0.0;
	double sxr2   = new_state[9]; // 1.0;
	double sxs    = new_state[10]; // 0.0;
	double ss     = new_state[11]; // 1.0;
	double sr     = new_state[12]; // 0.0;
	double sd     = new_state[13]; // 0.0;
	double sf     = new_state[14]; // 1.0;
	double sf2    = new_state[15]; // 1.0;
	double sfcass = new_state[16]; // 1.0;
	double sRR    = new_state[17]; // 1.0;
	double mNaL   = new_state[18];
	double hNaL   = new_state[19];
	double wt_C3  = new_state[20];  // 1.0;
	double wt_O   = new_state[21];  // 0.0;
	double wt_C1  = new_state[22];  // 0.0;
	double wt_C2  = new_state[23];  // 0.0;
	double wt_I   = new_state[24];  // 0.0;

	///Rice myofillament variables
	double N         = new_state[25];
	double XBprer    = new_state[26];
	double XBpostr   = new_state[27];
	double SL        = new_state[28];
	double xXBpostr  = new_state[29];
	double xXBprer   = new_state[30];
	double TRPNCaL   = new_state[31];
	double TRPNCaH   = new_state[32];
	double intf      = new_state[33];

	unsigned m_celltype = Cellstate.get_cell_type();

	double ABIndex;
	double IS_Index;
	double RV_Index;
	extract_black_box_parameters_ExplicitTNNPVentWithRice(ABIndex,IS_Index,RV_Index,Cellstate);
	double dt = Cellstate.get_dt();

	// std::cout << ABIndex << " " << IS_Index << " " << RV_Index << " " << dt << std::endl;

	// std::cout << ABIndex << "\t" << m_celltype << "\t" << IS_Index << "\t" << RV_Index << "\t" << dt << std::endl;

	// std::exit(0);

	/*  default values for ABIndex = 0.0 (Value set 0 == Apical Cells)
	    and ABIndex ranges from 0.0 -> 1.0 */
	/*
	*   0.0 -- Apical Cells
	*   1.0 -- Basal Cells
	*/

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

		if (m_celltype == 105 or m_celltype == 104 or m_celltype == 103)
		{
			A = 0.4;
		}

		if (m_celltype == 103 or m_celltype == 100) {
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

		if (m_celltype == 105 or m_celltype == 104 or m_celltype == 103)
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
	/*if (m_celltype == 105 or m_celltype == 104 or m_celltype == 103)
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

	if (m_celltype == 100) {
		GNaL = 0.0065;
		Gks = 0.392 * 1.4;  // was 0.392, changed in the same way as the ORd model did, to increase the APD difference in EPI and ENDO by haibo
		Gto = 0.294;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.*exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (m_celltype == 102) {
		GNaL = 0.0065;
		Gks = 0.392;
		Gto = 0.073;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 28 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 1000.*exp(-(Vm + 67) * (Vm + 67) / 1000.) + 8.;
	} else if (m_celltype == 101) {
		GNaL = 0.0095;
		Gks = 0.098 * Mfactor;
		Gto = 0.294;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.0 * exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (m_celltype == 105) {
		GNaL = 0.0065;
		Gks = 0.392 * (1.0 + 0.8 * RV_Index); // was 2.0
		Gto = 0.073 * (1.0 + 0.8 * RV_Index); // was 3.0
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 28 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 1000.*exp(-(Vm + 67) * (Vm + 67) / 1000.) + 8.;
	} else if (m_celltype == 104) {
		GNaL = 0.0095;
		Gks = 0.098 * Mfactor *  (1.0 + 0.8 * RV_Index); // was  0.098 * 2.0
		Gto = 0.294 * (1.0 + 0.8 * RV_Index); // was 3.0;
		R_INF = 1. / (1. + exp((20 - Vm) / 6.));
		S_INF = 1. / (1. + exp((Vm + 20 + InAcTo_Vhalf_ABh) / (5.0 * InAcTo_Vk_ABh)));
		TAU_R = 9.5 * exp(-(Vm + 40.) * (Vm + 40.) / 1800.) + 0.8;
		TAU_S = 85.0 * exp(-(Vm + 45.) * (Vm + 45.) / 320.) + 5. / (1. + exp((Vm - 20.) / 5.)) + 3.;
	} else if (m_celltype == 103) {
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
	ss = S_INF - (S_INF - ss) * exp(-dt / TAU_S);
	sr = R_INF - (R_INF - sr) * exp(-dt / TAU_R);
	

	/* INa */
	double GNa = 14.838;
	if (m_celltype == 104 or m_celltype == 101)
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

	sm = M_INF - (M_INF - sm) * exp(-dt / TAU_M);
	sh = H_INF - (H_INF - sh) * exp(-dt / TAU_H);
	sj = J_INF - (J_INF - sj) * exp(-dt / TAU_J);
	
	

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
	if (m_celltype == 102 or m_celltype == 105)
		Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * endo_factor;
	else if (m_celltype == 101 or m_celltype == 104)
		Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * mcell_factor;
	else
		Gkr = 0.0135 * pow(Ischemia_TTCell_Ko, 0.59) * epi_factor;

	double IKr = Gkr * wt_O * (Vm - Ek);

	

	{	
		double wt_dC3 = (wt_b * wt_C2) - (wt_a * wt_C3);
		double wt_dC2 = -((wt_b + wt_a1) * wt_C2) + (wt_a * wt_C3) + (wt_b1 * wt_C1);
		double wt_dC1 = -((wt_b1 + wt_a2 + wt_a2) * wt_C1) + (wt_a1 * wt_C2) + (wt_b2 * wt_O) + (wt_mu * wt_I);
		double wt_dO  =  -((wt_b2 + wt_bi) * wt_O) + (wt_a2 * wt_C1) + (wt_ai * wt_I);
		double wt_dI  = -((wt_mu + wt_ai) * wt_I) + (wt_a2 * wt_C1) + (wt_bi * wt_O);

		wt_O  += dt * wt_dO;
		wt_C1 += dt * wt_dC1;
		wt_C2 += dt * wt_dC2;
		wt_C3 += dt * wt_dC3;
		wt_I  += dt * wt_dI;
	}
	
	
	//An addition so we can output the currents
	// Cellstate.set_new_general_cell_model_data(IKr);
	//////////////////

	/* IKs */
	double IKs = GKs_ABh * Gks * sxs * sxs * (Vm - Eks);  // with apicalbasal heterogeneity wit direct scaling factor, by haibo
	double Xs_INF = 1. / (1. + exp((-5. - Vm) / 14.));
	// double Axs = (1400. / (sqrt(1. + exp((5. - Vm) / 6))));
	// double Bxs = (1. / (1. + exp((Vm - 35.) / 15.)));
	double TAU_Xs = (1400. / (sqrt(1. + exp((5. - Vm) / 6)))) * (1. / (1. + exp((Vm - 35.) / 15.))) + 80;
	TAU_Xs = TauKs_ABh * TAU_Xs;
	sxs = Xs_INF - (Xs_INF - sxs) * exp(-dt / TAU_Xs);

	
	
	//An addition so we can output the currents
	Cellstate.set_new_general_cell_model_data(IKs);
	//////////////////

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

	sd = D_INF - (D_INF - sd) * exp(-dt / TAU_D);
	sf = F_INF - (F_INF - sf) * exp(-dt / TAU_F);
	sf2 = F2_INF - (F2_INF - sf2) * exp(-dt / TAU_F2);

	double FCaSS_INF = 0.6 / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 0.4;
	double TAU_FCaSS = 80. / (1 + (CaSS / 0.05) * (CaSS / 0.05)) + 2.;
	sfcass = FCaSS_INF - (FCaSS_INF - sfcass) * exp(-dt / TAU_FCaSS);

	
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

	mNaL = TMP_INF - (TMP_INF - mNaL) * exp(-dt / TAU_M);  // double tmL = TAU_M;
	TMP_INF = 1.0 / (1.0 + exp((Vm + 87.61) / 7.488));

	hNaL = TMP_INF - (TMP_INF - hNaL) * exp(-dt / 200.0);  // double thL = 200.0;
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
	if (m_celltype == 100)
		GK1 *= 1.2; // according to the ORd model supple, added by haibo
	double IK1 = GK1 * rec_iK1 * (Vm - Ek);
	// double sItot = IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL + Istim + IKATP;

	/* update concentrations */
	double kCaSR = TTCell_maxsr - ((TTCell_maxsr - TTCell_minsr) / (1 + (TTCell_EC / CaSR) * (TTCell_EC / CaSR)));
	double k1 = TTCell_k1_ / kCaSR;
	double k2 = TTCell_k2_ * kCaSR;
	double dRR = TTCell_k4 * (1 - sRR) - k2 * CaSS * sRR;
	sRR += dt * dRR;
	double sOO = k1 * CaSS * CaSS * sRR / (TTCell_k3 + k1 * CaSS * CaSS);

	double Irel = TTCell_Vrel * sOO * (CaSR - CaSS);
	double Ileak = TTCell_Vleak * (CaSR - cai);
	double Iup = TTCell_Vmaxup / (1. + ((TTCell_Kup * TTCell_Kup) / (cai * cai)));
	double Ixfer = TTCell_Vxfer * (CaSS - cai);

	{	double CaCSQN = TTCell_Bufsr * CaSR / (CaSR + TTCell_Kbufsr);
		double dCaSR = dt * (Iup - Irel - Ileak);
		double bjsr = TTCell_Bufsr - CaCSQN - dCaSR - CaSR + TTCell_Kbufsr;
		double cjsr = TTCell_Kbufsr * (CaCSQN + dCaSR + CaSR);
		CaSR = (sqrt(bjsr * bjsr + 4 * cjsr) - bjsr) / 2;
	}


	{	double CaSSBuf = TTCell_Bufss * CaSS / (CaSS + TTCell_Kbufss);
		double dCaSS = dt * (-Ixfer * (TTCell_Vc / TTCell_Vss) + Irel * (TTCell_Vsr / TTCell_Vss) + (-ICaL * TTCell_inversevssF2 * TTCell_CAPACITANCE));
		double bcss = TTCell_Bufss - CaSSBuf - dCaSS - CaSS + TTCell_Kbufss;
		double ccss = TTCell_Kbufss * (CaSSBuf + dCaSS + CaSS);
		CaSS = (sqrt(bcss * bcss + 4 * ccss) - bcss) / 2;
	}


	double CaBuf = TTCell_Bufc * cai / (cai + TTCell_Kbufc);
	double dCai = dt * ((-(IbCa + IpCa - 2 * INaCa) * TTCell_inverseVcF2 * TTCell_CAPACITANCE) - (Iup - Ileak) * (TTCell_Vsr / TTCell_Vc) + Ixfer);
	double bc = TTCell_Bufc - CaBuf - dCai - cai + TTCell_Kbufc;
	double cc = TTCell_Kbufc * (CaBuf + dCai + cai);
	cai = (sqrt(bc * bc + 4 * cc) - bc) / 2;


	double dNai = -(INa + INaL + IbNa + 3 * INaK + 3 * INaCa) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	Nai += dt * dNai;

	double dKi = -(/*Istim +*/ IK1 + Ito + IKr + IKs - 2 * INaK + IpK + IKATP) * TTCell_inverseVcF * TTCell_CAPACITANCE;
	Ki += dt * dKi;





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
    
    N        += dt*dN;
    XBprer   += dt*dXBprer;
    XBpostr  += dt*dXBpostr;
    SL       += dt*dSL;
    xXBpostr += dt*dxXBpostr;
    xXBprer  += dt*dxXBprer;
    TRPNCaL  += dt*dTRPNCaL;
    TRPNCaH  += dt*dTRPNCaH;
    intf     += dt*dintf;

    double dCai_feeback = dTropToT;

    //End solve force model



	new_state[0]   = cai    ;
	new_state[1]   = CaSR   ;
	new_state[2]   = CaSS   ;
	new_state[3]   = Nai    ;
	new_state[4]   = Ki     ;
	new_state[5]   = sm     ;
	new_state[6]   = sh     ;
	new_state[7]   = sj     ;
	new_state[8]   = sxr1   ;
	new_state[9]  = sxr2   ;
	new_state[10]  = sxs    ;
	new_state[11]  = ss     ;
	new_state[12]  = sr     ;
	new_state[13]  = sd     ;
	new_state[14]  = sf     ;
	new_state[15]  = sf2    ;
	new_state[16]  = sfcass ;
	new_state[17]  = sRR    ;
	new_state[18]  = mNaL   ;
	new_state[19]  = hNaL   ;
	new_state[20]  = wt_C3  ;
	new_state[21]  = wt_O   ;
	new_state[22]  = wt_C1  ;
	new_state[23]  = wt_C2  ;
	new_state[24]  = wt_I   ;

	new_state[25] = N;
	new_state[26] = XBprer;
	new_state[27] = XBpostr;
	new_state[28] = SL;
	new_state[29] = xXBpostr;
	new_state[30] = xXBprer;
	new_state[31] = TRPNCaL;
	new_state[32] = TRPNCaH;
	new_state[33] = intf;

	Cellstate.set_active_strain(-(SL - SLset) / SLset);

	// std::cout << IKr << "\t" << IKs<< "\t" << IK1<< "\t" << Ito<< "\t" << INa<< "\t" << IbNa << "\t" <<ICaL << "\t" <<IbCa<< "\t" << INaK<< "\t" << INaCa << "\t" <<IpCa << "\t" <<IpK << "\t" <<INaL /*+ Istim*/<< "\t" << IKATP;
	// std::exit(0);
	Cellstate.set_membrane_current(IKr + IKs + IK1 + Ito + INa + IbNa + ICaL + IbCa + INaK + INaCa + IpCa + IpK + INaL /*+ Istim*/ + IKATP);
}

} //end namespace