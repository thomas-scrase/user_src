#include "Ismail_TNNP_06.h"


namespace oomph{

IsmailTNNP06::IsmailTNNP06()
{
	Cell_Model_Name = "IsmailTNNP06";
	
	//The active strain is the first (and only) other variable the model generates
	active_strain_index = 0;


	Ko=5.4;
	Cao=2.0;
	Nao=140.0;
	Vc=0.016404;
	Vsr=0.001094;
	Vss=0.00005468;
	Bufc=0.2;
	Kbufc=0.001;
	Bufsr=10.;
	Kbufsr=0.3;
	Bufss=0.4;
	Kbufss=0.00025;
	Vmaxup=0.006375;
	Kup=0.00025;
	Vrel=0.102;//40.8;
	k1_=0.15;
	k2_=0.045;
	k3=0.060;
	k4=0.005;//0.000015;
	EC=1.5;
	maxsr=2.5;
	minsr=1.;
	Vleak=0.00036;
	Vxfer=0.0038;
	R=8314.472;
	F=96485.3415;
	T=310.0;
	RTONF=(R*T)/F;
	CAPACITANCE=0.185;
	Gkr=0.153;
	pKNa=0.03;

	// if (m_celltype == 102)
	// 	Gks=0.392;
	// else if (m_celltype == 100)
	// 	Gks=0.392;
	// else
	// 	Gks=0.098;
	
	// GK1=5.405;
	// //Parameters for Ito
	// if (m_celltype == 0)
	// 	Gto=0.294;
	// else if (m_celltype == 100)
	// 	Gto=0.073;
	// else
	// 	Gto=0.294;
	
	GNa=14.838;
	GbNa=0.00029;
	KmK=1.0;
	KmNa=40.0;
	knak=2.724;
	GCaL=0.00003980;  
	GbCa=0.000592;
	knaca=1000;
	KmNai=87.5;
	KmCa=1.38;
	ksat=0.1;
	n=0.35;
	GpCa=0.1238;
	KpCa=0.0005;
	GpK=0.0146;
	// svolt=-86.2;
	// Cai=0.00007;
	// CaSR=1.3;
	// CaSS=0.00007;
	// Nai=7.67;
	// Ki=138.3;
	
	// stimduration=1.;
	// stimstrength=-52;
	// tbegin=0;
	// tend=tbegin+stimduration;
	// counter=1;
	// dia=5000;
	// basicperiod=1000.0;
	// basicapd=274;
	// repeats=10;
	
	//
	inverseVcF2=1/(2*Vc*F);
	inverseVcF=1./(Vc*F);
	inversevssF2=1/(2*Vss*F);

	// SQT1
	// sqt1_C3 = 1.0;
	sqt1_a1 = 2.172; sqt1_b1 = 1.077; 
	
	// INaL
	/*mNaL = hNaL =*/ alpha_mNaL = beta_mNaL = mNaL_INF = hNaL_INF = TAU_mNaL = 0.0;
	TAU_hNaL = 600;

	//===========
	Qkon  = 1.5;
	Qkoff = 1.2;
	Qkn_p = 1.6;
	Qkp_n = 1.6;
	Qfapp = 6.25;
	Qgapp = 2.5;
	Qhf	  = 6.25;
	Qhb   = 6.25;
	Qgxb  = 6.25;
	//===========

	/* fix min and max length */
	SLmax=2.4;
	SLmin=1.4;
	
	// Define lengths for thick and think filaments (units = um)
	len_thick=1.65;
	len_hbare=0.1;
	len_thin=1.2;
	
	/* Handle Ca binding to troponin here */
	// Rates for troponin Ca binding (unit = (uM*ms)^-1)
	kon=50.e-3;
	// Low affinty offrate (unit = ms^-1)
	koffL=0.9*250.e-3;
	// High affinty offrate (unit = ms^-1)
	koffH=0.9*25.e-3;
	
	/* Set perm variable that control n to p transiton in Hill-like fashion, No Units. */
	perm50=0.5;  
	nperm=15.0;
	
	/* Define rates for RU   Units are s^-1. add e-3 would be ms^-1 */
	/* Define base rates for N-P and P-N transitions. Units are s^-1. add e-3 would be ms^-1 */ 	
	kn_p=500e-3;
	kp_n=50e-3;
	
	/*  Set the XB modifier to convert rat to rabbit*/
	//xbspec=0.2;
	
	/* Set the XB off rate Units are s^-1. add e-3 would be ms^-1 */
	fapp=500.e-3*0.2;
	gapp=70.0e-3*0.2; 
	/* # Parameter to control SL effect on gapp, set between 0 and 1 */
	gslmod=6.;
	
	/* Set rate for rates between pre-force and force states  Units are s^-1. add e-3 would be ms^-1 */
	hf=2000.0e-3*0.2;
	hb=400.0e-3*0.2; 
	hbmdc=0;
	hfmdc=5;
	
	/* Set rate for ATP consuming transition   Units are s^-1. add e-3 would be ms^-1 */
	gxb=70.e-3*0.2;
	//gxmdc=60.0;
	
	/* Add term for distortion dependence of STP using transition gxb */
	sigmap=8.0;
	sigman=1.0;
	
	/* Compute SL,  Set to 0 for isometric, to 1 for active contraction */
	/* #define contr 1.0 */
	
	// Set strain induced by head rotation (units = um)
	x_0=0.007;
	
	// Set scaling factor that set balance of two competing effects on mean strain -
	// SL motion and XB cycling 
	xPsi=2;
	
	// Use the kxb scaling term to compute a force value in units of mN/mm^2
	// as presented in many papers.  This value also corresponds to the maximal force.
	kxb=120.; 
	
	// This term sets the passive viscocity in Fmax /(um/ms)  
	// Here is Fmax is normalized to 1. 
	// Value is taken from de Tomble & ter Keurs, 1991 
	visc=0.003e3;
	
	// Parameter for muscle mass
	massf=0.00005/0.2;	//Rat
	//massf=0.00025/0.2;		//Rabbit
	
	// Note that passive force is specified in terms of maximal active force.
	// Passive force is comnpute is two ways depending on of cell is skinned or intact.#define singcell=1
	
	// Define resting sarcomere length where passive force = 0 
	SLrest=1.9;
	
	// Parmeter for initial sarcomere length (um)
	// For Fig 8, please use this so initial SL is resting SL
	//SLset=1.9; 
	
	// For Fig 9, please use this so initial SL is 2.2 um
	SLset=2.2; 
	
	// These apply to trabeculae and single cells - assume to represent titin
	PCon_t = 0.002; 
	PExp_t  = 10; 
	
	// These apply to trabeculae only - assume to represent collagen
	SLcol = 2.25; 
	SLcol  = 2.25; 
	PCon_col = 0.02; 
	PExp_col = 70; 
	
	// Code for defining a series elastic element
	// Units are in Normalized Force/um 
	KSE = 1;
	
	Trop_conc = 70.0; /* uM, add e-3 for mM */
	
	TropToT_dt = 0.0;

	/*The follow section of code implements is used to simulate ejection 
	 from canine whole heart into a Windkessel model of the circulatory system.
	 The complete description is available in:
	 De Tombe, PP and Little, WC. Inotropic effects of ejection 
	 are myocardial properties. Am J Physiol. 1994 Mar;266(3 Pt 2):H1202-13. */
	
	myo = true;

	//Parameters for Windkessel model
	// units for Cwind are ml/mmHg
	Cwind = 0.4; 
	// units for Rchar and Ra are ml*ms/mmHg
	Rchar = 0.2e3;
	// Ra is in range 1.5 to 12.0 ml*ms/mmHg
	Ra = 1.5e3;
	
	// Parameters for spherical heart model from de Tombe and Little
	// units for volumes are ml
	Vref = 80.0;
	Vwall = 100.0;
	
	// units for SLref are um
	SLref = 2.4;
	
	// Units of SigFmx is maximum pressure corresponding to Sigma_es 
	// units are g/cm^2 after conversion from mN/mm^2
	SigFmx = 1224.0;
	
	// Define inital condition for arterial pressures
	// units for ARTpset are mmHg
	ARTpset = 50.0;
	ARTpress = ARTpset;
	
	// Define inital condition for  left venticular volume in ml
	// Note that initial LVV and SL should match (e.g., 80 ml = 2.4 um and 0 ml = 1.6 um)
	// LVV = std::pow((SLset/SLref),3.)*(Vref+0.333*Vwall)-0.333*Vwall;


	//By default let's assume it is a mutant cell type
	m_mutant = true;

	m_epiMidRatio = 1.0;


	//Assign the names of the variables used by this model
	Names_Of_Cell_Variables =
	{
		"sm",
		"sh",
		"sj",

		//For HH type
		// "sxr1",
		// "sxr2",

		//For MC type
		"sqt1_O",
		"sqt1_C1",
		"sqt1_C2",
		"sqt1_C3",
		"sqt1_I",

		"sxs",
		"ss",
		"sr",
		"sd",
		"sf",
		"sf2",
		"sfcass",
		"sRR",
		// "sOO",
		"Cai",
		"CaSR",
		"CaSS",
		"mNaL",
		"hNaL",
		"Nai",
		"Ki",
		"N_NoXB",
		"N",
		"P_NoXB",
		"P",
		"XBprer",
		"XBpostr",
		"SL",
		"xXBpostr",
		"xXBprer",
		"TropCaL",
		"TropCaH",
		"intf0"
		//We only want the single cell variables, not the left ventricle pressure, etc
		// "ARTpset"
		// "LVV"
	};
	Names_Of_Other_Parameters =
	{

	};
	Names_Of_Other_Variables =
	{

	};
	Names_Of_Output_Data =
	{

	};

	FinalizeConstruction();
}

IsmailTNNP06::~IsmailTNNP06()
{

}

double IsmailTNNP06::return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
{
	if ( m_mutant ) { // SQT1
		if ( cell_type == 100 ) {
			switch(v){
				case sm_tt: return 0.00758532;
				case sh_tt: return 0.498893;
				case sj_tt: return 0.500138;
				//For HH type
				// case sxr1_tt: return 0.00127169;
				// case sxr2_tt: return 0.39808;
				//For MC type
				case sqt1_O_tt: return 0.0;
				case sqt1_C1_tt: return 0.0;
				case sqt1_C2_tt: return 0.0;
				case sqt1_C3_tt: return 1.0;
				case sqt1_I_tt: return 0.0;

				case sxs_tt: return 0.00542449;
				case ss_tt: return 0.999991;
				case sr_tt: return 7.96401e-08;
				case sd_tt: return 8.75153e-05;
				case sf_tt: return 0.970237;
				case sf2_tt: return 0.998578;
				case sfcass_tt: return 0.999914;
				case sRR_tt: return 0.980806;
				// case sOO_tt: return 2.17258e-07;
				case Cai_tt: return 0.000144459;
				case CaSR_tt: return 4.44326;
				case CaSS_tt: return 0.000319686;
				case mNaL_tt: return 0.0;
				case hNaL_tt: return 0.0;
				case Nai_tt: return 14.4595;
				case Ki_tt: return 128.885;
				case N_NoXB_tt: return 0.99;
				case N_tt: return 0.97;
				case P_NoXB_tt: return 0.01;
				case P_tt: return 0.01;
				case XBprer_tt: return 0.01;
				case XBpostr_tt: return 0.01;
				case SL_tt: return 2.17149;
				case xXBpostr_tt: return 0.007;
				case xXBprer_tt: return 0;
				case TropCaL_tt: return 0.0144725;
				case TropCaH_tt: return 0.232095;
				case intf0_tt: return 0.0;
			}
		}
		else if ( cell_type == 101 ) {
			switch(v){
				case sm_tt: return 0.00782063;
				case sh_tt: return 0.493186;
				case sj_tt: return 0.495549;
				//For HH type
				// case sxr1_tt: return 0.00150735;
				// case sxr2_tt: return 0.396571;
				//For MC type
				case sqt1_O_tt: return 0.0;
				case sqt1_C1_tt: return 0.0;
				case sqt1_C2_tt: return 0.0;
				case sqt1_C3_tt: return 1.0;
				case sqt1_I_tt: return 0.0;

				case sxs_tt: return 0.00549386;
				case ss_tt: return 0.999991;
				case sr_tt: return 8.16837e-08;
				case sd_tt: return 8.92981e-05;
				case sf_tt: return 0.965796;
				case sf2_tt: return 0.998548;
				case sfcass_tt: return 0.999906;
				case sRR_tt: return 0.980156;
				// case sOO_tt: return 2.25888e-07;
				case Cai_tt: return 0.00014691;
				case CaSR_tt: return 4.5881;
				case CaSS_tt: return 0.000324862;
				case mNaL_tt: return 0.0;
				case hNaL_tt: return 0.0;
				case Nai_tt: return 14.2504;
				case Ki_tt: return 128.612;
				case N_NoXB_tt: return 0.99;
				case N_tt: return 0.97;
				case P_NoXB_tt: return 0.01;
				case P_tt: return 0.01;
				case XBprer_tt: return 0.01;
				case XBpostr_tt: return 0.01;
				case SL_tt: return 2.17152;
				case xXBpostr_tt: return 0.007;
				case xXBprer_tt: return 0;
				case TropCaL_tt: return 0.0144725;
				case TropCaH_tt: return 0.232095;
				case intf0_tt: return 0.0;
			}
		}
		else if ( cell_type == 102 ) {
			switch(v){
				case sm_tt: return 0.00752173;
				case sh_tt: return 0.500809;
				case sj_tt: return 0.502428;
				//For HH type
				// case sxr1_tt: return 0.00121345;
				// case sxr2_tt: return 0.398501;
				//For MC type
				case sqt1_O_tt: return 0.0;
				case sqt1_C1_tt: return 0.0;
				case sqt1_C2_tt: return 0.0;
				case sqt1_C3_tt: return 1.0;
				case sqt1_I_tt: return 0.0;

				case sxs_tt: return 0.00540245;
				case ss_tt: return 0.547071;
				case sr_tt: return 7.90571e-08;
				case sd_tt: return 8.70287e-05;
				case sf_tt: return 0.971538;
				case sf2_tt: return 0.998588;
				case sfcass_tt: return 0.999914;
				case sRR_tt: return 0.980858;
				// case sOO_tt: return 2.12427e-07;
				case Cai_tt: return 0.000142578;
				case CaSR_tt: return 4.35494;
				case CaSS_tt: return 0.000316873;
				case mNaL_tt: return 0.0;
				case hNaL_tt: return 0.0;
				case Nai_tt: return 14.3286;
				case Ki_tt: return 129.24;
				case N_NoXB_tt: return 0.99;
				case N_tt: return 0.97;
				case P_NoXB_tt: return 0.01;
				case P_tt: return 0.01;
				case XBprer_tt: return 0.01;
				case XBpostr_tt: return 0.01;
				case SL_tt: return 2.17157;
				case xXBpostr_tt: return 0.007;
				case xXBprer_tt: return 0;
				case TropCaL_tt: return 0.0144725;
				case TropCaH_tt: return 0.232095;
				case intf0_tt: return 0.0;
			}
		}
	}
	else { // WT
		if ( cell_type == 100 ) {
			switch(v){
				case sm_tt: return 0.00808475;
				case sh_tt: return 0.486634;
				case sj_tt: return 0.489046;
				//For HH type
				// case sxr1_tt: return 0.00260931;
				// case sxr2_tt: return 0.394936;
				//For MC type
				case sqt1_O_tt: return 0.0;
				case sqt1_C1_tt: return 0.0;
				case sqt1_C2_tt: return 0.0;
				case sqt1_C3_tt: return 1.0;
				case sqt1_I_tt: return 0.0;

				case sxs_tt: return 0.00567107;
				case ss_tt: return 0.99999;
				case sr_tt: return 8.39168e-08;
				case sd_tt: return 9.12767e-05;
				case sf_tt: return 0.948576;
				case sf2_tt: return 0.998513;
				case sfcass_tt: return 0.999874;
				case sRR_tt: return 0.978115;
				// case sOO_tt: return 2.47117e-07;
				case Cai_tt: return 0.000153134;
				case CaSR_tt: return 4.97334;
				case CaSS_tt: return 0.000337189;
				case mNaL_tt: return 0.0;
				case hNaL_tt: return 0.0;
				case Nai_tt: return 13.3024;
				case Ki_tt: return 129.765;
				case N_NoXB_tt: return 0.99;
				case N_tt: return 0.97;
				case P_NoXB_tt: return 0.01;
				case P_tt: return 0.01;
				case XBprer_tt: return 0.01;
				case XBpostr_tt: return 0.01;
				case SL_tt: return 2.17143;
				case xXBpostr_tt: return 0.007;
				case xXBprer_tt: return 0;
				case TropCaL_tt: return 0.0144725;
				case TropCaH_tt: return 0.232095;
				case intf0_tt: return 0.0;
			}
		}
		else if ( cell_type == 101 ) {
			switch(v){
				case sm_tt: return 0.00852496;
				case sh_tt: return 0.477014;
				case sj_tt: return 0.474622;
				//For HH type
				// case sxr1_tt: return 0.00741203;
				// case sxr2_tt: return 0.392327;
				//For MC type
				case sqt1_O_tt: return 0.0;
				case sqt1_C1_tt: return 0.0;
				case sqt1_C2_tt: return 0.0;
				case sqt1_C3_tt: return 1.0;
				case sqt1_I_tt: return 0.0;

				case sxs_tt: return 0.00646137;
				case ss_tt: return 0.99999;
				case sr_tt: return 8.76076e-08;
				case sd_tt: return 9.45337e-05;
				case sf_tt: return 0.90623;
				case sf2_tt: return 0.99846;
				case sfcass_tt: return 0.999812;
				case sRR_tt: return 0.974792;
				// case sOO_tt: return 2.80325e-07;
				case Cai_tt: return 0.000161959;
				case CaSR_tt: return 5.46463;
				case CaSS_tt: return 0.000356533;
				case mNaL_tt: return 0.0;
				case hNaL_tt: return 0.0;
				case Nai_tt: return 11.7376;
				case Ki_tt: return 131.076;
				case N_NoXB_tt: return 0.99;
				case N_tt: return 0.97;
				case P_NoXB_tt: return 0.01;
				case P_tt: return 0.01;
				case XBprer_tt: return 0.01;
				case XBpostr_tt: return 0.01;
				case SL_tt: return 2.17115;
				case xXBpostr_tt: return 0.007;
				case xXBprer_tt: return 0;
				case TropCaL_tt: return 0.0144725;
				case TropCaH_tt: return 0.232095;
				case intf0_tt: return 0.0;
			}
		}
		else if ( cell_type == 102 ) {
			switch(v){
				case sm_tt: return 0.0081098;
				case sh_tt: return 0.48616;
				case sj_tt: return 0.48855;
				//For HH type
				// case sxr1_tt: return 0.00264338;
				// case sxr2_tt: return 0.394787;
				//For MC type
				case sqt1_O_tt: return 0.0;
				case sqt1_C1_tt: return 0.0;
				case sqt1_C2_tt: return 0.0;
				case sqt1_C3_tt: return 1.0;
				case sqt1_I_tt: return 0.0;

				case sxs_tt: return 0.00568398;
				case ss_tt: return 0.501349;
				case sr_tt: return 8.41039e-08;
				case sd_tt: return 9.14621e-05;
				case sf_tt: return 0.948521;
				case sf2_tt: return 0.998511;
				case sfcass_tt: return 0.999868;
				case sRR_tt: return 0.977888;
				// case sOO_tt: return 2.44788e-07;
				case Cai_tt: return 0.000151786;
				case CaSR_tt: return 4.92014;
				case CaSS_tt: return 0.000336006;
				case mNaL_tt: return 0.0;
				case hNaL_tt: return 0.0;
				case Nai_tt: return 13.0618;
				case Ki_tt: return 129.98;
				case N_NoXB_tt: return 0.99;
				case N_tt: return 0.97;
				case P_NoXB_tt: return 0.01;
				case P_tt: return 0.01;
				case XBprer_tt: return 0.01;
				case XBpostr_tt: return 0.01;
				case SL_tt: return 2.17147;
				case xXBpostr_tt: return 0.007;
				case xXBprer_tt: return 0;
				case TropCaL_tt: return 0.0144725;
				case TropCaH_tt: return 0.232095;
				case intf0_tt: return 0.0;
			}
		}
	}
}

double IsmailTNNP06::return_initial_membrane_potential(const unsigned &cell_type)
{
	// SQT1
	if ( m_mutant ) {
		if ( cell_type == 100 ) {
			return -78.0773;
		}
		else if ( cell_type == 101 ) {
			return -77.9261;
		}
		else if ( cell_type == 102 ) {
			return -78.119;
		}
	}
	// WT
	else {
		if ( cell_type == 100 ) {
			return -77.7615;
		}
		else if ( cell_type == 101 ) {
			return -77.4981;
		}
		else if ( cell_type == 102 ) {
			return -77.7461;
		}
	}
}

// void IsmailTNNP06::Calculate_Derivatives(const double &Vm,
// 							const Vector<double> &CellVariables,
// 							const double &t,
// 							const unsigned &cell_type,
// 							const double &Istim,
// 							const Vector<double> &Other_Parameters,
// 							const Vector<double> &Other_Variables,
// 							Vector<double> &Variable_Derivatives,
// 							double &Iion)
void IsmailTNNP06::Calculate_Derivatives(const Boost_State_Type &Variables,
									const double &t,
									const unsigned &cell_type,
									const double &Istim,
									const Vector<double> &Other_Parameters,
									const Vector<double> &Other_Variables,
									Vector<double> &Variable_Derivatives,
									double &Iion)
{

	const double sm			= Variables[sm_tt];
	const double sh			= Variables[sh_tt];
	const double sj			= Variables[sj_tt];

	//For H-H type
	// const double sxr1		= Variables[sxr1_tt];
	// const double sxr2		= Variables[sxr2_tt];

	//For MC type
	const double sqt1_O		= Variables[sqt1_O_tt];
	const double sqt1_C1	= Variables[sqt1_C1_tt];
	const double sqt1_C2	= Variables[sqt1_C2_tt];
	const double sqt1_C3	= Variables[sqt1_C3_tt];
	const double sqt1_I		= Variables[sqt1_I_tt];

	const double sxs		= Variables[sxs_tt];
	const double ss			= Variables[ss_tt];
	const double sr			= Variables[sr_tt];
	const double sd			= Variables[sd_tt];
	const double sf			= Variables[sf_tt];
	const double sf2		= Variables[sf2_tt];
	const double sfcass		= Variables[sfcass_tt];
	const double sRR		= Variables[sRR_tt];
	// const double sOO		= Variables[sOO_tt];
	const double Cai		= Variables[Cai_tt];
	const double CaSR		= Variables[CaSR_tt];
	const double CaSS		= Variables[CaSS_tt];
	const double mNaL		= Variables[mNaL_tt];
	const double hNaL		= Variables[hNaL_tt];
	const double Nai		= Variables[Nai_tt];
	const double Ki			= Variables[Ki_tt];
	const double N_NoXB		= Variables[N_NoXB_tt];
	const double N			= Variables[N_tt];
	const double P_NoXB		= Variables[P_NoXB_tt];
	const double P			= Variables[P_tt];
	const double XBprer		= Variables[XBprer_tt];
	const double XBpostr	= Variables[XBpostr_tt];
	const double SL			= Variables[SL_tt];
	const double xXBpostr	= Variables[xXBpostr_tt];
	const double xXBprer	= Variables[xXBprer_tt];
	const double TropCaL	= Variables[TropCaL_tt];
	const double TropCaH	= Variables[TropCaH_tt];
	const double intf0 		= Variables[intf0_tt];

	const double Vm = Variables[Num_Cell_Vars];


	if (cell_type == 100)
	{
		Gks=0.392;
	}
	else if (cell_type == 102)
	{
		Gks=0.392;
	}
	else
	{
		Gks=0.098;
	}
	
	GK1=5.405;
	//Parameters for Ito
	if (cell_type == 100)
	{
		Gto=0.294;
	}
	else if (cell_type == 102)
	{
		Gto=0.073;
	}
	else
	{
		Gto=0.294;
	}

	//ORiginally cvode solved myocardium stuff
	// flag = CVode(cvode_mem, TOUT, y, &time, CV_NORMAL);	

	
	//N_NoXB = Cai * 1000;
	
	/* Handle Ca binding to troponin here */
	//  Adjust for temperature
	konT	= kon	*std::pow(Qkon,(T-273.0-37.0)/10.);
	koffLT	= koffL	*std::pow(Qkoff,(T-273.0-37.0)/10.);
	koffHT	= koffH	*std::pow(Qkoff,(T-273.0-37.0)/10.);
	
	/* Compute derived regulatory unit rates, Ca in in nM to prevent Ca overflow (max ~ 10 in auto), TropCaL and TropCaH are normailized so no Units. */
  	Variable_Derivatives[TropCaL_tt] = konT * Cai*1000 * (1.0 - TropCaL) - koffLT * TropCaL;                                              
  	Variable_Derivatives[TropCaH_tt] = konT * Cai*1000 * (1.0 - TropCaH) - koffHT * TropCaH; 
	/* Compute rates for N to P transitions */
	// Compute z-line end of single overlap region
	sovr_ze = std::min(len_thick/2.,SL/2.);
	
	// Compute centerline of end of single overlap region
    sovr_cle = std::max(SL/2.-(SL-len_thin),len_hbare/2.);
	
	// Compute length of the single overlap 
	len_sovr = sovr_ze-sovr_cle;
	
	// Compute overlap fraction for thick filament
	SOVFThick = len_sovr*2./(len_thick-len_hbare);
	
	// Compute overlap fraction for thin filament
	SOVFThin = len_sovr/len_thin;
	
	/* Compute combined Ca binding to high (w/XB) and low (w/o XB)  No Units. */
	perm = (1.0 - SOVFThin) * TropCaL + SOVFThin * TropCaH;
	
	/* Set perm variable that control n to p transiton in Hill-like fashion, No Units. */
	permtot = std::pow((1.0 / (1.0 + std::pow((perm50 / perm), nperm))), 0.5);
 	if (100.0 < 1.0/permtot)
		inprmt=100.0;
  	else 
		inprmt = 1.0/permtot;
	
	/* Adjust for Ca activation level and temp. */
	kn_pT = kn_p*permtot*std::pow(Qkn_p,(T-273.0-37.0)/10.);
	kp_nT = kp_n*inprmt*std::pow(Qkp_n,(T-273.0-37.0)/10.);
	
	/* Compute fapp, the attachment rate from weak to strong, pre-rotated state */
    fappT = fapp*std::pow(Qfapp,(T-273.0-37.0)/10.);
	
	/* Compute gapp, the detachment rate from strong, pre-rotated to the weakly-bound state */
	/* Compute SL modifier for gapp to increase rate st shorter SL */
  	gapslmd = 1.0 + (1.0 - SOVFThick)*gslmod;
    gappT = gapp*gapslmd*std::pow(Qgapp,(T-273.0-37.0)/10.);
	
	/* Set rate for rates between pre-force and force states*/
	/* Compute modifiers based on mean strain of states. */
  	hfmd=std::exp(-sign(xXBprer)*hfmdc*((xXBprer/x_0)*(xXBprer/x_0)));
	// oomph_info << "xXBprer " << xXBprer << std::endl;
	// oomph_info << "hfmdc " << hfmdc << std::endl;
	// oomph_info << "x_0 " << x_0 << std::endl;
  	hbmd=std::exp(sign((xXBpostr-x_0))*hbmdc*(((xXBpostr-x_0)/x_0)*((xXBpostr-x_0)/x_0)));
	
	/* Combine modifiers of hf and hb */
	hfT=hf*hfmd*std::pow(Qhf,(T-273.0-37.0)/10.);
	// oomph_info << "hf " << hf << std::endl;
	// oomph_info << "hfmd " << hfmd << std::endl;
	// oomph_info << "Qhf " << Qhf << std::endl;
	// oomph_info << "T " << T << std::endl;
	hbT=hb*hbmd*std::pow(Qhb,(T-273.0-37.0)/10.);
	
	/* Set rate for rates gxb, the ATP using XB transition */
	/* Add term for distortion dependence of gxb */
	if(x_0>xXBpostr)
	{
		gxbmd = std::exp(sigmap*((x_0-xXBpostr)/x_0)*((x_0-xXBpostr)/x_0))+(1.0-std::exp(sigman*(((xXBpostr-x_0)/x_0)*(xXBpostr-x_0)/x_0)));
	}
	else
	{
		gxbmd = 1.0;
	}
    // gxbmd = heav(x_0-xXBpostr)*std::exp(sigmap*((x_0-xXBpostr)/x_0)*((x_0-xXBpostr)/x_0))+(1.-heav(x_0-xXBpostr)*std::exp(sigman*(((xXBpostr-x_0)/x_0)*(xXBpostr-x_0)/x_0)));
	// oomph_info << "x_0 " << x_0 << std::endl;
	// oomph_info << "xXBpostr " << xXBpostr << std::endl;
	// oomph_info << "sigmap " << sigmap << std::endl;
	// oomph_info << "sigman " << sigman << std::endl;

    gxbT = gxb*gxbmd*std::pow(Qgxb,(T-273.0-37.0)/10.);

	// oomph_info << "gxb " << gxb << std::endl;
	// oomph_info << "gxbmd " << gxbmd << std::endl;
	// oomph_info << "Qgxb " << Qgxb << std::endl;
	// oomph_info << "T " << T << std::endl;

	
	/* update all RUs */
	Variable_Derivatives[N_NoXB_tt] = -kn_pT * N_NoXB + kp_nT * P_NoXB;
	Variable_Derivatives[P_NoXB_tt] = -kp_nT * P_NoXB + kn_pT * N_NoXB;
	
	Variable_Derivatives[N_tt] = -kn_pT * N + kp_nT * P; 
	Variable_Derivatives[P_tt] = -kp_nT * P + kn_pT * N - fappT * P + gappT * XBprer + gxbT * XBpostr; 
	Variable_Derivatives[XBprer_tt] = fappT * P - gappT * XBprer - hfT * XBprer + hbT * XBpostr;
	Variable_Derivatives[XBpostr_tt] = hfT * XBprer - hbT * XBpostr - gxbT * XBpostr;
	
	/* compute steady-state fractions in XBprer and XBpostr using King-Altman rule */
	SSXBprer = (hb*fapp+gxb*fapp)/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
	SSXBpostr = fapp*hf/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
	/* compute normalization for scaling active and passive force */
	Fnordv = kxb*x_0*SSXBpostr;
	
	/* compute force here */
	force = kxb*SOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer);
	
	/* Compute SL, Set to 0 for isometric, to 1 for active contraction */
	/* Note that passive force is specified in terms of maximal active force */
	/* read passive force from table file specified in terms of maximal active force */
	/* compyet derivative of SL */
	
	// Note that passive force is specified in terms of maximal active force.
	// Passive force is computed in two ways depending on if cell is skinned or intact.
	ppforce=sign((SL-SLrest))*PCon_t*(std::exp(PExp_t*std::abs((SL-SLrest)))-1);
    // if (singlecell == false) 
	ppforce += std::max(0.0,(PCon_col*std::exp(PExp_col*(SL-SLcol))));
	
	// Compute afterload either as a contant or from the series elastic element
 	afterload = 0.0;
	
 //    if (SEon == true) {
	// 	afterload = KSE*(SLset-SL);
	// 	if ((SEon_LengthClamp==true)&&(index==(pulse_number-1))) {
	// 		afterload=5000.0*(SLset-SL);
	// 	}
 //    } 
	// else {
	// 	afterload=0.0; 
 //    }
	
	//  Compute the integral of forces 
	Variable_Derivatives[intf0_tt]=(-ppforce+PreloadF+(-force/Fnordv)+afterload);
	
	// Compute the derivative of SL for an isolated muscle
  	dSLisolated = ((intf0 + (SLset-SL)*visc)/massf )*heav(SL-SLmin)*heav(SLmax-SL);
 	// if (contraction==false) 
		// dSLisolated = 0.;
	

	//We're not using the whole heart model
	

	// //========>>>>>>>>
	// /*  The follow section of code implements is used to simulate ejection 
	//  from canine whole heart into a Windkessel model of the circulatory system.
	//  The complete description is available in:
	//  De Tombe, PP and Little, WC. Inotropic effects of ejection 
	//  are myocardial properties. Am J Physiol. 1994 Mar;266(3 Pt 2):H1202-13. */
	
	// // Betavl is unitless
	// Betavl=M_PI*std::pow(0.75*M_PI,0.6667);
	
 //    // Compute A_0 which is the cross-sectional area of the ventricular wall
 //    // A_0 has unit of area (volume^0.6667) in cm^2
	// A_0 = Betavl*( std::pow(Vref + Vwall,0.6667) - std::pow(Vref,0.6667) );
	
 //    // Compute SigmaP from normalized force
 // 	SigPrs=SigFmx*(ppforce + force/Fnordv);
	
 //    // Compute left venticular pressure in mmHg
 // 	LVpress=SigPrs*A_0/(1.36*Betavl*std::pow(Ith(y,14),0.667));
	
	// // Compute Flows - assumes aortic valve so flow in unidirectional
	// // Units are l/s
	// iflow=heav(LVpress-Ith(y,13))*((LVpress-Ith(y,13))/Rchar);
	// // Allow active ejection on last beat only 
	// if (index<(pulse_number-1))
	// 	iflow=0.0;
	// if (Isovolume==true) 
	// 	iflow=0.0;
	
 //    // Compute the derivative of arterial pressure
	// Variable_Derivatives[ARTpset_tt]=(iflow/Cwind - (Ith(y,13)-ARTpset)/(Cwind*Ra));
 //    if (WholeHeart==false) 
	// 	Variable_Derivatives[ARTpset_tt]=0.0;
    
	// // Compute the  derivative of the change in Left Ventricular Volume
	// // Do not allow volume less than 0
 //    Variable_Derivatives[LVV_tt]=-(iflow)*heav(Ith(y,14));
 //    if (WholeHeart==false) 
	// 	Variable_Derivatives[LVV_tt]=0.0;
	
	// // Compute the  derivative of the change in SL
 //    dSLwholeheart=SLref*(1/std::pow(Vref+0.333*Vwall,0.333))*(Variable_Derivatives[LVV_tt]/3)/std::pow(Ith(y,14)+0.333*Vwall,0.667);
	

	// Compute the derivative of SL based using cell/tissue or simulated whole heart model
	// if (WholeHeart==true) {
	// 	Variable_Derivatives[SL_tt] = dSLwholeheart;
	// } else {
		Variable_Derivatives[SL_tt] = dSLisolated;
	// }
	
	//========>>>>>>>>
	
	
	if(myo == false)
	{
		Variable_Derivatives[N_NoXB_tt]=0.0;
		Variable_Derivatives[N_tt]=0.0;
		Variable_Derivatives[P_NoXB_tt]=0.0;
		Variable_Derivatives[P_tt]=0.0;
		Variable_Derivatives[XBprer_tt]=0.0;
		Variable_Derivatives[XBpostr_tt]=0.0;  
		Variable_Derivatives[SL_tt]=0.0; 
		Variable_Derivatives[TropCaL_tt]=0.0;
		Variable_Derivatives[TropCaH_tt]=0.0;
		Variable_Derivatives[intf0_tt]=0.0; 
		// Variable_Derivatives[ARTpset_tt]=0.0;
		// Variable_Derivatives[LVV_tt]=0.0;
	}
	
	dSL=Variable_Derivatives[SL_tt];
	
	// Compute the duty fractions using King-Alman Rule
	// Compute for states XBpref and XBf
	dtyf_prer=(hbT*fappT+gxbT*fappT)/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
	dtyf_postr=fappT*hfT/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
	
	// oomph_info << "fappT " << fappT << std::endl;
	// oomph_info << "hfT " << hfT << std::endl;
	// oomph_info << "gxbT " << gxbT << std::endl;
	// oomph_info << "gappT " << gappT << std::endl;
	// oomph_info << "hbT " << hbT << std::endl;

	// Compute mean strain by consider two competing effects - SL motion and XB cycling 
	if(myo == true)
	{
		Variable_Derivatives[xXBprer_tt] = dSL/2. + xPsi*(1/dtyf_prer)*(-xXBprer*fappT+(xXBpostr-x_0-xXBprer)*hbT);
		Variable_Derivatives[xXBpostr_tt] = dSL/2. + xPsi*(1/dtyf_postr)*(x_0+xXBprer-xXBpostr)*hfT;
		// if(!std::isfinite(Variable_Derivatives[xXBpostr_tt]))
		// {
		// 	Variable_Derivatives[xXBpostr_tt] = dSL/2. + xPsi*(x_0+xXBprer-xXBpostr)*(gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT)/fappT;
		// }
	}
	else
	{
		Variable_Derivatives[xXBprer_tt]=0.0;
		Variable_Derivatives[xXBpostr_tt]=0.0;
	}
	
	// oomph_info << "dSL " << dSL << std::endl;
	// oomph_info << "xPsi " << xPsi << std::endl;
	// oomph_info << "dtyf_postr " << dtyf_postr << std::endl;
	// oomph_info << "x_0 " << x_0 << std::endl;
	// oomph_info << "xXBprer " << xXBprer << std::endl;
	// oomph_info << "xXBpostr " << xXBpostr << std::endl;
	// oomph_info << "hfT " << hfT << std::endl;

	// Compute derivative of z-line end of single overlap region
	sovr_ze_dt=0.0;
	if(len_thick/2.>SL/2.) 
	{
		sovr_ze_dt=-dSL/2.;
	}
	
	// Compute derivative of centerline of end of single overlap region
    sovr_cle_dt=0.0;

    if (SL/2.0-(SL-len_thin)>len_hbare/2.0)
    {
		sovr_cle_dt=-dSL/2.0;
	}
	
	// Compute the derivative of the length of the single overlap 
	len_sovr_dt=sovr_ze_dt-sovr_cle_dt;
 	
	// Compute the derivative of the overlap fraction for thin filament

	SOVFThin_dt = len_sovr_dt/len_thin;
	
	// Compute the derivative of the overlap fraction for thick filament
	SOVFThick_dt = len_sovr_dt*2./(len_thick-len_hbare);
	
	/*Compute fraction of strongly bound for compute effective binding of Ca */
	//double FrSBXB=((XBpostr+XBprer)/(SSXBpostr + SSXBprer))*SOVFThick;
	//double FrSBXB_dt=((Variable_Derivatives[XBpostr_tt]+Variable_Derivatives[XBprer_tt])/(SSXBpostr + SSXBprer))*SOVFThick + ((XBpostr+XBprer)/(SSXBpostr + SSXBprer))*SOVFThick_dt;
	FrSBXB=((XBpostr+XBprer)/(SSXBpostr + SSXBprer));

	// oomph_info << "Variable_Derivatives[XBpostr_tt] " << Variable_Derivatives[XBpostr_tt] << std::endl;
	// oomph_info << "Variable_Derivatives[XBprer_tt] " << Variable_Derivatives[XBprer_tt] << std::endl;
	// oomph_info << "SSXBpostr " << SSXBpostr << std::endl;
	// oomph_info << "SSXBprer " << SSXBprer << std::endl;

	FrSBXB_dt=((Variable_Derivatives[XBpostr_tt]+Variable_Derivatives[XBprer_tt])/(SSXBpostr + SSXBprer));
	
	TropToT=Trop_conc*((1.0-SOVFThin)*TropCaL + SOVFThin*(FrSBXB*TropCaH+(1.0-FrSBXB)*TropCaL));
	
	// oomph_info << "TropCaL "<< TropCaL << std::endl;
	// oomph_info << "SOVFThin "<< SOVFThin << std::endl;
	// oomph_info << "Variable_Derivatives[TropCaL_tt] "<< Variable_Derivatives[TropCaL_tt] << std::endl;
	// oomph_info << "SOVFThin_dt "<< SOVFThin_dt << std::endl;
	// oomph_info << "FrSBXB "<< FrSBXB << std::endl;
	// oomph_info << "TropCaH "<< TropCaH << std::endl;
	// oomph_info << "FrSBXB "<< FrSBXB << std::endl;
	// oomph_info << "SOVFThin "<< SOVFThin << std::endl;
	// oomph_info << "FrSBXB_dt "<< FrSBXB_dt << std::endl;
	// oomph_info << "TropCaH "<< TropCaH << std::endl;
	// oomph_info << "FrSBXB "<< FrSBXB << std::endl;
	// oomph_info << "Variable_Derivatives[TropCaH_tt] "<< Variable_Derivatives[TropCaH_tt] << std::endl;
	// oomph_info << "FrSBXB_dt "<< FrSBXB_dt << std::endl;
	// oomph_info << "FrSBXB "<< FrSBXB << std::endl;
	// oomph_info << "Variable_Derivatives[TropCaL_tt] "<< Variable_Derivatives[TropCaL_tt] << std::endl;

	TropToT_dt=Trop_conc*(-1.0*SOVFThin_dt*TropCaL+(1.0-SOVFThin)*Variable_Derivatives[TropCaL_tt] + SOVFThin_dt*(FrSBXB*TropCaH+(1.0-FrSBXB)*TropCaL) + SOVFThin*(FrSBXB_dt*TropCaH+FrSBXB*Variable_Derivatives[TropCaH_tt]-FrSBXB_dt*TropCaL+(1.0-FrSBXB)*Variable_Derivatives[TropCaL_tt]));
	/*
	 if(myo == true)
	 Ith(ydot,26)=TropToT_dt; 
	 else 
	 Ith(ydot,89)=0.0; 
	 */
	//	/* Save non-state variable data */
	//	temp1 = force/Fnordv;         /* normalized force */ 
	//	temp1 = 120.0*force/Fnordv;   /* in units of mN/(mm^2) */
	//active_force = force/Fnordv;
	//active_force = force;
	active_force = 10*force/Fnordv;   /* in units of mN/(mm^2) */
	
	//	temp2 = SOVFThin;
	//	temp3 = SOVFThick;
	//	temp4 = SOVFThin_dt;
	//	temp5 = SOVFThick_dt;
	//	temp6 = perm;
	//	temp7 = permtot;
	//	temp8 = ppforce;
	//	temp9 = PreloadF;
	//	temp10 = FrSBXB; 
	//	temp11 = inprmt;
	//	temp12 = kn_pT;
	//	temp13 = kp_nT;
	//	temp14 = fappT;
	//	temp15 = gappT;
	//	temp16 = hfT;
	//	temp17 = hbT;
	//	temp18 = gxbT;
	//	temp19 = sovr_ze;
	//	temp20 = sovr_cle;
	//	temp21 = sovr_ze_dt;
	//	temp22 = sovr_cle_dt;
	//	temp23 = len_sovr_dt;
	//	temp24 = TropToT;
	//	temp25 = TropToT_dt;
	//	temp26 = dtyf_prer;
	//	temp27 = dtyf_postr;
	/****************************** myofilament ended *****************/
	//std::cout << "In feval..." << std::endl;










	Ek=RTONF*(log((Ko/Ki)));
	Ena=RTONF*(log((Nao/Nai)));
	Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
	
	
	Eca=0.5*RTONF*(log((Cao/Cai)));




	Ak1=0.1/(1.+exp(0.06*(Vm-Ek-200)));
	Bk1=(3.*exp(0.0002*(Vm-Ek+100))+exp(0.1*(Vm-Ek-10)))/(1.+exp(-0.5*(Vm-Ek)));
	rec_iK1=Ak1/(Ak1+Bk1);
	rec_iNaK=(1./(1.+0.1245*exp(-0.1*Vm*F/(R*T))+0.0353*exp(-Vm*F/(R*T))));
	rec_ipK=1./(1.+exp((25-Vm)/5.98));




	INa=GNa*sm*sm*sm*sh*sj*(Vm-Ena);



	ICaL=GCaL*sd*sf*sf2*sfcass*4.0*(Vm-15.0)*(F*F/(R*T))*(0.25*exp(2.0*(Vm-15.0)*F/(R*T))*CaSS-Cao)/(exp(2.0*(Vm-15.0)*F/(R*T))-1.0);




	Ito=Gto*sr*ss*(Vm-Ek);




	// IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(Vm-Ek);



	epi_factor      = 1.8	*	m_epiMidRatio;
	endo_factor     = 1.8; 
	mcell_factor    = 1.8; 
	
	if (m_mutant) {
		sqt1_a1 = 			2.172;
		sqt1_b1 = 0.5*		1.077;
		sqt1_a2 = 0.3*		0.00655*exp(0.05547153*(Vm-36. + 15.));
		sqt1_a =  			0.00555*exp(0.05547153*(Vm-12. + 15.));
		sqt1_b =  			0.002357*exp(-0.036588*(Vm));
		sqt1_b2 = 0.00077*	0.0029357*exp(1.3*	3.3*-0.02158*(Vm));
		sqt1_ai = 			0.439*exp(-0.02352*(Vm+25. + 15.))*(4.5/Ko);
		sqt1_bi = 0.025*		0.656*exp(0.000942*(Vm + 15.))*((pow((4.5/Ko),0.3)));
		sqt1_mu = 			(sqt1_ai*sqt1_b2*sqt1_a2)/(sqt1_a2*sqt1_bi);
	}
	else {
		sqt1_a1 = 			2.172;
		sqt1_b1 = 			1.077;
		sqt1_a2 = 			0.00655   * exp(0.5*			0.05547153*(Vm-36.));
		sqt1_a  =  			0.00555   * exp(				0.05547153*(Vm-12.));
		sqt1_b  =  			0.002357  * exp(				-0.036588*(Vm));
		sqt1_b2 = 0.65*		0.0029357 * exp(0.69*			-0.02158*(Vm));
		sqt1_ai = 0.11*		0.439     * exp(1.7*			-0.02352*(Vm+25.))*(4.5/Ko);
		sqt1_bi = 0.4*		0.656     * exp(				0.000942*(Vm))*((pow((4.5/Ko),0.3)));
		sqt1_mu = 			(sqt1_ai*sqt1_b2)/sqt1_bi;			
	}

	if (cell_type == 102)
		Gkr = 0.0135*pow(Ko, 0.59) * endo_factor;
	else if(cell_type == 101)
		Gkr = 0.0135*pow(Ko, 0.59) * mcell_factor;
	else			
		Gkr = 0.0135*pow(Ko,0.59) * epi_factor; 
	
    const double sqt1_dC3 = (sqt1_b * sqt1_C2)-(sqt1_a * sqt1_C3);
    const double sqt1_dC2 = -((sqt1_b + sqt1_a1) * sqt1_C2)+(sqt1_a * sqt1_C3)+(sqt1_b1 * sqt1_C1);
    const double sqt1_dC1 = -((sqt1_b1 + sqt1_a2 + sqt1_a2) * sqt1_C1) + (sqt1_a1 * sqt1_C2) + (sqt1_b2 * sqt1_O) + (sqt1_mu * sqt1_I);
    const double sqt1_dO  =  -((sqt1_b2 + sqt1_bi) * sqt1_O) + (sqt1_a2 * sqt1_C1) + (sqt1_ai * sqt1_I);
    const double sqt1_dI  = -((sqt1_mu + sqt1_ai) * sqt1_I) + (sqt1_a2 * sqt1_C1) + (sqt1_bi * sqt1_O);
	
    Variable_Derivatives[sqt1_O_tt]  = sqt1_dO;
    Variable_Derivatives[sqt1_C1_tt] = sqt1_dC1;
    Variable_Derivatives[sqt1_C2_tt] = sqt1_dC2;
    Variable_Derivatives[sqt1_C3_tt] = sqt1_dC3;
    Variable_Derivatives[sqt1_I_tt]  = sqt1_dI;

	
	IKr = Gkr*sqt1_O*(Vm-Ek);


	IKs=Gks*sxs*sxs*(Vm-Eks);



	IK1=GK1*rec_iK1*(Vm-Ek);


	INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
	(1./(1+ksat*exp((n-1)*Vm*F/(R*T))))*
	(exp(n*Vm*F/(R*T))*Nai*Nai*Nai*Cao-
	exp((n-1)*Vm*F/(R*T))*Nao*Nao*Nao*Cai*2.5);


	INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;


	IpCa=GpCa*Cai/(KpCa+Cai);


	IpK=GpK*rec_ipK*(Vm-Ek);


	IbNa=GbNa*(Vm-Ena);


	IbCa=GbCa*(Vm-Eca);


	if       (cell_type == 100)   GNaL = 0.0065;
	else if  (cell_type == 102)  GNaL = 0.0065;
	else if  (cell_type == 101) GNaL = 0.0095;
	
	INaL = GNaL*mNaL*mNaL*mNaL*hNaL*(Vm-Ena);

	Gsac = 1;//0.025/185;//110.0;
	Esac = 1.0;//-20.0;
	lambda_sac = 2;//1.10;
	
	strain = (SL - SLrest)/SLrest;
	strain_half_maximal = 0.233;
	slope_factor = 0.0203;//0.018;
	open_probability = 1.0/( 1 + std::exp( -( (strain-strain_half_maximal)/slope_factor ) ) );
	//Isac = Gsac*(lambda_sac-1)*(Vm-Esac);
	Isac = Gsac*open_probability*(Vm-Esac);
	//Isac = 0.0;
	
	//Calculate total ionic current
	Iion = -(IKr + IKs + IK1 + Ito +	INa + IbNa + ICaL +	IbCa + INaK + INaCa + IpCa + IpK + INaL + Istim + Isac);

	
	//Concentration derivatives	
	kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
	k1=k1_/kCaSR;
	k2=k2_*kCaSR;
	dRR=k4*(1-sRR)-k2*CaSS*sRR;
	Variable_Derivatives[sRR_tt]=dRR;
	const double sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);
	
	
	Irel=Vrel*sOO*(CaSR-CaSS);
	Ileak=Vleak*(CaSR-Cai);
	Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
	Ixfer=Vxfer*(CaSS-Cai);
	
	
	CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
	dCaSR=dt_rushlarsen*(Iup-Irel-Ileak);
	bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
	cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
	// CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	Variable_Derivatives[CaSR_tt] = (-CaSR + (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2)/dt_rushlarsen;
	// double dCaSR = (Iup - Irel - Ileak)/(1.0  + (Kbufsr*Bufsr)/pow(CaSR + Kbufsr, 2.0));
	// Variable_Derivatives[CaSR_tt] = (Iup - Irel - Ileak)/(1.0  + (Kbufsr*Bufsr)/pow(CaSR + Kbufsr, 2.0));
	Variable_Derivatives[CaSR_tt] = (Iup - Irel - Ileak)/(1.0  + (Kbufsr*Bufsr)/pow(CaSR + Kbufsr, 2.0));
	

	CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
	dCaSS=dt_rushlarsen*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
	bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
	ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
	// CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;
	Variable_Derivatives[CaSS_tt] = (-CaSS + (sqrt(bcss*bcss+4*ccss)-bcss)/2)/dt_rushlarsen;
	// double dCaSS =(-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE))/(1.0 + (Bufss*Kbufss)/pow(CaSS+Kbufss, 2.0));
	// Variable_Derivatives[CaSS_tt] = (-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE))/(1.0 + (Bufss*Kbufss)/pow(CaSS+Kbufss, 2.0));
	Variable_Derivatives[CaSS_tt] = (-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE))/(1.0 + (Bufss*Kbufss)/pow(CaSS+Kbufss, 2.0));



	CaBuf=Bufc*Cai/(Cai+Kbufc);

	// oomph_info << "IbCa " << IbCa << std::endl;
	// oomph_info << "IpCa " << IpCa << std::endl;
	// oomph_info << "INaCa " << INaCa << std::endl;
	// oomph_info << "Isac " << Isac << std::endl;
	// oomph_info << "inverseVcF2 " << inverseVcF2 << std::endl;
	// oomph_info << "CAPACITANCE " << CAPACITANCE << std::endl;
	// oomph_info << "Iup " << Iup << std::endl;
	// oomph_info << "Ileak " << Ileak << std::endl;
	// oomph_info << "Vsr " << Vsr << std::endl;
	// oomph_info << "Vc " << Vc << std::endl;
	// oomph_info << "Ixfer " << Ixfer << std::endl;
	// oomph_info << "TropToT_dt " << TropToT_dt << std::endl;


	dCai=dt_rushlarsen*(((-(IbCa+IpCa-2*INaCa+(Isac/3.0))*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer) - TropToT_dt/1000.0);


	// oomph_info << "Bufc " << Bufc << std::endl;
	// oomph_info << "CaBuf " << CaBuf << std::endl;
	// oomph_info << "dCai " << dCai << std::endl;
	// oomph_info << "Cai " << Cai << std::endl;
	// oomph_info << "Kbufc " << Kbufc << std::endl;

	const double bc=Bufc-CaBuf-dCai-Cai+Kbufc;
	cc=Kbufc*(CaBuf+dCai+Cai);
	// Cai=(sqrt(bc*bc+4*cc)-bc)/2;
	
	// oomph_info << Cai << std::endl;
	// oomph_info << bc << std::endl;
	// oomph_info << cc << std::endl;

	Variable_Derivatives[Cai_tt] = (-Cai + (sqrt(bc*bc+4*cc)-bc)/2)/dt_rushlarsen;
	// oomph_info << Variable_Derivatives[Cai_tt] << std::endl;
	// double dCai = ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer)/(1.0 + (Bufc*Kbufc)/pow(Cai + Kbufc, 2.0));
	// Variable_Derivatives[Cai_tt] = ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer)/(1.0 + (Bufc*Kbufc)/pow(Cai + Kbufc, 2.0));
    	
	Variable_Derivatives[Cai_tt] = ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer)/(1.0 + (Bufc*Kbufc)/pow(Cai + Kbufc, 2.0)) - TropToT_dt/1000.0;

	
	dNai=-(INa+IbNa+3*INaK+3*INaCa+(Isac/3.0))*inverseVcF*CAPACITANCE;
	Variable_Derivatives[Nai_tt]=dNai;
	
	dKi=-(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK+(Isac/3.0))*inverseVcF*CAPACITANCE;
	Variable_Derivatives[Ki_tt]=dKi;



			
	AM=1./(1.+exp((-60.-Vm)/5.));
	BM=0.1/(1.+exp((Vm+35.)/5.))+0.10/(1.+exp((Vm-50.)/200.));
	TAU_M=AM*BM;
	M_INF=1./((1.+exp((-56.86-Vm)/9.03))*(1.+exp((-56.86-Vm)/9.03)));
	if (Vm>=-40.)
	{
		AH_1=0.; 
		BH_1=(0.77/(0.13*(1.+exp(-(Vm+10.66)/11.1))));
		TAU_H= 1.0/(AH_1+BH_1);
	}
	else
	{
		AH_2=(0.057*exp(-(Vm+80.)/6.8));
		BH_2=(2.7*exp(0.079*Vm)+(3.1e5)*exp(0.3485*Vm));
		TAU_H=1.0/(AH_2+BH_2);
	}
	H_INF=1./((1.+exp((Vm+71.55)/7.43))*(1.+exp((Vm+71.55)/7.43)));
	if(Vm>=-40.)
	{
		AJ_1=0.;      
		BJ_1=(0.6*exp((0.057)*Vm)/(1.+exp(-0.1*(Vm+32.))));
		TAU_J= 1.0/(AJ_1+BJ_1);
	}
	else
	{
		AJ_2=(((-2.5428e4)*exp(0.2444*Vm)-(6.948e-6)*
			   exp(-0.04391*Vm))*(Vm+37.78)/
			  (1.+exp(0.311*(Vm+79.23))));    
		BJ_2=(0.02424*exp(-0.01052*Vm)/(1.+exp(-0.1378*(Vm+40.14))));
		TAU_J= 1.0/(AJ_2+BJ_2);
	}
	J_INF=H_INF;
	
	Xr1_INF=1./(1.+exp((-26.-Vm)/7.));
	axr1=450./(1.+exp((-45.-Vm)/10.));
	bxr1=6./(1.+exp((Vm-(-30.))/11.5));
	TAU_Xr1=axr1*bxr1;
	Xr2_INF=1./(1.+exp((Vm-(-88.))/24.));
	axr2=3./(1.+exp((-60.-Vm)/20.));
	bxr2=1.12/(1.+exp((Vm-60.)/20.));
	TAU_Xr2=axr2*bxr2;
	
	Xs_INF=1./(1.+exp((-5.-Vm)/14.));
	Axs=(1400./(sqrt(1.+exp((5.-Vm)/6))));
	Bxs=(1./(1.+exp((Vm-35.)/15.)));
	TAU_Xs=Axs*Bxs+80;
	
	if (cell_type == 100) {
		R_INF=1./(1.+exp((20-Vm)/6.));
		S_INF=1./(1.+exp((Vm+20)/5.));
		TAU_R=9.5*exp(-(Vm+40.)*(Vm+40.)/1800.)+0.8;
		TAU_S=85.*exp(-(Vm+45.)*(Vm+45.)/320.)+5./(1.+exp((Vm-20.)/5.))+3.;
	}
	else if (cell_type == 102) {
		R_INF=1./(1.+exp((20-Vm)/6.));
		S_INF=1./(1.+exp((Vm+28)/5.));
		TAU_R=9.5*exp(-(Vm+40.)*(Vm+40.)/1800.)+0.8;
		TAU_S=1000.*exp(-(Vm+67)*(Vm+67)/1000.)+8.;
	}
	else {
		R_INF=1./(1.+exp((20-Vm)/6.));
		S_INF=1./(1.+exp((Vm+20)/5.));
		TAU_R=9.5*exp(-(Vm+40.)*(Vm+40.)/1800.)+0.8;
		TAU_S=85.*exp(-(Vm+45.)*(Vm+45.)/320.)+5./(1.+exp((Vm-20.)/5.))+3.;
	}		
	
	D_INF=1./(1.+exp((-8-Vm)/7.5));
	Ad=1.4/(1.+exp((-35-Vm)/13))+0.25;
	Bd=1.4/(1.+exp((Vm+5)/5));
	Cd=1./(1.+exp((50-Vm)/20));
	TAU_D=Ad*Bd+Cd;
	F_INF=1./(1.+exp((Vm+20)/7));
	Af=1102.5*exp(-(Vm+27)*(Vm+27)/225);
	Bf=200./(1+exp((13-Vm)/10.));
	Cf=(180./(1+exp((Vm+30)/10)))+20;
	TAU_F=Af+Bf+Cf;
	F2_INF=0.67/(1.+exp((Vm+35)/7))+0.33;
	Af2=600*exp(-(Vm+25)*(Vm+25)/170);
	Bf2=31/(1.+exp((25-Vm)/10));
	Cf2=16/(1.+exp((Vm+30)/10));
	TAU_F2=Af2+Bf2+Cf2;
	FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
	TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;
	
	//INaL Area
	alpha_mNaL = (0.32*(Vm+47.13))/(1-std::exp(-0.1*(Vm+47.13)));
	beta_mNaL  = 0.08*std::exp(-Vm/11.0);
	mNaL_INF   = alpha_mNaL/(alpha_mNaL + beta_mNaL);
	TAU_mNaL   = 1.0/(alpha_mNaL + beta_mNaL); //TAU_hNaL initialised in constructor.
	hNaL_INF   = 1.0/(1+std::exp((Vm+91.0)/6.1));


	// Variable_Derivatives[sm_tt] 	= ( -sm			+	M_INF-(M_INF-sm)*exp(-dt_rushlarsen/TAU_M))/dt_rushlarsen;
	// Variable_Derivatives[sh_tt] 	= ( -sh			+	H_INF-(H_INF-sh)*exp(-dt_rushlarsen/TAU_H))/dt_rushlarsen;
	// Variable_Derivatives[sj_tt] 	= ( -sj			+	J_INF-(J_INF-sj)*exp(-dt_rushlarsen/TAU_J))/dt_rushlarsen;
	// Variable_Derivatives[sxr1_tt] 	=( -sxr1		+	Xr1_INF-(Xr1_INF-sxr1)*exp(-dt_rushlarsen/TAU_Xr1))/dt_rushlarsen;
	// Variable_Derivatives[sxr2_tt] 	=( -sxr2		+	Xr2_INF-(Xr2_INF-sxr2)*exp(-dt_rushlarsen/TAU_Xr2))/dt_rushlarsen;
	// Variable_Derivatives[sxs_tt] 	= ( -sxs		+	Xs_INF-(Xs_INF-sxs)*exp(-dt_rushlarsen/TAU_Xs))/dt_rushlarsen;
	// Variable_Derivatives[ss_tt] 	= ( -ss			+	S_INF-(S_INF-ss)*exp(-dt_rushlarsen/TAU_S))/dt_rushlarsen;
	// Variable_Derivatives[sr_tt] 	= ( -sr			+	R_INF-(R_INF-sr)*exp(-dt_rushlarsen/TAU_R))/dt_rushlarsen;
	// Variable_Derivatives[sd_tt] 	= ( -sd			+	D_INF-(D_INF-sd)*exp(-dt_rushlarsen/TAU_D))/dt_rushlarsen;
	// Variable_Derivatives[sf_tt] 	= ( -sf			+	F_INF-(F_INF-sf)*exp(-dt_rushlarsen/TAU_F))/dt_rushlarsen;
	// Variable_Derivatives[sf2_tt] 	= ( -sf2		+	F2_INF-(F2_INF-sf2)*exp(-dt_rushlarsen/TAU_F2))/dt_rushlarsen;
	// Variable_Derivatives[sfcass_tt] = ( -sfcass		+	FCaSS_INF-(FCaSS_INF-sfcass)*exp(-dt_rushlarsen/TAU_FCaSS))/dt_rushlarsen;
	// Variable_Derivatives[mNaL_tt] 	= (-mNaL		+	mNaL_INF-(mNaL_INF - mNaL)*exp(-dt_rushlarsen/TAU_mNaL))/dt_rushlarsen;
	// Variable_Derivatives[hNaL_tt] 	= (-hNaL		+	hNaL_INF-(hNaL_INF - hNaL)*exp(-dt_rushlarsen/TAU_hNaL))/dt_rushlarsen;

	Variable_Derivatives[sm_tt] 	= 	(M_INF-sm)/TAU_M;
	Variable_Derivatives[sh_tt] 	= 	(H_INF-sh)/TAU_H;
	Variable_Derivatives[sj_tt] 	= 	(J_INF-sj)/TAU_J;
	// Variable_Derivatives[sxr1_tt] 	=	(Xr1_INF-sxr1)/TAU_Xr1;
	// Variable_Derivatives[sxr2_tt] 	=	(Xr2_INF-sxr2)/TAU_Xr2;
	Variable_Derivatives[sxs_tt] 	= 	(Xs_INF-sxs)/TAU_Xs;
	Variable_Derivatives[ss_tt] 	= 	(S_INF-ss)/TAU_S;
	Variable_Derivatives[sr_tt] 	= 	(R_INF-sr)/TAU_R;
	Variable_Derivatives[sd_tt] 	= 	(D_INF-sd)/TAU_D;
	Variable_Derivatives[sf_tt] 	= 	(F_INF-sf)/TAU_F;
	Variable_Derivatives[sf2_tt] 	= 	(F2_INF-sf2)/TAU_F2;
	Variable_Derivatives[sfcass_tt] = 	(FCaSS_INF-sfcass)/TAU_FCaSS;
	Variable_Derivatives[mNaL_tt] 	=	(mNaL_INF - mNaL)/TAU_mNaL;
	Variable_Derivatives[hNaL_tt] 	=	(hNaL_INF - hNaL)/TAU_hNaL;


	//Not doing this here
	// UpdateVoltage();


	// oomph_info << "sm_tt " << Variable_Derivatives[sm_tt] << std::endl;
	// oomph_info << "sh_tt " << Variable_Derivatives[sh_tt] << std::endl;
	// oomph_info << "sj_tt " << Variable_Derivatives[sj_tt] << std::endl;
	// oomph_info << "sqt1_O_tt " << Variable_Derivatives[sqt1_O_tt] << std::endl;
	// oomph_info << "sqt1_C1_tt " << Variable_Derivatives[sqt1_C1_tt] << std::endl;
	// oomph_info << "sqt1_C2_tt " << Variable_Derivatives[sqt1_C2_tt] << std::endl;
	// oomph_info << "sqt1_C3_tt " << Variable_Derivatives[sqt1_C3_tt] << std::endl;
	// oomph_info << "sqt1_I_tt " << Variable_Derivatives[sqt1_I_tt] << std::endl;
	// oomph_info << "sxs_tt " << Variable_Derivatives[sxs_tt] << std::endl;
	// oomph_info << "ss_tt " << Variable_Derivatives[ss_tt] << std::endl;
	// oomph_info << "sr_tt " << Variable_Derivatives[sr_tt] << std::endl;
	// oomph_info << "sd_tt " << Variable_Derivatives[sd_tt] << std::endl;
	// oomph_info << "sf_tt " << Variable_Derivatives[sf_tt] << std::endl;
	// oomph_info << "sf2_tt " << Variable_Derivatives[sf2_tt] << std::endl;
	// oomph_info << "sfcass_tt " << Variable_Derivatives[sfcass_tt] << std::endl;
	// oomph_info << "sRR_tt " << Variable_Derivatives[sRR_tt] << std::endl;
	// oomph_info << "Cai_tt " << Variable_Derivatives[Cai_tt] << std::endl;
	// oomph_info << "CaSR_tt " << Variable_Derivatives[CaSR_tt] << std::endl;
	// oomph_info << "CaSS_tt " << Variable_Derivatives[CaSS_tt] << std::endl;
	// oomph_info << "mNaL_tt " << Variable_Derivatives[mNaL_tt] << std::endl;
	// oomph_info << "hNaL_tt " << Variable_Derivatives[hNaL_tt] << std::endl;
	// oomph_info << "Nai_tt " << Variable_Derivatives[Nai_tt] << std::endl;
	// oomph_info << "Ki_tt " << Variable_Derivatives[Ki_tt] << std::endl;
	// oomph_info << "N_NoXB_tt " << Variable_Derivatives[N_NoXB_tt] << std::endl;
	// oomph_info << "N_tt " << Variable_Derivatives[N_tt] << std::endl;
	// oomph_info << "P_NoXB_tt " << Variable_Derivatives[P_NoXB_tt] << std::endl;
	// oomph_info << "P_tt " << Variable_Derivatives[P_tt] << std::endl;
	// oomph_info << "XBprer_tt " << Variable_Derivatives[XBprer_tt] << std::endl;
	// oomph_info << "XBpostr_tt " << Variable_Derivatives[XBpostr_tt] << std::endl;
	// oomph_info << "SL_tt " << Variable_Derivatives[SL_tt] << std::endl;
	// oomph_info << "xXBpostr_tt " << Variable_Derivatives[xXBpostr_tt] << std::endl;
	// oomph_info << "xXBprer_tt " << Variable_Derivatives[xXBprer_tt] << std::endl;
	// oomph_info << "TropCaL_tt " << Variable_Derivatives[TropCaL_tt] << std::endl;
	// oomph_info << "TropCaH_tt " << Variable_Derivatives[TropCaH_tt] << std::endl;
	// oomph_info << "intf0_tt " << Variable_Derivatives[intf0_tt] << std::endl;
}



void IsmailTNNP06::get_optional_output(const Boost_State_Type &Variables,
								const double &t,
								const unsigned &cell_type,
								const double &Istim,
								const Vector<double> &Other_Parameters,
								const Vector<double> &Other_Variables,
								Vector<double> &Out)
{
	Out[0] = (Variables[SL_tt] - SLrest)/SLrest;
}

double IsmailTNNP06::get_active_strain()
{
	return (Variables[SL_tt] - SLrest)/SLrest;
}



}; //End namespace