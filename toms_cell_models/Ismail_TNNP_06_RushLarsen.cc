#include "Ismail_TNNP_06_RushLarsen.h"


namespace oomph{

IsmailTNNP06RushLarsen::IsmailTNNP06RushLarsen(const unsigned& number_of_backup_values) : CellModelBaseFullySegregated(number_of_backup_values)
{
	//Default cell type
	cell_type = 100;

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
	// state[sqt1_C3_tt] = 1.0;
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

	m_epiMidRatio = 1.6;


	//======
	// CVODE
	//======
	NEQ = 12+2;
	// time = 0.0;
	T1 = 2.0;			// first output time 
	TMULT = 10;			// output time factor
	NOUT = 4500;
	rtol = 1e-7;		// scalar relative tolerance
	atol = 1e-6;		// absolute tolerance low computational cost
	// TOUT = time + m_HT;
	IOUT = 0.0;
	// NOUT = static_cast<int>(0.1/m_HT);


	SUNContext_Create(NULL, &sunctx);


	abstol = N_VNew_Serial(12, sunctx);

	NV_Ith_S(abstol,0) = atol;
	NV_Ith_S(abstol,1) = atol;
	NV_Ith_S(abstol,2) = atol;
	NV_Ith_S(abstol,3) = atol;
	NV_Ith_S(abstol,4) = atol;
	NV_Ith_S(abstol,5) = atol;
	NV_Ith_S(abstol,6) = atol;
	NV_Ith_S(abstol,7) = atol;
	NV_Ith_S(abstol,8) = atol;
	NV_Ith_S(abstol,9) = atol;
	NV_Ith_S(abstol,10) = atol;
	NV_Ith_S(abstol,11) = atol;

	// // Call CVodeCreate to create the solver memory and specify the 
	// // Backward Differentiation Formula and the use of a Newton iteration
	// cvode_mem = CVodeCreate(CV_BDF, sunctx);
	// if (cvode_mem == 0) {
	// 	std::cerr << "Error in CVodeMalloc: could not allocate" << std::endl;
	// 	exit(0);
	// }
	
	// // Should the BDF stability limit detection algorithm be used
	// flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
	
	// // Call CVodeSetUserData to speciﬁes the user data block
	// flag = CVodeSetUserData(cvode_mem, static_cast<void*>(this));
	// if (flag == CV_MEM_NULL) {
	// 	std::cerr << "Error in attaching User Data" << std::endl;
	// 	exit(0);
	// }
	
	// /* Call CVodeInit to initialize the integrator memory and specify the
	//  * user's right hand side function in y'=f(t,y), the inital time T0, and
	//  * the initial dependent variable vector y. */

	y = N_VNew_Serial(12, sunctx);
	Jacobian = SUNDenseMatrix(12, 12, sunctx);

	LS = SUNLinSol_Dense(y, Jacobian, sunctx);

	// NV_Ith_S(y,0) =	0.99;
	// NV_Ith_S(y,1) =	0.97;
	// NV_Ith_S(y,2) =	0.01;
	// NV_Ith_S(y,3) =	0.01;
	// NV_Ith_S(y,4) =	0.01;
	// NV_Ith_S(y,5) =	0.01;
	// NV_Ith_S(y,6) =	2.17149;
	// NV_Ith_S(y,7) =	0.007;
	// NV_Ith_S(y,8) =	0;
	// NV_Ith_S(y,9) =	0.0144725;
	// NV_Ith_S(y,10) =	0.232095;
	// NV_Ith_S(y,11) =	0.0;

	// flag = CVodeInit(cvode_mem, rice_myofilament_ode, 0.0, y);
	// //		if (flag < 0) {
	// //			std::cerr << "Error in CVodeMalloc: could not allocate" << std::endl;
	// //			exit(0);
	// //		}
	
	// /* Call CVodeSVtolerances to specify the scalar relative tolerance
	//  * and vector absolute tolerances */
	// flag = CVodeSVtolerances(cvode_mem, rtol, abstol);
	// if (flag < 0) {
	// 	std::cerr << "Error in CVodeMalloc: could not allocate" << std::endl;
	// 	exit(0);
	// }

	
	// /* Call CVDense to specify the CVDENSE dense linear solver */
	// // flag = CVDense(cvode_mem, 12);
	// // flag = CVBand(cvode_mem, 12, 2,2);

	// LS = SUNLinSol_Dense(y, Jacobian, sunctx);

	// CVodeSetLinearSolver(cvode_mem, LS, Jacobian);



	//Assign the names of the variables used by this model
	Names_Of_Cell_Variables =
	{
		"vm",
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
		//For rice myofilament
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

	Names_Of_Output_Data =
	{
		"Strain"
	};

	FinalizeConstruction();
}

IsmailTNNP06RushLarsen::~IsmailTNNP06RushLarsen()
{
	SUNLinSolFree(LS);
	SUNContext_Free(&sunctx);
}

double IsmailTNNP06RushLarsen::get_initial_state_variable(const unsigned &v)
{
	if ( m_mutant ) { // SQT1
		if ( cell_type == EPI ) {
			switch(v){
				case vm_tt: return -78.0773;
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
		else if ( cell_type == MCELL ) {
			switch(v){
				case vm_tt: return -77.9261;
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
		else if ( cell_type == ENDO ) {
			switch(v){
				case vm_tt: return -78.119;
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
		else
		{
			std::ostringstream error_message;
			error_message << cell_type << " is not a valid cell type"<< std::endl;

			throw OomphLibError(
				        error_message.str(),
				        OOMPH_CURRENT_FUNCTION,
				        OOMPH_EXCEPTION_LOCATION);
			return 0.0; //to shut up compiler
		}
	}
	else { // WT
		if ( cell_type == EPI ) {
			switch(v){
				case vm_tt: return -77.7615;
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
		else if ( cell_type == MCELL ) {
			switch(v){
				case vm_tt: return -77.4981;
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
		else if ( cell_type == ENDO ) {
			switch(v){
				case vm_tt: return -77.7461;
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
		else
		{
			std::ostringstream error_message;
	           error_message
	            << cell_type << " is not a valid cell type"<< std::endl;

			throw OomphLibError(
				        error_message.str(),
				        OOMPH_CURRENT_FUNCTION,
				        OOMPH_EXCEPTION_LOCATION);
			return 0.0; //to shut up compiler
		}
	}
	std::ostringstream error_message;
	           error_message
	            << v << " is not a variable number"<< std::endl;

	throw OomphLibError(
		        error_message.str(),
		        OOMPH_CURRENT_FUNCTION,
		        OOMPH_EXCEPTION_LOCATION);
	return 0.0; //to shut up compiler
}

void IsmailTNNP06RushLarsen::TakeTimestep(const double& dt, const double& t, double* state)
{

	//We don't want to take a time-step larger than 0.01ms
	if(dt>0.01)
	{
		const unsigned N_solve_at_0p01 = unsigned(dt/0.01);
		const double remainder = dt - N_solve_at_0p01*0.01;
		double t_running = t;

		for(unsigned n=0; n<N_solve_at_0p01; n++)
		{
			TakeTimestep(0.01, t_running, state);
			t_running+=0.01;
		}
		if(remainder>0.0)
		{
			TakeTimestep(remainder, t_running, state);
		}

		return;
	}

	// oomph_info << state[vm_tt] << std::endl;
	// oomph_info << state[sm_tt] << std::endl;
	// oomph_info << state[sh_tt] << std::endl;
	// oomph_info << state[sj_tt] << std::endl;
	// oomph_info << state[sqt1_O_tt] << std::endl;
	// oomph_info << state[sqt1_C1_tt] << std::endl;
	// oomph_info << state[sqt1_C2_tt] << std::endl;
	// oomph_info << state[sqt1_C3_tt] << std::endl;
	// oomph_info << state[sqt1_I_tt] << std::endl;
	// oomph_info << state[sxs_tt] << std::endl;
	// oomph_info << state[ss_tt] << std::endl;
	// oomph_info << state[sr_tt] << std::endl;
	// oomph_info << state[sd_tt] << std::endl;
	// oomph_info << state[sf_tt] << std::endl;
	// oomph_info << state[sf2_tt] << std::endl;
	// oomph_info << state[sfcass_tt] << std::endl;
	// oomph_info << state[sRR_tt] << std::endl;
	// oomph_info << state[Cai_tt] << std::endl;
	// oomph_info << state[CaSR_tt] << std::endl;
	// oomph_info << state[CaSS_tt] << std::endl;
	// oomph_info << state[mNaL_tt] << std::endl;
	// oomph_info << state[hNaL_tt] << std::endl;
	// oomph_info << state[Nai_tt] << std::endl;
	// oomph_info << state[Ki_tt] << std::endl;
	// oomph_info << state[N_NoXB_tt] << std::endl;
	// oomph_info << state[N_tt] << std::endl;
	// oomph_info << state[P_NoXB_tt] << std::endl;
	// oomph_info << state[P_tt] << std::endl;
	// oomph_info << state[XBprer_tt] << std::endl;
	// oomph_info << state[XBpostr_tt] << std::endl;
	// oomph_info << state[SL_tt] << std::endl;
	// oomph_info << state[xXBpostr_tt] << std::endl;
	// oomph_info << state[xXBprer_tt] << std::endl;
	// oomph_info << state[TropCaL_tt] << std::endl;
	// oomph_info << state[TropCaH_tt] << std::endl;
	// oomph_info << state[intf0_tt] << std::endl;

	// const double Vm = state[vm_tt];
	// const double state[sm_tt]			= state[sm_tt];
	// const double state[sh_tt]			= state[sh_tt];
	// const double state[sj_tt]			= state[sj_tt];

	// //For H-H type
	// // const double sxr1		= state[sxr1_tt];
	// // const double sxr2		= state[sxr2_tt];

	// //For MC type
	// const double state[sqt1_O_tt]		= state[sqt1_O_tt];
	// const double state[sqt1_C1_tt]	= state[sqt1_C1_tt];
	// const double state[sqt1_C2_tt]	= state[sqt1_C2_tt];
	// const double state[sqt1_C3_tt]	= state[sqt1_C3_tt];
	// const double state[sqt1_I_tt]		= state[sqt1_I_tt];

	// const double state[sxs_tt]		= state[sxs_tt];
	// const double state[ss_tt]			= state[ss_tt];
	// const double state[sr_tt]			= state[sr_tt];
	// const double state[sd_tt]			= state[sd_tt];
	// const double state[sf_tt]			= state[sf_tt];
	// const double state[sf2_tt]		= state[sf2_tt];
	// const double state[sfcass_tt]		= state[sfcass_tt];
	// const double state[sRR_tt]		= state[sRR_tt];
	// const double sOO		= state[sOO_tt];
	// const double state[Cai_tt]		= state[Cai_tt];
	// const double state[CaSR_tt]		= state[CaSR_tt];
	// const double state[CaSS_tt]		= state[CaSS_tt];
	// const double state[mNaL_tt]		= state[mNaL_tt];
	// const double state[hNaL_tt]		= state[hNaL_tt];
	// const double state[Nai_tt]		= state[Nai_tt];
	// const double state[Ki_tt]			= state[Ki_tt];


	// const double state[N_NoXB_tt]		= state[N_NoXB_tt];
	// const double state[N_tt]			= state[N_tt];
	// const double state[P_NoXB_tt]		= state[P_NoXB_tt];
	// const double P			= state[P_tt];
	// const double XBprer		= state[XBprer_tt];
	// const double XBpostr	= state[XBpostr_tt];
	// double SL			= state[SL_tt];
	// const double xXBpostr	= state[xXBpostr_tt];
	// const double xXBprer	= state[xXBprer_tt];
	// const double TropCaL	= state[TropCaL_tt];
	// const double TropCaH	= state[TropCaH_tt];
	// const double intf0 		= state[intf0_tt];

	


	if (cell_type == EPI)
	{
		Gks=0.392;
	}
	else if (cell_type == ENDO)
	{
		Gks=0.392;
	}
	else
	{
		Gks=0.098;
	}
	
	GK1=5.405;
	//Parameters for Ito
	if (cell_type == EPI)
	{
		Gto=0.294;
	}
	else if (cell_type == ENDO)
	{
		Gto=0.073;
	}
	else
	{
		Gto=0.294;
	}

	// //ORiginally cvode solved myocardium stuff
	// flag = CVode(cvode_mem, TOUT, y, &time, CV_NORMAL);	

	
	// //state[N_NoXB_tt] = state[Cai_tt] * 1000;
	
	// /* Handle Ca binding to troponin here */
	// //  Adjust for temperature
	// konT	= kon	*std::pow(Qkon,(T-273.0-37.0)/10.);
	// koffLT	= koffL	*std::pow(Qkoff,(T-273.0-37.0)/10.);
	// koffHT	= koffH	*std::pow(Qkoff,(T-273.0-37.0)/10.);
	
	// /* Compute derived regulatory unit rates, Ca in in nM to prevent Ca overflow (max ~ 10 in auto), TropCaL and TropCaH are normailized so no Units. */
 //  	Variable_Derivatives[TropCaL_tt] = konT * state[Cai_tt]*1000 * (1.0 - TropCaL) - koffLT * TropCaL;                                              
 //  	Variable_Derivatives[TropCaH_tt] = konT * state[Cai_tt]*1000 * (1.0 - TropCaH) - koffHT * TropCaH; 
	// /* Compute rates for state[N_tt] to P transitions */
	// // Compute z-line end of single overlap region
	// sovr_ze = std::min(len_thick/2.,SL/2.);
	
	// // Compute centerline of end of single overlap region
 //    sovr_cle = std::max(SL/2.-(SL-len_thin),len_hbare/2.);
	
	// // Compute length of the single overlap 
	// len_sovr = sovr_ze-sovr_cle;
	
	// // Compute overlap fraction for thick filament
	// SOVFThick = len_sovr*2./(len_thick-len_hbare);
	
	// // Compute overlap fraction for thin filament
	// SOVFThin = len_sovr/len_thin;
	
	// /* Compute combined Ca binding to high (w/XB) and low (w/o XB)  No Units. */
	// perm = (1.0 - SOVFThin) * TropCaL + SOVFThin * TropCaH;
	
	// /* Set perm variable that control n to p transiton in Hill-like fashion, No Units. */
	// permtot = std::pow((1.0 / (1.0 + std::pow((perm50 / perm), nperm))), 0.5);
 // 	if (100.0 < 1.0/permtot)
	// 	inprmt=100.0;
 //  	else 
	// 	inprmt = 1.0/permtot;
	
	// /* Adjust for Ca activation level and temp. */
	// kn_pT = kn_p*permtot*std::pow(Qkn_p,(T-273.0-37.0)/10.);
	// kp_nT = kp_n*inprmt*std::pow(Qkp_n,(T-273.0-37.0)/10.);
	
	// /* Compute fapp, the attachment rate from weak to strong, pre-rotated state */
 //    fappT = fapp*std::pow(Qfapp,(T-273.0-37.0)/10.);
	
	// /* Compute gapp, the detachment rate from strong, pre-rotated to the weakly-bound state */
	// /* Compute SL modifier for gapp to increase rate st shorter SL */
 //  	gapslmd = 1.0 + (1.0 - SOVFThick)*gslmod;
 //    gappT = gapp*gapslmd*std::pow(Qgapp,(T-273.0-37.0)/10.);
	
	// /* Set rate for rates between pre-force and force states*/
	// /* Compute modifiers based on mean strain of states. */
	// oomph_info << "xXBprer " << xXBprer << std::endl;
	// oomph_info << "hfmdc " << hfmdc << std::endl;
	// oomph_info << "x_0 " << x_0 << std::endl;
 //  	hfmd=std::exp(-sign(xXBprer)*hfmdc*((xXBprer/x_0)*(xXBprer/x_0)));
	// // oomph_info << "xXBprer " << xXBprer << std::endl;
	// // oomph_info << "hfmdc " << hfmdc << std::endl;
	// // oomph_info << "x_0 " << x_0 << std::endl;
 //  	hbmd=std::exp(sign((xXBpostr-x_0))*hbmdc*(((xXBpostr-x_0)/x_0)*((xXBpostr-x_0)/x_0)));
	
	// /* Combine modifiers of hf and hb */
	// oomph_info << "hf " << hf << std::endl;
	// oomph_info << "hfmd " << hfmd << std::endl;
	// hfT=hf*hfmd*std::pow(Qhf,(T-273.0-37.0)/10.);
	// // oomph_info << "hf " << hf << std::endl;
	// // oomph_info << "hfmd " << hfmd << std::endl;
	// // oomph_info << "Qhf " << Qhf << std::endl;
	// // oomph_info << "T " << T << std::endl;
	// hbT=hb*hbmd*std::pow(Qhb,(T-273.0-37.0)/10.);
	
	// /* Set rate for rates gxb, the ATP using XB transition */
	// /* Add term for distortion dependence of gxb */
	// if(x_0>xXBpostr)
	// {
	// 	gxbmd = std::exp(sigmap*((x_0-xXBpostr)/x_0)*((x_0-xXBpostr)/x_0))+(1.0-std::exp(sigman*(((xXBpostr-x_0)/x_0)*(xXBpostr-x_0)/x_0)));
	// }
	// else
	// {
	// 	gxbmd = 1.0;
	// }
 //    // gxbmd = heav(x_0-xXBpostr)*std::exp(sigmap*((x_0-xXBpostr)/x_0)*((x_0-xXBpostr)/x_0))+(1.-heav(x_0-xXBpostr)*std::exp(sigman*(((xXBpostr-x_0)/x_0)*(xXBpostr-x_0)/x_0)));
	// // oomph_info << "x_0 " << x_0 << std::endl;
	// // oomph_info << "xXBpostr " << xXBpostr << std::endl;
	// // oomph_info << "sigmap " << sigmap << std::endl;
	// // oomph_info << "sigman " << sigman << std::endl;

 //    gxbT = gxb*gxbmd*std::pow(Qgxb,(T-273.0-37.0)/10.);

	// // oomph_info << "gxb " << gxb << std::endl;
	// // oomph_info << "gxbmd " << gxbmd << std::endl;
	// // oomph_info << "Qgxb " << Qgxb << std::endl;
	// // oomph_info << "T " << T << std::endl;

	
	// /* update all RUs */
	// Variable_Derivatives[N_NoXB_tt] = -kn_pT * state[N_NoXB_tt] + kp_nT * state[P_NoXB_tt];
	// Variable_Derivatives[P_NoXB_tt] = -kp_nT * state[P_NoXB_tt] + kn_pT * state[N_NoXB_tt];
	
	// Variable_Derivatives[N_tt] = -kn_pT * state[N_tt] + kp_nT * P; 
	// // oomph_info << "kp_nT " << kp_nT << std::endl;
	// // oomph_info << "P " << P << std::endl;
	// // oomph_info << "kn_pT " << kn_pT << std::endl;
	// // oomph_info << "state[N_tt] " << state[N_tt] << std::endl;
	// // oomph_info << "fappT " << fappT << std::endl;
	// // oomph_info << "gappT " << gappT << std::endl;
	// // oomph_info << "XBprer " << XBprer << std::endl;
	// // oomph_info << "gxbT " << gxbT << std::endl;
	// // oomph_info << "XBpostr " << XBpostr << std::endl;

	// oomph_info << "hfT " << hfT << std::endl;
	// oomph_info << "XBprer " << XBprer << std::endl;
	// oomph_info << "hbT " << hbT << std::endl;
	// oomph_info << "XBpostr " << XBpostr << std::endl;
	// oomph_info << "gxbT " << gxbT << std::endl;
	// oomph_info << "XBpostr " << XBpostr << std::endl;
	// Variable_Derivatives[P_tt] = -kp_nT * P + kn_pT * state[N_tt] - fappT * P + gappT * XBprer + gxbT * XBpostr; 
	// Variable_Derivatives[XBprer_tt] = fappT * P - gappT * XBprer - hfT * XBprer + hbT * XBpostr;
	// Variable_Derivatives[XBpostr_tt] = hfT * XBprer - hbT * XBpostr - gxbT * XBpostr;
	
	// /* compute steady-state fractions in XBprer and XBpostr using King-Altman rule */
	// SSXBprer = (hb*fapp+gxb*fapp)/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
	// SSXBpostr = fapp*hf/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
	// /* compute normalization for scaling active and passive force */
	// Fnordv = kxb*x_0*SSXBpostr;
	
	// /* compute force here */
	// force = kxb*SOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer);
	
	// /* Compute SL, Set to 0 for isometric, to 1 for active contraction */
	// /* Note that passive force is specified in terms of maximal active force */
	// /* read passive force from table file specified in terms of maximal active force */
	// /* compyet derivative of SL */
	
	// // Note that passive force is specified in terms of maximal active force.
	// // Passive force is computed in two ways depending on if cell is skinned or intact.
	// ppforce=sign((SL-SLrest))*PCon_t*(std::exp(PExp_t*std::abs((SL-SLrest)))-1);
 //    // if (singlecell == false) 
	// ppforce += std::max(0.0,(PCon_col*std::exp(PExp_col*(SL-SLcol))));
	
	// // Compute afterload either as a contant or from the series elastic element
 // 	afterload = 0.0;
	
 //  //   if (SEon == true) {
	// 	// afterload = KSE*(SLset-SL);
	// 	// if ((SEon_LengthClamp==true)&&(index==(pulse_number-1))) {
	// 	// 	afterload=5000.0*(SLset-SL);
	// 	// }
 // //    } 
	// // else {
	// // 	afterload=0.0; 
 // //    }
	
	// //  Compute the integral of forces 
	// Variable_Derivatives[intf0_tt]=(-ppforce+PreloadF+(-force/Fnordv)+afterload);
	
	// // Compute the derivative of SL for an isolated muscle
 //  	dSLisolated = ((intf0 + (SLset-SL)*visc)/massf )*heav(SL-SLmin)*heav(SLmax-SL);
 // 	// if (contraction==false) 
	// 	// dSLisolated = 0.;
	

	// //We're not using the whole heart model
	

	// // //========>>>>>>>>
	// // /*  The follow section of code implements is used to simulate ejection 
	// //  from canine whole heart into a Windkessel model of the circulatory system.
	// //  The complete description is available in:
	// //  De Tombe, PP and Little, WC. Inotropic effects of ejection 
	// //  are myocardial properties. Am J Physiol. 1994 Mar;266(3 Pt 2):H1202-13. */
	
	// // // Betavl is unitless
	// // Betavl=M_PI*std::pow(0.75*M_PI,0.6667);
	
 // //    // Compute A_0 which is the cross-sectional area of the ventricular wall
 // //    // A_0 has unit of area (volume^0.6667) in cm^2
	// // A_0 = Betavl*( std::pow(Vref + Vwall,0.6667) - std::pow(Vref,0.6667) );
	
 // //    // Compute SigmaP from normalized force
 // // 	SigPrs=SigFmx*(ppforce + force/Fnordv);
	
 // //    // Compute left venticular pressure in mmHg
 // // 	LVpress=SigPrs*A_0/(1.36*Betavl*std::pow(Ith(y,14),0.667));
	
	// // // Compute Flows - assumes aortic valve so flow in unidirectional
	// // // Units are l/s
	// // iflow=heav(LVpress-Ith(y,13))*((LVpress-Ith(y,13))/Rchar);
	// // // Allow active ejection on last beat only 
	// // if (index<(pulse_number-1))
	// // 	iflow=0.0;
	// // if (Isovolume==true) 
	// // 	iflow=0.0;
	
 // //    // Compute the derivative of arterial pressure
	// // Variable_Derivatives[ARTpset_tt]=(iflow/Cwind - (Ith(y,13)-ARTpset)/(Cwind*Ra));
 // //    if (WholeHeart==false) 
	// // 	Variable_Derivatives[ARTpset_tt]=0.0;
    
	// // // Compute the  derivative of the change in Left Ventricular Volume
	// // // Do not allow volume less than 0
 // //    Variable_Derivatives[LVV_tt]=-(iflow)*heav(Ith(y,14));
 // //    if (WholeHeart==false) 
	// // 	Variable_Derivatives[LVV_tt]=0.0;
	
	// // // Compute the  derivative of the change in SL
 // //    dSLwholeheart=SLref*(1/std::pow(Vref+0.333*Vwall,0.333))*(Variable_Derivatives[LVV_tt]/3)/std::pow(Ith(y,14)+0.333*Vwall,0.667);
	

	// // Compute the derivative of SL based using cell/tissue or simulated whole heart model
	// // if (WholeHeart==true) {
	// // 	Variable_Derivatives[SL_tt] = dSLwholeheart;
	// // } else {
	// 	Variable_Derivatives[SL_tt] = dSLisolated;
	// // }
	
	// //========>>>>>>>>
	
	
	// if(myo == false)
	// {
	// 	Variable_Derivatives[N_NoXB_tt]=0.0;
	// 	Variable_Derivatives[N_tt]=0.0;
	// 	Variable_Derivatives[P_NoXB_tt]=0.0;
	// 	Variable_Derivatives[P_tt]=0.0;
	// 	Variable_Derivatives[XBprer_tt]=0.0;
	// 	Variable_Derivatives[XBpostr_tt]=0.0;  
	// 	Variable_Derivatives[SL_tt]=0.0; 
	// 	Variable_Derivatives[TropCaL_tt]=0.0;
	// 	Variable_Derivatives[TropCaH_tt]=0.0;
	// 	Variable_Derivatives[intf0_tt]=0.0; 
	// 	// Variable_Derivatives[ARTpset_tt]=0.0;
	// 	// Variable_Derivatives[LVV_tt]=0.0;
	// }
	
	// dSL=Variable_Derivatives[SL_tt];
	
	// // Compute the duty fractions using King-Alman Rule
	// // Compute for states XBpref and XBf
	// dtyf_prer=(hbT*fappT+gxbT*fappT)/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
	// dtyf_postr=fappT*hfT/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
	
	// // oomph_info << "fappT " << fappT << std::endl;
	// // oomph_info << "hfT " << hfT << std::endl;
	// // oomph_info << "gxbT " << gxbT << std::endl;
	// // oomph_info << "gappT " << gappT << std::endl;
	// // oomph_info << "hbT " << hbT << std::endl;

	// // Compute mean strain by consider two competing effects - SL motion and XB cycling 
	// if(myo == true)
	// {
	// 	Variable_Derivatives[xXBprer_tt] = dSL/2. + xPsi*(1/dtyf_prer)*(-xXBprer*fappT+(xXBpostr-x_0-xXBprer)*hbT);
	// 	Variable_Derivatives[xXBpostr_tt] = dSL/2. + xPsi*(1/dtyf_postr)*(x_0+xXBprer-xXBpostr)*hfT;
	// 	if(std::abs(hfT)<1e-6)
	// 	{
	// 		Variable_Derivatives[xXBpostr_tt] = dSL/2. + xPsi*(x_0+xXBprer-xXBpostr)*(gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT)/fappT;
	// 	}

	// 	//Calculate by rush larsen
	// 	const double A = ((dSL/2.0)*dtyf_prer + xPsi*hbT*(xXBpostr-x_0))/(xPsi*(hbT+fappT));
	// 	const double B = -(xPsi/dtyf_prer)*(hbT + fappT);
	// 	Variable_Derivatives[xXBprer_tt] = (-xXBprer + A + (xXBprer - A)*exp(B*dt_rushlarsen) )/dt_rushlarsen;
	// }
	// else
	// {
	// 	Variable_Derivatives[xXBprer_tt]=0.0;
	// 	Variable_Derivatives[xXBpostr_tt]=0.0;
	// }
	
	// // oomph_info << "dSL " << dSL << std::endl;
	// // oomph_info << "xPsi " << xPsi << std::endl;
	// // oomph_info << "dtyf_postr " << dtyf_postr << std::endl;
	// // oomph_info << "x_0 " << x_0 << std::endl;
	// // oomph_info << "xXBprer " << xXBprer << std::endl;
	// // oomph_info << "xXBpostr " << xXBpostr << std::endl;
	// // oomph_info << "hfT " << hfT << std::endl;

	// // Compute derivative of z-line end of single overlap region
	// sovr_ze_dt=0.0;
	// if(len_thick/2.>SL/2.) 
	// {
	// 	sovr_ze_dt=-dSL/2.;
	// }
	
	// // Compute derivative of centerline of end of single overlap region
 //    sovr_cle_dt=0.0;

 //    if (SL/2.0-(SL-len_thin)>len_hbare/2.0)
 //    {
	// 	sovr_cle_dt=-dSL/2.0;
	// }
	
	// // Compute the derivative of the length of the single overlap 
	// len_sovr_dt=sovr_ze_dt-sovr_cle_dt;
 	
	// // Compute the derivative of the overlap fraction for thin filament

	// SOVFThin_dt = len_sovr_dt/len_thin;
	
	// // Compute the derivative of the overlap fraction for thick filament
	// SOVFThick_dt = len_sovr_dt*2./(len_thick-len_hbare);
	
	// /*Compute fraction of strongly bound for compute effective binding of Ca */
	// //double FrSBXB=((XBpostr+XBprer)/(SSXBpostr + SSXBprer))*SOVFThick;
	// //double FrSBXB_dt=((Variable_Derivatives[XBpostr_tt]+Variable_Derivatives[XBprer_tt])/(SSXBpostr + SSXBprer))*SOVFThick + ((XBpostr+XBprer)/(SSXBpostr + SSXBprer))*SOVFThick_dt;
	// FrSBXB=((XBpostr+XBprer)/(SSXBpostr + SSXBprer));

	// oomph_info << "Variable_Derivatives[XBpostr_tt] " << Variable_Derivatives[XBpostr_tt] << std::endl;
	// oomph_info << "Variable_Derivatives[XBprer_tt] " << Variable_Derivatives[XBprer_tt] << std::endl;
	// oomph_info << "SSXBpostr " << SSXBpostr << std::endl;
	// oomph_info << "SSXBprer " << SSXBprer << std::endl;

	// FrSBXB_dt=((Variable_Derivatives[XBpostr_tt]+Variable_Derivatives[XBprer_tt])/(SSXBpostr + SSXBprer));
	
	// TropToT=Trop_conc*((1.0-SOVFThin)*TropCaL + SOVFThin*(FrSBXB*TropCaH+(1.0-FrSBXB)*TropCaL));
	
	// oomph_info << "TropCaL "<< TropCaL << std::endl;
	// oomph_info << "Variable_Derivatives[TropCaL_tt] "<< Variable_Derivatives[TropCaL_tt] << std::endl;
	// oomph_info << "SOVFThin_dt "<< SOVFThin_dt << std::endl;
	// oomph_info << "SOVFThin "<< SOVFThin << std::endl;
	// oomph_info << "FrSBXB_dt "<< FrSBXB_dt << std::endl;
	// oomph_info << "TropCaH "<< TropCaH << std::endl;
	// oomph_info << "FrSBXB "<< FrSBXB << std::endl;
	// oomph_info << "Variable_Derivatives[TropCaH_tt] "<< Variable_Derivatives[TropCaH_tt] << std::endl;

	// TropToT_dt=Trop_conc*(-1.0*SOVFThin_dt*TropCaL+(1.0-SOVFThin)*Variable_Derivatives[TropCaL_tt] + SOVFThin_dt*(FrSBXB*TropCaH+(1.0-FrSBXB)*TropCaL) + SOVFThin*(FrSBXB_dt*TropCaH+FrSBXB*Variable_Derivatives[TropCaH_tt]-FrSBXB_dt*TropCaL+(1.0-FrSBXB)*Variable_Derivatives[TropCaL_tt]));
	// /*
	//  if(myo == true)
	//  Ith(ydot,26)=TropToT_dt; 
	//  else 
	//  Ith(ydot,89)=0.0; 
	//  */
	// //	/* Save non-state variable data */
	// //	temp1 = force/Fnordv;         /* normalized force */ 
	// //	temp1 = 120.0*force/Fnordv;   /* in units of mN/(mm^2) */
	// //active_force = force/Fnordv;
	// //active_force = force;
	// active_force = 10*force/Fnordv;   /* in units of mN/(mm^2) */
	
	// //	temp2 = SOVFThin;
	// //	temp3 = SOVFThick;
	// //	temp4 = SOVFThin_dt;
	// //	temp5 = SOVFThick_dt;
	// //	temp6 = perm;
	// //	temp7 = permtot;
	// //	temp8 = ppforce;
	// //	temp9 = PreloadF;
	// //	temp10 = FrSBXB; 
	// //	temp11 = inprmt;
	// //	temp12 = kn_pT;
	// //	temp13 = kp_nT;
	// //	temp14 = fappT;
	// //	temp15 = gappT;
	// //	temp16 = hfT;
	// //	temp17 = hbT;
	// //	temp18 = gxbT;
	// //	temp19 = sovr_ze;
	// //	temp20 = sovr_cle;
	// //	temp21 = sovr_ze_dt;
	// //	temp22 = sovr_cle_dt;
	// //	temp23 = len_sovr_dt;
	// //	temp24 = TropToT;
	// //	temp25 = TropToT_dt;
	// //	temp26 = dtyf_prer;
	// //	temp27 = dtyf_postr;
	// /****************************** myofilament ended *****************/
	// //std::cout << "In feval..." << std::endl;

if(myo)
	{
		//Assign the variables required by the filament model
		Cai_rice = state[Cai_tt];

		// Call CVodeCreate to create the solver memory and specify the 
		// Backward Differentiation Formula and the use of a Newton iteration
		cvode_mem = CVodeCreate(CV_BDF, sunctx);
		if (cvode_mem == 0) {
			oomph_info << "Error in CVodeMalloc: could not allocate" << std::endl;
			exit(0);
		}

		realtype t_cvode = t;

		NV_Ith_S(y,0) = state[N_NoXB_tt];
		NV_Ith_S(y,1) = state[N_tt];
		NV_Ith_S(y,2) = state[P_NoXB_tt];
		NV_Ith_S(y,3) = state[P_tt];
		NV_Ith_S(y,4) = state[XBprer_tt];
		NV_Ith_S(y,5) = state[XBpostr_tt];
		NV_Ith_S(y,6) = state[SL_tt];
		NV_Ith_S(y,7) = state[xXBpostr_tt];
		NV_Ith_S(y,8) = state[xXBprer_tt];
		NV_Ith_S(y,9) = state[TropCaL_tt];
		NV_Ith_S(y,10) = state[TropCaH_tt];
		NV_Ith_S(y,11) = state[intf0_tt];


		// oomph_info << "NV_Ith_S(y,0) " << NV_Ith_S(y,0) << std::endl;
		// oomph_info << "NV_Ith_S(y,1) " << NV_Ith_S(y,1) << std::endl;
		// oomph_info << "NV_Ith_S(y,2) " << NV_Ith_S(y,2) << std::endl;
		// oomph_info << "NV_Ith_S(y,3) " << NV_Ith_S(y,3) << std::endl;
		// oomph_info << "NV_Ith_S(y,4) " << NV_Ith_S(y,4) << std::endl;
		// oomph_info << "NV_Ith_S(y,5) " << NV_Ith_S(y,5) << std::endl;
		// oomph_info << "NV_Ith_S(y,6) " << NV_Ith_S(y,6) << std::endl;
		// oomph_info << "NV_Ith_S(y,7) " << NV_Ith_S(y,7) << std::endl;
		// oomph_info << "NV_Ith_S(y,8) " << NV_Ith_S(y,8) << std::endl;
		// oomph_info << "NV_Ith_S(y,9) " << NV_Ith_S(y,9) << std::endl;
		// oomph_info << "NV_Ith_S(y,10) " << NV_Ith_S(y,10) << std::endl;
		// oomph_info << "NV_Ith_S(y,11) " << NV_Ith_S(y,11) << std::endl;
		
		flag = CVodeInit(cvode_mem, rice_myofilament_ode, t_cvode, y);
		
		// flag = CVodeSVtolerances(cvode_mem, rtol, abstol);
		flag = CVodeSStolerances(cvode_mem, rtol, atol);
		if (flag < 0) {
			oomph_info << "Error in CVodeSStolerances: " << flag << std::endl;
			exit(0);
		}

		CVodeSetLinearSolver(cvode_mem, LS, Jacobian);

		// Should the BDF stability limit detection algorithm be used
		flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
		
		// Call CVodeSetUserData to speciﬁes the user data block
		flag = CVodeSetUserData(cvode_mem, static_cast<void*>(this));
		if (flag != CV_SUCCESS) {
			oomph_info << "Error in attaching User Data: " << flag << std::endl;
			exit(0);
		}
		
		/* Call CVodeInit to initialize the integrator memory and specify the
		 * user's right hand side function in y'=f(t,y), the inital time T0, and
		 * the initial dependent variable vector y. */

		// y = N_VNew_Serial(12, sunctx);
		// Jacobian = SUNDenseMatrix(12, 12, sunctx);

		



		

		//		if (flag < 0) {
		//			std::cerr << "Error in CVodeMalloc: could not allocate" << std::endl;
		//			exit(0);
		//		}
		/* Call CVodeSVtolerances to specify the scalar relative tolerance
		 * and vector absolute tolerances */
		

		/* Call CVDense to specify the CVDENSE dense linear solver */
		// flag = CVDense(cvode_mem, 12);
		// flag = CVBand(cvode_mem, 12, 2,2);
		// LS = SUNLinSol_Dense(y, Jacobian, sunctx);
		// CVodeSetLinearSolver(cvode_mem, LS, Jacobian);

		// N_Vector x;
		// y = N_VNew_Serial(12, sunctx);
		// IsmailTNNP06RushLarsen_ODE_Wrapper ode_wrapper(this);
		// integrate_adaptive( boost::numeric::odeint::make_controlled<controlled_error_stepper_type>(1e-7, 1e-7), ode_wrapper, x, t, t+dt_rushlarsen, dt_rushlarsen/100.0);

		// realtype T0 = t;
		// if (flag < 0) {
		// 	std::cerr << "Error in CVodeMalloc: could not allocate" << std::endl;
		// 	exit(0);
		// }

		
		// flag = CVodeInit(cvode_mem, rice_myofilament_ode, t_cvode, y);
		
		// CVodeSetStopTime(cvode_mem, t_cvode+dt);
		// oomph_info << "dt " << dt << std::endl;
		
		// oomph_info << NV_Ith_S(y,0) << std::endl;
		// oomph_info << NV_Ith_S(y,1) << std::endl;
		// oomph_info << NV_Ith_S(y,2) << std::endl;
		// oomph_info << NV_Ith_S(y,3) << std::endl;
		// oomph_info << NV_Ith_S(y,4) << std::endl;
		// oomph_info << NV_Ith_S(y,5) << std::endl;
		// oomph_info << NV_Ith_S(y,6) << std::endl;
		// oomph_info << NV_Ith_S(y,7) << std::endl;
		// oomph_info << NV_Ith_S(y,8) << std::endl;
		// oomph_info << NV_Ith_S(y,9) << std::endl;
		// oomph_info << NV_Ith_S(y,10) << std::endl;
		// oomph_info << NV_Ith_S(y,11) << std::endl;

		flag = CVode(cvode_mem, t_cvode+dt, y, &t_cvode, CV_NORMAL);

		if (flag != CV_SUCCESS) {
			oomph_info << "Error in CVode:" << flag << std::endl;
			exit(0);
		}

		state[N_NoXB_tt] = NV_Ith_S(y,0);
		state[N_tt] = NV_Ith_S(y,1);
		state[P_NoXB_tt] = NV_Ith_S(y,2);
		state[P_tt] = NV_Ith_S(y,3);
		state[XBprer_tt] = NV_Ith_S(y,4);
		state[XBpostr_tt] = NV_Ith_S(y,5);
		state[SL_tt] = NV_Ith_S(y,6);
		state[xXBpostr_tt] = NV_Ith_S(y,7);
		state[xXBprer_tt] = NV_Ith_S(y,8);
		state[TropCaL_tt] = NV_Ith_S(y,9);
		state[TropCaH_tt] = NV_Ith_S(y,10);
		state[intf0_tt] = NV_Ith_S(y,11);

		CVodeFree(&cvode_mem);
	}
	// else
	// {
	// 	state[N_NoXB_tt]	=	0.0;
	// 	state[N_tt]	=	0.0;
	// 	state[P_NoXB_tt]	=	0.0;
	// 	state[P_tt]	=	0.0;
	// 	state[XBprer_tt]	=	0.0;
	// 	state[XBpostr_tt]	=	0.0;
	// 	state[SL_tt]	=	0.0;
	// 	state[xXBpostr_tt]	=	0.0;
	// 	state[xXBprer_tt]	=	0.0;
	// 	state[TropCaL_tt]	=	0.0;
	// 	state[TropCaH_tt]	=	0.0;
	// 	state[intf0_tt]	=	0.0;
	// 	TropToT_dt = 0.0;
	// }







	Ek=RTONF*(log((Ko/state[Ki_tt])));
	Ena=RTONF*(log((Nao/state[Nai_tt])));
	Eks=RTONF*(log((Ko+pKNa*Nao)/(state[Ki_tt]+pKNa*state[Nai_tt])));
	
	
	Eca=0.5*RTONF*(log((Cao/state[Cai_tt])));




	Ak1=0.1/(1.+exp(0.06*(state[vm_tt]-Ek-200)));
	Bk1=(3.*exp(0.0002*(state[vm_tt]-Ek+100))+exp(0.1*(state[vm_tt]-Ek-10)))/(1.+exp(-0.5*(state[vm_tt]-Ek)));
	rec_iK1=Ak1/(Ak1+Bk1);
	rec_iNaK=(1./(1.+0.1245*exp(-0.1*state[vm_tt]*F/(R*T))+0.0353*exp(-state[vm_tt]*F/(R*T))));
	rec_ipK=1./(1.+exp((25-state[vm_tt])/5.98));




	INa=GNa*state[sm_tt]*state[sm_tt]*state[sm_tt]*state[sh_tt]*state[sj_tt]*(state[vm_tt]-Ena);



	ICaL=GCaL*state[sd_tt]*state[sf_tt]*state[sf2_tt]*state[sfcass_tt]*4.0*(state[vm_tt]-15.0)*(F*F/(R*T))*(0.25*exp(2.0*(state[vm_tt]-15.0)*F/(R*T))*state[CaSS_tt]-Cao)/(exp(2.0*(state[vm_tt]-15.0)*F/(R*T))-1.0);




	Ito=Gto*state[sr_tt]*state[ss_tt]*(state[vm_tt]-Ek);




	// IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(state[vm_tt]-Ek);



	epi_factor      = 1.8	*	m_epiMidRatio;
	endo_factor     = 1.8; 
	mcell_factor    = 1.8; 
	
	if (m_mutant) {
		sqt1_a1 = 			2.172;
		sqt1_b1 = 0.5*		1.077;
		sqt1_a2 = 0.3*		0.00655*exp(0.05547153*(state[vm_tt]-36. + 15.));
		sqt1_a =  			0.00555*exp(0.05547153*(state[vm_tt]-12. + 15.));
		sqt1_b =  			0.002357*exp(-0.036588*(state[vm_tt]));
		sqt1_b2 = 0.00077*	0.0029357*exp(1.3*	3.3*-0.02158*(state[vm_tt]));
		sqt1_ai = 			0.439*exp(-0.02352*(state[vm_tt]+25. + 15.))*(4.5/Ko);
		sqt1_bi = 0.025*		0.656*exp(0.000942*(state[vm_tt] + 15.))*((pow((4.5/Ko),0.3)));
		sqt1_mu = 			(sqt1_ai*sqt1_b2*sqt1_a2)/(sqt1_a2*sqt1_bi);
	}
	else {
		sqt1_a1 = 			2.172;
		sqt1_b1 = 			1.077;
		sqt1_a2 = 			0.00655   * exp(0.5*			0.05547153*(state[vm_tt]-36.));
		sqt1_a  =  			0.00555   * exp(				0.05547153*(state[vm_tt]-12.));
		sqt1_b  =  			0.002357  * exp(				-0.036588*(state[vm_tt]));
		sqt1_b2 = 0.65*		0.0029357 * exp(0.69*			-0.02158*(state[vm_tt]));
		sqt1_ai = 0.11*		0.439     * exp(1.7*			-0.02352*(state[vm_tt]+25.))*(4.5/Ko);
		sqt1_bi = 0.4*		0.656     * exp(				0.000942*(state[vm_tt]))*((pow((4.5/Ko),0.3)));
		sqt1_mu = 			(sqt1_ai*sqt1_b2)/sqt1_bi;			
	}

	if (cell_type == 102)
		Gkr = 0.0135*pow(Ko, 0.59) * endo_factor;
	else if(cell_type == 101)
		Gkr = 0.0135*pow(Ko, 0.59) * mcell_factor;
	else			
		Gkr = 0.0135*pow(Ko,0.59) * epi_factor; 
	
    const double sqt1_dC3 = (sqt1_b * state[sqt1_C2_tt])-(sqt1_a * state[sqt1_C3_tt]);
    const double sqt1_dC2 = -((sqt1_b + sqt1_a1) * state[sqt1_C2_tt])+(sqt1_a * state[sqt1_C3_tt])+(sqt1_b1 * state[sqt1_C1_tt]);
    const double sqt1_dC1 = -((sqt1_b1 + sqt1_a2 + sqt1_a2) * state[sqt1_C1_tt]) + (sqt1_a1 * state[sqt1_C2_tt]) + (sqt1_b2 * state[sqt1_O_tt]) + (sqt1_mu * state[sqt1_I_tt]);
    const double sqt1_dO  =  -((sqt1_b2 + sqt1_bi) * state[sqt1_O_tt]) + (sqt1_a2 * state[sqt1_C1_tt]) + (sqt1_ai * state[sqt1_I_tt]);
    const double sqt1_dI  = -((sqt1_mu + sqt1_ai) * state[sqt1_I_tt]) + (sqt1_a2 * state[sqt1_C1_tt]) + (sqt1_bi * state[sqt1_O_tt]);
	
    state[sqt1_O_tt]  +=dt* sqt1_dO;
    state[sqt1_C1_tt] +=dt* sqt1_dC1;
    state[sqt1_C2_tt] +=dt* sqt1_dC2;
    state[sqt1_C3_tt] +=dt* sqt1_dC3;
    state[sqt1_I_tt]  +=dt* sqt1_dI;

	
	IKr = Gkr*state[sqt1_O_tt]*(state[vm_tt]-Ek);


	IKs=Gks*state[sxs_tt]*state[sxs_tt]*(state[vm_tt]-Eks);



	IK1=GK1*rec_iK1*(state[vm_tt]-Ek);


	INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
	(1./(1+ksat*exp((n-1)*state[vm_tt]*F/(R*T))))*
	(exp(n*state[vm_tt]*F/(R*T))*state[Nai_tt]*state[Nai_tt]*state[Nai_tt]*Cao-
	exp((n-1)*state[vm_tt]*F/(R*T))*Nao*Nao*Nao*state[Cai_tt]*2.5);


	INaK=knak*(Ko/(Ko+KmK))*(state[Nai_tt]/(state[Nai_tt]+KmNa))*rec_iNaK;


	IpCa=GpCa*state[Cai_tt]/(KpCa+state[Cai_tt]);


	IpK=GpK*rec_ipK*(state[vm_tt]-Ek);


	IbNa=GbNa*(state[vm_tt]-Ena);


	IbCa=GbCa*(state[vm_tt]-Eca);


	if       (cell_type == EPI)   GNaL = 0.0065;
	else if  (cell_type == ENDO)  GNaL = 0.0065;
	else if  (cell_type == MCELL) GNaL = 0.0095;
	
	INaL = GNaL*state[mNaL_tt]*state[mNaL_tt]*state[mNaL_tt]*state[hNaL_tt]*(state[vm_tt]-Ena);

	Gsac = 1;//0.025/185;//110.0;
	Esac = 1.0;//-20.0;
	lambda_sac = 2;//1.10;
	
	strain = (state[SL_tt] - SLrest)/SLrest;
	strain_half_maximal = 0.233;
	slope_factor = 0.0203;//0.018;
	open_probability = 1.0/( 1 + exp( -( (strain-strain_half_maximal)/slope_factor ) ) );
	//Isac = Gsac*(lambda_sac-1)*(state[vm_tt]-Esac);
	Isac = Gsac*open_probability*(state[vm_tt]-Esac);
	//Isac = 0.0;
	
	//Concentration derivatives	
	kCaSR=maxsr-((maxsr-minsr)/(1+(EC/state[CaSR_tt])*(EC/state[CaSR_tt]))); 
	k1=k1_/kCaSR;
	k2=k2_*kCaSR;
	dRR=k4*(1-state[sRR_tt])-k2*state[CaSS_tt]*state[sRR_tt];
	state[sRR_tt] +=dt* dRR;
	const double sOO=k1*state[CaSS_tt]*state[CaSS_tt]*state[sRR_tt]/(k3+k1*state[CaSS_tt]*state[CaSS_tt]);
	
	
	Irel=Vrel*sOO*(state[CaSR_tt]-state[CaSS_tt]);
	Ileak=Vleak*(state[CaSR_tt]-state[Cai_tt]);
	Iup=Vmaxup/(1.+((Kup*Kup)/(state[Cai_tt]*state[Cai_tt])));
	Ixfer=Vxfer*(state[CaSS_tt]-state[Cai_tt]);
	
	
	CaCSQN=Bufsr*state[CaSR_tt]/(state[CaSR_tt]+Kbufsr);
	dCaSR=dt*(Iup-Irel-Ileak);
	bjsr=Bufsr-CaCSQN-dCaSR-state[CaSR_tt]+Kbufsr;
	cjsr=Kbufsr*(CaCSQN+dCaSR+state[CaSR_tt]);
	state[CaSR_tt]=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	// state[CaSR_tt] = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	// double dCaSR = (Iup - Irel - Ileak)/(1.0  + (Kbufsr*Bufsr)/pow(state[CaSR_tt] + Kbufsr, 2.0));
	// Variable_Derivatives[CaSR_tt] = (Iup - Irel - Ileak)/(1.0  + (Kbufsr*Bufsr)/pow(state[CaSR_tt] + Kbufsr, 2.0));
	
	CaSSBuf=Bufss*state[CaSS_tt]/(state[CaSS_tt]+Kbufss);
	dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
	bcss=Bufss-CaSSBuf-dCaSS-state[CaSS_tt]+Kbufss;
	ccss=Kbufss*(CaSSBuf+dCaSS+state[CaSS_tt]);
	state[CaSS_tt]=(sqrt(bcss*bcss+4*ccss)-bcss)/2;
	// state[CaSS_tt] +=dt* (-state[CaSS_tt] + (sqrt(bcss*bcss+4*ccss)-bcss)/2)/dt_rushlarsen;
	// double dCaSS =(-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE))/(1.0 + (Bufss*Kbufss)/pow(state[CaSS_tt]+Kbufss, 2.0));
	// Variable_Derivatives[CaSS_tt] = (-Ixfer * (Vc / Vss) + Irel * (Vsr / Vss) + (-ICaL * inversevssF2 * CAPACITANCE))/(1.0 + (Bufss*Kbufss)/pow(state[CaSS_tt]+Kbufss, 2.0));
	
	CaBuf=Bufc*state[Cai_tt]/(state[Cai_tt]+Kbufc);

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


	dCai=dt*(((-(IbCa+IpCa-2*INaCa+(Isac/3.0))*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer) - TropToT_dt/1000.0);

	// oomph_info << "Bufc " << Bufc << std::endl;
	// oomph_info << "CaBuf " << CaBuf << std::endl;
	// oomph_info << "dCai " << dCai << std::endl;
	// oomph_info << "state[Cai_tt] " << state[Cai_tt] << std::endl;
	// oomph_info << "Kbufc " << Kbufc << std::endl;

	const double bc=Bufc-CaBuf-dCai-state[Cai_tt]+Kbufc;
	cc=Kbufc*(CaBuf+dCai+state[Cai_tt]);
	// state[Cai_tt]=(sqrt(bc*bc+4*cc)-bc)/2;
	
	// oomph_info << state[Cai_tt] << std::endl;
	// oomph_info << bc << std::endl;
	// oomph_info << cc << std::endl;

	state[Cai_tt] = (sqrt(bc*bc+4*cc)-bc)/2;
	// state[Cai_tt] +=dt* (-state[Cai_tt] + (sqrt(bc*bc+4*cc)-bc)/2)/dt_rushlarsen;
	// oomph_info << Variable_Derivatives[Cai_tt] << std::endl;
	// double dCai = ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer)/(1.0 + (Bufc*Kbufc)/pow(state[Cai_tt] + Kbufc, 2.0));
	// Variable_Derivatives[Cai_tt] = ((-(IbCa + IpCa - 2 * INaCa) * inverseVcF2 * CAPACITANCE) - (Iup - Ileak) * (Vsr / Vc) + Ixfer)/(1.0 + (Bufc*Kbufc)/pow(state[Cai_tt] + Kbufc, 2.0));
    
	
	dNai=-(INa+IbNa+3*INaK+3*INaCa+(Isac/3.0))*inverseVcF*CAPACITANCE;
	state[Nai_tt] +=dt* dNai;
	
	dKi=-(this->get_stimulus(t)+IK1+Ito+IKr+IKs-2*INaK+IpK+(Isac/3.0))*inverseVcF*CAPACITANCE;
	state[Ki_tt] +=dt* dKi;

	if(myo)
	{
		state[SL_tt] = NV_Ith_S(y,6);
	}

			
	AM=1./(1.+exp((-60.-state[vm_tt])/5.));
	BM=0.1/(1.+exp((state[vm_tt]+35.)/5.))+0.10/(1.+exp((state[vm_tt]-50.)/200.));
	TAU_M=AM*BM;
	M_INF=1./((1.+exp((-56.86-state[vm_tt])/9.03))*(1.+exp((-56.86-state[vm_tt])/9.03)));
	if (state[vm_tt]>=-40.)
	{
		AH_1=0.; 
		BH_1=(0.77/(0.13*(1.+exp(-(state[vm_tt]+10.66)/11.1))));
		TAU_H= 1.0/(AH_1+BH_1);
	}
	else
	{
		AH_2=(0.057*exp(-(state[vm_tt]+80.)/6.8));
		BH_2=(2.7*exp(0.079*state[vm_tt])+(3.1e5)*exp(0.3485*state[vm_tt]));
		TAU_H=1.0/(AH_2+BH_2);
	}
	H_INF=1./((1.+exp((state[vm_tt]+71.55)/7.43))*(1.+exp((state[vm_tt]+71.55)/7.43)));
	if(state[vm_tt]>=-40.)
	{
		AJ_1=0.;      
		BJ_1=(0.6*exp((0.057)*state[vm_tt])/(1.+exp(-0.1*(state[vm_tt]+32.))));
		TAU_J= 1.0/(AJ_1+BJ_1);
	}
	else
	{
		AJ_2=(((-2.5428e4)*exp(0.2444*state[vm_tt])-(6.948e-6)*
			   exp(-0.04391*state[vm_tt]))*(state[vm_tt]+37.78)/
			  (1.+exp(0.311*(state[vm_tt]+79.23))));    
		BJ_2=(0.02424*exp(-0.01052*state[vm_tt])/(1.+exp(-0.1378*(state[vm_tt]+40.14))));
		TAU_J= 1.0/(AJ_2+BJ_2);
	}
	J_INF=H_INF;
	
	Xr1_INF=1./(1.+exp((-26.-state[vm_tt])/7.));
	axr1=450./(1.+exp((-45.-state[vm_tt])/10.));
	bxr1=6./(1.+exp((state[vm_tt]-(-30.))/11.5));
	TAU_Xr1=axr1*bxr1;
	Xr2_INF=1./(1.+exp((state[vm_tt]-(-88.))/24.));
	axr2=3./(1.+exp((-60.-state[vm_tt])/20.));
	bxr2=1.12/(1.+exp((state[vm_tt]-60.)/20.));
	TAU_Xr2=axr2*bxr2;
	
	Xs_INF=1./(1.+exp((-5.-state[vm_tt])/14.));
	Axs=(1400./(sqrt(1.+exp((5.-state[vm_tt])/6))));
	Bxs=(1./(1.+exp((state[vm_tt]-35.)/15.)));
	TAU_Xs=Axs*Bxs+80;
	
	if (cell_type == 100) {
		R_INF=1./(1.+exp((20-state[vm_tt])/6.));
		S_INF=1./(1.+exp((state[vm_tt]+20)/5.));
		TAU_R=9.5*exp(-(state[vm_tt]+40.)*(state[vm_tt]+40.)/1800.)+0.8;
		TAU_S=85.*exp(-(state[vm_tt]+45.)*(state[vm_tt]+45.)/320.)+5./(1.+exp((state[vm_tt]-20.)/5.))+3.;
	}
	else if (cell_type == 102) {
		R_INF=1./(1.+exp((20-state[vm_tt])/6.));
		S_INF=1./(1.+exp((state[vm_tt]+28)/5.));
		TAU_R=9.5*exp(-(state[vm_tt]+40.)*(state[vm_tt]+40.)/1800.)+0.8;
		TAU_S=1000.*exp(-(state[vm_tt]+67)*(state[vm_tt]+67)/1000.)+8.;
	}
	else {
		R_INF=1./(1.+exp((20-state[vm_tt])/6.));
		S_INF=1./(1.+exp((state[vm_tt]+20)/5.));
		TAU_R=9.5*exp(-(state[vm_tt]+40.)*(state[vm_tt]+40.)/1800.)+0.8;
		TAU_S=85.*exp(-(state[vm_tt]+45.)*(state[vm_tt]+45.)/320.)+5./(1.+exp((state[vm_tt]-20.)/5.))+3.;
	}		
	
	D_INF=1./(1.+exp((-8-state[vm_tt])/7.5));
	Ad=1.4/(1.+exp((-35-state[vm_tt])/13))+0.25;
	Bd=1.4/(1.+exp((state[vm_tt]+5)/5));
	Cd=1./(1.+exp((50-state[vm_tt])/20));
	TAU_D=Ad*Bd+Cd;
	F_INF=1./(1.+exp((state[vm_tt]+20)/7));
	Af=1102.5*exp(-(state[vm_tt]+27)*(state[vm_tt]+27)/225);
	Bf=200./(1+exp((13-state[vm_tt])/10.));
	Cf=(180./(1+exp((state[vm_tt]+30)/10)))+20;
	TAU_F=Af+Bf+Cf;
	F2_INF=0.67/(1.+exp((state[vm_tt]+35)/7))+0.33;
	Af2=600*exp(-(state[vm_tt]+25)*(state[vm_tt]+25)/170);
	Bf2=31/(1.+exp((25-state[vm_tt])/10));
	Cf2=16/(1.+exp((state[vm_tt]+30)/10));
	TAU_F2=Af2+Bf2+Cf2;
	FCaSS_INF=0.6/(1+(state[CaSS_tt]/0.05)*(state[CaSS_tt]/0.05))+0.4;
	TAU_FCaSS=80./(1+(state[CaSS_tt]/0.05)*(state[CaSS_tt]/0.05))+2.;
	
	//INaL Area
	alpha_mNaL = (0.32*(state[vm_tt]+47.13))/(1-std::exp(-0.1*(state[vm_tt]+47.13)));
	beta_mNaL  = 0.08*std::exp(-state[vm_tt]/11.0);
	mNaL_INF   = alpha_mNaL/(alpha_mNaL + beta_mNaL);
	TAU_mNaL   = 1.0/(alpha_mNaL + beta_mNaL); //TAU_hNaL initialised in constructor.
	hNaL_INF   = 1.0/(1+std::exp((state[vm_tt]+91.0)/6.1));


	// state[sm_tt] 	+=dt* ( -state[sm_tt]			+	M_INF-(M_INF-state[sm_tt])*exp(-dt_rushlarsen/TAU_M))/dt_rushlarsen;
	// state[sh_tt] 	+=dt* ( -state[sh_tt]			+	H_INF-(H_INF-state[sh_tt])*exp(-dt_rushlarsen/TAU_H))/dt_rushlarsen;
	// state[sj_tt] 	+=dt* ( -state[sj_tt]			+	J_INF-(J_INF-state[sj_tt])*exp(-dt_rushlarsen/TAU_J))/dt_rushlarsen;
	// // Variable_Derivatives[sxr1_tt] 	=( -sxr1		+	Xr1_INF-(Xr1_INF-sxr1)*exp(-dt_rushlarsen/TAU_Xr1))/dt_rushlarsen;
	// // Variable_Derivatives[sxr2_tt] 	=( -sxr2		+	Xr2_INF-(Xr2_INF-sxr2)*exp(-dt_rushlarsen/TAU_Xr2))/dt_rushlarsen;
	// state[sxs_tt] 	+=dt* ( -state[sxs_tt]		+	Xs_INF-(Xs_INF-state[sxs_tt])*exp(-dt_rushlarsen/TAU_Xs))/dt_rushlarsen;
	// state[ss_tt] 	+=dt* ( -state[ss_tt]			+	S_INF-(S_INF-state[ss_tt])*exp(-dt_rushlarsen/TAU_S))/dt_rushlarsen;
	// state[sr_tt] 	+=dt* ( -state[sr_tt]			+	R_INF-(R_INF-state[sr_tt])*exp(-dt_rushlarsen/TAU_R))/dt_rushlarsen;
	// state[sd_tt] 	+=dt* ( -state[sd_tt]			+	D_INF-(D_INF-state[sd_tt])*exp(-dt_rushlarsen/TAU_D))/dt_rushlarsen;
	// state[sf_tt] 	+=dt* ( -state[sf_tt]			+	F_INF-(F_INF-state[sf_tt])*exp(-dt_rushlarsen/TAU_F))/dt_rushlarsen;
	// state[sf2_tt] 	+=dt* ( -state[sf2_tt]		+	F2_INF-(F2_INF-state[sf2_tt])*exp(-dt_rushlarsen/TAU_F2))/dt_rushlarsen;
	// state[sfcass_tt] +=dt* ( -state[sfcass_tt]		+	FCaSS_INF-(FCaSS_INF-state[sfcass_tt])*exp(-dt_rushlarsen/TAU_FCaSS))/dt_rushlarsen;
	// state[mNaL_tt] 	+=dt* (-state[mNaL_tt]		+	mNaL_INF-(mNaL_INF - state[mNaL_tt])*exp(-dt_rushlarsen/TAU_mNaL))/dt_rushlarsen;
	// state[hNaL_tt] 	+=dt* (-state[hNaL_tt]		+	hNaL_INF-(hNaL_INF - state[hNaL_tt])*exp(-dt_rushlarsen/TAU_hNaL))/dt_rushlarsen;

	state[sm_tt] 	= M_INF-(M_INF-state[sm_tt])*exp(-dt/TAU_M);
	state[sh_tt] 	= H_INF-(H_INF-state[sh_tt])*exp(-dt/TAU_H);
	state[sj_tt] 	= J_INF-(J_INF-state[sj_tt])*exp(-dt/TAU_J);
	// Variable_Derivatives[sxr1_tt] 	=( -sxr1		+	Xr1_INF-(Xr1_INF-sxr1)*exp(-dt_rushlarsen/TAU_Xr1))/dt_rushlarsen;
	// Variable_Derivatives[sxr2_tt] 	=( -sxr2		+	Xr2_INF-(Xr2_INF-sxr2)*exp(-dt_rushlarsen/TAU_Xr2))/dt_rushlarsen;
	state[sxs_tt] 	= Xs_INF-(Xs_INF-state[sxs_tt])*exp(-dt/TAU_Xs);
	state[ss_tt] 	= S_INF-(S_INF-state[ss_tt])*exp(-dt/TAU_S);
	state[sr_tt] 	= R_INF-(R_INF-state[sr_tt])*exp(-dt/TAU_R);
	state[sd_tt] 	= D_INF-(D_INF-state[sd_tt])*exp(-dt/TAU_D);
	state[sf_tt] 	= F_INF-(F_INF-state[sf_tt])*exp(-dt/TAU_F);
	state[sf2_tt] 	= F2_INF-(F2_INF-state[sf2_tt])*exp(-dt/TAU_F2);
	state[sfcass_tt] = FCaSS_INF-(FCaSS_INF-state[sfcass_tt])*exp(-dt/TAU_FCaSS);
	state[mNaL_tt] 	= mNaL_INF-(mNaL_INF - state[mNaL_tt])*exp(-dt/TAU_mNaL);
	state[hNaL_tt] 	= hNaL_INF-(hNaL_INF - state[hNaL_tt])*exp(-dt/TAU_hNaL);


	//Calculate total ionic current
	
	// oomph_info << IKr << std::endl;
	// oomph_info << IKs << std::endl;
	// oomph_info << IK1 << std::endl;
	// oomph_info << Ito << std::endl;
	// oomph_info << INa << std::endl;
	// oomph_info << IbNa << std::endl;
	// oomph_info << ICaL << std::endl;
	// oomph_info << IbCa << std::endl;
	// oomph_info << INaK << std::endl;
	// oomph_info << INaCa << std::endl;
	// oomph_info << IpCa << std::endl;
	// oomph_info << IpK << std::endl;
	// oomph_info << INaL << std::endl;
	// oomph_info << this->get_stimulus(t) << std::endl;
	// oomph_info << Isac << std::endl;

	state[vm_tt] -= dt*(IKr + IKs + IK1 + Ito +	INa + IbNa + ICaL +	IbCa + INaK + INaCa + IpCa + IpK + INaL + this->get_stimulus(t) + Isac);


	// oomph_info << state[vm_tt] << std::endl;
	// oomph_info << state[sm_tt] << std::endl;
	// oomph_info << state[sh_tt] << std::endl;
	// oomph_info << state[sj_tt] << std::endl;
	// oomph_info << state[sqt1_O_tt] << std::endl;
	// oomph_info << state[sqt1_C1_tt] << std::endl;
	// oomph_info << state[sqt1_C2_tt] << std::endl;
	// oomph_info << state[sqt1_C3_tt] << std::endl;
	// oomph_info << state[sqt1_I_tt] << std::endl;
	// oomph_info << state[sxs_tt] << std::endl;
	// oomph_info << state[ss_tt] << std::endl;
	// oomph_info << state[sr_tt] << std::endl;
	// oomph_info << state[sd_tt] << std::endl;
	// oomph_info << state[sf_tt] << std::endl;
	// oomph_info << state[sf2_tt] << std::endl;
	// oomph_info << state[sfcass_tt] << std::endl;
	// oomph_info << state[sRR_tt] << std::endl;
	// oomph_info << state[Cai_tt] << std::endl;
	// oomph_info << state[CaSR_tt] << std::endl;
	// oomph_info << state[CaSS_tt] << std::endl;
	// oomph_info << state[mNaL_tt] << std::endl;
	// oomph_info << state[hNaL_tt] << std::endl;
	// oomph_info << state[Nai_tt] << std::endl;
	// oomph_info << state[Ki_tt] << std::endl;
	// oomph_info << state[N_NoXB_tt] << std::endl;
	// oomph_info << state[N_tt] << std::endl;
	// oomph_info << state[P_NoXB_tt] << std::endl;
	// oomph_info << state[P_tt] << std::endl;
	// oomph_info << state[XBprer_tt] << std::endl;
	// oomph_info << state[XBpostr_tt] << std::endl;
	// oomph_info << state[SL_tt] << std::endl;
	// oomph_info << state[xXBpostr_tt] << std::endl;
	// oomph_info << state[xXBprer_tt] << std::endl;
	// oomph_info << state[TropCaL_tt] << std::endl;
	// oomph_info << state[TropCaH_tt] << std::endl;
	// oomph_info << state[intf0_tt] << std::endl;
}



void IsmailTNNP06RushLarsen::get_output(double *state, double *out)
{
	out[0] = (state[SL_tt] - SLrest)/SLrest;
}

double IsmailTNNP06RushLarsen::GetActiveStrain(double *state)
{
	return (state[SL_tt] - SLrest)/SLrest;
}



}; //End namespace