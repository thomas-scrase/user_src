#ifndef OOMPH_ISMAIL_TNNP06_RUSHLARSEN_VENT_HEADER
#define OOMPH_ISMAIL_TNNP06_RUSHLARSEN_VENT_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif
	
#include "../cell_models/cell_model_base.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_band.h>        /* prototype for CVBand */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace oomph{

// namespace ISMAIL_TT
// {
// 	template <typename T> int sgn( T val )
// 	{
// 		return (val > T(0)) - (val < T(0));
// 	}
// };

class IsmailTNNP06RushLarsen : public CellModelBaseFullyPartitioned
{
public:
	IsmailTNNP06RushLarsen(const unsigned& number_of_backup_values);

	~IsmailTNNP06RushLarsen();

	std::string get_cell_model_name() {return "IsmailTNNP06RushLarsen";}

protected:
	double get_initial_state_variable(const unsigned &v);

	double GetActiveStrain(double *state);

	void TakeTimestep(const double& dt, const double& t, double* state);

	unsigned index_of_membrane_potential_in_cell_data(){return vm_tt;}

	virtual void get_output(double *state, double *out);
public:
	void set_cell_type(const unsigned &c)
	{
		if((c==EPI) || (c==MCELL) || (c==ENDO))
		{
			cell_type = c;
		}
		else
		{
			throw OomphLibError("Cell type is not compatible with model",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
		}
	}

	enum Cell_Variables_Enum : unsigned
	{
		vm_tt,
		sm_tt,
		sh_tt,
		sj_tt,
		
		// sxr1_tt,
		// sxr2_tt,
		
		sqt1_O_tt,
		sqt1_C1_tt,
		sqt1_C2_tt,
		sqt1_C3_tt,
		sqt1_I_tt,

		sxs_tt,
		ss_tt,
		sr_tt,
		sd_tt,
		sf_tt,
		sf2_tt,
		sfcass_tt,
		sRR_tt,
		// sOO_tt,
		Cai_tt,
		CaSR_tt,
		CaSS_tt,
		mNaL_tt,
		hNaL_tt,
		Nai_tt,
		Ki_tt,
		N_NoXB_tt,
		N_tt,
		P_NoXB_tt,
		P_tt,
		XBprer_tt,
		XBpostr_tt,
		SL_tt,
		xXBpostr_tt,
		xXBprer_tt,
		TropCaL_tt,
		TropCaH_tt,
		intf0_tt
	};

	enum TypeCell : unsigned
	{
		EPI,
		MCELL,
		ENDO
	};
	

	// void rice_myofilament_ode( const Boost_State_Type &x , Boost_State_Type &dxdt , const double t )
	// {
	// 	if(myo == false)
	// 	{
	// 			dxdt[0]	= 0.0;
	// 			dxdt[1]	= 0.0;
	// 			dxdt[2]	= 0.0;
	// 			dxdt[3]	= 0.0;
	// 			dxdt[4]	= 0.0;
	// 			dxdt[5]	= 0.0;
	// 			dxdt[6]	= 0.0;  
	// 			dxdt[7]	= 0.0; 
	// 			dxdt[8]	= 0.0;
	// 			dxdt[9]	= 0.0;
	// 			dxdt[10]	= 0.0;
	// 			dxdt[11]	= 0.0;
	// 			return;
	// 		}

	// 		/* Handle Ca binding to troponin here */
	// 		//  Adjust for temperature
	// 		konT		= kon	*std::pow(Qkon,(T-273.0-37.0)/10.);
	// 		koffLT	= koffL	*std::pow(Qkoff,(T-273.0-37.0)/10.);
	// 		koffHT	= koffH	*std::pow(Qkoff,(T-273.0-37.0)/10.);
			
	// 		/* Compute derived regulatory unit rates, Ca in in nM to prevent Ca overflow (max ~ 10 in auto), TropCaL and TropCaH are normailized so no Units. */
	// 	  	dxdt[9] = konT * Cai_rice*1000 * (1.0 - x[9]) - koffLT * x[9];                                              
	// 	  	dxdt[10] = konT * Cai_rice*1000 * (1.0 - x[10]) - koffHT * x[10]; 
			
	// 		/* Compute rates for N to P transitions */
	// 		// Compute z-line end of single overlap region
	// 		sovr_ze = std::min(len_thick/2.,x[6]/2.);
			
	// 		// Compute centerline of end of single overlap region
	// 	    sovr_cle = std::max(x[6]/2.-(x[6]-len_thin),len_hbare/2.);
			
	// 		// Compute length of the single overlap 
	// 		len_sovr = sovr_ze-sovr_cle;
			
	// 		// Compute overlap fraction for thick filament
	// 		SOVFThick = len_sovr*2./(len_thick-len_hbare);
			
	// 		// Compute overlap fraction for thin filament
	// 		SOVFThin = len_sovr/len_thin;
			
	// 		/* Compute combined Ca binding to high (w/XB) and low (w/o XB)  No Units. */
	// 		perm = (1.0 - SOVFThin) * x[9] + SOVFThin * x[10];
			
	// 		/* Set perm variable that control n to p transiton in Hill-like fashion, No Units. */
	// 		permtot = std::pow((1.0 / (1.0 + std::pow((perm50 / perm), nperm))), 0.5);
	// 	 	if (100.0 < 1.0/permtot)
	// 			inprmt=100.0;
	// 	  	else 
	// 			inprmt = 1.0/permtot;
			
	// 		/* Adjust for Ca activation level and temp. */
	// 		kn_pT = kn_p*permtot*std::pow(Qkn_p,(T-273.0-37.0)/10.);
	// 		kp_nT = kp_n*inprmt*std::pow(Qkp_n,(T-273.0-37.0)/10.);
			
	// 		/* Compute fapp, the attachment rate from weak to strong, pre-rotated state */
	// 	    fappT = fapp*std::pow(Qfapp,(T-273.0-37.0)/10.);
			
	// 		/* Compute gapp, the detachment rate from strong, pre-rotated to the weakly-bound state */
	// 		/* Compute SL modifier for gapp to increase rate st shorter SL */
	// 	  	gapslmd = 1.0 + (1.0 - SOVFThick)*gslmod;
	// 	    gappT = gapp*gapslmd*std::pow(Qgapp,(T-273.0-37.0)/10.);
			
	// 		/* Set rate for rates between pre-force and force states*/
	// 		/* Compute modifiers based on mean strain of states. */
	// 	  	hfmd=std::exp(-sign(x[8])*hfmdc*((x[8]/x_0)*(x[8]/x_0)));
	// 	  	hbmd=std::exp(sign((x[7]-x_0))*hbmdc*(((x[7]-x_0)/x_0)*((x[7]-x_0)/x_0)));
			
	// 		/* Combine modifiers of hf and hb */
	// 		hfT=hf*hfmd*std::pow(Qhf,(T-273.0-37.0)/10.);
	// 		hbT=hb*hbmd*std::pow(Qhb,(T-273.0-37.0)/10.); 
			
	// 		/* Set rate for rates gxb, the ATP using XB transition */
	// 		/* Add term for distortion dependence of gxb */
	// 	    gxbmd = heav(x_0-x[7])*std::exp(sigmap*((x_0-x[7])/x_0)*((x_0-x[7])/x_0))+(1.-heav(x_0-x[7])*std::exp(sigman*(((x[7]-x_0)/x_0)*(x[7]-x_0)/x_0)));
	// 	    gxbT = gxb*gxbmd*std::pow(Qgxb,(T-273.0-37.0)/10.);
			
	// 		/* update all RUs */
	// 		dxdt[0] = -kn_pT * x[0] + kp_nT * x[2];
	// 		dxdt[2] = -kp_nT * x[2] + kn_pT * x[0];
			
	// 		dxdt[1] = -kn_pT * x[1] + kp_nT * x[3]; 
	// 		dxdt[3] = -kp_nT * x[3] + kn_pT * x[1] - fappT * x[3] + gappT * x[4] + gxbT * x[5]; 
	// 		dxdt[4] = fappT * x[3] - gappT * x[4] - hfT * x[4] + hbT * x[5];
	// 		dxdt[5] = hfT * x[4] - hbT * x[5] - gxbT * x[5];
			
	// 		/* compute steady-state fractions in XBprer and XBpostr using King-Altman rule */
	// 		SSXBprer = (hb*fapp+gxb*fapp)/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
	// 		SSXBpostr = fapp*hf/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
	// 		/* compute normalization for scaling active and passive force */
	// 		Fnordv = kxb*x_0*SSXBpostr;
			
	// 		/* compute force here */
	// 		force = kxb*SOVFThick*(x[7]*x[5]+x[8]*x[4]);
			
	// 		/* Compute SL, Set to 0 for isometric, to 1 for active contraction */
	// 		/* Note that passive force is specified in terms of maximal active force */
	// 		/* read passive force from table file specified in terms of maximal active force */
	// 		/* compyet derivative of SL */
			
	// 		// Note that passive force is specified in terms of maximal active force.
	// 		// Passive force is computed in two ways depending on if cell is skinned or intact.
	// 		ppforce=sign((x[6]-SLrest))*PCon_t*(std::exp(PExp_t*std::abs((x[6]-SLrest)))-1);
	// 	  //   if (singlecell == false) 
	// 			// ppforce += std::max(0.0,(PCon_col*std::exp(PExp_col*(x[6]-SLcol))));
			
	// 		// Compute afterload either as a contant or from the series elastic element
	// 	 	afterload = 0.0;
			
	// 	 //    if (SEon == true) {
	// 		// 	afterload = KSE*(SLset-x[6]);
	// 		// 	if ((SEon_LengthClamp==true)&&(index==(pulse_number-1))) {
	// 		// 		afterload=5000.0*(SLset-x[6]);
	// 		// 	}
	// 	 //    } 
	// 		// else {
	// 			afterload=0.0; 
	// 	    // }
			
	// 		//  Compute the integral of forces 
	// 		dxdt[11]=(-ppforce+PreloadF+(-force/Fnordv)+afterload);
			
	// 		// Compute the derivative of SL for an isolated muscle
	// 	  	dSLisolated = ((x[11] + (SLset-x[6])*visc)/massf )*heav(x[6]-SLmin)*heav(SLmax-x[6]);
	// 	 	// if (contraction==false) 
	// 			// dSLisolated = 0.;
			
			
	// 		// //========>>>>>>>>
	// 		// /*  The follow section of code implements is used to simulate ejection 
	// 		//  from canine whole heart into a Windkessel model of the circulatory system.
	// 		//  The complete description is available in:
	// 		//  De Tombe, PP and Little, WC. Inotropic effects of ejection 
	// 		//  are myocardial properties. Am J Physiol. 1994 Mar;266(3 Pt 2):H1202-13. */
			
	// 		// // Betavl is unitless
	// 		// Betavl=M_PI*std::pow(0.75*M_PI,0.6667);
			
	// 	 //    // Compute A_0 which is the cross-sectional area of the ventricular wall
	// 	 //    // A_0 has unit of area (volume^0.6667) in cm^2
	// 		// A_0 = Betavl*( std::pow(Vref + Vwall,0.6667) - std::pow(Vref,0.6667) );
			
	// 	 //    // Compute SigmaP from normalized force
	// 	 // 	SigPrs=SigFmx*(ppforce + force/Fnordv);
			
	// 	 //    // Compute left venticular pressure in mmHg
	// 	 // 	LVpress=SigPrs*A_0/(1.36*Betavl*std::pow(x[13],0.667));
			
	// 		// // Compute Flows - assumes aortic valve so flow in unidirectional
	// 		// // Units are l/s
	// 		// iflow=heav(LVpress-x[12])*((LVpress-x[12])/Rchar);
	// 		// // Allow active ejection on last beat only 
	// 		// if (index<(pulse_number-1))
	// 		// 	iflow=0.0;
	// 		// if (Isovolume==true) 
	// 		// 	iflow=0.0;
			
	// 	 //    // Compute the derivative of arterial pressure
	// 		// dxdt[12]=(iflow/Cwind - (x[12]-ARTpset)/(Cwind*Ra));
	// 	 //    if (WholeHeart==false) 
	// 		// 	dxdt[12]=0.0;
		    
	// 		// // Compute the  derivative of the change in Left Ventricular Volume
	// 		// // Do not allow volume less than 0
	// 	 //    dxdt[13]=-(iflow)*heav(x[13]);
	// 	 //    if (WholeHeart==false) 
	// 		// 	dxdt[13]=0.0;
			
	// 		// // Compute the  derivative of the change in SL
	// 	 //    dSLwholeheart=SLref*(1/std::pow(Vref+0.333*Vwall,0.333))*(dxdt[13]/3)/std::pow(x[13]+0.333*Vwall,0.667);
			
	// 		// // Compute the derivative of SL based using cell/tissue or simulated whole heart model
	// 		// if (WholeHeart==true) {
	// 		// 	dxdt[6] = dSLwholeheart;
	// 		// } else {
	// 			dxdt[6] = dSLisolated;
	// 		// }
			
	// 		// //========>>>>>>>>
			
			
			
			
	// 		dSL=dxdt[6];
			
	// 		// Compute the duty fractions using King-Alman Rule
	// 		// Compute for states XBpref and XBf
	// 		dtyf_prer=(hbT*fappT+gxbT*fappT)/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
	// 		dtyf_postr=fappT*hfT/(fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
			
	// 		// Compute mean strain by consider two competing effects - SL motion and XB cycling 
	// 		if(myo == true) {
	// 			dxdt[8] = dSL/2. + xPsi*(1/dtyf_prer)*(-x[8]*fappT+(x[7]-x_0-x[8])*hbT);
	// 			dxdt[7] = dSL/2. + xPsi*(1/dtyf_postr)*(x_0+x[8]-x[7])*hfT;
	// 		} 
	// 		else {
	// 			dxdt[8]=0.0;
	// 			dxdt[7]=0.0;
	// 		}
			
	// 		// Compute derivative of z-line end of single overlap region
	// 		sovr_ze_dt=0.;
	// 		if(len_thick/2.>x[6]/2.) 
	// 			sovr_ze_dt=-dSL/2.;
			
	// 		// Compute derivative of centerline of end of single overlap region
	// 	    sovr_cle_dt=0.;
	// 	    if (x[6]/2.-(x[6]-len_thin)>len_hbare/2.)
	// 			sovr_cle_dt=-dSL/2;
			
	// 		// Compute the derivative of the length of the single overlap 
	// 		len_sovr_dt=sovr_ze_dt-sovr_cle_dt;
		 	
	// 		// Compute the derivative of the overlap fraction for thin filament
	// 		SOVFThin_dt = len_sovr_dt/len_thin;
			
	// 		// Compute the derivative of the overlap fraction for thick filament
	// 		SOVFThick_dt = len_sovr_dt*2./(len_thick-len_hbare);
			
	// 		/*Compute fraction of strongly bound for compute effective binding of Ca */
	// 		FrSBXB=((x[5]+x[4])/(SSXBpostr + SSXBprer));
	// 		FrSBXB_dt=((dxdt[5]+dxdt[4])/(SSXBpostr + SSXBprer));
			
	// 		TropToT=Trop_conc*((1.0-SOVFThin)*x[9] + SOVFThin*(FrSBXB*x[10]+(1.0-FrSBXB)*x[9]));
	// 		TropToT_dt=Trop_conc*(-1.0*SOVFThin_dt*x[9]+(1.0-SOVFThin)*dxdt[9] + SOVFThin_dt*(FrSBXB*x[10]+(1.0-FrSBXB)*x[9]) + SOVFThin*(FrSBXB_dt*x[10]+FrSBXB*dxdt[10]-FrSBXB_dt*x[9]+(1.0-FrSBXB)*dxdt[9]));

	// 		active_force = 10*force/Fnordv;   /* in units of mN/(mm^2) */
	// }


	// void rice_myofilament_ode( const Boost_State_Type &x , Boost_State_Type &dxdt , const double t )
	// {}
	
	// double sign(double a){ return ((a) < (0.) ? (-1.0) : (1.0));}
	// double heav(double a){ return ((a) < (0.) ? (0.0) : (1.0));}

	realtype sign(realtype a){ return ((a) < (0.) ? (-1.0) : (1.0));}
	realtype heav(realtype a){ return ((a) < (0.) ? (0.0) : (1.0));}


	double m_epiMidRatio;
	double m_mutant;
	// TypeCell m_celltype;
	unsigned cell_type;	

	//External concentrations
	double Ko;
	double Cao;
	double Nao;
	
	//Intracellular volumes
	double Vc;
	double Vsr;
	double Vss;
	
	//Calcium buffering dynamics
	double Bufc;
	double Kbufc;
	double Bufsr;
	double Kbufsr;
	double Bufss;
	double Kbufss;
	
	//Intracellular calcium flux dynamics
	double Vmaxup;
	double Kup;
	double Vrel;
	double k1_;
	double k2_;
	double k3;
	double k4;
	double EC;
	double maxsr;
	double minsr;
	double Vleak;
	double Vxfer;	
	
	//Constants
	double R;
	double F;
	double T;
	double RTONF;
	
	//Cellular capacitance         
	double CAPACITANCE;
	
	//Parameters for currents
	//Parameters for IKr
	double Gkr;
	
	//Parameters for Iks
	double pKNa;
	double Gks;
	
	//Parameters for Ik1
	double GK1;
	
	//Parameters for Ito
	double Gto;
	
	//Parameters for INa
	double GNa;
	
	//Parameters for IbNa
	double GbNa;
	
	//Parameters for INaK
	double KmK;
	double KmNa;
	double knak;
	
	//Parameters for ICaL
	double GCaL;
	
	//Parameters for IbCa
	double GbCa;
	
	//Parameters for INaCa
	double knaca;
	double KmNai;
	double KmCa;
	double ksat;
	double n;
	
	//Parameters for IpCa
	double GpCa;
	double KpCa;
	
	//Parameters for IpK;
	double GpK;

	double Cai_rice;
	
	
	//==========================
	// PARAMETER FOR INTEGRATION
	//==========================
	// double  m_HT; // time step
	
	//==================================
	// PARAMETERS FOR INITIAL CONDITIONS 
	//==================================
	//Initial values of state variables
	// double svolt;
	// double Cai;
	// double CaSR;
	// double CaSS;
	// double Nai;
	// double Ki;
	
	//==================================
	// PARAMETER FOR SIMULATION DURATION
	//==================================
	//duration of the simulation 
	//double STOPTIME=100000;
	
	//=====================================
	// PARAMETERS FOR STIMULATION PROTOCOLS 
	//=====================================	
	// double stimduration;
	// double stimstrength;
	// double tbegin;
	// double tend;
	// int counter;
	// double dia;
	// double basicperiod;
	// double basicapd;
	// int repeats;
	// double Istim;
	// double time;
	//double m_HT;
	
	//
	double IKr;
	double IKs;
	double IK1;
	double Ito;
	double INa;
	double IbNa;
	double ICaL;
	double IbCa;
	double INaCa;
	double IpCa;
	double IpK;
	double INaK;
	double Irel;
	double Ileak;
	double Iup;
	double Ixfer;
	double k1;
	double k2;
	double kCaSR;
	
	
	double dNai;
	double dKi;
	double dCai;
	double dCaSR;
	double dCaSS;
	double dRR;
	
	
	double Ek;
	double Ena;
	double Eks;
	double Eca;
	double CaCSQN;
	double bjsr;
	double cjsr;
	double CaSSBuf;
	double bcss;
	double ccss;
	double CaBuf;
	// double bc;
	double cc;
	double Ak1;
	double Bk1;
	double rec_iK1;
	double rec_ipK;
	double rec_iNaK;
	double AM;
	double BM;
	double AH_1;
	double BH_1;
	double AH_2;
	double BH_2;
	double AJ_1;
	double BJ_1;
	double AJ_2;
	double BJ_2;
	double M_INF;
	double H_INF;
	double J_INF;
	double TAU_M;
	double TAU_H;
	double TAU_J;
	double axr1;
	double bxr1;
	double axr2;
	double bxr2;
	double Xr1_INF;
	double Xr2_INF;
	double TAU_Xr1;
	double TAU_Xr2;
	double Axs;
	double Bxs;
	double Xs_INF;
	double TAU_Xs;
	double R_INF;
	double TAU_R;
	double S_INF;
	double TAU_S;
	double Ad;
	double Bd;
	double Cd;
	double Af;
	double Bf;
	double Cf;
	double Af2;
	double Bf2;
	double Cf2;
	double TAU_D;
	double D_INF;
	double TAU_F;
	double F_INF;
	double TAU_F2;
	double F2_INF;
	double TAU_FCaSS;
	double FCaSS_INF;
	
	
	double inverseVcF2;
	double inverseVcF;
	double inversevssF2;
	
	// double sm;
	// double sh;
	// double sj;
	// double sxr1;
	// double sxr2;
	// double sxs; 
	// double ss;  
	// double sr;
	// double sd;
	// double sf;
	// double sf2;
	// double sfcass;
	// double sRR;
	// double sOO;
	// double sItot;
	
	// CVODE
	int NEQ;
	double T0;		// Initial time
	double T1;		// First output time
	double TMULT;
	int NOUT;
	realtype rtol;	// Scalar relative tolerance
	realtype atol;	// Vector absolute tolerance
	double TOUT;
	double IOUT;
	SUNContext sunctx;
	SUNMatrix Jacobian;
	SUNLinearSolver LS;
	
	int flag, k;
	N_Vector y, abstol;
	void *cvode_mem;
	// double *data;
	// double *ropt;
	// long int *iopt;
	
	// Rice et al.
	//============
	// double N_NoXB;
	// double P_NoXB;
	// double N;
	// double P;
	// double XBprer;
	// double XBpostr;
	// double SL;
	// double xXBpostr;
	// double xXBprer;
	// double TropCaL;
	// double TropCaH; 
	// double intf0; 
	double kon, konT, koffL, koffLT, koffH, koffHT, Qkon, Qkoff;
	double sovr_ze, len_thick, sovr_cle, len_thin, len_hbare, len_sovr, SOVFThick, SOVFThin, perm, permtot, perm50, nperm, inprmt;
	double kn_p, kn_pT, Qkn_p, kp_n, kp_nT, Qkp_n, fapp, fappT, Qfapp, gapslmd, gslmod, gapp, gappT, Qgapp;
	double hfmd, hfmdc, x_0, hbmd, hbmdc, hfT, hf, Qhf, hbT, hb, Qhb;
	double gxbmd, gxbT, gxb, Qgxb, sigmap, sigman;
	double SSXBprer, SSXBpostr, Fnordv, kxb, force, active_force;
	double ppforce, PCon_t, PExp_t, SLrest, PCon_col, PExp_col, SLcol;
	bool /*singlecell,*//* SEon, SEon_LengthClamp,*/ myo /*,contraction, *//*Isovolume,*/ /*WholeHeart*/;
	double afterload, KSE, SLset, SLref, pulse_number, index, PreloadF;
	double dSLisolated, visc, massf, SLmin, SLmax, dSL, dtyf_prer, dtyf_postr, xPsi;
	double sovr_ze_dt, sovr_cle_dt, len_sovr_dt, SOVFThick_dt, SOVFThin_dt;
	double FrSBXB, FrSBXB_dt, TropToT, Trop_conc, TropToT_dt;
	double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13;
	double temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27;
	double Betavl, A_0, Vref, Vwall, SigPrs, SigFmx, LVpress, iflow, Rchar, Cwind, ARTpset, Ra, dSLwholeheart, ARTpress;
	double LVV;
	
	//==========================
	//Markov Formulation for IKr
	//==========================
	// double sqt1_O;									/* Markov Model - Open Probability for IKr */
	// double sqt1_C1;									/* Markov Model - C1 Probability for IKr */
	// double sqt1_C2;									/* Markov Model - C2 Probability for IKr */
	// double sqt1_C3;									/* Markov Model - C3 Probability for IKr */
	// double sqt1_I;									/* Markov Model - I Probability for IKr */
	
	// double sqt1_dO;									/* Markov Model - dOpen Probability for IKr */
	// double sqt1_dC1;								/* Markov Model - dC1 Probability for IKr */
	// double sqt1_dC2;								/* Markov Model - dC2 Probability for IKr */
	// double sqt1_dC3;								/* Markov Model - dC3 Probability for IKr */
	// double sqt1_dI;									/* Markov Model - dI Probability for IKr */
	
	double sqt1_a1;									/* C1->O or C1->I */
	double sqt1_a2;									/* C2->C1 */
	double sqt1_a;									/* C3->C2 */
	double sqt1_b;									/* C2->C3 */
	double sqt1_b1;									/* C1->C2 */
	double sqt1_b2;									/* O->C1 */
	double sqt1_ai;									/* I->O */
	double sqt1_bi;									/* O->I */
	double sqt1_mu;
	
	double epi_factor, endo_factor, mcell_factor;
	
	//=====
	// INaL
	//=====
	double INaL;
	double GNaL;
	// double mNaL;
	// double hNaL;
	double alpha_mNaL;
	double beta_mNaL;
	double mNaL_INF;
	double hNaL_INF;
	double TAU_mNaL;
	double TAU_hNaL;        //600ms
	
	//=========================
	//Stretch-activated Channel
	//=========================
	double Isac, Isac_Na, Isac_K, Isac_Ca;
	double Esac;
	double Gsac, gsac;
	double pNa, pK, pCa, zna, zk, zca;
	double Ksac, alpha_sac, lambda_sac;
	double strain, strain_half_maximal, open_probability, slope_factor;
};

// class IsmailTNNP06RushLarsen_ODE_Wrapper
// {
// private:
// 	IsmailTNNP06RushLarsen* m_obj_pt;
	
// public:
// 	IsmailTNNP06RushLarsen_ODE_Wrapper(IsmailTNNP06RushLarsen* M_Obj_pt)
// 	{
// 		m_obj_pt = M_Obj_pt;
// 	}	

// 	void operator()( const Boost_State_Type &x , Boost_State_Type &dxdt , const double t )
// 	{
// 		m_obj_pt->rice_myofilament_ode(x,dxdt,t);
// 	}
// };



static int rice_myofilament_ode(realtype t, N_Vector x, N_Vector dxdt, void *user_data)
	{
		IsmailTNNP06RushLarsen* p = static_cast<IsmailTNNP06RushLarsen*>(user_data);

		if(p->myo == false)
		{
				NV_Ith_S(dxdt,0)	= 0.0;
				NV_Ith_S(dxdt,1)	= 0.0;
				NV_Ith_S(dxdt,2)	= 0.0;
				NV_Ith_S(dxdt,3)	= 0.0;
				NV_Ith_S(dxdt,4)	= 0.0;
				NV_Ith_S(dxdt,5)	= 0.0;
				NV_Ith_S(dxdt,6)	= 0.0;  
				NV_Ith_S(dxdt,7)	= 0.0; 
				NV_Ith_S(dxdt,8)	= 0.0;
				NV_Ith_S(dxdt,9)	= 0.0;
				NV_Ith_S(dxdt,10)	= 0.0;
				NV_Ith_S(dxdt,11)	= 0.0;
				return 0;
			}

			/* Handle Ca binding to troponin here */
			//  Adjust for temperature
			p->konT		= p->kon	*std::pow(p->Qkon,(p->T-273.0-37.0)/10.);
			p->koffLT	= p->koffL	*std::pow(p->Qkoff,(p->T-273.0-37.0)/10.);
			p->koffHT	= p->koffH	*std::pow(p->Qkoff,(p->T-273.0-37.0)/10.);
			
			/* Compute derived regulatory unit rates, Ca in in nM to prevent Ca overflow (max ~ 10 in auto), TropCaL and TropCaH are normailized so no Units. */
		  	NV_Ith_S(dxdt,9) = p->konT * p->Cai_rice*1000 * (1.0 - NV_Ith_S(x,9)) - p->koffLT * NV_Ith_S(x,9);                                              
		  	NV_Ith_S(dxdt,10) = p->konT * p->Cai_rice*1000 * (1.0 - NV_Ith_S(x,10)) - p->koffHT * NV_Ith_S(x,10); 
			
			/* Compute rates for N to P transitions */
			// Compute z-line end of single overlap region
			p->sovr_ze = std::min(p->len_thick/2.,NV_Ith_S(x,6)/2.);
			
			// Compute centerline of end of single overlap region
		  p->sovr_cle = std::max(NV_Ith_S(x,6)/2.-(NV_Ith_S(x,6)-p->len_thin),p->len_hbare/2.);
			
			// Compute length of the single overlap 
			p->len_sovr = p->sovr_ze-p->sovr_cle;
			
			// Compute overlap fraction for thick filament
			p->SOVFThick = p->len_sovr*2./(p->len_thick-p->len_hbare);
			
			// Compute overlap fraction for thin filament
			p->SOVFThin = p->len_sovr/p->len_thin;
			
			/* Compute combined Ca binding to high (w/XB) and low (w/o XB)  No Units. */
			p->perm = (1.0 - p->SOVFThin) * NV_Ith_S(x,9) + p->SOVFThin * NV_Ith_S(x,10);
			
			/* Set perm variable that control n to p transiton in Hill-like fashion, No Units. */
			p->permtot = std::pow((1.0 / (1.0 + std::pow((p->perm50 / p->perm), p->nperm))), 0.5);
		 	if (100.0 < 1.0/p->permtot)
				p->inprmt=100.0;
		  	else 
				p->inprmt = 1.0/p->permtot;
			
			/* Adjust for Ca activation level and temp. */
			p->kn_pT = p->kn_p*p->permtot*std::pow(p->Qkn_p,(p->T-273.0-37.0)/10.);
			p->kp_nT = p->kp_n*p->inprmt*std::pow(p->Qkp_n,(p->T-273.0-37.0)/10.);
			
			/* Compute fapp, the attachment rate from weak to strong, pre-rotated state */
		    p->fappT = p->fapp*std::pow(p->Qfapp,(p->T-273.0-37.0)/10.);
			
			/* Compute gapp, the detachment rate from strong, pre-rotated to the weakly-bound state */
			/* Compute SL modifier for gapp to increase rate st shorter SL */
		  	p->gapslmd = 1.0 + (1.0 - p->SOVFThick)*p->gslmod;
		    p->gappT = p->gapp*p->gapslmd*std::pow(p->Qgapp,(p->T-273.0-37.0)/10.);
			
			/* Set rate for rates between pre-force and force states*/
			/* Compute modifiers based on mean strain of states. */
			// oomph_info << p->sign(NV_Ith_S(x,8)) << std::endl;
			// oomph_info << p->hfmdc << std::endl;
			// oomph_info << NV_Ith_S(x,8) << std::endl;
			// oomph_info << p->x_0 << std::endl;
			// oomph_info << NV_Ith_S(x,8) << std::endl;
			// oomph_info << p->x_0 << std::endl;

		  	p->hfmd=std::exp(-p->sign(NV_Ith_S(x,8))*p->hfmdc*((NV_Ith_S(x,8)/p->x_0)*(NV_Ith_S(x,8)/p->x_0)));
		  	p->hbmd=std::exp(p->sign((NV_Ith_S(x,7)-p->x_0))*p->hbmdc*(((NV_Ith_S(x,7)-p->x_0)/p->x_0)*((NV_Ith_S(x,7)-p->x_0)/p->x_0)));
			
			/* Combine modifiers of hf and hb */
			p->hfT=p->hf*p->hfmd*std::pow(p->Qhf,(p->T-273.0-37.0)/10.);
			p->hbT=p->hb*p->hbmd*std::pow(p->Qhb,(p->T-273.0-37.0)/10.); 
			
			/* Set rate for rates gxb, the ATP using XB transition */
			/* Add term for distortion dependence of gxb */
	    p->gxbmd = p->heav(p->x_0-NV_Ith_S(x,7))*std::exp(p->sigmap*((p->x_0-NV_Ith_S(x,7))/p->x_0)*((p->x_0-NV_Ith_S(x,7))/p->x_0))+(1.-p->heav(p->x_0-NV_Ith_S(x,7))*std::exp(p->sigman*(((NV_Ith_S(x,7)-p->x_0)/p->x_0)*(NV_Ith_S(x,7)-p->x_0)/p->x_0)));
	    p->gxbT = p->gxb*p->gxbmd*std::pow(p->Qgxb,(p->T-273.0-37.0)/10.);
			
			/* update all RUs */
			NV_Ith_S(dxdt,0) = -p->kn_pT * NV_Ith_S(x,0) + p->kp_nT * NV_Ith_S(x,2);
			NV_Ith_S(dxdt,2) = -p->kp_nT * NV_Ith_S(x,2) + p->kn_pT * NV_Ith_S(x,0);
			
			NV_Ith_S(dxdt,1) = -p->kn_pT * NV_Ith_S(x,1) + p->kp_nT * NV_Ith_S(x,3); 
			NV_Ith_S(dxdt,3) = -p->kp_nT * NV_Ith_S(x,3) + p->kn_pT * NV_Ith_S(x,1) - p->fappT * NV_Ith_S(x,3) + p->gappT * NV_Ith_S(x,4) + p->gxbT * NV_Ith_S(x,5); 
			NV_Ith_S(dxdt,4) = p->fappT * NV_Ith_S(x,3) - p->gappT * NV_Ith_S(x,4) - p->hfT * NV_Ith_S(x,4) + p->hbT * NV_Ith_S(x,5);
			NV_Ith_S(dxdt,5) = p->hfT * NV_Ith_S(x,4) - p->hbT * NV_Ith_S(x,5) - p->gxbT * NV_Ith_S(x,5);
			
			/* compute steady-state fractions in XBprer and XBpostr using King-Altman rule */
			p->SSXBprer = (p->hb*p->fapp+p->gxb*p->fapp)/(p->gxb*p->hf+p->fapp*p->hf+p->gxb*p->gapp+p->hb*p->fapp+p->hb*p->gapp+p->gxb*p->fapp);
			p->SSXBpostr = p->fapp*p->hf/(p->gxb*p->hf+p->fapp*p->hf+p->gxb*p->gapp+p->hb*p->fapp+p->hb*p->gapp+p->gxb*p->fapp);
			/* compute normalization for scaling active and passive force */
			p->Fnordv = p->kxb*p->x_0*p->SSXBpostr;
			
			/* compute force here */
			p->force = p->kxb*p->SOVFThick*(NV_Ith_S(x,7)*NV_Ith_S(x,5)+NV_Ith_S(x,8)*NV_Ith_S(x,4));
			
			/* Compute SL, Set to 0 for isometric, to 1 for active contraction */
			/* Note that passive force is specified in terms of maximal active force */
			/* read passive force from table file specified in terms of maximal active force */
			/* compyet derivative of SL */
			
			// Note that passive force is specified in terms of maximal active force.
			// Passive force is computed in two ways depending on if cell is skinned or intact.
			p->ppforce=p->sign((NV_Ith_S(x,6)-p->SLrest))*p->PCon_t*(std::exp(p->PExp_t*std::fabs((NV_Ith_S(x,6)-p->SLrest)))-1);
		  //   if (singlecell == false) 
				// ppforce += std::max(0.0,(PCon_col*std::exp(PExp_col*(NV_Ith_S(x,6)-SLcol))));
			
			// Compute afterload either as a contant or from the series elastic element
		 	p->afterload = 0.0;
		 	p->PreloadF = 0.0;
			
		 //    if (SEon == true) {
			// 	afterload = KSE*(SLset-NV_Ith_S(x,6));
			// 	if ((SEon_LengthClamp==true)&&(index==(pulse_number-1))) {
			// 		afterload=5000.0*(SLset-NV_Ith_S(x,6));
			// 	}
		 //    } 
			// else {
				p->afterload=0.0; 
		    // }
			
			//  Compute the integral of forces 
			NV_Ith_S(dxdt,11)=(-p->ppforce+p->PreloadF+(-p->force/p->Fnordv)+p->afterload);
			
			// Compute the derivative of SL for an isolated muscle
			// oomph_info << "NV_Ith_S(x,11) " <<  NV_Ith_S(x,11) << std::endl;
			// oomph_info << "p->SLset " <<  p->SLset << std::endl;
			// oomph_info << "NV_Ith_S(x,6) " <<  NV_Ith_S(x,6) << std::endl;
			// oomph_info << "p->visc " <<  p->visc << std::endl;
			// oomph_info << "p->massf " <<  p->massf << std::endl;
			// oomph_info << "p->SLmin " <<  p->SLmin << std::endl;
			// oomph_info << "p->SLmax " <<  p->SLmax << std::endl;

		  p->dSLisolated = ((NV_Ith_S(x,11) + (p->SLset-NV_Ith_S(x,6))*p->visc)/p->massf )*p->heav(NV_Ith_S(x,6)-p->SLmin)*p->heav(p->SLmax-NV_Ith_S(x,6));
		 	// if (contraction==false) 
				// dSLisolated = 0.;
			
			
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
		 // 	LVpress=SigPrs*A_0/(1.36*Betavl*std::pow(NV_Ith_S(x,13),0.667));
			
			// // Compute Flows - assumes aortic valve so flow in unidirectional
			// // Units are l/s
			// iflow=heav(LVpress-NV_Ith_S(x,12))*((LVpress-NV_Ith_S(x,12))/Rchar);
			// // Allow active ejection on last beat only 
			// if (index<(pulse_number-1))
			// 	iflow=0.0;
			// if (Isovolume==true) 
			// 	iflow=0.0;
			
		 //    // Compute the derivative of arterial pressure
			// NV_Ith_S(dxdt,12)=(iflow/Cwind - (NV_Ith_S(x,12)-ARTpset)/(Cwind*Ra));
		 //    if (WholeHeart==false) 
			// 	NV_Ith_S(dxdt,12)=0.0;
		    
			// // Compute the  derivative of the change in Left Ventricular Volume
			// // Do not allow volume less than 0
		 //    NV_Ith_S(dxdt,13)=-(iflow)*heav(NV_Ith_S(x,13));
		 //    if (WholeHeart==false) 
			// 	NV_Ith_S(dxdt,13)=0.0;
			
			// // Compute the  derivative of the change in SL
		 //    dSLwholeheart=SLref*(1/std::pow(Vref+0.333*Vwall,0.333))*(NV_Ith_S(dxdt,13)/3)/std::pow(NV_Ith_S(x,13)+0.333*Vwall,0.667);
			
			// // Compute the derivative of SL based using cell/tissue or simulated whole heart model
			// if (WholeHeart==true) {
			// 	NV_Ith_S(dxdt,6) = dSLwholeheart;
			// } else {
				NV_Ith_S(dxdt,6) = p->dSLisolated;
			// }
			
			// //========>>>>>>>>
			
			
			
			
			p->dSL=NV_Ith_S(dxdt,6);
			
			// Compute the duty fractions using King-Alman Rule
			// Compute for states XBpref and XBf
			p->dtyf_prer=(p->hbT*p->fappT+p->gxbT*p->fappT)/(p->fappT*p->hfT+p->gxbT*p->hfT+p->gxbT*p->gappT+p->hbT*p->fappT+p->hbT*p->gappT+p->gxbT*p->fappT);
			p->dtyf_postr=p->fappT*p->hfT/(p->fappT*p->hfT+p->gxbT*p->hfT+p->gxbT*p->gappT+p->hbT*p->fappT+p->hbT*p->gappT+p->gxbT*p->fappT);
			
			// Compute mean strain by consider two competing effects - SL motion and XB cycling 
			if(p->myo == true) {
				NV_Ith_S(dxdt,8) = p->dSL/2. + p->xPsi*(1/p->dtyf_prer)*(-NV_Ith_S(x,8)*p->fappT+(NV_Ith_S(x,7)-p->x_0-NV_Ith_S(x,8))*p->hbT);
				NV_Ith_S(dxdt,7) = p->dSL/2. + p->xPsi*(1/p->dtyf_postr)*(p->x_0+NV_Ith_S(x,8)-NV_Ith_S(x,7))*p->hfT;
			} 
			else {
				NV_Ith_S(dxdt,8)=0.0;
				NV_Ith_S(dxdt,7)=0.0;
			}
			
			// Compute derivative of z-line end of single overlap region
			p->sovr_ze_dt=0.;
			if(p->len_thick/2.>NV_Ith_S(x,6)/2.) 
				p->sovr_ze_dt=-p->dSL/2.;
			
			// Compute derivative of centerline of end of single overlap region
		    p->sovr_cle_dt=0.;
		    if (NV_Ith_S(x,6)/2.-(NV_Ith_S(x,6)-p->len_thin)>p->len_hbare/2.)
				p->sovr_cle_dt=-p->dSL/2;
			
			// Compute the derivative of the length of the single overlap 
			p->len_sovr_dt=p->sovr_ze_dt-p->sovr_cle_dt;
		 	
			// Compute the derivative of the overlap fraction for thin filament
			p->SOVFThin_dt = p->len_sovr_dt/p->len_thin;
			
			// Compute the derivative of the overlap fraction for thick filament
			p->SOVFThick_dt = p->len_sovr_dt*2./(p->len_thick-p->len_hbare);
			
			/*Compute fraction of strongly bound for compute effective binding of Ca */
			p->FrSBXB=((NV_Ith_S(x,5)+NV_Ith_S(x,4))/(p->SSXBpostr + p->SSXBprer));
			p->FrSBXB_dt=((NV_Ith_S(dxdt,5)+NV_Ith_S(dxdt,4))/(p->SSXBpostr + p->SSXBprer));
			
			p->TropToT=p->Trop_conc*((1.0-p->SOVFThin)*NV_Ith_S(x,9) + p->SOVFThin*(p->FrSBXB*NV_Ith_S(x,10)+(1.0-p->FrSBXB)*NV_Ith_S(x,9)));
			p->TropToT_dt=p->Trop_conc*(-1.0*p->SOVFThin_dt*NV_Ith_S(x,9)+(1.0-p->SOVFThin)*NV_Ith_S(dxdt,9) + p->SOVFThin_dt*(p->FrSBXB*NV_Ith_S(x,10)+(1.0-p->FrSBXB)*NV_Ith_S(x,9)) + p->SOVFThin*(p->FrSBXB_dt*NV_Ith_S(x,10)+p->FrSBXB*NV_Ith_S(dxdt,10)-p->FrSBXB_dt*NV_Ith_S(x,9)+(1.0-p->FrSBXB)*NV_Ith_S(dxdt,9)));

			p->active_force = 10*p->force/p->Fnordv;   /* in units of mN/(mm^2) */

			return 0;
	}



} //End namespace

#endif