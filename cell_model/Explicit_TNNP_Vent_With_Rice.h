#ifndef OOMPH_EXPLICIT_TNNP_VENT_WITH_RICE_HEADER
#define OOMPH_EXPLICIT_TNNP_VENT_WITH_RICE_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "cell_model_base.h"

namespace oomph{

class ExplicitTNNPVentWithRice : public CellModelBase
{
public:
	ExplicitTNNPVentWithRice();

	void explicit_timestep(CellState &Cellstate, Vector<double> &new_state);

	bool compatible_cell_types(const unsigned& cell_type);

	inline void return_initial_membrane_potential(double &v, const unsigned &cell_type=0);
	inline bool return_initial_state_variable(const unsigned &n, double &v, const unsigned &cell_type);

	//this can be overloaded in case the membrane capacitance is dependent on cell type or state parameters
	inline double cm(CellState &state) {return TTCell_CAPACITANCE;}

	inline void custom_output(CellState &state, Vector<double> &output) override {}

	inline bool model_calculates_jacobian_entries() {return true;}

	inline unsigned required_nodal_variables(const unsigned &cell_type=0){return 34;}
	//the model does not require derivatives to be provided, it's
	//	explicit and calculates them itself
	inline unsigned required_derivatives(){return 0;}
	//The model has three black box parameters
	inline unsigned required_black_box_parameters(){return 3;}
	inline unsigned required_external_data() {return 0;}


	//a virtual function to extract black box nodal parameters from the cell state
	//	we implement this as virtual so that in the case this class is used in a
	//	combined cell model class. Then it can be overloaded so that black box nodal
	//	parameters associated with each model have a unique index
	virtual inline void extract_black_box_parameters_ExplicitTNNPVentWithRice(double &abindex,
																double &isindex,
																double &rvindex,
																CellState &Cellstate);
protected:

	int sign(double a){ return ((a) < (0.) ? (-1.0) : (1.0));}
	double heav(double a){ return ((a) < (0.) ? (0.0) : (1.0));}
	// constants for Tent ventricle model
	//	from TNNP_MarkovIKr_function.h in haibo ni folders
	double TTCell_Ko;
	double TTCell_Cao;
	double TTCell_Nao;
	double TTCell_Vc;
	double TTCell_Vsr;
	double TTCell_Vss;
	double TTCell_Bufc;
	double TTCell_Kbufc;
	double TTCell_Bufsr;
	double TTCell_Kbufsr;
	double TTCell_Bufss;
	double TTCell_Kbufss;
	double TTCell_Vmaxup;
	double TTCell_Kup;
	double TTCell_Vrel;
	double TTCell_k1_;
	double TTCell_k2_;
	double TTCell_k3;
	double TTCell_k4;
	double TTCell_EC;
	double TTCell_maxsr;
	double TTCell_minsr;
	double TTCell_Vleak;
	double TTCell_Vxfer;
	double TTCell_R;
	double TTCell_F;
	double TTCell_T;
	double TTCell_CAPACITANCE;
	double TTCell_pKNa;
	double TTCell_GbNa;
	double TTCell_KmK;
	double TTCell_KmNa;
	double TTCell_knak;
	double TTCell_GCaL;
	double TTCell_GbCa;
	double TTCell_knaca;
	double TTCell_KmNai;
	double TTCell_KmCa;
	double TTCell_ksat;
	double TTCell_n;
	double TTCell_GpCa;
	double TTCell_KpCa;
	double TTCell_GpK;

	double TTCell_RTONF;
	double TTCell_inverseVcF2;
	double TTCell_inverseVcF;
	double TTCell_inversevssF2;





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
	double SLmax  ;
	double SLmin  ;
	double len_thin ;
	double len_thick;
	double len_hbare;
	//   Temperature Dependence
	double Qkon;
	double Qkoff;
	// double Qkoff = 1.4;
	double Qkn_p;
	double Qkp_n;
	double Qfapp;
	// double Qgapp = 6.25;
	double Qgapp;
	double Qhf;
	double Qhb;
	double Qgxb;
	//   Ca binding to troponin
	double kon    ;
	double koffL  ;
	double koffH  ;
	double perm50 ;
	double nperm  ;
	double kn_p   ;
	double kp_n   ;
	double koffmod;
	//   Thin filament regulation and crossbridge cycling
	double fapp   ;
	double gapp   ;
	double gslmod ;
	double hf     ;
	double hfmdc  ;
	double hb     ;
	double hbmdc  ;
	double gxb    ;
	double sigmap ;
	double sigman ;
	double xbmodsp;
	//   Mean strain of strongly-bound states
	double x_0    ;
	double xPsi   ;
	//   Normalized active a nd passive force
	double PCon_t ;
	double PExp_t ;
	double SL_c   ;
	double PCon_c ;
	double PExp_c ;
	//   Calculation of complete muscle response
	double massf  ;
	double visc   ;
	double KSE    ;
	double kxb    ;
	double Trop_conc;
	double Temp;

	double SLset ;         // initial length
    double SLrest;       //   (um) rest SL length for 0 passive force
    //End Rice Myofilaments constants
};


}

#endif