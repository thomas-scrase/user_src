#ifndef OOMPH_TNNP06_VENT_HEADER
#define OOMPH_TNNP06_VENT_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif
	
#include "../cell_models/cell_model_base.h"

namespace oomph{

class TNNP06 : public CellModelBaseFullyPartitioned
{
public:
	TNNP06();

	~TNNP06();

	double return_initial_state_variable(const unsigned &v, const unsigned &cell_type);

	double return_initial_membrane_potential(const unsigned &cell_type);


	void Calculate_Derivatives(const Boost_State_Type &Variables,
														const double &t,
														const unsigned &cell_type,
														const double &Istim,
														const Vector<double> &Other_Parameters,
														const Vector<double> &Other_Variables,
														Vector<double> &Variable_Derivatives,
														double &Iion);


	void get_optional_output(const Boost_State_Type &Variables,
													const double &t,
													const unsigned &cell_type,
													const double &Istim,
													const Vector<double> &Other_Parameters,
													const Vector<double> &Other_Variables,
													Vector<double> &Out);

	enum Cell_Variables_Enum : unsigned
	{
		Cai_TT,
		Nai_TT,
		Ki_TT,
		CaSRf_TT,
		CaSS_TT,
		sm_TT,
		sh_TT,
		sj_TT,
		sd_TT,
		sf_TT,
		sf2_TT,
		sfcass_TT,
		sr_TT,
		ss_TT,
		sxr1_TT,
		sxr2_TT,
		sxs_TT,
		y_TT,
		sRR_TT//,
		// sOO_TT
	};
	enum Other_Parameters_Enum : unsigned
	{
		ZIndex_TT
	};
		

protected:
		double Tent_R;
		double Tent_F;
		double Tent_T;
		double Tent_Ko;
		double Tent_Cao;
		double Tent_Nao;
		double Tent_Vc;
		double Tent_Vsr;
		double Tent_Vss;
		double Tent_Bufc;
		double Tent_Kbufc;
		double Tent_Bufsr;
		double Tent_Kbufsr;
		double Tent_Bufss;
		double Tent_Kbufss;
		double Tent_Vmaxup;
		double Tent_Kup;
		double Tent_Vrel;
		double Tent_k1_;
		double Tent_k2_;
		double Tent_k3;
		double Tent_k4;
		double Tent_EC;
		double Tent_maxsr;
		double Tent_minsr;
		double Tent_Vleak;
		double Tent_Vxfer;
		double Tent_CAPACITANCE;
		double Tent_pKNa;
		double Tent_GbNa;
		double Tent_KmK;
		double Tent_KmNa;
		double Tent_knak;
		double Tent_GCaL;
		double Tent_GbCa;
		double Tent_knaca;
		double Tent_KmNai;
		double Tent_KmCa;
		double Tent_ksat;
		double Tent_n;
		double Tent_GpCa;
		double Tent_KpCa;
		double Tent_GpK;
		double Tent_RTONF;
		double Tent_inverseVcF2;
		double Tent_inverseVcF;
		double Tent_inversevssF2;


};


} //End namespace

#endif