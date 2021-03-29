#ifndef OOMPH_EXPLICIT_TNNP_06_VENT_HEADER
#define OOMPH_EXPLICIT_TNNP_06_VENT_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "cell_model_base.h"

namespace oomph{

	class ExplicitTNNP06Vent : public CellModelBase
	{
	public:
		ExplicitTNNP06Vent();

		void explicit_timestep(CellState &Cellstate, Vector<double> &new_state);

		bool compatible_cell_types(const unsigned& cell_type);

		inline void return_initial_membrane_potential(double &v, const unsigned &cell_type=0);
		inline bool return_initial_state_variable(const unsigned &n, double &v, const unsigned &cell_type);

		//this can be overloaded in case the membrane capacitance is dependent on cell type or state parameters
		inline double cm(CellState &state) {return Tent_CAPACITANCE;}

		inline void custom_output(CellState &state, Vector<double> &output) override {}

		inline bool model_calculates_jacobian_entries() {return true;}

		inline unsigned required_nodal_variables(const unsigned &cell_type=0){return 20;}
		//the model does not require derivatives to be provided, it's
		//	explicit and calculates them itself
		inline unsigned required_derivatives(){return 0;}
		//The model has three black box parameters
		inline unsigned required_black_box_parameters(){return 1;}
		inline unsigned required_external_data() {return 0;}

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
}

#endif