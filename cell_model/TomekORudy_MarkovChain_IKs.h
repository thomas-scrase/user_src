// TorORd dynamic chloride model
// https://github.com/jtmff/torord/blob/master/matlab/model_Torord_dynCl.m

// Swapped out IKs formulation for the MC model I have optimised based on that of Rudy

#ifndef OOMPH_IMPLICIT_TOMEKORUDY_VENT_MARKOVCHAIN_IKS_HEADER
#define OOMPH_IMPLICIT_TOMEKORUDY_VENT_MARKOVCHAIN_IKS_HEADER
	
// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "cell_model_base.h"

namespace oomph{

	class TomekORudyVentMarkovChainIKs : public CellModelBase
	{
	public:
		TomekORudyVentMarkovChainIKs();

		void explicit_timestep(CellState &Cellstate, Vector<double> &new_state);

		bool compatible_cell_types(const unsigned& cell_type);

		inline void return_initial_membrane_potential(double &v, const unsigned &cell_type=0);
		inline bool return_initial_state_variable(const unsigned &n, double &v, const unsigned &cell_type);

		//this can be overloaded in case the membrane capacitance is dependent on cell type or state parameters
		inline double cm(CellState &state) {return 1.0;}

		inline void custom_output(CellState &state, Vector<double> &output) override {}

		inline bool model_calculates_jacobian_entries() {return true;}

		inline unsigned Num_Variables(){return 59;}

		//the model does not require derivatives to be provided, it's
		//	explicit and calculates them itself
		inline unsigned required_derivatives(){return 0;}
		//The model has three black box parameters
		inline unsigned required_black_box_parameters(){return 1;}
		inline unsigned required_external_data() {return 0;}


		//a virtual function to extract black box nodal parameters from the cell state
		//	we implement this as virtual so that in the case this class is used in a
		//	combined cell model class. Then it can be overloaded so that black box nodal
		//	parameters associated with each model have a unique index
		virtual inline void extract_black_box_parameters_TomekORudyVent(CellState &Cellstate);

		//Set the temperature of the model
		virtual void set_temp(const double& new_t){
			T = new_t;
		}

	protected:

		double nao;
		double ko;
		double cao;

		double R;
		double T;
		double F;

		double FNORT;

		double L;
		double rad;
		double vcell;
		double Ageo;
		double Acap;
		double vmyo;
		double vnsr;
		double vjsr;
		double vss;

	};

}


#endif