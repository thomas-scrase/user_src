#include "cell_model_base.h"

namespace oomph{
	CellModelBase::CellModelBase() : Mutation_pt(0)
	{
		//Most cell models are too complicated to have an analytic jacobian matrix, exceptions include explicit timestepping models
		Model_Calculates_Jacobian_Entries = false;

		//set the default values for requests - minimises data transfer for default values while
		//	keeping most commonly used data
		Required_Derivatives = 1; //number of required derivatives
		Requires_Vm = true; //does the model require membrane potential
		Requires_Strain = false; //does the model require strain
		Requires_Na_o = true; //does the model require external sodium concentration
		Requires_K_o = true; //does the model require external potassium concentration
		Requires_Ca_o = true; //does the model require external calcium concentration
		Requires_Cell_Type = true; //does the model require cell type
		Requires_Fibrosis = false; //does the model require fibrosis type?
		Requires_AB_Index = false; //does the model require ab index
		Requires_RV_Index = false; //does the model require rv index
		Requires_IS_Index = false; //does the model require is index

		Requires_dt = false; //does the model use an explicit time-stepping method and hence require the temporal increment
		Requires_previous_values = false; //does the model require the previous time values of variables? Generally used by explicit time stepping methods
	}
}