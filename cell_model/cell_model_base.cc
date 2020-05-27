#include "cell_model_base.h"

namespace oomph{
	CellModelBase::CellModelBase() : Mutation_pt(0)
	{
		//set the default values for requests
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
	}
}