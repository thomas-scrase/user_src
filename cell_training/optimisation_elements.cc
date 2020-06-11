#include "optimisation_elements.h"

namespace oomph{

	void OptimisationEquations::fill_in_generic_residual_contribution_optimisation(
	  			Vector<double> &residuals, DenseMatrix<double> &jacobian, 
	  			DenseMatrix<double> &mass_matrix, unsigned flag)
	{
		if(flag){
			throw OomphLibError("OptimisationEquations has no insight into the source of the residual vector, flag should never be set to 1.",
		                       OOMPH_CURRENT_FUNCTION,
		                       OOMPH_EXCEPTION_LOCATION);
		}
		//Declare the vector of resiodual contributions derrived from difference between experiment
		//	and simulation. This does not need to be the same size as the number of dof associated
		//	with this element.
		//This is where the sub problems are run
		Vector<double> contribution_from_external_source;
		external_source_contribution(contribution_from_external_source);

		// std::cout << "contribution_from_external_source.size() " << contribution_from_external_source.size() << std::endl;
		// std::cout << "residuals.size() " << residuals.size() << std::endl;
		//Loop over the external residual contributions and distribute them over the residual entires
		//	this is required for the case of the oomph-lib newton solver is used and the residual entries
		//	are required to have non-identical derivatives wrt the dofs
		for(unsigned i=0; i<contribution_from_external_source.size(); i++){
			// std::cout << "i " << i << std::endl;
			//local equations haven't been set up so maybe just get change all fill in generic residual things
			//	to other function names and remove all possibility of the elements being used as normal oomph
			//	elements
			// int local_eqn = internal_local_eqn(Internal_Data_Pt,i%N_Internal_Data);
			residuals[i%N_Internal_Data] += contribution_from_external_source[i];
		}

		//Get the residual contribution from variable values
		Vector<double> contribution_from_variables_values;
		variable_values_contribution(contribution_from_variables_values);

		//add them to their respective residual components
		//	contribution_from_variables_values should be of size N_Internal_Data
		//	but to be sure we loop over its size and take %N_Internal_Data
		for(unsigned i=0; i<contribution_from_variables_values.size(); i++){
			residuals[i%N_Internal_Data] += contribution_from_variables_values[i];
		}

	}


	//output function simple returns all internal data in braces
	void OptimisationEquations::output(std::ostream &outfile){
		outfile << "{";
		for(unsigned i=0; i<N_Internal_Data; i++){
			outfile << "\t" << get_internal_data(i);
		}
		outfile << "\t}\t";
	}


}