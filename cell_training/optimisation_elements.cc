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
		Vector<double> contribution_from_external_source;
		residual_from_comparrison_to_experiment(contribution_from_external_source);
	
		//Loop over the external residual contributions
		for(unsigned n=0; n<contribution_from_external_source.size(); n++){
			int local_eqn = internal_local_eqn(Internal_Data_Pt,n%N_Internal_Data);

			residuals[n] += contribution_from_external_source[n];
		}


		for(unsigned i=0; i<N_Internal_Data; i++){
			double x = this->internal_data_pt(Internal_Data_Pt)->value(i);
			//if the dof has a acceptable range
			if(Add_Internal_Data_Range[i]){
				double res_cont;
				residual_from_data_outside_of_range(x, Internal_Data_Mean[i], Internal_Data_Percentage_Range[i], res_cont);
				residuals[i] += res_cont;
			}
		}
	}

	//Let's define the default residual contribution from outside of range as linear ramps bounding a zero region
	inline void OptimisationEquations::residual_from_data_outside_of_range(const double &x, const double &x_mean, const double &relative_range, double &residual){
		residual = opt_heav(x-(x_mean+relative_range*std::abs(x_mean)))*(x-(x_mean+relative_range*std::abs(x_mean)))
				+	opt_heav(-x+(x_mean-relative_range*std::abs(x_mean)))*(-x+(x_mean-relative_range*std::abs(x_mean)));
	}

	//get the total contrubution to the residual from data being outside of permitted range
	inline double OptimisationEquations::total_residual_from_data_outside_of_range(){
		double residuals = 0.0;
		for(unsigned i=0; i<N_Internal_Data; i++){
			double x = this->internal_data_pt(Internal_Data_Pt)->value(i);
			//if the dof has a acceptable range
			if(Add_Internal_Data_Range[i]){
				double res_cont;
				residual_from_data_outside_of_range(x, Internal_Data_Mean[i], Internal_Data_Percentage_Range[i], res_cont);
				residuals += res_cont;
			}
		}
	}

}