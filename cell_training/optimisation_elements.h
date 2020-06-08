//Optimisation elements.
//In jacobian type meta problems a single Q optimisation element is used
//	to solve for the parameters using the usual oomph lib system

//In nelder mead type meta problems the optimisation equations are used only
// to store the values of the parameters at each node in the simplex, this is
// so a single trainable cell type class is required to be inherited from for
// a cell type to be trainable


#ifndef OOMPH_OPTIMISATION_ELEMENTS_HEADER
#define OOMPH_OPTIMISATION_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
	#include <oomph-lib-config.h>
#endif

//OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Qelements.h"
#include "../generic/Telements.h"
#include "../generic/error_estimator.h"

namespace oomph{

	#define opt_heav(a) ((a) < (0.) ? (0.0) : (1.0))

	class OptimisationEquations	:	public virtual FiniteElement
	{
	public:

		typedef void (*ResidualContributionFromExternalSourceFctPt) (Vector<double> &residual);

		OptimisationEquations()	{	}

		/// Broken copy constructor
		OptimisationEquations(const OptimisationEquations& dummy){ 
			BrokenCopy::broken_copy("OptimisationEquations");
		} 

		/// Broken assignment operator
		void operator=(const OptimisationEquations&) {
			BrokenCopy::broken_assign("OptimisationEquations");
		}

		//We don't need any values at the nodes
		inline unsigned required_nvalue(const unsigned &n) const {return 0;}

		void create_internal_data(const unsigned &n_data){
			Internal_Data_Pt = this->add_internal_data(new Data(n_data), true);
			N_Internal_Data = n_data;

			//Resize the bools and mean and relative errors used to calculate contribution to the residual from dof value
			Add_Internal_Data_Range.resize(N_Internal_Data);
			Internal_Data_Mean.resize(N_Internal_Data);
			Internal_Data_Percentage_Range.resize(N_Internal_Data);

			for(unsigned i=0; i<Add_Internal_Data_Range.size(); i++){
				//Default bools to false
				Add_Internal_Data_Range[i] = false;
				//Default means to zero
				Internal_Data_Mean[i] = 0.0;
				//Default percentage range to zero
				Internal_Data_Percentage_Range[i] = 0.0;
			}
		}

		//return the ith internal data value
		double get_internal_data(const unsigned &i){
			return this->internal_data_pt(Internal_Data_Pt)->value(i);
		}

		void set_internal_data(const unsigned &i, double &value){
			this->internal_data_pt(Internal_Data_Pt)->set_value(i, value);
		}

		inline void turn_on_relative_error_contribution(const unsigned &i){
			Add_Internal_Data_Range[i] = true;
		}

		inline void turn_off_relative_error_contribution(const unsigned &i){
			Add_Internal_Data_Range[i] = false;
		}

		inline void set_variable_mean(const unsigned &i, const double &mean){
			Internal_Data_Mean[i] = mean;
		}

		inline void set_variable_permitted_relative_error(const unsigned &i, const double &permitted_error){
			Internal_Data_Percentage_Range[i] = permitted_error;
		}

		//====================================================================
		//Residual and Jacobian functions
		//====================================================================
		/// Add the element's contribution to its residual vector (wrapper)
		void fill_in_contribution_to_residuals(Vector<double> &residuals)
		{
		   	//Call the generic residuals function with flag set to 0 and using
		   	//a dummy matrix
			fill_in_generic_residual_contribution_optimisation(
				residuals,GeneralisedElement::Dummy_matrix,
				GeneralisedElement::Dummy_matrix,0);
		}
		 
		/// \short Add the element's contribution to its residual vector and 
		/// the element Jacobian matrix (wrapper)
		void fill_in_contribution_to_jacobian(Vector<double> &residuals,
		                                   DenseMatrix<double> &jacobian)
		{
			//Perform fill in procedure using finite differencing
			FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
		}

		/// Add the element's contribution to its residuals vector,
		/// jacobian matrix and mass matrix
		void fill_in_contribution_to_jacobian_and_mass_matrix(
											Vector<double> &residuals, DenseMatrix<double> &jacobian,
											DenseMatrix<double> &mass_matrix)
		{
			//Call the generic routine with the flag set to 2
			// fill_in_generic_residual_contribution_cell_interface(residuals,jacobian,mass_matrix,2);
			FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,jacobian,mass_matrix);
		}

		//Access function to external source
		ResidualContributionFromExternalSourceFctPt& get_external_contribution_to_residual()
			{return Get_External_Contribution_To_Residual;}


		//overload to make element non-abstract
		void output(std::ostream &outfile) override {}
		void output(std::ostream &outfile, const unsigned &nplot) override {}
		void output(FILE* file_pt) override {}
		void output(FILE* file_pt, const unsigned &n_plot) override {}
		void output_fct(std::ostream &outfile, const unsigned &nplot, FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) override {}
		void output_fct(std::ostream &outfile, const unsigned &nplot, const double& time, FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt) override {}
		void shape(const Vector<double> &s, Shape &psi) const override {}

		inline double total_residual_from_data_outside_of_range();

	protected:
		//the fill in residual function for the element
		virtual void fill_in_generic_residual_contribution_optimisation(
  			Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  			DenseMatrix<double> &mass_matrix, unsigned flag);


	private:

		ResidualContributionFromExternalSourceFctPt Get_External_Contribution_To_Residual;

		virtual inline void residual_from_comparrison_to_experiment(Vector<double> &residual){
			//call the function pointer to an external source
			Get_External_Contribution_To_Residual(residual);
		}

		virtual inline void residual_from_data_outside_of_range(const double &x, const double &x_mean, const double &relative_range, double &residual);

		//Pointer to the internal datas
		unsigned Internal_Data_Pt;
		//The number of internal data
		unsigned N_Internal_Data;


		//Flag, add range of Internal_Data_Pt values to the residual
		std::vector<bool> Add_Internal_Data_Range;
		//The Mean of the internal data
		Vector<double> Internal_Data_Mean;
		//The permitted range aroung the mean value
		Vector<double> Internal_Data_Percentage_Range;
	};



	class QOptimisationElement :
	public virtual QElement<1, 2>,
	public virtual OptimisationEquations
	{
	public:
		//provide unique final overrider
		void output(std::ostream &outfile){OptimisationEquations::output(outfile);}
		void output(std::ostream &outfile, const unsigned &nplot){OptimisationEquations::output(outfile, nplot);}
		void output(FILE* file_pt){OptimisationEquations::output(file_pt);}
		void output(FILE* file_pt, const unsigned &n_plot){OptimisationEquations::output(file_pt, n_plot);}
		void output_fct(std::ostream &outfile, const unsigned &nplot, FiniteElement::SteadyExactSolutionFctPt exact_soln_pt){OptimisationEquations::output_fct(outfile, nplot, exact_soln_pt);}
		void output_fct(std::ostream &outfile, const unsigned &nplot, const double& time, FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt){OptimisationEquations::output_fct(outfile, nplot, time, exact_soln_pt);}
		void shape(const Vector<double> &s, Shape &psi) const override {}

	};


}
#endif