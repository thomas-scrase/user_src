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

	//forward declare OptimisationEquations so that MetaProblemBase can be defined
	class OptimisationEquations;

	//I exist solely to remove circular dependency, OptimisationEquations require these two functions but if
	//	they have a pointer to NelderMeadProblemClass then it's a pain to compile becase NelderMeadProblemClass
	//	contains references to OptimisationEquations
	//Implemented as borken to ensure these functions are overloaded
	class MetaProblemBase
	{
	public:

		MetaProblemBase(){}


		virtual void get_sub_problem_contribution_to_node(OptimisationEquations &node, Vector<double> &contribution){
			throw OomphLibError("get_sub_problem_contribution_to_node has been called but not defined, this should never happen",
						OOMPH_CURRENT_FUNCTION,
						OOMPH_EXCEPTION_LOCATION);
		}

		//get the cost incurred by the values of the variables in the optimisation element passed
		//used as VariableValuesContribution for the optimisation elements
		virtual void get_cost_of_variables(OptimisationEquations &optimisation_equations, Vector<double> &costs){
			throw OomphLibError("get_cost_of_variables has been called but not defined, this should never happen",
						OOMPH_CURRENT_FUNCTION,
						OOMPH_EXCEPTION_LOCATION);
		}
	};



	#define opt_heav(a) ((a) < (0.) ? (0.0) : (1.0))

	class OptimisationEquations	:	public virtual FiniteElement
	{
	public:

		OptimisationEquations(MetaProblemBase* owner_problem_pt) : 	OwnerProblemPt(owner_problem_pt),
																	Internal_Data_Pt(0),
																	N_Internal_Data(0)	{	}

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
			// std::cout << "create_internal_data" << std::endl;
			Internal_Data_Pt = this->add_internal_data(new Data(n_data), true);
			N_Internal_Data = n_data;
			// std::cout << this << N_Internal_Data << std::endl;
		}

		//return the ith internal data value
		double get_internal_data(const unsigned &i){
			return this->internal_data_pt(Internal_Data_Pt)->value(i);
		}
		//set the ith internal data value
		void set_internal_data(const unsigned &i, const double &value){
			this->internal_data_pt(Internal_Data_Pt)->set_value(i, value);
		}
		//pack all internal data into a vector
		void get_all_internal_data(Vector<double> &parameters){
			// std::cout << this << std::endl;
			// std::cout << N_Internal_Data << std::endl;
			// std::cout << "get_all_internal_data" << std::endl;
			// std::cout << N_Internal_Data << std::endl;
			parameters.resize(N_Internal_Data);
			// std::cout << "resized" << std::endl;
			for(unsigned i=0; i<N_Internal_Data; i++){
				// std::cout << "paramter " << i << std::endl;
				parameters[i] = get_internal_data(i);
			}
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


		//overload to make element non-abstract
		void output(std::ostream &outfile) override;
		//all of these throw errors, it makes no sense to call them in the context of these elements
		void output(std::ostream &outfile, const unsigned &nplot) override {
			throw OomphLibError("Broken output fct",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}
		void output(FILE* file_pt) override {
			throw OomphLibError("Broken output fct",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}
		void output(FILE* file_pt, const unsigned &n_plot) override {
			throw OomphLibError("Broken output fct",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}
		void output_fct(std::ostream &outfile, const unsigned &nplot, FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) override {
			throw OomphLibError("Broken output fct",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}
		void output_fct(std::ostream &outfile, const unsigned &nplot, const double& time, FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt) override {
			throw OomphLibError("Broken output fct",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}
		void shape(const Vector<double> &s, Shape &psi) const override {
			throw OomphLibError("Broken output fct",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}


		//NOTE: the following two functions are separate for the (admittedly rare) case in which the oomph-lib
		//	newton solver is used for the solving. In such a case the variable_value_contribution must be
		//	split out among the residual entries
		//ENDNOTE

		//function called to get external source contribution
		void external_source_contribution(Vector<double> &residual){
			if(OwnerProblemPt==NULL){
				throw OomphLibError("WARNING: OwnerProblemPt has not been set for this optimisation_equation\n",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
			}
			else{
				OwnerProblemPt->get_sub_problem_contribution_to_node(*this, residual);
			}
		}
		//function called to get contribution from values of variables
		void variable_values_contribution(Vector<double> &residual){
			if(OwnerProblemPt==NULL){
				throw OomphLibError("WARNING: OwnerProblemPt has not been set for this optimisation_equation\n",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
			}
			else{
				OwnerProblemPt->get_cost_of_variables(*this, residual);
			}
		}

	protected:
		//the fill in residual function for the element
		virtual void fill_in_generic_residual_contribution_optimisation(
  			Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  			DenseMatrix<double> &mass_matrix, unsigned flag);

	private:

		//Index of pointer to the internal datas
		unsigned Internal_Data_Pt;
		//The number of internal data
		unsigned N_Internal_Data;

		//The problem to which this optimisation element belongs
		MetaProblemBase* OwnerProblemPt;
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