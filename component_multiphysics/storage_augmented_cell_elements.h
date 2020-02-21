//LIC// ====================================================================
//LIC// This file contains component_multiphysics elements, elements which are
//LIC//	themselves useless but occupy a domain in a multi-discretised multi-
//LIC//	physics element.
//LIC//
//LIC//	The storage_augmented_cell_element allows for communication with a
//LIC//	cell model and storage of fibre orientation and diffusion coefficients.
//LIC//	Since it doesn't make physical sense to refine the cell mesh it suffices
//LIC// that these elements are non-refineable.
//LIC//
//LIC// ====================================================================

#ifndef OOMPH_STORAGE_AUGMENTED_CELL
#define OOMPH_STORAGE_AUGMENTED_CELL

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

//Cell interface elements (includes cell models)
#include "../cell_interface/cell_interface_elements.h"

//storage enrichment elements for fibre and diffusion ceofficient data
#include "../toms_utilities/vector_with_diffusion_storage_enrichment_elements.h"

//Monodomain elements
#include "../monodomain/monodomain_elements.h"

//Solid elements for the external solid element for geometric data
#include "../anisotropic_solid/anisotropic_solid_elements.h"


namespace oomph{
	//Storage convention: mono, fibre orientation, diffusion coefficients
	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class QStorageAugmentedCellElement	:	public virtual QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>,
											public virtual QVectorWithDiffusionStorageEnrichmentElement<DIM, NNODE_1D>
	{
	public:
		QStorageAugmentedCellElement()	:	QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>(),
											QVectorWithDiffusionStorageEnrichmentElement<DIM, NNODE_1D>()
		{
			//Pin the storage_enrichment_element dofs
			// unsigned n_node = this->nnode();
			// for(unsigned l=0; l < n_node; l++){
			// 	for(unsigned var = min_index_storage_enrichment(); var < max_index_storage_enrichment(); var++){
			// 		this->node_pt(l)->pin(var);
			// 	}
			// }
		}

		//Min and max variable indexes for output function and for ease of multiphysics elements
		virtual inline unsigned min_index_CellInterfaceEquations() const {return 0;}
		virtual inline unsigned max_index_CellInterfaceEquations() const {return NUM_VARS;}
		
		//identify the indexes of the diffusion matrix data
		virtual inline unsigned min_index_storage_enrichment() const {return NUM_VARS;}
		virtual inline unsigned max_index_storage_enrichment() const {return NUM_VARS + DIM*(DIM+1);}


		unsigned required_nvalue(const unsigned &n) const
		{return (QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>::required_nvalue(n) +
				QVectorWithDiffusionStorageEnrichmentElement<DIM, NNODE_1D>::required_nvalue(n));}

		//The output functions
		void output(std::ostream &outfile) {FiniteElement::output(outfile);}
		void output(std::ostream &outfile, const unsigned &nplot);	//Defined in .cc
		void output(FILE* file_pt){FiniteElement::output(file_pt);}
		void output(FILE* file_pt, const unsigned &n_plot){FiniteElement::output(file_pt,n_plot);}
		void output_fct(std::ostream &outfile, const unsigned &Nplot,FiniteElement::SteadyExactSolutionFctPt exact_soln_pt){FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}
		void output_fct(std::ostream &outfile, const unsigned &Nplot,const double& time,FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt){FiniteElement::output_fct(outfile,Nplot,time,exact_soln_pt);}
		// void compute_norm(double& el_norm){QUnsteadyHeatElement<DIM,3>::compute_norm(el_norm);}
		void compute_error(std::ostream &outfile,FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,const double& time,double& error, double& norm){FiniteElement::compute_error(outfile,exact_soln_pt,time,error,norm);}
		void compute_error(std::ostream &outfile,FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,double& error, double& norm){FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}

		//Residual and Jacobian functions
		void fill_in_contribution_to_residuals(Vector<double> &residuals)
		{
			CellInterfaceEquations<DIM>::fill_in_contribution_to_residuals(residuals);
			VectorWithDiffusionStorageEnrichmentEquations<DIM*(DIM+1)>::fill_in_contribution_to_residuals(residuals);
		}

		//No need to finite difference since there is no coupling between these elements
		void fill_in_contribution_to_jacobian(Vector<double> &residuals,DenseMatrix<double> &jacobian)
		{
			CellInterfaceEquations<DIM>::fill_in_contribution_to_jacobian(residuals,jacobian);
			VectorWithDiffusionStorageEnrichmentEquations<DIM*(DIM+1)>::fill_in_contribution_to_jacobian(residuals,jacobian);
		}

		void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals, DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
  		{
			CellInterfaceEquations<DIM>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,jacobian,mass_matrix);
			VectorWithDiffusionStorageEnrichmentEquations<DIM*(DIM+1)>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,jacobian,mass_matrix);
		}
	};



	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class TStorageAugmentedCellElement	:	public virtual TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>,
											public virtual TVectorWithDiffusionStorageEnrichmentElement<DIM, NNODE_1D>
	{
	public:
		TStorageAugmentedCellElement()	:	TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>(),
											TVectorWithDiffusionStorageEnrichmentElement<DIM, NNODE_1D>()
		{
			//Pin the storage_enrichment_element dofs
			// unsigned n_node = this->nnode();
			// for(unsigned l=0; l < n_node; l++){
			// 	for(unsigned var = min_index_storage_enrichment(); var < max_index_storage_enrichment(); var++){
			// 		this->node_pt(l)->pin(var);
			// 	}
			// }
		}

		//Min and max variable indexes for output function and for ease of multiphysics elements
		virtual inline unsigned min_index_CellInterfaceEquations() const {return 0;}
		virtual inline unsigned max_index_CellInterfaceEquations() const {return NUM_VARS;}
		
		//identify the indexes of the diffusion matrix data
		virtual inline unsigned min_index_storage_enrichment() const {return NUM_VARS;}
		virtual inline unsigned max_index_storage_enrichment() const {return NUM_VARS + DIM*(DIM+1);}

		unsigned required_nvalue(const unsigned &n) const
		{return (TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>::required_nvalue(n) +
				TVectorWithDiffusionStorageEnrichmentElement<DIM, NNODE_1D>::required_nvalue(n));}

		//The output functions
		void output(std::ostream &outfile) {FiniteElement::output(outfile);}
		void output(std::ostream &outfile, const unsigned &nplot); //Defined in .cc
		void output(FILE* file_pt){FiniteElement::output(file_pt);}
		void output(FILE* file_pt, const unsigned &n_plot){FiniteElement::output(file_pt,n_plot);}
		void output_fct(std::ostream &outfile, const unsigned &Nplot,FiniteElement::SteadyExactSolutionFctPt exact_soln_pt){FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}
		void output_fct(std::ostream &outfile, const unsigned &Nplot,const double& time,FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt){FiniteElement::output_fct(outfile,Nplot,time,exact_soln_pt);}
		// void compute_norm(double& el_norm){QUnsteadyHeatElement<DIM,3>::compute_norm(el_norm);}
		void compute_error(std::ostream &outfile,FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,const double& time,double& error, double& norm){FiniteElement::compute_error(outfile,exact_soln_pt,time,error,norm);}
		void compute_error(std::ostream &outfile,FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,double& error, double& norm){FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}

		//Residual and Jacobian functions: No need to finite difference since there is no coupling between these elements
		//	and storage enrichment data is pinned
		void fill_in_contribution_to_residuals(Vector<double> &residuals)
		{
			CellInterfaceEquations<DIM>::fill_in_contribution_to_residuals(residuals);
			VectorWithDiffusionStorageEnrichmentEquations<DIM*(DIM+1)>::fill_in_contribution_to_residuals(residuals);
		}

		void fill_in_contribution_to_jacobian(Vector<double> &residuals,DenseMatrix<double> &jacobian)
		{
			CellInterfaceEquations<DIM>::fill_in_contribution_to_jacobian(residuals,jacobian);
			VectorWithDiffusionStorageEnrichmentEquations<DIM*(DIM+1)>::fill_in_contribution_to_jacobian(residuals,jacobian);
		}

		void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals, DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
  		{
			CellInterfaceEquations<DIM>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,jacobian,mass_matrix);
			VectorWithDiffusionStorageEnrichmentEquations<DIM*(DIM+1)>::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,jacobian,mass_matrix);
		}
	};
}

#endif