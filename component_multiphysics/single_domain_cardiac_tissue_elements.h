//LIC// ====================================================================
//LIC// 
//LIC// ====================================================================


#ifndef OOMPH_SINGLE_DOMAIN_CARDIAC_TISSUE_SUB_ELEMENTS
#define OOMPH_SINGLE_DOMAIN_CARDIAC_TISSUE_SUB_ELEMENTS

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif


//Generic for the element with external element
// #include "generic.h"

//Include the storage_augmented_cell_elements header
#include "../toms_utilities/diff_augmented_cell_wrapper.h"

//Monodomain elements
#include "../cell_membrane_potential/cell_membrane_potential_elements.h"

//Solid elements for the external solid element for geometric data
#include "../anisotropic_solid/anisotropic_solid_elements.h"
#include "../anisotropic_solid/refineable_anisotropic_solid_elements.h"

//Cell interface elements (includes cell models)
#include "../cell_interface/cell_interface_elements.h"

#include "../solid/solid_traction_elements.h"



namespace oomph{

//We make this work for both monodomain and bidomain based conducting cell elements
template<unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
class QConductingCellElementWithIncompressibleAnisotropicSolidElement :
public virtual DiffAugmentedCell< QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >,
public virtual QAnisotropicPVDElementWithContinuousPressure<DIM>
{
public:
	QConductingCellElementWithIncompressibleAnisotropicSolidElement() :
	DiffAugmentedCell< QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >(),
	QAnisotropicPVDElementWithContinuousPressure<DIM>()
	{
	}


	unsigned required_nvalue(const unsigned &n) const
	{return QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::required_nvalue(n) +
			QAnisotropicPVDElementWithContinuousPressure<DIM>::required_nvalue(n);}



	//////////////OUTPUT FUNCTIONS/////////////////////////////////////////////
	unsigned int nscalar_paraview() const
	{
		return QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			nscalar_paraview();
	}


	void scalar_value_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot) const
	{
		QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			scalar_value_paraview(file_out, i, nplot);
	}

	void scalar_value_fct_paraview(std::ofstream& file_out,
										const unsigned& i,
										const unsigned& nplot,
										FiniteElement::SteadyExactSolutionFctPt
										exact_soln_pt) const
	{
		QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			scalar_value_fct_paraview(file_out, i, nplot, exact_soln_pt);
	}

	std::string scalar_name_paraview(const unsigned& i) const
	{
		return QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			scalar_name_paraview(i);
	}


	double get_interpolated_vm(const Vector<double>& s) const 
	{
		//Find number of nodes
		unsigned n_node = this->nnode();
		//Local shape function
		Shape psi(n_node);
		//Find values of shape function
		this->shape(s,psi);

		double interpolated_var = 0.0;
		//Loop over the local nodes and sum
		for(unsigned l=0;l<n_node;l++)
		{
			interpolated_var += QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::get_nodal_membrane_potential_BaseCellMembranePotential(l)*psi(l);
		}
		return interpolated_var;
	}

	///  Overload the standard output function with the broken default
 	void output(ostream &outfile) {FiniteElement::output(outfile);}

	/// \short Output function:  
	///  Output x, y, u, v, p, theta at Nplot^DIM plot points
	// Start of output function
	void output(ostream &outfile, const unsigned &nplot)
	{
		//vector of local coordinates
		Vector<double> s(DIM);
		Vector<double> xi(DIM);

		// Tecplot header info
		outfile << this->tecplot_zone_string(nplot);

		// Loop over plot points
		unsigned num_plot_points=this->nplot_points(nplot);
		for (unsigned iplot=0;iplot<num_plot_points;iplot++)
		{
			// Get local coordinates of plot point
			this->get_s_plot(iplot,nplot,s);

			// Get the Lagrangian coordinate
			this->interpolated_xi(s,xi);
			// Output the position of the plot point
			for(unsigned i=0;i<DIM;i++) 
			{
				outfile << this->interpolated_x(s,i) << " ";
			}
			// Output vm
			outfile << get_interpolated_vm(s) << std::endl;   
		}
		outfile << std::endl;

		// Write tecplot footer (e.g. FE connectivity lists)
		this->write_tecplot_zone_footer(outfile,nplot);
	} //End of output function

	/// \short C-style output function: Broken default
	void output(FILE* file_pt)
	{
		FiniteElement::output(file_pt);
	}
	
	///  \short C-style output function: Broken default
	void output(FILE* file_pt, const unsigned &n_plot)
	{
		FiniteElement::output(file_pt,n_plot);
	}
	
	/// \short Output function for an exact solution: Broken default
	void output_fct(ostream &outfile, const unsigned &Nplot,
	         FiniteElement::SteadyExactSolutionFctPt 
	         exact_soln_pt)
	{
		FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);
	}
	
	/// \short Output function for a time-dependent exact solution:
	/// Broken default.
	void output_fct(ostream &outfile, const unsigned &Nplot,
	         const double& time,
	         FiniteElement::UnsteadyExactSolutionFctPt 
	         exact_soln_pt)
	{
		FiniteElement::output_fct(outfile,Nplot,time,exact_soln_pt);
	}

	///////////////////////////////////////////////////////////////////////////




	//Inform cell element where to get diffusion tensor data from
	//Monodomain
	void get_diff_monodomain(const unsigned& ipt,
                            const Vector<double> &s,
                            const Vector<double>& x,
                            DenseMatrix<double>& D) const
	{
		QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			get_interpolated_diffusion_matrix(s, D);
	}
	//Bidomain
	void get_intracellular_conductivity_bidomain(const unsigned& ipt,
	                                            const Vector<double> &s,
	                                            const Vector<double>& x,
	                                            DenseMatrix<double>& G_i) const
	{
		QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			get_interpolated_diffusion_matrix(s, G_i);
	}
	void get_extracellular_conductivity_bidomain(const unsigned& ipt,
		                                            const Vector<double> &s,
		                                            const Vector<double>& x,
		                                            DenseMatrix<double>& G_e) const
	{
		QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			get_interpolated_diffusion_matrix(s, G_e);
	}


	//Inform solid element where to get diffusion tensor data from
	void anisotropic_matrix(const unsigned& ipt,
                           const Vector<double> &s,
                           const Vector<double>& xi,
                           DenseMatrix<double>& A)
    {
    	DiffAugmentedCell< QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
    		get_interpolated_preferential_vectors(s, A);
    }


	//Inform solid element where to get active strain from
	void driving_strain(const unsigned& ipt,
                       const Vector<double>& s,
                       const Vector<double>& xi,
                       Vector<double>& V)
    {
    	DiffAugmentedCell< QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
    		get_interpolated_active_strain(s);
    }




	//Tell the elements how to fill in residuals
	//Residuals fill in should be unchanged from just that of the underlying cell element
	void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
												DenseMatrix<double> &jacobian)
	{	
		//Call the residuals of the advection-diffusion eqautions
		DiffAugmentedCell< QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			fill_in_contribution_to_residuals(residuals);
		//Call the residuals of the Navier-Stokes equations
		QAnisotropicPVDElementWithContinuousPressure<DIM>::
			fill_in_contribution_to_residuals(residuals);
	}

	/// Add the element's contribution to its residuals vector,
	/// jacobian matrix and mass matrix
	void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
															DenseMatrix<double> &jacobian,
															DenseMatrix<double> &mass_matrix)
	{
		DiffAugmentedCell< QConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			fill_in_contribution_to_jacobian(residuals,jacobian);

		QAnisotropicPVDElementWithContinuousPressure<DIM>::
			fill_in_contribution_to_jacobian(residuals,jacobian);
	}

};



//============================================================================
/// FaceGeometry of a 3D QAnisotropicPVDElement element
//============================================================================
template<unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
class FaceGeometry<QConductingCellElementWithIncompressibleAnisotropicSolidElement<DIM, CELL_MODEL, CONDUCTANCE_MODEL> >:
public virtual SolidQElement<2,3>
{
public:
	FaceGeometry() : SolidQElement<2,3>() {}
};

//============================================================================
/// FaceGeometry of FaceGeometry of a 3D QAnisotropicPVDElement element
//============================================================================
template<unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
class FaceGeometry<FaceGeometry<QConductingCellElementWithIncompressibleAnisotropicSolidElement<DIM, CELL_MODEL, CONDUCTANCE_MODEL> > >:
public virtual SolidQElement<1,3>
{
public:
/// Constructor must call the constructor of the underlying solid element
	FaceGeometry() : SolidQElement<1,3>() {}
};














//We make this work for both monodomain and bidomain based conducting cell elements
template<unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
class TConductingCellElementWithIncompressibleAnisotropicSolidElement :
public virtual DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >,
public virtual TAnisotropicPVDElementWithContinuousPressure<DIM>
{
public:
	TConductingCellElementWithIncompressibleAnisotropicSolidElement() :
	DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >(),
	TAnisotropicPVDElementWithContinuousPressure<DIM>()
	{
	}


	unsigned required_nvalue(const unsigned &n) const
	{return TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::required_nvalue(n) +
			TAnisotropicPVDElementWithContinuousPressure<DIM>::required_nvalue(n);}


	// inline unsigned vm_index_BaseCellMembranePotential() const {return (TAnisotropicPVDElementWithContinuousPressure<DIM>::required_nvalue());}
	inline int solid_p_nodal_index() const {return TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::max_cell_variable_derivatives_index_plus_one_ConductingCellEquations();}

	//////////////OUTPUT FUNCTIONS/////////////////////////////////////////////

	unsigned nscalar_paraview() const
	{
		return TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			nscalar_paraview();
	}


	void scalar_value_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot) const
	{
		TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			scalar_value_paraview(file_out, i, nplot);
	}

	void scalar_value_fct_paraview(std::ofstream& file_out,
										const unsigned& i,
										const unsigned& nplot,
										FiniteElement::SteadyExactSolutionFctPt
										exact_soln_pt) const
	{
		TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			scalar_value_fct_paraview(file_out, i, nplot, exact_soln_pt);
	}

	std::string scalar_name_paraview(const unsigned& i) const
	{
		return TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::
			scalar_name_paraview(i);
	}


	double get_interpolated_vm(const Vector<double>& s) const 
	{
		//Find number of nodes
		unsigned n_node = this->nnode();
		//Local shape function
		Shape psi(n_node);
		//Find values of shape function
		this->shape(s,psi);

		double interpolated_var = 0.0;
		//Loop over the local nodes and sum
		for(unsigned l=0;l<n_node;l++)
		{
			interpolated_var += TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL >::get_nodal_membrane_potential_BaseCellMembranePotential(l)*psi(l);
		}
		return interpolated_var;
	}

	///  Overload the standard output function with the broken default
 	void output(ostream &outfile) {FiniteElement::output(outfile);}

	/// \short Output function:  
	///  Output x, y, u, v, p, theta at Nplot^DIM plot points
	// Start of output function
	void output(ostream &outfile, const unsigned &nplot)
	{
		//vector of local coordinates
		Vector<double> s(DIM);
		Vector<double> xi(DIM);

		// Tecplot header info
		outfile << this->tecplot_zone_string(nplot);

		// Loop over plot points
		unsigned num_plot_points=this->nplot_points(nplot);
		for (unsigned iplot=0;iplot<num_plot_points;iplot++)
		{
			// Get local coordinates of plot point
			this->get_s_plot(iplot,nplot,s);

			// Get the Lagrangian coordinate
			this->interpolated_xi(s,xi);
			// Output the position of the plot point
			for(unsigned i=0;i<DIM;i++) 
			{
				outfile << this->interpolated_x(s,i) << " ";
			}
			// Output vm
			outfile << get_interpolated_vm(s) << std::endl;   
		}
		outfile << std::endl;

		// Write tecplot footer (e.g. FE connectivity lists)
		this->write_tecplot_zone_footer(outfile,nplot);
	} //End of output function

	/// \short C-style output function: Broken default
	void output(FILE* file_pt)
	{
		FiniteElement::output(file_pt);
	}
	
	///  \short C-style output function: Broken default
	void output(FILE* file_pt, const unsigned &n_plot)
	{
		FiniteElement::output(file_pt,n_plot);
	}
	
	/// \short Output function for an exact solution: Broken default
	void output_fct(ostream &outfile, const unsigned &Nplot,
	         FiniteElement::SteadyExactSolutionFctPt 
	         exact_soln_pt)
	{
		FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);
	}
	
	/// \short Output function for a time-dependent exact solution:
	/// Broken default.
	void output_fct(ostream &outfile, const unsigned &Nplot,
	         const double& time,
	         FiniteElement::UnsteadyExactSolutionFctPt 
	         exact_soln_pt)
	{
		FiniteElement::output_fct(outfile,Nplot,time,exact_soln_pt);
	}

	///////////////////////////////////////////////////////////////////////////




	//Inform cell element where to get diffusion tensor data from
	//Monodomain
	void get_diff_monodomain(const unsigned& ipt,
                            const Vector<double> &s,
                            const Vector<double>& x,
                            DenseMatrix<double>& D) const
	{
		// oomph_info << "Using multielement version of get_diff_monodomain" << std::endl;
		DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			get_interpolated_diffusion_matrix(s, D);
	}
	//Bidomain
	void get_intracellular_conductivity_bidomain(const unsigned& ipt,
	                                            const Vector<double> &s,
	                                            const Vector<double>& x,
	                                            DenseMatrix<double>& G_i) const
	{
		DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			get_interpolated_diffusion_matrix(s, G_i);
	}
	void get_extracellular_conductivity_bidomain(const unsigned& ipt,
		                                            const Vector<double> &s,
		                                            const Vector<double>& x,
		                                            DenseMatrix<double>& G_e) const
	{
		DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			get_interpolated_diffusion_matrix(s, G_e);
	}


	//Inform solid element where to get diffusion tensor data from
	void anisotropic_matrix(const unsigned& ipt,
                           const Vector<double> &s,
                           const Vector<double>& xi,
                           DenseMatrix<double>& A) override
    {
    	// oomph_info << "Using multielement version of anisotropic_matrix" << std::endl;
    	DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
    		get_interpolated_preferential_vectors(s, A);
    }


	//Inform solid element where to get active strain from
	void driving_strain(const unsigned& ipt,
                       const Vector<double>& s,
                       const Vector<double>& xi,
                       Vector<double>& V) override
    {
    	// oomph_info << "Using multielement version of driving_strain" << std::endl;
    	// DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
    	// 	get_interpolated_active_strain(s);
    	V.resize(this->dim(), 0.0);
    }



    void fill_in_contribution_to_residuals(Vector<double> &residuals)
	{
		//Call the residuals of the advection-diffusion eqautions
		DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			fill_in_contribution_to_residuals(residuals);
		//Call the residuals of the Navier-Stokes equations
		TAnisotropicPVDElementWithContinuousPressure<DIM>::
			fill_in_contribution_to_residuals(residuals);
	}


	// void fill_in_jacobian_from_solid_position_by_fd(DenseMatrix<double> &jacobian)
	// {
	// 	QAnisotropicPVDElementWithContinuousPressure<DIM>::
	// 		fill_in_jacobian_from_solid_position_by_fd(jacobian);
	// }


	//Tell the elements how to fill in residuals
	//Residuals fill in should be unchanged from just that of the underlying cell element
	void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
												DenseMatrix<double> &jacobian)
	{	
		//Call the residuals of the advection-diffusion eqautions
		DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			fill_in_contribution_to_jacobian(residuals, jacobian);
		//Call the residuals of the Navier-Stokes equations
		TAnisotropicPVDElementWithContinuousPressure<DIM>::
			fill_in_contribution_to_jacobian(residuals, jacobian);
	}

	/// Add the element's contribution to its residuals vector,
	/// jacobian matrix and mass matrix
	void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
															DenseMatrix<double> &jacobian,
															DenseMatrix<double> &mass_matrix)
	{
		DiffAugmentedCell< TConductingCellElement<DIM, 3, CELL_MODEL, CONDUCTANCE_MODEL > >::
			fill_in_contribution_to_jacobian(residuals,jacobian);

		TAnisotropicPVDElementWithContinuousPressure<DIM>::
			fill_in_contribution_to_jacobian(residuals,jacobian);
	}

};



//============================================================================
/// FaceGeometry of a 3D QAnisotropicPVDElement element
//============================================================================
template<unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
class FaceGeometry<TConductingCellElementWithIncompressibleAnisotropicSolidElement<DIM, CELL_MODEL, CONDUCTANCE_MODEL> >:
public virtual SolidTElement<2,3>
{
public:
	FaceGeometry() : SolidTElement<2,3>() {}
};

//============================================================================
/// FaceGeometry of FaceGeometry of a 3D QAnisotropicPVDElement element
//============================================================================
template<unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
class FaceGeometry<FaceGeometry<TConductingCellElementWithIncompressibleAnisotropicSolidElement<DIM, CELL_MODEL, CONDUCTANCE_MODEL> > >:
public virtual SolidTElement<1,3>
{
public:
/// Constructor must call the constructor of the underlying solid element
	FaceGeometry() : SolidTElement<1,3>() {}
};











}//End namespace

#endif