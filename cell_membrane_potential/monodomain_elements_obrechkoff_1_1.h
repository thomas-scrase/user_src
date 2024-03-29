//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations - use strang splitting method
//LIC//====================================================================

//IMPLEMENTS HERMITE SINGLE STEP FORMULATION OF DIFFUSION EQUATION,
// TAKES THE VM PREDICTED BY CELL SOLVE AS THE INITIAL POINT AND FORMULATES
// INTEGRAL OF THE HERMITE INTERPOLATION OF THE FUNCTION OVER THE
// TIME INTERVAL GIVEN BY DT. HERMITE BASIS FUNCTIONS ARE INTEGRATED ANALYTICALLY
//THIS METHOD IS PROBABLY STABLE BECAUSE IT IS IMPLICIT
//3RD (OR 2ND) ORDER ACCURACY, BUMPED UP TO 5TH (OR 4TH) ORDER IF A QUINTIC SPLINE INTERPOLATION
// IS USED.
// THIS HIGH ORDER ACCURACY MAKES IT SUITABLE FOR USE WITH HIGH ORDER OPERATOR
// SPLITTING TECHNIQUES, AN IMPROVEMENT OVER SIMPLE CRANK-NICOLSON/IMPLICIT EULER
// TYPE IMPLEMENTATION

#ifndef OOMPH_MONODOMAIN_HERMITE_APPROXIMATION
#define OOMPH_MONODOMAIN_HERMITE_APPROXIMATION

#include "cell_membrane_potential_elements.h"

namespace oomph{

//Monodomain Equations
template <unsigned DIM>
class MonodomainEquationsHermiteApproximation :
public BaseCellMembranePotentialEquations<DIM>
{
public:

	//change this to take s instead of x? (no functional change, just notation)
	/// \short Funciton pointer to a diffusivity function
	typedef void (*MonodomainEquationsHermiteApproximationDiffFctPt)
	(const Vector<double> &x, DenseMatrix<double> &D);

	//change this to take s instead of x? (no functional change, just notation)
	/// \short Funciton pointer to a diffusivity function
	typedef void (*MonodomainEquationsHermiteApproximationDivDiffFctPt)
	(const Vector<double> &x, Vector<double> &D);


    MonodomainEquationsHermiteApproximation()	:	/*Diff_fct_pt(0),*/ Div_diff_fct_pt(0)
    {
    	this->disable_ALE();
    }


    virtual inline void update_nodal_membrane_potential_BaseCellMembranePotential(const unsigned &l, const double& vm)
	{
		this->node_pt(l)->set_value(this->vm_index_BaseCellMembranePotential(), vm);
	}

	//The membrane potential is just stored at the nodes
	virtual inline double get_nodal_membrane_potential_BaseCellMembranePotential(const unsigned &n) const
	{
		return this->get_nodal_membrane_potential_BaseCellMembranePotential(0, n);
	}

	//Get the t-th history value of membrane potential at the nth node
	virtual inline double get_nodal_membrane_potential_BaseCellMembranePotential(const unsigned &t, const unsigned &n) const
	{
		return this->node_pt(n)->value(t, this->vm_index_BaseCellMembranePotential());
	}

	//Number of second derivatives for the element according to its dimension
	const unsigned d2_dof[3]={1,3,6};
	const unsigned get_d2_dof(){return d2_dof[DIM-1];}


	//Overload the residual for the monodomain equations
	void fill_in_generic_residual_contribution_BaseCellMembranePotential(
	Vector<double> &residuals, DenseMatrix<double> &jacobian, 
	DenseMatrix<double> &mass_matrix, unsigned flag)
	{
		//Find out how many nodes there are
		const unsigned n_node = this->nnode();

		//Get the nodal index at which the unknown is stored
		const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();

		//Set up memory for the shape and test functions
		Shape psi(n_node), test(n_node);
		DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);


		//The number of 
		DShape d2psidx(n_node, 1, get_d2_dof()), d2testdx(n_node, 1, get_d2_dof());

		//Set the value of n_intpt
		const unsigned n_intpt = this->integral_pt()->nweight();

		//Set the Vector to hold local coordinates
		Vector<double> s(DIM);

		//Integers used to store the local equation number and local unknown
		//indices for the residuals and jacobians
		int local_eqn=0, local_unknown=0;
		// oomph_info << n_intpt << std::endl;
		//Loop over the integration points
		for(unsigned ipt=0;ipt<n_intpt;ipt++)
		{
			//Assign values of s
			for(unsigned i=0;i<DIM;i++){s[i] = this->integral_pt()->knot(ipt,i);}

			//Get the integral weight
			double w = this->integral_pt()->weight(ipt);

			if(w<1e-9){continue;}

			//Call the derivatives of the shape and test functions
			double J = this->d2shape_and_d2test_eulerian_at_knot_BaseCellMembranePotential(ipt, 
																							psi, dpsidx, d2psidx,
																							test, dtestdx, d2testdx);

			//Put the 2nd derivatives in matrix form because I can't figure out how to do it otherwise
			Vector<DenseMatrix<double>> m_d2psidx(n_node);
			Vector<DenseMatrix<double>> m_d2testdx(n_node);
			// Loop over nodes
			for(unsigned l=0;l<n_node;l++) 
			{
				m_d2psidx[l].resize(DIM, DIM, 0.0);
				m_d2testdx[l].resize(DIM, DIM, 0.0);

				// Loop over directions
				for(unsigned j=0;j<DIM;j++)
				{
					//Do the diagonals
					m_d2psidx[l](j,j) = d2psidx(l, j);
					m_d2testdx[l](j,j) = d2testdx(l, j);

					for(unsigned i=j+1;i<DIM;i++)
					{
						//Do the upper triangular
						m_d2psidx[l](i,j) = d2psidx(l, i+j);
						m_d2testdx[l](i,j) = d2psidx(l, i+j);
					}

					for(unsigned i=0;i<j;i++)
					{
						//Copy to the lower triangular
						m_d2psidx[l](i,j) = m_d2psidx[l](j,i);
						m_d2testdx[l](i,j) = m_d2testdx[l](j,i);
					}
				}
			}


			//Premultiply the weights and the Jacobian
			double W = w*J;

			//Calculate local values of the solution and its derivatives

			//Coordinate of the integral point
			Vector<double> interpolated_x(DIM,0.0);

			//Interpolated membrane potential at nodes and that given by cells
			double interpolatd_cell_vm = 0.0;
			double interpolated_vm = 0.0;

			//Interpolated gradient of membrane potential at nodes and that given by cells
			Vector<double> interpolated_dcell_vm_dx(DIM, 0.0);
			Vector<double> interpolated_dvm_dx(DIM, 0.0);

			//Interpolated second derivatives of the membrane potential at the nodes
			DenseMatrix<double> interpolated_d2cell_vm_dx(DIM, DIM, 0.0);
			DenseMatrix<double> interpolated_d2vm_dx(DIM, DIM, 0.0);


			//We assume that dt is constant over the nodes
			const double dt = this->node_pt(0)->time_stepper_pt()->time_pt()->dt();

			//Calculate function value and derivatives:
			//-----------------------------------------
			// Loop over nodes
			for(unsigned l=0;l<n_node;l++) 
			{
				//If we're feeling paranoid then we check to ensure dt is infact the same for all the nodes
				#ifdef PARANOID
				if(std::fabs(dt-this->node_pt(0)->time_stepper_pt()->time_pt()->dt())>1e-9)
				{
					throw OomphLibError(
									"dt is not the same over all the nodes, but it has to be for this implementation",
									OOMPH_CURRENT_FUNCTION,
									OOMPH_EXCEPTION_LOCATION);
				}
				#endif

				//Get the interpolated values at the nodes
				double vm_value = this->raw_nodal_value(l,vm_nodal_index);
				interpolated_vm += vm_value*psi(l);

				//Get the interpolated membrane potential from the cells
				double cell_vm_value = this->get_nodal_predicted_vm_BaseCellMembranePotential(l);
				interpolatd_cell_vm += cell_vm_value*psi(l);

				// Loop over directions
				for(unsigned j=0;j<DIM;j++)
				{
					//Get the interpolated global coordinate from the nodes
					interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);

					//Get the interpolated gradients of the membrane potential from the nodes and the cells
					interpolated_dvm_dx[j] += vm_value*dpsidx(l,j);
					interpolated_dcell_vm_dx[j] += cell_vm_value*dpsidx(l,j);

					// for(unsigned i=j+1;i<DIM;i++)
					for(unsigned i=0;i<DIM;i++)
					{
						interpolated_d2vm_dx(j,i) += vm_value*m_d2psidx[l](i,j);

						interpolated_d2cell_vm_dx(j,i) += cell_vm_value*m_d2psidx[l](i,j);
					}
				}
			}

			// Mesh velocity?
			if (!this->ALE_is_disabled)
			{
				// for(unsigned l=0;l<n_node;l++) 
				// {
				// 	for(unsigned j=0;j<DIM;j++)
				// 	{
				// 		mesh_velocity[j] += this->raw_dnodal_position_dt(l,j)*psi(l);
				// 	}
				// }
				throw OomphLibError("Monodomain elements with hermite approximation are not implemented with ALE yet",
									OOMPH_CURRENT_FUNCTION,
									OOMPH_EXCEPTION_LOCATION);
			}

			//Get diffusivity tensor
			DenseMatrix<double> D(DIM, DIM, 0.0);
			this->get_diff_BaseCellMembranePotential(ipt,s,interpolated_x,D);

			//Get the divergence of the diffusivity tensor
			Vector<double> divD(DIM, 0.0);
			this->get_div_diff_monodomain(ipt,s,interpolated_x,divD);

			// for(unsigned i=0;i<DIM;i++)
			// {
			// 	oomph_info << divD[i] << std::endl;
			// }
			// exit(0);



			// Assemble residuals and Jacobian
			//--------------------------------

			// Loop over the test functions
			for(unsigned l=0;l<n_node;l++)
			{
				//Fill in residual and jacobian contribution for u
				//Set the local equation number
				local_eqn = this->nodal_local_eqn(l,this->vm_index_BaseCellMembranePotential());

				/*IF it's not a boundary condition*/
				if(local_eqn >= 0)
				{
					//Add contribution due to time derivative
					residuals[local_eqn] += -(interpolated_vm - interpolatd_cell_vm)*test(l)*W;


					double tmp2 = 0.0;
					double tmp3 = 0.0;

					for(unsigned i=0;i<DIM;i++)
					{
						double tmp1 = 0.0;

						tmp2 += divD[i]*(interpolated_dcell_vm_dx[i] - interpolated_dvm_dx[i]);

						tmp3 += divD[i]*dtestdx(l,i);

						for(unsigned j=0;j<DIM;j++)
						{
							//tmp1[i] but we don't bother making it a vector
							tmp1 += D(i,j)*(interpolated_dcell_vm_dx[j] + interpolated_dvm_dx[j]);

							tmp2 += D(i,j)*(interpolated_d2cell_vm_dx(i,j) - interpolated_d2vm_dx(i,j));

							tmp3 += D(i,j)*m_d2testdx[l](i,j);
						}
						residuals[local_eqn] += (dt/2.0)*tmp1*dtestdx(l, i)*W;
					}

					residuals[local_eqn] += (dt*dt/12.0)*tmp2*tmp3*W;

					// Calculate the jacobian
					//-----------------------
					if(flag)
					{
						//Loop over the velocity shape functions again
						for(unsigned l2=0;l2<n_node;l2++)
						{ 
							//Set the number of the unknown
							local_unknown = this->nodal_local_eqn(l2,this->vm_index_BaseCellMembranePotential());

							//If at a non-zero degree of freedom add in the entry
							if(local_unknown >= 0)
							{
								jacobian(local_eqn, local_unknown) += -test(l)*psi(l2)*W;

								//Add the mass matrix term
								if(flag==2)
								{
									mass_matrix(local_eqn,local_unknown) += -test(l)*psi(l2)*W;
								}

								double tmp2 = 0.0;
								double tmp3 = 0.0;

								//Add contribution to Elemental Matrix
								for(unsigned i=0;i<DIM;i++)
								{
									double tmp1 = 0.0;
								
									tmp2 += -divD[i]*dpsidx(l2, i);

									tmp3 += divD[i]*dtestdx(l, i);

									for(unsigned j=0;j<DIM;j++)
									{
										tmp1 += D(i,j)*dpsidx(l2,j);

										tmp2 += -D(i,j)*m_d2psidx[l2](i, j);

										tmp3 += D(i,j)*m_d2testdx[l](i,j);
									}

									jacobian(local_eqn, local_unknown) += (dt/2.0)*tmp1*dtestdx(l,i)*W;
								}

								jacobian(local_eqn, local_unknown) += (dt*dt/12.0)*tmp2*tmp3*W;
							}
						}
					}
				}
			}
		} // End of loop over integration points
	}
	
	inline unsigned required_nvalue(const unsigned &n) const 
  		{return 1;}

	// /// Access function: Pointer to diffusion  function
	// MonodomainEquationsHermiteApproximationDiffFctPt& diff_fct_pt() 
	// {return Diff_fct_pt;}

	// /// Access function: Pointer to diffusion function. Const version
	// MonodomainEquationsHermiteApproximationDiffFctPt diff_fct_pt() const 
	// {return Diff_fct_pt;}


	/// Access function: Pointer to diffusion  function
	MonodomainEquationsHermiteApproximationDivDiffFctPt& div_diff_fct_pt() 
	{return Div_diff_fct_pt;}

	/// Access function: Pointer to diffusion function. Const version
	MonodomainEquationsHermiteApproximationDivDiffFctPt div_diff_fct_pt() const 
	{return Div_diff_fct_pt;}


	void assign_additional_initial_conditions(const unsigned &l)
	{
		//Do nothing - we have no additional variables to assign values to
	}

	void assign_initial_conditions_consistent_with_cell_model(const unsigned &l, const double& vm)
	{
		this->node_pt(l)->set_value(this->vm_index_BaseCellMembranePotential(), vm);
	}


	//Calculate dvmdt but instead of using the stored value for the previous timestep, use the value predicted by the cell model
	double dvm_dt_Strang_Split(const unsigned &n) const
	{
		// Get the data's timestepper
		TimeStepper* time_stepper_pt= this->node_pt(n)->time_stepper_pt();

		//Initialise dudt
		double dvmdt=0.0;
		//Loop over the timesteps, if there is a non Steady timestepper
		if (!time_stepper_pt->is_steady())
		{
		 //Find the index at which the variable is stored
		 const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();

		 // Number of timsteps (past & present)
		 const unsigned n_time = time_stepper_pt->ntstorage();
		 
		 //Start at the one before the previous timestep
		 for(unsigned t=2;t<n_time;t++)
		  {
		   dvmdt += time_stepper_pt->weight(1,t)*this->nodal_value(t,n,vm_nodal_index);
		  }
		  //Add the previous timestep and the current timestep
		  dvmdt += time_stepper_pt->weight(1,1)*this->get_nodal_predicted_vm_BaseCellMembranePotential(n);

		  dvmdt += time_stepper_pt->weight(1,0)*this->nodal_value(0,n,vm_nodal_index);
		}
		return dvmdt;
	}


	// /// \short Get diffusivity tensor at (Eulerian) position 
	// /// x and/or local coordinate s. 
	// /// This function is
	// /// virtual to allow overloading in multi-physics problems where
	// /// the diff function might be determined by
	// /// another system of equations 
	// inline virtual void get_diff_monodomain(const unsigned& ipt,
 //                                            const Vector<double> &s,
 //                                            const Vector<double>& x,
 //                                            DenseMatrix<double>& D) const
	// {
	// 	//If no diff function has been set, return identity
	// 	if(Diff_fct_pt==0){
	// 		// oomph_info << "Not using fct pt" << std::endl;
	// 		for(unsigned i=0; i<DIM; i++){
	// 			for(unsigned j=0; j<DIM; j++){
	// 				D(i,j) =  0.0;
	// 			}
	// 			D(i,i)  = 1.0;
	// 		}
	// 	}
	// 	else{
	// 		// oomph_info << "Using fct pt" << std::endl;
	// 		// Get diffusivity tensor from function
	// 		(*Diff_fct_pt)(x,D);
	// 	}
	// }


	//Get the 'divergence' of the diffusion tensor divD[i] = dD(j,i)/dx_j
	//By default it finite differences the diffusion tensor
	inline virtual void get_div_diff_monodomain(const unsigned& ipt,
                                            const Vector<double> &s,
                                            const Vector<double>& x,
                                            Vector<double>& divD) const
	{
		//If no diff function has been set, return the zero vector
		if(this->diff_fct_pt()==0)
		{
			for(unsigned i=0; i<DIM; i++)
			{
				divD[i] =  0.0;
			}
			return;
		}
		
		if(Div_diff_fct_pt==0)//otherwise if no div function is set perform finite differencing of the diffusion tensor function
		{
			//Zero the vector
			for(unsigned i=0; i<DIM; i++)
			{
				divD[i] =  0.0;
			}

			Vector<double> x_fd = x;

			DenseMatrix<double> D0(DIM);
			(*this->diff_fct_pt())(x,D0);

			for(unsigned i=0; i<DIM; i++)
			{

				x_fd[i] += this->Default_fd_jacobian_step;

				DenseMatrix<double> D1(DIM);
				(*this->diff_fct_pt())(x_fd,D1);

				for(unsigned j=0; j<DIM; j++)
				{
					divD[i] += (D1(j,i)-D0(j,i))/this->Default_fd_jacobian_step;
				}
				x_fd[i] -= this->Default_fd_jacobian_step;
			}
		}
		else//Otherwise, if the div function is set, use it
		{
			(*Div_diff_fct_pt)(x,divD);
		}
	}


	/// Get flux: \f$\mbox{flux}[i] = \mbox{d}u / \mbox{d}x_i \f$
	void get_total_flux_monodomain(const Vector<double>& s, 
	                			Vector<double>& total_flux) const
	{
		//Find out how many nodes there are in the element
		const unsigned n_node = this->nnode();

		//Get the nodal index at which the unknown is stored
		const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();

		//Set up memory for the shape and test functions
		Shape psi(n_node);
		DShape dpsidx(n_node,DIM);

		//Call the derivatives of the shape and test functions
		this->dshape_eulerian(s,psi,dpsidx);

		//Storage for the Eulerian position
		Vector<double> interpolated_x(DIM,0.0);
		//Storage for the concentration
		double interpolated_vm = 0.0;
		//Storage for the derivatives of the concentration
		Vector<double> interpolated_dvmdx(DIM,0.0);

		// Loop over nodes
		for(unsigned l=0;l<n_node;l++) 
		{
		 //Get the value at the node
		 const double vm_value = this->nodal_value(l,vm_nodal_index);
		 interpolated_vm += vm_value*psi(l);
		 // Loop over directions
		 for(unsigned j=0;j<DIM;j++)
		  {
		   interpolated_x[j] += this->nodal_position(l,j)*psi(l);
		   interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
		  }
		}

		//Dummy integration point 
		unsigned ipt=0;

		//Get diffusivity tensor
		DenseMatrix<double> D(DIM,DIM);
		this->get_diff_BaseCellMembranePotential(ipt,s,interpolated_x,D);

		//Calculate the total flux made up of the diffusive flux
		//and the conserved wind
		for(unsigned i=0;i<DIM;i++)
		{                               
		 total_flux[i] = 0.0;
		 for(unsigned j=0;j<DIM;j++)
		  {
		   total_flux[i] +=  D(i,j)*interpolated_dvmdx[j];
		  }
		}
	}

	inline std::vector<std::string> get_variable_names_BaseCellMembranePotentialEquations() const override
	{
		std::vector<std::string> v = {"Vm"};
		return v;
	}


	unsigned nscalar_paraview() const
	{
		return 1;
	}

	void scalar_value_paraview(std::ofstream& file_out,
								const unsigned& i,
								const unsigned& nplot) const
	{
		if(i!=0)
		{
			throw OomphLibError(
				"Monodomain element has only transmembrane potential to report",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}
		//Vector of local coordinates
		Vector<double> s(DIM);

		const unsigned n_node = this->nnode();
		// const unsigned vm_index = this->vm_index_BaseCellMembranePotential();
		Shape psi(n_node);
		DShape dpsidx(n_node,DIM);

		// Loop over plot points
		unsigned num_plot_points=this->nplot_points_paraview(nplot);
		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
			// Get local coordinates of plot point
			this->get_s_plot(iplot,nplot,s);

			file_out << this->get_interpolated_membrane_potential_BaseCellMembranePotential(s) << std::endl;
		}
	}

	void scalar_value_fct_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot,
									FiniteElement::SteadyExactSolutionFctPt
									exact_soln_pt) const
	{
		scalar_value_paraview(file_out,i,nplot);
	}

	std::string scalar_name_paraview(const unsigned& i) const
	{
		return "Transmembrane potential";
	}


	/////////////////////////////////////////////////////////////////////////////////
	//Residual and Jacobian functions
	// The conductance model is the only thing that does actual oomph lib solving
	/////////////////////////////////////////////////////////////////////////////////
	/// Add the element's contribution to its residual vector (wrapper)
	void fill_in_contribution_to_residuals(Vector<double> &residuals)
	{
		//Call the generic residuals function with flag set to 0 and using
		//a dummy matrix
		fill_in_generic_residual_contribution_BaseCellMembranePotential(
		    residuals, GeneralisedElement::Dummy_matrix, GeneralisedElement::Dummy_matrix, 0);
	}

	/// \short Add the element's contribution to its residual vector and 
	/// the element Jacobian matrix (wrapper)
	void fill_in_contribution_to_jacobian(Vector<double> &residuals,
	DenseMatrix<double> &jacobian)
	{
		fill_in_generic_residual_contribution_BaseCellMembranePotential(
		    residuals, jacobian, GeneralisedElement::Dummy_matrix, 1);

		// FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
	}

	/// Add the element's contribution to its residuals vector,
	/// jacobian matrix and mass matrix
	void fill_in_contribution_to_jacobian_and_mass_matrix(
	Vector<double> &residuals, DenseMatrix<double> &jacobian,
	DenseMatrix<double> &mass_matrix)
	{
		//Call the generic routine with the flag set to 2
		// fill_in_generic_residual_contribution_cell_interface(residuals,jacobian,mass_matrix,2);
		fill_in_generic_residual_contribution_BaseCellMembranePotential(
		    residuals, jacobian, mass_matrix, 2);
	}



protected:
	/// Pointer to diffusivity function:
	///		function typedef is given in BaseCellMembranePotentialEquations
	// MonodomainEquationsHermiteApproximationDiffFctPt Diff_fct_pt;

	MonodomainEquationsHermiteApproximationDivDiffFctPt Div_diff_fct_pt;	
};




//Monodomain Elements

//======================================================================
/// \short QMonodomainElementHermiteApproximation elements are 
/// linear/quadrilateral/brick-shaped Advection Diffusion elements with 
/// isoparametric interpolation for the function.
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
 class QMonodomainElementHermiteApproximation : 
 public virtual QElement<DIM,NNODE_1D>,
 public virtual MonodomainEquationsHermiteApproximation<DIM>
{

private:

 /// \short Static array of ints to hold number of variables at 
 /// nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue = 1;
 
  public:


 ///\short  Constructor: Call constructors for QElement and 
 /// Advection Diffusion equations
 QMonodomainElementHermiteApproximation() : QElement<DIM,NNODE_1D>(), 
  MonodomainEquationsHermiteApproximation<DIM>()
  { }

 /// Broken copy constructor
 QMonodomainElementHermiteApproximation(
  const QMonodomainElementHermiteApproximation<DIM,NNODE_1D>&  dummy) 
  { 
   BrokenCopy::broken_copy("QMonodomainElementHermiteApproximation");
  } 
 
 /// Broken assignment operator
 void operator=(const QMonodomainElementHermiteApproximation<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("QMonodomainElementHermiteApproximation");
  }

 /// \short  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {BaseCellMembranePotentialEquations<DIM>::output(outfile);}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &nplot) override
 {
 	// std::cout << "BOOM\nBOOM\nBOOM\nBOOM\nBOOM" << std::endl;
	//Vector of local coordinates
	Vector<double> s(DIM);

	// Tecplot header info
	outfile << this->tecplot_zone_string(nplot);

	const unsigned n_node = this->nnode();
	const unsigned vm_index = this->vm_index_BaseCellMembranePotential();
	// std::cout << "n_node " << n_node << std::endl;
	Shape psi(n_node);
	DShape dpsidx(n_node,DIM);

	// Loop over plot points
	unsigned num_plot_points=this->nplot_points(nplot);
	for (unsigned iplot=0;iplot<num_plot_points;iplot++){
		// Get local coordinates of plot point
		this->get_s_plot(iplot,nplot,s);

		// Get Eulerian coordinate of plot point
		Vector<double> x(DIM);
		this->interpolated_x(s,x);

		for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}
		outfile << this->get_interpolated_membrane_potential_BaseCellMembranePotential(s) << " ";

		//Get the gradients
		(void)this->dshape_eulerian(s,psi,dpsidx);
		Vector<double> interpolated_dvmdx(DIM,0.0);
		double dvmdt = 0.0;
		for(unsigned n=0;n<n_node;n++){
			const double vm_ = this->nodal_value(n,vm_index);
			dvmdt += this->dvm_dt_BaseCellMembranePotential(n)*psi(n);
			for(unsigned i=0;i<DIM;i++){interpolated_dvmdx[i] += vm_*dpsidx(n,i);}
		}

		outfile << dvmdt << " ";

		for(unsigned i=0;i<DIM;i++){outfile << interpolated_dvmdx[i]  << " ";}

		//Get diffusivity tensor
		DenseMatrix<double> D(DIM,DIM,0.0);
		this->get_diff_BaseCellMembranePotential(iplot,s,x,D);
		for(unsigned i=0; i<DIM; i++){
			for(unsigned j=0; j<=i; j++){
				outfile << D(i,j) << " ";
			}
		}

		outfile  << std::endl;
	}
	// Write tecplot footer (e.g. FE connectivity lists)
	this->write_tecplot_zone_footer(outfile,nplot);
}


 /// \short C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {
   BaseCellMembranePotentialEquations<DIM>::output(file_pt);
  }

 ///  \short C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {
   BaseCellMembranePotentialEquations<DIM>::output(file_pt,n_plot);
  }

 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
  {BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,exact_soln_pt);}


 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt 
                 exact_soln_pt)
  {
   BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,time,exact_soln_pt);
  }


protected:

 	inline double dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s, 
                                               Shape &psi, 
                                               DShape &dpsidx, 
                                               Shape &test, 
                                               DShape &dtestdx) const;

 	inline double dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt,
												Shape &psi, 
												DShape &dpsidx,
												Shape &test, 
												DShape &dtestdx) const;

 	// inline double d2shape_and_d2test_eulerian_BaseCellMembranePotential(const Vector<double> &s, 
  //                                              Shape &psi, 
  //                                              DShape &dpsidx, 
  //                                              Shape &test, 
  //                                              DShape &dtestdx) const;

 	inline double d2shape_and_d2test_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt, 
																				Shape &psi, 
																				DShape &dpsidx, 
																				DShape &d2psidx,
																				Shape &test, 
																				DShape &dtestdx, 
																				DShape &d2testdx) const;

};

//Inline functions:


//======================================================================
/// \short Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double QMonodomainElementHermiteApproximation<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s,
                                         Shape &psi, 
                                         DShape &dpsidx,
                                         Shape &test, 
                                         DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian(s,psi,dpsidx);

 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<NNODE_1D;i++)
  {
   test[i] = psi[i]; 
   for(unsigned j=0;j<DIM;j++)
    {
     dtestdx(i,j) = dpsidx(i,j);
    }
  }
 
 //Return the jacobian
 return J;
}



//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double QMonodomainElementHermiteApproximation<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(
 const unsigned &ipt,
 Shape &psi, 
 DShape &dpsidx,
 Shape &test, 
 DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Set the test functions equal to the shape functions (pointer copy)
 test = psi;
 dtestdx = dpsidx;

 //Return the jacobian
 return J;
}
// template<unsigned DIM, unsigned NNODE_1D>
// double QMonodomainElementHermiteApproximation<DIM, NNODE_1D>::
// d2shape_and_d2test_eulerian_BaseCellMembranePotential(const Vector<double> &s, 
//                                                Shape &psi, 
//                                                DShape &dpsidx, 
//                                                Shape &test, 
//                                                DShape &dtestdx) const
// {

// }
template<unsigned DIM, unsigned NNODE_1D>
double QMonodomainElementHermiteApproximation<DIM, NNODE_1D>::
d2shape_and_d2test_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt, 
															Shape &psi, 
															DShape &dpsidx, 
															DShape &d2psidx,
															Shape &test, 
															DShape &dtestdx, 
															DShape &d2testdx) const
{
	const double J = this->d2shape_eulerian_at_knot(ipt, psi, dpsidx, d2psidx);
	test = psi;
	dtestdx = dpsidx;
	d2testdx = d2psidx;

	return J;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======================================================================
/// \short Face geometry for the QMonodomainElementHermiteApproximation elements: 
/// The spatial dimension of the face elements is one lower than that 
/// of the bulk element but they have the same number of points along 
/// their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<QMonodomainElementHermiteApproximation<DIM,NNODE_1D> >: 
 public virtual QElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the 1D QMonodomain elements: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<QMonodomainElementHermiteApproximation<1,NNODE_1D> >: 
 public virtual PointElement
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : PointElement() {}

};


//Override functions in specific implementations of the diff augmented wrapper
template<unsigned DIM, unsigned NNODE_1D>
class DiffAugmentedCell<QMonodomainElementHermiteApproximation<DIM, NNODE_1D>>:
	public TMonodomainElement<DIM, NNODE_1D>
{
public:

	inline void get_diff_BaseCellMembranePotential(const unsigned& ipt,
                                    const Vector<double> &s,
                                    const Vector<double>& x,
                                    DenseMatrix<double>& D) const
	{
		this->get_interpolated_diffusion_matrix(s, D);
	}

};





/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// TMonodomainElementHermiteApproximation
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//======================================================================
/// \short TMonodomainElementHermiteApproximation elements are isoparametric triangular 
/// DIM-dimensional General Advection Diffusion Equations with  NNODE_1D nodal points along each
/// element edge. Inherits from TElement and MonodomainEquationsHermiteApproximation
//======================================================================
template <unsigned DIM, unsigned NNODE_1D>
 class TMonodomainElementHermiteApproximation : 
 public virtual TElement<DIM,NNODE_1D>,
 public virtual MonodomainEquationsHermiteApproximation<DIM>
{

private:

 /// \short Static array of ints to hold number of variables at 
 /// nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue = 1;
 
  public:


 ///\short  Constructor: Call constructors for TElement and 
 /// Advection Diffusion equations
 TMonodomainElementHermiteApproximation() : TElement<DIM,NNODE_1D>(), 
  MonodomainEquationsHermiteApproximation<DIM>()
  { }

 /// Broken copy constructor
 TMonodomainElementHermiteApproximation(
  const TMonodomainElementHermiteApproximation<DIM,NNODE_1D>&  dummy) 
  { 
   BrokenCopy::broken_copy("TMonodomainElementHermiteApproximation");
  } 
 
 /// Broken assignment operator
 void operator=(const TMonodomainElementHermiteApproximation<DIM,NNODE_1D>&) 
  {
   BrokenCopy::broken_assign("TMonodomainElementHermiteApproximation");
  }

 /// \short  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {BaseCellMembranePotentialEquations<DIM>::output(outfile);}

 /// \short Output function:  
 ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {BaseCellMembranePotentialEquations<DIM>::output(outfile,n_plot);}


 /// \short C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {
   BaseCellMembranePotentialEquations<DIM>::output(file_pt);
  }

 ///  \short C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {
   BaseCellMembranePotentialEquations<DIM>::output(file_pt,n_plot);
  }

 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
  {BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,exact_soln_pt);}


 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt 
                 exact_soln_pt)
  {
   BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,time,exact_soln_pt);
  }


protected:

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_and_dtest_eulerian_BaseCellMembranePotential(
  const Vector<double> &s, 
  Shape &psi, 
  DShape &dpsidx, 
  Shape &test, 
  DShape &dtestdx) const;
 
 /// \short Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(
  const unsigned& ipt,
  Shape &psi, 
  DShape &dpsidx, 
  Shape &test,
  DShape &dtestdx) 
  const;


  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 // inline double d2shape_and_d2test_eulerian_BaseCellMembranePotential(
 //  const Vector<double> &s, 
 //  Shape &psi, 
 //  DShape &dpsidx, 
 //  Shape &test, 
 //  DShape &dtestdx) const;
 
 /// \short Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt, 
																			Shape &psi, 
																			DShape &dpsidx, 
																			DShape &d2psidx,
																			Shape &test, 
																			DShape &dtestdx, 
																			DShape &d2testdx) const;

};

//Inline functions:


//======================================================================
/// \short Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double TMonodomainElementHermiteApproximation<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s,
                                         Shape &psi, 
                                         DShape &dpsidx,
                                         Shape &test, 
                                         DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian(s,psi,dpsidx);

 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<NNODE_1D;i++)
  {
   test[i] = psi[i]; 
   for(unsigned j=0;j<DIM;j++)
    {
     dtestdx(i,j) = dpsidx(i,j);
    }
  }
 
 //Return the jacobian
 return J;
}



//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
double TMonodomainElementHermiteApproximation<DIM,NNODE_1D>::
 dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(
 const unsigned &ipt,
 Shape &psi, 
 DShape &dpsidx,
 Shape &test, 
 DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Set the test functions equal to the shape functions (pointer copy)
 test = psi;
 dtestdx = dpsidx;

 //Return the jacobian
 return J;
}

// template<unsigned DIM, unsigned NNODE_1D>
// double QMonodomainElementHermiteApproximation<DIM, NNODE_1D>::
// d2shape_and_d2test_eulerian_BaseCellMembranePotential(const Vector<double> &s, 
//                                                Shape &psi, 
//                                                DShape &dpsidx, 
//                                                Shape &test, 
//                                                DShape &dtestdx) const
// {

// }

template<unsigned DIM, unsigned NNODE_1D>
double TMonodomainElementHermiteApproximation<DIM, NNODE_1D>::
d2shape_and_d2test_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt, 
															  Shape &psi, 
															  DShape &dpsidx, 
															  DShape &d2psidx,
															  Shape &test, 
															  DShape &dtestdx, 
															  DShape &d2testdx) const
{
	const double J = this->d2shape_eulerian_at_knot(ipt, psi, dpsidx, d2psidx);

	test = psi;
	dtestdx = dpsidx;
	d2testdx = d2psidx;

	return J;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======================================================================
/// \short Face geometry for the TMonodomainElementHermiteApproximation elements: 
/// The spatial dimension of the face elements is one lower than that 
/// of the bulk element but they have the same number of points along 
/// their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<TMonodomainElementHermiteApproximation<DIM,NNODE_1D> >: 
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the 1D TMonodomainElementHermiteApproximation: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<TMonodomainElementHermiteApproximation<1,NNODE_1D> >: 
 public virtual PointElement
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : PointElement() {}

};


//Override functions in specific implementations of the diff augmented wrapper
template<unsigned DIM, unsigned NNODE_1D>
class DiffAugmentedCell<TMonodomainElementHermiteApproximation<DIM, NNODE_1D>>:
	public TMonodomainElement<DIM, NNODE_1D>
{
public:

	inline void get_diff_BaseCellMembranePotential(const unsigned& ipt,
                                    const Vector<double> &s,
                                    const Vector<double>& x,
                                    DenseMatrix<double>& D) const
	{
		this->get_interpolated_diffusion_matrix(s, D);
	}

};



}


#endif