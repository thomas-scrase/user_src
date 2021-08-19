//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations
//LIC//====================================================================

#ifndef OOMPH_MONODOMAIN_INTEGRAL_SPLITTING_METHOD
#define OOMPH_MONODOMAIN_INTEGRAL_SPLITTING_METHOD

#include "cell_membrane_potential_elements.h"

namespace oomph{

//Monodomain Equations
	template <unsigned DIM>
	class MonodomainEquationsIntegralSplittingMethod :
	public BaseCellMembranePotentialEquations<DIM>
	{
	public:

		//change this to take s instead of x? (no functional change, just notation)
		/// \short Funciton pointer to a diffusivity function
		typedef void (*MonodomainEquationsIntegralSplittingMethodDiffFctPt)
		(const Vector<double> &x, DenseMatrix<double> &D);


        MonodomainEquationsIntegralSplittingMethod()	:	Diff_fct_pt(0)
        {

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
		   
		   //Set the value of n_intpt
		   const unsigned n_intpt = this->integral_pt()->nweight();

		   //Set the Vector to hold local coordinates
		   Vector<double> s(DIM);

		   //Integers used to store the local equation number and local unknown
		   //indices for the residuals and jacobians
		   int local_eqn=0, local_unknown=0;

		   //Loop over the integration points
		   for(unsigned ipt=0;ipt<n_intpt;ipt++)
		    {

		     //Assign values of s
		     for(unsigned i=0;i<DIM;i++) s[i] = this->integral_pt()->knot(ipt,i);

		     //Get the integral weight
		     double w = this->integral_pt()->weight(ipt);

		      if(w<=1e-9){/* oomph_info << ipt << " continue."<<std::endl;*/ continue;}
		      
		     //Call the derivatives of the shape and test functions
		     double J = 
		      this->dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(ipt,psi,dpsidx,test,dtestdx);
		         
		     //Premultiply the weights and the Jacobian
		     double W = w*J;

		     //Calculate local values of the solution and its derivatives
		     //Allocate
		     double interpolated_vm=0.0;
		     double dvmdt=0.0;

		     Vector<double> interpolated_x(DIM,0.0);
		     Vector<double> interpolated_dvmdx(DIM,0.0);
		     Vector<double> mesh_velocity(DIM,0.0);

		     //Get source function
		     //-------------------
		     double source = 0.0;

		     //Calculate function value and derivatives:
		     //-----------------------------------------
		     // Loop over nodes
		     for(unsigned l=0;l<n_node;l++) 
		      {
		       //Get the value at the node
		       double vm_value = this->raw_nodal_value(l,vm_nodal_index);
		       interpolated_vm += vm_value*psi(l);

		       //The average source over the timestep
		       source -= this->get_nodal_integral_iion_BaseCellMembranePotential(l)*psi(l)/(this->node_pt(l)->time_stepper_pt()->time_pt()->dt());


		       dvmdt += this->dvm_dt_BaseCellMembranePotential(l)*psi(l);

		       // Loop over directions
		       for(unsigned j=0;j<DIM;j++)
		        {
		         interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
		         interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
		        }
		      }
		     
		     // Mesh velocity?
		     if (!this->ALE_is_disabled)
		      {
		       for(unsigned l=0;l<n_node;l++) 
		        {
		         for(unsigned j=0;j<DIM;j++)
		          {
		           mesh_velocity[j] += this->raw_dnodal_position_dt(l,j)*psi(l);
		          }
		        }
		      }

		     //Get source function
		     //-------------------
		     // double source;
		     // this->get_source_BaseCellMembranePotential(ipt,s,interpolated_x,source);


		     //Get diffusivity tensor
		     DenseMatrix<double> D(DIM,DIM,0.0);
		     this->get_diff_monodomain(ipt,s,interpolated_x,D);

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
		          // Add body force/source term and time derivative
		          residuals[local_eqn] -= (dvmdt + source)*test(l)*W;
		         
		          // The Generalised Advection Diffusion bit itself
		          for(unsigned k=0;k<DIM;k++)
		            {
		             //Terms that multiply the test function 
		              double tmp = 0.0;
		             // //If the mesh is moving need to subtract the mesh velocity
		             if(!this->ALE_is_disabled) {tmp -= mesh_velocity[k];}
		             tmp *= interpolated_dvmdx[k];

		             //Terms that multiply the derivative of the test function
		              double tmp2 = 0.0;
		             //Now the diuffusive term
		             for(unsigned j=0;j<DIM;j++)
		              {
		               tmp2 += interpolated_dvmdx[j]*D(k,j);
		              }
		             //Now construct the contribution to the residuals
		             residuals[local_eqn] -= (tmp*test(l) + tmp2*dtestdx(l,k))*W;
		            }

		         
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
		                //Mass matrix term
		                jacobian(local_eqn,local_unknown) 
		                  -= test(l)*psi(l2)*
		                  this->node_pt(l2)->time_stepper_pt()->weight(1,0)*W;

		                //Add the mass matrix term
		                if(flag==2)
		                  {
		                  mass_matrix(local_eqn,local_unknown)
		                    += test(l)*psi(l2)*W;
		                  }

		                //Add contribution to Elemental Matrix
		                for(unsigned k=0;k<DIM;k++){
		                  //Temporary term used in assembly
		                  double tmp = 0.0;
		                  if(!this->ALE_is_disabled)
		                   {tmp -= mesh_velocity[k];}
		                  tmp *= dpsidx(l2,k);

		                  double tmp2 = 0.0;
		                  //Now the diffusive term
		                  for(unsigned j=0;j<DIM;j++)
		                    {
		                    tmp2 += D(k,j)*dpsidx(l2,j);
		                    }
		                 
		                  //Now assemble Jacobian term
		                  jacobian(local_eqn,local_unknown) 
		                    -= (tmp*test(l) + tmp2*dtestdx(l,k))*W;
		                  }
		                }
		              }
		            }
		          }
		        }
		      } // End of loop over integration points
		    }
	
		inline unsigned required_nvalue(const unsigned &n) const 
	  		{return 1;}

		/// Access function: Pointer to diffusion  function
		MonodomainEquationsIntegralSplittingMethodDiffFctPt& diff_fct_pt() 
		{return Diff_fct_pt;}

		/// Access function: Pointer to diffusion function. Const version
		MonodomainEquationsIntegralSplittingMethodDiffFctPt diff_fct_pt() const 
		{return Diff_fct_pt;}


		void assign_additional_initial_conditions()
		{
			//Do nothing - we have no additional variables to assign values to
		}

		void assign_initial_conditions_consistent_with_cell_model(const unsigned &l, const double& vm)
		{
			this->node_pt(l)->set_value(this->vm_index_BaseCellMembranePotential(), vm);
		}


		/// \short Get diffusivity tensor at (Eulerian) position 
		/// x and/or local coordinate s. 
		/// This function is
		/// virtual to allow overloading in multi-physics problems where
		/// the diff function might be determined by
		/// another system of equations 
		inline virtual void get_diff_monodomain(const unsigned& ipt,
	                                            const Vector<double> &s,
	                                            const Vector<double>& x,
	                                            DenseMatrix<double>& D) const
		{
			//If no diff function has been set, return identity
			if(Diff_fct_pt==0){
				// oomph_info << "Not using fct pt" << std::endl;
				for(unsigned i=0; i<DIM; i++){
					for(unsigned j=0; j<DIM; j++){
						D(i,j) =  0.0;
					}
					D(i,i)  = 1.0;
				}
			}
			else{
				// oomph_info << "Using fct pt" << std::endl;
				// Get diffusivity tensor from function
				(*Diff_fct_pt)(x,D);
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
			get_diff_monodomain(ipt,s,interpolated_x,D);

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

		inline std::vector<std::string> get_variable_names() const override
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
		MonodomainEquationsIntegralSplittingMethodDiffFctPt Diff_fct_pt;

};


//Monodomain Elements

	//======================================================================
	/// \short QMonodomainElementIntegralSplittingMethod elements are 
	/// linear/quadrilateral/brick-shaped Advection Diffusion elements with 
	/// isoparametric interpolation for the function.
	//======================================================================
	template <unsigned DIM, unsigned NNODE_1D>
	 class QMonodomainElementIntegralSplittingMethod : 
	 public virtual QElement<DIM,NNODE_1D>,
	 public virtual MonodomainEquationsIntegralSplittingMethod<DIM>
	{

	private:

	 /// \short Static array of ints to hold number of variables at 
	 /// nodes: Initial_Nvalue[n]
	 static const unsigned Initial_Nvalue = 1;
	 
	  public:


	 ///\short  Constructor: Call constructors for QElement and 
	 /// Advection Diffusion equations
	 QMonodomainElementIntegralSplittingMethod() : QElement<DIM,NNODE_1D>(), 
	  MonodomainEquationsIntegralSplittingMethod<DIM>()
	  { }

	 /// Broken copy constructor
	 QMonodomainElementIntegralSplittingMethod(
	  const QMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>&  dummy) 
	  { 
	   BrokenCopy::broken_copy("QMonodomainElementIntegralSplittingMethod");
	  } 
	 
	 /// Broken assignment operator
	 void operator=(const QMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>&) 
	  {
	   BrokenCopy::broken_assign("QMonodomainElementIntegralSplittingMethod");
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
			outfile << this->interpolated_vm_BaseCellMembranePotential(s) << " ";

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
			this->get_diff_monodomain(iplot,s,x,D);
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

	};

	//Inline functions:


	//======================================================================
	/// \short Define the shape functions and test functions and derivatives
	/// w.r.t. global coordinates and return Jacobian of mapping.
	///
	/// Galerkin: Test functions = shape functions
	//======================================================================
	template<unsigned DIM, unsigned NNODE_1D>
	double QMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>::
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
	double QMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>::
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


	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////



	//=======================================================================
	/// \short Face geometry for the QMonodomainElementIntegralSplittingMethod elements: 
	/// The spatial dimension of the face elements is one lower than that 
	/// of the bulk element but they have the same number of points along 
	/// their 1D edges.
	//=======================================================================
	template<unsigned DIM, unsigned NNODE_1D>
	class FaceGeometry<QMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D> >: 
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
	class FaceGeometry<QMonodomainElementIntegralSplittingMethod<1,NNODE_1D> >: 
	 public virtual PointElement
	{

	  public:
	 
	 /// \short Constructor: Call the constructor for the
	 /// appropriate lower-dimensional QElement
	 FaceGeometry() : PointElement() {}

	};



	//Override functions in specific implementations of the diff augmented wrapper
	template<unsigned DIM, unsigned NNODE_1D>
	class DiffAugmentedCell<QMonodomainElementIntegralSplittingMethod<DIM, NNODE_1D>>:
		public QMonodomainElementIntegralSplittingMethod<DIM, NNODE_1D>
	{
	public:

		inline void get_diff_monodomain(const unsigned& ipt,
                                        const Vector<double> &s,
                                        const Vector<double>& x,
                                        DenseMatrix<double>& D) const
		{
			this->get_interpolated_diffusion_matrix(s, D);
		}

	};





	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	// TMonodomainElementIntegralSplittingMethod
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////



	//======================================================================
	/// \short TMonodomainElementIntegralSplittingMethod elements are isoparametric triangular 
	/// DIM-dimensional General Advection Diffusion Equations with  NNODE_1D nodal points along each
	/// element edge. Inherits from TElement and MonodomainEquationsIntegralSplittingMethod
	//======================================================================
	template <unsigned DIM, unsigned NNODE_1D>
	 class TMonodomainElementIntegralSplittingMethod : 
	 public virtual TElement<DIM,NNODE_1D>,
	 public virtual MonodomainEquationsIntegralSplittingMethod<DIM>
	{

	private:

	 /// \short Static array of ints to hold number of variables at 
	 /// nodes: Initial_Nvalue[n]
	 static const unsigned Initial_Nvalue = 1;
	 
	  public:


	 ///\short  Constructor: Call constructors for TElement and 
	 /// Advection Diffusion equations
	 TMonodomainElementIntegralSplittingMethod() : TElement<DIM,NNODE_1D>(), 
	  MonodomainEquationsIntegralSplittingMethod<DIM>()
	  { }

	 /// Broken copy constructor
	 TMonodomainElementIntegralSplittingMethod(
	  const TMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>&  dummy) 
	  { 
	   BrokenCopy::broken_copy("TMonodomainElementIntegralSplittingMethod");
	  } 
	 
	 /// Broken assignment operator
	 void operator=(const TMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>&) 
	  {
	   BrokenCopy::broken_assign("TMonodomainElementIntegralSplittingMethod");
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

	};

	//Inline functions:


	//======================================================================
	/// \short Define the shape functions and test functions and derivatives
	/// w.r.t. global coordinates and return Jacobian of mapping.
	///
	/// Galerkin: Test functions = shape functions
	//======================================================================
	template<unsigned DIM, unsigned NNODE_1D>
	double TMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>::
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
	double TMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D>::
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


	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////



	//=======================================================================
	/// \short Face geometry for the TMonodomainElementIntegralSplittingMethod elements: 
	/// The spatial dimension of the face elements is one lower than that 
	/// of the bulk element but they have the same number of points along 
	/// their 1D edges.
	//=======================================================================
	template<unsigned DIM, unsigned NNODE_1D>
	class FaceGeometry<TMonodomainElementIntegralSplittingMethod<DIM,NNODE_1D> >: 
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
	/// Face geometry for the 1D TMonodomainElementIntegralSplittingMethod: Point elements
	//=======================================================================
	template<unsigned NNODE_1D>
	class FaceGeometry<TMonodomainElementIntegralSplittingMethod<1,NNODE_1D> >: 
	 public virtual PointElement
	{

	  public:
	 
	 /// \short Constructor: Call the constructor for the
	 /// appropriate lower-dimensional TElement
	 FaceGeometry() : PointElement() {}

	};


	//Override functions in specific implementations of the diff augmented wrapper
	template<unsigned DIM, unsigned NNODE_1D>
	class DiffAugmentedCell<TMonodomainElementIntegralSplittingMethod<DIM, NNODE_1D>>:
		public TMonodomainElementIntegralSplittingMethod<DIM, NNODE_1D>
	{
	public:

		inline void get_diff_monodomain(const unsigned& ipt,
                                        const Vector<double> &s,
                                        const Vector<double>& x,
                                        DenseMatrix<double>& D) const
		{
			this->get_interpolated_diffusion_matrix(s, D);
		}

	};


}


#endif