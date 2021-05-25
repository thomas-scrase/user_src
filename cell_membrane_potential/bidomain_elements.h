//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations
//LIC//====================================================================

#ifndef OOMPH_BIDOMAIN_ELEMENTS_HEADER
#define OOMPH_BIDOMAIN_ELEMENTS_HEADER

#include "cell_membrane_potential_elements.h"

namespace oomph{

//Bidomain Equations
	template <unsigned DIM>
	class BidomainEquations :
	public virtual BaseCellMembranePotentialEquations<DIM>
	{
	public:

		//change this to take s instead of x? (no functional change, just notation)
		/// \short Funciton pointer to a diffusivity function
		typedef void (*BidomainEquationsConductFctPt)
		(const Vector<double> &x, DenseMatrix<double> &G);

		/// \short Function pointer to source function fct(x,f(x)) -- 
		/// x is a Vector! 
		typedef void (*BidomainEquationsSourceFctPt)
		(const Vector<double>& x, double& f);

        BidomainEquations()	:	Intracellular_Conductivity_Fct_Pt(0),
        						Extracellular_Conductivity_Fct_Pt(0),
        						Extracellular_Source_fct_pt(0)
        {

        }

        //======================================================================
		/// \short Compute element residual Vector and/or element Jacobian matrix 
		/// 
		/// flag=1: compute both
		/// flag=0: compute only residual Vector
		///
		/// Pure version without hanging nodes
		//======================================================================
		void fill_in_generic_residual_contribution_BaseCellMembranePotential(Vector<double> &residuals, 
																			DenseMatrix<double> &jacobian, 
																			DenseMatrix<double> &mass_matrix,
																			unsigned flag) 
		{
		//Find out how many nodes there are
		const unsigned n_node = this->nnode();

		//Get the nodal index at which the unknown is stored
		const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();
		const unsigned phie_nodal_index = this->phie_index_Bidomain();
		 
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
		  // std::cout << "integration point " << ipt << std::endl;

		 //Assign values of s
		 for(unsigned i=0;i<DIM;i++) s[i] = this->integral_pt()->knot(ipt,i);

		 //Get the integral weight
		 double w = this->integral_pt()->weight(ipt);

		  if(w==0.0){continue;}
		  
		 //Call the derivatives of the shape and test functions
		 double J = 
		  this->dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(ipt,psi,dpsidx,test,dtestdx);
		     
		 //Premultiply the weights and the Jacobian
		 double W = w*J;

		 //Calculate local values of the solution and its derivatives
		 //Allocate
		 double interpolated_vm=0.0;
		 double interpolated_phie=0.0;
		 double dvmdt=0.0;
		 double dphiedt=0.0;

		 Vector<double> interpolated_x(DIM,0.0);
		 Vector<double> interpolated_dvmdx(DIM,0.0);
		 Vector<double> interpolated_phiedx(DIM,0.0);
		 Vector<double> mesh_velocity(DIM,0.0);

		 //Calculate function value and derivatives:
		 //-----------------------------------------
		 // Loop over nodes
		 for(unsigned l=0;l<n_node;l++) 
		  {
		   //Get the value at the node
		   double vm_value = this->raw_nodal_value(l,vm_nodal_index);
		   double phie_value = this->raw_nodal_value(l, phie_nodal_index);
		   interpolated_vm += vm_value*psi(l);
		   interpolated_phie += phie_value*psi(l);
		   dvmdt += this->dvm_dt_BaseCellMembranePotential(l)*psi(l);
		   dphiedt += this->dphie_dt_Bidomain(l)*psi(l);
		   // Loop over directions
		   for(unsigned j=0;j<DIM;j++)
		    {
		     interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
		     interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
		     interpolated_phiedx[j] += phie_value*dpsidx(l,j);
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

		 //Get source functions
		 //-------------------
		 double intracellular_source;
		 this->get_source_BaseCellMembranePotential(ipt,s,interpolated_x,intracellular_source);

		 double extracellular_source;
		 this->get_extracellular_source_BidomainEquations(ipt,s,interpolated_x,extracellular_source);

		 //Get conductivity tensors
		 DenseMatrix<double> Gi(DIM,DIM,0.0);
		 this->get_intracellular_conductivity_bidomain(ipt,s,interpolated_x,Gi);

		 DenseMatrix<double> Ge(DIM,DIM,0.0);
		 this->get_extracellular_conductivity_bidomain(ipt,s,interpolated_x,Ge);

		 // Assemble residuals and Jacobian
		 //--------------------------------
		     
		 // Loop over the test functions
		 for(unsigned l=0;l<n_node;l++)
		 { 


		    //Residual and jacobian for transmembrane potential dof


		    //Fill in residual and jacobian contribution for u
		    //Set the local equation number
		    local_eqn = this->nodal_local_eqn(l,vm_nodal_index);

		    /*IF it's not a boundary condition*/
		    if(local_eqn >= 0)
		    {
		      // Add body force/source term and time derivative 
		      residuals[local_eqn] -= (dvmdt + intracellular_source)*test(l)*W;
		     
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
		         tmp2 += (interpolated_dvmdx[j] + interpolated_phiedx[j])*Gi(k,j);
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
		          local_unknown = this->nodal_local_eqn(l2,vm_nodal_index);
		         
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
		                tmp2 += Gi(k,j)*dpsidx(l2,j);
		                }
		             
		              //Now assemble Jacobian term
		              jacobian(local_eqn,local_unknown) 
		                -= (tmp*test(l) + tmp2*dtestdx(l,k))*W;
		              }
		            }

		            local_unknown = this->nodal_local_eqn(l2, this->phie_index_Bidomain());
		            //If at a non-zero degree of freedom add in the entry
		            if(local_unknown >= 0)
		            {

		              //Add the mass matrix term
		              if(flag==2)
		              {
		                mass_matrix(local_eqn,local_unknown)
		                  += test(l)*psi(l2)*W;
		              }

		              //Add contribution to Elemental Matrix
		              for(unsigned k=0;k<DIM;k++){
		                double tmp2 = 0.0;
		                //Now the diffusive term
		                for(unsigned j=0;j<DIM;j++)
		                {
		                  tmp2 += Gi(k,j)*dpsidx(l2,j);
		                }
		               
		                //Now assemble Jacobian term
		                jacobian(local_eqn,local_unknown) 
		                  -= (tmp2*dtestdx(l,k))*W;
		              }
		            
		            }

		          }
		        }
		      }


		      //Residual and jacobian for extracellular potential dof


		      //Fill in residual and jacobian contribution for u
		      //Set the local equation number
		      local_eqn = this->nodal_local_eqn(l,this->phie_index_Bidomain());

		      /*IF it's not a boundary condition*/
		      if(local_eqn >= 0)
		      {
		        residuals[local_eqn] -= (extracellular_source)*test(l)*W;

		        // The Generalised Advection Diffusion bit itself
		        for(unsigned k=0;k<DIM;k++)
		          {
		           //Terms that multiply the derivative of the test function
		            double tmp2 = 0.0;
		           //Now the diuffusive term
		           for(unsigned j=0;j<DIM;j++)
		            {
		             tmp2 += (interpolated_dvmdx[j] + interpolated_phiedx[j])*Gi(k,j) +
		                      interpolated_phiedx[j]*Ge(k,j);
		            }
		           //Now construct the contribution to the residuals
		           residuals[local_eqn] -= (tmp2*dtestdx(l,k))*W;
		          }
		       
		        // Calculate the jacobian
		        //-----------------------
		        if(flag)
		        {
		          //Loop over the velocity shape functions again
		          for(unsigned l2=0;l2<n_node;l2++)
		          { 
		            //Set the number of the unknown
		            local_unknown = this->nodal_local_eqn(l2,vm_nodal_index);
		           
		            //If at a non-zero degree of freedom add in the entry
		            if(local_unknown >= 0)
		            {
		              //Add the mass matrix term
		              if(flag==2)
		              {
		              mass_matrix(local_eqn,local_unknown)
		                += test(l)*psi(l2)*W;
		              }

		              //Add contribution to Elemental Matrix
		              for(unsigned k=0;k<DIM;k++){
		                double tmp2 = 0.0;
		                //Now the diffusive term
		                for(unsigned j=0;j<DIM;j++)
		                {
		                  tmp2 += Gi(k,j)*dpsidx(l2,j);
		                }
		               
		                //Now assemble Jacobian term
		                jacobian(local_eqn,local_unknown) 
		                  -= (tmp2*dtestdx(l,k))*W;
		              }
		            }

		            local_unknown = this->nodal_local_eqn(l2, this->phie_index_Bidomain());
		            //If at a non-zero degree of freedom add in the entry
		            if(local_unknown >= 0)
		            {
		              //Add the mass matrix term
		              if(flag==2)
		                {
		                mass_matrix(local_eqn,local_unknown)
		                  += test(l)*psi(l2)*W;
		                }

		              //Add contribution to Elemental Matrix
		              for(unsigned k=0;k<DIM;k++){

		                double tmp2 = 0.0;
		                //Now the diffusive term
		                for(unsigned j=0;j<DIM;j++)
		                  {
		                  tmp2 += (Gi(k,j) + Ge(k,j))*dpsidx(l2,j);
		                  }
		               
		                //Now assemble Jacobian term
		                jacobian(local_eqn,local_unknown)
		                  -= (tmp2*dtestdx(l,k))*W;
		              }
		            }
		          }

		        }
		    }

		    
		  }
		} // End of loop over integration points
		}

		inline unsigned required_nvalue(const unsigned &n) const {return 2;}

        virtual inline unsigned phie_index_Bidomain() const {return this->vm_index_BaseCellMembranePotential()+1;}

        virtual inline unsigned max_index_plus_one_BaseCellMembranePotential() const {return phie_index_Bidomain()+1;}

        double dphie_dt_Bidomain(const unsigned &n) const
        { 	
			// Get the data's timestepper
			TimeStepper* time_stepper_pt= this->node_pt(n)->time_stepper_pt();

			//Initialise dudt
			double dphiedt=0.0;
			//Loop over the timesteps, if there is a non Steady timestepper
			if (!time_stepper_pt->is_steady())
			{
				//Find the index at which the variable is stored
				const unsigned phie_nodal_index = phie_index_Bidomain();

				// Number of timsteps (past & present)
				const unsigned n_time = time_stepper_pt->ntstorage();

				for(unsigned t=0;t<n_time;t++)
				{
					dphiedt += time_stepper_pt->weight(1,t)*this->nodal_value(t,n,phie_nodal_index);
				}
			}
			return dphiedt;
        }

		/// Access function: Pointer to diffusion  function
		BidomainEquationsConductFctPt& intracellular_conductivity_fct_pt() 
		{return Intracellular_Conductivity_Fct_Pt;}

		/// Access function: Pointer to diffusion function. Const version
		BidomainEquationsConductFctPt intracellular_conductivity_fct_pt() const 
		{return Intracellular_Conductivity_Fct_Pt;}

		/// Access function: Pointer to diffusion  function
		BidomainEquationsConductFctPt& extracellular_conductivity_fct_pt() 
		{return Extracellular_Conductivity_Fct_Pt;}

		/// Access function: Pointer to diffusion function. Const version
		BidomainEquationsConductFctPt extracellular_conductivity_fct_pt() const 
		{return Extracellular_Conductivity_Fct_Pt;}


		void assign_additional_initial_conditions()
		{
			for(unsigned n=0; n<this->nnode(); n++){
				//Set the initial value of the extracellular potential to zero, then intracellular potential is implicitly equal to transmembrane potential
				this->node_pt(n)->set_value(phie_index_Bidomain(), 0.0);
			}
		}


		/// Access function: Pointer to source function
		BidomainEquationsSourceFctPt& extracellular_source_fct_pt() 
		{return Extracellular_Source_fct_pt;}

		/// Access function: Pointer to source function. Const version
		BidomainEquationsSourceFctPt extracellular_source_fct_pt() const 
		{return Extracellular_Source_fct_pt;}

		/// \short Get source term at (Eulerian) position x. This function is
		/// virtual to allow overloading in multi-physics problems where
		/// the strength of the source function might be determined by
		/// another system of equations 
		inline virtual void get_extracellular_source_BidomainEquations(const unsigned& ipt,
																		const Vector<double>& s,
																		const Vector<double>& x,
																		double& source) const
		{
			//If no source function has been set, return zero
			if(Extracellular_Source_fct_pt==0)
			{
				source = 0.0;
			}
			else
			{
				// Get source strength
				(*Extracellular_Source_fct_pt)(x,source);
			}
		}

		/// \short Get diffusivity tensor at (Eulerian) position 
		/// x and/or local coordinate s. 
		/// This function is
		/// virtual to allow overloading in multi-physics problems where
		/// the diff function might be determined by
		/// another system of equations 
		inline virtual void get_intracellular_conductivity_bidomain(const unsigned& ipt,
		                                            const Vector<double> &s,
		                                            const Vector<double>& x,
		                                            DenseMatrix<double>& G_i) const
		{
			//If no diff function has been set, return identity
			if(Intracellular_Conductivity_Fct_Pt==0){
				for(unsigned i=0; i<DIM; i++){
					for(unsigned j=0; j<DIM; j++){
						G_i(i,j) =  0.0;
					}
					G_i(i,i)  = 1.0;
				}
			}
			else{
				// oomph_info << "Using Intracellular conductivity fct pt" << std::endl;
				// Get diffusivity tensor from function
				(*Intracellular_Conductivity_Fct_Pt)(x,G_i);
			}
		}

		inline virtual void get_extracellular_conductivity_bidomain(const unsigned& ipt,
		                                            const Vector<double> &s,
		                                            const Vector<double>& x,
		                                            DenseMatrix<double>& G_e) const
		{
			//If no diff function has been set, return identity
			if(Extracellular_Conductivity_Fct_Pt==0){
				for(unsigned i=0; i<DIM; i++){
					for(unsigned j=0; j<DIM; j++){
						G_e(i,j) =  0.0;
					}
					G_e(i,i)  = 1.0;
				}
			}
			else{
				// oomph_info << "Using Extracellular conductivity fct pt" << std::endl;
				// Get diffusivity tensor from function
				(*Extracellular_Conductivity_Fct_Pt)(x,G_e);
			}
		}


		/// Return FE representation of function value u(s) at local coordinate s
		inline double interpolated_phie_BaseCellMembranePotential(const Vector<double> &s) const
		{
			//Find number of nodes
			unsigned n_node = this->nnode();

			//Get the nodal index at which the unknown is stored
			unsigned phie_nodal_index = this->phie_index_Bidomain();

			//Local shape function
			Shape psi(n_node);

			//Find values of shape function
			this->shape(s,psi);

			//Initialise value of u
			double interpolated_phie = 0.0;

			//Loop over the local nodes and sum
			for(unsigned l=0;l<n_node;l++) 
			{
			 interpolated_phie += this->nodal_value(l,phie_nodal_index)*psi[l];
			}

			return interpolated_phie;
		}



		//!!!!!this is mostly used for error checking. What form should it take,
		//	perhaps flux of the transmembrane + intracellular potential?
		// Maybe get G_i.grad(V+phi_e) + G_e.grad(phi_e) ?
		/// Get flux: \f$\mbox{flux}[i] = \mbox{d}u / \mbox{d}x_i \f$
		void get_total_flux_bidomain(const Vector<double>& s,
		                			Vector<double>& total_flux) const
		{
			//Find out how many nodes there are in the element
			const unsigned n_node = this->nnode();

			//Get the nodal index at which the unknown is stored
			const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();
			const unsigned phie_nodal_index = this->phie_index_Bidomain();

			//Set up memory for the shape and test functions
			Shape psi(n_node);
			DShape dpsidx(n_node,DIM);

			//Call the derivatives of the shape and test functions
			this->dshape_eulerian(s,psi,dpsidx);

			//Storage for the Eulerian position
			Vector<double> interpolated_x(DIM,0.0);
			//Storage for the derivatives of the concentration
			Vector<double> interpolated_dvmdx(DIM,0.0);
			Vector<double> interpolated_dphiedx(DIM,0.0);
			// Loop over nodes
			for(unsigned l=0;l<n_node;l++) 
			{
			 //Get the value at the node
			 const double vm_value = this->nodal_value(l,vm_nodal_index);
			 const double phie_value = this->nodal_value(l,phie_nodal_index);
			 // Loop over directions
			 for(unsigned j=0;j<DIM;j++)
			  {
			   interpolated_x[j] += this->nodal_position(l,j)*psi(l);
			   interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
			   interpolated_dphiedx[j] += phie_value*dpsidx(l,j);
			  }
			}

			//Dummy integration point 
			unsigned ipt=0;

			//Get conductivity tensors
			DenseMatrix<double> Gi(DIM,DIM,0.0);
			this->get_intracellular_conductivity_bidomain(ipt,s,interpolated_x,Gi);

			DenseMatrix<double> Ge(DIM,DIM,0.0);
			this->get_extracellular_conductivity_bidomain(ipt,s,interpolated_x,Ge);

			//Calculate the total flux made up of the diffusive flux
			//and the conserved wind
			for(unsigned i=0;i<DIM;i++)
			{                               
			 total_flux[i] = 0.0;
			 for(unsigned j=0;j<DIM;j++)
			  {
			   total_flux[i] += (interpolated_dvmdx[j] + interpolated_dphiedx[j])*Gi(i,j) + interpolated_dphiedx[j]*Ge(i,j);
			  }
			}
		}

		inline std::vector<std::string> get_variable_names() const override
		{
			std::vector<std::string> v = {"Vm", "phie"};
			return v;
		}


		unsigned nscalar_paraview() const
		{
			return 2;
		}

		void scalar_value_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot) const
		{
			//Vector of local coordinates
			Vector<double> s(DIM);

			const unsigned n_node = this->nnode();
			// const unsigned vm_index = this->vm_index_BaseCellMembranePotential();
			Shape psi(n_node);
			DShape dpsidx(n_node,DIM);

			// Loop over plot points
			unsigned num_plot_points=this->nplot_points_paraview(nplot);
			for (unsigned iplot=0;iplot<num_plot_points;iplot++){

				switch(i){
					case 0 :
						// Get local coordinates of plot point
						this->get_s_plot(iplot,nplot,s);

						file_out << this->interpolated_vm_BaseCellMembranePotential(s) << " ";
					case 1 :
						// Get local coordinates of plot point
						this->get_s_plot(iplot,nplot,s);

						file_out << this->interpolated_phie_BaseCellMembranePotential(s) << " ";

				}
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

	protected:
		/// Pointer to diffusivity function:
		///		function typedef is given in BaseCellMembranePotentialEquations
 		BidomainEquationsConductFctPt Intracellular_Conductivity_Fct_Pt;
 		BidomainEquationsConductFctPt Extracellular_Conductivity_Fct_Pt;

 		BidomainEquationsSourceFctPt Extracellular_Source_fct_pt;

	};




//Bidomain Elements

	//======================================================================
	/// \short QBidomainElement elements are 
	/// linear/quadrilateral/brick-shaped Advection Diffusion elements with 
	/// isoparametric interpolation for the function.
	//======================================================================
	template <unsigned DIM, unsigned NNODE_1D>
	 class QBidomainElement : 
	 public virtual QElement<DIM,NNODE_1D>,
	 public virtual BidomainEquations<DIM>
	{

	private:

	 /// \short Static array of ints to hold number of variables at 
	 /// nodes: Initial_Nvalue[n]
	 static const unsigned Initial_Nvalue = 2;
	 
	  public:


	 ///\short  Constructor: Call constructors for QElement and 
	 /// Advection Diffusion equations
	 QBidomainElement() : QElement<DIM,NNODE_1D>(), 
	  BidomainEquations<DIM>()
	  { }

	 /// Broken copy constructor
	 QBidomainElement(
	  const QBidomainElement<DIM,NNODE_1D>&  dummy) 
	  { 
	   BrokenCopy::broken_copy("QBidomainElement");
	  } 
	 
	 /// Broken assignment operator
	 void operator=(const QBidomainElement<DIM,NNODE_1D>&) 
	  {
	   BrokenCopy::broken_assign("QBidomainElement");
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

			// //Get diffusivity tensor
			// DenseMatrix<double> D(DIM,DIM,0.0);
			// this->get_diff_monodomain(iplot,s,x,D);
			// for(unsigned i=0; i<DIM; i++){
			// 	for(unsigned j=0; j<=i; j++){
			// 		outfile << D(i,j) << " ";
			// 	}
			// }

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
	double QBidomainElement<DIM,NNODE_1D>::
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
	double QBidomainElement<DIM,NNODE_1D>::
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
	/// \short Face geometry for the QBidomainElement elements: 
	/// The spatial dimension of the face elements is one lower than that 
	/// of the bulk element but they have the same number of points along 
	/// their 1D edges.
	//=======================================================================
	template<unsigned DIM, unsigned NNODE_1D>
	class FaceGeometry<QBidomainElement<DIM,NNODE_1D> >: 
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
	/// Face geometry for the 1D QBidomain elements: Point elements
	//=======================================================================
	template<unsigned NNODE_1D>
	class FaceGeometry<QBidomainElement<1,NNODE_1D> >: 
	 public virtual PointElement
	{

	  public:
	 
	 /// \short Constructor: Call the constructor for the
	 /// appropriate lower-dimensional QElement
	 FaceGeometry() : PointElement() {}

	};





	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	// TBidomainElement
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////



	//======================================================================
	/// \short TBidomainElement elements are isoparametric triangular 
	/// DIM-dimensional General Advection Diffusion Equations with  NNODE_1D nodal points along each
	/// element edge. Inherits from TElement and BidomainEquations
	//======================================================================
	template <unsigned DIM, unsigned NNODE_1D>
	 class TBidomainElement : 
	 public virtual TElement<DIM,NNODE_1D>,
	 public virtual BidomainEquations<DIM>
	{

	private:

	 /// \short Static array of ints to hold number of variables at 
	 /// nodes: Initial_Nvalue[n]
	 static const unsigned Initial_Nvalue = 2;
	 
	  public:


	 ///\short  Constructor: Call constructors for TElement and 
	 /// Advection Diffusion equations
	 TBidomainElement() : TElement<DIM,NNODE_1D>(), 
	  BidomainEquations<DIM>()
	  { }

	 /// Broken copy constructor
	 TBidomainElement(
	  const TBidomainElement<DIM,NNODE_1D>&  dummy) 
	  { 
	   BrokenCopy::broken_copy("TBidomainElement");
	  } 
	 
	 /// Broken assignment operator
	 void operator=(const TBidomainElement<DIM,NNODE_1D>&) 
	  {
	   BrokenCopy::broken_assign("TBidomainElement");
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
	double TBidomainElement<DIM,NNODE_1D>::
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
	double TBidomainElement<DIM,NNODE_1D>::
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
	/// \short Face geometry for the TBidomainElement elements: 
	/// The spatial dimension of the face elements is one lower than that 
	/// of the bulk element but they have the same number of points along 
	/// their 1D edges.
	//=======================================================================
	template<unsigned DIM, unsigned NNODE_1D>
	class FaceGeometry<TBidomainElement<DIM,NNODE_1D> >: 
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
	/// Face geometry for the 1D TBidomainElement: Point elements
	//=======================================================================
	template<unsigned NNODE_1D>
	class FaceGeometry<TBidomainElement<1,NNODE_1D> >: 
	 public virtual PointElement
	{

	  public:
	 
	 /// \short Constructor: Call the constructor for the
	 /// appropriate lower-dimensional TElement
	 FaceGeometry() : PointElement() {}

	};
}


#endif