//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations - use strang splitting method
//LIC//====================================================================

//IMPLEMENTS NOTHING. USED TO BUILD CELL MESHES WHERE WE DO NOT WANT ANY DIFFUSION. USED FOR TESTING
// AND FOR IMPLEMENTING NON-REFINEABLE CELL MESH EXTERNAL TO A REFINEABLE DIFFUSION MESH

#ifndef OOMPH_NO_DIFFUSION_ELEMENTS_HEADER
#define OOMPH_NO_DIFFUSION_ELEMENTS_HEADER

#include "cell_membrane_potential_elements.h"

namespace oomph{

//Monodomain Equations
	template <unsigned DIM>
	class NoDiffusionEquations :
	public BaseCellMembranePotentialEquations<DIM>
	{
	public:

        NoDiffusionEquations() : ipt_not_at_nodes(0)
        {
        }

        //Do nothing, we have no locally stored membrane potential to update
        virtual inline void update_nodal_membrane_potential_BaseCellMembranePotential(const unsigned &l, const double& vm)
		{

		}

		//The membrane potential is just stored at the nodes
		virtual inline double get_nodal_membrane_potential_BaseCellMembranePotential(const unsigned &n) const
		{

			throw OomphLibError(
				"Not using overriden version",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
			return 0.0;
		}

		//Get the t-th history value of membrane potential at the nth node
		virtual inline double get_nodal_membrane_potential_BaseCellMembranePotential(const unsigned &t, const unsigned &n) const
		{
			throw OomphLibError(
				"Not using overriden version",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
			return 0.0;
		}
		

		//There are no equations so there are no residuals
		void fill_in_generic_residual_contribution_BaseCellMembranePotential(
	    Vector<double> &residuals, DenseMatrix<double> &jacobian, 
	    DenseMatrix<double> &mass_matrix, unsigned flag)
		{
		}
		
		//The element has no degrees of freedom
		inline unsigned required_nvalue(const unsigned &n) const {return 0;}


		void assign_additional_initial_conditions(const unsigned &l)
		{
			//Do nothing - we have no additional variables to assign values to
		}

		void assign_initial_conditions_consistent_with_cell_model(const unsigned &l, const double& vm)
		{
		}


		//Calculate dvmdt but instead of using the stored value for the previous timestep, use the value predicted by the cell model
		double dvm_dt_Strang_Split(const unsigned &n) const
		{
			throw OomphLibError(
				"Not using overriden version",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
			return 0.0;
		}


		/// Get flux: \f$\mbox{flux}[i] = \mbox{d}u / \mbox{d}x_i \f$
		void get_total_flux_monodomain(const Vector<double>& s, 
		                			Vector<double>& total_flux) const
		{
		}

		inline std::vector<std::string> get_variable_names_BaseCellMembranePotentialEquations() const override
		{
			std::vector<std::string> v = {""};
			return v;
		}


		unsigned nscalar_paraview() const
		{
			return 0;
		}

		void scalar_value_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot) const
		{
			throw OomphLibError(
				"I have nothing to output",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
			
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
			return "";
		}


		/////////////////////////////////////////////////////////////////////////////////
		//Residual and Jacobian functions
		// The conductance model is the only thing that does actual oomph lib solving
		/////////////////////////////////////////////////////////////////////////////////
		/// Add the element's contribution to its residual vector (wrapper)
		void fill_in_contribution_to_residuals(Vector<double> &residuals)
		{
		}

		/// \short Add the element's contribution to its residual vector and 
		/// the element Jacobian matrix (wrapper)
		void fill_in_contribution_to_jacobian(Vector<double> &residuals,
		DenseMatrix<double> &jacobian)
		{
		}

		/// Add the element's contribution to its residuals vector,
		/// jacobian matrix and mass matrix
		void fill_in_contribution_to_jacobian_and_mass_matrix(
		Vector<double> &residuals, DenseMatrix<double> &jacobian,
		DenseMatrix<double> &mass_matrix)
		{
		}

	protected:
		//The number of integral points which are not additional ones placed at the nodes
		unsigned ipt_not_at_nodes;

};


//WE DON'T NEED Q OR T ELEMENTS BECAUSE THESE EQUATIONS ARE ONLY BUILD WITHIN THE CONTEXT OF A CELL MESH, WHICH IS HANDLED VIA A Q OR T WRAPPER AUTOMATICALLY


// //Monodomain Elements

// 	//======================================================================
// 	/// \short QNoDiffusionElement elements are 
// 	/// linear/quadrilateral/brick-shaped Advection Diffusion elements with 
// 	/// isoparametric interpolation for the function.
// 	//======================================================================
// 	template <unsigned DIM, unsigned NNODE_1D>
// 	 class QNoDiffusionElement : 
// 	 public virtual QElement<DIM,NNODE_1D>,
// 	 public virtual NoDiffusionEquations<DIM>
// 	{

// 	private:

// 	 /// \short Static array of ints to hold number of variables at 
// 	 /// nodes: Initial_Nvalue[n]
// 	 static const unsigned Initial_Nvalue = 1;
	 
// 	  public:


// 	 ///\short  Constructor: Call constructors for QElement and 
// 	 /// Advection Diffusion equations
// 	 QNoDiffusionElement() : QElement<DIM,NNODE_1D>(), 
// 	  NoDiffusionEquations<DIM>()
// 	  { }

// 	 /// Broken copy constructor
// 	 QNoDiffusionElement(
// 	  const QNoDiffusionElement<DIM,NNODE_1D>&  dummy) 
// 	  { 
// 	   BrokenCopy::broken_copy("QNoDiffusionElement");
// 	  } 
	 
// 	 /// Broken assignment operator
// 	 void operator=(const QNoDiffusionElement<DIM,NNODE_1D>&) 
// 	  {
// 	   BrokenCopy::broken_assign("QNoDiffusionElement");
// 	  }

// 	 /// \short  Required  # of `values' (pinned or dofs) 
// 	 /// at node n
// 	 inline unsigned required_nvalue(const unsigned &n) const 
// 	  {return Initial_Nvalue;}

// 	 /// \short Output function:  
// 	 ///  x,y,u   or    x,y,z,u
// 	 void output(std::ostream &outfile)
// 	  {BaseCellMembranePotentialEquations<DIM>::output(outfile);}

// 	 /// \short Output function:  
// 	 ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
// 	 void output(std::ostream &outfile, const unsigned &nplot) override
// 	 {
// 	 	// std::cout << "BOOM\nBOOM\nBOOM\nBOOM\nBOOM" << std::endl;
// 		//Vector of local coordinates
// 		Vector<double> s(DIM);

// 		// Tecplot header info
// 		outfile << this->tecplot_zone_string(nplot);

// 		const unsigned n_node = this->nnode();
// 		const unsigned vm_index = this->vm_index_BaseCellMembranePotential();
// 		// std::cout << "n_node " << n_node << std::endl;
// 		Shape psi(n_node);
// 		DShape dpsidx(n_node,DIM);

// 		// Loop over plot points
// 		unsigned num_plot_points=this->nplot_points(nplot);
// 		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
// 			// Get local coordinates of plot point
// 			this->get_s_plot(iplot,nplot,s);

// 			// Get Eulerian coordinate of plot point
// 			Vector<double> x(DIM);
// 			this->interpolated_x(s,x);

// 			for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}
// 			outfile << this->interpolated_vm_BaseCellMembranePotential(s) << " ";

// 			//Get the gradients
// 			(void)this->dshape_eulerian(s,psi,dpsidx);
// 			Vector<double> interpolated_dvmdx(DIM,0.0);
// 			double dvmdt = 0.0;
// 			for(unsigned n=0;n<n_node;n++){
// 				const double vm_ = this->nodal_value(n,vm_index);
// 				dvmdt += this->dvm_dt_BaseCellMembranePotential(n)*psi(n);
// 				for(unsigned i=0;i<DIM;i++){interpolated_dvmdx[i] += vm_*dpsidx(n,i);}
// 			}

// 			outfile << dvmdt << " ";

// 			for(unsigned i=0;i<DIM;i++){outfile << interpolated_dvmdx[i]  << " ";}

// 			//Get diffusivity tensor
// 			DenseMatrix<double> D(DIM,DIM,0.0);
// 			this->get_diff_monodomain(iplot,s,x,D);
// 			for(unsigned i=0; i<DIM; i++){
// 				for(unsigned j=0; j<=i; j++){
// 					outfile << D(i,j) << " ";
// 				}
// 			}

// 			outfile  << std::endl;
// 		}
// 		// Write tecplot footer (e.g. FE connectivity lists)
// 		this->write_tecplot_zone_footer(outfile,nplot);
// 	}


// 	 /// \short C-style output function:  
// 	 ///  x,y,u   or    x,y,z,u
// 	 void output(FILE* file_pt)
// 	  {
// 	   BaseCellMembranePotentialEquations<DIM>::output(file_pt);
// 	  }

// 	 ///  \short C-style output function:  
// 	 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
// 	 void output(FILE* file_pt, const unsigned &n_plot)
// 	  {
// 	   BaseCellMembranePotentialEquations<DIM>::output(file_pt,n_plot);
// 	  }

// 	 /// \short Output function for an exact solution:
// 	 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
// 	 void output_fct(std::ostream &outfile, const unsigned &n_plot,
// 	                 FiniteElement::SteadyExactSolutionFctPt 
// 	                 exact_soln_pt)
// 	  {BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,exact_soln_pt);}


// 	 /// \short Output function for a time-dependent exact solution.
// 	 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
// 	 /// (Calls the steady version)
// 	 void output_fct(std::ostream &outfile, const unsigned &n_plot,
// 	                 const double& time,
// 	                 FiniteElement::UnsteadyExactSolutionFctPt 
// 	                 exact_soln_pt)
// 	  {
// 	   BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,time,exact_soln_pt);
// 	  }


// 	protected:

// 	 	inline double dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s, 
// 	                                               Shape &psi, 
// 	                                               DShape &dpsidx, 
// 	                                               Shape &test, 
// 	                                               DShape &dtestdx) const;

// 	 	inline double dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt,
// 													Shape &psi, 
// 													DShape &dpsidx,
// 													Shape &test, 
// 													DShape &dtestdx) const;

// 	};

// 	//Inline functions:


// 	//======================================================================
// 	/// \short Define the shape functions and test functions and derivatives
// 	/// w.r.t. global coordinates and return Jacobian of mapping.
// 	///
// 	/// Galerkin: Test functions = shape functions
// 	//======================================================================
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	double QNoDiffusionElement<DIM,NNODE_1D>::
// 	 dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s,
// 	                                         Shape &psi, 
// 	                                         DShape &dpsidx,
// 	                                         Shape &test, 
// 	                                         DShape &dtestdx) const
// 	{
// 	 //Call the geometrical shape functions and derivatives  
// 	 double J = this->dshape_eulerian(s,psi,dpsidx);

// 	 //Loop over the test functions and derivatives and set them equal to the
// 	 //shape functions
// 	 for(unsigned i=0;i<NNODE_1D;i++)
// 	  {
// 	   test[i] = psi[i]; 
// 	   for(unsigned j=0;j<DIM;j++)
// 	    {
// 	     dtestdx(i,j) = dpsidx(i,j);
// 	    }
// 	  }
	 
// 	 //Return the jacobian
// 	 return J;
// 	}



// 	//======================================================================
// 	/// Define the shape functions and test functions and derivatives
// 	/// w.r.t. global coordinates and return Jacobian of mapping.
// 	///
// 	/// Galerkin: Test functions = shape functions
// 	//======================================================================
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	double QNoDiffusionElement<DIM,NNODE_1D>::
// 	 dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(
// 	 const unsigned &ipt,
// 	 Shape &psi, 
// 	 DShape &dpsidx,
// 	 Shape &test, 
// 	 DShape &dtestdx) const
// 	{
// 	 //Call the geometrical shape functions and derivatives  
// 	 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

// 	 //Set the test functions equal to the shape functions (pointer copy)
// 	 test = psi;
// 	 dtestdx = dpsidx;

// 	 //Return the jacobian
// 	 return J;
// 	}


// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////



// 	//=======================================================================
// 	/// \short Face geometry for the QNoDiffusionElement elements: 
// 	/// The spatial dimension of the face elements is one lower than that 
// 	/// of the bulk element but they have the same number of points along 
// 	/// their 1D edges.
// 	//=======================================================================
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	class FaceGeometry<QNoDiffusionElement<DIM,NNODE_1D> >: 
// 	 public virtual QElement<DIM-1,NNODE_1D>
// 	{

// 	  public:
	 
// 	 /// \short Constructor: Call the constructor for the
// 	 /// appropriate lower-dimensional QElement
// 	 FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

// 	};



// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////


// 	//=======================================================================
// 	/// Face geometry for the 1D QMonodomain elements: Point elements
// 	//=======================================================================
// 	template<unsigned NNODE_1D>
// 	class FaceGeometry<QNoDiffusionElement<1,NNODE_1D> >: 
// 	 public virtual PointElement
// 	{

// 	  public:
	 
// 	 /// \short Constructor: Call the constructor for the
// 	 /// appropriate lower-dimensional QElement
// 	 FaceGeometry() : PointElement() {}

// 	};


// 	//Override functions in specific implementations of the diff augmented wrapper
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	class DiffAugmentedCell<QNoDiffusionElement<DIM, NNODE_1D>>:
// 		public TMonodomainElement<DIM, NNODE_1D>
// 	{
// 	public:

// 		inline void get_diff_monodomain(const unsigned& ipt,
//                                         const Vector<double> &s,
//                                         const Vector<double>& x,
//                                         DenseMatrix<double>& D) const
// 		{
// 			this->get_interpolated_diffusion_matrix(s, D);
// 		}

// 	};





// 	/////////////////////////////////////////////////////////////////////////
// 	/////////////////////////////////////////////////////////////////////////
// 	// TNoDiffusionElement
// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////



// 	//======================================================================
// 	/// \short TNoDiffusionElement elements are isoparametric triangular 
// 	/// DIM-dimensional General Advection Diffusion Equations with  NNODE_1D nodal points along each
// 	/// element edge. Inherits from TElement and NoDiffusionEquations
// 	//======================================================================
// 	template <unsigned DIM, unsigned NNODE_1D>
// 	 class TNoDiffusionElement : 
// 	 public virtual TElement<DIM,NNODE_1D>,
// 	 public virtual NoDiffusionEquations<DIM>
// 	{

// 	private:

// 	 /// \short Static array of ints to hold number of variables at 
// 	 /// nodes: Initial_Nvalue[n]
// 	 static const unsigned Initial_Nvalue = 1;
	 
// 	  public:


// 	 ///\short  Constructor: Call constructors for TElement and 
// 	 /// Advection Diffusion equations
// 	 TNoDiffusionElement() : TElement<DIM,NNODE_1D>(), 
// 	  NoDiffusionEquations<DIM>()
// 	  { }

// 	 /// Broken copy constructor
// 	 TNoDiffusionElement(
// 	  const TNoDiffusionElement<DIM,NNODE_1D>&  dummy) 
// 	  { 
// 	   BrokenCopy::broken_copy("TNoDiffusionElement");
// 	  } 
	 
// 	 /// Broken assignment operator
// 	 void operator=(const TNoDiffusionElement<DIM,NNODE_1D>&) 
// 	  {
// 	   BrokenCopy::broken_assign("TNoDiffusionElement");
// 	  }

// 	 /// \short  Required  # of `values' (pinned or dofs) 
// 	 /// at node n
// 	 inline unsigned required_nvalue(const unsigned &n) const 
// 	  {return Initial_Nvalue;}

// 	 /// \short Output function:  
// 	 ///  x,y,u   or    x,y,z,u
// 	 void output(std::ostream &outfile)
// 	  {BaseCellMembranePotentialEquations<DIM>::output(outfile);}

// 	 /// \short Output function:  
// 	 ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
// 	 void output(std::ostream &outfile, const unsigned &n_plot)
// 	  {BaseCellMembranePotentialEquations<DIM>::output(outfile,n_plot);}


// 	 /// \short C-style output function:  
// 	 ///  x,y,u   or    x,y,z,u
// 	 void output(FILE* file_pt)
// 	  {
// 	   BaseCellMembranePotentialEquations<DIM>::output(file_pt);
// 	  }

// 	 ///  \short C-style output function:  
// 	 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
// 	 void output(FILE* file_pt, const unsigned &n_plot)
// 	  {
// 	   BaseCellMembranePotentialEquations<DIM>::output(file_pt,n_plot);
// 	  }

// 	 /// \short Output function for an exact solution:
// 	 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
// 	 void output_fct(std::ostream &outfile, const unsigned &n_plot,
// 	                 FiniteElement::SteadyExactSolutionFctPt 
// 	                 exact_soln_pt)
// 	  {BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,exact_soln_pt);}


// 	 /// \short Output function for a time-dependent exact solution.
// 	 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
// 	 /// (Calls the steady version)
// 	 void output_fct(std::ostream &outfile, const unsigned &n_plot,
// 	                 const double& time,
// 	                 FiniteElement::UnsteadyExactSolutionFctPt 
// 	                 exact_soln_pt)
// 	  {
// 	   BaseCellMembranePotentialEquations<DIM>::output_fct(outfile,n_plot,time,exact_soln_pt);
// 	  }


// 	protected:

// 	 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
// 	 inline double dshape_and_dtest_eulerian_BaseCellMembranePotential(
// 	  const Vector<double> &s, 
// 	  Shape &psi, 
// 	  DShape &dpsidx, 
// 	  Shape &test, 
// 	  DShape &dtestdx) const;
	 
// 	 /// \short Shape, test functions & derivs. w.r.t. to global coords. at
// 	 /// integration point ipt. Return Jacobian.
// 	 inline double dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(
// 	  const unsigned& ipt,
// 	  Shape &psi, 
// 	  DShape &dpsidx, 
// 	  Shape &test,
// 	  DShape &dtestdx) 
// 	  const;

// 	};

// 	//Inline functions:


// 	//======================================================================
// 	/// \short Define the shape functions and test functions and derivatives
// 	/// w.r.t. global coordinates and return Jacobian of mapping.
// 	///
// 	/// Galerkin: Test functions = shape functions
// 	//======================================================================
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	double TNoDiffusionElement<DIM,NNODE_1D>::
// 	 dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s,
// 	                                         Shape &psi, 
// 	                                         DShape &dpsidx,
// 	                                         Shape &test, 
// 	                                         DShape &dtestdx) const
// 	{
// 	 //Call the geometrical shape functions and derivatives  
// 	 double J = this->dshape_eulerian(s,psi,dpsidx);

// 	 //Loop over the test functions and derivatives and set them equal to the
// 	 //shape functions
// 	 for(unsigned i=0;i<NNODE_1D;i++)
// 	  {
// 	   test[i] = psi[i]; 
// 	   for(unsigned j=0;j<DIM;j++)
// 	    {
// 	     dtestdx(i,j) = dpsidx(i,j);
// 	    }
// 	  }
	 
// 	 //Return the jacobian
// 	 return J;
// 	}



// 	//======================================================================
// 	/// Define the shape functions and test functions and derivatives
// 	/// w.r.t. global coordinates and return Jacobian of mapping.
// 	///
// 	/// Galerkin: Test functions = shape functions
// 	//======================================================================
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	double TNoDiffusionElement<DIM,NNODE_1D>::
// 	 dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(
// 	 const unsigned &ipt,
// 	 Shape &psi, 
// 	 DShape &dpsidx,
// 	 Shape &test, 
// 	 DShape &dtestdx) const
// 	{
// 	 //Call the geometrical shape functions and derivatives  
// 	 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

// 	 //Set the test functions equal to the shape functions (pointer copy)
// 	 test = psi;
// 	 dtestdx = dpsidx;

// 	 //Return the jacobian
// 	 return J;
// 	}


// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////



// 	//=======================================================================
// 	/// \short Face geometry for the TNoDiffusionElement elements: 
// 	/// The spatial dimension of the face elements is one lower than that 
// 	/// of the bulk element but they have the same number of points along 
// 	/// their 1D edges.
// 	//=======================================================================
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	class FaceGeometry<TNoDiffusionElement<DIM,NNODE_1D> >: 
// 	 public virtual TElement<DIM-1,NNODE_1D>
// 	{

// 	  public:
	 
// 	 /// \short Constructor: Call the constructor for the
// 	 /// appropriate lower-dimensional TElement
// 	 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

// 	};



// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////
// 	////////////////////////////////////////////////////////////////////////


// 	//=======================================================================
// 	/// Face geometry for the 1D TNoDiffusionElement: Point elements
// 	//=======================================================================
// 	template<unsigned NNODE_1D>
// 	class FaceGeometry<TNoDiffusionElement<1,NNODE_1D> >: 
// 	 public virtual PointElement
// 	{

// 	  public:
	 
// 	 /// \short Constructor: Call the constructor for the
// 	 /// appropriate lower-dimensional TElement
// 	 FaceGeometry() : PointElement() {}

// 	};


// 	//Override functions in specific implementations of the diff augmented wrapper
// 	template<unsigned DIM, unsigned NNODE_1D>
// 	class DiffAugmentedCell<TNoDiffusionElement<DIM, NNODE_1D>>:
// 		public TMonodomainElement<DIM, NNODE_1D>
// 	{
// 	public:

// 		inline void get_diff_monodomain(const unsigned& ipt,
//                                         const Vector<double> &s,
//                                         const Vector<double>& x,
//                                         DenseMatrix<double>& D) const
// 		{
// 			this->get_interpolated_diffusion_matrix(s, D);
// 		}

// 	};



}


#endif