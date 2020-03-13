//LIC// ====================================================================
//LIC// This file contains the CellInterface cell model from Haibo Ni re-written as
//LIC// oomph-lib equations and elements. This is a MASSIVE undertaking...
//LIC//		+	Contains all cell functionality EXCEPT FOR MEMBRANE POTENTIAL
//LIC//			this is provided through combining with a monodomain element
//LIC//			in order to avoid overlapping functionality
//LIC//		+	Must be combined with Monodomain Element in order to function
//LIC//====================================================================

//!!!!!
//REQUIRED ALTERATIONS
//	Add required_concs		-	the required number of concentrations, this might change if
//								e.g. the conc of ATP is required by a channel formulation
//
//	Add no_repeated_cells	-	if the entries in the residual corresponding to the current
//								nodes cell variables, don't run the cell code, just skip it

//!!!!! PERHAPS IMPLEMENT THE DIFFUSION COEFFICIENTS IN THE CELL INTERFACE ELEMENTS

//Header file for CellInterface elements
#ifndef OOMPH_CELL_INTERFACE
#define OOMPH_CELL_INTERFACE

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

//include the cell model library
#include "../cell_model/cell_model_elements.h"

namespace oomph
{
	template <unsigned DIM>
	class CellInterfaceEquations : public virtual FiniteElement
	{
	public:

		//Define the function template used for forcing terms and stuff
		typedef void (*CellInterfaceScalarFctPt)
		(const double& t, const unsigned& ipt, const Vector<double>& s, const Vector<double>& x, double& Scal);

		CellInterfaceEquations() : 	Mutation_pt(0),
									Cell_model_pt(0),
									Membrane_potential_fct_pt_CellInterface(0),
									Strain_fct_pt(0),
									External_Na_conc_fct_pt_CellInterface(0),
									External_Ca_conc_fct_pt_CellInterface(0),
									External_K_conc_fct_pt_CellInterface(0),
									Ignore_Repeated_Cells(true)
		{	}

		CellInterfaceEquations(const CellInterfaceEquations& dummy){BrokenCopy::broken_copy("CellInterfaceEquations");}

		void operator=(const CellInterfaceEquations&){BrokenCopy::broken_assign("CellInterfaceEquations");}

		//Min and max variable indexes for output function and for ease of multiphysics elements
		virtual inline unsigned min_index_CellInterfaceEquations() const {return 0;}
		virtual inline unsigned max_index_CellInterfaceEquations() const {return (min_index_CellInterfaceEquations() + cell_model_pt()->Required_storage());}

		// Access functions to ignore repeated cells variable
		void ignore_repeated_cells(){Ignore_Repeated_Cells = true;}
		void do_not_ignore_repeated_cells(){Ignore_Repeated_Cells = false;}

		//Return the cell model pt
		CellModelBase* &cell_model_pt()
			{
				// std::cout << "in cell_model_pt()" << std::endl;
			#ifdef PARANOID
				if(Cell_model_pt == 0){
					//throw an error			    
					throw OomphLibError("No Cell model assigned to element Cell_interface_element",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
				}
			#endif
				return Cell_model_pt;
			}

		CellModelBase* const & cell_model_pt() const
			{
				// std::cout << "in const cell_model_pt()" << std::endl;
			#ifdef PARANOID
				if(Cell_model_pt == 0){
					//throw an error			    
					throw OomphLibError("No Cell model assigned to element Cell_interface_element",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
				}
			#endif
				return Cell_model_pt;
			}


		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			unsigned nplot=5;
			output(outfile,nplot);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot);
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			unsigned n_plot=5;
			output(file_pt,n_plot);
		}
		/// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// n_plot^DIM plot points
		void output(FILE* file_pt, const unsigned &n_plot);



		//====================================================================
		//====================================================================
		//====================================================================
		//Access functions
		//====================================================================
		//====================================================================
		//====================================================================


		//====================================================================
		//Membrane Potential Access function
		//(In the functional, multiphysics element this will be overwritten
		//	to interpolate the potential from the other parent element)
		//====================================================================
		inline virtual void get_membrane_potential_CellInterface(	const unsigned& ipt,
																	const Vector<double>& s,
																	const Vector<double>& x,
																	double& V) const
		{
			if(Membrane_potential_fct_pt_CellInterface!=0){
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*Membrane_potential_fct_pt_CellInterface)(time, ipt, s, x, V);
			}
			else{
				V=0.0;		//A default value
			}
		}

		//====================================================================
		//Strain access function
		//====================================================================
		//!!!!!
		// Might remove Strain_fct_pt, it serves little purpose
		//  this fct could be default implemented to just get the cell
		//  active strain, but overloaded in a tissue element to use the
		//  interpolated strain in the tissue
		inline virtual void get_strain_CellInterface(const unsigned& ipt,
											const Vector<double>& s,
											const Vector<double>& x,
											double& strain) const
		{
			if(Strain_fct_pt!=0){
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*Strain_fct_pt)(time, ipt, s, x, strain);
			}
			else{
				strain = get_interpolated_cell_active_strain(s);
			}
		}


		//====================================================================
		//External ionic concentrations
		//====================================================================
		inline virtual double get_external_Na_conc_CellInterface(const unsigned& ipt,
													const Vector<double>& s,
													const Vector<double>& x) const
		{
			if(External_Na_conc_fct_pt_CellInterface==0){
				return 140;		//Constant from unedited CNZ cell model from Haibo-Ni
			}
			else{
				double temp_conc;
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*External_Na_conc_fct_pt_CellInterface)(time, ipt, s, x, temp_conc);
				return temp_conc;
			}
		}
		inline virtual double get_external_Ca_conc_CellInterface(const unsigned& ipt,
													const Vector<double>& s,
													const Vector<double>& x) const
		{
			if(External_Ca_conc_fct_pt_CellInterface==0){
				return 1.8;		//Constant from unedited CNZ cell model from Haibo-Ni
			}
			else{
				double temp_conc;
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*External_Ca_conc_fct_pt_CellInterface)(time, ipt, s, x, temp_conc);
				return temp_conc;
			}
		}
		inline virtual double get_external_K_conc_CellInterface(const unsigned& ipt,
													const Vector<double>& s,
													const Vector<double>& x) const
		{
			if(External_K_conc_fct_pt_CellInterface==0){
				return 5.4;		//Constant from unedited CNZ cell model from Haibo-Ni
			}
			else{
				double temp_conc;
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*External_K_conc_fct_pt_CellInterface)(time, ipt, s, x, temp_conc);
				return temp_conc;
			}
		}

		//====================================================================
		//Cell type at node
		//	Used in switch function with the following correspondence:
		//		0 RA, 1 PM, 2 CT, 3 RAA, 4 AS, 5 AVR, 6 BB, 7 LA, 8 LAA, 9 PV, 10 SAN_C, 11 SAN_P
		//====================================================================
		virtual unsigned get_cell_type_at_node_CellInterface(const unsigned &n) const = 0;

		virtual unsigned get_fibrosis_type_at_node_CellInterface(const unsigned &n) const =0;

		//====================================================================
		//Drug Action
		//Default to zero
		//(Might be worth implementing these node-wise so diffusion of drugs could be investigated in future)
		//====================================================================


		//====================================================================
		//Access functions for mutation type
		//(Element-wise since we assume all cells are affected by the same mutation)
		//====================================================================
		const unsigned &mutation_CellInterface() const {return *Mutation_pt;}
		unsigned* &mutation_pt_CellInterface() {return Mutation_pt;}

		void set_mutation_pt(unsigned* mutation_pt_){Mutation_pt = mutation_pt_;}

		//====================================================================
		//Access functions to function pointers
		//====================================================================
		CellInterfaceScalarFctPt& membrane_potential_fct_pt_CellInterface()
			{return Membrane_potential_fct_pt_CellInterface;}

		CellInterfaceScalarFctPt& strain_fct_pt_CellInterface()
			{return Strain_fct_pt;}

		CellInterfaceScalarFctPt& external_Na_conc_fct_pt_CellInterface()
			{return External_Na_conc_fct_pt_CellInterface;}

		CellInterfaceScalarFctPt& external_Ca_conc_fct_pt_CellInterface()
			{return External_Ca_conc_fct_pt_CellInterface;}

		CellInterfaceScalarFctPt& external_K_conc_fct_pt_CellInterface()
			{return External_K_conc_fct_pt_CellInterface;}


		//====================================================================
		//====================================================================
		//Residual and Jacobian functions
		//====================================================================
		//====================================================================
		/// Add the element's contribution to its residual vector (wrapper)
		void fill_in_contribution_to_residuals(Vector<double> &residuals)
		{
		   	//Call the generic residuals function with flag set to 0 and using
		   	//a dummy matrix
			fill_in_generic_residual_contribution_cell_interface(
				residuals,GeneralisedElement::Dummy_matrix,
				GeneralisedElement::Dummy_matrix,0);
		}
		 
		/// \short Add the element's contribution to its residual vector and 
		/// the element Jacobian matrix (wrapper)
		void fill_in_contribution_to_jacobian(Vector<double> &residuals,
		                                   DenseMatrix<double> &jacobian)
		{
			//Call the generic routine with the flag set to 1
			// fill_in_generic_residual_contribution_cell_interface(residuals,jacobian,GeneralisedElement::Dummy_matrix,1);
			// DenseMatrix<double> temp_jacobian(jacobian.nrow());
			// FiniteElement::fill_in_contribution_to_jacobian(residuals,temp_jacobian);

			// const unsigned n_node = nnode();
			// for(unsigned l=0;l<n_node;l++){
			// 	Vector<int> local_eqn(cell_model_pt()->Required_storage());
			// 	for(unsigned var=0; var<cell_model_pt()->Required_storage(); var++){
			// 		local_eqn[var] = nodal_local_eqn(l, min_index_CellInterfaceEquations() + var);

			// 		if(local_eqn[var]>=0){
			// 			jacobian(local_eqn[var], local_eqn[var]) += temp_jacobian(local_eqn[var], local_eqn[var]);
			// 		}
			// 	}
			// }
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


		//====================================================================
		//Interpolated total membrane current
		//====================================================================
		inline double get_interpolated_cell_active_strain(const Vector<double> &s) const
		{
			//number of nodes in the element
			unsigned n_node = nnode();
			//The local and global coordinates of the node being considered
			Vector<double> s_node(DIM);
			Vector<double> x_node(DIM);
			//The values of the shape functions at the position interpolation is being calculated at
			Shape psi(n_node);
			shape(s,psi);
			//running total of the interpolated active strain
			double interpolated_active_strain = 0.0;
			// The local locations of the cell variables
			Vector<unsigned> local_ind(cell_model_pt()->Required_storage());

			//loop over nodes in the element and add their contributions
			for(unsigned n=0;n<n_node;n++){
				//Compile the local locations of the cell variables
				for(unsigned var=0; var<cell_model_pt()->Required_storage(); var++){local_ind[var] = min_index_CellInterfaceEquations() + var;}

				interpolated_active_strain += cell_model_pt()->active_strain(this->node_pt(n),local_ind)
															*psi[n];
			}

			return interpolated_active_strain;
		}

		//====================================================================
		//Interpolated total membrane current
		//====================================================================
		inline double interpolated_membrane_current_CellInterface(const Vector<double> &s) const
		{	
			//number of nodes in the element
			unsigned n_node = nnode();
			//The local and global coordinates of the node being considered
			Vector<double> s_node(DIM);
			Vector<double> x_node(DIM);
			//The values of the shape functions at the position interpolation is being calculated at
			Shape psi(n_node);
			shape(s,psi);

			//The thing we're calculating
			double interpolated_membrane_current=0.0;

			//loop over the nodes and add up their contributions
			for(unsigned n=0;n<n_node;n++){

				//Calculate the local and global coordinate of node n
				local_coordinate_of_node(n,s_node);
				for(unsigned j=0;j<DIM;j++){
					x_node[j] = raw_nodal_position(n,j);
				}

				//Calculate membrane potential
				double Vm;
				get_membrane_potential_CellInterface(0, s_node, x_node, Vm);
				
				//Calculate strain
				double strain;
				get_strain_CellInterface(0, s, x_node, strain);

				//Calculate external concentrations
				//Na, Ca, K
				Vector<double> Ext_conc(3);
				double temp_conc;
				
				Ext_conc[0] = get_external_Na_conc_CellInterface(0, s, x_node);
				Ext_conc[1] = get_external_Ca_conc_CellInterface(0, s, x_node);
				Ext_conc[2] = get_external_K_conc_CellInterface(0, s, x_node);

				//Compile the local locations of the cell variables
				Vector<unsigned> local_ind(cell_model_pt()->Required_storage());
				for(unsigned var=0; var<cell_model_pt()->Required_storage(); var++){
					local_ind[var] = min_index_CellInterfaceEquations() + var;
				}

				//Add nodal contribution to interpolated current
				interpolated_membrane_current += cell_model_pt()->membrane_current(	this->node_pt(n),		
																					Vm,
																					strain,
																					Ext_conc,
																					local_ind,
																					get_cell_type_at_node_CellInterface(n),
																					mutation_CellInterface(),
																					get_fibrosis_type_at_node_CellInterface(n)	)*psi[n];
			}
			if(interpolated_membrane_current==0.0){
				throw OomphLibError("did not get interpolated_membrane_current",
									OOMPH_CURRENT_FUNCTION,
									OOMPH_EXCEPTION_LOCATION);
			}
			return interpolated_membrane_current;
		}

	protected:
		//====================================================================
		//====================================================================
		//Residual and Jacobian functions
		//====================================================================
		//====================================================================
		/// \short Add the element's contribution to its residual vector only 
 		/// (if flag=and/or element  Jacobian matrix 
 		virtual void fill_in_generic_residual_contribution_cell_interface(
  			Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  			DenseMatrix<double> &mass_matrix, unsigned flag);


		//====================================================================
		//Definition of function pointers
		//====================================================================
		CellInterfaceScalarFctPt Membrane_potential_fct_pt_CellInterface;
		CellInterfaceScalarFctPt Strain_fct_pt;

		//External concentration function pointers
		CellInterfaceScalarFctPt External_Na_conc_fct_pt_CellInterface;
		CellInterfaceScalarFctPt External_Ca_conc_fct_pt_CellInterface;
		CellInterfaceScalarFctPt External_K_conc_fct_pt_CellInterface;


		//====================================================================
		//Drug action function pointers
		//====================================================================


		//====================================================================
		//Mutation type
		//	Used in switch function with the following correspondence:
		//		0 WT, 1 D322H, 2 E48G, 3 A305T, 4 Y155C, 5 D469E, 6 P488S
		//====================================================================
		unsigned *Mutation_pt;

		//Pointer to the cell model
		CellModelBase *Cell_model_pt;

		// If true, before calculating single cell for a node, check if residual
		//	entries corresponding to that node are zero. If any are not zero, do
		//	not calculate single cell.
		bool Ignore_Repeated_Cells;


	private:
	};


	//====================================================================
	//====================================================================
	//Q Element
	//====================================================================
	//====================================================================

	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class QCellInterfaceElement	:
		public virtual QElement<DIM, NNODE_1D>,
		public virtual CellInterfaceEquations<DIM>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		QCellInterfaceElement()	:	QElement<DIM, NNODE_1D>(),
									CellInterfaceEquations<DIM>()
		{
			//Create data for cell type and fibrosis type and pin them immediately
			Cell_type_internal_index = this->add_internal_data(new Data(this->nnode()), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Cell_type_internal_index)->pin(l);
			}
			Fibrosis_type_internal_index = this->add_internal_data(new Data(this->nnode()), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Fibrosis_type_internal_index)->pin(l);
			}
		}

		QCellInterfaceElement(const QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>& dummy){BrokenCopy::broken_copy("QCellInterfaceElement");}

		void operator=(const QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>&){BrokenCopy::broken_assign("QCellInterfaceElement");}

		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			CellInterfaceEquations<DIM>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			CellInterfaceEquations<DIM>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			CellInterfaceEquations<DIM>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	CellInterfaceEquations<DIM>::output(file_pt, n_plot);
		}




		inline unsigned required_nvalue(const unsigned &n) const {return NUM_VARS;}//this->cell_model_pt()->Required_storage();}

		void set_cell_model_pt(CellModelBase* cell_model_pt_){
			//Check if the required number of values from cell_model_pt is the same as that passed to the element constructor
			if(NUM_VARS!=cell_model_pt_->Required_storage()){
				//throw an error
				std::string error_message =
						"The number of variables passed to the QCellInterfaceElement constructor (";
			    error_message += std::to_string(NUM_VARS);
			    error_message += ") does not match\n\tthe number defined by the Cell_model_pt (";
			    error_message += std::to_string(cell_model_pt_->Required_storage());
			    error_message += ".";
			    
			   	throw OomphLibError(error_message,
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			else{
				//set the cell_model_pt
				this->Cell_model_pt = cell_model_pt_;	
			}
		}

		//====================================================================
		//Cell type set and access functions
		//====================================================================
		unsigned  get_cell_type_at_node_CellInterface(const unsigned &n) const
		{return this->internal_data_pt(Cell_type_internal_index)->value(n);}

		void set_cell_type(const unsigned &n, const unsigned &cell_type){this->internal_data_pt(Cell_type_internal_index)->set_value(n, cell_type);}
		unsigned get_cell_type_at_node(const unsigned &n){return this->internal_data_pt(Cell_type_internal_index)->value(n);}

		//====================================================================
		//Fibrosis type set and access functions
		//====================================================================
		unsigned  get_fibrosis_type_at_node_CellInterface(const unsigned &n) const {return this->internal_data_pt(Fibrosis_type_internal_index)->value(n);}
		
		void set_fibrosis_type(const unsigned &n, const unsigned &fibrosis_type){this->internal_data_pt(Fibrosis_type_internal_index)->set_value(n, fibrosis_type);}
		unsigned get_fibrosis_type_at_node(const unsigned &n){return this->internal_data_pt(Fibrosis_type_internal_index)->value(n);}


	protected:
		unsigned Cell_type_internal_index;

		unsigned Fibrosis_type_internal_index;
	};

	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class FaceGeometry<QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D> >:
		public virtual QElement<DIM-1, NNODE_1D>
	{
	public:
		FaceGeometry()	:	QElement<DIM-1, NNODE_1D>()	{}
	};	

	template<unsigned NUM_VARS, unsigned NNODE_1D>
	class FaceGeometry<QCellInterfaceElement<1, NUM_VARS, NNODE_1D> >:
		public virtual PointElement
	{
	public:
		FaceGeometry()	:	PointElement()	{}
	};


	//====================================================================
	//====================================================================
	//T Element
	//====================================================================
	//====================================================================

	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class TCellInterfaceElement	:
		public virtual TElement<DIM, NNODE_1D>,
		public virtual CellInterfaceEquations<DIM>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		TCellInterfaceElement()	:	TElement<DIM, NNODE_1D>(),
									CellInterfaceEquations<DIM>()
		{
			//Create data for cell type and fibrosis type and pin them immediately
			Cell_type_internal_index = this->add_internal_data(new Data(this->nnode()), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Cell_type_internal_index)->pin(l);
			}
			Fibrosis_type_internal_index = this->add_internal_data(new Data(this->nnode()), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Fibrosis_type_internal_index)->pin(l);
			}
		}

		TCellInterfaceElement(const TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>& dummy){BrokenCopy::broken_copy("TCellInterfaceElement");}

		void operator=(const TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>&){BrokenCopy::broken_assign("TCellInterfaceElement");}

		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			CellInterfaceEquations<DIM>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			CellInterfaceEquations<DIM>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			CellInterfaceEquations<DIM>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	CellInterfaceEquations<DIM>::output(file_pt, n_plot);
		}




		inline unsigned required_nvalue(const unsigned &n) const {return NUM_VARS;}//this->cell_model_pt()->Required_storage();}

		void set_cell_model_pt(CellModelBase* cell_model_pt_){
			//Check if the required number of values from cell_model_pt is the same as that passed to the element constructor
			if(NUM_VARS!=cell_model_pt_->Required_storage()){
				//throw an error
				std::string error_message =
						"The number of variables passed to the TCellInterfaceElement constructor (";
			    error_message += NUM_VARS;
			    error_message += ") does not match\n\tthe number defined by the Cell_model_pt (";
			    error_message += cell_model_pt_->Required_storage();
			    error_message += ".";
			    
			   	throw OomphLibError(error_message,
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			else{
				//set the cell_model_pt
				this->Cell_model_pt = cell_model_pt_;
			}
		}

		//====================================================================
		//Cell type set and access functions
		//====================================================================
		unsigned  get_cell_type_at_node_CellInterface(const unsigned &n) const
		{return this->internal_data_pt(Cell_type_internal_index)->value(n);}

		void set_cell_type(const unsigned &n, const unsigned &cell_type){this->internal_data_pt(Cell_type_internal_index)->set_value(n, cell_type);}
		unsigned get_cell_type_at_node(const unsigned &n){return this->internal_data_pt(Cell_type_internal_index)->value(n);}

		//====================================================================
		//Fibrosis type set and access functions
		//====================================================================
		unsigned  get_fibrosis_type_at_node_CellInterface(const unsigned &n) const {return this->internal_data_pt(Fibrosis_type_internal_index)->value(n);}
		
		void set_fibrosis_type(const unsigned &n, const unsigned &fibrosis_type){this->internal_data_pt(Fibrosis_type_internal_index)->set_value(n, fibrosis_type);}
		unsigned get_fibrosis_type_at_node(const unsigned &n){return this->internal_data_pt(Fibrosis_type_internal_index)->value(n);}


	protected:
		unsigned Cell_type_internal_index;

		unsigned Fibrosis_type_internal_index;
	};

	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class FaceGeometry<TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D> >:
		public virtual TElement<DIM-1, NNODE_1D>
	{
	public:
		FaceGeometry()	:	TElement<DIM-1, NNODE_1D>()	{}
	};	

	template<unsigned NUM_VARS, unsigned NNODE_1D>
	class FaceGeometry<TCellInterfaceElement<1, NUM_VARS, NNODE_1D> >:
		public virtual PointElement
	{
	public:
		FaceGeometry()	:	PointElement()	{}
	};




	//====================================================================
	//====================================================================
	//Point Element
	//How many dimensions should I assign to this element?
	//====================================================================
	//====================================================================

	template<unsigned DIM, unsigned NUM_VARS>
	class PointCellInterfaceElement	:
		public virtual PointElement,
		public virtual CellInterfaceEquations<DIM>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		PointCellInterfaceElement()	:	PointElement(),
										CellInterfaceEquations<DIM>()
		{
			//Create data for cell type and fibrosis type and pin them immediately
			Cell_type_internal_index = this->add_internal_data(new Data(1), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Cell_type_internal_index)->pin(l);
			}
			Fibrosis_type_internal_index = this->add_internal_data(new Data(1), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Fibrosis_type_internal_index)->pin(l);
			}
		}
		PointCellInterfaceElement(const PointCellInterfaceElement<DIM,NUM_VARS>& dummy){BrokenCopy::broken_copy("PointCellInterfaceElement");}

		void operator=(const PointCellInterfaceElement<DIM,NUM_VARS>&){BrokenCopy::broken_assign("PointCellInterfaceElement");}

		//====================================================================
		//Output functions
		//====================================================================
		// /// Output with default number of plot points
		// void output(std::ostream &outfile){
		// 	PointCellInterfaceElement<DIM, NUM_VARS>::output(outfile, 0);
		// }
		// /// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		// /// nplot^DIM plot points
		// void output(std::ostream &outfile, const unsigned &nplot);
		// /// C_style output with default number of plot points
		// void output(FILE* file_pt){
		// 	PointCellInterfaceElement<DIM, NUM_VARS>::output(file_pt, 0);
		// }
		//  /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		//  /// n_plot^DIM plot points
		//  void output(FILE* file_pt, const unsigned &n_plot);

		 /// Output with default number of plot points
		void output(std::ostream &outfile){
			CellInterfaceEquations<DIM>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			CellInterfaceEquations<DIM>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			CellInterfaceEquations<DIM>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	CellInterfaceEquations<DIM>::output(file_pt, n_plot);
		}


		inline unsigned required_nvalue(const unsigned &n) const {return NUM_VARS;}//this->cell_model_pt()->Required_storage();}

		void set_cell_model_pt(CellModelBase* cell_model_pt_){
			//Check if the required number of values from cell_model_pt is the same as that passed to the element constructor
			if(NUM_VARS!=cell_model_pt_->Required_storage()){
				//throw an error
				std::string error_message =
			    "The number of variables passed to the PointCellInterfaceElement constructor (";
			    error_message += NUM_VARS;
			    error_message += ") does not match\n\tthe number defined by the Cell_model_pt (";
			    error_message += cell_model_pt_->Required_storage();
			    error_message += ".";
			    
			   	throw OomphLibError(error_message,
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			else{
				//set the cell_model_pt
				this->Cell_model_pt = cell_model_pt_;
			}
		}

		//====================================================================
		//Cell type set and access functions
		//====================================================================
		unsigned  get_cell_type_at_node_CellInterface(const unsigned &n) const
		{return this->internal_data_pt(Cell_type_internal_index)->value(0);}

		void set_cell_type(const unsigned &n, const unsigned &cell_type){this->internal_data_pt(Cell_type_internal_index)->set_value(0, cell_type);}
		unsigned get_cell_type_at_node(const unsigned &n){return this->internal_data_pt(Cell_type_internal_index)->value(0);}

		//====================================================================
		//Fibrosis type set and access functions
		//====================================================================
		unsigned  get_fibrosis_type_at_node_CellInterface(const unsigned &n) const {return this->internal_data_pt(Fibrosis_type_internal_index)->value(0);}
		
		void set_fibrosis_type(const unsigned &n, const unsigned &fibrosis_type){this->internal_data_pt(Fibrosis_type_internal_index)->set_value(0, fibrosis_type);}
		unsigned get_fibrosis_type_at_node(const unsigned &n){return this->internal_data_pt(Fibrosis_type_internal_index)->value(0);}


	protected:
		unsigned Cell_type_internal_index;

		unsigned Fibrosis_type_internal_index;
	};
}

#endif