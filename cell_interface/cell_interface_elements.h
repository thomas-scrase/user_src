//LIC// ====================================================================
//LIC// This file contains the CellInterface cell model from Haibo Ni re-written as
//LIC// oomph-lib equations and elements. This is a MASSIVE undertaking...
//LIC//		+	Contains all cell functionality EXCEPT FOR MEMBRANE POTENTIAL
//LIC//			this is provided through combining with a monodomain element
//LIC//			in order to avoid overlapping functionality
//LIC//		+	Must be combined with Monodomain Element in order to function
//LIC//			- that is unless you want to prescribe membrane potential
//LIC//====================================================================

//!!!!!
//REQUIRED ALTERATIONS
//	Add no_repeated_cells

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
#include "../cell_model/cell_model_base.h"
#include "../cell_model/cell_state_container.h"

//For the custom integration scheme
#include "../toms_utilities/toms_integral.h"


namespace oomph
{
	template <unsigned DIM, unsigned NUM_VARS>
	class CellInterfaceEquations : public virtual FiniteElement
	{
	public:

		//Define the function template used for forcing terms and stuff
		typedef void (*CellInterfaceScalarFctPt)
		(const double& t, const unsigned& ipt, const Vector<double>& s, const Vector<double>& x, double& Scal);

		// \short function pointer to boundary source function fct(bounds, f(bounds)) --
		// bounds_of_node is a vector of the bounds the node exists on
		typedef void (* CellInterfaceBoundarySourceFctPt)
		(std::set<unsigned>* &boundaries_pt, double& bound_source);

		CellInterfaceEquations() :	Cell_model_pt(0),
									Membrane_potential_fct_pt_CellInterface(0),
									Strain_fct_pt(0),
									Boundary_source_fct_pt(0),
									Ignore_Repeated_Cells(true)
		{	}

		CellInterfaceEquations(const CellInterfaceEquations& dummy){BrokenCopy::broken_copy("CellInterfaceEquations");}

		void operator=(const CellInterfaceEquations&){BrokenCopy::broken_assign("CellInterfaceEquations");}

		//Min and max variable indexes for output function and for ease of multiphysics elements
		virtual inline unsigned min_index_CellInterfaceEquations() const {return 0;}
		virtual inline unsigned max_index_CellInterfaceEquations() const {return (min_index_CellInterfaceEquations() + cell_model_pt()->required_nodal_variables());}

		// Access functions to ignore repeated cells variable
		void ignore_repeated_cells(){Ignore_Repeated_Cells = true;}
		void do_not_ignore_repeated_cells(){Ignore_Repeated_Cells = false;}

		//====================================================================
		//Handle cell model pt
		//====================================================================
		//Return the cell model pt
		// CellModelBase* &cell_model_pt(){
		// 	std::cout << "boom" << std::endl;
		// 	#ifdef PARANOID
		// 	if(Cell_model_pt == 0){
		// 		//throw an error			    
		// 		throw OomphLibError("No Cell model assigned to element Cell_interface_element",
		// 		OOMPH_CURRENT_FUNCTION,
		// 		OOMPH_EXCEPTION_LOCATION);
		// 	}
		// 	#endif
		// 	return Cell_model_pt;
		// }

		CellModelBase* const & cell_model_pt() const{
			// std::cout << "BOOM" << std::endl;
			// std::cout << "BOOM" << ipt_not_at_nodes << std::endl;
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

		void set_cell_model_pt(CellModelBase* cell_model_pt_){
			// std::cout << "BOOM" << std::endl;
			//Check if the required number of values from cell_model_pt is the same as that passed to the element constructor
			if(NUM_VARS!=cell_model_pt_->required_nodal_variables()){
				//throw an error
				std::string error_message =
						"The number of variables passed to the QCellInterfaceElement constructor (";
			    error_message += std::to_string(NUM_VARS);
			    error_message += ") does not match\n\tthe number defined by the Cell_model_pt (";
			    error_message += std::to_string(cell_model_pt_->required_nodal_variables());
			    error_message += ".";
			    
			   	throw OomphLibError(error_message,
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}

			//set the cell_model_pt
			this->Cell_model_pt = cell_model_pt_;

			//build the required nodal parameters

			//Create data for cell type and pin them immediately
			Cell_type_internal_index = this->add_internal_data(new Data(this->nnode()), false);
			for(unsigned l=0;l<this->nnode();l++)
			{
				this->internal_data_pt(Cell_type_internal_index)->pin(l);
			}

			//Create data for black-box nodal parameters and pin them immediately
			Black_box_nodal_parameters_internal_index = this->add_internal_data(new Data(this->nnode()*cell_model_pt_->required_black_box_parameters()), false);
			for(unsigned l=0;l<this->nnode()*cell_model_pt_->required_black_box_parameters();l++)
			{
				this->internal_data_pt(Black_box_nodal_parameters_internal_index)->set_value(l,0.0);
				this->internal_data_pt(Black_box_nodal_parameters_internal_index)->pin(l);
			}
		
		}

		//dependent on cell type to allow for distinction between cell models
		//	in a combined cell model class.
		inline unsigned required_nvalue(const unsigned &n) const {
			// return cell_model_pt()->required_nodal_variables(get_cell_type_at_node(n));
			return NUM_VARS;
		}

		/////////////////////////////////////////////////////////////////////////////////
		//Get and Set locally stored data which is to be passed to the cell model
		/////////////////////////////////////////////////////////////////////////////////

		//get the t-th history value of the v-th cell variable associated with the l-th node in the element
		inline double get_nodal_cell_variable(const unsigned &t, const unsigned &l, const unsigned &v) const {
			// return get_nodal_value(t, l, min_index_CellInterfaceEquations()+v);
			return node_pt(l)->value(t,v);
		}

		//get the d-th time derivative of the v-th cell variable associated with the l-th node in the element
		inline double get_nodal_cell_variable_derivative(const unsigned &l, const unsigned &d, const unsigned &v) const {
			// Get the data's timestepper
			TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();
			double dvdt=0.0;
			//Loop over the timesteps, if there is a non Steady timestepper
			if (!time_stepper_pt->is_steady()){
				//Initialise dudt
				const unsigned var_index = min_index_CellInterfaceEquations()+v;
				// Number of timsteps (past & present)
				const unsigned n_time = time_stepper_pt->ntstorage();
				for(unsigned t=0;t<n_time;t++){
					dvdt += time_stepper_pt->weight(d,t)*get_nodal_cell_variable(t,l,var_index);
				}
			}
			return dvdt;
		}

		//Cell type at node
		//	Used in switch function with the following correspondence:
		//		ATRIA 0 - 99,
		//		VENTS 100 - 199,
		//		OTHER 200 - 299 (?)
		//		CNZAtria
		//		0 RA, 1 PM, 2 CT, 3 RAA, 4 AS, 5 AVR, 6 BB, 7 LA, 8 LAA, 9 PV, 10 SAN_C, 11 SAN_P
		//		TNNPVentricle
		//		100 LVEPI 101 LVMCELL 102 LVENDO 103 RVEPI 104 RVMCELL 105 RVENDO 106 PFI 107 PFMB 108 PF
		void set_cell_type(const unsigned &n, const unsigned &cell_type){
			this->internal_data_pt(Cell_type_internal_index)->set_value(n, cell_type);
		}
		unsigned get_cell_type_at_node(const unsigned &n) const {
			return this->internal_data_pt(Cell_type_internal_index)->value(n);
		}

		//set the var-th black box nodal parameter associated with the l-th node to value
		inline void set_black_box_nodal_parameters(const unsigned &l, const unsigned &var, const double &value) const{
			this->internal_data_pt(Black_box_nodal_parameters_internal_index)
					->set_value(l*cell_model_pt()->required_black_box_parameters()+var, value);
		}

		// //get the black box nodal parameters for the l-th node
		// inline Vector<double> get_black_box_nodal_parameters(const unsigned &l) const {
		// 	Vector<double> black_box_parameters(cell_model_pt()->required_black_box_parameters());
		// 	for(unsigned var=0; var<cell_model_pt()->required_black_box_parameters(); var++){
		// 		black_box_parameters[var] = this->internal_data_pt(Black_box_nodal_parameters_internal_index)
		// 										->value(l*cell_model_pt()->required_black_box_parameters()+var);
		// 	}
		// }

		//====================================================================
		//====================================================================
		//Fill in data of state container
		//====================================================================
		//====================================================================
		void fill_state_container_at_node(CellState &state, const unsigned &l) const {
			//Loop through the variables which are requested by the cell model and
			//	grab the data from the interface element

			//Fill in current state and derivative matrix
			DenseMatrix<double> new_vars;
			new_vars.resize(cell_model_pt()->required_derivatives()+1, cell_model_pt()->required_nodal_variables());
			for(unsigned v = 0; v < cell_model_pt()->required_nodal_variables(); v++){
				new_vars(0,v) = get_nodal_cell_variable(0,l,min_index_CellInterfaceEquations() + v);
				for(unsigned d = 1; d <= cell_model_pt()->required_derivatives(); d++){
					new_vars(d,v) = get_nodal_cell_variable_derivative(l, d, min_index_CellInterfaceEquations()+v);
				}
			}
			state.set_vars(new_vars);

			//FIll in time_stepper_weights
			DenseMatrix<double> new_weights;
			new_weights.resize(cell_model_pt()->required_derivatives()+1, cell_model_pt()->required_nodal_variables());
			for(unsigned v = 0; v < cell_model_pt()->required_nodal_variables(); v++){
				new_weights(0,v) = 1.0;
				for(unsigned d = 1; d <= cell_model_pt()->required_derivatives(); d++){
					new_weights(d,v) = this->node_pt(l)->time_stepper_pt()->weight(d, 0);
				}
			}
			state.set_time_stepper_weights(new_weights);

			//set the cell type
			state.set_cell_type(get_cell_type_at_node(l));

			//fill in black box parameters
			Vector<double> black_box_parameters(cell_model_pt()->required_black_box_parameters());
			for(unsigned var=0; var<cell_model_pt()->required_black_box_parameters(); var++){
				black_box_parameters[var] = this->internal_data_pt(Black_box_nodal_parameters_internal_index)
												->value(l*cell_model_pt()->required_black_box_parameters()+var);
			}
			state.set_black_box_nodal_parameters(black_box_parameters);

			//Fill in transmembrane potential
			state.set_vm(get_nodal_membrane_potential(l));
			
			//Fill in mechanical stress
			state.set_stress(get_nodal_mechanical_strain(l));

			//Fill in black box external data
			Vector<double> black_box_external_data;
			get_nodal_black_box_external_data(l, black_box_external_data);
			state.set_black_box_external_data(black_box_external_data);
			
			//Fill in dt
			state.set_dt(node_pt(l)->time_stepper_pt()->time_pt()->dt(0));

			//Fill in previous time variables
			Vector<double> New_previous_values;
			New_previous_values.resize(cell_model_pt()->required_nodal_variables());
			for(unsigned v=0; v<cell_model_pt()->required_nodal_variables();v++){
				New_previous_values[v] = get_nodal_cell_variable(1, l, min_index_CellInterfaceEquations()+v);
			}
			state.set_previous_variables(New_previous_values);
			

		}



		/////////////////////////////////////////////////////////////////////////////////
		//Get Data from the cell model which is to be passed to external sources
		/////////////////////////////////////////////////////////////////////////////////

		//get the membrane current from the cell model at the l-th node
		inline double get_nodal_membrane_current(const unsigned &l) const
		{
			double nodal_membrane_current = 0.0;

			//Construct the state container
			CellState state;
			fill_state_container_at_node(state, l);

			//Add nodal contribution to interpolated current
			nodal_membrane_current += cell_model_pt()->membrane_current(state);

			//If Boundary_source_fct_pt has been set, get the contribution from the node
			//  This check prevents bulk non boundary elements from contributing
			//  unnecessary overhead
			if(Boundary_source_fct_pt!=0){
				// Preallocate boundaries the node is on
				std::set<unsigned>* boundaries_pt;
				// Get the pointer to set of boundaries node lies on
				node_pt(l)->get_boundaries_pt(boundaries_pt);
				// If the set is non-zero, get a contribution to nodal_membrane_current
				if(boundaries_pt!=0){
					double bound_source = 0.0;
					Boundary_source_fct_pt(boundaries_pt, bound_source);
					nodal_membrane_current += bound_source;
				}
			}

			return nodal_membrane_current;
		}
		//get the interpolated membrane current from the cell model at the local coordinate s
		//	Uses the get_nodal function because we want to ensure we account for any boundary
		//	source contributions from the nodes
		inline double get_interpolated_membrane_current_CellInterface(const Vector<double> &s) const
		{	
			//number of nodes in the element
			unsigned n_node = nnode();
			
			//The values of the shape functions at the position interpolation is being calculated at
			Shape psi(n_node);
			shape(s,psi);

			//The thing we're calculating
			double interpolated_membrane_current=0.0;

			//loop over the nodes and add up their contributions
			for(unsigned n=0;n<n_node;n++){
				interpolated_membrane_current += get_nodal_membrane_current(n)*psi[n];
			}

			return interpolated_membrane_current;
		}


		//get the active stress from the cell model at the l-th node
		inline double get_nodal_active_stress(const unsigned &l) const {			
			CellState state;
			fill_state_container_at_node(state, l);

			return cell_model_pt()->active_strain(state);
		}
		//get the interpolated active stress from the cell model at the local coordinate s
		inline double get_interpolated_active_stress(const Vector<double> &s) const
		{
			//number of nodes in the element
			unsigned n_node = nnode();
			//The local and global coordinates of the node being considered
			Vector<double> s_node(DIM);
			//The values of the shape functions at the position interpolation is being calculated at
			Shape psi(n_node);
			shape(s,psi);
			//running total of the interpolated active stress
			double interpolated_active_stress = 0.0;

			CellState state;

			//loop over nodes in the element and add their contributions
			for(unsigned n=0;n<n_node;n++){
				fill_state_container_at_node(state, n);
				interpolated_active_stress += cell_model_pt()->active_strain(state)*psi[n];
			}

			return interpolated_active_stress;
		}

		//get the membrane capacitance from the cell model at the l-th node
		inline double get_nodal_membrane_capacitance(const unsigned &l){
			CellState state;
			fill_state_container_at_node(state, l);

			return cell_model_pt()->cm(state);
			// return 1.0;
		}

		//get the interpolated membrane capacitance from the cell model at the local coordinate s
		inline double get_interpolated_membrane_capacitance(const unsigned &l){
			//number of nodes in the element
			unsigned n_node = nnode();
			//The local and global coordinates of the node being considered
			Vector<double> s_node(DIM);
			//The values of the shape functions at the position interpolation is being calculated at
			Shape psi(n_node);
			shape(s_node,psi);
			//running total of the interpolated active stress
			double interpolated_membrane_capacitance = 0.0;

			CellState state;

			//loop over nodes in the element and add their contributions
			for(unsigned n=0;n<n_node;n++){
				fill_state_container_at_node(state,n);
				interpolated_membrane_capacitance += cell_model_pt()->cm(state)*psi[n];
			}
			// interpolated_membrane_capacitance=1.0;
			return interpolated_membrane_capacitance;
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
		//predefined external sources
		//====================================================================

		//get the membrane potential at the coordinate of the l-th node
		inline double get_nodal_membrane_potential(const unsigned &l) const {
			//get the local and global coordinates of the node
			unsigned ipt_node = ipt_at_node(l);
			Vector<double> s_node(DIM);
			Vector<double> x_node(DIM);
			local_coordinate_of_node(l,s_node);
			for(unsigned j=0;j<DIM;j++){x_node[j] = raw_nodal_position(l,j);}
			double Vm;
			//call the external communicator function
			get_membrane_potential_CellInterface(ipt_node, s_node, x_node, Vm);
			return Vm;
		}

		//Membrane Potential Access function
		//(In the functional, multiphysics element this will be overwritten
		//	to interpolate the potential from the other parent element)
		inline virtual void get_membrane_potential_CellInterface(const unsigned& ipt,
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

		//get the mechanical strain at the coordinate of the l-th node
		inline double get_nodal_mechanical_strain(const unsigned &l) const {
			//get the local and global coordinates of the node
			unsigned ipt_node = ipt_at_node(l);
			Vector<double> s_node(DIM);
			Vector<double> x_node(DIM);
			local_coordinate_of_node(l,s_node);
			for(unsigned j=0;j<DIM;j++){x_node[j] = raw_nodal_position(l,j);}
			double strain;
			//call the external communicator function
			get_mechanical_strain_CellInterface(ipt_node, s_node, x_node, strain);
			return strain;
		}

		//Mechanical strain access functions
		inline virtual void get_mechanical_strain_CellInterface(const unsigned& ipt,
																const Vector<double>& s,
																const Vector<double>& x,
																double& strain) const
		{
			if(Strain_fct_pt!=0){
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*Strain_fct_pt)(time, ipt, s, x, strain);
			}
			else{
				strain = 0.0;
			}
		}

		// //get nodal black-box external data
		// inline Vector<double> get_nodal_black_box_external_data(const unsigned &l) const {
		// 	//get the local and global coordinates of the node
		// 	unsigned ipt_node = ipt_at_node(l);
		// 	Vector<double> s_node(DIM);
		// 	Vector<double> x_node(DIM);
		// 	local_coordinate_of_node(l,s_node);
		// 	for(unsigned j=0;j<DIM;j++){x_node[j] = raw_nodal_position(l,j);}
		// 	Vector<double> Ext_data;
		// 	get_black_box_external_data_CellInterface(ipt_node, s_node, x_node, Ext_data);
		// 	return Ext_data;
		// }
		inline void get_nodal_black_box_external_data(const unsigned &l, Vector<double> &data) const {
			//get the local and global coordinates of the node
			unsigned ipt_node = ipt_at_node(l);
			Vector<double> s_node(DIM);
			Vector<double> x_node(DIM);
			local_coordinate_of_node(l,s_node);
			for(unsigned j=0;j<DIM;j++){x_node[j] = raw_nodal_position(l,j);}
			get_black_box_external_data_CellInterface(ipt_node, s_node, x_node, data);
		}

		//get black-box external data
		//by default this function runs over the external function pointer which are provided to the element
		//virtual so it can be overloaded in multi-physics problems where the black-box external
		// sources may be defined by an external set of equations
		inline virtual void get_black_box_external_data_CellInterface(const unsigned& ipt,
																		const Vector<double>& s,
																		const Vector<double>& x,
																		Vector<double>& Ext_Data) const
		{
			Ext_Data.resize(Black_box_external_fct_pts.size());
			for(unsigned i=0; i<Black_box_external_fct_pts.size(); i++){
				double val = 0.0;
				double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
				(*Black_box_external_fct_pts[i])(time, ipt, s, x, val);
				Ext_Data[i] = val;
			}
		}

	

		//====================================================================
		//Access functions to predefined external function pointers
		//====================================================================
		CellInterfaceScalarFctPt& membrane_potential_fct_pt_CellInterface()
			{return Membrane_potential_fct_pt_CellInterface;}

		CellInterfaceScalarFctPt& strain_fct_pt_CellInterface()
			{return Strain_fct_pt;}

		/// Access function: Pointer to boundary source function
		CellInterfaceBoundarySourceFctPt& boundary_source_fct_pt() 
			{return Boundary_source_fct_pt;}

		/// Access function: Pointer to boundary source function
		CellInterfaceBoundarySourceFctPt boundary_source_fct_pt() const
			{return Boundary_source_fct_pt;}




		//====================================================================
		//====================================================================
		//Assign initial conditions at nodes from cell model
		//====================================================================
		//====================================================================
		inline void assign_initial_conditions_from_cell_model(){
			double current_var;
			for(unsigned n=0; n < this->nnode(); n++){
				for(unsigned v=0; v<NUM_VARS;v++){
					//try and get a default value for the v-th variable for node n,
					if(cell_model_pt()->return_initial_state_variable(v,current_var,get_cell_type_at_node(n))){
						this->node_pt(n)->set_value(min_index_CellInterfaceEquations() + v, current_var);
					}
					//if a suitable value cannot be returned by the cell model then the value is pinned
					else{
						this->node_pt(n)->pin(min_index_CellInterfaceEquations() + v);
					}
				}
			}
		}

		//====================================================================
		//Nodal cell model custom output
		//====================================================================
		inline void get_nodal_cell_custom_output(const unsigned &l, Vector<double> &output) const
		{
			CellState state;
			fill_state_container_at_node(state, l);

			cell_model_pt()->custom_output(state,output);
		}

		//return the Gauss point associated with node n
		inline unsigned ipt_at_node(const unsigned &n) const
			{return ipt_not_at_nodes + n;}


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

			if(cell_model_pt()->model_calculates_jacobian_entries()){
				//if the cell model is capable, allow it to calculate the jacobian
				fill_in_generic_residual_contribution_cell_interface(residuals,jacobian,GeneralisedElement::Dummy_matrix,1);
			}
			else{
				//Otherwise perform fill in procedure using finite differencing
				FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
			}
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

 		CellInterfaceBoundarySourceFctPt Boundary_source_fct_pt;

		CellInterfaceScalarFctPt Membrane_potential_fct_pt_CellInterface;
		CellInterfaceScalarFctPt Strain_fct_pt;

		//External black box function pointers
		Vector<CellInterfaceScalarFctPt> Black_box_external_fct_pts;


		//Pointer to the cell model
		CellModelBase *Cell_model_pt;

		// If true, before calculating single cell for a node, check if residual
		//	entries corresponding to that node are zero. If any are not zero, do
		//	not calculate single cell.
		bool Ignore_Repeated_Cells;

		unsigned ipt_not_at_nodes;




	private:
		unsigned Cell_type_internal_index;
		unsigned Black_box_nodal_parameters_internal_index;
	};








	//====================================================================
	//====================================================================
	//Q Element
	//====================================================================
	//====================================================================

	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class QCellInterfaceElement	:
		public virtual QElement<DIM, NNODE_1D>,
		public virtual CellInterfaceEquations<DIM, NUM_VARS>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		QCellInterfaceElement()	:	QElement<DIM, NNODE_1D>(),
									CellInterfaceEquations<DIM, NUM_VARS>()
		{
			//set the integration scheme to one with integral points aligned with the nodes
			this->set_integration_scheme(new GaussWithNodes<DIM, NNODE_1D>);
			//set the number of integral points which are not aligned with nodes
			this->ipt_not_at_nodes = this->integral_pt()->nweight() - this->nnode();
		}

		QCellInterfaceElement(const QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>& dummy){BrokenCopy::broken_copy("QCellInterfaceElement");}

		void operator=(const QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>&){BrokenCopy::broken_assign("QCellInterfaceElement");}

		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			CellInterfaceEquations<DIM, NUM_VARS>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			CellInterfaceEquations<DIM, NUM_VARS>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			CellInterfaceEquations<DIM, NUM_VARS>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	CellInterfaceEquations<DIM, NUM_VARS>::output(file_pt, n_plot);
		}
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
		public virtual CellInterfaceEquations<DIM, NUM_VARS>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		TCellInterfaceElement()	:	TElement<DIM, NNODE_1D>(),
									CellInterfaceEquations<DIM, NUM_VARS>()
		{
			//set the integration scheme to one with integral points aligned with the nodes
			this->set_integration_scheme(new GaussWithNodes<DIM, NNODE_1D>);
			//set the number of integral points which are not aligned with nodes
			this->ipt_not_at_nodes = this->integral_pt()->nweight() - this->nnode();
		}

		TCellInterfaceElement(const TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>& dummy){BrokenCopy::broken_copy("TCellInterfaceElement");}

		void operator=(const TCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>&){BrokenCopy::broken_assign("TCellInterfaceElement");}

		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			CellInterfaceEquations<DIM, NUM_VARS>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			CellInterfaceEquations<DIM, NUM_VARS>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			CellInterfaceEquations<DIM, NUM_VARS>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	CellInterfaceEquations<DIM, NUM_VARS>::output(file_pt, n_plot);
		}
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
		public virtual CellInterfaceEquations<DIM, NUM_VARS>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		PointCellInterfaceElement()	:	PointElement(),
										CellInterfaceEquations<DIM, NUM_VARS>()
		{
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
			CellInterfaceEquations<DIM, NUM_VARS>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			CellInterfaceEquations<DIM, NUM_VARS>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			CellInterfaceEquations<DIM, NUM_VARS>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	CellInterfaceEquations<DIM, NUM_VARS>::output(file_pt, n_plot);
		}
	};
}

#endif