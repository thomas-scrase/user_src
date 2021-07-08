#ifndef CONDUCTING_CELL_SOLVES_FOR_VM_HEADER
#define CONDUCTING_CELL_SOLVES_FOR_VM_HEADER


#include "conducting_cell_elements.h"

namespace oomph{

template <class CELL_MODEL, class CONDUCTANCE_MODEL>
	class ConductingCellEquationsSolvesforVm : public virtual FiniteElement,
												public CELL_MODEL,
												public virtual CONDUCTANCE_MODEL,
												public virtual ConductingCellFunctionsBase
	{
	public:

		// \short function pointer to boundary source function fct(bounds, f(bounds)) --
		// bounds_of_node is a vector of the bounds the node exists on
		typedef void (* CellInterfaceBoundarySourceFctPt)
		(std::set<unsigned>* &boundaries_pt, double& bound_source);

		ConductingCellEquationsSolvesforVm() : CELL_MODEL(), CONDUCTANCE_MODEL(),
									Interpolated_Vm_Solver_Flag(Implicit),
									Cell_Variables_Are_Intentionally_Blanket_Pinned(false)
		{
			//Initialise the vector of nodes which are computed by this element
			//	by default we compute them all.
			Cell_Inds_To_Compute.resize(this->nnode(), 1);

			//Resize the cell type data suitably
			Cell_Type_Data.resize(this->nnode(), -1);

			BoostExplicitTimestepMethodPt = &BoostSolve;

			Other_Nodal_Parameters.resize(this->nnode());
			//unfortunately since the node positional data is not set up at this point we can't generate the containers, but we can allocate enough memory to hold them
			SolverContainers.resize(this->nnode());

			for(unsigned l=0; l < this->nnode(); l++){
				//Resize the nodal parameters suitably
				Other_Nodal_Parameters[l].resize(CELL_MODEL::Num_Other_Params, 0.0);

				//Create the solver containers
				SolverContainers[l] = new CellSourcesPackagedWithLocationData(dynamic_cast<ConductingCellFunctionsBase*>(this), dynamic_cast<CellModelBaseUpdated*>(this));
				SolverContainers[l]->set_timestepper_solves_for_vm();
			}
		}

		~ConductingCellEquationsSolvesforVm()
		{
			for(unsigned i=0; i < this->nnode(); i++){
				delete SolverContainers[i];
				SolverContainers[i] = 0;
			}
		}

		ConductingCellEquationsSolvesforVm(const ConductingCellEquationsSolvesforVm& dummy){BrokenCopy::broken_copy("ConductingCellEquationsSolvesforVm");}

		void operator=(const ConductingCellEquationsSolvesforVm&){BrokenCopy::broken_assign("ConductingCellEquationsSolvesforVm");}

		//How much data do we need to store in this element
		//In the order in which they are stored - conduction variables, forcing term away from decoupled solution due to diffusion, cell variables
		inline unsigned required_nvalue(const unsigned &n) const {
			return (CONDUCTANCE_MODEL::required_nvalue(n) + CELL_MODEL::Num_Cell_Vars + CELL_MODEL::Num_Output_Data + 1 + CELL_MODEL::Num_Cell_Vars);
		}



		//Indexes of variables storage in this element

		//Cell variables
		inline unsigned min_cell_variable_index_ConductingCellEquationsSolvesforVm() const {return CONDUCTANCE_MODEL::max_index_plus_one_BaseCellMembranePotential();}
		inline unsigned max_cell_variable_index_plus_one_ConductingCellEquationsSolvesforVm() const {return min_cell_variable_index_ConductingCellEquationsSolvesforVm() + CELL_MODEL::Num_Cell_Vars;}

		//Cell model output for writing to file AND communicating with other elements, i.e. active strain
		inline unsigned min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() const {return max_cell_variable_index_plus_one_ConductingCellEquationsSolvesforVm();}
		inline unsigned max_cell_model_output_data_index_plus_one_ConductingCellEquationsSolvesforVm() const {return min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() + CELL_MODEL::Num_Output_Data;}

		//vm predictd by decoupled cell solve from model
		inline unsigned cell_predicted_vm_index_ConductingCellEquationsSolvesforVm() const {return max_cell_model_output_data_index_plus_one_ConductingCellEquationsSolvesforVm();}

		//Derivatives of cell variables according to cell model, so that they can be calculated once in mpi and copied to all other procs
		inline unsigned min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() const {return cell_predicted_vm_index_ConductingCellEquationsSolvesforVm() + 1;}
		inline unsigned max_cell_variable_derivatives_index_plus_one_ConductingCellEquationsSolvesforVm() const {return min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + CELL_MODEL::Num_Cell_Vars;}

		//End indexes of variables storage in this element
		


		/////////////////////////////////////////////////////////////////////////////////
		//Get node-wise and interpolated variables
		/////////////////////////////////////////////////////////////////////////////////
		//Cell variables
		//Nodal, only
		inline double get_nodal_cell_variable(const unsigned &l, const unsigned &v) const {
			return get_nodal_cell_variable(0, l, v);
		}
		inline double get_nodal_cell_variable(const unsigned &t, const unsigned &l, const unsigned &v) const {
			return CONDUCTANCE_MODEL::node_pt(l)->value(t,min_cell_variable_index_ConductingCellEquationsSolvesforVm() + v);
		}
		inline void get_nodal_cell_variables(const unsigned &t, const unsigned &l, Vector<double> &v) const {
			for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
				v[i] = get_nodal_cell_variable(t, l, i);
			}
			if(VectorHelpers::magnitude(v)<1e-12){
				throw OomphLibError("vars are all zero.",
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
		}
		double get_dnodal_cell_variable_dt(const unsigned &l, const unsigned& i) const
		{
			// Get the data's timestepper
			TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();

			//Initialise dudt
			double dvardt=0.0;
			//Loop over the timesteps, if there is a non Steady timestepper
			if (!time_stepper_pt->is_steady())
			{
				// Number of timsteps (past & present)
				const unsigned n_time = time_stepper_pt->ntstorage();

				for(unsigned t=0;t<n_time;t++)
				{
					dvardt += time_stepper_pt->weight(1,t)*get_nodal_cell_variable(t,l,i);
				}
			}
			return dvardt;
		}

		//Output data
		//nodal
		inline double get_nodal_output_variable(const unsigned &l, const unsigned &v) const {
			return get_nodal_output_variable(0, l, v);
		}
		inline double get_nodal_output_variable(const unsigned &t, const unsigned &l, const unsigned &v) const {
			return CONDUCTANCE_MODEL::node_pt(l)->value(t,min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() + v);
		}
		inline void get_nodal_output_variables(const unsigned &t, const unsigned &l, Vector<double> &v) const {
			for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
				v[i] = get_nodal_output_variable(t, l, i);
			}
		}
		double get_dnodal_output_variable_dt(const unsigned& i, const unsigned &l) const
		{
			// Get the data's timestepper
			TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();

			//Initialise dudt
			double dvardt=0.0;
			//Loop over the timesteps, if there is a non Steady timestepper
			if (!time_stepper_pt->is_steady())
			{
				// Number of timsteps (past & present)
				const unsigned n_time = time_stepper_pt->ntstorage();

				for(unsigned t=0;t<n_time;t++)
				{
					dvardt += time_stepper_pt->weight(1,t)*get_nodal_output_variable(t,l,i);
				}
			}
			return dvardt;
		}
		inline double get_nodal_active_strain(const unsigned &t, const unsigned &l) const {
			return get_nodal_output_variable(t, l, CELL_MODEL::get_index_of_active_strain());
		}
		inline double get_interpolated_active_strain(const Vector<double> &s){
			return get_interpolated_active_strain(0, s);	
		}
		inline double get_interpolated_active_strain(const unsigned &t, const Vector<double> &s){
			//Find number of nodes
			unsigned n_node = nnode();
			//Local shape function
			Shape psi(n_node);
			//Find values of shape function
			shape(s,psi);
			double interpolated_var = 0.0;
			//Loop over the local nodes and sum
			for(unsigned l=0;l<n_node;l++)
			{
				interpolated_var += get_nodal_active_strain(t,l)*psi(l);
			}
			return interpolated_var;
		}

		//interpolatd
		inline double get_interpolated_output_variable(const Vector<double> &s, const unsigned &v) const {
			return get_interpolated_output_variable(0, s, v);
		}
		inline double get_interpolated_output_variable(const unsigned &t, const Vector<double> &s, const unsigned &v) const {
			//Find number of nodes
			unsigned n_node = nnode();
			//Local shape function
			Shape psi(n_node);
			//Find values of shape function
			shape(s,psi);
			double interpolated_var = 0.0;
			//Loop over the local nodes and sum
			for(unsigned l=0;l<n_node;l++)
			{
				interpolated_var += get_nodal_output_variable(t,l,v)*psi(l);
			}
			return interpolated_var;
		}
		inline void get_interpolated_output_variables(const unsigned &t, const Vector<double> &s, Vector<double> &v)
		{
			for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
				v[i] = get_interpolated_output_variable(t, s, i);
			}
		}
		double get_dinterpolated_output_variable_dt(const unsigned& i, const unsigned &l) const
		{
			//Find number of nodes
			unsigned n_node = nnode();
			//Local shape function
			Shape psi(n_node);
			
			Vector<double> s(CONDUCTANCE_MODEL::dim());
			//Find values of shape function
			shape(s,psi);
			double interpolated_var = 0.0;
			//Loop over the local nodes and sum
			for(unsigned l=0;l<n_node;l++)
			{
				interpolated_var += get_dnodal_output_variable_dt(i,l)*psi(l);
			}
			return interpolated_var;
		}


		//predicted vm
		//Nodal
		inline double get_nodal_predicted_vm(const unsigned &l) const {
			return get_nodal_predicted_vm(0, l);
		}
		inline double get_nodal_predicted_vm(const unsigned &t, const unsigned &l) const {
			return CONDUCTANCE_MODEL::node_pt(l)->value(t,cell_predicted_vm_index_ConductingCellEquationsSolvesforVm());
		}
		double get_dnodal_predicted_vm_dt(const unsigned& i, const unsigned &l) const
		{
			// Get the data's timestepper
			TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();

			//Initialise dudt
			double dvardt=0.0;
			//Loop over the timesteps, if there is a non Steady timestepper
			if (!time_stepper_pt->is_steady())
			{
				//Find the index at which the variable is stored
				const unsigned var_nodal_index = cell_predicted_vm_index_ConductingCellEquationsSolvesforVm();

				// Number of timsteps (past & present)
				const unsigned n_time = time_stepper_pt->ntstorage();

				for(unsigned t=0;t<n_time;t++)
				{
					dvardt += time_stepper_pt->weight(1,t)*nodal_value(t,l,var_nodal_index);
				}
			}
			return dvardt;
		}
		//Interpolated
		inline double get_interpolated_predicted_vm(const Vector<double> &s) const {
			return get_interpolated_predicted_vm(0, s);
		}
		inline double get_interpolated_predicted_vm(const unsigned &t, const Vector<double> &s) const {
			//Find number of nodes
			unsigned n_node = nnode();
			//Local shape function
			Shape psi(n_node);
			//Find values of shape function
			shape(s,psi);
			double interpolated_var = 0.0;
			//Loop over the local nodes and sum
			for(unsigned l=0;l<n_node;l++)
			{
				interpolated_var += get_nodal_predicted_vm(t,l)*psi(l);
			}

			// oomph_info << interpolated_var << std::endl;
			return interpolated_var;
		}
		double get_dinterpolated_predicted_vm_dt(const unsigned& i, const unsigned &l) const
		{
			//Find number of nodes
			unsigned n_node = nnode();
			//Local shape function
			Shape psi(n_node);

			Vector<double> s(CONDUCTANCE_MODEL::dim());
			//Find values of shape function
			shape(s,psi);
			double interpolated_var = 0.0;
			//Loop over the local nodes and sum
			for(unsigned l=0;l<n_node;l++)
			{
				interpolated_var += get_dnodal_predicted_vm_dt(i,l)*psi(l);
			}
			return interpolated_var;
		}

		//Cell variable derivatives as from cell model
		//Nodal, only
		inline double get_nodal_cell_variable_derivative(const unsigned &l, const unsigned &v) const {
			return get_nodal_cell_variable_derivative(0, l, v);
		}
		inline double get_nodal_cell_variable_derivative(const unsigned &t, const unsigned &l, const unsigned &v) const {
			return CONDUCTANCE_MODEL::node_pt(l)->value(t,min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + v);
		}
		inline void get_nodal_cell_variable_derivatives(const unsigned &t, const unsigned &l, Vector<double> &v) const {
			for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
				v[i] = get_nodal_cell_variable_derivative(t, l, i);
			}
		}


		///End get nodal values of variables stored by this element



		/////////////////////////////////////////////////////////////////////////////////
		//Override the conducting element to take transmembrane current from this element
		/////////////////////////////////////////////////////////////////////////////////
		inline double get_nodal_predicted_vm_BaseCellMembranePotential(const unsigned& l) const
		{
			return this->get_nodal_predicted_vm(l);
		}

		//Get the membrane potential at the node l, with suitable solver method
		inline double get_nodal_vm_Segregated(const unsigned& l) const
		{
			switch(Interpolated_Vm_Solver_Flag){
				case Implicit :
					return CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(l);
				case Explicit :
					return CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(1, l);
				case CrankNicolson : 
					return 0.5*(CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(l)
								+ CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(1,l));
			}

			throw OomphLibError("that is not a solver flag",
	                       	OOMPH_CURRENT_FUNCTION,
	                       	OOMPH_EXCEPTION_LOCATION);
			//To make compiler shut up
			return 0.0;
		}

		inline double get_nodal_cell_variable_derivative_Segregated(const unsigned& l, const unsigned& v){
			switch(Interpolated_Vm_Solver_Flag){
				case Implicit :
					return get_nodal_cell_variable_derivative(l, v);
				case Explicit :
					return get_nodal_cell_variable_derivative(1, l, v);
				case CrankNicolson : 
					return 0.5*(get_nodal_cell_variable_derivative(l, v)
								+ get_nodal_cell_variable_derivative(1, l, v));
			}

			throw OomphLibError("that is not a solver flag",
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			//To make compiler shut up
			return 0.0;
		}


		/////////////////////////////////////////////////////////////////////////////////
		//Abuse of the oomph lib integral schemes to allow us to grab data from external
		//elements at the location of the nodes of this element
		/////////////////////////////////////////////////////////////////////////////////
		//return the Gauss point associated with node n
		inline unsigned ipt_at_node(const unsigned &n) const
			{return ipt_not_at_nodes + n;}



		/////////////////////////////////////////////////////////////////////////////////
		//For distributing the cell timestep solve properly over the elements in a mesh
		/////////////////////////////////////////////////////////////////////////////////
		unsigned n_computed_node(){
			unsigned temp = 0;
			for(unsigned i=0; i<this->nnode(); i++){
				temp += Cell_Inds_To_Compute[i];
			}
			return temp;
		}
		void do_not_compute_node(const unsigned &n){
			#ifdef PARANOID
			if(n>this->nnode()){
				throw OomphLibError("n > nnode",
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			#endif
			Cell_Inds_To_Compute[n] = 0;
		}
		void do_compute_node(const unsigned &n){
			#ifdef PARANOID
			if(n>this->nnode()){
				throw OomphLibError("n > nnode",
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			#endif
			Cell_Inds_To_Compute[n] = 1;
		}
		unsigned is_node_computed(const unsigned &n) const {
			#ifdef PARANOID
			if(n>this->nnode()){
				throw OomphLibError("n > nnode",
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			#endif
			return Cell_Inds_To_Compute[n];
		}


		


		/////////////////////////////////////////////////////////////////////////////////
		//Get and Set cell type at node n
		/////////////////////////////////////////////////////////////////////////////////

		//Cell type at node
		//	Used in switch function with the following correspondence:
		//		ATRIA 0 - 99,
		//		VENTS 100 - 199,
		//		OTHER 200 - 299 (?)
		//		CNZAtria
		//		0 RA, 1 PM, 2 CT, 3 RAA, 4 AS, 5 AVR, 6 BB, 7 LA, 8 LAA, 9 PV, 10 SAN_C, 11 SAN_P
		//		TNNPVentricle
		//		100 LVEPI 101 LVMCELL 102 LVENDO 103 RVEPI 104 RVMCELL 105 RVENDO 106 PFI 107 PFMB 108 PF
		void set_cell_type(const unsigned &l, const unsigned &cell_type){
			Cell_Type_Data[l] = cell_type;

			SolverContainers[l]->set_cell_type(Cell_Type_Data[l]);
		}


		/////////////////////////////////////////////////////////////////////////////////
		//Set other node parameters at node l,
		//	these could be, for example:
		//	 apico-basal ratio,
		//	 other heterogeneity factors,
		//	 ion channel block factors,
		//	Parameters which do not change and are assigned cell (node) wise
		/////////////////////////////////////////////////////////////////////////////////
		//set the var-th black box nodal parameter associated with the l-th node to value
		inline void set_Other_Nodal_Parameter(const unsigned &l, const unsigned &var, const double &value){
			Other_Nodal_Parameters[l][var] = value;

			SolverContainers[l]->set_other_nodal_parameters(Other_Nodal_Parameters[l]);
		}



		void probe_eqn_numbers(){
			std::ofstream outfile;
  			outfile.open("MULTI/probe.dat", std::ios_base::app);

			for(unsigned n=0; n<this->nnode(); n++){
				outfile << this->node_pt(n)->nvalue() << std::endl;
				for(unsigned i=min_cell_variable_index_ConductingCellEquationsSolvesforVm(); i<min_cell_variable_index_ConductingCellEquationsSolvesforVm() + CELL_MODEL::Num_Cell_Vars; i++){
					// const unsigned local_eqn = this->nodal_local_eqn(n,min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i);
					int eqn_number = this->node_pt(n)->eqn_number(i);
					if(eqn_number>=0){
						// outfile << this << " " << this->node_pt(n) << " " << n << " cell variable " << i << " is unpinned" << std::endl;
						outfile << "Eqn: "<<eqn_number<< " | is cell variable " << i << std::endl;
					}
					// else{
					// 	outfile << this << " " << this->node_pt(n) << " " << n << " cell variable " << i << " is pinned" << std::endl;	
					// }
				}

				for(unsigned i=min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm(); i<min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + CELL_MODEL::Num_Cell_Vars; i++){
					// const unsigned local_eqn = this->nodal_local_eqn(n,min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + i);
					int eqn_number = this->node_pt(n)->eqn_number(i);
					if(eqn_number>=0){
						// outfile << this << " " << this->node_pt(n) << " " << n << " cell variable deriv " << i << " is unpinned" << std::endl;
						outfile << "Eqn: "<<eqn_number<< " | is cell variable derivative " << i << std::endl;
					}
					// else{
					// 	outfile << this << " " << this->node_pt(n) << " " << n << " cell variable deriv " << i << " is pinned" << std::endl;	
					// }
				}
				


			}
		}
		/////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////
		//Assign initial conditions at nodes from cell model
		/////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////
		inline void assign_initial_conditions(){
			for(unsigned n=0; n<this->nnode(); n++){
				unsigned counter = 0;
				//Assign initial conditions of the membrane potential

				//Set the initial value of the membrane potential
				this->node_pt(n)->set_value(CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential(),
												CELL_MODEL::return_initial_membrane_potential(Cell_Type_Data[n]));
				
				//How to handle bidomain? it requires an extra variable to be assigned
				//Assign intitial conditions of cell model variables
				for(unsigned v=0; v<CELL_MODEL::Num_Cell_Vars; v++){
					this->node_pt(n)->set_value(min_cell_variable_index_ConductingCellEquationsSolvesforVm() + v,
												CELL_MODEL::return_initial_state_variable(v,Cell_Type_Data[n]));
					if(this->node_pt(n)->value(min_cell_variable_index_ConductingCellEquationsSolvesforVm() + v)>1e-12){
						counter++;
					}
				}

				for(unsigned v=0; v<CELL_MODEL::Num_Output_Data; v++){
					this->node_pt(n)->set_value(min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() + v, 0.0);
				}

				//Set the initial value of the membrane potential
				this->node_pt(n)->set_value(cell_predicted_vm_index_ConductingCellEquationsSolvesforVm(), this->node_pt(n)->value(CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential()));

				for(unsigned v=0; v<CELL_MODEL::Num_Cell_Vars; v++){
					this->node_pt(n)->set_value(min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + v, 0.0);
				}

				if(counter==0){
					oomph_info << "Node " << n << " cell vars were not initialized." << std::endl;
				}
			}

			//Let the conductance model assign any other initial conditions it needs.
			//	The monodomain model does not do anything, however the bidomain model
			//	Needs to be able to set initial values for extracellular membrane potential
			CONDUCTANCE_MODEL::assign_additional_initial_conditions();
		}


		/////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////
		//Access the explict timestepping method used by the solver, used to set the pointer
		/////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////
		TomsExplicitTimestepMethodsFctPt &Boost_Explicit_Timestep_Method_Pt()
		{
			return BoostExplicitTimestepMethodPt;
		}




		/////////////////////////////////////////////////////////////////////////////////
		//Pinning and unpinning of specific groups of local variables for the purpose
		//	of segregated solving.
		/////////////////////////////////////////////////////////////////////////////////

		//Vm
		void pin_all_vm(){
			//Ask that the underlying conduction element pins all of it's dofs, note this is not just vm but any other variables required, e.g. bidomain variables
			CONDUCTANCE_MODEL::pin_all_vars();
		}
		void unpin_all_vm(){
			//Ask that the underlying conduction element pins all of it's dofs, note this is not just vm but any other variables required, e.g. bidomain variables
			CONDUCTANCE_MODEL::unpin_all_vars();
		}

		//Cell vars
		void pin_all_cell_vars(){
			for(unsigned l=0; l<this->nnode(); l++){
				//Pin all the cell variables and vm predicted by decoupled solver
				for(unsigned i=min_cell_variable_index_ConductingCellEquationsSolvesforVm(); i<max_cell_variable_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
					this->node_pt(l)->pin(i);
				}
			}
			Cell_Variables_Are_Intentionally_Blanket_Pinned = true;
		}
		void unpin_all_cell_vars(){
			// oomph_info << "Unpin all cell vars has been called" << std::endl;
			for(unsigned l=0; l<this->nnode(); l++){
				for(unsigned i=min_cell_variable_index_ConductingCellEquationsSolvesforVm(); i<max_cell_variable_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
					this->node_pt(l)->unpin(i);
				}
			}
			Cell_Variables_Are_Intentionally_Blanket_Pinned = false;
		}

		//Output vars
		void pin_all_cell_model_output_vars(){
			for(unsigned l=0; l<this->nnode(); l++){
				//Pin all the cell variables and vm predicted by decoupled solver
				for(unsigned i=min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm(); i<max_cell_model_output_data_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
					this->node_pt(l)->pin(i);
				}
			}
		}
		void unpin_all_cell_model_output_vars(){
			for(unsigned l=0; l<this->nnode(); l++){
				for(unsigned i=min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm(); i<max_cell_model_output_data_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
					this->node_pt(l)->unpin(i);
				}
			}
		}

		//Ionic current
		void pin_cell_model_predicted_vm(){
			for(unsigned l=0; l<this->nnode(); l++){
				this->node_pt(l)->pin(cell_predicted_vm_index_ConductingCellEquationsSolvesforVm());
				
			}
		}
		void unpin_cell_model_predicted_vm(){
			for(unsigned l=0; l<this->nnode(); l++){
				this->node_pt(l)->unpin(cell_predicted_vm_index_ConductingCellEquationsSolvesforVm());
				
			}
		}

		//Output vars
		void pin_all_cell_model_vars_derivatives(){
			for(unsigned l=0; l<this->nnode(); l++){
				//Pin all the cell variables and vm predicted by decoupled solver
				for(unsigned i=min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm(); i<max_cell_variable_derivatives_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
					this->node_pt(l)->pin(i);
				}
			}
		}
		void unpin_all_cell_model_vars_derivatives(){
			for(unsigned l=0; l<this->nnode(); l++){
				for(unsigned i=min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm(); i<max_cell_variable_derivatives_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
					this->node_pt(l)->unpin(i);
				}
			}
		}


		void pin_all_cell_related_storage(){
			pin_all_cell_vars();
			pin_all_cell_model_output_vars();
			pin_cell_model_predicted_vm();
			pin_all_cell_model_vars_derivatives();
		}

		void unpin_all_cell_related_storage(){
			unpin_all_cell_vars();
			unpin_all_cell_model_output_vars();
			unpin_cell_model_predicted_vm();
			unpin_all_cell_model_vars_derivatives();
		}


		void Extern_Has_Intentionally_Blanket_Pinned_All_Cell_Variables(){Cell_Variables_Are_Intentionally_Blanket_Pinned = true;}

		void Extern_Has_Intentionally_Blanket_UnPinned_All_Cell_Variables(){Cell_Variables_Are_Intentionally_Blanket_Pinned = false;}


		void set_interpolated_vm_solver_flag(const unsigned& flag){
			switch(flag){
				case Implicit:
					Interpolated_Vm_Solver_Flag = flag;
					return;
				case Explicit:
					Interpolated_Vm_Solver_Flag = flag;
					return;
				case CrankNicolson:
					Interpolated_Vm_Solver_Flag = flag;
					return;
			};
		}

		//This is the decoupled solver, assuming cells are all fully decoupled from the diffusion problem. It takes a timestep
		// from the previous time values of cell variables and vm and returns predicted values for the current timestep.
		//If the user sets the segregated solve flag then the calculated value of vm is inserted into the residual for 
		//	vm predicted by the decoupled solve, otherwise it is just put into that of vm
		//THIS of course assumes that the user has properly set the segregated solve flag. Segregated solve flag should only
		// be set by the segregated solve problem class since it then properly pins redundant dofs and rebuilds the global mesh
		// accordingly ensuring that
		void perform_decoupled_solve(const double& dt, Vector<double>& residuals, const bool &use_current_as_initial){
			int local_eqn = 0;

			//Allocate memory for the New Variables
			Vector<double> New_Variables(CELL_MODEL::Num_Cell_Vars, 0.0);
			//Allocate New_Vm
			double New_Vm = 0.0;

			//loop over the cells
			for(unsigned l=0; l<this->nnode(); l++){
				// oomph_info << "Node " << l << " is pinned? " << this->is_node_computed(l) << std::endl;
				//if the cell is not to be computed then don't compute it
				if(this->is_node_computed(l)){
					// oomph_info << "Solving node " << l << " of " << this->nnode() << std::endl;
					//Zero the results of the derivatives calculation
					std::fill(New_Variables.begin(), New_Variables.end(), 0.0); 
					New_Vm = 0.0;

					//The variables from which we perform the timestepping
					//Get vm
					double Vm;// = get_nodal_vm_Segregated(l);
					if(use_current_as_initial)
					{
						Vm = CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(0, l);
					}
					else
					{
						Vm = CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(1, l);
					}

					//Get t
					double t;
					if(use_current_as_initial)
					{
						t = this->node_pt(l)->time_stepper_pt()->time_pt()->time();
					}
					else
					{
						t = this->node_pt(l)->time_stepper_pt()->time_pt()->time(1);
					}					

					//Get cell variables
					Vector<double> CellVariables(CELL_MODEL::Num_Cell_Vars, 0.0);
					if(use_current_as_initial)
					{
						get_nodal_cell_variables(0, l , CellVariables);
					}
					else
					{
						get_nodal_cell_variables(1, l , CellVariables);
					}
					
					//Get Cell type

					//Set the location data of the solver container: ideally this wouldn't have to be performed every time we take a timestep but I'm not sure where
					//	else to put it

					//Get coordinate of the node
					unsigned ipt = ipt_at_node(l);
					//Get the local coordinate in the element
					Vector<double> s(CONDUCTANCE_MODEL::dim());
	 				//Get the global coordinate
	 				Vector<double> x(CONDUCTANCE_MODEL::dim());

					for(unsigned j=0;j<CONDUCTANCE_MODEL::dim();j++)
					{
						s[j] = this->integral_pt()->knot(ipt,j);
						x[j] = this->raw_nodal_position(l,j);
					}
					SolverContainers[l]->set_location_data(ipt, s, x, l);				

					//Perform the explicit timestep
					(*BoostExplicitTimestepMethodPt)(Vm,
													CellVariables,
													t,
													dt,
													*(SolverContainers[l]),

													New_Variables,
													New_Vm);

					//update the nodal variables and the predicted membrane potential
					for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
						local_eqn = this->nodal_local_eqn(l,min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i);
						if(local_eqn>=0){
							residuals[local_eqn] = New_Variables[i];
						}
					}

					//Update the vm predicted by the decoupled cell solve
					local_eqn = this->nodal_local_eqn(l, cell_predicted_vm_index_ConductingCellEquationsSolvesforVm());
					if(local_eqn>=0){
						residuals[local_eqn] = New_Vm;
					}

					//Update the vm predicted by the decoupled cell solve
					local_eqn = this->nodal_local_eqn(l, CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential());
					if(local_eqn>=0){
						// oomph_info << "Predvm " << New_Vm << std::endl;
						residuals[local_eqn] = New_Vm;
					}
				}
			}
		}

		//Get the new cell variables from the residuals vector - 
		//Residuals vector contains ALL dofs - therefore we need to pull out the global eqn number index
		void update_cell_values_from_assembled_vector(const Vector<double>& residuals)
		{
			int local_eqn = 0;
			//loop over the cells
			for(unsigned l=0; l<this->nnode(); l++){
				// oomph_info << "Residuals pertaining to node " << l << std::endl;
				unsigned vars_updated = 0;
				double magnitude_of_update = 0;
				//update the nodal variables and the predicted membrane potential
				for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
					local_eqn = this->nodal_local_eqn(l,min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i);
					if(local_eqn>=0){
						vars_updated++;
						magnitude_of_update += residuals[this->eqn_number(local_eqn)];
						// oomph_info << this->eqn_number(local_eqn) << " " << residuals[this->eqn_number(local_eqn)] << std::endl;
						this->node_pt(l)->set_value(min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i, residuals[this->eqn_number(local_eqn)]);
					}
				}

				local_eqn = this->nodal_local_eqn(l,cell_predicted_vm_index_ConductingCellEquationsSolvesforVm());
				if(local_eqn>=0){
					// oomph_info << this->eqn_number(local_eqn) << " " << New_Vm << std::endl;
					this->node_pt(l)->set_value(cell_predicted_vm_index_ConductingCellEquationsSolvesforVm(), residuals[this->eqn_number(local_eqn)]);
				}

				//Update the vm predicted by the decoupled cell solve
				local_eqn = this->nodal_local_eqn(l, CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential());
				if(local_eqn>=0){
					// oomph_info << "Predvm " << New_Vm << std::endl;
					// residuals[local_eqn] = New_Vm;
					this->node_pt(l)->set_value(CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential(), residuals[this->eqn_number(local_eqn)]);
				}

				// oomph_info << "Vars updated on node " << l << ": " << vars_updated << ", with magnitude " << magnitude_of_update << std::endl;


				// Vector<double> Vars(CELL_MODEL::Num_Cell_Vars, 0.0);
				// get_nodal_cell_variables(0, l, Vars);


			}
		}

		//Take a timestep using the segregated solve from the previous history values of vm and cell variables
		// and assign the new values to vm and cell variables
		void explicit_decoupled_timestep(const double &dt){
			Vector<double> residuals(this->ndof(), 0.0);

			//Perform the decoupled timestep
			perform_decoupled_solve(dt, residuals);

			unsigned n_node = this->nnode();

			int local_eqn=0;
			//loop over the cells
			for(unsigned l=0; l<n_node; l++){

				//if the cell is not to be computed then don't compute it
				if(this->is_node_computed(l)){
					for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
						local_eqn = this->nodal_local_eqn(l,min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i);
						if(local_eqn>=0){
							this->node_pt(l)->set_value(min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i, residuals[local_eqn]);
						}
					}

					//Set the new value of the actual membrane potential if it is not pinned
					local_eqn = this->nodal_local_eqn(l,CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential());
					if(local_eqn>=0){
						this->node_pt(l)->set_value(CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential(), residuals[local_eqn]);
					}
				}
			}
		}



		//Get the current values of the derivatives, iion, and output data from the cell model
		// The results are placed into a residual vector so that this function can be reused
		//  in mpi problems for efficient calculation of these values
		void get_data_from_cell_model(Vector<double>& residuals)
		{
			int local_eqn = 0;

			unsigned n_node = this->nnode();

			double Vm = 0.0;

			Vector<double> CellVariables(CELL_MODEL::Num_Cell_Vars, 0.0);

			//The current nodal time
			double t = 0.0;
			//Get coordinate of the node
			unsigned ipt = 0;
			//Get the local coordinate in the element
			Vector<double> s(CONDUCTANCE_MODEL::dim());
			//Get the global coordinate
			Vector<double> x(CONDUCTANCE_MODEL::dim());
			//Other nodal variables - the time dependent ones
			Vector<double> OtherVariables(CELL_MODEL::Num_Other_Vars, 0.0);

			//What we are solving for
			Vector<double> Variable_Derivatives(CELL_MODEL::Num_Cell_Vars, 0.0);
			double Iion = 0.0;
			Vector<double> Cell_Output(CELL_MODEL::Num_Output_Data, 0.0);

			for(unsigned l=0;l<n_node;l++){
				// //if the cell is not to be computed then don't compute it
				if(this->is_node_computed(l)){
					//Get node-specific coordinates
					t = this->node_pt(l)->time_stepper_pt()->time_pt()->dt();
					
					ipt = ipt_at_node(l);

					for(unsigned j=0;j<CONDUCTANCE_MODEL::dim();j++)
					{
						s[j] = this->integral_pt()->knot(ipt,j);
						x[j] = this->raw_nodal_position(l,j);
					}
					//Get vm at the l-th node
					Vm = CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(l);

					//Get the cell variables
					get_nodal_cell_variables(0, l , CellVariables);

					//Get time dependent variables
					get_other_variables(ipt, s, x, l, t, OtherVariables);
					
					//Calculate derivatives and Iion
					CELL_MODEL::Calculate_Derivatives(Vm,
													CellVariables,
													t,
													Cell_Type_Data[l],
													get_stimulus(ipt, s, x, t),
													Other_Nodal_Parameters[l],
													OtherVariables,
													Variable_Derivatives,
													Iion);

					//Calculate ouptut data
					CELL_MODEL::get_optional_output(Vm,
													CellVariables,
													t,
													Cell_Type_Data[l],
													get_stimulus(ipt, s, x, t),
													Other_Nodal_Parameters[l],
													OtherVariables,
													Cell_Output);

					
					//Suitably fill in residual vector
					for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
						local_eqn = this->nodal_local_eqn(l,min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + i);
						if(local_eqn>=0){
							residuals[local_eqn] = Variable_Derivatives[i];
						}
					}
					
					for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
						local_eqn = this->nodal_local_eqn(l,min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() + i);
						if(local_eqn>=0){
							residuals[local_eqn] = Cell_Output[i];
						}
					}
				}
			}
		}

		//Residuals vector contains ALL dofs - therefore we need to pull out the global eqn number index
		void update_cell_model_data_from_assembled_vector(const Vector<double>& residuals)
		{
			int local_eqn = 0;

			unsigned ipt = 0;
			//loop over the cells
			for(unsigned l=0; l<this->nnode(); l++){
				//Update output data
				for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
					local_eqn = this->nodal_local_eqn(l,min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() + i);
					if(local_eqn>=0){
						this->node_pt(l)->set_value(min_cell_model_output_data_index_ConductingCellEquationsSolvesforVm() + i, residuals[this->eqn_number(local_eqn)]);
					}
				}
				
				//update the nodal variable derivative data
				for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
					local_eqn = this->nodal_local_eqn(l,min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + i);
					if(local_eqn>=0){
						this->node_pt(l)->set_value(min_cell_variable_derivatives_index_ConductingCellEquationsSolvesforVm() + i, residuals[this->eqn_number(local_eqn)]);
					}
				}
			}
		}

		//Assemble the residuals and jacobians for the cell variables using the oomph-lib machinery
		void fill_in_generic_residual_contribution_ConductingCellElements(
			Vector<double> &residuals, DenseMatrix<double> &jacobian, 
		    DenseMatrix<double> &mass_matrix, unsigned flag)
		{
			unsigned n_node = this->nnode();

			int local_eqn = 0;

			for(unsigned l=0;l<n_node;l++){
				//Cell variable residuals - these are all just node-wise and use the stored values of cell variable derivatives
				for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
					local_eqn = this->nodal_local_eqn(l,min_cell_variable_index_ConductingCellEquationsSolvesforVm() + i);
					if(local_eqn>=0){
						residuals[local_eqn] -= (get_dnodal_cell_variable_dt(l, i) - get_nodal_cell_variable_derivative_Segregated(l, i));
					}
				}//End fill in cell variable residuals
			}
			
		}//End get residuals





		/////////////////////////////////////////////////////////////////////////////////
		//Residual and Jacobian functions
		// The conductance model is the only thing that does actual oomph lib solving
		/////////////////////////////////////////////////////////////////////////////////
		/// Add the element's contribution to its residual vector (wrapper)
		void fill_in_contribution_to_residuals(Vector<double> &residuals)
		{	
			// for(unsigned l=0;l<this->nnode();l++){
			// 	for(unsigned i=min_cell_variable_index_ConductingCellEquationsSolvesforVm(); i<max_cell_variable_derivatives_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
			// 		if(this->nodal_local_eqn(l, i)>=0){
			// 			oomph_info << "local variable " << i << " is unpinned going into fill in residuals" << std::endl;
			// 		}
			// 	}
			// }

			// oomph_info << "Cell_Variables_Are_Intentionally_Blanket_Pinned" << std::endl;
			CONDUCTANCE_MODEL::fill_in_generic_residual_contribution_BaseCellMembranePotential(residuals, GeneralisedElement::Dummy_matrix,
																				GeneralisedElement::Dummy_matrix, 0);
			//If we have blanket pinned all cell variables then don't bother filling in the jacobian for them
			if(!Cell_Variables_Are_Intentionally_Blanket_Pinned){
				// oomph_info << "Filling in cell variable residuals with oomph lib machinery" << std::endl;
				fill_in_generic_residual_contribution_ConductingCellElements(residuals, GeneralisedElement::Dummy_matrix,
																	GeneralisedElement::Dummy_matrix, 0);
				
			}			
		}
		
		/// \short Add the element's contribution to its residual vector and 
		/// the element Jacobian matrix (wrapper)
		void fill_in_contribution_to_jacobian(Vector<double> &residuals,
		                                   DenseMatrix<double> &jacobian)
		{	
			// for(unsigned l=0;l<this->nnode();l++){
			// 	for(unsigned i=min_cell_variable_index_ConductingCellEquationsSolvesforVm(); i<max_cell_variable_derivatives_index_plus_one_ConductingCellEquationsSolvesforVm(); i++){
			// 		if(this->nodal_local_eqn(l, i)>=0){
			// 			oomph_info << "local variable " << i << " is unpinned going into fill in jacobian" << std::endl;
			// 		}
			// 	}
			// }

			CONDUCTANCE_MODEL::fill_in_contribution_to_jacobian(residuals, jacobian);

			// //If we have blanket pinned all cell variables then don't bother filling in the jacobian for them
			// if(Cell_Variables_Are_Intentionally_Blanket_Pinned){
			// 	CONDUCTANCE_MODEL::fill_in_contribution_to_jacobian(residuals, jacobian);
			// 	// oomph_info << "Filling in conductance jacobian analytically" << std::endl;
			// }
			// else{//Otherwise we have to do the prohibitively expensive task of finite differencing them all
			// 	FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
			// 	// std::cout << "Filling in all variables, possibly including cell variables, by finite differencing" << std::endl;
			// }
		}

		/// Add the element's contribution to its residuals vector,
		/// jacobian matrix and mass matrix
		void fill_in_contribution_to_jacobian_and_mass_matrix(
											Vector<double> &residuals, DenseMatrix<double> &jacobian,
											DenseMatrix<double> &mass_matrix)
		{
			//Broken - how can this be filled in for the cell variables? - Not obvious
			FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);
		}





		/////////////////////////////////////////////////////////////////////////////////
		//Get the output data from the cell
		/////////////////////////////////////////////////////////////////////////////////
		void get_optional_cell_output(const unsigned &l, Vector<double> &Out) const
		{

			//Get the variables required by the cell model
			//get ipt
			const unsigned ipt = ipt_at_node(l);
			//Get vm
			const double Vm = CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(l);
			//Get t
			const double t = this->node_pt(l)->time_stepper_pt()->time();
			//Get dt, the model presumably only uses dt to calculate derivatives, so we'll just give it a very small value 
			// const double dt = 1e-9;
			//////////
			//Get cell variables
			Vector<double> CellVariables(CELL_MODEL::Num_Cell_Vars, 0.0);
			get_nodal_cell_variables(0, l , CellVariables);
			//Get Cell type
			const unsigned CellType = Cell_Type_Data[l];
			//Get Other Parameters
			// const Vector<double> OtherParameters = Other_Nodal_Parameters[l];
			const Vector<double> OtherParameters = Other_Nodal_Parameters[l];

			
			Vector<double> s(CONDUCTANCE_MODEL::dim());
			//Assign values of s
			for(unsigned i=0;i<CONDUCTANCE_MODEL::dim();i++){
				s[i] = this->integral_pt()->knot(ipt,i);
			}
			Vector<double> x(CONDUCTANCE_MODEL::dim());
			for(unsigned j=0;j<CONDUCTANCE_MODEL::dim();j++){
				x[j] += this->raw_nodal_position(l,j);
			}

			Vector<double> OtherVariables(CELL_MODEL::Num_Other_Vars, 0.0);
			ConductingCellFunctionsBase::get_other_variables(ipt,
															s,
															x,
															l,
															t,
															OtherVariables);

			//Call the function from the cell model
			CELL_MODEL::get_optional_output(Vm,
											CellVariables,
											t,
											CellType,
											ConductingCellFunctionsBase::get_stimulus(ipt,s,x,t),
											OtherParameters,
											OtherVariables,

											Out);

		}




		//Get the interpolated optional cell output for outputting purposes
		void get_interpolated_optional_cell_output(const Vector<double>& s, Vector<double> &Out) const
		{
			Out.resize(CELL_MODEL::Num_Output_Data, 0.0);

			const unsigned n_node = this->nnode();
			Shape psi(n_node);
			shape(s,psi);

			for(unsigned l=0; l<n_node; l++){
				Vector<double> node_out(CELL_MODEL::Num_Output_Data, 0.0);
				get_optional_cell_output(l, node_out);
				for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
					Out[i] += psi[l]*node_out[i];
				}
			}
		}
		

		/////////////////////////////////////////////////////////////////////////////////
		//Oomph lib Output functions
		/////////////////////////////////////////////////////////////////////////////////
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			unsigned nplot=5;
			output(outfile,nplot);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			//Vector of local coordinates
			Vector<double> s(CONDUCTANCE_MODEL::dim());
			// Tecplot header info
			outfile << tecplot_zone_string(2);

			//Get the shape functions
			const unsigned n_node = this->nnode();
			Shape psi(n_node);
			shape(s,psi);

			// Loop over plot points
			unsigned num_plot_points=nplot_points(2);
			for (unsigned iplot=0;iplot<num_plot_points;iplot++)
			{
				// Get local coordinates of plot point
				get_s_plot(iplot,2,s);

				// Get Eulerian coordinate of plot point
				Vector<double> x(CONDUCTANCE_MODEL::dim());
				interpolated_x(s,x);
				for(unsigned i=0;i<CONDUCTANCE_MODEL::dim();i++) {outfile << x[i] << " ";}

				//Get the membrane potential to output
				double vm_iplot = 0.0;
				//Get the cell variables to output
				Vector<double> cell_var_iplot(CELL_MODEL::Num_Cell_Vars, 0.0);
				//Get the optional cell output
				Vector<double> opt_out_iplot(CELL_MODEL::Num_Output_Data, 0.0);

				//Loop over the nodes in the element
				for(unsigned l=0; l<n_node; l++){
					vm_iplot += CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(l)*psi[l];

					Vector<double> nodal_cell_vars(CELL_MODEL::Num_Cell_Vars, 0.0);
					get_nodal_cell_variables(0, l, nodal_cell_vars);
					for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
						cell_var_iplot[i] += nodal_cell_vars[i]*psi[l];
					}

					Vector<double> nodal_output_vars(CELL_MODEL::Num_Output_Data, 0.0);
					get_nodal_output_variables(0, l, nodal_output_vars);
					for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
						opt_out_iplot[i] += nodal_output_vars[i]*psi[l];
					}
				}

				outfile << vm_iplot << " ";

				for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
					outfile << cell_var_iplot[i] << " ";
				}
				for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
					outfile << opt_out_iplot[i] << " ";
				}

				outfile << std::endl;
			}
		};

		//Get the output for iplot and put it into a vector
		void output_to_vect(Vector<double> &outvect, const unsigned &iplot){

			outvect.resize(CELL_MODEL::Num_Cell_Vars+CELL_MODEL::Num_Output_Data, 0.0);

			//Vector of local coordinates
			Vector<double> s(CONDUCTANCE_MODEL::dim());

			//Get the shape functions
			const unsigned n_node = this->nnode();
			Shape psi(n_node);
			shape(s,psi);

			// Get local coordinates of plot point
			get_s_plot(iplot,2,s);

			// Get Eulerian coordinate of plot point
			Vector<double> x(CONDUCTANCE_MODEL::dim());
			interpolated_x(s,x);
			// for(unsigned i=0;i<CONDUCTANCE_MODEL::dim();i++) {outvect[i] = x[i];}

			//Loop over the nodes in the element
			for(unsigned l=0; l<n_node; l++){
				outvect[0] += CONDUCTANCE_MODEL::get_nodal_membrane_potential_BaseCellMembranePotential(l)*psi[l];

				Vector<double> nodal_cell_vars(CELL_MODEL::Num_Cell_Vars, 0.0);
				get_nodal_cell_variables(0, l, nodal_cell_vars);
				for(unsigned i=0; i<CELL_MODEL::Num_Cell_Vars; i++){
					outvect[1+i] += nodal_cell_vars[i]*psi[l];
				}

				Vector<double> nodal_output_vars(CELL_MODEL::Num_Output_Data, 0.0);
				get_nodal_output_variables(0, l, nodal_output_vars);
				for(unsigned i=0; i<CELL_MODEL::Num_Output_Data; i++){
					outvect[1+CELL_MODEL::Num_Cell_Vars+i] += nodal_output_vars[i]*psi[l];
				}
			}
		};

		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			unsigned n_plot=5;
			output(file_pt,n_plot);
		}
		/// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// n_plot^DIM plot points
		void output(FILE* file_pt, const unsigned &n_plot){};



		
		/////////////////////////////////////////////////////////////////////////////////
		//Paraview Output
		/////////////////////////////////////////////////////////////////////////////////
		unsigned nscalar_paraview() const
		{	
			return (CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars + CELL_MODEL::Num_Output_Data);
		}

		void scalar_value_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot) const
		{
			//Vector of local coordinates
	 		Vector<double> s(CONDUCTANCE_MODEL::dim());

	 		//Get the number of nodes
	 		const unsigned n_node = this->nnode();

	 		//Preallocate the shape function
	 		Shape psi(n_node);

	 		for(unsigned l=0; l<n_node; l++){
				if(i<CONDUCTANCE_MODEL::get_variable_names().size())
				{
					file_out << this->nodal_value(l, CONDUCTANCE_MODEL::vm_index_BaseCellMembranePotential()+i) << std::endl;
					continue;
				}

				if(i<(CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars)){
					file_out << get_nodal_cell_variable(0, l, i-CONDUCTANCE_MODEL::get_variable_names().size()) << std::endl;
					continue;
				}

				if(i < (CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars + CELL_MODEL::Num_Cell_Vars)){
					Vector<double> Out(CELL_MODEL::Num_Output_Data, 1e300);
					get_optional_cell_output(l, Out);
					file_out << Out[i-(CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars)] << std::endl;
					continue;
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
			if(i<CONDUCTANCE_MODEL::get_variable_names().size())
			{
				return CONDUCTANCE_MODEL::get_variable_names()[i];
			}

			if(i<(CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars)){
				return CELL_MODEL::Names_Of_Cell_Variables[i-CONDUCTANCE_MODEL::get_variable_names().size()];
			}

			if(i < (CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars + CELL_MODEL::Num_Cell_Vars)){
				return CELL_MODEL::Names_Of_Output_Data[ i - (CONDUCTANCE_MODEL::get_variable_names().size() + CELL_MODEL::Num_Cell_Vars) ];
			}

			throw OomphLibError(
					"A variable index which was too large was requested when getting names in paraview output",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
		}




	protected:
		//The number of integral points which are not additional ones placed at the nodes
		unsigned ipt_not_at_nodes;
			



	private:

		std::vector<CellSourcesPackagedWithLocationData*> SolverContainers;

		// unsigned Cell_type_internal_index;
		Vector<unsigned> Cell_Type_Data;

		std::vector<Vector<double>> Other_Nodal_Parameters;

		// list of local node indexes for which the jacobian and
		//	residual entries are computed. The fill in is
		//	exactly the same as that of neighbouring elements
		//	so we can safely ignore repeated nodes. By default all
		//	nodes are computed. 1 indicates that a node is computed,
		//	0 indicates that it is not computed
		Vector<unsigned> Cell_Inds_To_Compute;

		//The explicit timestepping method we use to solve for the cell variables and also
		//	the spatially independent portion of the membrane potential eqauation
		TomsExplicitTimestepMethodsFctPt BoostExplicitTimestepMethodPt;


		//Flags to inform functions which time values to use for getting data,
		// they are used so that segregated solvers can correctly operate
		//By default we just use the implicit, curent timestep, method
		unsigned Interpolated_Vm_Solver_Flag; //Which timestep do we get Vm from


		bool Cell_Variables_Are_Intentionally_Blanket_Pinned;
	};





































	//====================================================================
	//====================================================================
	//Q Element
	//====================================================================
	//====================================================================
	template<unsigned DIM, unsigned NNODE_1D, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
	class QConductingCellElementSolvesforVm :
		public virtual QElement<DIM, NNODE_1D>,
		public virtual ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		QConductingCellElementSolvesforVm()	:	QElement<DIM, NNODE_1D>(),
										ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>()
		{
			//set the integration scheme to one with integral points aligned with the nodes
			GaussWithNodes<DIM, NNODE_1D>* new_integral_pt = new GaussWithNodes<DIM, NNODE_1D>;
			this->set_integration_scheme(new_integral_pt);
			//set the number of integral points which are not aligned with nodes
			this->ipt_not_at_nodes = this->integral_pt()->nweight() - this->nnode();
		}

		QConductingCellElementSolvesforVm(const QConductingCellElementSolvesforVm<DIM, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL>& dummy){BrokenCopy::broken_copy("QConductingCellElementSolvesforVm");}

		void operator=(const QConductingCellElementSolvesforVm<DIM, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL>&){BrokenCopy::broken_assign("QConductingCellElementSolvesforVm");}



		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(file_pt, n_plot);
		}

		//======================================================================
		/// \short Define the shape functions and test functions and derivatives
		/// w.r.t. global coordinates and return Jacobian of mapping.
		///
		/// Galerkin: Test functions = shape functions
		//======================================================================
		double dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s,
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
		double dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt,
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
	};



	template<unsigned NNODE_1D, unsigned DIM, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
	class FaceGeometry<QConductingCellElementSolvesforVm<DIM, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL> >:
		public virtual QElement<DIM-1, NNODE_1D>
	{
	public:
		FaceGeometry()	:	QElement<DIM-1, NNODE_1D>()	{}
	};	

	template<unsigned NNODE_1D, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
	class FaceGeometry<QConductingCellElementSolvesforVm<1, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL> >:
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
	template<unsigned DIM, unsigned NNODE_1D, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
	class TConductingCellElementSolvesforVm :
		public virtual TElement<DIM, NNODE_1D>,
		public virtual ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>
	{
	private:

	public:
		//====================================================================
		//Constructors
		//====================================================================
		TConductingCellElementSolvesforVm()	:	TElement<DIM, NNODE_1D>(),
										ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>()
		{
			//set the integration scheme to one with integral points aligned with the nodes
			TGaussWithNodes<DIM, NNODE_1D>* new_integral_pt = new TGaussWithNodes<DIM, NNODE_1D>;
			this->set_integration_scheme(new_integral_pt);
			//set the number of integral points which are not aligned with nodes
			this->ipt_not_at_nodes = this->integral_pt()->nweight() - this->nnode();
		}

		TConductingCellElementSolvesforVm(const TConductingCellElementSolvesforVm<DIM, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL>& dummy){BrokenCopy::broken_copy("TConductingCellElementSolvesforVm");}

		void operator=(const TConductingCellElementSolvesforVm<DIM, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL>&){BrokenCopy::broken_assign("TConductingCellElementSolvesforVm");}



		//====================================================================
		//Output functions
		//====================================================================
		/// Output with default number of plot points
		void output(std::ostream &outfile){
			ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(outfile);
		}
		/// \short Output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		/// nplot^DIM plot points
		void output(std::ostream &outfile, const unsigned &nplot){
			ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(outfile, nplot);
		}
		/// C_style output with default number of plot points
		void output(FILE* file_pt){
			ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(file_pt);
		}
		 /// \short C-style output FE representation of soln: x,y,V_fct,[vars] or x,y,z,V_fct,[vars] at 
		 /// n_plot^DIM plot points
		 void output(FILE* file_pt, const unsigned &n_plot){
		 	ConductingCellEquationsSolvesforVm<CELL_MODEL, CONDUCTANCE_MODEL<DIM>>::output(file_pt, n_plot);
		}

		//======================================================================
		/// \short Define the shape functions and test functions and derivatives
		/// w.r.t. global coordinates and return Jacobian of mapping.
		///
		/// Galerkin: Test functions = shape functions
		//======================================================================
		double dshape_and_dtest_eulerian_BaseCellMembranePotential(const Vector<double> &s,
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
		double dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt,
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

	};

	template<unsigned DIM, unsigned NNODE_1D, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
	class FaceGeometry<TConductingCellElementSolvesforVm<DIM, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL> >:
		public virtual TElement<DIM-1, NNODE_1D>
	{
	public:
		FaceGeometry()	:	TElement<DIM-1, NNODE_1D>()	{}
	};	

	template<unsigned NNODE_1D, class CELL_MODEL, template<unsigned> class CONDUCTANCE_MODEL>
	class FaceGeometry<TConductingCellElementSolvesforVm<1, NNODE_1D, CELL_MODEL, CONDUCTANCE_MODEL> >:
		public virtual PointElement
	{
	public:
		FaceGeometry()	:	PointElement()	{}
	};








}//end namespace

#endif