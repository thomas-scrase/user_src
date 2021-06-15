//LIC// ====================================================================
//LIC// Contains the base cell model class along with two simple cell models
//LIC// ZeroCell - purely for testing
//LIC// FitzHighNagumo - an old redundant model also for testing, albeit more
//LIC//		interesting
//LIC// ====================================================================



#ifndef OOMPH_CELL_MODEL_UPDATED_HEADER
#define OOMPH_CELL_MODEL_UPDATED_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

//OOMPH-LIB includes
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/matrices.h"
#include "../generic/Qelements.h"
#include "../generic/Telements.h"


namespace oomph{


	class CellModelBaseUpdated
	{
	public:
		CellModelBaseUpdated() : active_strain_index(-1)
		{
			//Resize the variable names vectors to zero
			Names_Of_Cell_Variables.resize(0);
			Names_Of_Other_Parameters.resize(0);
			Names_Of_Other_Variables.resize(0);
			Names_Of_Output_Data.resize(0);

			// std::cout << "Names of variables size at base construction " << Names_Of_Cell_Variables.size() << std::endl;
		}
		virtual ~CellModelBaseUpdated(){}
		

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//BEGIN functions which need to be overridden by a user

		//Get the v-th variables initial condition for cell_type,
		//	these MUST be written in the same order that you write the Names_Of_Cell_Variables vector entries.
		virtual double return_initial_state_variable(const unsigned &v, const unsigned &cell_type)
		{
			//Broken by default
			throw OomphLibError("return_initial_state_variable: This function has not been implemented yet",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}


		//Get the initial membrane potential for the cell_type-th cell type
		virtual double return_initial_membrane_potential(const unsigned &cell_type)
		{
			//Broken by default
			throw OomphLibError("return_initial_membrane_potential: This function has not been implemented yet",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}


		//Calculate the derivatives of the cell variables, also calculate the trans membrane current
		//	We provide the membrane current, the contemporary cell variable valyues, the time (just in case)
		//	the timestep (for rush larson methods), the cell type,
		//	other parameters (which represets time independent variables, e.g. a drug concentration or ion channel blockade),
		//	other variables (generally time dependent variables, e.g. strain in the tissue)
		//
		//	The function is to calculate the time derivatives of the cell variables, and populate the Variable_Derivatives,
		//	 the user must use the same keyword names as the variable, e.g. CellVariables[sxs], Variable_Derivatives[sxs] = (...)
		//	The user also provides the Iion which is the total transmembrane ionic current, which must be in pA/pF, i.e. current per capacitance
		//	 this way the membrane capacitance doesn't need to be passed to any other methods and the computations are neater.
		virtual void Calculate_Derivatives(const double &Vm,
											const Vector<double> &CellVariables,
											const double &t,
											const unsigned &cell_type,
											const double &Istim,
											const Vector<double> &Other_Parameters,
											const Vector<double> &Other_Variables,

											Vector<double> &Variable_Derivatives,
											double &Iion)
		{
			//Broken by default
			throw OomphLibError("Calculate_Derivatives: This function has not been implemented yet",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}


		//Generate an unordered map of data required by things outside of the cell model, this could be
		//	active strain, or perhaps ion channel currents, or ionic species concentrations, we want to keep track of
		virtual void get_optional_output(const double &Vm,
									const Vector<double> &CellVariables,
									const double &t,
									const unsigned &cell_type,
									const double &Istim,
									const Vector<double> &Other_Parameters,
									const Vector<double> &Other_Variables,
									Vector<double> &Out) const
		{
			//Broken by default
			throw OomphLibError("get_optional_output: This function has not been implemented yet",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}

		//END functions which need to be overridden by a user
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		//The names given to the cell variables, we use this because it's more user friendly,
		// it means a user does not need to also provide the number of cell variables, and they
		// can just call the variables by name rather than worrying about index
		const std::vector<std::string>& names_of_cell_variables() const
		{
			return Names_Of_Cell_Variables;
		}

		//The names of the other parameters passed to the function, we provide this for error checking,
		//	if the library is built with PARANOID then we can automatically check to ensure that the
		//	user is correctly creating the other parameters data.
		const std::vector<std::string>& names_of_other_parameters() const
		{
			return Names_Of_Other_Parameters;
		}

		//Same as above but for other variables which are in general time dependent.
		const std::vector<std::string>& names_of_other_variables() const
		{
			return Names_Of_Other_Variables;
		}

		//Same as above but for the data the cell model returns for other things to use
		const std::vector<std::string>& names_of_output_data() const 
		{
			return Names_Of_Output_Data;
		}


		const unsigned get_Num_Cell_Vars() const{
			return Num_Cell_Vars;
		}
		const unsigned get_Num_Other_Params() const{
			return 	Num_Other_Params;
		}
		const unsigned get_Num_Other_Vars() const{
			return 	Num_Other_Vars;
		}
		const unsigned get_Num_Output_Data() const{
			return 	Num_Output_Data;
		}

		inline unsigned get_index_of_active_strain() const{
			#ifdef PARANOID
			// oomph_info << "In fct which should die" << std::endl;
			//If the user has failed to specify at which index the active strain is stored at then we throw an error 
			if(active_strain_index<0){
				throw OomphLibError("You appear to be attempting to access active strain stored in a cell element\nbut you haven't set the index at which it is stored. Please make sure that you override\nget_index_of_active_strain() to return the correct index you have stored it at.",
			                       	OOMPH_CURRENT_FUNCTION,
			                       	OOMPH_EXCEPTION_LOCATION);
			}
			#endif
			return active_strain_index;
		}

	protected:
		//The vectors of variable names, these need to be assigned at construction of the cell model.
		//	these are used for the headers of output files
		std::vector<std::string> Names_Of_Cell_Variables;
		std::vector<std::string> Names_Of_Other_Parameters;
		std::vector<std::string> Names_Of_Other_Variables;
		std::vector<std::string> Names_Of_Output_Data;

		//The number of variables, these avoid unnecessary overhead by constantly calling Names_Of_Cell_Variables.size() etc
		//	They are set automatically by the following function which must be called by the user at the end of their cell
		//	class
		unsigned Num_Cell_Vars;
		unsigned Num_Other_Params;
		unsigned Num_Other_Vars;
		unsigned Num_Output_Data;

		void FinalizeConstruction()
		{
			Num_Cell_Vars = 	Names_Of_Cell_Variables.size();
			Num_Other_Params = 	Names_Of_Other_Parameters.size();
			Num_Other_Vars = 	Names_Of_Other_Variables.size();
			Num_Output_Data = 	Names_Of_Output_Data.size();

			// std::cout << Num_Cell_Vars << std::endl;
			// std::cout << Num_Other_Params << std::endl;
			// std::cout << Num_Other_Vars << std::endl;
			// std::cout << Num_Output_Data << std::endl;
		}


		int active_strain_index;



		//Enums which contain the names of the variables, these are used for several reasons
		//	1. Make it easier to index vectors containing variables and variable derivatives
		//	2. Make removing and adding variables from and to a model easier. You don't need
		//		to shift the indexing of vectors, just use the correct enumerator from the
		//		enums below
		// enum Cell_Variables_Enum	: unsigned;
		// enum Other_Parameters_Enum	: unsigned;
		// enum Other_Variables_Enum	: unsigned;
		// enum Output_Data_Enum		: unsigned;


	private:

	};







	// //====================================================================
	// //====================================================================
	// //Combined cell model, useful for computations over an entire organ
	// //	i.e. atria computed by CNZ and ventricles by TNNP
	// //	This way multiple cell meshes are not required.
	// //====================================================================
	// //====================================================================
	// template<class CELL_MODEL_1, class CELL_MODEL_2>
	// class CombinedCellModel :	public CELL_MODEL_1,
	// 							public CELL_MODEL_2
	// {
	// public:
	// 	CombinedCellModel() :	CELL_MODEL_1(), CELL_MODEL_2(){ }

	// 	//Identify the cell model which is to be used in the subsequent computation:
	// 	//	by default this function will simply check if the passed cell type compatible
	// 	//	with both cell models. If it is only compatible with one then return true for
	// 	//	CELL_MODEL_1, or false for CELL_MODEL_2. If it is compatible with both then
	// 	//	in PARANOID throw an error.
	// 	//	Implemented as virtual so that user defined cell type distributions over
	// 	//	the cell models can be used
	// 	//		i.e. both models are of atria type but right atrium is to be computed
	// 	//		by CELL_MODEL_2 and the rest are to be computed by CELL_MODEL_1
	// 	virtual inline bool Identify_Correct_Cell_Model(const unsigned& cell_type){
	// 		#ifdef PARANOID
	// 		if(CELL_MODEL_1::compatible_cell_types(cell_type) and CELL_MODEL_2::compatible_cell_types(cell_type)){
	// 			throw OomphLibError("Cell type is compatible with both cell models, since paranoid is enabled I\n"
	// 								"am killing the process. Please use non-overlapping cell models or alternatively redefine this function to\n"
	// 								"a user defined cell type distribution",
	// 	                       OOMPH_CURRENT_FUNCTION,
	// 	                       OOMPH_EXCEPTION_LOCATION);
	// 		}
	// 		#endif
			
	// 		if(CELL_MODEL_1::compatible_cell_types(cell_type)){
	// 			return true;
	// 		}
	// 		if(CELL_MODEL_2::compatible_cell_types(cell_type)){
	// 			return false;
	// 		}

	// 		//if we get here then neither cell model can handle the passed cell type
	// 		throw OomphLibError("Cell type is not compatible with either cell model",
	// 	                       OOMPH_CURRENT_FUNCTION,
	// 	                       OOMPH_EXCEPTION_LOCATION);

	// 	}
	// 	//return the correct model's membrane current
	// 	double membrane_current(CellState &state)
	// 	{
	// 		if(Identify_Correct_Cell_Model(state.get_cell_type())){
	// 			return CELL_MODEL_1::membrane_current(state);
	// 		}
	// 		else{
	// 			return CELL_MODEL_2::membrane_current(state);
	// 		}
	// 	}
	// 	//return the correct model's active strain
	// 	double active_strain(CellState &state)
	// 	{
	// 		if(Identify_Correct_Cell_Model(state.get_cell_type())){
	// 			return CELL_MODEL_1::active_strain(state);
	// 		}
	// 		else{
	// 			return CELL_MODEL_2::active_strain(state);
	// 		}
	// 	}

	// 	//return the correct model's membrane capacitance
	// 	// double cm(CellState &state) {
	// 	// 	if(Identify_Correct_Cell_Model(state.get_cell_type())){
	// 	// 		return CELL_MODEL_1::cm(state);
	// 	// 	}
	// 	// 	else{
	// 	// 		return CELL_MODEL_2::cm(state);
	// 	// 	}
	// 	// }

	// 	inline void custom_output(CellState &state, Vector<double> &output)
	// 	{
	// 		if(Identify_Correct_Cell_Model(state.get_cell_type())){
	// 			CELL_MODEL_1::custom_output(state, output);
	// 		}
	// 		else{
	// 			CELL_MODEL_2::custom_output(state, output);
	// 		}
	// 	}

	// 	// Calculate the sub residual and sub jacobian objects
	// 	void fill_in_generic_residual_contribution_cell_base(CellState &state,
	// 														Vector<double> &residuals,
	// 														DenseMatrix<double> &jacobian,
	// 														unsigned flag)
	// 	{
	// 		if(Identify_Correct_Cell_Model(state.get_cell_type())){
	// 			CELL_MODEL_1::fill_in_generic_residual_contribution_cell_base(state,residuals, jacobian, flag);
	// 		}
	// 		else{
	// 			CELL_MODEL_2::fill_in_generic_residual_contribution_cell_base(state,residuals, jacobian, flag);
	// 		}
	// 	}

	// 	inline void return_initial_membrane_potential(double &v, const unsigned &cell_type=0){
	// 		if(Identify_Correct_Cell_Model(cell_type)){
	// 			CELL_MODEL_1::return_initial_membrane_potential(v, cell_type);
	// 		}
	// 		else{
	// 			CELL_MODEL_2::return_initial_membrane_potential(v, cell_type);
	// 		}
	// 	}

	// 	inline bool return_initial_state_variable(const unsigned &n, double &v, const unsigned &cell_type=0){
	// 		if(Identify_Correct_Cell_Model(cell_type)){
	// 			return CELL_MODEL_1::return_initial_state_variable(n, v, cell_type);
	// 		}
	// 		else{
	// 			return CELL_MODEL_2::return_initial_state_variable(n, v, cell_type);
	// 		}
	// 	}

		
	// 	//It would be very difficult to write a cell model which performs automatic Jacobian fill in
	// 	//	via finite differencing for some nodes and not for others. For this reason this flag is
	// 	//	only set to true if both models can calculate residuals
	// 	inline bool model_calculates_jacobian_entries() {
	// 		if(CELL_MODEL_1::model_calculates_jacobian_entries() && CELL_MODEL_2::model_calculates_jacobian_entries()){
	// 			return true;
	// 		}
	// 	 	else{
	// 	 		return false;
	// 	 	}
	// 	}
	// 	//Return the maximum required storage
	// 	inline unsigned required_nodal_variables(){
	// 		return std::max(CELL_MODEL_1::Num_Variables(),CELL_MODEL_2::Num_Variables());
	// 	}
	// 	//return the maximum number of required derivatives
	// 	inline unsigned required_derivatives(const unsigned &cell_type=0){
	// 		if(Identify_Correct_Cell_Model(cell_type)){
	// 			return CELL_MODEL_1::required_derivatives();
	// 		}
	// 		else{
	// 			return CELL_MODEL_2::required_derivatives();
	// 		}
	// 	}
	// 	//return the maximum number of required black box parameters
	// 	inline unsigned required_black_box_parameters(){
	// 		return std::max(CELL_MODEL_1::required_black_box_parameters(),CELL_MODEL_2::required_black_box_parameters());
	// 	}
	// };





}

#endif