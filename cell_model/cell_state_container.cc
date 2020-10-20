#include "cell_state_container.h"

namespace oomph{
	//Empty constructor and destructor
	CellState::CellState(){	}
	
	//Access functions for the base class member data
	void CellState::set_vars(const DenseMatrix<double> &new_vars){
		Cell_Variables_And_Derivatives.resize(new_vars.nrow(), new_vars.ncol());
		for(unsigned i = 0; i < new_vars.nrow(); i++){
			for(unsigned j = 0; j < new_vars.ncol(); j++){
				Cell_Variables_And_Derivatives(i,j) = new_vars(i,j);
			}
		}
	}
	//return the dth derivative of the nth cell state variable
	double CellState::get_var(const unsigned &d, const unsigned &n) const {return Cell_Variables_And_Derivatives(d,n);}

	//time_stepper_weights
	void CellState::set_time_stepper_weights(const DenseMatrix<double> &new_time_stepper_weights){
		Time_Stepper_Weights.resize(new_time_stepper_weights.nrow(), new_time_stepper_weights.ncol());
		for(unsigned i = 0; i < new_time_stepper_weights.nrow(); i++){
			for(unsigned j = 0; j < new_time_stepper_weights.ncol(); j++){
				Time_Stepper_Weights(i,j) = new_time_stepper_weights(i,j);
			}
		}
	}
	//get the dth derivative weight of the nth cell state variable
	double CellState::get_time_stepper_weight(const unsigned &d, const unsigned &n) const {return Time_Stepper_Weights(d,n);}

	//cell_type
	void CellState::set_cell_type(const unsigned &new_cell_type){Cell_Type = new_cell_type;}
	unsigned CellState::get_cell_type() const {return Cell_Type;}

	//Black-box nodal parameters
	void CellState::set_black_box_nodal_parameters(const Vector<double> &new_nodal_parameters){
		Black_Box_Nodal_Parameters = new_nodal_parameters;
	}
	double CellState::get_black_box_nodal_parameters(const unsigned &paramter_index) const {return Black_Box_Nodal_Parameters[paramter_index];}

	//Trans-Membrane potential
	void CellState::set_vm(const double &new_vm){Membrane_Potential = new_vm;}
	double CellState::get_vm() const {return Membrane_Potential;}

	//Mechanical strain
	void CellState::set_stress(const double &new_stress){Mechanical_Stress = new_stress;}
	double CellState::get_stress() const {return Mechanical_Stress;}

	//Black-box external data
	void CellState::set_black_box_external_data(const Vector<double> &new_external_data){
		Black_Box_External_Data.resize(new_external_data.size());
		for(unsigned i=0; i<new_external_data.size();i++){
			Black_Box_External_Data[i] = new_external_data[i];
		}
	}
	double CellState::get_black_box_external_data(const unsigned &data_index) const {
		return Black_Box_External_Data[data_index];
	}

	/////////Explicit cell model////
	//dt Value
	void CellState::set_dt(const double &new_dt){Dt = new_dt;}
	double CellState::get_dt() const {return Dt;}

	//Previous values
	void CellState::set_previous_variables(const Vector<double> &new_previous_variables){
		Previous_Variables.resize(new_previous_variables.size());
		for(unsigned i = 0; i < new_previous_variables.size(); i++){
			Previous_Variables[i] = new_previous_variables[i];
		}
	}
	double CellState::get_previous_variables(const unsigned &n) const {return Previous_Variables[n];}
	//////////////////////////////

	//Total cell membrane current
	void CellState::set_membrane_current(const double &new_cell_membrane_current){Membrane_Current = new_cell_membrane_current;}
	double CellState::get_membrane_current() const {return Membrane_Current;}

	//Active strain
	void CellState::set_active_strain(const double &new_cell_model_strain){Active_Strain = new_cell_model_strain;}
	double CellState::get_active_strain() const {return Active_Strain;}
	
	//General cell model data - generated and changed by the cell model class during a run sequence
	//	allows for compartmentalisation of the model into several member functions
	void CellState::set_new_general_cell_model_data(const double &new_data){ //appends to the general cell model data
		General_cell_model_data.push_back(new_data);
	}
	void CellState::set_general_cell_model_data_index(const unsigned &data_index, const double &new_data){ //change the value of a specific index general cell model data
		General_cell_model_data[data_index] = new_data;
	}
	void CellState::resize_general_cell_model_data(const unsigned &new_size){	//resize the vector representing general cell model data
		General_cell_model_data.resize(new_size,0.0);
	}
	double CellState::get_general_cell_model_data(const unsigned &data_index) const { //get the data at data_index
		return General_cell_model_data[data_index];
	}


	
	//Force build CellState
	class CellState;
}