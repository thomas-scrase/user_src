#include "meta_problem_class.h"

namespace oomph{

	void JacobianMetaProblem::setup_optimisation(const unsigned &n_values){
		//set the number of variables used in problem
		N_Variables = n_values;

		optimisation_equation_pt()->create_internal_data(N_Variables);

		//fill in the mean and permitted range of the variables
		send_add_internal_data_range_to_optimisation_elements();
		send_internal_data_mean_to_optimisation_elements();
		send_internal_data_percentage_range_to_optimisation_elements();
	}


	//Communicator with the optimisation element
	void JacobianMetaProblem::communicate_sub_problem_error_to_optimisation_element(Vector<double> &residual){
		residual.resize(Number_Of_Sub_Problems);
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			residual[i] = Sub_Problem_Error[i];
		}
	}


	//Add a sub problem 
	void JacobianMetaProblem::add_sub_problem_pt(SubProblem* sub_prob_pt){
		//Check if the new pointer is already in the vector
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			if(sub_prob_pt==Sub_Problem_Pts[i]){
				std::ostringstream warn_message;
				warn_message << "WARNING: Sub problem has already been added\n";
				OomphLibWarning(warn_message.str(),
								"JacobianMetaProblem::add_sub_problem_pt(SubProblem* sub_prob_pt)",
								OOMPH_EXCEPTION_LOCATION);
				return;
			}
		}

		//If the new problem has not been added yet, resize the number of sub problems
		Number_Of_Sub_Problems++;
		//resize the sub problem errors vector
		Sub_Problem_Error.resize(Number_Of_Sub_Problems);
		//Make the new vector
		Vector<SubProblem*> New_Sub_Problem_Pts(Number_Of_Sub_Problems);
		//Loop over the old vector and copy accross
		for(unsigned i=0; i<Number_Of_Sub_Problems-1; i++){
			New_Sub_Problem_Pts[i] = Sub_Problem_Pts[i];
		}
		//Add the new problem pointer
		New_Sub_Problem_Pts[Number_Of_Sub_Problems-1] = Sub_Problem_Pts[Number_Of_Sub_Problems-1];
		//resize the old vector
		Sub_Problem_Pts.resize(Number_Of_Sub_Problems);
		//Copy accross from the New vector
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			Sub_Problem_Pts[i] = New_Sub_Problem_Pts[i];
		}
	}

	void JacobianMetaProblem::pass_cell_model_pt_to_sub_problems(){
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			Sub_Problem_Pts[i]->assign_cell_model_pt(cell_model_pt());
		}
	}


	//Actions before newton solve
	void JacobianMetaProblem::actions_before_newton_solve() {
		//run the sub problems and assign their outputs to the error vector
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			Sub_Problem_Error[i] = Sub_Problem_Pts[i]->run();
		}
	}

	//send inform the optimisation elements if and how to calculate the cost due to
	//	variables lying outside of a range
	//Flag, add range of Internal_Data_Pt values to the residual
	void JacobianMetaProblem::send_add_internal_data_range_to_optimisation_elements(){
		//loop over variables
		for(unsigned var=0; var<N_Variables; var++){
			//check if true
			if(Add_Internal_Data_Range[var]){
				//turn on relative error contribution
				optimisation_equation_pt()->turn_on_relative_error_contribution(var);
			}
			else{
				//turn off relative error contribuution
				optimisation_equation_pt()->turn_off_relative_error_contribution(var);	
			}
		}
	}
	//The Mean of the internal data
	void JacobianMetaProblem::send_internal_data_mean_to_optimisation_elements(){
		//loop over variables
		for(unsigned var=0; var<N_Variables; var++){
			//set_variable_mean
			optimisation_equation_pt()->set_variable_mean(var, Internal_Data_Mean[var]);
		}
	}
	//The permitted range aroung the mean value
	void JacobianMetaProblem::send_internal_data_percentage_range_to_optimisation_elements(){
		//loop over variables
		for(unsigned var=0; var<N_Variables; var++){
			//set_variable_permitted_relative_error
			optimisation_equation_pt()->set_variable_permitted_relative_error(var, Internal_Data_Percentage_Range[var]);
		}
	}










	///Functions for the Nelder-Mead simplex meta problem

	//Destructor
	NelderMeadMetaProblemClass::~NelderMeadMetaProblemClass(){
		//delete the sub problems
		for(unsigned i=0; i<N_Variables+1; i++){
			delete Simplex_nodes[i];
		}

		//delete the simplex nodes
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			delete Sub_Problem_Pts[i];
		}
	}

	void NelderMeadMetaProblemClass::setup_optimisation(const unsigned &n_values){
		//can only be called once
		if(Simplex_nodes.size()>0){
			throw OomphLibError("simplex can only be setup once",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}
		//set the number of variables used in problem
		N_Variables = n_values;
		//make the simplex
		Simplex_nodes.resize(N_Variables+1);
		for(unsigned i=0; i<N_Variables+1;i++){
			//make the nodes
			Simplex_nodes[i] = new OptimisationEquations;
			//give the nodes the number of variables so they can properly set themselves up
			Simplex_nodes[i]->create_internal_data(N_Variables);
		}

		//fill in the mean and permitted range of the variables
		send_add_internal_data_range_to_optimisation_elements();
		send_internal_data_mean_to_optimisation_elements();
		send_internal_data_percentage_range_to_optimisation_elements();
	}

	//Add a sub problem 
	void NelderMeadMetaProblemClass::add_sub_problem_pt(SubProblem* sub_prob_pt){
		//Check if the new pointer is already in the vector
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			if(sub_prob_pt==Sub_Problem_Pts[i]){
				std::ostringstream warn_message;
				warn_message << "WARNING: Sub problem has already been added\n";
				OomphLibWarning(warn_message.str(),
								"JacobianMetaProblem::add_sub_problem_pt(SubProblem* sub_prob_pt)",
								OOMPH_EXCEPTION_LOCATION);
				return;
			}
		}

		//If the new problem has not been added yet, resize the number of sub problems
		Number_Of_Sub_Problems++;
		//resize the sub problem errors vector
		Sub_Problem_Error.resize(Number_Of_Sub_Problems);
		//Make the new vector
		Vector<SubProblem*> New_Sub_Problem_Pts(Number_Of_Sub_Problems);
		//Loop over the old vector and copy accross
		for(unsigned i=0; i<Number_Of_Sub_Problems-1; i++){
			New_Sub_Problem_Pts[i] = Sub_Problem_Pts[i];
		}
		//Add the new problem pointer
		New_Sub_Problem_Pts[Number_Of_Sub_Problems-1] = Sub_Problem_Pts[Number_Of_Sub_Problems-1];
		//resize the old vector
		Sub_Problem_Pts.resize(Number_Of_Sub_Problems);
		//Copy accross from the New vector
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			Sub_Problem_Pts[i] = New_Sub_Problem_Pts[i];
		}
	}

	void NelderMeadMetaProblemClass::pass_cell_model_pt_to_sub_problems(){
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			Sub_Problem_Pts[i]->assign_cell_model_pt(cell_model_pt());
		}
	}



	//send inform the optimisation elements if and how to calculate the cost due to
	//	variables lying outside of a range
	//Flag, add range of Internal_Data_Pt values to the residual
	void NelderMeadMetaProblemClass::send_add_internal_data_range_to_optimisation_elements(){
		//loop over the elements
		for(unsigned node=0; node<N_Variables+1; node++){
			//loop over variables
			for(unsigned var=0; var<N_Variables; var++){
				//check if true
				if(Add_Internal_Data_Range[var]){
					//turn on relative error contribution
					Simplex_nodes[node]->turn_on_relative_error_contribution(var);
				}
				else{
					//turn off relative error contribuution
					Simplex_nodes[node]->turn_off_relative_error_contribution(var);	
				}
			}
		}
	}
	//The Mean of the internal data
	void NelderMeadMetaProblemClass::send_internal_data_mean_to_optimisation_elements(){
		//loop over the elements
		for(unsigned node=0; node<N_Variables+1; node++){
			//loop over variables
			for(unsigned var=0; var<N_Variables; var++){
				//set_variable_mean
				Simplex_nodes[node]->set_variable_mean(var, Internal_Data_Mean[var]);
			}
		}
	}
	//The permitted range aroung the mean value
	void NelderMeadMetaProblemClass::send_internal_data_percentage_range_to_optimisation_elements(){
		//loop over the elements
		for(unsigned node=0; node<N_Variables+1; node++){
			//loop over variables
			for(unsigned var=0; var<N_Variables; var++){
				//set_variable_permitted_relative_error
				Simplex_nodes[node]->set_variable_permitted_relative_error(var, Internal_Data_Percentage_Range[var]);
			}
		}
	}




	void NelderMeadMetaProblemClass::run_algorithm(){
		Vector<double> Error; //error at each timestep
		Vector<double> Node_Qualities(N_Variables+1);
		unsigned best_node_index;
		unsigned worst_node_index;
		unsigned second_worst_node_index;
		Vector<unsigned> Sorted_node_indexes(N_Variables+1);

		//the extra nodes required for the algorithm to operate
		OptimisationEquations x0;
		OptimisationEquations xr;
		OptimisationEquations xe;
		OptimisationEquations xc;

		//evaluate entire simplex
		evaluate_quality_of_simplex(Node_Qualities);

		//rely on calculations from previous run, this way only nodes which are changed
		//	need to be reevaluated

		while(!terminate(Node_Qualities, Error)){
			//sort indexes from best to worst
			sort_nodes(Node_Qualities, Sorted_node_indexes);

			best_node_index = Sorted_node_indexes[0];
			worst_node_index = Sorted_node_indexes[N_Variables];
			second_worst_node_index = Sorted_node_indexes[N_Variables-1];

			//make x0
			fill_in_x0(x0, Sorted_node_indexes);

			//make xr
			fill_in_xr(xr, x0, *Simplex_nodes[second_worst_node_index]);
			//evaluate quality of xr
			double xr_quality;
			evaluate_quality_of_node(xr, xr_quality);

			if(Node_Qualities[best_node_index] <= xr_quality && xr_quality <= Node_Qualities[second_worst_node_index]){
				replace_node(*Simplex_nodes[worst_node_index], xr);
				//calculate the quality of the new node
				double new_worst_node_quality;
				evaluate_quality_of_node(*Simplex_nodes[worst_node_index],new_worst_node_quality);
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = new_worst_node_quality;
				//return to the beginning of the loop
				continue;
			}

			//if the reflected node quality is better than that of the best node
			if(xr_quality < Node_Qualities[best_node_index]){
				//make xe
				fill_in_xe(xe, x0, xr);
				//evaluate quality of xe
				double xe_quality;
				evaluate_quality_of_node(xe, xe_quality);

				//if the expansion node is better than the reflected node replace the worst node with the expansion node
				if(xe_quality < xr_quality){replace_node(*Simplex_nodes[worst_node_index], xe);}
				//else replace the worst node with the reflected node
				else{replace_node(*Simplex_nodes[worst_node_index], xr);}
				//calculate the quality of the new node
				double new_worst_node_quality;
				evaluate_quality_of_node(*Simplex_nodes[worst_node_index],new_worst_node_quality);
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = new_worst_node_quality;
				//return to the beginning of the loop
				continue;
			}


			//make xc
			fill_in_xc(xc, x0, *Simplex_nodes[worst_node_index]);
			//evaluate quality of xc
			double xc_quality;
			evaluate_quality_of_node(xc, xc_quality);

			//if the contraction node is better than the worst node replace then worst node with the contracted node
			if(xc_quality < Node_Qualities[worst_node_index]){
				replace_node(*Simplex_nodes[worst_node_index], xc);
				//calculate the quality of the new node
				double new_worst_node_quality;
				evaluate_quality_of_node(*Simplex_nodes[worst_node_index],new_worst_node_quality);
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = new_worst_node_quality;
				continue;
			}


			//shrink all points towards best
			for(unsigned i=1; i<N_Variables+1; i++){
				unsigned node_index = Sorted_node_indexes[i];
				shrink_node(*Simplex_nodes[node_index],*Simplex_nodes[best_node_index]);
				//calculate the quality of the new node
				double new_node_quality;
				evaluate_quality_of_node(*Simplex_nodes[node_index],new_node_quality);
				//record quality for the new node
				Node_Qualities[node_index] = new_node_quality;
			}

		}

		//terminate the program


	}



	void NelderMeadMetaProblemClass::evaluate_quality_of_simplex(Vector<double>& node_qualities){
		double quality;
		for(unsigned n=0; n < N_Variables; n++){
			evaluate_quality_of_node(*Simplex_nodes[n], quality);
			node_qualities[n] = quality;
		}
	}

	void NelderMeadMetaProblemClass::evaluate_quality_of_node(OptimisationEquations &node, double &quality){
		//pass the current node to the cell model
		dynamic_cast<TrainableCellModelBase*>(cell_model_pt())->set_optimisation_equations_pt(node);
		// cell_model_pt()->set_optimisation_equations_pt(node);

		//initialise quality to zero
		quality = 0.0;
		//run through the sub problems
		for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
			//add their output to the value of quality
			quality += Sub_Problem_Pts[i]->run();
		}

		//add the cost due to variable values existing outside of their permitted range
		quality += node.total_residual_from_data_outside_of_range();
	}


	bool NelderMeadMetaProblemClass::terminate(Vector<double> &node_qualities, Vector<double> &error_time_series){

		double new_error = error_from_node_qualities(node_qualities);

		Vector<double> new_error_time_series;
		new_error_time_series.resize(error_time_series.size()+1);

		for(unsigned i=0; i<error_time_series.size(); i++){
			new_error_time_series[i] = error_time_series[i];
		}

		new_error_time_series[error_time_series.size()] = new_error;

		if(new_error < 1e-9){
			return true;
		}
	}

	//the total error, default to the minimum node quality
	double NelderMeadMetaProblemClass::error_from_node_qualities(Vector<double> &node_qualities){
		return *std::min_element(node_qualities.begin(), node_qualities.end());
	}


	void NelderMeadMetaProblemClass::fill_in_x0(OptimisationEquations &node, Vector<unsigned> &sorted_node_indexes){
		//loop over the variables
		for(unsigned var=0; var<N_Variables; var++){
			double centroid_var = 0.0;
			//loop over the nodes in the simplex, skipping the worst one
			for(unsigned n=0; n<N_Variables; n++){
				//get the index of the current node
				unsigned index = sorted_node_indexes[n];
				//add the value from that node
				centroid_var += Simplex_nodes[index]->get_internal_data(var);
			}
			//take average
			centroid_var/=N_Variables;
			//assign the internal data for the node
			node.set_internal_data(var, centroid_var);
		}
	}

	void NelderMeadMetaProblemClass::fill_in_xr(OptimisationEquations &node, OptimisationEquations &x0, OptimisationEquations &xnp1){
		for(unsigned var=0; var<N_Variables; var++){
			double val = x0.get_internal_data(var) + alpha*(x0.get_internal_data(var) - xnp1.get_internal_data(var));
			node.set_internal_data(var, val);
		}
	}

	void NelderMeadMetaProblemClass::fill_in_xe(OptimisationEquations &node, OptimisationEquations &x0, OptimisationEquations &xr){
		for(unsigned var=0; var<N_Variables; var++){
			double val = x0.get_internal_data(var) + gamma*(xr.get_internal_data(var) - x0.get_internal_data(var));
			node.set_internal_data(var, val);
		}
	}

	void NelderMeadMetaProblemClass::fill_in_xc(OptimisationEquations &node, OptimisationEquations &x0, OptimisationEquations &xnp1){
		for(unsigned var=0; var<N_Variables; var++){
			double val = x0.get_internal_data(var) + rho*(xnp1.get_internal_data(var) - x0.get_internal_data(var));
			node.set_internal_data(var, val);
		}
	}

	void NelderMeadMetaProblemClass::shrink_node(OptimisationEquations &node, OptimisationEquations &x0){
		for(unsigned var=0; var<N_Variables; var++){
			double val = x0.get_internal_data(var) + sigma*(node.get_internal_data(var) - x0.get_internal_data(var));
			node.set_internal_data(var, val);
		}
	}

	void NelderMeadMetaProblemClass::replace_node(OptimisationEquations &node, OptimisationEquations &replacement_node){
		for(unsigned var=0; var<N_Variables; var++){
			double val = replacement_node.get_internal_data(var);
			node.set_internal_data(var, val);
		}
	}

}	//end namespace