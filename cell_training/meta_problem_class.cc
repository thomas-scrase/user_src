#include "meta_problem_class.h"

namespace oomph{

	// void JacobianMetaProblem::setup_optimisation(const unsigned &n_values){
	// 	//set the number of variables used in problem
	// 	N_Variables = n_values;

	// 	optimisation_equation_pt()->create_internal_data(N_Variables);
	// }


	// //Communicator with the optimisation element
	// void JacobianMetaProblem::communicate_sub_problem_error_to_optimisation_element(Vector<double> &residual){
	// 	residual.resize(Number_Of_Sub_Problems);
	// 	for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
	// 		residual[i] = Sub_Problem_Error[i];
	// 	}
	// }


	// //Add a sub problem 
	// void JacobianMetaProblem::add_sub_problem_pt(SubProblem* sub_prob_pt){
	// 	//Check if the new pointer is already in the vector
	// 	for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
	// 		if(sub_prob_pt==Sub_Problem_Pts[i]){
	// 			std::ostringstream warn_message;
	// 			warn_message << "WARNING: Sub problem has already been added\n";
	// 			OomphLibWarning(warn_message.str(),
	// 							"JacobianMetaProblem::add_sub_problem_pt(SubProblem* sub_prob_pt)",
	// 							OOMPH_EXCEPTION_LOCATION);
	// 			return;
	// 		}
	// 	}

	// 	//If the new problem has not been added yet, resize the number of sub problems
	// 	Number_Of_Sub_Problems++;
	// 	//resize the sub problem errors vector
	// 	Sub_Problem_Error.resize(Number_Of_Sub_Problems);
	// 	//Make the new vector
	// 	Vector<SubProblem*> New_Sub_Problem_Pts(Number_Of_Sub_Problems);
	// 	//Loop over the old vector and copy accross
	// 	for(unsigned i=0; i<Number_Of_Sub_Problems-1; i++){
	// 		New_Sub_Problem_Pts[i] = Sub_Problem_Pts[i];
	// 	}
	// 	//Add the new problem pointer
	// 	New_Sub_Problem_Pts[Number_Of_Sub_Problems-1] = Sub_Problem_Pts[Number_Of_Sub_Problems-1];
	// 	//resize the old vector
	// 	Sub_Problem_Pts.resize(Number_Of_Sub_Problems);
	// 	//Copy accross from the New vector
	// 	for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
	// 		Sub_Problem_Pts[i] = New_Sub_Problem_Pts[i];
	// 	}
	// }


	// //Actions before newton solve
	// void JacobianMetaProblem::actions_before_newton_solve() {
	// 	//run the sub problems and assign their outputs to the error vector
	// 	for(unsigned i=0; i<Number_Of_Sub_Problems; i++){
	// 		Sub_Problem_Error[i] = Sub_Problem_Pts[i]->run();
	// 	}
	// }

	// //After a solve we want to record the current variables and performance etc so we do that
	// void JacobianMetaProblem::actions_after_newton_solve(){

	// }







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
			Simplex_nodes[i] = new OptimisationEquations(this);
			//give the nodes the number of variables so they can properly set themselves up
			Simplex_nodes[i]->create_internal_data(N_Variables);
		}

		add_all_trainable_elements_from_sub_problems_as_dependents();
	}

	void NelderMeadMetaProblemClass::set_initial_simplex(const Vector<double>& min_values, const Vector<double>& max_values){
		//loop over the nodes
		for(unsigned node=0; node<N_Variables+1;node++){
			//loop over the variables
			for(unsigned var=0; var<N_Variables;var++){
				//set all variables to their default
				Simplex_nodes[node]->set_internal_data(var, min_values[var]);		
			}
			if(node>0){
				//set the nodeth value to the max value
				Simplex_nodes[node]->set_internal_data(node-1, max_values[node-1]);
			}
		}
	}

	//Add a sub problem 
	void NelderMeadMetaProblemClass::add_sub_problem_pt(SubProblem* sub_prob_pt){
		//if there are no sub problems yet
		if(Number_Of_Sub_Problems == 0){
			//just add it
			Number_Of_Sub_Problems=1;
			Sub_Problem_Error.resize(1);
			Sub_Problem_Pts.resize(1);
			Sub_Problem_Pts[0] = sub_prob_pt;
		}
		else{
			//more care is needed to ensure it hasn't been added yet and
			//	to make sure we don't lose any prevously added sub problems
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
	}

	//report on the nodes
	void NelderMeadMetaProblemClass::output(const unsigned &iteration, std::ostream &outfile){
		// std::cout << "NelderMeadMetaProblemClass::output" << std::endl;
		outfile << "Results from iteration " << iteration << std::endl;
		for(unsigned i=0; i<N_Variables+1; i++){
			Simplex_nodes[i]->output(outfile);
			double node_quality;
			evaluate_quality_of_node(*Simplex_nodes[i], node_quality);
			outfile << "\t" << node_quality << std::endl;
		}
	}


	void NelderMeadMetaProblemClass::run_algorithm(std::ostream &outfile){
		Vector<double> Error; //error at each timestep
		unsigned iteration = 0;
		Vector<double> Node_Qualities(N_Variables+1);
		unsigned best_node_index;
		unsigned worst_node_index;
		unsigned second_worst_node_index;
		Vector<unsigned> Sorted_node_indexes(N_Variables+1);

		

		//evaluate entire simplex
		evaluate_quality_of_simplex(Node_Qualities);
		std::cout << "node qualities\t";
		for(unsigned i=0; i<Node_Qualities.size(); i++){
			std::cout << Node_Qualities[i] << "\t";
		}
		std::cout << std::endl;

		//rely on calculations from previous run, this way only nodes which are changed
		//	need to be reevaluated

		while(!terminate(Node_Qualities, Error)){

			//report on the state of the simplex
			output(iteration, outfile);


			//sort indexes from best to worst
			sort_nodes(Node_Qualities, Sorted_node_indexes);

			best_node_index = Sorted_node_indexes[0];
			worst_node_index = Sorted_node_indexes[N_Variables];
			second_worst_node_index = Sorted_node_indexes[N_Variables-1];

			//make x0
			OptimisationEquations x0(this);
			x0.create_internal_data(N_Variables);
			fill_in_x0(x0, Sorted_node_indexes);

			//make xr
			OptimisationEquations xr(this);
			xr.create_internal_data(N_Variables);
			fill_in_xr(xr, x0, *Simplex_nodes[second_worst_node_index]);
			//evaluate quality of xr
			double xr_quality;
			evaluate_quality_of_node(xr, xr_quality);

			std::cout << "Node_Qualities[best_node_index] " << Node_Qualities[best_node_index] << std::endl;
			std::cout << "xr_quality " << xr_quality << std::endl;
			std::cout << "Node_Qualities[second_worst_node_index] " << Node_Qualities[second_worst_node_index] << std::endl;

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
				OptimisationEquations xe(this);
				xe.create_internal_data(N_Variables);
				fill_in_xe(xe, x0, xr);
				//evaluate quality of xe
				double xe_quality;
				evaluate_quality_of_node(xe, xe_quality);

				std::cout << "xe_quality " << xe_quality << std::endl;

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
			OptimisationEquations xc(this);
			xc.create_internal_data(N_Variables);
			fill_in_xc(xc, x0, *Simplex_nodes[worst_node_index]);
			//evaluate quality of xc
			double xc_quality;
			evaluate_quality_of_node(xc, xc_quality);

			std::cout << "xc_quality " << xc_quality << std::endl;

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

		//report on the state of the simplex
		output(iteration, outfile);

	}

	//add dependent element
	void NelderMeadMetaProblemClass::add_trainable_element_as_dependent(TrainableElement* dependent_element){
		// std::cout << "add_trainable_element_as_dependent" << std::endl;
		//make a backup
		// Vector<TrainableElement*> temp_dependent_elements(Dependent_Elements.size());
		// for(unsigned i=0; i<Dependent_Elements.size(); i++){
		// 	temp_dependent_elements[i] = Dependent_Elements[i];
		// }
		// //resize the Dependent_Elements vector to fit the new element
		// Dependent_Elements.resize(Dependent_Elements.size()+1);
		// //copy from the backup to the original
		// for(unsigned i=0; i<temp_dependent_elements.size(); i++){
		// 	Dependent_Elements[i] = temp_dependent_elements[i];
		// }
		// //add the new element
		// Dependent_Elements[temp_dependent_elements.size()] = dependent_element;
		Dependent_Elements.push_back(dependent_element);
		// std::cout << Dependent_Elements.size() << std::endl;
	}

	//add all trainable elements from all sub problems
	void NelderMeadMetaProblemClass::add_all_trainable_elements_from_sub_problems_as_dependents(){
		// std::cout << "add_all_trainable_elements_from_sub_problems_as_dependents" << std::endl;
		//loop over the sub problems
		for(unsigned i=0; i< Sub_Problem_Pts.size(); i++){
			//a temporary vector
			Vector<TrainableElement*> trainable_elements_from_sub_problem;
			//get the list of trainable elements from the problem
			// std::cout << Sub_Problem_Pts[i] << std::endl;
			Sub_Problem_Pts[i]->get_all_trainable_elements(trainable_elements_from_sub_problem);
			// std::cout << "trainable_elements_from_sub_problem " << trainable_elements_from_sub_problem.size() << std::endl;
			//loop over the elements
			for(unsigned j=0; j<trainable_elements_from_sub_problem.size(); j++){
				//add them all
				add_trainable_element_as_dependent(trainable_elements_from_sub_problem[j]);
			}
		}
	}


	//send a particular optimisation element to all dependent trainable elements
	//	performed before sub problems are run
	void NelderMeadMetaProblemClass::push_optimisation_element_to_dependent_elements(OptimisationEquations &node){
		// std::cout << "push_optimisation_element_to_dependent_elements" << std::endl;
		for(unsigned i=0; i<Dependent_Elements.size(); i++){
			// std::cout << "\t" << i << std::endl;
			Dependent_Elements[i]->set_parameter_source_pt(node);
		}
	}


	void NelderMeadMetaProblemClass::evaluate_quality_of_simplex(Vector<double>& node_qualities){
		// std::cout << "evaluate_quality_of_simplex" << std::endl;
		double quality;
		for(unsigned n=0; n < N_Variables+1; n++){
			evaluate_quality_of_node(*Simplex_nodes[n], quality);
			node_qualities[n] = quality;
		}
	}
	

	void NelderMeadMetaProblemClass::evaluate_quality_of_node(OptimisationEquations &node, double &quality){
		// std::cout << "evaluate_quality_of_node" << std::endl;
		//declare the residal vector we will be filling in using the optimisation equations element
		Vector<double> residuals(N_Variables);
		//populate the residuals vector
		node.fill_in_contribution_to_residuals(residuals);
		//condense the residuals into a single number and pass it to quality
		for(unsigned i=0; i<residuals.size(); i++){quality += residuals[i];}
	}


	//push the optimisation element to and run the sub problems and get their norms
	// used as ExternalSourceContribution of the optimisation elements
	void NelderMeadMetaProblemClass::get_sub_problem_contribution_to_node(OptimisationEquations &node, Vector<double> &contribution){
		// std::cout << "get_sub_problem_contribution_to_node" << std::endl;
		//push the node to all dependent elements
		push_optimisation_element_to_dependent_elements(node);
		//resize the contribution vector to the number of sub problems
		contribution.resize(Sub_Problem_Pts.size());
		// std::cout << "contribution.size() " << contribution.size() << std::endl;
		//run each of the problems and add their norm to the contribution vector
		for(unsigned i=0; i<Sub_Problem_Pts.size(); i++){
			contribution[i] = Sub_Problem_Pts[i]->run();
		}

		// std::cout << "contribution.size() " << contribution.size() << std::endl;
	}

	//get the cost incurred by the values of the variables in the optimisation element passed
	//used as VariableValuesContribution for the optimisation elements
	void NelderMeadMetaProblemClass::get_cost_of_variables(OptimisationEquations &optimisation_equations, Vector<double> &costs){
		costs.resize(N_Variables, 0.0);
	}


	bool NelderMeadMetaProblemClass::terminate(Vector<double> &node_qualities, Vector<double> &error_time_series){
		// std::cout << "NelderMeadMetaProblemClass::terminate" << std::endl;
		// for(unsigned i=0; i<node_qualities.size(); i++){
		// 	std::cout << node_qualities[i] << std::endl;
		// }
		//get the measure we are using as the error, by default this is the minimum node quality
		double new_error = error_from_node_qualities(node_qualities);
		// std::cout << "error_from_node_qualities " << new_error << std::endl;

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
		// std::cout << "NelderMeadMetaProblemClass::error_from_node_qualities" << std::endl;
		// return *std::min_element(node_qualities.begin(), node_qualities.end());
		double min_val = node_qualities[0];
		for(unsigned i=1; i<node_qualities.size(); i++){
			if(node_qualities[i] < min_val){min_val = node_qualities[i];}
		}
		// std::cout << min_val << std::endl;
		return min_val;
	}

	void NelderMeadMetaProblemClass::sort_nodes(const Vector<double> &node_qualities, Vector<unsigned> &sorted_node_indexes){
		//index sorting algorithm
		sorted_node_indexes.resize(node_qualities.size());
		for(unsigned i=0; i<node_qualities.size(); i++){sorted_node_indexes[i] = 0;}
		stable_sort(sorted_node_indexes.begin(), sorted_node_indexes.end(), [node_qualities](size_t i1, size_t i2) {return node_qualities[i1] < node_qualities[i2];});
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