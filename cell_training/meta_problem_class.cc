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





	//Default variables used by the class
	double NelderMeadMetaProblemClass::Default_Acceptable_Edge_Length = 1e-4;




	///Functions for the Nelder-Mead simplex meta problem

	//Destructor
	NelderMeadMetaProblemClass::~NelderMeadMetaProblemClass(){
		//delete the simplex
		for(unsigned i=0; i<N_Variables+1; i++){
			delete Simplex_nodes[i];
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
			// Simplex_nodes[i]->create_internal_data(N_Variables);
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
	void NelderMeadMetaProblemClass::output(const unsigned &iteration, Vector<double> &Node_Qualities,  std::ostream &outfile){
		outfile << "Results from iteration " << iteration << std::endl;
		for(unsigned i=0; i<N_Variables+1; i++){
			Simplex_nodes[i]->output(outfile);
			outfile << "\t" << Node_Qualities[i] << std::endl;
		}
	}


	//add dependent element
	void NelderMeadMetaProblemClass::add_trainable_element_as_dependent(TrainableElement* dependent_element){
		Dependent_Elements.push_back(dependent_element);
	}

	//add all trainable elements from all sub problems
	void NelderMeadMetaProblemClass::add_all_trainable_elements_from_sub_problems_as_dependents(){
		//loop over the sub problems
		for(unsigned i=0; i< Sub_Problem_Pts.size(); i++){
			//a temporary vector
			Vector<TrainableElement*> trainable_elements_from_sub_problem;
			//get the list of trainable elements from the problem
			Sub_Problem_Pts[i]->get_all_trainable_elements(trainable_elements_from_sub_problem);
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
		for(unsigned i=0; i<Dependent_Elements.size(); i++){
			Dependent_Elements[i]->set_parameter_source_pt(node);
		}
	}


	void NelderMeadMetaProblemClass::evaluate_quality_of_simplex(Vector<double>& node_qualities){
		double quality=0.0;
		for(unsigned n=0; n < N_Variables+1; n++){
			quality=0.0;
			evaluate_quality_of_node(*Simplex_nodes[n], quality);
			node_qualities[n] = quality;
		}
	}
	

	void NelderMeadMetaProblemClass::evaluate_quality_of_node(OptimisationEquations &node, double &quality){
		//push the node to all dependent elements
		push_optimisation_element_to_dependent_elements(node);
		//declare the residal vector we will be filling in using the optimisation equations element
		Vector<double> residuals(N_Variables, 0.0);
		//populate the residuals vector
		node.fill_in_contribution_to_residuals(residuals);
		//condense the residuals into a single number and pass it to quality
		for(unsigned i=0; i<residuals.size(); i++){
			quality += residuals[i];
		}
	}


	//push the optimisation element to and run the sub problems and get their norms
	// used as ExternalSourceContribution of the optimisation elements
	void NelderMeadMetaProblemClass::get_sub_problem_contribution_to_node(OptimisationEquations &node, Vector<double> &contribution){
		//resize the contribution vector to the number of sub problems
		contribution.resize(Sub_Problem_Pts.size());
		//run each of the problems and add their norm to the contribution vector
		for(unsigned i=0; i<Sub_Problem_Pts.size(); i++){
			contribution[i] = Sub_Problem_Pts[i]->run();
		}
	}

	//get the cost incurred by the values of the variables in the optimisation element passed
	//used as VariableValuesContribution for the optimisation elements
	void NelderMeadMetaProblemClass::get_cost_of_variables(OptimisationEquations &optimisation_equations, Vector<double> &costs){
		costs.resize(N_Variables, 0.0);
	}


	bool NelderMeadMetaProblemClass::terminate(Vector<double> &node_qualities){
		//get the measure we are using as the error, by default this is the minimum node quality
		// double simplex_maximum_edge_length = error_from_node_qualities(node_qualities);
		
		double Simplex_Maximum_Edge_Length = simplex_maximum_edge_length();
		
		if(std::abs(Simplex_Maximum_Edge_Length) < termination_tolerance()){
			return true;
		}
		else{
			return false;
		}
	}

	//the total error, default to the minimum node quality
	double NelderMeadMetaProblemClass::error_from_node_qualities(Vector<double> &node_qualities){
		double min_val = node_qualities[0];
		for(unsigned i=1; i<node_qualities.size(); i++){
			if(std::abs(node_qualities[i]) < std::abs(min_val)){
				min_val = node_qualities[i];
			}
		}
		return min_val;
	}


	double NelderMeadMetaProblemClass::simplex_maximum_edge_length(){
		double max_length = 0.0;
		//loop over all edges and get their lengths
		for(unsigned node1=0; node1<N_Variables+1; node1++){
			for(unsigned node2=node1+1; node2<N_Variables+1; node2++){
				double length = 0.0;
				for(unsigned index=0; index<N_Variables; index++){
					length += Simplex_nodes[node1]->get_internal_data(index) - Simplex_nodes[node2]->get_internal_data(index);
				}
				if(std::abs(length) > max_length){max_length = std::abs(length);}
			}
		}
		return max_length;
	}



	void NelderMeadMetaProblemClass::sort_nodes(const Vector<double> &node_qualities, Vector<unsigned> &sorted_node_indexes){
		//index sorting algorithm
		sorted_node_indexes.resize(node_qualities.size());
		for(unsigned i=0; i<node_qualities.size(); i++){sorted_node_indexes[i] = i;}
		stable_sort(sorted_node_indexes.begin(), sorted_node_indexes.end(), [node_qualities](size_t i1, size_t i2) {return std::abs(node_qualities[i1]) < std::abs(node_qualities[i2]);});
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

	//run the algorithm to completion
	void NelderMeadMetaProblemClass::run_algorithm(std::ostream &outfile){
		unsigned iteration = 0;
		Vector<double> Node_Qualities(N_Variables+1);
		unsigned best_node_index;
		unsigned worst_node_index;
		unsigned second_worst_node_index;
		Vector<unsigned> Sorted_node_indexes(N_Variables+1);

		//evaluate entire simplex
		evaluate_quality_of_simplex(Node_Qualities);

		//rely on calculations from previous run, this way only nodes which are changed
		//	need to be reevaluated

		while(!terminate(Node_Qualities)){
			iteration++;

			//report on the state of the simplex
			output(iteration, Node_Qualities, outfile);


			//sort indexes from best to worst
			sort_nodes(Node_Qualities, Sorted_node_indexes);
			outfile << std::endl << "Sorted indices: ";
			for(unsigned i=0; i<N_Variables+1; i++){
				outfile << "\t" << Sorted_node_indexes[i] << " (" << Node_Qualities[Sorted_node_indexes[i]] << ")\t";
			}
			outfile << std::endl;

			best_node_index = Sorted_node_indexes[0];
			worst_node_index = Sorted_node_indexes[N_Variables];
			second_worst_node_index = Sorted_node_indexes[N_Variables-1];

			//make x0
			OptimisationEquations x0(this);
			// x0.create_internal_data(N_Variables);
			fill_in_x0(x0, Sorted_node_indexes);

			//make xr
			OptimisationEquations xr(this);
			// xr.create_internal_data(N_Variables);
			fill_in_xr(xr, x0, *Simplex_nodes[worst_node_index]);
			//evaluate quality of xr
			double xr_quality=0.0;
			evaluate_quality_of_node(xr, xr_quality);

			if(Node_Qualities[best_node_index] <= xr_quality && xr_quality <= Node_Qualities[second_worst_node_index]){
				replace_node(*Simplex_nodes[worst_node_index], xr);
				//calculate the quality of the new node
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = xr_quality;
				//return to the beginning of the loop
				outfile << "Replaced node " << worst_node_index << " for xr" << std::endl;
				continue;
			}

			//if the reflected node quality is better than that of the best node
			if(xr_quality < Node_Qualities[best_node_index]){
				//make xe
				OptimisationEquations xe(this);
				// xe.create_internal_data(N_Variables);
				fill_in_xe(xe, x0, xr);
				//evaluate quality of xe
				double xe_quality=0.0;
				evaluate_quality_of_node(xe, xe_quality);

				//if the expansion node is better than the reflected node replace the worst node with the expansion node
				if(xe_quality < xr_quality){
					replace_node(*Simplex_nodes[worst_node_index], xe);
					outfile << "Replaced node " << worst_node_index << " for xe" << std::endl;
					Node_Qualities[worst_node_index] = xe_quality;
				}
				//else replace the worst node with the reflected node
				else{
					replace_node(*Simplex_nodes[worst_node_index], xr);
					outfile << "Replaced node " << worst_node_index << " for xr" << std::endl;
					Node_Qualities[worst_node_index] = xr_quality;
				}
				//return to the beginning of the loop
				continue;
			}


			//make xc
			OptimisationEquations xc(this);
			// xc.create_internal_data(N_Variables);
			fill_in_xc(xc, x0, *Simplex_nodes[worst_node_index]);
			//evaluate quality of xc
			double xc_quality=0.0;
			evaluate_quality_of_node(xc, xc_quality);

			//if the contraction node is better than the worst node replace then worst node with the contracted node
			if(xc_quality < Node_Qualities[worst_node_index]){
				replace_node(*Simplex_nodes[worst_node_index], xc);
				outfile << "Replaced node " << worst_node_index << " for xc" << std::endl;
				//calculate the quality of the new node
				// double new_worst_node_quality;
				// evaluate_quality_of_node(*Simplex_nodes[worst_node_index],new_worst_node_quality);
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = xc_quality;
				continue;
			}


			//shrink all points towards best
			for(unsigned i=1; i<N_Variables+1; i++){
				unsigned node_index = Sorted_node_indexes[i];
				shrink_node(*Simplex_nodes[node_index],*Simplex_nodes[best_node_index]);
				//calculate the quality of the new node
				double new_node_quality=0.0;
				evaluate_quality_of_node(*Simplex_nodes[node_index],new_node_quality);
				//record quality for the new node
				Node_Qualities[node_index] = new_node_quality;
			}
			outfile << "Shrank all nodes towards " << best_node_index << std::endl;

		}
		//terminate the program

		//report on the state of the simplex
		outfile << "Simplex at end of program" << std::endl;
		output(iteration, Node_Qualities, outfile);

	}




		//run the algorithm with probabilistic restarts
		void NelderMeadMetaProblemClass::run_algorithm_with_restarts(std::ostream &outfile){
			
		}

}	//end namespace