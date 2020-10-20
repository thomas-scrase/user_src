#include "meta_problem_class.h"

namespace oomph{
	//Default variables used by the class
	double NelderMeadOptimisation::Default_Acceptable_Edge_Length = 1e-4;




	///Functions for the Nelder-Mead simplex meta problem

	//Destructor
	NelderMeadOptimisation::~NelderMeadOptimisation(){
	}

	void NelderMeadOptimisation::setup_optimisation(const unsigned &n_values){
		//can only be called once
		if(Simplex.size()>0){
			throw OomphLibError("simplex can only be setup once",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}
		//set the number of variables used in problem
		N_Variables = n_values;
		//make the simplex
		Simplex.resize(n_values+1);
		for(unsigned i=0; i< n_values+1; i++){
			Simplex[i].resize(n_values);
		}

		Node_Qualities.resize(N_Variables+1);

		//Get all elements which depend on parameters supplied by the simplex
		add_all_trainable_elements_from_trainable_problems_as_dependents();
	}

	void NelderMeadOptimisation::set_simplex_value(const unsigned& node, const unsigned& var, const double& value){
		Simplex[node][var] = value;
	}

	//Add a sub problem 
	void NelderMeadOptimisation::add_trainable_problem(TrainableProblem& trainable_problem){
		Trainable_Problem_Pts.push_back(&trainable_problem);
	}

	//report on the nodes
	void NelderMeadOptimisation::output(const unsigned &iteration, Vector<double> &Node_Qualities, std::ostream &outfile){
		for(unsigned i=0; i<N_Variables+1; i++){
			outfile << "Node " << i << ": ( ";
			for(unsigned var=0; var<N_Variables; var++){
				outfile << Simplex[i][var] << " ";
			}
			outfile << ") -> " << Node_Qualities[i] << std::endl;
		}
	}


	//add dependent element
	void NelderMeadOptimisation::add_trainable_element_as_dependent(TrainableElement* dependent_element){
		Dependent_Elements.push_back(dependent_element);
	}

	//add all trainable elements from all sub problems
	void NelderMeadOptimisation::add_all_trainable_elements_from_trainable_problems_as_dependents(){
		//loop over the sub problems
		for(unsigned i=0; i< Trainable_Problem_Pts.size(); i++){
			//a temporary vector
			Vector<TrainableElement*> trainable_elements_from_sub_problem;
			//get the list of trainable elements from the problem
			Trainable_Problem_Pts[i]->get_all_trainable_elements(trainable_elements_from_sub_problem);
			//loop over the elements
			for(unsigned j=0; j<trainable_elements_from_sub_problem.size(); j++){
				//add them all
				add_trainable_element_as_dependent(trainable_elements_from_sub_problem[j]);
			}
		}
	}


	//send a particular optimisation element to all dependent trainable elements
	//	performed before sub problems are run
	void NelderMeadOptimisation::push_optimisation_element_to_dependent_elements(Vector<double> &node){
		for(unsigned i=0; i<Dependent_Elements.size(); i++){
			Dependent_Elements[i]->set_parameter_source_pt(node);
		}
	}


	void NelderMeadOptimisation::evaluate_quality_of_simplex(){
		double quality;
		for(unsigned n=0; n < N_Variables+1; n++){
			quality=0.0;
			evaluate_quality_of_node(Simplex[n], quality);
			Node_Qualities[n] = quality;
		}
	}
	

	void NelderMeadOptimisation::evaluate_quality_of_node(Vector<double> &node, double &quality){
		//push the node to all dependent elements
		push_optimisation_element_to_dependent_elements(node);
		//condense the residuals from each trainable problem into
		//	a single number and pass it to quality
		for(unsigned i=0; i<Trainable_Problem_Pts.size(); i++){
			quality += Trainable_Problem_Pts[i]->run();
		}
	}

	double NelderMeadOptimisation::simplex_maximum_edge_length(){
		double max_length = 0.0;
		//loop over all edges and get their lengths
		for(unsigned node1=0; node1<N_Variables+1; node1++){
			for(unsigned node2=node1+1; node2<N_Variables+1; node2++){
				double length = 0.0;
				for(unsigned index=0; index<N_Variables; index++){
					length += Simplex[node1][index] - Simplex[node2][index];
				}
				if(std::abs(length) > max_length){max_length = std::abs(length);}
			}
		}
		return max_length;
	}



	void NelderMeadOptimisation::sort_nodes(Vector<unsigned> &sorted_node_indexes){
		//index sorting algorithm
		sorted_node_indexes.resize(Node_Qualities.size());
		for(unsigned i=0; i<N_Variables+1; i++){sorted_node_indexes[i] = i;}
		
		// std::sort(sorted_node_indexes.begin(), sorted_node_indexes.end(), Node_Qualities);;



		for(unsigned i=0; i < N_Variables;i++){
			for(unsigned j=0; j < N_Variables-i; j++){
				if(abs(Node_Qualities[sorted_node_indexes[j]]) > abs(Node_Qualities[sorted_node_indexes[j+1]])){
					unsigned tmp = sorted_node_indexes[j];
					sorted_node_indexes[j] = sorted_node_indexes[j+1];
					sorted_node_indexes[j+1] = tmp;
				}
			}
		}
	}


	void NelderMeadOptimisation::fill_in_x0(Vector<double> &node, Vector<unsigned> &sorted_node_indexes){
		//loop over the variables
		for(unsigned var=0; var<N_Variables; var++){
			node[var] = 0.0;
			//loop over the nodes in the simplex, skipping the worst one
			for(unsigned n=0; n<N_Variables; n++){
				//get the index of the current node
				unsigned index = sorted_node_indexes[n];
				//add the value from that node
				node[var] += Simplex[index][var];
			}
			//take average
			node[var]/=N_Variables;
		}
	}

	void NelderMeadOptimisation::fill_in_xr(Vector<double> &node, Vector<double> &x0, Vector<double> &xnp1){
		for(unsigned var=0; var<N_Variables; var++){
			node[var] = x0[var] + alpha*(x0[var] - xnp1[var]);
		}
	}

	void NelderMeadOptimisation::fill_in_xe(Vector<double> &node, Vector<double> &x0, Vector<double> &xr){
		for(unsigned var=0; var<N_Variables; var++){
			node[var] = x0[var] + gamma*(xr[var] - x0[var]);
		}
	}

	void NelderMeadOptimisation::fill_in_xc(Vector<double> &node, Vector<double> &x0, Vector<double> &xnp1){
		for(unsigned var=0; var<N_Variables; var++){
			node[var] = x0[var] + rho*(xnp1[var] - x0[var]);
		}
	}

	void NelderMeadOptimisation::shrink_node(Vector<double> &node, Vector<double> &x0){
		for(unsigned var=0; var<N_Variables; var++){
			double val = 
			node[var] = x0[var] + sigma*(node[var] - x0[var]);
		}
	}

	void NelderMeadOptimisation::replace_node(Vector<double> &node, Vector<double> &replacement_node){
		node = replacement_node;
	}

	//run one step of the nelder mead algorithm
	void NelderMeadOptimisation::run_algorithm(std::ostream &outfile){
		outfile << "This is the Nelder-Mead Coordinator." << std::endl;
		unsigned iteration = 0;
		unsigned best_node_index;
		unsigned worst_node_index;
		unsigned second_worst_node_index;
		Vector<unsigned> Sorted_node_indexes(N_Variables+1);

		//evaluate entire simplex
		evaluate_quality_of_simplex();

		//rely on calculations from previous run, this way only nodes which are changed
		//	need to be reevaluated

		while(simplex_maximum_edge_length() > Acceptable_Edge_Length){
			iteration++;

			//report on the state of the simplex
			output(iteration, Node_Qualities, outfile);


			//sort indexes from best to worst
			sort_nodes(Sorted_node_indexes);
			outfile << std::endl << "Sorted indices: ";
			for(unsigned i=0; i<N_Variables+1; i++){
				outfile << "\t" << Sorted_node_indexes[i] << " (" << Node_Qualities[Sorted_node_indexes[i]] << ")\t";
			}
			outfile << std::endl;

			best_node_index = Sorted_node_indexes[0];
			worst_node_index = Sorted_node_indexes[N_Variables];
			second_worst_node_index = Sorted_node_indexes[N_Variables-1];

			//make x0
			Vector<double> x0(N_Variables);
			// x0.create_internal_data(N_Variables);
			fill_in_x0(x0, Sorted_node_indexes);

			//make xr
			Vector<double> xr(N_Variables);
			// xr.create_internal_data(N_Variables);
			fill_in_xr(xr, x0, Simplex[worst_node_index]);
			//evaluate quality of xr
			double xr_quality=0.0;
			evaluate_quality_of_node(xr, xr_quality);

			if(Node_Qualities[best_node_index] <= xr_quality && xr_quality <= Node_Qualities[second_worst_node_index]){
				replace_node(Simplex[worst_node_index], xr);
				//calculate the quality of the new node
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = xr_quality;
				//return to the beginning of the loop
				outfile << "Replaced node " << worst_node_index << " for xr:\t";
				output(xr, outfile);
				outfile << std::endl;
				continue;
			}

			//if the reflected node quality is better than that of the best node
			if(xr_quality < Node_Qualities[best_node_index]){
				//make xe
				Vector<double> xe(N_Variables);
				// xe.create_internal_data(N_Variables);
				fill_in_xe(xe, x0, xr);
				//evaluate quality of xe
				double xe_quality=0.0;
				evaluate_quality_of_node(xe, xe_quality);

				//if the expansion node is better than the reflected node replace the worst node with the expansion node
				if(xe_quality < xr_quality){
					replace_node(Simplex[worst_node_index], xe);
					outfile << "Replaced node " << worst_node_index << " for xe:\t";
					output(xe, outfile);
					outfile << std::endl;
					Node_Qualities[worst_node_index] = xe_quality;
				}
				//else replace the worst node with the reflected node
				else{
					replace_node(Simplex[worst_node_index], xr);
					outfile << "Replaced node " << worst_node_index << " for xr:\t";
					output(xr, outfile);
					outfile << std::endl;
					Node_Qualities[worst_node_index] = xr_quality;
				}
				//return to the beginning of the loop
				continue;
			}


			//make xc
			Vector<double> xc(N_Variables);
			// xc.create_internal_data(N_Variables);
			fill_in_xc(xc, x0, Simplex[worst_node_index]);
			//evaluate quality of xc
			double xc_quality=0.0;
			evaluate_quality_of_node(xc, xc_quality);

			//if the contraction node is better than the worst node replace then worst node with the contracted node
			if(xc_quality < Node_Qualities[worst_node_index]){
				replace_node(Simplex[worst_node_index], xc);
				outfile << "Replaced node " << worst_node_index << " for xc:\t";
				output(xc, outfile);
				outfile << std::endl;
				//calculate the quality of the new node
				// double new_worst_node_quality;
				// evaluate_quality_of_node(Simplex[worst_node_index],new_worst_node_quality);
				//record quality for the new worst node
				Node_Qualities[worst_node_index] = xc_quality;
				continue;
			}


			//shrink all points towards best
			for(unsigned i=1; i<N_Variables+1; i++){
				unsigned node_index = Sorted_node_indexes[i];
				shrink_node(Simplex[node_index],Simplex[best_node_index]);
				//calculate the quality of the new node
				double new_node_quality=0.0;
				evaluate_quality_of_node(Simplex[node_index],new_node_quality);
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

}	//end namespace