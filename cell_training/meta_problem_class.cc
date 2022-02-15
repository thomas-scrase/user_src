#include "meta_problem_class.h"
#include <random>


namespace oomph{
	//Default variables used by the class
	double NelderMeadOptimisation::Default_Acceptable_Edge_Length = 1e-4;
	double NelderMeadOptimisation::Default_Acceptable_Homogenous_Fitness = 1e-4;




	///Functions for the Nelder-Mead simplex meta problem

	//Destructor
	NelderMeadOptimisation::~NelderMeadOptimisation(){}

	void NelderMeadOptimisation::setup_optimisation(const unsigned &n_values){
		//can only be called once
		if(Simplex.size()>0){
			throw OomphLibError("simplex can only be setup once",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}

		if(n_values==1){
			throw OomphLibError("Nelder Mead cannot solve for a single variable",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}
		//set the number of variables used in problem
		this->n_variables() = n_values;
		//make the simplex
		Simplex.resize(n_values+1);
		for(unsigned i=0; i< n_values+1; i++){
			Simplex[i].resize(n_values);
		}

		Node_Fitnesses.resize(this->n_variables()+1);

		//Get all elements which depend on parameters supplied by the simplex
		add_all_trainable_elements_from_trainable_problems_as_dependents();

		//The record of number of node evaluations performed
		Evaluations_performed = 0;

		//Set these to a big value
		Default_Minimum_Permitted_Value = -1e9;
		Default_Maximum_Permitted_Value = 1e9;
		Minimum_Permitted_Values.resize(this->n_variables());
		Maximum_Permitted_Values.resize(this->n_variables());
		for(unsigned i=0; i<this->n_variables(); i++){
			Minimum_Permitted_Values[i] = Default_Minimum_Permitted_Value;
			Maximum_Permitted_Values[i] = Default_Maximum_Permitted_Value;
		}
	}

	void NelderMeadOptimisation::set_simplex_value(const unsigned& node, const unsigned& var, const double& value){
		Simplex[node][var] = value;
	}

	//report on the nodes
	void NelderMeadOptimisation::output(const unsigned &iteration, Vector<double> &Node_Fitnesses, std::ostream &outfile){
		for(unsigned i=0; i<this->n_variables()+1; i++){
			outfile << "Node " << i << ": ( ";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << Simplex[i][var] << " ";
			}
			outfile << ") -> " << Node_Fitnesses[i] << std::endl;
		}
	}

	void NelderMeadOptimisation::check_node_for_forbidden_values(Vector<double> &node){
		for(unsigned i=0; i<this->n_variables(); i++){
			if(node[i]<Minimum_Permitted_Values[i]){node[i] += 2.0*(Minimum_Permitted_Values[i]-node[i]);}
			if(node[i]>Maximum_Permitted_Values[i]){node[i] -= 2.0*(node[i]-Maximum_Permitted_Values[i]);}
		}
	}

	void NelderMeadOptimisation::evaluate_fitness_of_simplex(std::ostream &raw_data_outfile){
		double fitness;
		for(unsigned n=0; n < this->n_variables()+1; n++){
			fitness=0.0;
			evaluate_fitness_of_node(Simplex[n], fitness, raw_data_outfile);
			Node_Fitnesses[n] = fitness;
		}
	}
	

	void NelderMeadOptimisation::evaluate_fitness_of_node(Vector<double> &node,
														double &fitness,
														std::ostream &raw_data_outfile){
		if(Num_Concurrent_Node_Evaluation>Max_Num_Concurrent_Node_Evaluation){
			Vector<unsigned> Sorted_node_indexes(this->n_variables()+1);
			sort_nodes(Sorted_node_indexes);
			node = Simplex[Sorted_node_indexes[0]]; //just set it to be the best node if we cannot find another non-nan node
		}
		//push the node to all dependent elements
		push_optimisation_parameter_values_to_dependent_elements(node);
		//condense the residuals from each trainable problem into
		//	a single number and pass it to fitness
		for(unsigned i=0; i<Trainable_Problem_Pts.size(); i++){
			fitness += Trainable_Problem_Pts[i]->run();
		}
		if(std::isnan(fitness)){
			fitness = 1e300;

			//So it's easy to match run number with output from sub problems
			// raw_data_outfile << Evaluations_performed << "\t\t";
			// for(unsigned var=0; var<this->n_variables(); var++){
			// 	raw_data_outfile << node[var] << " ";
			// }
			// raw_data_outfile << std::to_string(std::nan("")) << std::endl;

			// Evaluations_performed++;




			// //Zero the faulty node
			// for(unsigned v=0; v<this->n_variables(); v++){
			// 	node[v] = 0.0;
			// }

			// //Create a new random node from the rest of the simplex				
			// for(unsigned i=0; i<this->n_variables()+1; i++){
			// 	std::srand(std::time(0));
			// 	double ran_val = 2.0*((double)std::rand()/(RAND_MAX) - 0.5); //random number between -1 and 1
			// 	for(unsigned v=0; v<this->n_variables(); v++){
			// 		if(std::fabs(Simplex[i][v] - node[v])<1e-12){continue;}//If the node vertex is identical to the one creating a nan then skip it
			// 		node[v] += Simplex[i][v]/(this->n_variables()); //Take random weighted average over all nodes in the simplex, except the one we are changing
			// 	}
			// }

			// //Call evaluate fitness again
			// fitness = 0.0;
			// Num_Concurrent_Node_Evaluation++;
			// return evaluate_fitness_of_node(node, fitness, raw_data_outfile);
		}

		//So it's easy to match run number with output from sub problems
		raw_data_outfile << Evaluations_performed << "\t\t";
		for(unsigned var=0; var<this->n_variables(); var++){
			raw_data_outfile << node[var] << " ";
		}
		raw_data_outfile << fitness << std::endl;

		Num_Concurrent_Node_Evaluation = 0;
		Evaluations_performed++;
	}

	double NelderMeadOptimisation::simplex_maximum_edge_length(){
		double max_length = 0.0;
		//loop over all edges and get their lengths
		for(unsigned node1=0; node1<this->n_variables()+1; node1++){
			for(unsigned node2=node1+1; node2<this->n_variables()+1; node2++){
				double length = 0.0;
				for(unsigned index=0; index<this->n_variables(); index++){
					length += Simplex[node1][index] - Simplex[node2][index];
				}
				if(std::fabs(length) > max_length){max_length = std::fabs(length);}
			}
		}
		return max_length;
	}



	void NelderMeadOptimisation::sort_nodes(Vector<unsigned> &sorted_node_indexes){
		//index sorting algorithm
		sorted_node_indexes.resize(Node_Fitnesses.size());
		for(unsigned i=0; i<this->n_variables()+1; i++){sorted_node_indexes[i] = i;}
	
		for(unsigned i=0; i < this->n_variables();i++){
			for(unsigned j=0; j < this->n_variables()-i; j++){
				if(std::fabs(Node_Fitnesses[sorted_node_indexes[j]]) > std::fabs(Node_Fitnesses[sorted_node_indexes[j+1]])){
					unsigned tmp = sorted_node_indexes[j];
					sorted_node_indexes[j] = sorted_node_indexes[j+1];
					sorted_node_indexes[j+1] = tmp;
				}
			}
		}
	}


	void NelderMeadOptimisation::fill_in_x0(Vector<double> &node, Vector<unsigned> &sorted_node_indexes){
		//loop over the variables
		for(unsigned var=0; var<this->n_variables(); var++){
			node[var] = 0.0;
			//loop over the nodes in the simplex, skipping the worst one
			for(unsigned n=0; n<this->n_variables(); n++){
				//get the index of the current node
				unsigned index = sorted_node_indexes[n];
				//add the value from that node
				node[var] += Simplex[index][var];
			}
			//take average
			node[var]/=this->n_variables();
		}
	}

	void NelderMeadOptimisation::fill_in_xr(Vector<double> &node, Vector<double> &x0, Vector<double> &xnp1){
		for(unsigned var=0; var<this->n_variables(); var++){
			node[var] = x0[var] + alpha*(x0[var] - xnp1[var]);
		}
	}

	void NelderMeadOptimisation::fill_in_xe(Vector<double> &node, Vector<double> &x0, Vector<double> &xr){
		for(unsigned var=0; var<this->n_variables(); var++){
			node[var] = x0[var] + gamma*(xr[var] - x0[var]);
		}
	}

	void NelderMeadOptimisation::fill_in_xc(Vector<double> &node, Vector<double> &x0, Vector<double> &xnp1){
		for(unsigned var=0; var<this->n_variables(); var++){
			node[var] = x0[var] + rho*(xnp1[var] - x0[var]);
		}
	}

	void NelderMeadOptimisation::shrink_node(Vector<double> &node, Vector<double> &x0){
		for(unsigned var=0; var<this->n_variables(); var++){
			double val = 
			node[var] = x0[var] + sigma*(node[var] - x0[var]);
		}
	}

	void NelderMeadOptimisation::replace_node(Vector<double> &node, Vector<double> &replacement_node){
		node = replacement_node;
	}

	//run one step of the nelder mead algorithm
	Vector<double> NelderMeadOptimisation::run_algorithm(std::ostream &outfile,
												std::ostream &raw_data_outfile){
		outfile << "This is the Nelder-Mead Coordinator." << std::endl;
		Iterations = 0;
		unsigned best_node_index;
		unsigned worst_node_index;
		unsigned second_worst_node_index;
		Vector<unsigned> Sorted_node_indexes(this->n_variables()+1);

		//evaluate entire simplex
		evaluate_fitness_of_simplex(raw_data_outfile);

		//rely on calculations from previous run, this way only nodes which are changed
		//	need to be reevaluated
		double Maximum_Edge_Length = simplex_maximum_edge_length();
		while(Maximum_Edge_Length > Acceptable_Edge_Length && Iterations < Max_Iters){

			//Check if all fitnesses are identical, if they are then terminate
			bool simplex_is_uniformly_fit = true;
			for(unsigned i=1; i<this->n_variables()+1; i++){
				if(std::fabs(Node_Fitnesses[i]-Node_Fitnesses[i-1])>Acceptable_Homogenous_Fitness){simplex_is_uniformly_fit = false;break;}
			}
			if(simplex_is_uniformly_fit){
				outfile << "Fitness of all nodes is identical and we don't know how to proceed. Terminating." << std::endl;
				sort_nodes(Sorted_node_indexes);
				best_node_index = Sorted_node_indexes[0];
				//Return the best calculated node
				return Simplex[best_node_index];
			}

			//Calculate maximum edge length
			Maximum_Edge_Length = simplex_maximum_edge_length();
			Iterations++;

			//sort indexes from best to worst
			sort_nodes(Sorted_node_indexes);
			outfile << "Sorted indices: ";
			for(unsigned i=0; i<this->n_variables()+1; i++){
				outfile << "\t" << Sorted_node_indexes[i] << " (" << Node_Fitnesses[Sorted_node_indexes[i]] << ")\t";
			}
			outfile << std::endl << std::endl;

			best_node_index = Sorted_node_indexes[0];
			worst_node_index = Sorted_node_indexes[this->n_variables()];
			second_worst_node_index = Sorted_node_indexes[this->n_variables()-1];


			//Detailed output
			outfile << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
			outfile << "Begin detailed rundown, pre Iterations " << Iterations << "." << std::endl;
			//report on the state of the simplex
			output(Iterations, Node_Fitnesses, outfile);
			outfile << "Start of Iterations detailed rundown: " << std::endl;
			outfile << "best node: " << best_node_index << ": (";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << Simplex[best_node_index][var] << " ";
			}
			outfile << ") -> " << Node_Fitnesses[best_node_index] << std::endl;
			outfile << "worst node: " << worst_node_index << ": (";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << Simplex[worst_node_index][var] << " ";
			}
			outfile << ") -> " << Node_Fitnesses[worst_node_index] << std::endl;
			outfile << "second worst node: " << second_worst_node_index << ": (";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << Simplex[second_worst_node_index][var] << " ";
			}
			outfile << ") -> " << Node_Fitnesses[second_worst_node_index] << std::endl;
			outfile << "End detailed rundown" << std::endl;
			outfile << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

			

			//make x0
			Vector<double> x0(this->n_variables());
			fill_in_x0(x0, Sorted_node_indexes);
			check_node_for_forbidden_values(x0);
			outfile << "Centroid: ( ";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << x0[var] << " ";
			}
			outfile << ")" << std::endl;

			//make xr
			Vector<double> xr(this->n_variables());
			fill_in_xr(xr, x0, Simplex[worst_node_index]);
			check_node_for_forbidden_values(xr);
			//evaluate fitness of xr
			double xr_fitness=0.0;
			evaluate_fitness_of_node(xr, xr_fitness, raw_data_outfile);
			outfile << "xr: ( ";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << xr[var] << " ";
			}
			outfile << ") -> " << xr_fitness <<  std::endl;

			if(std::fabs(Node_Fitnesses[best_node_index]) <= std::fabs(xr_fitness) && std::fabs(xr_fitness) <= std::fabs(Node_Fitnesses[second_worst_node_index])){
				replace_node(Simplex[worst_node_index], xr);
				//calculate the fitness of the new node
				//record fitness for the new worst node
				Node_Fitnesses[worst_node_index] = xr_fitness;
				//return to the beginning of the loop
				outfile << "Replaced node " << worst_node_index << " for xr:\t";
				output(xr, outfile);
				outfile << std::endl;
				continue;
			}

			//if the reflected node fitness is better than that of the best node
			if(std::fabs(xr_fitness) < std::fabs(Node_Fitnesses[best_node_index])){
				//make xe
				Vector<double> xe(this->n_variables());
				fill_in_xe(xe, x0, xr);
				check_node_for_forbidden_values(xe);
				//evaluate fitness of xe
				double xe_fitness=0.0;
				evaluate_fitness_of_node(xe, xe_fitness, raw_data_outfile);
				outfile << "xe: ( ";
				for(unsigned var=0; var<this->n_variables(); var++){
					outfile << xe[var] << " ";
				}
				outfile << ") -> " << xe_fitness <<  std::endl;

				//if the expansion node is better than the reflected node replace the worst node with the expansion node
				if(std::fabs(xe_fitness) < std::fabs(xr_fitness)){
					replace_node(Simplex[worst_node_index], xe);
					outfile << "Replaced node " << worst_node_index << " for xe:\t";
					output(xe, outfile);
					outfile << std::endl;
					Node_Fitnesses[worst_node_index] = xe_fitness;
				}
				//else replace the worst node with the reflected node
				else{
					replace_node(Simplex[worst_node_index], xr);
					outfile << "Replaced node " << worst_node_index << " for xr:\t";
					output(xr, outfile);
					outfile << std::endl;
					Node_Fitnesses[worst_node_index] = xr_fitness;
				}
				//return to the beginning of the loop
				continue;
			}


			//make xc
			Vector<double> xc(this->n_variables());
			fill_in_xc(xc, x0, Simplex[worst_node_index]);
			check_node_for_forbidden_values(xc);
			//evaluate fitness of xc
			double xc_fitness=0.0;
			evaluate_fitness_of_node(xc, xc_fitness, raw_data_outfile);
			outfile << "xc: ( ";
			for(unsigned var=0; var<this->n_variables(); var++){
				outfile << xc[var] << " ";
			}
			outfile << ") -> " << xc_fitness <<  std::endl;

			//if the contraction node is better than the worst node replace then worst node with the contracted node
			if(std::fabs(xc_fitness) < std::fabs(Node_Fitnesses[worst_node_index])){
				replace_node(Simplex[worst_node_index], xc);
				outfile << "Replaced node " << worst_node_index << " for xc:\t";
				output(xc, outfile);
				outfile << std::endl;
				//record fitness for the new worst node
				Node_Fitnesses[worst_node_index] = xc_fitness;
				continue;
			}


			//shrink all points towards best
			for(unsigned i=1; i<this->n_variables()+1; i++){
				unsigned node_index = Sorted_node_indexes[i];
				shrink_node(Simplex[node_index],Simplex[best_node_index]);
				//calculate the fitness of the new node
				double new_node_fitness=0.0;
				evaluate_fitness_of_node(Simplex[node_index],new_node_fitness, raw_data_outfile);
				//record fitness for the new node
				Node_Fitnesses[node_index] = new_node_fitness;

				outfile << "shrank node " << node_index << ": ( ";
				for(unsigned var=0; var<this->n_variables(); var++){
					outfile << Simplex[node_index][var] << " ";
				}
				outfile << ") -> " << new_node_fitness <<  std::endl;
			}
			outfile << "Shrank all nodes towards " << best_node_index << std::endl;

		}
		//terminate the program

		//report on the state of the simplex
		outfile << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		outfile << "Simplex at end of program, took " << Iterations << " Iterationss." << std::endl;
		output(Iterations, Node_Fitnesses, outfile);
		outfile << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

		sort_nodes(Sorted_node_indexes);
		best_node_index = Sorted_node_indexes[0];
		//Return the best calculated node
		return Simplex[best_node_index];

	}






	///Gradient descent
	void GradientDescentOptimisation::setup_optimisation(const unsigned &n_values){
		//can only be called once
		if(this->n_variables()>0){
			throw OomphLibError("Gradient descent can only be setup once",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}
		//set the number of variables used in problem
		this->n_variables() = n_values;
		//make the current point
		Current_Point.resize(n_values);

		//Get all elements which depend on parameters supplied by the simplex
		add_all_trainable_elements_from_trainable_problems_as_dependents();
	}

	void GradientDescentOptimisation::set_initial_position(const Vector<double>& v_0){
		Current_Point = v_0;
		History_Points.resize(1);
		History_Points[0] = v_0;
	}

	void GradientDescentOptimisation::evaluate_fitness_of_point(Vector<double>& point,
																double& fitness,
																std::ostream &raw_data_outfile){
		// for(unsigned i=0; i<this->n_variables(); i++){
		// 	std::cout << point[i] << std::endl;	
		// }
		
		//push the node to all dependent elements
		push_optimisation_parameter_values_to_dependent_elements(point);
		//condense the residuals from each trainable problem into
		//	a single number and pass it to fitness
		fitness = 0.0;
		for(unsigned i=0; i<Trainable_Problem_Pts.size(); i++){
			fitness += Trainable_Problem_Pts[i]->run();
		}

		//Output the evaluation to the raw data output file
		for(unsigned var=0; var<this->n_variables(); var++){
			raw_data_outfile << point[var] << " ";
		}
		raw_data_outfile << fitness << std::endl;
	}

	void GradientDescentOptimisation::evaluate_fitness_of_current_point(std::ostream &raw_data_outfile){
		evaluate_fitness_of_point(Current_Point,
								Current_Fitness,
								raw_data_outfile);
	}

	//run the gradient descent algorithm until convergence
	Vector<double> GradientDescentOptimisation::run_algorithm(std::ostream &outfile,
													std::ostream &raw_data_outfile){

		outfile << "This is the gradient descent algorithm coordinator...Hello." << std::endl;
		unsigned iterations = 0;

		while(true){
			outfile << std::endl << "Entering iteration " << iterations << std::endl;

			outfile << "Current point: {";
			for(unsigned j=0; j<this->n_variables(); j++){
				outfile << std::setprecision(15) << " " << Current_Point[j];
			}
			outfile << " }" << std::endl;


			//Evaluate current point
			evaluate_fitness_of_current_point(raw_data_outfile);
			History_Fitness.push_back(Current_Fitness);
			outfile << "Current Position fitness " << Current_Fitness << std::endl;
			//Check termination
			if(iterations>1){
				if(std::fabs(History_Fitness[iterations] - History_Fitness[iterations-1])
						<Convergence_Test_Constant){
					outfile << "Converged." << std::endl;
					break;
				}
			}

			//Perform finite difference in all dimensions to find Grad.F
			Vector<double> GradF(this->n_variables(), 0.0);
			outfile << "Calculating gradient:";
			for(unsigned i=0; i<this->n_variables(); i++){
				outfile << " " << i+1 << " of " << this->n_variables();
				GradF[i] = 0.0;
				double perturbed_fitness = 0.0;
				Current_Point[i] += Finite_Difference_Step;

				outfile << "{";
				for(unsigned j=0; j<this->n_variables(); j++){
					outfile << std::setprecision(15) << " " << Current_Point[j];
				}
				outfile << " }";

				evaluate_fitness_of_point(Current_Point,
										perturbed_fitness,
										raw_data_outfile);
				outfile << ". Perturbed fitness: " << perturbed_fitness << std::endl;
				//Calculate the gradient wrt the perturbed quantity;
				GradF[i] = (perturbed_fitness - Current_Fitness)/Finite_Difference_Step;
				Current_Point[i] -= Finite_Difference_Step;
			}
			History_Gradients.push_back(GradF);
			outfile << "Gradient is {";
			for(unsigned i=0; i<this->n_variables(); i++){
				outfile << std::setprecision(15) << " " << GradF[i];
			}
			outfile << " }"<< std::endl;

			///Calculate Gamma
			///	If it is the first iteration use a very small value of gamma,
			///	Maybe the finite difference step length?
			outfile << "Calculating gamma..." << std::endl;
			double gamma;
			// if(iterations==0){
				// gamma = Finite_Difference_Step;
				gamma = 1e-9;
				// outfile << "This is the first iteration so we are using the value " << gamma << std::endl;
			// }
			//Otherwise calculate it assuming the Wolfe conditions apply
			// (dubious, I know...but fun!)
			// else{
			// 	gamma = 0.0;
			// 	double grad_norm = 0.0;
			// 	for(unsigned i=0; i<this->n_variables(); i++){
			// 		gamma += std::fabs((Current_Point[i] - History_Points[iterations-1][i])*
			// 					(GradF[i] - History_Gradients[iterations-1][i]));

			// 		grad_norm += pow(GradF[i] - History_Gradients[iterations-1][i], 2.0);
			// 	}

			// 	gamma /= pow(grad_norm, 0.5);
			// 	outfile << "This is not the first iteration so we are using the value " << gamma << std::endl;
			// }

			//Update the current point
			History_Points.push_back(Current_Point); //Record it first
			for(unsigned i=0; i<this->n_variables(); i++){
				double new_point = 0.0;
				new_point = Current_Point[i] - gamma * GradF[i];
				Current_Point[i] = new_point;
			}

			//Finish the iteration
			iterations++;
		}

		//Not finished, just so it doesn't complain when copmpiling
		Vector<double> DummyVector(1,1.0);
		return DummyVector;
	}
















	//Setup the simplex, construct n+1 optimisation equations each with n internal data
	void ParticleSwarmOptimisation::setup_optimisation(const unsigned &n_variables){
		//can only be called once
		if(this->n_variables()>0){
			throw OomphLibError("Particle Swarm can only be setup once",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}

		//Display a warning
		throw OomphLibWarning("The particle swarm optimisation algorithm requires several parameters for which defaults are provided,\nmake sure you understand and have set them appropriately.",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);

		//If number of particles is 0
		if(N_Swarm == 0)
		{
			throw OomphLibError("N_Swarm is zero",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}

		//If the bounding ball is of radius 0
		if(N_Ball_Radius<1e-12)
		{
			throw OomphLibError("Bounding hypersphere has radius zero",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}

		//If the bounding ball is of radius 0
		if(Max_Iters==0)
		{
			throw OomphLibError("Max iterations is zero",
		                       	OOMPH_CURRENT_FUNCTION,
		                       	OOMPH_EXCEPTION_LOCATION);
		}


		//Get all elements which depend on parameters supplied by the simplex
		add_all_trainable_elements_from_trainable_problems_as_dependents();


		//The record of number of fitness evaluations performed
		Evaluations_performed = 0;

		this->n_variables() = n_variables;

		Current_Swarm.resize(N_Swarm);
		Previous_Swarm.resize(N_Swarm);
		Previous_Fitness.resize(N_Swarm);

		for(unsigned i=0; i<N_Swarm; i++){
			Current_Swarm[i].resize(this->n_variables());
			Previous_Swarm[i].resize(Max_Iters);
			for(unsigned t=0; t<Max_Iters; t++){
				Previous_Swarm[i][t].resize(this->n_variables());
			}
			Previous_Fitness[i].resize(Max_Iters);
		}
	}

	//run the nelder-mead simplex algorithm until convergence
	Vector<double> ParticleSwarmOptimisation::run_algorithm(std::ostream &outfile,
													std::ostream &raw_data_outfile){

		outfile << "This is the Particle Swarm Coordinator..." << std::endl;

		//Set up the swarm, we do this here because we need a random number generator and we may as well use this one
		
		//Setup a mersenne prime rng
		std::mt19937_64 rng;
		// initialize the random number generator with random seed
		std::seed_seq ss{uint32_t(Random_Seed & 0xffffffff), uint32_t(Random_Seed>>32)};
		rng.seed(ss);
		// initialize a uniform distribution between 0 and 1
		std::uniform_real_distribution<double> unif(0, 1);

		//Loop over the particles
		for(unsigned p=0; p<N_Swarm; p++){
			//Construct a random unit vector
			Vector<double> unit_vector(this->n_variables());
			double Vect_Len = 0.0;
			for(unsigned i=0; i< this->n_variables(); i++){
				unit_vector[i] = unif(rng)-0.5;
				Vect_Len += pow(unit_vector[i],2.0);
			}
			//Normalize our vector
			for(unsigned i=0; i< this->n_variables(); i++){
				unit_vector[i]/=std::sqrt(Vect_Len);
			}

			//Now set the particle coordinate to be N_Ball_Centre perturbed by a vector of
			//	length between 0 and N_Ball_Radius
			double R = unif(rng)*N_Ball_Radius;
			for(unsigned i=0; i<this->n_variables(); i++){
				Current_Swarm[p][i] = N_Ball_Centre[i] + R*unit_vector[i];
			}
		}
		//Finish set up swarm


		//The fitness of the current swarm
		Vector<double> Current_Fitness(N_Swarm);

		//The best known position so far for each particle
		std::vector<Vector<double>> Best_Position(N_Swarm);

		Vector<double> Best_Particle_Fitness(N_Swarm);

		//The velocity of each particle in the swarm
		std::vector<Vector<double>> Swarm_Velocity(N_Swarm);
		for(unsigned p=0; p<N_Swarm; p++){
			Swarm_Velocity[p].resize(this->n_variables(),0.0);
		}

		//The best fitness achieved so far
		double Best_Fitness = 1e12;
		//The best known position found so far
		Vector<double> Best_Known_Position;

		//We first evaluate the fitness of all of the particles
		for(unsigned p=0; p<N_Swarm; p++){
			double fitness = 0.0;
			evaluate_fitness_of_particle(p, fitness, raw_data_outfile);
			Current_Fitness[p] = fitness;
			//And we initialize the best known position to this one
			Best_Position[p] = Current_Swarm[p];

			Best_Particle_Fitness[p] = fitness;

			//Get the best known fitness of the swarm
			if(fitness < Best_Fitness){
				Best_Known_Position = Current_Swarm[p];
				Best_Fitness = fitness;
			}
		}

		//Initialize the velocity of the swarm
		for(unsigned p=0; p<N_Swarm; p++){
			for(unsigned v=0; v<this->n_variables(); v++){
				//get two random numbers
				double r1 = unif(rng);
				double r2 = unif(rng);
				Swarm_Velocity[p][v] = W*Swarm_Velocity[p][v]
									 + C1*r1*(Best_Position[p][v] - Current_Swarm[p][v])
									 + C2*r2*(Best_Known_Position[v] - Current_Swarm[p][v]);
			}
		}

		//Perform the loop
		unsigned ITER = 0;
		while(ITER<Max_Iters){
			//Change the parameter values every iteration to aid in convergence
			C1 = (0.5-2.025)*double(ITER)/double(Max_Iters) + 2.025;
			C2 = (0.5-2.025)*double(ITER)/double(Max_Iters) + 2.025;
			W = 1.0 + double(ITER)/double(Max_Iters)*(0.5);

			for(unsigned p=0; p<N_Swarm; p++){
				//Record the position
				Previous_Swarm[p][ITER] = Current_Swarm[p];

				//The new coordinate of the particle	
				Vector<double> New_Position(this->n_variables(), 0.0);

				//get two random numbers
				double r1 = unif(rng);
				double r2 = unif(rng);
				for(unsigned v=0; v<this->n_variables(); v++){
					//Update the particles velocity
					Swarm_Velocity[p][v] = W*Swarm_Velocity[p][v]
									 + C1*r1*(Best_Position[p][v] - Current_Swarm[p][v])
									 + C2*r2*(Best_Known_Position[v] - Current_Swarm[p][v]);	 

					//Calculate what the new position of the particle will be
					New_Position[v] = Current_Swarm[p][v] + Swarm_Velocity[p][v];
				}

				//Update the position of the particle
				//Check to ensure we haven't intersected the N_Ball, doing here for readability
				double R = 0.0;
				for(unsigned v=0; v<this->n_variables(); v++){
					R += pow(New_Position[v] - N_Ball_Centre[v] ,2.0);
				}

				//If the point is outside of the hypersphere
				if(R>N_Ball_Radius*N_Ball_Radius){

					//Move the particle to the intersection point and do nothing more
					//Stuff we need to calculate alpha
					Vector<double> P_k = Current_Swarm[p];
					//The current velocity
					Vector<double> V_k = Swarm_Velocity[p];
					Vector<double> pmxc(this->n_variables(), 0.0); //P-x_c
					double vdotv = 0.0; //V.V
					double pmxcdotpmxc = 0.0; //(P-x_c).(P-x_c)
					double pmxcdotv = 0.0;
					for(unsigned v=0; v<this->n_variables(); v++){
						pmxc[v] = P_k[v] - N_Ball_Centre[v];
						vdotv += V_k[v]*V_k[v];
						pmxcdotpmxc += pmxc[v]*pmxc[v];
						pmxcdotv += pmxc[v]*V_k[v];
					}
					//Fraction of the velocity which takes the particle from the current point to the intersection point
					double alpha = 1.0;
					//Calculate alpha
					alpha = ( -(pmxcdotv) + std::sqrt( pmxcdotv*pmxcdotv - vdotv*(pmxcdotpmxc-pow(N_Ball_Radius,2.0))) )/vdotv;
					if(alpha > 1.0){
						alpha = ( -(pmxcdotv) - std::sqrt( pmxcdotv*pmxcdotv - vdotv*(pmxcdotpmxc-pow(N_Ball_Radius,2.0))) )/vdotv;
					}
					if(alpha < 0.0){
						throw OomphLibError("An error occurred when calculating a particles reflected position",
	                       	OOMPH_CURRENT_FUNCTION,
	                       	OOMPH_EXCEPTION_LOCATION);
					}

					for(unsigned v=0; v<this->n_variables(); v++){
						New_Position[v] = P_k[v]+alpha*V_k[v];
						Swarm_Velocity[p][v] = 0.0;
					}




					//Perform bounces until the particle is placed within the hemisphere
					//The current position
					// Vector<double> P_k = Current_Swarm[p];
					// //The current velocity
					// Vector<double> V_k = Swarm_Velocity[p];
					// //Fraction of the velocity which takes the particle from the current point to the intersection point
					// double alpha = 1.0;
					// //How many times we've gone through this, for debugging
					// unsigned num_bounces = 0.0;
					// do{
					// 	std::cout << "Performing bounce number " << num_bounces << std::endl;
					// 	num_bounces++;
					// 	//Stuff we need to calculate alpha
					// 	Vector<double> pmxc(this->n_variables(), 0.0); //P-x_c
					// 	double vdotv = 0.0; //V.V
					// 	double pmxcdotpmxc = 0.0; //(P-x_c).(P-x_c)
					// 	double pmxcdotv = 0.0;
					// 	for(unsigned v=0; v<this->n_variables(); v++){
					// 		pmxc[v] = P_k[v] - N_Ball_Centre[v];
					// 		vdotv += V_k[v]*V_k[v];
					// 		pmxcdotpmxc += pmxc[v]*pmxc[v];
					// 		pmxcdotv += pmxc[v]*V_k[v];
					// 	}
					// 	//Calculate alpha
					// 	alpha = ( -(pmxcdotv) + std::sqrt( pmxcdotv*pmxcdotv - vdotv*(pmxcdotpmxc-pow(N_Ball_Radius,2.0))) )/vdotv;
					// 	if(alpha > 1.0){
					// 		alpha = ( -(pmxcdotv) - std::sqrt( pmxcdotv*pmxcdotv - vdotv*(pmxcdotpmxc-pow(N_Ball_Radius,2.0))) )/vdotv;
					// 	}
					// 	if(alpha < 0.0){
					// 		throw OomphLibError("An error occurred when calculating a particles reflected position",
		   //                     	OOMPH_CURRENT_FUNCTION,
		   //                     	OOMPH_EXCEPTION_LOCATION);
					// 	}
					// 	//Calculate the new P_k+1 = P_k + alpha*V_k
					// 	for(unsigned v=0; v<this->n_variables(); v++){
					// 		P_k[v] += alpha*V_k[v];
					// 	}
					// 	//Stuff we need for calculating the new velocity vector
					// 	double vdotxcmpk = 0.0;
					// 	double normxcmpk = 0.0;
					// 	for(unsigned v=0; v<this->n_variables(); v++){
					// 		vdotxcmpk += V_k[v]*(N_Ball_Centre[v] - P_k[v]);
					// 		normxcmpk += (N_Ball_Centre[v] - P_k[v])*(N_Ball_Centre[v] - P_k[v]);
					// 	}
					// 	double C = 2*vdotxcmpk/normxcmpk;
					// 	//Calculate the new velocity vector r = (1-alpha)V_k - 2*(1-alpha)*( V_k . (x_c-P_k) ) / |(x_c-P_k)|
					// 	for(unsigned v=0; v<this->n_variables(); v++){
					// 		V_k[v] = (1-alpha)*(V_k[v] - C*(N_Ball_Centre[v] - P_k[v]));
					// 	}
					// 	//Calculate the new position
					// 	R = 0.0;
					// 	for(unsigned v=0; v<this->n_variables(); v++){
					// 		New_Position[v] = P_k[v]+V_k[v];
					// 		R += pow(New_Position[v] - N_Ball_Centre[v] ,2.0);
					// 		Swarm_Velocity[v] = V_k[v];
					// 	}
					// }while(R>N_Ball_Radius*N_Ball_Radius);
				

				}


				//Recalculate the radius of the new position from the centre of the ball					
				R = 0.0;
				for(unsigned v=0; v<this->n_variables(); v++){
					R += pow(New_Position[v] - N_Ball_Centre[v] 
							,2.0);
				}
				if(R>N_Ball_Radius*N_Ball_Radius+1e-6){
					std::cout << "About to die, Nball centre:";
					for(unsigned v=0; v<this->n_variables(); v++){
						std::cout << " " << N_Ball_Centre[v];
					}
					std::cout << std::endl << "N_Ball_Radius: " << N_Ball_Radius << std::endl;
					std::cout << "attempted particle position ";
					for(unsigned v=0; v<this->n_variables(); v++){
						std::cout << " " <<New_Position[v];
					}
					std::cout << std::endl << "Radius of attempted position " << std::sqrt(R) << std::endl;
					throw OomphLibError("A particle has escaped from it's hypersphere cage!",
	                       	OOMPH_CURRENT_FUNCTION,
	                       	OOMPH_EXCEPTION_LOCATION);
				}



				//Set the position to be the new calculated position
				for(unsigned v=0; v<this->n_variables(); v++){
					Current_Swarm[p][v] = New_Position[v];
				}


				//Evaluate the fitness of the particle
				Previous_Fitness[p][ITER] = Current_Fitness[p];
				double fitness = 0.0;
				evaluate_fitness_of_particle(p, fitness, raw_data_outfile);
				Current_Fitness[p] = fitness;
				//Is it the best?
				if(fitness < Best_Particle_Fitness[p]){
					Best_Particle_Fitness[p] = fitness;
					Best_Position[p] = Current_Swarm[p];
				}
			}

			//check for the new best particle
			for(unsigned p=0; p<N_Swarm; p++){
				if(Best_Particle_Fitness[p] < Best_Fitness){
					Best_Fitness = Best_Particle_Fitness[p];
					Best_Known_Position = Current_Swarm[p];
				}
			}

			//Record the best known position
			outfile << ITER << ": ";
			for(unsigned v=0; v<this->n_variables(); v++){
				outfile << Best_Known_Position[v] << " ";
			}
			outfile << ":     " << Best_Fitness << std::endl;

			//Increment the iteration counter
			ITER++;

		}

		//Not finished, just so it doesn't complain when copmpiling
		Vector<double> DummyVector(1,1.0);
		return DummyVector;
		
	}





}	//end namespace