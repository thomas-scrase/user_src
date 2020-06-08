#include "cell_interface_elements.h"

namespace oomph
{
	//Only ever output 2 nplot points
	template <unsigned DIM, unsigned NUM_VARS>
	void CellInterfaceEquations<DIM, NUM_VARS>::output(std::ostream &outfile, const unsigned &nplot){
		// std::cout << "BOOM" << std::endl;
		//Vector of local coordinates
 		Vector<double> s(DIM);

		//Tecplot header info
 		outfile << tecplot_zone_string(2);

 		//Get the number of nodes
 		const unsigned n_node = this->nnode();

 		//Preallocate the shape function
 		Shape psi(n_node);

 		for(unsigned l=0; l<n_node; l++){
 			//Get local coordinates of plot point
 			local_coordinate_of_node(l,s);

 			//Get Eulerian coordinate of plot point
			Vector<double> x(DIM);
			interpolated_x(s,x);

			//output the Eulerian coordinate
			for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}

			//Get the shape function
			(void)this->shape(s,psi);

			// outfile << nodal_vm(l) << " ";

			outfile << interpolated_membrane_current_CellInterface(s) << " ";

			Vector<double> custom_output;
			get_nodal_cell_custom_output(l, custom_output);

			for(unsigned i=0; i<custom_output.size(); i++){
				outfile << custom_output[i] << " ";
			}

			//output active strain
			outfile << get_interpolated_cell_active_strain(s) << " ";
			
			//Loop over the variables
			for(unsigned var=min_index_CellInterfaceEquations();var<max_index_CellInterfaceEquations();var++){
				//loop over the nodes
				double interp_var = 0.0;
				for(unsigned n=0;n<n_node;n++){
					//interpolate
					interp_var += nodal_value(n,var)*psi(n);
				}
				//output the variable
				outfile << interp_var << " ";
			}

			outfile << std::endl;
 		}
 		write_tecplot_zone_footer(outfile,nplot);
	}

	template <unsigned DIM, unsigned NUM_VARS>
	void CellInterfaceEquations<DIM, NUM_VARS>::output(FILE* file_pt, const unsigned &n_plot){	}


	//====================================================================
	//====================================================================
	//Generic fill in of residual contribution
	//====================================================================
	//====================================================================
	template <unsigned DIM, unsigned NUM_VARS>
	void CellInterfaceEquations<DIM, NUM_VARS>::fill_in_generic_residual_contribution_cell_interface(Vector<double> &residuals, 
											                                               DenseMatrix<double> &jacobian, 
											                                               DenseMatrix<double> &mass_matrix,
											                                               unsigned flag)
	{
		const unsigned n_node = nnode();
		
		unsigned ipt_node;

 	// 	//should the function include contribution from this cell
 	// 	bool compute_this_node;

 		//Storage vectors for the local_index and local_eqn numbers for the single cell data
 		Vector<int> local_eqn(cell_model_pt()->required_storage());
		//Construct the state container
		CellState state;

		//Loop over the nodes
		for(unsigned l=0;l<n_node;l++){
			ipt_node = this->ipt_at_node(l);
			// Get the local ind and local eqn of the cell data
			for(unsigned var=0; var<cell_model_pt()->required_storage(); var++){
				local_eqn[var] = nodal_local_eqn(l, min_index_CellInterfaceEquations() + var);
			}

			// compute_this_node = true;
			// // If Ignore_Repeated_Cells is true then check if the cell has already been computed
			// if(Ignore_Repeated_Cells){
			// 	//loop over the single cell data
			// 	for(unsigned var=0; var<cell_model_pt()->Required_storage(); var++){

			// 		if(local_eqn[var]>=0){
			// 			// Check the value of the residuals corresponding to the single cell data
			// 			if(residuals[local_eqn[var]]!=0){
			// 				// If it's not zero then the cell has already been computed
			// 				// Do not compute this cell
			// 				compute_this_node = false;
			// 				// Stop checking
			// 				break;
			// 			}
			// 		}
			// 	}
			// }
			// // If the single cell for this node has already been computed then go to the next node
			// if(!compute_this_node){
			// 	std::cout << "+!+!+!+!+!+!+!+!+!+!+!+!+!+! SKIPPED A NODE +!+!+!+!+!+!+!+!+!+!+!+!+!+!" << std::endl;
			// 	for(unsigned printy=0;printy<10;printy++){
			// 		std::cout << std::endl;
			// 	}
			// 	continue;
			// }


			//Preallocate memory for the residual and jacobian sub objects
			Vector<double> residual_sub(cell_model_pt()->required_storage(), 0.0);
			DenseMatrix<double> jacobian_sub(cell_model_pt()->required_storage(),cell_model_pt()->required_storage(),0.0);

			fill_state_container_at_node(state, l);
			//get the cell model to update the residual and jacobian entries
			cell_model_pt()->fill_in_generic_residual_contribution_cell_base(state,
																			residual_sub,
																			jacobian_sub,
																			flag);
			// Loop over the entries in the residual and jacobian
			for(unsigned var=0; var<cell_model_pt()->required_storage(); var++){
				if(local_eqn[var]>=0){
					residuals[local_eqn[var]] += residual_sub[var];


					if(flag){
						for(unsigned var1=0; var1<cell_model_pt()->required_storage(); var1++){
							jacobian(local_eqn[var], local_eqn[var1]) += jacobian_sub(var, var1);
						}
						// jacobian(local_eqn[var], local_eqn[var]) += jacobian_sub(var, var);
					}
				}
			}
		}
	}







	//====================================================================
	//Force build of templates
	//====================================================================
	template class CellInterfaceEquations<1,45>;
	template class CellInterfaceEquations<2,45>;
	template class CellInterfaceEquations<3,45>;

	template class CellInterfaceEquations<1,25>;
	template class CellInterfaceEquations<2,25>;
	template class CellInterfaceEquations<3,25>;

	template class CellInterfaceEquations<1,1>;
	template class CellInterfaceEquations<2,1>;
	template class CellInterfaceEquations<3,1>;

	//====================================================================
	//Force build of Q - templates
	//====================================================================

	//!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+
	//IF YOU ADD A NEW CELL MODEL:
	//		Copy and paste the commented-out commands below, edit the last number (0)
	//		in the <x,y,z> to reflect the Required_Storage of the new cell model
	//!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+

	// template class QCellInterfaceElement<1,2>;
	// template class QCellInterfaceElement<1,0,3>;
	// template class QCellInterfaceElement<1,0,4>;

	// template class QCellInterfaceElement<2,0,2>;
	// template class QCellInterfaceElement<2,0,3>;
	// template class QCellInterfaceElement<2,0,4>;

	// template class QCellInterfaceElement<3,0,2>;
	// template class QCellInterfaceElement<3,0,3>;
	// template class QCellInterfaceElement<3,0,4>;



	
	//Force build for the CNZ cell model
	template class QCellInterfaceElement<1,45,2>;
	template class QCellInterfaceElement<1,45,3>;
	// template class QCellInterfaceElement<1,45,4>;

	template class QCellInterfaceElement<2,45,2>;
	template class QCellInterfaceElement<2,45,3>;
	// template class QCellInterfaceElement<2,45,4>;

	template class QCellInterfaceElement<3,45,2>;
	template class QCellInterfaceElement<3,45,3>;
	// template class QCellInterfaceElement<3,45,4>;


	template class TCellInterfaceElement<1,45,2>;
	template class TCellInterfaceElement<1,45,3>;
	// template class TCellInterfaceElement<1,45,4>;

	template class TCellInterfaceElement<2,45,2>;
	template class TCellInterfaceElement<2,45,3>;
	// template class TCellInterfaceElement<2,45,4>;

	template class TCellInterfaceElement<3,45,2>;
	template class TCellInterfaceElement<3,45,3>;
	// template class TCellInterfaceElement<3,45,4>;



	//Force build for the TNNP cell model
	template class QCellInterfaceElement<1,25,2>;
	template class QCellInterfaceElement<1,25,3>;
	// template class QCellInterfaceElement<1,25,4>;

	template class QCellInterfaceElement<2,25,2>;
	template class QCellInterfaceElement<2,25,3>;
	// template class QCellInterfaceElement<2,25,4>;

	template class QCellInterfaceElement<3,25,2>;
	template class QCellInterfaceElement<3,25,3>;
	// template class QCellInterfaceElement<3,25,4>;


	template class TCellInterfaceElement<1,25,2>;
	template class TCellInterfaceElement<1,25,3>;
	// template class TCellInterfaceElement<1,25,4>;

	template class TCellInterfaceElement<2,25,2>;
	template class TCellInterfaceElement<2,25,3>;
	// template class TCellInterfaceElement<2,25,4>;

	template class TCellInterfaceElement<3,25,2>;
	template class TCellInterfaceElement<3,25,3>;
	// template class TCellInterfaceElement<3,25,4>;




	//Force build for the 1D variable cell models
	template class QCellInterfaceElement<1,1,2>;
	template class QCellInterfaceElement<1,1,3>;
	// template class QCellInterfaceElement<1,1,4>;

	template class QCellInterfaceElement<2,1,2>;
	template class QCellInterfaceElement<2,1,3>;
	// template class QCellInterfaceElement<2,1,4>;

	template class QCellInterfaceElement<3,1,2>;
	template class QCellInterfaceElement<3,1,3>;
	// template class QCellInterfaceElement<3,1,4>;


	template class TCellInterfaceElement<1,1,2>;
	template class TCellInterfaceElement<1,1,3>;
	// template class TCellInterfaceElement<1,1,4>;

	template class TCellInterfaceElement<2,1,2>;
	template class TCellInterfaceElement<2,1,3>;
	// template class TCellInterfaceElement<2,1,4>;

	template class TCellInterfaceElement<3,1,2>;
	template class TCellInterfaceElement<3,1,3>;
	// template class TCellInterfaceElement<3,1,4>;
}