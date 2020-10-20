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

			outfile << get_nodal_membrane_potential(l) << " ";

			outfile << get_nodal_membrane_current(l) << " ";

			//output active strain
			outfile << get_nodal_active_stress(l) << " ";

			
			Vector<double> custom_output;
			get_nodal_cell_custom_output(l, custom_output);
			for(unsigned i=0; i<custom_output.size(); i++){
				outfile << custom_output[i] << " ";
			}
			
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
		//Get the number of nodes in the element
		const unsigned n_node = nnode();
 		//Preallocate storage vectors for the local_index and local_eqn numbers for the single cell data
 		Vector<int> local_eqn(cell_model_pt()->required_nodal_variables());
		//Preallocate the state container
		CellState state;

		//Loop over the nodes: we don't loop over integration points
		//	because cells are the nodes. If the cells are not points then
		//	this is reflected by how data passed to it is calculated.
		//Note that the cell model is still a point, these cell models are
		//	always points, but by setting cells are not points means any
		//	external "stimulus" to the cell model is calculated from integrating
		//	over the element.
		for(unsigned l=0;l<n_node;l++){

			//If cells are represented by a point and this cell should not
			//	be computed then skip it
			if(Cells_Are_Points && !Cell_Inds_To_Compute[l]){continue;}

			// Get the local ind and local eqn of the cell data
			for(unsigned var=0; var<cell_model_pt()->required_nodal_variables(); var++){
				local_eqn[var] = nodal_local_eqn(l, min_index_CellInterfaceEquations() + var);
			}

			//Preallocate memory for the residual and jacobian sub objects
			Vector<double> residual_sub(cell_model_pt()->required_nodal_variables(), 0.0);
			DenseMatrix<double> jacobian_sub(cell_model_pt()->required_nodal_variables(),cell_model_pt()->required_nodal_variables(),0.0);

			fill_state_container_at_node(l, state);
			//get the cell model to update the residual and jacobian entries
			cell_model_pt()->fill_in_generic_residual_contribution_cell_base(state,
																			residual_sub,
																			jacobian_sub,
																			flag);
			// Loop over the entries in the residual and jacobian
			for(unsigned var=0; var<cell_model_pt()->required_nodal_variables(); var++){
				if(local_eqn[var]>=0){
					residuals[local_eqn[var]] += residual_sub[var];
					//Compute Jacobian?
					if(flag){
						for(unsigned var1=0; var1<cell_model_pt()->required_nodal_variables(); var1++){
							jacobian(local_eqn[var], local_eqn[var1]) += jacobian_sub(var, var1);
						}
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

	template class CellInterfaceEquations<1,49>;
	template class CellInterfaceEquations<2,49>;
	template class CellInterfaceEquations<3,49>;

	template class CellInterfaceEquations<1,25>;
	template class CellInterfaceEquations<2,25>;
	template class CellInterfaceEquations<3,25>;

	template class CellInterfaceEquations<1,34>;
	template class CellInterfaceEquations<2,34>;
	template class CellInterfaceEquations<3,34>;

	template class CellInterfaceEquations<1,1>;
	template class CellInterfaceEquations<2,1>;
	template class CellInterfaceEquations<3,1>;

	//====================================================================
	//Force build of Q - templates
	//====================================================================

	//!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+
	//IF YOU ADD A NEW CELL MODEL:
	//		Copy and paste the commented-out commands below, edit the last number (0)
	//		in the <x,y,z> to reflect the required_nodal_variables of the new cell model
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

	template class PointCellInterfaceElement<1,45>;
	template class PointCellInterfaceElement<2,45>;
	template class PointCellInterfaceElement<3,45>;




	//Force build for the Explicit CNZ cell model
	template class QCellInterfaceElement<1,49,2>;
	template class QCellInterfaceElement<1,49,3>;
	// template class QCellInterfaceElement<1,49,4>;

	template class QCellInterfaceElement<2,49,2>;
	template class QCellInterfaceElement<2,49,3>;
	// template class QCellInterfaceElement<2,49,4>;

	template class QCellInterfaceElement<3,49,2>;
	template class QCellInterfaceElement<3,49,3>;
	// template class QCellInterfaceElement<3,49,4>;


	template class TCellInterfaceElement<1,49,2>;
	template class TCellInterfaceElement<1,49,3>;
	// template class TCellInterfaceElement<1,49,4>;

	template class TCellInterfaceElement<2,49,2>;
	template class TCellInterfaceElement<2,49,3>;
	// template class TCellInterfaceElement<2,49,4>;

	template class TCellInterfaceElement<3,49,2>;
	template class TCellInterfaceElement<3,49,3>;
	// template class TCellInterfaceElement<3,49,4>;

	template class PointCellInterfaceElement<1,49>;
	template class PointCellInterfaceElement<2,49>;
	template class PointCellInterfaceElement<3,49>;




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

	template class PointCellInterfaceElement<1,25>;
	template class PointCellInterfaceElement<2,25>;
	template class PointCellInterfaceElement<3,25>;


	//Force build for the TNNP with Rice myofilament cell model
	template class QCellInterfaceElement<1,34,2>;
	template class QCellInterfaceElement<1,34,3>;
	// template class QCellInterfaceElement<1,34,4>;

	template class QCellInterfaceElement<2,34,2>;
	template class QCellInterfaceElement<2,34,3>;
	// template class QCellInterfaceElement<2,34,4>;

	template class QCellInterfaceElement<3,34,2>;
	template class QCellInterfaceElement<3,34,3>;
	// template class QCellInterfaceElement<3,34,4>;


	template class TCellInterfaceElement<1,34,2>;
	template class TCellInterfaceElement<1,34,3>;
	// template class TCellInterfaceElement<1,34,4>;

	template class TCellInterfaceElement<2,34,2>;
	template class TCellInterfaceElement<2,34,3>;
	// template class TCellInterfaceElement<2,34,4>;

	template class TCellInterfaceElement<3,34,2>;
	template class TCellInterfaceElement<3,34,3>;
	// template class TCellInterfaceElement<3,34,4>;

	template class PointCellInterfaceElement<1,34>;
	template class PointCellInterfaceElement<2,34>;
	template class PointCellInterfaceElement<3,34>;




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

	template class PointCellInterfaceElement<1,1>;
	template class PointCellInterfaceElement<2,1>;
	template class PointCellInterfaceElement<3,1>;


	template class MonodomainSingleCellElement<1>;
	template class MonodomainSingleCellElement<25>;
	template class MonodomainSingleCellElement<34>;
	template class MonodomainSingleCellElement<45>;
	template class MonodomainSingleCellElement<49>;
}