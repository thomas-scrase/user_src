#include "cell_interface_elements.h"

namespace oomph
{

	template <unsigned DIM>
	void CellInterfaceEquations<DIM>::output(std::ostream &outfile, const unsigned &nplot){
		// std::cout << "Entered output" << std::endl;
		//Vector of local coordinates
 		Vector<double> s(DIM);

		//Tecplot header info
 		outfile << tecplot_zone_string(nplot);

 		//Get the number of nodes
 		const unsigned n_node = this->nnode();

 		//Preallocate the shape function
 		Shape psi(n_node);

		unsigned num_plot_points=nplot_points(nplot);
 		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
 			// std::cout << "\tOutputting point " << iplot << std::endl;
 			//Get local coordinates of plot point
 			get_s_plot(iplot,nplot,s);

 			//Get Eulerian coordinate of plot point
			Vector<double> x(DIM);
			interpolated_x(s,x);

			//output the Eulerian coordinate
			for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}

			//Get the shape function
			(void)this->shape(s,psi);

			outfile << interpolated_membrane_current_CellInterface(s) << " ";

			//Loop over the variables
			for(unsigned var=min_index_CellInterfaceEquations();var<max_index_CellInterfaceEquations();var++){
		// unsigned var = 8;
				// std::cout << "\t\tvar " << var << std::endl;
				//loop over the nodes
				double interp_var = 0.0;
				for(unsigned n=0;n<n_node;n++){
					//interpolate
					interp_var += nodal_value(n,var)*psi(n);
				}
				//output the variable
				outfile << interp_var << " ";
			}

			//output active strain
			outfile << get_interpolated_cell_active_strain(s) << " ";

			outfile << std::endl;
 		}
 		write_tecplot_zone_footer(outfile,nplot);
	}

	template <unsigned DIM>
	void CellInterfaceEquations<DIM>::output(FILE* file_pt, const unsigned &n_plot){	}


	//====================================================================
	//====================================================================
	//Generic fill in of residual contribution
	//====================================================================
	//====================================================================
	template <unsigned DIM>
	void CellInterfaceEquations<DIM>::fill_in_generic_residual_contribution_cell_interface(Vector<double> &residuals, 
											                                               DenseMatrix<double> &jacobian, 
											                                               DenseMatrix<double> &mass_matrix,
											                                               unsigned flag)
	{
		// std::cout << "Entered cell interface generic residual" << std::endl;
		// std::cout << "entered fill_in_generic_residual_contribution_cell_interface" << std::endl;
		const unsigned n_node = nnode();
		//Preallocate data for the local and global coordinate of the node
		Vector<double> s_node(DIM, 0.0);
		Vector<double> x_node(DIM, 0.0);
		//Preallocate memory for the membrane potential
		double Vm;
		//Preallocate memory for the cell type at the node
		unsigned cell_type;
		//Preallocate memory for the external concentrations
		Vector<double> Ext_conc(3,0.0);
		//Preallocate memory for the strain
		double strain;
		//Preallocate memory for the fibrosis type of the cell
		unsigned fibrosis;
		//Set the value of n_intpt
 		const unsigned n_intpt = integral_pt()->nweight();
		//Set the Vector to hold local coordinates
 		Vector<double> s(DIM);
 		//should the function include contribution from this cell
 		bool compute_this_node;
 		//Storage vectors for the local_index and local_eqn numbers for the single cell data
 		Vector<unsigned> local_eqn(Cell_model_pt->Required_storage());
		Vector<unsigned> local_ind(Cell_model_pt->Required_storage());

		// std::cout << "Finished get variables, entering node loop" << std::endl;
		//Loop over the nodes
		for(unsigned l=0;l<n_node;l++){
			// std::cout << "node:\t" << l << std::endl << "getting local ind and eqn numbers" << std::endl;

			// Get the local ind and local eqn of the cell data
			for(unsigned var=0; var<Cell_model_pt->Required_storage(); var++){
				// std::cout << "var:\t" << var << std::endl;
				// std::cout << "got to var " << var << std::endl;
				local_eqn[var] = nodal_local_eqn(l, min_index_CellInterfaceEquations() + var);
				local_ind[var] = min_index_CellInterfaceEquations() + var;
			}

			// std::cout << "got indexes and eqn numbers storage" << std::endl;
			compute_this_node = true;
			// If Ignore_Repeated_Cells is true then check if the cell has already been computed
			if(Ignore_Repeated_Cells){
				// std::cout << "ignoring repeated cells" << std::endl;
				//loop over the single cell data
				for(unsigned var=0; var<Cell_model_pt->Required_storage(); var++){
					// std::cout << "var:\t" << var << std::endl;
					// Check the value of the residuals corresponding to the single cell data
					if(residuals[local_eqn[var]]!=0){
						// If it's not zero then the cell has already been computed
						// Do not compute this cell
						compute_this_node = false;
						// Stop checking
						break;
					}
				}
			}
			// std::cout << "done ignore cells" << std::endl;
			// If the single cell for this node has already been computed then go to the next node
			if(!compute_this_node){
				std::cout << "+!+!+!+!+!+!+!+!+!+!+!+!+!+! SKIPPED A NODE +!+!+!+!+!+!+!+!+!+!+!+!+!+!" << std::endl;
				for(unsigned printy=0;printy<10;printy++){
					std::cout << std::endl;
				}
				continue;
			}

			// std::cout << "done compute this node" << std::endl;

			//Preallocate memory for the residual and jacobian sub objects
			Vector<double> residual_sub(Cell_model_pt->Required_storage(), 0.0);
			DenseMatrix<double> jacobian_sub(Cell_model_pt->Required_storage(),Cell_model_pt->Required_storage(),0.0);

			// std::cout << "allocated residuals and jacobian" << std::endl;

			local_coordinate_of_node(l,s_node);

			// std::cout << "getting raw nodal position" << std::endl;
			for(unsigned j=0;j<DIM;j++){
				x_node[j] = raw_nodal_position(l,j);
			}
			// std::cout << "got node coords" << std::endl;

			//Get the membrane potential
			// std::cout << "getting potential..." << std::endl;
			get_membrane_potential_CellInterface(0, s_node, x_node, Vm);

			//Get the strain
			// std::cout << "getting strain..." << std::endl;
    		get_strain_CellInterface(0, s_node, x_node, strain);
    		// std::cout << "got potential and strain" << std::endl;

    		//Calculate the external concentrations
			Ext_conc[0] = get_external_Na_conc_CellInterface(0,s_node,x_node);
    		Ext_conc[1] = get_external_Ca_conc_CellInterface(0,s_node,x_node);
    		Ext_conc[2] = get_external_K_conc_CellInterface(0,s_node,x_node);
    		// std::cout << "got concs" << std::endl;

    		//Get the cell type at the node
			cell_type = get_cell_type_at_node_CellInterface(l);

			//get the fibrosis type
			fibrosis = get_fibrosis_type_at_node_CellInterface(l);

			// std::cout << "got cell and fibrosis type" << std::endl;

			//get the cell model to update the residual and jacobian entries
			Cell_model_pt->fill_in_generic_residual_contribution_cell_base(	node_pt(l),
																			Vm,
																			strain,
																			Ext_conc,
																			local_ind,
																			cell_type,
																			mutation_CellInterface(),
																			fibrosis,
																			residual_sub,
																			jacobian_sub,
																			flag);
			// std::cout << "successfully called Cell_model_pt generic residual" << std::endl;
			if(flag){
				throw OomphLibError(
				    "Populate jacobian not implemented yet",
				    OOMPH_CURRENT_FUNCTION,
				    OOMPH_EXCEPTION_LOCATION);
				// for(unsigned test=0; test<Cell_model_pt->Required_storage(); test++){
				// 	jacobian(local_eqn[var], local_eqn[test]) += jacobian_sub(var, test)*W*psi(l);
				// }
			}

			//Loop over the entries in the residual and jacobian
			for(unsigned var=0; var<Cell_model_pt->Required_storage(); var++){
				residuals[local_eqn[var]] += residual_sub[var];//*W*psi(l);
			}
		}
	}





	//====================================================================
	//Force build of templates
	//====================================================================
	template class CellInterfaceEquations<1>;
	template class CellInterfaceEquations<2>;
	template class CellInterfaceEquations<3>;

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
	template class QCellInterfaceElement<1,40,2>;
	template class QCellInterfaceElement<1,40,3>;
	// template class QCellInterfaceElement<1,40,4>;

	template class QCellInterfaceElement<2,40,2>;
	template class QCellInterfaceElement<2,40,3>;
	// template class QCellInterfaceElement<2,40,4>;

	template class QCellInterfaceElement<3,40,2>;
	template class QCellInterfaceElement<3,40,3>;
	// template class QCellInterfaceElement<3,40,4>;
}