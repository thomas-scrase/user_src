#include "conducting_cell_elements.h"

#include "../cell_model_updated/TNNP_06_Updated.h"
#include "../cell_model_updated/validate_cell.h"
#include "../cell_model_updated/validate_NSpecies_Lotka_Volterra.h"
#include "../cell_model_updated/TomekORudy_Updated.h"
#include "../cell_model_updated/TNNP_06_MC_IKr.h"


namespace oomph
{
	//Only ever output 2 nplot points
	// template <class CELL_MODEL, class CONDUCTANCE_MODEL>
	// void ConductingCellEquations<DIM, CELL_MODEL>::output(std::ostream &outfile, const unsigned &nplot){
		// // std::cout << "BOOM" << std::endl;
		// //Vector of local coordinates
 	// 	Vector<double> s(DIM);

		// //Tecplot header info
 	// 	outfile << tecplot_zone_string(2);

 	// 	//Get the number of nodes
 	// 	const unsigned n_node = this->nnode();

 	// 	//Preallocate the shape function
 	// 	Shape psi(n_node);

 	// 	for(unsigned l=0; l<n_node; l++){
 	// 		//Get local coordinates of plot point
 	// 		local_coordinate_of_node(l,s);

 	// 		//Get Eulerian coordinate of plot point
		// 	Vector<double> x(DIM);
		// 	interpolated_x(s,x);

		// 	//output the Eulerian coordinate
		// 	for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}

		// 	//Get the shape function
		// 	(void)this->shape(s,psi);

		// 	outfile << get_nodal_membrane_potential(l) << " ";

		// 	outfile << get_nodal_membrane_current(l) << " ";

		// 	//output active strain
		// 	outfile << get_nodal_active_stress(l) << " ";

			
		// 	Vector<double> custom_output;
		// 	get_nodal_cell_custom_output(l, custom_output);
		// 	for(unsigned i=0; i<custom_output.size(); i++){
		// 		outfile << custom_output[i] << " ";
		// 	}
			
		// 	//Loop over the variables
		// 	for(unsigned var=min_index_CellInterfaceEquations();var<max_index_CellInterfaceEquations();var++){
		// 		//loop over the nodes
		// 		double interp_var = 0.0;
		// 		for(unsigned n=0;n<n_node;n++){
		// 			//interpolate
		// 			interp_var += nodal_value(n,var)*psi(n);
		// 		}
		// 		//output the variable
		// 		outfile << interp_var << " ";
		// 	}

		// 	outfile << std::endl;
 	// 	}
 	// 	write_tecplot_zone_footer(outfile,nplot);
	// }

	// template <unsigned DIM, class CELL_MODEL>
	// void CellInterfaceEquations<DIM, CELL_MODEL>::output(FILE* file_pt, const unsigned &n_plot){	}

	// template<class CLASS>
	// std::unordered_map<std::string, double> forwarder(void* context, const double& t)
	// {
	// 	return (static_cast<CLASS*>(context)->get_other_variables(t));
	// }

	// template class FastSingleCellUpdated<TNNP06MCIKr>;

	// template class FastSingleCellUpdated<ExplicitTNNP06VentUpdated>;

	// template class FastSingleCellUpdated<TomekORudyVentUpdated>;

	// template class FastSingleCellUpdated<ValidateCell>;

	// template class FastSingleCellUpdated<ValidateCellLotkaVolterra>;


	// template<class CELL_MODEL>
	// template class ConductingCellEquations<CELL_MODEL, MonodomainEquations<1>>;


	// template<class CELL_MODEL>
	// template class ConductingCellEquations<CELL_MODEL, MonodomainEquations<2>>;


	// template<class CELL_MODEL>
	// template class ConductingCellEquations<CELL_MODEL, MonodomainEquations<3>>;
}