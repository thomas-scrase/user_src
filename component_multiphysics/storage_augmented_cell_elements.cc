#include "storage_augmented_cell_elements.h"

namespace oomph{

	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	void QStorageAugmentedCellElement<DIM, NUM_VARS, NNODE_1D>::output(std::ostream &outfile, const unsigned &nplot){
		//Tecplot header info
 		outfile << this->tecplot_zone_string(nplot);

 		//Vector of local coordinates
 		Vector<double> s(DIM);

 		//Get the number of nodes
 		const unsigned n_node = this->nnode();

 		//Preallocate the shape function
 		Shape psi(n_node);

 		unsigned num_plot_points=this->nplot_points(nplot);
 		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
 			//Get local coordinates of plot point
 			this->get_s_plot(iplot,nplot,s);

 			//Get Eulerian coordinate of plot point
			Vector<double> x(DIM);
			this->interpolated_x(s,x);

			//output the Eulerian coordinate
			for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}

			//Get the shape function
			(void)this->shape(s,psi);

			outfile << this->interpolated_membrane_current_CellInterface(s) << " ";


			DenseMatrix<double> diffusion_matrix = this->get_interpolated_diffusion_matrix_augmented_cell(s);

			for(unsigned i = 0; i < DIM; i++){
				for(unsigned j = i; j < DIM; j++){
					outfile << diffusion_matrix(i,j) << " ";
				}
			}

			outfile << std::endl;
 		}
 		this->write_tecplot_zone_footer(outfile,nplot);
	}


	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	void TStorageAugmentedCellElement<DIM, NUM_VARS, NNODE_1D>::output(std::ostream &outfile, const unsigned &nplot){
		//Tecplot header info
 		outfile << this->tecplot_zone_string(nplot);

 		//Vector of local coordinates
 		Vector<double> s(DIM);

 		//Get the number of nodes
 		const unsigned n_node = this->nnode();

 		//Preallocate the shape function
 		Shape psi(n_node);

 		unsigned num_plot_points=this->nplot_points(nplot);
 		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
 			//Get local coordinates of plot point
 			this->get_s_plot(iplot,nplot,s);

 			//Get Eulerian coordinate of plot point
			Vector<double> x(DIM);
			this->interpolated_x(s,x);

			//output the Eulerian coordinate
			for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}

			//Get the shape function
			(void)this->shape(s,psi);

			outfile << this->interpolated_membrane_current_CellInterface(s) << " ";


			DenseMatrix<double> diffusion_matrix = this->get_interpolated_diffusion_matrix_augmented_cell(s);

			for(unsigned i = 0; i < DIM; i++){
				for(unsigned j = i; j < DIM; j++){
					outfile << diffusion_matrix(i,j) << " ";
				}
			}

			outfile << std::endl;
 		}
 		this->write_tecplot_zone_footer(outfile,nplot);
	}

	template<unsigned DIM, unsigned NUM_VARS>
	void PointStorageAugmentedCellElement<DIM, NUM_VARS>::output(std::ostream &outfile, const unsigned &nplot){
		//Tecplot header info
 		outfile << this->tecplot_zone_string(nplot);

 		//Vector of local coordinates
 		Vector<double> s(DIM);

 		//Get the number of nodes
 		const unsigned n_node = this->nnode();

 		//Preallocate the shape function
 		Shape psi(n_node);

 		unsigned num_plot_points=this->nplot_points(nplot);
 		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
 			//Get local coordinates of plot point
 			this->get_s_plot(iplot,nplot,s);

 			//Get Eulerian coordinate of plot point
			Vector<double> x(DIM);
			this->interpolated_x(s,x);

			//output the Eulerian coordinate
			for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}

			//Get the shape function
			(void)this->shape(s,psi);

			outfile << this->interpolated_membrane_current_CellInterface(s) << " ";


			DenseMatrix<double> diffusion_matrix = this->get_interpolated_diffusion_matrix_augmented_cell(s);

			for(unsigned i = 0; i < DIM; i++){
				for(unsigned j = i; j < DIM; j++){
					outfile << diffusion_matrix(i,j) << " ";
				}
			}

			outfile << std::endl;
 		}
 		this->write_tecplot_zone_footer(outfile,nplot);
	}


























	template class QStorageAugmentedCellElement<1,1,2>;
	template class QStorageAugmentedCellElement<1,1,3>;
	template class QStorageAugmentedCellElement<2,1,2>;
	template class QStorageAugmentedCellElement<2,1,3>;
	template class QStorageAugmentedCellElement<3,1,2>;
	template class QStorageAugmentedCellElement<3,1,3>;

	template class TStorageAugmentedCellElement<1,1,2>;
	template class TStorageAugmentedCellElement<1,1,3>;
	template class TStorageAugmentedCellElement<2,1,2>;
	template class TStorageAugmentedCellElement<2,1,3>;
	template class TStorageAugmentedCellElement<3,1,2>;
	template class TStorageAugmentedCellElement<3,1,3>;

	template class PointStorageAugmentedCellElement<1,1>;
	template class PointStorageAugmentedCellElement<2,1>;
	template class PointStorageAugmentedCellElement<3,1>;


	template class QStorageAugmentedCellElement<1,40,2>;
	template class QStorageAugmentedCellElement<1,40,3>;
	template class QStorageAugmentedCellElement<2,40,2>;
	template class QStorageAugmentedCellElement<2,40,3>;
	template class QStorageAugmentedCellElement<3,40,2>;
	template class QStorageAugmentedCellElement<3,40,3>;

	template class TStorageAugmentedCellElement<1,40,2>;
	template class TStorageAugmentedCellElement<1,40,3>;
	template class TStorageAugmentedCellElement<2,40,2>;
	template class TStorageAugmentedCellElement<2,40,3>;
	template class TStorageAugmentedCellElement<3,40,2>;
	template class TStorageAugmentedCellElement<3,40,3>;

	template class PointStorageAugmentedCellElement<1,40>;
	template class PointStorageAugmentedCellElement<2,40>;
	template class PointStorageAugmentedCellElement<3,40>;
}