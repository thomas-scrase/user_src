#ifndef OOMPH_FULL_MONODOMAIN_WITH_CELL
#define OOMPH_FULL_MONODOMAIN_WITH_CELL

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

//generic oomph header
#include "generic.h"

//Monodomain elements
#include "../monodomain/monodomain_elements.h"
//DIffusion tensor elements
#include "../monodomain/monodomain_d_expansion_elements.h"
//Cell interface elements (includes cell models)
#include "../cell_interface/cell_interface_elements.h"

namespace oomph
{
	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class QFullMonodomainWithCellElement	:	public virtual QMonodomainElement<DIM, NNODE_1D>,
                                  				public virtual QMonodomainDExpansionElement<DIM, NNODE_1D>,
                                  				public virtual QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>
    {
   		public:
   			QFullMonodomainWithCellElement() : 	QMonodomainElement<DIM, NNODE_1D>(),
												QMonodomainDExpansionElement<DIM, NNODE_1D>(),
												QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>()
												{	}

			//Shuffle the local location of the data
			unsigned diff_min_index_monodomain() const {return (this->u_index_monodomain()+1);}
			unsigned min_index_CellInterfaceEquations() const {return (this->diff_max_index_monodomain());}

			//The number of variables stored at each node
			unsigned required_nvalue(const unsigned &n) const
			{return QMonodomainElement<DIM, NNODE_1D>::required_nvalue(n) + 
					QMonodomainDExpansionElement<DIM, NNODE_1D>::required_nvalue(n) +
					QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>::required_nvalue(n);}


			//the get functions

			//The diffusion matrix
			void get_diff_monodomain(const unsigned& ipt,
									const Vector<double> &s,
									const Vector<double>& x,
									DenseMatrix<double>& D) const
			{
				QMonodomainDExpansionElement<DIM, NNODE_1D>::get_diff_monodomain(ipt, s, x, D);
			}
			//Source
			void get_source_monodomain(const unsigned& ipt,
										const Vector<double>& x,
										double& source) const
			{
				//Set the Vector to hold local coordinates
				Vector<double> s(DIM);
				//Assign values of s
   				for(unsigned i=0;i<DIM;i++) s[i] = this->integral_pt()->knot(ipt,i);

   				//If a source function has been set, use it
				if(QMonodomainElement<DIM, NNODE_1D>::Source_fct_pt!=0){
					//Get source strength
					(*QMonodomainElement<DIM, NNODE_1D>::Source_fct_pt)(x,source);
				}
				//Add the membrane current from the cell model
				source += QCellInterfaceElement<DIM, NUM_VARS, NNODE_1D>::interpolated_membrane_current_CellInterface(s);
			}
			//Membrane potential
			void get_membrane_potential_CellInterface(const unsigned& ipt,
														const Vector<double>& s,
														const Vector<double>& x,
														double& V) const
			{
				V = QMonodomainElement<DIM, NNODE_1D>::interpolated_u_monodomain(s);
			}



			//The output functions
			void output(ostream &outfile) {FiniteElement::output(outfile);}
			void output(ostream &outfile, const unsigned &nplot){
				//vector of local coordinates
				Vector<double> s(DIM);
				Vector<double> xi(DIM);

				// Tecplot header info
				outfile << this->tecplot_zone_string(nplot);

				const unsigned n_node = this->nnode();
				Shape psi(n_node);
				DShape dpsidx(n_node,DIM);

				// Loop over plot points
				unsigned num_plot_points=this->nplot_points(nplot);
				for (unsigned iplot=0;iplot<num_plot_points;iplot++){
					// Get local coordinates of plot point
					this->get_s_plot(iplot,nplot,s);

					// Get the Lagrangian coordinate
					this->interpolated_x(s,xi);

					// Output the position of the plot point
					for(unsigned i=0;i<DIM;i++) {outfile << xi[i] << " ";}

					outfile << this->interpolated_u_monodomain(s) << " ";

					(void)this->dshape_eulerian(s,psi,dpsidx);
					Vector<double> interpolated_dudx(DIM,0.0);
					double dudt = 0.0;
					for(unsigned n=0;n<n_node;n++){
						const double u_ = this->nodal_value(n,this->u_index_monodomain());

						dudt += this->du_dt_monodomain(n)*psi(n);

						for(unsigned i=0;i<DIM;i++){
							interpolated_dudx[i] += u_*dpsidx(n,i);
						}
					}

					outfile << dudt << " ";
					for(unsigned i=0;i<DIM;i++){
						outfile << interpolated_dudx[i]  << " ";
					}


					DenseMatrix<double> Dprint(DIM);
					get_diff_monodomain(0, s, xi, Dprint);
					for(unsigned i=0; i<DIM; i++){
						for(unsigned j=i; j<DIM; j++){
							outfile << Dprint(i,j) << " ";
						}
					}
					outfile << std::endl;
				}

					this->write_tecplot_zone_footer(outfile, nplot);
			}
			void output(FILE* file_pt){FiniteElement::output(file_pt);}
			void output(FILE* file_pt, const unsigned &n_plot){FiniteElement::output(file_pt,n_plot);}
			void output_fct(ostream &outfile, const unsigned &Nplot,FiniteElement::SteadyExactSolutionFctPt exact_soln_pt){FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}
			void output_fct(ostream &outfile, const unsigned &Nplot,const double& time,FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt){FiniteElement::output_fct(outfile,Nplot,time,exact_soln_pt);}
			// void compute_norm(double& el_norm){QUnsteadyHeatElement<DIM,3>::compute_norm(el_norm);}
			void compute_error(ostream &outfile,FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,const double& time,double& error, double& norm){FiniteElement::compute_error(outfile,exact_soln_pt,time,error,norm);}
			void compute_error(ostream &outfile,FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,double& error, double& norm){FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}


			//Residual and Jacobian functions
			void fill_in_contribution_to_residuals(Vector<double> &residuals)
			{
				QMonodomainElement<DIM,NNODE_1D>::fill_in_contribution_to_residuals(residuals);
				QMonodomainDExpansionElement<DIM,NNODE_1D>::fill_in_contribution_to_residuals(residuals);
				QCellInterfaceElement<DIM,NUM_VARS,NNODE_1D>::fill_in_contribution_to_residuals(residuals);
			}
			void fill_in_contribution_to_jacobian(Vector<double> &residuals,DenseMatrix<double> &jacobian)
			{
				FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
			}
			void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals, DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
	  		{
			   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(residuals,jacobian,mass_matrix);
			}
    };

	//GEOMETRIES
	template<unsigned DIM, unsigned NUM_VARS, unsigned NNODE_1D>
	class FaceGeometry<QFullMonodomainWithCellElement<DIM, NUM_VARS, NNODE_1D> > :	public virtual QElement<DIM-1, NNODE_1D>
	{
	public:
		FaceGeometry() : QElement<DIM-1, NNODE_1D>() {}
	};

	template<unsigned NUM_VARS, unsigned NNODE_1D>
	class FaceGeometry<QFullMonodomainWithCellElement<1, NUM_VARS, NNODE_1D> >	:	public virtual PointElement
	{
	public:
		FaceGeometry() : PointElement() {}
	};
}

#endif