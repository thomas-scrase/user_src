//LIC// ====================================================================
//LIC// This file contains a custom equation and element type for oomph-lib
//LIC// Solves the monodomain equation, based on the gen_advection_diffusion_elements
//LIC//
//LIC// Main changes:
//LIC// + Convection terms are removed in their entirety for performance 
//LIC// + A parameter theta is added to allow for Crank-Nicholson in computation of
//LIC//     Jacobian and residuals
//LIC// + Number of storage values at each node is increased from 1 to 10
//LIC//     to accommodate for the 3, 3D unit vectors of fibre orientation
//LIC// 
//LIC//====================================================================


//Non-inline functions for monodomain Diffusion elements
#include "diffusion_coefficient_expansion_elements.h"

namespace oomph
{                                                                                                                           //0.5 Crank Nicolson

//======================================================================
/// \short Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template <unsigned DIM, unsigned NVECT>
void  DiffusionCoefficientExpansionEquations<DIM, NVECT>::
fill_in_generic_residual_contribution_diffusion_coefficient_expansion(Vector<double> &residuals, 
                                               DenseMatrix<double> &jacobian, 
                                               DenseMatrix<double> 
                                               &mass_matrix,
                                               unsigned flag) 
{
  // Loop over the test functions
  const unsigned n_node = nnode();
  int local_eqn, local_unknown = 0;
    for(unsigned l=0;l<n_node;l++){
        for(unsigned index = diffusion_coefficient_min_index_expansion(); index < diffusion_coefficient_max_index_expansion(); index++){
          local_eqn = nodal_local_eqn(l, index);
          //IF it's not a boundary condition
          if(local_eqn >= 0){
              //residual is just the fibre vector component
              residuals[local_eqn] = 0;

              //Jacobian
              if(flag){
                for(unsigned l2=0;l2<n_node;l2++){
                    local_unknown = nodal_local_eqn(l2,index);
                    if(local_eqn==local_unknown){
                      jacobian(local_eqn,local_unknown) = 1;
                    }
                    else{
                      jacobian(local_eqn,local_unknown) = 0;
                    }
                }
              }
              
          }
        }//End handle diffusion tensor variables
    }
}


//!!!!!!REMOVE???
//======================================================================
/// \short Output function:
///
///   x,y,u,w_x,w_y   or    x,y,z,u,w_x,w_y,w_z
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NVECT>
void  DiffusionCoefficientExpansionEquations<DIM, NVECT>::output(std::ostream &outfile, 
                                               const unsigned &nplot)
{ 
 //  // std::cout << "I am plotting with this plot" << std::endl;
 // //Vector of local coordinates
 // Vector<double> s(DIM);

 
 // // Tecplot header info
 // outfile << tecplot_zone_string(nplot);
 
 // const unsigned n_node = this->nnode();
 // Shape psi(n_node);
 // DShape dpsidx(n_node,DIM);

 // // Loop over plot points
 // unsigned num_plot_points=nplot_points(nplot);
 // for (unsigned iplot=0;iplot<num_plot_points;iplot++)
 //  {
 //   // Get local coordinates of plot point
 //   get_s_plot(iplot,nplot,s);
   
 //   // Get Eulerian coordinate of plot point
 //   Vector<double> x(DIM);
 //   interpolated_x(s,x);
   
 //   for(unsigned i=0;i<DIM;i++) 
 //    {
 //     outfile << x[i] << " ";
 //    }
   
 //   //Get the gradients
 //   (void)this->dshape_eulerian(s,psi,dpsidx);

 //    DenseMatrix<double> PrintD(DIM, DIM);
 //    get_diff_monodomain(0, s, x, PrintD);

 //    for(unsigned i=0; i<DIM; i++){
 //      for(unsigned j=0; j<DIM; j++){
 //        outfile << PrintD(i,j) << " ";
 //      }
 //    }

 //    outfile  << std::endl;
   
 //  }

 // // Write tecplot footer (e.g. FE connectivity lists)
 // write_tecplot_zone_footer(outfile,nplot);

}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NVECT>
void DiffusionCoefficientExpansionEquations<DIM, NVECT>::output(FILE* file_pt,
                                              const unsigned &n_plot)
{ }



//======================================================================
 /// \short  Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM, unsigned NVECT>
void DiffusionCoefficientExpansionEquations<DIM, NVECT>::output_fct(std::ostream &outfile, 
             const unsigned &nplot, 
             FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {  }

//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<unsigned DIM, unsigned NVECT, unsigned NNODE_1D>
const unsigned QDiffusionCoefficientExpansionElement<DIM, NVECT,NNODE_1D>::Initial_Nvalue = NVECT;
                                                                  //D

//====================================================================
// Force build of templates with number of vectors up to the dimension
//====================================================================
template class DiffusionCoefficientExpansionEquations<1,1>;

template class DiffusionCoefficientExpansionEquations<2,1>;
template class DiffusionCoefficientExpansionEquations<2,2>;

template class DiffusionCoefficientExpansionEquations<3,1>;
template class DiffusionCoefficientExpansionEquations<3,2>;
template class DiffusionCoefficientExpansionEquations<3,3>;


template class QDiffusionCoefficientExpansionElement<1,1,2>;
template class QDiffusionCoefficientExpansionElement<1,1,3>;
template class QDiffusionCoefficientExpansionElement<1,1,4>;


template class QDiffusionCoefficientExpansionElement<2,1,2>;
template class QDiffusionCoefficientExpansionElement<2,2,2>;
template class QDiffusionCoefficientExpansionElement<2,1,3>;
template class QDiffusionCoefficientExpansionElement<2,2,3>;
template class QDiffusionCoefficientExpansionElement<2,1,4>;
template class QDiffusionCoefficientExpansionElement<2,2,4>;


template class QDiffusionCoefficientExpansionElement<3,1,2>;
template class QDiffusionCoefficientExpansionElement<3,2,2>;
template class QDiffusionCoefficientExpansionElement<3,3,2>;
template class QDiffusionCoefficientExpansionElement<3,1,3>;
template class QDiffusionCoefficientExpansionElement<3,2,3>;
template class QDiffusionCoefficientExpansionElement<3,3,3>;
template class QDiffusionCoefficientExpansionElement<3,1,4>;
template class QDiffusionCoefficientExpansionElement<3,2,4>;
template class QDiffusionCoefficientExpansionElement<3,3,4>;


/////////////////////////////////////////////////////////////////////////
// TDiffusionCoefficientExpansionElement
/////////////////////////////////////////////////////////////////////////



//======================================================================
// Set the data for the number of Variables at each node, always 1
//======================================================================
template<unsigned DIM, unsigned NVECT, unsigned NNODE_1D>
const unsigned TDiffusionCoefficientExpansionElement<DIM, NVECT,NNODE_1D>::Initial_Nvalue = NVECT;
                                                                  //D

//====================================================================
// Force build of templates with number of vectors up to the dimension
//====================================================================
template class TDiffusionCoefficientExpansionElement<1,1,2>;
template class TDiffusionCoefficientExpansionElement<1,1,3>;
template class TDiffusionCoefficientExpansionElement<1,1,4>;

template class TDiffusionCoefficientExpansionElement<2,1,2>;
template class TDiffusionCoefficientExpansionElement<2,2,2>;
template class TDiffusionCoefficientExpansionElement<2,1,3>;
template class TDiffusionCoefficientExpansionElement<2,2,3>;
template class TDiffusionCoefficientExpansionElement<2,1,4>;
template class TDiffusionCoefficientExpansionElement<2,2,4>;

template class TDiffusionCoefficientExpansionElement<3,1,2>;
template class TDiffusionCoefficientExpansionElement<3,2,2>;
template class TDiffusionCoefficientExpansionElement<3,3,2>;
template class TDiffusionCoefficientExpansionElement<3,1,3>;
template class TDiffusionCoefficientExpansionElement<3,2,3>;
template class TDiffusionCoefficientExpansionElement<3,3,3>;
}
