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
#include "monodomain_d_expansion_elements.h"

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
template <unsigned DIM>
void  MonodomainDExpansionEquations<DIM>::
fill_in_generic_residual_contribution_monodomain_d_expansion(Vector<double> &residuals, 
                                               DenseMatrix<double> &jacobian, 
                                               DenseMatrix<double> 
                                               &mass_matrix,
                                               unsigned flag) 
{
  // Loop over the test functions
  const unsigned n_node = nnode();
  int local_eqn, local_unknown = 0;
    for(unsigned l=0;l<n_node;l++){
      // std::cout << "node " << l << std::endl;
        //fill in residual and jacobian contribution for the diffusion tensor
        //The monodomain equation is incapable of changing fibre orientation,
        //  once coupled with a solid element the right cauchy-green deformation tensor
        //  is used to push forward the interpolated fibre direction which is handled
        //  by overloading the diff_fct_pt
        // std::cout << "diff mat" << std::endl;
        for(unsigned index = diff_min_index_monodomain(); index < diff_max_index_monodomain(); index++){
          local_eqn = nodal_local_eqn(l, index);
          //IF it's not a boundary condition
          if(local_eqn >= 0){
              //residual is just the fibre vector component
              residuals[local_eqn] = 0;

              //Jacobian
              if(flag){
                for(unsigned l2=0;l2<n_node;l2++){
                    local_unknown = nodal_local_eqn(l2,index);        ///THIS WAS u_nodal_index, but surely I should be looping over the diffusion indexes?
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
template <unsigned DIM>
void  MonodomainDExpansionEquations<DIM>::output(std::ostream &outfile, 
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
template <unsigned DIM>
void MonodomainDExpansionEquations<DIM>::output(FILE* file_pt,
                                              const unsigned &nplot)
{ }



//======================================================================
 /// \short  Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM>
void MonodomainDExpansionEquations<DIM>::output_fct(std::ostream &outfile, 
             const unsigned &nplot, 
             FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {  }

//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
const unsigned QMonodomainDExpansionElement<DIM,NNODE_1D>::Initial_Nvalue = DIM*DIM;
                                                                  //D

//====================================================================
// Force build of templates
//====================================================================
template class MonodomainDExpansionEquations<1>;
template class MonodomainDExpansionEquations<2>;
template class MonodomainDExpansionEquations<3>;

template class QMonodomainDExpansionElement<1,2>;
template class QMonodomainDExpansionElement<1,3>;
template class QMonodomainDExpansionElement<1,4>;

template class QMonodomainDExpansionElement<2,2>;
template class QMonodomainDExpansionElement<2,3>;
template class QMonodomainDExpansionElement<2,4>;

template class QMonodomainDExpansionElement<3,2>;
template class QMonodomainDExpansionElement<3,3>;
template class QMonodomainDExpansionElement<3,4>;


/////////////////////////////////////////////////////////////////////////
// TMonodomainDExpansionElement
/////////////////////////////////////////////////////////////////////////



//======================================================================
// Set the data for the number of Variables at each node, always 1
//======================================================================
template<unsigned DIM, unsigned NNODE_1D>
const unsigned TMonodomainDExpansionElement<DIM,NNODE_1D>::Initial_Nvalue = DIM*DIM;
                                                                  //D

//====================================================================
// Force build of templates
//====================================================================
template class TMonodomainDExpansionElement<1,2>;
template class TMonodomainDExpansionElement<1,3>;
template class TMonodomainDExpansionElement<1,4>;

template class TMonodomainDExpansionElement<2,2>;
template class TMonodomainDExpansionElement<2,3>;
template class TMonodomainDExpansionElement<2,4>;

template class TMonodomainDExpansionElement<3,2>;
template class TMonodomainDExpansionElement<3,3>;


}
