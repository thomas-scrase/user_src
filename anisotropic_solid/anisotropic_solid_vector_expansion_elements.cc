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
#include "anisotropic_solid_vector_expansion_elements.h"

namespace oomph
{

//======================================================================
/// \short Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
template <unsigned DIM, unsigned NVECT>
void  AnisotropicSolidVectorExpansionEquations<DIM, NVECT>::
fill_in_generic_residual_contribution_anisotropic_solid_vector_expansion(Vector<double> &residuals, 
                                               DenseMatrix<double> &jacobian, 
                                               DenseMatrix<double> 
                                               &mass_matrix,
                                               unsigned flag) 
{
  // Loop over the test functions
  const unsigned n_node = nnode();
  int local_eqn, local_unknown = 0;
  Vector<double> dvectdt(NVECT, 0.0);
  for(unsigned l=0;l<n_node;l++){

    dvectors_dt(l, dvectdt);

    //loop over the vectors and their components
    for(unsigned index = vect_min_index_anisotropic_solid_vector_expansion(); index < vect_max_index_anisotropic_solid_vector_expansion(); index++){
      local_eqn = nodal_local_eqn(l, index);
      //IF it's not a boundary condition
      if(local_eqn >= 0){

        //!!!!! INCLUDE ALE????

        residuals[local_eqn] = dvectdt[index - vect_min_index_anisotropic_solid_vector_expansion()];

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
        } //End jacobian

      } //End if not boundary
    } //End loop over vector variables
  } //End loop over nodes
} //End generic residual etc


//======================================================================
/// \short Output function:
///
///   x,y,u,w_x,w_y   or    x,y,z,u,w_x,w_y,w_z
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NVECT>
void  AnisotropicSolidVectorExpansionEquations<DIM, NVECT>::output(std::ostream &outfile, 
                                               const unsigned &nplot)
{ }

//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM, unsigned NVECT>
void AnisotropicSolidVectorExpansionEquations<DIM, NVECT>::output(FILE* file_pt,
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
template <unsigned DIM, unsigned NVECT>
void AnisotropicSolidVectorExpansionEquations<DIM, NVECT>::output_fct(std::ostream &outfile, 
                                                                      const unsigned &nplot, 
                                                                      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{  }






















//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<unsigned DIM, unsigned NVECT, unsigned NNODE_1D>
const unsigned QAnisotropicSolidVectorExpansionElement<DIM, NVECT,NNODE_1D>::Initial_Nvalue = DIM*NVECT;
                                                                                              //dim*number of vectors


//====================================================================
// Force build of templates with number of vectors up to the dimension
//====================================================================


//Force build eqations
template class AnisotropicSolidVectorExpansionEquations<1,1>;

template class AnisotropicSolidVectorExpansionEquations<2,1>;
template class AnisotropicSolidVectorExpansionEquations<2,2>;

template class AnisotropicSolidVectorExpansionEquations<3,1>;
template class AnisotropicSolidVectorExpansionEquations<3,2>;
template class AnisotropicSolidVectorExpansionEquations<3,3>;


//Force build Q elements

//Dim 1
template class QAnisotropicSolidVectorExpansionElement<1,1,2>;
template class QAnisotropicSolidVectorExpansionElement<1,1,3>;
template class QAnisotropicSolidVectorExpansionElement<1,1,4>;

//Dim 2
  //2 nnode_1d
template class QAnisotropicSolidVectorExpansionElement<2,1,2>;
template class QAnisotropicSolidVectorExpansionElement<2,2,2>;
  //3 nnode_1d
template class QAnisotropicSolidVectorExpansionElement<2,1,3>;
template class QAnisotropicSolidVectorExpansionElement<2,2,3>;
  //4 nnode_1d
template class QAnisotropicSolidVectorExpansionElement<2,1,4>;
template class QAnisotropicSolidVectorExpansionElement<2,2,4>;

//Dim 3
  //2 nnode_1d
template class QAnisotropicSolidVectorExpansionElement<3,1,2>;
template class QAnisotropicSolidVectorExpansionElement<3,2,2>;
template class QAnisotropicSolidVectorExpansionElement<3,3,2>;
  //3 nnode_1d
template class QAnisotropicSolidVectorExpansionElement<3,1,3>;
template class QAnisotropicSolidVectorExpansionElement<3,2,3>;
template class QAnisotropicSolidVectorExpansionElement<3,3,3>;
  //4 nnode_1d
template class QAnisotropicSolidVectorExpansionElement<3,1,4>;
template class QAnisotropicSolidVectorExpansionElement<3,2,4>;
template class QAnisotropicSolidVectorExpansionElement<3,3,4>;


/////////////////////////////////////////////////////////////////////////
// TAnisotropicSolidVectorExpansionElement
/////////////////////////////////////////////////////////////////////////



//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<unsigned DIM, unsigned NVECT, unsigned NNODE_1D>
const unsigned TAnisotropicSolidVectorExpansionElement<DIM, NVECT,NNODE_1D>::Initial_Nvalue = DIM*NVECT;
                                                                  //dim*number of vectors

//====================================================================
// Force build of templates with number of vectors up to the dimension
//====================================================================
//dim 1
template class TAnisotropicSolidVectorExpansionElement<1,1,2>;
template class TAnisotropicSolidVectorExpansionElement<1,1,3>;
template class TAnisotropicSolidVectorExpansionElement<1,1,4>;

//dim 2
  //2 nnode_1d
template class TAnisotropicSolidVectorExpansionElement<2,1,2>;
template class TAnisotropicSolidVectorExpansionElement<2,2,2>;
  //3 nnode_1d
template class TAnisotropicSolidVectorExpansionElement<2,1,3>;
template class TAnisotropicSolidVectorExpansionElement<2,2,3>;
  //4 nnode_1d
template class TAnisotropicSolidVectorExpansionElement<2,1,4>;
template class TAnisotropicSolidVectorExpansionElement<2,2,4>;

//dim 3
  //2 nnode_1d
template class TAnisotropicSolidVectorExpansionElement<3,1,2>;
template class TAnisotropicSolidVectorExpansionElement<3,2,2>;
template class TAnisotropicSolidVectorExpansionElement<3,3,2>;
  //3 nnode_1d
template class TAnisotropicSolidVectorExpansionElement<3,1,3>;
template class TAnisotropicSolidVectorExpansionElement<3,2,3>;
template class TAnisotropicSolidVectorExpansionElement<3,3,3>;
}
