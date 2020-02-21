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
#include "storage_enrichment_elements.h"

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
template <unsigned NUM>
void  StorageEnrichmentEquations<NUM>::
fill_in_generic_residual_contribution_storage_enrichment(Vector<double> &residuals, 
                                               DenseMatrix<double> &jacobian, 
                                               DenseMatrix<double> 
                                               &mass_matrix,
                                               unsigned flag) 
{
  // Loop over the test functions
  const unsigned n_node = nnode();
  int local_eqn = 0;
  for(unsigned l=0;l<n_node;l++){

    // Get the data's timestepper
    TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();
    const unsigned n_time = time_stepper_pt->ntstorage();

    //loop over the vectors and their components
    for(unsigned index = min_index_storage_enrichment(); index < max_index_storage_enrichment(); index++){
      //get the local equation number
      local_eqn = nodal_local_eqn(l, index);

      //IF it's not a boundary condition
      if(local_eqn >= 0){

        for(unsigned t=0;t<n_time;t++){

          //Since the data is purely passive we DON'T want it to change with the movement of the mesh
          residuals[local_eqn] -= time_stepper_pt->weight(1,t)*nodal_value(t,l,index);

        }

        //Jacobian, just fill in the diagonal entry for this variable with the suitable weight
        if(flag){
          //Add to the jacobian
          jacobian(local_eqn,local_eqn) -= time_stepper_pt->weight(1,0);

        } //End jacobian

      } //End if not boundary

    } //End loop over variables

  } //End loop over nodes

} //End generic residual etc


//======================================================================
/// \short Output function:
///
///   x,y,u,w_x,w_y   or    x,y,z,u,w_x,w_y,w_z
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned NUM>
void  StorageEnrichmentEquations<NUM>::output(std::ostream &outfile, 
                                             const unsigned &nplot)
{ }

//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned NUM>
void StorageEnrichmentEquations<NUM>::output(FILE* file_pt,
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
template <unsigned NUM>
void StorageEnrichmentEquations<NUM>::output_fct(std::ostream &outfile, 
                                                const unsigned &nplot, 
                                                FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{  }








































template<unsigned NUM>
const unsigned StorageEnrichmentEquations<NUM>::Initial_Nvalue = NUM;


//Do the force build for
//  DIM = 1,2,3
//  1 <= NUM <= 100, if you need to store more than this number you'll need to add it yourself
//  2 <= NNODE_1D <= 3
template class StorageEnrichmentEquations<1>;
template class StorageEnrichmentEquations<2>;
template class StorageEnrichmentEquations<3>;
template class StorageEnrichmentEquations<4>;
template class StorageEnrichmentEquations<5>;
template class StorageEnrichmentEquations<6>;
template class StorageEnrichmentEquations<7>;
template class StorageEnrichmentEquations<8>;
template class StorageEnrichmentEquations<9>;
template class StorageEnrichmentEquations<10>;
template class StorageEnrichmentEquations<11>;
template class StorageEnrichmentEquations<12>;
template class StorageEnrichmentEquations<13>;
template class StorageEnrichmentEquations<14>;
template class StorageEnrichmentEquations<15>;
template class StorageEnrichmentEquations<16>;
template class StorageEnrichmentEquations<17>;
template class StorageEnrichmentEquations<18>;
template class StorageEnrichmentEquations<19>;
template class StorageEnrichmentEquations<20>;
template class StorageEnrichmentEquations<21>;
template class StorageEnrichmentEquations<22>;
template class StorageEnrichmentEquations<23>;
template class StorageEnrichmentEquations<24>;
template class StorageEnrichmentEquations<25>;
template class StorageEnrichmentEquations<26>;
template class StorageEnrichmentEquations<27>;
template class StorageEnrichmentEquations<28>;
template class StorageEnrichmentEquations<29>;
template class StorageEnrichmentEquations<30>;
template class StorageEnrichmentEquations<31>;
template class StorageEnrichmentEquations<32>;
template class StorageEnrichmentEquations<33>;
template class StorageEnrichmentEquations<34>;
template class StorageEnrichmentEquations<35>;
template class StorageEnrichmentEquations<36>;
template class StorageEnrichmentEquations<37>;
template class StorageEnrichmentEquations<38>;
template class StorageEnrichmentEquations<39>;
template class StorageEnrichmentEquations<40>;
template class StorageEnrichmentEquations<41>;
template class StorageEnrichmentEquations<42>;
template class StorageEnrichmentEquations<43>;
template class StorageEnrichmentEquations<44>;
template class StorageEnrichmentEquations<45>;
template class StorageEnrichmentEquations<46>;
template class StorageEnrichmentEquations<47>;
template class StorageEnrichmentEquations<48>;
template class StorageEnrichmentEquations<49>;
template class StorageEnrichmentEquations<50>;
template class StorageEnrichmentEquations<51>;
template class StorageEnrichmentEquations<52>;
template class StorageEnrichmentEquations<53>;
template class StorageEnrichmentEquations<54>;
template class StorageEnrichmentEquations<55>;
template class StorageEnrichmentEquations<56>;
template class StorageEnrichmentEquations<57>;
template class StorageEnrichmentEquations<58>;
template class StorageEnrichmentEquations<59>;
template class StorageEnrichmentEquations<60>;
template class StorageEnrichmentEquations<61>;
template class StorageEnrichmentEquations<62>;
template class StorageEnrichmentEquations<63>;
template class StorageEnrichmentEquations<64>;
template class StorageEnrichmentEquations<65>;
template class StorageEnrichmentEquations<66>;
template class StorageEnrichmentEquations<67>;
template class StorageEnrichmentEquations<68>;
template class StorageEnrichmentEquations<69>;
template class StorageEnrichmentEquations<70>;
template class StorageEnrichmentEquations<71>;
template class StorageEnrichmentEquations<72>;
template class StorageEnrichmentEquations<73>;
template class StorageEnrichmentEquations<74>;
template class StorageEnrichmentEquations<75>;
template class StorageEnrichmentEquations<76>;
template class StorageEnrichmentEquations<77>;
template class StorageEnrichmentEquations<78>;
template class StorageEnrichmentEquations<79>;
template class StorageEnrichmentEquations<80>;
template class StorageEnrichmentEquations<81>;
template class StorageEnrichmentEquations<82>;
template class StorageEnrichmentEquations<83>;
template class StorageEnrichmentEquations<84>;
template class StorageEnrichmentEquations<85>;
template class StorageEnrichmentEquations<86>;
template class StorageEnrichmentEquations<87>;
template class StorageEnrichmentEquations<88>;
template class StorageEnrichmentEquations<89>;
template class StorageEnrichmentEquations<90>;
template class StorageEnrichmentEquations<91>;
template class StorageEnrichmentEquations<92>;
template class StorageEnrichmentEquations<93>;
template class StorageEnrichmentEquations<94>;
template class StorageEnrichmentEquations<95>;
template class StorageEnrichmentEquations<96>;
template class StorageEnrichmentEquations<97>;
template class StorageEnrichmentEquations<98>;
template class StorageEnrichmentEquations<99>;
template class StorageEnrichmentEquations<100>;

template class QStorageEnrichmentElement<1,1,2>;
template class QStorageEnrichmentElement<1,1,3>;
template class QStorageEnrichmentElement<1,2,2>;
template class QStorageEnrichmentElement<1,2,3>;
template class QStorageEnrichmentElement<1,3,2>;
template class QStorageEnrichmentElement<1,3,3>;
template class QStorageEnrichmentElement<1,4,2>;
template class QStorageEnrichmentElement<1,4,3>;
template class QStorageEnrichmentElement<1,5,2>;
template class QStorageEnrichmentElement<1,5,3>;
template class QStorageEnrichmentElement<1,6,2>;
template class QStorageEnrichmentElement<1,6,3>;
template class QStorageEnrichmentElement<1,7,2>;
template class QStorageEnrichmentElement<1,7,3>;
template class QStorageEnrichmentElement<1,8,2>;
template class QStorageEnrichmentElement<1,8,3>;
template class QStorageEnrichmentElement<1,9,2>;
template class QStorageEnrichmentElement<1,9,3>;
template class QStorageEnrichmentElement<1,10,2>;
template class QStorageEnrichmentElement<1,10,3>;
template class QStorageEnrichmentElement<1,11,2>;
template class QStorageEnrichmentElement<1,11,3>;
template class QStorageEnrichmentElement<1,12,2>;
template class QStorageEnrichmentElement<1,12,3>;
template class QStorageEnrichmentElement<1,13,2>;
template class QStorageEnrichmentElement<1,13,3>;
template class QStorageEnrichmentElement<1,14,2>;
template class QStorageEnrichmentElement<1,14,3>;
template class QStorageEnrichmentElement<1,15,2>;
template class QStorageEnrichmentElement<1,15,3>;
template class QStorageEnrichmentElement<1,16,2>;
template class QStorageEnrichmentElement<1,16,3>;
template class QStorageEnrichmentElement<1,17,2>;
template class QStorageEnrichmentElement<1,17,3>;
template class QStorageEnrichmentElement<1,18,2>;
template class QStorageEnrichmentElement<1,18,3>;
template class QStorageEnrichmentElement<1,19,2>;
template class QStorageEnrichmentElement<1,19,3>;
template class QStorageEnrichmentElement<1,20,2>;
template class QStorageEnrichmentElement<1,20,3>;
template class QStorageEnrichmentElement<1,21,2>;
template class QStorageEnrichmentElement<1,21,3>;
template class QStorageEnrichmentElement<1,22,2>;
template class QStorageEnrichmentElement<1,22,3>;
template class QStorageEnrichmentElement<1,23,2>;
template class QStorageEnrichmentElement<1,23,3>;
template class QStorageEnrichmentElement<1,24,2>;
template class QStorageEnrichmentElement<1,24,3>;
template class QStorageEnrichmentElement<1,25,2>;
template class QStorageEnrichmentElement<1,25,3>;
template class QStorageEnrichmentElement<1,26,2>;
template class QStorageEnrichmentElement<1,26,3>;
template class QStorageEnrichmentElement<1,27,2>;
template class QStorageEnrichmentElement<1,27,3>;
template class QStorageEnrichmentElement<1,28,2>;
template class QStorageEnrichmentElement<1,28,3>;
template class QStorageEnrichmentElement<1,29,2>;
template class QStorageEnrichmentElement<1,29,3>;
template class QStorageEnrichmentElement<1,30,2>;
template class QStorageEnrichmentElement<1,30,3>;
template class QStorageEnrichmentElement<1,31,2>;
template class QStorageEnrichmentElement<1,31,3>;
template class QStorageEnrichmentElement<1,32,2>;
template class QStorageEnrichmentElement<1,32,3>;
template class QStorageEnrichmentElement<1,33,2>;
template class QStorageEnrichmentElement<1,33,3>;
template class QStorageEnrichmentElement<1,34,2>;
template class QStorageEnrichmentElement<1,34,3>;
template class QStorageEnrichmentElement<1,35,2>;
template class QStorageEnrichmentElement<1,35,3>;
template class QStorageEnrichmentElement<1,36,2>;
template class QStorageEnrichmentElement<1,36,3>;
template class QStorageEnrichmentElement<1,37,2>;
template class QStorageEnrichmentElement<1,37,3>;
template class QStorageEnrichmentElement<1,38,2>;
template class QStorageEnrichmentElement<1,38,3>;
template class QStorageEnrichmentElement<1,39,2>;
template class QStorageEnrichmentElement<1,39,3>;
template class QStorageEnrichmentElement<1,40,2>;
template class QStorageEnrichmentElement<1,40,3>;
template class QStorageEnrichmentElement<1,41,2>;
template class QStorageEnrichmentElement<1,41,3>;
template class QStorageEnrichmentElement<1,42,2>;
template class QStorageEnrichmentElement<1,42,3>;
template class QStorageEnrichmentElement<1,43,2>;
template class QStorageEnrichmentElement<1,43,3>;
template class QStorageEnrichmentElement<1,44,2>;
template class QStorageEnrichmentElement<1,44,3>;
template class QStorageEnrichmentElement<1,45,2>;
template class QStorageEnrichmentElement<1,45,3>;
template class QStorageEnrichmentElement<1,46,2>;
template class QStorageEnrichmentElement<1,46,3>;
template class QStorageEnrichmentElement<1,47,2>;
template class QStorageEnrichmentElement<1,47,3>;
template class QStorageEnrichmentElement<1,48,2>;
template class QStorageEnrichmentElement<1,48,3>;
template class QStorageEnrichmentElement<1,49,2>;
template class QStorageEnrichmentElement<1,49,3>;
template class QStorageEnrichmentElement<1,50,2>;
template class QStorageEnrichmentElement<1,50,3>;
template class QStorageEnrichmentElement<1,51,2>;
template class QStorageEnrichmentElement<1,51,3>;
template class QStorageEnrichmentElement<1,52,2>;
template class QStorageEnrichmentElement<1,52,3>;
template class QStorageEnrichmentElement<1,53,2>;
template class QStorageEnrichmentElement<1,53,3>;
template class QStorageEnrichmentElement<1,54,2>;
template class QStorageEnrichmentElement<1,54,3>;
template class QStorageEnrichmentElement<1,55,2>;
template class QStorageEnrichmentElement<1,55,3>;
template class QStorageEnrichmentElement<1,56,2>;
template class QStorageEnrichmentElement<1,56,3>;
template class QStorageEnrichmentElement<1,57,2>;
template class QStorageEnrichmentElement<1,57,3>;
template class QStorageEnrichmentElement<1,58,2>;
template class QStorageEnrichmentElement<1,58,3>;
template class QStorageEnrichmentElement<1,59,2>;
template class QStorageEnrichmentElement<1,59,3>;
template class QStorageEnrichmentElement<1,60,2>;
template class QStorageEnrichmentElement<1,60,3>;
template class QStorageEnrichmentElement<1,61,2>;
template class QStorageEnrichmentElement<1,61,3>;
template class QStorageEnrichmentElement<1,62,2>;
template class QStorageEnrichmentElement<1,62,3>;
template class QStorageEnrichmentElement<1,63,2>;
template class QStorageEnrichmentElement<1,63,3>;
template class QStorageEnrichmentElement<1,64,2>;
template class QStorageEnrichmentElement<1,64,3>;
template class QStorageEnrichmentElement<1,65,2>;
template class QStorageEnrichmentElement<1,65,3>;
template class QStorageEnrichmentElement<1,66,2>;
template class QStorageEnrichmentElement<1,66,3>;
template class QStorageEnrichmentElement<1,67,2>;
template class QStorageEnrichmentElement<1,67,3>;
template class QStorageEnrichmentElement<1,68,2>;
template class QStorageEnrichmentElement<1,68,3>;
template class QStorageEnrichmentElement<1,69,2>;
template class QStorageEnrichmentElement<1,69,3>;
template class QStorageEnrichmentElement<1,70,2>;
template class QStorageEnrichmentElement<1,70,3>;
template class QStorageEnrichmentElement<1,71,2>;
template class QStorageEnrichmentElement<1,71,3>;
template class QStorageEnrichmentElement<1,72,2>;
template class QStorageEnrichmentElement<1,72,3>;
template class QStorageEnrichmentElement<1,73,2>;
template class QStorageEnrichmentElement<1,73,3>;
template class QStorageEnrichmentElement<1,74,2>;
template class QStorageEnrichmentElement<1,74,3>;
template class QStorageEnrichmentElement<1,75,2>;
template class QStorageEnrichmentElement<1,75,3>;
template class QStorageEnrichmentElement<1,76,2>;
template class QStorageEnrichmentElement<1,76,3>;
template class QStorageEnrichmentElement<1,77,2>;
template class QStorageEnrichmentElement<1,77,3>;
template class QStorageEnrichmentElement<1,78,2>;
template class QStorageEnrichmentElement<1,78,3>;
template class QStorageEnrichmentElement<1,79,2>;
template class QStorageEnrichmentElement<1,79,3>;
template class QStorageEnrichmentElement<1,80,2>;
template class QStorageEnrichmentElement<1,80,3>;
template class QStorageEnrichmentElement<1,81,2>;
template class QStorageEnrichmentElement<1,81,3>;
template class QStorageEnrichmentElement<1,82,2>;
template class QStorageEnrichmentElement<1,82,3>;
template class QStorageEnrichmentElement<1,83,2>;
template class QStorageEnrichmentElement<1,83,3>;
template class QStorageEnrichmentElement<1,84,2>;
template class QStorageEnrichmentElement<1,84,3>;
template class QStorageEnrichmentElement<1,85,2>;
template class QStorageEnrichmentElement<1,85,3>;
template class QStorageEnrichmentElement<1,86,2>;
template class QStorageEnrichmentElement<1,86,3>;
template class QStorageEnrichmentElement<1,87,2>;
template class QStorageEnrichmentElement<1,87,3>;
template class QStorageEnrichmentElement<1,88,2>;
template class QStorageEnrichmentElement<1,88,3>;
template class QStorageEnrichmentElement<1,89,2>;
template class QStorageEnrichmentElement<1,89,3>;
template class QStorageEnrichmentElement<1,90,2>;
template class QStorageEnrichmentElement<1,90,3>;
template class QStorageEnrichmentElement<1,91,2>;
template class QStorageEnrichmentElement<1,91,3>;
template class QStorageEnrichmentElement<1,92,2>;
template class QStorageEnrichmentElement<1,92,3>;
template class QStorageEnrichmentElement<1,93,2>;
template class QStorageEnrichmentElement<1,93,3>;
template class QStorageEnrichmentElement<1,94,2>;
template class QStorageEnrichmentElement<1,94,3>;
template class QStorageEnrichmentElement<1,95,2>;
template class QStorageEnrichmentElement<1,95,3>;
template class QStorageEnrichmentElement<1,96,2>;
template class QStorageEnrichmentElement<1,96,3>;
template class QStorageEnrichmentElement<1,97,2>;
template class QStorageEnrichmentElement<1,97,3>;
template class QStorageEnrichmentElement<1,98,2>;
template class QStorageEnrichmentElement<1,98,3>;
template class QStorageEnrichmentElement<1,99,2>;
template class QStorageEnrichmentElement<1,99,3>;
template class QStorageEnrichmentElement<1,100,2>;
template class QStorageEnrichmentElement<1,100,3>;
template class QStorageEnrichmentElement<2,1,2>;
template class QStorageEnrichmentElement<2,1,3>;
template class QStorageEnrichmentElement<2,2,2>;
template class QStorageEnrichmentElement<2,2,3>;
template class QStorageEnrichmentElement<2,3,2>;
template class QStorageEnrichmentElement<2,3,3>;
template class QStorageEnrichmentElement<2,4,2>;
template class QStorageEnrichmentElement<2,4,3>;
template class QStorageEnrichmentElement<2,5,2>;
template class QStorageEnrichmentElement<2,5,3>;
template class QStorageEnrichmentElement<2,6,2>;
template class QStorageEnrichmentElement<2,6,3>;
template class QStorageEnrichmentElement<2,7,2>;
template class QStorageEnrichmentElement<2,7,3>;
template class QStorageEnrichmentElement<2,8,2>;
template class QStorageEnrichmentElement<2,8,3>;
template class QStorageEnrichmentElement<2,9,2>;
template class QStorageEnrichmentElement<2,9,3>;
template class QStorageEnrichmentElement<2,10,2>;
template class QStorageEnrichmentElement<2,10,3>;
template class QStorageEnrichmentElement<2,11,2>;
template class QStorageEnrichmentElement<2,11,3>;
template class QStorageEnrichmentElement<2,12,2>;
template class QStorageEnrichmentElement<2,12,3>;
template class QStorageEnrichmentElement<2,13,2>;
template class QStorageEnrichmentElement<2,13,3>;
template class QStorageEnrichmentElement<2,14,2>;
template class QStorageEnrichmentElement<2,14,3>;
template class QStorageEnrichmentElement<2,15,2>;
template class QStorageEnrichmentElement<2,15,3>;
template class QStorageEnrichmentElement<2,16,2>;
template class QStorageEnrichmentElement<2,16,3>;
template class QStorageEnrichmentElement<2,17,2>;
template class QStorageEnrichmentElement<2,17,3>;
template class QStorageEnrichmentElement<2,18,2>;
template class QStorageEnrichmentElement<2,18,3>;
template class QStorageEnrichmentElement<2,19,2>;
template class QStorageEnrichmentElement<2,19,3>;
template class QStorageEnrichmentElement<2,20,2>;
template class QStorageEnrichmentElement<2,20,3>;
template class QStorageEnrichmentElement<2,21,2>;
template class QStorageEnrichmentElement<2,21,3>;
template class QStorageEnrichmentElement<2,22,2>;
template class QStorageEnrichmentElement<2,22,3>;
template class QStorageEnrichmentElement<2,23,2>;
template class QStorageEnrichmentElement<2,23,3>;
template class QStorageEnrichmentElement<2,24,2>;
template class QStorageEnrichmentElement<2,24,3>;
template class QStorageEnrichmentElement<2,25,2>;
template class QStorageEnrichmentElement<2,25,3>;
template class QStorageEnrichmentElement<2,26,2>;
template class QStorageEnrichmentElement<2,26,3>;
template class QStorageEnrichmentElement<2,27,2>;
template class QStorageEnrichmentElement<2,27,3>;
template class QStorageEnrichmentElement<2,28,2>;
template class QStorageEnrichmentElement<2,28,3>;
template class QStorageEnrichmentElement<2,29,2>;
template class QStorageEnrichmentElement<2,29,3>;
template class QStorageEnrichmentElement<2,30,2>;
template class QStorageEnrichmentElement<2,30,3>;
template class QStorageEnrichmentElement<2,31,2>;
template class QStorageEnrichmentElement<2,31,3>;
template class QStorageEnrichmentElement<2,32,2>;
template class QStorageEnrichmentElement<2,32,3>;
template class QStorageEnrichmentElement<2,33,2>;
template class QStorageEnrichmentElement<2,33,3>;
template class QStorageEnrichmentElement<2,34,2>;
template class QStorageEnrichmentElement<2,34,3>;
template class QStorageEnrichmentElement<2,35,2>;
template class QStorageEnrichmentElement<2,35,3>;
template class QStorageEnrichmentElement<2,36,2>;
template class QStorageEnrichmentElement<2,36,3>;
template class QStorageEnrichmentElement<2,37,2>;
template class QStorageEnrichmentElement<2,37,3>;
template class QStorageEnrichmentElement<2,38,2>;
template class QStorageEnrichmentElement<2,38,3>;
template class QStorageEnrichmentElement<2,39,2>;
template class QStorageEnrichmentElement<2,39,3>;
template class QStorageEnrichmentElement<2,40,2>;
template class QStorageEnrichmentElement<2,40,3>;
template class QStorageEnrichmentElement<2,41,2>;
template class QStorageEnrichmentElement<2,41,3>;
template class QStorageEnrichmentElement<2,42,2>;
template class QStorageEnrichmentElement<2,42,3>;
template class QStorageEnrichmentElement<2,43,2>;
template class QStorageEnrichmentElement<2,43,3>;
template class QStorageEnrichmentElement<2,44,2>;
template class QStorageEnrichmentElement<2,44,3>;
template class QStorageEnrichmentElement<2,45,2>;
template class QStorageEnrichmentElement<2,45,3>;
template class QStorageEnrichmentElement<2,46,2>;
template class QStorageEnrichmentElement<2,46,3>;
template class QStorageEnrichmentElement<2,47,2>;
template class QStorageEnrichmentElement<2,47,3>;
template class QStorageEnrichmentElement<2,48,2>;
template class QStorageEnrichmentElement<2,48,3>;
template class QStorageEnrichmentElement<2,49,2>;
template class QStorageEnrichmentElement<2,49,3>;
template class QStorageEnrichmentElement<2,50,2>;
template class QStorageEnrichmentElement<2,50,3>;
template class QStorageEnrichmentElement<2,51,2>;
template class QStorageEnrichmentElement<2,51,3>;
template class QStorageEnrichmentElement<2,52,2>;
template class QStorageEnrichmentElement<2,52,3>;
template class QStorageEnrichmentElement<2,53,2>;
template class QStorageEnrichmentElement<2,53,3>;
template class QStorageEnrichmentElement<2,54,2>;
template class QStorageEnrichmentElement<2,54,3>;
template class QStorageEnrichmentElement<2,55,2>;
template class QStorageEnrichmentElement<2,55,3>;
template class QStorageEnrichmentElement<2,56,2>;
template class QStorageEnrichmentElement<2,56,3>;
template class QStorageEnrichmentElement<2,57,2>;
template class QStorageEnrichmentElement<2,57,3>;
template class QStorageEnrichmentElement<2,58,2>;
template class QStorageEnrichmentElement<2,58,3>;
template class QStorageEnrichmentElement<2,59,2>;
template class QStorageEnrichmentElement<2,59,3>;
template class QStorageEnrichmentElement<2,60,2>;
template class QStorageEnrichmentElement<2,60,3>;
template class QStorageEnrichmentElement<2,61,2>;
template class QStorageEnrichmentElement<2,61,3>;
template class QStorageEnrichmentElement<2,62,2>;
template class QStorageEnrichmentElement<2,62,3>;
template class QStorageEnrichmentElement<2,63,2>;
template class QStorageEnrichmentElement<2,63,3>;
template class QStorageEnrichmentElement<2,64,2>;
template class QStorageEnrichmentElement<2,64,3>;
template class QStorageEnrichmentElement<2,65,2>;
template class QStorageEnrichmentElement<2,65,3>;
template class QStorageEnrichmentElement<2,66,2>;
template class QStorageEnrichmentElement<2,66,3>;
template class QStorageEnrichmentElement<2,67,2>;
template class QStorageEnrichmentElement<2,67,3>;
template class QStorageEnrichmentElement<2,68,2>;
template class QStorageEnrichmentElement<2,68,3>;
template class QStorageEnrichmentElement<2,69,2>;
template class QStorageEnrichmentElement<2,69,3>;
template class QStorageEnrichmentElement<2,70,2>;
template class QStorageEnrichmentElement<2,70,3>;
template class QStorageEnrichmentElement<2,71,2>;
template class QStorageEnrichmentElement<2,71,3>;
template class QStorageEnrichmentElement<2,72,2>;
template class QStorageEnrichmentElement<2,72,3>;
template class QStorageEnrichmentElement<2,73,2>;
template class QStorageEnrichmentElement<2,73,3>;
template class QStorageEnrichmentElement<2,74,2>;
template class QStorageEnrichmentElement<2,74,3>;
template class QStorageEnrichmentElement<2,75,2>;
template class QStorageEnrichmentElement<2,75,3>;
template class QStorageEnrichmentElement<2,76,2>;
template class QStorageEnrichmentElement<2,76,3>;
template class QStorageEnrichmentElement<2,77,2>;
template class QStorageEnrichmentElement<2,77,3>;
template class QStorageEnrichmentElement<2,78,2>;
template class QStorageEnrichmentElement<2,78,3>;
template class QStorageEnrichmentElement<2,79,2>;
template class QStorageEnrichmentElement<2,79,3>;
template class QStorageEnrichmentElement<2,80,2>;
template class QStorageEnrichmentElement<2,80,3>;
template class QStorageEnrichmentElement<2,81,2>;
template class QStorageEnrichmentElement<2,81,3>;
template class QStorageEnrichmentElement<2,82,2>;
template class QStorageEnrichmentElement<2,82,3>;
template class QStorageEnrichmentElement<2,83,2>;
template class QStorageEnrichmentElement<2,83,3>;
template class QStorageEnrichmentElement<2,84,2>;
template class QStorageEnrichmentElement<2,84,3>;
template class QStorageEnrichmentElement<2,85,2>;
template class QStorageEnrichmentElement<2,85,3>;
template class QStorageEnrichmentElement<2,86,2>;
template class QStorageEnrichmentElement<2,86,3>;
template class QStorageEnrichmentElement<2,87,2>;
template class QStorageEnrichmentElement<2,87,3>;
template class QStorageEnrichmentElement<2,88,2>;
template class QStorageEnrichmentElement<2,88,3>;
template class QStorageEnrichmentElement<2,89,2>;
template class QStorageEnrichmentElement<2,89,3>;
template class QStorageEnrichmentElement<2,90,2>;
template class QStorageEnrichmentElement<2,90,3>;
template class QStorageEnrichmentElement<2,91,2>;
template class QStorageEnrichmentElement<2,91,3>;
template class QStorageEnrichmentElement<2,92,2>;
template class QStorageEnrichmentElement<2,92,3>;
template class QStorageEnrichmentElement<2,93,2>;
template class QStorageEnrichmentElement<2,93,3>;
template class QStorageEnrichmentElement<2,94,2>;
template class QStorageEnrichmentElement<2,94,3>;
template class QStorageEnrichmentElement<2,95,2>;
template class QStorageEnrichmentElement<2,95,3>;
template class QStorageEnrichmentElement<2,96,2>;
template class QStorageEnrichmentElement<2,96,3>;
template class QStorageEnrichmentElement<2,97,2>;
template class QStorageEnrichmentElement<2,97,3>;
template class QStorageEnrichmentElement<2,98,2>;
template class QStorageEnrichmentElement<2,98,3>;
template class QStorageEnrichmentElement<2,99,2>;
template class QStorageEnrichmentElement<2,99,3>;
template class QStorageEnrichmentElement<2,100,2>;
template class QStorageEnrichmentElement<2,100,3>;
template class QStorageEnrichmentElement<3,1,2>;
template class QStorageEnrichmentElement<3,1,3>;
template class QStorageEnrichmentElement<3,2,2>;
template class QStorageEnrichmentElement<3,2,3>;
template class QStorageEnrichmentElement<3,3,2>;
template class QStorageEnrichmentElement<3,3,3>;
template class QStorageEnrichmentElement<3,4,2>;
template class QStorageEnrichmentElement<3,4,3>;
template class QStorageEnrichmentElement<3,5,2>;
template class QStorageEnrichmentElement<3,5,3>;
template class QStorageEnrichmentElement<3,6,2>;
template class QStorageEnrichmentElement<3,6,3>;
template class QStorageEnrichmentElement<3,7,2>;
template class QStorageEnrichmentElement<3,7,3>;
template class QStorageEnrichmentElement<3,8,2>;
template class QStorageEnrichmentElement<3,8,3>;
template class QStorageEnrichmentElement<3,9,2>;
template class QStorageEnrichmentElement<3,9,3>;
template class QStorageEnrichmentElement<3,10,2>;
template class QStorageEnrichmentElement<3,10,3>;
template class QStorageEnrichmentElement<3,11,2>;
template class QStorageEnrichmentElement<3,11,3>;
template class QStorageEnrichmentElement<3,12,2>;
template class QStorageEnrichmentElement<3,12,3>;
template class QStorageEnrichmentElement<3,13,2>;
template class QStorageEnrichmentElement<3,13,3>;
template class QStorageEnrichmentElement<3,14,2>;
template class QStorageEnrichmentElement<3,14,3>;
template class QStorageEnrichmentElement<3,15,2>;
template class QStorageEnrichmentElement<3,15,3>;
template class QStorageEnrichmentElement<3,16,2>;
template class QStorageEnrichmentElement<3,16,3>;
template class QStorageEnrichmentElement<3,17,2>;
template class QStorageEnrichmentElement<3,17,3>;
template class QStorageEnrichmentElement<3,18,2>;
template class QStorageEnrichmentElement<3,18,3>;
template class QStorageEnrichmentElement<3,19,2>;
template class QStorageEnrichmentElement<3,19,3>;
template class QStorageEnrichmentElement<3,20,2>;
template class QStorageEnrichmentElement<3,20,3>;
template class QStorageEnrichmentElement<3,21,2>;
template class QStorageEnrichmentElement<3,21,3>;
template class QStorageEnrichmentElement<3,22,2>;
template class QStorageEnrichmentElement<3,22,3>;
template class QStorageEnrichmentElement<3,23,2>;
template class QStorageEnrichmentElement<3,23,3>;
template class QStorageEnrichmentElement<3,24,2>;
template class QStorageEnrichmentElement<3,24,3>;
template class QStorageEnrichmentElement<3,25,2>;
template class QStorageEnrichmentElement<3,25,3>;
template class QStorageEnrichmentElement<3,26,2>;
template class QStorageEnrichmentElement<3,26,3>;
template class QStorageEnrichmentElement<3,27,2>;
template class QStorageEnrichmentElement<3,27,3>;
template class QStorageEnrichmentElement<3,28,2>;
template class QStorageEnrichmentElement<3,28,3>;
template class QStorageEnrichmentElement<3,29,2>;
template class QStorageEnrichmentElement<3,29,3>;
template class QStorageEnrichmentElement<3,30,2>;
template class QStorageEnrichmentElement<3,30,3>;
template class QStorageEnrichmentElement<3,31,2>;
template class QStorageEnrichmentElement<3,31,3>;
template class QStorageEnrichmentElement<3,32,2>;
template class QStorageEnrichmentElement<3,32,3>;
template class QStorageEnrichmentElement<3,33,2>;
template class QStorageEnrichmentElement<3,33,3>;
template class QStorageEnrichmentElement<3,34,2>;
template class QStorageEnrichmentElement<3,34,3>;
template class QStorageEnrichmentElement<3,35,2>;
template class QStorageEnrichmentElement<3,35,3>;
template class QStorageEnrichmentElement<3,36,2>;
template class QStorageEnrichmentElement<3,36,3>;
template class QStorageEnrichmentElement<3,37,2>;
template class QStorageEnrichmentElement<3,37,3>;
template class QStorageEnrichmentElement<3,38,2>;
template class QStorageEnrichmentElement<3,38,3>;
template class QStorageEnrichmentElement<3,39,2>;
template class QStorageEnrichmentElement<3,39,3>;
template class QStorageEnrichmentElement<3,40,2>;
template class QStorageEnrichmentElement<3,40,3>;
template class QStorageEnrichmentElement<3,41,2>;
template class QStorageEnrichmentElement<3,41,3>;
template class QStorageEnrichmentElement<3,42,2>;
template class QStorageEnrichmentElement<3,42,3>;
template class QStorageEnrichmentElement<3,43,2>;
template class QStorageEnrichmentElement<3,43,3>;
template class QStorageEnrichmentElement<3,44,2>;
template class QStorageEnrichmentElement<3,44,3>;
template class QStorageEnrichmentElement<3,45,2>;
template class QStorageEnrichmentElement<3,45,3>;
template class QStorageEnrichmentElement<3,46,2>;
template class QStorageEnrichmentElement<3,46,3>;
template class QStorageEnrichmentElement<3,47,2>;
template class QStorageEnrichmentElement<3,47,3>;
template class QStorageEnrichmentElement<3,48,2>;
template class QStorageEnrichmentElement<3,48,3>;
template class QStorageEnrichmentElement<3,49,2>;
template class QStorageEnrichmentElement<3,49,3>;
template class QStorageEnrichmentElement<3,50,2>;
template class QStorageEnrichmentElement<3,50,3>;
template class QStorageEnrichmentElement<3,51,2>;
template class QStorageEnrichmentElement<3,51,3>;
template class QStorageEnrichmentElement<3,52,2>;
template class QStorageEnrichmentElement<3,52,3>;
template class QStorageEnrichmentElement<3,53,2>;
template class QStorageEnrichmentElement<3,53,3>;
template class QStorageEnrichmentElement<3,54,2>;
template class QStorageEnrichmentElement<3,54,3>;
template class QStorageEnrichmentElement<3,55,2>;
template class QStorageEnrichmentElement<3,55,3>;
template class QStorageEnrichmentElement<3,56,2>;
template class QStorageEnrichmentElement<3,56,3>;
template class QStorageEnrichmentElement<3,57,2>;
template class QStorageEnrichmentElement<3,57,3>;
template class QStorageEnrichmentElement<3,58,2>;
template class QStorageEnrichmentElement<3,58,3>;
template class QStorageEnrichmentElement<3,59,2>;
template class QStorageEnrichmentElement<3,59,3>;
template class QStorageEnrichmentElement<3,60,2>;
template class QStorageEnrichmentElement<3,60,3>;
template class QStorageEnrichmentElement<3,61,2>;
template class QStorageEnrichmentElement<3,61,3>;
template class QStorageEnrichmentElement<3,62,2>;
template class QStorageEnrichmentElement<3,62,3>;
template class QStorageEnrichmentElement<3,63,2>;
template class QStorageEnrichmentElement<3,63,3>;
template class QStorageEnrichmentElement<3,64,2>;
template class QStorageEnrichmentElement<3,64,3>;
template class QStorageEnrichmentElement<3,65,2>;
template class QStorageEnrichmentElement<3,65,3>;
template class QStorageEnrichmentElement<3,66,2>;
template class QStorageEnrichmentElement<3,66,3>;
template class QStorageEnrichmentElement<3,67,2>;
template class QStorageEnrichmentElement<3,67,3>;
template class QStorageEnrichmentElement<3,68,2>;
template class QStorageEnrichmentElement<3,68,3>;
template class QStorageEnrichmentElement<3,69,2>;
template class QStorageEnrichmentElement<3,69,3>;
template class QStorageEnrichmentElement<3,70,2>;
template class QStorageEnrichmentElement<3,70,3>;
template class QStorageEnrichmentElement<3,71,2>;
template class QStorageEnrichmentElement<3,71,3>;
template class QStorageEnrichmentElement<3,72,2>;
template class QStorageEnrichmentElement<3,72,3>;
template class QStorageEnrichmentElement<3,73,2>;
template class QStorageEnrichmentElement<3,73,3>;
template class QStorageEnrichmentElement<3,74,2>;
template class QStorageEnrichmentElement<3,74,3>;
template class QStorageEnrichmentElement<3,75,2>;
template class QStorageEnrichmentElement<3,75,3>;
template class QStorageEnrichmentElement<3,76,2>;
template class QStorageEnrichmentElement<3,76,3>;
template class QStorageEnrichmentElement<3,77,2>;
template class QStorageEnrichmentElement<3,77,3>;
template class QStorageEnrichmentElement<3,78,2>;
template class QStorageEnrichmentElement<3,78,3>;
template class QStorageEnrichmentElement<3,79,2>;
template class QStorageEnrichmentElement<3,79,3>;
template class QStorageEnrichmentElement<3,80,2>;
template class QStorageEnrichmentElement<3,80,3>;
template class QStorageEnrichmentElement<3,81,2>;
template class QStorageEnrichmentElement<3,81,3>;
template class QStorageEnrichmentElement<3,82,2>;
template class QStorageEnrichmentElement<3,82,3>;
template class QStorageEnrichmentElement<3,83,2>;
template class QStorageEnrichmentElement<3,83,3>;
template class QStorageEnrichmentElement<3,84,2>;
template class QStorageEnrichmentElement<3,84,3>;
template class QStorageEnrichmentElement<3,85,2>;
template class QStorageEnrichmentElement<3,85,3>;
template class QStorageEnrichmentElement<3,86,2>;
template class QStorageEnrichmentElement<3,86,3>;
template class QStorageEnrichmentElement<3,87,2>;
template class QStorageEnrichmentElement<3,87,3>;
template class QStorageEnrichmentElement<3,88,2>;
template class QStorageEnrichmentElement<3,88,3>;
template class QStorageEnrichmentElement<3,89,2>;
template class QStorageEnrichmentElement<3,89,3>;
template class QStorageEnrichmentElement<3,90,2>;
template class QStorageEnrichmentElement<3,90,3>;
template class QStorageEnrichmentElement<3,91,2>;
template class QStorageEnrichmentElement<3,91,3>;
template class QStorageEnrichmentElement<3,92,2>;
template class QStorageEnrichmentElement<3,92,3>;
template class QStorageEnrichmentElement<3,93,2>;
template class QStorageEnrichmentElement<3,93,3>;
template class QStorageEnrichmentElement<3,94,2>;
template class QStorageEnrichmentElement<3,94,3>;
template class QStorageEnrichmentElement<3,95,2>;
template class QStorageEnrichmentElement<3,95,3>;
template class QStorageEnrichmentElement<3,96,2>;
template class QStorageEnrichmentElement<3,96,3>;
template class QStorageEnrichmentElement<3,97,2>;
template class QStorageEnrichmentElement<3,97,3>;
template class QStorageEnrichmentElement<3,98,2>;
template class QStorageEnrichmentElement<3,98,3>;
template class QStorageEnrichmentElement<3,99,2>;
template class QStorageEnrichmentElement<3,99,3>;
template class QStorageEnrichmentElement<3,100,2>;
template class QStorageEnrichmentElement<3,100,3>;

template class TStorageEnrichmentElement<1,1,2>;
template class TStorageEnrichmentElement<1,1,3>;
template class TStorageEnrichmentElement<1,2,2>;
template class TStorageEnrichmentElement<1,2,3>;
template class TStorageEnrichmentElement<1,3,2>;
template class TStorageEnrichmentElement<1,3,3>;
template class TStorageEnrichmentElement<1,4,2>;
template class TStorageEnrichmentElement<1,4,3>;
template class TStorageEnrichmentElement<1,5,2>;
template class TStorageEnrichmentElement<1,5,3>;
template class TStorageEnrichmentElement<1,6,2>;
template class TStorageEnrichmentElement<1,6,3>;
template class TStorageEnrichmentElement<1,7,2>;
template class TStorageEnrichmentElement<1,7,3>;
template class TStorageEnrichmentElement<1,8,2>;
template class TStorageEnrichmentElement<1,8,3>;
template class TStorageEnrichmentElement<1,9,2>;
template class TStorageEnrichmentElement<1,9,3>;
template class TStorageEnrichmentElement<1,10,2>;
template class TStorageEnrichmentElement<1,10,3>;
template class TStorageEnrichmentElement<1,11,2>;
template class TStorageEnrichmentElement<1,11,3>;
template class TStorageEnrichmentElement<1,12,2>;
template class TStorageEnrichmentElement<1,12,3>;
template class TStorageEnrichmentElement<1,13,2>;
template class TStorageEnrichmentElement<1,13,3>;
template class TStorageEnrichmentElement<1,14,2>;
template class TStorageEnrichmentElement<1,14,3>;
template class TStorageEnrichmentElement<1,15,2>;
template class TStorageEnrichmentElement<1,15,3>;
template class TStorageEnrichmentElement<1,16,2>;
template class TStorageEnrichmentElement<1,16,3>;
template class TStorageEnrichmentElement<1,17,2>;
template class TStorageEnrichmentElement<1,17,3>;
template class TStorageEnrichmentElement<1,18,2>;
template class TStorageEnrichmentElement<1,18,3>;
template class TStorageEnrichmentElement<1,19,2>;
template class TStorageEnrichmentElement<1,19,3>;
template class TStorageEnrichmentElement<1,20,2>;
template class TStorageEnrichmentElement<1,20,3>;
template class TStorageEnrichmentElement<1,21,2>;
template class TStorageEnrichmentElement<1,21,3>;
template class TStorageEnrichmentElement<1,22,2>;
template class TStorageEnrichmentElement<1,22,3>;
template class TStorageEnrichmentElement<1,23,2>;
template class TStorageEnrichmentElement<1,23,3>;
template class TStorageEnrichmentElement<1,24,2>;
template class TStorageEnrichmentElement<1,24,3>;
template class TStorageEnrichmentElement<1,25,2>;
template class TStorageEnrichmentElement<1,25,3>;
template class TStorageEnrichmentElement<1,26,2>;
template class TStorageEnrichmentElement<1,26,3>;
template class TStorageEnrichmentElement<1,27,2>;
template class TStorageEnrichmentElement<1,27,3>;
template class TStorageEnrichmentElement<1,28,2>;
template class TStorageEnrichmentElement<1,28,3>;
template class TStorageEnrichmentElement<1,29,2>;
template class TStorageEnrichmentElement<1,29,3>;
template class TStorageEnrichmentElement<1,30,2>;
template class TStorageEnrichmentElement<1,30,3>;
template class TStorageEnrichmentElement<1,31,2>;
template class TStorageEnrichmentElement<1,31,3>;
template class TStorageEnrichmentElement<1,32,2>;
template class TStorageEnrichmentElement<1,32,3>;
template class TStorageEnrichmentElement<1,33,2>;
template class TStorageEnrichmentElement<1,33,3>;
template class TStorageEnrichmentElement<1,34,2>;
template class TStorageEnrichmentElement<1,34,3>;
template class TStorageEnrichmentElement<1,35,2>;
template class TStorageEnrichmentElement<1,35,3>;
template class TStorageEnrichmentElement<1,36,2>;
template class TStorageEnrichmentElement<1,36,3>;
template class TStorageEnrichmentElement<1,37,2>;
template class TStorageEnrichmentElement<1,37,3>;
template class TStorageEnrichmentElement<1,38,2>;
template class TStorageEnrichmentElement<1,38,3>;
template class TStorageEnrichmentElement<1,39,2>;
template class TStorageEnrichmentElement<1,39,3>;
template class TStorageEnrichmentElement<1,40,2>;
template class TStorageEnrichmentElement<1,40,3>;
template class TStorageEnrichmentElement<1,41,2>;
template class TStorageEnrichmentElement<1,41,3>;
template class TStorageEnrichmentElement<1,42,2>;
template class TStorageEnrichmentElement<1,42,3>;
template class TStorageEnrichmentElement<1,43,2>;
template class TStorageEnrichmentElement<1,43,3>;
template class TStorageEnrichmentElement<1,44,2>;
template class TStorageEnrichmentElement<1,44,3>;
template class TStorageEnrichmentElement<1,45,2>;
template class TStorageEnrichmentElement<1,45,3>;
template class TStorageEnrichmentElement<1,46,2>;
template class TStorageEnrichmentElement<1,46,3>;
template class TStorageEnrichmentElement<1,47,2>;
template class TStorageEnrichmentElement<1,47,3>;
template class TStorageEnrichmentElement<1,48,2>;
template class TStorageEnrichmentElement<1,48,3>;
template class TStorageEnrichmentElement<1,49,2>;
template class TStorageEnrichmentElement<1,49,3>;
template class TStorageEnrichmentElement<1,50,2>;
template class TStorageEnrichmentElement<1,50,3>;
template class TStorageEnrichmentElement<1,51,2>;
template class TStorageEnrichmentElement<1,51,3>;
template class TStorageEnrichmentElement<1,52,2>;
template class TStorageEnrichmentElement<1,52,3>;
template class TStorageEnrichmentElement<1,53,2>;
template class TStorageEnrichmentElement<1,53,3>;
template class TStorageEnrichmentElement<1,54,2>;
template class TStorageEnrichmentElement<1,54,3>;
template class TStorageEnrichmentElement<1,55,2>;
template class TStorageEnrichmentElement<1,55,3>;
template class TStorageEnrichmentElement<1,56,2>;
template class TStorageEnrichmentElement<1,56,3>;
template class TStorageEnrichmentElement<1,57,2>;
template class TStorageEnrichmentElement<1,57,3>;
template class TStorageEnrichmentElement<1,58,2>;
template class TStorageEnrichmentElement<1,58,3>;
template class TStorageEnrichmentElement<1,59,2>;
template class TStorageEnrichmentElement<1,59,3>;
template class TStorageEnrichmentElement<1,60,2>;
template class TStorageEnrichmentElement<1,60,3>;
template class TStorageEnrichmentElement<1,61,2>;
template class TStorageEnrichmentElement<1,61,3>;
template class TStorageEnrichmentElement<1,62,2>;
template class TStorageEnrichmentElement<1,62,3>;
template class TStorageEnrichmentElement<1,63,2>;
template class TStorageEnrichmentElement<1,63,3>;
template class TStorageEnrichmentElement<1,64,2>;
template class TStorageEnrichmentElement<1,64,3>;
template class TStorageEnrichmentElement<1,65,2>;
template class TStorageEnrichmentElement<1,65,3>;
template class TStorageEnrichmentElement<1,66,2>;
template class TStorageEnrichmentElement<1,66,3>;
template class TStorageEnrichmentElement<1,67,2>;
template class TStorageEnrichmentElement<1,67,3>;
template class TStorageEnrichmentElement<1,68,2>;
template class TStorageEnrichmentElement<1,68,3>;
template class TStorageEnrichmentElement<1,69,2>;
template class TStorageEnrichmentElement<1,69,3>;
template class TStorageEnrichmentElement<1,70,2>;
template class TStorageEnrichmentElement<1,70,3>;
template class TStorageEnrichmentElement<1,71,2>;
template class TStorageEnrichmentElement<1,71,3>;
template class TStorageEnrichmentElement<1,72,2>;
template class TStorageEnrichmentElement<1,72,3>;
template class TStorageEnrichmentElement<1,73,2>;
template class TStorageEnrichmentElement<1,73,3>;
template class TStorageEnrichmentElement<1,74,2>;
template class TStorageEnrichmentElement<1,74,3>;
template class TStorageEnrichmentElement<1,75,2>;
template class TStorageEnrichmentElement<1,75,3>;
template class TStorageEnrichmentElement<1,76,2>;
template class TStorageEnrichmentElement<1,76,3>;
template class TStorageEnrichmentElement<1,77,2>;
template class TStorageEnrichmentElement<1,77,3>;
template class TStorageEnrichmentElement<1,78,2>;
template class TStorageEnrichmentElement<1,78,3>;
template class TStorageEnrichmentElement<1,79,2>;
template class TStorageEnrichmentElement<1,79,3>;
template class TStorageEnrichmentElement<1,80,2>;
template class TStorageEnrichmentElement<1,80,3>;
template class TStorageEnrichmentElement<1,81,2>;
template class TStorageEnrichmentElement<1,81,3>;
template class TStorageEnrichmentElement<1,82,2>;
template class TStorageEnrichmentElement<1,82,3>;
template class TStorageEnrichmentElement<1,83,2>;
template class TStorageEnrichmentElement<1,83,3>;
template class TStorageEnrichmentElement<1,84,2>;
template class TStorageEnrichmentElement<1,84,3>;
template class TStorageEnrichmentElement<1,85,2>;
template class TStorageEnrichmentElement<1,85,3>;
template class TStorageEnrichmentElement<1,86,2>;
template class TStorageEnrichmentElement<1,86,3>;
template class TStorageEnrichmentElement<1,87,2>;
template class TStorageEnrichmentElement<1,87,3>;
template class TStorageEnrichmentElement<1,88,2>;
template class TStorageEnrichmentElement<1,88,3>;
template class TStorageEnrichmentElement<1,89,2>;
template class TStorageEnrichmentElement<1,89,3>;
template class TStorageEnrichmentElement<1,90,2>;
template class TStorageEnrichmentElement<1,90,3>;
template class TStorageEnrichmentElement<1,91,2>;
template class TStorageEnrichmentElement<1,91,3>;
template class TStorageEnrichmentElement<1,92,2>;
template class TStorageEnrichmentElement<1,92,3>;
template class TStorageEnrichmentElement<1,93,2>;
template class TStorageEnrichmentElement<1,93,3>;
template class TStorageEnrichmentElement<1,94,2>;
template class TStorageEnrichmentElement<1,94,3>;
template class TStorageEnrichmentElement<1,95,2>;
template class TStorageEnrichmentElement<1,95,3>;
template class TStorageEnrichmentElement<1,96,2>;
template class TStorageEnrichmentElement<1,96,3>;
template class TStorageEnrichmentElement<1,97,2>;
template class TStorageEnrichmentElement<1,97,3>;
template class TStorageEnrichmentElement<1,98,2>;
template class TStorageEnrichmentElement<1,98,3>;
template class TStorageEnrichmentElement<1,99,2>;
template class TStorageEnrichmentElement<1,99,3>;
template class TStorageEnrichmentElement<1,100,2>;
template class TStorageEnrichmentElement<1,100,3>;
template class TStorageEnrichmentElement<2,1,2>;
template class TStorageEnrichmentElement<2,1,3>;
template class TStorageEnrichmentElement<2,2,2>;
template class TStorageEnrichmentElement<2,2,3>;
template class TStorageEnrichmentElement<2,3,2>;
template class TStorageEnrichmentElement<2,3,3>;
template class TStorageEnrichmentElement<2,4,2>;
template class TStorageEnrichmentElement<2,4,3>;
template class TStorageEnrichmentElement<2,5,2>;
template class TStorageEnrichmentElement<2,5,3>;
template class TStorageEnrichmentElement<2,6,2>;
template class TStorageEnrichmentElement<2,6,3>;
template class TStorageEnrichmentElement<2,7,2>;
template class TStorageEnrichmentElement<2,7,3>;
template class TStorageEnrichmentElement<2,8,2>;
template class TStorageEnrichmentElement<2,8,3>;
template class TStorageEnrichmentElement<2,9,2>;
template class TStorageEnrichmentElement<2,9,3>;
template class TStorageEnrichmentElement<2,10,2>;
template class TStorageEnrichmentElement<2,10,3>;
template class TStorageEnrichmentElement<2,11,2>;
template class TStorageEnrichmentElement<2,11,3>;
template class TStorageEnrichmentElement<2,12,2>;
template class TStorageEnrichmentElement<2,12,3>;
template class TStorageEnrichmentElement<2,13,2>;
template class TStorageEnrichmentElement<2,13,3>;
template class TStorageEnrichmentElement<2,14,2>;
template class TStorageEnrichmentElement<2,14,3>;
template class TStorageEnrichmentElement<2,15,2>;
template class TStorageEnrichmentElement<2,15,3>;
template class TStorageEnrichmentElement<2,16,2>;
template class TStorageEnrichmentElement<2,16,3>;
template class TStorageEnrichmentElement<2,17,2>;
template class TStorageEnrichmentElement<2,17,3>;
template class TStorageEnrichmentElement<2,18,2>;
template class TStorageEnrichmentElement<2,18,3>;
template class TStorageEnrichmentElement<2,19,2>;
template class TStorageEnrichmentElement<2,19,3>;
template class TStorageEnrichmentElement<2,20,2>;
template class TStorageEnrichmentElement<2,20,3>;
template class TStorageEnrichmentElement<2,21,2>;
template class TStorageEnrichmentElement<2,21,3>;
template class TStorageEnrichmentElement<2,22,2>;
template class TStorageEnrichmentElement<2,22,3>;
template class TStorageEnrichmentElement<2,23,2>;
template class TStorageEnrichmentElement<2,23,3>;
template class TStorageEnrichmentElement<2,24,2>;
template class TStorageEnrichmentElement<2,24,3>;
template class TStorageEnrichmentElement<2,25,2>;
template class TStorageEnrichmentElement<2,25,3>;
template class TStorageEnrichmentElement<2,26,2>;
template class TStorageEnrichmentElement<2,26,3>;
template class TStorageEnrichmentElement<2,27,2>;
template class TStorageEnrichmentElement<2,27,3>;
template class TStorageEnrichmentElement<2,28,2>;
template class TStorageEnrichmentElement<2,28,3>;
template class TStorageEnrichmentElement<2,29,2>;
template class TStorageEnrichmentElement<2,29,3>;
template class TStorageEnrichmentElement<2,30,2>;
template class TStorageEnrichmentElement<2,30,3>;
template class TStorageEnrichmentElement<2,31,2>;
template class TStorageEnrichmentElement<2,31,3>;
template class TStorageEnrichmentElement<2,32,2>;
template class TStorageEnrichmentElement<2,32,3>;
template class TStorageEnrichmentElement<2,33,2>;
template class TStorageEnrichmentElement<2,33,3>;
template class TStorageEnrichmentElement<2,34,2>;
template class TStorageEnrichmentElement<2,34,3>;
template class TStorageEnrichmentElement<2,35,2>;
template class TStorageEnrichmentElement<2,35,3>;
template class TStorageEnrichmentElement<2,36,2>;
template class TStorageEnrichmentElement<2,36,3>;
template class TStorageEnrichmentElement<2,37,2>;
template class TStorageEnrichmentElement<2,37,3>;
template class TStorageEnrichmentElement<2,38,2>;
template class TStorageEnrichmentElement<2,38,3>;
template class TStorageEnrichmentElement<2,39,2>;
template class TStorageEnrichmentElement<2,39,3>;
template class TStorageEnrichmentElement<2,40,2>;
template class TStorageEnrichmentElement<2,40,3>;
template class TStorageEnrichmentElement<2,41,2>;
template class TStorageEnrichmentElement<2,41,3>;
template class TStorageEnrichmentElement<2,42,2>;
template class TStorageEnrichmentElement<2,42,3>;
template class TStorageEnrichmentElement<2,43,2>;
template class TStorageEnrichmentElement<2,43,3>;
template class TStorageEnrichmentElement<2,44,2>;
template class TStorageEnrichmentElement<2,44,3>;
template class TStorageEnrichmentElement<2,45,2>;
template class TStorageEnrichmentElement<2,45,3>;
template class TStorageEnrichmentElement<2,46,2>;
template class TStorageEnrichmentElement<2,46,3>;
template class TStorageEnrichmentElement<2,47,2>;
template class TStorageEnrichmentElement<2,47,3>;
template class TStorageEnrichmentElement<2,48,2>;
template class TStorageEnrichmentElement<2,48,3>;
template class TStorageEnrichmentElement<2,49,2>;
template class TStorageEnrichmentElement<2,49,3>;
template class TStorageEnrichmentElement<2,50,2>;
template class TStorageEnrichmentElement<2,50,3>;
template class TStorageEnrichmentElement<2,51,2>;
template class TStorageEnrichmentElement<2,51,3>;
template class TStorageEnrichmentElement<2,52,2>;
template class TStorageEnrichmentElement<2,52,3>;
template class TStorageEnrichmentElement<2,53,2>;
template class TStorageEnrichmentElement<2,53,3>;
template class TStorageEnrichmentElement<2,54,2>;
template class TStorageEnrichmentElement<2,54,3>;
template class TStorageEnrichmentElement<2,55,2>;
template class TStorageEnrichmentElement<2,55,3>;
template class TStorageEnrichmentElement<2,56,2>;
template class TStorageEnrichmentElement<2,56,3>;
template class TStorageEnrichmentElement<2,57,2>;
template class TStorageEnrichmentElement<2,57,3>;
template class TStorageEnrichmentElement<2,58,2>;
template class TStorageEnrichmentElement<2,58,3>;
template class TStorageEnrichmentElement<2,59,2>;
template class TStorageEnrichmentElement<2,59,3>;
template class TStorageEnrichmentElement<2,60,2>;
template class TStorageEnrichmentElement<2,60,3>;
template class TStorageEnrichmentElement<2,61,2>;
template class TStorageEnrichmentElement<2,61,3>;
template class TStorageEnrichmentElement<2,62,2>;
template class TStorageEnrichmentElement<2,62,3>;
template class TStorageEnrichmentElement<2,63,2>;
template class TStorageEnrichmentElement<2,63,3>;
template class TStorageEnrichmentElement<2,64,2>;
template class TStorageEnrichmentElement<2,64,3>;
template class TStorageEnrichmentElement<2,65,2>;
template class TStorageEnrichmentElement<2,65,3>;
template class TStorageEnrichmentElement<2,66,2>;
template class TStorageEnrichmentElement<2,66,3>;
template class TStorageEnrichmentElement<2,67,2>;
template class TStorageEnrichmentElement<2,67,3>;
template class TStorageEnrichmentElement<2,68,2>;
template class TStorageEnrichmentElement<2,68,3>;
template class TStorageEnrichmentElement<2,69,2>;
template class TStorageEnrichmentElement<2,69,3>;
template class TStorageEnrichmentElement<2,70,2>;
template class TStorageEnrichmentElement<2,70,3>;
template class TStorageEnrichmentElement<2,71,2>;
template class TStorageEnrichmentElement<2,71,3>;
template class TStorageEnrichmentElement<2,72,2>;
template class TStorageEnrichmentElement<2,72,3>;
template class TStorageEnrichmentElement<2,73,2>;
template class TStorageEnrichmentElement<2,73,3>;
template class TStorageEnrichmentElement<2,74,2>;
template class TStorageEnrichmentElement<2,74,3>;
template class TStorageEnrichmentElement<2,75,2>;
template class TStorageEnrichmentElement<2,75,3>;
template class TStorageEnrichmentElement<2,76,2>;
template class TStorageEnrichmentElement<2,76,3>;
template class TStorageEnrichmentElement<2,77,2>;
template class TStorageEnrichmentElement<2,77,3>;
template class TStorageEnrichmentElement<2,78,2>;
template class TStorageEnrichmentElement<2,78,3>;
template class TStorageEnrichmentElement<2,79,2>;
template class TStorageEnrichmentElement<2,79,3>;
template class TStorageEnrichmentElement<2,80,2>;
template class TStorageEnrichmentElement<2,80,3>;
template class TStorageEnrichmentElement<2,81,2>;
template class TStorageEnrichmentElement<2,81,3>;
template class TStorageEnrichmentElement<2,82,2>;
template class TStorageEnrichmentElement<2,82,3>;
template class TStorageEnrichmentElement<2,83,2>;
template class TStorageEnrichmentElement<2,83,3>;
template class TStorageEnrichmentElement<2,84,2>;
template class TStorageEnrichmentElement<2,84,3>;
template class TStorageEnrichmentElement<2,85,2>;
template class TStorageEnrichmentElement<2,85,3>;
template class TStorageEnrichmentElement<2,86,2>;
template class TStorageEnrichmentElement<2,86,3>;
template class TStorageEnrichmentElement<2,87,2>;
template class TStorageEnrichmentElement<2,87,3>;
template class TStorageEnrichmentElement<2,88,2>;
template class TStorageEnrichmentElement<2,88,3>;
template class TStorageEnrichmentElement<2,89,2>;
template class TStorageEnrichmentElement<2,89,3>;
template class TStorageEnrichmentElement<2,90,2>;
template class TStorageEnrichmentElement<2,90,3>;
template class TStorageEnrichmentElement<2,91,2>;
template class TStorageEnrichmentElement<2,91,3>;
template class TStorageEnrichmentElement<2,92,2>;
template class TStorageEnrichmentElement<2,92,3>;
template class TStorageEnrichmentElement<2,93,2>;
template class TStorageEnrichmentElement<2,93,3>;
template class TStorageEnrichmentElement<2,94,2>;
template class TStorageEnrichmentElement<2,94,3>;
template class TStorageEnrichmentElement<2,95,2>;
template class TStorageEnrichmentElement<2,95,3>;
template class TStorageEnrichmentElement<2,96,2>;
template class TStorageEnrichmentElement<2,96,3>;
template class TStorageEnrichmentElement<2,97,2>;
template class TStorageEnrichmentElement<2,97,3>;
template class TStorageEnrichmentElement<2,98,2>;
template class TStorageEnrichmentElement<2,98,3>;
template class TStorageEnrichmentElement<2,99,2>;
template class TStorageEnrichmentElement<2,99,3>;
template class TStorageEnrichmentElement<2,100,2>;
template class TStorageEnrichmentElement<2,100,3>;
template class TStorageEnrichmentElement<3,1,2>;
template class TStorageEnrichmentElement<3,1,3>;
template class TStorageEnrichmentElement<3,2,2>;
template class TStorageEnrichmentElement<3,2,3>;
template class TStorageEnrichmentElement<3,3,2>;
template class TStorageEnrichmentElement<3,3,3>;
template class TStorageEnrichmentElement<3,4,2>;
template class TStorageEnrichmentElement<3,4,3>;
template class TStorageEnrichmentElement<3,5,2>;
template class TStorageEnrichmentElement<3,5,3>;
template class TStorageEnrichmentElement<3,6,2>;
template class TStorageEnrichmentElement<3,6,3>;
template class TStorageEnrichmentElement<3,7,2>;
template class TStorageEnrichmentElement<3,7,3>;
template class TStorageEnrichmentElement<3,8,2>;
template class TStorageEnrichmentElement<3,8,3>;
template class TStorageEnrichmentElement<3,9,2>;
template class TStorageEnrichmentElement<3,9,3>;
template class TStorageEnrichmentElement<3,10,2>;
template class TStorageEnrichmentElement<3,10,3>;
template class TStorageEnrichmentElement<3,11,2>;
template class TStorageEnrichmentElement<3,11,3>;
template class TStorageEnrichmentElement<3,12,2>;
template class TStorageEnrichmentElement<3,12,3>;
template class TStorageEnrichmentElement<3,13,2>;
template class TStorageEnrichmentElement<3,13,3>;
template class TStorageEnrichmentElement<3,14,2>;
template class TStorageEnrichmentElement<3,14,3>;
template class TStorageEnrichmentElement<3,15,2>;
template class TStorageEnrichmentElement<3,15,3>;
template class TStorageEnrichmentElement<3,16,2>;
template class TStorageEnrichmentElement<3,16,3>;
template class TStorageEnrichmentElement<3,17,2>;
template class TStorageEnrichmentElement<3,17,3>;
template class TStorageEnrichmentElement<3,18,2>;
template class TStorageEnrichmentElement<3,18,3>;
template class TStorageEnrichmentElement<3,19,2>;
template class TStorageEnrichmentElement<3,19,3>;
template class TStorageEnrichmentElement<3,20,2>;
template class TStorageEnrichmentElement<3,20,3>;
template class TStorageEnrichmentElement<3,21,2>;
template class TStorageEnrichmentElement<3,21,3>;
template class TStorageEnrichmentElement<3,22,2>;
template class TStorageEnrichmentElement<3,22,3>;
template class TStorageEnrichmentElement<3,23,2>;
template class TStorageEnrichmentElement<3,23,3>;
template class TStorageEnrichmentElement<3,24,2>;
template class TStorageEnrichmentElement<3,24,3>;
template class TStorageEnrichmentElement<3,25,2>;
template class TStorageEnrichmentElement<3,25,3>;
template class TStorageEnrichmentElement<3,26,2>;
template class TStorageEnrichmentElement<3,26,3>;
template class TStorageEnrichmentElement<3,27,2>;
template class TStorageEnrichmentElement<3,27,3>;
template class TStorageEnrichmentElement<3,28,2>;
template class TStorageEnrichmentElement<3,28,3>;
template class TStorageEnrichmentElement<3,29,2>;
template class TStorageEnrichmentElement<3,29,3>;
template class TStorageEnrichmentElement<3,30,2>;
template class TStorageEnrichmentElement<3,30,3>;
template class TStorageEnrichmentElement<3,31,2>;
template class TStorageEnrichmentElement<3,31,3>;
template class TStorageEnrichmentElement<3,32,2>;
template class TStorageEnrichmentElement<3,32,3>;
template class TStorageEnrichmentElement<3,33,2>;
template class TStorageEnrichmentElement<3,33,3>;
template class TStorageEnrichmentElement<3,34,2>;
template class TStorageEnrichmentElement<3,34,3>;
template class TStorageEnrichmentElement<3,35,2>;
template class TStorageEnrichmentElement<3,35,3>;
template class TStorageEnrichmentElement<3,36,2>;
template class TStorageEnrichmentElement<3,36,3>;
template class TStorageEnrichmentElement<3,37,2>;
template class TStorageEnrichmentElement<3,37,3>;
template class TStorageEnrichmentElement<3,38,2>;
template class TStorageEnrichmentElement<3,38,3>;
template class TStorageEnrichmentElement<3,39,2>;
template class TStorageEnrichmentElement<3,39,3>;
template class TStorageEnrichmentElement<3,40,2>;
template class TStorageEnrichmentElement<3,40,3>;
template class TStorageEnrichmentElement<3,41,2>;
template class TStorageEnrichmentElement<3,41,3>;
template class TStorageEnrichmentElement<3,42,2>;
template class TStorageEnrichmentElement<3,42,3>;
template class TStorageEnrichmentElement<3,43,2>;
template class TStorageEnrichmentElement<3,43,3>;
template class TStorageEnrichmentElement<3,44,2>;
template class TStorageEnrichmentElement<3,44,3>;
template class TStorageEnrichmentElement<3,45,2>;
template class TStorageEnrichmentElement<3,45,3>;
template class TStorageEnrichmentElement<3,46,2>;
template class TStorageEnrichmentElement<3,46,3>;
template class TStorageEnrichmentElement<3,47,2>;
template class TStorageEnrichmentElement<3,47,3>;
template class TStorageEnrichmentElement<3,48,2>;
template class TStorageEnrichmentElement<3,48,3>;
template class TStorageEnrichmentElement<3,49,2>;
template class TStorageEnrichmentElement<3,49,3>;
template class TStorageEnrichmentElement<3,50,2>;
template class TStorageEnrichmentElement<3,50,3>;
template class TStorageEnrichmentElement<3,51,2>;
template class TStorageEnrichmentElement<3,51,3>;
template class TStorageEnrichmentElement<3,52,2>;
template class TStorageEnrichmentElement<3,52,3>;
template class TStorageEnrichmentElement<3,53,2>;
template class TStorageEnrichmentElement<3,53,3>;
template class TStorageEnrichmentElement<3,54,2>;
template class TStorageEnrichmentElement<3,54,3>;
template class TStorageEnrichmentElement<3,55,2>;
template class TStorageEnrichmentElement<3,55,3>;
template class TStorageEnrichmentElement<3,56,2>;
template class TStorageEnrichmentElement<3,56,3>;
template class TStorageEnrichmentElement<3,57,2>;
template class TStorageEnrichmentElement<3,57,3>;
template class TStorageEnrichmentElement<3,58,2>;
template class TStorageEnrichmentElement<3,58,3>;
template class TStorageEnrichmentElement<3,59,2>;
template class TStorageEnrichmentElement<3,59,3>;
template class TStorageEnrichmentElement<3,60,2>;
template class TStorageEnrichmentElement<3,60,3>;
template class TStorageEnrichmentElement<3,61,2>;
template class TStorageEnrichmentElement<3,61,3>;
template class TStorageEnrichmentElement<3,62,2>;
template class TStorageEnrichmentElement<3,62,3>;
template class TStorageEnrichmentElement<3,63,2>;
template class TStorageEnrichmentElement<3,63,3>;
template class TStorageEnrichmentElement<3,64,2>;
template class TStorageEnrichmentElement<3,64,3>;
template class TStorageEnrichmentElement<3,65,2>;
template class TStorageEnrichmentElement<3,65,3>;
template class TStorageEnrichmentElement<3,66,2>;
template class TStorageEnrichmentElement<3,66,3>;
template class TStorageEnrichmentElement<3,67,2>;
template class TStorageEnrichmentElement<3,67,3>;
template class TStorageEnrichmentElement<3,68,2>;
template class TStorageEnrichmentElement<3,68,3>;
template class TStorageEnrichmentElement<3,69,2>;
template class TStorageEnrichmentElement<3,69,3>;
template class TStorageEnrichmentElement<3,70,2>;
template class TStorageEnrichmentElement<3,70,3>;
template class TStorageEnrichmentElement<3,71,2>;
template class TStorageEnrichmentElement<3,71,3>;
template class TStorageEnrichmentElement<3,72,2>;
template class TStorageEnrichmentElement<3,72,3>;
template class TStorageEnrichmentElement<3,73,2>;
template class TStorageEnrichmentElement<3,73,3>;
template class TStorageEnrichmentElement<3,74,2>;
template class TStorageEnrichmentElement<3,74,3>;
template class TStorageEnrichmentElement<3,75,2>;
template class TStorageEnrichmentElement<3,75,3>;
template class TStorageEnrichmentElement<3,76,2>;
template class TStorageEnrichmentElement<3,76,3>;
template class TStorageEnrichmentElement<3,77,2>;
template class TStorageEnrichmentElement<3,77,3>;
template class TStorageEnrichmentElement<3,78,2>;
template class TStorageEnrichmentElement<3,78,3>;
template class TStorageEnrichmentElement<3,79,2>;
template class TStorageEnrichmentElement<3,79,3>;
template class TStorageEnrichmentElement<3,80,2>;
template class TStorageEnrichmentElement<3,80,3>;
template class TStorageEnrichmentElement<3,81,2>;
template class TStorageEnrichmentElement<3,81,3>;
template class TStorageEnrichmentElement<3,82,2>;
template class TStorageEnrichmentElement<3,82,3>;
template class TStorageEnrichmentElement<3,83,2>;
template class TStorageEnrichmentElement<3,83,3>;
template class TStorageEnrichmentElement<3,84,2>;
template class TStorageEnrichmentElement<3,84,3>;
template class TStorageEnrichmentElement<3,85,2>;
template class TStorageEnrichmentElement<3,85,3>;
template class TStorageEnrichmentElement<3,86,2>;
template class TStorageEnrichmentElement<3,86,3>;
template class TStorageEnrichmentElement<3,87,2>;
template class TStorageEnrichmentElement<3,87,3>;
template class TStorageEnrichmentElement<3,88,2>;
template class TStorageEnrichmentElement<3,88,3>;
template class TStorageEnrichmentElement<3,89,2>;
template class TStorageEnrichmentElement<3,89,3>;
template class TStorageEnrichmentElement<3,90,2>;
template class TStorageEnrichmentElement<3,90,3>;
template class TStorageEnrichmentElement<3,91,2>;
template class TStorageEnrichmentElement<3,91,3>;
template class TStorageEnrichmentElement<3,92,2>;
template class TStorageEnrichmentElement<3,92,3>;
template class TStorageEnrichmentElement<3,93,2>;
template class TStorageEnrichmentElement<3,93,3>;
template class TStorageEnrichmentElement<3,94,2>;
template class TStorageEnrichmentElement<3,94,3>;
template class TStorageEnrichmentElement<3,95,2>;
template class TStorageEnrichmentElement<3,95,3>;
template class TStorageEnrichmentElement<3,96,2>;
template class TStorageEnrichmentElement<3,96,3>;
template class TStorageEnrichmentElement<3,97,2>;
template class TStorageEnrichmentElement<3,97,3>;
template class TStorageEnrichmentElement<3,98,2>;
template class TStorageEnrichmentElement<3,98,3>;
template class TStorageEnrichmentElement<3,99,2>;
template class TStorageEnrichmentElement<3,99,3>;
template class TStorageEnrichmentElement<3,100,2>;
template class TStorageEnrichmentElement<3,100,3>;
}
