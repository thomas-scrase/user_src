//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations
//LIC//====================================================================

#include "monodomain_elements.h"

namespace oomph{
  //======================================================================
  /// \short Compute element residual Vector and/or element Jacobian matrix 
  /// 
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template <unsigned DIM>
  void  MonodomainEquations<DIM>::fill_in_generic_residual_contribution_BaseCellMembranePotential(Vector<double> &residuals, 
                                                                                   DenseMatrix<double> &jacobian, 
                                                                                   DenseMatrix<double> &mass_matrix,
                                                                                   unsigned flag) 
  {
   //Find out how many nodes there are
   const unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();
     
   //Set up memory for the shape and test functions
   Shape psi(n_node), test(n_node);
   DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
   
   //Set the value of n_intpt
   const unsigned n_intpt = this->integral_pt()->nweight();
     
   //Set the Vector to hold local coordinates
   Vector<double> s(DIM);

   //Get the membrane capacitance
   const double cm = this->cm();

   //Integers used to store the local equation number and local unknown
   //indices for the residuals and jacobians
   int local_eqn=0, local_unknown=0;

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
      // std::cout << "integration point " << ipt << std::endl;

     //Assign values of s
     for(unsigned i=0;i<DIM;i++) s[i] = this->integral_pt()->knot(ipt,i);

     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);

      if(w==0.0){continue;}
      
     //Call the derivatives of the shape and test functions
     double J = 
      this->dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(ipt,psi,dpsidx,test,dtestdx);
         
     //Premultiply the weights and the Jacobian
     double W = w*J;

     //Calculate local values of the solution and its derivatives
     //Allocate
     double interpolated_vm=0.0;
     double dvmdt=0.0;

     double interpolated_boundary_source=0.0;

     Vector<double> interpolated_x(DIM,0.0);
     Vector<double> interpolated_dvmdx(DIM,0.0);
     Vector<double> mesh_velocity(DIM,0.0);

     //Calculate function value and derivatives:
     //-----------------------------------------
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Get the value at the node
       double vm_value = this->raw_nodal_value(l,vm_nodal_index);
       interpolated_vm += vm_value*psi(l);
       dvmdt += this->dvm_dt_BaseCellMembranePotential(l)*psi(l);
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
         interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
        }

        //If Boundary_source_fct_pt has been set, get the contribution from the node
        //  This check prevents bulk non boundary elements from contributing
        //  unnecessary overhead
       if(this->Boundary_source_fct_pt){
        // Preallocate boundaries the node is on
        std::set<unsigned>* boundaries_pt;
        // Get the pointer to set of boundaries node lies on
        this->node_pt(l)->get_boundaries_pt(boundaries_pt);
        // If the set is non-zero, get a contribution to interpolated_boundary_source
        if(boundaries_pt!=0){
          double bound_source = 0.0;
          this->Boundary_source_fct_pt(boundaries_pt, bound_source);
          interpolated_boundary_source += bound_source*psi(l);
        }
       }

      }
     
     // Mesh velocity?
     if (!this->ALE_is_disabled)
      {
       for(unsigned l=0;l<n_node;l++) 
        {
         for(unsigned j=0;j<DIM;j++)
          {
           mesh_velocity[j] += this->raw_dnodal_position_dt(l,j)*psi(l);
          }
        }
      }
    

     //Get source function
     //-------------------
     double source;
     this->get_source_BaseCellMembranePotential(ipt,interpolated_x,source);

     //add the boundary source
     if(this->Boundary_source_fct_pt!=0){
      // std::cout << "Source from bound " << interpolated_boundary_source << std::endl;
      source += interpolated_boundary_source;
     }


     //Get diffusivity tensor
     DenseMatrix<double> D(DIM,DIM,0.0);
     this->get_diff_monodomain(ipt,s,interpolated_x,D);

     // Assemble residuals and Jacobian
     //--------------------------------
         
     // Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
        //Fill in residual and jacobian contribution for u
        //Set the local equation number
        // std::cout << "Vm" << std::endl;
        local_eqn = this->nodal_local_eqn(l,this->vm_index_BaseCellMembranePotential());

        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
          {
          // Add body force/source term and time derivative 
          residuals[local_eqn] -= (cm*dvmdt + source)*test(l)*W;
         
          // The Generalised Advection Diffusion bit itself
          for(unsigned k=0;k<DIM;k++)
            {
             //Terms that multiply the test function 
              double tmp = 0.0;
             // //If the mesh is moving need to subtract the mesh velocity
             if(!this->ALE_is_disabled) {tmp -= cm*mesh_velocity[k];}
             tmp *= interpolated_dvmdx[k];

             //Terms that multiply the derivative of the test function
              double tmp2 = 0.0;
             //Now the diuffusive term
             for(unsigned j=0;j<DIM;j++)
              {
               tmp2 += interpolated_dvmdx[j]*D(k,j);
              }
             //Now construct the contribution to the residuals
             residuals[local_eqn] -= (tmp*test(l) + tmp2*dtestdx(l,k))*W;
            }
         
          // Calculate the jacobian
          //-----------------------
          if(flag)
            {
            //Loop over the velocity shape functions again
            for(unsigned l2=0;l2<n_node;l2++)
              { 
              //Set the number of the unknown
              local_unknown = this->nodal_local_eqn(l2,this->vm_index_BaseCellMembranePotential());
             
              //If at a non-zero degree of freedom add in the entry
              if(local_unknown >= 0)
                {
                //Mass matrix term
                jacobian(local_eqn,local_unknown) 
                  -= cm*test(l)*psi(l2)*
                  this->node_pt(l2)->time_stepper_pt()->weight(1,0)*W;

                //Add the mass matrix term
                if(flag==2)
                  {
                  mass_matrix(local_eqn,local_unknown)
                    += cm*test(l)*psi(l2)*W;
                  }

                //Add contribution to Elemental Matrix
                for(unsigned k=0;k<DIM;k++){
                  //Temporary term used in assembly
                  double tmp = 0.0;
                  if(!this->ALE_is_disabled)
                   {tmp -= cm*mesh_velocity[k];}
                  tmp *= dpsidx(l2,k);

                  double tmp2 = 0.0;
                  //Now the diffusive term
                  for(unsigned j=0;j<DIM;j++)
                    {
                    tmp2 += D(k,j)*dpsidx(l2,j);
                    }
                 
                  //Now assemble Jacobian term
                  jacobian(local_eqn,local_unknown) 
                    -= (tmp*test(l) + tmp2*dtestdx(l,k))*W;
                  }
                }
              }
            }
          }
        }
      } // End of loop over integration points
    }





  template class MonodomainEquations<1>;
  template class MonodomainEquations<2>;
  template class MonodomainEquations<3>;


  template class QMonodomainElement<1,2>;
  template class QMonodomainElement<1,3>;
  template class QMonodomainElement<1,4>;

  template class QMonodomainElement<2,2>;
  template class QMonodomainElement<2,3>;
  template class QMonodomainElement<2,4>;

  template class QMonodomainElement<3,2>;
  template class QMonodomainElement<3,3>;
  template class QMonodomainElement<3,4>;


  /////////////////////////////////////////////////////////////////////////
  // TMonodomainElement
  /////////////////////////////////////////////////////////////////////////


  //====================================================================
  // Force build of templates
  //====================================================================
  template class TMonodomainElement<1,2>;
  template class TMonodomainElement<1,3>;
  template class TMonodomainElement<1,4>;

  template class TMonodomainElement<2,2>;
  template class TMonodomainElement<2,3>;
  template class TMonodomainElement<2,4>;

  template class TMonodomainElement<3,2>;
  template class TMonodomainElement<3,3>;



  /////////////////////////////////////////////////////////////////////////
  // PointMonodomainElement
  /////////////////////////////////////////////////////////////////////////


  //====================================================================
  // Force build of templates
  //====================================================================
  // template class PointMonodomainElement<1>;
  // template class PointMonodomainElement<2>;
  // template class PointMonodomainElement<3>;


}