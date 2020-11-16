//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations
//LIC//====================================================================

#include "bidomain_elements.h"

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
  void  BidomainEquations<DIM>::fill_in_generic_residual_contribution_BaseCellMembranePotential(Vector<double> &residuals, 
                                                                                   DenseMatrix<double> &jacobian, 
                                                                                   DenseMatrix<double> &mass_matrix,
                                                                                   unsigned flag) 
  {
   //Find out how many nodes there are
   const unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();
   const unsigned phie_nodal_index = this->phie_index_Bidomain();
     
   //Set up memory for the shape and test functions
   Shape psi(n_node), test(n_node);
   DShape dpsidx(n_node,DIM), dtestdx(n_node,DIM);
   
   //Set the value of n_intpt
   const unsigned n_intpt = this->integral_pt()->nweight();
     
   //Set the Vector to hold local coordinates
   Vector<double> s(DIM);

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
     double interpolated_phie=0.0;
     double dvmdt=0.0;
     double dphiedt=0.0;

     Vector<double> interpolated_x(DIM,0.0);
     Vector<double> interpolated_dvmdx(DIM,0.0);
     Vector<double> interpolated_phiedx(DIM,0.0);
     Vector<double> mesh_velocity(DIM,0.0);

     //Calculate function value and derivatives:
     //-----------------------------------------
     // Loop over nodes
     for(unsigned l=0;l<n_node;l++) 
      {
       //Get the value at the node
       double vm_value = this->raw_nodal_value(l,vm_nodal_index);
       double phie_value = this->raw_nodal_value(l, phie_nodal_index);
       interpolated_vm += vm_value*psi(l);
       interpolated_phie += phie_value*psi(l);
       dvmdt += this->dvm_dt_BaseCellMembranePotential(l)*psi(l);
       dphiedt += this->dphie_dt_Bidomain(l)*psi(l);
       // Loop over directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
         interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
         interpolated_phiedx[j] += phie_value*dpsidx(l,j);
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

     //Get conductivity tensors
     DenseMatrix<double> Gi(DIM,DIM,0.0);
     this->get_intracellular_conductivity_bidomain(ipt,s,interpolated_x,Gi);

     DenseMatrix<double> Ge(DIM,DIM,0.0);
     this->get_extracellular_conductivity_bidomain(ipt,s,interpolated_x,Ge);

     // Assemble residuals and Jacobian
     //--------------------------------
         
     // Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
     { 


        //Residual and jacobian for transmembrane potential dof


        //Fill in residual and jacobian contribution for u
        //Set the local equation number
        local_eqn = this->nodal_local_eqn(l,vm_nodal_index);

        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
        {
          // Add body force/source term and time derivative 
          residuals[local_eqn] -= (dvmdt + source)*test(l)*W;
         
          // The Generalised Advection Diffusion bit itself
          for(unsigned k=0;k<DIM;k++)
          {
           //Terms that multiply the test function 
            double tmp = 0.0;
           // //If the mesh is moving need to subtract the mesh velocity
           if(!this->ALE_is_disabled) {tmp -= mesh_velocity[k];}
           tmp *= interpolated_dvmdx[k];

           //Terms that multiply the derivative of the test function
            double tmp2 = 0.0;
           //Now the diuffusive term
           for(unsigned j=0;j<DIM;j++)
            {
             tmp2 += (interpolated_dvmdx[j] + interpolated_phiedx[j])*Gi(k,j);
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
              local_unknown = this->nodal_local_eqn(l2,vm_nodal_index);
             
              //If at a non-zero degree of freedom add in the entry
              if(local_unknown >= 0)
                {
                //Mass matrix term
                jacobian(local_eqn,local_unknown) 
                  -= test(l)*psi(l2)*
                  this->node_pt(l2)->time_stepper_pt()->weight(1,0)*W;

                //Add the mass matrix term
                if(flag==2)
                  {
                  mass_matrix(local_eqn,local_unknown)
                    += test(l)*psi(l2)*W;
                  }

                //Add contribution to Elemental Matrix
                for(unsigned k=0;k<DIM;k++){
                  //Temporary term used in assembly
                  double tmp = 0.0;
                  if(!this->ALE_is_disabled)
                   {tmp -= mesh_velocity[k];}
                  tmp *= dpsidx(l2,k);

                  double tmp2 = 0.0;
                  //Now the diffusive term
                  for(unsigned j=0;j<DIM;j++)
                    {
                    tmp2 += Gi(k,j)*dpsidx(l2,j);
                    }
                 
                  //Now assemble Jacobian term
                  jacobian(local_eqn,local_unknown) 
                    -= (tmp*test(l) + tmp2*dtestdx(l,k))*W;
                  }
                }

                local_unknown = this->nodal_local_eqn(l2, this->phie_index_Bidomain());
                //If at a non-zero degree of freedom add in the entry
                if(local_unknown >= 0)
                {

                  //Add the mass matrix term
                  if(flag==2)
                  {
                    mass_matrix(local_eqn,local_unknown)
                      += test(l)*psi(l2)*W;
                  }

                  //Add contribution to Elemental Matrix
                  for(unsigned k=0;k<DIM;k++){
                    double tmp2 = 0.0;
                    //Now the diffusive term
                    for(unsigned j=0;j<DIM;j++)
                    {
                      tmp2 += Gi(k,j)*dpsidx(l2,j);
                    }
                   
                    //Now assemble Jacobian term
                    jacobian(local_eqn,local_unknown) 
                      -= (tmp2*dtestdx(l,k))*W;
                  }
                
                }

              }
            }
          }


          //Residual and jacobian for extracellular potential dof


          //Fill in residual and jacobian contribution for u
          //Set the local equation number
          local_eqn = this->nodal_local_eqn(l,this->phie_index_Bidomain());

          /*IF it's not a boundary condition*/
          if(local_eqn >= 0)
          {
            // The Generalised Advection Diffusion bit itself
            for(unsigned k=0;k<DIM;k++)
              {
               //Terms that multiply the derivative of the test function
                double tmp2 = 0.0;
               //Now the diuffusive term
               for(unsigned j=0;j<DIM;j++)
                {
                 tmp2 += (interpolated_dvmdx[j] + interpolated_phiedx[j])*Gi(k,j) +
                          interpolated_phiedx[j]*Ge(k,j);
                }
               //Now construct the contribution to the residuals
               residuals[local_eqn] -= (tmp2*dtestdx(l,k))*W;
              }
           
            // Calculate the jacobian
            //-----------------------
            if(flag)
            {
              //Loop over the velocity shape functions again
              for(unsigned l2=0;l2<n_node;l2++)
              { 
                //Set the number of the unknown
                local_unknown = this->nodal_local_eqn(l2,vm_nodal_index);
               
                //If at a non-zero degree of freedom add in the entry
                if(local_unknown >= 0)
                {
                  //Add the mass matrix term
                  if(flag==2)
                  {
                  mass_matrix(local_eqn,local_unknown)
                    += test(l)*psi(l2)*W;
                  }

                  //Add contribution to Elemental Matrix
                  for(unsigned k=0;k<DIM;k++){
                    double tmp2 = 0.0;
                    //Now the diffusive term
                    for(unsigned j=0;j<DIM;j++)
                    {
                      tmp2 += Gi(k,j)*dpsidx(l2,j);
                    }
                   
                    //Now assemble Jacobian term
                    jacobian(local_eqn,local_unknown) 
                      -= (tmp2*dtestdx(l,k))*W;
                  }
                }

                local_unknown = this->nodal_local_eqn(l2, this->phie_index_Bidomain());
                //If at a non-zero degree of freedom add in the entry
                if(local_unknown >= 0)
                {
                  //Add the mass matrix term
                  if(flag==2)
                    {
                    mass_matrix(local_eqn,local_unknown)
                      += test(l)*psi(l2)*W;
                    }

                  //Add contribution to Elemental Matrix
                  for(unsigned k=0;k<DIM;k++){

                    double tmp2 = 0.0;
                    //Now the diffusive term
                    for(unsigned j=0;j<DIM;j++)
                      {
                      tmp2 += (Gi(k,j) + Ge(k,j))*dpsidx(l2,j);
                      }
                   
                    //Now assemble Jacobian term
                    jacobian(local_eqn,local_unknown)
                      -= (tmp2*dtestdx(l,k))*W;
                  }
                }
              }

            }
        }

        
      }
    } // End of loop over integration points
  }





  template class BidomainEquations<1>;
  template class BidomainEquations<2>;
  template class BidomainEquations<3>;


  template class QBidomainElement<1,2>;
  template class QBidomainElement<1,3>;
  template class QBidomainElement<1,4>;

  template class QBidomainElement<2,2>;
  template class QBidomainElement<2,3>;
  template class QBidomainElement<2,4>;

  template class QBidomainElement<3,2>;
  template class QBidomainElement<3,3>;
  template class QBidomainElement<3,4>;


  /////////////////////////////////////////////////////////////////////////
  // TBidomainElement
  /////////////////////////////////////////////////////////////////////////


  //====================================================================
  // Force build of templates
  //====================================================================
  template class TBidomainElement<1,2>;
  template class TBidomainElement<1,3>;
  template class TBidomainElement<1,4>;

  template class TBidomainElement<2,2>;
  template class TBidomainElement<2,3>;
  template class TBidomainElement<2,4>;

  template class TBidomainElement<3,2>;
  template class TBidomainElement<3,3>;



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