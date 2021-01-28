#include "refineable_bidomain_elements.h"

namespace oomph
{
	template <unsigned DIM>
	void RefineableBidomainEquations<DIM>::
	fill_in_generic_residual_contribution_BaseCellMembranePotential(Vector<double> &residuals,
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
		const unsigned n_intpt = integral_pt()->nweight();
 
		//Set the Vector to hold local coordinates
		Vector<double> s(DIM);

		//Integers used to store the local equation number and local unknown
		//indices for the residuals and jacobians
		int local_eqn=0, local_unknown=0;

		// Local storage for pointers to hang_info objects
		HangInfo *hang_info_pt=0, *hang_info2_pt=0;

		//Local variable to determine the ALE stuff
		bool ALE_is_disabled_flag = this->ALE_is_disabled;

		//Loop over the integration points
		for(unsigned ipt=0;ipt<n_intpt;ipt++){
			// std::cout << "integration point " << ipt << std::endl;

			//Assign values of s
			for(unsigned i=0;i<DIM;i++) s[i] = integral_pt()->knot(ipt,i);

			//Get the integral weight
			double w = integral_pt()->weight(ipt);
		
			if(w==0.0){continue;}

			//Call the derivatives of the shape and test functions
			double J = this->dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(ipt,psi,dpsidx,test,dtestdx);

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
			for(unsigned l=0;l<n_node;l++){
				//Get the value at the node
				double vm_value = this->nodal_value(l,vm_nodal_index);
				double phie_value = this->raw_nodal_value(l, phie_nodal_index);
				interpolated_vm += vm_value*psi(l);
				interpolated_phie += phie_value*psi(l);
				dvmdt += this->dvm_dt_BaseCellMembranePotential(l)*psi(l);
				dphiedt += this->dphie_dt_Bidomain(l)*psi(l);
				// Loop over directions
				for(unsigned j=0;j<DIM;j++){
					interpolated_x[j] = this->nodal_position(l,j)*psi(l);
					interpolated_dvmdx[j] += vm_value*dpsidx(l,j);
					interpolated_phiedx[j] += phie_value*dpsidx(l,j);
				}
			}

			// Mesh velocity?
			if (!ALE_is_disabled_flag){
				for(unsigned l=0;l<n_node;l++){
					for(unsigned j=0;j<DIM;j++){
						mesh_velocity[j] += this->dnodal_position_dt(l,j)*psi(l);
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
			for(unsigned l=0;l<n_node;l++){
				//Local variables to store the number of master nodes and 
				//the weight associated with the shape function if the node is hanging
				unsigned n_master=1; double hang_weight=1.0;
				//Local bool (is the node hanging)
				bool is_node_hanging = this->node_pt(l)->is_hanging();
				//If the node is hanging, get the number of master nodes
				if(is_node_hanging){
					hang_info_pt = this->node_pt(l)->hanging_pt();
					n_master = hang_info_pt->nmaster();
				}
				//Otherwise there is just one master node, the node itself
				else {
					n_master = 1;
				}

				//Loop over the number of master nodes
				for(unsigned m=0;m<n_master;m++){


					//Residual and jacobian for transmembrane potential dof


					//Get the local equation number and hang_weight
					//If the node is hanging
					if(is_node_hanging){
						//Read out the local equation from the master node
						local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),vm_nodal_index);
						//Read out the weight from the master node
						hang_weight = hang_info_pt->master_weight(m);
					}
					//If the node is not hanging
					else{
						//The local equation number comes from the node itself
						local_eqn = this->nodal_local_eqn(l,vm_nodal_index);
						//The hang weight is one
						hang_weight = 1.0;
					}

					//IF it's not a boundary condition
					if(local_eqn >= 0){
						// Add body force/source term and time derivative 
						residuals[local_eqn] -= (dvmdt + source)*test(l)*W*hang_weight;
						// The Generalised Advection Diffusion bit itself
						for(unsigned k=0;k<DIM;k++){
							//Terms that multiply the test function 
							double tmp = 0.0;
							// //If the mesh is moving need to subtract the mesh velocity
							if(!ALE_is_disabled_flag){
								tmp -= mesh_velocity[k];
							}
							tmp *= interpolated_dvmdx[k];
							//Terms that multiply the derivative of the test function
							double tmp2 = 0.0;
							//Now the diuffusive term
							for(unsigned j=0;j<DIM;j++){
								tmp2 += (interpolated_dvmdx[j] + interpolated_phiedx[j])*Gi(k,j);
							}
							//Now construct the contribution to the residuals
							residuals[local_eqn] -= (tmp*test(l) + tmp2*dtestdx(l,k))*W*hang_weight;
						}
						// Calculate the jacobian
						if(flag){
							//Local variables to store the number of master nodes
							//and the weights associated with each hanging node
							unsigned n_master2=1; double hang_weight2=1.0;
							//Loop over the velocity shape functions again
							for(unsigned l2=0;l2<n_node;l2++){
								//Local bool (is the node hanging)
								bool is_node2_hanging = this->node_pt(l2)->is_hanging();
								//If the node is hanging, get the number of master nodes
								if(is_node2_hanging){
									hang_info2_pt = this->node_pt(l2)->hanging_pt();
									n_master2 = hang_info2_pt->nmaster();
								}
								//Otherwise there is one master node, the node itself
								else{
									n_master2 = 1;
								}

								//Loop over the master nodes
								for(unsigned m2=0;m2<n_master2;m2++){
									//Get the local unknown and weight
									//If the node is hanging
									if(is_node2_hanging){
										//Read out the local unknown from the master node
										local_unknown = 
										this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
										vm_nodal_index);
										//Read out the hanging weight from the master node
										hang_weight2 = hang_info2_pt->master_weight(m2);
									}
									//If the node is not hanging
									else{
										//The local unknown number comes from the node
										local_unknown = this->nodal_local_eqn(l2,vm_nodal_index);
										//The hang weight is one
										hang_weight2 = 1.0;
									}

									//If at a non-zero degree of freedom add in the entry
									if(local_unknown >= 0){
										//Mass matrix term
										jacobian(local_eqn,local_unknown) 
										-= test(l)*psi(l2)*
										node_pt(l2)->time_stepper_pt()->weight(1,0)
										*W*hang_weight*hang_weight2;

										//Add the mass matrix term
										if(flag==2){
											mass_matrix(local_eqn,local_unknown)
											+= test(l)*psi(l2)*W*hang_weight*hang_weight2;
										}

										//Add contribution to Elemental Matrix
										for(unsigned k=0;k<DIM;k++){
											//Temporary term used in assembly
											double tmp = 0.0;
											if(!ALE_is_disabled_flag){tmp -= mesh_velocity[k];}
											tmp *= dpsidx(l2,k);

											double tmp2 = 0.0;
											//Now the diffusive term
											for(unsigned j=0;j<DIM;j++){
												tmp2 += Gi(k,j)*dpsidx(l2,j);
											}

											//Now assemble Jacobian term
											jacobian(local_eqn,local_unknown) 
											-= (tmp*test(l) + tmp2*dtestdx(l,k))
											*W*hang_weight*hang_weight2;
										}
									}//End if off diagonal is not fixed



									//Get the local unknown and weight
									//If the node is hanging
									if(is_node2_hanging){
										//Read out the local unknown from the master node
										local_unknown = 
										this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
										phie_nodal_index);
										//Read out the hanging weight from the master node
										hang_weight2 = hang_info2_pt->master_weight(m2);
									}
									//If the node is not hanging
									else{
										//The local unknown number comes from the node
										local_unknown = this->nodal_local_eqn(l2,phie_nodal_index);
										//The hang weight is one
										hang_weight2 = 1.0;
									}

									//If at a non-zero degree of freedom add in the entry
									if(local_unknown >= 0){

										//Add the mass matrix term
										if(flag==2){
											mass_matrix(local_eqn,local_unknown)
											+= test(l)*psi(l2)*W*hang_weight*hang_weight2;
										}

										//Add contribution to Elemental Matrix
										for(unsigned k=0;k<DIM;k++){
											double tmp2 = 0.0;
											//Now the diffusive term
											for(unsigned j=0;j<DIM;j++){
												tmp2 += Gi(k,j)*dpsidx(l2,j);
											}

											//Now assemble Jacobian term
											jacobian(local_eqn,local_unknown) 
											-= (tmp2*dtestdx(l,k))*W*hang_weight*hang_weight2;
										}
									}//End if off diagonal is not fixed

								}//End Loop over master nodes (jacobian)
							}//End loop over nodes
						}//End of jacobian calculation
					}//End of non-zero equation




					//Residual and jacobian for extracellular potential dof

					//Get the local equation number and hang_weight
					//If the node is hanging
					if(is_node_hanging){
						//Read out the local equation from the master node
						local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),phie_nodal_index);
						//Read out the weight from the master node
						hang_weight = hang_info_pt->master_weight(m);
					}
					//If the node is not hanging
					else{
						//The local equation number comes from the node itself
						local_eqn = this->nodal_local_eqn(l,phie_nodal_index);
						//The hang weight is one
						hang_weight = 1.0;
					}

					//IF it's not a boundary condition
					if(local_eqn >= 0){
						// The Generalised Advection Diffusion bit itself
						for(unsigned k=0;k<DIM;k++){
							//Terms that multiply the derivative of the test function
							double tmp2 = 0.0;
							//Now the diuffusive term
							for(unsigned j=0;j<DIM;j++){
								tmp2 += (interpolated_dvmdx[j] + interpolated_phiedx[j])*Gi(k,j) +
                          				interpolated_phiedx[j]*Ge(k,j);
							}
							//Now construct the contribution to the residuals
							residuals[local_eqn] -= (tmp2*dtestdx(l,k))*W*hang_weight;
						}
						// Calculate the jacobian
						if(flag){
							//Local variables to store the number of master nodes
							//and the weights associated with each hanging node
							unsigned n_master2=1; double hang_weight2=1.0;
							//Loop over the velocity shape functions again
							for(unsigned l2=0;l2<n_node;l2++){
								//Local bool (is the node hanging)
								bool is_node2_hanging = this->node_pt(l2)->is_hanging();
								//If the node is hanging, get the number of master nodes
								if(is_node2_hanging){
									hang_info2_pt = this->node_pt(l2)->hanging_pt();
									n_master2 = hang_info2_pt->nmaster();
								}
								//Otherwise there is one master node, the node itself
								else{
									n_master2 = 1;
								}

								//Loop over the master nodes
								for(unsigned m2=0;m2<n_master2;m2++){
									//Get the local unknown and weight
									//If the node is hanging
									if(is_node2_hanging){
										//Read out the local unknown from the master node
										local_unknown = 
										this->local_hang_eqn(hang_info2_pt->master_node_pt(m2), vm_nodal_index);
										//Read out the hanging weight from the master node
										hang_weight2 = hang_info2_pt->master_weight(m2);
									}
									//If the node is not hanging
									else{
										//The local unknown number comes from the node
										local_unknown = this->nodal_local_eqn(l2,vm_nodal_index);
										//The hang weight is one
										hang_weight2 = 1.0;
									}

									//If at a non-zero degree of freedom add in the entry
									if(local_unknown >= 0){

										//Add the mass matrix term
										if(flag==2){
											mass_matrix(local_eqn,local_unknown)
											+= test(l)*psi(l2)*W*hang_weight*hang_weight2;
										}

										//Add contribution to Elemental Matrix
										for(unsigned k=0;k<DIM;k++){

											double tmp2 = 0.0;
											//Now the diffusive term
											for(unsigned j=0;j<DIM;j++){
												tmp2 += Gi(k,j)*dpsidx(l2,j);
											}

											//Now assemble Jacobian term
											jacobian(local_eqn,local_unknown) 
											-= (tmp2*dtestdx(l,k))
											*W*hang_weight*hang_weight2;
										}
									}//End if off diagonal is not fixed



									//Get the local unknown and weight
									//If the node is hanging
									if(is_node2_hanging){
										//Read out the local unknown from the master node
										local_unknown = 
										this->local_hang_eqn(hang_info2_pt->master_node_pt(m2), phie_nodal_index);
										//Read out the hanging weight from the master node
										hang_weight2 = hang_info2_pt->master_weight(m2);
									}
									//If the node is not hanging
									else{
										//The local unknown number comes from the node
										local_unknown = this->nodal_local_eqn(l2,phie_nodal_index);
										//The hang weight is one
										hang_weight2 = 1.0;
									}

									//If at a non-zero degree of freedom add in the entry
									if(local_unknown >= 0){

										//Add the mass matrix term
										if(flag==2){
											mass_matrix(local_eqn,local_unknown)
											+= test(l)*psi(l2)*W*hang_weight*hang_weight2;
										}

										//Add contribution to Elemental Matrix
										for(unsigned k=0;k<DIM;k++){
											double tmp2 = 0.0;
											//Now the diffusive term
											for(unsigned j=0;j<DIM;j++){
												tmp2 += (Gi(k,j) + Ge(k,j))*dpsidx(l2,j);
											}

											//Now assemble Jacobian term
											jacobian(local_eqn,local_unknown) 
											-= (tmp2*dtestdx(l,k))*W*hang_weight*hang_weight2;
										}
									}//End if off diagonal is not fixed

								}//End Loop over master nodes (jacobian)
							}//End loop over nodes
						}//End of jacobian calculation
					}//End of non-zero equation




				}//End loop over master nodes (residual)
			}//End loop over nodes
		}//End loop over integration points
	}//End fill in generic contribution



//====================================================================
// Force build of templates
//====================================================================
template class RefineableQBidomainElement<1,2>;
template class RefineableQBidomainElement<1,3>;
template class RefineableQBidomainElement<1,4>;

template class RefineableQBidomainElement<2,2>;
template class RefineableQBidomainElement<2,3>;
template class RefineableQBidomainElement<2,4>;

template class RefineableQBidomainElement<3,2>;
template class RefineableQBidomainElement<3,3>;
template class RefineableQBidomainElement<3,4>;

}