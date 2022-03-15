//LIC// ====================================================================
//LIC// This file contains the monodomain equations and elements, derived
//LIC// from the base cell membrane potential equations - use strang splitting method
//LIC//====================================================================

//IMPLEMENTS CRANK-NICOLSON TYPE FORMULATION OF DIFFUSION ALA ISMAIL THESIS
// SECOND ORDER ACCURATE AND STABLE SINCE IT USES BACKWARD EULER METHOD.
// 


#ifndef OOMPH_REFINEABLE_MONODOMAIN_CRANK_NICOLSON_HEADER
#define OOMPH_REFINEABLE_MONODOMAIN_CRANK_NICOLSON_HEADER

//oomph-lib headers
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"

#include "monodomain_elements_crank_nicolson.h"


namespace oomph
{

//Monodomain Equations
template <unsigned DIM>
class RefineableMonodomainEquationsCrankNicolson : 	public virtual MonodomainEquationsCrankNicolson<DIM>,
													public virtual RefineableElement,
													public virtual ElementWithZ2ErrorEstimator
{
public:

	RefineableMonodomainEquationsCrankNicolson() :	MonodomainEquationsCrankNicolson<DIM>(),
													RefineableElement(),
													ElementWithZ2ErrorEstimator()
    {

    }

    /// Broken copy constructor
	RefineableMonodomainEquationsCrankNicolson(
		const RefineableMonodomainEquationsCrankNicolson<DIM>& dummy) 
	{ 
		BrokenCopy::broken_copy("RefineableMonodomainEquationsCrankNicolson");
	}

	/// Broken assignment operator
	void operator=(const RefineableMonodomainEquationsCrankNicolson<DIM>&) 
	{
		BrokenCopy::broken_assign(
		"RefineableMonodomainEquationsCrankNicolson");
	}

	/// Number of 'flux' terms for Z2 error estimation 
	unsigned num_Z2_flux_terms() {return DIM;}

	/// \short Get 'flux' for Z2 error recovery:
	/// Standard flux.from GeneralisedAdvectionDiffusion equations
	void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
	{this->get_flux(s,flux);}

	///  Further build: Copy source function pointer from father element
	void further_build()
	{
		RefineableMonodomainEquationsCrankNicolson<DIM>* cast_father_element_pt 
			= dynamic_cast<RefineableMonodomainEquationsCrankNicolson<DIM>*>(this->father_element_pt());

		//Set the values of the pointers from the father
		this->Source_fct_pt = cast_father_element_pt->source_fct_pt();
		this->Predicted_vm_fct_pt = cast_father_element_pt->predicted_vm_fct_pt();
		this->Diff_fct_pt = cast_father_element_pt->diff_fct_pt();
		this->Cm_pt = cast_father_element_pt->cm_pt();

		//Set the ALE status
		this->ALE_is_disabled = cast_father_element_pt->ALE_is_disabled;
	}

	void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
	{
		// Set size of Vector: u
		values.resize(1);

		//Find number of nodes
		const unsigned n_node = this->nnode();

		//Find the index at which the unknown is stored
		const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();

		//Local shape function
		Shape psi(n_node);

		//Find values of shape function
		this->shape(s,psi);

		//Initialise value of u
		values[0] = 0.0;

		//Loop over the local nodes and sum
		for(unsigned l=0;l<n_node;l++)
		{
			values[0] += this->nodal_value(l,vm_nodal_index)*psi[l];
		}
	}


protected:

	//Override the fill is residuals function to do correct interpolation from hanging nodes
	//Overload the residual for the monodomain equations
	void fill_in_generic_residual_contribution_BaseCellMembranePotential(Vector<double> &residuals,
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

		//Integers used to store the local equation number and local unknown
		//indices for the residuals and jacobians
		int local_eqn=0, local_unknown=0;
		
		// Local storage for pointers to hang_info objects
		HangInfo *hang_info_pt=0, *hang_info2_pt=0;

		//Loop over the integration points
		for(unsigned ipt=0;ipt<n_intpt;ipt++)
		{

			//Assign values of s
			for(unsigned i=0;i<DIM;i++) s[i] = this->integral_pt()->knot(ipt,i);

			//Get the integral weight
			double w = this->integral_pt()->weight(ipt);

			if(w<1e-9){continue;}

			//Call the derivatives of the shape and test functions
			double J = 
			this->dshape_and_dtest_eulerian_at_knot_BaseCellMembranePotential(ipt,psi,dpsidx,test,dtestdx);

			//Premultiply the weights and the Jacobian
			double W = w*J;

			//Calculate local values of the solution and its derivatives
			//Allocate
			double interpolated_vm=0.0;
			// double dvmdt=0.0;

			double interpolated_predvm = 0.0;
			Vector<double> interpolated_x(DIM,0.0);
			Vector<double> interpolated_dvmdx(DIM,0.0);
			Vector<double> interpolated_dpredvmdx(DIM, 0.0);
			Vector<double> mesh_velocity(DIM,0.0);

			//We assume that dt is constant over the nodes
			const double dt = this->node_pt(0)->time_stepper_pt()->time_pt()->dt();

			//Calculate function value and derivatives:
			//-----------------------------------------
			// Loop over nodes
			for(unsigned l=0;l<n_node;l++) 
			{	
				//If we're feeling paranoid then we check to ensure dt is infact the same for all the nodes
				#ifdef PARANOID
				if(std::fabs(dt-this->node_pt(0)->time_stepper_pt()->time_pt()->dt())>1e-9)
				{
					throw OomphLibError(
						"dt is not the same over all the nodes, but it has to be for this implementation",
						OOMPH_CURRENT_FUNCTION,
						OOMPH_EXCEPTION_LOCATION);
				}		      			
				#endif

				//Get the value at the node
				double vm_value = this->nodal_value(l,vm_nodal_index);
				interpolated_vm += vm_value*psi(l);

				double predicted_vm_value = this->get_nodal_predicted_vm_BaseCellMembranePotential(l);
				interpolated_predvm += predicted_vm_value*psi(l);

				// Loop over directions
				for(unsigned j=0;j<DIM;j++)
				{
					interpolated_x[j] += this->nodal_position(l,j)*psi(l);
					interpolated_dvmdx[j] += vm_value*dpsidx(l,j);

					interpolated_dpredvmdx[j] += predicted_vm_value*dpsidx(l,j);
				}
			}

			//Get the predicted membrane potential and its derivatives from the cells
			//Doing it this way is more efficient then interpolating from the nodes because predicted Vm may come
			// from an external element which would require interpolation anyway. This function is implemented 
			// to interpolate the source function of this element by default but:
			// in the case of this element belonging to a cell mesh interpolates the values at the nodes given by the cells
			// in the case of multiphysics is overridden to interpolate the value in the external element
			// this->get_interpolated_predicted_vm_and_dpredicted_vm_dx(interpolated_predvm, interpolated_dpredvmdx, s, ipt, psi, dpsidx);

			// Mesh velocity?
			if (!this->ALE_is_disabled)
			{
				for(unsigned l=0;l<n_node;l++) 
				{
					for(unsigned j=0;j<DIM;j++)
					{
						mesh_velocity[j] += this->dnodal_position_dt(l,j)*psi(l);
					}
				}
			}

			//Get diffusivity tensor
			DenseMatrix<double> D(DIM,DIM,0.0);
			this->get_diff_monodomain(ipt,s,interpolated_x,D);

			// Assemble residuals and Jacobian
			//--------------------------------

			// Loop over the test functions
			for(unsigned l=0;l<n_node;l++)
			{
				//Local variables to store the number of master nodes and 
				//the weight associated with the shape function if the node is hanging
				unsigned n_master=1; double hang_weight=1.0;
				//Local bool (is the node hanging)
				bool is_node_hanging = this->node_pt(l)->is_hanging();
				//If the node is hanging, get the number of master nodes
				if(is_node_hanging)
				{
					hang_info_pt = this->node_pt(l)->hanging_pt();
					n_master = hang_info_pt->nmaster();
				}
				//Otherwise there is just one master node, the node itself
				else
				{
					n_master = 1;
				}


				//Loop over the number of master nodes
				for(unsigned m=0;m<n_master;m++)
				{
					//Get the local equation number and hang_weight
					//If the node is hanging
					if(is_node_hanging)
					{
						//Fill in residual and jacobian contribution for u
						//Set the local equation number
						local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),vm_nodal_index);
						//Read out the weight from the master node
						hang_weight = hang_info_pt->master_weight(m);
					}
					//If the node is not hanging
					else
					{
						//The local equation number comes from the node itself
						local_eqn = this->nodal_local_eqn(l,vm_nodal_index);
						//The hang weight is one
						hang_weight = 1.0;
					}

					/*IF it's not a boundary condition*/
					if(local_eqn >= 0)
					{
						residuals[local_eqn] -= (interpolated_vm - interpolated_predvm)*test(l)*W*hang_weight;

						// The Generalised Advection Diffusion bit itself
						for(unsigned k=0;k<DIM;k++)
						{
							//Terms that multiply the test function 
							double tmp = 0.0;
							// //If the mesh is moving need to subtract the mesh velocity
							if(!this->ALE_is_disabled)
							{
								tmp -= mesh_velocity[k];
							}
							tmp *= interpolated_dvmdx[k];

							//Terms that multiply the derivative of the test function
							double tmp2 = 0.0;
							//Now the diuffusive term
							for(unsigned j=0;j<DIM;j++)
							{
								tmp2 += 0.5*dt*(interpolated_dvmdx[j] + interpolated_dpredvmdx[j])*D(k,j);
							}
							//Now construct the contribution to the residuals
							residuals[local_eqn] -= (tmp*test(l) + tmp2*dtestdx(l,k))*W*hang_weight;
						}
						// Calculate the jacobian
						//-----------------------
						if(flag)
						{
							//Local variables to store the number of master nodes
							//and the weights associated with each hanging node
							unsigned n_master2=1; double hang_weight2=1.0;
							//Loop over the velocity shape functions again
							for(unsigned l2=0;l2<n_node;l2++)
							{
								//Local bool (is the node hanging)
								bool is_node2_hanging = this->node_pt(l2)->is_hanging();
								//If the node is hanging, get the number of master nodes
								if(is_node2_hanging)
								{
									hang_info2_pt = this->node_pt(l2)->hanging_pt();
									n_master2 = hang_info2_pt->nmaster();
								}
								//Otherwise there is one master node, the node itself
								else
								{
									n_master2 = 1;
								}

								//Loop over the master nodes
								for(unsigned m2=0;m2<n_master2;m2++)
								{
									//Get the local unknown and weight
									//If the node is hanging
									if(is_node2_hanging)
									{
										//Read out the local unknown from the master node
										local_unknown = this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),vm_nodal_index);
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
									if(local_unknown >= 0)
									{
										//Mass matrix term
										jacobian(local_eqn,local_unknown) 
										-= test(l)*psi(l2)*W*hang_weight*hang_weight2;

										//Add the mass matrix term
										if(flag==2){
											mass_matrix(local_eqn,local_unknown)
											+= test(l)*psi(l2)*W*hang_weight*hang_weight2;
										}

										//Add contribution to Elemental Matrix
										for(unsigned k=0;k<DIM;k++)
										{
											// jacobian(local_eqn, local_unknown) -= test(l)*psi(l2)*W*hang_weight*hang_weight2;

											//Add the mass matrix term
											if(flag==2)
											{
												mass_matrix(local_eqn,local_unknown) += test(l)*psi(l2)*W*hang_weight*hang_weight2;
											}

											//Add contribution to Elemental Matrix
											for(unsigned k=0;k<DIM;k++)
											{
												//Temporary term used in assembly
												double tmp = 0.0;
												if(!this->ALE_is_disabled)
												{
													tmp -= mesh_velocity[k];
												}
												tmp *= dpsidx(l2,k);

												double tmp2 = 0.0;
												//Now the diffusive term
												for(unsigned j=0;j<DIM;j++)
												{
													tmp2 += 0.5*dt*D(k,j)*dpsidx(l2,j);
												}

												//Now assemble Jacobian term
												jacobian(local_eqn,local_unknown) -= (tmp*test(l) + tmp2*dtestdx(l,k))*W*hang_weight*hang_weight2;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};



template<unsigned DIM, unsigned NNODE_1D>
class RefineableQMonodomainElementCrankNicolson	:	public	QMonodomainElementCrankNicolson<DIM, NNODE_1D>,
													public virtual RefineableMonodomainEquationsCrankNicolson<DIM>,
													public virtual RefineableQElement<DIM>
{
public:
	RefineableQMonodomainElementCrankNicolson()	: RefineableElement(),
												RefineableQElement<DIM>(),
												QMonodomainElementCrankNicolson<DIM, NNODE_1D>(),
												RefineableMonodomainEquationsCrankNicolson<DIM>()
	{	}

	/// Broken copy constructor
	RefineableQMonodomainElementCrankNicolson(
		const RefineableQMonodomainElementCrankNicolson<DIM,NNODE_1D>& dummy){ 
		BrokenCopy::broken_copy(
		"RefineableQMonodomainElementCrankNicolson");
	} 

	/// Broken assignment operator
	void operator=(const RefineableQMonodomainElementCrankNicolson<DIM,NNODE_1D>&){
		BrokenCopy::broken_assign(
		"RefineableQMonodomainElementCrankNicolson");
	}

	/// Number of continuously interpolated values: 1
	unsigned ncont_interpolated_values() const {return 1;}

	/// \short Number of vertex nodes in the element
	unsigned nvertex_node() const
	{
		return QMonodomainElementCrankNicolson<DIM,NNODE_1D>::nvertex_node();
	}

	/// \short Pointer to the j-th vertex node in the element
	Node* vertex_node_pt(const unsigned& j) const
	{
		return QMonodomainElementCrankNicolson<DIM,NNODE_1D>::vertex_node_pt(j);
	}

	/// Rebuild from sons: empty
	void rebuild_from_sons(Mesh* &mesh_pt) {}

	/// \short Order of recovery shape functions for Z2 error estimation:
	/// Same order as shape functions.
	unsigned nrecovery_order()
	{
		return (NNODE_1D-1);
	}

	///  \short Perform additional hanging node procedures for variables
	/// that are not interpolated by all nodes. Empty.
	void further_setup_hanging_nodes(){}
protected:
	double d2shape_and_d2test_eulerian_at_knot_BaseCellMembranePotential(const unsigned &ipt, 
																		Shape &psi, 
																		DShape &dpsidx, 
																		DShape &d2psidx,
																		Shape &test, 
																		DShape &dtestdx, 
																		DShape &d2testdx) const override
	{
		const double J = this->d2shape_eulerian_at_knot(ipt, psi, dpsidx, d2psidx);
		test = psi;
		dtestdx = dpsidx;
		d2testdx = d2psidx;

		return J;
	}
};

//=======================================================================
/// Face geometry for the 
/// RefineableQuadGeneralisedAdvectionDiffusionElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<RefineableQMonodomainElementCrankNicolson<DIM,NNODE_1D> >: 
public virtual QElement<DIM-1,NNODE_1D>
{

public:

	/// \short Constructor: Call the constructor for the
	/// appropriate lower-dimensional QElement
	FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}

};

//=======================================================================
/// Face geometry for the 1D QGeneralisedAdvectionDiffusion elements: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableQMonodomainElementCrankNicolson<1,NNODE_1D> >: 
public virtual PointElement
{

public:

	/// \short Constructor: Call the constructor for the
	/// appropriate lower-dimensional QElement
	FaceGeometry() : PointElement() {}

};


}


#endif