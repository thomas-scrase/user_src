#ifndef	OOMPH_ELECTRO_MECHANICAL_MULTIPHYSICS_ELEMENTS_HEADER
#define OOMPH_ELECTRO_MECHANICAL_MULTIPHYSICS_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif



#include "../toms_utilities/diff_augmented_cell_wrapper.h"

//For matrix inverse
#include "../anisotropic_constitutive/anisotropic_constitutive_laws.h"

//Solid elements for the external solid element for geometric data
#include "../anisotropic_solid/anisotropic_solid_elements.h"


namespace oomph
{
	//Cell solver with external solid element

	template<class CELL_SOLVER_ELEMENT, class EXT_SOLID_ELEMENT>
	class CellSolverWithExternalSolidElement :
		public virtual CELL_SOLVER_ELEMENT,
		public virtual ElementWithExternalElement
	{
	public:
		CellSolverWithExternalSolidElement() :
		CELL_SOLVER_ELEMENT(),
		ElementWithExternalElement()
		{
			ElementWithExternalElement::set_ninteraction(1);
		}

		void fill_in_contribution_to_residuals(oomph::Vector<double>& residuals, oomph::DenseMatrix<double>& jacobian)
		{
			CELL_SOLVER_ELEMENT::fill_in_contribution_to_residuals(residuals, jacobian);
		}

		void fill_in_contribution_to_jacobian(oomph::Vector<double>& residuals, oomph::DenseMatrix<double>& jacobian)
		{
			CELL_SOLVER_ELEMENT::fill_in_contribution_to_jacobian(residuals, jacobian);
		}

		void identify_field_data_for_interactions(std::set<std::pair<Data*,unsigned> > &paired_field_data) override
		{
			FiniteElement::identify_field_data_for_interactions(paired_field_data);
		}

		void describe_local_dofs(std::ostream &out, const std::string &current_string) const
		{
			CELL_SOLVER_ELEMENT::describe_local_dofs(out, current_string);
		}

		//Get the deformed metric tensor of the external solid element at the local integral point ipt.
		// we assume that G is suitably sized and zeroed before the function is called
		void get_external_interpolated_deformed_metric_tensor(const unsigned& ipt, DenseMatrix<double>& G)
		{
			//Pointer to the external solid element
			EXT_SOLID_ELEMENT* ext_elem_pt = dynamic_cast<EXT_SOLID_ELEMENT*>(this->external_element_pt(0, ipt));

			//Local coordinate of the integral point within the external element
			Vector<double> s(this->dim(), 0.0);
			s = this->external_element_local_coord(0, ipt);


			//Find out how many nodes there are
			unsigned n_node = ext_elem_pt->nnode();

			//Find out how many positional dofs there are
			unsigned n_position_type = ext_elem_pt->nnodal_position_type();

			//Set up memory for the shape functions
			Shape psi(n_node,n_position_type);
			DShape dpsidxi(n_node,n_position_type,this->dim());

			//Call the derivatives of the shape functions (ignore Jacobian)
			(void) ext_elem_pt->dshape_lagrangian(s,psi,dpsidxi);


			//Calculate interpolated values of the derivative of global position
			//wrt lagrangian coordinates
			DenseMatrix<double> interpolated_G(this->dim());

			//Initialise to zero
			for(unsigned i=0;i<this->dim();i++){for(unsigned j=0;j<this->dim();j++) {interpolated_G(i,j) = 0.0;}}

			//Calculate displacements and derivatives
			for(unsigned l=0;l<n_node;l++)
			{
				//Loop over positional dofs
				for(unsigned k=0;k<n_position_type;k++)
				{
					//Loop over displacement components (deformed position)
					for(unsigned i=0;i<this->dim();i++)
					{
						//Loop over derivative directions
						for(unsigned j=0;j<this->dim();j++)
						{
							interpolated_G(j,i) += ext_elem_pt->nodal_position_gen(l,k,i)*dpsidxi(l,k,j);
						}
					}
				}
			}

			//Assign values of G
			for(unsigned i=0;i<this->dim();i++)
			{
				//Do upper half of matrix
				//Note that j must be signed here for the comparison test to work
				//Also i must be cast to an int
				for(int j=(this->dim()-1);j>=static_cast<int>(i);j--)
				{
					//Initialise G(i,j) to zero
					G(i,j) = 0.0;
					//Now calculate the dot product
					for(unsigned k=0;k<this->dim();k++)
					{
						G(i,j) += interpolated_G(i,k)*interpolated_G(j,k);
					}
				}
				//Matrix is symmetric so just copy lower half
				for(int j=(i-1);j>=0;j--)
				{
					G(i,j) = G(j,i);
				}
			}
		}
	};

	//Specific implementations

	//If the diffusion equations are MonodomainEquations
	template<unsigned DIM, template<class> class CELL_SOLVER_ELEMENT, class EXT_CELL_ELEMENT>
	class CellSolverWithExternalSolidElement<CELL_SOLVER_ELEMENT<MonodomainEquations<DIM>>, EXT_CELL_ELEMENT> :
		public virtual CELL_SOLVER_ELEMENT<MonodomainEquations<DIM>>,
		public virtual ElementWithExternalElement
	{
		//Overwrite the get diff function to take it instead from the values stored in the diffusion augmentation wrapper
		inline void get_diff_monodomain(const unsigned& ipt,
	                                            const Vector<double> &s,
	                                            const Vector<double>& x,
	                                            DenseMatrix<double>& D) const
		{
			DenseMatrix<double> Dundef(DIM, DIM, 0.0);
			DiffAugmentedCell<MonodomainEquations<DIM>>::get_interpolated_diffusion_matrix(s, Dundef);

			DenseMatrix<double> G(DIM, DIM, 0.0);
			this->get_external_interpolated_deformed_metric_tensor(ipt, G);

			//Post multiply D by the inverse of G
			DenseMatrix<double> Ginv(DIM, DIM, 0.0);
			AnisotropicConstitutiveLaw::calculate_contravariant(G, Ginv);

			//Postmultiply D by Ginv
			for(unsigned i=0; i<DIM; i++)
			{
				for(unsigned j=0; j<DIM; j++)
				{
					//Make sure it is zeroed
					D(i,j) = 0.0;
					//Calculate the entry
					for(unsigned k=0; k<DIM; k++)
					{
						D(i,j) += Dundef(i,k)*Ginv(k,j);
					}
				}
			}
		}
	};


	//If the diffusion equations are Monodomain, implicit simpsons-cubic hermite spline method
	template<unsigned DIM, template<class> class CELL_SOLVER_ELEMENT, class EXT_CELL_ELEMENT>
	class CellSolverWithExternalSolidElement<CELL_SOLVER_ELEMENT<MonodomainEquationsSimpsonsApproximation<DIM>>, EXT_CELL_ELEMENT> :
		public virtual CELL_SOLVER_ELEMENT<MonodomainEquationsSimpsonsApproximation<DIM>>,
		public virtual ElementWithExternalElement
	{
		//Overwrite the get diff function to take it instead from the values stored in the diffusion augmentation wrapper
		inline void get_diff_monodomain(const unsigned& ipt,
	                                            const Vector<double> &s,
	                                            const Vector<double>& x,
	                                            DenseMatrix<double>& D) const
		{
			DenseMatrix<double> Dundef(DIM, DIM, 0.0);
			DiffAugmentedCell<MonodomainEquationsSimpsonsApproximation<DIM>>::get_interpolated_diffusion_matrix(s, Dundef);

			DenseMatrix<double> G(DIM, DIM, 0.0);
			this->get_external_interpolated_deformed_metric_tensor(ipt, G);

			//Post multiply D by the inverse of G
			DenseMatrix<double> Ginv(DIM, DIM, 0.0);
			AnisotropicConstitutiveLaw::calculate_contravariant(G, Ginv);

			//Postmultiply D by Ginv
			for(unsigned i=0; i<DIM; i++)
			{
				for(unsigned j=0; j<DIM; j++)
				{
					//Make sure it is zeroed
					D(i,j) = 0.0;
					//Calculate the entry
					for(unsigned k=0; k<DIM; k++)
					{
						D(i,j) += Dundef(i,k)*Ginv(k,j);
					}
				}
			}
		}
	};


	//If the diffusion equations are Monodomain, implicit simpsons-cubic hermite spline method
	template<unsigned DIM, template<class> class CELL_SOLVER_ELEMENT, class EXT_CELL_ELEMENT>
	class CellSolverWithExternalSolidElement<CELL_SOLVER_ELEMENT<MonodomainEquationsStrangSplitting<DIM>>, EXT_CELL_ELEMENT> :
		public virtual CELL_SOLVER_ELEMENT<MonodomainEquationsStrangSplitting<DIM>>,
		public virtual ElementWithExternalElement
	{
		//Overwrite the get diff function to take it instead from the values stored in the diffusion augmentation wrapper
		inline void get_diff_monodomain(const unsigned& ipt,
	                                            const Vector<double> &s,
	                                            const Vector<double>& x,
	                                            DenseMatrix<double>& D) const
		{
			DenseMatrix<double> Dundef(DIM, DIM, 0.0);
			DiffAugmentedCell<MonodomainEquationsStrangSplitting<DIM>>::get_interpolated_diffusion_matrix(s, Dundef);

			DenseMatrix<double> G(DIM, DIM, 0.0);
			this->get_external_interpolated_deformed_metric_tensor(ipt, G);

			//Post multiply D by the inverse of G
			DenseMatrix<double> Ginv(DIM, DIM, 0.0);
			AnisotropicConstitutiveLaw::calculate_contravariant(G, Ginv);

			//Postmultiply D by Ginv
			for(unsigned i=0; i<DIM; i++)
			{
				for(unsigned j=0; j<DIM; j++)
				{
					//Make sure it is zeroed
					D(i,j) = 0.0;
					//Calculate the entry
					for(unsigned k=0; k<DIM; k++)
					{
						D(i,j) += Dundef(i,k)*Ginv(k,j);
					}
				}
			}
		}
	};












	//Solid with external cell solver element

	template<class SOLID_ELEMENT, class EXT_CELL_SOLVER_ELEMENT>
	class SolidElementWithExternalCellElement :
		public virtual SOLID_ELEMENT,
		public virtual ElementWithExternalElement
	{
	public:
		SolidElementWithExternalCellElement() :
		SOLID_ELEMENT(),
		ElementWithExternalElement()
		{
			ElementWithExternalElement::set_ninteraction(1);
		}

		void fill_in_contribution_to_residuals(oomph::Vector<double>& residuals, oomph::DenseMatrix<double>& jacobian)
		{
			SOLID_ELEMENT::fill_in_contribution_to_residuals(residuals, jacobian);
		}

		void fill_in_contribution_to_jacobian(oomph::Vector<double>& residuals, oomph::DenseMatrix<double>& jacobian)
		{
			SOLID_ELEMENT::fill_in_contribution_to_jacobian(residuals, jacobian);
		}

		void identify_field_data_for_interactions(std::set<std::pair<Data*,unsigned> > &paired_field_data) override
		{
			FiniteElement::identify_field_data_for_interactions(paired_field_data);
		}

		void describe_local_dofs(std::ostream &out, const std::string &current_string) const
		{
			SOLID_ELEMENT::describe_local_dofs(out, current_string);
		}

		//Get the anisotropic matrix - f,s,n
		inline void anisotropic_matrix(const unsigned& ipt,
										const Vector<double> &s,
										const Vector<double>& xi,
										DenseMatrix<double>& A)
		{
			throw OomphLibError(
					"I cannot get anisotropic matrix if the external element is not wrapped in DiffAugmentedCell",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
		}

		//Get the active strain
		inline void driving_strain(const unsigned& ipt,
											const Vector<double>& s,
											const Vector<double>& xi,
											Vector<double>& V)
		{
			//Although this function doesn't explicitly depend on the DiffAugmentedCell wrapper, I am leaving it blank
		}
	};



	//Specific implementations

	//If the cell solver is wrapped in DiffAugmeneted solver
	template<class SOLID_ELEMENT, class EXT_CELL_SOLVER_ELEMENT>
	class SolidElementWithExternalCellElement<SOLID_ELEMENT, DiffAugmentedCell<EXT_CELL_SOLVER_ELEMENT>>  :
	public virtual SOLID_ELEMENT,
	public virtual ElementWithExternalElement
	{
	public:
		SolidElementWithExternalCellElement() :
		SOLID_ELEMENT(),
		ElementWithExternalElement()
		{
			ElementWithExternalElement::set_ninteraction(1);

			//We require that the external element is derived from DiffAugmentedCell
		}

		void fill_in_contribution_to_residuals(oomph::Vector<double>& residuals, oomph::DenseMatrix<double>& jacobian)
		{
			SOLID_ELEMENT::fill_in_contribution_to_residuals(residuals, jacobian);
		}

		void fill_in_contribution_to_jacobian(oomph::Vector<double>& residuals, oomph::DenseMatrix<double>& jacobian)
		{
			SOLID_ELEMENT::fill_in_contribution_to_jacobian(residuals, jacobian);
		}

		void identify_field_data_for_interactions(std::set<std::pair<Data*,unsigned> > &paired_field_data) override
		{
			FiniteElement::identify_field_data_for_interactions(paired_field_data);
		}

		void describe_local_dofs(std::ostream &out, const std::string &current_string) const
		{
			SOLID_ELEMENT::describe_local_dofs(out, current_string);
		}

		//Get the anisotropic matrix - f,s,n
		inline void anisotropic_matrix(const unsigned& ipt,
										const Vector<double> &s,
										const Vector<double>& xi,
										DenseMatrix<double>& A)
		{
			//Pointer to the external solid element
			DiffAugmentedCell<EXT_CELL_SOLVER_ELEMENT>* ext_elem_pt = dynamic_cast<DiffAugmentedCell<EXT_CELL_SOLVER_ELEMENT>*>(this->external_element_pt(0, ipt));

			//Local coordinate of the integral point within the external element
			Vector<double> s_ext(this->dim(), 0.0);
			s_ext = this->external_element_local_coord(0, ipt);

			//Get the vectors from the external element
			ext_elem_pt->get_interpolated_preferential_vectors(s_ext, A);
		}

		//Get the active strain
		inline void driving_strain(const unsigned& ipt,
											const Vector<double>& s,
											const Vector<double>& xi,
											Vector<double>& V)
		{
			//Pointer to the external solid element
			DiffAugmentedCell<EXT_CELL_SOLVER_ELEMENT>* ext_elem_pt = dynamic_cast<DiffAugmentedCell<EXT_CELL_SOLVER_ELEMENT>*>(this->external_element_pt(0, ipt));

			//Local coordinate of the integral point within the external element
			Vector<double> s_ext(this->dim(), 0.0);
			s_ext = this->external_element_local_coord(0, ipt);

			//Zero the active strains
			V.resize(this->dim(), 0.0);

			//Get the active strain from the external element
			V[0] = ext_elem_pt->get_interpolated_active_strain_from_cell_model(s_ext);
		}


	};


}

#endif