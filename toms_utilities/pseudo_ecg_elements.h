#ifndef OOMPH_PSEUDO_ECG
#define OOMPH_PSEUDO_ECG

#include "../cell_membrane_potential/cell_membrane_potential_elements.h"
#include "../generic/mesh.h"

namespace oomph{

//the default value of alpha^2*sigma_i/ 4*sigma_e
static double Default_ECG_parameter = 1.0;

	
//We provide an element-free version so we can use them free of meshes
template<unsigned DIM>
class PseudoECGEquationsBase
{
public:
	PseudoECGEquationsBase()
	{
		ECG_parameter_pt = &Default_ECG_parameter;
		Mesh_pts.resize(0);
	}

	///ECG parameter
	const double &ecg_parameter() const {return *ECG_parameter_pt;}

	/// Pointer to the ECG_parameter_pt
	double* &ecg_parameter_pt() {return ECG_parameter_pt;}

	

	//Add a mesh to the list of meshes we get ecg contribution from
	void add_tissue_mesh(Mesh* mesh_pt)
	{
		Mesh_pts.push_back(mesh_pt);
	}
protected:
	//calculate the pseudo ecg, this assumes constant intracellular and extracellular conduction in the formulation,
	// if an element is not derived from DimensionlessMembranePotentialEquationsBase this will cause an error if PARANOID is
	// set, otherwise it will cause a horrible error.
	double get_ecg_at_global_coordinate(const Vector<double>& x) const
	{
		double ecg_val = 0.0;
		for(unsigned m=0; m<Mesh_pts.size(); m++)
		{
			const unsigned n_elements = Mesh_pts[m]->nelement();
			for(unsigned e=0; e<n_elements; e++)
			{
				DimensionlessMembranePotentialEquationsBase* elem_pt=dynamic_cast<DimensionlessMembranePotentialEquationsBase*>(Mesh_pts[m]->element_pt(e));
				#ifdef PARANOID
				if(elem_pt==nullptr)
				{
					throw OomphLibError(
						    "An element in a mesh is not derived from DimensionlessMembranePotentialEquationsBase",
						    OOMPH_CURRENT_FUNCTION,
						    OOMPH_EXCEPTION_LOCATION);
				}
				#endif
				ecg_val += elem_pt->integrate_ecg(x);
			}
		}

		return (*ECG_parameter_pt)*ecg_val;
	}


private:
	//Vector of pointers to meshes which ecg is constructed from
	Vector<Mesh*> Mesh_pts;

	double* ECG_parameter_pt;


};

template<unsigned DIM>
class PseudoECGEquations : public virtual FiniteElement,
							public virtual PseudoECGEquationsBase<DIM>
{
public:
	PseudoECGEquations() : PseudoECGEquationsBase<DIM>()
	{
		
	}

	//We don't need any dofs
	inline unsigned required_nvalue(const unsigned &n) const {return 0;}
	
	//Output
	unsigned nscalar_paraview() const
	{
		return 1;
	}

	void scalar_value_paraview(std::ofstream& file_out,
								const unsigned& i,
								const unsigned& nplot) const
	{
		if(i!=0)
		{
			throw OomphLibError(
				"PseudoECGEquations element has only ecg to output",
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}

		// Loop over plot points
		unsigned num_plot_points=this->nplot_points_paraview(nplot);
		for (unsigned iplot=0;iplot<num_plot_points;iplot++){
			// Get local coordinates of plot point
			Vector<double> s(DIM);
			this->get_s_plot(iplot,nplot,s);
			// Get Eulerian coordinate of plot point
			Vector<double> x(DIM);
			this->interpolated_x(s,x);
			//return the ecg
			file_out << this->get_ecg_at_global_coordinate(x) << std::endl;
		}
	}

	void scalar_value_fct_paraview(std::ofstream& file_out,
									const unsigned& i,
									const unsigned& nplot,
									FiniteElement::SteadyExactSolutionFctPt
									exact_soln_pt) const
	{
		scalar_value_paraview(file_out,i,nplot);
	}

	std::string scalar_name_paraview(const unsigned& i) const
	{
		return "PseudoECG";
	}



	void fill_in_contribution_to_residuals(Vector<double> &residuals)
	{
	}

	void fill_in_contribution_to_jacobian(Vector<double> &residuals,
											DenseMatrix<double> &jacobian)
	{
	}

	void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
															DenseMatrix<double> &jacobian, 
															DenseMatrix<double> &mass_matrix)
	{
	}


};


template<unsigned DIM, unsigned NNODE_1D>
class QPseudoECGElement : public virtual QElement<DIM, NNODE_1D>,
							public virtual PseudoECGEquations<DIM>
{
public:
	QPseudoECGElement() : QElement<DIM, NNODE_1D>(),
							PseudoECGEquations<DIM>()
	{ }
};

template<unsigned DIM, unsigned NNODE_1D>
class TPseudoECGElement : public virtual TElement<DIM, NNODE_1D>,
							public virtual PseudoECGEquations<DIM>
{
public:
	TPseudoECGElement() : TElement<DIM, NNODE_1D>(),
							PseudoECGEquations<DIM>()
	{ }
};





template class PseudoECGEquationsBase<1>;
template class PseudoECGEquationsBase<2>;
template class PseudoECGEquationsBase<3>;


template class PseudoECGEquations<1>;
template class PseudoECGEquations<2>;
template class PseudoECGEquations<3>;


template class QPseudoECGElement<1,2>;
template class QPseudoECGElement<1,3>;
template class QPseudoECGElement<1,4>;

template class QPseudoECGElement<2,2>;
template class QPseudoECGElement<2,3>;
template class QPseudoECGElement<2,4>;

template class QPseudoECGElement<3,2>;
template class QPseudoECGElement<3,3>;
template class QPseudoECGElement<3,4>;


template class TPseudoECGElement<1,2>;
template class TPseudoECGElement<1,3>;
template class TPseudoECGElement<1,4>;

template class TPseudoECGElement<2,2>;
template class TPseudoECGElement<2,3>;
template class TPseudoECGElement<2,4>;

template class TPseudoECGElement<3,2>;
template class TPseudoECGElement<3,3>;
// template class TPseudoECGElement<3,4>;

}

#endif
