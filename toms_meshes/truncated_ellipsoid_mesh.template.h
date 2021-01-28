#ifndef OOMPH_TRUNCATED_ELLIPSOID_MESH_HEADER
#define OOMPH_TRUNCATED_ELLIPSOID_MESH_HEADER

// Headers
#include "../generic/refineable_brick_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/brick_mesh.h"
#include "../generic/macro_element_node_update_element.h"

//For 2D ventricle
#include "rectangular_quadmesh.template.h"
#include "rectangular_quadmesh.template.cc"

#include "truncated_ellipsoid_domain.h"

namespace oomph
{
	template<class ELEMENT>
	class TwoDTruncatedEllipseMesh : public virtual RectangularQuadMesh<ELEMENT>
	{
	public:

		TwoDTruncatedEllipseMesh(	const unsigned& transmural_layers,
									const unsigned& azimuthal_elements,
									const double& internal_volume,
									const double& truncated_thickness,
									const double& apex_thickness,
									const double& truncated_opening,
									const double& element_spread_at_apex,
									TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)	:
		RectangularQuadMesh<ELEMENT>(azimuthal_elements,transmural_layers,1.0,1.0,time_stepper_pt)
		{
			// Mesh can only be built with 2D Qelements.
			MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

			//Set the variables describing the endo and epi ellipses
			A_endo = 2.0*internal_volume/(MathematicalConstants::Pi*truncated_opening);
			if(A_endo<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain A_endo is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}
			B_endo = truncated_opening;
			if(B_endo<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain B_endo is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}
			A_epi  = apex_thickness + 2*internal_volume/(MathematicalConstants::Pi*truncated_opening);
			if(A_epi<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain A_epi is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}
			B_epi  = truncated_thickness + truncated_opening;
			if(B_epi<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain B_epi is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}


			if(abs(A_endo - A_epi)<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain A_endo is equal to A_epi",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}

			if(abs(B_endo - B_epi)<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain B_endo is equal to B_epi",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}


			wrap_into_ellipse_shape(A_endo, B_endo, A_epi, B_epi, element_spread_at_apex);
		}

		void recover_ratio_from_x(const Vector<double>& x,
									double& sigma);

		void recover_theta_from_x(const Vector<double>& x,
									double& theta);

		void recover_r_from_x(const Vector<double>& x,
									double& r);



	private:
		void wrap_into_ellipse_shape(const double& a_endo,
									const double& b_endo,
									const double& a_epi,
									const double& b_epi,
									const double& element_spread_at_apex);

		double A_endo;
		double B_endo;
		double A_epi;
		double B_epi;

	};




	template<class ELEMENT>
	class RefineableTwoDTruncatedEllipseMesh : public virtual TwoDTruncatedEllipseMesh<ELEMENT>,
												public virtual RefineableQuadMesh<ELEMENT>
	{
	public:
		RefineableTwoDTruncatedEllipseMesh(const unsigned& transmural_layers,
											const unsigned& azimuthal_elements,
											const double& internal_volume,
											const double& truncated_thickness,
											const double& apex_thickness,
											const double& truncated_opening,
											const double& element_spread_at_apex,
											TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)	:
		RectangularQuadMesh<ELEMENT>(azimuthal_elements,transmural_layers,1.0,1.0,time_stepper_pt),
		TwoDTruncatedEllipseMesh<ELEMENT>(transmural_layers,azimuthal_elements,internal_volume,truncated_thickness,apex_thickness,truncated_opening,element_spread_at_apex,time_stepper_pt)
		{
			// Mesh can only be built with 2D Qelements.
    		MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(2);

    		this->A_endo = 2.0*internal_volume/(MathematicalConstants::Pi*truncated_opening);
			this->B_endo = truncated_opening;
			this->A_epi  = apex_thickness + 2*internal_volume/(MathematicalConstants::Pi*truncated_opening);
			this->B_epi  = truncated_thickness + truncated_opening;

    		Domain_pt = new TruncatedEllipseDomain(transmural_layers,
													azimuthal_elements,
													internal_volume,
													truncated_thickness,
													apex_thickness,
													truncated_opening,
													element_spread_at_apex);

			// Loop over all elements and set macro element pointer
			unsigned nel=this->nelement();
			for (unsigned ielem=0;ielem<nel;ielem++){
				dynamic_cast<RefineableQElement<2>*>(this->element_pt(ielem))->
					set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
			}

			// Update nodal positions based on macro-element representation
			this->node_update();

			// Nodal positions etc. were created in constructor for
			// RectangularMesh<...>. Only need to setup quadtree forest
			this->setup_quadtree_forest();
		}

	private:
		//The domain
		TruncatedEllipseDomain* Domain_pt;
	};




	template<class ELEMENT>
	class ElasticTwoDTruncatedEllipseMesh :
	public virtual TwoDTruncatedEllipseMesh<ELEMENT>,
	public virtual SolidMesh
	{
	public:
		ElasticTwoDTruncatedEllipseMesh<ELEMENT>(const unsigned& transmural_layers,
											const unsigned& azimuthal_elements,
											const double& internal_volume,
											const double& truncated_thickness,
											const double& apex_thickness,
											const double& truncated_opening,
											const double& element_spread_at_apex,
											TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper) : 
			RectangularQuadMesh<ELEMENT>(azimuthal_elements,transmural_layers,1.0,1.0,time_stepper_pt),
			TwoDTruncatedEllipseMesh<ELEMENT>(transmural_layers,azimuthal_elements,internal_volume,truncated_thickness,apex_thickness,truncated_opening,element_spread_at_apex,time_stepper_pt)
		{
			set_lagrangian_nodal_coordinates();
		}
		virtual ~ElasticTwoDTruncatedEllipseMesh()	{}
	};


	template<class ELEMENT>
	class ElasticRefineableTwoDTruncatedEllipseMesh :
	public virtual ElasticTwoDTruncatedEllipseMesh<ELEMENT>,
	public virtual RefineableQuadMesh<ELEMENT>
	{
	public:
		ElasticRefineableTwoDTruncatedEllipseMesh<ELEMENT>(const unsigned& transmural_layers,
											const unsigned& azimuthal_elements,
											const double& internal_volume,
											const double& truncated_thickness,
											const double& apex_thickness,
											const double& truncated_opening,
											const double& element_spread_at_apex,
											TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)	: 
			RectangularQuadMesh<ELEMENT>(azimuthal_elements,transmural_layers,1.0,1.0,time_stepper_pt),
			TwoDTruncatedEllipseMesh<ELEMENT>(transmural_layers,azimuthal_elements,internal_volume,truncated_thickness,apex_thickness,truncated_opening,element_spread_at_apex,time_stepper_pt),
			ElasticTwoDTruncatedEllipseMesh<ELEMENT>(transmural_layers,azimuthal_elements,internal_volume,truncated_thickness,apex_thickness,truncated_opening,element_spread_at_apex,time_stepper_pt)
		{
			this->setup_quadtree_forest();
		}
	};

}

#endif