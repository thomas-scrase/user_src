#ifndef OOMPH_PROLATE_SPHEROIDAL_LEFT_VENTRICLE_HEADER
#define OOMPH_PROLATE_SPHEROIDAL_LEFT_VENTRICLE_HEADER

// Headers
#include "../generic/refineable_brick_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/brick_mesh.h"
#include "../generic/macro_element_node_update_element.h"

#include "simple_cubic_mesh.template.h"
#include "prolate_spheroidal_left_ventricle_domain.h"

namespace oomph{

template<class ELEMENT>
class ProlateSpheroidalLeftVentricleMesh : public virtual SimpleCubicMesh<ELEMENT>
{
public:
	ProlateSpheroidalLeftVentricleMesh(const unsigned& Transmural_elements_in_endo_region,
									const unsigned& Transmural_elements_in_mid_region,
									const unsigned& Transmural_elements_in_epi_region,
									const unsigned& Apico_basal_elements,
									const double& Phi_en,
									const double& F_en,
									const double& F_ep,
									const double& L_en,
									const double& L_ep,
									const bool& Top_is_flat=true,
									TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)	:
	SimpleCubicMesh<ELEMENT>(2*Apico_basal_elements,
							2*Apico_basal_elements,
							Transmural_elements_in_endo_region+Transmural_elements_in_mid_region+Transmural_elements_in_epi_region,
							2.0,2.0,2.0,
							time_stepper_pt)
	{
		// Mesh can only be built with 3D Qelements.
		MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(3);

		//Assign the parameters which define the epithelial, enothelial, and basal boundaries
		Phi_En = Phi_en;
		if(Top_Is_Flat){
			Phi_Ep = acos(F_En*cosh(L_En)*cos(Phi_En)/(F_Ep*cosh(L_Ep)));
		}
		else{
			Phi_Ep = Phi_En;
		}
		F_En = F_en;
		F_Ep = F_ep;
		L_En = L_en;
		L_Ep = L_ep;
		Top_Is_Flat = Top_is_flat;

		wrap_into_prolate_spheroid();
	}

private:
	double Phi_En;
	double Phi_Ep;
	double F_En;
	double F_Ep;
	double L_En;
	double L_Ep;

	bool Top_Is_Flat;

	void wrap_into_prolate_spheroid();
};

}


#endif