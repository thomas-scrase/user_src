#ifndef OOMPH_PROLATE_SPHEROIDAL_LEFT_VENTRICLE_DOMAIN_HEADER
#define OOMPH_PROLATE_SPHEROIDAL_LEFT_VENTRICLE_DOMAIN_HEADER


// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph{

	class ProlateSpheroidalLeftVentricleDomain : public Domain
	{
	public:
		ProlateSpheroidalLeftVentricleDomain(const unsigned& Transmural_elements_in_endo_region,
											const unsigned& Transmural_elements_in_mid_region,
											const unsigned& Transmural_elements_in_epi_region,
											const unsigned& Apico_basal_elements,
											const double& Phi_en,
											const double& F_en,
											const double& F_ep,
											const double& L_en,
											const double& L_ep,
											const bool& Top_is_flat=true) :
			Transmural_Elements_In_Endo_Region(Transmural_elements_in_endo_region),
			Transmural_Elements_In_Mid_Region(Transmural_elements_in_mid_region),
			Transmural_Elements_In_Epi_Region(Transmural_elements_in_epi_region),
			Apico_Basal_Elements(Apico_basal_elements),
			Circumferential_Elements(4*Apico_basal_elements),
			Phi_En(Phi_en),
			F_En(F_en),
			F_Ep(F_ep),
			L_En(L_en),
			L_Ep(L_ep),
			Top_Is_Flat(Top_is_flat)
		{
			if(Top_Is_Flat){
				Phi_Ep = acos(F_En*cosh(L_En)*cos(Phi_En)/(F_Ep*cosh(L_Ep)));
			}
			else{
				Phi_Ep = Phi_En;
			}

			unsigned nmacro = (Transmural_Elements_In_Endo_Region +
							Transmural_Elements_In_Mid_Region +
							Transmural_Elements_In_Epi_Region)*
							Apico_Basal_Elements*
							Apico_Basal_Elements*4;

			Macro_element_pt.resize(nmacro);
			// Create macro elements
			for (unsigned i=0;i<nmacro;i++)
			{
				Macro_element_pt[i]=new QMacroElement<3>(this,i);
			}
		}

		ProlateSpheroidalLeftVentricleDomain(const ProlateSpheroidalLeftVentricleDomain&) { 
			BrokenCopy::broken_copy("ProlateSpheroidalLeftVentricleDomain");
		}

		void operator=(const ProlateSpheroidalLeftVentricleDomain&) {
			BrokenCopy::broken_assign("ProlateSpheroidalLeftVentricleDomain");
		}

		~ProlateSpheroidalLeftVentricleDomain(){
			unsigned n=nmacro_element();
			for(unsigned i=0; i<n; i++){
				delete Macro_element_pt[i];
			}
		}

		void macro_element_boundary(const unsigned& t,
									const unsigned& i_macro,
									const unsigned& i_direct,
									const Vector<double>& s,
									Vector<double>& f)
		{

		}


	private:

		unsigned Transmural_Elements_In_Endo_Region;
		unsigned Transmural_Elements_In_Mid_Region;
		unsigned Transmural_Elements_In_Epi_Region;
		unsigned Apico_Basal_Elements;
		unsigned Circumferential_Elements;
		double Phi_En;
		double Phi_Ep;
		double F_En;
		double F_Ep;
		double L_En;
		double L_Ep;

		bool Top_Is_Flat;

	};

}// End namespace

#endif