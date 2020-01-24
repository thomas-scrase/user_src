#include "anisotropic_solid_with_vector_and_monodomain_elements.h"

namespace oomph
{
	//force build the multiphysics cardiac elements
	// template class QAnisotropicWithVectAndMonoPVDElement<1,1,2>;
	// template class QAnisotropicWithVectAndMonoPVDElement<1,1,3>;
	// template class QAnisotropicWithVectAndMonoPVDElement<1,1,4>;

	template class QAnisotropicWithVectAndMonoPVDElement<2,1,2>;
	template class QAnisotropicWithVectAndMonoPVDElement<2,2,2>;
	template class QAnisotropicWithVectAndMonoPVDElement<2,1,3>;
	template class QAnisotropicWithVectAndMonoPVDElement<2,2,3>;
	// template class QAnisotropicWithVectAndMonoPVDElement<2,1,4>;
	// template class QAnisotropicWithVectAndMonoPVDElement<2,2,4>;

	template class QAnisotropicWithVectAndMonoPVDElement<3,1,2>;
	template class QAnisotropicWithVectAndMonoPVDElement<3,2,2>;
	template class QAnisotropicWithVectAndMonoPVDElement<3,3,2>;
	template class QAnisotropicWithVectAndMonoPVDElement<3,1,3>;
	template class QAnisotropicWithVectAndMonoPVDElement<3,2,3>;
	template class QAnisotropicWithVectAndMonoPVDElement<3,3,3>;
	// template class QAnisotropicWithVectAndMonoPVDElement<3,1,4>;
	// template class QAnisotropicWithVectAndMonoPVDElement<3,2,4>;
	// template class QAnisotropicWithVectAndMonoPVDElement<3,3,4>;



	// template class QAnisotropicWithVectAndMonoPVDElementWithPressure<1,1>;

	template class QAnisotropicWithVectAndMonoPVDElementWithPressure<2,1>;
	template class QAnisotropicWithVectAndMonoPVDElementWithPressure<2,2>;

	template class QAnisotropicWithVectAndMonoPVDElementWithPressure<3,1>;
	template class QAnisotropicWithVectAndMonoPVDElementWithPressure<3,2>;
	template class QAnisotropicWithVectAndMonoPVDElementWithPressure<3,3>;

	// template class QAnisotropicWithVectAndMonoPVDElementWithContinuousPressure<2,1>;
	// template class QAnisotropicWithVectAndMonoPVDElementWithContinuousPressure<2,2>;

	// template class QAnisotropicWithVectAndMonoPVDElementWithContinuousPressure<3,1>;
	// template class QAnisotropicWithVectAndMonoPVDElementWithContinuousPressure<3,2>;
	// template class QAnisotropicWithVectAndMonoPVDElementWithContinuousPressure<3,3>;



	template class TAnisotropicWithVectAndMonoPVDElement<2,1,2>;
	template class TAnisotropicWithVectAndMonoPVDElement<2,2,2>;
	template class TAnisotropicWithVectAndMonoPVDElement<2,1,3>;
	template class TAnisotropicWithVectAndMonoPVDElement<2,2,3>;

}