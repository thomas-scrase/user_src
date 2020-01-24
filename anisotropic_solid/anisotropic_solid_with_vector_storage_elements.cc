#include "anisotropic_solid_with_vector_storage_elements.h"

namespace oomph
{
	// template class QAnisotropicWithVectPVDElement<1,1,2>;
	// template class QAnisotropicWithVectPVDElement<1,1,3>;
	// template class QAnisotropicWithVectPVDElement<1,1,4>;

	template class QAnisotropicWithVectPVDElement<2,1,2>;
	template class QAnisotropicWithVectPVDElement<2,2,2>;
	template class QAnisotropicWithVectPVDElement<2,1,3>;
	template class QAnisotropicWithVectPVDElement<2,2,3>;
	// template class QAnisotropicWithVectPVDElement<2,1,4>;
	// template class QAnisotropicWithVectPVDElement<2,2,4>;

	template class QAnisotropicWithVectPVDElement<3,1,2>;
	template class QAnisotropicWithVectPVDElement<3,2,2>;
	template class QAnisotropicWithVectPVDElement<3,3,2>;
	template class QAnisotropicWithVectPVDElement<3,1,3>;
	template class QAnisotropicWithVectPVDElement<3,2,3>;
	template class QAnisotropicWithVectPVDElement<3,3,3>;
	// template class QAnisotropicWithVectPVDElement<3,1,4>;
	// template class QAnisotropicWithVectPVDElement<3,2,4>;
	// template class QAnisotropicWithVectPVDElement<3,3,4>;



	// template class QAnisotropicWithVectPVDElementWithPressure<1,1>;

	template class QAnisotropicWithVectPVDElementWithPressure<2,1>;
	template class QAnisotropicWithVectPVDElementWithPressure<2,2>;

	template class QAnisotropicWithVectPVDElementWithPressure<3,1>;
	template class QAnisotropicWithVectPVDElementWithPressure<3,2>;
	template class QAnisotropicWithVectPVDElementWithPressure<3,3>;


	// template class QAnisotropicWithVectPVDElementWithContinuousPressure<1,1>;

	template class QAnisotropicWithVectPVDElementWithContinuousPressure<2,1>;
	template class QAnisotropicWithVectPVDElementWithContinuousPressure<2,2>;

	template class QAnisotropicWithVectPVDElementWithContinuousPressure<3,1>;
	template class QAnisotropicWithVectPVDElementWithContinuousPressure<3,2>;
	template class QAnisotropicWithVectPVDElementWithContinuousPressure<3,3>;



	// template class TAnisotropicWithVectPVDElement<1,1,2>;
	// template class TAnisotropicWithVectPVDElement<1,1,3>;
	// template class TAnisotropicWithVectPVDElement<1,1,4>;

	template class TAnisotropicWithVectPVDElement<2,1,2>;
	template class TAnisotropicWithVectPVDElement<2,2,2>;
	template class TAnisotropicWithVectPVDElement<2,1,3>;
	template class TAnisotropicWithVectPVDElement<2,2,3>;
	template class TAnisotropicWithVectPVDElement<2,1,4>;
	template class TAnisotropicWithVectPVDElement<2,2,4>;

	template class TAnisotropicWithVectPVDElement<3,1,2>;
	template class TAnisotropicWithVectPVDElement<3,2,2>;
	template class TAnisotropicWithVectPVDElement<3,3,2>;
	template class TAnisotropicWithVectPVDElement<3,1,3>;
	template class TAnisotropicWithVectPVDElement<3,2,3>;
	template class TAnisotropicWithVectPVDElement<3,3,3>;

	// template class TAnisotropicWithVectPVDElementWithContinuousPressure<1,1>;

	template class TAnisotropicWithVectPVDElementWithContinuousPressure<2,1>;
	template class TAnisotropicWithVectPVDElementWithContinuousPressure<2,2>;

	template class TAnisotropicWithVectPVDElementWithContinuousPressure<3,1>;
	template class TAnisotropicWithVectPVDElementWithContinuousPressure<3,2>;
	template class TAnisotropicWithVectPVDElementWithContinuousPressure<3,3>;
}