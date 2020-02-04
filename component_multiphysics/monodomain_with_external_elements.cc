#include "monodomain_with_external_elements.h"

namespace oomph{

	template class QMonodomainElementWithExternalCellAndSolidElements<1,2>;
	template class QMonodomainElementWithExternalCellAndSolidElements<1,3>;
	template class QMonodomainElementWithExternalCellAndSolidElements<2,2>;
	template class QMonodomainElementWithExternalCellAndSolidElements<2,3>;
	template class QMonodomainElementWithExternalCellAndSolidElements<3,2>;
	template class QMonodomainElementWithExternalCellAndSolidElements<3,3>;

	template class RefineableQMonodomainElementWithExternalCellAndSolidElements<1,2>;
	template class RefineableQMonodomainElementWithExternalCellAndSolidElements<1,3>;
	template class RefineableQMonodomainElementWithExternalCellAndSolidElements<2,2>;
	template class RefineableQMonodomainElementWithExternalCellAndSolidElements<2,3>;
	template class RefineableQMonodomainElementWithExternalCellAndSolidElements<3,2>;
	template class RefineableQMonodomainElementWithExternalCellAndSolidElements<3,3>;
	
}