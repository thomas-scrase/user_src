#include "monodomain_with_d_elements.h"

namespace oomph
{
	template class QMonodomainWithDELement<1,2>;
	template class QMonodomainWithDELement<1,3>;
	template class QMonodomainWithDELement<2,2>;
	template class QMonodomainWithDELement<2,3>;
	template class QMonodomainWithDELement<3,2>;
	template class QMonodomainWithDELement<3,3>;

	template class TMonodomainWithDELement<1,2>;
	template class TMonodomainWithDELement<1,3>;
	template class TMonodomainWithDELement<2,2>;
	template class TMonodomainWithDELement<2,3>;
	template class TMonodomainWithDELement<3,2>;
	template class TMonodomainWithDELement<3,3>;
}