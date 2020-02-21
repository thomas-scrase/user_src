#include "monodomain_with_cell_elements.h"

namespace oomph
{
	template class QMonodomainWithCellElement<1, 40, 2>;
	template class QMonodomainWithCellElement<1, 40, 3>;

	template class QMonodomainWithCellElement<2, 40, 2>;
	template class QMonodomainWithCellElement<2, 40, 3>;

	template class QMonodomainWithCellElement<3, 40, 2>;
	template class QMonodomainWithCellElement<3, 40, 3>;
}