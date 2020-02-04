#include "storage_augmented_cell_elements.h"

namespace oomph{

	template class QStorageAugmentedCellElement<1,40,2>;
	template class QStorageAugmentedCellElement<1,40,3>;
	template class QStorageAugmentedCellElement<2,40,2>;
	template class QStorageAugmentedCellElement<2,40,3>;
	template class QStorageAugmentedCellElement<3,40,2>;
	template class QStorageAugmentedCellElement<3,40,3>;

	template class TStorageAugmentedCellElement<1,40,2>;
	template class TStorageAugmentedCellElement<1,40,3>;
	template class TStorageAugmentedCellElement<2,40,2>;
	template class TStorageAugmentedCellElement<2,40,3>;
	template class TStorageAugmentedCellElement<3,40,2>;
	template class TStorageAugmentedCellElement<3,40,3>;

	

	template class QStorageAugmentedCellElementWithExternalMonoAndSolidElements<1,40,2>;
	template class QStorageAugmentedCellElementWithExternalMonoAndSolidElements<1,40,3>;
	template class QStorageAugmentedCellElementWithExternalMonoAndSolidElements<2,40,2>;
	template class QStorageAugmentedCellElementWithExternalMonoAndSolidElements<2,40,3>;
	template class QStorageAugmentedCellElementWithExternalMonoAndSolidElements<3,40,2>;
	template class QStorageAugmentedCellElementWithExternalMonoAndSolidElements<3,40,3>;
}