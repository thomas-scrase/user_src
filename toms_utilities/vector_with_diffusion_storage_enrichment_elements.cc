#include"vector_with_diffusion_storage_enrichment_elements.h"

namespace oomph
{
	template class VectorWithDiffusionStorageEnrichmentEquations<2>;
	template class VectorWithDiffusionStorageEnrichmentEquations<6>;
	template class VectorWithDiffusionStorageEnrichmentEquations<12>;

	template class QVectorWithDiffusionStorageEnrichmentElement<1,2>;
	template class QVectorWithDiffusionStorageEnrichmentElement<1,3>;
	template class QVectorWithDiffusionStorageEnrichmentElement<2,2>;
	template class QVectorWithDiffusionStorageEnrichmentElement<2,3>;
	template class QVectorWithDiffusionStorageEnrichmentElement<3,2>;
	template class QVectorWithDiffusionStorageEnrichmentElement<3,3>;

	template class TVectorWithDiffusionStorageEnrichmentElement<1,2>;
	template class TVectorWithDiffusionStorageEnrichmentElement<1,3>;
	template class TVectorWithDiffusionStorageEnrichmentElement<2,2>;
	template class TVectorWithDiffusionStorageEnrichmentElement<2,3>;
	template class TVectorWithDiffusionStorageEnrichmentElement<3,2>;
	template class TVectorWithDiffusionStorageEnrichmentElement<3,3>;

	template class PointVectorWithDiffusionStorageEnrichmentElement<1>;
	template class PointVectorWithDiffusionStorageEnrichmentElement<2>;
	template class PointVectorWithDiffusionStorageEnrichmentElement<3>;
}