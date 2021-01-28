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
#include "simple_cubic_mesh.template.cc"
#include "prolate_spheroidal_left_ventricle_domain.h"

namespace oomph{

template<class ELEMENT>
class ProlateSpheroidalLeftVentricleMesh : public virtual SimpleCubicMesh<ELEMENT>
{
public:
	ProlateSpheroidalLeftVentricleMesh(const unsigned& n_elements_per_edge,
									const double& sphere_radius,
									const double& lambda,
									TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)	:
	SimpleCubicMesh<ELEMENT>(n_elements_per_edge,
							n_elements_per_edge,
							n_elements_per_edge,
							-1.0, 1.0,
							-1.0, 1.0,
							-1.0, 1.0,
							time_stepper_pt)
	{
		// Mesh can only be built with 3D Qelements.
		MeshChecker::assert_geometric_element<QElementGeometricBase,ELEMENT>(3);

		N_Elements_Per_Edge = n_elements_per_edge;
		Sphere_Radius = sphere_radius;
		Lambda = lambda;

		// wrap_into_prolate_spheroid();
	}

private:
	unsigned N_Elements_Per_Edge;
	double Sphere_Radius;
	double Lambda;

public:
	void wrap_into_prolate_spheroid();
};


template<class ELEMENT>
class ElasticProlateSpheroidalLeftVentricleMesh :
public virtual ProlateSpheroidalLeftVentricleMesh<ELEMENT>,
// public virtual SimpleCubicMesh<ELEMENT>,
public virtual SolidMesh
{
public:
	ElasticProlateSpheroidalLeftVentricleMesh<ELEMENT>(const unsigned& n_elements_per_edge,
													const double& sphere_radius,
													const double& lambda,
													TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)	:
		SimpleCubicMesh<ELEMENT>(n_elements_per_edge,
								n_elements_per_edge,
								n_elements_per_edge,
								-1.0, 1.0,
								-1.0, 1.0,
								-1.0, 1.0,
								time_stepper_pt),
		ProlateSpheroidalLeftVentricleMesh<ELEMENT>(n_elements_per_edge,
													sphere_radius,
													lambda,
													time_stepper_pt)
	{
		set_lagrangian_nodal_coordinates();
	}
	virtual ~ElasticProlateSpheroidalLeftVentricleMesh() {}
	
};

}


#endif