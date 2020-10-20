#include "point_mesh.template.h"

namespace oomph{

	template<class ELEMENT>
	PointMesh<ELEMENT>::PointMesh(TimeStepper* time_stepper_pt){
		//This mesh has no boundaries
		set_nboundary(0);
		//This mesh has 1 element
		Element_pt.resize(1);
		//Build the element
		Element_pt[0] = new ELEMENT;
		//This mesh has 1 node
		Node_pt.resize(1);
		//Build the node
		Node_pt[0] = finite_element_pt(0)->construct_node(0,time_stepper_pt);
		node_pt(0)->x(0) = 0.0;
	}


	template<class ELEMENT>
	void PointMesh<ELEMENT>::setup_boundary_element_info(std::ostream &outfile){
		bool doc = false;
		if(outfile) { doc = true;}

		Boundary_element_pt.clear();
		Face_index_at_boundary.clear();
		Boundary_element_pt.resize(0);
		Face_index_at_boundary.resize(0);

		if(doc)
		{
			outfile << "PointMesh has no boundaries" << std::endl;
		}

		Lookup_for_elements_next_boundary_is_setup = true;
	}

}//End of namespace