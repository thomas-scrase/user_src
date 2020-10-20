#include "prolate_spheroidal_left_ventricle_mesh.template.h"

namespace oomph{
	
template<class ELEMENT>
void ProlateSpheroidalLeftVentricleMesh<ELEMENT>::wrap_into_prolate_spheroid(){
	//Loop over all the nodes
	const unsigned n_node = this->nnode();
	for(unsigned n=0;n<n_node;n++){
		// Pointer to node
		Node* nod_pt = this->node_pt(n);

		//First we shift the nodes such that the origin lies at the centre of the cube
		nod_pt->x(0) -= 1.0;
		nod_pt->x(1) -= 1.0;

		Vector<double> xi(3);
		xi[0] = nod_pt->x(0);
		xi[1] = nod_pt->x(1);
		xi[2] = nod_pt->x(2);


		//Preallocate memory for the radius, and circumferential, and azimuthal anglges
		double r, theta, phi;
		//Calculate the suitable new radius of the point
		r = std::max(std::max(-xi[0], xi[0]),std::max(-xi[1], xi[1]));
		//Calculate current (and final ) value of theta so we can shift x and y
		if(r>0.0){theta = arcsin(xi[1]/r);}
		else{theta = 0.0;}

		//Shift the x and y positions
		xi[0] = r*cos(theta);
		xi[1] = r*sin(theta);
	}
}

}