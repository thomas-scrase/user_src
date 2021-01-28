#ifndef OOMPH_PROLATE_SPHEROIDAL_LEFT_VENTRICLE_TEMPLATE_CC
#define OOMPH_PROLATE_SPHEROIDAL_LEFT_VENTRICLE_TEMPLATE_CC

#include "prolate_spheroidal_left_ventricle_mesh.template.h"
#include "simple_cubic_mesh.template.cc"

namespace oomph{
	
template<class ELEMENT>
void ProlateSpheroidalLeftVentricleMesh<ELEMENT>::wrap_into_prolate_spheroid(){
	//Loop over all the nodes
	const unsigned n_node = this->nnode();
	for(unsigned n=0;n<n_node;n++){
		// Pointer to node
		Node* nod_pt = this->node_pt(n);

		//First we shift the nodes such that the origin lies at the centre of the cube
		// nod_pt->x(0) -= 1.0;
		// nod_pt->x(1) -= 1.0;
		// nod_pt->x(2) -= 1.0;

		Vector<double> xi(3);
		xi[0] = nod_pt->x(0);
		xi[1] = nod_pt->x(1);
		xi[2] = nod_pt->x(2);

		double r_unrotated = sqrt(xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2]);
		//Skip the point if it lies at {0,0,0}
		if(/*r_unrotated>0.0*/ false){

			//Detail the node
			oomph_info << "Node " << n << ":" << std::endl;
			oomph_info << "{" << xi[0] << "," << xi[1] << "," << xi[2] << "}" << std::endl;
			//Get the cube segment id the point lies in
			unsigned hexant_unrotated=10;
			//Corresponds to unrotated x>=y,z
			if(xi[0]>=std::fabs(xi[1]) && xi[0]>=std::fabs(xi[2])){
				hexant_unrotated = 0;
			}
			//Corresponds to unrotated -x>=y,z
			else if(-xi[0]>=std::fabs(xi[1]) && -xi[0]>=std::fabs(xi[2])){
				hexant_unrotated = 1;
			}
			//Corresponds to unrotated y>=x,z
			else if(xi[1]>=std::fabs(xi[0]) && xi[1]>=std::fabs(xi[2])){
				hexant_unrotated = 2;
			}
			//Corresponds to unrotated -y>=x,z
			else if(-xi[1]>=std::fabs(xi[0]) && -xi[1]>=std::fabs(xi[2])){
				hexant_unrotated = 3;
			}
			//Corresponds to unrotated z>=x,y
			else if(xi[2]>=std::fabs(xi[0]) && xi[2]>=std::fabs(xi[1])){
				hexant_unrotated = 4;
			}
			//Corresponds to unrotated -z>=x,y
			else if(-xi[2]>=std::fabs(xi[0]) && -xi[2]>=std::fabs(xi[1])){
				hexant_unrotated = 5;
			}

			//Undeformed face information, point and unit normal
			Vector<double> r_0(3, 0.0);
			Vector<double> n_0(3, 0.0);

			//Populate them suitably
			switch(hexant_unrotated){
				case 0:
					r_0[0] = 1.0; r_0[1] = 0.0; r_0[2] = 0.0;
					n_0[0] = 1.0; n_0[1] = 0.0; n_0[2] = 0.0;
					break;
				case 1:
					r_0[0] = -1.0; r_0[1] = 0.0; r_0[2] = 0.0;
					n_0[0] = -1.0; n_0[1] = 0.0; n_0[2] = 0.0;
					break;
				case 2:
					r_0[1] = 1.0; r_0[0] = 0.0; r_0[2] = 0.0;
					n_0[1] = 1.0; n_0[0] = 0.0; n_0[2] = 0.0;
					break;
				case 3:
					r_0[1] = -1.0; r_0[0] = 0.0; r_0[2] = 0.0;
					n_0[1] = -1.0; n_0[0] = 0.0; n_0[2] = 0.0;
					break;
				case 4:
					r_0[2] = 1.0; r_0[0] = 0.0; r_0[1] = 0.0;
					n_0[2] = 1.0; n_0[0] = 0.0; n_0[1] = 0.0;
					break;
				case 5:
					r_0[2] = -1.0; r_0[0] = 0.0; r_0[1] = 0.0;
					n_0[2] = -1.0; n_0[0] = 0.0; n_0[1] = 0.0;
					break;
				default:
					oomph_info << "Could not attribute original point to a hexant of the initial cube, hexant: " << hexant_unrotated << std::endl; exit(0);
			}

			//Get the distance to the face defined by the haxant we are in, along the directional vector to the unrotated node
			double R_surface_unrotated = ((r_0[0]*r_0[0] + r_0[1]*r_0[1] + r_0[2]*r_0[2])/
										(r_0[0]*xi[0] + r_0[1]*xi[1] + r_0[2]*xi[2]))*r_unrotated;


			//Rotate the cube so that a corner is pointing upwards
			// nod_pt->x(0) = (1.0/sqrt(2.0))*xi[0] 							- (1.0/sqrt(2.0))*xi[2];
			// nod_pt->x(1) =-(1.0/sqrt(6.0))*xi[0]  +  sqrt(2.0/3.0)*xi[1] 	- (1.0/sqrt(6.0))*xi[2];
			// nod_pt->x(2) = (1.0/sqrt(3.0))*xi[0]  + (1.0/sqrt(3.0))*xi[1]	+ (1.0/sqrt(3.0))*xi[2];


			//Get the rotated coordinates
			Vector<double> Rxi(3);
			Rxi[0] = (1.0/sqrt(2.0))*xi[0] 							- (1.0/sqrt(2.0))*xi[2];
			Rxi[1] =-(1.0/sqrt(6.0))*xi[0]  +  sqrt(2.0/3.0)*xi[1] 	- (1.0/sqrt(6.0))*xi[2];
			Rxi[2] = (1.0/sqrt(3.0))*xi[0]  + (1.0/sqrt(3.0))*xi[1]	+ (1.0/sqrt(3.0))*xi[2];

			//Get the radial distance of the rotated point
			double r_rotated = sqrt(Rxi[0]*Rxi[0] + Rxi[1]*Rxi[1] + Rxi[2]*Rxi[2]);

			//Get the rotated cube segment id the point rotated lies in
			unsigned hexant_rotated = 10;
			//Corresponds to unrotated x>=y,z
			if(Rxi[0]>=std::fabs(Rxi[1]) && Rxi[0]>=std::fabs(Rxi[2])){
				hexant_rotated = 0;
			}
			//Corresponds to unrotated -x>=y,z
			else if(-Rxi[0]>=std::fabs(Rxi[1]) && -Rxi[0]>=std::fabs(Rxi[2])){
				hexant_rotated = 1;
			}
			//Corresponds to unrotated y>=x,z
			else if(Rxi[1]>=std::fabs(Rxi[0]) && Rxi[1]>=std::fabs(Rxi[2])){
				hexant_rotated = 2;
			}
			//Corresponds to unrotated -y>=x,z
			else if(-Rxi[1]>=std::fabs(Rxi[0]) && -Rxi[1]>=std::fabs(Rxi[2])){
				hexant_rotated = 3;
			}
			//Corresponds to unrotated z>=x,y
			else if(Rxi[2]>=std::fabs(Rxi[0]) && Rxi[2]>=std::fabs(Rxi[1])){
				hexant_rotated = 4;
			}
			//Corresponds to unrotated -z>=x,y
			else if(-Rxi[2]>=std::fabs(Rxi[0]) && -Rxi[2]>=std::fabs(Rxi[1])){
				hexant_rotated = 5;
			}
			Vector<double> r_0_rotated(3, 0.0);
			Vector<double> n_0_rotated(3, 0.0);
			//Populate them suitably
			switch(hexant_rotated){
				case 0:
					r_0_rotated[0] = 1.0; r_0_rotated[1] = 0.0; r_0_rotated[2] = 0.0;
					n_0_rotated[0] = 1.0; n_0_rotated[1] = 0.0; n_0_rotated[2] = 0.0;
					break;
				case 1:
					r_0_rotated[0] = -1.0; r_0_rotated[1] = 0.0; r_0_rotated[2] = 0.0;
					n_0_rotated[0] = -1.0; n_0_rotated[1] = 0.0; n_0_rotated[2] = 0.0;
					break;
				case 2:
					r_0_rotated[1] = 1.0; r_0_rotated[0] = 0.0; r_0_rotated[2] = 0.0;
					n_0_rotated[1] = 1.0; n_0_rotated[0] = 0.0; n_0_rotated[2] = 0.0;
					break;
				case 3:
					r_0_rotated[1] = -1.0; r_0_rotated[0] = 0.0; r_0_rotated[2] = 0.0;
					n_0_rotated[1] = -1.0; n_0_rotated[0] = 0.0; n_0_rotated[2] = 0.0;
					break;
				case 4:
					r_0_rotated[2] = 1.0; r_0_rotated[0] = 0.0; r_0_rotated[1] = 0.0;
					n_0_rotated[2] = 1.0; n_0_rotated[0] = 0.0; n_0_rotated[1] = 0.0;
					break;
				case 5:
					r_0_rotated[2] = -1.0; r_0_rotated[0] = 0.0; r_0_rotated[1] = 0.0;
					n_0_rotated[2] = -1.0; n_0_rotated[0] = 0.0; n_0_rotated[1] = 0.0;
					break;
				default:
					oomph_info << "Could not attribute rotated point to a hexant of the target cube: " << hexant_rotated << std::endl; exit(0);
			}

			//Get the radial distance to the target cube along the position vector to Rxi
			double R_surface_rotated = ((r_0_rotated[0]*r_0_rotated[0] + r_0_rotated[1]*r_0_rotated[1] + r_0_rotated[2]*r_0_rotated[2])/
										(r_0_rotated[0]*Rxi[0] + r_0_rotated[1]*Rxi[1] + r_0_rotated[2]*Rxi[2]))*r_rotated;


			Vector<double> Xi_new(3, 0.0);

			Xi_new[0] = (1.0/r_rotated)*Rxi[0]*R_surface_rotated*(r_unrotated/R_surface_unrotated);
			Xi_new[1] = (1.0/r_rotated)*Rxi[1]*R_surface_rotated*(r_unrotated/R_surface_unrotated);
			Xi_new[2] = (1.0/r_rotated)*Rxi[2]*R_surface_rotated*(r_unrotated/R_surface_unrotated);

			nod_pt->x(0) = Xi_new[0];
			nod_pt->x(1) = Xi_new[1];
			nod_pt->x(2) = Xi_new[2];


			//Get the new coords
			xi[0] = nod_pt->x(0);
			xi[1] = nod_pt->x(1);
			xi[2] = nod_pt->x(2);

		}

		//THIS WORKS, DON'T CHANGE.
		// Transform the cube into a cylinder and then into a hemispherical shell
		// Deform the cube into a cylinder
		double r_xy = sqrt(xi[0]*xi[0] + xi[1]*xi[1]);
		double r_new = std::max(std::max(-xi[0],xi[0]), std::max(-xi[1],xi[1]));
		if(r_xy>0){
			nod_pt->x(0) = (r_new/r_xy)*xi[0];
			nod_pt->x(1) = (r_new/r_xy)*xi[1];
		}
		//Get the new coords
		xi[0] = nod_pt->x(0);
		xi[1] = nod_pt->x(1);
		xi[2] = nod_pt->x(2);
		r_xy = sqrt(xi[0]*xi[0] + xi[1]*xi[1]);
		//Deform the cylinder into a hemispherical shell
		double R = ((3.0-xi[2])/2.0);
		if(R<1e-12){oomph_info << "Mapping to {0,0,0} in shellify" << std::endl; exit(0);}
		double Phi = MathematicalConstants::Pi - MathematicalConstants::Pi/2.0 * r_xy;
		//Move the node to it's new location on the hemispherical shell
		if(r_xy>0.0){
			nod_pt->x(0) = R*(xi[0]/r_xy)*sin(Phi);
			nod_pt->x(1) = R*(xi[1]/r_xy)*sin(Phi);
			nod_pt->x(2) = R*cos(Phi);
		}
		else{
			nod_pt->x(0) = 0.0;
			nod_pt->x(1) = 0.0;
			nod_pt->x(2) = R*cos(Phi);
		}

		//Checking because it breaks
		//Get the new coords
		xi[0] = nod_pt->x(0);
		xi[1] = nod_pt->x(1);
		xi[2] = nod_pt->x(2);
		if(xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2] < 1e-12){oomph_info << "Mapped to {0,0,0} in shellify" << std::endl; exit(0);}
		/////////////////////////////////

	}
}

}

#endif