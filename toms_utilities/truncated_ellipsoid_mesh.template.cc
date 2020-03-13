#ifndef OOMPH_TRUNCATED_ELLIPSOID_MESH_CC
#define OOMPH_TRUNCATED_ELLIPSOID_MESH_CC

#include "truncated_ellipsoid_mesh.template.h"

namespace oomph{

	template<class ELEMENT>
	void TwoDTruncatedEllipseMesh<ELEMENT>::wrap_into_ellipse_shape(
		const double& a_endo,
		const double& b_endo,
		const double& a_epi,
		const double& b_epi,
		const double& element_spread_at_apex)
	{
		//For the local coordinate, angle and sigma (sigma used instead of r, easier for internal layers)
		//sigma = 0 coincides with endo layer
		//sigma = 1 coincides with epi layer
		Vector<double> xi(2);

		const unsigned n_node = this->nnode();
		for(unsigned n=0;n<n_node;n++){
			// Pointer to node
			Node* nod_pt = this->node_pt(n);

			double deriv_at_apex = element_spread_at_apex;
			double dderiv_at_apex = 0.0;

			double a = MathematicalConstants::Pi;
			double b = -2.0*deriv_at_apex - dderiv_at_apex/2.0 + 3.0*MathematicalConstants::Pi;
			double c = 6.0*deriv_at_apex + dderiv_at_apex/2.0 - 6.0*MathematicalConstants::Pi;
			double d = 4.0*MathematicalConstants::Pi - 4.0*deriv_at_apex;

			//Assign the local coordiantes
			//theta
			xi[0] = a + 
					b*nod_pt->x(0) + 
					c*pow(nod_pt->x(0), 2.0) + 
					d*pow(nod_pt->x(0), 3.0);

			//sigma
			xi[1] = nod_pt->x(1);


			double a_sigmasquared = xi[1]*a_endo*a_endo + (1 - xi[1])*a_epi*a_epi;
			double b_sigmasquared = xi[1]*b_endo*b_endo + (1 - xi[1])*b_epi*b_epi;
			double r_squared = a_sigmasquared*b_sigmasquared/( b_sigmasquared*pow(sin(xi[0]),2) + a_sigmasquared*pow(cos(xi[0]),2) );
			double r = sqrt(r_squared);


			//set the new coordinate of the cell
			nod_pt->x(0) = r*sin(xi[0]);
			nod_pt->x(1) = r*cos(xi[0]);

			// Set boundary coordinates
    		Vector<double> xi_bound(1);

			// Polar angle for boundary coordinate on boundary 0
			if(nod_pt->is_on_boundary(0))
			{       
			xi_bound[0]=xi[0];
			nod_pt->set_coordinates_on_boundary(0,xi_bound);
			}
			// Radius for boundary coordinate on boundary 1
			if(nod_pt->is_on_boundary(1))
			{      
			xi_bound[0]=r;
			nod_pt->set_coordinates_on_boundary(1,xi_bound);
			}
			// Polar angle for boundary coordinate on boundary 2
			if(nod_pt->is_on_boundary(2))
			{
			xi_bound[0]=xi[0];
			nod_pt->set_coordinates_on_boundary(2,xi_bound);
			}
			// Radius for boundary coordinate on boundary 3
			if(nod_pt->is_on_boundary(3))
			{      
			xi_bound[0]=r;
			nod_pt->set_coordinates_on_boundary(3,xi_bound);
			}

		}

		this->Boundary_coordinate_exists[0]=true;
		this->Boundary_coordinate_exists[1]=true;
		this->Boundary_coordinate_exists[2]=true;
		this->Boundary_coordinate_exists[3]=true;

	}

	template<class ELEMENT>
	void TwoDTruncatedEllipseMesh<ELEMENT>::recover_ratio_from_x(const Vector<double>& x,
																double& sigma)
	{
		if(x[0]==0){sigma = (x[1]*x[1] - B_epi*B_epi)/(B_endo*B_endo - B_epi*B_epi);}
		if(x[1]==0){sigma = (x[0]*x[0] - A_epi*A_epi)/(A_endo*A_endo - A_epi*A_epi);}

		if(x[0]!=0 && x[1]!=0){
			double a = (A_endo*A_endo-A_epi*A_epi)*(B_endo*B_endo-B_epi*B_epi);
			double b = (A_endo*A_endo-A_epi*A_epi)*(B_epi*B_epi - x[1]*x[1]) + (B_endo*B_endo-B_epi*B_epi)*(A_epi*A_epi - x[0]*x[0]);
			double c = A_epi*A_epi*B_epi*B_epi - B_epi*B_epi*x[0]*x[0] - A_epi*A_epi*x[1]*x[1];

			sigma = (-b - sqrt(b*b - 4*a*c))/(2*a);

			//if value is illegal (for whatever reason) try the other value
			if(sigma<(0-1.0e-12) || sigma>(1.0+1.0e-12)){
				sigma = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);
			}

			//if value is still illegal
			if(sigma<(0-1.0e-12) || sigma>(1.0+1.0e-12)){
				//Throw error, not allowed
				std::string error_message = "Sigma is horribly broken. Values ";
				error_message += std::to_string((-b - sqrt(b*b - 4.0*a*c))/(2.0*a));
				error_message += ", and ";
				error_message += std::to_string((-b + sqrt(b*b - 4.0*a*c))/(2.0*a));
				error_message += " are both not permitted";
				throw OomphLibError(error_message,
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
			}


		}

		if(sigma<(0-1.0e-12) || sigma>(1.0+1.0e-12)){
			//Throw error, not allowed
			std::string error_message = "Sigma is horribly broken. Values ";
			error_message += std::to_string(sigma);
			error_message += " is both not permitted";
			throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
		}

		// sigma = (-sqrt(
		// 			pow((B_endo*B_endo*A_epi*A_epi - B_endo*B_endo*x[0]*x[0] - 2*B_epi*B_epi*A_epi*A_epi + B_epi*B_epi*A_endo*A_endo + B_epi*B_epi*x[0]*x[0] + A_epi*A_epi*x[1]*x[1] - A_endo*A_endo*x[1]*x[1]),2) -
		// 			4*(-B_endo*B_endo*A_epi*A_epi + B_endo*B_endo*A_endo*A_endo + B_epi*B_epi*A_epi*A_epi - B_epi*B_epi*A_endo*A_endo)*(B_epi*B_epi*A_epi*A_epi - B_epi*B_epi*x[0]*x[0] - A_epi*A_epi*x[1]*x[1])
		// 			) +
		// 			B_endo*B_endo*A_epi*A_epi -
		// 			B_endo*B_endo*x[0]*x[0] -
		// 			2*B_epi*B_epi*A_epi*A_epi +
		// 			B_epi*B_epi*A_endo*A_endo +
		// 			B_epi*B_epi*x[0]*x[0] +
		// 			A_epi*A_epi*x[1]*x[1] -
		// 			A_endo*A_endo*x[1]*x[1]
		// 		)/
		// 		(
		// 			2*(B_endo*B_endo - B_epi*B_epi)*(A_epi*A_epi - A_endo*A_endo)
		// 		);


	}

	template<class ELEMENT>
	void TwoDTruncatedEllipseMesh<ELEMENT>::recover_theta_from_x(const Vector<double>& x,
																double& theta)
	{
		theta = asin( sqrt( x[0]*x[0]/(x[0]*x[0]+x[1]*x[1]) ) );
	}

	template<class ELEMENT>
	void TwoDTruncatedEllipseMesh<ELEMENT>::recover_r_from_x(const Vector<double>& x,
															double& r)
	{
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
	}


}

#endif