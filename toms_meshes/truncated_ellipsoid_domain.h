#ifndef OOMPH_TRUNCATED_ELLIPSOID_DOMAIN_HEADER
#define OOMPH_TRUNCATED_ELLIPSOID_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph{
	
	class TruncatedEllipseDomain : public Domain
	{
	public:

		// TruncatedEllipseDomain(	const unsigned& transmural_layers,
		// 						const unsigned& azimuthal_elements,
		// 						const double& a_endo,
		// 						const double& b_endo,
		// 						const double& a_epi,
		// 						const double& b_epi,
		// 						const double& element_spread_at_apex) :
		// Transmural_layers(transmural_layers),
		// Azimuthal_elements(azimuthal_elements),
		// A_endo(a_endo),
		// B_endo(b_endo),
		// A_epi(a_epi),
		// B_epi(b_epi),
		// Element_spread_at_apex(element_spread_at_apex)
		// {
		// 	const unsigned n_macro = transmural_layers*azimuthal_elements;
		// 	Macro_element_pt.resize(n_macro);

		// 	for(unsigned i=0; i<n_macro; i++){
		// 		Macro_element_pt[i] = new QMacroElement<2>(this,i);
		// 	}
		// }

		//Second constructor
		TruncatedEllipseDomain( const unsigned& transmural_layers,
								const unsigned& azimuthal_elements,
								const double& internal_volume,
								const double& truncated_thickness,
								const double& apex_thickness,
								const double& truncated_opening,
								const double& element_spread_at_apex) :
		Transmural_layers(transmural_layers),
		Azimuthal_elements(azimuthal_elements),
		Element_spread_at_apex(element_spread_at_apex)
		{
			A_endo = 2.0*internal_volume/(MathematicalConstants::Pi*truncated_opening);
			if(A_endo<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain A_endo is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}
			B_endo = truncated_opening;
			if(B_endo<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain B_endo is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}
			A_epi  = apex_thickness + 2*internal_volume/(MathematicalConstants::Pi*truncated_opening);
			if(A_epi<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain A_epi is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}
			B_epi  = truncated_thickness + truncated_opening;
			if(B_epi<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain B_epi is zero",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}


			if(abs(A_endo - A_epi)<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain A_endo is equal to A_epi",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}

			if(abs(B_endo - B_epi)<1e-12){
				throw OomphLibError(
					"TruncatedEllipseDomain B_endo is equal to B_epi",
					OOMPH_CURRENT_FUNCTION,
					OOMPH_EXCEPTION_LOCATION);
			}


			const unsigned n_macro = transmural_layers*azimuthal_elements;
			Macro_element_pt.resize(n_macro);

			for(unsigned i=0; i<n_macro; i++){
				Macro_element_pt[i] = new QMacroElement<2>(this,i);
			}
		}



		TruncatedEllipseDomain(const TruncatedEllipseDomain&) { 
			BrokenCopy::broken_copy("TruncatedEllipseDomain");
		}

		void operator=(const TruncatedEllipseDomain&) {
			BrokenCopy::broken_assign("TruncatedEllipseDomain");
		}

		~TruncatedEllipseDomain(){
			unsigned n=nmacro_element();
			for(unsigned i=0; i<n; i++){
				delete Macro_element_pt[i];
			}
		}

		void macro_element_boundary(const unsigned& t,
									const unsigned& i_macro,
									const unsigned& i_direct,
									const Vector<double>& s,
									Vector<double>& f);


		void recover_ratio_from_x(const Vector<double>& x,
									double& sigma);

		void recover_theta_from_x(const Vector<double>& x,
									double& theta);

		void recover_r_from_x(const Vector<double>& x,
									double& r);


	private:

	unsigned Transmural_layers;
	unsigned Azimuthal_elements;
	double A_endo;
	double B_endo;
	double A_epi;
	double B_epi;
	double Element_spread_at_apex;

	};



	void TruncatedEllipseDomain::macro_element_boundary(const unsigned& t,
								const unsigned& imacro,
								const unsigned& idirect,
								const Vector<double>& s,
								Vector<double>& f)
	{
		using namespace QuadTreeNames;

		// Get coordinates of macro element
		unsigned i_theta=imacro%Azimuthal_elements;
		unsigned i_r=(imacro-i_theta)/Azimuthal_elements;


		double deriv_at_apex = Element_spread_at_apex;
		double dderiv_at_apex = 0.0;
		double a = MathematicalConstants::Pi;
		double b = -2.0*deriv_at_apex - dderiv_at_apex/2.0 + 3.0*MathematicalConstants::Pi;
		double c = 6.0*deriv_at_apex + dderiv_at_apex/2.0 - 6.0*MathematicalConstants::Pi;
		double d = 4.0*MathematicalConstants::Pi - 4.0*deriv_at_apex;

		// Angle limits
		double lo_theta_frac_of_el = double(i_theta)/double(Azimuthal_elements);
		double theta_lo=a +
						b*lo_theta_frac_of_el +
						c*pow(lo_theta_frac_of_el, 2.0) + 
						d*pow(lo_theta_frac_of_el, 3.0);

		double hi_theta_frac_of_el = double(i_theta+1)/double(Azimuthal_elements);
		double theta_hi=a +
						b*hi_theta_frac_of_el +
						c*pow(hi_theta_frac_of_el, 2.0) + 
						d*pow(hi_theta_frac_of_el, 3.0);


		// Radius limits

		double r_lo_lo_squared, r_lo_hi_squared;
		double r_hi_lo_squared, r_hi_hi_squared;

		double lo_radius_frac_of_el = double(i_r)/double(Transmural_layers);
		double hi_radius_frac_of_el = double(i_r+1)/double(Transmural_layers);


		//r_lo_lo
		double a_sigmasquared = lo_radius_frac_of_el*A_endo*A_endo + (1 - lo_radius_frac_of_el)*A_epi*A_epi;
		double b_sigmasquared = lo_radius_frac_of_el*B_endo*B_endo + (1 - lo_radius_frac_of_el)*B_epi*B_epi;
		r_lo_lo_squared = a_sigmasquared*b_sigmasquared/( b_sigmasquared*pow(sin(theta_lo),2) + a_sigmasquared*pow(cos(theta_lo),2) );
		double r_lo_lo = sqrt(r_lo_lo_squared);

		//r_lo_hi
		a_sigmasquared = hi_radius_frac_of_el*A_endo*A_endo + (1 - hi_radius_frac_of_el)*A_epi*A_epi;
		b_sigmasquared = hi_radius_frac_of_el*B_endo*B_endo + (1 - hi_radius_frac_of_el)*B_epi*B_epi;
		r_lo_hi_squared = a_sigmasquared*b_sigmasquared/( b_sigmasquared*pow(sin(theta_lo),2) + a_sigmasquared*pow(cos(theta_lo),2) );
		double r_lo_hi = sqrt(r_lo_hi_squared);

		//r_hi_lo
		a_sigmasquared = lo_radius_frac_of_el*A_endo*A_endo + (1 - lo_radius_frac_of_el)*A_epi*A_epi;
		b_sigmasquared = lo_radius_frac_of_el*B_endo*B_endo + (1 - lo_radius_frac_of_el)*B_epi*B_epi;
		r_hi_lo_squared = a_sigmasquared*b_sigmasquared/( b_sigmasquared*pow(sin(theta_hi),2) + a_sigmasquared*pow(cos(theta_hi),2) );
		double r_hi_lo = sqrt(r_hi_lo_squared);

		//r_hi_hi
		a_sigmasquared = hi_radius_frac_of_el*A_endo*A_endo + (1 - hi_radius_frac_of_el)*A_epi*A_epi;
		b_sigmasquared = hi_radius_frac_of_el*B_endo*B_endo + (1 - hi_radius_frac_of_el)*B_epi*B_epi;
		r_hi_hi_squared = a_sigmasquared*b_sigmasquared/( b_sigmasquared*pow(sin(theta_hi),2) + a_sigmasquared*pow(cos(theta_hi),2) );
		double r_hi_hi = sqrt(r_hi_hi_squared);


		// Actual radius and angle
		double r=0.0;
		double theta=0.0;

		double a_squared;
		double b_squared;

		// Which direction?
		switch(idirect){
			case N:

				theta=theta_lo+0.5*(s[0]+1.0)*(theta_hi-theta_lo);
				// r=r_lo_hi+0.5*(s[0]+1.0)*(r_hi_hi-r_lo_hi);
				a_squared = hi_radius_frac_of_el*A_endo*A_endo + (1 - hi_radius_frac_of_el)*A_epi*A_epi;
				b_squared = hi_radius_frac_of_el*B_endo*B_endo + (1 - hi_radius_frac_of_el)*B_epi*B_epi;
				r = sqrt( a_squared*b_squared/( b_squared*pow(sin(theta),2) + a_squared*pow(cos(theta),2) ) );

				break;

			case S:

				theta=theta_lo+0.5*(s[0]+1.0)*(theta_hi-theta_lo);
				// r=r_lo_lo+0.5*(s[0]+1.0)*(r_hi_lo-r_lo_lo);
				a_squared = lo_radius_frac_of_el*A_endo*A_endo + (1 - lo_radius_frac_of_el)*A_epi*A_epi;
				b_squared = lo_radius_frac_of_el*B_endo*B_endo + (1 - lo_radius_frac_of_el)*B_epi*B_epi;
				r = sqrt( a_squared*b_squared/( b_squared*pow(sin(theta),2) + a_squared*pow(cos(theta),2) ) );

				break;

			case W:

				theta=theta_lo;
				r=r_lo_lo+0.5*(s[0]+1.0)*(r_lo_hi-r_lo_lo);

				break;

			case E:

				theta=theta_hi;
				r=r_hi_lo+0.5*(s[0]+1.0)*(r_hi_hi-r_hi_lo);

				break;

			default:
				std::ostringstream error_stream;
				error_stream << "idirect is " << idirect << " not one of N, S, W, E" <<  std::endl;
				throw OomphLibError(
				error_stream.str(),
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
		}

		f[0]=r*sin(theta);
		f[1]=r*cos(theta);

	}


	void TruncatedEllipseDomain::recover_ratio_from_x(const Vector<double>& x,
																double& sigma)
	{
		if(x[0]==0){sigma = (x[1]*x[1] - B_epi*B_epi)/(B_endo*B_endo - B_epi*B_epi);}
		if(x[1]==0){sigma = (x[0]*x[0] - B_epi*B_epi)/(B_endo*B_endo - B_epi*B_epi);}

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
			std::string error_message = "Sigma is horribly broken. Value ";
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

	void TruncatedEllipseDomain::recover_theta_from_x(const Vector<double>& x,
																double& theta)
	{
		theta = asin( sqrt( x[0]*x[0]/(x[0]*x[0]+x[1]*x[1]) ) );
	}

	void TruncatedEllipseDomain::recover_r_from_x(const Vector<double>& x,
															double& r)
	{
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
	}


}

#endif