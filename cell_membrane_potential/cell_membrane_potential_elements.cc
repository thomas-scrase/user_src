//LIC// ====================================================================
//LIC// This file contains the base equations for all Cell membrane potential
//LIC// equations and elements.
//LIC//  Contains all members common to both monodomain and bidomain equations
//LIC//====================================================================


//Non-inline functions for monodomain Diffusion elements
#include "cell_membrane_potential_elements.h"

namespace oomph
{

/// Default value for Peclet number
template<unsigned DIM>
double BaseCellMembranePotentialEquations<DIM>::Default_membrane_capacitance=1.0;


//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
template <unsigned DIM>
unsigned  BaseCellMembranePotentialEquations<DIM>::self_test()
{

 bool passed=true;

 // Check lower-level stuff
 if (FiniteElement::self_test()!=0)
  {
   passed=false;
  }

 // Return verdict
 if (passed)
  {
   return 0;
  }
 else
  {
   return 1;
  }
   
}


//======================================================================
/// \short Output function:
///
///   x,y,u,w_x,w_y   or    x,y,z,u,w_x,w_y,w_z
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void  BaseCellMembranePotentialEquations<DIM>::output(std::ostream &outfile, 
                                               const unsigned &nplot)
{
  // std::cout << "BOOM" << std::endl;
  //Vector of local coordinates
  Vector<double> s(DIM);

  // Tecplot header info
  outfile << tecplot_zone_string(nplot);

  const unsigned n_node = this->nnode();
  const unsigned vm_index = vm_index_BaseCellMembranePotential();
  // std::cout << "n_node " << n_node << std::endl;
  Shape psi(n_node);
  DShape dpsidx(n_node,DIM);

  // Loop over plot points
  unsigned num_plot_points=nplot_points(nplot);
  // std::cout << "Begin loop over ipt" << std::endl;
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
    // Get local coordinates of plot point
    get_s_plot(iplot,nplot,s);

    // Get Eulerian coordinate of plot point
    Vector<double> x(DIM);
    interpolated_x(s,x);

    for(unsigned i=0;i<DIM;i++) {outfile << x[i] << " ";}
    // std::cout << "outputting interpolated_u_monodomain" << std::endl;
    outfile << interpolated_vm_BaseCellMembranePotential(s) << " ";

    //Get the gradients
    (void)this->dshape_eulerian(s,psi,dpsidx);
    Vector<double> interpolated_dvmdx(DIM,0.0);
    double dvmdt = 0.0;
    for(unsigned n=0;n<n_node;n++){
      const double vm_ = this->nodal_value(n,vm_index);
      dvmdt += dvm_dt_BaseCellMembranePotential(n)*psi(n);
      for(unsigned i=0;i<DIM;i++){interpolated_dvmdx[i] += vm_*dpsidx(n,i);}
    }

    outfile << dvmdt << " ";

    for(unsigned i=0;i<DIM;i++){outfile << interpolated_dvmdx[i]  << " ";}

    outfile  << std::endl;
  }
  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile,nplot);
  // std::cout << "ENDBOOM" << std::endl;
}


//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
template <unsigned DIM>
void BaseCellMembranePotentialEquations<DIM>::output(FILE* file_pt,
                                              const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(DIM);
 
 // Tecplot header info
 fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<DIM;i++) 
    {
     fprintf(file_pt,"%g ",interpolated_x(s,i));

    }
   fprintf(file_pt,"%g \n",interpolated_vm_BaseCellMembranePotential(s));
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(file_pt,nplot);

}



//======================================================================
 /// \short  Output exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
template <unsigned DIM>
void BaseCellMembranePotentialEquations<DIM>::output_fct(std::ostream &outfile, 
             const unsigned &nplot, 
             FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {

   //Vector of local coordinates
   Vector<double> s(DIM);

   // Vector for coordintes
   Vector<double> x(DIM);

   // Tecplot header info
   outfile << tecplot_zone_string(nplot);
   
   // Exact solution Vector (here a scalar)
   Vector<double> exact_soln(1);

   // Loop over plot points
   unsigned num_plot_points=nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     
     // Get local coordinates of plot point
     get_s_plot(iplot,nplot,s);

     // Get x position as Vector
     interpolated_x(s,x);

     // Get exact solution at this point
     (*exact_soln_pt)(x,exact_soln);

     //Output x,y,...,u_exact
     for(unsigned i=0;i<DIM;i++)
      {
       outfile << x[i] << " ";
      }
     outfile << exact_soln[0] << std::endl;  
    }

   // Write tecplot footer (e.g. FE connectivity lists)
   write_tecplot_zone_footer(outfile,nplot);
   
  }




//======================================================================
 /// \short Validate against exact solution
 /// 
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
//======================================================================
template <unsigned DIM>
void BaseCellMembranePotentialEquations<DIM>::compute_error(std::ostream &outfile, 
                                           FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                                           double& error, double& norm)
{ 

 // Initialise
 error=0.0;
 norm=0.0;

 //Vector of local coordinates
 Vector<double> s(DIM);

 // Vector for coordintes
 Vector<double> x(DIM);

 //Find out how many nodes there are in the element
 unsigned n_node = nnode();

 Shape psi(n_node);

 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();
   
 // Tecplot header info
 outfile << "ZONE" << std::endl;
   
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<DIM;i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J=J_eulerian(s);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   // Get x position as Vector
   interpolated_x(s,x);

   // Get FE function value
   double vm_fe=interpolated_vm_BaseCellMembranePotential(s);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);

   //Output x,y,...,error
   for(unsigned i=0;i<DIM;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << " " << exact_soln[0]-vm_fe << std::endl;  

   // Add to error and norm
   norm+=exact_soln[0]*exact_soln[0]*W;
   error+=(exact_soln[0]-vm_fe)*(exact_soln[0]-vm_fe)*W;

  }

}

//======================================================================
/// \short Calculate the integrated value of the unknown over the element
///
//======================================================================
template <unsigned DIM>
double BaseCellMembranePotentialEquations<DIM>::integrate_vm()
{ 
 // Initialise
 double sum = 0.0;

 //Vector of local coordinates
 Vector<double> s(DIM);

 //Find out how many nodes there are in the element
 const unsigned n_node = nnode();

 //Find the index at which the concentration is stored
 const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();

 //Allocate memory for the shape functions
 Shape psi(n_node);

 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   const double w = integral_pt()->weight(ipt);
   
   //Get the shape functions
   this->shape_at_knot(ipt,psi);

   //Calculate the concentration
   double interpolated_vm = 0.0;
   for(unsigned l=0;l<n_node;l++) 
    {interpolated_vm += this->nodal_value(l,vm_nodal_index)*psi(l);}

   // Get jacobian of mapping
   const double J=J_eulerian_at_knot(ipt);

   //Add the values to the sum
   sum += interpolated_vm*w*J;
  }

 //return the sum
 return sum;
}
//const version
template <unsigned DIM>
double BaseCellMembranePotentialEquations<DIM>::integrate_vm() const{
  // Initialise
  double sum = 0.0;

  //Vector of local coordinates
  Vector<double> s(DIM);

  //Find out how many nodes there are in the element
  const unsigned n_node = nnode();

  //Find the index at which the concentration is stored
  const unsigned vm_nodal_index = this->vm_index_BaseCellMembranePotential();

  //Allocate memory for the shape functions
  Shape psi(n_node);

  //Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
    //Get the integral weight
    const double w = integral_pt()->weight(ipt);

    //Get the shape functions
    this->shape_at_knot(ipt,psi);

    //Calculate the concentration
    double interpolated_vm = 0.0;
    for(unsigned l=0;l<n_node;l++) {
      interpolated_vm += this->nodal_value(l,vm_nodal_index)*psi(l);
    }

    // Get jacobian of mapping
    const double J=J_eulerian_at_knot(ipt);

    //Add the values to the sum
    sum += interpolated_vm*w*J;
  }

 //return the sum
 return sum;
}




//====================================================================
// Force build of templates
//====================================================================
template class BaseCellMembranePotentialEquations<1>;
template class BaseCellMembranePotentialEquations<2>;
template class BaseCellMembranePotentialEquations<3>;

}