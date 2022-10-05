#ifndef OOMPH_ANISOTROPIC_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_ANISOTROPIC_ELASTICITY_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "../solid/solid_elements.h"
#include "../anisotropic_constitutive/anisotropic_constitutive_laws.h"


namespace oomph
{


  template <unsigned DIM>
  class AnisotropicPVDEquationsBase : public virtual PVDEquationsBase<DIM>
  {
  public:

    //Constructor
    AnisotropicPVDEquationsBase() : PVDEquationsBase<DIM>()
    {}

    //We use an anisotropic constitutive law instead of an isotropic default Oomph-lib one
    AnisotropicConstitutiveLaw* &anisotropic_constitutive_law_pt() {return Anisotropic_constitutive_law_pt;}

    //Add some function templates:
    //  preferential vectors  
    typedef void (*PreferentialVectorsFctPt)(const unsigned& ipt,
                                            const Vector<double>& s,
                                            const Vector<double>& x,
                                            Vector<Vector<double>>& A);

    //Get strains in the DIM directions, it is assumed that these are in the direction of
    // the preferential vectors
    typedef void (*DrivingStrainFctPt)(const unsigned& ipt,
                                      const Vector<double>& s,
                                      const Vector<double>& x,
                                      Vector<double>& V);

    //Access functions for preferential vectors function pointers
    PreferentialVectorsFctPt& preferential_vectors_fct_pt() {return Preferential_Vectors_fct_pt;}
    PreferentialVectorsFctPt preferential_vectors_fct_pt() const {return Preferential_Vectors_fct_pt;}

    //Access functions for driving strain function pointers
    DrivingStrainFctPt& driving_strain_fct_pt() {return Driving_strain_fct_pt;}
    DrivingStrainFctPt driving_strain_fct_pt() const {return Driving_strain_fct_pt;}

    virtual inline void preferential_vectors(const unsigned& ipt,
                                            const Vector<double> &s,
                                            const Vector<double>& xi,
                                            Vector<Vector<double>>& A)
    {
      if(Preferential_Vectors_fct_pt==0)
      {
        A.resize(DIM);
        for(unsigned i=0; i<DIM; i++)
        {
          A[i].resize(DIM);
          for(unsigned j=0; j<DIM; j++)
          {
            (A[i])[j] = 0.0;
          }
        }
      }
      else
      {
        (*Preferential_Vectors_fct_pt)(ipt, s, xi, A);

        #ifdef PARANOID
        for(unsigned i=0; i<A.size(); i++)
        {
          double len = 0;
          for(unsigned j=0; j<A[i].size(); j++)
          {
            len += (A[i])[j]*(A[i])[j];
          }
          if(std::fabs(len-1.0)>1e-9)
          {
            throw OomphLibError("Preferential vector in anisotropic solid is not normal",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
        #endif
      }
    }


    virtual inline void driving_strain(const unsigned& ipt,
                                        const Vector<double>& s,
                                        const Vector<double>& xi,
                                        Vector<double>& V)
    {
      if(Driving_strain_fct_pt==0)
      {
        V.resize(DIM);
        for(unsigned i=0; i<DIM; i++)
        {
          V[i] = 0.0;
        }
      }
      else
      {
        (*Driving_strain_fct_pt)(ipt, s, xi, V);
      }
    }


  protected:
    AnisotropicConstitutiveLaw *Anisotropic_constitutive_law_pt;

  private:
    PreferentialVectorsFctPt Preferential_Vectors_fct_pt;

    DrivingStrainFctPt Driving_strain_fct_pt; 
  };

  template<unsigned DIM>
  class AnisotropicPVDEquations : public virtual PVDEquations<DIM>,
                                  public virtual AnisotropicPVDEquationsBase<DIM>
  {
  public:
    AnisotropicPVDEquations() : AnisotropicPVDEquationsBase<DIM>()
    {}

    //We need to re-implement get_stress since we need it to also compute
    // and pass to the anisotropic constitutive law, the preferential vectors
    // and active strains
    void get_stress(const Vector<double> &s, DenseMatrix<double> &sigma);


    unsigned nscalar_paraview() const
    {
      return (DIM+1+(DIM*DIM));
    }

    void scalar_value_paraview(std::ofstream& file_out,
                              const unsigned& i,
                              const unsigned& nplot) const
    {
      Vector<double> x(DIM);
      Vector<double> xi(DIM);
      Vector<double> s(DIM);
      DenseMatrix<double> stress_or_strain(DIM,DIM);

      // Loop over plot points
      unsigned num_plot_points=this->nplot_points_paraview(nplot);
      for (unsigned iplot=0;iplot<num_plot_points;iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot,nplot,s);

        // Get Eulerian and Lagrangian coordinates
        this->interpolated_x(s,x);
        this->interpolated_xi(s,xi);

        //xi
        if(i<DIM)
        {
          file_out << xi[i] << std::endl;
          continue;
        }
        //isotropic growth
        if(i<DIM+1)
        {
          // Get isotropic growth
          double gamma;
          // Dummy integration point
          unsigned ipt=0;
          this->get_isotropic_growth(ipt,s,xi,gamma);
          file_out << gamma << std::endl;
          continue;
        }
        //strain
        if(i<(DIM+1 + (DIM*DIM)))
        {
          this->get_strain(s,stress_or_strain);
          file_out << stress_or_strain((i-DIM-1)%DIM,unsigned((i-DIM-1)/DIM)) << std::endl;
          continue;
        }
        //Can't include stress because get_stress is not a const function
        // //stress
        // if(i<(2*DIM+1 + (DIM*DIM)))
        // {
        //   this->get_stress(s,stress_or_strain);
        //   file_out << stress_or_strain((i-(DIM+1 + (DIM*DIM)))%DIM,unsigned((i-(DIM+1 + (DIM*DIM)))/DIM)) << std::endl;
        //   return;
        // }

        //Never get here
        #ifdef PARANOID
        std::stringstream error_stream;
        error_stream << "paraview anisotropic solid output died" << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        #endif
      }
    }

    void scalar_value_fct_paraview(std::ofstream& file_out,
    const unsigned& i,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt
    exact_soln_pt) const
    {
      scalar_value_paraview(file_out,i,nplot);
    }

    std::string scalar_name_paraview(const unsigned& i) const
    {
      //xi
      if(i<DIM)
      {
        return ("xi"+std::to_string(i));
      }
      //isotropic growth
      if(i<DIM+1)
      {
        return "Isotropic growth";
      }
      //strain
      if(i<(DIM+1+(DIM*DIM)))
      {
        return ("Strain"+std::to_string(i-(DIM+1)));
      }

      // if(i<(DIM+1+2*(DIM*DIM)))
      // {
      //   return ("Stress"+std::to_string(i-(DIM+1 + (DIM*DIM))));
      // }

      throw OomphLibError("Requested variable index is too large",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    }

    

  protected:
    // //We need to override this since get_stress(g,G,sigma) does not have positional
    // // information as an argument so cannot be overridden itself to get the
    // // preferential vectors or active strain. Overriding but only 4 lines are
    // // changed in any significant way
    void fill_in_generic_contribution_to_residuals_pvd(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag) override;

    //We have to reimplement get_stress(g,G,sigma (,etc)) to also accept the
    // preferential vectors and active strains
    inline void get_stress(const DenseMatrix<double> &g, 
                          const DenseMatrix<double> &G,
                          const Vector<Vector<double>>& A,
                          const Vector<double> &V,
                          DenseMatrix<double> &sigma)
    {
    #ifdef PARANOID
      //If the pointer to the constitutive law hasn't been set, issue an error
      if(this->Anisotropic_constitutive_law_pt==0)
      {
        //Write an error message
        std::string error_message = "Elements derived from AnisotropicPVDEquations must have a constitutive law:\n";
        error_message += "set one using the anisotropic_constitutive_law_pt() member function";
        //Throw the error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    #endif
    this->Anisotropic_constitutive_law_pt->calculate_second_piola_kirchhoff_stress(g,G,A,V,sigma);
    } 

    inline void get_d_stress_dG_upper(const DenseMatrix<double> &g, 
                                      const DenseMatrix<double> &G,
                                      const Vector<Vector<double>>& A,
                                      const Vector<double> &V,
                                      const DenseMatrix<double> &sigma,
                                      RankFourTensor<double> &d_sigma_dG)
    {
    #ifdef PARANOID
      //If the pointer to the constitutive law hasn't been set, issue an error
      if(this->Anisotropic_constitutive_law_pt==0)
      {
        //Write an error message
        std::string error_message = "Elements derived from AnisotropicPVDEquations must have a constitutive law:\n";
        error_message += "set one using the anisotropic_constitutive_law_pt() member function";
        //Throw the error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    #endif
      //Only bother with the symmetric part by passing false as last entry
      this->Anisotropic_constitutive_law_pt->calculate_d_second_piola_kirchhoff_stress_dG(g,G,A,V,sigma,d_sigma_dG,false);
    }
  };

  template<unsigned DIM, unsigned NNODE_1D>
  class QAnisotropicPVDElement : public virtual QPVDElement<DIM,NNODE_1D>,
                                 public virtual AnisotropicPVDEquations<DIM>
  {
  public:
    /// Constructor, there are no internal data points
    QAnisotropicPVDElement() : QPVDElement<DIM,NNODE_1D>(), AnisotropicPVDEquations<DIM>() { }

    //Final overrider
    void get_stress(const Vector<double> &s, DenseMatrix<double> &sigma)
    {
      AnisotropicPVDEquations<DIM>::get_stress(s, sigma);
    }

    /// Output function
    void output(std::ostream &outfile) {AnisotropicPVDEquations<DIM>::output(outfile);}

    /// Output function
    void output(std::ostream &outfile, const unsigned &n_plot){AnisotropicPVDEquations<DIM>::output(outfile,n_plot);}


    /// C-style output function
    void output(FILE* file_pt) {AnisotropicPVDEquations<DIM>::output(file_pt);}

    /// C-style output function
    void output(FILE* file_pt, const unsigned &n_plot){AnisotropicPVDEquations<DIM>::output(file_pt,n_plot);}
  };

  //Face geometries

  template<unsigned NNODE_1D>
  class FaceGeometry<QAnisotropicPVDElement<2,NNODE_1D> > :
  public virtual SolidQElement<1,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<1,NNODE_1D>() {}
  };

  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<QAnisotropicPVDElement<2,NNODE_1D> > >:
  public virtual PointElement
  {
  public:
    //Make sure that we call the constructor of the SolidQElement
    //Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };


  template<unsigned NNODE_1D>
  class FaceGeometry<QAnisotropicPVDElement<3,NNODE_1D> > :
  public virtual SolidQElement<2,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<2,NNODE_1D>() {}
  };

  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<QAnisotropicPVDElement<3,NNODE_1D> > > :
  public virtual SolidQElement<1,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<1,NNODE_1D>() {}
  };


  template<unsigned DIM>
  class HermiteAnisotropicPVDElement : public virtual HermitePVDElement<DIM>, 
                                        public virtual AnisotropicPVDEquations<DIM>
  {
  public:
    /// Constructor, there are no internal data points
    HermiteAnisotropicPVDElement() : HermitePVDElement<DIM>(), 
    AnisotropicPVDEquations<DIM>() { }

    /// SolidQHermiteElement output function
    void output(std::ostream &outfile)
    {HermitePVDElement<DIM>::output(outfile);}

    /// SolidQHermiteElement output function
    void output(std::ostream &outfile, const unsigned &n_plot)
    {HermitePVDElement<DIM>::output(outfile,n_plot);}

    /// C-style SolidQHermiteElement output function
    void output(FILE* file_pt) {HermitePVDElement<DIM>::output(file_pt);}

    /// C-style SolidQHermiteElement output function
    void output(FILE* file_pt, const unsigned &n_plot)
    {HermitePVDElement<DIM>::output(file_pt,n_plot);}
  };


//No need to override ProjectablePVDElement

  template <unsigned DIM>
  class AnisotropicPVDEquationsWithPressure : public virtual AnisotropicPVDEquationsBase<DIM>,
                                              public virtual PVDEquationsWithPressure<DIM>,
                                              public virtual SolidElementWithDiagonalMassMatrix
  {
  public:
    AnisotropicPVDEquationsWithPressure() : AnisotropicPVDEquationsBase<DIM>(),
                                            PVDEquationsWithPressure<DIM>()
    {}

    //Need to override this one to get preferential vectors and active strain
    // so they can be passed to the anisotropic constitutive law
    void get_stress(const Vector<double> &s, DenseMatrix<double> &sigma);


    //Const version for output
    double interpolated_solid_p(const Vector<double> &s) const 
    {
      //Find number of nodes
      unsigned n_solid_pres = this->npres_solid();
      //Local shape function
      Shape psisp(n_solid_pres);
      //Find values of shape function
      this->solid_pshape(s,psisp);

      //Initialise value of solid_p
      double interpolated_solid_p = 0.0;
      //Loop over the local nodes and sum
      for(unsigned l=0;l<n_solid_pres;l++) 
      {
        interpolated_solid_p += this->solid_p(l)*psisp[l];
      }

      return(interpolated_solid_p);
    }

    //const version
    virtual double solid_p(const unsigned &l) const = 0;

    unsigned nscalar_paraview() const
    {
      return (DIM + 1 + (DIM*DIM) + 1);
    }

    void scalar_value_paraview(std::ofstream& file_out,
                              const unsigned& i,
                              const unsigned& nplot) const
    {
      Vector<double> x(DIM);
      Vector<double> xi(DIM);
      Vector<double> s(DIM);
      DenseMatrix<double> stress_or_strain(DIM,DIM);

      // Loop over plot points
      unsigned num_plot_points=this->nplot_points_paraview(nplot);
      for (unsigned iplot=0;iplot<num_plot_points;iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot,nplot,s);

        // Get Eulerian and Lagrangian coordinates
        this->interpolated_x(s,x);
        this->interpolated_xi(s,xi);

        //xi
        if(i<DIM)
        {
          file_out << xi[i] << std::endl;
          continue;
        }
        //isotropic growth
        if(i<DIM+1)
        {
          // Get isotropic growth
          double gamma;
          // Dummy integration point
          unsigned ipt=0;
          this->get_isotropic_growth(ipt,s,xi,gamma);
          file_out << gamma << std::endl;
          continue;
        }
        //strain
        if(i<(DIM+1 + (DIM*DIM)))
        {
          this->get_strain(s,stress_or_strain);
          file_out << stress_or_strain((i-DIM-1)%DIM,unsigned((i-DIM-1)/DIM)) << std::endl;
          continue;
        }
        if(i<(DIM + 1 + (DIM*DIM)) + 1)
        {
          file_out << this->interpolated_solid_p(s) << std::endl;
          continue;
        }
        //Can't include stress because get_stress is not a const function
        // //stress
        // if(i<(2*DIM+1 + (DIM*DIM)))
        // {
        //   this->get_stress(s,stress_or_strain);
        //   file_out << stress_or_strain((i-(DIM+1 + (DIM*DIM)))%DIM,unsigned((i-(DIM+1 + (DIM*DIM)))/DIM)) << std::endl;
        //   return;
        // }

        //Never get here
        #ifdef PARANOID
        std::stringstream error_stream;
        error_stream << "paraview anisotropic solid output died" << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        #endif
      }
    }

    void scalar_value_fct_paraview(std::ofstream& file_out,
    const unsigned& i,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt
    exact_soln_pt) const
    {
      scalar_value_paraview(file_out,i,nplot);
    }

    std::string scalar_name_paraview(const unsigned& i) const
    {
      //xi
      if(i<DIM)
      {
        return ("xi"+std::to_string(i));
      }
      //isotropic growth
      if(i<DIM+1)
      {
        return "Isotropic growth";
      }
      //strain
      if(i<(DIM+1+(DIM*DIM)))
      {
        return ("Strain"+std::to_string(i-(DIM+1)));
      }

      if(i<(DIM + 1 + (DIM*DIM)) + 1)
      {
        return "Pressure";
      }

      // if(i<(DIM+1+2*(DIM*DIM)))
      // {
      //   return ("Stress"+std::to_string(i-(DIM+1 + (DIM*DIM))));
      // }

      throw OomphLibError("Requested variable index is too large",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    }



  protected:
    //We need to override this since get_stress(g,G,sigma) does not have positional
    // information as an argument so cannot be overridden itself to get the
    // preferential vectors or active strain. Overriding but only 4 lines are
    // changed in any significant way
    void fill_in_generic_residual_contribution_pvd_with_pressure(Vector<double> &residuals, DenseMatrix<double> &jacobian,
                                                                        DenseMatrix<double> &mass_matrix,
                                                                        const unsigned& flag) override;

    //We have to reimplement get_stress(g,G,sigma (,etc)) to also accept the
    // preferential vectors and active strainss
    inline void get_stress(const DenseMatrix<double> &g, 
                          const DenseMatrix<double> &G,
                          const Vector<Vector<double>>& A,
                          const Vector<double> &V,
                          DenseMatrix<double> &sigma_dev, 
                          DenseMatrix<double> &Gcontra, 
                          double &gen_dil, double &inv_kappa) 
    {
    #ifdef PARANOID
      //If the pointer to the constitutive law hasn't been set, issue an error
      if(this->Anisotropic_constitutive_law_pt == 0)
      {
        //Write an error message
        std::string error_message =
        "Elements derived from AnisotropicPVDEquationsWithPressure \n";
        error_message += "must have a constitutive law:\n";
        error_message +=
        "set one using the anisotropic_constitutive_law_pt() member function";
        //Throw the error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    #endif
      this->Anisotropic_constitutive_law_pt->
      calculate_second_piola_kirchhoff_stress(g,G,A,V,sigma_dev,Gcontra,gen_dil,inv_kappa);
    }

    /// \short Return the derivative of the 
    /// deviatoric part of the 2nd Piola Kirchhoff stress 
    /// tensor, as calculated from the constitutive law in the nearly 
    /// incompresible formulation. Also return the derivative of the 
    /// generalised dilatation.
    inline void get_d_stress_dG_upper(const DenseMatrix<double> &g, 
                                      const DenseMatrix<double> &G,
                                      const Vector<Vector<double>>& A,
                                      const Vector<double> &V,
                                      const DenseMatrix<double> &sigma,
                                      const double &gen_dil,                    
                                      const double &inv_kappa,           
                                      const double &interpolated_solid_p,
                                      RankFourTensor<double> &d_sigma_dG,
                                      DenseMatrix<double> &d_gen_dil_dG)
    {
    #ifdef PARANOID
      //If the pointer to the constitutive law hasn't been set, issue an error
      if(this->Anisotropic_constitutive_law_pt == 0)
      {
        //Write an error message
        std::string error_message =
        "Elements derived from AnisotropicPVDEquationsWithPressure \n";
        error_message += "must have a constitutive law:\n";
        error_message +=
        "set one using the anisotropic_constitutive_law_pt() member function";
        //Throw the error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    #endif
      //Only bother with the symmetric part by passing false as last entry
      this->Anisotropic_constitutive_law_pt->
      calculate_d_second_piola_kirchhoff_stress_dG(g,G,A,V,sigma,gen_dil,inv_kappa,interpolated_solid_p,d_sigma_dG,d_gen_dil_dG,false);
    }

    inline void get_stress(const DenseMatrix<double> &g, 
                          const DenseMatrix<double> &G,
                          const Vector<Vector<double>>& A,
                          const Vector<double> &V,
                          DenseMatrix<double> &sigma_dev, 
                          DenseMatrix<double> &Gcontra, 
                          double &detG)
    {
    #ifdef PARANOID
      //If the pointer to the constitutive law hasn't been set, issue an error
      if(this->Anisotropic_constitutive_law_pt == 0)
      {
        //Write an error message
        std::string error_message =
        "Elements derived from AnisotropicPVDEquationsWithPressure \n";
        error_message += "must have a constitutive law:\n";
        error_message +=
        "set one using the anisotropic_constitutive_law_pt() member function";
        //Throw the error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    #endif
      this->Anisotropic_constitutive_law_pt->calculate_second_piola_kirchhoff_stress(g,G,A,V,sigma_dev,Gcontra,detG);
    }

    /// \short  Return the derivative of the 2nd Piola Kirchhoff stress 
    /// tensor, as calculated from the constitutive law in the 
    /// incompresible formulation. Also return 
    /// derivative of the determinant of the deformed covariant metric tensor 
    /// (likely to be needed in the incompressibility constraint)
    inline void get_d_stress_dG_upper(const DenseMatrix<double> &g, 
                                      const DenseMatrix<double> &G,
                                      const Vector<Vector<double>>& A,
                                      const Vector<double> &V,
                                      const DenseMatrix<double> &sigma,
                                      const double &detG,                    
                                      const double &interpolated_solid_p,
                                      RankFourTensor<double> &d_sigma_dG,
                                      DenseMatrix<double> &d_detG_dG)
    {
    #ifdef PARANOID
      //If the pointer to the constitutive law hasn't been set, issue an error
      if(this->Anisotropic_constitutive_law_pt == 0)
      {
        //Write an error message
        std::string error_message =
        "Elements derived from AnisotropicPVDEquationsWithPressure \n";
        error_message += "must have a constitutive law:\n";
        error_message +=
        "set one using the anisotropic_constitutive_law_pt() member function";
        //Throw the error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    #endif
      //Only bother with the symmetric part by passing false as last entry
      this->Anisotropic_constitutive_law_pt->calculate_d_second_piola_kirchhoff_stress_dG(g,G,A,V,sigma,detG,
                                                                                          interpolated_solid_p,
                                                                                          d_sigma_dG,d_detG_dG,false);
    }
  };

  template<unsigned DIM>
  class QAnisotropicPVDElementWithPressure : public virtual QPVDElementWithPressure<DIM>,
                                            public virtual AnisotropicPVDEquationsWithPressure<DIM>
  {
  public:
    QAnisotropicPVDElementWithPressure() : QPVDElementWithPressure<DIM>(),
                                            AnisotropicPVDEquationsWithPressure<DIM>()
    {}

    double solid_p(const unsigned &l) const
    {
      return this->internal_data_pt(this->P_solid_internal_index)->value(l);
    }
 
  };

  //======================================================================
  /// FaceGeometry of 2D QAnisotropicPVDElementWithPressure
  //======================================================================
  template<>
  class FaceGeometry<QAnisotropicPVDElementWithPressure<2> >: 
  public virtual SolidQElement<1,3>
  {
  public:
    /// Constructor must call constructor of underlying solid element
    FaceGeometry() : SolidQElement<1,3>() {}
  };


  //======================================================================
  /// FaceGeometry of FaceGeometry of 2D QAnisotropicPVDElementWithPressure
  //======================================================================
  template<>
  class FaceGeometry<FaceGeometry<QAnisotropicPVDElementWithPressure<2> > >: 
  public virtual PointElement
  {
  public:
    /// Constructor must call constructor of underlying solid element
    FaceGeometry() : PointElement() {}
  };

  //======================================================================
  /// FaceGeometry of 3D QAnisotropicPVDElementWithPressure
  //======================================================================
  template<>
  class FaceGeometry<QAnisotropicPVDElementWithPressure<3> >: 
  public virtual SolidQElement<2,3>
  {
  public:
    /// Constructor must call constructor of underlying solid element
    FaceGeometry() : SolidQElement<2,3>() {}
  };


  //======================================================================
  /// FaceGeometry of FaceGeometry of 3D QAnisotropicPVDElementWithPressure
  //======================================================================
  template<>
  class FaceGeometry<FaceGeometry<QAnisotropicPVDElementWithPressure<3> > >: 
  public virtual SolidQElement<1,3>
  {
  public:
    /// Constructor must call constructor of underlying solid element
    FaceGeometry() : SolidQElement<1,3>() {}
  };

  template<unsigned DIM>
  class QAnisotropicPVDElementWithContinuousPressure : public virtual QPVDElementWithContinuousPressure<DIM>,
                                                      public virtual AnisotropicPVDEquationsWithPressure<DIM>
  {
  public:
    QAnisotropicPVDElementWithContinuousPressure() : QPVDElementWithContinuousPressure<DIM>(),
                                                    AnisotropicPVDEquationsWithPressure<DIM>()
    {}

    double solid_p(const unsigned &l) const 
    {
      return this->nodal_value(this->Pconv[l],this->solid_p_nodal_index());
    }

  };

  //===============================================================
  /// FaceGeometry for 2D QAnisotropicPVDElementWithContinuousPressure element
  //===============================================================
  template<>
  class FaceGeometry<QAnisotropicPVDElementWithContinuousPressure<2> >: 
  public virtual SolidQElement<1,3>
  {
  public:
    /// Constructor must call constructor of the underlying Solid element
    FaceGeometry() : SolidQElement<1,3>() {}
  };



  //===============================================================
  /// FaceGeometry of FaceGeometry 
  /// for 2D QAnisotropicPVDElementWithContinuousPressure element
  //===============================================================
  template<>
  class FaceGeometry<FaceGeometry<QAnisotropicPVDElementWithContinuousPressure<2> > >: 
  public virtual PointElement
  {
  public:
    /// Constructor must call constructor of the underlying Point element
    FaceGeometry() : PointElement() {}
  };


  //===============================================================
  /// FaceGeometry for 3D QAnisotropicPVDElementWithContinuousPressure element
  //===============================================================
  template<>
  class FaceGeometry<QAnisotropicPVDElementWithContinuousPressure<3> >: 
  public virtual SolidQElement<2,3>
  {
  public:
    /// Constructor must call constructor of the underlying Solid element
    FaceGeometry() : SolidQElement<2,3>() {}
  };



  //===============================================================
  /// FaceGeometry of FaceGeometry 
  /// for 3D QAnisotropicPVDElementWithContinuousPressure element
  //===============================================================
  template<>
  class FaceGeometry<FaceGeometry<QAnisotropicPVDElementWithContinuousPressure<3> > >: 
  public virtual SolidQElement<1,3>
  {
  public:
    /// Constructor must call constructor of the underlying element
    FaceGeometry() : SolidQElement<1,3>() {}
  };


  template<unsigned DIM, unsigned NNODE_1D>
  class TAnisotropicPVDElement : public virtual TPVDElement<DIM, NNODE_1D>,
                                  public virtual AnisotropicPVDEquations<DIM>
  {
  public:

    //Final overrider
    void get_stress(const Vector<double> &s, DenseMatrix<double> &sigma)
    {
      AnisotropicPVDEquations<DIM>::get_stress(s, sigma);
    }

    TAnisotropicPVDElement() : TPVDElement<DIM, NNODE_1D>(), AnisotropicPVDEquations<DIM>()
    {}
  };

  //============================================================================
  /// FaceGeometry of a 2D TAnisotropicPVDElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TAnisotropicPVDElement<2,NNODE_1D> > :
  public virtual SolidTElement<1,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidTElement<1,NNODE_1D>() {}
  };


  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 2D TAnisotropicPVDElement 
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<TAnisotropicPVDElement<2,NNODE_1D> > >:
  public virtual PointElement
  {
  public:
    //Make sure that we call the constructor of the SolidQElement
    //Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };


  //============================================================================
  /// FaceGeometry of a 3D TAnisotropicPVDElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TAnisotropicPVDElement<3,NNODE_1D> > :
  public virtual SolidTElement<2,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidTElement<2,NNODE_1D>() {}
  };

  //============================================================================
  /// FaceGeometry of FaceGeometry of a 3D TAnisotropicPVDElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<TAnisotropicPVDElement<3,NNODE_1D> > > :
  public virtual SolidTElement<1,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidTElement<1,NNODE_1D>() {}
  };


  template<unsigned DIM, unsigned NNODE_1D>
  class TAnisotropicPVDBubbleEnrichedElement : TPVDBubbleEnrichedElement<DIM, NNODE_1D>,
                                                AnisotropicPVDEquations<DIM>
  {
  public:
    TAnisotropicPVDBubbleEnrichedElement() : TPVDBubbleEnrichedElement<DIM, NNODE_1D>(),
                                              AnisotropicPVDEquations<DIM>()
    {}
  };

  //============================================================================
  /// FaceGeometry of a 2D TAnisotropicPVDBubbleEnrichedElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TAnisotropicPVDBubbleEnrichedElement<2,NNODE_1D> > :
  public virtual SolidTElement<1,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidTElement<1,NNODE_1D>() {}
  };


  //==============================================================
  /// FaceGeometry of the FaceGeometry of the 2D TAnisotropicPVDBubbleEnrichedElement 
  //==============================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<TAnisotropicPVDBubbleEnrichedElement<2,NNODE_1D> > >:
  public virtual PointElement
  {
  public:
    //Make sure that we call the constructor of the SolidQElement
    //Only the Intel compiler seems to need this!
    FaceGeometry() : PointElement() {}
  };


  //============================================================================
  /// FaceGeometry of a 3D TAnisotropicPVDBubbleEnrichedElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<TAnisotropicPVDBubbleEnrichedElement<3,NNODE_1D> > :
  public virtual SolidTBubbleEnrichedElement<2,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidTBubbleEnrichedElement<2,NNODE_1D>() {}
  };

  //============================================================================
  /// FaceGeometry of FaceGeometry of a 3D TAnisotropicPVDElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FaceGeometry<TAnisotropicPVDBubbleEnrichedElement<3,NNODE_1D> > > :
  public virtual SolidTElement<1,NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidTElement<1,NNODE_1D>() {}
  };


  template <unsigned DIM>
  class TAnisotropicPVDElementWithContinuousPressure : public virtual TPVDElementWithContinuousPressure<DIM>,
                                                        public virtual AnisotropicPVDEquationsWithPressure<DIM>
  {
    TAnisotropicPVDElementWithContinuousPressure() : TPVDElementWithContinuousPressure<DIM>(),
                                                      AnisotropicPVDEquationsWithPressure<DIM>()
    {}

    /// Broken copy constructor
    TAnisotropicPVDElementWithContinuousPressure(const TAnisotropicPVDElementWithContinuousPressure<DIM>& dummy) 
    { 
      BrokenCopy::broken_copy("TAnisotropicPVDElementWithContinuousPressure");
    }
  };

  //=======================================================================
  /// Face geometry of the 2D Taylor_Hood elements
  //=======================================================================
  template<>
  class FaceGeometry<TAnisotropicPVDElementWithContinuousPressure<2> >:
  public virtual SolidTElement<1,3>
  {
  public:
    /// Constructor: Call constructor of base
    FaceGeometry() : SolidTElement<1,3>() {}
  };

  //=======================================================================
  /// Face geometry of the 3D Taylor_Hood elements
  //=======================================================================
  template<>
  class FaceGeometry<TAnisotropicPVDElementWithContinuousPressure<3> >: 
  public virtual SolidTElement<2,3>
  {
  public:
    /// Constructor: Call constructor of base
    FaceGeometry() : SolidTElement<2,3>() {}
  };










} //end namespace

#endif