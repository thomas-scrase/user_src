//Bidomain Equations
// template <unsigned DIM>
// class BidomainEquations :
// public virtual BaseCellMembranePotentialEquations<DIM>
// {
// private:
//   /// \short Static int to hold number of variables at nodes
//   /// transmembrane potential and extracellular potential
//   static const unsigned Initial_Nvalue = 2;

//   /// Add Phi_e_index_Bidomain

//   /// Add interpolate extracellular potential

//   /// Add interpolate d extracellular potential / dt

//   /// Add calculate intracellular potential

//   //Overload the residual for the bidomain equations
//   void fill_in_generic_residual_contribution_BaseCellMembranePotential(
//     Vector<double> &residuals, DenseMatrix<double> &jacobian, 
//     DenseMatrix<double> &mass_matrix, unsigned flag)
//   {
//   throw OomphLibError(
//     "fill_in_generic_residual_contribution_BaseCellMembranePotential not implemented for bidomain equations yet",
//     OOMPH_CURRENT_FUNCTION,
//     OOMPH_EXCEPTION_LOCATION);
//  }

// };