//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef OOMPH_PARTITIONED_CELL_SOLVER
#define OOMPH_PARTITIONED_CELL_SOLVER


#include "../generic/problem.h"
#include "../generic/geom_objects.h"
#include "../generic/mesh.h"
//We use the objects defined in partitioned fsi solver
// For handling the picard convergence
#include "../multi_physics/partitioned_fsi_solver.h"

#include "../conducting_cell/conducting_cell_elements.h"

namespace oomph
{

  enum PartitionedSolverTypes{
    StrangSplitting,
    Whiteley2006,
    Whiteley2006WithCrankNicolson,
    ExplicitCellImplicitDiff,
    IterativeStrangSplitting,
    CrankNicolsonCrankNicolson
  };

//===============================================================
/// Base class for problems that can be solved by partitioned 
/// FSI solver
//===============================================================
  template<class ELEMENT_TYPE>
 class SegregatableCellProblem : public virtual Problem
  {
    protected:

   /// \short This function is called once at the start of each
   /// partitioned solve.
   virtual void actions_before_partitioned_solve() {}

   /// \short This function is called once at the end of each
   /// partitioned solve.
   virtual void actions_after_partitioned_solve() {}

   /// \short This function is to be filled with actions that take place
   /// before the check for convergence of the entire partitioned solve
   virtual void actions_before_partitioned_convergence_check() {}

    public:

   /// \short Constructor. Set default values for solver parameters:
   /// - Don't use pointwise Aitken extrapolation but if it's used at 
   ///   all, it's used immediately.
   /// - No under-relaxation at all (neither classical nor Irons&Tuck)
   /// - Base convergence on max. residual of coupled system of eqns
   /// - Convergence tolerance = tolerance for Newton solver
   ///   defined in Problem base class
   /// - Max. 50 Picard iterations
   SegregatableCellProblem() 
    {
     // Use pointwise Aitken extrapolation?
     Use_pointwise_aitken=false;

     // Default: No under-relaxation
     Omega_relax=1.0;
 
     // Don't use of Irons and Tuck's extrapolation for vm values
     Use_irons_and_tuck_extrapolation=false;

     // Start using pointwise Aitken immediately
     Pointwise_aitken_start=0;

     // By default we don't recheck convergence
     Recheck_convergence_after_pointwise_aitken=false;

     // Convergence criterion
     Convergence_criterion=Assess_convergence_based_on_max_global_residual;

     // Convergence tolerance (as in global Newton solver)
     Convergence_tolerance=Problem::Newton_solver_tolerance;

     // Doc max. global residual during iteration?
     Doc_max_global_residual=false;

     // Max. number of Picard iterations
     Max_picard=50;

     // Pointer to Mesh containing only cell elements -- the elements in this
     // Mesh will be excluded from the assembly process when
     // the vm problem is solved
     Cell_mesh_pt=0;


     // Initialise timer that allows doc of iteration/cpu time
     T_ref = clock();
     T_spent_on_actual_solve=0.0;

     /// \short boolean flag to indicate if timer has been halted
     Timer_has_been_halted=false;

     //By default we use whiteley 2006 which solves for vm(tnp1) with Iion(tn)
     // Then solves for cell vars(tnp1) with vm(tnp1)
     PartitionedSolveType = Whiteley2006;

     // enable_under_relaxation();
     // enable_pointwise_aitken();
    }

   /// Empty virtual destructor
   virtual ~SegregatableCellProblem(){}

   void set_partitioned_solve_type(const unsigned &type){
    switch(type){
      case StrangSplitting :
        PartitionedSolveType = type;
        set_cells_interpolate_flags(Implicit,Implicit);
        break;
      case Whiteley2006 : 
        PartitionedSolveType = type;
        set_cells_interpolate_flags(Explicit, Implicit);
        break;
      case Whiteley2006WithCrankNicolson : 
        PartitionedSolveType = type;
        set_cells_interpolate_flags(CrankNicolson,CrankNicolson);
        break;
      case ExplicitCellImplicitDiff :
        PartitionedSolveType = type;
        set_cells_interpolate_flags(Implicit, Explicit);
        break;
      case IterativeStrangSplitting :
        PartitionedSolveType = type;
        set_cells_interpolate_flags(Implicit,Implicit);
        break;
      case CrankNicolsonCrankNicolson :
        PartitionedSolveType = type;
        set_cells_interpolate_flags(CrankNicolson,CrankNicolson);
        break;

    }
   }

   void set_cells_interpolate_flags(const unsigned &Iion_Flag,
                                    const unsigned &Vm_Flag)
   {
      for(unsigned e=0; e<Cell_mesh_pt->nelement(); e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e));
        elem_pt->set_interpolated_iion_solver_flag(Iion_Flag);
        elem_pt->set_interpolated_vm_solver_flag(Vm_Flag);
      }
   }

   //Calculate derivatives, output data, and Iion from the cell model,
   // this is generally performed after cell variables are solved for
   // and it's results are inserted into the current time values
   void update_cell_model_data()
   {
    Oomph_Cell_Interface_Helpers::Get_And_Update_Cell_Model_Data_Over_Mesh<ELEMENT_TYPE>(this, Cell_mesh_pt);
   }

   void explicit_cell_model_timestep(const double &dt)
   {
    Oomph_Cell_Interface_Helpers::Partitioned_Cell_Timestep_Solve_Over_Mesh<ELEMENT_TYPE>(dt, this, Cell_mesh_pt, false);
   }

   void explicit_cell_model_timestep(const double &dt, const bool &use_current_as_initial)
   {
    Oomph_Cell_Interface_Helpers::Partitioned_Cell_Timestep_Solve_Over_Mesh<ELEMENT_TYPE>(dt, this, Cell_mesh_pt, use_current_as_initial);
   }


   /// \short Setup the partitioned solver: Backup the pinned status of
   /// the cell and vm dofs and allocate the internal storage
   /// based on the input provided by identify_cell_and_vm_dofs(...)
   /// In addition, reset storage associated with convergence acceleration
   /// techniques.
   /// If the problem and degrees of freedom has not changed between
   /// calls to the solver then it is wasteful to call 
   /// identify_cell_and_vm_dofs(...) again and again. If the optional
   /// boolean flag is set to false then the storage for convergence
   /// acceleration techniques is reset, but the cell and vm dofs
   /// are not altered.
   void setup_partitioned_solver()
   {      
    if (Cell_mesh_pt==0)
    {
      oomph_info 
       << std::endl << std::endl 
       << "Warning: Your implementation of the pure virtual\n"
       << "         function identify_cell_and_vm_dofs(...)\n"
       << "         returned a NULL pointer for Cell_mesh_pt.\n"
       << "         --> The cell elements will remain activated\n"
       << "         during the vm solve. This is inefficient!\n"
       << "         You should combine all cell elements into a combined\n"
       << "         mesh and specify this mesh in your\n"
       << "         implementation of \n\n"
       << "         SegregatableCellProblem<ELEMENT_TYPE>::identify_cell_and_vm_dofs(...)"
       << std::endl << std::endl;
    }

    //Pin all the nodal variables we use simply for storage. We will unpin them
    // for fill in with unpin_all_dofs
    //This is okay since we only unpin these when we call on the cell model
    // to calculate the values, under which circumstances it is okay if all
    // the other dofs are unpinned sincen we are not assembling any residuals.
    for(unsigned e=0; e<Cell_mesh_pt->nelement(); e++){
      ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e));
      elem_pt->pin_all_cell_model_output_vars();
      elem_pt->pin_all_cell_model_vars_derivatives();
      elem_pt->pin_cell_model_iion();
    }

        
    //The total number of variables
    unsigned n_values = 0;
    
    // Back up the pointers to the submeshes in the original problem
    // so we can restore the problem when we're done
    unsigned orig_n_sub_mesh=nsub_mesh();

    Value_is_pinned.resize(orig_n_sub_mesh);

    Orig_sub_mesh_pt.resize(orig_n_sub_mesh);

    for (unsigned i=0;i<orig_n_sub_mesh;i++)
    {
      Orig_sub_mesh_pt[i]=mesh_pt(i);

      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();

      Value_is_pinned[i].resize(n_elem);

      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        Value_is_pinned[i][e].resize(n_node);
        for(unsigned l=0; l<n_node; l++){
          Value_is_pinned[i][e][l].resize(elem_pt->required_nvalue(l));
          //Loop over cell variables and fill
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            Value_is_pinned[i][e][l][k]=elem_pt->node_pt(l)->is_pinned(k);
            n_values++;
          }
        }
      }
    }


    //Resize storage for the previous values
    Previous_value.resize(n_values);

    // Allocate storage and initialise Irons and Tuck extrapolation
    if (Use_irons_and_tuck_extrapolation)
    {
      Del_irons_and_tuck.resize(n_values);
    }
    
    // Make space for pointwise Aitken extrapolation
    if (Use_pointwise_aitken)
    {
      Pointwise_aitken_value.resize(n_values);
      for (unsigned i=0;i<n_values;i++)
      {
        Pointwise_aitken_value[i].resize(3);
      }
    }

    // Initialise Irons and Tuck relaxation factor
    R_irons_and_tuck=1.0-Omega_relax;

    //Initialise Irons and Tuck delta values
    unsigned n_del = Del_irons_and_tuck.size();
    for (unsigned i=0;i<n_del;i++) {Del_irons_and_tuck[i]=1.0e20;}

    // Initialise counter for the number of pointwise Aitken values stored
    Pointwise_aitken_counter=0;

    //Set the partitioned solve type to the partitioned solve type
    // so that element get value flags are setup correctly
    set_partitioned_solve_type(PartitionedSolveType);
  }
   
   /// \short Partitioned solver. Peform a partitioned step from
   /// the present state of the system.
   /// Returns PicardConvergenceData object that contains the vital
   /// stats of the iteration
   PicardConvergenceData partitioned_solve()
    {
    // Initialise timer for essential bits of code
    reset_timer();
    // Start timer for overall solve
    clock_t t_total_start = clock();

    // If convergence is based on max. global residual we may as well
    // document it...
    bool get_max_global_residual=Doc_max_global_residual;
    if(Convergence_criterion==Assess_convergence_based_on_max_global_residual)
     {
      get_max_global_residual=true;
     }

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;
     
    // Initialise total time for computation of global residuals
    double cpu_for_global_residual=0.0;
     
    //Update anything that needs updating
    actions_before_partitioned_solve();
     
    // Set flags to values that are appropriate if Picard iteration
    // does not converge with Max_picard iterations
    bool converged=false;
    unsigned iter_taken=Max_picard;

    // This one will be overwritten during the convergence checks
    double tol_achieved=0.0;
           
    // Store the current values of the vm dofs as reference values
    // and for the pointwise Aitken extrapolation 
    store_dofs();

    // Loop over Picard iterations
    //----------------------------
    for (unsigned iter=1;iter<=Max_picard;iter++)
     {

      //Calculate the initial residuals
      if(iter==1)
       {
        //Problem is always non-linear?
        //Perform any actions before the convergence check
        actions_before_partitioned_convergence_check();

        //Update the cell model data, i.e. Iion, and active strain etc.
        unpin_all_dofs();
        assign_eqn_numbers();
        update_cell_model_data();
         
        double max_res = 0.0;
        if (get_max_global_residual)
         {
          clock_t t_start = clock(); 
          // restore_cell_variable_dofs();
          //We only need the vm and diffusion forcing term dofs to be unpinned
          restore_dofs();
          rebuild_monolithic_mesh();
          assign_eqn_numbers();
          DoubleVector residual;
          get_residuals(residual);
          //Get maximum residuals using our own abscmp function
          max_res = residual.max();
          clock_t t_end = clock();
          cpu_for_global_residual+=double(t_end-t_start)/CLOCKS_PER_SEC;
         }
         
        oomph_info << "==================================================\n";
        oomph_info <<   "Initial iteration     : " 
                   << 0 << std::endl;
        oomph_info <<   "RMS  change           : "
                   << 0 << std::endl;
        oomph_info <<   "Max. change           : "
                   << 0 << std::endl;
        oomph_info <<   "RMS norm              : "
                   << 0   << std::endl;
        if (get_max_global_residual)
         {
          oomph_info << "Max. global residual  : "
                     << max_res   << std::endl;
         }
        oomph_info << "==================================================\n\n";
         
        // Check for convergence, but this only makes sense
        // for the maximum (rather than relative case)
        if((Convergence_criterion==
            Assess_convergence_based_on_max_global_residual) &&
           (max_res < Convergence_tolerance))
         {
          oomph_info
           << "\n\n\n////////////////////////////////////////////////////////"
           << "\nPicard iteration converged after " 
           << 0 << " steps!" << std::endl
           << "Convergence was based on max. residual of coupled eqns \n"
           << "being less than " << Convergence_tolerance << "." << std::endl
           << "////////////////////////////////////////////////////////\n\n\n"
           << std::endl;
           converged = true;

          //Converged, so bail out
          // Number of iterations taken
          iter_taken=0;
             
          // We have converged (overwrites default of false)
          conv_data.set_solver_converged();
             
          // Break loop using a GO TO! This is the first (and hopefully 
          // the last) one in the whole of oomph-lib. Here it's 
          // it's actually the cleanest way to exit from these
          // nested control structures
          goto jump_out_of_picard;
         }
       }
       


      // Stage 1: Perform a decoupled solve of the single cell variables
      //------------------------------------------------------------
      // and re-assign the equation numbers
      oomph_info <<"\n\nDOING OOMPH-LIB CELL SOLVE\n\n";
      pin_all_dofs();
      unpin_cell_variable_dofs();
      assign_eqn_numbers();
      newton_solve();


      //Update the cell model data, i.e. Iion, and active strain etc.
      unpin_all_dofs();
      assign_eqn_numbers();
      update_cell_model_data();


      // Stage 2: Pin cell variables and predictive membrane potential and solve for Vm using oomph lib
      //------------------------------------------------------
      // Solve the vm problem for the wall solution
      oomph_info <<"\n\nDOING OOMPH-LIB VM SOLVE\n\n";
      restore_dofs();
      pin_cell_variable_dofs();
      assign_eqn_numbers();
      newton_solve();


      // if(PartitionedSolveType==Whiteley2006){goto jump_out_of_picard;}



      // Under-relax, either by classical relaxation or by Irons & Tuck
      // Note that we are under-relaxing before the convergence check
      // because under-relaxtion may be required to acheieve any
      // kind of convergence. If the convergence check is on the RELATIVE
      // change of the residuals, however, then a small under-relaxation
      // parameter will cause a false convergence. You have been warned!
      under_relax_data();
       
      // Stage 3: Convergence check (possibly again after pointwise Aitken 
      //------------------------------------------------------------------
      // extrapolation)
      //---------------
      //Assume that we are not doing Aitken acceleration
      Recheck_convergence_after_pointwise_aitken=false;
      do
       {
        //Perform any actions before the convergence check
        actions_before_partitioned_convergence_check();
         
        //Get the change in the vm variables
        double rms_change;
        double max_change;
        double rms_norm;
        double max_res=0.0;
        get_change(rms_change,max_change,rms_norm);
         
         
        //If we are computing the maximum global residual, do so
        if (get_max_global_residual)
         {
          clock_t t_start = clock();

          // restore_cell_variable_dofs();
          restore_dofs();
          rebuild_monolithic_mesh();
          assign_eqn_numbers();

          //Get the residuals
          DoubleVector residual;
          get_residuals(residual);
           
          //Get maximum residuals, using our own abscmp function
          max_res =  residual.max();
           
          clock_t t_end = clock();
          cpu_for_global_residual+=double(t_end-t_start)/CLOCKS_PER_SEC;
         }
         
        oomph_info << "==================================================\n";
        oomph_info <<   "Iteration             : " 
                   << iter << std::endl;
        oomph_info <<   "RMS  change           : "
                   << rms_change << std::endl;
        oomph_info <<   "Max. change           : "
                   << max_change << std::endl;
        oomph_info <<   "RMS norm              : "
                   << rms_norm   << std::endl;
        if (get_max_global_residual)
         {
          oomph_info << "Max. global residual  : "
                     << max_res   << std::endl;
         }
        oomph_info << "==================================================\n\n";
       
        // Check for convergence
        switch (Convergence_criterion)
         {
         
         case Assess_convergence_based_on_absolute_vm_change:
          tol_achieved=max_change;
          if (tol_achieved<Convergence_tolerance)
           {
            oomph_info 
             << "\n\n\n////////////////////////////////////////////////////////"
             << "\nPicard iteration converged after " 
             << iter << " steps!" << std::endl
             << "Convergence was based on absolute change in vm dofs \n"
             << "being less than " << Convergence_tolerance << "." << std::endl
             << "////////////////////////////////////////////////////////\n\n\n"
             << std::endl;
            converged=true;
           }
          break;
           
           
         case Assess_convergence_based_on_relative_vm_change:
          tol_achieved=std::fabs(rms_change/rms_norm);
          if (tol_achieved<Convergence_tolerance)
           {
            oomph_info
             << "\n\n\n///////////////////////////////////////////////////////" 
             << "\nPicard iteration converged after " 
             << iter << " steps!" << std::endl
             << "Convergence was based on relative change in vm dofs \n"
             << "being less than " << Convergence_tolerance << "." << std::endl
             << "////////////////////////////////////////////////////////\n\n\n"
             << std::endl;
            converged=true;
           }
          break;
           
         case Assess_convergence_based_on_max_global_residual:
          tol_achieved=max_res;
          if (tol_achieved<Convergence_tolerance)
           {
            oomph_info
             << "\n\n\n////////////////////////////////////////////////////////"
             << "\nPicard iteration converged after " 
             << iter << " steps!" << std::endl
             << "Convergence was based on max. residual of coupled eqns \n"
             << "being less than " << Convergence_tolerance << "." << std::endl
             << "////////////////////////////////////////////////////////\n\n\n"
             << std::endl;
            converged=true;
           }
          break;

         }

        // If converged bail out
        if (converged)
         {
          // Number of iterations taken
          iter_taken=iter;

          // We have converged (overwrites default of false)
          conv_data.set_solver_converged();

          // Break loop using a GO TO! This is the first (and hopefully 
          // the last) one in the whole of oomph-lib. Here it's 
          // it's actually the cleanest way to exit from these
          // nested control structures
          goto jump_out_of_picard;
         }
        
        // Store the current values of the vm dofs as reference values
        // and for the pointwise Aitken extrapolation 
        store_dofs();
       
        // Correct wall shape by pointwise Aitken extrapolation if possible
        //-----------------------------------------------------------------
        // and desired
        //------------
        //This is an acceleration method for the (possibly under-relaxed)
        //sequence of iterates.
        if ((3==Pointwise_aitken_counter)&&(Use_pointwise_aitken)&&
            (iter>Pointwise_aitken_start))
         {
          pointwise_aitken_extrapolate();
           
          // Repeat the convergence check
          Recheck_convergence_after_pointwise_aitken=true;
         }
        else
         {
          // Didn't change anything: Don't repeat the convergence check
          Recheck_convergence_after_pointwise_aitken=false;
         } 
       }
      //Repeat convergence while we are doing aitken extrapolation
      while(Recheck_convergence_after_pointwise_aitken);
     
     } // End of loop over iterations
     
     

    // Jump address for converged or diverged iteration
     jump_out_of_picard:

    // Reset everything
    restore_dofs();
    rebuild_monolithic_mesh();
    assign_eqn_numbers();

    // Do any updates that are required 
    actions_after_partitioned_solve();
     
    // Number of iterations (either this is still Max_iter from 
    // the initialisation or it's been overwritten on convergence)
    conv_data.niter()=iter_taken;
     
    // Total cpu time
    clock_t t_total_end = clock();
    conv_data.cpu_total()=double(t_total_end-t_total_start)/
     CLOCKS_PER_SEC;
     
    // Total essential cpu time (exluding doc etc)
    conv_data.essential_cpu_total()=t_spent_on_actual_solve();

    // cpu time for check/doc of global residual
    conv_data.cpu_for_global_residual()=cpu_for_global_residual;
     
    // Final tolerance achieved by the iteration
    conv_data.tol_achieved()=tol_achieved;

    // Doc non-convergence
    if (!converged)
     {
      switch (Convergence_criterion)
       {
         
       case Assess_convergence_based_on_absolute_vm_change:
        oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nPicard iteration did not converge after " 
         << iter_taken << " steps!" << std::endl
         << "Convergence was based on absolute change in vm dofs \n"
         << "being less than " << Convergence_tolerance << " " << std::endl
         << "but we achieved only " << tol_achieved << "." << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;
        //Throw an error indicating if we ran out of iterations
        if (iter_taken==Max_picard)
         {
          throw PartitionedSolverError(true);
         }
        else 
         {
          throw PartitionedSolverError(false);
         }
        break;
         
         
       case Assess_convergence_based_on_relative_vm_change:
        oomph_info
         << "\n\n\n///////////////////////////////////////////////////////" 
         << "\nPicard iteration did not converge after " 
         << iter_taken << " steps!" << std::endl
         << "Convergence was based on relative change in vm dofs \n"
         << "being less than " << Convergence_tolerance << " " << std::endl
         << "but we achieved only " << tol_achieved << "." << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;
        //Throw an error indicating if we ran out of iterations
        if (iter_taken==Max_picard)
         {
          throw PartitionedSolverError(true);
         }
        else 
         {
          throw PartitionedSolverError(false);
         }
        break;
         
       case Assess_convergence_based_on_max_global_residual:
        oomph_info
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nPicard iteration did not converge after " 
         << iter_taken << " steps!" << std::endl
         << "Convergence was based on max. residual of coupled eqns \n"
         << "being less than " << Convergence_tolerance << " " << std::endl
         << "but we achieved only " << tol_achieved << "." << std::endl
          
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;
        //Throw an error indicating if we ran out of iterations
        if (iter_taken==Max_picard)
         {
          throw PartitionedSolverError(true);
         }
        else 
         {
          throw PartitionedSolverError(false);
         }
        break;
         
       }
     }

    return conv_data;
   }
   
   /// \short Steady version of partitioned solver. Makes all
   /// timesteppers steady before solving.
   /// Returns PicardConvergenceData object that contains the
   /// vital stats of the iteration. 
   PicardConvergenceData steady_partitioned_solve()
   {
    //Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();
     
    // Vector of bools to store the is_steady status of the various
    // timesteppers when we came in here
    std::vector<bool> was_steady(n_time_steppers);

    //Loop over them all and make them (temporarily) static
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      was_steady[i]=time_stepper_pt(i)->is_steady();
      time_stepper_pt(i)->make_steady();
     }

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;
     
    //Solve the non-linear problem by the partitioned solver
    try
     {
      conv_data = partitioned_solve();
     }
    //Catch any exceptions thrown in the partitioned solver
    catch(PartitionedSolverError &error)
     {
      if (!error.Ran_out_of_iterations)
       {
        std::ostringstream error_stream;
        error_stream << "Error occured in Partitioned solver. "
                     << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
       }
      else
       {
        oomph_info << "Note: Ran out of iterations but continuing anyway"
                   << std::endl;
       }
     }
     
    // Reset the is_steady status of all timesteppers that
    // weren't already steady when we came in here and reset their
    // weights
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      if (!was_steady[i])
       {
        time_stepper_pt(i)->undo_make_steady();
       }
     }
     
    // Since we performed a steady solve, the history values 
    // now have to be set as if we had performed an impulsive start from
    // the current solution. This ensures that the time-derivatives
    // evaluate to zero even now that the timesteppers have been
    // reactivated. 
    assign_initial_values_impulsive();
     
    //Return the convergence data
    return conv_data;
   }
   

   /// \short Unsteady partitioned solver, advance time by dt and solve
   /// by the partitioned solver. The time values are always shifted by
   /// this function.
   /// Returns PicardConvergenceData object that contains the
   /// vital stats of the iteration. 
   PicardConvergenceData unsteady_partitioned_solve(const double& dt)
   {
    //We shift the values, so shift_values is true
    return unsteady_partitioned_solve(dt,true);
   }


   /// \short Unsteady partitioned solver. Advance time by dt and solve
   /// the system by a partitioned method. The boolean flag is used to
   /// control whether the time values should be shifted. If it is true 
   /// the current data values will be shifted (stored as previous 
   /// timesteps) before the solution step.
   /// Returns PicardConvergenceData object that contains the
   /// vital stats of the iteration. 
   PicardConvergenceData unsteady_partitioned_solve(const double& dt,
                                                   const bool &shift_values)
   {
    //Shift the time values and the dts according to the control flag
    if(shift_values) {shift_time_values();}

    // Advance global time and set current value of dt 
    time_pt()->time()+=dt;
    time_pt()->dt()=dt;
     
    //Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();
     
    //Loop over them all and set the weights
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      time_stepper_pt(i)->set_weights();
     }
     
    //Now update anything that needs updating before the timestep
    //This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();
     
    // Extrapolate the vm data and then update fluid mesh during unsteady run
    extrapolate_data();
     
    // Create object to doc convergence stats
    PicardConvergenceData conv_data;
     
    try
     {
      //Solve the non-linear problem for this timestep with Newton's method
      conv_data = partitioned_solve();
     }
    //Catch any exceptions thrown in the partitioned solver
    catch(PartitionedSolverError &error)
     {
      if (!error.Ran_out_of_iterations)
       {
        std::ostringstream error_stream;
        error_stream << "Error occured in Partitioned solver. "
                     << std::endl;
        throw OomphLibError(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
       }
      else
       {
        oomph_info << "Note: Ran out of iterations but continuing anyway"
                   << std::endl;
       }
     }
     
    //Now update anything that needs updating after the timestep
    actions_after_implicit_timestep();

    return conv_data;
   }








   PicardConvergenceData unsteady_semi_implicit_solve(const double& dt,
                                                   const bool &shift_values)
   {
    clock_t t_start = clock();

    //Shift the time values and the dts according to the control flag
    if(shift_values) {shift_time_values();}

    // Advance global time and set current value of dt 
    time_pt()->time()+=dt;
    time_pt()->dt()=dt;

    //Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();
     
    //Loop over them all and set the weights
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      time_stepper_pt(i)->set_weights();
     }
     
    //Now update anything that needs updating before the timestep
    //This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();
    


    // Extrapolate the vm data and then update fluid mesh during unsteady run
    extrapolate_data();

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;

    //Prepare for explicit cell solving - we just unpin everything
    unpin_all_dofs();
    assign_eqn_numbers();

    //Perform explicit cell solve
    explicit_cell_model_timestep(dt);
    
    //Get updated current timestep cell model data
    update_cell_model_data();  

    //Prepare for implicit diffusion model solving
    restore_dofs();
    pin_cell_variable_dofs();
    assign_eqn_numbers();

    //Perform solve diffusion model
    newton_solve();


    //Get updated current timestep cell model data
    unpin_all_dofs();
    assign_eqn_numbers();
    update_cell_model_data();




    actions_after_implicit_timestep();

    clock_t t_end = clock();

    oomph_info << "Time taken to complete semi implicit partitioned timestep " << double(t_end - t_start)/CLOCKS_PER_SEC << std::endl;

    return conv_data;
   }


   //Perform a timestep with strang splitting
   PicardConvergenceData unsteady_strang_splitting_solve(const double& dt,
                                                                  const bool &shift_values)
   {

    // oomph_info << std::endl << "Attempting to solve using implicit strang splitting method over a timestep " << dt << std::endl;

    oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nAttempting to solve using implicit strang splitting method over a timestep " << dt << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;

    clock_t t_start = clock();

    //Shift the time values and the dts according to the control flag
    if(shift_values) {shift_time_values();}


    
     
    //Now update anything that needs updating before the timestep
    //This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();

    
    // Advance global time and set current value of dt/2.0 and set dt = dt/2
    time_pt()->time()+=dt;
    time_pt()->dt()=dt;
    //Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();
    //Loop over them all and set the weights
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      time_stepper_pt(i)->set_weights();
     }




    // Store the current values of the vm dofs as reference values
    // and for the pointwise Aitken extrapolation 
    store_dofs();

    // Extrapolate the vm data and then update fluid mesh during unsteady run
    extrapolate_data();

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;
    
    unsigned iter;

    for(iter=1;iter<=Max_picard;iter++){
      ///////////////////////////////////////////////////////////////////////////////////
      //UPDATE CELL MODEL DATA
      ////////////////////////////////////////////////////////////////////////////////////
      //prepare for update cell model data
      unpin_all_dofs();
      assign_eqn_numbers();
      //Get updated current timestep cell model data
      update_cell_model_data();

      ////////////////////////////////////////////////////////////////////////////////////
      //CELL SOLVE OVER DT/2.0
      ////////////////////////////////////////////////////////////////////////////////////
      //Perform explicit cell solve
      explicit_cell_model_timestep(dt/2.0, false);

      ////////////////////////////////////////////////////////////////////////////////////
      //UPDATE CELL MODEL DATA
      ////////////////////////////////////////////////////////////////////////////////////
      //prepare for update cell model data
      unpin_all_dofs();
      assign_eqn_numbers();
      //Get updated current timestep cell model data
      update_cell_model_data();



      ////////////////////////////////////////////////////////////////////////////////////
      //DIFFUSION SOLVE OVER DT
      ////////////////////////////////////////////////////////////////////////////////////
      //Prepare for implicit diffusion model solving
      restore_dofs();
      pin_cell_variable_dofs();
      assign_eqn_numbers();
      //Perform solve diffusion model
      newton_solve();


      ////////////////////////////////////////////////////////////////////////////////////
      //UPDATE CELL MODEL DATA
      ////////////////////////////////////////////////////////////////////////////////////
      //prepare for update cell model data
      unpin_all_dofs();
      assign_eqn_numbers();
      //Get updated current timestep cell model data
      update_cell_model_data();


      ////////////////////////////////////////////////////////////////////////////////////
      //CELL SOLVE OVER DT/2.0
      ////////////////////////////////////////////////////////////////////////////////////
      //Perform explicit cell solve
      explicit_cell_model_timestep(dt/2.0, true);


      //Check for convergence
      double rms_change;
      double max_change;
      double rms_norm;
      get_change(rms_change,max_change,rms_norm);

      oomph_info << "==================================================\n";
      oomph_info <<   "Iteration             : " 
                 << iter << std::endl;
      oomph_info <<   "RMS  change           : "
                 << rms_change << std::endl;
      oomph_info <<   "Max. change           : "
                 << max_change << std::endl;
      oomph_info <<   "RMS norm              : "
                 << rms_norm   << std::endl;
      oomph_info << "==================================================\n\n";

      //If relative change is smaller than the tolerance for convergence
      if(std::fabs(rms_change/rms_norm)<Convergence_tolerance){
        break;
      }
      else{
        // Store the current values of the vm dofs as reference values
        // and for the pointwise Aitken extrapolation 
        store_dofs();
      }
    }//End of iteration loop

    if(iter==(Max_picard-1)){
      throw OomphLibError("Implcit iterative strang-splitting could not converge",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }



    actions_after_implicit_timestep();

    clock_t t_end = clock();

    // oomph_info << "Time taken to complete implicit strang-splitting timestep " << double(t_end - t_start)/CLOCKS_PER_SEC << std::endl;

    oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nTime taken to complete implicit strang-splitting timestep " << double(t_end - t_start)/CLOCKS_PER_SEC << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;

    return conv_data;
   }






   //Perform a timestep with strang splitting
   PicardConvergenceData unsteady_single_step_strang_splitting_solve(const double& dt,
                                                                  const bool &shift_values)
   {

    // oomph_info << std::endl << "Attempting to solve using implicit strang splitting method over a timestep " << dt << std::endl;

    oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nAttempting to solve using implicit strang splitting method over a timestep " << dt << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;

    clock_t t_start = clock();

    //Shift the time values and the dts according to the control flag
    if(shift_values) {shift_time_values();}


    
     
    //Now update anything that needs updating before the timestep
    //This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();

    
    // Advance global time and set current value of dt/2.0 and set dt = dt/2
    time_pt()->time()+=dt;
    time_pt()->dt()=dt;
    //Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();
    //Loop over them all and set the weights
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      time_stepper_pt(i)->set_weights();
     }




    // Extrapolate the vm data and then update fluid mesh during unsteady run
    extrapolate_data();

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;
    
    ///////////////////////////////////////////////////////////////////////////////////
    //UPDATE CELL MODEL DATA
    ////////////////////////////////////////////////////////////////////////////////////
    //prepare for update cell model data
    unpin_all_dofs();
    assign_eqn_numbers();
    //Get updated current timestep cell model data
    update_cell_model_data();

    ////////////////////////////////////////////////////////////////////////////////////
    //CELL SOLVE OVER DT/2.0
    ////////////////////////////////////////////////////////////////////////////////////
    //Perform explicit cell solve
    oomph_info << "Solving explicit cell half step" << std::endl;
    explicit_cell_model_timestep(dt/2.0, false);

    ////////////////////////////////////////////////////////////////////////////////////
    //UPDATE CELL MODEL DATA
    ////////////////////////////////////////////////////////////////////////////////////
    //prepare for update cell model data
    unpin_all_dofs();
    assign_eqn_numbers();
    //Get updated current timestep cell model data
    update_cell_model_data();



    ////////////////////////////////////////////////////////////////////////////////////
    //DIFFUSION SOLVE OVER DT
    ////////////////////////////////////////////////////////////////////////////////////
    //Prepare for implicit diffusion model solving
    restore_dofs();
    pin_cell_variable_dofs();
    assign_eqn_numbers();
    //Perform solve diffusion model
    newton_solve();


    ////////////////////////////////////////////////////////////////////////////////////
    //UPDATE CELL MODEL DATA
    ////////////////////////////////////////////////////////////////////////////////////
    //prepare for update cell model data
    unpin_all_dofs();
    assign_eqn_numbers();
    //Get updated current timestep cell model data
    update_cell_model_data();


    ////////////////////////////////////////////////////////////////////////////////////
    //CELL SOLVE OVER DT/2.0
    ////////////////////////////////////////////////////////////////////////////////////
    //Perform explicit cell solve
    oomph_info << "Solving 'implicit' cell half step" << std::endl;
    explicit_cell_model_timestep(dt/2.0, true);



    actions_after_implicit_timestep();

    clock_t t_end = clock();

    // oomph_info << "Time taken to complete implicit strang-splitting timestep " << double(t_end - t_start)/CLOCKS_PER_SEC << std::endl;

    oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nTime taken to complete single step semi-implicit strang-splitting timestep " << double(t_end - t_start)/CLOCKS_PER_SEC << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;

    return conv_data;
   }

   







   //Perform a timestep with whiteley 2006 - solve diffusion model with Iion of previous timestep,
   // solve for cell variables implicitely using current timestep vm (implicit boost solve is not implemented yet so we "implicit" solve, uising an adaptive methods)
   // This is SUPER unstable
   PicardConvergenceData unsteady_whitely_2006_solve(const double& dt,
                                                    const bool &shift_values)
   {


    oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nAttempting to solve using implicit strang splitting method over a timestep " << dt << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;

    clock_t t_start = clock();

    //We get the Iion value from the previous timestep and the Vm value from the current timestep
    set_cells_interpolate_flags(Explicit, Implicit);

    //Shift the time values and the dts according to the control flag
    if(shift_values) {shift_time_values();}
     
    //Now update anything that needs updating before the timestep
    //This could be time-dependent boundary conditions, for example.
    actions_before_implicit_timestep();

    
    // Advance global time and set current value of dt and set dt = dt
    time_pt()->time()+=dt;
    time_pt()->dt()=dt;
    //Find out how many timesteppers there are
    unsigned n_time_steppers = ntime_stepper();
    //Loop over them all and set the weights
    for(unsigned i=0;i<n_time_steppers;i++)
     {
      time_stepper_pt(i)->set_weights();
     }




    // Extrapolate the vm data and then update fluid mesh during unsteady run
    extrapolate_data();

    // Create object to doc convergence stats
    PicardConvergenceData conv_data;



    
    ////////////////////////////////////////////////////////////////////////////////////
    //DIFFUSION SOLVE OVER DT
    ////////////////////////////////////////////////////////////////////////////////////
    //Prepare for implicit diffusion model solving
    restore_dofs();
    pin_cell_variable_dofs();
    assign_eqn_numbers();
    //Perform solve diffusion model
    newton_solve();


    ////////////////////////////////////////////////////////////////////////////////////
    //UPDATE CELL MODEL DATA
    ////////////////////////////////////////////////////////////////////////////////////
    //prepare for update cell model data
    unpin_all_dofs();
    assign_eqn_numbers();
    //Get updated current timestep cell model data
    update_cell_model_data();


    ////////////////////////////////////////////////////////////////////////////////////
    //CELL SOLVE OVER DT
    ////////////////////////////////////////////////////////////////////////////////////
    //Perform explicit cell solve
    oomph_info << "Solving 'implicit' cell half step" << std::endl;
    explicit_cell_model_timestep(dt, true);


    ////////////////////////////////////////////////////////////////////////////////////
    //UPDATE CELL MODEL DATA - Get Iion at the current timestep for the next step solve
    ////////////////////////////////////////////////////////////////////////////////////
    //prepare for update cell model data
    unpin_all_dofs();
    assign_eqn_numbers();
    //Get updated current timestep cell model data
    update_cell_model_data();


    actions_after_implicit_timestep();

    clock_t t_end = clock();

    oomph_info 
         << "\n\n\n////////////////////////////////////////////////////////"
         << "\nTime taken to complete single step semi-implicit strang-splitting timestep " << double(t_end - t_start)/CLOCKS_PER_SEC << std::endl
         << "////////////////////////////////////////////////////////\n\n\n"
         << std::endl;

    return conv_data;
   }




   //Perform a timestep with modified whiteley 2006 - solve diffusion model with Iion of previous timestep,
   // solve for cell variables implicitely using current timestep vm (implicit boost solve is not implemented yet so we "implicit" solve, uising an adaptive methods).
   //Original whiteley 2006 method is a single step partitioned solver - instead force both diffusion and cell solver to use crank-nicolson method and iteratively
   // seek a solution
   









   /// \short Assess convergence based on max. residual of coupled system of 
   /// eqns. The argument specifies the convergence tolerance.
   void assess_convergence_based_on_max_global_residual(
    const double& tol)
    {
     Convergence_criterion=Assess_convergence_based_on_max_global_residual;
     Convergence_tolerance=tol;
    }

   /// \short Assess convergence based on max. residuals of coupled
   /// system of eqns. This interface has no argument
   /// and the default convergence tolerance
   /// for the Newton solver, Problem::Newton_solver_tolerance is used.
   void assess_convergence_based_on_max_global_residual()
    {
     assess_convergence_based_on_max_global_residual(
      Problem::Newton_solver_tolerance);
    }

   /// \short Assess convergence based on max. absolute change of vm
   /// dofs. The argument specifies the convergence tolerance.
   void assess_convergence_based_on_absolute_vm_change(
    const double& tol)
    {
     Convergence_criterion=Assess_convergence_based_on_absolute_vm_change;
     Convergence_tolerance=tol;
    }

   /// \short Assess convergence based on max. absolute change of vm
   /// dofs. This interface has no argument and the default 
   /// convergence tolerance
   /// for the Newton solver, Problem::Newton_solver_tolerance is used.
   void assess_convergence_based_on_absolute_vm_change()
    {
     assess_convergence_based_on_absolute_vm_change(
      Problem::Newton_solver_tolerance);
    }
   
   /// \short Assess convergence based on max. relative change of vm
   /// dofs. The argument specifies the convergence tolerance.
   void assess_convergence_based_on_relative_vm_change(
    const double& tol)
    {
     Convergence_criterion=Assess_convergence_based_on_relative_vm_change;
     Convergence_tolerance=tol;
    }
      
   /// \short Assess convergence based on max. relative change of vm
   /// dofs. This interface has no argument and the default 
   /// convergence tolerance
   /// for the Newton solver, Problem::Newton_solver_tolerance is used.
   void assess_convergence_based_on_relative_vm_change()
    {
     assess_convergence_based_on_relative_vm_change(
      Problem::Newton_solver_tolerance);
    }


   /// \short Use pointwise Aitken extrapolation. The argument is used to 
   /// specify the Picard iteration after which pointwise Aitken extrapolation 
   /// is to be used for the first time. 
   void enable_pointwise_aitken(const unsigned& pointwise_aitken_start)
    {
     Pointwise_aitken_start = pointwise_aitken_start;
     Use_pointwise_aitken = true;
    }

   /// \short Use pointwise Aitken extrapolation. This interface has
   /// no argument and the current value of Pointwise_aitken_start will
   /// be used. The default is zero, extrapolation starts immediately
   void enable_pointwise_aitken() {Use_pointwise_aitken = true;}

   /// \short Disable the use of pointwise Aitken extrapolation
   void disable_pointwise_aitken() {Use_pointwise_aitken = false;}

   ///\short Use under-relaxation and (optionally) specify under-relaxation 
   /// parameter. Default: omega=1.0 (i.e. no actual under-relaxation; 
   /// Other extreme: omega=0.0 (freeze wall shape). Under-relaxation
   /// parameter can also be computed dynamically by setting
   /// use_irons_and_tuck_extrapolation()
   void enable_under_relaxation(const double& omega=1.0)
    {Omega_relax=omega;}

   ///\short Use Irons and Tuck extrapolation for vm dofs
   void enable_irons_and_tuck_extrapolation()
   {Use_irons_and_tuck_extrapolation = true;}

   ///\short Do not use Irons and Tuck extrapolation for vm dofs
   void disable_irons_and_tuck_extrapolation()
   {Use_irons_and_tuck_extrapolation = false;}
   
   /// Enumerated flags for convergence criteria
   enum convergence_criteria{Assess_convergence_based_on_absolute_vm_change,
                             Assess_convergence_based_on_relative_vm_change,
                             Assess_convergence_based_on_max_global_residual};


   /// \short Get rms of change in the vm dofs; the max. change of the
   /// vm dofs and the rms norm of the vm dofs themselves.
   /// Change is computed relative to the reference values stored when
   /// store_vm_dofs() was last called.
   void get_change(double& rms_change, 
                   double& max_change,
                   double& rms_norm)
   {
    // Initialise
    rms_change=0.0;
    max_change=0.0;
    rms_norm=0.0;

    // Counter for the number of values:
    unsigned value_count=0;

    //Loop over meshes
    unsigned orig_n_sub_mesh=Orig_sub_mesh_pt.size();
    for (unsigned i=0;i<orig_n_sub_mesh;i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            // Change
            double change=Previous_value[value_count]-
            elem_pt->node_pt(l)->value(k);

            // Max change?
            if (std::fabs(change)>max_change) max_change=std::fabs(change);

            //Add square of change relative to previous value
            rms_change+=pow(change,2);

            //Add square of previous value
            rms_norm+=pow(Previous_value[value_count],2);

            // Increment counter
            value_count++;
          }
        }
      }
    }

    // Turn into rms:
    rms_change=sqrt(rms_change/double(value_count));
    rms_norm=sqrt(rms_norm/double(value_count));
   }

   /// \short Store the current dof values as reference values for 
   /// future convergence check. Also add another entry to pointwise 
   /// Aitken history if required.
   void store_dofs()
   {
    oomph_info << Pointwise_aitken_counter << std::endl;
    // Counter for the number of values:
    unsigned value_count=0;

    //Loop over meshes
    unsigned orig_n_sub_mesh=Orig_sub_mesh_pt.size();
    for(unsigned i=0;i<orig_n_sub_mesh;i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            // Make backup
            Previous_value[value_count]=elem_pt->node_pt(l)->value(k);

            // Store in pointwise Aitken history
            if (Use_pointwise_aitken&&(Pointwise_aitken_counter>=0)){
              Pointwise_aitken_value[value_count][Pointwise_aitken_counter]=elem_pt->node_pt(l)->value(k);
            }
            // Increment counter
            value_count++;
          }
        }
      }
    }
    // we stored another level of Aitken history values
    Pointwise_aitken_counter++;
  }

   /// \short Reset timer
   void reset_timer()
    {
     T_spent_on_actual_solve=0.0;
     T_ref=clock();
     Timer_has_been_halted=false;
    }


   /// \short (Re-)start timer (e.g. after completing non-essential
   /// parts of the code such as documentation of the iteration's
   /// progress)
   void restart_timer()
    {
     T_ref=clock();
     Timer_has_been_halted=false;
    }


   /// \short Halt timer (e.g. before performing non-essential
   /// parts of the code such as documentation of the iteration's
   /// progress)
   void halt_timer()
    {
     if (!Timer_has_been_halted)
      {
       T_spent_on_actual_solve+=double(clock()-T_ref)/CLOCKS_PER_SEC;
       Timer_has_been_halted=true;
      }
    }


   /// \short Total elapsed time since start of solve
   double t_spent_on_actual_solve()
    {
     halt_timer();
     double time=T_spent_on_actual_solve;
     restart_timer();
     return time;
    }


    protected:

   /// Rebuild global mesh for monolithic discretisation
   void rebuild_monolithic_mesh()
   {
    // Get rid of the previous submeshes
    flush_sub_meshes();

    // Add original submeshes
    unsigned orig_n_sub_mesh=Orig_sub_mesh_pt.size();
    for (unsigned i=0;i<orig_n_sub_mesh;i++)
     {
      add_sub_mesh(Orig_sub_mesh_pt[i]);
     }

    // Rebuild global mesh
    rebuild_global_mesh();
   }

   /// \short Number of Aitken histories available (int because after
   /// extrapolation it's re-initialised to -1 to force the computation
   /// of three new genuine iterates).
   int Pointwise_aitken_counter;

   /// Use pointwise Aitken extrapolation?
   bool Use_pointwise_aitken;

   /// \short Start pointwise Aitken extrpolation after specified number 
   /// of Picard iterations
   unsigned Pointwise_aitken_start;
   
   /// Convergence tolerance for Picard iteration
   double Convergence_tolerance;

   /// Max. number of Picard iterations
   unsigned Max_picard;

   /// Doc maximum global residual during iteration? (default: false)
   bool Doc_max_global_residual;

   /// Pin/unpin cell dofs
   void pin_cell_variable_dofs()
   {
    //Loop over elements
    unsigned n_elem = Cell_mesh_pt->nelement();
    for(unsigned e=0; e<n_elem; e++){
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->pin_all_cell_vars();
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->Extern_Has_Intentionally_Blanket_Pinned_All_Cell_Variables();
    }
   }
   void unpin_cell_variable_dofs()
   {
    //Loop over elements
    unsigned n_elem = Cell_mesh_pt->nelement();
    for(unsigned e=0; e<n_elem; e++){
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->unpin_all_cell_vars();
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->Extern_Has_Intentionally_Blanket_UnPinned_All_Cell_Variables();
    }
   }

   void pin_all_dofs()
   {
    //Loop over meshes
    for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            elem_pt->node_pt(l)->pin(k);
          }
        }
      }
    }
    //Loop over elements
    unsigned n_elem = Cell_mesh_pt->nelement();
    for(unsigned e=0; e<n_elem; e++){
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->Extern_Has_Intentionally_Blanket_Pinned_All_Cell_Variables();
    }
   }

   void unpin_all_dofs()
   {
    //Loop over meshes
    for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            elem_pt->node_pt(l)->unpin(k);
          }
        }
      }
    }
    //Loop over elements
    unsigned n_elem = Cell_mesh_pt->nelement();
    for(unsigned e=0; e<n_elem; e++){
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->Extern_Has_Intentionally_Blanket_UnPinned_All_Cell_Variables();
    }
   }

   /// Restore pinned status of all dofs
   void restore_dofs()
   {
    //Loop over meshes
    for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            //grab suitable status from vector
            if(Value_is_pinned[i][e][l][k]){
              elem_pt->node_pt(l)->pin(k);
            }   
            else{
              elem_pt->node_pt(l)->unpin(k);
            }
          }
        }
      }
    }

    //We cannot necessarily say that we have unpinned/pinned all the cell variables in all elements but we can't say the opposite either
    // be safe and assume they have all been unpinned
    unsigned n_elem = Cell_mesh_pt->nelement();
    for(unsigned e=0; e<n_elem; e++){
      dynamic_cast<ELEMENT_TYPE*>(Cell_mesh_pt->element_pt(e))->Extern_Has_Intentionally_Blanket_UnPinned_All_Cell_Variables();
    }
  }






   /// Do pointwise Aitken extrapolation for vm
   void pointwise_aitken_extrapolate()
   {
    // Counter for the number of values:
    unsigned value_count=0;

    //Loop over meshes
    for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            // Shorthand
            double v0=Pointwise_aitken_value[value_count][0];
            double v1=Pointwise_aitken_value[value_count][1];
            double v2=Pointwise_aitken_value[value_count][2];

            double new_value=v2;

            double max_diff=std::max(std::fabs(v1-v0),std::fabs(v2-v1));
            if (max_diff>1.0e-10)
            {
              new_value=v2-std::pow((v2-v1),int(2))/(v2-2.0*v1+v0);
            }
            elem_pt->node_pt(l)->set_value(k,new_value);
            // Increment counter
            value_count++;
          }
        }
      }
    }

    // Reset the counter for the Aitken convergence check
    // (setting counter to -1 forces three new genuine
    // iterates to be computed). 
    Pointwise_aitken_counter=-1;
   
   }

   //Store the pinned status of each of the meshes, data elements, nodes, variable
   // so that after we unpin everything we can re-pin boundary conditions
   Vector<Vector<Vector<std::vector<bool> > > > Value_is_pinned;

   /// \short Vector storing the previous vm values -- used for
   /// convergence check
   Vector<double> Previous_value;

   /// \short Mesh containing only cell elements -- the elements in this
   /// Mesh will be excluded from the assembly process when
   /// the vm problem is solved
   Mesh* Cell_mesh_pt;
   
   /// \short Backup for the pointers to the submeshes in the original problem
   Vector<Mesh*> Orig_sub_mesh_pt;

   /// Vector of changes in Irons and Tuck under-relaxation
   Vector<double> Del_irons_and_tuck;

   /// Irons and Tuck relaxation factor
   double R_irons_and_tuck;

   /// \short Vector of Vectors containing up to three previous
   /// iterates for the vm dofs; used for pointwise Aitken extrapolation
   Vector<Vector<double> > Pointwise_aitken_value;

   /// Have we just done a pointwise Aitken step
   bool Recheck_convergence_after_pointwise_aitken;

    private:

   /// Extrapolate vm data and update cell mesh during unsteady run
   void extrapolate_data()
   {
    //Loop over meshes
    for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
      //Loop over elements
      unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
      for(unsigned e=0; e<n_elem; e++){
        ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
        unsigned n_node = elem_pt->nnode();
        //Loop over nodes
        for(unsigned l=0; l<n_node; l++){
          //Loop over the variables
          for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
            // Linear extrapolation based on previous two values,
            // assuming constant timestep.
            double new_value=2.0*elem_pt->node_pt(l)->value(1,k)-
             elem_pt->node_pt(l)->value(2,k);
            *(elem_pt->node_pt(l)->value_pt(0,k))=new_value;
          }
        }
      }
    }
  }

  /// \short Under-relax the most recently computed vm variables, either
  /// by classical relaxation or by Irons & Tuck
  void under_relax_data()
  {
    // Irons and Tuck extrapolation/relaxation; an extension of Aitken's method
    //-------------------------------------------------------------------------
    if (Use_irons_and_tuck_extrapolation)
    {
      double top=0.0;
      double den=0.0;
      double crit=0.0;

      // Counter for the number of values:
      unsigned value_count=0;
       
      //Loop over meshes
      for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
        //Loop over elements
        unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
        for(unsigned e=0; e<n_elem; e++){
          ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
          unsigned n_node = elem_pt->nnode();
          //Loop over nodes
          for(unsigned l=0; l<n_node; l++){
            //Loop over the variables
            for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){


              // Prediction from vm solver
              double new_vm_value=elem_pt->node_pt(l)->value(k);
               
              // Previus value
              double old_vm_value=Previous_value[value_count];
               
              // Change
              double change=old_vm_value-new_vm_value;

              // Change of change
              double del2=Del_irons_and_tuck[value_count]-change;

              // Update change
              Del_irons_and_tuck[value_count]=change;

              // Update top
              top+=del2*change;
               
              // Update denominator
              den+=del2*del2;

              // Update convergence criterion
              crit+=std::fabs(change);

              // Increment counter
              value_count++;
            }
          }
        }
      }

      // Update relaxation factor. The if buffers the case in which
      // we haven't realised that we've converged (so that den=0).
      // This can happen, e.g. if the convergence assessment is based on the
      // global residual or during validation. In that case we 
      // obviously don't want any changes to the iterates.
      if (den!=0.0)
       {
        double new_r=R_irons_and_tuck+(R_irons_and_tuck-1.0)*top/den;
        R_irons_and_tuck=new_r;
       }
      else
       {
        R_irons_and_tuck=0.0;
       }

      //Loop over meshes
      for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
        //Loop over elements
        unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
        for(unsigned e=0; e<n_elem; e++){
          ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
          unsigned n_node = elem_pt->nnode();
          //Loop over nodes
          for(unsigned l=0; l<n_node; l++){
            //Loop over the variables
            for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
              // Compute relaxed/extrapolated value
              double new_value=elem_pt->node_pt(l)->value(k)+
               R_irons_and_tuck*Del_irons_and_tuck[value_count];

              // Assign 
              elem_pt->node_pt(l)->set_value(k,new_value);
           
              // Increment counter
              value_count++;
            }
          }
        }
      }
      return;
    }

    // Standard relaxation
    //--------------------
    else
    {
       
       
      // No relaxation: Can return immediately
      if (Omega_relax==1.0) return;
       
      // Counter for the number of values:
      unsigned value_count=0;
       
      // Number of vm Data items:
      // unsigned n_data=Vm_data_pt.size();
       
      //Loop over meshes
      for(unsigned i=0;i<Orig_sub_mesh_pt.size();i++){
        //Loop over elements
        unsigned n_elem = Orig_sub_mesh_pt[i]->nelement();
        for(unsigned e=0; e<n_elem; e++){
          ELEMENT_TYPE* elem_pt = dynamic_cast<ELEMENT_TYPE*>(Orig_sub_mesh_pt[i]->element_pt(e));
          unsigned n_node = elem_pt->nnode();
          //Loop over nodes
          for(unsigned l=0; l<n_node; l++){
            //Loop over the variables
            for(unsigned k=0; k<elem_pt->required_nvalue(l); k++){
           
              // Prediction from vm solver
              double new_vm_value=elem_pt->node_pt(l)->value(k);
           
              // Previus value
              double old_vm_value=Previous_value[value_count];
           
              // Relax
              elem_pt->node_pt(l)->set_value(k,Omega_relax*new_vm_value+
                                          (1.0-Omega_relax)*old_vm_value);
                    
              // Increment counter
              value_count++;
            }
          }
        }
      }
    }

  }



   /// \short Under-relaxation parameter. (1.0: no under-relaxation; 
   /// 0.0: Freeze wall shape)
   double Omega_relax; 

   /// \short Boolean flag to indicate use of Irons and Tuck's extrapolation 
   /// for vm values
   bool Use_irons_and_tuck_extrapolation;

   /// Convergence criterion (enumerated flag)
   int Convergence_criterion;
  
   /// \short Reference time for partitioned solve. Can be 
   /// re-initialised whenever total elapsed time has been stored
   /// (before entering non-essential doc sections of the code)
   clock_t T_ref;

   /// \short Total elapsed time since start of solve, can be
   /// accumulated by adding bits of time spent in relevant parts of
   /// code (bypassing sections that only document the progress)
   double T_spent_on_actual_solve;

   /// \short boolean flag to indicate if timer has been halted
   bool Timer_has_been_halted;

   bool PartitionedSolveType;

  };

}



#endif
