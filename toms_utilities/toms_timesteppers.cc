#include "../generic/timesteppers.h"

namespace oomph
{

//====================================================================
 ///Assign the values of the weights; pass the value of the timestep
//====================================================================
 template<>
 void  BDF<5>::set_weights()
  {
#ifdef PARANOID
   double dt0=Time_pt->dt(0);
   for (unsigned i=0;i<Time_pt->ndt();i++)
    {
     if (dt0!=Time_pt->dt(i))
      {
       throw OomphLibError(
        "BDF4 currently only works for fixed timesteps \n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
   double dt=Time_pt->dt(0);
   Weight(1,0) =  137.0/60.0/dt;
   Weight(1,1) = -300.0/60.0/dt;
   Weight(1,2) =  300/60.0/dt;
   Weight(1,3) = -200.0/60.0/dt;
   Weight(1,4) =   75.0/60.0/dt;
   Weight(1,5) =  -12.0/60.0/dt;
  }


//======================================================================
///Calculate the predictor weights
//======================================================================
template<>
void BDF<5>::set_predictor_weights()
{
 //Read the value of the previous timesteps
 //double dt=Time_pt->dt(0);
 //double dtprev=Time_pt->dt(1);

 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
}

//=======================================================================
///Calculate the predicted values and store them at the appropriate
///location in the data structure
///This function must be called after the time-values have been shifted!
//=======================================================================
template<>
void BDF<5>::calculate_predicted_values(Data* const &data_pt)
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
 
}

///Calculate predictions for the positions
template<>
void BDF<5>::calculate_predicted_positions(Node* const &node_pt)
  {
   throw OomphLibError("Not implemented yet",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

///Function that sets the error weights
template<>
void BDF<5>::set_error_weights()
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
}

//===================================================================
///Function to compute the error in position i at node
//===================================================================
template<>
double BDF<5>::temporal_error_in_position(Node* const &node_pt, 
                                          const unsigned &i)
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
 return 0.0;
}
   

//=========================================================================
///Function to calculate the error in the data value i
//=========================================================================
template<>
double BDF<5>::temporal_error_in_value(Data* const &data_pt, const unsigned &i)
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
 return 0.0;
}








//====================================================================
 ///Assign the values of the weights; pass the value of the timestep
//====================================================================
 template<>
 void  BDF<6>::set_weights()
  {
#ifdef PARANOID
   double dt0=Time_pt->dt(0);
   for (unsigned i=0;i<Time_pt->ndt();i++)
    {
     if (dt0!=Time_pt->dt(i))
      {
       throw OomphLibError(
        "BDF4 currently only works for fixed timesteps \n",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
   double dt=Time_pt->dt(0);
   Weight(1,0) =  147.0/60.0/dt;
   Weight(1,1) = -360.0/60.0/dt;
   Weight(1,2) =  450/60.0/dt;
   Weight(1,3) = -400.0/60.0/dt;
   Weight(1,4) =  225.0/60.0/dt;
   Weight(1,5) = -72.0/60.0/dt;
   Weight(1,6) =  10.0/60.0/dt;
  }


//======================================================================
///Calculate the predictor weights
//======================================================================
template<>
void BDF<6>::set_predictor_weights()
{
 //Read the value of the previous timesteps
 //double dt=Time_pt->dt(0);
 //double dtprev=Time_pt->dt(1);

 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
}

//=======================================================================
///Calculate the predicted values and store them at the appropriate
///location in the data structure
///This function must be called after the time-values have been shifted!
//=======================================================================
template<>
void BDF<6>::calculate_predicted_values(Data* const &data_pt)
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
 
}

///Calculate predictions for the positions
template<>
void BDF<6>::calculate_predicted_positions(Node* const &node_pt)
  {
   throw OomphLibError("Not implemented yet",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

///Function that sets the error weights
template<>
void BDF<6>::set_error_weights()
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
}

//===================================================================
///Function to compute the error in position i at node
//===================================================================
template<>
double BDF<6>::temporal_error_in_position(Node* const &node_pt, 
                                          const unsigned &i)
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
 return 0.0;
}
   

//=========================================================================
///Function to calculate the error in the data value i
//=========================================================================
template<>
double BDF<6>::temporal_error_in_value(Data* const &data_pt, const unsigned &i)
{
 throw OomphLibError("Not implemented yet",
                     OOMPH_CURRENT_FUNCTION,
                     OOMPH_EXCEPTION_LOCATION);
 return 0.0;
}


}//End namespace