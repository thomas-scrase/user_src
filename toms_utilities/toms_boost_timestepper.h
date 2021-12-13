#ifndef BOOST_NUMERIC_ODEINT_STEPPER_CELL_HEUN_EULER_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_CELL_HEUN_EULER_HPP_INCLUDED


#include <boost/numeric/odeint/stepper/base/explicit_stepper_base.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

namespace boost {
namespace numeric {
namespace odeint {

	



	//Implementation of an explicit euler-heun method. Checks for non-finite values after initial step to avoid unnecessary
	// computation after corrupted data has been generated

	template<typename State, typename Value = double, typename Deriv = State, 
			typename Time = Value, 
			typename Algebra = typename algebra_dispatcher< State >::algebra_type, 
			typename Operations = typename operations_dispatcher< State >::operations_type, 
			typename Resizer = initially_resizer>
	class cell_heun_euler : public explicit_stepper_base {
	public:
		// types
		typedef explicit_stepper_base< euler< ... >,... > stepper_base_type;
		typedef stepper_base_type::state_type             state_type;
		typedef stepper_base_type::value_type             value_type;
		typedef stepper_base_type::deriv_type             deriv_type;
		typedef stepper_base_type::time_type              time_type;
		typedef stepper_base_type::algebra_type           algebra_type;
		typedef stepper_base_type::operations_type        operations_type;
		typedef stepper_base_type::resizer_type           resizer_type;

		cell_heun_euler(const algebra_type & algebra = algebra_type()) : stepper_base_type(algebra)
		{	}

		// public member functions
		template<typename System, typename StateIn> 
		void do_step(System system, const StateIn &x, time_type t, time_type dt) const
		{
			//Declare the storage we need
			StateOut x_np1_1;
			StateOut x_np1_2;

			deriv_type K1;
			deriv_type K2;
			system.first(x, K1);

			//Calculate derivative at the start
			norm_result_type<StateOut> norm_val;

			//The error in the solution
			value_type tau;

			//The end time
			time_type t_end = t+dt;
			while(t<t_end)
			{
				//Calculate the initial guess with explict euler
				stepper_base_type::m_algebra.for_each3( x_np1_1 , x , K1 ,
	                typename operations_type::template scale_sum2< value_type , time_type >( 1.0 , dt ) );

				norm_val = norm_inf<StateOut>(x_np1_1)

				//If the value is non-finite halve the timestep and try again
				if(norm == static_cast<norm_result_type<StateOut>>(std::nanf)){dt *= 0.5; continue;}

				//Calculate derivative at the end point
				system.first(x_np1_1, K2);

				//Calculate the second guess at the end point

				//Calculate the initial guess with explict euler
				stepper_base_type::m_algebra.for_each3( x_np1_1 , x , K1 , K2,
	                typename operations_type::template scale_sum3< value_type , time_type >( 1.0 , dt*0.5, dt*0.5 ) );

				


				//Get the error
				if()
			}
		}
		
		template<typename StateOut, typename StateIn1, typename StateIn2> 
		void calc_state(StateOut &, time_type, const StateIn1 &, time_type, const StateIn2 &, time_type) const;
		
		template<typename StateType>
		void adjust_size(const StateType &);
	};
}
}
}

#endif