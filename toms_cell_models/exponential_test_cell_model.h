#ifndef OOMPH_EXPONENTIAL_TEST_CELL_HEADER
#define OOMPH_EXPONENTIAL_TEST_CELL_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "../cell_models/cell_model_base.h"

namespace oomph
{
	//dydt = V
	//dVdt = -y+I_stim

	//Oscillatory solutions:
	//d2ydt2 = -y + I_stim
	//d2Vdt2 = -V + dI_stimdt (?)
	static double K_Default_ExponentialTestCellModel = 1.0;


	class ExponentialTestCellModel : public CellModelBaseFullyPartitioned
	{
	public:
		ExponentialTestCellModel(const unsigned& number_of_backup_values) : CellModelBaseFullyPartitioned(number_of_backup_values), K(K_Default_ExponentialTestCellModel)
		{
			Names_Of_Cell_Variables = { "Vm", "y" };

			Names_Of_Output_Data = { "ActiveStrain" };

			FinalizeConstruction();
		}

		~ExponentialTestCellModel(){ }

		std::string get_cell_model_name(){return "ExponentialTestCellModel";}
		
	protected:
		double get_initial_state_variable(const unsigned &v)
		{
			if(v==0)
			{
				return 0.0;
			}
			if(v==1)
			{
				return 1.0;
			}
			
			return 0.0;
		}

		unsigned index_of_membrane_potential_in_cell_data(){return 0;}

		double return_initial_membrane_potential(const unsigned &cell_type){return 0.0;}

		void TakeTimestep(const double& dt, const double&t, double* state)
		{
			// oomph_info << "Taking timestep " << t << ", " << dt << std::endl;

			//We don't want to take a time-step larger than 0.01ms
			if(dt>0.01)
			{
				const unsigned N_solve_at_0p01 = unsigned(dt/0.01);
				const double remainder = dt - N_solve_at_0p01*0.01;
				double t_running = t;

				for(unsigned n=0; n<N_solve_at_0p01; n++)
				{
					TakeTimestep(0.01, t_running, state);
					t_running+=0.01;
				}
				if(remainder>0.0)
				{
					TakeTimestep(remainder, t_running, state);
					t_running+=remainder;
				}

				return;
			}

			state[0] += -dt*(state[1]+this->get_stimulus(t));
			state[1] += dt*K*state[0];
		}

		void get_output(double *state, double *out)
		{
		}

		void set_K(const double& new_k)
		{
			K = new_k;
		}

		double GetActiveStrain(double *state)
		{
			const double min = 0.7;
			// oomph_info << 0.5*(min-1.0)*(1.0-cos(sqrt(K)*this->time()))+1 << std::endl;
			return 0.5*(min-1.0)*(1.0-cos(sqrt(K)*this->time()));
		}

	private:
		double K;
	};

}

#endif