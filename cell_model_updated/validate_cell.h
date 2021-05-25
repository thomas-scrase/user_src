#ifndef OOMPH_VALIDATE_CELL_HEADER
#define OOMPH_VALIDATE_CELL_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "cell_model_base_updated.h"
	
namespace oomph{

	class ValidateCell : public CellModelBaseUpdated
	{
	public:
		ValidateCell()
		{
			//Lotka-Volterra function of population model
			alpha = 1.1;
			beta = 0.4;
			delta = 0.1;
			gamma = 0.4;

			//Fitzhugh Nagumo Model
			tau = 12.5;
			a = 0.7;
			b = 0.8;

			//Assign the names of the variables used by this model
			Names_Of_Cell_Variables =
			{
				"x",
				"y",
				"w"
			};
			Names_Of_Other_Parameters =
			{
				
			};
			Names_Of_Other_Variables =
			{

			};
			Names_Of_Output_Data =
			{

			};

			FinalizeConstruction();
		}

		//Assign the names of the variables used by this model
		enum Cell_Variables_Enum : unsigned
		{
			x,
			y,
			w
		};
		enum Other_Parameters_Enum : unsigned
		{

		};
		enum Other_Variables_Enum : unsigned
		{

		};
		enum Output_Data_Enum : unsigned
		{

		};

		double return_initial_state_variable(const unsigned &v, const unsigned &cell_type) override
		{
			switch(v){
				case 0: //x
					return 1.0;
					break;
				case 1: //y
					return 0.0;
					break;
				case 2: //w
					return 1.0;
					break;
			}	
		}

		double return_initial_membrane_potential(const unsigned &cell_type) override
		{
			return 0.0;
		}


		void Calculate_Derivatives(const double &Vm,
									const Vector<double> &CellVariables,
									const double &t,
									const unsigned &cell_type,
									const double &Istim,
									const Vector<double> &Other_Parameters,
									const Vector<double> &Other_Variables,
									
									Vector<double> &Variable_Derivatives,
									double &Iion) override
		{
			//Lotka-Volterra function of population growth
			// Variable_Derivatives[x] = alpha*CellVariables[x] - beta*CellVariables[x]*CellVariables[y];
			// Variable_Derivatives[y] = delta*CellVariables[x]*CellVariables[y] - gamma*CellVariables[y];

			Variable_Derivatives[x] = CellVariables[y];
			// Variable_Derivatives[x] = sin(t);
			Variable_Derivatives[y] = -CellVariables[x];


			// //Fitzhugh Nagumo Model
			Variable_Derivatives[w] = (1.0/tau)*(Vm + a - b*CellVariables[w]);
			Iion = (Vm - Vm*Vm*Vm/3.0 - CellVariables[w] + Istim);
			// Variable_Derivatives[w] = 0.0;
			// Iion = 0.0;
		}

		void get_optional_output(const double &Vm,
							const Vector<double> &CellVariables,
							const double &t,
							const unsigned &cell_type,
							const double &Istim,
							const Vector<double> &Other_Parameters,
							const Vector<double> &Other_Variables,

							Vector<double> &Out) const override
		{

		}

		

	protected:

		//Lotka-Volterra function of population model
		double alpha;
		double beta;
		double delta;
		double gamma;

		//Fitzhugh Nagumo Model
		double tau;
		double a;
		double b;


	};
}

#endif