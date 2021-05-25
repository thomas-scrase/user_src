#ifndef OOMPH_VALIDATE_LOTKA_VOLTERRA_N_SPECIES_HEADER
#define OOMPH_VALIDATE_LOTKA_VOLTERRA_N_SPECIES_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "cell_model_base_updated.h"
	
namespace oomph{

	template<unsigned N_SPECIES>
	class ValidateCellLotkaVolterra : public CellModelBaseUpdated
	{
	public:
		ValidateCellLotkaVolterra()
		{
			//Lotka-Volterra function of population model
			r_i.resize(N_SPECIES);
			alpha_i.resize(N_SPECIES,N_SPECIES);


			double Max_r = 5.0;
			double Max_alpha = 2.0;
			std::srand(pow(N_SPECIES,4)+N_SPECIES+1);
			for(unsigned i=0; i<N_SPECIES; i++){
				r_i[i] = ((double)std::rand()/(RAND_MAX))*Max_r;
				std::cout << "r["<<i<<"] "<<r_i[i]<<std::endl;
				for(unsigned j=0; j<N_SPECIES; j++){
					alpha_i(i,j) = ((double)std::rand()/(RAND_MAX))*Max_alpha;
					std::cout << "a["<<i<<","<<j<<"] "<<alpha_i(i,j)<<std::endl;
				}
				alpha_i(i,i) = 1.0;
				// std::cout << "a["<<i<<","<<i<<"] "<<alpha_i(i,i)<<std::endl;
			}


			//Assign the names of the variables used by this model
			Names_Of_Cell_Variables.resize(N_SPECIES);
			for(unsigned i=0; i<N_SPECIES; i++){
				Names_Of_Cell_Variables[i] = ("x" + std::to_string(i));
			}

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
			std::srand(pow(N_SPECIES,3)+N_SPECIES+v);
			return ((double)std::rand()/(RAND_MAX));
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
			for(unsigned i=0; i<N_SPECIES; i++){
				Variable_Derivatives[i] = r_i[i]*CellVariables[i];
				for(unsigned j=0; j<N_SPECIES; j++){
					Variable_Derivatives[i] -= r_i[i]*CellVariables[i]*alpha_i(i,j)*CellVariables[j];
				}
			}

			//We don't use the membrane potential for this validation
			Iion = 0.0;
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

		Vector<double> r_i;
		DenseMatrix<double> alpha_i;

	};
}

#endif