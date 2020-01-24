#include "simulation_config.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

FileSys::FileSys() {
	/*Geometry_file = new char [100];
	Fibre_theta  = new char [100];
	Fibre_phi    = new char [100];
	Pacemap      = new char [100];
	Stim_file    = new char [100];
	SAN_type = new char[100];
	Stim_time_file = new char[100];
	Stim_amp_file = new char[100];
	Fcell_map = new char [100];
	FB_map = new char [100];
	Apicobasal_file = new char [100];*/
}
FileSys::~FileSys() {

	/* delete [] Geometry_file;
	 delete [] Fibre_theta ;
	 delete [] Fibre_phi   ;
	 delete [] Pacemap     ;
	 delete [] Stim_file   ;
	 delete [] SAN_type;
	 delete [] Stim_time_file;
	 delete [] Stim_amp_file;
	 delete [] Fcell_map;
	 delete [] FB_map;
	 delete [] Apicobasal_file;*/
}


Simulation_Config::Simulation_Config() {
	Initilise();
}

void Simulation_Config::Initilise() {

	BCL        = 700;
	S1_number  = 10;
	S2         = BCL;
	dt         = 0.005;

	Total_time = 2000.0;
	Model_type = 2;
	IKur_type  = 0;
	model_out  = 0;
	region     = 3;
	region_3D  = 202;
	AF_model   = 0;
	mutation   = 0;
	FB_type    = 0;
	FB_number  = 0;
	Ggap       = 3.0;
	tau_type   = 1;
	Diff_Scale = 1.0;

	region_char     = "RAA";
	Model_type_char = "Colman";
	tau_type_char   = "Fast";
	ICs             = "function_defined";
	Pacemap         = "wholeHeart_091148011-10.pace.gz";
	Stim_type       = "Paced";
	SAN_type        = "Heterogeneous";
	mutation_char   = "None";
	ICs             = "Default"; // Vent_seeman13.theta.gz
	Stim_time_file  = "Geometry/Vent_seeman.time_2.bin";
	Stim_amp_file   = "Geometry/Vent_seeman.volt_2.bin";
	Geometry_file   =  "Geometry/Vent_seeman_Add_RV_Stim_RG.mat.gz";
	Fibre_theta     = "Geometry/Vent_seeman13.theta.gz";
	Fibre_phi       = "Geometry/Vent_seeman13.phi.gz";
	Apicobasal_file = "Geometry/Vent_seeman_Apicobasal_file.bin";
	RVIndex_file    = "Geometry/RV_Gradient_Index.bin";
}
Simulation_Config::Simulation_Config(int argc, char *argv[]) {
	Initilise();
	// Simulation_Config();
	Config_handling(argc, argv);

}
Simulation_Config::~Simulation_Config() {
	/*delete [] Model_type_char;
	delete [] region_char    ;
	delete [] tau_type_char  ;
	delete [] mutation_char  ;
	delete [] ICs            ;
	delete [] Stim_type      ;*/
}


void Simulation_Config::Config_handling(int argc, char *argv[]) {

	int arg_counter;

	arg_counter = 1;
	while (arg_counter < argc) {
		std::string command = argv[arg_counter];
		if (command == "BCL") {
			BCL = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Total_time") {
			Total_time = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Model_type") {
			Model_type_char = argv[arg_counter + 1];
			if (Model_type_char == "Colman") {
				IKur_type = 0;
				model_out = 0;
			}
			else if (Model_type_char == "CNZ") {
				IKur_type = 1;
				model_out = 1;
			}
			else if (Model_type_char == "Colman_v") Model_type = 4;
			else {
				std::cerr << "Invalid Model type\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		} // end Model_type if
		else if (command == "Region") {
			region_char = argv[arg_counter + 1];
			if (region_char == "PM") {
				region = 1;
				region_3D = 12;
			}
			else if (region_char == "CT") {
				region = 2;
				region_3D = 11;
			}
			else if (region_char == "RAA") {
				region = 3;
				region_3D = 202;
			}
			else if (region_char == "AVR") {
				region = 4;
				region_3D = 18;
			}
			else if (region_char == "BB") {
				region = 5;
				region_3D = 15;
			}
			else if (region_char == "LA") {
				region = 6;
				region_3D = 16;
			}
			else if (region_char == "AS") {
				region = 7;
				region_3D = 17;
			}
			else if (region_char == "LAA") {
				region = 8;
				region_3D = 201;
			}
			else if (region_char == "PV") {
				region = 11;
				region_3D = 101;
			}
			else if (region_char == "PV_jones") region = 9;
			else if (region_char == "SAN_C") {
				region = 14;
				region_3D = 10;
			}
			else if (region_char == "SAN_P") {
				region = 15;
				region_3D = 1011;
			}
			arg_counter++;
		}
		else if (command == "AF") {
			AF_model = atoi(argv[arg_counter + 1]);
			if (AF_model > 4) {
				std::cerr << "AF model must be 0, 1, 2, 3 or 4\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}
		else if (command == "Mutation") {
			mutation_char = argv[arg_counter + 1];

			if (mutation_char == "D322H") mutation = 1;
			else if (mutation_char == "E48G") mutation = 2;
			else if (mutation_char == "A305T") mutation = 3;
			else if (mutation_char == "Y155C") mutation = 4;
			else if (mutation_char == "D469E") mutation = 5;
			else if (mutation_char == "P488S") mutation = 6;
			else if (mutation_char == "None") mutation  = 0;
			else if (mutation_char == "A545P") mutation  = 10;
			else {
				std::cerr << "Invalid mutation\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}
		else if (command == "FB_type") {
			FB_type = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "FB_number") {
			FB_number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Ggap") {
			Ggap = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Tau_type") {
			tau_type_char = argv[arg_counter + 1];
			if (tau_type_char ==  "Slow") tau_type = 0;
			else if (tau_type_char ==  "Fast") tau_type = 1;
			else {
				std::cerr << "invalid tau type";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}
		else if (command == "BCL") {
			BCL = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "S1") {
			BCL = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "S2") {
			S2 = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "SAN_type") {
			SAN_type = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Stim_type") {
			Stim_type = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "S1_number") {
			S1_number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Diffusion_Scale") {
			Diff_Scale = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Pacemap") {
			Pacemap = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Stim_Time_File") {
			Stim_time_file = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Stim_Amp_File") {
			Stim_amp_file = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Apicobasal_File") {
			Apicobasal_file = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "RVIndex_File") {
			RVIndex_file = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Geometry_File") {
			Geometry_file = argv[arg_counter + 1];
			arg_counter++;
		}
		else {
			std::cerr << argv[arg_counter] << " is not a valid argument. Valid arguments are as follows:\n\t "
			          << "BCL\n\tTotal_time\n\tModel_type\n\tRegion\n\tAF\n\tMutation\n\tFB_type\n\tFB_number\n\tGgap\n\tTau_type\n";
			std::exit(EXIT_FAILURE);
		}
		arg_counter ++;
	} // end while

} // end Argument_handling

void Simulation_Config::Report_Config() {

	std::cerr << "BCL\t\t = " << BCL << std::endl;
	std::cerr << "Total_time\t = " << Total_time << std::endl;
	std::cerr << "dt\t\t = " << dt << std::endl;
	std::cerr << "Model_type\t = " << Model_type_char << std::endl;
	std::cerr << "Region\t\t = " << region << std::endl;
	std::cerr << "Region_char\t = " << region_char << std::endl;
	std::cerr << "AF_model\t = " << AF_model << std::endl;
	std::cerr << "tau_type_char\t = " << tau_type_char << std::endl;
	std::cerr << "Diff_Scale\t = " << Diff_Scale << std::endl;
	if (FB_number != 0)  {
		std::cerr << "Fibroblasts on:\n\t"
		          << "FB_type = " << FB_type
		          << "\n\tFB_number = " << FB_number
		          << "\n\tGgap = " << Ggap << std::endl;
	}

	std::cerr << "Pacemap\t\t = " << Pacemap << std::endl
	          << "Geometry_file\t = " << Geometry_file << std::endl
	          << "Apicobasal_file\t = " << Apicobasal_file << std::endl
	          << "Fibre_theta\t = "  << Fibre_theta << std::endl
	          << "Fibre_phi\t = "  << Fibre_phi << std::endl
	          <<  "Stim_type\t= " << Stim_type << std::endl
	          <<  "Stim_Time_File\t= " << Stim_time_file << std::endl
	          <<  "Stim_Amp_File\t= " << Stim_amp_file << std::endl;
}


void Simulation_Config::Report_All() {
	std::cerr << "BCL\t\t = " << BCL << std::endl;
	std::cerr << "Total_time\t = " << Total_time << std::endl;
	std::cerr << "dt\t\t = " << dt << std::endl;
	std::cerr << "Model_type\t = " << Model_type_char << std::endl;
	std::cerr << "Region\t\t = " << region << std::endl;
	std::cerr << "Region_char\t = " << region_char << std::endl;
	std::cerr << "AF_model\t = " << AF_model << std::endl;
	std::cerr << "tau_type_char\t = " << tau_type_char << std::endl;
	std::cerr << "Diff_Scale\t = " << Diff_Scale << std::endl;
}
