#ifndef OOMPH_CELL_PROBE_UTILITIES_HEADER
#define OOMPH_CELL_PROBE_UTILITIES_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
	#include <oomph-lib-config.h>
#endif

#include "../cell_models/cell_model_base.h"



namespace oomph{

namespace CELLPROBE_PARAMS{
	//dvdt above which we assume we are in an upstroke
	double Upstroke_dvdt = 50.0;

	//Magic number indicating that a duration was not measured
	double Duration_Not_Measured = -1.0;

	//Magic number indicating that a minimum measurement (e.g. vmin) was not measured
	double Min_Not_Measured = 1e9;
}

class APDMeasurer
{
protected:

	//Open the file we write apd information to
	void OpenOutFile(std::string OutFileName)
	{
		OutFile.open(OutFileName);

		OutputFileHeaderInformation();
	}

	void InitializeMeasurements(){
		Num_Upstrokes = 0;

		Post_Upstroke = false;

		dvdt_upstroke = 0.0;
		t_upstroke = CELLPROBE_PARAMS::Duration_Not_Measured;

		Candidate_dvdt_max.clear();
		Candidate_t_dvdt_max.clear();
		
		APD20 = CELLPROBE_PARAMS::Duration_Not_Measured;
		APD50 = CELLPROBE_PARAMS::Duration_Not_Measured;
		APD75 = CELLPROBE_PARAMS::Duration_Not_Measured;
		APD90 = CELLPROBE_PARAMS::Duration_Not_Measured;
		v20 = 0.0;
		v50 = 0.0;
		v75 = 0.0;
		v90 = 0.0;

		DI20 = CELLPROBE_PARAMS::Duration_Not_Measured;
		DI50 = CELLPROBE_PARAMS::Duration_Not_Measured;
		DI75 = CELLPROBE_PARAMS::Duration_Not_Measured;
		DI90 = CELLPROBE_PARAMS::Duration_Not_Measured;

		vmin = CELLPROBE_PARAMS::Min_Not_Measured;

		resting_potential = vmin;
	}

	void OutputFileHeaderInformation(){
		OutFile
		<< "Num_Upstrokes "
		<< "t_upstroke "
		<< "dvdt_upstroke "
		<< "vmax_current "
		<< "v_plateau "

		<< "APD20 "
		<< "APD50 "
		<< "APD75 "
		<< "APD90 "

		<< "v20 "
		<< "v50 "
		<< "v75 "
		<< "v90 "

		<< "Repolarization_Rate_Between_Between_APD50_and_APD90 "

		<< "DI20 "
		<< "DI50 "
		<< "DI75 "
		<< "DI90 "

		<< "vmin "
		<< "resting_potential"
		<< std::endl;
	}

	void MeasureAPD(const double& t, const double& v, const double& dvdt)
	{
		if(dvdt>CELLPROBE_PARAMS::Upstroke_dvdt){
			// if(Post_Upstroke==false){
				
			// }
			//Post upstroke
			Post_Upstroke = true;

			//Record the time and dvdt so we can calculate the correct time of dvdt_max
			Candidate_dvdt_max.push_back(dvdt);
			Candidate_t_dvdt_max.push_back(t);
		}
		else{
			if(Post_Upstroke == true){
				//Calculate the true time of upstroke and max dvdt
				unsigned index_of_max = 0;
				for(unsigned ind = 0; ind < Candidate_dvdt_max.size(); ind++){
					if(Candidate_dvdt_max[ind]>Candidate_dvdt_max[index_of_max]){index_of_max = ind;}
				}
				dvdt_upstroke = Candidate_dvdt_max[index_of_max];
				t_upstroke = Candidate_t_dvdt_max[index_of_max];


				//Calculate durations
				DI20_next = t_upstroke - t_APD20;
				DI50_next = t_upstroke - t_APD50;
				DI75_next = t_upstroke - t_APD75;
				DI90_next = t_upstroke - t_APD90;

				//Output the data
				//////////////////////////
				//////////////////////////
				//////////////////////////
				OutFile
				<< Num_Upstrokes << " "
				<< t_upstroke << " "
				<< dvdt_upstroke << " "
				<< vmax_current << " "
				<< v_plateau/(50.0-10.0) << " "

				<< APD20 << " "
				<< APD50 << " "
				<< APD75 << " "
				<< APD90 << " "

				<< v20 << " "
				<< v50 << " "
				<< v75 << " "
				<< v90 << " "

				<< (v90 - v50)/(APD90-APD50) << " "

				<< DI20 << " "
				<< DI50 << " "
				<< DI75 << " "
				<< DI90 << " "

				<< vmin << " "
				<< resting_potential << std::endl;


				//Set the next durations to be output to be those we just measured
				// DI20 = t_upstroke - t_APD20;
				// DI50 = t_upstroke - t_APD50;
				// DI75 = t_upstroke - t_APD75;
				// DI90 = t_upstroke - t_APD90;
				DI20 = DI20_next;
				DI50 = DI50_next;
				DI75 = DI75_next;
				DI90 = DI90_next;




				//Reset the value of tmax_current, represents the maximum membrane potential achieved in this action potential
				vmax_current = 0;

				//Reset plateau potential counter and value
				v_plateau = 0.0;

				//Reset the action potential durations to a negative number,
				//	this way we know if they are not measured, for example in EAD
				// APD20 = CELLPROBE_PARAMS::Duration_Not_Measured;
				// APD50 = CELLPROBE_PARAMS::Duration_Not_Measured;
				// APD75 = CELLPROBE_PARAMS::Duration_Not_Measured;
				// APD90 = CELLPROBE_PARAMS::Duration_Not_Measured;
				APD20 = APD20_next;
				APD50 = APD50_next;
				APD75 = APD75_next;
				APD90 = APD90_next;

				//Do the same for the diastolic interval measurements
				// DI20_next = CELLPROBE_PARAMS::Duration_Not_Measured;
				// DI50_next = CELLPROBE_PARAMS::Duration_Not_Measured;
				// DI75_next = CELLPROBE_PARAMS::Duration_Not_Measured;
				// DI90_next = CELLPROBE_PARAMS::Duration_Not_Measured;

				//Calculate the APD membrane potentials
				// v20 = vmax_current - 0.20 * (vmax_current - vmin);
				// v50 = vmax_current - 0.50 * (vmax_current - vmin);
				// v75 = vmax_current - 0.75 * (vmax_current - vmin);
				// v90 = vmax_current - 0.90 * (vmax_current - vmin);

				v20 = vmax_current - 0.20 * (vmax_current - resting_potential);
				v50 = vmax_current - 0.50 * (vmax_current - resting_potential);
				v75 = vmax_current - 0.75 * (vmax_current - resting_potential);
				v90 = vmax_current - 0.90 * (vmax_current - resting_potential);

				resting_potential = 0.0;

				//Record t_upstroke as the time the APD begins
				// timeAPDstart = t_upstroke;

				//Clear the candidate data
				Candidate_dvdt_max.clear();
				Candidate_t_dvdt_max.clear();

				//We have encountered a new upstroke
				Num_Upstrokes++;

				//We have passed the upstroke so we reset the flag
				Post_Upstroke = false;
			}
			else{
				//Calculate plateau potential
				if(t>t_upstroke+10.0 && t<t_upstroke+50){
					v_plateau += (t-tprev)*(v+vprev)*0.5;
				}

				if ((vprev >= v20) && (v <= v20) ) {
	           		// APD20 = t - t_upstroke ;
	           		APD20_next = t - t_upstroke ;
	           		t_APD20 = t;
		        }
		        else if ((vprev >= v50) && (v <= v50 )) {
		            // APD50 = t - t_upstroke ;
		            APD50_next = t - t_upstroke ;
		            t_APD50 = t;
		        }
		        else if (vprev >= v75 && v <= v75 ) {
		            // APD75 = t - t_upstroke ;
		            APD75_next = t - t_upstroke ;
		            t_APD75 = t;
		        }
		        else if (vprev >= v90 && v <= v90 ) {
		            // APD90 = t - t_upstroke ;
		            APD90_next = t - t_upstroke ;
		            t_APD90 = t;
		        }
			}
		}

		//If we are post APD90, measure the resting potential as the minimum potential until the next upstroke
		if(t_APD90>=0.0){
			if(v<resting_potential){
				resting_potential = v;
			}
		}

		//Check for maximum v, change if necessary
		if(v>vmax_current){
			vmax_current = v;

			// //Calculate the APD membrane potentials
			// v20 = vmax_current - 0.20 * (vmax_current - vmin);
			// v50 = vmax_current - 0.50 * (vmax_current - vmin);
			// v75 = vmax_current - 0.75 * (vmax_current - vmin);
			// v90 = vmax_current - 0.90 * (vmax_current - vmin);
		}

		//Check for Vmin, change if necessary
		if(v<vmin){
			vmin = v;

			// //Calculate the APD membrane potentials
			// v20 = vmax_current - 0.20 * (vmax_current - vmin);
			// v50 = vmax_current - 0.50 * (vmax_current - vmin);
			// v75 = vmax_current - 0.75 * (vmax_current - vmin);
			// v90 = vmax_current - 0.90 * (vmax_current - vmin);
		}


		vprev = v;
		tprev = t;
	}

	//Get the time which the previous apd90 was measured
	double get_tAPD90()
	{
		return t_APD90;
	}



private:
	//The file we output data to
	std::ofstream OutFile;

	//The number of times the cell has experienced an upstroke
	unsigned Num_Upstrokes;

	bool Post_Upstroke;
	double t_upstroke;
	double dvdt_upstroke;
	// double timeAPDstart;

	Vector<double> Candidate_dvdt_max;
	Vector<double> Candidate_t_dvdt_max;
	double vmax_current;
	double vmin;

	double resting_potential;

	double v_plateau;

	double APD20;
	double APD50;
	double APD75;
	double APD90;

	double APD20_next;
	double APD50_next;
	double APD75_next;
	double APD90_next;

	double t_APD20;
	double t_APD50;
	double t_APD75;
	double t_APD90;

	double v20;
	double v50;
	double v75;
	double v90;

	double DI20;
	double DI50;
	double DI75;
	double DI90;
	
	double DI20_next;
	double DI50_next;
	double DI75_next;
	double DI90_next;

	double tprev;
	double vprev;

//Outputting functions
public:
	void get_APD90_and_DI90(double& apd90, double& di90)
	{
		apd90 = APD90;
		di90 = DI90;
	}
};



//Specific implementations of the cell probe:
//It's too difficult to implement one single class which works for both,
//I also think it's useful for a user to be sure which probe type they are using
//	because of the argument differences in the Record(...) function

//For Cell Interface calculated cells
class CellProbe : public virtual APDMeasurer
{
public:

	CellProbe(CellModelBaseFullySegregated* cell_pt, std::string OutFileName, const unsigned& l = 0)
	{
		Cell_pt = cell_pt;

		this->OpenOutFile(OutFileName);

		//Grab the current value of membrane potential as vmin
		InitializeMeasurements();
	}

	//Perform recordings if we are an oomph lib element, grab the data straight from the element
	void Record()
	{
		const double t = Cell_pt->get_time();
		const double V = Cell_pt->get_predicted_vm();
		const double dVdt = Cell_pt->get_predicted_dvmdt();

		this->MeasureAPD(t, V, dVdt);
	}

	//Reset the probe, useful if we are using a single probe for measuring multiple different conditions
	// e.g. applying varying dose of a drug to a cell model.
	void Reset_Running_Values()
	{
		InitializeMeasurements();
	}


	double get_tAPD90()
	{
		return APDMeasurer::get_tAPD90();
	}

private:

	//The cell element
	CellModelBaseFullySegregated* Cell_pt;
};


}

#endif