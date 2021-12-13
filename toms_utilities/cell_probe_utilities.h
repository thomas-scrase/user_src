#ifndef OOMPH_CELL_PROBE_UTILITIES_HEADER
#define OOMPH_CELL_PROBE_UTILITIES_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
	#include <oomph-lib-config.h>
#endif

#include "../cell_models/cell_model_base.h"



namespace oomph{

namespace CELLPROBE_PARAMS{
	//dvdt above which we assume we are in an upstroke, mV/ms
	const double Upstroke_dvdt = 100.0;

	//Magic number indicating that a duration was not measured
	const double Duration_Not_Measured = -1.0;

	//Magic number indicating that a minimum measurement (e.g. vmin) was not measured
	const double Value_Not_Measured = 1e9;

	const unsigned Maximum_measurements = 10000;

	const double T_After_Upstroke_Start_Measure_Plateau = 10.0;

	const double T_After_Upstroke_Stop_Measure_Plateau = 50.0;
}

class APDMeasurer
{
protected:

	//Open the file we write apd information to
	void OpenOutFile(const std::string OutFileName)
	{
		if(!OutFileName.empty())
		{
			OutFile.open(OutFileName);
			OutputFileHeaderInformation();
		}
	}

	void InitializeMeasurements(){
		Num_Upstrokes = 0;

		Post_Upstroke = false;
		Within_Upstroke = false;

		Candidate_dvdt_max.clear();
		Candidate_t_dvdt_max.clear();
		Candidate_vmax.clear();

		
		t_upstroke.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		dvdt_upstroke.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);

		

		vmax.resize(CELLPROBE_PARAMS::Maximum_measurements, -CELLPROBE_PARAMS::Value_Not_Measured);
		vmin.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);

		resting_potential.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);
		v_plateau.resize(CELLPROBE_PARAMS::Maximum_measurements, 0.0);//because this is measured through repeated addition we have to set it to zero as a default value
		
		APD20.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		APD50.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		APD75.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		APD90.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);

		t_APD20.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		t_APD50.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		t_APD75.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		t_APD90.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);

		CompletedAPD20 = false;
		CompletedAPD50 = false;
		CompletedAPD75 = false;
		CompletedAPD90 = false;


		v20.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);
		v50.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);
		v75.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);
		v90.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Value_Not_Measured);

		DI20.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		DI50.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		DI75.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);
		DI90.resize(CELLPROBE_PARAMS::Maximum_measurements, CELLPROBE_PARAMS::Duration_Not_Measured);


	}

	void OutputFileHeaderInformation(){
		if(!OutFile.is_open())
		{
			OomphLibWarning(("A cell probe did not open an output file but is attempting to output."),
				"CellProbe::OutputFileHeaderInformation()",
				OOMPH_EXCEPTION_LOCATION);
		}
		OutFile
		<< "Num_Upstrokes "
		<< "t_upstroke "
		<< "dvdt_upstroke "
		<< "vmax "
		<< "v_plateau "

		<< "APD20 "
		<< "APD50 "
		<< "APD75 "
		<< "APD90 "

		<< "t_APD20 "
		<< "t_APD50 "
		<< "t_APD75 "
		<< "t_APD90 "

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

	void OutputData(){
		if(!OutFile.is_open())
		{
			OomphLibWarning(("A cell probe did not open an output file but is attempting to output."),
				"CellProbe::OutputData()",
				OOMPH_EXCEPTION_LOCATION);
		}
		//Output the data
		//////////////////////////
		//////////////////////////
		//////////////////////////
		OutFile
		<< Num_Upstrokes << " "
		<< t_upstroke[Num_Upstrokes] << " "
		<< dvdt_upstroke[Num_Upstrokes] << " "
		<< vmax[Num_Upstrokes] << " "
		<< v_plateau[Num_Upstrokes] << " "

		<< APD20[Num_Upstrokes] << " "
		<< APD50[Num_Upstrokes] << " "
		<< APD75[Num_Upstrokes] << " "
		<< APD90[Num_Upstrokes] << " "

		<< t_APD20[Num_Upstrokes] << " "
		<< t_APD50[Num_Upstrokes] << " "
		<< t_APD75[Num_Upstrokes] << " "
		<< t_APD90[Num_Upstrokes] << " "

		<< v20[Num_Upstrokes] << " "
		<< v50[Num_Upstrokes] << " "
		<< v75[Num_Upstrokes] << " "
		<< v90[Num_Upstrokes] << " "

		<< (v90[Num_Upstrokes] - v50[Num_Upstrokes])/(APD90[Num_Upstrokes]-APD50[Num_Upstrokes]) << " "

		<< DI20[Num_Upstrokes] << " "
		<< DI50[Num_Upstrokes] << " "
		<< DI75[Num_Upstrokes] << " "
		<< DI90[Num_Upstrokes] << " "

		<< vmin[Num_Upstrokes] << " "
		<< resting_potential[Num_Upstrokes] << std::endl;
	}



	void MeasureAPD(const double& t, const double& v, const double& dvdt)
	{	
		//If we are encountering an upstroke
		if(dvdt>CELLPROBE_PARAMS::Upstroke_dvdt)
		{
			Within_Upstroke = true;
		}
		if(Within_Upstroke)
		{
			//Record the time and dvdt so we can calculate the correct time of dvdt_max
			Candidate_dvdt_max.push_back(dvdt);
			Candidate_t_dvdt_max.push_back(t);
			Candidate_vmax.push_back(v);

			if(dvdt<0.0)//if the membrane potential starts to drop again we are no longer in an upstroke
			{
				Within_Upstroke = false;
				Post_Upstroke = true;
			}
		}

		//if we are not within an upstroke
		if(!Within_Upstroke)
		{
			//and if we have just encountered one
			if(Post_Upstroke == true)
			{
				Num_Upstrokes++;

				//Calculate the true time of upstroke and the maximum dvdt
				unsigned index_of_max = 0;
				for(unsigned ind = 0; ind < Candidate_dvdt_max.size(); ind++)
				{
					if(Candidate_dvdt_max[ind]>Candidate_dvdt_max[index_of_max])
					{
						index_of_max = ind;
					}

					if(Candidate_vmax[ind]>vmax[Num_Upstrokes])
					{
						vmax[Num_Upstrokes] = Candidate_vmax[ind];
					}
				}
				dvdt_upstroke[Num_Upstrokes] = Candidate_dvdt_max[index_of_max];
				t_upstroke[Num_Upstrokes] = Candidate_t_dvdt_max[index_of_max];
				

				//Calculate the diastolic intervals, these are of the next set of measurements since they
				// are what illicit the next measurements
				DI20[Num_Upstrokes] = t_upstroke[Num_Upstrokes] - t_APD20[Num_Upstrokes-1];
				DI50[Num_Upstrokes] = t_upstroke[Num_Upstrokes] - t_APD50[Num_Upstrokes-1];
				DI75[Num_Upstrokes] = t_upstroke[Num_Upstrokes] - t_APD75[Num_Upstrokes-1];
				DI90[Num_Upstrokes] = t_upstroke[Num_Upstrokes] - t_APD90[Num_Upstrokes-1];

				//Determine what the membrane potential values at which APDs are measured for the next AP
				v20[Num_Upstrokes] = vmax[Num_Upstrokes] - 0.20 * (vmax[Num_Upstrokes] - resting_potential[Num_Upstrokes-1]);
				v50[Num_Upstrokes] = vmax[Num_Upstrokes] - 0.50 * (vmax[Num_Upstrokes] - resting_potential[Num_Upstrokes-1]);
				v75[Num_Upstrokes] = vmax[Num_Upstrokes] - 0.75 * (vmax[Num_Upstrokes] - resting_potential[Num_Upstrokes-1]);
				v90[Num_Upstrokes] = vmax[Num_Upstrokes] - 0.90 * (vmax[Num_Upstrokes] - resting_potential[Num_Upstrokes-1]);

				//Finish calculating the plateau potential of this APD
				v_plateau[Num_Upstrokes-1]/=(CELLPROBE_PARAMS::T_After_Upstroke_Stop_Measure_Plateau - CELLPROBE_PARAMS::T_After_Upstroke_Start_Measure_Plateau);


				//Clear the candidate data
				Candidate_dvdt_max.clear();
				Candidate_t_dvdt_max.clear();
				Candidate_vmax.clear();


				

				//We have encountered the next upstroke
				Post_Upstroke = false;

				//We need to remember to record the action potential durations for this upstroke
				CompletedAPD20=false;
				CompletedAPD50=false;
				CompletedAPD75=false;
				CompletedAPD90=false;
			}
			else
			{
				//Calculate plateau potential
				if(t>t_upstroke[Num_Upstrokes]+CELLPROBE_PARAMS::T_After_Upstroke_Start_Measure_Plateau && t<t_upstroke[Num_Upstrokes]+CELLPROBE_PARAMS::T_After_Upstroke_Stop_Measure_Plateau){
					v_plateau[Num_Upstrokes] += (t-tprev)*(v+vprev)*0.5;
				}


				if (!CompletedAPD20)
				{
					if((vprev >= v20[Num_Upstrokes]) && (v <= v20[Num_Upstrokes]) )
					{
			       		APD20[Num_Upstrokes] = t - t_upstroke[Num_Upstrokes];
			       		t_APD20[Num_Upstrokes] = t;
			       		CompletedAPD20 = true;
	       			}
		        }
		        if (!CompletedAPD50)
	        	{
	        		if((vprev >= v50[Num_Upstrokes]) && (v <= v50[Num_Upstrokes] ))
	        		{
			            APD50[Num_Upstrokes] = t - t_upstroke[Num_Upstrokes];
			            t_APD50[Num_Upstrokes] = t;
			            CompletedAPD50 = true;
	        		}
		        }
		        if (!CompletedAPD75)
		        {
	        		if(vprev >= v75[Num_Upstrokes] && v <= v75[Num_Upstrokes] )
	        		{
			            APD75[Num_Upstrokes] = t - t_upstroke[Num_Upstrokes];
			            t_APD75[Num_Upstrokes] = t;
			            CompletedAPD75 = true;
		        	}
		        }
		        if (!CompletedAPD90)
		        {
		        	if(vprev >= v90[Num_Upstrokes] && v <= v90[Num_Upstrokes] )
		        	{
			            APD90[Num_Upstrokes] = t - t_upstroke[Num_Upstrokes];
			            t_APD90[Num_Upstrokes] = t;
			            CompletedAPD90 = true;
		        	}
		        }

		        //If we are post APD90, measure the resting potential as the minimum potential until the next upstroke
				// if(t_APD90[Num_Upstrokes]>=0.0){
					if(v<resting_potential[Num_Upstrokes]){
						resting_potential[Num_Upstrokes] = v;
					}
				// }

				//Check for Vmin, change if necessary
				if(v<vmin[Num_Upstrokes]){
					vmin[Num_Upstrokes] = v;
				}
			}

			vprev = v;
			tprev = t;
		}

	}


	//Get the values from the last completed APD cycle. Since some measurements can only be made certain at an upstroke we
	// define an APD from one upstroke to the next. Therefore, we index all of these values at Num_Upstrokes-1

	//Get the time which the previous apd90 was measured
	double get_tAPD90()
	{
		return t_APD90[Num_Upstrokes];
	}



private:
	//The file we output data to
	std::ofstream OutFile;

	//The number of times the cell has experienced an upstroke
	unsigned Num_Upstrokes;

	bool Post_Upstroke;
	bool Within_Upstroke;

	Vector<double> Candidate_dvdt_max;
	Vector<double> Candidate_t_dvdt_max;
	Vector<double> Candidate_vmax;
	Vector<double> vmax;
	Vector<double> vmin;

	Vector<double> t_upstroke;
	Vector<double> dvdt_upstroke;	

	Vector<double> resting_potential;

	Vector<double> v_plateau;

	Vector<double> APD20;
	Vector<double> APD50;
	Vector<double> APD75;
	Vector<double> APD90;

	Vector<double> t_APD20;
	Vector<double> t_APD50;
	Vector<double> t_APD75;
	Vector<double> t_APD90;

	bool CompletedAPD20;
	bool CompletedAPD50;
	bool CompletedAPD75;
	bool CompletedAPD90;

	Vector<double> v20;
	Vector<double> v50;
	Vector<double> v75;
	Vector<double> v90;

	Vector<double> DI20;
	Vector<double> DI50;
	Vector<double> DI75;
	Vector<double> DI90;

	double tprev;
	double vprev;

//Outputting functions
public:
	void get_APD90_and_DI90(double& apd90, double& di90)
	{
		apd90 = APD90[Num_Upstrokes];
		di90 = DI90[Num_Upstrokes];
	}

	//Get the time of the upstroke of the last complete AP the cell undertook, note that since the number of upstrokes
	// is incremented when an upstroke is encountered, taking the latest upstroke encountered will always return a
	// default (nonsensical) value
	double get_t_upstroke()
	{
		return t_upstroke[Num_Upstrokes];
	}
	void get_t_upstroke(double& t)
	{
		t = t_upstroke[Num_Upstrokes];
	}

	//Get the time of the upstroke of the last complete AP the cell undertook, note that since the number of upstrokes
	// is incremented when an upstroke is encountered, taking the latest upstroke encountered will always return a
	// default (nonsensical) value
	double get_latest_t_upstroke()
	{
		return t_upstroke[Num_Upstrokes-1];
	}
	void get_latest_t_upstroke(double& t)
	{
		t = t_upstroke[Num_Upstrokes-1];
	}


	//Get the latest reliable measurements of apd90 and di90
	void get_latest_APD90_and_DI90(double& apd90, double& di90)
	{
		apd90 = APD90[Num_Upstrokes-1];
		di90 = DI90[Num_Upstrokes-1];
	}

	void get_latest_resting_potential(double& rp)
	{
		rp = resting_potential[Num_Upstrokes-1];
	}

	double get_latest_resting_potential()
	{
		return resting_potential[Num_Upstrokes-1];
	}

	void get_latest_maximum_potential(double& max_vm)
	{
		max_vm = vmax[Num_Upstrokes-1]; //This records the maximum potential of the previous completed cycle, not the last upstroke encountered.
	}

	double get_latest_upstroke_velocity()
	{
		return dvdt_upstroke[Num_Upstrokes-1];
	}




	unsigned get_n_upstroke()
	{
		return Num_Upstrokes;
	}


	void dump_measurements_to_file()
	{
		if(!OutFile.is_open())
		{
			OomphLibWarning(("A cell probe did not open an output file but is attempting to output."),
				"CellProbe::dump_measurements_to_file()",
				OOMPH_EXCEPTION_LOCATION);
		}

		for(unsigned i=0; i<Num_Upstrokes+1; i++)
		{
			//Output the data
			//////////////////////////
			//////////////////////////
			//////////////////////////
			OutFile
			<< i << " "
			<< t_upstroke[i] << " "
			<< dvdt_upstroke[i] << " "
			<< vmax[i] << " "
			<< v_plateau[i] << " "

			<< APD20[i] << " "
			<< APD50[i] << " "
			<< APD75[i] << " "
			<< APD90[i] << " "

			<< t_APD20[i] << " "
			<< t_APD50[i] << " "
			<< t_APD75[i] << " "
			<< t_APD90[i] << " "

			<< v20[i] << " "
			<< v50[i] << " "
			<< v75[i] << " "
			<< v90[i] << " "

			<< (v90[i] - v50[i])/(APD90[i]-APD50[i]) << " "

			<< DI20[i] << " "
			<< DI50[i] << " "
			<< DI75[i] << " "
			<< DI90[i] << " "

			<< vmin[i] << " "
			<< resting_potential[i] << std::endl;
		}
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

	CellProbe(CellModelBaseFullySegregated* cell_pt, const std::string OutFileName = "", const unsigned& l = 0)
	{
		Cell_pt = cell_pt;

		this->OpenOutFile(OutFileName);

		//Grab the current value of membrane potential as vmin
		InitializeMeasurements();

		Record();
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