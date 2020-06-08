//LIC/////////////////////////////////////////////////////////////////// 
//LIC// The sub problem class, used in training cell models.
//LIC// Sub problems are ownewd by a meta problem. A sub problem can
//LIC//  perform any simulation required but must inherit from this class.
//LIC// A sub problem is required to run it's simulation when the function
//LIC//  "run" is called. This function must return a double, generally
//LIC//  describing the discrepency between the simulation just computed
//LIC//  and the experiment it is intended to replicate.
//LIC///////////////////////////////////////////////////////////////////

#ifndef OOMPH_SUB_PROBLEM_HEADER
#define OOMPH_SUB_PROBLEM_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
	#include <oomph-lib-config.h>
#endif

#include "../generic/problem.h"

#include "../cell_model/cell_model.h"

namespace oomph{
	class SubProblem	:	public Problem
	{
	public:
		SubProblem(){	}
		~SubProblem(){	}

		//Perform the simulation the problem is designed to do,
		// compare result to experimental_data and return
		// comparrison
		double run();

		void assign_cell_model_pt(CellModelBase* new_cell_model_pt){
			Cell_Model_Pt = new_cell_model_pt;
		}

	protected:

		CellModelBase* Cell_Model_Pt;

	private:

	};

}//end namespace

#endif