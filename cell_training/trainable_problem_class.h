//LIC/////////////////////////////////////////////////////////////////// 
//LIC// The sub problem class, used in training cell models.
//LIC// Sub problems are ownewd by a meta problem. A sub problem can
//LIC//  perform any simulation required but must inherit from this class.
//LIC// A sub problem is required to run it's simulation when the function
//LIC//  "run" is called. This function must return a double, generally
//LIC//  describing the discrepency between the simulation just computed
//LIC//  and the experiment it is intended to replicate.
//LIC///////////////////////////////////////////////////////////////////

#ifndef OOMPH_TRAINABLE_PROBLEM_HEADER
#define OOMPH_TRAINABLE_PROBLEM_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
	#include <oomph-lib-config.h>
#endif


#include "trainable_element.h"

#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Qelements.h"
#include "../generic/Telements.h"
#include "../generic/error_estimator.h"
#include "../generic/problem.h"

namespace oomph{

	//To be a base class of problems used in optimisation problems, allows for any trainable elements
	//	(elements which get parameters from a nelder mead optimisation class) to be handed to the
	//	nelder mead class so that those parameters can be handled
	class TrainableProblem : public Problem
	{
	public:
		TrainableProblem(){	}
		~TrainableProblem(){	}

		//Perform the simulation the problem is designed to do,
		// compare result to experimental_data and return
		// comparrison
		//This function must reset the problem to it's pre simulation
		//	state upon completion
		virtual double run(){
			throw OomphLibError("run has not been implemented for sub problem yet\nYour implementation must ensure the problem is reset completely before beign re-run.\nIt is also stressed that find_all_trainable_elements() should be called for this element before\nthe coordinator problem locates all dependent elements",
								OOMPH_CURRENT_FUNCTION,
								OOMPH_EXCEPTION_LOCATION);
		}

		//return a vector of pointers to elements within this problem with the trainable element
		void get_all_trainable_elements(Vector<TrainableElement*> &trainable_elements){
				
			//add all non-global mesh trainable elements to the vector
			for(unsigned i=0; i<Problem_Trainable_Elements.size(); i++){
				trainable_elements.push_back(Problem_Trainable_Elements[i]);
			}

			//get the trainable elements from the global mesh
			Vector<TrainableElement*> trainable_elements_from_global_mesh;
			find_all_trainable_elements_in_global_mesh(trainable_elements_from_global_mesh);
			for(unsigned i=0; i<trainable_elements_from_global_mesh.size(); i++){
				trainable_elements.push_back(trainable_elements_from_global_mesh[i]);
			}
		}

	protected:

		//If you have added a trainable element to this problem which doesn't live in the
		//	global mesh then you can make the machinery aware of it here
		void add_trainable_element(TrainableElement* candidate_element_pt){
			Problem_Trainable_Elements.push_back(candidate_element_pt);
		}

	private:

		//a helper function intended to be used within find_all_trainable_elements() to
		//	automatically get all trainable elements found inside the global mesh
		void find_all_trainable_elements_in_global_mesh(Vector<TrainableElement*> &trainable_elements_from_global_mesh){
			//get the global mesh
			// Mesh* mesh_pt = this->mesh_pt();
			//if there is no global mesh then return
			if(this->mesh_pt()==NULL){return;}
			//get the number of elements in the mesh
			unsigned n_element = this->mesh_pt()->nelement();
			//loop over them
			for(unsigned e=0; e<n_element; e++){
				//attempt to cast the element to a trainable element
				TrainableElement* candidate_element_pt = dynamic_cast<TrainableElement*>(this->mesh_pt()->element_pt(e));
				//if the cast was successfull
				if(candidate_element_pt!=NULL){
					trainable_elements_from_global_mesh.push_back(candidate_element_pt);
				}
			}
		}

		//all trainable elemnts which are not in the global mesh
		Vector<TrainableElement*> Problem_Trainable_Elements;		
	};

}//end namespace

#endif