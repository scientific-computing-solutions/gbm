//-----------------------------------
//
// File: varsplitter.cpp
//
// Description: class that implements the splitting of a node on a variable.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "varsplitter.h"

//---------------------
// Public Functions
//---------------------
VarSplitter::VarSplitter(unsigned long minNumObs):bestSplit(), proposedSplit()
{
	adGroupSumZ.resize(1024);
	adGroupW.resize(1024);
	acGroupN.resize(1024);
	adGroupMean.resize(1024);
	aiCurrentCategory.resize(1024);
	aiBestCategory.resize(1024);
	cMinObsInNode = minNumObs;

	/*InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;
	fIsSplit = false;
	dLastXValue = -HUGE_VAL;
	*/
}

VarSplitter::~VarSplitter()
{
}

void VarSplitter::IncorporateObs
(
    double dX,
    double dZ,
    double dW,
    long lMonotone
)
{

	static double dWZ = 0.0;

	//if(fIsSplit) return;

	dWZ = dW*dZ;

	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dWZ, dW);

	}
	else if(cCurrentVarClasses == 0)   // variable is continuous
	{
		if(dLastXValue > dX)
		{
			throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
		}

		// Evaluate the current split
		// the newest observation is still in the right child
		proposedSplit.SplitValue = 0.5*(dLastXValue + dX);

		if((dLastXValue != dX) &&
			proposedSplit.HasMinNumOfObs(cMinObsInNode) &&
			proposedSplit.SplitIsCorrMonotonic(lMonotone))
		{
			proposedSplit.NodeGradResiduals();
			if(proposedSplit.HasMinNumOfObs(cMinObsInNode) &&
								(proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
			{
				bestSplit = proposedSplit;
				//std::cout << "Current Best Var: " << iBestSplitVar << " " << "best split var: " << bestSplit.SplitVar << endl;
			}
		}

		// now move the new observation to the left
		// if another observation arrives we will evaluate this
		proposedSplit.UpdateLeftNode(dWZ, dW);
		dLastXValue = dX;
	}
	else // variable is categorical, evaluates later
	{
		adGroupSumZ[(unsigned long)dX] += dWZ;
		adGroupW[(unsigned long)dX] += dW;
		acGroupN[(unsigned long)dX] ++;
	}
	/*if(fIsSplit) return;
	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dW*dZ, dW);
	}
	else if(proposedSplit.SplitClass == 0)
	{
		if(dLastXValue > dX)
		{
			throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
		}

		// Evaluate the current split
		// the newest observation is still in the right child
		proposedSplit.SplitValue = 0.5*(dLastXValue + dX);

		if((dLastXValue != dX) &&
			proposedSplit.HasMinNumOfObs(minObsInNode) &&
			proposedSplit.SplitIsCorrMonotonic(lMonotone))
		{
			proposedSplit.NodeGradResiduals();

			if(proposedSplit.HasMinNumOfObs(minObsInNode) &&
					(proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
			{
				bestSplit = proposedSplit;
				WrapUpSplit();

			}

		}

		// now move the new observation to the left
		// if another observation arrives we will evaluate this
		proposedSplit.UpdateLeftNode(dW*dZ, dW);
		dLastXValue = dX;
	}
	else // variable is categorical, evaluates later
	{
		proposedSplit.IncrementCategories((unsigned long) dX, dW*dZ, dW);
	}*/
}


void VarSplitter::EvaluateCategoricalSplit()
{
	 long i=0;
	  unsigned long cFiniteMeans = 0;

	  //if(fIsSplit) return;

	  if(cCurrentVarClasses == 0)
	    {
	      throw GBM::invalid_argument();
	    }

	  cFiniteMeans = 0;
	  for(i=0; i<cCurrentVarClasses; i++)
	    {
	      aiCurrentCategory[i] = i;
	      if(adGroupW[i] != 0.0)
	        {
		  adGroupMean[i] = adGroupSumZ[i]/adGroupW[i];
		  cFiniteMeans++;
	        }
	      else
	        {
		  adGroupMean[i] = HUGE_VAL;
	        }
	    }

	  rsort_with_index(&adGroupMean[0],&aiCurrentCategory[0],cCurrentVarClasses);

	  // if only one group has a finite mean it will not consider
	  // might be all are missing so no categories enter here
	  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
	    {
	   

	      proposedSplit.SplitValue = (double) i;
	      proposedSplit.UpdateLeftNode(adGroupSumZ[aiCurrentCategory[i]], adGroupW[aiCurrentCategory[i]],
	    		  	  	  	  	  	   acGroupN[aiCurrentCategory[i]]);
	      proposedSplit.NodeGradResiduals();

		  if(proposedSplit.HasMinNumOfObs(cMinObsInNode)
		      		  && (proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
		{

		  if(bestSplit.SplitVar!= proposedSplit.SplitVar)
		  {
		      
		      cBestVarClasses = cCurrentVarClasses;
		      std::copy(aiCurrentCategory.begin(),
				aiCurrentCategory.end(),
				aiBestCategory.begin());
		  }
		  bestSplit = proposedSplit;

	        }
	    }
 /* long i=0;
  unsigned long cFiniteMeans = 0;
  if(fIsSplit) return;
  if(proposedSplit.SplitClass == 0)
	{
	  throw GBM::invalid_argument("Evaluate Categorical Split - Split variable is not categorical");
	}

  cFiniteMeans = proposedSplit.SetAndReturnNumGroupMeans();

  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
      proposedSplit.SplitValue = (double)i;
      proposedSplit.UpdateLeftNodeWithCat(i);
      proposedSplit.NodeGradResiduals();
      proposedSplit.setBestCategory();

      if(proposedSplit.HasMinNumOfObs(minObsInNode)
    		  && (proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
      {
    	  bestSplit = proposedSplit;
		  WrapUpSplit();

      }

    }*/
}

void VarSplitter::Set(CNode& nodeToSplit)
{
	dInitSumZ = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
	dInitTotalW = nodeToSplit.dTrainW;
	cInitN = nodeToSplit.cN;

	bestSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN);
	proposedSplit.ResetSplitProperties(0, dInitTotalW, cInitN);


	/*
	InitWeightResiduals = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
	InitTotalWeight = nodeToSplit.dTrainW;
	InitNumObs = nodeToSplit.cN;

	bestSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs);
	dLastXValue = -HUGE_VAL;*/

}

/*void VarSplitter::SetForVariable(unsigned long iWhichVar, long cVarClasses)
{
	if(fIsSplit) return;
	//bestSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs);
	proposedSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs,
		  proposedSplit.SplitValue,	cVarClasses, iWhichVar);



}*/

void VarSplitter::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
  //if(fIsSplit) return;

  if (int(cCurrentVarClasses) > adGroupSumZ.size()) {
    throw GBM::failure("too many variable classes");
  }

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cCurrentVarClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + cCurrentVarClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + cCurrentVarClasses, 0);
  this->cCurrentVarClasses = cCurrentVarClasses;


  proposedSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN,
  		  proposedSplit.SplitValue,	cCurrentVarClasses, iWhichVar);
  dLastXValue = -HUGE_VAL;
}

void VarSplitter::WrapUpCurrentVariable()
{
  if(proposedSplit.SplitVar == bestSplit.SplitVar)
    {
      if(proposedSplit.MissingNumObs > 0)
        {
    	  bestSplit.MissingWeightResiduals = proposedSplit.MissingWeightResiduals;
		bestSplit.MissingTotalWeight = proposedSplit.MissingTotalWeight;
		bestSplit.MissingNumObs = proposedSplit.MissingNumObs;

        }
      else // DEBUG: consider a weighted average with parent node?
        {
		bestSplit.MissingWeightResiduals   = dInitSumZ;
		bestSplit.MissingTotalWeight = dInitTotalW;
		bestSplit.MissingNumObs      = 0;
        }
    }
}

/*void VarSplitter::Reset()
{
	// Reset the splitter for new searching
	InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;

	dLastXValue = -HUGE_VAL;

	// Reset best split
	bestSplit.ResetSplitProperties(0.0, 0.0, 0);
	proposedSplit.ResetSplitProperties(0.0, 0.0, 0);

}*/

//---------------------
// Private Functions
//---------------------

/*
void VarSplitter::WrapUpSplit()
{
	if(proposedSplit.MissingNumObs <= 0)
	{
		bestSplit.MissingWeightResiduals   = InitWeightResiduals;
		bestSplit.MissingTotalWeight = InitTotalWeight;
		bestSplit.MissingNumObs      = 0;
	}
	else
	{
		bestSplit.MissingWeightResiduals = proposedSplit.MissingWeightResiduals;
		bestSplit.MissingTotalWeight = proposedSplit.MissingTotalWeight;
		bestSplit.MissingNumObs = proposedSplit.MissingNumObs;

	}
}

*/


