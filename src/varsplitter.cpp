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
VarSplitter::VarSplitter(unsigned long minNumObs)//:bestSplit(), proposedSplit()
{

	iBestSplitVar = 0;

	dBestSplitValue = 0.0;
	fIsSplit = false;

	dBestMissingTotalW = 0.0;
	dCurrentMissingTotalW = 0.0;
	dBestMissingSumZ = 0.0;
	dCurrentMissingSumZ = 0.0;

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

	if(fIsSplit) return;

	dWZ = dW*dZ;

	if(ISNA(dX))
	{
		dCurrentMissingSumZ += dWZ;
		dCurrentMissingTotalW += dW;
		cCurrentMissingN++;
		dCurrentRightSumZ -= dWZ;
		dCurrentRightTotalW -= dW;
		cCurrentRightN--;
	}
	else if(cCurrentVarClasses == 0)   // variable is continuous
	{
		if(dLastXValue > dX)
		{
			throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
		}

		// Evaluate the current split
		// the newest observation is still in the right child
		dCurrentSplitValue = 0.5*(dLastXValue + dX);
		if((dLastXValue != dX) &&
			(cCurrentLeftN >= cMinObsInNode) &&
			(cCurrentRightN >= cMinObsInNode) &&
			((lMonotone==0) ||
			(lMonotone*(dCurrentRightSumZ*dCurrentLeftTotalW -
						dCurrentLeftSumZ*dCurrentRightTotalW) > 0)))
		{
			dCurrentImprovement =
			  Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
									dCurrentMissingTotalW,
									dCurrentLeftSumZ,dCurrentRightSumZ,
									dCurrentMissingSumZ);
			if(dCurrentImprovement > dBestImprovement)
			{
				iBestSplitVar = iCurrentSplitVar;
				dBestSplitValue = dCurrentSplitValue;
				cBestVarClasses = 0;

				dBestLeftSumZ    = dCurrentLeftSumZ;
				dBestLeftTotalW  = dCurrentLeftTotalW;
				cBestLeftN       = cCurrentLeftN;
				dBestRightSumZ   = dCurrentRightSumZ;
				dBestRightTotalW = dCurrentRightTotalW;
				cBestRightN      = cCurrentRightN;
				dBestImprovement = dCurrentImprovement;
			}
		}

		// now move the new observation to the left
		// if another observation arrives we will evaluate this
		dCurrentLeftSumZ += dWZ;
		dCurrentLeftTotalW += dW;
		cCurrentLeftN++;
		dCurrentRightSumZ -= dWZ;
		dCurrentRightTotalW -= dW;
		cCurrentRightN--;

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

	  if(fIsSplit) return;

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
	      dCurrentSplitValue = (double)i;

	      dCurrentLeftSumZ    += adGroupSumZ[aiCurrentCategory[i]];
	      dCurrentLeftTotalW  += adGroupW[aiCurrentCategory[i]];
	      cCurrentLeftN       += acGroupN[aiCurrentCategory[i]];
	      dCurrentRightSumZ   -= adGroupSumZ[aiCurrentCategory[i]];
	      dCurrentRightTotalW -= adGroupW[aiCurrentCategory[i]];
	      cCurrentRightN      -= acGroupN[aiCurrentCategory[i]];

	      dCurrentImprovement =
		Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
				   dCurrentMissingTotalW,
				   dCurrentLeftSumZ,dCurrentRightSumZ,
				   dCurrentMissingSumZ);
	      if((cCurrentLeftN >= cMinObsInNode) &&
		 (cCurrentRightN >= cMinObsInNode) &&
		 (dCurrentImprovement > dBestImprovement))
	        {
		  dBestSplitValue = dCurrentSplitValue;
		  if(iBestSplitVar != iCurrentSplitVar)
	            {
		      iBestSplitVar = iCurrentSplitVar;
		      cBestVarClasses = cCurrentVarClasses;
		      std::copy(aiCurrentCategory.begin(),
				aiCurrentCategory.end(),
				aiBestCategory.begin());
	            }

		  dBestLeftSumZ      = dCurrentLeftSumZ;
		  dBestLeftTotalW    = dCurrentLeftTotalW;
		  cBestLeftN         = cCurrentLeftN;
		  dBestRightSumZ     = dCurrentRightSumZ;
		  dBestRightTotalW   = dCurrentRightTotalW;
		  cBestRightN        = cCurrentRightN;
		  dBestImprovement   = dCurrentImprovement;
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

	dBestLeftSumZ       = 0.0;
	dBestLeftTotalW     = 0.0;
	cBestLeftN          = 0;
	dCurrentLeftSumZ    = 0.0;
	dCurrentLeftTotalW  = 0.0;
	cCurrentLeftN       = 0;

	dBestRightSumZ      = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
	dBestRightTotalW    = nodeToSplit.dTrainW;
	cBestRightN         = nodeToSplit.cN;
	dCurrentRightSumZ   = 0.0;
	dCurrentRightTotalW = nodeToSplit.dTrainW;
	cCurrentRightN      = nodeToSplit.cN;

	dBestMissingSumZ      = 0.0;
	dBestMissingTotalW    = 0.0;
	cBestMissingN         = 0;
	dCurrentMissingSumZ   = 0.0;
	dCurrentMissingTotalW = 0.0;
	cCurrentMissingN      = 0;

	dBestImprovement    = 0.0;
	iBestSplitVar       = UINT_MAX;

	dCurrentImprovement = 0.0;
	iCurrentSplitVar    = UINT_MAX;
	dCurrentSplitValue  = -HUGE_VAL;

	fIsSplit = false;

	//this->pThisNode = pThisNode;
	//this->ppParentPointerToThisNode = ppParentPointerToThisNode;

	/*fIsSplit = false;
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
  if(fIsSplit) return;

  if (int(cCurrentVarClasses) > adGroupSumZ.size()) {
    throw GBM::failure("too many variable classes");
  }

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cCurrentVarClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + cCurrentVarClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + cCurrentVarClasses, 0);

  iCurrentSplitVar = iWhichVar;
  this->cCurrentVarClasses = cCurrentVarClasses;

  dCurrentLeftSumZ      = 0.0;
  dCurrentLeftTotalW    = 0.0;
  cCurrentLeftN         = 0;
  dCurrentRightSumZ     = dInitSumZ;
  dCurrentRightTotalW   = dInitTotalW;
  cCurrentRightN        = cInitN;
  dCurrentMissingSumZ   = 0.0;
  dCurrentMissingTotalW = 0.0;
  cCurrentMissingN      = 0;

  dCurrentImprovement = 0.0;

  dLastXValue = -HUGE_VAL;
}

void VarSplitter::WrapUpCurrentVariable()
{
  if(iCurrentSplitVar == iBestSplitVar)
    {
      if(cCurrentMissingN > 0)
        {
	  dBestMissingSumZ   = dCurrentMissingSumZ;
	  dBestMissingTotalW = dCurrentMissingTotalW;
	  cBestMissingN      = cCurrentMissingN;
        }
      else // DEBUG: consider a weighted average with parent node?
        {
	  dBestMissingSumZ   = dInitSumZ;
	  dBestMissingTotalW = dInitTotalW;
	  cBestMissingN      = 0;
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


