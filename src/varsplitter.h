//------------------------------------------------------------------------------
//
//  File:       varsplitter.h
//
//  Description: header for class that splits a node on a particular variable.
//
//------------------------------------------------------------------------------

#ifndef __varsplitter_h__
#define __varsplitter_h__

//------------------------------
// Includes
//------------------------------
#include "node.h"
#include "nodeParameters.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class VarSplitter
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	VarSplitter(unsigned long minNumObs);

	//---------------------
	// Public destructor
	//---------------------
	~VarSplitter();

	//---------------------
	// Public Functions
	//---------------------
	/*void SetForNode(CNode& nodeToSet);
	void SetForVariable(unsigned long iWhichVar, long cVarClasses);

	inline double GetBestImprovement() { return bestSplit.GetImprovement(); };
	void IncorporateObs(double dX,
			double dZ,
			double dW,
			long lMonotone);
	void EvaluateCategoricalSplit();
	void SetToSplit()
	{
		fIsSplit = true;
	};
	inline NodeParams GetBestSplit() { return bestSplit;};
	void Reset();*/
	 void IncorporateObs(double dX,
				double dZ,
				double dW,
				long lMonotone);

	void Set(CNode& nodeToSplit);
	void ResetForNewVar(unsigned long iWhichVar,
			long cVarClasses);

	static double Improvement
		(
			double dLeftW,
			double dRightW,
			double dMissingW,
			double dLeftSum,
			double dRightSum,
			double dMissingSum
		)
		{
			double dTemp = 0.0;
			double dResult = 0.0;

			if(dMissingW == 0.0)
			{
				dTemp = dLeftSum/dLeftW - dRightSum/dRightW;
				dResult = dLeftW*dRightW*dTemp*dTemp/(dLeftW+dRightW);
			}
			else
			{
				dTemp = dLeftSum/dLeftW - dRightSum/dRightW;
				dResult += dLeftW*dRightW*dTemp*dTemp;
				dTemp = dLeftSum/dLeftW - dMissingSum/dMissingW;
				dResult += dLeftW*dMissingW*dTemp*dTemp;
				dTemp = dRightSum/dRightW - dMissingSum/dMissingW;
				dResult += dRightW*dMissingW*dTemp*dTemp;
				dResult /= (dLeftW + dRightW + dMissingW);
			}

			return dResult;
		}

	double BestImprovement() { return dBestImprovement; }
	void SetupNewNodes(CNode& nodeToSplit)
	{
		nodeToSplit.SplitNode(iBestSplitVar,
		cBestVarClasses,
		 dBestSplitValue,
		dBestLeftSumZ,
		 dBestLeftTotalW,
		 cBestLeftN,
		 dBestRightSumZ,
		 dBestRightTotalW,
		cBestRightN,
		dBestMissingSumZ,
		dBestMissingTotalW,
		cBestMissingN,
		dBestImprovement,
		aiBestCategory);
	}

	void EvaluateCategoricalSplit();
	void WrapUpCurrentVariable();

	unsigned long iBestSplitVar;
	double dBestSplitValue;

	double dBestLeftSumZ;
	double dBestLeftTotalW;
	unsigned long cBestLeftN;

	double dBestRightSumZ;
	double dBestRightTotalW;
	unsigned long cBestRightN;

	double dBestMissingSumZ;
	double dBestMissingTotalW;
	unsigned long cBestMissingN;

	double dCurrentMissingSumZ;
	double dCurrentMissingTotalW;
	unsigned long cCurrentMissingN;

	long cCurrentVarClasses;
	long cBestVarClasses;
	double dInitTotalW;
	double dInitSumZ;
	unsigned long cInitN;
	double dBestImprovement;

private:

	unsigned long cMinObsInNode;



	double dCurrentLeftSumZ;
	double dCurrentLeftTotalW;
	unsigned long cCurrentLeftN;
	double dCurrentRightSumZ;
	double dCurrentRightTotalW;
	unsigned long cCurrentRightN;
	double dCurrentImprovement;
	unsigned long iCurrentSplitVar;
	double dCurrentSplitValue;

	double dLastXValue;

	std::vector<double> adGroupSumZ;
	std::vector<double> adGroupW;
	std::vector<unsigned long> acGroupN;
	std::vector<double> adGroupMean;
	// this is an int to fit in with R API
	// it's probably best not to ask.
	std::vector<int> aiCurrentCategory;
	std::vector<unsigned long> aiBestCategory;

	NodeParams bestSplit, proposedSplit;
	/*//---------------------
	// Private Functions
	//---------------------
	void WrapUpSplit();
	
	//---------------------
	// Private Variables
	//---------------------
	double InitTotalWeight, InitWeightResiduals, dLastXValue;
	unsigned long InitNumObs;
	unsigned long minObsInNode;

	bool fIsSplit;*/



};
#endif // __varplitter_h__
