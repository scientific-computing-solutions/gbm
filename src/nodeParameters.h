//------------------------------------------------------------------------------
//
//  File:       nodeParameters.h
//
//  Description:  header  contains the parameters used to split a node
//
//------------------------------------------------------------------------------

#ifndef NODEPARAMETERS_H
#define NODEPARAMETERS_H

//------------------------------
// Includes
//------------------------------
#include <Rcpp.h>

//------------------------------
// Struct Definition
//------------------------------
struct NodeDef
{
  //----------------------
  // Public Constructors
  //----------------------
  NodeDef() : numObs(0), weightResid(0), totalWeight(0) {};
  
  NodeDef(double weightResid, double totalWeight, long numObs) :
  numObs(numObs), weightResid(weightResid), totalWeight(totalWeight) {};
  
  //---------------------
  // Public Functions
  //---------------------
  void clear()
  {
    numObs = 0;
    weightResid = totalWeight = 0;
  };
  
  void increment(const double pred, const double trainWeight, long num)
  {
    weightResid += pred;
    totalWeight += trainWeight;
    numObs += num;
  };
  
  double prediction() const
  {
    return weightResid / totalWeight;
  };
  
  double unweightedGradient(const NodeDef& other) const
  {

    return weightResid * other.totalWeight - other.weightResid * totalWeight;
  };
  
  bool hasMinObs(long minObsInNode) const
  {
    return (numObs >= minObsInNode);
  }
  
  bool hasObs() const {
    return numObs;
  }

  double weightSum(const NodeDef& a) const {
    return totalWeight + a.totalWeight;
  }
  
  double weightSum(const NodeDef& a, const NodeDef& b) const {
    return totalWeight + a.weightSum(b);
  }

  double getTotalWeight() const {
    return totalWeight;
  }

  long getNumObs() const {
    return numObs;
  }
  
private:
  
  long numObs;
  double weightResid;
  double totalWeight;

};

//------------------------------
// Class Definition
//------------------------------
class NodeParams
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	NodeParams() : aiBestCategory(1024) {};

	//---------------------
	// Public destructor
	//---------------------
        ~NodeParams();

	//---------------------
	// Public Functions
	//---------------------
	void ResetSplitProperties(double weightedResiduals, double trainingWeight, unsigned long numObs,
				  double splitValue = -HUGE_VAL, unsigned long variableClasses=1, unsigned long splitVar = UINT_MAX);
	void UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement = 1)
	{
	  // Move data point from right node to missing
	  missing.increment(predIncrement, trainWIncrement, numIncrement);
	  right.increment(-predIncrement, -trainWIncrement, -numIncrement);
	}
	void UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement = 1)
	{
		// Move data point from right node to left node
		left.increment(predIncrement, trainWIncrement, numIncrement);
		right.increment(-predIncrement, -trainWIncrement, -numIncrement);
	}
	inline double GetImprovement() { return ImprovedResiduals;};
	bool SplitIsCorrMonotonic(long specifyMonotone)
	{
		return (
			(specifyMonotone == 0) ||
			((specifyMonotone * right.unweightedGradient(left)) > 0)
			);
	}
	void NodeGradResiduals()
	{
	  // Only need to look at left and right
	  if(!missing.hasObs())
	    {
	      ImprovedResiduals = left.unweightedGradient(right) /
		left.weightSum(right);
	    }
	  else
	    {
	      // Grad - left/right
	      ImprovedResiduals =
		(left.unweightedGradient(right) +
		 left.unweightedGradient(missing) +
		 right.unweightedGradient(missing)) /
		left.weightSum(right, missing);
	    }
	};
	
	bool HasMinNumOfObs(long minObsInNode)
	{
		return (left.hasMinObs(minObsInNode) &&
			  right.hasMinObs(minObsInNode));
	}
	inline void setBestCategory(std::vector<std::pair<double, int> >& groupMeanAndCat)
	{

	  int count = 0;
	  aiBestCategory.resize(groupMeanAndCat.size());
	  for(std::vector<std::pair<double, int> >::const_iterator it = groupMeanAndCat.begin();
	      it != groupMeanAndCat.end();
	      ++it)
	    {
	      aiBestCategory[count] = it->second;
	      count++;
	    }
	};

	NodeParams& operator=(const NodeParams& rhs)
	{
	  right = rhs.right;
	  left = rhs.left;
	  missing = rhs.missing;

	  SplitValue = rhs.SplitValue;
	  SplitVar = rhs.SplitVar;
	  SplitClass = rhs.SplitClass;
	  ImprovedResiduals = rhs.ImprovedResiduals;
	  
	  // Copy best category
	  aiBestCategory = rhs.aiBestCategory;
	  return *this;
	}
	bool hasMissing() const
	{
	  return missing.hasObs();
	};
	//---------------------
	// Public Variables
	//---------------------
	// Left Node Definition
	NodeDef left, right, missing;

	// Splitting values
	double SplitValue; // Continuous Split Value
	unsigned long SplitVar; // Which feature to split on
	unsigned long SplitClass; // Categorical Split Value
	std::vector<int> aiBestCategory; // Vector of levels ordering
	double ImprovedResiduals;
};

#endif // NODEPARAMETERS_H
