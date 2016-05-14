//------------------------------------------------------------------------------
//
//  File:       origCoxState.h
//
//  Description: original CoxPH methods
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __origCoxState_h__
#define __origCoxState_h__

//------------------------------
// Includes
//------------------------------
#include "coxph.h"
#include "genericCoxState.h"
#include <Rcpp.h>


//------------------------------
// Class Definition
//------------------------------
class OrigCoxState: public GenericCoxState
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	OrigCoxState(CCoxPH* coxPhPtr): coxPh(coxPhPtr){};

	//---------------------
	// Public destructor
	//---------------------
	~OrigCoxState(){coxPh = NULL;};

	//---------------------
	// Public Functions
	//---------------------
	void ComputeWorkingResponse
	(
		const CDataset& data,
	    const double *adF,
	    double *adZ
	)
	{

		double dF = 0.0;
		double dTot = 0.0;
		double dRiskTot = 0.0;

		vector<double> vecdRiskTot(data.get_trainSize(), 0.0);
		dRiskTot = 0.0;
		for(unsigned long i=0; i<data.get_trainSize(); i++)
		{
			if(data.GetBagElem(i))
			{
				dF = adF[i] +  data.offset_ptr()[i];
				dRiskTot += data.weight_ptr()[i]*std::exp(dF);
				vecdRiskTot[i] = dRiskTot;
			}
		}
		dTot = 0.0;
		for(long i= data.get_trainSize()-1; i != -1; i--)
		{
			if(data.GetBagElem(i))
			{
				if(coxPh->StatusVec()[i]==1.0)
				{
					dTot += data.weight_ptr()[i]/vecdRiskTot[i];
				}
				dF = adF[i] +  data.offset_ptr()[i];
				adZ[i] = coxPh->StatusVec()[i] - std::exp(dF)*dTot;
			}
		}


	}

	void FitBestConstant
	(
		const CDataset& data,
	    const double *adF,
	    unsigned long cTermNodes,
	    double* adZ,
	    CTreeComps& treeComps
	)
	{
	    double dF = 0.0;
	    double dRiskTot = 0.0;
	    unsigned long i = 0;
	    unsigned long k = 0;
	    unsigned long m = 0;

	    double dTemp = 0.0;
	    bool fTemp = false;
	    unsigned long K = 0;

	    vector<double> vecdP;
	    vector<double> vecdG;
	    vector<unsigned long> veciK2Node(cTermNodes, 0);
	    vector<unsigned long> veciNode2K(cTermNodes, 0);

	    matrix<double> matH;
	    matrix<double> matHinv;

	    for(i=0; i<cTermNodes; i++)
	    {
	        veciNode2K[i] = 0;
	        if(treeComps.GetTermNodes()[i]->cN >= treeComps.GetMinNodeObs())
	        {
	            veciK2Node[K] = i;
	            veciNode2K[i] = K;
	            K++;
	        }
	    }

	    vecdP.resize(K);

	    matH.setactualsize(K-1);
	    vecdG.resize(K-1);
	    vecdG.assign(K-1,0.0);

	    // zero the Hessian
	    for(k=0; k<K-1; k++)
	    {
	        for(m=0; m<K-1; m++)
	        {
	            matH.setvalue(k,m,0.0);
	        }
	    }

	    // get the gradient & Hessian, Ridgeway (1999) pp. 100-101
	    // correction from Ridgeway (1999): fix terminal node K-1 prediction to 0.0
	    //      for identifiability
	    dRiskTot = 0.0;
	    vecdP.assign(K,0.0);
	    for(i=0; i<data.get_trainSize(); i++)
	    {
	        if(data.GetBagElem(i) && (treeComps.GetTermNodes()[treeComps.GetNodeAssign()[i]]->cN >= treeComps.GetMinNodeObs()))
	        {
	            dF = adF[i] + data.offset_ptr()[i];
	            vecdP[veciNode2K[treeComps.GetNodeAssign()[i]]] += data.weight_ptr()[i]*std::exp(dF);
	            dRiskTot += data.weight_ptr()[i]*std::exp(dF);

	            if(coxPh->StatusVec()[i]==1.0)
	            {
	                // compute g and H
	                for(k=0; k<K-1; k++)
	                {
	                    vecdG[k] +=
	                        data.weight_ptr()[i]*((treeComps.GetNodeAssign()[i]==veciK2Node[k]) - vecdP[k]/dRiskTot);

	                    matH.getvalue(k,k,dTemp,fTemp);
	                    matH.setvalue(k,k,dTemp -
	                        data.weight_ptr()[i]*vecdP[k]/dRiskTot*(1-vecdP[k]/dRiskTot));
	                    for(m=0; m<k; m++)
	                    {
	                        matH.getvalue(k,m,dTemp,fTemp);
	                        dTemp += data.weight_ptr()[i]*vecdP[k]/dRiskTot*vecdP[m]/dRiskTot;
	                        matH.setvalue(k,m,dTemp);
	                        matH.setvalue(m,k,dTemp);
	                    }
	                }
	            }
	        }
	    }

	    /*
	    for(k=0; k<K-1; k++)
	    {
	        for(m=0; m<K-1; m++)
	        {
	            matH.getvalue(k,m,dTemp,fTemp);
	            Rprintf("%f ",dTemp);
	        }
	        Rprintf("\n");
	    }
	    */

	    // one step to get leaf predictions
	    matH.invert();

	    for(k=0; k<cTermNodes; k++)
	    {
	        treeComps.GetTermNodes()[k]->dPrediction = 0.0;
	    }
	    for(m=0; m<K-1; m++)
	    {
	        for(k=0; k<K-1; k++)
	        {
	            matH.getvalue(k,m,dTemp,fTemp);
	            if(!R_FINITE(dTemp)) // occurs if matH was not invertible
	            {
	                treeComps.GetTermNodes()[veciK2Node[k]]->dPrediction = 0.0;
	                break;
	            }
	            else
	            {
	                treeComps.GetTermNodes()[veciK2Node[k]]->dPrediction -= dTemp*vecdG[m];
	            }
	          }
	    }
	    // vecpTermNodes[veciK2Node[K-1]]->dPrediction = 0.0; // already set to 0.0
	}

	double Deviance
	(
		const long cLength,
		const CDataset& data,
	    const double *adF
	)
	{
	    unsigned long i=0;
	    double dL = 0.0;
	    double dF = 0.0;
	    double dW = 0.0;
	    double dTotalAtRisk = 0.0;



		// Original CoxPH implementation
	    dTotalAtRisk = 0.0;
	    for(i=0; i!=cLength; i++)
	    {
	        dF = adF[i] +  data.offset_ptr()[i];
	        dTotalAtRisk += data.weight_ptr()[i]*std::exp(dF);
	        if(coxPh->StatusVec()[i]==1.0)
	        {
	            dL += data.weight_ptr()[i]*(dF - std::log(dTotalAtRisk));
	            dW += data.weight_ptr()[i];
	        }
	    }

	    //TODO: Check if weights are all zero for validation set
	   if((dW == 0.0) && (dL == 0.0))
	   {
		   return nan("");
	   }
	   else if(dW == 0.0)
	   {
		   return copysign(HUGE_VAL, -dL);
	   }

	    return -2*dL/dW;
	}

	double BagImprovement
	(
		const CDataset& data,
		const double *adF,
		const bag& afInBag,
	  const double shrinkage,
	  const double* adFadj
	)
	{
	    double dReturnValue = 0.0;
	    double dNum = 0.0;
	    double dDen = 0.0;
	    double dF = 0.0;
	    double dW = 0.0;
	    unsigned long i = 0;

	    dNum = 0.0;
	    dDen = 0.0;
	    for(i=0; i< data.get_trainSize(); i++)
	    {
	        if(!data.GetBagElem(i))
	        {
	            dNum += data.weight_ptr()[i]*std::exp(dF + shrinkage*adFadj[i]);
	            dDen += data.weight_ptr()[i]*std::exp(dF);
	            if(coxPh->StatusVec()[i]==1.0)
	            {
	                dReturnValue +=
	                    data.weight_ptr()[i]*(shrinkage*adFadj[i] - std::log(dNum) + log(dDen));
	                dW += data.weight_ptr()[i];
	            }
	        }
	    }

	    return dReturnValue/dW;
	}

private:
	CCoxPH* coxPh;
};
#endif //__origCoxState_h__
