
// $Id: computePijComponent.cpp 9253 2011-01-31 01:37:21Z rubi $

#include "definitions.h"
#include "treeIt.h"
#include "computePijComponent.h"
#include "logFile.h"

void computePijHomSpec::fillPij(
	const MDOUBLE dis,
	const stochasticProcess& sp,
	int derivationOrder,
	bool isReversible)
{
	if (!(isReversible && sp.isReversible())) {// if one is false 
		isReversible = false;
	}
    resize(sp.alphabetSize());
    size_t i;
    for (i=0; i<sp.alphabetSize(); i++) {
	switch (derivationOrder) {
	case 0:
	  _V[i][i] = sp.Pij_t(i,i,dis);
	 // if (_V[i][i] > 1 || _V[i][i] < 0) {
	//	  errorMsg::reportError("error in function computePijHomSpec::fillPij", 1);
	 // }
	  break;
	case 1:
	  _V[i][i] = sp.dPij_dt(i,i,dis);
	  break;
	case 2:
	  _V[i][i] = sp.d2Pij_dt2(i,i,dis);
	  break;
	default:
	  errorMsg::reportError("error in function fillPij - derivationOrder must be 0, 1 or 2");		
	}

	size_t j;
	  
	for (j=i+1; j<sp.alphabetSize(); j++) {
	  switch (derivationOrder) {
	  case 0:
	      _V[i][j] = sp.Pij_t(i,j,dis);
	      if ((_V[i][j] == 0 )&& (dis !=0)){
		  
		  _V[i][j] = EPSILON;
	      }
		 // if (_V[i][j] > 1 || _V[i][j] < 0) {
		//	  errorMsg::reportError("error in function computePijHomSpec::fillPij", 1);
		 // }
	      break;
	  case 1:
	      _V[i][j] = sp.dPij_dt(i,j,dis);
	      break;
	  case 2:
	      _V[i][j] = sp.d2Pij_dt2(i,j,dis);
	      break;
	  default:
	      errorMsg::reportError("error in function fillPij - derivationOrder must be 0, 1 or 2");		
	  }
	  if (sp.freq(j) == 0.0) {
	      if (isReversible) {
		  errorMsg::reportError("error in function fillPij");
	      }
	      
	  }
//	  else {
	      if (isReversible){
			  _V[j][i] = _V[i][j]* sp.freq(i)/sp.freq(j);
			  //if (_V[j][i] > 1 || _V[j][i] < 0) {
				//  cout << "i= " << i << " freq(i)= " << sp.freq(i) << " j= " << j << " freq(j)= " << sp.freq(j) << endl;
				//  errorMsg::reportError("error in function computePijHomSpec::fillPij", 1);
			  //}
	      }
	      else {
		  switch (derivationOrder) {
		  case 0:
		      _V[j][i] = sp.Pij_t(j,i,dis);
		      if ((_V[j][i] == 0 )&& (dis !=0))
			  _V[j][i] = EPSILON;
		      break;
		  case 1:
		      _V[j][i] = sp.dPij_dt(j,i,dis);
		      break;
		  case 2:
		      _V[j][i] = sp.d2Pij_dt2(j,i,dis);
		      break;
		  default:
		      errorMsg::reportError("error in function fillPij - derivationOrder must be 0, 1 or 2");		
		  }
	      }
//	  }
	}
    }
}


void computePijHom::fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder, bool isReversible) {
	_V.resize(et.getNodesNum());
	treeIterTopDownConst tIt(et);
	tree::nodeP myNode = tIt.first();
	{// skipping the root, but allocating place for the root pij even if they are not use
	 // to maintain that all arrays have the same size.
		_V[myNode->id()].resize(sp.alphabetSize());
	}
	LOGDO(50,et.output(myLog::LogFile(),tree::ANCESTOR));
	LOGDO(50,et.output(myLog::LogFile(),tree::PHYLIP));
	for (; myNode != tIt.end(); myNode = tIt.next()) {
	  if (!(myNode->isRoot()))
		  _V[myNode->id()].fillPij(myNode->dis2father()*sp.getGlobalRate(),sp,derivationOrder,isReversible);
//	  else
//	    myLog::LogFile()<<"ROOT IS "<<myNode->name()<<endl;
	}
}


void computePijGam::fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder, bool isReversible)
{
	_V.resize(sp.categories());
	for (size_t i=0; i < _V.size(); ++i) {
		tree cp = et;
		cp.multipleAllBranchesByFactor(sp.rates(i)/sp.getGlobalRate());// the global rate is taken care of in the hom pij.
		_V[i].fillPij(cp,sp,derivationOrder,isReversible);
	}
}

void PijGamMultipleStochasticProcesses::fillPij(const tree& et, const vector<stochasticProcess>& V_sp)
{
	_V.resize(V_sp.size());
	for (size_t i = 0; i < V_sp.size(); i++)
	{
		_V[i].fillPij(et, V_sp[i]);
	}
}

