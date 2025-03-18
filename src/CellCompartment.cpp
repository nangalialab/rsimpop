/*
 * CellCompartment.cpp
 *
 *  Created on: 20 Mar 2020
 *      Author: nw14
 */

#include "CellCompartment.h"
#include "CellSimulation.h"
//#include <R.h>
#include <Rmath.h>
#include <string>
#include <stack>
#include "Event.h"
using namespace std;
const int MAX_NUM_ADDED_DRIVERS=100000;
CellCompartment::CellCompartment(int id,int targetPopSize,double divisionRate,double deathRate,std::vector<std::pair<double,int>> fitnessID): 
                                 id(id),mTargetPopSize(targetPopSize),mDivisionRate(divisionRate),mDeathRate(deathRate) {  
  nsub=fitnessID.size();
  for(int i=0;i<nsub;i++){
    idxByID[fitnessID[i].second]=i;
    mFitness.push_back(fitnessID[i].first);
    mID.push_back(fitnessID[i].second);
    vector<shared_ptr<PhyloNode>> tmp;
    subCompartments.push_back(tmp);
    bSubActive.push_back(false);
  }
  totalpop=0;
  rndGen=RandomNumberGenerator::getInstance();
  //Get random number generator
  prob=new double[nsub];
  printf("Cellcomparment : %d : divRate=%3.2f deathRate=%3.2f : targetpop=%d\n",id,mDivisionRate,mDeathRate,mTargetPopSize);
  if(mDivisionRate < -0.001){
    active=false;
  }else{
    active=true;
  }
}

CellCompartment::~CellCompartment() {

	//printf("destroying %d",id);

}

int CellCompartment::getSub(int ID){
	return idxByID[ID];
}


void CellCompartment::addNode(shared_ptr<PhyloNode> node,int subid){
	subCompartments[subid].push_back(node);
	totalpop++;
	bSubActive[subid]=true;//not ideal dirty way of setting bSubActive flag if we start with non-empty subcompartment
	if(totalpop>=mTargetPopSize){
		atEquilibrium=true;
	}
}

void CellCompartment::setNumNonEmptyIndices(){
	//printf("In setNumNonEmptyIndices\n");
	nonEmptyCompartmentIndices.clear();
	emptyCompartmentIndices.clear();
	for(int i=0;i<subCompartments.size();i++){
		if(subCompartments[i].size()>0){
			nonEmptyCompartmentIndices.push_back(i);
		}else{
			emptyCompartmentIndices.push_back(i);
		}
	}
	numNonEmptyCompartments=nonEmptyCompartmentIndices.size();
	numEmptyCompartments=emptyCompartmentIndices.size();
	delete [] prob;
	prob=new double[numNonEmptyCompartments];
}

void CellCompartment::checkPop(){
	int chk=0;
	int k=0;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
		// printf("%d=%d,",k++,subCompartment.size());
			 chk+=subCompartment.size();
	}
    //printf("::totals=%d %d\n",chk,totalpop);
	if(chk!=totalpop){
		printf("pop mismatch:%d != %d\n",chk,totalpop);
		throw "CellCompartment::checkPop: inconsistency in population tracking";
	}
}

void CellCompartment::printInfo(){
	printf("id=%d\nrate=%5.4f\ntargetPop=%d\nmTotalDeathRate=%7.6f\nmTotalDivRate=%7.6f\n%s\n",id,mDivisionRate,mTargetPopSize,mTotalDeathRate,mTotalDivRate,active ? "ACTIVE" : "INACTIVE");
	for(int k=0;k<nsub;k++){
		printf("fitness[%d]=%5.4f,count[%d]=%lu\n",k,mFitness[k],k,subCompartments[k].size());
	}
}


void CellCompartment::clear(){
	for(auto & sub:subCompartments){
		sub.clear();
	}
	totalpop=0;
}

vector<shared_ptr<PhyloNode>> CellCompartment::getNodes(){
	vector<shared_ptr<PhyloNode>> nodes;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
				 nodes.insert(nodes.end(),subCompartment.begin(),subCompartment.end());
	}
	return nodes;
}

vector<pair<bool,int>> CellCompartment::getSubCounts(){
	vector<pair<bool,int>> counts;
	int k=0;
	for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
				counts.push_back(pair<bool,int>(bSubActive[k],subCompartment.size()));
				k++;
	}
	return counts;
}

vector<tuple<int,int,double>> CellCompartment::getSubInfo(){
  vector<tuple<int,int,double>> infov;
  int k=0;
  double fitness=0;
  int popCount=0;
  int idDriver=0;
  for(const vector<shared_ptr<PhyloNode>> & subCompartment : subCompartments){
    popCount=subCompartment.size();
    if(popCount>0){
      //printf("getSubInfo k=%d  n=%d\n",k,popCount);
      fitness=mFitness[k];
      //printf("getSubInfo f=%3.2f\n",fitness);
      idDriver=mID[k];
      //printf("getSubInfo id=%d\n",idDriver);
      infov.push_back(tuple<int,int,double>(idDriver,popCount,fitness));
      //printf("getSubInfo done push id=%d\n",idDriver);
    }
    k++;
  }
  //printf("getSubInfo:done\n");
  return infov;
}


double CellCompartment::getTotalDivisionRate(){
	if(mDivisionRate<0){
		return 0.0;
	}
	double tot=0.0;
	int i;
	int ii;
	for(i=0;i<numNonEmptyCompartments;i++){
		ii=nonEmptyCompartmentIndices[i];
		tot+=mDivisionRate*(1+mFitness[ii])*subCompartments[ii].size();
	}
	return tot;
}

double CellCompartment::getTotalDeathRate(){
  if(mDivisionRate<0){
    return 0.0;
  }
  double tot=0.0;
  int i;
  int ii;
  for(i=0;i<numNonEmptyCompartments;i++){
    ii=nonEmptyCompartmentIndices[i];
    tot+=mDeathRate*subCompartments[ii].size();
  }
  return tot;
}

double CellCompartment::getTotalRate(){
	double tot=mTotalDivRate+mTotalDeathRate+mTotalAsymmetricDivRate+mTotalSymmetricDivRate;
	if(tot<0.0){
		printInfo();
		throw "getTotalRate: UnexpectedNegative";

	}
//	if(atEquilibrium){
		return tot;//mTotalDivRate+mTotalDeathRate;
	//}else{
	  //Following means we have pure exponential growth until equilibrium
		//return mTotalDivRate+mTotalDeathRate;
	//}
}

void CellCompartment::setRates(){
  double totalDeathRate=0.0;
  double tol=1e-6;
  int sz=0;
  int ii=0;
	//We want 20% (alpha) of the deviation to go per day - so rate_delta needs to be rate_delta=0.05*pop_delta
	//printf("In set rates\n");
	mTotalDivRate=0.0;
	totalDeathRate=0.0;
	mTotalAsymmetricDivRate=0;
	mTotalSymmetricDivRate=0;
	if(!active || mDivisionRate<0.0){
	
		// Do nothing
	}else{

		//mTotalDivRate=getTotalDivisionRate();
	  //totalDeathRate=getTotalDeathRate();
	  
	  //
	  mTotalAsymmetricDivRate=0;
	  for(int i=0;i<numNonEmptyCompartments;i++){
	    ii=nonEmptyCompartmentIndices[i];
	    sz=subCompartments[ii].size();
	    mTotalDivRate+=mDivisionRate*(1+mFitness[ii])*sz;
	    totalDeathRate=mDeathRate*sz;
	    sz=subCompartments[ii].size();
	    for(pair<shared_ptr<CellCompartment>,double> otherCompartment : asymmetricDifferentiation){
	      mTotalAsymmetricDivRate+=sz*otherCompartment.second;
	    }
	    for(pair<shared_ptr<CellCompartment>,double> otherCompartment : symmetricDifferentiation){
	      mTotalSymmetricDivRate+=sz*otherCompartment.second;
	    }
	  }
	  //TODO Review this. e.g. actual contraction rate is also governed by flows coming in (not tracked here..)
	  if(totalDeathRate>mTotalDivRate){
	    throw "deathRate: Should be <= divrate!";
	  }
	  
	  // if we are within 1% of the target population size or above the target then actively change death rate
	  // to target the population size
	  double targetTolerance=0.01;
	  if(totalDeathRate<tol && !atEquilibrium){
	    mTotalDeathRate=0.0;
	    return;
	  }
	  if(atEquilibrium && totalpop>(1-targetTolerance)*mTargetPopSize){
	    double incoming=0;
	    // Following only works if the feeding compartments only only feed this compartment..
	    for(pair<shared_ptr<CellCompartment>,double> otherCompartment : incomingAsymmetricDifferentiation){
	      incoming+=otherCompartment.first->mTotalAsymmetricDivRate*otherCompartment.second;
	    }
	    for(pair<shared_ptr<CellCompartment>,double> otherCompartment : incomingSymmetricDifferentiation){
	      incoming+=2*otherCompartment.first->mTotalSymmetricDivRate*otherCompartment.second;
	    }
	    
      // Actively target a population size
	    double tmp=(mTotalDivRate+incoming-mTotalSymmetricDivRate-alpha*(mTargetPopSize-totalpop));
	    mTotalDeathRate=tmp>0?tmp:0.0;
	    //printf("id=%d incoming=%7.6f divrate=%7.6f sdivrate=%7.6f mTargetPopSize=%d totalpop=%d  deathrate=%7.6f\n",id,incoming,mTotalDivRate,mTotalSymmetricDivRate,mTargetPopSize,totalpop,mTotalDeathRate);
	    return;
	  }
	  mTotalDeathRate=totalDeathRate;
	}
}


/* Moved to CellSimulation.cpp
void CellCompartment::doEvent(CellSimulation & sim);
*/

/**
 * This method randomly selects a cell in this compartment (equal probability) and reassigns it to a new sub-compartment with
 * a new fitness according to an additive model (new fitness=old fitness + fitness).
 */
double CellCompartment::addDriver(CellSimulation & sim,double ts){   //,double fitness
		
		pair<int,double> driver=sim.getNextDriver();
		int driverid=driver.first;
		double fitnessValue=driver.second;
		setNumNonEmptyIndices();
		//double prob[nsub]; //what does nsub represent? what is the difference between compartment and sub-compartment?
 		for(int i=0;i<numNonEmptyCompartments;i++){
			prob[i]=subCompartments[nonEmptyCompartmentIndices[i]].size();
		}
		int ii=rndGen->sample(numNonEmptyCompartments,prob,false); //random selection of non-empty subcompartment
		int i=nonEmptyCompartmentIndices[ii];
		int sz=subCompartments[i].size();
		int k=rndGen->sample(sz);
		shared_ptr<PhyloNode> node=subCompartments[i][k];  //random selection of a cell from the chosen subcompartment
		//Node ID of event only matters when provided externally.
		Event event(-1,ts,id,driverid,driverid); 
		shared_ptr<Event> thisEvent=std::make_shared<Event>(event);
		node->addEvent(thisEvent);
		if(fitnessValue>0){
			double fitness = fitnessValue+mFitness[i];
			mFitnessOut.push_back(fitness);
			if(numEmptyCompartments==0){
				subCompartments.push_back(vector<shared_ptr<PhyloNode>>());
				subCompartments[nsub].push_back(node);
				mFitness.push_back(fitness);
				mID.push_back(driverid);
				idxByID[driverid]=nsub;    //sz;  //
				///subcompIndexInAddDriver.push_back(nsub);
				numOfAddedDrivers++;
				bSubActive.push_back(true);
				nsub++;
				//printf("adding new compartment: %d fit=%3.2f\n",nsub,fitness);
			}else{
				int l=emptyCompartmentIndices[0];
				subCompartments[l].push_back(node);
				mFitness[l]=fitness; 
				mID[l]=driverid;
				bSubActive[l]=true;
				idxByID[driverid]=l;
				//subcompIndexInAddDriver.push_back(l);
				numOfAddedDrivers++;
				//printf("recycling compartment: %d fit=%3.2f\n",nsub,fitness);
			}
    		//Remove from vector in old compartment
			if(k<sz-1){
				subCompartments[i][k]=subCompartments[i][sz-1];
			}
			subCompartments[i].pop_back();
			//checkPop();
			//Finally reset nonEmptyCompartments etc
			setNumNonEmptyIndices();
			return fitness; //fitness;
		}
		else
		{
		  throw "addDriver Unexpected 0 or negative driver!";
			return fitnessValue;
		  
		}
	
	//Create a new sub-compartment for this node and remove from old compartment..

	
	
	
}

void CellCompartment::addSymmetricDifferentiationRate(shared_ptr<CellCompartment> otherCompartment,double rate){
  printf("Adding symmetricDifferentiation %d -> %d : rate = %5.4f\n",this->id,otherCompartment->id,rate);
  symmetricDifferentiation.push_back(pair<shared_ptr<CellCompartment>,double>(otherCompartment,rate));
  
}
void CellCompartment::addAsymmetricDifferentiationRate(shared_ptr<CellCompartment> otherCompartment,double rate){
  printf("Adding aSymmetricDifferentiation %d -> %d : rate = %5.4f\n",this->id,otherCompartment->id,rate);
  asymmetricDifferentiation.push_back(pair<shared_ptr<CellCompartment>,double>(otherCompartment,rate));
}

void CellCompartment::addIncomingSymmetricDifferentiationRate(shared_ptr<CellCompartment> otherCompartment,double rate){
  printf("Adding incoming symmetricDifferentiation %d -> %d : probability this compartment = %3.2f\n",otherCompartment->id,this->id,rate);
  incomingSymmetricDifferentiation.push_back(pair<shared_ptr<CellCompartment>,double>(otherCompartment,rate));
  
}
void CellCompartment::addIncomingAsymmetricDifferentiationRate(shared_ptr<CellCompartment> otherCompartment,double rate){
  printf("Adding incoming aSymmetricDifferentiation %d -> %d : probability this compartment = %3.2f\n",otherCompartment->id,this->id,rate);
  incomingAsymmetricDifferentiation.push_back(pair<shared_ptr<CellCompartment>,double>(otherCompartment,rate));
}
