/*
 * CellSimulation.cpp
 *
 *  Created on: 24 Mar 2020
 *      Author: nw14
 */

#include "CellSimulation.h"
#include "PhyloNode.h"
#include "Event.h"
#include "PopulationTimePoint.h"
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <algorithm>
#include <numeric>

using namespace std;


CellSimulation::CellSimulation(int * edges,
		int * ndivs,
		double * tBirth,
		int nedge,
		int ntips,
		vector<Event> events,
		vector<shared_ptr<CellCompartment>> cellCompartments,
		double startTime,
		double driverAcquisitionRate,
		double * trajectoryTs,
		int * trajectoryPop,
		double * trajectoryDivRate,
		double * trajectoryDeathRate, // Todo replace with trajectory class
		int * trajectoryCompartment,
		int trajectorySize,
		int maxSize,
		double * fitnessDistribution, 
		int maxPreviousDriverID,
		int driverFitnessSize,
		int bVerbose){ 
	
	//cout << "nrow trajectory: " << trajectorySize << endl;
	if(trajectorySize>=1 && trajectoryTs[0]!=0){
		for(int i=0;i<trajectorySize;i++){
			trajectoryStack.push(tuple<double,int,double,int,double>(trajectoryTs[trajectorySize-1-i],trajectoryPop[trajectorySize-1-i],trajectoryDivRate[trajectorySize-1-i],
			trajectoryCompartment[trajectorySize-1-i],trajectoryDeathRate[trajectorySize-1-i]));
  		} 
	}
	//cout << "trajectory size in C++: " << trajectoryStack.size() << endl;
	

	rndGen=RandomNumberGenerator::getInstance();
	int i,inode,iparent;
	int popsize=1;//Add a ROOT compartment

	for(const shared_ptr<CellCompartment> & compartment : cellCompartments){
		popsize+=compartment->mTargetPopSize;
		compartments.push_back(compartment);
		//compartment->printInfo();
	}

	int MAX_SIZE=maxSize;//4*max(popsize,nedge);
	for(int i=0;i<MAX_SIZE;i++){
		shared_ptr<PhyloNode> p=std::make_shared<PhyloNode>();
		nodeStack.push(p);
	}
	vector<shared_ptr<PhyloNode>> nodeArray;
	nodeArray.reserve(nedge+1);
	for(i=0;i<nedge+1;i++){
		nodeArray.push_back(getNewPhyloNode());
	}
	shared_ptr<PhyloNode> cnode;
	root=nodeArray[ntips];
	root->id=ntips+1;
	this->ntips=ntips;
	//printf("full pop size=%d  and max_size=%d  nedge=%d",popsize,MAX_SIZE,nedge);
	for(i=0;i<nedge;i++){
		inode=edges[i+nedge]-1;

		iparent=edges[i]-1;
		if(inode>nedge+1 || iparent> nedge+1 || inode<0 || iparent<0){
			throw "CellSimulation:Edge Matrix Values out of Range";
		}
		cnode=nodeArray[inode];
		cnode->ndiv=ndivs[i];
		cnode->tbirth=tBirth[i];
		cnode->id=inode+1;//ID same as in edge space...
		cnode->parent=nodeArray[iparent];
		auto pp=cnode->parent.lock();
		pp->addChild(cnode);
	}
	//printf("done allocating now doing cell compartments...\n");
	//cout << "root: born " << root->tbirth << endl;
	//now implement the cell compartments...
	for(const Event & event : events){
		shared_ptr<Event> thisEvent=std::make_shared<Event>(event);
		nodeArray[event.node-1]->addEvent(thisEvent);
		//cout << "event:" << thisEvent->node << ":" << thisEvent->value << endl;
	}
	currentTime=startTime;
	driverRate=driverAcquisitionRate;
	this->bVerbose=bVerbose;
	for(int i=0;i<driverFitnessSize;i++){
	  //Add backwars because stack is lifo
	  driverFitnessStack.push(fitnessDistribution[driverFitnessSize-1-i]); 
	}
	lastDriverID=maxPreviousDriverID;
	
	
	//printf("bVerbose=%d\n",bVerbose);
	//The following actually populates the CellCompartments by copying pointers to PhyloNodes that 
	//correspond to extant cells into their compartments and sub-compartments.
	setCompartmentInfo();
	
	resetNtips(false);
}

CellSimulation::~CellSimulation() {
	//printf("Destroying cellSimulation\n");
}

shared_ptr<PhyloNode>  CellSimulation::getNewPhyloNode(){
	if(nodeStack.empty()){
		throw "nodeStack is unexpectedly empty!";
	}
	shared_ptr<PhyloNode> out=nodeStack.top();
	nodeStack.pop();
	out->clear();
	return out;
}
void CellSimulation::recycle(shared_ptr<PhyloNode> node){
	node->clear();
	nodeStack.push(node);

}
//void evolveBirthDeath()

int CellSimulation::getNumberOfAddedDrivers(){  
	return nAddedDrivers;
}


int CellSimulation::run(double stopTime,bool bStopAtEquilibrium,bool bStopIfEmpty, int maxDriverCount){  //,  , double * nCompFitnessOut
	int i;
	double totrate;
	int ncomp=compartments.size();
	double rates[ncomp];
	double popbycomp[ncomp];
	///collect daily timestamps
	double lastSnap=currentTime;
	double lastYear=floor(currentTime/365.0);
	//cout << "current time is " << currentTime << endl;
	//cout << "last year is " << lastYear << endl;
	double alpha;
	bool allAtEquilibrium=true;
	bool bStop=false;
	int kk=0;
	int k;
	//printf("stop if @ equi=%d\n" ,bStopAtEquilibrium );
	int status=0;
	//printf("init params:t=%3.2f:till=%3.2f:rate=%3.2f:poptarget=%d\n",currentTime,stopTime,compartments[1]->mDivisionRate,compartments[1]->mTargetPopSize);
	while(currentTime<stopTime){  //stopTime is the number of simulation days
		//kk++;
		totrate=0.0;
		i=0;
		allAtEquilibrium=true;
		int pop=0;
		for(const shared_ptr<CellCompartment> & compartment : compartments){
			if(compartment->active && compartment->totalpop>0){
				if(bStopAtEquilibrium){
					allAtEquilibrium=allAtEquilibrium && compartment->atEquilibrium;
				}
				if(bStopIfEmpty){
					auto counts=compartment->getSubCounts();
					k=0;
					for(pair<bool,int> & cnt: counts){
						if(cnt.first && cnt.second==0 && k>0){
							printf("Compartment %d: One of the active driver sub-compartments is empty.. Stopping sim\n",compartment->id);
														bStop=true;
														break;
						}
						k++;
					}
				}
			}
			if(compartment->active){
			  popbycomp[i]=compartment->totalpop;
			}else{
			  popbycomp[i]=0.0;
			}
			rates[i]=compartment->getTotalRate();
			//printf("total rates %d = %7.6f\n",i,rates[i]);
			totrate+=rates[i++];
			pop+=compartment->totalpop;
			
		}
		if(bStop){
					break;
		}
		if(bStopAtEquilibrium && allAtEquilibrium){
			printf("All active compartments are either empty or at target popsize.. Stopping sim @ %7.6f\n",currentTime);
			status=1;
			break;
			
		}
		
		currentTime+=rndGen->getExponential(totrate+pop*driverRate);
		//printf("\n");
		for(i=0;i<ncomp;i++){
			rates[i]/=totrate;
		  //popbycomp[i]/=pop;
		  //printf("p[%d]=%3.2f",i,popbycomp[i]);
		}
		//printf("\n");
		//if(compartments.size()==4){
		//  printf("compartment 3 is active? %d:  total rate =%3.2g  %3.2g\n",compartments[3]->active,compartments[3]->getTotalDivisionRate(),compartments[3]->getTotalRate());
		//}
		
		if(driverRate>1e-10){
		  if(rndGen->getUniform()<pop*driverRate/(totrate+pop*driverRate)){
		    status=2;
		    i=rndGen->sample(ncomp,popbycomp,false);
		    //i=1;//Force all drivers into compartment 1
		    double dv=compartments[i]->addDriver(*this,currentTime);
		    if( dv > 0){
		      nAddedDrivers++; 
		      mfitnessValues.push_back(dv); 
		      continue;
		    }
		  }
		}
		i=rndGen->sample(ncomp,rates,true);
		compartments[i]->doEvent(*this);
		if(currentTime-lastSnap>=1.0){
			snap();
			lastSnap=currentTime;
			if(maxDriverCount>0){
			  ndrivers=populationTrace.back().getTotalDriverCount();// std::get<2>(populationTrace.back());
			  if(ndrivers>=maxDriverCount){
			    return status;
			  }
			}
		}
		if(currentTime/365.0 -lastYear > 1 ){
			//printf("PROGRESS: T=%3.2f year.  population=%d",currentTime/365.0,ntips);
			if(bVerbose){
			printf("PROGRESS: T=%d year.  population=%d",int(currentTime/365.0),ntips);
			if (!populationTrace.empty()){
			//	printf(":drivers=%d\n",std::get<2>(populationTrace.back()));
			  printf(":drivers=%d\n",populationTrace.back().getTotalDriverCount());// report the main compartment fornow
			}else{
				printf("\n");
			}
			}
			lastYear=currentTime/365.0;
		}
		
		while(currentTime>stopTime && !trajectoryStack.empty()){
				std::tuple<double,int,double,int,double> newParam = trajectoryStack.top();
				int compartmentID = std::get<3>(newParam);
				stopTime=std::get<0>(newParam);
				//Logic is currently broken for multiple compartments - need a per compartment trajectory stack
				if(currentTime<=stopTime){
				  compartments[compartmentID]->mDivisionRate = std::get<2>(newParam);
				  compartments[compartmentID]->mTargetPopSize = std::get<1>(newParam);
				  compartments[compartmentID]->mDeathRate = std::get<4>(newParam);
				  //stopTime = std::get<0>(newParam);
				  if(bVerbose>1){
				    printf("Updated params:t=%3.2f:till=%3.2f:comp=%d;rate=%5.6f:dr==%5.6f;poptarget=%d\n",currentTime,stopTime,compartmentID,compartments[compartmentID]->mDivisionRate,compartments[compartmentID]->mDeathRate,compartments[compartmentID]->mTargetPopSize);
				  }
				  trajectoryStack.pop();
				}else{
				  //setting is in past so just pop..
				  trajectoryStack.pop();
				}
		}

	}
	printf("done!");
	return status;

}



/**
 * Only valid immediately after construction
 */
void CellSimulation::deleteTips(std::set<int> tipsToDelete){

	//printf("deleting %d tips in tree with %d tips\n",tipsToDelete.size(),ntips);
	//now loop through all nodes and delete as required...

	for(auto & compartment: compartments){
		auto nodes=compartment->getNodes();
		for( shared_ptr<PhyloNode> & node: nodes){
			auto it=tipsToDelete.find(node->id);
			if(it != tipsToDelete.end()){
				die(node);
			}else{
				//printf("%d,",node->id);
			}
		}
	}
	//Can now reconstruct compartment info..
	setCompartmentInfo();
	resetNtips();
}

void CellSimulation::snap(){
  double tol=1e-8;
	resetNtips();
  if(!populationTrace.empty()){
    if(currentTime>(populationTrace.back().timeStamp+tol)){
	  //populationTrace.push_back(tuple<double,int,int>(currentTime,ntips,ndrivers));
	    populationTrace.push_back(PopulationTimePoint( currentTime,compartmentIDs,populationCounts,driverCounts));
    }
  }else{
    populationTrace.push_back(PopulationTimePoint( currentTime,compartmentIDs,populationCounts,driverCounts));
  }
}


vector<PopulationTimePoint> CellSimulation::getPopulationTrace(){
	return populationTrace;
}
//void CellSimulation::assignCellCompartments();

void CellSimulation::resetNtips(bool bCheckMatch){
	int totalpop=0;
	int ndr=0;
	driverCounts.resize(compartments.size(),0);
	compartmentIDs.resize(compartments.size(),0);
	populationCounts.resize(compartments.size(),0);
	int j=0;
	int cndr=0;
	//printf("in resettips\n");
	for(const shared_ptr<CellCompartment> & compartment : compartments){
				vector<pair<bool,int>> subc=compartment->getSubCounts();
				int n=subc.size();
				cndr=0;
				for(int i=1;i<n;i++){
					if(subc[i].first){
					  cndr+=subc[i].second;
					}
				}
				ndr+=cndr;
				compartmentIDs[j]=compartment->id;
				populationCounts[j]=compartment->totalpop;
				driverCounts[j]=cndr;
				totalpop+=compartment->totalpop;
				compartment->checkPop();
				j++;
	}
	//printf("done resettips\n");
	if(bCheckMatch && totalpop!=ntips){
		printf("Mismatch in tip count %d != %d!\n",totalpop,ntips);
		throw "CellSimulation: Inconsistency - mismatch between sum of compartment pops and CellSimulation pop";
	}
	ntips=totalpop;
	ndrivers=ndr;

}

void CellSimulation::resetNtips(){
	resetNtips(true);
}

/**
 * Populates preallocated arrays with tree info.
 */
void CellSimulation::populate_edge_info(int * edgesOut,int * ndivsOut,int * stateOut,int * driverIDOut,double * tBirthOut,int max_nedge,int * nedgeOut,int *nInternalNodes,int * numTip,vector<Event> & eventsOut){
	resetNtips();
	//int ntip=ntips;
	int root_id=ntips+1;
	int tip_counter=1;
	int row_counter=0;
	int internal_counter=root_id+1;
	update_edge_matrix(edgesOut,ndivsOut,stateOut,driverIDOut,tBirthOut,root_id,&internal_counter,&tip_counter,&row_counter,root,max_nedge,eventsOut,0,0);
	if(tip_counter-1 != ntips){
		printf("Mismatch in tip count! %d %d\n",tip_counter,ntips);
		throw "CellSimulation::populate_edge_info: tip count mismatch";
	}
	*numTip=ntips;
	*nedgeOut=row_counter;
	*nInternalNodes=internal_counter-root_id;
}

//TODO : Fix export of events!
void CellSimulation::update_edge_matrix(int * edges,int * ndivsOut,int * stateOut,int * driverIDOut,
		double *tBirthOut,int parent_counter,int * internal_counter,
		int * tip_counter,int * row_counter,shared_ptr<PhyloNode> thisNode,int max_size,vector<Event> & eventsOut,int parentState,int parentDriverID){
	if(*row_counter>max_size){
		printf("update_edge_matrix:Exceeded max size %d:%d\n",*row_counter,max_size);
	}

	std::list<shared_ptr<PhyloNode>>::iterator it;
	int currentState;
	int currentDriverID;
	for (it = thisNode->children.begin(); it != thisNode->children.end(); it++)
	{
		shared_ptr<PhyloNode> currentChild=*it;
		edges[*row_counter]=parent_counter;
		currentState=parentState;
		currentDriverID=parentDriverID;
		if(currentChild->children.size()==0){
			//terminal node
			for(const shared_ptr<Event> & pevent: currentChild->events){
				Event tmp=*pevent;
				tmp.node=*tip_counter;
				eventsOut.push_back(tmp);
				currentState=tmp.value;
				currentDriverID=tmp.driverid;
				//currentUid=tmp.uid;
			}
			edges[*row_counter+max_size]=(*tip_counter)++;
			ndivsOut[*row_counter]=currentChild->ndiv;
			tBirthOut[*row_counter]=currentChild->tbirth;
			stateOut[*row_counter]=currentState;
			driverIDOut[*row_counter]=currentDriverID;
			//uidOut[*row_counter]=currentUid;
			(*row_counter)++;
		}else{
			//internal node
			int child_counter=(*internal_counter)++;
			edges[*row_counter+max_size]=child_counter;
			for(const shared_ptr<Event> & pevent: currentChild->events){
							//printf("event:%d %d %d\n",pevent->node,pevent->value,pevent->uid);
							Event tmp=*pevent;
							tmp.node=child_counter;
							///printf("event:%d %d\n",tmp.node,tmp.value);
							eventsOut.push_back(tmp);
							currentState=tmp.value;
							currentDriverID=tmp.driverid;
							//currentUid=tmp.uid;
			}
			ndivsOut[*row_counter]=currentChild->ndiv;
			tBirthOut[*row_counter]=currentChild->tbirth;
			stateOut[*row_counter]=currentState;
			driverIDOut[*row_counter]=currentDriverID;
			//uidOut[*row_counter]=currentUid;
			(*row_counter)++;
			update_edge_matrix(edges,ndivsOut,stateOut,driverIDOut,tBirthOut,child_counter,internal_counter,tip_counter,row_counter,currentChild,max_size,eventsOut,currentState,currentDriverID);
		}
	}
}

void CellSimulation::setCompartmentInfo(){
	int ntip=ntips;
	int idx_tip=0;
	//printf("setCompartmentInfo: %d\n",ntip);
	//First clear the compartment...
	for(auto & compartment:compartments){
		compartment->clear();
	}
	setCompartmentInfoRecursively(root,0,0,&idx_tip,ntip);
	for(auto & compartment:compartments){
		//printf("#Cells in compartment %d = %d\n",compartment->id,compartment->totalpop);
		//compartment->printInfo();
		compartment->setNumNonEmptyIndices();
		compartment->setRates();

	}
}

/**
 * Assigns nodes to compartments
 */
void CellSimulation::setCompartmentInfoRecursively(shared_ptr<PhyloNode> thisNode,int compartment,int driverid,int * tip_idx,int ntip){
	int val;
	int nsub;
	if(thisNode->events.size()>0){
		for(const shared_ptr<Event> & event: thisNode->events){
			val=event->value;
			if(event->driverid!=0){
			  driverid=event->driverid;
			}
			compartment=val;
			//Could reset driverid if compartment changes..
		}
	}
	list<shared_ptr<PhyloNode>> childList=thisNode->children;
	if(childList.empty()){
		if(*tip_idx>=ntip){
			throw "CellSimulaton:setCompartmentInfoRecursively: Found too many tips!";
		}
		nsub=compartments[compartment]->nsub;
		//Note that we only addNode to a compartment if it is a tip
		//compartments[compartment]->addNode(thisNode,driverid);
		compartments[compartment]->addNode(thisNode,compartments[compartment]->getSub(driverid));
		(*tip_idx)++;
		return;
	}
	for(const shared_ptr<PhyloNode> & currentChild: childList )
	{
		setCompartmentInfoRecursively(currentChild,compartment,driverid,tip_idx,ntip);
	}
}

void CellSimulation::die(shared_ptr<PhyloNode> node){
		shared_ptr<PhyloNode> parent=(node->parent).lock();
		if(parent  && parent->children.size()==2){
			shared_ptr<PhyloNode> otherSibling=node->getSibling();
			//other_sibling->nmut+=parent->nmut;
			otherSibling->ndiv+=parent->ndiv;
			otherSibling->tbirth=parent->tbirth;
			otherSibling->parent=parent->parent;
			parent->events.insert(parent->events.end(), otherSibling->events.begin(), otherSibling->events.end());
			otherSibling->events=parent->events;
			shared_ptr<PhyloNode> grandParent=(parent->parent).lock();
			grandParent->removeChild(parent);
			grandParent->addChild(otherSibling);
			recycle(parent);
		}
		ntips--;
}

pair<shared_ptr<PhyloNode>,shared_ptr<PhyloNode>> CellSimulation::divide(shared_ptr<PhyloNode> node){
	shared_ptr<PhyloNode> & parent=node;//rename
	shared_ptr<PhyloNode> child1=getNewPhyloNode();
	shared_ptr<PhyloNode> child2=getNewPhyloNode();
	child1->parent=parent;
	child2->parent=parent;
	parent->addChild(child1);
	parent->addChild(child2);

	child1->ndiv=1;
	child2->ndiv=1;
			//Sort out class relationship here
	child1->tbirth=currentTime;
	child2->tbirth=currentTime;
	ntips++;
	return pair<shared_ptr<PhyloNode>,shared_ptr<PhyloNode>>(child1,child2);
}

/*void CellSimulation::setCurrentDriverID(int id){
	driverID=id;
}*/

/*
int CellSimulation::incrementCurrentDriverID(){
	int tmp=driverID;
	driverID++;
	return tmp;
}
*/

pair<int,double> CellSimulation::getNextDriver(){
  if(driverFitnessStack.empty()){
    //TODO throw catchable error so that the simulator just passes control back to R..
    throw "Run out of drivers... is unexpectedly empty!";
  }
  double fitness=driverFitnessStack.top();
  driverFitnessStack.pop();
  lastDriverID++;
  return pair<int,double>(lastDriverID,fitness);
}

void CellCompartment::doEvent(CellSimulation & sim){
  if(numNonEmptyCompartments==0){
    throw "Attempting to doEvent on empty compartment!";
  }
  
  /* int ct=round(sim.currentTime);
   if( abs(sim.currentTime-ct)<0.001 && ct % 20 == 0){
   printf("totals ts=%3.2f divr=%4.3f deathr=%4.3f dprob=%4.2f\n",sim.currentTime,mTotalDivRate,mTotalDeathRate,mTotalDeathRate/(mTotalDeathRate+mTotalDivRate));
   }*/
  double totalRate=mTotalDeathRate+mTotalDivRate+mTotalAsymmetricDivRate+mTotalSymmetricDivRate;
  //printf("totalrate=%7.6f..\n",totalRate);
  double rnd=totalRate*rndGen->getUniform();
  int ii=0;
  if(rnd < mTotalDeathRate){
    // death is the same for each subcomparment
    for(int i=0;i<numNonEmptyCompartments;i++){
      prob[i]=subCompartments[nonEmptyCompartmentIndices[i]].size();
    }
    int ii=rndGen->sample(numNonEmptyCompartments,prob,false);
    int i=nonEmptyCompartmentIndices[ii];
    int sz=subCompartments[i].size();
    int k=rndGen->sample(sz);
    shared_ptr<PhyloNode> node=subCompartments[i][k];
    sim.die(node);
    //remove it from the compartment...
    if(k==sz-1){
      subCompartments[i].pop_back();
    }else{
      subCompartments[i][k]=subCompartments[i][sz-1];
      subCompartments[i].pop_back();
    }
    totalpop--;
    sim.recycle(node);
    
  }else if(rnd<mTotalDeathRate+mTotalDivRate){
    //Divide
    for(int i=0;i<numNonEmptyCompartments;i++){
      ii=nonEmptyCompartmentIndices[i];
      prob[i]=(1+mFitness[ii])*subCompartments[ii].size();
    }
    ii=rndGen->sample(numNonEmptyCompartments,prob,false);
    int i=nonEmptyCompartmentIndices[ii];
    
    //now choose which sample to divide
    int k=rndGen->sample(subCompartments[i].size());
    shared_ptr<PhyloNode> parent=subCompartments[i][k];
    auto children=sim.divide(parent);
    subCompartments[i][k]=children.first;
    
    // [START]Simulate 2 way movement 
    /*		if(rndGen->getUniform()< 0.1 && sim.compartments.size()>3 && (this->id)==2){
     //int compid=this->id==2?3:2;
     int compid=3;
     sim.compartments[compid]->active=true;
     sim.compartments[compid]->addNode(children.second,0);
     //printf("Doing transition %d:%d\n",this->id,compid);
     Event event(-1,sim.currentTime,compid,0,0); 
     shared_ptr<Event> thisEvent=std::make_shared<Event>(event);
     children.second->addEvent(thisEvent);
     //subCompartments[i].push_back(children.second);
     sim.compartments[compid]->setNumNonEmptyIndices();
     sim.compartments[compid]->setRates();
    }else{
     subCompartments[i].push_back(children.second);
     totalpop++;
    }
     */
    subCompartments[i].push_back(children.second);
    totalpop++;
    //printf("totpop=%d..\n",totalpop);
    
  }else{
    //Divide
    int ii=0;
    for(int i=0;i<numNonEmptyCompartments;i++){
      ii=nonEmptyCompartmentIndices[i];
      prob[i]=subCompartments[ii].size();
    }
    ii=rndGen->sample(numNonEmptyCompartments,prob,false);
    int i=nonEmptyCompartmentIndices[ii];
    //now choose which sample to divide
    int k=rndGen->sample(subCompartments[i].size());
    shared_ptr<PhyloNode> parent=subCompartments[i][k];
    auto children=sim.divide(parent);
    if((rnd<mTotalDeathRate+mTotalDivRate+mTotalAsymmetricDivRate)){
      subCompartments[i][k]=children.first;
      int compnum=0;
      int nOut=asymmetricDifferentiation.size();
      if(nOut>1){
        //double rnd=totalRate*rndGen->getUniform();
        double pp[nOut];
        for(int j=0;j<nOut;j++){
          pp[j]=asymmetricDifferentiation[j].second;
        }
        compnum=rndGen->sample(nOut,pp,false);
      }
        int compid=asymmetricDifferentiation[compnum].first->id;
        sim.compartments[compid]->active=true;
        sim.compartments[compid]->addNode(children.second,0);
        //printf("Doing transition %d:%d\n",this->id,compid);
        Event event(-1,sim.currentTime,compid,0,0); 
        shared_ptr<Event> thisEvent=std::make_shared<Event>(event);
        children.second->addEvent(thisEvent);
        sim.compartments[compid]->setNumNonEmptyIndices();
        sim.compartments[compid]->setRates();
      //}else{
      //  throw "need to implement support for differentiation to multiple comparments.";
      //}
    }else{
        int sz=subCompartments[i].size();
        //remove it from the compartment...
        if(k==sz-1){
          subCompartments[i].pop_back();
        }else{
          subCompartments[i][k]=subCompartments[i][sz-1];
          subCompartments[i].pop_back();
        }
        
        int compnum=0;
        int nOut=symmetricDifferentiation.size();
        if(nOut>1){
          //double rnd=totalRate*rndGen->getUniform();
          double pp[nOut];
          for(int j=0;j<nOut;j++){
            pp[j]=symmetricDifferentiation[j].second;
          }
          compnum=rndGen->sample(nOut,pp,false);
        }
        int compid=symmetricDifferentiation[compnum].first->id;
        sim.compartments[compid]->active=true;
        // We're not really supporting differentiation to driver compartments.
        sim.compartments[compid]->addNode(children.first,0);
        sim.compartments[compid]->addNode(children.second,0);
        //printf("Doing transition %d:%d\n",this->id,compid);
        Event event(-1,sim.currentTime,compid,0,0);// TOD change ID 
        shared_ptr<Event> thisEvent=std::make_shared<Event>(event);
        children.first->addEvent(thisEvent);
        children.second->addEvent(thisEvent);// Potential issue here!
        sim.compartments[compid]->setNumNonEmptyIndices();
        sim.compartments[compid]->setRates();
        totalpop--;
    }
  }
  if(!atEquilibrium && totalpop>=mTargetPopSize){
      //printf("Compartment %d At Equilibrium\n",id);
    atEquilibrium=true;
  }
  setRates();
}
