/*
 * rsimpop.cpp
 *
 *  Created on: 25 Mar 2020
 *      Author: nw14
 */

#include <math.h>
#include <map>
#include <list>
#include <set>
#include <utility>
#include "CellSimulation.h"
#include "CellCompartment.h"
#include "PhyloNode.h"
#include "Event.h"
//#include "TrajectoryTimePoint.h"
extern "C" {
/**
 * Functions called by R.
 */

/*
void createTrajectoryStack(std::stack<TrajectoryTimePoint> & trajectoryStack,
                                double * trajectoryTs,  //trajectoryYear
                                int * trajectoryPop,
                                double* * trajectoryDivRate,
                                double * trajectoryDeathRate,
                                int * trajectoryCompartment,
                                int * trajectoryCompartment2,
                                double * trajectorySDivTo2,
                                double * trajectoryADivTo2,
                                int * trajectorySize
                                ){
  double pts=-1.0f;
  double tol=1e-6;
  TrajectoryTimePoint * currentTrajectory;
  for(int i=0;i<*trajectorySize;i++){
    if(trajectoryTs[i]<pts-tol){
      throw "Trajectory time stamps should be in order..";
    }
    if(trajectoryTs[i]>pts+tol){
      pts=trajectoryTs[i];
      currentTrajectory=new TrajectoryTimePoint(pts,
                                            vector<int> compartment,
                                            vector<int> pop,
                                            vector<double> divRate,
                                            vector<double> deathRate,
                                            vector<int> compartment2,
                                            vector<double> sDivRate2,
                                            vector<double> aDivRate2
                                              );    
    }
    currentTrajectory->compartment.push_back(compartment[i]);
    
  }
  
  
  
}
 */

/**
 * Creates and populates events and creates cellCompartments with the correct number of sub-compartments with the specified fitness.
 *
 * Note that the sub-compartments are actually populated in the CellSimulation constructor from the edge matrix in conjunction with events.
 */
void setSimData(vector<Event> & events,vector<shared_ptr<CellCompartment>> & cellCompartments,
		int * eventnode,double * eventts,int * eventval,int * eventdriverid,int * eventuid,
		int * nevents,int * compinfoval,double * compinfofitness,
		int * ncomp,int* compartmentsize,double * compartmentrate,double * compartmentdeathrate,
		int * compartmentval,int * ncompartment,int * driverid, double * fitnessDistribution, int * ndriverevents,int * driverFitnessSize
		){   //
  
  //printf("checkpoint 1\n");
	for(int i=0;i<*nevents;i++){
		events.push_back(Event(eventnode[i],eventts[i],eventval[i],eventdriverid[i],eventuid[i]));
	}
	//printf("checkpoint 2\n");
	std::map<int, vector<std::pair<double,int>>> fitnessByCompartment;
	for(int i=0;i<*ncomp;i++){
		auto it=fitnessByCompartment.find(compinfoval[i]);
		if(it==fitnessByCompartment.end()){
			vector<std::pair<double,int>> tmp;
			fitnessByCompartment[compinfoval[i]]=tmp;
		}
		//printf("adding %d:%d:%3.2f\n",compinfoval[i],driverid[i],compinfofitness[i]);
		fitnessByCompartment[compinfoval[i]].push_back(std::pair<double,int>(compinfofitness[i],driverid[i]));

	}
	//printf("checkpoint 3 %d\n",*ncompartment);
	for(int i=0;i<*ncompartment;i++){
	  //printf("comp %d\n",i);
		if(compartmentval[i]!=i){
			printf("Should be contiguous %d %d\n",compartmentval[i],i);
			throw "compartment should have contiguous ordered value 0..n";
		}
		
		cellCompartments.push_back(std::make_shared<CellCompartment>(CellCompartment(i,
		compartmentsize[i],
		compartmentrate[i],
    compartmentdeathrate[i],      
		fitnessByCompartment[compartmentval[i]])));
	}
	
}
  
void setMigrations(vector<shared_ptr<CellCompartment>> & cellCompartments,
                  int * c1,
                  int * c2,
                  double * arate,
                  double * srate,
                  int nmigration){
  // Add migration config to compartments
  for(int i=0;i<nmigration;i++){
  //  cellCompartments[c1[i]]->addSymmetricDifferentiationRate(cellCompartments[c2[i]],srate[i]);
  //  cellCompartments[c1[i]]->addAsymmetricDifferentiationRate(cellCompartments[c2[i]],arate[i]);
  }
}

void populateEvents(const vector<Event> & eventsOut,int * neventOut,int * eventvalOut,int * eventdriveridOut,double * eventtsOut,int * eventuidOut,int * eventnodeOut){
	int k=0;
	*neventOut=eventsOut.size();
	printf("Number of events: %d\n",*neventOut);
	for(const Event & evento : eventsOut){
				eventvalOut[k]=evento.value;
				eventdriveridOut[k]=evento.driverid;
				eventtsOut[k]=evento.timeStamp;
				eventuidOut[k]=evento.uid;
				eventnodeOut[k++]=evento.node;
	}
}

void populateCompartmentInfo(const CellSimulation & sim, int * nCompID,int * nCompPops, int * nSubCompDriverID, double * subCompFitness, int * nformedPops){
	int k=0;
  int compid;
  //printf("in popC\n");
  for(const auto & compartment : sim.compartments){
    //tuple: driverid, populationCount,fitness
    compid=compartment->id;
    vector<tuple<int,int,double>> subs=compartment->getSubInfo();
    for(const auto & sub: subs ){
      nCompID[k]=compid;
      nSubCompDriverID[k]=std::get<0>(sub);
      nCompPops[k]=std::get<1>(sub);
      subCompFitness[k++]=std::get<2>(sub);
      
    }
  }
  *nformedPops=k;
}


void initRSimPop(int * seed){
	RandomNumberGenerator::createInstance(*seed);
}



void sim_pop2(
		int * edges,
		int * ndivs,
		double * tBirth,
		int * ntips,
		int * nedge,
		int * eventnode,
		int * eventval,
		int * eventdriverid,
		double * eventts,
		int * eventuid,
		int * nevents,
		int * compartmentval,
		double *  compartmentrate,
		double * compartmentdeathrate,
		int * compartmentsize,
		int * ncompartment,
		int * compinfoval,
		double * compinfofitness,
		int * driverid,
		int *  ncomp,
		double * params,
		//int * nparams,
		double * fitnessDistribution, 
		int * nMaxPreviousDriverID, 
		int * max_size,
		int * max_events,
		double * trajectoryTs,
		int * trajectoryPop,
		double * trajectoryDivRate,
		double * trajectoryDeathRate,
		int * trajectoryCompartment,
		int * trajectorySize,
		int * c1,
		int * c2,
		double * arate,
		double * srate,
		int * driverFitnessSize,
		int * bVerbose,
		int * edgesOut,
		int * nDivsOut,
		int * statusOut,
		int * driverIDOut,
		double * tBirthOut,
		int * nedgeOut,
		int * nInternalNodeOut,
		int * eventvalOut,
		int * eventdriveridOut,
		double * eventtsOut,
		int * eventuidOut,
		int * eventnodeOut,
		int * neventOut,
		double * timestampOut,
		int * nPopSizeOut,
		int * nDriverOut,
		int * nEventsCount,
		int * ntipsOut,
		int * nCompPop,
		int * nMaxDriverIDOut,
		double * nCompFitnessOut,
		int * nCompPopsOut,
		int * nCompIDOut,
		int * nSubCompIDOut,
		int * nformedPops,
		int * ndrivereventsOut,
		int * status)
{
	try{
		//TODO:
		//RandomNumberGenerator::createInstance(-1);
		int n_days=ceil(params[0]);
		double startTime=params[3];
		int b_stop_at_pop_size=round(params[2]);
		int b_stop_if_empty=round(params[1]);
		double driverAcquisitionRate=params[4];
		if(driverAcquisitionRate>1e-5){
			throw "driverAcquisitionRate too high!";
		}
		int maxDriverCount=round(params[5]);
		int nmigration=0;//params[6];
		printf("nmigration=%d\n",nmigration);
		vector<Event> events;
		vector<shared_ptr<CellCompartment>> cellCompartments;
		//Initialise Trajectories
		
		
		setSimData(events,cellCompartments,eventnode,eventts,eventval,eventdriverid,eventuid,nevents,compinfoval,compinfofitness,
				ncomp,compartmentsize,compartmentrate,compartmentdeathrate,compartmentval,ncompartment,driverid,fitnessDistribution,nMaxPreviousDriverID,driverFitnessSize);  //
		//setMigrations(cellCompartments,c1,c2,arate,srate,nmigration);
		CellSimulation sim(edges,
				ndivs,
				tBirth,
				*nedge,
				*ntips,
				events,
				cellCompartments,
				startTime,
				driverAcquisitionRate,
				trajectoryTs,  //trajectoryYear
				trajectoryPop,
				trajectoryDivRate,
				trajectoryDeathRate,
				trajectoryCompartment,
				*trajectorySize,
				*max_size,
				fitnessDistribution,
				*nMaxPreviousDriverID,
				*driverFitnessSize,
				*bVerbose
				);
		*status=sim.run((double) n_days,(bool) b_stop_at_pop_size, (bool) b_stop_if_empty,maxDriverCount);  // (double *) nCompFitnessOut
		int nad = sim.getNumberOfAddedDrivers(); 
		//printf("nad=%d\n",nad);
		*nMaxDriverIDOut=*nMaxPreviousDriverID+nad;
		*ndrivereventsOut=*nMaxPreviousDriverID+nad;
		sim.snap();
		vector<Event> eventsOut;
		sim.populate_edge_info(edgesOut,nDivsOut,statusOut,driverIDOut,tBirthOut,*max_size,nedgeOut,nInternalNodeOut,ntipsOut,eventsOut);
		populateEvents(eventsOut,neventOut,eventvalOut,eventdriveridOut,eventtsOut,eventuidOut,eventnodeOut);
		populateCompartmentInfo(sim, nCompIDOut,nCompPopsOut, nSubCompIDOut, nCompFitnessOut, nformedPops);
		
		auto popTrace=sim.getPopulationTrace();

		int k=0;
		for(tuple<double,int,int> timepoint : popTrace){
			timestampOut[k]=std::get<0>(timepoint);//.first;
			nPopSizeOut[k]=std::get<1>(timepoint);//timepoint.second;
			nDriverOut[k++]=std::get<2>(timepoint);//
		}
		*nEventsCount=popTrace.size();
		//printf("status=%d\n",*status);
		//*status=0;

	} catch (const char* msg) {
		//cerr << msg << endl;
		printf("%s\n",msg);
		*status=-1;
	} catch (const std::exception& ex) {
		printf("%s\n",ex.what());
		*status=-1;
	}
}

/**
 * Restricts the tree (here specified by edges,muts and drivers and times) to nodes connected to the specified subtips.
 */
void sub_sample(
		int * edges,
		int * ndivs,
		double * tBirth,
		int * ntips,
		int * nedge,
		int * eventnode,
		int * eventval,
		int * eventdriverid,
		double * eventts,
		int * eventuid,
		int * nevents,
		int * compartmentval,
		double *  compartmentrate,
		double* compartmentdeathrate,
		int * compartmentsize,
		int * ncompartment,
		int * compinfoval,
		double * compinfofitness,
		int * driverid,
		int *  ncomp,
		int * subtips,
		int * nsubtips,
		double * fitnessDistribution, 
		int * ndriverevents, 
		int * max_size,
		double * trajectoryTs,
		int * trajectoryPop,
		double * trajectoryDivRate,
		double * trajectoryDeathRate,
		int * trajectoryCompartment,
		int * trajectorySize,
		int * driverFitnessSize,
		int * edgesOut,
		int * nDivsOut,
		int * stateOut,
		int * driverIDOut,
		double * tBirthOut,
		int * nedgeOut,
		int * nInternalNodeOut,
		double * eventtsOut,
		int * eventvalOut,
		int * eventdriveridOut,
		int * eventuidOut,
		int * eventnodeOut,
		int * neventOut,
		int * ntipsOut,
		double * nCompFitnessOut,
		int * nCompPopsOut,
		int * nCompIDOut,
		int * nSubCompIDOut,
		int * nformedPops,
		int * status)
{

	try{
		///RandomNumberGenerator::createInstance(1234567);
		vector<Event> events;
		vector<shared_ptr<CellCompartment>> cellCompartments;
		setSimData(events,cellCompartments,eventnode,eventts,eventval,eventdriverid,eventuid,nevents,compinfoval,compinfofitness,
				ncomp,compartmentsize,compartmentrate,compartmentdeathrate,compartmentval,ncompartment,driverid,fitnessDistribution,ndriverevents,driverFitnessSize);   //,fitnessDistribution
		
		
		
		printf("Specifying %d tips to keep!\n",*nsubtips);
		std::set<int> tipsToDelete;
		int i;
		for(i=0;i<*ntips;i++){
			tipsToDelete.insert(i+1);
		}
		int nkeep=*nsubtips;
		//printf("keeping:");
		for(i=0;i<nkeep;i++){
			tipsToDelete.erase(subtips[i]);
		}

		CellSimulation sim(edges,
				ndivs,
				tBirth,
				*nedge,
				*ntips,
				events,
				cellCompartments,
				0.0,
				0.0,
				trajectoryTs,  //trajectoryYear
				trajectoryPop,
				trajectoryDivRate,
				trajectoryDeathRate,
				trajectoryCompartment,
				*trajectorySize,
				*max_size,
				fitnessDistribution,
				*ndriverevents,
				*driverFitnessSize,
				0);
		sim.deleteTips(tipsToDelete);
		vector<Event> eventsOut;
		sim.populate_edge_info(edgesOut,nDivsOut,stateOut,driverIDOut,tBirthOut,*max_size,nedgeOut,nInternalNodeOut,ntipsOut,eventsOut);
		populateEvents(eventsOut,neventOut,eventvalOut,eventdriveridOut,eventtsOut,eventuidOut,eventnodeOut);
		populateCompartmentInfo(sim, nCompIDOut,nCompPopsOut, nSubCompIDOut, nCompFitnessOut, nformedPops);
		*status=0;
	} catch (const char* msg) {
		//cerr << msg << endl;
		printf("%s\n",msg);
		*status=1;
	} catch (const std::exception& ex) {
		printf("%s\n",ex.what());
		*status=1;
	}

}





}
