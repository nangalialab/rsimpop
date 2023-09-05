/*
 * CellCompartment.h
 *
 *  Created on: 20 Mar 2020
 *      Author: nw14
 */

#include <vector>
#include <iostream>
#include <map>
#include <stack>
#include "PhyloNode.h"
#include "RandomNumberGenerator.h"

#ifndef CELLCOMPARTMENT_H_
#define CELLCOMPARTMENT_H_
using namespace std;
class CellSimulation;
class CellCompartment {

private:
	vector<vector<shared_ptr<PhyloNode>>> subCompartments;
	vector<bool> bSubActive;
	RandomNumberGenerator * rndGen;
	double * prob;
	vector<int> nonEmptyCompartmentIndices;
	vector<int> emptyCompartmentIndices;
	int numNonEmptyCompartments;
	int numEmptyCompartments;
	std::map<int, int> idxByID;
	double * mfitnessDistribution=NULL;
	int numOfAddedDrivers=0;
	std::vector<double> mFitnessOut;
	std::vector<int> mID;
	



public:
	CellCompartment(int id,int targetPopSize,double divisionRate,double deathRate,std::vector<std::pair<double,int>> fitnessID);  //,
	virtual ~CellCompartment();
	int id;
	int mTargetPopSize;
	double mDivisionRate;
	double mDeathRate;
	double mTotalDivRate=0.0;
	double mTotalDeathRate=0.0;
	int mindexInFitnessDistro; 
	double alpha=0.2;
	std::vector<double> mFitness;
	int nsub;
	int totalpop;
	bool active=true;
	bool atEquilibrium=false;
	void addNode(shared_ptr<PhyloNode> node,int subid);
	double getTotalRate();
	double getTotalDivisionRate();
	double getTotalDeathRate();
	void setRates();
	void setNumNonEmptyIndices();
	int getSub(int ID);
	

	void doEvent(CellSimulation & sim);
	void checkPop();
	int die(shared_ptr<PhyloNode> node,CellSimulation & sim );
	vector<shared_ptr<PhyloNode>> getNodes();
	void clear();
	vector<pair<bool,int>> getSubCounts();
	//driverid, populationCount,fitness
	vector<tuple<int,int,double>> getSubInfo();
	void printInfo();
    //friend ostream & operator << (ostream &out, const CellCompartment &c);
	double addDriver(CellSimulation & sim,double ts);   
//	std::vector<int> subcompIndexInAddDriver; 
	
};
//ostream & operator << (ostream &out, const CellCompartment &c);
#endif /* CELLCOMPARTMENT_H_ */
