/*
 * PopulationTimePoint.cpp
 *
 *  Created on: 28 Sep 2023
 *      Author: nw14
 */

#include "PopulationTimePoint.h"

PopulationTimePoint::PopulationTimePoint(double ts,
                               vector<int> compartment,
                               vector<int> pop,
                               vector<int> driverCount
){
  timeStamp=ts;
  int n=compartment.size();
  for(int i=0;i<n;i++){
    counts[compartment[i]]=pair<int,int>(pop[i],driverCount[i]);
  }
}

int PopulationTimePoint::getDriverCount(int compartment){
  return counts[compartment].second;
}
int PopulationTimePoint::getPopulation(int compartment){
  return counts[compartment].first;
}
int PopulationTimePoint::getTotalDriverCount(){
  int tot=0;
  for(pair<int,pair<int,int>> count : counts){
    tot+=count.second.second;
  }
  return tot;
}
int PopulationTimePoint::getTotalPopulation(){
  int tot=0;
  for(pair<int,pair<int,int>> count : counts){
    tot+=count.second.first;
  }
  return tot;
}