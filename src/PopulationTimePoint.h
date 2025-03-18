/*
 * PopulationTimePoint.h
 *
 *  Created on: 1 Apr 2020
 *      Author: nw14
 */

#include <string>
#include <map>
#include <vector>

#ifndef POPULATIONTIMEPOINT_H_
#define POPULATIONTIMEPOINT_H_
using namespace std;
class PopulationTimePoint {
  map<int,pair<int,int>> counts;
  
public:
  PopulationTimePoint(double ts,
                      vector<int> compartment,
                      vector<int> pop,
                      vector<int> driverCount
                      );
  double timeStamp;
  int getDriverCount(int compartment);
  int getPopulation(int compartment);
  int getTotalDriverCount();
  int getTotalPopulation();
};

#endif /* POPULATIONTIMEPOINT_H_ */
