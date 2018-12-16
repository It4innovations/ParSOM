#pragma once
#include <string>
#include <cstring>
#include <stdlib.h>

enum TypOfHybrid
{
	normal, incHybrid, decHybrid
};
class ParameterParser
{
public:
	ParameterParser(void);
	~ParameterParser(void);
	bool Parse(int numberOfitems, char *charArray[]);
	char *nameOfOutputFile;
	char *nameOfConfigFile;
	char* nameOfInputDataFile;
	bool debug;
	bool writeAllresults;
	bool info;
	bool gng;
	double gngTreshHold;
	int gngMaximumParts;
	int numberOfThread;
	int janPrint; //0- vypnuto, 1-Smery VZSJ, 2- vsech 8 smeru
	double MinkowskehoNorma;
	int versionOfParallelizationThreads; //0
	bool completeThreads;
	bool cosin;
	bool gngCosin;
	bool onlyGNG;
	bool gngTestGroups;
	double gngMinkow;
	double version;
	int somHybrid;
	TypOfHybrid hybridState;
	double hybridStep;
	bool gngTree;
	bool saveOutput;
};

