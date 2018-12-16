#include "ParameterParser.h"


ParameterParser::ParameterParser(void)
{
	debug=false;
	writeAllresults=false;
	info=false;
	numberOfThread=1;
	gng=false;
	MinkowskehoNorma=2;
	versionOfParallelizationThreads=0;
	completeThreads=false;
	cosin=false;
	gngTreshHold=0;
	gngMaximumParts=-1;
	 gngCosin =false;
	 onlyGNG=false;
	 gngMinkow=2;
	 version = 2;
	 gngTestGroups = false;
	 somHybrid = -1;
	 hybridState = normal;
	 hybridStep = 10;
	 gngTree = false;
	 saveOutput = true;
}


ParameterParser::~ParameterParser(void)
{
}


bool ParameterParser::Parse(int numberOfitems, char *charArray[])
{
	for (int i = 0; i < numberOfitems; i++)
	{
		if(strcmp(charArray[i],"-c")==0)
		{
			nameOfConfigFile =charArray[i+1];
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-o")==0)
		{
			nameOfOutputFile =charArray[i+1];
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-i")==0)
		{
			nameOfInputDataFile =charArray[i+1];
			i++;
			continue;
		}

		if(strcmp(charArray[i],"-debug")==0)
		{
			debug=true;

		}
		if(strcmp(charArray[i],"-all")==0)
		{
			writeAllresults=true;

		}
		if(strcmp(charArray[i],"-info")==0)
		{
			info=true;

		}
		if (strcmp(charArray[i], "-gngTree") == 0)
		{
			gngTree = true;

		}
		if(strcmp(charArray[i],"-t")==0)
		{
			numberOfThread=atoi(charArray[i+1]);
			if(numberOfThread<1)
				numberOfThread=1;
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-j")==0)
		{
			janPrint=atoi(charArray[i+1]);
			if(janPrint<1)
				janPrint=0;
			i++;
			continue;
		}

		if(strcmp(charArray[i],"-gng")==0)
		{
			gng=true;
			gngTreshHold=atof(charArray[i+1]);
			if(gngTreshHold<0)
				gngTreshHold=0;
			i++;
			continue;

		}
		if(strcmp(charArray[i],"-gngMaximumParts")==0)
		{
 			gngMaximumParts=atof(charArray[i+1]);
			if(gngMaximumParts<0)
				gngMaximumParts=-1;
			i++;
			continue;

		}
		if (strcmp(charArray[i], "-somHybrid") == 0)
		{
			somHybrid = atoi(charArray[i + 1]);
			if (somHybrid<-1)
				somHybrid = 0;

			if (somHybrid>100)
				somHybrid = 100;
			i++;
			continue;
		}

		if (strcmp(charArray[i], "-hybridInc") == 0)
		{
			hybridState = incHybrid;
		}

		if (strcmp(charArray[i], "-hybridDec") == 0)
		{
			hybridState = decHybrid;
		}

		if (strcmp(charArray[i], "-hybridStep") == 0)
		{
			hybridStep = atof(charArray[i + 1]);
			if (hybridStep<=0)
				hybridStep = 1;

			if (hybridStep>100)
				hybridStep = 100;


			i++;
			continue;
		}


		if(strcmp(charArray[i],"-cosin")==0)
		{
			cosin=true;
		}
		if (strcmp(charArray[i], "-gngcosin") == 0)
		{
			gngCosin = true;
		}
		if (strcmp(charArray[i], "-gngOnly") == 0)
		{
			onlyGNG = true;
		}
		if (strcmp(charArray[i], "-gngTest") == 0)
		{
			gngTestGroups = true;
		}
		if (strcmp(charArray[i], "-gngminkow") == 0)
		{
			gngMinkow = atof(charArray[i + 1]);
			if (gngMinkow<0)
				gngMinkow = 2;
			i++;
			continue;
		}
		if(strcmp(charArray[i],"-TC")==0)
		{
			completeThreads=true;
		}

		if(strcmp(charArray[i],"-minkow")==0)
		{
			MinkowskehoNorma=atof(charArray[i+1]);
			if(MinkowskehoNorma<0)
				MinkowskehoNorma=2;
			i++;
			continue;
		}
		if (strcmp(charArray[i], "-version") == 0)
		{
			version = atof(charArray[i + 1]);
			if (version<0)
				version = 2;
			i++;
			continue;
		}


		if(strcmp(charArray[i],"-update")==0)
		{
			versionOfParallelizationThreads=atoi(charArray[i+1]);
			if(versionOfParallelizationThreads<1)
				versionOfParallelizationThreads=0;
			i++;
			continue;
		}
		if (strcmp(charArray[i], "-noOutput") == 0)
		{
			saveOutput = false;
		}

		//return true;
	}
	return false;
}
