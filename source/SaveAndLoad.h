#pragma once
#include "ParameterParser.h"
#include <mpi.h>
//#include <list>

struct InputData
{
	int numberOfRecords;
	double **data;
	int **noZero;
	int *numberOfDataInRows;
	double *SumOfInput;
	int dimensionX;
	int dimensionY;
	int numberOfRepeating;
	int numberOfInput;
	std::string nameOfFile;
	int problemWithFile;   //0- OK, 1-missing cconfig, 2 - missing data

};

class SaveAndLoad
{
private:
	ParameterParser *inputParameters;

public:
	SaveAndLoad(ParameterParser *Parameters);
	~SaveAndLoad(void);

	InputData Load(void);
	void Load(InputData *dataCollection);
};

