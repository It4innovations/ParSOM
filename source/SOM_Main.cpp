#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include "ParameterParser.h"
#include "SaveAndLoad.h"
#include "SOMnetwork.h"

using namespace std;
//On GCC compile under -std=c++0x -pthread to get it to work.

int main(int argc, char *argv[])
{
 	  MPI_Init(&argc,&argv);

	ParameterParser inputParameters;
	inputParameters.Parse(argc,argv);
	
	if (inputParameters.onlyGNG)
	{
		SaveAndLoad fileFunction(&inputParameters);
		//InputData inputData = fileFunction.Load();
		InputData inputData;
		fileFunction.Load(&inputData);
		if (inputData.problemWithFile == 0)
		{
			SOMnetwork som(&inputData, &inputParameters);
			som.Inicialization();
			if (som.Load())
			{

				cout << "Load SOM ok" << endl;
				som.SaveGNG();
			}
			else{
				cout << "Error" << endl;
			}
		}
	}
	else
	{

		SaveAndLoad fileFunction(&inputParameters);
		InputData inputData;
		fileFunction.Load(&inputData);
		if (inputData.problemWithFile == 0)
		{
			SOMnetwork som(&inputData, &inputParameters);
			som.Inicialization();
			MPI_Barrier(MPI_COMM_WORLD);

			if (inputParameters.somHybrid == -1)
			{

				if (inputParameters.completeThreads)
					som.Run2();
				else
					som.Run();
			}
			else
			{
				som.RunHybrid();
			}
			if (inputParameters.saveOutput)
				som.Save();

			if (inputParameters.gng)
			{
				cout << "gng" << endl;
				som.SaveGNG();
			}

			MPI_Barrier(MPI_COMM_WORLD);		
			if (som.Rank==0)
			cout << "ALL OK";
		}
		else
		{
			if (inputData.problemWithFile == 1)
				cout << "Problem with Files - missing config file" << endl;
			if (inputData.problemWithFile == 2)
				cout << "Problem with Files - missing data file" << endl;
		}
	}
	MPI_Finalize();
	return 0;
}
