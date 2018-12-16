#pragma once
#include <mpi.h>
#include "SaveAndLoad.h"
#include "ParameterParser.h"
#include "Neuron.h"

typedef struct Point {
        int xval, yval;
} Point;

class SOMnetwork
{
public:
	SOMnetwork(InputData *inputDatatemp,ParameterParser *inputParameter);
	~SOMnetwork(void);
	bool Inicialization(void);
private:
	InputData *inputData;
	ParameterParser *inputParameters;
    int Size;
    
	MPI_Status stat;
	Neuron ***poleNeuronu;
	int NUMBEROFTHREAD;
	int *startPralel;
	int* endParallel;
	int* NeNullPrvkyX;
	int* NeNullPrvkyY;
	double* pomocna;
	double* resultMinParallelFor;
	Point* resultPoinParallelFor;
	double *poleVysledku;
	int * poleProUpdate;	//Pole obsahuji ID neuronu, kterou budou nasledne aktualizovany
	int pocetZaznamuPoleProUpdate; //Pocet zaznamu v poli o radek vyse
	double parametrUceni;                    
	double fiNa2;
	int pocetNeuronu1 ;
	bool InitOk;
public:
	int Rank;
	void Run(void);
	void Run2(void);
	void RunHybrid(void);
private:
	void FindBMUTask(int ID, int citacDat);
public:
	void UpdateTask(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID);
	void UpdateTask2(Point resWin,int ID);
	void UpdateTask3(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID);
	void UpdateTask4(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID);
	void UpdateTask5(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID);
	void FillArrayForUpdate(Point LevyHorniBod, Point PravyDolniBod);
	void Save();
	void HonzaExportToFigi(int** MapaPrirazeni,bool JANVYPIS);
	void HonzaExportToFigiMPI(int** MapaPrirazeni,bool JANVYPIS);
	bool SendAndReciveWeight(int Row, int Column,double * buffer);
	Point * TestFinal(void);
	double MQECount(void);
	double MQECount2(int ID);
	int * FindWinnerNeuron(void);
	bool SaveGNG(void);
	void MargeTrainingDataSpanningTree(double treshHold, int * resultTrainingSet);
	void MargeTrainingData(double treshHold,int * resultTrainingSet);
	bool Load();
};

