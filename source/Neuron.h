#pragma once
#include <vector>

class Neuron
{
public:
	Neuron(int x, int y, int pocetVstupu,int seedNeuron,int poziceVNeuronovesiti, double Mink);
	Neuron(int x, int y);
	~Neuron(void);
	double Euklid2(void);
	double Minkowskeho(void);
	double Cosin2(void);
	double MinkowskehoNumber;
	bool UpdateWeigh(int coordinatesOfWinnerX,int coordinatesOfWinnerY, double fiNa2, double learningRatio);
	int Compute(int counter, int maxCounte, double omega, int count);
	int GetNumberOfZero(void);
	int threadID;
	double RidkostVypocet(int x);
	bool cosin;
private:
	 double ZERO;
public:
	double *vahy;
public:
	int *nonZeroVahy;
	bool *seznamNenulAll;
	int nonZeroCount;
public:
	int souradnice[2];
	int poziceVRikosti;
	double VstupY;
	double *vstupy;
	int * seznamNenul;
	int  numberOfseznamNenul;
	double r2(void);
public:
		double VahyX;
		int numberOfdimension;
public:
	int LIMITITERACE;
    std::vector<int> pocetNulvEpochach;
	
	double Euklid(void);
	double Cosin(void);
	double EuklidMikowsky(void);
	void RecalculateWeight(void);
	bool gngCosin;
	double gngMinkowskehoNumber;
};

