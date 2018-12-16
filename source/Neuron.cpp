#include "Neuron.h"
#include <stdlib.h>  
#include <time.h>   
#include <math.h>

Neuron::Neuron(int x, int y, int pocetVstupu, int seedNeuron,int poziceVNeuronovesiti, double Mink)
{
	srand(seedNeuron);
	
	 LIMITITERACE = 50;
	 ZERO= 0.0000000001;
	vahy=new double[pocetVstupu];
	nonZeroVahy = new int[pocetVstupu];
    seznamNenulAll = new bool[pocetVstupu];
	//0=x, 1=y
	souradnice[0]=x;
	souradnice[1]=y;
	poziceVRikosti=poziceVNeuronovesiti;
	MinkowskehoNumber=Mink;

	numberOfdimension=pocetVstupu;
	numberOfseznamNenul=0;

	VahyX = 0;
	nonZeroCount = 0;
	VstupY = 0.0;
	
	  for (int i = 0; i < pocetVstupu; i++)
	  {
		  vahy[i]=r2();
		  VahyX += vahy[i] * vahy[i];
		  nonZeroVahy[i] = i;
		  seznamNenulAll[i]=false;
	  }
	  nonZeroCount = pocetVstupu;
}

Neuron::Neuron(int x, int y)
{
		 LIMITITERACE = 50;
		 ZERO= 0.0000000001;
		 souradnice[0]=x;
		 souradnice[1]=y;
		 numberOfseznamNenul=0;
		 VahyX = 0;
		 nonZeroCount = 0;
	     VstupY = 0.0;
}



Neuron::~Neuron(void)
{	
	delete vahy;
	delete nonZeroVahy;
	delete seznamNenulAll;

}

double Neuron::Minkowskeho(void)
{
		    double suma = 0;
            double pom = 0;
            int coutner = 0;

            for (int i = 0; i < numberOfdimension; i++)
            {
				 if (i == seznamNenul[coutner])
				 {
					 pom = vstupy[coutner] - vahy[i];
					suma += pow(abs(pom),MinkowskehoNumber);
					coutner++;
					if (coutner == numberOfseznamNenul)
                        coutner = 0;
                }
                else
                {
					 suma += pow(vahy[i],MinkowskehoNumber);
				}
            }
			return pow(suma,1/MinkowskehoNumber);// sqrt(suma);

}
double Neuron::Euklid2(void)
{
	if(cosin)
		return Cosin2();
	else
	{
	 double suma = 0;

            int pocet = numberOfseznamNenul;

                for (int i = 0; i < pocet; i++)
                {
                    suma -= 2 * vstupy[i] * vahy[seznamNenul[i]];
                }
				
            suma += VahyX;
            suma += VstupY;					

            return sqrt(suma);
	}
}
double Neuron::Cosin2()
{
		 double suma = 0;

            int pocet = numberOfseznamNenul;

                for (int i = 0; i < pocet; i++)
                {
                    suma += vstupy[i] * vahy[seznamNenul[i]];
                }
				           			

			return 1-(suma/(sqrt(VahyX)*sqrt(VstupY)));
}


bool Neuron::UpdateWeigh(int coordinatesOfWinnerX,int coordinatesOfWinnerY, double fiNa2, double learningRatio)
{
	bool jeVetsiVahaNez0 = false;
            //Potrebuji vzdalenost mezi prvky
            int rozdilX = souradnice[0] - coordinatesOfWinnerX;
            int rozdilY = souradnice[1] - coordinatesOfWinnerY;
            //Nedelam odmocninu protoze pro dalsi vypocty potrebuji mit euklida ^2
            float euklid = rozdilX * rozdilX + rozdilY * rozdilY;

            if (euklid < fiNa2)
			{
				int counter = 0;
                int maxCounte = numberOfseznamNenul;
                double omega = exp(-(double)euklid / (2 * fiNa2));
                omega *= learningRatio;     //zmena
                int pocet = numberOfdimension;
          //      double pomocnaNasob = 0;

				counter = Compute(counter, maxCounte, omega, pocet);
			}
	return jeVetsiVahaNez0;
}


int Neuron::Compute(int counter, int maxCounte, double omega, int count)
{
	VahyX = 0;

                for (int i = 0; i < numberOfseznamNenul; i++)
                {
                    seznamNenulAll[seznamNenul[i]] = true;
                }

                int tmpNonZero = 0;
                for (int i = 0; i < nonZeroCount; i++)
                {
                    int index = nonZeroVahy[i];
                    
                    if (seznamNenulAll[index]) continue;

                    double vahyI = vahy[index];

                    vahyI -= omega * vahyI;
                    if (vahyI < ZERO)
                    {
                        vahy[index] = 0;
                    }
                    else
                    {
                        VahyX += vahyI * vahyI;
                        nonZeroVahy[tmpNonZero] = index;
                        tmpNonZero++;
                        vahy[index] = vahyI;
                    }
                }


				for (int i = 0; i < numberOfseznamNenul; i++)
                {
                    int index = seznamNenul[i];
                    double vahyI = vahy[index];
                    vahyI += omega * (vstupy[i] - vahyI);
                    VahyX += vahyI * vahyI;
                    vahy[index] = vahyI;
                    seznamNenulAll[index] = false;
                    nonZeroVahy[tmpNonZero] = index;
                    tmpNonZero++;
                }
                nonZeroCount = tmpNonZero;

            
            return counter;
}

void Neuron::RecalculateWeight(void)
{
	VahyX = 0;
	for (int i = 0; i < numberOfdimension; i++)
	{
		VahyX += vahy[i] * vahy[i];
	}
}
int Neuron::GetNumberOfZero(void)
{
	 int result = 0;
            for (int i = 0; i < numberOfdimension; i++)
            {
                if (vahy[i] < ZERO)
                    result++;
            }
            return result;
}


double Neuron::r2(void)
{
	return (double)0.6+(((double)rand() / (double)RAND_MAX)/100) ;

}

double  Neuron::RidkostVypocet(int x)
{
      return ((double)pocetNulvEpochach[x]) / numberOfdimension;
}

double Neuron::EuklidMikowsky(void)
{
	if (gngCosin)
		return Cosin();
	            double suma = 0;
            double pom = 0;
           
            for (int i = 0; i < numberOfdimension; i++)
            {
                pom = vstupy[i] - vahy[i];
				suma += pow(abs(pom), gngMinkowskehoNumber);
            }
			return  pow(suma, 1 / gngMinkowskehoNumber);// sqrt(suma);
}
double Neuron::Euklid(void)
{
	double suma = 0;
	double pom = 0;

	for (int i = 0; i < numberOfdimension; i++)
	{
		pom = vstupy[i] - vahy[i];
		suma +=pom*pom;
	}
	return sqrt(suma);
}
double Neuron::Cosin()
{
		 double suma = 0;

        //    int pocet = numberOfseznamNenul;
			double vs=0;
			double va=0;
                for (int i = 0; i < numberOfdimension; i++)
                {
                    suma += vstupy[i] * vahy[i];
					va+=vahy[i]*vahy[i];
					vs+=vstupy[i] *vstupy[i];
                }

			return 1-(suma/(sqrt(va)*sqrt(vs)));
}