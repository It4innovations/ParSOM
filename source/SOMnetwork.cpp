#include "SOMnetwork.h"
#include <time.h>       /* time */
#include <iostream>
#include <math.h>
#include <thread>
#include <fstream>
#include<sstream>
#include <omp.h>
#include <map>


using namespace std;

template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

SOMnetwork::SOMnetwork(InputData *inputDatatemp,ParameterParser *inputParameter)
{
	inputData=inputDatatemp;
	inputParameters =inputParameter;
	MPI_Comm_size(MPI_COMM_WORLD,&Size);
   /* and this processes' rank is */
	MPI_Comm_rank(MPI_COMM_WORLD,&Rank);
	NUMBEROFTHREAD=inputParameter->numberOfThread;
	InitOk=false;

}


SOMnetwork::~SOMnetwork(void)
{
	if(InitOk)
	{
	delete startPralel;
	delete endParallel;
	delete NeNullPrvkyX;
	delete NeNullPrvkyY;
	delete pomocna;
	delete resultMinParallelFor;
	delete resultPoinParallelFor;
	if(Rank==0)
	delete poleVysledku;

	for (int i = 0; i < inputData->dimensionY; i++)
	{
		for (int j = 0; j < inputData->dimensionX; j++)
		{
			if(poleNeuronu[i][j]!=NULL)
				delete poleNeuronu[i][j];
		}
		delete poleNeuronu[i];
	}
	delete poleNeuronu;
	}
			
}


bool SOMnetwork::Inicialization(void)
{
				int celkemNeuronu = inputData->dimensionX  * inputData->dimensionY;
//                Random rand = new Random(50);//unchecked((int)(DateTime.Now.Ticks >> (comm.Rank))));  //Je zaruceno ze kazdy procesor vygeneruje nahodna jina data

                /* Cilem je rozdelit vsechny neurony do jednotlivych bloku
                 * Celkovy pocet neuronu vydelim poctem procesoru (je pravdepodobne ze bude zbytek po deleni) a tak
                 * jdu od procesoru s rankem 0 a postupne pridavam procesorum jeden neuron tak abych celkovy soucet neuronu na vsech
                 * procesorech byl shodny s velikosti NS*/

                int obecnyPocetNeuronu = celkemNeuronu / Size;
                int aktualniPocetNeuronu = obecnyPocetNeuronu;
                int prebytekNeuronu = celkemNeuronu % Size;
                if (prebytekNeuronu > Rank)
                    aktualniPocetNeuronu++;
                int pocatecniPozice = 0;

                for (int i = 1; i <= Rank; i++)
                {
                    pocatecniPozice += obecnyPocetNeuronu;

                    if (prebytekNeuronu > i - 1)
                    {
                        pocatecniPozice++;
                    }
                }
           //     int celkemOpakovani = pocatecniPozice + aktualniPocetNeuronu;


                poleNeuronu = new Neuron**[inputData->dimensionY];
                for (int i = 0; i < inputData->dimensionY; i++)
                {
                    poleNeuronu[i] = new Neuron*[inputData->dimensionX ];                    
                }

                int pocetProcesoru = Size;
                int pomocnaRadek=0;
                int pocetNeuroniku = 0;

                 NeNullPrvkyX = new int[aktualniPocetNeuronu];
                 NeNullPrvkyY = new int[aktualniPocetNeuronu];
				 poleProUpdate = new int[aktualniPocetNeuronu];
				 pocetZaznamuPoleProUpdate=0;
                for (int i = 0; i < celkemNeuronu; i++)
                {
                    pomocnaRadek=i%pocetProcesoru;
                    if(pomocnaRadek==Rank)
                    {
                    
						poleNeuronu[i / inputData->dimensionX ][i %inputData->dimensionX ] = new Neuron(i % inputData->dimensionX , i / inputData->dimensionX , inputData->numberOfInput,i,pocetNeuroniku,inputParameters->MinkowskehoNorma);
						poleNeuronu[i / inputData->dimensionX ][i %inputData->dimensionX ]->threadID=pocetNeuroniku%NUMBEROFTHREAD;
						poleNeuronu[i / inputData->dimensionX ][i %inputData->dimensionX ]->cosin=inputParameters->cosin;
                        NeNullPrvkyX[pocetNeuroniku] = i % inputData->dimensionX ;
                        NeNullPrvkyY[pocetNeuroniku]=i / inputData->dimensionX ;
                        pocetNeuroniku++;
                    }else
					{
						poleNeuronu[i / inputData->dimensionX ][i %inputData->dimensionX ]=NULL;
					}
                }
                if (aktualniPocetNeuronu != pocetNeuroniku)
                {
					cout<<"error"<<aktualniPocetNeuronu <<"!="<<pocetNeuroniku<<endl;
					MPI_Abort(MPI_COMM_WORLD,1);
					return false;
                }

				  /* Samotna sit
                 * Vypocitam si vsechny konstanty ktere budu potrebovat pro vypocitani vzdalenosti a jejich dalsi upravu
                 * tyto konstatny uz pak dale jen pouzivam a neni treba je pocitat znova.*/
                 pocetNeuronu1 = pocetNeuroniku;


          //      Point LevyHorniBod = new Point(0, 0);
           //     Point PravyDolniBod = new Point(0, 0);
          

              //  bool testik=false;
                /* Zde zacina samotne prochazeni jednotlivych iteraci
                 * kazdy procesor pocita vzdalenost od vstupu k jednotlivym neuronum
                 * a zaroven hleda minimum*/
               
				MPI_Barrier(MPI_COMM_WORLD);
               

                                           
                // new double[comm.Size];
                if (Rank == 0)
                {
                    int pom = 4 * Size;
                    poleVysledku = new double[pom];
                }

               // double pomoc=0;
//                int polexy[2];// = new int[2];
//                double Poleposilani[2];// = new double[2];       
				inputData->numberOfRepeating*=inputData->numberOfRecords;
           //     int PocetEpoch = inputData->numberOfRepeating / inputData->numberOfRecords;
                
                
                //double[][] trenovaciData1 = new double[trenovaciData.Count][];
                //for (int j = 0; j < trenovaciData.Count; j++)
                //{
                //    trenovaciData1[j] = trenovaciData[j].ToArray();
                //}

                //int[][] seznamNenul1 = new int[trenovaciData.Count][];
                //for (int j = 0; j < trenovaciData.Count; j++)
                //{
                //    seznamNenul1[j] = senzamNenulData[j].ToArray();
                //}
                //double[] sumaVstupu1 = sumaVstupu.ToArray(); 
                
		
               // double  pomocnaA=0;
               // bool first = true;
              
                


                //Priprava pro paralelFor
                pomocna = new double[NUMBEROFTHREAD];
                int rangeParalel = pocetNeuronu1 / NUMBEROFTHREAD;
                int modParalel = pocetNeuronu1 % NUMBEROFTHREAD;

                startPralel = new int[NUMBEROFTHREAD];
                endParallel = new int[NUMBEROFTHREAD];

                int help = 0;
                for (int i = 0; i < NUMBEROFTHREAD; i++)
                {
                    startPralel[i] = help;
                    endParallel[i] = rangeParalel +startPralel[i];
                    if (i < modParalel)
                        endParallel[i]++;

                    help = endParallel[i];
                }

                resultMinParallelFor = new double[NUMBEROFTHREAD];
                resultPoinParallelFor = new Point[NUMBEROFTHREAD];
                for (int i = 0; i < NUMBEROFTHREAD; i++)
                {
                    resultMinParallelFor[i] = 0;
                }

                //Konec Pripravy
                
				MPI_Barrier(MPI_COMM_WORLD);
                if (Rank == 0)
                {
					cout<<"Run verze:"<<inputParameters->completeThreads<<endl;
					cout<<"Update verze:"<<inputParameters->versionOfParallelizationThreads<<endl;
					cout<<"Pocet vlaken:"<<NUMBEROFTHREAD<<endl;
					cout<<"Pocet procesu:"<<Size<<endl;
					cout<<"Velikost double:"<<sizeof(double)<<endl;
					cout<<"Velikost int:"<<sizeof(int)<<endl;
					cout<<"Velikost pameti pro jeden neuron:"<<sizeof(double)*inputData->numberOfInput+sizeof(int)*inputData->numberOfInput<<endl;
					cout<<"Mapa vytvorena. Celkem je " << inputData->numberOfRecords <<" tenovacich vzoru. Velikost mapy je:"<<inputData->dimensionX<<"x"<<inputData->dimensionY<<endl;
                }
				InitOk=true;
	return true;
}




void SOMnetwork::Run(void)
{
	if(Rank==0)
		cout<<"Spusteni funkce Run"<<endl;
	double nextIteraceEu=0;
                double couterNenul=0;
				double min = 0;
				int citacDat = 0;				
                double  pole[4];// = new double[4];
                int vyslednePole[3];// = new int[3]; 
				double lokalniMin = 0;
                int lokalniRank = 0;
				int pocetProcesoruPomocnych = Size * 4;
				int citacEpoch = 0;
			//	 bool posledniEpocha = true;
				 double MQE = 0;     
				 double completTime=0;
				
				Point LevyHorniBod;
                Point PravyDolniBod;

                bool ZacatekEpochy = true;

			//	int pocetTrenovacichDat = inputData->numberOfRecords;// trenovaciData.Count;
                
				double pocetEpoch =  inputData->numberOfRepeating / inputData->numberOfRecords;// pocetOpakovani / trenovaciData.Count;
                double PocatecniPolomer = (inputData->dimensionX < inputData->dimensionY) ? (double)inputData->dimensionX / 2 : (double)inputData->dimensionY / 2;
               // PocatecniPolomer = PocatecniPolomer / 2;
                double CasovaKonstanta = pocetEpoch / log(PocatecniPolomer);
                 parametrUceni = 1;
                double fi=0;
                 fiNa2=0;
				 double uceni = 0.5;

				Point *pozice;
			//	thread **poolOfThread=new thread *[NUMBEROFTHREAD];
								Point resTemp;
                double Start = MPI_Wtime();
				double wmtimestart = MPI_Wtime();
				double wmtimekonec = 0;
            //    Stopwatch stopky = new Stopwatch() ;
           //     stopky.Start();
			//	omp_set_num_threads(NUMBEROFTHREAD);

                for (int j = 0; j < inputData->numberOfRepeating; j++)
                {
 
                    min = 0; //pomocna = 0;                    
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]=new thread(&SOMnetwork::FindBMUTask,this,i,citacDat);
					//}
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]->join();
					//	delete poolOfThread[i];
					//}						
					
						  #pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id = omp_get_thread_num();
							if(th_id<NUMBEROFTHREAD)
							{
								//if(inputParameters->versionOfParallelizationThreads!=5)								
							//	FindBMUTask(th_id,citacDat);						

							//	int nenulX=0;
							//	int neunlY=0;
								double pomocnaTh=0;
								Point *tempPoint=&resultPoinParallelFor[th_id];
								double *tempMin=&resultMinParallelFor[th_id];
								int tempstart=startPralel[th_id];
								int tempend=endParallel[th_id];
								Neuron *tempNeuron;
											for (int i = tempstart; i < tempend; i++)
										{
										//	neunlY=NeNullPrvkyY[i];
										//	nenulX=NeNullPrvkyX[i];
											tempNeuron=poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]];
											//poleNeuronu[neunlY][nenulX];
											tempNeuron->vstupy =inputData->data[citacDat];	//trenovaciData1[citacDat];
											tempNeuron->VstupY =inputData->SumOfInput[citacDat]; //sumaVstupu1[citacDat];
											tempNeuron->seznamNenul =inputData->noZero[citacDat]; //seznamNenul1[citacDat];
											tempNeuron->numberOfseznamNenul =inputData->numberOfDataInRows[citacDat]; //seznamNenul1[citacDat];
                      
											
													if(inputParameters->MinkowskehoNorma==2)
														pomocnaTh = tempNeuron->Euklid2();
													else
														pomocnaTh = tempNeuron->Minkowskeho();
											

											//     if(j>=10)
											//       Console.WriteLine("X:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.x + " Y:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.y + " min:" + pomocna);                                  

											if (i == tempstart)
											{
												*tempMin = pomocnaTh;
												tempPoint->xval = tempNeuron->souradnice[0];
												tempPoint->yval = tempNeuron->souradnice[1];

											}
											else
											{
												if (*tempMin > pomocnaTh)
												{
													*tempMin = pomocnaTh;
												//    resultPoinParallelFor[ID] = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice;
													tempPoint->xval = tempNeuron->souradnice[0];
													tempPoint->yval = tempNeuron->souradnice[1];
												}
											}
                           
										}
							}
								
						  }



                    for (int i = 0; i < NUMBEROFTHREAD; i++)
                    {
                        if (i == 0)
                        {
                            min = resultMinParallelFor[i];
                            pozice = &resultPoinParallelFor[i];
                        }
                        else
                        {
                            if (min > resultMinParallelFor[i])
                            {
                                min = resultMinParallelFor[i];
                                pozice = &resultPoinParallelFor[i];
                            }
                            else if (min == resultMinParallelFor[i])
                            {
							//	if(inputParameters->versionOfParallelizationThreads!=5)
							//	cout<<"Rovno data:"<<citacDat<<"\tepocha:"<<j<<endl;
                            }
                            
                        }
                    }


                    /* Fukce ktera ze vsech procesoru najde nejmensi vzdalenost a tu pak odesle vsem procesorum*/
                   // Prenos vysledek = comm.Allreduce(new Prenos(min, pozice), Operace.Minimum);					
                    pole[0] = min;
					pole[1] =  pozice->xval;
					pole[2] =  pozice->yval;
                    pole[3]=nextIteraceEu;

					MPI_Gather(pole,4,MPI_DOUBLE,poleVysledku,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
					
                    //comm.GatherFlattened(pole, 0, ref poleVysledku);
                    if (Rank == 0)
                    {
                        lokalniMin = poleVysledku[0];
                        lokalniRank = 0;
                        couterNenul=poleVysledku[3];
                        for (int i = 4; i < pocetProcesoruPomocnych; i += 4)
                        {
                            if (lokalniMin > poleVysledku[i])
                            {
                                lokalniMin = poleVysledku[i];
                                lokalniRank = i;
                            }
                            couterNenul+=poleVysledku[i+3];
                        }
                        vyslednePole[0] = (int) poleVysledku[lokalniRank + 1];
                        vyslednePole[1] = (int) poleVysledku[lokalniRank + 2];
                        vyslednePole[2]= (int) couterNenul;
                        //if (couterNenul != 0)
                        //    Console.WriteLine("vetsi nez 1");
                    }
                    //comm.Broadcast(ref vyslednePole, 0);
					MPI_Bcast(vyslednePole,3,MPI_INT,0,MPI_COMM_WORLD);
                //    vysledek.souradnice.x = vyslednePole[0];
                 //   vysledek.souradnice.y = vyslednePole[1];



                    if (ZacatekEpochy)
                    {
                        ZacatekEpochy = false;
                        /*S kazdou iteraci je nutne snizit ucici pomer a take oblast ve ktere se bude projevovat zmena vahy*/
                        parametrUceni = uceni * exp(-(double)citacEpoch / CasovaKonstanta);
                        fi = PocatecniPolomer * exp(-(double)citacEpoch / CasovaKonstanta);
                        fiNa2 = fi * fi;

                        /*test*/
                        fi++;
                    }
     //               vysledek.souradnice = new Point(0, 0);
     //               Console.WriteLine("Polomer:" + fi+" Souradnice:X"+vysledek.souradnice.x+" Y:"+vysledek.souradnice.y);

                    /* Vypocet ctverce
                     * Jednodusi reseni na naprogramovani je prochazet vsechny neurony a testovat jestli nahodou nejsou ve spravne 
                     * vzdalenosti od vitezneho neuronu.
                     * Efektivnejsi reseni je je prochazet jen ctverec, ktery ohranicuje danou ktuhovou oblast.
                     * TAdy jsem mel akorat problem, ze neexistuje 2D pole  ale ze mam jen jedno rozmerne pole ve kterem kazdy neuron
                     * ma svoje souradnicem,tudiz je bylo nutne prepocitat*/
       //             vysledek.souradnice.x = 50;
        //            vysledek.souradnice.y = 50;
                    
					LevyHorniBod.xval =  vyslednePole[0]- (int)fi;
                    if (LevyHorniBod.xval < 0) LevyHorniBod.xval = 0;
                    LevyHorniBod.yval =  vyslednePole[1]- (int)fi;
                    if (LevyHorniBod.yval < 0) LevyHorniBod.yval = 0;
                    PravyDolniBod.xval =  vyslednePole[0] + (int)fi;
					if (PravyDolniBod.xval >= inputData->dimensionX ) PravyDolniBod.xval = inputData->dimensionX ;
                    PravyDolniBod.yval =  vyslednePole[1] + (int)fi;
                    if (PravyDolniBod.yval >= inputData->dimensionY) PravyDolniBod.yval = inputData->dimensionY ;

                    nextIteraceEu = 0;

                   //         for (int i = 0; i < velikostVystupu.y; i++)
                   //             for (int k = 0; k < velikostVystupu.x; k++)

                    //for (int i = LevyHorniBod.y; i < PravyDolniBod.y; i++)
                    //Parallel.For(LevyHorniBod.y,PravyDolniBod.y, i=>
                    //{
                    //    for (int k = LevyHorniBod.x; k < PravyDolniBod.x; k++)
                    //    {
                    //        if (poleNeuronu[i][k] != null)
                    //        {

                    //            poleNeuronu[i][k].AktualizujVahy(vysledek.souradnice, fiNa2, parametrUceni, true);
                    //            //                poleNeuronu[i][k].AktualizujVahyOLD(vysledek.souradnice, fiNa2, parametrUceni, true);

                    //        }
                    //    }
                    //
                    //});

					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	//poolOfThread[i]=new thread(&SOMnetwork::UpdateTask,LevyHorniBod.yval,PravyDolniBod.yval,LevyHorniBod.xval,PravyDolniBod.xval,vyslednePole[0],vyslednePole[1]);
					//	Point resTemp;
					//	resTemp.xval=vyslednePole[0];
					//	resTemp.yval=vyslednePole[1];
					//	poolOfThread[i]=new thread(&SOMnetwork::UpdateTask,this,LevyHorniBod,PravyDolniBod,resTemp,i);
					//}
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]->join();
					//	delete poolOfThread[i];
					//}
					resTemp.xval=vyslednePole[0];
					resTemp.yval=vyslednePole[1];
					if(inputParameters->versionOfParallelizationThreads==0)
					{
						#pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id2 = omp_get_thread_num();
							if(th_id2<NUMBEROFTHREAD)
							{								
								UpdateTask(LevyHorniBod,PravyDolniBod,resTemp,th_id2);
							}						
						  }
					}
					else
					if(inputParameters->versionOfParallelizationThreads==1)
					{
						FillArrayForUpdate(LevyHorniBod,PravyDolniBod);
					    #pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id2 = omp_get_thread_num();
							if(th_id2<NUMBEROFTHREAD)
							{								
								UpdateTask2(resTemp,th_id2);
							}
							//#pragma omp barrier
						  }
					}
					else if(inputParameters->versionOfParallelizationThreads==2)
					{

						  #pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id2 = omp_get_thread_num();
							if(th_id2<NUMBEROFTHREAD)
							{								
								UpdateTask3(LevyHorniBod,PravyDolniBod,resTemp,th_id2);
							}						
						  }
					}
					else if(inputParameters->versionOfParallelizationThreads==3)
					{
						 #pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id2 = omp_get_thread_num();
							if(th_id2<NUMBEROFTHREAD)
							{								
								UpdateTask4(LevyHorniBod,PravyDolniBod,resTemp,th_id2);
							}						
						  }
					}
					else if(inputParameters->versionOfParallelizationThreads==4)
					{
						#pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id2 = omp_get_thread_num();
							if(th_id2<NUMBEROFTHREAD)
							{								
								UpdateTask5(LevyHorniBod,PravyDolniBod,resTemp,th_id2);
							}						
						  }
					}else
					{

					}


                    
       //         VYHODNOCENI:
					if (((j + 1) % inputData->numberOfRecords) == 0)
                    {
         
						wmtimekonec = MPI_Wtime();
                        citacEpoch++;
                        ZacatekEpochy = true;
                      //  stopky.Stop();
 //                       posledniEpocha = false;
                    
 //                       if((citacEpoch%10==0)||(citacEpoch==pocetEpoch-1))
  //                          posledniEpocha = true;                        

                           // Console.WriteLine(comm.Rank);
                  if(inputParameters->info)
				  {
					MQE=MQECount(); 
					
				  }else
				  {
					  MQE=0;
				  }

                      //  if((citacEpoch==1))
                      //      MQE = Test(trenovaciData1,sumaVstupu1,seznamNenul1, comm, poleNeuronu);                            
                        if (Rank == 0)
                        {
                                                     
                             //   Console.WriteLine(citacEpoch + "\t" + (MQE / ((double)inputData->numberOfRecords)) + "\t" + (wmtimekonec - wmtimestart));
								cout<<citacEpoch<< "\t" << (MQE ) << "\t" << (wmtimekonec - wmtimestart)<<endl;
								completTime+=wmtimekonec - wmtimestart;
                               
                        }
                   //     PocetPurumernychNull(poleNeuronu, velikostVystupu, comm);
                   //     VypocetRidkostiDat(poleNeuronu, NeNullPrvkyX, NeNullPrvkyY, pocetNeuronu1);
                        
                        MQE = 0;                  

                       // if ((citacEpoch==10)||(citacEpoch==20)||(citacEpoch==99))
                        //if(false)
                        //{
                        //    Uloz(comm,NazevSouboru,velikostVystupu,pocetVstupu,pocetOpakovani,poleNeuronu,pocetNeuronu,((j + 1) / trenovaciData.Count));
                        //}
                      //  stopky.Start();
						wmtimestart = MPI_Wtime();
                    }
  
                    /*Prochazeni vsech neuronu- pro overeni ze funguje ctverec, pokud se ctverec zakomentuje a toto odkomentuje je to funkcni*/
            /*           for (int i = 0; i < pocetNeuronu; i++)
                       {
                           poleNeuronu[i].AktualizujVahy(vysledek.souradnice,fiNa2,parametrUceni);
                       }
               */      
                    /* overuji jestli jsem neprojel uz vsechna trenovaci data jinak jdi opet na zacatek*/
                    if (citacDat < (inputData->numberOfRecords - 1))
                        citacDat++;
                    else
                        citacDat = 0;
                }
                
				MPI_Barrier(MPI_COMM_WORLD);
				double konecT = MPI_Wtime();
 //               stopky.Stop();

                MQE =0;//= Test(trenovaciData1, sumaVstupu1, seznamNenul1, comm, poleNeuronu);			
					MQE=MQECount();
                if (Rank == 0)
                {
					double CallTime = konecT - Start;
                    cout<<"Celkovy cas :" << CallTime<<endl;
					cout<<"Cas bez testovani a vypisu :" << completTime<<endl;
					cout<<"PocetCyklu" << inputData->numberOfRepeating<<endl;
                    double casZaJedenCyklus = CallTime / inputData->numberOfRepeating;
					cout<<"Cas za jeden cyklus(sekundy):" << casZaJedenCyklus<<endl;
					cout<<"MQE posledni prumerne je:" << MQE<<endl;
                }
}

void SOMnetwork::RunHybrid(void)
{

	int HybridNumberFromPercent = ((double)inputData->numberOfRecords / 100)*inputParameters->somHybrid;
	if (HybridNumberFromPercent == 0)
		HybridNumberFromPercent = 1;
	if (Rank == 0)
		cout << "Spusteni funkce RunHybid\t Number of Hybrid:" << HybridNumberFromPercent;
	double nextIteraceEu = 0;
	double couterNenul = 0;
	double min = 0;
	int citacDat = 0;
	//int HybridNumberFromPercent = ((double)inputData->numberOfRecords / 100)*inputParameters->somHybrid;

	double * pole=new double[4 * HybridNumberFromPercent];// = new double[4];
	if (Rank == 0)
	{

		delete poleVysledku;
		poleVysledku = new double[4 * HybridNumberFromPercent*Size];
	}

	int *vyslednePole = new int[3 * HybridNumberFromPercent];// = new int[3]; 
	int hybNumberOfvyslednePole = 3 * HybridNumberFromPercent;
	double lokalniMin = 0;
	int lokalniRank = 0;
	int pocetProcesoruPomocnych = Size * 4 * HybridNumberFromPercent;
	int citacEpoch = 0;
	//	 bool posledniEpocha = true;
	double MQE = 0;
	double completTime = 0;

	Point LevyHorniBod;
	Point PravyDolniBod;

	bool ZacatekEpochy = true;

	//	int pocetTrenovacichDat = inputData->numberOfRecords;// trenovaciData.Count;

	double pocetEpoch = inputData->numberOfRepeating / inputData->numberOfRecords;// pocetOpakovani / trenovaciData.Count;
	double PocatecniPolomer = (inputData->dimensionX < inputData->dimensionY) ? (double)inputData->dimensionX / 2 : (double)inputData->dimensionY / 2;
	// PocatecniPolomer = PocatecniPolomer / 2;
	double CasovaKonstanta = pocetEpoch / log(PocatecniPolomer);
	parametrUceni = 1;
	double fi = 0;
	fiNa2 = 0;
	double uceni = 0.5;

	Point *pozice;
	//	thread **poolOfThread=new thread *[NUMBEROFTHREAD];
	Point resTemp;
	double Start = MPI_Wtime();
	double wmtimestart = MPI_Wtime();
	double wmtimekonec = 0;
	//    Stopwatch stopky = new Stopwatch() ;
	//     stopky.Start();
	//	omp_set_num_threads(NUMBEROFTHREAD);

	int counterHybrid = 0;
	int tempcounterHybrid = 0;
	bool newEpochHybrid = false;
	int actualHybridNumber = 1;
	if (inputParameters->hybridState == normal || inputParameters->hybridState == decHybrid)
		actualHybridNumber = HybridNumberFromPercent;
	double StepWithHybridNewTemp=((double)inputData->numberOfRecords / 100)*inputParameters->hybridStep;
	int StepWithHybridNew = round(StepWithHybridNewTemp);
	if (StepWithHybridNew == 0)
		StepWithHybridNew++;
	if (Rank == 0)
		cout << ", Step:" << StepWithHybridNew << endl;
	for (int j = 0; j < inputData->numberOfRepeating; j++)
	{

		min = 0; //pomocna = 0;           
		newEpochHybrid = false;
		//for (int i = 0; i < NUMBEROFTHREAD; i++)
		//{
		//	poolOfThread[i]=new thread(&SOMnetwork::FindBMUTask,this,i,citacDat);
		//}
		//for (int i = 0; i < NUMBEROFTHREAD; i++)
		//{
		//	poolOfThread[i]->join();
		//	delete poolOfThread[i];
		//}						

#pragma omp parallel num_threads(NUMBEROFTHREAD)
		{
			int th_id = omp_get_thread_num();
			if (th_id<NUMBEROFTHREAD)
			{
				//if(inputParameters->versionOfParallelizationThreads!=5)								
				//	FindBMUTask(th_id,citacDat);						

				//	int nenulX=0;
				//	int neunlY=0;
				double pomocnaTh = 0;
				Point *tempPoint = &resultPoinParallelFor[th_id];
				double *tempMin = &resultMinParallelFor[th_id];
				int tempstart = startPralel[th_id];
				int tempend = endParallel[th_id];
				Neuron *tempNeuron;
				for (int i = tempstart; i < tempend; i++)
				{
					//	neunlY=NeNullPrvkyY[i];
					//	nenulX=NeNullPrvkyX[i];
					tempNeuron = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]];
					//poleNeuronu[neunlY][nenulX];
					tempNeuron->vstupy = inputData->data[citacDat];	//trenovaciData1[citacDat];
					tempNeuron->VstupY = inputData->SumOfInput[citacDat]; //sumaVstupu1[citacDat];
					tempNeuron->seznamNenul = inputData->noZero[citacDat]; //seznamNenul1[citacDat];
					tempNeuron->numberOfseznamNenul = inputData->numberOfDataInRows[citacDat]; //seznamNenul1[citacDat];


					if (inputParameters->MinkowskehoNorma == 2)
						pomocnaTh = tempNeuron->Euklid2();
					else
						pomocnaTh = tempNeuron->Minkowskeho();


					//     if(j>=10)
					//       Console.WriteLine("X:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.x + " Y:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.y + " min:" + pomocna);                                  

					if (i == tempstart)
					{
						*tempMin = pomocnaTh;
						tempPoint->xval = tempNeuron->souradnice[0];
						tempPoint->yval = tempNeuron->souradnice[1];

					}
					else
					{
						if (*tempMin > pomocnaTh)
						{
							*tempMin = pomocnaTh;
							//    resultPoinParallelFor[ID] = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice;
							tempPoint->xval = tempNeuron->souradnice[0];
							tempPoint->yval = tempNeuron->souradnice[1];
						}
					}

				}
			}

		}



		for (int i = 0; i < NUMBEROFTHREAD; i++)
		{
			if (i == 0)
			{
				min = resultMinParallelFor[i];
				pozice = &resultPoinParallelFor[i];
			}
			else
			{
				if (min > resultMinParallelFor[i])
				{
					min = resultMinParallelFor[i];
					pozice = &resultPoinParallelFor[i];
				}
				else if (min == resultMinParallelFor[i])
				{
					//	if(inputParameters->versionOfParallelizationThreads!=5)
					//	cout<<"Rovno data:"<<citacDat<<"\tepocha:"<<j<<endl;
				}

			}
		}


		/* Fukce ktera ze vsech procesoru najde nejmensi vzdalenost a tu pak odesle vsem procesorum*/
		// Prenos vysledek = comm.Allreduce(new Prenos(min, pozice), Operace.Minimum);
		tempcounterHybrid = 4 * counterHybrid;
		pole[tempcounterHybrid] = min;
		pole[tempcounterHybrid + 1] = pozice->xval;
		pole[tempcounterHybrid + 2] = pozice->yval;
		pole[tempcounterHybrid + 3] = nextIteraceEu;

	
		counterHybrid++;
		if (counterHybrid == actualHybridNumber)
		{
		hybridRun:
			hybNumberOfvyslednePole = 3 * counterHybrid;
			pocetProcesoruPomocnych =  4 * counterHybrid;

			newEpochHybrid = true;
			MPI_Gather(pole, pocetProcesoruPomocnych, MPI_DOUBLE, poleVysledku, pocetProcesoruPomocnych, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			pocetProcesoruPomocnych *= Size;


			//comm.GatherFlattened(pole, 0, ref poleVysledku);
			if (Rank == 0)
			{
				int tempHyb = 0;
				int tempTimes = 4 * counterHybrid;
				for (int g = 0; g < counterHybrid; g++)
				{
					tempHyb = g * 4;

					lokalniMin = poleVysledku[tempHyb + 0];
					lokalniRank = 0;
					couterNenul = poleVysledku[tempHyb + 3];
					for (int i = tempTimes + tempHyb; i < pocetProcesoruPomocnych; i += tempTimes)
					{
						if (lokalniMin > poleVysledku[i])
						{
							lokalniMin = poleVysledku[i];
							lokalniRank = i;
						}
						couterNenul += poleVysledku[i + 3];
					}
					vyslednePole[g + 0] = (int)poleVysledku[lokalniRank + 1];
					vyslednePole[g + 1] = (int)poleVysledku[lokalniRank + 2];
					vyslednePole[g + 2] = (int)couterNenul;
					//if (couterNenul != 0)
					//    Console.WriteLine("vetsi nez 1");
				}
			}

			//comm.Broadcast(ref vyslednePole, 0);
			MPI_Bcast(vyslednePole, hybNumberOfvyslednePole, MPI_INT, 0, MPI_COMM_WORLD);
			//    vysledek.souradnice.x = vyslednePole[0];
			//   vysledek.souradnice.y = vyslednePole[1];


			if (ZacatekEpochy)
			{
				ZacatekEpochy = false;
				/*S kazdou iteraci je nutne snizit ucici pomer a take oblast ve ktere se bude projevovat zmena vahy*/
				parametrUceni = uceni * exp(-(double)citacEpoch / CasovaKonstanta);
				fi = PocatecniPolomer * exp(-(double)citacEpoch / CasovaKonstanta);
				fiNa2 = fi * fi;

				/*test*/
				fi++;
			}
			//               vysledek.souradnice = new Point(0, 0);
			//               Console.WriteLine("Polomer:" + fi+" Souradnice:X"+vysledek.souradnice.x+" Y:"+vysledek.souradnice.y);

			/* Vypocet ctverce
			* Jednodusi reseni na naprogramovani je prochazet vsechny neurony a testovat jestli nahodou nejsou ve spravne
			* vzdalenosti od vitezneho neuronu.
			* Efektivnejsi reseni je je prochazet jen ctverec, ktery ohranicuje danou ktuhovou oblast.
			* TAdy jsem mel akorat problem, ze neexistuje 2D pole  ale ze mam jen jedno rozmerne pole ve kterem kazdy neuron
			* ma svoje souradnicem,tudiz je bylo nutne prepocitat*/
			//             vysledek.souradnice.x = 50;
			//            vysledek.souradnice.y = 50;

			for (int i = 0; i < hybNumberOfvyslednePole; i += 3)
			{


				LevyHorniBod.xval = vyslednePole[i + 0] - (int)fi;
				if (LevyHorniBod.xval < 0) LevyHorniBod.xval = 0;
				LevyHorniBod.yval = vyslednePole[i + 1] - (int)fi;
				if (LevyHorniBod.yval < 0) LevyHorniBod.yval = 0;
				PravyDolniBod.xval = vyslednePole[i + 0] + (int)fi;
				if (PravyDolniBod.xval >= inputData->dimensionX) PravyDolniBod.xval = inputData->dimensionX;
				PravyDolniBod.yval = vyslednePole[i + 1] + (int)fi;
				if (PravyDolniBod.yval >= inputData->dimensionY) PravyDolniBod.yval = inputData->dimensionY;

				nextIteraceEu = 0;

				//         for (int i = 0; i < velikostVystupu.y; i++)
				//             for (int k = 0; k < velikostVystupu.x; k++)

				//for (int i = LevyHorniBod.y; i < PravyDolniBod.y; i++)
				//Parallel.For(LevyHorniBod.y,PravyDolniBod.y, i=>
				//{
				//    for (int k = LevyHorniBod.x; k < PravyDolniBod.x; k++)
				//    {
				//        if (poleNeuronu[i][k] != null)
				//        {

				//            poleNeuronu[i][k].AktualizujVahy(vysledek.souradnice, fiNa2, parametrUceni, true);
				//            //                poleNeuronu[i][k].AktualizujVahyOLD(vysledek.souradnice, fiNa2, parametrUceni, true);

				//        }
				//    }
				//
				//});

				//for (int i = 0; i < NUMBEROFTHREAD; i++)
				//{
				//	//poolOfThread[i]=new thread(&SOMnetwork::UpdateTask,LevyHorniBod.yval,PravyDolniBod.yval,LevyHorniBod.xval,PravyDolniBod.xval,vyslednePole[0],vyslednePole[1]);
				//	Point resTemp;
				//	resTemp.xval=vyslednePole[0];
				//	resTemp.yval=vyslednePole[1];
				//	poolOfThread[i]=new thread(&SOMnetwork::UpdateTask,this,LevyHorniBod,PravyDolniBod,resTemp,i);
				//}
				//for (int i = 0; i < NUMBEROFTHREAD; i++)
				//{
				//	poolOfThread[i]->join();
				//	delete poolOfThread[i];
				//}
				resTemp.xval = vyslednePole[i + 0];
				resTemp.yval = vyslednePole[i + 1];
				if (inputParameters->versionOfParallelizationThreads == 0)
				{
#pragma omp parallel num_threads(NUMBEROFTHREAD)
					{
						int th_id2 = omp_get_thread_num();
						if (th_id2 < NUMBEROFTHREAD)
						{
							UpdateTask(LevyHorniBod, PravyDolniBod, resTemp, th_id2);
						}
					}
				}
				else
				if (inputParameters->versionOfParallelizationThreads == 1)
				{
					FillArrayForUpdate(LevyHorniBod, PravyDolniBod);
#pragma omp parallel num_threads(NUMBEROFTHREAD)
					{
						int th_id2 = omp_get_thread_num();
						if (th_id2 < NUMBEROFTHREAD)
						{
							UpdateTask2(resTemp, th_id2);
						}
						//#pragma omp barrier
					}
				}
				else if (inputParameters->versionOfParallelizationThreads == 2)
				{

#pragma omp parallel num_threads(NUMBEROFTHREAD)
					{
						int th_id2 = omp_get_thread_num();
						if (th_id2 < NUMBEROFTHREAD)
						{
							UpdateTask3(LevyHorniBod, PravyDolniBod, resTemp, th_id2);
						}
					}
				}
				else if (inputParameters->versionOfParallelizationThreads == 3)
				{
#pragma omp parallel num_threads(NUMBEROFTHREAD)
					{
						int th_id2 = omp_get_thread_num();
						if (th_id2 < NUMBEROFTHREAD)
						{
							UpdateTask4(LevyHorniBod, PravyDolniBod, resTemp, th_id2);
						}
					}
				}
				else if (inputParameters->versionOfParallelizationThreads == 4)
				{
#pragma omp parallel num_threads(NUMBEROFTHREAD)
					{
						int th_id2 = omp_get_thread_num();
						if (th_id2 < NUMBEROFTHREAD)
						{
							UpdateTask5(LevyHorniBod, PravyDolniBod, resTemp, th_id2);
						}
					}
				}
				else
				{

				}

			}
			counterHybrid = 0;
		}
		//         VYHODNOCENI:
		if (((j + 1) % inputData->numberOfRecords) == 0)
		{
			if (newEpochHybrid==false)
				goto hybridRun;

			if (inputParameters->hybridState == incHybrid)
			{
				actualHybridNumber += StepWithHybridNew;// ((double)inputData->numberOfRecords / 100)*inputParameters->hybridStep;

				if (actualHybridNumber > HybridNumberFromPercent)
					actualHybridNumber = HybridNumberFromPercent;
			}
			else if (inputParameters->hybridState == decHybrid)
			{
				actualHybridNumber -= StepWithHybridNew;//((double)inputData->numberOfRecords / 100)*inputParameters->hybridStep;

				if (actualHybridNumber < 1)
					actualHybridNumber = 1;
			}

			wmtimekonec = MPI_Wtime();
			citacEpoch++;
			ZacatekEpochy = true;
			//  stopky.Stop();
			//                       posledniEpocha = false;

			//                       if((citacEpoch%10==0)||(citacEpoch==pocetEpoch-1))
			//                          posledniEpocha = true;                        

			// Console.WriteLine(comm.Rank);
			if (inputParameters->info)
			{
				MQE = MQECount();

			}
			else
			{
				MQE = 0;
			}

			//  if((citacEpoch==1))
			//      MQE = Test(trenovaciData1,sumaVstupu1,seznamNenul1, comm, poleNeuronu);                            
			if (Rank == 0)
			{

				//   Console.WriteLine(citacEpoch + "\t" + (MQE / ((double)inputData->numberOfRecords)) + "\t" + (wmtimekonec - wmtimestart));
				cout << citacEpoch << "\t" << (MQE) << "\t" << (wmtimekonec - wmtimestart) <<"\t"<< actualHybridNumber << endl;
				completTime += wmtimekonec - wmtimestart;

			}
			//     PocetPurumernychNull(poleNeuronu, velikostVystupu, comm);
			//     VypocetRidkostiDat(poleNeuronu, NeNullPrvkyX, NeNullPrvkyY, pocetNeuronu1);

			MQE = 0;

			// if ((citacEpoch==10)||(citacEpoch==20)||(citacEpoch==99))
			//if(false)
			//{
			//    Uloz(comm,NazevSouboru,velikostVystupu,pocetVstupu,pocetOpakovani,poleNeuronu,pocetNeuronu,((j + 1) / trenovaciData.Count));
			//}
			//  stopky.Start();
			wmtimestart = MPI_Wtime();
		}

		/*Prochazeni vsech neuronu- pro overeni ze funguje ctverec, pokud se ctverec zakomentuje a toto odkomentuje je to funkcni*/
		/*           for (int i = 0; i < pocetNeuronu; i++)
		{
		poleNeuronu[i].AktualizujVahy(vysledek.souradnice,fiNa2,parametrUceni);
		}
		*/
		/* overuji jestli jsem neprojel uz vsechna trenovaci data jinak jdi opet na zacatek*/
		if (citacDat < (inputData->numberOfRecords - 1))
			citacDat++;
		else
			citacDat = 0;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double konecT = MPI_Wtime();
	//               stopky.Stop();

	MQE = 0;//= Test(trenovaciData1, sumaVstupu1, seznamNenul1, comm, poleNeuronu);			
	MQE = MQECount();
	if (Rank == 0)
	{
		double CallTime = konecT - Start;
		cout << "Celkovy cas :" << CallTime << endl;
		cout << "Cas bez testovani a vypisu :" << completTime << endl;
		cout << "PocetCyklu" << inputData->numberOfRepeating << endl;
		double casZaJedenCyklus = CallTime / inputData->numberOfRepeating;
		cout << "Cas za jeden cyklus(sekundy):" << casZaJedenCyklus << endl;
		cout << "MQE posledni prumerne je:" << MQE << endl;
	}
	delete pole;
	delete vyslednePole;
}
void SOMnetwork::Run2(void)
{
	if(Rank==0)
		cout<<"Spusteni funkce Run"<<endl;
	double nextIteraceEu=0;
                double couterNenul=0;
				double min = 0;
				int citacDat = 0;				
                double  pole[4];// = new double[4];
                int vyslednePole[3];// = new int[3]; 
				double lokalniMin = 0;
                int lokalniRank = 0;
				int pocetProcesoruPomocnych = Size * 4;
				int citacEpoch = 0;
			//	 bool posledniEpocha = true;
				 double MQE = 0;     
				 double completTime=0;
				
				Point LevyHorniBod;
                Point PravyDolniBod;

                bool ZacatekEpochy = true;

			//	int pocetTrenovacichDat = inputData->numberOfRecords;// trenovaciData.Count;
                
				double pocetEpoch =  inputData->numberOfRepeating / inputData->numberOfRecords;// pocetOpakovani / trenovaciData.Count;
                double PocatecniPolomer = (inputData->dimensionX < inputData->dimensionY) ? (double)inputData->dimensionX / 2 : (double)inputData->dimensionY / 2;
               // PocatecniPolomer = PocatecniPolomer / 2;
                double CasovaKonstanta = pocetEpoch / log(PocatecniPolomer);
                 parametrUceni = 1;
                double fi=0;
                 fiNa2=0;
				 double uceni = 0.5;

				Point *pozice;
			//	thread **poolOfThread=new thread *[NUMBEROFTHREAD];
								Point resTemp;
                double Start = MPI_Wtime();
				double wmtimestart = MPI_Wtime();
				double wmtimekonec = 0;
            //    Stopwatch stopky = new Stopwatch() ;
           //     stopky.Start();
			//	omp_set_num_threads(NUMBEROFTHREAD);
			#pragma omp parallel num_threads(NUMBEROFTHREAD) private(citacDat,min)
				{
					int th_id = omp_get_thread_num();
					citacDat=0;
						double pomocnaTh=0;
								Point *tempPoint=&resultPoinParallelFor[th_id];
								double *tempMin=&resultMinParallelFor[th_id];
								int tempstart=startPralel[th_id];
								int tempend=endParallel[th_id];
								Neuron *tempNeuron;
                for (int j = 0; j < inputData->numberOfRepeating; j++)
                {
 
                    min = 0; //pomocna = 0;                    
							
							
							if(th_id<NUMBEROFTHREAD)
							{
								//FindBMUTask(th_id,citacDat);	

											for (int i = tempstart; i < tempend; i++)
										{
										//	neunlY=NeNullPrvkyY[i];
										//	nenulX=NeNullPrvkyX[i];
											tempNeuron=poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]];
											//poleNeuronu[neunlY][nenulX];
											tempNeuron->vstupy =inputData->data[citacDat];	//trenovaciData1[citacDat];
											tempNeuron->VstupY =inputData->SumOfInput[citacDat]; //sumaVstupu1[citacDat];
											tempNeuron->seznamNenul =inputData->noZero[citacDat]; //seznamNenul1[citacDat];
											tempNeuron->numberOfseznamNenul =inputData->numberOfDataInRows[citacDat]; //seznamNenul1[citacDat];
                      
											if(inputParameters->MinkowskehoNorma==2)
												pomocnaTh = tempNeuron->Euklid2();
											else
												pomocnaTh = tempNeuron->Minkowskeho();
					

											//     if(j>=10)
											//       Console.WriteLine("X:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.x + " Y:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.y + " min:" + pomocna);                                  

											if (i == tempstart)
											{
												*tempMin = pomocnaTh;
												tempPoint->xval = tempNeuron->souradnice[0];
												tempPoint->yval = tempNeuron->souradnice[1];

											}
											else
											{
												if (*tempMin > pomocnaTh)
												{
													*tempMin = pomocnaTh;
												//    resultPoinParallelFor[ID] = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice;
													tempPoint->xval = tempNeuron->souradnice[0];
													tempPoint->yval = tempNeuron->souradnice[1];
												}
											}
                           
										}
							}
				#pragma omp barrier
				#pragma omp master
							{
									for (int i = 0; i < NUMBEROFTHREAD; i++)
									{
										if (i == 0)
										{
											min = resultMinParallelFor[i];
											pozice = &resultPoinParallelFor[i];
										}
										else
										{
											if (min > resultMinParallelFor[i])
											{
												min = resultMinParallelFor[i];
												pozice = &resultPoinParallelFor[i];
											}
											else if (min == resultMinParallelFor[i])
											{
												//cout<<"Rovno data:"<<citacDat<<"\tepocha:"<<j<<endl;
											}
                            
										}
									}


									/* Fukce ktera ze vsech procesoru najde nejmensi vzdalenost a tu pak odesle vsem procesorum*/
								   // Prenos vysledek = comm.Allreduce(new Prenos(min, pozice), Operace.Minimum);					
									pole[0] = min;
									pole[1] =  pozice->xval;
									pole[2] =  pozice->yval;
									pole[3]=nextIteraceEu;

									MPI_Gather(pole,4,MPI_DOUBLE,poleVysledku,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
					
									//comm.GatherFlattened(pole, 0, ref poleVysledku);
									if (Rank == 0)
									{
										lokalniMin = poleVysledku[0];
										lokalniRank = 0;
										couterNenul=poleVysledku[3];
										for (int i = 4; i < pocetProcesoruPomocnych; i += 4)
										{
											if (lokalniMin > poleVysledku[i])
											{
												lokalniMin = poleVysledku[i];
												lokalniRank = i;
											}
											couterNenul+=poleVysledku[i+3];
										}
										vyslednePole[0] = (int) poleVysledku[lokalniRank + 1];
										vyslednePole[1] = (int) poleVysledku[lokalniRank + 2];
										vyslednePole[2]= (int) couterNenul;
										//if (couterNenul != 0)
										//    Console.WriteLine("vetsi nez 1");
									}
									//comm.Broadcast(ref vyslednePole, 0);
									MPI_Bcast(vyslednePole,3,MPI_INT,0,MPI_COMM_WORLD);


									if (ZacatekEpochy)
									{
										ZacatekEpochy = false;
										/*S kazdou iteraci je nutne snizit ucici pomer a take oblast ve ktere se bude projevovat zmena vahy*/
										parametrUceni = uceni * exp(-(double)citacEpoch / CasovaKonstanta);
										fi = PocatecniPolomer * exp(-(double)citacEpoch / CasovaKonstanta);
										fiNa2 = fi * fi;

										/*test*/
										fi++;
									}

                    
									LevyHorniBod.xval =  vyslednePole[0]- (int)fi;
									if (LevyHorniBod.xval < 0) LevyHorniBod.xval = 0;
									LevyHorniBod.yval =  vyslednePole[1]- (int)fi;
									if (LevyHorniBod.yval < 0) LevyHorniBod.yval = 0;
									PravyDolniBod.xval =  vyslednePole[0] + (int)fi;
									if (PravyDolniBod.xval >= inputData->dimensionX ) PravyDolniBod.xval = inputData->dimensionX ;
									PravyDolniBod.yval =  vyslednePole[1] + (int)fi;
									if (PravyDolniBod.yval >= inputData->dimensionY) PravyDolniBod.yval = inputData->dimensionY ;

									nextIteraceEu = 0;

									resTemp.xval=vyslednePole[0];
									resTemp.yval=vyslednePole[1];
							} //END MASTER
				#pragma omp barrier


					if(inputParameters->versionOfParallelizationThreads==0)
					{

							if(th_id<NUMBEROFTHREAD)
							{								
								UpdateTask(LevyHorniBod,PravyDolniBod,resTemp,th_id);
							}						
						  
					}
					else
					if(inputParameters->versionOfParallelizationThreads==1)
					{
						 #pragma omp master
						FillArrayForUpdate(LevyHorniBod,PravyDolniBod);
						#pragma omp barrier

					   
							
							if(th_id<NUMBEROFTHREAD)
							{								
								UpdateTask2(resTemp,th_id);
							}
							//#pragma omp barrier
						  
					}
					else if(inputParameters->versionOfParallelizationThreads==2)
					{


							if(th_id<NUMBEROFTHREAD)
							{								
								UpdateTask3(LevyHorniBod,PravyDolniBod,resTemp,th_id);
							}						
						  
					}
					else if(inputParameters->versionOfParallelizationThreads==3)
					{

							if(th_id<NUMBEROFTHREAD)
							{								
								UpdateTask4(LevyHorniBod,PravyDolniBod,resTemp,th_id);
							}						
						  
					}
					else if(inputParameters->versionOfParallelizationThreads==4)
					{

							if(th_id<NUMBEROFTHREAD)
							{								
								UpdateTask5(LevyHorniBod,PravyDolniBod,resTemp,th_id);
							}						
						  
					}else
					{

					}
					#pragma omp barrier
					
                    
						   //         VYHODNOCENI:
										if (((j + 1) % inputData->numberOfRecords) == 0)
										{
											#pragma omp master
											{
												wmtimekonec = MPI_Wtime();
												citacEpoch++;
												ZacatekEpochy = true;
											}

											   // Console.WriteLine(comm.Rank);
									  if(inputParameters->info)
									  {
										double locMQE=MQECount2(th_id); 
										if(th_id==0)
											MQE=locMQE;
					
									  }else
									  {
										  MQE=0;
									  }

										  //  if((citacEpoch==1))
										  //      MQE = Test(trenovaciData1,sumaVstupu1,seznamNenul1, comm, poleNeuronu);        
									   #pragma omp master
										{
											if (Rank == 0)
											{
                                                     
												 //   Console.WriteLine(citacEpoch + "\t" + (MQE / ((double)inputData->numberOfRecords)) + "\t" + (wmtimekonec - wmtimestart));
													cout<<citacEpoch<< "\t" << (MQE ) << "\t" << (wmtimekonec - wmtimestart)<<endl;
													completTime+=wmtimekonec - wmtimestart;													
											}
                        
											MQE = 0;                  
											wmtimestart = MPI_Wtime();
										}
							}
	
                    /* overuji jestli jsem neprojel uz vsechna trenovaci data jinak jdi opet na zacatek*/				
						if (citacDat < (inputData->numberOfRecords - 1))
							citacDat++;
						else
							citacDat = 0;
					
                }
				}
				MPI_Barrier(MPI_COMM_WORLD);
				double konecT = MPI_Wtime();
 //               stopky.Stop();

                MQE =0;//= Test(trenovaciData1, sumaVstupu1, seznamNenul1, comm, poleNeuronu);			
					MQE=MQECount();
                if (Rank == 0)
                {
					double CallTime = konecT - Start;
                    cout<<"Celkovy cas :" << CallTime<<endl;
					cout<<"Cas bez testovani a vypisu :" << completTime<<endl;
					cout<<"PocetCyklu" << inputData->numberOfRepeating<<endl;
                    double casZaJedenCyklus = CallTime / inputData->numberOfRepeating;
					cout<<"Cas za jeden cyklus(sekundy):" << casZaJedenCyklus<<endl;
					cout<<"MQE posledni prumerne je:" << MQE<<endl;
                }
}

void SOMnetwork::FindBMUTask(int ID, int citacDat)
{
	 for (int i = startPralel[ID]; i < endParallel[ID]; i++)
                        {
							poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->vstupy =inputData->data[citacDat];	//trenovaciData1[citacDat];
							poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->VstupY =inputData->SumOfInput[citacDat]; //sumaVstupu1[citacDat];
							poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->seznamNenul =inputData->noZero[citacDat]; //seznamNenul1[citacDat];
							poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->numberOfseznamNenul =inputData->numberOfDataInRows[citacDat]; //seznamNenul1[citacDat];
                      
							if(inputParameters->MinkowskehoNorma==2)
	                            pomocna[ID] = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->Euklid2();
							else
								pomocna[ID] = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->Minkowskeho();
					

                            //     if(j>=10)
                            //       Console.WriteLine("X:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.x + " Y:" + poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice.y + " min:" + pomocna);                                  

                            if (i == startPralel[ID])
                            {
                                resultMinParallelFor[ID] = pomocna[ID];
								resultPoinParallelFor[ID].xval = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->souradnice[0];
								resultPoinParallelFor[ID].yval = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->souradnice[1];

                            }
                            else
                            {
                                if (resultMinParallelFor[ID] > pomocna[ID])
                                {
                                    resultMinParallelFor[ID] = pomocna[ID];
                                //    resultPoinParallelFor[ID] = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]].souradnice;
									resultPoinParallelFor[ID].xval = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->souradnice[0];
									resultPoinParallelFor[ID].yval = poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->souradnice[1];
                                }
                            }
                           
                        }
}


void SOMnetwork::FillArrayForUpdate(Point LevyHorniBod, Point PravyDolniBod)
{
	pocetZaznamuPoleProUpdate=0;
	for (int Row = LevyHorniBod.yval; Row < PravyDolniBod.yval; Row++)
	{	
		for (int k = LevyHorniBod.xval; k < PravyDolniBod.xval; k++)
		{
			if (poleNeuronu[Row][k] != NULL)
			{
				poleProUpdate[pocetZaznamuPoleProUpdate]=poleNeuronu[Row][k]->poziceVRikosti;
				pocetZaznamuPoleProUpdate++;
			}
		}
	}
}
void SOMnetwork::UpdateTask2(Point resWin,int ID)
{
	int temp=0;
	for (int i = 0+ID; i < pocetZaznamuPoleProUpdate; i+=NUMBEROFTHREAD)
	{
		temp=poleProUpdate[i];
		poleNeuronu[NeNullPrvkyY[temp]][NeNullPrvkyX[temp]]->UpdateWeigh(resWin.xval,resWin.yval, fiNa2, parametrUceni);
	}
}


void SOMnetwork::UpdateTask(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID)
{
	for (int Row = LevyHorniBod.yval+ID; Row < PravyDolniBod.yval; Row+=NUMBEROFTHREAD)
	{	
		for (int k = LevyHorniBod.xval; k < PravyDolniBod.xval; k++)
		{
			if (poleNeuronu[Row][k] != NULL)
			{
				poleNeuronu[Row][k]->UpdateWeigh(resWin.xval,resWin.yval, fiNa2, parametrUceni);                                
			}
		}
	}
}

void SOMnetwork::UpdateTask3(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID)
{
	 for (int i = startPralel[ID]; i < endParallel[ID]; i++)
	 {
		 if(NeNullPrvkyY[i]>=LevyHorniBod.yval&&NeNullPrvkyY[i]<=PravyDolniBod.yval)
			 if(NeNullPrvkyX[i]>=LevyHorniBod.xval&&NeNullPrvkyX[i]<=PravyDolniBod.xval)
			 {
				 poleNeuronu[NeNullPrvkyY[i]][NeNullPrvkyX[i]]->UpdateWeigh(resWin.xval,resWin.yval, fiNa2, parametrUceni);
			 }

	 }
}

void SOMnetwork::UpdateTask4(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID)
{
	int position=0;
	int counter=0;
	for (int Row = LevyHorniBod.yval; Row < PravyDolniBod.yval; Row++)
	{	
		for (int k = LevyHorniBod.xval; k < PravyDolniBod.xval; k++)
		{
			position=((Row*inputData->dimensionX)+k);
			if(position%Size==Rank)
			{
				if(counter%NUMBEROFTHREAD==ID)
				{
					position=(position-Rank)/Size;
					
					if (poleNeuronu[NeNullPrvkyY[position]][NeNullPrvkyX[position]] != NULL)
					{
						poleNeuronu[NeNullPrvkyY[position]][NeNullPrvkyX[position]] ->UpdateWeigh(resWin.xval,resWin.yval, fiNa2, parametrUceni);                                
					}else
					{
						cout<<"Chyba pri updatu"<<endl;
					}
				}
				counter++;
			}
		}
	}
}
void SOMnetwork::UpdateTask5(Point LevyHorniBod, Point PravyDolniBod, Point resWin,int ID)
{
		for (int Row = LevyHorniBod.yval; Row < PravyDolniBod.yval; Row++)
		{	
			for (int k = LevyHorniBod.xval; k < PravyDolniBod.xval; k++)
			{
				if (poleNeuronu[Row][k] != NULL)
					if (poleNeuronu[Row][k]->threadID==ID)	
						poleNeuronu[Row][k]->UpdateWeigh(resWin.xval,resWin.yval, fiNa2, parametrUceni);   
			}
		}
}

void SOMnetwork::Save()
{
	string nameOfOutputFileCon(inputParameters->nameOfOutputFile);
	string temp=nameOfOutputFileCon;// inputData->nameOfFile;
						temp.append(NumberToString<int>(Size));
	                if(Rank==0)
                    {
						char nullChar=nameOfOutputFileCon.size();
						string localFile=temp;
						localFile.append(".Bbin");
						ofstream myFile (localFile, ios::out | ios::binary);
						myFile.write(reinterpret_cast<const char*>(&inputData->dimensionX), sizeof inputData->dimensionX);
						myFile.write(reinterpret_cast<const char*>(&inputData->dimensionY), sizeof inputData->dimensionY);
						myFile.write(reinterpret_cast<const char*>(&inputData->numberOfInput), sizeof inputData->numberOfInput);
						myFile.write(reinterpret_cast<const char*>(&inputData->numberOfRepeating), sizeof inputData->numberOfRepeating);
						myFile.write(reinterpret_cast<const char*>(&Size), sizeof Size);
						myFile.write(reinterpret_cast<const char*>(&nullChar), sizeof nullChar);
					//	cout<<nullChar<<endl;
						myFile.write(nameOfOutputFileCon.c_str(),nameOfOutputFileCon.size());	
						nullChar = strlen(inputParameters->nameOfInputDataFile);// sizeof inputParameters->nameOfInputDataFile;
						myFile.write(reinterpret_cast<const char*>(&nullChar), sizeof nullChar);
					//	cout<<nullChar<<endl;
						myFile.write(inputParameters->nameOfInputDataFile, strlen(inputParameters->nameOfInputDataFile));// (sizeof inputParameters->nameOfInputDataFile));
						myFile.write(reinterpret_cast<const char*>(&inputParameters->version), sizeof inputParameters->version);
						myFile.close();                                                
                    }

					  if (0==0)
                    {
						string localFile=temp;
						localFile.append(NumberToString<int>(Rank));
						localFile.append(".1Bbin");
						ofstream myFile (localFile, ios::out | ios::binary);
						myFile.write(reinterpret_cast<const char*>(&pocetNeuronu1), sizeof pocetNeuronu1);

                        //BinaryWriter bw = new BinaryWriter(File.Create(NazevSouboru + comm.Size + comm.Rank + ".1Bbin"));
                        //bw.Write(pocetNeuronu1);

                        for (int i = 0; i < inputData->dimensionY; i++)
                            for (int k = 0; k < inputData->dimensionX; k++)
                            {
								if (inputParameters->version >= 2)
								{
									if (poleNeuronu[i][k] == NULL)
										continue;
									myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->souradnice[0]), sizeof poleNeuronu[i][k]->souradnice[0]);
									myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->souradnice[1]), sizeof poleNeuronu[i][k]->souradnice[1]);

								}

								for (int j = 0; j < inputData->numberOfInput; j++)
                                {
									
									if (inputParameters->version < 2)
									{
										if (poleNeuronu[i][k] == NULL)
											continue;
										myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->souradnice[0]), sizeof poleNeuronu[i][k]->souradnice[0]);
										myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->souradnice[1]), sizeof poleNeuronu[i][k]->souradnice[1]);
									}
									myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->vahy[j]),sizeof poleNeuronu[i][k]->vahy[j]);                                    
                                    
                                }
                            }
							myFile.close();
                    }

					  
                //ZAPISE I TEXTAK
					  if (inputParameters->writeAllresults)
                    {
						cout<<Rank<<endl;
                        Point* ress = TestFinal();
                        int** MapaPrirazeni = new int*[inputData->dimensionY]; //Test zda to bude fungovat, jinak zakomentovat toto a odkomentovat o prad radku nize, cca 14
                        if (Rank == 0)
                        {
							string localFile=temp;
							localFile.append(".txt");
							ofstream myfile;
							myfile.open(localFile);
                            //StreamWriter tw1 = File.CreateText(NazevSouboru + comm.Size + ".txt");
							myfile<<inputData->dimensionX<<endl;
							myfile<<inputData->dimensionY<<endl;
							myfile<<inputData->numberOfInput<<endl;
							myfile<<inputData->numberOfRepeating<<endl;
							myfile<<Size<<endl;
							myfile<<nameOfOutputFileCon<<endl;
							myfile<<inputParameters->nameOfInputDataFile<<endl;                            
							myfile.close();


                            //int[][] MapaPrirazeni = new int[velikostVystupu.y][];
                            for (int i = 0; i < inputData->dimensionY; i++)
                            {
                                MapaPrirazeni[i] = new int[inputData->dimensionX];
								for (int j = 0; j < inputData->dimensionX; j++)
								{
									 MapaPrirazeni[i][j]=0;
								}
                            }

							string localFile2=temp;
							localFile2.append(" Prirazeni.txt");

							ofstream myfile2;
							myfile2.open(localFile2);

                           // StreamWriter tw2 = File.CreateText(NazevSouboru + comm.Size + " Prirazeni.txt");

							myfile2<<"ID\tY\tX"<<endl;
                            
                            for (int i = 0; i < inputData->numberOfRecords; i++)
                            {
                               // if(numberPartRecord==null)
                                    
								myfile2<<i<<"\t" << ress[i].yval << "\t" << ress[i].xval<<endl;
                             //   else
                               //     tw2.WriteLine(numberPartRecord[i] + "\t" + ress[i].y + "\t" + ress[i].x);
                                MapaPrirazeni[ress[i].yval][ress[i].xval]++;
                            }
                            myfile2.close();

							string localFile3=temp;
							localFile3.append(" PrirazeniCelkem.txt");

							ofstream myfile3;
							myfile3.open(localFile3);

                            //StreamWriter tw4 = File.CreateText(NazevSouboru + comm.Size + " PrirazeniCelkem.txt");


                            for (int i = 0; i < inputData->dimensionY; i++)
                            {
                                for (int j = 0; j < inputData->dimensionX; j++)
                                {
                                    myfile3<<MapaPrirazeni[i][j] << "\t";
                                }
                                myfile3<<endl;
                            }
                            myfile3.close();

                           // HonzaExportToFigi(poleNeuronu, velikostVystupu, comm, MapaPrirazeni,JANVYPISFIGI); //Pro HONZU a je to mozne pouzit pouze bez paralelni pocitani. Pokud i u neho, je nutne odstranit MapuPrizrazeni i z funkce:D
                        }

                        if (0 == 0)
                        {
							string localFile4=temp;
							localFile4.append(NumberToString<int>(Rank));
							localFile4.append(".txt");

							ofstream myfile4;
							myfile4.open(localFile4);
                           // StreamWriter tw = File.CreateText(NazevSouboru + comm.Size + comm.Rank + ".txt");
							myfile4<<pocetNeuronu1<<endl;
                            

                            for (int i = 0; i < inputData->dimensionY; i++)
                                for (int k = 0; k < inputData->dimensionX; k++)
                                {
                                    if (poleNeuronu[i][k] == NULL)
                                        continue;
									myfile4<<poleNeuronu[i][k]->souradnice[0];
									myfile4<<"\t";
									myfile4<<poleNeuronu[i][k]->souradnice[1];
									myfile4<<"\t";

									for (int j = 0; j < inputData->numberOfInput; j++)
                                    {
                                        myfile4<<poleNeuronu[i][k]->vahy[j] << "\t";
                                    }
                                    myfile4<<endl;
                                }
                            myfile4.close();
                        }

						string localFile5=temp;
							localFile5.append(NumberToString<int>(Rank));
							localFile5.append("-nuloveHodnoty.txt");

							ofstream myfile5;
							myfile5.open(localFile5);

                       // StreamWriter tw3 = File.CreateText(NazevSouboru + comm.Size + comm.Rank + "-nuloveHodnoty.txt");
                        int celkemPok = 0;
                        for (int i = 0; i < inputData->dimensionY; i++)
                            for (int k = 0; k < inputData->dimensionX; k++)
                            {
                                if (poleNeuronu[i][k] == NULL)
                                    continue;
									myfile5<<poleNeuronu[i][k]->souradnice[0];
									myfile5<<"\t";
									myfile5<<poleNeuronu[i][k]->souradnice[1];
									myfile5<<"\t";

									celkemPok = poleNeuronu[i][k]->pocetNulvEpochach.size();
                                for (int j = 0; j < celkemPok; j++)
                                {
                                    myfile5<<poleNeuronu[i][k]->RidkostVypocet(j) << "\t";
                                }
								myfile5<<endl;
                                
                            }
                        myfile5.close();

						if(inputParameters->janPrint>0)
						{
							if(inputParameters->janPrint>1)
								HonzaExportToFigiMPI( MapaPrirazeni, true);    
							else
								HonzaExportToFigiMPI( MapaPrirazeni, false);    


						}
						if(Rank==0)
						{
						delete ress;
						delete MapaPrirazeni;
						}
                    }
}


Point * SOMnetwork::TestFinal(void)
{
	//double MQE = 0;
            double min = 0, pomocna = 0;
			int count = inputData->numberOfRecords;
			int countY = inputData->dimensionY;
			int countX = inputData->dimensionX;
            double pole[3];
			double* poleVysledkuLoc =NULL;
			if(Rank==0)
			{
				poleVysledkuLoc=	new double[3 * Size];
			}


         
            Point* result = new Point[count];
            bool first = true;
            double lokalniMin;
            int lokalniRank;
            for (int j = 0; j < count; j++)
            //    for (int j = PRVEKTFORTEST; j < count; j+=PRVEKTFORTEST)
            {
                min = 10000; pomocna = 0;

                first = true;
               

                for (int i = 0; i < countY; i++)
                {
                    for (int k = 0; k < countX; k++)
                    {
						if (poleNeuronu[i][k] == NULL)
                            continue;

						poleNeuronu[i][k]->vstupy = inputData->data[j];// testovaciData[j];
						poleNeuronu[i][k]->VstupY = inputData->SumOfInput[j];//sumaVstupu1[j];
                        poleNeuronu[i][k]->seznamNenul = inputData->noZero[j];
						poleNeuronu[i][k]->numberOfseznamNenul =inputData->numberOfDataInRows[j]; //seznamNenul1[citacDat];
                        //            pomocna = poleNeuronu[i][k].Euklid3();

						if(inputParameters->MinkowskehoNorma==2)
							pomocna = poleNeuronu[i][k]->Euklid2();
						else
							pomocna = poleNeuronu[i][k]->Minkowskeho();

                        //if ((i == 0)&&(k==0))
                        if (first)
                        {
                            min = pomocna;
                            first = false;
                            //pozice = poleNeuronu[i][k].souradnice;
							pole[1]=poleNeuronu[i][k]->souradnice[0];
							pole[2]=poleNeuronu[i][k]->souradnice[1];
                        }
                        else
                        {
                            if (min > pomocna)
                            {
                                min = pomocna;
    							pole[1]=poleNeuronu[i][k]->souradnice[0];
								pole[2]=poleNeuronu[i][k]->souradnice[1];

                            }
                        }

                    }
                }


                pole[0] = min;// pole[1] = pozice.x; pole[2] = pozice.y;
				//MPI_Barrier(MPI_COMM_WORLD);
               // comm.GatherFlattened(pole, 0, ref poleVysledku);

				MPI_Gather(pole,3,MPI_DOUBLE,poleVysledkuLoc,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
				
              
                if (Rank == 0)
                {
                    lokalniMin = poleVysledkuLoc[0];
                    lokalniRank = 0;    
					int tempik=Size*3;
					for (int i = 3; i < tempik; i += 3)
                    {
                        if (lokalniMin > poleVysledkuLoc[i])
                        {
                            lokalniMin = poleVysledkuLoc[i];
                            lokalniRank = i;
                        }                        
                    }
                    //result[j] = new Point((int)poleVysledkuLoc[lokalniRank + 1],(int)poleVysledkuLoc[lokalniRank + 2]);
					result[j].xval=(int)poleVysledkuLoc[lokalniRank + 1];
					result[j].yval=(int)poleVysledkuLoc[lokalniRank + 2];

                    //result[j].x = (int)poleVysledku[lokalniRank + 1];
                    //result[j].y = (int)poleVysledku[lokalniRank + 2];
                    
                    //if (couterNenul != 0)
                    //    Console.WriteLine("vetsi nez 1");
                }

                

            }

	//		delete poleVysledkuLoc;
            return result;
	
}

  void SOMnetwork::HonzaExportToFigi(int** MapaPrirazeni,bool JANVYPIS)
       {
           if (Size > 1)
           {
               if (Rank == 0)
               {
                   int counter = 1;
                   int y, x,intTemp,sendCode;
              
				   MPI_Status mpistat;
				   sendCode=1;
				   MPI_Send(&sendCode,1,MPI_INT,1,111,MPI_COMM_WORLD);
                   while (counter < Size )
                   {
                      // comm.Receive(Intracommunicator.anySource, counter + 1, out x);
					   MPI_Recv(&x,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
                       if (x == -1)
                       {
                           counter++;
						   
                           continue;
                       }
                     //  comm.Receive(Intracommunicator.anySource, counter + 1, out y);
					   MPI_Recv(&y,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
                       
					   //comm.Receive(Intracommunicator.anySource, counter + 1, out helpik);
					   
					   MPI_Recv(&intTemp,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
					   poleNeuronu[y][x] = new Neuron(x,y,intTemp,1,0,2);					   
					   MPI_Recv(poleNeuronu[y][x]->vahy,intTemp,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);					   					   
                   }
               }
               else
               {
				   int tempRecv=0;
				   MPI_Status mpistat;
				   MPI_Recv(&tempRecv,1,MPI_INT,Rank-1,111,MPI_COMM_WORLD,&mpistat);
                   for (int i = 0; i < inputData->dimensionY; i++)
                       for (int k = 0; k < inputData->dimensionX; k++)
                       {
						   if (poleNeuronu[i][k] != NULL)
                           {
                               //comm.Send(k, 0, comm.Rank);
                               //comm.Send(i, 0, comm.Rank);
							   MPI_Send(&k,1,MPI_INT,0,0,MPI_COMM_WORLD);
							   MPI_Send(&i,1,MPI_INT,0,0,MPI_COMM_WORLD);
                               //comm.Send(poleNeuronu[i][k], 0, comm.Rank);

							//   MPI_Send(&(poleNeuronu[i][k]->nonZeroCount),1,MPI_INT,0,0,MPI_COMM_WORLD);
							   MPI_Send(&(poleNeuronu[i][k]->numberOfdimension),1,MPI_INT,0,0,MPI_COMM_WORLD);
							//   MPI_Send(&(poleNeuronu[i][k]->nonZeroVahy),poleNeuronu[i][k]->numberOfdimension,MPI_INT,0,0,MPI_COMM_WORLD);
							//   MPI_Send(&(poleNeuronu[i][k]->numberOfseznamNenul),1,MPI_INT,0,0,MPI_COMM_WORLD);
							//   MPI_Send(&(poleNeuronu[i][k]->seznamNenulAll),poleNeuronu[i][k]->numberOfdimension,MPI_BYTE,0,0,MPI_COMM_WORLD);
							   MPI_Send(poleNeuronu[i][k]->vahy,poleNeuronu[i][k]->numberOfdimension,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
							//   MPI_Send(&(poleNeuronu[i][k]->VahyX),1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
							  
                           }
                       }
                   //comm.Send(-1, 0, comm.Rank);
					   int intT=-1;
					   MPI_Send(&intT,1,MPI_INT,0,0,MPI_COMM_WORLD);
					   if(Rank!=Size-1)
						MPI_Send(&intT,1,MPI_INT,Rank+1,111,MPI_COMM_WORLD);

                   return;
               }
           }

           //StreamWriter tw = File.CreateText("ProFigiho.txt");
		   				string nameOfOutputFileCon1(inputParameters->nameOfOutputFile);
				string localFile=nameOfOutputFileCon1;

		   ofstream myfile5;
		   localFile.append("ProFigiho.txt");
		   myfile5.open(localFile);
           double err = 0;
           double min = -1, max = -1;
           for (int i = 0; i < inputData->dimensionY; i++)
               for (int k = 0; k < inputData->dimensionX; k++)
               {
                   if (MapaPrirazeni[i][k] == 0)
                       continue;
                   if (JANVYPIS)
                   if (((k + 1) <  inputData->dimensionX) && (i - 1) >= 0)
                   {
                       if (MapaPrirazeni[i-1][k + 1] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i - 1][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }

                   
                   if (((k + 1) < inputData->dimensionX))
                   {
                       if (MapaPrirazeni[i][k + 1] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }

                   if (JANVYPIS)
                   if (((k + 1) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i+1][k + 1] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }
                   
                   
                   if (((k) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i + 1][k] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }

               }
           string text = "";
           for (int i = 0; i < inputData->dimensionY; i++)
               for (int k = 0; k < inputData->dimensionX; k++)
               {
                   bool nasel = false;
                   if (MapaPrirazeni[i][k] == 0)
                       continue;

                   if (JANVYPIS)
                   if (((k + 1) < inputData->dimensionX) && (i - 1) >= 0)
                   {
                       if (MapaPrirazeni[i - 1][k + 1] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i - 1][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();

                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if(text.compare("0")!=0)
						   {
                            //tw.WriteLine(i + "-" + k + " " + (i - 1) + "-" + (k + 1) + " " + text);
							myfile5<<i << "-" << k << " " << (i - 1) << "-" << (k + 1) << " " << text<<endl;
						   }
                  
                           nasel = true;
                       }
                   }

                   
                   if (((k + 1) < inputData->dimensionX))
                   {
                       if (MapaPrirazeni[i][k + 1] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if (text.compare("0") != 0)
                           {
							   //tw.WriteLine(i + "-" + k + " " + (i) + "-" + (k + 1) + " " + text);
							   myfile5<<i << "-" << k << " " << (i) << "-" << (k + 1) << " " << text<<endl;
						   }                   
                           nasel = true;
                       }
                   }
                   if (JANVYPIS)
                   if (((k + 1) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i + 1][k + 1] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if (text.compare("0") != 0)
						   {
								//tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k + 1) + " " + text);
							   myfile5<<i << "-" << k << " " << (i + 1) << "-" << (k + 1) << " " << text<<endl;
						   }
                        //   else
                        //       tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k + 1));
                           nasel = true;
                       }
                   }

                   
                   if (((k) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i + 1][k] != 0)
                       {
                           poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if (text.compare("0") != 0)
						   {
								//tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k) + " " + text);
							   myfile5<<i << "-" << k << " " << (i + 1) << "-" << (k) << " " << text<<endl;
						   }
                       //    else
                       //        tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k));
                           nasel = true;
                       }
                   }
                   if (nasel == false)
                   {
                       //tw.WriteLine(i + "-" + k);
					   myfile5<<i << "-" << k<<endl;
                   }
               }
           //tw.Close();
			   myfile5.close();

	   }

	   double SOMnetwork::MQECount(void)
	   {

			//	thread **poolOfThread=new thread *[NUMBEROFTHREAD];
				double min,temp;
				double MQE=0;
				for (int j = 0; j < inputData->numberOfRecords; j++)
                {
 

                    min = 0; //pomocna = 0;                    
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]=new thread(&SOMnetwork::FindBMUTask,this,i,j);
					//}
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]->join();
					//	delete poolOfThread[i];
					//}
						
						  #pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id = omp_get_thread_num();
							if(th_id<NUMBEROFTHREAD)
								FindBMUTask(th_id,j);						
						  }


                    for (int i = 0; i < NUMBEROFTHREAD; i++)
                    {
                        if (i == 0)
                        {
                            min = resultMinParallelFor[i];                            
                        }
                        else
                        {
                            if (min > resultMinParallelFor[i])
                            {
                                min = resultMinParallelFor[i];
                            }
                            else if (min == resultMinParallelFor[i])
                            {
								cout<<"Rovno"<<endl;
                            }
                            
                        }
                    }

					MPI_Reduce(&min,&temp,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
					if(Rank==0)
					{
						//cout<<j<<"\t"<<MQE<<"\t"<<temp<<endl;
						MQE+=temp;
						
						//return NULL;
					}
				}
					
				return MQE/inputData->numberOfRecords;
	   }
	    double SOMnetwork::MQECount2(int ID)
	   {

			//	thread **poolOfThread=new thread *[NUMBEROFTHREAD];
				double min,temp;
				double MQE=0;				
				for (int j = 0; j < inputData->numberOfRecords; j++)
                {
 

                    min = 0; //pomocna = 0;                    
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]=new thread(&SOMnetwork::FindBMUTask,this,i,j);
					//}
					//for (int i = 0; i < NUMBEROFTHREAD; i++)
					//{
					//	poolOfThread[i]->join();
					//	delete poolOfThread[i];
					//}
						
							if(ID<NUMBEROFTHREAD)
								FindBMUTask(ID,j);						
						  
			#pragma omp barrier
			#pragma omp master
							{
								for (int i = 0; i < NUMBEROFTHREAD; i++)
								{
									if (i == 0)
									{
										min = resultMinParallelFor[i];      
										
									}
									else
									{
										if (min > resultMinParallelFor[i])
										{
											min = resultMinParallelFor[i];
										}
										else if (min == resultMinParallelFor[i])
										{
											cout<<"Rovno"<<endl;
										}
										resultMinParallelFor[i]=-1;
									}
								}

								MPI_Reduce(&min,&temp,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
								if(Rank==0)
								{
									//cout<<j<<"\t"<<MQE<<"\t"<<temp<<endl;
									MQE+=temp;
						
									//return NULL;
								}
							}
							#pragma omp barrier
				}
				#pragma omp barrier	
				return MQE/inputData->numberOfRecords;
	   }


	   int * SOMnetwork::FindWinnerNeuron(void)
	   {
		   double min;
				int VyslednyXY[2];
				double ResultMin=0;
				MPI_Status status;
				int* resultTraningSet = new int[inputData->numberOfRecords];
				for (int j = 0; j < inputData->numberOfRecords; j++)
                {
                    min = 0; //pomocna = 0;                    

						  #pragma omp parallel num_threads(NUMBEROFTHREAD)
						  {
							int th_id = omp_get_thread_num();
							if(th_id<NUMBEROFTHREAD)
								FindBMUTask(th_id,j);						
						  }

			     for (int i = 0; i < NUMBEROFTHREAD; i++)
                    {
                        if (i == 0)
                        {
                            min = resultMinParallelFor[i];
                         //   pozice = &resultPoinParallelFor[i];
							VyslednyXY[0]=resultPoinParallelFor[i].yval;
							VyslednyXY[1]=resultPoinParallelFor[i].xval;
                        }
                        else
                        {
                            if (min > resultMinParallelFor[i])
                            {
                                min = resultMinParallelFor[i];
                                //pozice = &resultPoinParallelFor[i];
								VyslednyXY[0]=resultPoinParallelFor[i].yval;
								VyslednyXY[1]=resultPoinParallelFor[i].xval;
                            }
                            else if (min == resultMinParallelFor[i])
                            {
								cout<<"Rovno:"<<j<<endl;
                            }
                            
                        }
                    }

				 MPI_Allreduce(&min,&ResultMin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
				 if (ResultMin == min)
                {
                    if (Rank != 0)
						MPI_Send(&VyslednyXY,2,MPI_INT, 0, 10,MPI_COMM_WORLD);

                }
                else if (Rank == 0)
                {
					MPI_Recv(&VyslednyXY,2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                }
				 resultTraningSet[j] = VyslednyXY[0] * inputData->dimensionX + VyslednyXY[1];

		   
				}
				return resultTraningSet;
	   }

	   bool SOMnetwork::SaveGNG(void)
	   {
		   if (inputParameters->gngCosin)
			   cout << "GNG Cosin" << endl;
		   else
		   {
			   cout << "GNG " << inputParameters->gngMinkow<<endl;
		   }
			   

		    int *VyslednyResult = FindWinnerNeuron();
			double StartT = MPI_Wtime();
			if (inputParameters->gngTree)
			{
				cout << "GNG Tree" << endl;
				MargeTrainingDataSpanningTree(inputParameters->gngTreshHold, VyslednyResult);
			}
			else
			{
				cout << "GNG Trash" << endl;
				MargeTrainingData(inputParameters->gngTreshHold, VyslednyResult);
			}
			double konecT = MPI_Wtime();

			if(Rank==0)
			{
				cout << "GNGGroupsTime:" << (konecT - StartT) << "GNGEndTime"<<endl;
				int tempCounter = 0;
				std::map<int,int> temp;
								//Dictionary<int, int> temp = new Dictionary<int, int>();
								for (int i = 0; i < inputData->numberOfRecords; i++)
								{
									//if (!temp.TryGetValue(VyslednyResult[i], out res))
									//int x=0;
									if(temp.count(VyslednyResult[i])==0)
									{
										
										temp.insert(std::pair<int,int>(VyslednyResult[i], tempCounter));
										tempCounter++;
									}
								}
				string nameOfOutputFileCon(inputParameters->nameOfOutputFile);
				string localFile=nameOfOutputFileCon;
				localFile.append(".conf");
				ofstream myfile;
				myfile.open(localFile);

				myfile<<inputData->numberOfRecords<<endl;
				myfile<<inputData->numberOfInput<<endl;
				myfile<<inputParameters->nameOfInputDataFile<<endl;

				localFile=nameOfOutputFileCon;
				localFile.append(".VRbin");
				myfile<<localFile<<endl;

				localFile=nameOfOutputFileCon;
				myfile<<localFile<<Size<<endl;

				myfile<<temp.size()<<endl;
				myfile.close();

				localFile=nameOfOutputFileCon;
				localFile.append(".VRbin");
				ofstream myFile1 (localFile, ios::out | ios::binary);
				for (int i = 0; i < inputData->numberOfRecords; i++)
				{
					myFile1.write(reinterpret_cast<const char*>(&(temp[VyslednyResult[i]])),sizeof (temp[VyslednyResult[i]]));
				}
				myFile1.close();

				if (inputParameters->gngTestGroups)
				{
					myfile.open(nameOfOutputFileCon + "-Groups-" + NumberToString(inputParameters->gngMaximumParts)+".txt");
					for (int i = 0; i < inputData->numberOfRecords; i++)
					{
						myfile << VyslednyResult[i] << endl;
					}
					myfile.close();
				}

			}
			delete VyslednyResult;
		   return false;
	   }

	   void SOMnetwork::HonzaExportToFigiMPI(int** MapaPrirazeni,bool JANVYPIS)
       {
           if (Size > 1)
           {
               if (Rank != 0)
               {
				   int temprecv[2]; // prvni y druhe x
				   
				   MPI_Status mpistat;
				   while(true)
				   {
					   MPI_Recv(&temprecv,2,MPI_INT,0,0,MPI_COMM_WORLD,&mpistat);
					   if(temprecv[0]==-1)
						   break;

							   if (poleNeuronu[temprecv[0]][temprecv[1]] != NULL)
							   {
								   MPI_Send(poleNeuronu[temprecv[0]][temprecv[1]]->vahy,poleNeuronu[temprecv[0]][temprecv[1]]->numberOfdimension,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
								   }
							   else
							   {
								   cout<<"Error in function HonzaExporttoFigiMPI"<<endl;
							   }
						   
				   }

                   return;
               }
           }
		   else
		   {
			   HonzaExportToFigi(MapaPrirazeni,JANVYPIS);
			   return;
		   }

			double *temparray=new double[inputData->numberOfInput];
			Neuron * tempNeuron= new Neuron(-1,-1,inputData->numberOfInput,0,0,0);
		   				string nameOfOutputFileCon1(inputParameters->nameOfOutputFile);
				string localFile=nameOfOutputFileCon1;

		   ofstream myfile5;
		   localFile.append("ProFigiho.txt");
		   myfile5.open(localFile);
           double err = 0;
           double min = -1, max = -1;
           for (int i = 0; i < inputData->dimensionY; i++)
               for (int k = 0; k < inputData->dimensionX; k++)
               {
                   if (MapaPrirazeni[i][k] == 0)
                       continue;
				   if( poleNeuronu[i][k]==NULL)
				   {
					   SendAndReciveWeight(i,k,tempNeuron->vahy);
					   poleNeuronu[i][k]=tempNeuron;
				   }

                   if (JANVYPIS)
                   if (((k + 1) <  inputData->dimensionX) && (i - 1) >= 0)
                   {
                       if (MapaPrirazeni[i-1][k + 1] != 0)
                       {
						   if(SendAndReciveWeight(i-1,k + 1,temparray))
							poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i - 1][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i - 1][k + 1]->vahy;

                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }

                   
                   if (((k + 1) < inputData->dimensionX))
                   {
                       if (MapaPrirazeni[i][k + 1] != 0)
                       {
                           if(SendAndReciveWeight(i,k + 1,temparray))

							  poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i ][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }

                   if (JANVYPIS)
                   if (((k + 1) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i+1][k + 1] != 0)
                       {
                           if(SendAndReciveWeight(i+1,k + 1,temparray))
							    poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i + 1][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }
                   
                   
                   if (((k) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i + 1][k] != 0)
                       {
                          if( SendAndReciveWeight(i+1,k,temparray))
							poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i + 1][k]->vahy;
						  else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k]->vahy;
                           err = poleNeuronu[i][k]->Euklid();
                           if (min == -1)
                               min = err;
                           if (max == -1)
                               max = err;
                           if (min > err)
                               min = err;
                           if (max < err)
                               max = err;
                       }
                   }

				   if( poleNeuronu[i][k]->souradnice[0]==-1)
					     poleNeuronu[i][k]=NULL;
               }
		//	   cout<<"max:"<<max<<" min:"<<min<<endl;
           string text = "";
           for (int i = 0; i < inputData->dimensionY; i++)
               for (int k = 0; k < inputData->dimensionX; k++)
               {
                   bool nasel = false;
                   if (MapaPrirazeni[i][k] == 0)
                       continue;

				if( poleNeuronu[i][k]==NULL)
				   {
					   SendAndReciveWeight(i,k,tempNeuron->vahy);
					   poleNeuronu[i][k]=tempNeuron;
				   }

                   if (JANVYPIS)
                   if (((k + 1) < inputData->dimensionX) && (i - 1) >= 0)
                   {
                       if (MapaPrirazeni[i - 1][k + 1] != 0)
                       {
                           if(SendAndReciveWeight(i-1,k + 1,temparray))
								poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i - 1][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i - 1][k + 1]->vahy;
                           err = poleNeuronu[i][k]->Euklid();

                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if(text.compare("0")!=0)
						   {
                            //tw.WriteLine(i + "-" + k + " " + (i - 1) + "-" + (k + 1) + " " + text);
							myfile5<<i << "-" << k << " " << (i - 1) << "-" << (k + 1) << " " << text<<endl;
						   }
                  
                           nasel = true;
                       }
                   }

                   
                   if (((k + 1) < inputData->dimensionX))
                   {
                       if (MapaPrirazeni[i][k + 1] != 0)
                       {
                           if(SendAndReciveWeight(i,k + 1,temparray))

								 poleNeuronu[i][k]->vstupy =temparray;//poleNeuronu[i][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i][k + 1]->vahy;

                           err = poleNeuronu[i][k]->Euklid();
                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if (text.compare("0") != 0)
                           {
							   //tw.WriteLine(i + "-" + k + " " + (i) + "-" + (k + 1) + " " + text);
							   myfile5<<i << "-" << k << " " << (i) << "-" << (k + 1) << " " << text<<endl;
						   }                   
                           nasel = true;
                       }
                   }
                   if (JANVYPIS)
                   if (((k + 1) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i + 1][k + 1] != 0)
                       {
                           if(SendAndReciveWeight(i+1,k + 1,temparray))
							   poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i + 1][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k + 1]->vahy;

                           err = poleNeuronu[i][k]->Euklid();
                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if (text.compare("0") != 0)
						   {
								//tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k + 1) + " " + text);
							   myfile5<<i << "-" << k << " " << (i + 1) << "-" << (k + 1) << " " << text<<endl;
						   }
                        //   else
                        //       tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k + 1));
                           nasel = true;
                       }
                   }

                   
                   if (((k) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       if (MapaPrirazeni[i + 1][k] != 0)
                       {
                           if(SendAndReciveWeight(i+1,k,temparray))

							 poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i + 1][k]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k]->vahy;

                           err = poleNeuronu[i][k]->Euklid();
                           text = to_string(1 - ((err - min) / (double)(max - min)));//.ToString().Replace(',', '.');
                           if (text.compare("0") != 0)
						   {
								//tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k) + " " + text);
							   myfile5<<i << "-" << k << " " << (i + 1) << "-" << (k) << " " << text<<endl;
						   }
                       //    else
                       //        tw.WriteLine(i + "-" + k + " " + (i + 1) + "-" + (k));
                           nasel = true;
                       }
                   }
                   if (nasel == false)
                   {
                       //tw.WriteLine(i + "-" + k);
					   myfile5<<i << "-" << k<<endl;
                   }
				   				   if( poleNeuronu[i][k]->souradnice[0]==-1)
					     poleNeuronu[i][k]=NULL;
               }
           //tw.Close();
			   myfile5.close();
			   cout<<"File close"<<endl;
			   int arrayForPosition[2];
			   arrayForPosition[0]=-1;
		   arrayForPosition[1]=-1;
			   for (int i = 1; i < Size; i++)
			   {
				   MPI_Send(arrayForPosition,2,MPI_INT,i,0, MPI_COMM_WORLD);
			   }
			   delete temparray;
			delete tempNeuron;
	   }
	   


	   bool SOMnetwork::SendAndReciveWeight(int Row, int Column,double * buffer)
	   {
		   MPI_Status mpistat;
		   int  procesForSend= (((Row)*inputData->dimensionX)+(Column))%Size;
		//   cout<<"Row:"<<Row<<" Colum:"<<Column<<endl;
		   if(procesForSend==0)
		   {
			   return false;
		   }else
		   {
				   int arrayForPosition[2];
				   arrayForPosition[0]=Row;
				   arrayForPosition[1]=Column;
				   MPI_Send(arrayForPosition,2,MPI_INT,procesForSend,0, MPI_COMM_WORLD);

				   MPI_Recv(buffer,inputData->numberOfInput,MPI_DOUBLE,procesForSend,0,MPI_COMM_WORLD,&mpistat);
		   }
		   return true;
	   }

	   void SOMnetwork::MargeTrainingDataSpanningTree(double treshHold, int * resultTrainingSet)
	   {
		   if (Size > 1)
		   {
			   if (Rank != 0)
			   {
				   int temprecv[2]; // prvni y druhe x

				   MPI_Status mpistat;
				   while (true)
				   {
					   MPI_Recv(&temprecv, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistat);
					   if (temprecv[0] == -1)
						   break;

					   if (poleNeuronu[temprecv[0]][temprecv[1]] != NULL)
					   {
						   MPI_Send(poleNeuronu[temprecv[0]][temprecv[1]]->vahy, poleNeuronu[temprecv[0]][temprecv[1]]->numberOfdimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
					   }
					   else
					   {
						   cout << "Error in the function MargeTrainingData" << endl;
					   }

				   }

				   return;
			   }
		   }

		   double *temparray = new double[inputData->numberOfInput];
		   Neuron * tempNeuron = new Neuron(-1, -1, inputData->numberOfInput, 0, 0, 0);
		   string nameOfOutputFileCon1(inputParameters->nameOfOutputFile);
		   string localFile = nameOfOutputFileCon1;

		   double **arrayResult = new double*[inputData->dimensionY* inputData->dimensionX];


		   double min = -1, max = -1;

		   string text = "";
		   int position = 0;
		   for (int i = 0; i < inputData->dimensionY; i++)
		   for (int k = 0; k < inputData->dimensionX; k++)
		   {
			   position = i * inputData->dimensionX + k;
			   arrayResult[position] = new double[8];
			   for (int j = 0; j < 8; j++)
			   {
				   arrayResult[position][j] = 10;
			   }

			   if (poleNeuronu[i][k] == NULL)
			   {
				   SendAndReciveWeight(i, k, tempNeuron->vahy);
				   poleNeuronu[i][k] = tempNeuron;
			   }

			   poleNeuronu[i][k]->gngCosin = inputParameters->gngCosin;
			   poleNeuronu[i][k]->gngMinkowskehoNumber = inputParameters->gngMinkow;


			   if (((k + 1) < inputData->dimensionX) && (i - 1) >= 0)
			   {

				   if (SendAndReciveWeight(i - 1, k + 1, temparray))
					   poleNeuronu[i][k]->vstupy = temparray;// poleNeuronu[i - 1][k + 1]->vahy;
				   else
					   poleNeuronu[i][k]->vstupy = poleNeuronu[i - 1][k + 1]->vahy;

				   arrayResult[position][0] = poleNeuronu[i][k]->EuklidMikowsky();

				   if (min == -1)
					   min = arrayResult[position][0];
				   if (max == -1)
					   max = arrayResult[position][0];
				   if (min >  arrayResult[position][0])
					   min = arrayResult[position][0];
				   if (max <  arrayResult[position][0])
					   max = arrayResult[position][0];

			   }


			   if (((k + 1) < inputData->dimensionX))
			   {

				   if (SendAndReciveWeight(i, k + 1, temparray))

					   poleNeuronu[i][k]->vstupy = temparray;//poleNeuronu[i][k + 1]->vahy;
				   else
					   poleNeuronu[i][k]->vstupy = poleNeuronu[i][k + 1]->vahy;

				   arrayResult[position][1] = poleNeuronu[i][k]->EuklidMikowsky();
				   if (min == -1)
					   min = arrayResult[position][1];
				   if (max == -1)
					   max = arrayResult[position][1];
				   if (min >  arrayResult[position][1])
					   min = arrayResult[position][1];
				   if (max <  arrayResult[position][1])
					   max = arrayResult[position][1];
			   }

			   if (((k + 1) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
			   {

				   if (SendAndReciveWeight(i + 1, k + 1, temparray))
					   poleNeuronu[i][k]->vstupy = temparray;// poleNeuronu[i + 1][k + 1]->vahy;
				   else
					   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k + 1]->vahy;

				   arrayResult[position][2] = poleNeuronu[i][k]->EuklidMikowsky();
				   if (min == -1)
					   min = arrayResult[position][2];
				   if (max == -1)
					   max = arrayResult[position][2];
				   if (min >  arrayResult[position][2])
					   min = arrayResult[position][2];
				   if (max <  arrayResult[position][2])
					   max = arrayResult[position][2];
			   }


			   if (((k) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
			   {

				   if (SendAndReciveWeight(i + 1, k, temparray))

					   poleNeuronu[i][k]->vstupy = temparray;// poleNeuronu[i + 1][k]->vahy;
				   else
					   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k]->vahy;

				   arrayResult[position][3] = poleNeuronu[i][k]->EuklidMikowsky();
				   if (min == -1)
					   min = arrayResult[position][3];
				   if (max == -1)
					   max = arrayResult[position][3];
				   if (min >  arrayResult[position][3])
					   min = arrayResult[position][3];
				   if (max <  arrayResult[position][3])
					   max = arrayResult[position][3];
			   }

			   if (poleNeuronu[i][k]->souradnice[0] == -1)
				   poleNeuronu[i][k] = NULL;
		   }
		   //tw.Close();

		   if (Size != 1)
		   {
			   int arrayForPosition[2];
			   arrayForPosition[0] = -1;
			   arrayForPosition[1] = -1;
			   for (int i = 1; i < Size; i++)
			   {
				   MPI_Send(arrayForPosition, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
			   }
		   }
		   //normalizace
		   int celkemNeurons = inputData->dimensionY* inputData->dimensionX;
		   for (int i = celkemNeurons - 1; i >= 0; i--)
		   {
			   for (int j = 0; j < 4; j++)
			   {
				   if (((int)(arrayResult[i][j])) != 10)
					   arrayResult[i][j] = ((arrayResult[i][j] - min) / (double)(max - min));
			   }
		   }


		   int tempXY = 0;
		   for (int p = inputData->dimensionY - 1; p >= 0; p--)
		   for (int k = inputData->dimensionX - 1; k >= 0; k--)
		   {
			   int i = p * inputData->dimensionX + k;

			   if ((k != 0) && (p != 0))
			   {
				   tempXY = i - 1 - inputData->dimensionX;
				   arrayResult[i][6] = arrayResult[tempXY][2];
			   }

			   if ((k != 0))
			   {
				   tempXY = i - 1;
				   arrayResult[i][5] = arrayResult[tempXY][1];
			   }

			   if ((p != 0))
			   {
				   tempXY = i - inputData->dimensionX;
				   arrayResult[i][7] = arrayResult[tempXY][3];
			   }
			   if ((k != 0) && (p != (inputData->dimensionY - 1)))
			   {
				   tempXY = i - 1 + inputData->dimensionX;
				   arrayResult[i][4] = arrayResult[tempXY][0];
			   }

		   }

		   int positionXY;
		   int *mapNeurons = new int[inputData->dimensionY* inputData->dimensionX];
		   int tempSizeMapNeuorns = inputData->dimensionY* inputData->dimensionX;

		   std::map<double, int*> sortedDistance;
		   int *posSorDis;
		   for (int p = inputData->dimensionY - 1; p >= 0; p--)
		   for (int k = inputData->dimensionX - 1; k >= 0; k--)
		   {
			   int i = p * inputData->dimensionX + k;
			   for (int j = 0; j < 4; j++)
			   {
				   if (arrayResult[i][j] == 10)
					   continue;
				   switch (j)
				   {
				   case 0:positionXY = i - inputData->dimensionX + 1;
					   break;
				   case 1:positionXY = i + 1;
					   break;
				   case 2:positionXY = i + inputData->dimensionX + 1;
					   break;
				   case 3:positionXY = i + inputData->dimensionX;
					   break;
				   default:
					   break;
				   }
				
				   if (sortedDistance.count(arrayResult[i][j]) == 0)
				   {
					   posSorDis = new int[2];
					   posSorDis[0] = i;
					   posSorDis[1] = positionXY;
					   sortedDistance.insert(std::pair<double, int*>(arrayResult[i][j],posSorDis));
				   }

			   }
			   
		   }

		   for (int u = 0; u < tempSizeMapNeuorns; u++)
		   {
			   mapNeurons[u] = u;
		   }

		   double errTemp = treshHold;
		   int numberOfIncreaseTreshold = 0;;
		   //for (int i = celkemNeurons-1; i  >=0; i --)
		   double trassAdd = 0.005;
		   std::map<int, int> tempCounter;
	    int counter = tempSizeMapNeuorns;

		int tempMapa = 0;
		bool completeGo;
		   for (std::map<double, int *>::const_iterator it = sortedDistance.begin(); it != sortedDistance.end(); ++it)
		   {
			   completeGo = true;
			   if (inputParameters->debug)
				cout << it->first << "\t" << it->second[0] << "\t" << it->second[1]<<endl;
			   if (mapNeurons[it->second[0]] == mapNeurons[it->second[1]])
				   continue;

			   tempMapa = mapNeurons[it->second[0]];
			   for (int u = 0; u < tempSizeMapNeuorns; u++)
			   {
				   if (mapNeurons[u] == tempMapa)
				   {
					   mapNeurons[u] = mapNeurons[it->second[1]];
				   }
			   }
			   mapNeurons[it->second[0]] = mapNeurons[it->second[1]];

			   tempCounter.clear();
			   for (int u = 0; u < tempSizeMapNeuorns; u++)
			   {
				   if (tempCounter.count(mapNeurons[u]) == 0)
				   {
					   tempCounter.insert(std::pair<int, int>(mapNeurons[u], 1));
				   }
			   }
			   counter = tempCounter.size();
			   if (inputParameters->debug)
				cout << "\tBumber of separate parts:" << counter << "" << inputParameters->gngMaximumParts << endl;
			   if ((counter <= inputParameters->gngMaximumParts) && (inputParameters->gngMaximumParts != -1))
			   {
				   completeGo = false;
				   
				   break;
			   }
			   

		   }
		   
		   if (completeGo)
		   {
			   cout << "Full walked"<<"\tBumber of separate parts:" << counter << endl;
		   }


		   for (int u = 0; u < tempSizeMapNeuorns; u++)
		   {
			   cout << mapNeurons[u] << endl;
		   }

		   std::map<int, int> numberintogroup;

		   for (int i = 0; i < inputData->numberOfRecords; i++)
		   {
			   resultTrainingSet[i] = mapNeurons[resultTrainingSet[i]];
			   //	   cout << i << "\t" << resultTrainingSet[i]<<endl;
			   if (numberintogroup.count(resultTrainingSet[i]) == 0)
			   {

				   numberintogroup.insert(std::pair<int, int>(resultTrainingSet[i], 1));
			   }
			   else
			   {
				   numberintogroup[resultTrainingSet[i]]++;
			   }
		   }
		   int checkSum = 0;
		   int internalCounter = 0;
		   for (int i = 0; i < celkemNeurons; i++)
		   {
			   if (numberintogroup[i] != 0)
			   {
				   cout << "Group:" << internalCounter << " has:" << numberintogroup[i] << " records" << endl;
				   checkSum += numberintogroup[i];
				   internalCounter++;
			   }
		   }
		   cout << "Complete sum is:" << checkSum << endl;

		   delete mapNeurons;
		   delete arrayResult;
		   delete temparray;
		   delete tempNeuron;
	   }
	   void SOMnetwork::MargeTrainingData(double treshHold,int * resultTrainingSet)
	   {
		   if (Size > 1)
		   {
			   if (Rank != 0)
			   {
				   int temprecv[2]; // prvni y druhe x

				   MPI_Status mpistat;
				   while(true)
				   {
					   MPI_Recv(&temprecv,2,MPI_INT,0,0,MPI_COMM_WORLD,&mpistat);
					   if(temprecv[0]==-1)
						   break;

					   if (poleNeuronu[temprecv[0]][temprecv[1]] != NULL)
					   {
						   MPI_Send(poleNeuronu[temprecv[0]][temprecv[1]]->vahy,poleNeuronu[temprecv[0]][temprecv[1]]->numberOfdimension,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
					   }
					   else
					   {
						   cout<<"Error in the function MargeTrainingData"<<endl;
					   }

				   }

				   return;
			   }
		   }

			double *temparray=new double[inputData->numberOfInput];
			Neuron * tempNeuron= new Neuron(-1,-1,inputData->numberOfInput,0,0,0);
		   	string nameOfOutputFileCon1(inputParameters->nameOfOutputFile);
			string localFile=nameOfOutputFileCon1;

			double **arrayResult = new double*[inputData->dimensionY* inputData->dimensionX];
		 

           double min = -1, max = -1;
 
           string text = "";
		   int position=0;
           for (int i = 0; i < inputData->dimensionY; i++)
               for (int k = 0; k < inputData->dimensionX; k++)
               {
				   position=i * inputData->dimensionX + k;
				   arrayResult[position]=new double[8];
				   for (int j = 0; j < 8; j++)
				   {
					   arrayResult[position][j]=10;
				   }

				if( poleNeuronu[i][k]==NULL)
				   {
					   SendAndReciveWeight(i,k,tempNeuron->vahy);
					   poleNeuronu[i][k]=tempNeuron;
				   }

				poleNeuronu[i][k]->gngCosin = inputParameters->gngCosin;
				poleNeuronu[i][k]->gngMinkowskehoNumber = inputParameters->gngMinkow;


                   if (((k + 1) < inputData->dimensionX) && (i - 1) >= 0)
                   {

                           if(SendAndReciveWeight(i-1,k + 1,temparray))
								poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i - 1][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i - 1][k + 1]->vahy;
						   
						   arrayResult[position][0] = poleNeuronu[i][k]->EuklidMikowsky();

						      if (min == -1)
                               min =  arrayResult[position][0] ;
                           if (max == -1)
                               max =  arrayResult[position][0] ;
                           if (min >  arrayResult[position][0] )
                               min =  arrayResult[position][0] ;
                           if (max <  arrayResult[position][0] )
                               max =  arrayResult[position][0] ;
                    
                   }

                   
                   if (((k + 1) < inputData->dimensionX))
                   {
                       
                           if(SendAndReciveWeight(i,k + 1,temparray))

								 poleNeuronu[i][k]->vstupy =temparray;//poleNeuronu[i][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i][k + 1]->vahy;

						   arrayResult[position][1] = poleNeuronu[i][k]->EuklidMikowsky();
                            if (min == -1)
                               min =  arrayResult[position][1] ;
                           if (max == -1)
                               max =  arrayResult[position][1] ;
                           if (min >  arrayResult[position][1] )
                               min =  arrayResult[position][1] ;
                           if (max <  arrayResult[position][1] )
                               max =  arrayResult[position][1] ;
                   }

                   if (((k + 1) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       
                           if(SendAndReciveWeight(i+1,k + 1,temparray))
							   poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i + 1][k + 1]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k + 1]->vahy;

						   arrayResult[position][2] = poleNeuronu[i][k]->EuklidMikowsky();
                              if (min == -1)
                               min =  arrayResult[position][2] ;
                           if (max == -1)
                               max =  arrayResult[position][2] ;
                           if (min >  arrayResult[position][2] )
                               min =  arrayResult[position][2] ;
                           if (max <  arrayResult[position][2] )
                               max =  arrayResult[position][2] ;
                   }

                   
                   if (((k) < inputData->dimensionX) && (i + 1) < inputData->dimensionY)
                   {
                       
                           if(SendAndReciveWeight(i+1,k,temparray))

							 poleNeuronu[i][k]->vstupy =temparray;// poleNeuronu[i + 1][k]->vahy;
						   else
							   poleNeuronu[i][k]->vstupy = poleNeuronu[i + 1][k]->vahy;

						   arrayResult[position][3] = poleNeuronu[i][k]->EuklidMikowsky();
                              if (min == -1)
                               min =  arrayResult[position][3] ;
                           if (max == -1)
                               max =  arrayResult[position][3] ;
                           if (min >  arrayResult[position][3] )
                               min =  arrayResult[position][3] ;
                           if (max <  arrayResult[position][3] )
                               max =  arrayResult[position][3] ;
                   }

				   	if( poleNeuronu[i][k]->souradnice[0]==-1)
					     poleNeuronu[i][k]=NULL;
               }
           //tw.Close();

			   if(Size!=1)
			   {
				   int arrayForPosition[2];
				   arrayForPosition[0]=-1;
				   arrayForPosition[1]=-1;
				   for (int i = 1; i < Size; i++)
				   {
					   MPI_Send(arrayForPosition,2,MPI_INT,i,0, MPI_COMM_WORLD);
				   }
			   }
			   //normalizace
			   int celkemNeurons= inputData->dimensionY* inputData->dimensionX;
			   for (int i = celkemNeurons-1; i  >=0; i --)
			   {
				   for (int j = 0; j < 4; j++)
					   {
						   if(((int)(arrayResult[i][j]))!=10)
						   arrayResult[i][j]= (( arrayResult[i][j] - min) / (double)(max - min));
				   }
			   }

			  
			   //double **TempArrayResult = new double*[inputData->dimensionY];
			   //for (int i = 0; i < inputData->dimensionY; i++)
			   //{
				  // TempArrayResult[i] = new double[inputData->dimensionX];
			   //}
			   //int posx = 0;
			   //int posy = 0;
			   //for (int p = inputData->dimensionY - 1; p >= 0; p--)
			   //for (int k = inputData->dimensionX - 1; k >= 0; k--)			   
			   //{
				  // int i = p * inputData->dimensionX + k;
				  // for (int j = 0; j < 4; j++)
				  // {
					 //  if (arrayResult[i][j] == 10)
						//   continue;
					 //  switch (j)
					 //  {
					 //  case 0://positionXY = i - inputData->dimensionX + 1;
						//   posx = k + 1;
						//   posy = p - 1;
						//   break;
					 //  case 1://positionXY = i + 1;
						//   posx = k + 1;
						//   posy = p ;
						//   break;
					 //  case 2://positionXY = i + inputData->dimensionX + 1;
						//   posx = k + 1;
						//   posy = p + 1;
						//   break;
					 //  case 3://positionXY = i + inputData->dimensionX;
						//   posx = k ;
						//   posy = p + 1;
						//   break;
					 //  default:
						//   break;
					 //  }
					 //  TempArrayResult[posy][posx] = arrayResult[i][j];
					 //  TempArrayResult[posx][posy] = arrayResult[i][j];
				  // }
			   //}

			   //for (int p = inputData->dimensionY - 1; p >= 0; p--)
			   //{

				  // for (int k = inputData->dimensionX - 1; k >= 0; k--)
				  // {
					 //  cout << TempArrayResult[p][k] << "\t";
				  // }
				  // cout << endl;
			   //}
			   int tempXY = 0;
			   for (int p = inputData->dimensionY - 1; p >= 0; p--)
			   for (int k = inputData->dimensionX - 1; k >= 0; k--)			   
			   {
				   int i = p * inputData->dimensionX + k;
				   
				   if ((k != 0) && (p != 0))
				   {
					   tempXY = i - 1 - inputData->dimensionX;
					   arrayResult[i][6] = arrayResult[tempXY][2];
				   }

				   if ((k != 0))
				   {
					   tempXY = i - 1;
					   arrayResult[i][5] = arrayResult[tempXY][1];
				   }

				   if ((p != 0))
				   {
					   tempXY = i - inputData->dimensionX;
					   arrayResult[i][7] = arrayResult[tempXY][3];
				   }
				   if ((k != 0) && (p != (inputData->dimensionY - 1)))
				   {
					   tempXY = i - 1 + inputData->dimensionX;
					   arrayResult[i][4] = arrayResult[tempXY][0];
				   }
				   
			   }

			   int positionXY;
			   int *mapNeurons= new int[inputData->dimensionY* inputData->dimensionX];
			   int tempSizeMapNeuorns = inputData->dimensionY* inputData->dimensionX;

			   for (int u = 0; u < tempSizeMapNeuorns; u++)
			   {
			    mapNeurons[u] = u;
			   }

			 //  treshHold = 0.7;
			   double errTemp=treshHold;
			   int numberOfIncreaseTreshold=0;;
			   //for (int i = celkemNeurons-1; i  >=0; i --)
			   double trassAdd = 0.005;
			   int posXYtemp = 0;
			   std::map<int, int> tempCounter;
			   int tempMapa = 0;
		   incTreshold:			   int counter = tempSizeMapNeuorns;

			   //for (int u = 0; u < tempSizeMapNeuorns; u++)
			   //{
				  // mapNeurons[u] = u;
			   //}

			   for(int k=inputData->dimensionX-1;k>=0;k--)
				   for(int p=inputData->dimensionY-1;p>=0;p--)
				   {
					   bool find=false;
					   int i= p * inputData->dimensionX + k;
					   errTemp=treshHold;
					   for (int j = 0; j < 8; j++)
					   {
						   if(errTemp>=arrayResult[i][j])
						   {
							   switch (j)
							   {
							   case 0:positionXY=i-inputData->dimensionX+1;
								   break;
							   case 1:positionXY=i+1;
								   break;
							   case 2:positionXY=i+inputData->dimensionX+1;
								   break;
							   case 3:positionXY=i+inputData->dimensionX;
								   break;
							   case 4:positionXY = i + inputData->dimensionX - 1;
								   break;
							   case 5:positionXY = i - 1;
								   break;
							   case 6:positionXY = i - inputData->dimensionX - 1;
								   break;
							   case 7:positionXY = i - inputData->dimensionX;
								   break;
							   default:
								   break;
							   }
							   if (mapNeurons[i] == mapNeurons[positionXY])
								   continue;
							   errTemp = arrayResult[i][j];
							   posXYtemp = positionXY;
							   find = true;
						   }
					   }

					   if(find)
					   {
						   if (mapNeurons[i] == mapNeurons[posXYtemp])
							   continue;
						   tempMapa = mapNeurons[i];
						   for (int u = 0; u < tempSizeMapNeuorns; u++)
						   {
							   if (mapNeurons[u] == tempMapa)
							   {
								   mapNeurons[u] = mapNeurons[posXYtemp];
							   }
						   }
						   mapNeurons[i] = mapNeurons[posXYtemp];
						   counter--;
					   }
					   else
					   {
						 //  mapNeurons[i]=i;
						 //  counter++;
					   }
				   }
				   tempCounter.clear();
				   for (int u = 0; u < tempSizeMapNeuorns; u++)
				   {
					   if (tempCounter.count(mapNeurons[u]) == 0)
					   {
						   tempCounter.insert(std::pair<int, int>(mapNeurons[u], 1));
					   }
				   }
				   counter = tempCounter.size();
				   
				   if((counter>inputParameters->gngMaximumParts)&& (inputParameters->gngMaximumParts!=-1) )
				   {
					   if (trassAdd > 0.00001)
						   trassAdd = trassAdd*0.95;
					  treshHold += trassAdd;					   
					   numberOfIncreaseTreshold++;
					   cout << "Trashold:" << treshHold << "\tTrasAdd:" << trassAdd << "\tBumber of separate parts:" << counter << endl;
					   goto incTreshold;
				   }
				   cout<<"Number of separate parts:"<<counter<<endl;
				   cout<<"Set Treshold:"<<inputParameters->gngTreshHold<<"\tReal Treshold:"<<treshHold<<endl;
				   //for (int i = 0; i < celkemNeurons; i++)
				   //{
					  // cout<<i<<"\t"<<mapNeurons[i]<<endl;
				   //}

				   for (int u = 0; u < tempSizeMapNeuorns; u++)
				   {
					   cout << mapNeurons[u] << endl;
				   }

				   std::map<int,int> numberintogroup;

				   for (int i = 0; i < inputData->numberOfRecords; i++)
				   {
					   resultTrainingSet[i]=mapNeurons[resultTrainingSet[i]];
				//	   cout << i << "\t" << resultTrainingSet[i]<<endl;
					   	if(numberintogroup.count(resultTrainingSet[i])==0)
						{

							numberintogroup.insert(std::pair<int,int>(resultTrainingSet[i], 1));
						}else
						{
						numberintogroup[resultTrainingSet[i]]++;
						}
				   }
				   int checkSum=0;
				   int internalCounter=0;
				   for (int i = 0; i < celkemNeurons; i++)
				   {
					   if(numberintogroup[i]!=0)
					   {
						   cout<<"Group:"<<internalCounter<<" has:"<<numberintogroup[i] <<" records"<<endl;
						   checkSum+=numberintogroup[i];
						   internalCounter++;
					   }					   
				   }
				   cout<<"Complete sum is:" <<checkSum<<endl;

			   delete mapNeurons;
			   delete arrayResult;
			   delete temparray;
			   delete tempNeuron;
	   }

	   bool SOMnetwork::Load()
	   {
		   string nameOfOutputFileCon(inputParameters->nameOfOutputFile);
		   string temp = nameOfOutputFileCon;// inputData->nameOfFile;
		   temp.append(NumberToString<int>(Size));
		   if (Rank == 0)
		   {
			   char nullChar = nameOfOutputFileCon.size();
			   string localFile = temp;
			   localFile.append(".Bbin");
			   ifstream myFile(localFile, ios::out | ios::binary);
			   int tempLoad;
			   myFile.read((char*)&tempLoad, sizeof tempLoad);
			   if (tempLoad != inputData->dimensionX)
				   return false;

			   myFile.read((char*)&tempLoad, sizeof tempLoad);
			   if (tempLoad != inputData->dimensionY)
				   return false;
			   myFile.read((char*)&tempLoad, sizeof tempLoad);
			   if (tempLoad != inputData->numberOfInput)
				   return false;

			   myFile.read((char*)&tempLoad, sizeof tempLoad);
			   if (tempLoad != inputData->numberOfRepeating)
				   return false;

			   myFile.read((char*)&tempLoad, sizeof tempLoad);
			   if (tempLoad != Size)
				   return false;

			   char tempNull;

			   myFile.read((char*)&tempNull, sizeof tempNull);
			   if (tempNull != nullChar)
				   return false;

			   char * tempNameofOutputFile = new char[((int)nullChar)+1];
			   myFile.read(tempNameofOutputFile,  ((int)nullChar));
			   tempNameofOutputFile[((int)nullChar)] = '\0';
			   if (nameOfOutputFileCon.compare(tempNameofOutputFile)!=0)
				   return false;

			   nullChar = strlen(inputParameters->nameOfInputDataFile);
			   myFile.read((char*)&tempNull, sizeof tempNull);
			   if (tempNull != nullChar)
				   return false;

			   delete tempNameofOutputFile;



			   tempNameofOutputFile = new char[((int)nullChar)+1];
			   myFile.read(tempNameofOutputFile, (int)nullChar);
			   tempNameofOutputFile[((int)nullChar)] = '\0';
			   if (strcmp(inputParameters->nameOfInputDataFile,tempNameofOutputFile) != 0)
				   return false;
			   delete tempNameofOutputFile;

			   if (myFile.peek() == EOF)
			   {

				   std::cout << "all characters read successfully." << endl;
				   inputParameters->version = 1;
			   }
			   else
			   {
				   double tempLoad2;
				 //  cout << "error: only " << myFile.gcount() << " could be read" << endl;
				   cout << "New version" << endl;
				   myFile.read((char*)&tempLoad2, sizeof tempLoad2);
				   inputParameters->version = tempLoad2;
			   }
				   
			   //myFile.read(reinterpret_cast<const char*>(&inputData->dimensionY), sizeof inputData->dimensionY);
			   //myFile.read(reinterpret_cast<const char*>(&inputData->numberOfInput), sizeof inputData->numberOfInput);
			   //myFile.read(reinterpret_cast<const char*>(&inputData->numberOfRepeating), sizeof inputData->numberOfRepeating);
			   //myFile.read(reinterpret_cast<const char*>(&Size), sizeof Size);
			   //myFile.read(reinterpret_cast<const char*>(&nullChar), sizeof nullChar);
			   //	cout<<nullChar<<endl;
			//   myFile.read(nameOfOutputFileCon.c_str(), nameOfOutputFileCon.size());
			   //nullChar = sizeof inputParameters->nameOfInputDataFile;
			   //myFile.read(reinterpret_cast<const char*>(&nullChar), sizeof nullChar);
			   ////	cout<<nullChar<<endl;
			   //myFile.read(inputParameters->nameOfInputDataFile, (sizeof inputParameters->nameOfInputDataFile));

			   myFile.close();
			   cout << "Version:" << inputParameters->version << endl;
		   }
		   MPI_Bcast(&inputParameters->version, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		   if (0 == 0)
		   {
			   string localFile = temp;
			   localFile.append(NumberToString<int>(Rank));
			   localFile.append(".1Bbin");
			   ifstream myFile(localFile, ios::out | ios::binary);

			   int tempLoad;
			   myFile.read((char*)&tempLoad, sizeof tempLoad);
			   if (tempLoad != pocetNeuronu1)
				   return false;

			//   myFile.write(reinterpret_cast<const char*>(&pocetNeuronu1), sizeof pocetNeuronu1);

			   //BinaryWriter bw = new BinaryWriter(File.Create(NazevSouboru + comm.Size + comm.Rank + ".1Bbin"));
			   //bw.Write(pocetNeuronu1);
			   double tempLoad2;
			   for (int i = 0; i < inputData->dimensionY; i++)
			   for (int k = 0; k < inputData->dimensionX; k++)
			   {
				   if (inputParameters->version >= 2)
				   {

					   if (poleNeuronu[i][k] == NULL)
						   continue;
					   myFile.read((char*)&tempLoad, sizeof tempLoad);
					   poleNeuronu[i][k]->souradnice[0] = tempLoad;

					   myFile.read((char*)&tempLoad, sizeof tempLoad);
					   poleNeuronu[i][k]->souradnice[1] = tempLoad;
				   }
				   for (int j = 0; j < inputData->numberOfInput; j++)
				   {
					   if (inputParameters->version < 2)
					   {

						   if (poleNeuronu[i][k] == NULL)
							   continue;
						   myFile.read((char*)&tempLoad, sizeof tempLoad);
						   poleNeuronu[i][k]->souradnice[0] = tempLoad;

						   myFile.read((char*)&tempLoad, sizeof tempLoad);
						   poleNeuronu[i][k]->souradnice[1] = tempLoad;
					   }

					   myFile.read((char*)&tempLoad2, sizeof tempLoad2);
					   poleNeuronu[i][k]->vahy[j] = tempLoad2;
					  

					//   myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->souradnice[0]), sizeof poleNeuronu[i][k]->souradnice[0]);
					//   myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->souradnice[1]), sizeof poleNeuronu[i][k]->souradnice[1]);
					//   myFile.write(reinterpret_cast<const char*>(&poleNeuronu[i][k]->vahy[j]), sizeof poleNeuronu[i][k]->vahy[j]);

				   }
				   poleNeuronu[i][k]->RecalculateWeight();
			   }
			   myFile.close();
		   }
		   return true;
	   }