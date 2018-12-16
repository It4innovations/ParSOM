#include "SaveAndLoad.h"
#include <iostream>
#include<sstream>
#include <fstream>
#include <string>
#include <vector>


using namespace std;
SaveAndLoad::SaveAndLoad(ParameterParser *Parameters)
{
	inputParameters=Parameters;
}


SaveAndLoad::~SaveAndLoad(void)
{
}

int ConvertStringToInt(string text)
{
	int Result;
	stringstream convert(text); // stringstream used for the conversion initialized with the contents of Text

	if ( !(convert >> Result) )//give the value to Result using the characters in the string
	 Result = -1;
	return Result;
}

template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

template <typename T>
T StringToNumber ( const string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

void split(string &s1, string &s2, char delim){
    size_t i;

    i=s1.find_first_of(delim);
    s2.append(s1, i+1, s1.size());
    s1.erase(i, s1.size());
  }



InputData SaveAndLoad::Load(void)
{
	double wmtimestart = MPI_Wtime();
	InputData dataCollection;
	dataCollection.problemWithFile=0;
	ifstream myConfigFile(inputParameters->nameOfConfigFile);
	string line;
	if (myConfigFile.is_open())
	{		
		string sub;
			getline (myConfigFile,line);				
			istringstream iss(line);

    	    iss >> sub;        
			dataCollection.dimensionX= ConvertStringToInt(sub);
			iss >> sub;        
			dataCollection.dimensionY= ConvertStringToInt(sub);
			iss >> sub;        
			dataCollection.numberOfRepeating= ConvertStringToInt(sub);
			iss >> sub;        
			dataCollection.numberOfInput= ConvertStringToInt(sub);
			iss >> sub;        
			dataCollection.nameOfFile= sub;
			myConfigFile.close();
  }else
	{
		dataCollection.problemWithFile=1;
  }

	vector< vector<double> > tempData;
	vector< vector<int> > tempPosition;
	  int maximalniCall = -1;
	ifstream myDataFile(inputParameters->nameOfInputDataFile);
	if(myDataFile.is_open())
	{
		    while ( myDataFile.good() )
			{
				vector<double> tempOneLineData;
				vector<int> tempOneLinePosition;
			  getline (myDataFile,line);

			  istringstream iss(line);
					do
					{
						string sub;
						string secondPart;
						iss >> sub;
						if(sub=="")
							break;
						split(sub,secondPart,':');
						tempOneLinePosition.push_back(StringToNumber<int>(sub));
						tempOneLineData.push_back(StringToNumber<double>(secondPart));

						 int ook = StringToNumber<int>(sub);

						if (ook > maximalniCall)
                           maximalniCall=ook;

					//	cout << "Substring: " << sub << endl;
					} while (iss);
					if(tempOneLineData.size()!=0)
					{
					tempData.push_back(tempOneLineData);
					tempPosition.push_back(tempOneLinePosition);
					}
			}
				myDataFile.close();

		// convert dynamic data to usefull:)
 				dataCollection.numberOfRecords=tempData.size();
				dataCollection.data=new double*[dataCollection.numberOfRecords];
				dataCollection.noZero=new int*[dataCollection.numberOfRecords];
				dataCollection.numberOfDataInRows = new int[dataCollection.numberOfRecords];
				dataCollection.SumOfInput=new double[dataCollection.numberOfRecords];
				for (int i = 0; i < dataCollection.numberOfRecords; i++)
				{
					int tempNumber=tempData[i].size();
					dataCollection.data[i]=new double[tempNumber];
					dataCollection.noZero[i]=new int[tempNumber];
					dataCollection.numberOfDataInRows[i]=tempNumber;
					double sum=0;
					for (int j = 0; j < tempNumber; j++)
					{
						dataCollection.data[i][j]=tempData[i][j];
						dataCollection.noZero[i][j]=tempPosition[i][j];

						sum+=tempData[i][j]*tempData[i][j];
					}
					dataCollection.SumOfInput[i]=sum;
				}

				if(dataCollection.numberOfInput !=(maximalniCall+1))
				{
					cout<<"Divny pocet atributu, skutecny je:"<<(maximalniCall+1)<<endl;				
				}

               dataCollection.numberOfInput=maximalniCall+1;

	}
	else
	{
		dataCollection.problemWithFile=2;
		return dataCollection;
	}
	double wmtimekonec = MPI_Wtime();
	cout<<"Time:"<<wmtimekonec - wmtimestart<<endl;

	return dataCollection;
}

void SaveAndLoad::Load(InputData *dataCollection)
{
	double wmtimestart = MPI_Wtime();
	dataCollection->problemWithFile = 0;
	ifstream myConfigFile(inputParameters->nameOfConfigFile);
	string line;
	if (myConfigFile.is_open())
	{
		string sub;
		getline(myConfigFile, line);
		istringstream iss(line);

		iss >> sub;
		dataCollection->dimensionX = ConvertStringToInt(sub);
		iss >> sub;
		dataCollection->dimensionY = ConvertStringToInt(sub);
		iss >> sub;
		dataCollection->numberOfRepeating = ConvertStringToInt(sub);
		iss >> sub;
		dataCollection->numberOfInput = ConvertStringToInt(sub);
		iss >> sub;
		dataCollection->nameOfFile = sub;
		myConfigFile.close();
	}
	else
	{
		dataCollection->problemWithFile = 1;
	}

	vector< vector<double> > tempData;
	vector< vector<int> > tempPosition;
	int maximalniCall = -1;
	ifstream myDataFile(inputParameters->nameOfInputDataFile);
	if (myDataFile.is_open())
	{
		while (myDataFile.good())
		{
			vector<double> tempOneLineData;
			vector<int> tempOneLinePosition;
			getline(myDataFile, line);

			istringstream iss(line);
			do
			{
				string sub;
				string secondPart;
				iss >> sub;
				if (sub == "")
					break;
				split(sub, secondPart, ':');
				tempOneLinePosition.push_back(StringToNumber<int>(sub));
				tempOneLineData.push_back(StringToNumber<double>(secondPart));

				int ook = StringToNumber<int>(sub);

				if (ook > maximalniCall)
					maximalniCall = ook;

				//	cout << "Substring: " << sub << endl;
			} while (iss);
			if (tempOneLineData.size() != 0)
			{
				tempData.push_back(tempOneLineData);
				tempPosition.push_back(tempOneLinePosition);
			}
		}
		myDataFile.close();

		// convert dynamic data to usefull:)
		dataCollection->numberOfRecords = tempData.size();
		dataCollection->data = new double*[dataCollection->numberOfRecords];
		dataCollection->noZero = new int*[dataCollection->numberOfRecords];
		dataCollection->numberOfDataInRows = new int[dataCollection->numberOfRecords];
		dataCollection->SumOfInput = new double[dataCollection->numberOfRecords];
		for (int i = 0; i < dataCollection->numberOfRecords; i++)
		{
			int tempNumber = tempData[i].size();
			dataCollection->data[i] = new double[tempNumber];
			dataCollection->noZero[i] = new int[tempNumber];
			dataCollection->numberOfDataInRows[i] = tempNumber;
			double sum = 0;
			for (int j = 0; j < tempNumber; j++)
			{
				dataCollection->data[i][j] = tempData[i][j];
				dataCollection->noZero[i][j] = tempPosition[i][j];

				sum += tempData[i][j] * tempData[i][j];
			}
			dataCollection->SumOfInput[i] = sum;
		}

	//	if (dataCollection->numberOfInput != (maximalniCall + 1))
//		{
			cout << "Divny pocet atributu, skutecny je:" << (maximalniCall + 1) << endl;
	//	}

		dataCollection->numberOfInput = maximalniCall + 1;

	}
	else
	{
		dataCollection->problemWithFile = 2;
		return ;
	}
	double wmtimekonec = MPI_Wtime();
	cout << "Time:" << wmtimekonec - wmtimestart << endl;

	return ;
}