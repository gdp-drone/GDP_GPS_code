"""
Author: Lukasz Sliwinski luki3141@gmail.com
date: June 2019
project: ADD GDP
This code transforms the GPS data output from the rosbag files into data in ENU frame of reference using GeographicLib
"""

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> 
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>


using namespace std;
using namespace GeographicLib;

int main(int argc, char *argv[])
{
	//Initializing GeographicLib
	try 
	{
		Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
		fstream infile;
		fstream outfile;
		cout << argv[1]<<endl;
		//opening the input and output files
		string inname="data/gpsData";
		inname.append(argv[1]);
		string outname="data/reduced";
		outname.append(argv[1]);
		string line;
		infile.open(inname, fstream::in);
		outfile.open(outname, fstream::out | fstream::trunc);
		int i = 0;
		int numData = 0;
		size_t found,found2;
		long double lati,longi,alti;
		long int secs,nsecs;
		// Looking through the input file for latitude, longitude, altitude and time data to put into output file
		while (!(infile.eof()))
		{
			infile >> line;
			found = line.find("latitude");
			if (found!=string::npos)
			{
				infile >>lati;
				i++;
			}
			found = line.find("longitude");
			if (found!=string::npos)
			{
				infile>>longi;
				i++;
			}
			found = line.find("altitude");
			if (found!=string::npos)
			{
				infile>>alti;
				i++;
			}
			found = line.find("secs:");
			if (found!=string::npos)
			{
				found2 = line.find("nsecs:");
				if (found2!=string::npos)
				{
					infile>>nsecs;
					i++;
				}
				else
				{
					infile>>secs;
					i++;
				}
			}

			found = line.find("---");
			//If managed to collect every part of given data point the numbers are put into output file
			if (found!=string::npos)
			{
				if (i==5)
				{
					outfile<<setprecision(10)<<lati<<" "<<setprecision(10)<<longi<<" "<<setprecision(10)<<alti<<" "<<secs<<" "<<nsecs<<endl;
				}
				i=0; 
				numData++;
			}


		}
		outfile.close();
		infile.close();
		//opening the file to write transformed data in
		infile.open(outname, fstream::in);
		outname = "data/localized";
		outname.append(argv[1]);
		outfile.open(outname, fstream::out | fstream::trunc);
		double* longAr = new double[numData];
		double* latiAr = new double[numData];
		double* altiAr = new double[numData];
		long int* secsAr = new long int[numData];
		long int* nsecsAr = new long int[numData];
		double latiMean, longMean;
		latiMean = 0;
		longMean = 0;
		long int startT;
		//reading the data
		for (int j=0; j<(numData-1);j++)
		{
			infile>>latiAr[j];
			infile>>longAr[j];
			infile>>altiAr[j];
			infile>>secsAr[j];
			infile>>nsecsAr[j];
			if (j==0)
				startT=secsAr[0];
			secsAr[j]=secsAr[j]-startT;
			latiMean+=latiAr[j];
			longMean+=longAr[j];
		}

		latiMean = latiMean/(numData-1);
		longMean = longMean/(numData-1);

		cout << numData << endl;
		//Every data point is transformed to ENU and saved into localized<name.txt> file
		const double lat0 = 51.48244505, lon0 = -0.2180541661; //some point in testing location
    	LocalCartesian proj(lat0, lon0, 0, earth);
    	double x, y, z;
    	double h=0;
    	double xMean,yMean;
    	xMean=0;
    	yMean=0;
    	for (int j=0; j<numData-1;j++)
		{
			proj.Forward(latiAr[j], longAr[j], h, x, y, z);
			outfile<<setprecision(10)<<x<<" "<<setprecision(10)<<y<<" "<<setprecision(10)<<z<<" "<<secsAr[j]<<" "<<nsecsAr[j]<<endl;
			xMean+=x;
			yMean+=y;
			cout<<secsAr[j]<<endl;
		}
		cout <<setprecision(10)<<latiMean<<" "<<setprecision(10)<<longMean<<endl;
		cout <<setprecision(10)<<xMean<<" "<<setprecision(10)<<yMean<<endl;
		infile.close();
		outfile.close();
	}
  catch (const exception& e) 
  	{
   		cerr << "Caught exception: " << e.what() << "\n";
  	  	return 1;
 	}
	
	return 0;
}
