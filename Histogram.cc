/*
 *  Histogram.cc
 *  
 *
 *  Created by Christopher on 20/08/2010.
 *  Copyright 2010 University of Cambridge. All rights reserved.
 *
 */


using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>

double z(double t1, double t2, double a) {  //Evaluates length of line between points
	
	double z = a*sqrt(1/(cos(t1)*cos(t1) + a*sin(t1)*a*sin(t1)) + 1/(cos(t2)*cos(t2) +a*sin(t2)*a*sin(t2)) - 2.0*cos(t2-t1)/sqrt((cos(t1)*cos(t1) + a*sin(t1)*a*sin(t1))*(cos(t2)*cos(t2) + a*sin(t2)*a*sin(t2))));
	return z;
}


int main() {
	
	int n, m;		//n is number of theta values to test, m is no. of bins in histogram
	double a;		//a is the length of the semi major axis divided by that of the semi minor axis
	
	cout << "Enter aspect ration a" << endl;
	cin >> a;
	cout << "Enter number of theta values" << endl;
	cin >> n; 
	cout << "Enter the number of histogram bins required" << endl;
	cin >> m;
		
	double dtheta = 2.0*3.1415926/(double)n;	//Increment in theta
	
	double dz = 2.0*a/(double)m;		//This is used to calculate which histogram bin a z value should go into. It is essentially the width of each bin in terms of z
	
	int output[m];		//index of a bin in the array correponds to the lower bound of the histogram bin divided by dz
	
	for(int x=0; x<=m; x++) {
		output[x]=0;
	}
	
	int x = 0;						//Counts number of calculations performed
	
	for(int i=0; i<n; i++) {         //i and j correspond to theta1 and theta2
		for(int j=0; j<n; j++) {
			int k=0;
			x++;
			for(double zs = z(((double)i*dtheta), ((double)j*dtheta) ,a); zs>dz; zs = zs - dz) {	//Bin allocating procedure. Counts backwards, incrementing index
				k++;
			}
			output[k] = output[k] +1;
		}
	}
	cout << x << endl;
	ofstream fout;
	fout.open("HistogramOutput.txt");		//Prints data to pre existing file HistogramOutput.txt
	for(int i=0; i<m; i++) {
		fout << output[i] << endl;
	}
}
	
	
		
		

