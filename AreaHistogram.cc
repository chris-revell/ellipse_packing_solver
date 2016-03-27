/*
 *  AreaHistogram.cc
 *  
 *
 *  Created by Christopher on 04/09/2010.
 *  Copyright 2010 University of Cambridge. All rights reserved.
 *
 */


using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>

double z(double t1, double t2, double a) {		//Evaluates length of line between points
	
	double z = a*sqrt(1/(cos(t1)*cos(t1) + a*sin(t1)*a*sin(t1)) + 1/(cos(t2)*cos(t2) +a*sin(t2)*a*sin(t2)) - 2.0*cos(t2-t1)/sqrt((cos(t1)*cos(t1) + a*sin(t1)*a*sin(t1))*(cos(t2)*cos(t2) + a*sin(t2)*a*sin(t2))));
	return z;
}



double A(double z1, double z2, double z3) {

	double A = 0.5*sqrt(z1*z1*z2*z2-(z1*z1+z2*z2-z3*z3)*(z1*z1+z2*z2-z3*z3)/4.0); //Calculates the area of the triangle for given parameters
	return A;
	
}


int main() {
	
	int n, m;		//n is number of theta values to test, m is no. of bins in histogram
	double a;		//a is the length of the semi major axis divided by that of the semi minor axis
	
	cout << "Enter aspect ratio a" << endl;
	cin >> a;
	cout << "Enter number of theta values" << endl;
	cin >> n; 
	cout << "Enter the number of histogram bins required" << endl;
	cin >> m;
	
	double dtheta = 2.0*3.1415926/(double)n;			//Increment in theta
	
	double dA = 3.1415926*a/(double)m;				//This is used to calculate which histogram bin an A value should go into. It is essentially the width of each bin in terms of z, taking the max possible value of the area, and hence max bin value, as the area of the ellipse itself, pi*a
		
	int output[m];									//index of a bin in the array correponds to the lower bound of the histogram bin divided by dz
	
	for(int x=0; x<m; x++) {
		output[x]=0;
	}
		
	int x = 0;
	
	for(int i=0; i<n; i++) {						//i and j correspond to theta1 and theta2
		for(int j=i+1; j<n; j++) {
			for(int k=j+1; k<n; k++) {
				int l=0;
				double t1=(double)i*dtheta;
				double t2=(double)j*dtheta;
				double t3=(double)k*dtheta;
				double z1 = z(t1, t2, a);
				double z2 = z(t2, t3, a);
				double z3 = z(t3, t1, a);
				double As = A(z1, z2, z3);
				x++;
				for(double Ass = As; Ass>dA; Ass = Ass - dA) {
					l++;
				}
				output[l] = output[l] +1;
			}
		}	
	}
	cout << x << endl;
	ofstream fout;
	fout.open("AreaHistogramOutput.txt");			//Prints data to pre existing file HistogramOutput.txt
	for(int i=0; i<m; i++) {
		fout << output[i] << endl;
	}
}
	

