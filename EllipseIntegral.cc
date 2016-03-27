/*
 *  EllipseIntegral.cc
 *  
 *
 *  Created by Christopher on 23/08/2010.
 *  Copyright 2010 University of Cambridge. All rights reserved.
 *
 */

using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>

double zfun(double r1, double r2, double a) {	//Evaluates length of line between points
	
	double zfun = sqrt(r1*r1 + r2*r2 - 2.0*(a*a*sqrt(fabs((r2*r2-1.0)*(r1*r1-1.0))) + sqrt(fabs((a*a-r2*r2)*(a*a-r1*r1)))));
	return zfun;
}

double pfun(double a, double r1, double r2) {		//Evaluates P(rho1)*P(rho2)
	
	double pfun = a*a/(4.0*9.8696044*r1*r2*sqrt(fabs((a*a-r1*r1)*(a*a-r2*r2)*(r1*r1-1.0)*(r2*r2-1.0))));
	return pfun;
	
}

int main() {
	
	int n,m;			//n is number of rho values to test, m is the number of z values to test
	double a, e;		//a is the length of the semi major axis divided by that of the semi minor axis
	
	cout << "Enter aspect ration a" << endl;
	cin >> a;
	cout << "Enter number of z values" << endl;
	cin >> n;
	cout << "Enter the number of rho values in grid" << endl;
	cin >> m;
	cout << "Enter acceptable error on delta function" << endl;
	cin >> e;
		
	double drho = (a-1.0)/(double)m;			//Increment in rho
	
	double output[n+1];						//Array to store output data for plotting
	
	double dz = 2.0*a/(double)n;				//Increment in z
	
	for(int i=0; i<=n; i++) {
		
		double z = (double)i*dz;			//Convert i into corresponding value of z by multiplication by increment of z
		
		double p = 0.0;						//Value for P(z) at this value of z, to be updated to correct value by procedure below
		
		for(int j=0; j<=m; j++) {
						
			for(int k=0; k<=m; k++) {					//i and j correspond to rho1 and rho2 by multiplication by increment of rho
				
				double x = (double)j*drho;
				double y = (double)k*drho;
				
				cout << pfun(a,x,y) << "	" << (zfun(x,y,a)) << endl;
				
				if(fabs((z-zfun(x, y, a)))<=e) {				//If delta function is approximately zero then add element of integral sum
					p = p + pfun(a, x, y)*drho*drho;
				}
				
				else {											//otherwise add nothing to integral sum
				}	
			}
		}
		
		output[i]=p;
	}
	
	ofstream fout;
	fout.open("IntegralOutput.txt");		//Prints data to pre existing file IntegralOutput.txt
	for(int i=0; i<n; i++) {
		fout << output[i] << endl;
	}
}
		
		
		
		
		
		

