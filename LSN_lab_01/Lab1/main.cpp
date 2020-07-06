#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

double error(double* AV, double* AV2, int n) {
	if (n == 0)
		return 0;
	else
		return sqrt((AV2[n] - AV[n] * AV[n]) / n);
}

int main() {

	//Procedure for random number generation 
	Random rnd;
	int seed[4];
	int p1;
	int p2;
	ifstream Primes("Primes");
	if (Primes.is_open())
		Primes >> p1 >> p2;
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	ifstream input("seed.in");
	string property;
	if (input.is_open()) {
		while (!input.eof()) {
			input >> property;
			if (property == "RANDOMSEED") {
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed, p1, p2);
			}
		}
		input.close();
	}
	else cerr << "PROBLEM: Unable to open seed.in" << endl;

	// EXERCISE 1.1.a ----------------------------------------------------------------------------------------
	int M = 10000 ;              // Total number of throws
	int N = 100;                 // Number of blocks
	int L = int(M / N);          // Number of throws in each block
	double* r = new double [M] ; // Array for measures
	for (int i = 0; i < M; i++) 
		r[i] = rnd.Rannyu();
	double* x = new double [N];  // Variables for blocking method
	for (int j = 0; j < N; j++)
		x[j] = j;
	double* ave = new double[N] {0};
	double* av2 = new double[N] {0};
	double* sum_prog = new double[N] {0};
	double* su2_prog = new double[N] {0};
	double* err_prog = new double[N] {0};

	for (int i = 0; i < N; i++) { // Blocking method
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + r[k];
		}
		ave[i] = sum / L;				// r_i
		av2[i] = (ave[i]) * (ave[i]);	 // (r_i) ^ 2
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];	// SUM_{ j = 0,i } r_j
			su2_prog[i] = su2_prog[i] + av2[j];	// SUM_{ j = 0,i } (r_j) ^ 2
		}
		sum_prog[i] = sum_prog[i]/(i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i]/(i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream data1("data1.dat");
	for (int i = 0; i < N; i++)
		data1 << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
	data1.close();

	// EXERCISE 1.1.b ---------------------------------------------------------------------------------------------
	for (int i = 0; i < N; i++) { // Blocking method
		ave[i] = 0;
		av2[i] = 0;
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		err_prog[i] = 0;
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + (r[k]-0.5)* (r[k] - 0.5);
		}
		ave[i] = sum / L;				// r_i
		av2[i] = (ave[i]) * (ave[i]);	 // (r_i) ^ 2
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];	// SUM_{ j = 0,i } r_j
			su2_prog[i] = su2_prog[i] + av2[j];	// SUM_{ j = 0,i } (r_j) ^ 2
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream data2("data2.dat");
	for (int i = 0; i < N; i++)
		data2 << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
	data1.close();


	// EXERCISE 1.1.c ---------------------------------------------------------------------------------------------
	int step;
	int k = 0;
	int Ncycles = 100;
	int Nstep = 100;
	double* chi = new double[Nstep] {0};
	double* counter = new double[Nstep] {0}; // Array with the number of events in each interval
	while (k < Ncycles) {
		for (int y = 0;y < Nstep;y++)
			counter[y] = 0;        // Initialized to 0 at the beginning of each cycle
		for (int i = 0; i < M; i++) {
			r[i] = rnd.Rannyu();   // Redo the experiment in order to have indipendent results
			step = 0;
			while ((double)step / 100. < 1) {
				if (r[i] > (double)step / 100. && r[i] < ((double)step + 1.) / 100.) {
					counter[step]++;
					break;
				}
				else
					step++;
			}
		}
		for (int i = 0;i < 100;i++)
			chi[k] = chi[k] + pow((counter[i] - 100), 2) / 100;    // Calculate chi
		k++;
	}

	ofstream data3("data3.dat");
	for (int i = 0;i < N;i++)
		data3 << x[i] + 1 << "\t" << chi[i] << endl;
	data3.close();


	// EXERCISE 1.2 -------------------------------------------------------------------------------------------------------
	double Lambda = 1;	// Parameter exp
	double Mu = 0;		// Parameter lorentzian
	double Gamma = 1;
	M = pow(10, 4);		// # realizations
	double* s1 = new double[M] {0};
	double* s2 = new double[M] {0};
	double* s10 = new double[M] {0};
	double* s100 = new double[M] {0};

	// Standard dice
	for (int i = 0; i < M; i++) { 
		s1[i] = floor(rnd.Rannyu(1., 7.));
		s2[i] = floor(rnd.Rannyu(1., 7.)) + floor(rnd.Rannyu(1., 7.));
		for (int j = 0; j < 10; j++) {
			s10[i] += floor(rnd.Rannyu(1., 7.));
			for (int k = 0; k < 10; k++)
				s100[i] += floor(rnd.Rannyu(1., 7.));
		}
		s2[i] /= 2;
		s10[i] /= 10;
		s100[i] /= 100;
	} 
	ofstream ndice("ndice.dat");
	for (int i = 0; i < M; i++)
		ndice << s1[i] << "\t" << s2[i] << " \t" << s10[i] << " \t" << s100[i] << endl;
	ndice.close();

	// Exponential dice
	for (int i = 0; i < M; i++) {
		s1[i] = rnd.Exp(Lambda);
		s2[i] = rnd.Exp(Lambda) + rnd.Exp(Lambda);
		s10[i] = 0;
		s100[i] = 0;
		for (int j = 0; j < 10; j++) {
			s10[i] += rnd.Exp(Lambda);
			for (int k = 0; k < 10; k++)
				s100[i] += rnd.Exp(Lambda);
		}
		s2[i] /= 2;
		s10[i] /= 10;
		s100[i] /= 100;
	}
	ofstream edice("edice.dat");
	for (int i = 0; i < M; i++)
		edice << s1[i] << "\t" << s2[i] << " \t" << s10[i] << " \t" << s100[i] << endl;
	edice.close();

	// Lorentzian dice
	for (int i = 0; i < M; i++) {
		s1[i] = rnd.Lorentz(Gamma, Mu);
		s2[i] = rnd.Lorentz(Gamma, Mu) + rnd.Lorentz(Gamma, Mu);
		s10[i] = 0;
		s100[i] = 0;
		for (int j = 0; j < 10; j++) {
			s10[i] += rnd.Lorentz(Gamma, Mu);
			for (int k = 0; k < 10; k++)
				s100[i] += rnd.Lorentz(Gamma, Mu);
		}
		s2[i] /= 2;
		s10[i] /= 10;
		s100[i] /= 100;
	}
	ofstream ldice("ldice.dat");
	for (int i = 0; i < M; i++)
		ldice << s1[i] << "\t" << s2[i] << " \t" << s10[i] << " \t" << s100[i] << endl;
	ldice.close();


	// EXERCISE 1.3 ------------------------------------------------------------------------------------------------------
	M = 10000;					// Number of experiments
	double Nthrows = 10000;     // Number of throws in each experiment
	double Lgt = (double)1/12;	// Length of the needle (we need L<1)
	double d = (double)1/8;     // Space between vertical lines (we need d>L)
	int Hit;
	double X1, theta, PX, X2;
	int j;						// Number of lines: tot 8 lines
	double* pi = new double[M]; // Extimation of pi
	for (int k = 0;k < M;k++) {
		Hit = 0;
		for (int i = 0; i < Nthrows; i++) {
			X1 = rnd.Rannyu();				  // Starter point 
			theta = rnd.Rannyu(0, 2 * M_PI);  // Random angle
			PX = Lgt * cos(theta);			  // Projection on x axes
			X2 = X1 + PX;					  // End point
			j = 0;
			while (j * d <= 1) {
				if (X1 <= j * d && X2 >= j * d)       // Hit-or-Miss test
					Hit++;
				else if (X1 >= j * d && X2 <= j * d)
					Hit++;
				j++;
			}
		}
		pi[k] = 2 * Lgt * Nthrows / (Hit * d);  
	}

	// Blocking method
	N = 100;
	L = M / N;
	x = new double[N];
	for (int i = 0; i < N; i++)
		x[i] = L * i;
	ave = new double[N] {0};
	av2 = new double[N] {0};
	sum_prog = new double[N] {0};
	su2_prog = new double[N] {0};
	err_prog = new double[N] {0};
	for (int i = 0;i < N; i++) {
		double sum = 0;
		for (int j = 0;j < L;j++) {
			int k = j + i * L;
			sum = sum + pi[k];
		}
		ave[i] = sum / (double)L;
		av2[i] = pow(ave[i], 2);
	}

	for (int i = 0;i < N;i++) {
		for (int j = 0;j < i + 1;j++) {
			sum_prog[i] = sum_prog[i] + ave[j];
			su2_prog[i] = su2_prog[i] + av2[j];
		}
		sum_prog[i] = sum_prog[i] / ((double)(i)+1);
		su2_prog[i] = su2_prog[i] / ((double)(i)+1);
		err_prog[i] = error(sum_prog, su2_prog, i);
	}

	ofstream buffon("buffon.dat");
	for (int i = 0;i < N;i++)
		buffon << x[i] << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
	buffon.close();


	rnd.SaveSeed();
	return 0;
}