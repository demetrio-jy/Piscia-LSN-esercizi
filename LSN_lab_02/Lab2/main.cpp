#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "position.h"

using namespace std;

double error(double* AV, double* AV2, int n) {
	if (n == 0)
		return 0;
	else
		return sqrt((AV2[n] - AV[n] * AV[n]) / n);
}

double error(double AV, double AV2, int n) {
	if (n == 0)
		return 0;
	else
		return sqrt((AV2 - AV * AV) / n);
}

double f1(double x) {		//function for ex. 2.1.1
	return M_PI / 2 * cos(x * M_PI / 2);
}

double f2(double x) {		//function for ex. 2.1.2
	return M_PI / 2 * cos(x * M_PI / 2) / (2 - 2 * x);
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

	// EXERCISE 2.1.1 ----------------------------------------------------------------------------------------
	int M = 100000;             // Total number of throws
	int N = 100;                 // Number of blocks
	int L = int(M / N);          // Number of throws in each block
	double* x = new double[N];
	for (int j = 0; j < N; j++)
		x[j] = j;
	double* ave = new double[N] {0};
	double* av2 = new double[N] {0};
	double* sum_prog = new double[N] {0};
	double* su2_prog = new double[N] {0};
	double* err_prog = new double[N] {0};

	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + f1(rnd.Rannyu());
		}
		ave[i] = sum / L;				// function integrated 
		av2[i] = (ave[i]) * (ave[i]);
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];	// SUM_{ j = 0,i } f1
			su2_prog[i] = su2_prog[i] + av2[j];	// SUM_{ j = 0,i } f1 ^ 2
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream integral("integral.dat");
	for (int i = 0; i < N; i++)
		integral << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
	integral.close();

	// EXERCISE 2.1.2 ----------------------------------------------------------------------------------------
	for (int i = 0; i < N; i++) {
		ave[i] = 0;
		av2[i] = 0;
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + f2(rnd.Linear());
		}
		ave[i] = sum / L;				// f2 with prob. distr. 2(1-x)
		av2[i] = (ave[i]) * (ave[i]);
	}

	for (int i = 0; i < N; i++) {
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];	// SUM_{ j = 0,i } f2
			su2_prog[i] = su2_prog[i] + av2[j];	// SUM_{ j = 0,i } f2 ^ 2
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream importance("importance.dat");
	for (int i = 0; i < N; i++)
		importance << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
	importance.close();

	// EXERCISE 2.2.1 ----------------------------------------------------------------------------------------
	M = 10000;	//# experiments
	N = 100;	//# blocks
	L = M / N;	//# rw x block
	int Nstep = 100;
	int dice;
	Position pos;
	double* aveRW = new double[N * Nstep]{ 0 };
	double* av2RW = new double[N * Nstep]{ 0 };
	double* Rsquared_prog = new double[Nstep] {0};
	double* Rsquared2_prog = new double[Nstep] {0};
	double* err_progRW = new double[Nstep] {0};

	for(int k = 0; k < N; k++){ //sum over block number
		for (int i = 0; i < L; i++) { //sum over experiments
			pos.SetX(0);
			pos.SetY(0);
			pos.SetZ(0);
			for (int j = 0; j < Nstep; j++) { //sum over steps
				dice = floor(rnd.Rannyu(1, 7));
				if (dice == 1)
					pos.SetX(pos.GetX() + 1);
				else if (dice == 2)
					pos.SetX(pos.GetX() - 1);
				else if (dice == 3)
					pos.SetY(pos.GetY() + 1);
				else if (dice == 4)
					pos.SetY(pos.GetY() - 1);
				else if (dice == 5)
					pos.SetZ(pos.GetZ() + 1);
				else if (dice == 6)
					pos.SetZ(pos.GetZ() - 1);
				aveRW[j + k * Nstep] += pow(pos.GetR(), 2) / L;
				av2RW[j + k * Nstep] += pow(pos.GetR(), 4) / L;
			 }
		}
	} 

	ofstream RW1("RW1.dat");

	for (int k = 0; k < Nstep; k++) {
		Rsquared_prog[k] = 0.;
		Rsquared2_prog[k] = 0.;
		for (int i = 0; i < N; i++) {
			Rsquared_prog[k] += aveRW[i * Nstep + k];
			Rsquared2_prog[k] += av2RW[i * Nstep + k];
		}
		Rsquared_prog[k] /= N;
		Rsquared2_prog[k] /= N;

		err_progRW[k] = error(Rsquared_prog[k], Rsquared2_prog[k], N);

		RW1 << Rsquared_prog[k] << "\t" << err_progRW[k] << endl;
	}

	RW1.close();

	// EXERCISE 2.2.2 ----------------------------------------------------------------------------------------
	double t = 0;
	double f = 0;

	for (int i = 0; i < Nstep * N; i++) {
		aveRW[i] = 0.;
		av2RW[i] = 0.;
	}

	for (int k = 0; k < N; k++) { //sum over block number
		for (int i = 0; i < L; i++) { //sum over experiments
			pos.SetX(0);
			pos.SetY(0);
			pos.SetZ(0);
			for (int j = 0; j < Nstep; j++) { //sum over steps
				t = rnd.Theta();
				f = rnd.Phi();
				pos.SetX(pos.GetX() + sin(f) * cos(t));
				pos.SetY(pos.GetY() + sin(f) * sin(t));
				pos.SetZ(pos.GetZ() + cos(f));
				aveRW[j + k * Nstep] += pow(pos.GetR(), 2) / L;
				av2RW[j + k * Nstep] += pow(pos.GetR(), 4) / L;
			}
		}
	}

	ofstream RW2("RW2.dat");

	for (int k = 0; k < Nstep; k++) {
		Rsquared_prog[k] = 0.;
		Rsquared2_prog[k] = 0.;
		for (int i = 0; i < N; i++) {
			Rsquared_prog[k] += aveRW[i * Nstep + k];
			Rsquared2_prog[k] += av2RW[i * Nstep + k];
		}
		Rsquared_prog[k] /= N;
		Rsquared2_prog[k] /= N;

		err_progRW[k] = error(Rsquared_prog[k], Rsquared2_prog[k], N);

		RW2 << Rsquared_prog[k] << "\t" << err_progRW[k] << endl;
	}

	RW2.close();


	rnd.SaveSeed();
	return 0;
}