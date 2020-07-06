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

	// EXERCISE 3.1.1 ----------------------------------------------------------------------------------------
	double T = 1;		//final time
	double K = 100;		//strike price
	double r = 0.1;		//interest
	double sig = 0.25;	//volatility
	int M = 10000;
	double S_0 = 100;
	double S_T = 0;
	double z = 0;
	double* C = new double[M] {0};
	double* P = new double[M] {0};

	int N = 100;
	int L = M / N;
	double* sum_C = new double[N] {0};
	double* su2_C = new double[N] {0};
	double* err_C = new double[N] {0};
	double* sum_P = new double[N] {0};
	double* su2_P = new double[N] {0};
	double* err_P = new double[N] {0};
	double* avC = new double[N] {0};
	double* avC2 = new double[N] {0};
	double* avP = new double[N] {0};
	double* avP2 = new double[N] {0};
	int* x = new int[N] {0};

	for (int i = 0; i < M; i++) {
		z = rnd.Gauss(0, 1);
		S_T = S_0 * exp((r - pow(sig, 2) / 2) * T + sig * z * sqrt(T));
		if (S_T - K > 0) 
			C[i] = exp(-r * T) * (S_T - K);
		if (K - S_T > 0)
			P[i] = exp(-r * T) * (K - S_T);
	}

	for (int i = 0;i < N;i++)
		x[i] = i + 1;

	for (int i = 0;i < N; i++) {
		double sumC = 0;
		double sumP = 0;
		for (int j = 0;j < L;j++) {
			int k = j + i * L;
			sumC = sumC + C[k];
			sumP = sumP + P[k];
		}
		avC[i] = sumC / (double)L;
		avC2[i] = pow(avC[i], 2);
		avP[i] = sumP / (double)L;
		avP2[i] = pow(avP[i], 2);
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i + 1; j++) {
			sum_C[i] = sum_C[i] + avC[j];	
			su2_C[i] = su2_C[i] + avC2[j];	
			sum_P[i] = sum_P[i] + avP[j];
			su2_P[i] = su2_P[i] + avP2[j];
		}
		sum_C[i] = sum_C[i] / (i + 1); 
		su2_C[i] = su2_C[i] / (i + 1); 
		err_C[i] = error(sum_C, su2_C, i); 
		sum_P[i] = sum_P[i] / (i + 1);
		su2_P[i] = su2_P[i] / (i + 1);
		err_P[i] = error(sum_P, su2_P, i);
	}

	ofstream data1("data1.dat");
	for (int i = 0; i < N; i++)
		data1 << x[i] << "  \t" << sum_C[i] << "   \t" << err_C[i] << "  \t" << sum_P[i] << "  \t" << err_P[i] << endl;
	data1.close();



	// EXERCISE 3.1.2 ----------------------------------------------------------------------------------------
	int B = 100; 
	double* S_t = new double[B] {0};
	S_t[0] = S_0;

	for (int i = 0; i < M; i++) {
		C[i] = 0;
		P[i] = 0;
		for (int j = 1; j < B; j++) {
			z = rnd.Gauss(0, 1);
			S_t[j] = S_t[j-1] * exp((r - pow(sig, 2) / 2) * T / B + sig * z * sqrt(T / B));
		}
		if (S_t[B-1] - K > 0)
			C[i] = exp(-r * T) * (S_t[B-1] - K);
		if (K - S_t[B - 1] > 0)
			P[i] = exp(-r * T) * (K - S_t[B - 1]);
	}


	for (int i = 0;i < N; i++) {
		double sumC = 0;
		double sumP = 0;
		for (int j = 0;j < L;j++) {
			int k = j + i * L;
			sumC = sumC + C[k];
			sumP = sumP + P[k];
		}
		avC[i] = sumC / (double)L;
		avC2[i] = pow(avC[i], 2);
		avP[i] = sumP / (double)L;
		avP2[i] = pow(avP[i], 2);
	}

	for (int i = 0; i < N; i++) {
		sum_C[i] = 0;
		su2_C[i] = 0;
		sum_P[i] = 0;
		su2_P[i] = 0;
		for (int j = 0; j < i + 1; j++) {
			sum_C[i] = sum_C[i] + avC[j];
			su2_C[i] = su2_C[i] + avC2[j];
			sum_P[i] = sum_P[i] + avP[j];
			su2_P[i] = su2_P[i] + avP2[j];
		}
		sum_C[i] = sum_C[i] / (i + 1);
		su2_C[i] = su2_C[i] / (i + 1);
		err_C[i] = error(sum_C, su2_C, i);
		sum_P[i] = sum_P[i] / (i + 1);
		su2_P[i] = su2_P[i] / (i + 1);
		err_P[i] = error(sum_P, su2_P, i);
	}

	ofstream data2("data2.dat");
	for (int i = 0; i < N; i++)
		data2 << x[i] << "  \t" << sum_C[i] << "   \t" << err_C[i] << "  \t" << sum_P[i] << "  \t" << err_P[i] << endl;
	data2.close();

	rnd.SaveSeed();
	return 0;
}