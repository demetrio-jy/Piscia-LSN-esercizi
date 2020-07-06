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

double psi(double x, double sigma2, double mu) {
	return exp(-pow(x - mu, 2) / (2 * sigma2)) + exp(-pow(x + mu, 2) / (2 * sigma2));
}

double V(double x) {
	return pow(x, 4) - 2.5 * pow(x, 2);
}

double Hpsi(double x, double sigma2, double mu) {
	double e1 = exp(-pow(x - mu, 2) / (2 * sigma2));
	double e2 = exp(-pow(x + mu, 2) / (2 * sigma2));
	return (e1 / sigma2 - pow(x-mu,2) / pow(sigma2,2) * e1 + e2 / sigma2 - pow(x + mu, 2) / pow(sigma2, 2) * e2)/2 + V(x) * psi(x, sigma2, mu);
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

	// EXERCISE 8.1 and 8.2 -------------------------------------------------------------------------
	int M = 100000;		//# steps in metropolis
	int N = 100;		//# blocks
	int L = int(M / N);	//# elements in each block

	double x_old;		//variables for metropolis
	double x_new;
	double jump_max = 2;
	double alpha = 0;
	double R = 0;
	double accept = 0;
	double total = 0;
	double conf = 0;

	double mu_i = 0.52;	//variables for finding optimal mu and sigma2
	double sigma2_i = 0.22;
	double mu = 0;
	double sigma2 = 0;
	double min = 0;
	double best_mu = 0; 
	double best_sigma2 = 0;

	for (int l = 0; l < 45; l++) {
		mu = mu_i + l * 0.01;
		for (int m = 0; m < 45; m++) {
			sigma2 = sigma2_i + m * 0.01;
			x_old = 2.5;
			double* x = new double[M] {0};
			double* ave = new double[N] {0};
			double* av2 = new double[N] {0};
			double* sum_prog = new double[N] {0};
			double* su2_prog = new double[N] {0};
			double* err_prog = new double[N] {0};
			for (int j = 0; j < M; j++) {
				x_new = x_old + rnd.Rannyu(-jump_max, jump_max);
				alpha = pow(psi(x_new, sigma2, mu), 2) / pow(psi(x_old, sigma2, mu), 2);
				if (alpha > 1)
					alpha = 1;
				R = rnd.Rannyu();
				if (R <= alpha) {
					x_old = x_new;
					accept++;
				}
				total++;
				x[j] = x_old;
			}
			conf = accept / total;
			cout << "The accepted jumps over the total proposed moves are " << conf << endl;
			for (int i = 0; i < N; i++) {
				double sum = 0;
				for (int j = 0; j < L; j++) {
					int k = j + i * L;
					sum = sum + Hpsi(x[k], sigma2, mu) / psi(x[k], sigma2, mu);
				}
				ave[i] = sum / L;
				av2[i] = (ave[i]) * (ave[i]);
			}
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < i + 1; j++) {
					sum_prog[i] = sum_prog[i] + ave[j];
					su2_prog[i] = su2_prog[i] + av2[j];
				}
				sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
				su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
				err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
			}
			if (sum_prog[N - 1] < min) {
				min = sum_prog[N - 1];
				best_mu = mu;
				best_sigma2 = sigma2;
			}
			accept = 0; 
			total = 0;
		}
	}

	cout << endl << "The optimal parameters for the wave function are mu = " << best_mu << ", sigma2 = " << best_sigma2 << endl << endl;

	//The metropolis in now repeated with optimal parameters
	double* x = new double[M] {0};
	double* ave = new double[N] {0};
	double* av2 = new double[N] {0};
	double* sum_prog = new double[N] {0};
	double* su2_prog = new double[N] {0};
	double* err_prog = new double[N] {0};
	ofstream data1("mean.dat");
	ofstream data2("configuration.dat");
	for (int j = 0; j < M; j++) {
		x_new = x_old + rnd.Rannyu(-jump_max, jump_max);
		alpha = pow(psi(x_new, best_sigma2, best_mu), 2) / pow(psi(x_old, best_sigma2, best_mu), 2);
		if (alpha > 1)
			alpha = 1;
		R = rnd.Rannyu();
		if (R <= alpha) {
			x_old = x_new;
			accept++;
		}
		total++;
		x[j] = x_old;
		data2 << x[j] << endl;
	}
	conf = accept / total;
	cout << "The accepted jumps over the total proposed moves are " << conf << endl;
	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + Hpsi(x[k], best_sigma2, best_mu) / psi(x[k], best_sigma2, best_mu);
		}
		ave[i] = sum / L;
		av2[i] = (ave[i]) * (ave[i]);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];
			su2_prog[i] = su2_prog[i] + av2[j];
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
		data1 << sum_prog[i] << " \t" << err_prog[i] << endl;
	}
	data1.close();
	data2.close();

	rnd.SaveSeed();
	return 0;
}