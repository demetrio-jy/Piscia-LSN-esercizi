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

double psq_100(Position t) {
	return exp(-2*t.GetR());
}

double psq_210(Position t) {
	return pow(t.GetR() * cos(t.GetT()), 2) * exp(-t.GetR());
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

	// EXERCISE 5.1 : psi_100 uniform -------------------------------------------------------------------------
	int M = 100000;		//# steps in metropolis
	int N = 1000;		//# blocks
	int L = int(M / N);	//# elements in each block

	Position p;			
	Position jump;
	double jump_max = 1.223;
	double alpha = 0;
	double R = 0;
	double accept = 0;
	double reject = 0;
	double conf = 0;

	double* r = new double[M] {0};
	double* ave = new double[N] {0};
	double* av2 = new double[N] {0};
	double* sum_prog = new double[N] {0};
	double* su2_prog = new double[N] {0};
	double* err_prog = new double[N] {0};

	p.SetX(10);
	p.SetY(10);
	p.SetZ(10);

	for (int j = 0; j < M; j++) {
		jump.SetX(p.GetX() + rnd.Rannyu(-jump_max, jump_max));
		jump.SetY(p.GetY() + rnd.Rannyu(-jump_max, jump_max));
		jump.SetZ(p.GetZ() + rnd.Rannyu(-jump_max, jump_max));
		alpha = psq_100(jump) / psq_100(p);
		if (alpha > 1)
			alpha = 1;
		R = rnd.Rannyu();
		if (R <= alpha) {
			p = jump;
			accept++;
		}
		else {
			reject++;
		}
		r[j] = p.GetR();
	}

	conf = accept / reject;
	cout << "The ratio between accepted and rejected jumps in 'Psi_100 uniform' is " << conf << endl;

	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + r[k];
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

	ofstream data1("data1.dat");
	for (int i = 0; i < N; i++)
		data1 << sum_prog[i] << " \t" << err_prog[i] << endl;
	data1.close();



	// EXERCISE 5.1 : psi_210 uniform -------------------------------------------------------------------------
	accept = 0;
	reject = 0;
	jump_max = 3.18;

	p.SetX(10);
	p.SetY(10);
	p.SetZ(10);

	for (int j = 0; j < M; j++) {
		jump.SetX(p.GetX() + rnd.Rannyu(-jump_max, jump_max));
		jump.SetY(p.GetY() + rnd.Rannyu(-jump_max, jump_max));
		jump.SetZ(p.GetZ() + rnd.Rannyu(-jump_max, jump_max));
		alpha = psq_210(jump) / psq_210(p);
		if (alpha > 1)
			alpha = 1;
		R = rnd.Rannyu();
		if (R <= alpha) {
			p = jump;
			accept++;
		}
		else {
			reject++;
		}
		r[j] = p.GetR();
	}

	conf = accept / reject;
	cout << "The ratio between accepted and rejected jumps in 'Psi_210 uniform' is " << conf << endl;

	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + r[k];
		}
		ave[i] = sum / L;
		av2[i] = (ave[i]) * (ave[i]);
	}

	for (int i = 0; i < N; i++) {
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];
			su2_prog[i] = su2_prog[i] + av2[j];
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream data2("data2.dat");
	for (int i = 0; i < N; i++)
		data2 << sum_prog[i] << " \t" << err_prog[i] << endl;
	data2.close();



	// EXERCISE 5.1 : psi_100 gaussian -------------------------------------------------------------------------
	accept = 0;
	reject = 0;
	double sigma = 0.764;

	p.SetX(10);
	p.SetY(10);
	p.SetZ(10);

	for (int j = 0; j < M; j++) {
		jump.SetX(rnd.Gauss(p.GetX(), sigma));
		jump.SetY(rnd.Gauss(p.GetY(), sigma));
		jump.SetZ(rnd.Gauss(p.GetZ(), sigma));
		alpha = psq_100(jump) / psq_100(p);
		if (alpha > 1)
			alpha = 1;
		R = rnd.Rannyu();
		if (R <= alpha) {
			p = jump;
			accept++;
		}
		else {
			reject++;
		}
		r[j] = p.GetR();
	}

	conf = accept / reject;
	cout << "The ratio between accepted and rejected jumps in 'Psi_100 gaussian' is " << conf << endl;

	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + r[k];
		}
		ave[i] = sum / L;
		av2[i] = (ave[i]) * (ave[i]);
	}

	for (int i = 0; i < N; i++) {
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];
			su2_prog[i] = su2_prog[i] + av2[j];
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream data3("data3.dat");
	for (int i = 0; i < N; i++)
		data3 << sum_prog[i] << " \t" << err_prog[i] << endl;
	data3.close();



	// EXERCISE 5.1 : psi_210 gaussian -------------------------------------------------------------------------
	accept = 0;
	reject = 0;
	sigma = 1.996;

	p.SetX(10);
	p.SetY(10);
	p.SetZ(10);

	for (int j = 0; j < M; j++) {
		jump.SetX(rnd.Gauss(p.GetX(), sigma));
		jump.SetY(rnd.Gauss(p.GetY(), sigma));
		jump.SetZ(rnd.Gauss(p.GetZ(), sigma));
		alpha = psq_210(jump) / psq_210(p);
		if (alpha > 1)
			alpha = 1;
		R = rnd.Rannyu();
		if (R <= alpha) {
			p = jump;
			accept++;
		}
		else {
			reject++;
		}
		r[j] = p.GetR();
	}

	conf = accept / reject;
	cout << "The ratio between accepted and rejected jumps in 'Psi_210 gaussian' is " << conf << endl;

	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			int k = j + i * L;
			sum = sum + r[k];
		}
		ave[i] = sum / L;
		av2[i] = (ave[i]) * (ave[i]);
	}

	for (int i = 0; i < N; i++) {
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		for (int j = 0; j < i + 1; j++) {
			sum_prog[i] = sum_prog[i] + ave[j];
			su2_prog[i] = su2_prog[i] + av2[j];
		}
		sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
		su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
		err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
	}

	ofstream data4("data4.dat");
	for (int i = 0; i < N; i++)
		data4 << sum_prog[i] << " \t" << err_prog[i] << endl;
	data4.close();

	rnd.SaveSeed();
	return 0;
}