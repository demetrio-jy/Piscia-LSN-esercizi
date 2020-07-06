#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "individuo.h"
#include "functions.h"

using namespace std;

int main() {
	RandomInitializer(); //Procedure for random number generation 
	PositionGenerator(); //Procedure for creating the positions of the cities
	PopulationGenerator(); //Procedure for creating a random population

    Individual* nextGen = new Individual[Npopulation];

    for (int y = 0; y < iterations ; y++) {
        FitnessCalculator();
        Quicksort(population, 0, Npopulation - 1);
        Minimum(y);
        if (y < iterations / 100)
            prob = 0.14;
        else prob = 0.02;
        int check1 = 0;
        int check2 = 0;
        for (int w = 0; w < Npopulation; w = w++) {
            Individual parent1 = Selector();
            Individual parent2 = Selector();
            double r = rnd.Rannyu();
            if (r <= prob) {
                PairMutation(parent1);
            } else if (r <= prob*2) {
                ShiftMutation(parent1);
            } else if (r <= prob*3) {
                PermutationMutation(parent1);
            } else if (r <= prob*4) {
                InversionMutation(parent1);
            } else {
                CrossingOver(parent1, parent2);
            }
            check1 = CheckFunction(parent1);
            check2 = CheckFunction(parent2);
            if (check1 == 1 || check2 == 1) {
                cout << "Error in the creation of an individual" << endl;
            }
            nextGen[w].SetCities(parent1.GetCities());
            if (r > prob*4) {
                w++;
                nextGen[w].SetCities(parent2.GetCities());
            }
        }
        for (int k = 0; k < Npopulation; k++)
            population[k].SetCities(nextGen[k].GetCities());
    }
    Quicksort(population, 0, Npopulation - 1);
    Minimum(iterations);
    ofstream pos;
    pos.open("data_circ_GA.dat");
    for (int k = 0; k < Ncities; k++)
        pos << positions[best.GetCity(k)-1][0] << "\t" << positions[best.GetCity(k)-1][1] << "\t" << best.GetCity(k) << endl;
    pos.close();
    ofstream L;
    L.open("L_circ_GA.dat");
    for (int k = 0; k <= iterations; k++)
        L << Lbest[k] << "\t" << Lmean[k] << "\t" << Lerr[k] << endl;
    L.close();

    cout << "The minimum L(1) found is " << Lmin << endl;

	rnd.SaveSeed();
	return 0;
}

void RandomInitializer() { //procedure for using random class
	ifstream Primes("Primes");
	if (Primes.is_open()) {
		Primes >> p1 >> p2;
	}
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
	return;
}

void PositionGenerator() {
	if (method == 0) { //random positions on a circumference
		double theta; 
		for (int i = 0; i < Ncities; i++) {
			theta = rnd.Theta();
			positions[i][0] = cos(theta);
			positions[i][1] = sin(theta);
		}
	}
	else { //random positions in a square
		for (int i = 0; i < Ncities; i++) {
			positions[i][0] = rnd.Rannyu();
			positions[i][1] = rnd.Rannyu();
		}
	}
	return;
}

int CheckFunction(Individual I){ //return 1 if I doesn't fulfil the bonds
    int flag = 0;
    for (int i = 1; i < Ncities; i++) {
        for (int j = 1; j < i; j++) {
            if (I.GetCity(j) == I.GetCity(i))
                flag = 1;
        }
    }
    return flag;
}

void PopulationGenerator() { //random creation of a population of individuals
	for (int i = 0; i < Npopulation; i++) { //cycle over the population
		int flag = 1; //1 enter the while
		Individual appo;
		for (int j = 1; j < Ncities; j++) {
			while (flag == 1) {
				appo.SetCity(j, floor(rnd.Rannyu(2., 33.)));
				flag = 0;
				for (int k = 1; k < j; k++) {
					if (appo.GetCity(j) == appo.GetCity(k))
						flag = 1;
				}
			}
			flag = 1;
		}
		population[i] = appo;
	}
	return;
}

double L(Individual I) { //compute the distance between adiacent cities: I use L(1)
    double L = 0;
    int j;
    for (int i = 0; i < Ncities; i++) {
        if (i == Ncities - 1)
            j = 0;
        else
            j = i + 1;
        L += sqrt(pow(positions[I.GetCity(i)-1][0] - positions[I.GetCity(j)-1][0], 2) + pow(positions[I.GetCity(i)-1][1] - positions[I.GetCity(j)-1][1], 2));
    }
	return L;
}

void FitnessCalculator() { //compute L for the entire population 
    for (int k = 0; k < Npopulation; k++) {
        double appo = L(population[k]);
        population[k].SetL(appo);
    }
    return;
}

void Quicksort(Individual* I, int first, int last) {
    int i, j, pivot;
    Individual appo;
    if (first < last) {
        pivot = first;
        i = first;
        j = last;
        while (i < j) {
            while (I[i].GetL() <= I[pivot].GetL() && i < last)
                i++;
            while (I[j].GetL() > I[pivot].GetL())
                j--;
            if (i < j) {
                appo = I[i];
                I[i] = I[j];
                I[j] = appo;
            }
        }
        appo = I[pivot];
        I[pivot] = I[j];
        I[j] = appo;
        Quicksort(I, first, j - 1);
        Quicksort(I, j + 1, last);
    }
}

Individual Selector() { //use only the best half of the population
    double r = rnd.Rannyu();
    int p = 5.;
    int i = (int)floor(Npopulation/2. * pow(r, p));
    return population[i];
}

void PairMutation(Individual I) { //switch two random positions
    int pos1 = (int)floor(rnd.Rannyu(1., 32.));
    int pos2 = (int)floor(rnd.Rannyu(1., 32.));
    while (pos2 == pos1)
        pos2 = (int)floor(rnd.Rannyu(1., 32.));
    int appo1, appo2;
    appo1 = I.GetCity(pos1);
    appo2 = I.GetCity(pos2);
    I.SetCity(pos1, appo2);
    I.SetCity(pos2, appo1);
    return;
}

void ShiftMutation(Individual I) { //shift of +n positions for m contiguous cities
    int pos = (int)floor(rnd.Rannyu(1., 32.));
    int n = (int)floor(rnd.Rannyu(1., 16.));
    int m = (int)floor(rnd.Rannyu(1., 16.));
    int* appoM = new int[m];
    int* appoN = new int[n];
    for (int k = 0; k < m; k++)
        appoM[k] = I.GetCity(PBC(pos + k));
    for (int k = 0; k < n; k++)
        appoN[k] = I.GetCity(PBC(pos + m + k));
    for (int k = 0; k < m; k++)
        I.SetCity(PBC(pos + k + n), appoM[k]);
    for (int k = 0; k < n; k++)
        I.SetCity(PBC(pos + k), appoN[k]);
    return;
}

void PermutationMutation(Individual I) { //permutation among m contiguous cities with other m contiguous cities
    int n = (int)floor(rnd.Rannyu(1., 16.));
    int pos1 = (int)floor(rnd.Rannyu(1., 32.));
    int pos2 = PBC(pos1 + n + (int)floor(rnd.Rannyu(0,31-2*n)));
    //cout << pos1 << " " << pos2 << " " << n << endl;
    int appo1;
    int appo2;
    for (int k = 0; k < n; k++) {
        appo1 = I.GetCity(PBC(pos1 + k));
        appo2 = I.GetCity(PBC(pos2 + k));
        I.SetCity(PBC(pos1 + k), appo2);
        I.SetCity(PBC(pos2 + k), appo1);
    }
    return;
}

void InversionMutation(Individual I) { //invert the order in a block of m cities
    int pos = (int)floor(rnd.Rannyu(1., 32.));
    int m = (int)floor(rnd.Rannyu(1., 32.));
    int* appo = new int[m];
    for (int k = 0; k < m; k++)
        appo[m - k - 1] = I.GetCity(PBC(pos + k));
    for (int k = pos; k < pos + m; k++)
        I.SetCity(PBC(k), appo[k - pos]);
    return;
}

void CrossingOver(Individual I1, Individual I2) { //mix two individuals 
    int cut = (int)floor(rnd.Rannyu(2., 31.));  //choose a random point for the cut
    int* tail1 = new int[Ncities - cut]; //arrays for the last  
    int* tail2 = new int[Ncities - cut];
    for (int i = 0; i < Ncities - cut; i++) { //at this point are fulled with current positions
        tail1[i] = I1.GetCity(PBC(cut + i));
        tail2[i] = I2.GetCity(PBC(cut + i));
    }
    int* newtail1 = new int[Ncities - cut]; //arrays for the exchange
    int* newtail2 = new int[Ncities - cut];
    int k = 0;
    while (k < Ncities - cut) {
        for (int i = 1; i < Ncities; i++) {
            for (int j = 0; j < Ncities - cut; j++) {
                if (I2.GetCity(i) == tail1[j]) {
                    newtail1[k] = I2.GetCity(i);
                    k++;
                }
            }
        }
    }
    k = 0;
    while (k < Ncities - cut) {
        for (int j = 1;j < Ncities;j++) {
            for (int l = 0;l < Ncities - cut;l++) {
                if (I1.GetCity(j) == tail2[l]) {
                    newtail2[k] = I1.GetCity(j);
                    k++;
                }
            }
        }
    }
    for (int i = 0; i < Ncities - cut; i++) {
        I1.SetCity(cut + i, newtail1[i]);
        I2.SetCity(cut + i, newtail2[i]);
    }
    return;
}

int PBC(int k) { //do not consider the first city
    return k - 31 * rint(k / 32);
}

double Mean() {
    double sum = 0;
    for (int k = 0; k < Npopulation / 2; k++)
        sum += population[k].GetL();
    return sum * 2 / Npopulation;
}

double StdDev() {
    double m = Mean();
    double sum2 = 0;
    for (int k = 0; k < Npopulation / 2; k++)
        sum2 += pow(population[k].GetL()-m,2);
    sum2 = sum2 / ((int)Npopulation/2 -1);
    return sqrt(sum2);
}

void Minimum(int i) { //do generic stuff and save the minimum path
    Lbest[i] = population[0].GetL();
    Lmean[i] = Mean();
    Lerr[i] = StdDev();
    if (Lmin > population[0].GetL()) {
        Lmin = population[0].GetL();
        best = population[0];
    }
}

