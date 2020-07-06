#include "individuo.h"

Individual::Individual() {
	for (int i = 0;i < Nelem;i++)
		order[i] = i + 1;
}

Individual::~Individual() {}

int Individual::GetCity(int i) {
	return order[i];
}

void Individual::SetCity(int i, int argument) {
	order[i] = argument;
	return;
}

int * Individual::GetCities() {
	return order;
}

void Individual::SetCities(int * w) {
	for (int k = 0; k < Nelem; k++)
		order[k] = w[k];
	return;
}

double Individual::GetL(void) {
	return L;
}

void Individual::SetL(double fitness) {
	L = fitness;
	return;
}
