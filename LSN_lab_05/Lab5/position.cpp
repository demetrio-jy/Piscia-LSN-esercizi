#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "position.h"

using namespace std;

Position :: Position(){
	p_x = 0;
	p_y = 0;
	p_z = 0;
}

Position :: ~Position(){}

void Position :: SetX(double x) {
	p_x = x;
	return;
}

void Position :: SetY(double y) {
	p_y = y;
	return;
}

void Position :: SetZ(double z) {
	p_z = z;
	return;
}

double Position :: GetX() {
	return p_x;
}

double Position :: GetY() {
	return p_y;
}

double Position :: GetZ() {
	return p_z;
}

double Position :: GetR() {
	return sqrt(p_x * p_x + p_y * p_y + p_z * p_z);
}

double Position::GetT() {
	return atan(p_y / p_x);
}