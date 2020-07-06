/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#define _USE_MATH_DEFINES
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main() {
    Input(); //Inizialization
    for (int iblk = 1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for (int istep = 1; istep <= nstep; ++istep)
        {
            Move();
            Measure();
            Accumulate(); //Update block averages
            if (istep % 100 == 0)
                cout << "step:" << istep << endl;
            if (istep == (nstep - 1) && iblk == nblk) { OldFinal(); }
        }
        Averages(iblk);   //Print results for current block

    }
    ConfFinal(); //Write final configuration

    return 0;
}


void Input(void) { //Prepare all stuff for the simulation
    ifstream ReadInput, ReadConf, ReadOld;
    double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator

    ReadInput.open("input.dat"); //Read input

    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart / rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol, 1.0 / 3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;

    //Tail corrections for potential energy and pressure
    vtail = (8.0 * M_PI * rho) / (9.0 * pow(rcut, 9)) - (8.0 * M_PI * rho) / (3.0 * pow(rcut, 3));
    ptail = (32.0 * M_PI * rho) / (9.0 * pow(rcut, 9)) - (16.0 * M_PI * rho) / (3.0 * pow(rcut, 3));
    cout << "Tail correction for the potential energy = " << vtail << endl;
    cout << "Tail correction for the virial           = " << ptail << endl;

    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> nblk;
    ReadInput >> method;     // remember: 0 is used for random method, 1 for 2 configuration method

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();

    //Prepare array for measurements
    iv = 0; //Potential energy
    iw = 1; //Kinetic energy
    n_props = 2; //Number of observables

    //measurement of g(r)
    igofr = 2;
    nbins = 100;
    n_props = n_props + nbins;
    bin_size = (box / 2.0) / (double)nbins;

    if (method == 0) {
        //Read initial configuration: random method
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i = 0; i < npart; ++i) {
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();

        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < npart; ++i) {
            vx[i] = rand() / double(RAND_MAX) - 0.5;
            vy[i] = rand() / double(RAND_MAX) - 0.5;
            vz[i] = rand() / double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim = 0; idim < 3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs = 0.0;
        for (int i = 0; i < npart; ++i) {
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];
            sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
        for (int i = 0; i < npart; ++i) {
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }
    else {
        //Read initial and old configuration
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i = 0; i < npart; ++i) {
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();

        cout << "Read old configuration from file old.0 " << endl << endl;
        ReadOld.open("old.0");
        for (int i = 0; i < npart; ++i) {
            ReadOld >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadOld.close();

        Move();

        double sumv2 = 0.0, T = 0.0, fs = 0.0;
        for (int i = 0; i < npart; ++i) {
            vx[i] = Pbc(x[i] - xold[i]) / delta;
            vy[i] = Pbc(y[i] - yold[i]) / delta;
            vz[i] = Pbc(z[i] - zold[i]) / delta;
            sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        }
        sumv2 /= (double)npart;
        T = sumv2 / 3.0;  
        fs = temp / T;
        for (int i = 0; i < npart; ++i) {
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }
    return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, w, t, vij, wij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  w = 0.0;
  t = 0.0;

  //reset the hystogram of g(r)
  for (int k = igofr; k < igofr + nbins; ++k) walker[k] = 0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     int k = igofr;
     while (dr > (k - 1)* bin_size)
         k++;
     walker[k] = walker[k] + 2;

     if(dr < rcut){
         vij = 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
         wij = 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6);

         // contribution to energy and virial
         v += vij;
         w += wij;
     }
    }          
  }
    walker[iv] = 4.0 * v;
    walker[iw] = 48.0 * w / 3.0;
    return;
}

void Reset(int iblk) //Reset block averages
{

    if (iblk == 1)
    {
        for (int i = 0; i < n_props; ++i)
        {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < n_props; ++i)
    {
        blk_av[i] = 0;
    }
    blk_norm = 0;
}

void Accumulate(void) //Update block averages
{

    for (int i = 0; i < n_props; ++i)
    {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{

    double r, gdir;
    ofstream Gofr, Gave, Epot, Pres;
    const int wd = 12;

    cout << "Block number " << iblk << endl;

    Epot.open("output.epot.0", ios::app);
    Pres.open("output.pres.0", ios::app);
    Gofr.open("output.gofr.0", ios::app);
    Gave.open("output.gave.0", ios::app);

    stima_pot = blk_av[iv] / blk_norm / (double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot * stima_pot;
    err_pot = Error(glob_av[iv], glob_av2[iv], iblk);

    stima_pres = rho * temp + (blk_av[iw] / blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres * stima_pres;
    err_press = Error(glob_av[iw], glob_av2[iw], iblk);

    //Potential energy per particle
    Epot << setw(wd) << iblk << setw(wd) << stima_pot << setw(wd) << glob_av[iv] / (double)iblk << setw(wd) << err_pot << endl;
    //Pressure
    Pres << setw(wd) << iblk << setw(wd) << stima_pres << setw(wd) << glob_av[iw] / (double)iblk << setw(wd) << err_press << endl;

    //g(r)
    Gofr << iblk << setw(wd);
    for (int k = 0;k < nbins;k++) {
        double DV = 4 / 3 * M_PI * (pow((k + 1) * bin_size, 3) - pow(k * bin_size, 3));
        blk_av[igofr + k] /= (rho * npart * DV);
        Gofr << blk_av[igofr + k] / blk_norm << setw(wd);
    }
    Gofr << endl;

    for (int k = 0;k < nbins;k++) {
        glob_av[igofr + k] += blk_av[igofr + k] / blk_norm;
        glob_av2[igofr + k] += (blk_av[igofr + k] / blk_norm) * (blk_av[igofr + k] / blk_norm);
    }

    if (iblk == nblk) {
        for (int k = 0;k < nbins;k++) {
            err_gdir = Error(glob_av[igofr + k], glob_av2[igofr + k], iblk);
            Gave << k * bin_size << setw(wd) << glob_av[igofr + k] / iblk << setw(wd) << err_gdir << endl;
        }

    }

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void OldFinal(void) { //Write final configuration
    ofstream WriteOld;

    cout << "Print old configuration to file old.final " << endl << endl;
    WriteOld.open("old.final");

    for (int i = 0; i < npart; ++i) {
        WriteOld << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
    }
    WriteOld.close();
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if (iblk == 1) return 0.0;
    else return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
