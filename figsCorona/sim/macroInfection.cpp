
/*
  Template fuer main() Programm; Ort:
  
  ~/versionedProjects/lib/templates/templateMain.cpp

  makefile dazu:
  
  ~/versionedProjects/lib/templates/makefile
 
  Achtung! Auch ohne .h File muss man bei $OBJECTS immer auch das 
  File mit der Main-Methode dazunehmen!
  (sonst "ld: undefined reference to main"
*/


// c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//alternatively there is <cstdio> which declares
//everything in namespace std
//but the explicit "using namespace std;" puts
//everything in global namespace


// c++ 
#include <iostream>
using namespace std;

// own
#include "general.h"

#include "Statistics.h"
#include "RandomUtils.h" // contains, e.g.,  myRand()
#include "InOut.h"


// constants

static const int NDATA_MAX=5000;// max. number of data points
static const int MAXSTR=500;// max. number of data points


//#####################################################
//#####################################################
int main(int argc, char* argv[]) {


  //#####################################################
  // SI model
  //#####################################################

  double beta=0.1;  // infection rate #infected per day if everybody is S

  double initI=0.001; //initial infection percentage
  double dt=0.1;   // internal time step
  int tmax=200;


  // simulation

  int ndt=int(1./dt+0.5);
  int nData=tmax;
  double initS=1.-initI;  // initial percentage of susceptibles 
  double arrI[nData];
  double arrS[nData];
  double arrTime[nData];
  arrS[0]=initS;
  arrI[0]=initI;
  arrTime[0]=0;
  double S=initS;
  double I=initI;

  for(int i=1; i<ndt*tmax; i++){
    double Iold=I;
    double Sold=S;

    double Ipred=Iold+dt*beta*Iold*Sold;
    double Spred=1-Ipred;
    I+=dt*beta*0.5*(Iold*Sold+Ipred*Spred);
    S=1-I;

    if(i%ndt==0){
      arrTime[i/ndt]=i/ndt;
      arrS[i/ndt]=S;
      arrI[i/ndt]=I;
      cout <<"t="<<i*dt<<" I="<<I<<endl;
    }
  }

  // output

  InOut inout;
  char   fnameOut[MAXSTR];
  char   titleStr[MAXSTR];

  sprintf(fnameOut,"%s","modelSI.dat");
  sprintf(titleStr,"#SI model beta=%2.4f  produced by macroInfection[.cpp]\n#time[days]\tS\t\tI", beta);

  cout <<" writing to "<<fnameOut<<" ..."<<endl;
  inout.write_array(fnameOut, nData, arrTime, arrS, arrI, titleStr);


  //#####################################################
  // SIR model
  //#####################################################

  beta=0.2;          // infection rate #infected per day if everybody is S
  double gamma=0.1;  // inverse lifetime of the "I" status (contagious)


  double arrR[nData];  // percentage of recovered/removed
  double R;
  S=arrS[0]=initS;
  I=arrI[0]=initI;
  R=arrR[0]=0;
  arrTime[0]=0;

  for(int i=1; i<ndt*tmax; i++){
    double Iold=I;
    double Sold=S;

    double Spred=Sold - dt*beta*Iold*Sold;
    double Ipred=Iold + dt*beta*Iold*Sold - dt*gamma*Iold;

    S+= -dt*beta*0.5*(Iold*Sold+Ipred*Spred);
    I+=  dt*beta*0.5*(Iold*Sold+Ipred*Spred) - dt*gamma*0.5*(Iold+Ipred);
    R=1-S-I;

    if(i%ndt==0){
      arrTime[i/ndt]=i/ndt;
      arrS[i/ndt]=S;
      arrI[i/ndt]=I;
      arrR[i/ndt]=R;
    }
  }

  // output

  sprintf(fnameOut,"%s","modelSIR.dat");
  sprintf(titleStr,"#SIR model beta=%2.4f gamma=%2.4f  produced by macroInfection[.cpp]\n#time[days]\tS\tI\tR", beta, gamma);

  cout <<" writing to "<<fnameOut<<" ..."<<endl;
  inout.write_array(fnameOut, nData, arrTime, arrS, arrI, arrR, titleStr);


  //#####################################################
  // SEIR model susceptible-exposed-infectious-ecovered/removed
  //#####################################################

  // S ->beta-> E ->alpha-> I ->gamma->R
  // R0=beta/gamma

  beta=0.2;          // infection rate #infected per day if everybody is S
  gamma=0.1;         // inverse lifetime of the "I" status (contagious)
  double alpha=0.2;  // inverse incubation time if exponentially distributed

  double arrE[nData];  // percentage of recovered/removed
  double E;
  S=arrS[0]=initS;
  E=arrE[0]=1-initS;
  I=arrI[0]=0;
  R=arrR[0]=0;

  for(int i=1; i<ndt*tmax; i++){
    double Sold=S;
    double Eold=E;
    double Iold=I;

    double Spred=Sold - dt*beta*Iold*Sold;
    double Epred=Eold + dt*beta*Iold*Sold - dt*alpha*Eold;
    double Ipred=Iold + dt*alpha*Eold     - dt*gamma*Iold;

    S+= -dt*beta*0.5*(Iold*Sold+Ipred*Spred);
    E+=  dt*beta*0.5*(Iold*Sold+Ipred*Spred) - dt*alpha*0.5*(Eold+Epred);
    I+=  dt*alpha*0.5*(Eold+Epred) - dt*gamma*0.5*(Iold+Ipred);
    R=1-S-E-I;

    if(i%ndt==0){
      arrTime[i/ndt]=i/ndt;
      arrS[i/ndt]=S;
      arrE[i/ndt]=E;
      arrI[i/ndt]=I;
      arrR[i/ndt]=R;
    }
  }

  // output

  sprintf(fnameOut,"%s","modelSEIR.dat");
  sprintf(titleStr,"#SEIR model beta=%2.4f gamma=%2.4f alpha=%2.4f  produced by macroInfection[.cpp]\n#time[days]\tS\tE\tI\tR", beta, gamma, alpha);

  cout <<" writing to "<<fnameOut<<" ..."<<endl;
  inout.write_array(fnameOut, nData, arrTime, arrS, arrE, arrI, arrR, 
		    titleStr);


 return(0);
}

// #########################################################







