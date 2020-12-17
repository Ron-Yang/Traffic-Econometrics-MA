
//######################################
/**
 parameter-linear deterministic utility function V for alternative i for a given person
 - NI=number of alternatives
 - NC=number of characteristics (e.g. NC=2 if times and costs)
 - NS=number of socioeconomic+external vars (e.g.; NS=2 if gender and weather)
 - Mparam=number of model parameters
 - Cdata=characteristica in form C_00, C_01,...C_ji ... C_{NC-1,NI-1}
   e.g., C[i]=Traveltime(i), C[1*NI+i]=Cost(i)
 - Sdata=Socioeconomic vars, 
   e.g. S[0]=gender (0=male,1=female), S[1]=age, ...
 - beta = parameters
*/
//######################################

//#############################################################
// 4 alternatives, Alt0=ped, Alt1=bike, Alt2=pTransp, Alt3=car
//#############################################################

 int NC=0; // 
 int NS=1; // # only distance, must be treated as sociodemographic var 
 int NI=4; // # of alternatives
 int Mparam=6; // # of parameters

 double Vfunc(int i, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
   double V=
     beta[0]*delta(i,0)
     +beta[1]*delta(i,1)
     +beta[2]*delta(i,2)
     +beta[3]*SdataPers[0]*delta(i,0) 
     +beta[4]*SdataPers[0]*delta(i,1) //!!! don't forget to remove semicolon!!
     +beta[5]*SdataPers[0]*delta(i,2); //!!! don't forget to remove semicolon!!

return V;
}
