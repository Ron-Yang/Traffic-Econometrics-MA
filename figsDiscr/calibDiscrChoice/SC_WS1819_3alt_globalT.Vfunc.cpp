
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
// 3 alternatives, Alt0=ped, alt1=bike, alt2=publicTransport/car
//#############################################################

 int NC=2; // # characteristica in .inputData (times and costs)
 int NS=1; // # sociodemographic or external vars in .inputData (weather)
 int NI=3; // # of alternatives
 int Mparam=4; // # of parameters

 double Vfunc(int i, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  if(false){
    //if(CdataPers[NI+i]-CdataPers[NI]!=0){
    cout <<"i="<<i<<" Cost_i-Cost_0="<<CdataPers[NI+i]-CdataPers[NI]<<endl;
  }
  double V=
    +beta[0]*delta(i,0)        //!!! don't forget to remove semicolon!!
    +beta[1]*delta(i,1)        // start with AC to have same params for extensions!
    +beta[2]*CdataPers[NI+i]   // cost (only alt 2)
    +beta[3]*CdataPers[i];     // time sensitivity

 return V;
}
