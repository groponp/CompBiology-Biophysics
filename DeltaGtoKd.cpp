/* This script compute Dissociation constant from Gibbs binding energy
Writte by: Ropón-Palacios G.
E-mail: georopon@gmail.com 
Filiation: 
Computational biophysicist, Research
Condensed matter physics Lab, Departament of physics 
Universidade Federal de Alfenas, Minas Gerais, Brasil. 
date: 30 Jan 2020. 
usage: ./DeltaGtoKd -11.5 nM 
*/
#include<iostream>
#include<stdio.h> //
#include<stdlib.h> // pass char to float 
#include<math.h> // math date library
#include<cstdlib> // pass char to int date

using namespace std; 

int main(int argc, char *argv[]){
        /*
        argc: get total numbers of arguments
        argv: array with the command pass in command line
        */       
        if(argc >= 3){   
                // Variable definition
                float delg = atof(argv[1]); 
                //string conc  = argv[2]; 
                double kd_nM = 0, kd_mM = 0, kd_M = 0; 
                /*kd_nM = exp((delg*1000)/(1.98*298.15))*1000000000;
                cout<< kd_nM<<endl;*/
                // Math processing using conditional syntaxis 
               
                if (std::string(argv[2]) == "M"){
                        kd_M = exp((delg*1000)/(1.98*298.15)); 
                        cout<< kd_M<<endl; 
                }
                else if (std::string(argv[2]) == "mM"){
                        kd_mM = exp((delg*1000)/(1.98*298.15))*1000000; 
                        cout<< kd_mM<<endl;
                }
                else if (std::string(argv[2]) == "nM"){
                        kd_nM = exp((delg*1000)/(1.98*298.15))*1000000000;
                        cout<< kd_nM<<endl;
                }       
                else{
                        cout<< "process not sopported"<<endl;
                }             
        }
        else{
         cout<<"no command pass"; 
        }                               
        return 0;  
}
