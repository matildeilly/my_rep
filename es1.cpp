/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <vector>

using namespace std;

struct position {
    double x, y, z;

    position() : x(0), y(0), z(0) {}

    position(double x1, double x2, double x3) {
        x = x1;
        y = x2;
        z = x3;
    }

    double dist() {
        return sqrt(x * x + y * y + z * z);
    }
};

double min (double x, double y){
      if (x<y) {return x;}
      else {return y;}
}

/*int*/void avanza_unif_gs (vector<position> &V, int index, Random* rnd, double larghezza){
   position old = V[index - 1]; 

   position nuovo;
   nuovo.x = rnd->Rannyu(old.x - larghezza/2., old.x + larghezza/2.);
   nuovo.y = rnd->Rannyu(old.y - larghezza/2., old.y + larghezza/2.); 
   nuovo.z = rnd->Rannyu(old.z - larghezza/2., old.z + larghezza/2.); 

   double A = min(1., exp(2*(old.dist() - nuovo.dist())));

   double soglia = rnd->Rannyu();
   if (soglia <= A) {V[index] = nuovo;
   /*return 1;*/}
   else {V[index] = old;
   /*return 0;*/}   
}

// Parti commentate con return 0,1 sono per vedere che A sia in media 0.5

/*int*/void avanza_unif_2p (vector<position> &V, int index, Random* rnd, double larghezza){
   position old = V[index - 1]; 

   position nuovo;
   nuovo.x = rnd->Rannyu(old.x - larghezza/2., old.x + larghezza/2.); 
   nuovo.y = rnd->Rannyu(old.y - larghezza/2., old.y + larghezza/2.); 
   nuovo.z = rnd->Rannyu(old.z - larghezza/2., old.z + larghezza/2.); 

   double x1 = nuovo.dist();
   double x2 = old.dist();

   double cos_theta1 = nuovo.z/x1; // (Angolo polare in sferiche)
   double cos_theta2 = old.z/x2;

   double A = min(1., exp(x2 - x1) * x1*x1/(x2*x2) * cos_theta1*cos_theta1/(cos_theta2*cos_theta2));

   double soglia = rnd->Rannyu();
   if (soglia <= A) {V[index] = nuovo;
      /*return 1;*/}
   else {V[index] = old;
      /*return 0;*/}
}

/*int*/void avanza_gauss_gs (vector<position> &V, int index, Random* rnd, double sigma){
   position old = V[index - 1]; 

   position nuovo;
   nuovo.x = rnd->Gauss(old.x,sigma); 
   nuovo.y = rnd->Gauss(old.y,sigma);
   nuovo.z = rnd->Gauss(old.z,sigma);

   double A = min(1., exp(2*(old.dist() - nuovo.dist())));

   double soglia = rnd->Rannyu();
   if (soglia <= A) {V[index] = nuovo;
   /*return 1;*/}
   else {V[index] = old;
   /*return 0;*/}
   
}

/*int*/void avanza_gauss_2p (vector<position> &V, int index, Random* rnd, double sigma){
   position old = V[index - 1]; 

   position nuovo;
   nuovo.x = rnd->Gauss(old.x,sigma); 
   nuovo.y = rnd->Gauss(old.y,sigma);
   nuovo.z = rnd->Gauss(old.z,sigma);

   double x1 = nuovo.dist();
   double x2 = old.dist();

   double cos_theta1 = nuovo.z/x1;
   double cos_theta2 = old.z/x2;

   double A = min(1., exp(x2 - x1) * x1*x1/(x2*x2) * cos_theta1*cos_theta1/(cos_theta2*cos_theta2));

   double soglia = rnd->Rannyu();
   if (soglia <= A) {V[index] = nuovo;
      /*return 1;*/}
   else {V[index] = old;
      /*return 0;*/}
}
 
double Vector_mean (vector<double> V, int dim){
   double sum=0;
   for (int i=0; i<dim; i++){
      sum = sum + V[i];
   }
   return sum/dim;
}

double Vector_variance (vector<double> V, int dim){
   if (dim==1){return 0;}
   else{
   double sum=0;
   double sum2=0; 
   for (int i=0; i<dim; i++){
      sum = sum + V[i];
      sum2 = sum2 + pow(V[i],2);
   }

   return sqrt(1.0/dim * (1.0/dim*sum2 - pow(1.0/dim*sum,2)));

}
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();

   int M = pow(10,7); //numero di tiri
   int N = 50; //numero blocchi
   int L = M/N; //tiri per blocco

   double lunghezza_unif_gs = 2.4;
   double lunghezza_unif_2p = 5.9;
   double sigma_gauss_gs = 0.75;
   double sigma_gauss_2p = 1.85;

   vector<double> D_unif_gs(N); // Contiene le N distanze
   vector<position> V_unif_gs(M); // Contiene le M posizioni nel caso g.s. con distr uniforme
   double sum_unif_gs =0.;

   vector<double> D_unif_2p(N);
   vector<position> V_unif_2p(M);
   double sum_unif_2p =0.;

   vector<double> D_gauss_gs(N);
   vector<position> V_gauss_gs(M);
   double sum_gauss_gs =0.;

   vector<double> D_gauss_2p(N);
   vector<position> V_gauss_2p(M);
   double sum_gauss_2p =0.;

   int pos_in_vec =0;
// int test=0;

   ofstream out("metro_unif_gs_out.txt");
   ofstream pout("metro_unif_2p_out.txt");
   ofstream gout("metro_gauss_gs_out.txt");
   ofstream gpout("metro_gauss_2p_out.txt");

   for (int i=0; i<M; i++){
      
     if (i==0){
         V_unif_gs[i].x=1.5;
         V_unif_gs[i].y=0;
         V_unif_gs[i].z=0;
         
         V_unif_2p[i].x=5.;
         V_unif_2p[i].y=0;
         V_unif_2p[i].z=0;

         V_gauss_gs[i].x=1.5;
         V_gauss_gs[i].y=0;
         V_gauss_gs[i].z=0;

         V_gauss_2p[i].x=5.;
         V_gauss_2p[i].y=0;
         V_gauss_2p[i].z=0;
      }

      if (i!=0){
      /*test +=*/ avanza_unif_gs(V_unif_gs,i,&rnd,lunghezza_unif_gs);
      /*test +=*/ avanza_unif_2p(V_unif_2p,i,&rnd,lunghezza_unif_2p);
      /*test +=*/ avanza_gauss_gs(V_gauss_gs,i,&rnd,sigma_gauss_gs);
      /*test +=*/ avanza_gauss_2p(V_gauss_2p,i,&rnd,sigma_gauss_2p);
      }
      sum_unif_gs += V_unif_gs[i].dist();
      sum_unif_2p += V_unif_2p[i].dist();
      sum_gauss_gs += V_gauss_gs[i].dist();
      sum_gauss_2p += V_gauss_2p[i].dist();


      if ((i+1)%L==0){ // i.e. ho finito un blocco
         D_unif_gs[pos_in_vec] = sum_unif_gs/L;
         D_unif_2p[pos_in_vec] = sum_unif_2p/L;
         D_gauss_gs[pos_in_vec] = sum_gauss_gs/L;
         D_gauss_2p[pos_in_vec] = sum_gauss_2p/L;

         out<<pos_in_vec + 1<<"\t"<<Vector_mean(D_unif_gs,pos_in_vec + 1)<<"\t"<<Vector_variance(D_unif_gs,pos_in_vec + 1)<<endl;          
         pout<<pos_in_vec + 1<<"\t"<<Vector_mean(D_unif_2p,pos_in_vec + 1)<<"\t"<<Vector_variance(D_unif_2p,pos_in_vec + 1)<<endl; 
         gout<<pos_in_vec + 1<<"\t"<<Vector_mean(D_gauss_gs,pos_in_vec + 1)<<"\t"<<Vector_variance(D_gauss_gs,pos_in_vec + 1)<<endl; 
         gpout<<pos_in_vec + 1<<"\t"<<Vector_mean(D_gauss_2p,pos_in_vec + 1)<<"\t"<<Vector_variance(D_gauss_2p,pos_in_vec + 1)<<endl; 
                 
         sum_unif_gs =0.;
         sum_unif_2p =0.;
         sum_gauss_gs =0.;
         sum_gauss_2p =0.;

         pos_in_vec ++;

         /*cout<<double(test)/L<<endl;
         test = 0;*/
        }
   }


      out.close();
      pout.close();
      gout.close();
      gpout.close();

      //SVUTO E COMPRIMO I VECTORS
      V_unif_gs.clear();
      V_unif_gs.shrink_to_fit();
      D_unif_gs.clear();
      D_unif_gs.shrink_to_fit();

      V_unif_2p.clear();
      V_unif_2p.shrink_to_fit();
      D_unif_2p.clear();
      D_unif_2p.shrink_to_fit();

      V_gauss_gs.clear();
      V_gauss_gs.shrink_to_fit();
      D_gauss_gs.clear();
      D_gauss_gs.shrink_to_fit();

      V_gauss_2p.clear();
      V_gauss_2p.shrink_to_fit();
      D_gauss_2p.clear();
      D_gauss_2p.shrink_to_fit();

      return 0;
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
