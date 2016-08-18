//g++ -O3 -Wall spm.c++ -o spm



#include <openssl/rand.h>
#include <iostream>
#include <cmath>
#include "spm.h"

#define K_INF 0.1


long int Rand(long int N) //Return an integer in the range {0,...,N-1}
{
  return (long int)(N*(double)rand()/(1+(double)RAND_MAX));
};

void SPM::jump_coagulation()
{
  long int y=x+z[x_to_the_minus_two_thirds->select()];
  if (y>0) {x=y;} else printf("TOO LARGE!\n");
  coagulation_events++;
};

void SPM::jump_M()
{
  M=theta_of(x)*Q;//printf("Th(%f)=%f\n",(double)x,theta_of(x));
  M_events++;
};

void SPM::jump_Phi()
{
  Phi+=(M-Phi)/N; // ajs224 spinBox->value();
  Phi_events++;
};

void SPM::jump_swap()
{
  long int i=Rand(N); // ajs224
  x_to_the_minus_two_thirds->set(i,pow(x,0.0));
  z[i]=x; //printf("%f\n",(double)x);
  swap_events++;
};

void SPM::jump_sampling()
{
  nt+=nt/((double) N);
  sample_events++;
};

void SPM::jump_source()
{
  x=x_in_distribution->select()+1;
  source_events++;
};

void SPM::transition()
{
  double coag_rate, M_rate, Phi_rate, swap_rate, sample_rate, source_rate, rho;
  
  coag_rate = K_INF*pow(x,1.0)*Phi*x_to_the_minus_two_thirds->total()/N;
  sample_rate = ((double) N)/nt;
  M_rate = 1/theta_of(x);
  Phi_rate = ((double)N)/((double)tau*nt);
  source_rate = Q/M;
  swap_rate = ((double)N)*M/Phi/((double)tau*nt);
  rho = coag_rate+M_rate+Phi_rate+swap_rate+sample_rate+source_rate;
  double d=rho*((double)rand())/(1+(double)RAND_MAX)-coag_rate;
  
  //printf("(C=%f,M=%f,Phi=%f,Sw=%f,Sa=%f,Sc=%f,Q=%f)\n",coag_rate,M_rate,Phi_rate,swap_rate,sample_rate,source_rate,Q);
  
  if (d<0) {jump_coagulation();}
  else
    {
	  d-=M_rate;
      if (d<0){jump_M();}
      else
		{
		  d-=Phi_rate;
		  if (d<0) {jump_Phi();}
		  else
			{
			  d-=swap_rate;
			  if (d<0) {jump_swap();}
			  else
				{
				  d-=sample_rate;
				  if (d<0) {jump_sampling();}
				  else
					jump_source();
				}
			}
		}
    } 
  //printf("M=%f, Phi=%f, x=%f\n",M,Phi,(double)x);
};

unnormalized_distribution::unnormalized_distribution(int p){
  tree=new double[size=1<<(p+1)] ;
  for (long int i=0; i<size; i++)
    {
      tree[i]=0;
    }
};

unnormalized_distribution::~unnormalized_distribution()
{
  delete [] tree;
};

void unnormalized_distribution::set(long int i,double d)
{
  tree[i+=size>>1]=d; 
  for(i>>=1;i!=0;i>>=1) tree[i]=tree[i<<1]+tree[(i<<1)|1];
};

void unnormalized_distribution::zero()
{
  for(long int i=0; i<size; i++)
    tree[i]=0;
};

double unnormalized_distribution::total()
{
  return tree[1];
};

double unnormalized_distribution::get(long int i)
{
  return i<(size>>2)? tree[i+(size>>1)]:0;
};

long int unnormalized_distribution::select()
{
  long int i=1;
  double d=((double) rand())/((double)RAND_MAX)*tree[1];
  for(;i< (size>>1);)
    if (d>tree[i<<=1])
      d-=tree[i++];
  return i-(size>>1);
};

double SPM::theta_of(long int x0)
{

  //return 0.0e0;
  
  //return 1e-1;
  return 1e1;
  
  /*
    if (x0<=120)
    return theta_distribution[x0-1];
  else
    return doubleSpinBox->value()*pow((double)x0/(1000/10.3+1),doubleSpinBox_2->value()); 
  */
};

int main()
{

  using namespace std;
  
  SPM sim;
  sim.N = 1024;
  sim.tau=1e0;
  sim.Q = 1e-1;
  sim.x_in_distribution = new unnormalized_distribution(7);

  for (long int i=0; i<sim.N; i++)
    {
	  sim.x_in_distribution->set(i,1.0);
	}
 
 
  sim.coagulation_events = sim.swap_events = sim.Phi_events = sim.M_events = sim.sample_events = sim.source_events = 0;
  sim.z=new long int[sim.N];
  
  sim.x=1;
  sim.Phi=0;
  int Np=sim.N-1;
  sim.nt=1.0;
  int p=0;
  while(Np!=0)
    {
      Np>>=1;p++;
    }
  
  sim.x_to_the_minus_two_thirds = new unnormalized_distribution(p);

  for (int j=0; j<1000;j++)
    {
	  
	  
	  for (long int i=0; i<sim.N; i++)
		{
		  sim.z[i]=1+sim.x_in_distribution->select();
		  sim.x_to_the_minus_two_thirds->set(i,pow(sim.z[i],0.0));
		  sim.Phi += sim.Q*sim.theta_of(sim.z[i]) / sim.N;
		  //cout << sim.x << "\t" << sim.Phi<< "\t" << sim.M << "\t" << sim.z[i] <<endl;

		  sim.transition();
		}
	  sim.M=sim.Phi;


	  

	  
	  // calculate moments
	  double m0,m1,m2,m3;
	  m0=0.0;
	  m1=0.0;
	  m2=0.0;
	  m3=0.0;
  
	  for (long int i=0; i<sim.N; i++)
		{
		  m0 += 1.0;
		  m1 += sim.z[i];
		  m2 += sim.z[i]*sim.z[i];
		  m3 += sim.z[i]*sim.z[i]*sim.z[i];
		}
  
	  m0 /= sim.N;
	  m1 /= sim.N;
	  m2 /= sim.N;
	  m3 /= sim.N;
	  cout << m0 << "\t"
		   << m1 << "\t"
		   << m2 << "\t"
		   << m3 << "\t"
		   << sim.x << "\t"
		   << sim.Phi<< "\t"
		   << sim.M << endl;
  
	  
    


	  
	}
  
  cout << sim.coagulation_events << "\t"
	   << sim.swap_events << "\t"
	   << sim.Phi_events << "\t"
	   << sim.M_events << "\t"
	   << sim.sample_events << "\t"
	   << sim.source_events << endl;
  
  return 0; 
}

