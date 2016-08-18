#include<stdio.h>

class unnormalized_distribution
{
 private:
     long int size;
     double *tree;
  public:
     unnormalized_distribution(int p=0);
     ~unnormalized_distribution();
     void set(long int i,double d);
     void zero();
     double get(long int i);
     long int select();
     double total();
};

class SPM
{
 public:
  enum {Setting_up, Running, Stopped, Paused, Quitting} state;
  long int coagulation_events, swap_events, Phi_events, M_events, sample_events, source_events;    
  double Phi;	
  double M;
  long int *z;
  long int x;
  double nt;
  unnormalized_distribution *x_to_the_minus_two_thirds;
  //	unnormalized_distribution *reciprocal_of_theta;
  unnormalized_distribution *x_in_distribution;
  int *slider_xin, *slider_theta;
  double theta_distribution[120];
  double Q;
  long int N;
  double tau;
  
  

  //void recompute_gaussian_xin();
  //void recompute_linear_xin();
  //void recompute_gaussian_theta();
  // void recompute_linear_theta();
  double theta_of(long int);
  void jump_coagulation();
  void jump_M();
  void jump_Phi();
  void jump_swap();
  void jump_sampling();
  void jump_source();
  void transition();
  
};

