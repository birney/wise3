#include "probability.h"

int main(int argc,char ** argv)
{
  Probability p;
  double zscore;

  zscore = strtold(argv[1],NULL);

  p = crude_normal_probability(5,1,zscore);

  printf("From %f got %f\n",zscore,p);
  
}

