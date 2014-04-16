
#include "snplocus.h"
#include <stdio.h>



int main(int argc,char ** argv) 
{
  GenoVarSet * s;

  s = read_sanger_genotype_file(stdin);
  write_sanger_GenoVarSet(s,stdout);


  return(0);

}
