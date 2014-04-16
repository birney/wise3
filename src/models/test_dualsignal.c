
#include "dualsignal.h"



int main(int argc,char ** argv)
{
  SignalMap * map;
  int i;
  FILE * sm;
  Sequence * seq;
  SignalSeq * ss;
  SignalSim * sim;

  sm = openfile(argv[1],"r");

  map = read_SignalMap_tsv(5,sm);

  fprintf(stderr,"Read map\n");

  convert_normal_to_absolute_SignalMap(map,500);

  prepare_SignalMap(map,1.0);


  seq = read_fasta_Sequence(stdin);

  sim = new_SignalSim(0,0,0);

  ss = simulated_SignalSeq_from_Seq(seq,map,sim);


  write_SignalSeq(ss,stdout);

}
