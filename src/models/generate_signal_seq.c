
#include "dualsignal.h"



int main(int argc,char ** argv)
{
  SignalMap * map;
  int i;
  FILE * sm;
  Sequence * seq;
  SignalSeq * ss;
  SignalSim * sim;

  sim = new_SignalSim(0.0000000001,0.00000000001,0.000000000001);

  srandomdev();

  sm = openfile(argv[1],"r");


  map = read_SignalMap(sm);
  seq = read_fasta_file_Sequence(argv[2]);

  assert(map != NULL);
  assert(seq != NULL);
  assert(sim != NULL);

  for(i=0;i < 10 ;i++) {
    ss = simulated_SignalSeq_from_Seq(seq,map,sim);
    write_SignalSeq(ss,stdout);
  }

}
