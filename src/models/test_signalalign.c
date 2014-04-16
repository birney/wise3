#include "signalalign.h"


int main(int argc,char ** argv)
{
  SignalMap * map;
  int i;
  FILE * sm;
  Sequence * seq;
  SignalSeq * ss;
  FILE * sseqf;

  PackAln * pal;
  AlnBlock * alb;
  DPRunImpl * dpri;


  srandomdev();

  sm = openfile(argv[1],"r");

  dpri = new_DPRunImpl_from_argv(&argc,argv);


  map = read_SignalMap(sm);
  seq = read_fasta_file_Sequence(argv[2]);
  sseqf = openfile(argv[3],"r");

  ss = read_SignalSeq(sseqf);


  assert(map != NULL);
  assert(seq != NULL);
  assert(ss != NULL);

  pal = PackAln_bestmemory_SimpleSignalMat(ss,seq,map,-10,-2,NULL,dpri);

  alb = convert_PackAln_to_AlnBlock_SimpleSignalMat(pal);

  show_alignment_with_fit_SimpleSignalMat(alb,ss,seq,map,stdout);


}
