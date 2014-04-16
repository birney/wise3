#include "dualsignal.h"
#include "signalalign.h"
#include "version.h"




char * program_name = "porewise";


void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) EMBL and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@ebi.ac.uk>\n");
  exit(63);   
}


void show_help(FILE * ofp)
{
  fprintf(ofp,"%s signal-file dna-fasta\n",program_name);

  fprintf(ofp,"  -sigrev    signal is reverse strand\n");
  fprintf(ofp,"Model Parameter options\n");
  fprintf(ofp,"  -signalmap signal map tab delimited file\n");
  fprintf(ofp,"  -weight    weight of counts from normal distributions in SignalMap [500.0]\n");
  fprintf(ofp,"  -pscount   pseudocount on weights [1.0]\n");
  fprintf(ofp,"  -seqdiff <prob> probability of a sequence difference\n");
  fprintf(ofp,"  -seqdiffext <prob> extension probability across differences\n");
  fprintf(ofp,"  -signaldiff <prob> open for signal-only regions\n");
  fprintf(ofp,"  -signalext  <prob> extension for signal-ext regions\n");
  fprintf(ofp,"Display options\n");
  fprintf(ofp,"  -pal       show raw alignment\n");
  fprintf(ofp,"  -alb       show text (alb) alignment\n");
  fprintf(ofp,"  -fit       show fit alignment\n");


  show_help_DPRunImpl(ofp);

  show_standard_options(ofp);
}

int main(int argc,char ** argv)
{
  DPRunImpl * dpri;

  double dis_weight  = 500.0;
  double pseudocount = 1.0;

  SignalMap * sm;
  int kmer = 5;

  SignalSeq * sseq;
  SignalEventList * sel;
  
  Sequence * comp;

  char * signal_map_file;
  FILE * ifp;
  FILE * sfp;

  PackAln  * pal;
  AlnBlock * alb;

  Probability signal_prob = 0.2;
  Probability signal_ext  = 0.2;

  Probability seqdiff = 0.05;
  Probability seqdiff_ext = 0.01;

  boolean show_pal = FALSE;
  boolean show_alb = FALSE;
  boolean show_fit = FALSE;
  
  boolean sig_reverse  = FALSE;

  boolean print_model = TRUE;
  
  strip_out_float_argument(&argc,argv,"weight",&dis_weight);
  strip_out_float_argument(&argc,argv,"pscount",&pseudocount);

  strip_out_float_argument(&argc,argv,"seqdiff",&seqdiff);
  strip_out_float_argument(&argc,argv,"seqdiffext",&seqdiff_ext);

  strip_out_float_argument(&argc,argv,"signaldiff",&signal_prob);
  strip_out_float_argument(&argc,argv,"signalext",&signal_ext);
  

  signal_map_file = strip_out_assigned_argument(&argc,argv,"signalmap");

  dpri = new_DPRunImpl_from_argv(&argc,argv);

  show_pal = strip_out_boolean_argument(&argc,argv,"pal");
  show_alb = strip_out_boolean_argument(&argc,argv,"alb");
  show_fit = strip_out_boolean_argument(&argc,argv,"fit");

  sig_reverse = strip_out_boolean_argument(&argc,argv,"sigrev");

  strip_out_standard_options(&argc,argv,show_help,show_version);
  if( argc != 3 ) {
    show_help(stdout);
    exit(12);
  }

  sfp = openfile(signal_map_file,"r");
  assert(sfp != NULL);

  sm = read_SignalMap_tsv(kmer,sfp);

  if( validate_SignalMap_with_warnings(sm) == FALSE ) {
    warn("Bad signal Map File. Cannot continue");
    exit(-2);
  }


  convert_normal_to_absolute_SignalMap(sm,dis_weight);

  prepare_SignalMap(sm,pseudocount);
  

  ifp = openfile(argv[1],"r");

  sel = read_SignalEventList(ifp);

  sseq = SignalSeq_from_SignalEventList(sel);

  
  assert(sseq != NULL);

  comp = read_fasta_file_Sequence(argv[2]);

  assert(comp);

  pal = PackAln_bestmemory_SimpleSignalMat(sseq,comp,sm,Probability2Score(signal_prob),Probability2Score(signal_ext),Probability2Score(seqdiff),Probability2Score(seqdiff_ext),NULL,dpri);

  alb = convert_PackAln_to_AlnBlock_SimpleSignalMat(pal);

  if( show_pal == TRUE ) {
    show_simple_PackAln(pal,stdout);
  }

  if( show_alb == TRUE ) {
    show_flat_AlnBlock(alb,stdout);
  }
  

  if( show_fit == TRUE ) {
    show_alignment_with_fit_SimpleSignalMat(alb,sel,comp,sm,stdout);
  }

}
