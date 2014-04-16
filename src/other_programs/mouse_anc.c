
#include "mousethreestate.h"

#include "version.h"

char * program_name = "mouse_anc";


void show_help(FILE * ofp)
{

  show_help_DPRunImpl(ofp);
  show_standard_options(ofp);
  exit(63);
}

void show_version(FILE * ofp)
{
  fprintf(ofp,"%s\nVersion: %s\nReleased: %s\nCompiled: %s\n",program_name,VERSION_NUMBER,RELEASE_DAY,COMPILE_DATE);
  fprintf(ofp,"\nThis program is freely distributed under a Gnu Public License\n");
  fprintf(ofp,"The source code is copyright (c) GRL 1998 and others\n");
  fprintf(ofp,"There is no warranty, implied or otherwise on the performance of this program\n");
  fprintf(ofp,"For more information read the GNULICENSE file in the distribution\n\n");
  fprintf(ofp,"Credits: Ewan Birney <birney@sanger.ac.uk> wrote the core code.\n");
  exit(63);   
}


int main(int argc,char ** argv)
{
  GenoVarSet * s;
  GenoVarSet * sclean;
  PackAln * pal;
  AlnBlock * alb;
  
  DPRunImpl * dpri;
  int i;

  GenomePara * p;
  FrameSet * fs;

  FILE * ifp;
  FILE * ofp;


  AncestralVarSet * avs;
  AncestralBlockSet * abs;
  
  dpri =  new_DPRunImpl_from_argv(&argc,argv);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  
  ifp = openfile(argv[1],"r");
  assert(ifp != NULL);
  
  s = read_sanger_genotype_file(ifp);

  sclean = only_simple_snp_loci_GenoVarSet(s);
 
  p = new_GenomePara(0.99,0.01);

  avs = ancestral_prediction_MouseThreeState(sclean,"PWK","CAST","WSB",p,dpri);


  fs =  AncestralVarChr_as_FrameSet(avs,avs->chr[0]);

  abs = AncestralBlockSet_from_AncestralVarSet(avs,sclean);

  ofp = openfile("ancestor_single.txt","w");
  dump_AncestralBlockSet(abs,ofp);
  
  ofp = openfile("ancestor.ps","w");
  flat_no_frame_postscript_FrameSet(fs,ofp);
  
  ofp = openfile("ancestor.txt","w");
  write_simple_AncestralVarSet(avs,ofp);

  ofp = openfile("ancestor_pairwise.stat","w");
  dump_as_pairwise_AncestralVarSet(avs,sclean,ofp);

}
