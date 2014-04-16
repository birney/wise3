
#include "ibd_model.h"

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
  PairMatch * pm;

  DPRunImpl * dpri;

  int i;
  int j;

  GenomePairPara * p;

  FILE * ifp;
  FILE * ofp;



  dpri =  new_DPRunImpl_from_argv(&argc,argv);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  
  ifp = openfile(argv[1],"r");
  assert(ifp != NULL);
  
  s = read_sanger_genotype_file(ifp);

  sclean = only_simple_snp_loci_GenoVarSet(s);
 
  p = new_GenomePairPara_null(0.001,0.1,0.01,0.01);



  /*
  for(i=0;i<pm->len;i++) {
    fprintf(stdout,"Got %d at %d\n",pm->match[i],i);
  }

  */

  for(i=0;i<s->ind_len;i++) {
    for(j=i+1;j<s->ind_len;j++) {

      pm = new_PairMatch(sclean,sclean->chr[0],s->ind[i]->id,s->ind[j]->id);
      pal = PackAln_bestmemory_MousePairMatch(p,pm,NULL,dpri);

      alb = convert_PackAln_to_AlnBlock_MousePairMatch(pal);

      /*dump_ascii_AlnBlock(alb,stdout);*/
      dump_naive_AlnBlock_IBD(s->ind[i]->id,s->ind[j]->id,sclean->chr[0],alb,stdout);
      
      free_PairMatch(pm);
      free_AlnBlock(alb);
      free_PackAln(pal);
    }
  }

}
