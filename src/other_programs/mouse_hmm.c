
#include "mousethreestate.h"

#include "version.h"

char * program_name = "mouse_hmm";


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


#define MouseSNPMatch_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*3)+STATE][i+0]   
#define MouseSNPMatch_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]  
#define MouseSNPMatch_READ_OFF_ERROR -2
 



int main(int argc,char ** argv)
{
  GenoVarSet * s;
  PackAln * pal;
  AlnBlock * alb;
  
  DPRunImpl * dpri;
  SnpMatch * snpm;
  int i;

  GenomePara * p;

  MouseSNPMatch * forward;
  MouseSNPMatch * backward;
  

  FILE * ifp;

  
  dpri =  new_DPRunImpl_from_argv(&argc,argv);

  strip_out_standard_options(&argc,argv,show_help,show_version);
  
  ifp = openfile(argv[1],"r");
  assert(ifp != NULL);
  
  s = read_sanger_genotype_file(ifp);


  
  p = new_GenomePara(0.95,0.01);


  fprintf(stdout,"Score for match %d, score for mismatch %d\n",p->match,p->mismatch);

  for(i=0;i < s->len ;i++) {
    snpm = new_SnpMatch(s->chr[i],s,"DBA","PWK","CAST","WSB");

    /*    show_SnpMatchStats(snpm,stdout); */


    backward = backward_logsum_MouseSNPMatch(p,snpm,dpri);
    forward = forward_logsum_MouseSNPMatch(p,snpm,dpri); 


    fprintf(stdout,"Backward score is %d\n",
	    MouseSNPMatch_EXPL_SPECIAL(backward,0,-1,0));

    fprintf(stdout,"Forward score is %d Backward score is %d\n",MouseSNPMatch_EXPL_SPECIAL(forward,0,forward->lenj-1,1),
	    MouseSNPMatch_EXPL_SPECIAL(backward,0,-1,0));
    

  }
  
  
  
}
