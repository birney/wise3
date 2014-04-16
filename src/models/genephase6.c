#ifdef _cplusplus
extern "C" {
#endif
#include "genephase6.h"


# line 289 "genephase6.dy"
AlnBlock * AlnBlock_from_phased_protein_wrap(Protein * pro,ThreeStateModel * tsm,Genomic * gen,CodonMapper * cm,RandomModel * rm,CompMat * mat,PhasedProteinPara * ppp,GeneParameter21 * gpara,DPRunImpl * dpri,PackAln ** palret)
{
  PackAln * pal;
  AlnBlock * alb;
  GeneParser21Score * gp21s;
  GeneParser4Score * gp4s;
  ComplexSequence * cs;
  ThreeStateModel * tsmin;

  GeneralGeneModelScore * ggms;
  
  GenePhaseModel * gpm;
  GenePhaseScore * gps;

  int i;
  
  assert(gen != NULL);
  assert(ppp != NULL);

  ggms = vanilla_GeneralGeneModelScore(gpara->ct,Bits2Probability(10),1.0,Bits2Probability(10));

  if( tsm == NULL ) {
    if( pro == NULL ) {
      warn("Major problem - neither a protein nor a tsm passed into phase model");
      return NULL;
    }
    tsmin = ThreeStateModel_alloc_len(pro->baseseq->len);
    
    if( pro->baseseq->name != NULL )
      tsmin->name = stringalloc(pro->baseseq->name);
    else tsmin->name = stringalloc("NoName");
  
    tsmin->rm = hard_link_RandomModel(rm);

    for(i=0;i<pro->baseseq->len;i++) { 
      add_ThreeStateModel(tsmin,ThreeStateUnit_from_half_bit_aminoacid(pro->baseseq->seq[i],mat,rm,ppp->gap,ppp->ext));
    }
  } else {
    tsmin = hard_link_ThreeStateModel(tsm);
  }


  gpm = GenePhaseModel_from_ThreeStateModel(tsmin,cm,rm,mat,ppp);

  GeneWise_fold_in_synchronised_RandomModel(gpm->gw,rm,gpara->cm,gpara->ct,0.5);

  gps = GenePhaseScore_from_GenePhaseModel(gpm);
  
  if( (gp21s = GeneParser21Score_from_GeneParser21(gpara->gp)) == NULL) {
    warn("Unable to make GeneParserScore model");
    return NULL;
  }
  
  gp4s = GeneParser4Score_from_GeneParser21Score(gp21s);
  
  if( (cs=evaluate_ComplexSequence_Genomic(gen,gpara->cses,0,Probability2Score(0.01))) == FALSE ) {
    warn("Unable to make ComplexSequence in TMS2DNA wrap");
  }
  
  pal = PackAln_bestmemory_GenePhase6(gps,cs,gp4s,ggms,NULL,dpri);

  if( palret != NULL ) {
    *palret = pal;
  }

  alb = convert_PackAln_to_AlnBlock_GenePhase6(pal);

  free_GeneralGeneModelScore(ggms);
  free_GenePhaseModel(gpm);
  free_GenePhaseScore(gps);
  free_ComplexSequence(cs);
  free_GeneParser21Score(gp21s);
  free_GeneParser4Score(gp4s);

  free_ThreeStateModel(tsm);

  return alb;

}

# line 87 "genephase6.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sun Mar 18 12:37:13 2007           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define MATCH 0  
#define INSERT 1 
#define DELETE 2 
#define INTRON_0 3   
#define INTRON_1 4   
#define INTRON_2 5   


#define LOOP 0   
#define START 1  
#define END 2    
#define BEFORE_MATCH_CODING 3    
#define BEFORE_MATCH_INTRON_0 4  
#define BEFORE_MATCH_INTRON_1 5  
#define BEFORE_MATCH_INTRON_2 6  
#define AFTER_MATCH_CODING 7 


#define GenePhase6_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+10)*6)+STATE][i+1] 
#define GenePhase6_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+10]    
#define GenePhase6_READ_OFF_ERROR -12
   


#define GenePhase6_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+11)%11][((i+1)*6)+STATE]  
#define GenePhase6_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+11)%11][STATE]   




#define GenePhase6_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]  
#define GenePhase6_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)  


/* Function:  PackAln_read_Shatter_GenePhase6(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_GenePhase6(GenePhase6 * mat) 
{
    GenePhase6_access_func_holder holder;    


    holder.access_main    = GenePhase6_shatter_access_main;  
    holder.access_special = GenePhase6_shatter_access_special;   
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_GenePhase6(mat,holder);  
}    


/* Function:  GenePhase6_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenePhase6_shatter_access_main(GenePhase6 * mat,int i,int j,int state) 
{
    return GenePhase6_SHATTER_MATRIX(mat,i,j,state); 
}    


/* Function:  GenePhase6_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenePhase6_shatter_access_special(GenePhase6 * mat,int i,int j,int state) 
{
    return GenePhase6_SHATTER_SPECIAL(mat,i,j,state);    
}    


/* Function:  calculate_shatter_GenePhase6(mat,dpenv)
 *
 * Descrip:    This function calculates the GenePhase6 matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [GenePhase6 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_GenePhase6(GenePhase6 * mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int should_calc;     
    int leni;    
    int lenj;    
    int tot; 
    int num; 
    int starti;  
    int startj;  
    int endi;    
    int endj;    


    int * SIG_0_0;   
    int * SIG_1_3;   
    int * SIG_1_6;   
    int * SIG_1_5;   
    int * SIG_1_4;   
    int * SIG_1_2;   
    int * SIG_1_1;   
    int * SIG_0_3;   
    int * SIG_0_6;   
    int * SIG_0_5;   
    int * SIG_0_4;   
    int * SIG_0_2;   
    int * SIG_0_1;   
    int * SIG_1_0;   
    int * SIG_0_8;   
    int * SIG_0_9;   
    int * SIG_0_10;  


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,6,lenj,8);    
    prepare_DPEnvelope(dpenv);   
    starti = dpenv->starti;  
    if( starti < 0 ) 
      starti = 0;    
    startj = dpenv->startj;  
    if( startj < 0 ) 
      startj = 0;    
    endi = dpenv->endi;  
    if( endi > mat->leni )   
      endi = mat->leni;  
    endj = dpenv->endj;  
    if( endj > mat->lenj )   
      endj = mat->lenj;  
    tot = (endi-starti) * (endj-startj); 
    num = 0; 


    start_reporting("GenePhase6 Matrix calculation: ");  
    for(j=startj;j<endj;j++) {  
      auto int score;    
      auto int temp;     
      for(i=starti;i<endi;i++)   {  
        /* Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed */ 
        should_calc = 0; 
        for(k=0;k<dpenv->len;k++)    {  
          auto DPUnit * u;   
          u = dpenv->dpu[k]; 
          switch(u->type)    {  
            case DPENV_RECT :    
              if( i >= u->starti && j >= u->startj && i < (u->starti+u->height) && j < (u->startj+u->length))    
                should_calc = 1;     
              break; 
            case DPENV_DIAG :    
              if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)   
                should_calc = 1;     
              break; 
            }  
          if( should_calc == 1 ) 
            break;   
          }  
        if( should_calc == 0)    
          continue;  


        SIG_0_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i,j);   
        SIG_1_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-3);   
        SIG_1_6 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-6);   
        SIG_1_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-5);   
        SIG_1_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-4);   
        SIG_1_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-2);   
        SIG_1_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-1);   
        SIG_0_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-3);   
        SIG_0_6 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-6);   
        SIG_0_5 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-5);   
        SIG_0_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-4);   
        SIG_0_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-2);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   
        SIG_0_8 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-8);   
        SIG_0_9 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-9);   
        SIG_0_10 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-10); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_3[MATCH] + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_3[INSERT] + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_3[DELETE] + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = GenePhase6_SHATTER_SPECIAL(mat,i-1,j-3,START) + ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_MATCH_CODING to state MATCH */ 
        temp = GenePhase6_SHATTER_SPECIAL(mat,i-1,j-3,BEFORE_MATCH_CODING) + (mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = SIG_1_6[INTRON_0] + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = SIG_1_5[INTRON_1] + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = SIG_1_4[INTRON_2] + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_2[MATCH] + mat->gp->transition[GP4_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_1[MATCH] + mat->gp->transition[GP4_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_4[MATCH] + mat->gp->transition[GP4_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_1_5[MATCH] + mat->gp->transition[GP4_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_SHATTER_SPECIAL(mat,i,j,LOOP) )    {  
          GenePhase6_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;   
          }  




        /* state MATCH is a source for special AFTER_MATCH_CODING */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_SHATTER_SPECIAL(mat,i,j,AFTER_MATCH_CODING) )  {  
          GenePhase6_SHATTER_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = SIG_0_3[MATCH] + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_3[INSERT] + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = SIG_0_3[DELETE] + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = GenePhase6_SHATTER_SPECIAL(mat,i-0,j-3,LOOP) + ((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = SIG_0_6[INTRON_0] + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = SIG_0_5[INTRON_1] + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = SIG_0_4[INTRON_2] + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_2[INSERT] + mat->gp->transition[GP4_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_1[INSERT] + mat->gp->transition[GP4_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         SIG_0_0[INSERT] = score;    


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_INSERT2END]) + (0) ;  
        if( temp > GenePhase6_SHATTER_SPECIAL(mat,i,j,LOOP) )    {  
          GenePhase6_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;   
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH] + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = SIG_1_0[INSERT] + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = SIG_1_0[DELETE] + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[DELETE] = score;    


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_DELETE2END]) + (0) ;  
        if( temp > GenePhase6_SHATTER_SPECIAL(mat,i,j,LOOP) )    {  
          GenePhase6_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;   
          }  




        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = SIG_0_8[MATCH] + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);   
        /* From state INSERT to state INTRON_0 */ 
        temp = SIG_0_8[INSERT] + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = SIG_0_1[INTRON_0] + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INTRON_0] = score;  


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = SIG_0_9[MATCH] + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);   
        /* From state INSERT to state INTRON_1 */ 
        temp = SIG_0_9[INSERT] + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = SIG_0_1[INTRON_1] + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INTRON_1] = score;  


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = SIG_0_10[MATCH] + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);  
        /* From state INSERT to state INTRON_2 */ 
        temp = SIG_0_10[INSERT] + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = SIG_0_1[INTRON_2] + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INTRON_2] = score;  


        /* Finished calculating state INTRON_2 */ 
        }  


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,LOOP);  


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,LOOP) + (mat->gp->transition[GP4_LOOP2LOOP]) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,LOOP) = score;  
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,END);   


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Source AFTER_MATCH_CODING is a special source for END */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,END) = score;   
      /* Finished updating state END */ 




      /* Special state BEFORE_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_CODING);   


      /* Source START is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 3,START) + (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 3,BEFORE_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 6,BEFORE_MATCH_INTRON_0) + (CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 5,BEFORE_MATCH_INTRON_1) + (CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 4,BEFORE_MATCH_INTRON_2) + (CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_CODING) = score;   
      /* Finished updating state BEFORE_MATCH_CODING */ 




      /* Special state BEFORE_MATCH_INTRON_0 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0); 


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_0) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0) = score; 
      /* Finished updating state BEFORE_MATCH_INTRON_0 */ 




      /* Special state BEFORE_MATCH_INTRON_1 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1); 


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 9,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_1) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1) = score; 
      /* Finished updating state BEFORE_MATCH_INTRON_1 */ 




      /* Special state BEFORE_MATCH_INTRON_2 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2); 


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 10,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_2) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2) = score; 
      /* Finished updating state BEFORE_MATCH_INTRON_2 */ 




      /* Special state AFTER_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_SHATTER_SPECIAL(mat,0,j,AFTER_MATCH_CODING);    


      /* Source MATCH for state AFTER_MATCH_CODING is not special... already calculated */ 
      /* Source AFTER_MATCH_CODING is a special source for AFTER_MATCH_CODING */ 
      temp = GenePhase6_SHATTER_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_SHATTER_SPECIAL(mat,0,j,AFTER_MATCH_CODING) = score;    
      /* Finished updating state AFTER_MATCH_CODING */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_GenePhase6(dbsi,out,query,targetdb,gp,general_model)
 *
 * Descrip:    This function makes a database search of GenePhase6
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:                 dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:                  out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                query [UNKN ] Undocumented argument [GenePhaseScore*]
 * Arg:             targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:                   gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Undocumented argument [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_GenePhase6(DBSearchImpl * dbsi,Hscore * out,GenePhaseScore* query,GenomicDB* targetdb ,GeneParser4Score * gp,GeneralGeneModelScore * general_model) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_GenePhase6 * holder;   
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_GenePhase6. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_GenePhase6. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for GenePhase6, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);    
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_GenePhase6(out,query, targetdb ,gp,general_model);  
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_GenePhase6 *) ckalloc(sizeof(struct thread_pool_holder_GenePhase6)); 
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->query = query;   
        holder->targetdb = targetdb; 
        holder->gp = gp; 
        holder->general_model = general_model;   
        if( pthread_mutex_init(&(holder->input_lock),NULL) != 0 )    
        fatal("Unable to iniated input mutex lock"); 
        if( pthread_mutex_init(&(holder->output_lock),NULL) != 0 )   
        fatal("Unable to iniated output mutex lock");    
        /* Let us rock! */ 
        thr_no = number_of_threads_DBSearchImpl(dbsi);   
        holder->pool = ckcalloc (thr_no,sizeof(pthread_t));  
        if( holder->pool == NULL )   {  
          warn("Unable to allocated thread pools");  
          return SEARCH_ERROR;   
          }  
        /* Build a thread attribute to make sure we get the most out of SMP boxes */ 
        pthread_attr_init(&pat);     
        /* Give thread libraries a hint that threads should be kernel threads */ 
#ifndef __sgi /* SGI can't set system scope ... */   
#ifdef  HAS_PTHREAD_SETSCOPE 
        pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   
#endif /* sgi */ 
        /* Give thread libraries a hint that there are num of threads to run */ 
#ifdef HAS_PTHREAD_SETCONCURRENCY    
        pthread_setconcurrency(thr_no+1);    
#endif /* set concurrency */ 
        for(i=0;i<thr_no;i++)    {  
          if( pthread_create(holder->pool+i,&pat,thread_loop_GenePhase6,(void *)holder) )    
            fatal("Unable to create a thread!"); 
          }  
        /* Now - wait for all the threads to exit */ 
        for(i=0;i<thr_no;i++)    {  
          if( pthread_join(holder->pool[i],NULL) != 0 )  
            fatal("Unable to join a thread!");   
          }  
        /* Deallocate the thread structures */ 
        ckfree(holder->pool);    
        ckfree(holder);  
        return SEARCH_OK;    
#else /* not compiled with threads */    
        warn("You did not specifiy the PTHREAD compile when compiled the C code for GenePhase6");    
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from GenePhase6",impl_string_DBSearchImpl(dbsi)); 
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_GenePhase6(ptr)
 *
 * Descrip:    dummy loop code foreach thread for GenePhase6
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * thread_loop_GenePhase6(void * ptr) 
{
    fatal("dummy thread loop function"); 
}    


/* Function:  serial_search_GenePhase6(out,query,targetdb,gp,general_model)
 *
 * Descrip:    This function makes a database search of GenePhase6
 *             It is a single processor implementation
 *
 *
 * Arg:                  out [UNKN ] Undocumented argument [Hscore *]
 * Arg:                query [UNKN ] Undocumented argument [GenePhaseScore*]
 * Arg:             targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:                   gp [UNKN ] Undocumented argument [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Undocumented argument [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_GenePhase6(Hscore * out,GenePhaseScore* query,GenomicDB* targetdb ,GeneParser4Score * gp,GeneralGeneModelScore * general_model) 
{
    ComplexSequence* target;     
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  


    target = init_GenomicDB(targetdb,&db_status);    
    if( db_status == DB_RETURN_ERROR )   {  
      warn("In searching GenePhase6, got a database init error on the target [target] database");    
      return SEARCH_ERROR;   
      }  
    for(;;)  { /*For all target entries*/ 


      /* No maximum length - allocated on-the-fly */ 
      score = score_only_GenePhase6(query, target , gp, general_model);  
      if( should_store_Hscore(out,score) == TRUE )   { /*if storing datascore*/ 
        ds = new_DataScore_from_storage(out);    
        if( ds == NULL )     {  
          warn("GenePhase6 search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");   
          return SEARCH_ERROR;   
          }  
        /* Now: add query/target information to the entry */ 
        dataentry_add_GenomicDB(ds->target,target,targetdb);     
        ds->score = score;   
        add_Hscore(out,ds);  
        } /* end of if storing datascore */ 
      pop_errormsg_stack();  
      push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);  


       target = reload_GenomicDB(target,targetdb,&db_status);    
      if( db_status == DB_RETURN_ERROR ) {  
        warn("In searching GenePhase6, Reload error on database target, position %d,%d",query_pos,target_pos);   
        return SEARCH_ERROR; 
        }  
      if( db_status == DB_RETURN_END )   
        break;  /* Out of target loop */ 
      target_pos++;  
      } /* end of For all target entries */ 
    close_GenomicDB(target,targetdb);    
    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_GenePhase6(query,target,gp,general_model)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_GenePhase6_only
 *
 *
 * Arg:                query [UNKN ] query data structure [GenePhaseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_GenePhase6(GenePhaseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    GenePhase6 * mat;    


    mat = allocate_GenePhase6_only(query, target , gp, general_model);   
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(11,(mat->leni + 1) * 6,11,8)) == NULL)    {  
      warn("Score only matrix for GenePhase6 cannot be allocated, (asking for 10  by %d  cells)",mat->leni*6);   
      mat = free_GenePhase6(mat);    
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<12;j++)    {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<6;k++) 
          GenePhase6_VSMALL_MATRIX(mat,i,j,k) = NEGI;    
        }  
      GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) = NEGI;    
      GenePhase6_VSMALL_SPECIAL(mat,i,j,START) = 0;  
      GenePhase6_VSMALL_SPECIAL(mat,i,j,END) = NEGI; 
      GenePhase6_VSMALL_SPECIAL(mat,i,j,BEFORE_MATCH_CODING) = NEGI; 
      GenePhase6_VSMALL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_0) = NEGI;   
      GenePhase6_VSMALL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_1) = NEGI;   
      GenePhase6_VSMALL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_2) = NEGI;   
      GenePhase6_VSMALL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = NEGI;  
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_VSMALL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = GenePhase6_VSMALL_SPECIAL(mat,i-1,j-3,START) + ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_MATCH_CODING to state MATCH */ 
        temp = GenePhase6_VSMALL_SPECIAL(mat,i-1,j-3,BEFORE_MATCH_CODING) + (mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_VSMALL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) )     {  
          GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* state MATCH is a source for special AFTER_MATCH_CODING */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_VSMALL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) )   {  
          GenePhase6_VSMALL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = temp;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_VSMALL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = GenePhase6_VSMALL_SPECIAL(mat,i-0,j-3,LOOP) + ((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_VSMALL_MATRIX(mat,i,j,INSERT) = score;   


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_INSERT2END]) + (0) ;  
        if( temp > GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) )     {  
          GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];  
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_VSMALL_MATRIX(mat,i,j,DELETE) = score;   


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_DELETE2END]) + (0) ;  
        if( temp > GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) )     {  
          GenePhase6_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;    
          }  




        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_VSMALL_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);  
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_VSMALL_MATRIX(mat,i,j,INTRON_0) = score; 


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_VSMALL_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);  
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_VSMALL_MATRIX(mat,i,j,INTRON_1) = score; 


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_VSMALL_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);     
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_VSMALL_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_VSMALL_MATRIX(mat,i,j,INTRON_2) = score; 


        /* Finished calculating state INTRON_2 */ 
        } /* end of for all query positions */ 




      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,LOOP);   


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,START) + (0) + (0);     
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,LOOP) + (mat->gp->transition[GP4_LOOP2LOOP]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,LOOP) = score;   
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,END);    


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);  
        if( temp > score )   
          score = temp;  
        }  


      /* Source AFTER_MATCH_CODING is a special source for END */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,END) = score;    
      /* Finished updating state END */ 




      /* Special state BEFORE_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_CODING);    


      /* Source START is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 3,START) + (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 3,BEFORE_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 6,BEFORE_MATCH_INTRON_0) + (CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 5,BEFORE_MATCH_INTRON_1) + (CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 4,BEFORE_MATCH_INTRON_2) + (CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_CODING) = score;    
      /* Finished updating state BEFORE_MATCH_CODING */ 




      /* Special state BEFORE_MATCH_INTRON_0 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0);  


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_0) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0) = score;  
      /* Finished updating state BEFORE_MATCH_INTRON_0 */ 




      /* Special state BEFORE_MATCH_INTRON_1 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1);  


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 9,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_1) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1) = score;  
      /* Finished updating state BEFORE_MATCH_INTRON_1 */ 




      /* Special state BEFORE_MATCH_INTRON_2 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2);  


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 10,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_2) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2) = score;  
      /* Finished updating state BEFORE_MATCH_INTRON_2 */ 




      /* Special state AFTER_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_VSMALL_SPECIAL(mat,0,j,AFTER_MATCH_CODING); 


      /* Source MATCH for state AFTER_MATCH_CODING is not special... already calculated */ 
      /* Source AFTER_MATCH_CODING is a special source for AFTER_MATCH_CODING */ 
      temp = GenePhase6_VSMALL_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_VSMALL_SPECIAL(mat,0,j,AFTER_MATCH_CODING) = score; 
      /* Finished updating state AFTER_MATCH_CODING */ 


      if( bestscore < GenePhase6_VSMALL_SPECIAL(mat,0,j,END) )   
        bestscore = GenePhase6_VSMALL_SPECIAL(mat,0,j,END);  
      } /* end of for all target positions */ 


    mat = free_GenePhase6(mat);  
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_GenePhase6(query,target,gp,general_model,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GenePhase6
 *
 *
 * Arg:                query [UNKN ] query data structure [GenePhaseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 * Arg:                dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:                 dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_GenePhase6(GenePhaseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    GenePhase6 * mat;    
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->len * target->seq->len;   
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 6 * sizeof(int)) > 1000*dpri->kbyte_size) {  
        strategy = DPIM_Linear;  
        }  
      else   {  
        strategy = DPIM_Explicit;    
        }  
      }  
    else {  
      strategy = dpri->memory;   
      }  


    if( dpenv != NULL )  {  
      if( strategy == DPIM_Explicit) {  
        if( (mat=allocate_Expl_GenePhase6(query, target , gp, general_model,dpri)) == NULL ) {  
          warn("Unable to allocate large GenePhase6 version");   
          return NULL;   
          }  
        calculate_dpenv_GenePhase6(mat,dpenv);   
        out =  PackAln_read_Expl_GenePhase6(mat);    
        }  
      else   {  
        mat = allocate_GenePhase6_only(query, target , gp, general_model);   
        calculate_shatter_GenePhase6(mat,dpenv);     
        out = PackAln_read_Shatter_GenePhase6(mat);  
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_GenePhase6(query, target , gp, general_model)) == NULL ) {  
          warn("Unable to allocate small GenePhase6 version");   
          return NULL;   
          }  
        out = PackAln_calculate_Small_GenePhase6(mat,dpenv);     
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_GenePhase6(query, target , gp, general_model,dpri)) == NULL ) {  
          warn("Unable to allocate large GenePhase6 version");   
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_GenePhase6(mat);     
          out =  PackAln_read_Expl_GenePhase6(mat);  
          }  
        }  
      }  


    mat = free_GenePhase6(mat);  
    return out;  
}    


/* Function:  allocate_GenePhase6_only(query,target,gp,general_model)
 *
 * Descrip:    This function only allocates the GenePhase6 structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:                query [UNKN ] query data structure [GenePhaseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhase6 *]
 *
 */
GenePhase6 * allocate_GenePhase6_only(GenePhaseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model) 
{
    GenePhase6 * out;    


    if((out= GenePhase6_alloc()) == NULL)    {  
      warn("Allocation of basic GenePhase6 structure failed...");    
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->gp = gp;    
    out->general_model = general_model;  
    out->leni = query->len;  
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_GenePhase6(query,target,gp,general_model,dpri)
 *
 * Descrip:    This function allocates the GenePhase6 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GenePhase6_only
 *
 *
 * Arg:                query [UNKN ] query data structure [GenePhaseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 * Arg:                 dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhase6 *]
 *
 */
GenePhase6 * allocate_Expl_GenePhase6(GenePhaseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model,DPRunImpl * dpri) 
{
    GenePhase6 * out;    


    out = allocate_GenePhase6_only(query, target , gp, general_model);   
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+10)*6 && dpri->cache->maxlenj >= (out->leni+1))   
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+10)*6,(out->leni+1),8,out->lenj+10)) == NULL)   {  
        warn("Explicit matrix GenePhase6 cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);    
        free_GenePhase6(out);    
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_GenePhase6(out);    
    return out;  
}    


/* Function:  init_GenePhase6(mat)
 *
 * Descrip:    This function initates GenePhase6 matrix when in explicit mode
 *             Called in /allocate_Expl_GenePhase6
 *
 *
 * Arg:        mat [UNKN ] GenePhase6 which contains explicit basematrix memory [GenePhase6 *]
 *
 */
void init_GenePhase6(GenePhase6 * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->len;i++)   {  
      for(j= (-10);j<11;j++) {  
        GenePhase6_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
        GenePhase6_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
        GenePhase6_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        }  
      }  
    for(j= (-10);j<mat->target->seq->len;j++)    {  
      for(i= (-1);i<2;i++)   {  
        GenePhase6_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
        GenePhase6_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
        GenePhase6_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        }  
      GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = NEGI;  
      GenePhase6_EXPL_SPECIAL(mat,i,j,START) = 0;    
      GenePhase6_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_CODING) = NEGI;   
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_0) = NEGI; 
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_1) = NEGI; 
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_2) = NEGI; 
      GenePhase6_EXPL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = NEGI;    
      }  
    return;  
}    


/* Function:  recalculate_PackAln_GenePhase6(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GenePhase6
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 */
void recalculate_PackAln_GenePhase6(PackAln * pal,GenePhase6 * mat) 
{
    int i,j,k,offi,offj; 
    PackAlnUnit * prev;  
    PackAlnUnit * pau;   


    for(k=1,prev=pal->pau[0];k < pal->len;k++,prev=pau)  {  
      pau = pal->pau[k]; 
      i = pau->i;    
      j = pau->j;    
      offi = pau->i - prev->i;   
      offj = pau->j - prev->j;   
      switch(pau->state) {  
        case MATCH :     
          if( offi == 1 && offj == 3 && prev->state == MATCH )   {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == INSERT )  {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offi == 1 && offj == 3 && prev->state == DELETE )  {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
            continue;    
            }  
          if( offj == 3 && prev->state == (START+6) )    {  
            pau->score = ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offj == 3 && prev->state == (BEFORE_MATCH_CODING+6) )  {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 6 && prev->state == INTRON_0 )    {  
            pau->score = (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == INTRON_1 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-2)) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == INTRON_2 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-1)) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 1 && offj == 2 && prev->state == MATCH )   {  
            pau->score = mat->gp->transition[GP4_DELETE_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->gp->transition[GP4_DELETE_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 1 && offj == 4 && prev->state == MATCH )   {  
            pau->score = mat->gp->transition[GP4_INSERT_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 1 && offj == 5 && prev->state == MATCH )   {  
            pau->score = mat->gp->transition[GP4_INSERT_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case INSERT :    
          if( offi == 0 && offj == 3 && prev->state == MATCH )   {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == INSERT )  {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == DELETE )  {  
            pau->score = (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offj == 3 && prev->state == (LOOP+6) ) {  
            pau->score = ((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 0 && offj == 6 && prev->state == INTRON_0 )    {  
            pau->score = (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 0 && offj == 5 && prev->state == INTRON_1 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-2)) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == INTRON_2 )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-1)) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
            continue;    
            }  
          if( offi == 0 && offj == 2 && prev->state == INSERT )  {  
            pau->score = mat->gp->transition[GP4_DELETE_1_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INSERT )  {  
            pau->score = mat->gp->transition[GP4_DELETE_2_BASE] + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
            continue;    
            }  
          warn("In recaluclating PackAln with state INSERT, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case DELETE :    
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = mat->query->gws->seg[i]->transition[GW_MATCH2DELETE] + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == INSERT )  {  
            pau->score = mat->query->gws->seg[i]->transition[GW_INSERT2DELETE] + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == DELETE )  {  
            pau->score = mat->query->gws->seg[i]->transition[GW_DELETE2DELETE] + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state DELETE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case INTRON_0 :  
          if( offi == 0 && offj == 8 && prev->state == MATCH )   {  
            pau->score = ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 8 && prev->state == INSERT )  {  
            pau->score = ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_0 )    {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_0, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case INTRON_1 :  
          if( offi == 0 && offj == 9 && prev->state == MATCH )   {  
            pau->score = ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 9 && prev->state == INSERT )  {  
            pau->score = ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_1 )    {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_1, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case INTRON_2 :  
          if( offi == 0 && offj == 10 && prev->state == MATCH )  {  
            pau->score = ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 10 && prev->state == INSERT ) {  
            pau->score = ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron) + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INTRON_2 )    {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state INTRON_2, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case (LOOP+6) :  
          if( offj == 1 && prev->state == (START+6) )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->gws->seg[i]->transition[GW_MATCH2END] + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == INSERT )   {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->gws->seg[i]->transition[GW_INSERT2END] + (0);   
            continue;    
            }  
          if( offj == 0 && prev->state == DELETE )   {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->gws->seg[i]->transition[GW_DELETE2END] + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (LOOP+6) ) {  
            pau->score = mat->gp->transition[GP4_LOOP2LOOP] + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state LOOP, got a bad source state. Error!");  
          break; 
        case (START+6) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+6) :   
          if( offj == 1 && prev->state == (LOOP+6) ) {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (AFTER_MATCH_CODING+6) )   {  
            pau->score = mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)] + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state END, got a bad source state. Error!");   
          break; 
        case (BEFORE_MATCH_CODING+6) :   
          if( offj == 3 && prev->state == (START+6) )    {  
            pau->score = mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)] + (0);  
            continue;    
            }  
          if( offj == 3 && prev->state == (BEFORE_MATCH_CODING+6) )  {  
            pau->score = mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)] + (0);    
            continue;    
            }  
          if( offj == 6 && prev->state == (BEFORE_MATCH_INTRON_0+6) )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-3)) + (0);  
            continue;    
            }  
          if( offj == 5 && prev->state == (BEFORE_MATCH_INTRON_1+6) )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-2)) + (0);  
            continue;    
            }  
          if( offj == 4 && prev->state == (BEFORE_MATCH_INTRON_2+6) )    {  
            pau->score = CSEQ_GENOMIC_3SS(mat->target,(j-1)) + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_MATCH_CODING, got a bad source state. Error!");   
          break; 
        case (BEFORE_MATCH_INTRON_0+6) :     
          if( offj == 1 && prev->state == (BEFORE_MATCH_CODING+6) )  {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j)) + (0);  
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_MATCH_INTRON_0+6) )    {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_MATCH_INTRON_0, got a bad source state. Error!"); 
          break; 
        case (BEFORE_MATCH_INTRON_1+6) :     
          if( offj == 9 && prev->state == (BEFORE_MATCH_CODING+6) )  {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);  
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_MATCH_INTRON_1+6) )    {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_MATCH_INTRON_1, got a bad source state. Error!"); 
          break; 
        case (BEFORE_MATCH_INTRON_2+6) :     
          if( offj == 10 && prev->state == (BEFORE_MATCH_CODING+6) ) {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7))) + (0);  
            continue;    
            }  
          if( offj == 1 && prev->state == (BEFORE_MATCH_INTRON_2+6) )    {  
            pau->score = (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state BEFORE_MATCH_INTRON_2, got a bad source state. Error!"); 
          break; 
        case (AFTER_MATCH_CODING+6) :    
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->query->gws->seg[i]->transition[GW_MATCH2END] + (0);    
            continue;    
            }  
          if( offj == 3 && prev->state == (AFTER_MATCH_CODING+6) )   {  
            pau->score = mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)] + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state AFTER_MATCH_CODING, got a bad source state. Error!");    
          break; 
        default :    
          warn("In recaluclating PackAln got a bad recipient state. Error!");    
        }  
      prev = pau;    
      }  
    return;  
}    
/* divide and conquor macros are next */ 
#define GenePhase6_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+10)][(i+1)*6+state])   
#define GenePhase6_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+11)*8) % 88][(i+1)*6+state])   
#define GenePhase6_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+10)])  
#define GenePhase6_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)]) 
#define GenePhase6_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+11)*8)+(shadow+1)) % 88)][(i+1)*6 + state])    
#define GenePhase6_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)]) 
#define GenePhase6_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+10)% 10) * (leni+1) * 6) + ((i+1) * 6) + (state)])    
#define GenePhase6_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+10)% 10) * (leni+1) * 48) + ((i+1) * 48) + (state * 8) + shadow+1])    
#define GenePhase6_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+10)]) 
/* Function:  allocate_Small_GenePhase6(query,target,gp,general_model)
 *
 * Descrip:    This function allocates the GenePhase6 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GenePhase6_only
 *
 *
 * Arg:                query [UNKN ] query data structure [GenePhaseScore*]
 * Arg:               target [UNKN ] target data structure [ComplexSequence*]
 * Arg:                   gp [UNKN ] Resource [GeneParser4Score *]
 * Arg:        general_model [UNKN ] Resource [GeneralGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhase6 *]
 *
 */
#define GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+10)]) 
GenePhase6 * allocate_Small_GenePhase6(GenePhaseScore* query,ComplexSequence* target ,GeneParser4Score * gp,GeneralGeneModelScore * general_model) 
{
    GenePhase6 * out;    


    out = allocate_GenePhase6_only(query, target , gp, general_model);   
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(88,(out->leni + 1) * 6,64,out->lenj+10);  
    if(out == NULL)  {  
      warn("Small shadow matrix GenePhase6 cannot be allocated, (asking for 11 by %d main cells)",out->leni+2);  
      free_GenePhase6(out);  
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_GenePhase6(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GenePhase6 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GenePhase6 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GenePhase6 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_GenePhase6(GenePhase6 * mat,DPEnvelope * dpenv) 
{
    int endj;    
    int score;   
    PackAln * out;   
    PackAlnUnit * pau;   
    int starti;  
    int startj;  
    int startstate;  
    int stopi;   
    int stopj;   
    int stopstate;   
    int temp;    
    int donej;  /* This is for reporting, will be passed as a & arg in */ 
    int totalj; /* This also is for reporting, but as is not changed, can be passed by value */ 


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW )    {  
      warn("Could not calculate packaln small for GenePhase6 due to wrong type of matrix");  
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_GenePhase6(mat,dpenv);   
    score = start_end_find_end_GenePhase6(mat,&endj);    
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_GenePhase6(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )   {  
        warn("Problem in reading off special state system... going to return partial alignment");    
        break;   
        }  
      if( startstate == START || endj <= 0)  {  
        log_full_error(REPORT,0,"Recovered complete alignment"); 
        stop_reporting();    
        break;   
        }  


      log_full_error(REPORT,0,"Finished to %d",endj);    
      stop_reporting();  


      /* Ok... have to eat up another piece of matrix <sigh> */ 
      temp = startstate; 
      starti = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);   
      startj = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);   
      startstate = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);   
      stopi = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3);    
      stopj = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4);    
      stopstate = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5);    


      /* Get out the score of this block. V. important! */ 
      temp = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6); 
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE)    {  
        warn("In the alignment GenePhase6 [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj);  
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_GenePhase6(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);   
      if( stopi == GenePhase6_READ_OFF_ERROR)    {  
        warn("In GenePhase6 read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj);  
        invert_PackAln(out); 
        recalculate_PackAln_GenePhase6(out,mat); 
        return out;  
        }  
      /* if at start, break, otherwise, back to eat another strip */ 
      if( stopstate == START)    {  
        log_full_error(REPORT,0,"Recovered complete alignment      ");   
        stop_reporting();    
        break;   
        }  
      log_full_error(REPORT,0,"Finished  alignment to %d           ",startj);    
      stop_reporting();  
      endj = stopj;  
      /* stopstate is correct as it is */ 
      } /* end of while there are more special bits to recover */ 
    invert_PackAln(out); 
    recalculate_PackAln_GenePhase6(out,mat); 
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_GenePhase6(mat)
 *
 * Descrip:    This function calculates an alignment for GenePhase6 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GenePhase6 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GenePhase6
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GenePhase6 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_GenePhase6(GenePhase6 * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_GenePhase6(mat,NULL);    
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_GenePhase6(mat);  
    return out;  
}    


/* Function:  AlnRangeSet_from_GenePhase6(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GenePhase6 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GenePhase6 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GenePhase6
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_GenePhase6(GenePhase6 * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_GenePhase6");    
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_GenePhase6(mat,&jpos);   
    state = END; 


    while( (temp = AlnRange_build_GenePhase6(mat,jpos,state,&jpos,&state)) != NULL)  
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_GenePhase6(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_GenePhase6
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_GenePhase6(GenePhase6 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_GenePhase6");    
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_GenePhase6(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE) {  
      warn("In AlnRanger_build_GenePhase6 alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj); 
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);   
    out->startj = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);   
    out->startstate = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);   
    out->stopi = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);    
    out->stopj = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);    
    out->stopstate = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);    
    out->startscore = GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);   
    out->stopscore = GenePhase6_DC_SHADOW_SPECIAL(mat,0,jpos,state); 


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_GenePhase6(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);   
    if( jpos == GenePhase6_READ_OFF_ERROR)   {  
      warn("In AlnRange_build_GenePhase6 alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);  
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = GenePhase6_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);   
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:               out [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean read_hidden_GenePhase6(GenePhase6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
{
    int i;   
    int j;   
    int state;   
    int cellscore;   
    int isspecial;   
    /* We don't need hiddenj here, 'cause matrix access handled by max funcs */ 
    PackAlnUnit * pau;   


    /* stop position is on the path */ 
    i = stopi;   
    j = stopj;   
    state= stopstate;    
    isspecial = FALSE;   


    while( i >= starti && j >= startj)   {  
      /* Put away current i,j,state */ 
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state;   
      add_PackAln(out,pau);  


      max_hidden_GenePhase6(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  


      if( i == GenePhase6_READ_OFF_ERROR)    {  
        warn("In GenePhase6 hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);    
        return FALSE;    
        }  


      if( i == starti && j == startj && state == startstate) {  
/* Put away final state (start of this block) */ 
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state; 
        add_PackAln(out,pau);    
          return TRUE;   
        }  
      if( i == starti && j == startj)    {  
        warn("In GenePhase6 hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);  
        return FALSE;    
        }  
      }  
    warn("In GenePhase6 hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);   
    return FALSE;    
}    


/* Function:  max_hidden_GenePhase6(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:           hiddenj [UNKN ] Undocumented argument [int]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_hidden_GenePhase6(GenePhase6 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GenePhase6_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In GenePhase6 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = GenePhase6_HIDDEN_MATRIX(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->gp->transition[GP4_INSERT_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-5,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 5,MATCH);    
          }  
        temp = cscore - (mat->gp->transition[GP4_INSERT_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-4,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 4,MATCH);    
          }  
        temp = cscore - (mat->gp->transition[GP4_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);    
          }  
        temp = cscore - (mat->gp->transition[GP4_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-2,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 2,MATCH);    
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-1))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 4,INTRON_2) ) {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-4,INTRON_2);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 4,INTRON_2);     
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-2))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 5,INTRON_1) ) {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-5,INTRON_1);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 5,INTRON_1);     
          }  
        temp = cscore - ((mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 6,INTRON_0) ) {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-6,INTRON_0);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 6,INTRON_0);     
          }  
        /* Not allowing special sources.. skipping BEFORE_MATCH_CODING */ 
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-3,DELETE);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-3,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-3,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->gp->transition[GP4_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT);   
          }  
        temp = cscore - (mat->gp->transition[GP4_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 2,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-2,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 2,INSERT);   
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-1))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 4,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-4,INTRON_2);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 4,INTRON_2);     
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-2))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 5,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-5,INTRON_1);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 5,INTRON_1);     
          }  
        temp = cscore - ((mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 6,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-6,INTRON_0);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 6,INTRON_0);     
          }  
        /* Not allowing special sources.. skipping LOOP */ 
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-3,DELETE);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-3,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-3,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->gws->seg[i]->transition[GW_DELETE2DELETE]) -  (0);  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-0,DELETE);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE);   
          }  
        temp = cscore - (mat->query->gws->seg[i]->transition[GW_INSERT2DELETE]) -  (0);  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-0,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 0,INSERT);   
          }  
        temp = cscore - (mat->query->gws->seg[i]->transition[GW_MATCH2DELETE]) -  (0);   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-1,j-0,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_0 :    
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) -  (0); 
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_0);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_0);     
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron)) -  (0);  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 8,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-8,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 8,INSERT);   
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0)) -  (0);   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 8,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-8,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 8,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_1 :    
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) -  (0); 
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_1);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_1);     
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron)) -  (0);  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 9,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-9,INSERT);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 9,INSERT);   
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1)) -  (0);   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 9,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-9,MATCH);   
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 9,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_2 :    
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) -  (0); 
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_2);    
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 1,INTRON_2);     
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron)) -  (0);  
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-10,INSERT); 
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2)) -  (0);   
        if( temp == GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 10,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_HIDDEN_MATRIX(mat,i-0,j-10,MATCH);  
            }  
          return GenePhase6_HIDDEN_MATRIX(mat,i - 0,j - 10,MATCH);   
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_GenePhase6(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int *]
 * Arg:        startstate [UNKN ] Undocumented argument [int *]
 * Arg:               out [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean read_special_strip_GenePhase6(GenePhase6 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
{
    int i;   
    int j;   
    int state;   
    int cellscore;   
    int isspecial;   
    PackAlnUnit * pau;   


    /* stop position is on the path */ 
    i = stopi;   
    j = stopj;   
    state= stopstate;    
    isspecial = TRUE;    


    /* Loop until state has the same j as its stop in shadow pointers */ 
    /* This will be the state is came out from, OR it has hit !start */ 
    /* We may not want to get the alignment, in which case out will be NULL */ 
    while( j > GenePhase6_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)   { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 6; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_GenePhase6(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  
      if( i == GenePhase6_READ_OFF_ERROR)    {  
        warn("In special strip read GenePhase6, got a bad read off error. Sorry!");  
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < GenePhase6_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))  {  
      warn("In special strip read GenePhase6, at special [%d] state [%d] overshot!",j,state);    
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 6;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_GenePhase6(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_special_strip_GenePhase6(GenePhase6 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = GenePhase6_READ_OFF_ERROR;   
    if( isspecial == FALSE ) {  
      warn("In special strip max function for GenePhase6, got a non special start point. Problem! (bad!)");  
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In GenePhase6 matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    cscore = GenePhase6_DC_SHADOW_SPECIAL(mat,i,j,state);    
    switch(state)    { /*switch on special states*/ 
      case LOOP :    
        /* source LOOP is a special */ 
        temp = cscore - (mat->gp->transition[GP4_LOOP2LOOP]) - (0);  
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,LOOP) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,LOOP);    
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,LOOP) ;     
          }  
        /* Source DELETE is not a special */ 
        /* Source INSERT is not a special */ 
        /* Source MATCH is not a special */ 
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,START) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,START);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,START) ;    
          }  
      case START :   
      case END :     
        /* source AFTER_MATCH_CODING is a special */ 
        temp = cscore - (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);  
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,AFTER_MATCH_CODING) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_MATCH_CODING;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-3,AFTER_MATCH_CODING);  
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 3,AFTER_MATCH_CODING) ;   
          }  
        /* source LOOP is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,LOOP) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,LOOP);    
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,LOOP) ;     
          }  
      case BEFORE_MATCH_CODING :     
        /* source BEFORE_MATCH_INTRON_2 is a special */ 
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-1))) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 4,BEFORE_MATCH_INTRON_2) )    {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = BEFORE_MATCH_INTRON_2; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-4,BEFORE_MATCH_INTRON_2);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 4,BEFORE_MATCH_INTRON_2) ;    
          }  
        /* source BEFORE_MATCH_INTRON_1 is a special */ 
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-2))) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 5,BEFORE_MATCH_INTRON_1) )    {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = BEFORE_MATCH_INTRON_1; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-5,BEFORE_MATCH_INTRON_1);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 5,BEFORE_MATCH_INTRON_1) ;    
          }  
        /* source BEFORE_MATCH_INTRON_0 is a special */ 
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-3))) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 6,BEFORE_MATCH_INTRON_0) )    {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = BEFORE_MATCH_INTRON_0; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-6,BEFORE_MATCH_INTRON_0);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 6,BEFORE_MATCH_INTRON_0) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);   
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,BEFORE_MATCH_CODING) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-3,BEFORE_MATCH_CODING);     
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 3,BEFORE_MATCH_CODING) ;  
          }  
        /* source START is a special */ 
        temp = cscore - (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,START) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-3,START);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 3,START) ;    
          }  
      case BEFORE_MATCH_INTRON_0 :   
        /* source BEFORE_MATCH_INTRON_0 is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) - (0);  
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_0) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_INTRON_0; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_MATCH_INTRON_0);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_0) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_MATCH_CODING) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_MATCH_CODING);     
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_MATCH_CODING) ;  
          }  
      case BEFORE_MATCH_INTRON_1 :   
        /* source BEFORE_MATCH_INTRON_1 is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) - (0);  
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_1) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_INTRON_1; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_MATCH_INTRON_1);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_1) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 9,BEFORE_MATCH_CODING) )  {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-9,BEFORE_MATCH_CODING);     
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 9,BEFORE_MATCH_CODING) ;  
          }  
      case BEFORE_MATCH_INTRON_2 :   
        /* source BEFORE_MATCH_INTRON_2 is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) - (0);  
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_2) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_INTRON_2; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-1,BEFORE_MATCH_INTRON_2);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_2) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 10,BEFORE_MATCH_CODING) ) {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-10,BEFORE_MATCH_CODING);    
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 10,BEFORE_MATCH_CODING) ;     
          }  
      case AFTER_MATCH_CODING :  
        /* source AFTER_MATCH_CODING is a special */ 
        temp = cscore - (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);   
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,AFTER_MATCH_CODING) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_MATCH_CODING;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-3,AFTER_MATCH_CODING);  
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 3,AFTER_MATCH_CODING) ;   
          }  
        /* Source MATCH is not a special */ 
      default:   
        warn("Major problem (!) - in GenePhase6 special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_GenePhase6(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:            cscore [UNKN ] Undocumented argument [int]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_matrix_to_special_GenePhase6(GenePhase6 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = GenePhase6_READ_OFF_ERROR;   


    if( j < 0 || j > mat->lenj)  {  
      warn("In GenePhase6 matrix to special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source INTRON_0 is not a special, should not get here! */ 
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,BEFORE_MATCH_CODING) )  {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-1,j-3,BEFORE_MATCH_CODING);     
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 1,j - 3,BEFORE_MATCH_CODING) ;  
          }  
        temp = cscore - (((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 1,j - 3,START) )    {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-1,j-3,START);   
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 1,j - 3,START) ;    
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GenePhase6 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INSERT :  
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source INTRON_0 is not a special, should not get here! */ 
        temp = cscore - (((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == GenePhase6_DC_SHADOW_SPECIAL(mat,i - 0,j - 3,LOOP) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - GenePhase6_DC_SHADOW_SPECIAL(mat,i-0,j-3,LOOP);    
            }  
          return GenePhase6_DC_SHADOW_MATRIX(mat,i - 0,j - 3,LOOP) ;     
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GenePhase6 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DELETE :  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GenePhase6 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INTRON_0 :    
        /* Source INTRON_0 is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GenePhase6 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INTRON_1 :    
        /* Source INTRON_1 is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GenePhase6 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case INTRON_2 :    
        /* Source INTRON_2 is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in GenePhase6 matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_GenePhase6(GenePhase6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_GenePhase6(mat,starti,startj,stopi,stopj);   


    GenePhase6_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;  


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenePhase6_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;    
          GenePhase6_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;   
          GenePhase6_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;   
          GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_0) = NEGI;     
          GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_1) = NEGI;     
          GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_2) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_HIDDEN_MATRIX(mat,i,j,MATCH) = score;    
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_HIDDEN_MATRIX(mat,i,j,INSERT) = score;   
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];  
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_HIDDEN_MATRIX(mat,i,j,DELETE) = score;   
        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);  
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_0) = score; 
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);  
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_1) = score; 
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);     
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_HIDDEN_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_2) = score; 
        /* Finished calculating state INTRON_2 */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_GenePhase6(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_GenePhase6(GenePhase6 * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-10);j<=stopj;j++)  {  
      for(i=(starti-1);i<=stopi;i++) {  
        GenePhase6_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
 
        GenePhase6_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;
    
        GenePhase6_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;
    
        GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_0) = NEGI;
  
        GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_1) = NEGI;
  
        GenePhase6_HIDDEN_MATRIX(mat,i,j,INTRON_2) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_GenePhase6
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_GenePhase6 to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [GenePhase6 *]
 * Arg:            starti [UNKN ] Start position in i [int]
 * Arg:            startj [UNKN ] Start position in j [int]
 * Arg:        startstate [UNKN ] Start position state number [int]
 * Arg:             stopi [UNKN ] Stop position in i [int]
 * Arg:             stopj [UNKN ] Stop position in j [int]
 * Arg:         stopstate [UNKN ] Stop position state number [int]
 * Arg:               out [UNKN ] PackAln structure to put alignment into [PackAln *]
 * Arg:             donej [UNKN ] pointer to a number with the amount of alignment done [int *]
 * Arg:            totalj [UNKN ] total amount of alignment to do (in j coordinates) [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean full_dc_GenePhase6(GenePhase6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_GenePhase6");  
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 50) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);     
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)   {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)    {  
      warn("In divide and conquor for GenePhase6, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);   
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= GenePhase6_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);    
    lstartj= GenePhase6_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);    
    lstate = GenePhase6_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);    


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_GenePhase6(mat,GenePhase6_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),GenePhase6_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),GenePhase6_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_GenePhase6(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE) {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_GenePhase6(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         perc_done [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean do_dc_single_pass_GenePhase6(GenePhase6 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_GenePhase6(mat); 


    GenePhase6_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;   
    run_up_dc_GenePhase6(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);   
    push_dc_at_merge_GenePhase6(mat,starti,stopi,halfj,&halfj,dpenv);    
    follow_on_dc_GenePhase6(mat,starti,stopi,halfj,stopj,dpenv,perc_done);   
    return TRUE; 
}    


/* Function:  push_dc_at_merge_GenePhase6(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_GenePhase6(GenePhase6 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<10;count++,j++)   {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);  
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-100); 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = (-100); 
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-100); 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = (-100); 
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = (-100);   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = (-100);   
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = (-100);   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = (-100);   
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = (-100);   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( j - 3 <= mergej) {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1; 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
          }  
        else {  
          for(k=0;k<7;k++)   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);   
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k);    
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k);    
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-6;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INTRON_0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,INTRON_0,k);  
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INTRON_1;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INTRON_1,k);  
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INTRON_2;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INTRON_2,k);  
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-5;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += CSEQ_GENOMIC_CDSPOT(mat->target,j);     
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished with state MATCH */ 


        /* For state INSERT, pushing when j - offj <= mergej */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( j - 3 <= mergej) {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k);   
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-3;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = DELETE;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k);   
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp > score)    {  
          score = temp;  


          if( j - 6 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-6;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INTRON_0; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,INTRON_0,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp > score)    {  
          score = temp;  


          if( j - 5 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-5;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INTRON_1; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INTRON_1,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-4;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INTRON_2; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INTRON_2,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-2;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,k);   
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += CSEQ_GENOMIC_CDSPOT(mat->target,j);     
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;     
        /* Finished with state INSERT */ 


        /* For state DELETE, pushing when j - offj <= mergej */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];   
        if( j - 0 <= mergej) {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = MATCH;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = INSERT;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);   
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = DELETE;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;     
        /* Finished with state DELETE */ 


        /* For state INTRON_0, pushing when j - offj <= mergej */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);   
        if( j - 8 <= mergej) {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-8;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = MATCH;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,k);    
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp > score)    {  
          score = temp;  


          if( j - 8 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-8;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = INSERT; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = i-0;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,1) = j-1;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,2) = INTRON_0;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,3) = i;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,4) = j;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,5) = INTRON_0;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;   
        /* Finished with state INTRON_0 */ 


        /* For state INTRON_1, pushing when j - offj <= mergej */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);   
        if( j - 9 <= mergej) {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-9;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = MATCH;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,k);    
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp > score)    {  
          score = temp;  


          if( j - 9 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-9;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = INSERT; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,k); 
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = i-0;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,1) = j-1;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,2) = INTRON_1;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,3) = i;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,4) = j;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,5) = INTRON_1;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;   
        /* Finished with state INTRON_1 */ 


        /* For state INTRON_2, pushing when j - offj <= mergej */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);  
        if( j - 10 <= mergej)    {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0;  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-10; 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = MATCH;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j;    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2; 
          }  
        else {  
          for(k=0;k<7;k++)   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,k);   
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp > score)    {  
          score = temp;  


          if( j - 10 <= mergej)  {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-10;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = INSERT; 
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);    
            }  
          }  


        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = i-0;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,1) = j-1;    
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,2) = INTRON_2;   
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,3) = i;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,4) = j;  
            GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,5) = INTRON_2;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k);   
            }  
          }  
        /* Add any movement independant score */ 
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;   
        /* Finished with state INTRON_2 */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_GenePhase6(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_GenePhase6(GenePhase6 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
{
    int i;   
    int j;   
    int k;   
    int score;   
    int temp;    
    int localshadow[7];  
    long int total;  
    long int num;    


    total = (stopi - starti+1) * (stopj - startj+1); 
    num = 0;     


    for(j=startj;j<=stopj;j++)   { /*for each valid j column*/ 
      for(i=starti;i<=stopi;i++) { /*this is strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,k);  
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,k);   
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,k);   
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 6,INTRON_0,k); 
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INTRON_1,k); 
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INTRON_2,k); 
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,k);    
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score; 
        for(k=0;k<7;k++) 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,k);  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,k);   
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,k);   
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 6,INTRON_0,k); 
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INTRON_1,k); 
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INTRON_2,k); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,k);   
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k);   
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        for(k=0;k<7;k++) 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);  
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,k);   
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);   
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        for(k=0;k<7;k++) 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,k);  
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,k);   
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,k); 
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;  
        for(k=0;k<7;k++) 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,k);  
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,k);   
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,k); 
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;  
        for(k=0;k<7;k++) 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,k); 
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,k);  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = GenePhase6_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,k); 
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;  
        for(k=0;k<7;k++) 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INTRON_2 */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_GenePhase6(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_GenePhase6(GenePhase6 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    long int total;  
    long int num;    


    total = (stopi - starti+1) * (stopj - startj+1); 
    if( total <= 0 ) 
      total = 1; 
    num = 0;     


    for(j=startj;j<=stopj;j++)   { /*for each valid j column*/ 
      for(i=starti;i<=stopi;i++) { /*this is strip*/ 
        if( j == startj && i == starti)  
          continue;  
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score; 
        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];   
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);   
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;  
        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);   
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;  
        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);  
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_DC_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;  
        /* Finished calculating state INTRON_2 */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_GenePhase6(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 */
}    
void init_dc_GenePhase6(GenePhase6 * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<12;j++)    {  
      for(i=(-1);i<mat->query->len;i++)  {  
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;    
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;    
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;    
        for(k=0;k<7;k++) {  
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1);    
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = (-1);   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = (-1);   
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = (-1); 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = (-1); 
          GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_GenePhase6(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [GenePhase6 *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_GenePhase6(GenePhase6 * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = GenePhase6_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END);   
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,END) > max )  {  
        max = GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,END); 
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_GenePhase6(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [GenePhase6]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_GenePhase6(GenePhase6 *mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int score;   
    int temp;    
    int leni;    
    int lenj;    
    int localshadow[7];  
    long int total;  
    long int num=0;  
    int * score_pointers;    
    int * shadow_pointers;   
    int * localsp;   
    leni = mat->query->len;  
    lenj = mat->target->seq->len;    
    total = leni * lenj; 


    score_pointers = (int *) calloc (10 * (leni + 1) * 6,sizeof(int));   
    shadow_pointers = (int *) calloc (10 * (leni + 1) * 6 * 8,sizeof(int));  


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;     
          GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;    
          GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;    
          GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;  
          GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;  
          GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        /* assign local shadown pointer */ 
        localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,MATCH,0));    
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,INSERT,0)); 
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 3,DELETE,0)); 
          }  
        /* From state START to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,START) + ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state BEFORE_MATCH_CODING to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-3,BEFORE_MATCH_CODING) + (mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* This state [BEFORE_MATCH_CODING] is a special for MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 6,INTRON_0,0));   
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2)) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,INTRON_1,0));   
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1)) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,INTRON_2,0));   
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 2,MATCH,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 4,MATCH,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 5,MATCH,0));  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score; 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )  {  
          GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6); 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i;   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j;   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = MATCH;   
          }  




        /* state MATCH is a source for special AFTER_MATCH_CODING */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_MATCH_CODING) )    {  
          GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k); 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_MATCH_CODING,6) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_MATCH_CODING,3) = i; 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_MATCH_CODING,4) = j; 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_MATCH_CODING,5) = MATCH; 
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        /* assign local shadown pointer */ 
        localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,MATCH,0));    
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,INSERT,0)); 
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,DELETE,0)); 
          }  
        /* From state LOOP to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-3,LOOP) + ((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL]) + (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp  > score )  {  
          score = temp;  
          /* This state [LOOP] is a special for INSERT... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= INSERT;    
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3))) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 6,INTRON_0,0));   
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2)) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 5,INTRON_1,0));   
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1)) +(CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,INTRON_2,0));   
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,INSERT,0)); 
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE] +(CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,0)); 
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_INSERT2END]) + (0) ;  
        if( temp > GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )  {  
          GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k);  
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,6);    
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i;   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j;   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = INSERT;  
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE] + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));    
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE] +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,INSERT,0)); 
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE] +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,0)); 
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_DELETE2END]) + (0) ;  
        if( temp > GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )  {  
          GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k);  
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,6);    
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i;   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j;   
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = DELETE;  
          }  




        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,MATCH,0));    
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 8,INSERT,0)); 
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_0,0));   
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_0) = score;  
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,MATCH,0));    
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 9,INSERT,0)); 
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_1,0));   
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_1) = score;  
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,MATCH,0));   
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron) +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 10,INSERT,0));    
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INTRON_2,0));   
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         GenePhase6_DC_OPT_SHADOW_MATRIX(mat,i,j,INTRON_2) = score;  
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INTRON_2 */ 


        } /* end of for each i position in strip */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,LOOP);    


      /* Source START is a special source for LOOP */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,START) + (0) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,START,k);    
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,LOOP) + (mat->gp->transition[GP4_LOOP2LOOP]) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,LOOP,k); 
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,LOOP) = score;    
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,END); 


      /* Source LOOP is a special source for END */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);     
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,LOOP,k);  
        }  


      /* Source AFTER_MATCH_CODING is a special source for END */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,AFTER_MATCH_CODING,k);    
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,END) = score; 
      /* Finished updating state END */ 




      /* Special state BEFORE_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_CODING); 


      /* Source START is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,START) + (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,START,k); 
        }  


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,BEFORE_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,BEFORE_MATCH_CODING,k);   
        }  


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 6,BEFORE_MATCH_INTRON_0) + (CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 6,BEFORE_MATCH_INTRON_0,k); 
        }  


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 5,BEFORE_MATCH_INTRON_1) + (CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 5,BEFORE_MATCH_INTRON_1,k); 
        }  


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 4,BEFORE_MATCH_INTRON_2) + (CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);  
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 4,BEFORE_MATCH_INTRON_2,k); 
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_CODING) = score; 
      /* Finished updating state BEFORE_MATCH_CODING */ 




      /* Special state BEFORE_MATCH_INTRON_0 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0);   


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_INTRON_0,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_MATCH_CODING,k); 
        }  


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_0) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_INTRON_0,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_0,k);   
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0) = score;   
      /* Finished updating state BEFORE_MATCH_INTRON_0 */ 




      /* Special state BEFORE_MATCH_INTRON_1 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1);   


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 9,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_INTRON_1,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 9,BEFORE_MATCH_CODING,k); 
        }  


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_1) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_INTRON_1,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_1,k);   
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1) = score;   
      /* Finished updating state BEFORE_MATCH_INTRON_1 */ 




      /* Special state BEFORE_MATCH_INTRON_2 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2);   


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 10,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_INTRON_2,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 10,BEFORE_MATCH_CODING,k);    
        }  


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_2) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,BEFORE_MATCH_INTRON_2,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_2,k);   
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2) = score;   
      /* Finished updating state BEFORE_MATCH_INTRON_2 */ 




      /* Special state AFTER_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_MATCH_CODING);  


      /* Source MATCH for state AFTER_MATCH_CODING is not special... already calculated */ 
      /* Source AFTER_MATCH_CODING is a special source for AFTER_MATCH_CODING */ 
      temp = GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,AFTER_MATCH_CODING,k) = GenePhase6_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 3,AFTER_MATCH_CODING,k); 
        }  


      /* Put back score... (now updated!) */ 
      GenePhase6_DC_OPT_SHADOW_SPECIAL(mat,0,j,AFTER_MATCH_CODING) = score;  
      /* Finished updating state AFTER_MATCH_CODING */ 


      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_GenePhase6(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 */
void init_start_end_linear_GenePhase6(GenePhase6 * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<12;j++)    {  
      for(i=(-1);i<mat->query->len;i++)  {  
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
        GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);  
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;  
        GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-1); 
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;  
        GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-1); 
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_0) = NEGI;    
        GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_0,0) = (-1);   
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_1) = NEGI;    
        GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_1,0) = (-1);   
        GenePhase6_DC_SHADOW_MATRIX(mat,i,j,INTRON_2) = NEGI;    
        GenePhase6_DC_SHADOW_MATRIX_SP(mat,i,j,INTRON_2,0) = (-1);   
        }  
      }  


    for(j=(-10);j<mat->target->seq->len;j++) {  
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,LOOP) = NEGI; 
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,LOOP,0) = (-1);    
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;   
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;  
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;  
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1); 
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_CODING) = NEGI;  
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_MATCH_CODING,0) = (-1); 
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0) = NEGI;    
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_MATCH_INTRON_0,0) = (-1);   
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1) = NEGI;    
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_MATCH_INTRON_1,0) = (-1);   
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2) = NEGI;    
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,BEFORE_MATCH_INTRON_2,0) = (-1);   
      GenePhase6_DC_SHADOW_SPECIAL(mat,0,j,AFTER_MATCH_CODING) = NEGI;   
      GenePhase6_DC_SHADOW_SPECIAL_SP(mat,0,j,AFTER_MATCH_CODING,0) = (-1);  
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_GenePhase6(pal)
 *
 * Descrip:    Converts a path alignment to a label alignment
 *             The label alignment is probably much more useful than the path
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * convert_PackAln_to_AlnBlock_GenePhase6(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_GenePhase6();    
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "MATCH_STATE","INSERT_STATE","DELETE_STATE","INTRON_STATE","LOOP_STATE","END","BEFORE_MATCH","BEFORE_MATCH_INTRON","AFTER_MATCH" }; 
/* Function:  AlnConvertSet_GenePhase6(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "CODON","3SS_PHASE_0","3SS_PHASE_1","3SS_PHASE_2","SEQUENCE_DELETION","SEQUENCE_INSERTION","INSERT","5SS_PHASE_0","CENTRAL_INTRON","5SS_PHASE_1","5SS_PHASE_2","RANDOM_SEQUENCE","END" };  
AlnConvertSet * AlnConvertSet_GenePhase6(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_CODING + 6;   
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_0;  
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 6;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_1;  
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_2;  
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 5;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = LOOP + 6;  
    acu->is_from_special = TRUE; 
    acu->state2 = INSERT;    
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_0;  
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 6;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_1;  
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 5;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_2;  
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 2;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[6];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INTRON_0;  
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INTRON_0;  
    acu->offi = 0;   
    acu->offj = 8;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_0;  
    acu->state2 = INTRON_0;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INTRON_1;  
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INTRON_1;  
    acu->offi = 0;   
    acu->offj = 9;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_1;  
    acu->state2 = INTRON_1;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INTRON_2;  
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[3];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INTRON_2;  
    acu->offi = 0;   
    acu->offj = 10;  
    acu->label1 = query_label[3];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INTRON_2;  
    acu->state2 = INTRON_2;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->state2 = LOOP + 6;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = LOOP + 6;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = LOOP + 6;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = LOOP + 6;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = LOOP + 6;  
    acu->state2 = LOOP + 6;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[11];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = LOOP + 6;  
    acu->state2 = END + 6;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_MATCH_CODING + 6;    
    acu->state2 = END + 6;   
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[12];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 6; 
    acu->state2 = BEFORE_MATCH_CODING + 6;   
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_CODING + 6;   
    acu->state2 = BEFORE_MATCH_CODING + 6;   
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_INTRON_0 + 6; 
    acu->state2 = BEFORE_MATCH_CODING + 6;   
    acu->offi = (-1);    
    acu->offj = 6;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_INTRON_1 + 6; 
    acu->state2 = BEFORE_MATCH_CODING + 6;   
    acu->offi = (-1);    
    acu->offj = 5;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_INTRON_2 + 6; 
    acu->state2 = BEFORE_MATCH_CODING + 6;   
    acu->offi = (-1);    
    acu->offj = 4;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_CODING + 6;   
    acu->state2 = BEFORE_MATCH_INTRON_0 + 6;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[7];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_INTRON_0 + 6; 
    acu->state2 = BEFORE_MATCH_INTRON_0 + 6;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_CODING + 6;   
    acu->state2 = BEFORE_MATCH_INTRON_1 + 6;     
    acu->offi = (-1);    
    acu->offj = 9;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[9];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_INTRON_1 + 6; 
    acu->state2 = BEFORE_MATCH_INTRON_1 + 6;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_CODING + 6;   
    acu->state2 = BEFORE_MATCH_INTRON_2 + 6;     
    acu->offi = (-1);    
    acu->offj = 10;  
    acu->label1 = query_label[7];    
    acu->label2 = target_label[10];  
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BEFORE_MATCH_INTRON_2 + 6; 
    acu->state2 = BEFORE_MATCH_INTRON_2 + 6;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[7];    
    acu->label2 = target_label[8];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = AFTER_MATCH_CODING + 6;    
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = AFTER_MATCH_CODING + 6;    
    acu->state2 = AFTER_MATCH_CODING + 6;    
    acu->offi = (-1);    
    acu->offj = 3;   
    acu->label1 = query_label[8];    
    acu->label2 = target_label[0];   
    add_collapse_label_AlnConvertSet(out,"INTRON_STATE","CENTRAL_INTRON");   
    add_collapse_label_AlnConvertSet(out,"BEFORE_MATCH_INTRON","CENTRAL_INTRON");    
    add_collapse_label_AlnConvertSet(out,"AFTER_MATCH_INTRON","CENTRAL_INTRON"); 
    add_collapse_label_AlnConvertSet(out,"LOOP_STATE","RANDOM_SEQUENCE");    
    return out;  
}    


/* Function:  PackAln_read_Expl_GenePhase6(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_GenePhase6(GenePhase6 * mat) 
{
    GenePhase6_access_func_holder holder;    


    holder.access_main    = GenePhase6_explicit_access_main; 
    holder.access_special = GenePhase6_explicit_access_special;  
    return PackAln_read_generic_GenePhase6(mat,holder);  
}    


/* Function:  GenePhase6_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenePhase6_explicit_access_main(GenePhase6 * mat,int i,int j,int state) 
{
    return GenePhase6_EXPL_MATRIX(mat,i,j,state);    
}    


/* Function:  GenePhase6_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int GenePhase6_explicit_access_special(GenePhase6 * mat,int i,int j,int state) 
{
    return GenePhase6_EXPL_SPECIAL(mat,i,j,state);   
}    


/* Function:  PackAln_read_generic_GenePhase6(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:          h [UNKN ] Undocumented argument [GenePhase6_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_GenePhase6(GenePhase6 * mat,GenePhase6_access_func_holder h) 
{
    register PackAln * out;  
    int i;   
    int j;   
    int state;   
    int cellscore = (-1);    
    boolean isspecial;   
    PackAlnUnit * pau = NULL;    
    PackAlnUnit * prev = NULL;   


    assert(mat);     
    assert(h.access_main);   
    assert(h.access_special);    


    out = PackAln_alloc_std();   
    if( out == NULL )    
      return NULL;   


    out->score =  find_end_GenePhase6(mat,&i,&j,&state,&isspecial,h);    


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in GenePhase6_basic_read, returning a mess.(Sorry!)",out->len);   
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 6;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_GenePhase6(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      else   
        max_calc_GenePhase6(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      if(i == GenePhase6_READ_OFF_ERROR || j == GenePhase6_READ_OFF_ERROR || state == GenePhase6_READ_OFF_ERROR )    {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in GenePhase6_basic_read, returning partial alignment",out->len);   
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 6;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_GenePhase6(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [GenePhase6_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_GenePhase6(GenePhase6 * mat,int * ri,int * rj,int * state,boolean * isspecial,GenePhase6_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->target->seq->len-1,END);    
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( (temp =(*h.access_special)(mat,0,j,END)) > max )   {  
        max = temp;  
        maxj = j;    
        }  
      }  


    if( ri != NULL)  
       *ri = 0;  
    if( rj != NULL)  
       *rj = maxj;   
    if( state != NULL)   
       *state = END; 
    if( isspecial != NULL)   
       *isspecial = TRUE;    


    return max;  
}    


/* Function:  GenePhase6_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void GenePhase6_debug_show_matrix(GenePhase6 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->len;i++) {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",GenePhase6_EXPL_MATRIX(mat,i,j,MATCH));   
        fprintf(ofp,"State INSERT %d\n",GenePhase6_EXPL_MATRIX(mat,i,j,INSERT)); 
        fprintf(ofp,"State DELETE %d\n",GenePhase6_EXPL_MATRIX(mat,i,j,DELETE)); 
        fprintf(ofp,"State INTRON_0 %d\n",GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0)); 
        fprintf(ofp,"State INTRON_1 %d\n",GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1)); 
        fprintf(ofp,"State INTRON_2 %d\n",GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2)); 
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_GenePhase6(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [GenePhase6_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_GenePhase6(GenePhase6 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GenePhase6_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GenePhase6_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->seq->len)  {  
      warn("In GenePhase6 matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (mat->gp->transition[GP4_INSERT_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 1,j - 5,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 5,MATCH);    
          }  
        temp = cscore - (mat->gp->transition[GP4_INSERT_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 1,j - 4,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 4,MATCH);    
          }  
        temp = cscore - (mat->gp->transition[GP4_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH);    
          }  
        temp = cscore - (mat->gp->transition[GP4_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 1,j - 2,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 2; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-2,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 2,MATCH);    
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-1))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 1,j - 4,INTRON_2) ) {  
          *reti = i - 1; 
          *retj = j - 4; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-4,INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 1,j - 4,INTRON_2);     
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-2))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 1,j - 5,INTRON_1) ) {  
          *reti = i - 1; 
          *retj = j - 5; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-5,INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 1,j - 5,INTRON_1);     
          }  
        temp = cscore - ((mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 1,j - 6,INTRON_0) ) {  
          *reti = i - 1; 
          *retj = j - 6; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-6,INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 1,j - 6,INTRON_0);     
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_special)(mat,i - 1,j - 3,BEFORE_MATCH_CODING) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-3,BEFORE_MATCH_CODING);  
            }  
          return (*h.access_main)(mat,i - 1,j - 3,BEFORE_MATCH_CODING);  
          }  
        temp = cscore - (((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == (*h.access_special)(mat,i - 1,j - 3,START) ) {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-3,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 3,START);    
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 3,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-3,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));  
        if( temp == (*h.access_main)(mat,i - 1,j - 3,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-3,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 1,j - 3,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-3,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->gp->transition[GP4_DELETE_2_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INSERT);   
          }  
        temp = cscore - (mat->gp->transition[GP4_DELETE_1_BASE]) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 0,j - 2,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-2,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 2,INSERT);   
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-1))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 4,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 4,INTRON_2);     
          }  
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-2))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 5,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-5,INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 5,INTRON_1);     
          }  
        temp = cscore - ((mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)))) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 6,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-6,INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 6,INTRON_0);     
          }  
        temp = cscore - (((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,LOOP); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,LOOP);     
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 0,j - 3,DELETE) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,DELETE);  
            }  
          return (*h.access_main)(mat,i - 0,j - 3,DELETE);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j));    
        if( temp == (*h.access_main)(mat,i - 0,j - 3,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 3,INSERT);   
          }  
        temp = cscore - ((mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])) -  (CSEQ_GENOMIC_CDSPOT(mat->target,j)); 
        if( temp == (*h.access_main)(mat,i - 0,j - 3,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 3,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->query->gws->seg[i]->transition[GW_DELETE2DELETE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 0,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,DELETE);   
          }  
        temp = cscore - (mat->query->gws->seg[i]->transition[GW_INSERT2DELETE]) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 0,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,INSERT);   
          }  
        temp = cscore - (mat->query->gws->seg[i]->transition[GW_MATCH2DELETE]) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_0 :    
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_0;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INTRON_0);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTRON_0);     
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 8,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 8,INSERT);   
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 8,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 8; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-8,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 8,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_1 :    
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_1;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INTRON_1);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTRON_1);     
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 9,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 9,INSERT);   
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 9,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-9,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 9,MATCH);    
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case INTRON_2 :    
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INTRON_2;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INTRON_2);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INTRON_2);     
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 10,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,INSERT); 
            }  
          return (*h.access_main)(mat,i - 0,j - 10,INSERT);  
          }  
        temp = cscore - (((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 10,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-10,MATCH);  
            }  
          return (*h.access_main)(mat,i - 0,j - 10,MATCH);   
          }  
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_GenePhase6(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [GenePhase6 *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [GenePhase6_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_GenePhase6(GenePhase6 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,GenePhase6_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = GenePhase6_READ_OFF_ERROR;   


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In GenePhase6 matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case LOOP :    
        /* source LOOP is a special */ 
        temp = cscore - (mat->gp->transition[GP4_LOOP2LOOP]) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,LOOP);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,LOOP) ;     
          }  
        /* source DELETE is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->gws->seg[i]->transition[GW_DELETE2END]) - (0);    
          if( temp == (*h.access_main)(mat,i - 0,j - 0,DELETE) ) {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = DELETE;  
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,DELETE);    
              }  
            return (*h.access_main)(mat,i - 0,j - 0,DELETE) ;    
            }  
          } /* end of for i >= 0 */ 
        /* source INSERT is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->gws->seg[i]->transition[GW_INSERT2END]) - (0);    
          if( temp == (*h.access_main)(mat,i - 0,j - 0,INSERT) ) {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = INSERT;  
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,INSERT);    
              }  
            return (*h.access_main)(mat,i - 0,j - 0,INSERT) ;    
            }  
          } /* end of for i >= 0 */ 
        /* source MATCH is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->gws->seg[i]->transition[GW_MATCH2END]) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH) ;     
            }  
          } /* end of for i >= 0 */ 
        /* source START is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,START) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,START);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,START) ;    
          }  
      case START :   
      case END :     
        /* source AFTER_MATCH_CODING is a special */ 
        temp = cscore - (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 3,AFTER_MATCH_CODING) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_MATCH_CODING;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,AFTER_MATCH_CODING);   
            }  
          return (*h.access_special)(mat,i - 0,j - 3,AFTER_MATCH_CODING) ;   
          }  
        /* source LOOP is a special */ 
        temp = cscore - (0) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,LOOP);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,LOOP) ;     
          }  
      case BEFORE_MATCH_CODING :     
        /* source BEFORE_MATCH_INTRON_2 is a special */ 
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-1))) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 4,BEFORE_MATCH_INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = BEFORE_MATCH_INTRON_2; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-4,BEFORE_MATCH_INTRON_2);    
            }  
          return (*h.access_special)(mat,i - 0,j - 4,BEFORE_MATCH_INTRON_2) ;    
          }  
        /* source BEFORE_MATCH_INTRON_1 is a special */ 
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-2))) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 5,BEFORE_MATCH_INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 5; 
          *retstate = BEFORE_MATCH_INTRON_1; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-5,BEFORE_MATCH_INTRON_1);    
            }  
          return (*h.access_special)(mat,i - 0,j - 5,BEFORE_MATCH_INTRON_1) ;    
          }  
        /* source BEFORE_MATCH_INTRON_0 is a special */ 
        temp = cscore - (CSEQ_GENOMIC_3SS(mat->target,(j-3))) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 6,BEFORE_MATCH_INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 6; 
          *retstate = BEFORE_MATCH_INTRON_0; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-6,BEFORE_MATCH_INTRON_0);    
            }  
          return (*h.access_special)(mat,i - 0,j - 6,BEFORE_MATCH_INTRON_0) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,BEFORE_MATCH_CODING) )   {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,BEFORE_MATCH_CODING);  
            }  
          return (*h.access_special)(mat,i - 0,j - 3,BEFORE_MATCH_CODING) ;  
          }  
        /* source START is a special */ 
        temp = cscore - (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 3,START) ) {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,START);    
            }  
          return (*h.access_special)(mat,i - 0,j - 3,START) ;    
          }  
      case BEFORE_MATCH_INTRON_0 :   
        /* source BEFORE_MATCH_INTRON_0 is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_0) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_INTRON_0; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_MATCH_INTRON_0);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_0) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_CODING) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_MATCH_CODING);  
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_CODING) ;  
          }  
      case BEFORE_MATCH_INTRON_1 :   
        /* source BEFORE_MATCH_INTRON_1 is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_1) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_INTRON_1; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_MATCH_INTRON_1);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_1) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 9,BEFORE_MATCH_CODING) )   {  
          *reti = i - 0; 
          *retj = j - 9; 
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-9,BEFORE_MATCH_CODING);  
            }  
          return (*h.access_special)(mat,i - 0,j - 9,BEFORE_MATCH_CODING) ;  
          }  
      case BEFORE_MATCH_INTRON_2 :   
        /* source BEFORE_MATCH_INTRON_2 is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) - (0);  
        if( temp == (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_2) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BEFORE_MATCH_INTRON_2; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,BEFORE_MATCH_INTRON_2);    
            }  
          return (*h.access_special)(mat,i - 0,j - 1,BEFORE_MATCH_INTRON_2) ;    
          }  
        /* source BEFORE_MATCH_CODING is a special */ 
        temp = cscore - ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) - (0);     
        if( temp == (*h.access_special)(mat,i - 0,j - 10,BEFORE_MATCH_CODING) )  {  
          *reti = i - 0; 
          *retj = j - 10;    
          *retstate = BEFORE_MATCH_CODING;   
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-10,BEFORE_MATCH_CODING);     
            }  
          return (*h.access_special)(mat,i - 0,j - 10,BEFORE_MATCH_CODING) ;     
          }  
      case AFTER_MATCH_CODING :  
        /* source AFTER_MATCH_CODING is a special */ 
        temp = cscore - (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 3,AFTER_MATCH_CODING) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = AFTER_MATCH_CODING;    
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-3,AFTER_MATCH_CODING);   
            }  
          return (*h.access_special)(mat,i - 0,j - 3,AFTER_MATCH_CODING) ;   
          }  
        /* source MATCH is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (mat->query->gws->seg[i]->transition[GW_MATCH2END]) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH) ;     
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in GenePhase6 read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_GenePhase6(mat)
 *
 * Descrip:    This function calculates the GenePhase6 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GenePhase6
 *
 *
 * Arg:        mat [UNKN ] GenePhase6 which contains explicit basematrix memory [GenePhase6 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_GenePhase6(GenePhase6 * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_GenePhase6, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("GenePhase6 Matrix calculation: ");  
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = GenePhase6_EXPL_SPECIAL(mat,i-1,j-3,START) + ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_MATCH_CODING to state MATCH */ 
        temp = GenePhase6_EXPL_SPECIAL(mat,i-1,j-3,BEFORE_MATCH_CODING) + (mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_EXPL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) )   {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* state MATCH is a source for special AFTER_MATCH_CODING */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) )     {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = GenePhase6_EXPL_SPECIAL(mat,i-0,j-3,LOOP) + ((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_EXPL_MATRIX(mat,i,j,INSERT) = score; 


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_INSERT2END]) + (0) ;  
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) )   {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];    
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,DELETE) = score; 


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_DELETE2END]) + (0) ;  
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) )   {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);    
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0) = score;   


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);    
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1) = score;   


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);   
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2) = score;   


        /* Finished calculating state INTRON_2 */ 
        }  


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,LOOP); 


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (mat->gp->transition[GP4_LOOP2LOOP]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,LOOP) = score; 
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,END);  


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Source AFTER_MATCH_CODING is a special source for END */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 




      /* Special state BEFORE_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_CODING);  


      /* Source START is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,START) + (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,BEFORE_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 6,BEFORE_MATCH_INTRON_0) + (CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 5,BEFORE_MATCH_INTRON_1) + (CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 4,BEFORE_MATCH_INTRON_2) + (CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_CODING) = score;  
      /* Finished updating state BEFORE_MATCH_CODING */ 




      /* Special state BEFORE_MATCH_INTRON_0 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0);    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_0) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0) = score;    
      /* Finished updating state BEFORE_MATCH_INTRON_0 */ 




      /* Special state BEFORE_MATCH_INTRON_1 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1);    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 9,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_1) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1) = score;    
      /* Finished updating state BEFORE_MATCH_INTRON_1 */ 




      /* Special state BEFORE_MATCH_INTRON_2 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2);    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 10,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_2) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2) = score;    
      /* Finished updating state BEFORE_MATCH_INTRON_2 */ 




      /* Special state AFTER_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,AFTER_MATCH_CODING);   


      /* Source MATCH for state AFTER_MATCH_CODING is not special... already calculated */ 
      /* Source AFTER_MATCH_CODING is a special source for AFTER_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,AFTER_MATCH_CODING) = score;   
      /* Finished updating state AFTER_MATCH_CODING */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_GenePhase6(mat,dpenv)
 *
 * Descrip:    This function calculates the GenePhase6 matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] GenePhase6 which contains explicit basematrix memory [GenePhase6 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_GenePhase6(GenePhase6 * mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int starti;  
    int startj;  
    int endi;    
    int endj;    
    int tot; 
    int num; 
    int should_calc; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_GenePhase6, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    prepare_DPEnvelope(dpenv);   
    starti = dpenv->starti;  
    if( starti < 0 ) 
      starti = 0;    
    startj = dpenv->startj;  
    if( startj < 0 ) 
      startj = 0;    
    endi = dpenv->endi;  
    if( endi > mat->leni )   
      endi = mat->leni;  
    endj = dpenv->endj;  
    if( endj > mat->lenj )   
      endj = mat->lenj;  
    tot = (endi-starti) * (endj-startj); 
    num = 0; 


    for(j=startj-10;j<endj;j++)  {  
      for(i=1;i<mat->leni;i++)   {  
        GenePhase6_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
        GenePhase6_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;   
        GenePhase6_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;   
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI; 
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI; 
        GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI; 
        }  
      }  
    for(j=-10;j<mat->lenj;j++)   {  
      GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = NEGI;  
      GenePhase6_EXPL_SPECIAL(mat,i,j,START) = 0;    
      GenePhase6_EXPL_SPECIAL(mat,i,j,END) = NEGI;   
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_CODING) = NEGI;   
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_0) = NEGI; 
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_1) = NEGI; 
      GenePhase6_EXPL_SPECIAL(mat,i,j,BEFORE_MATCH_INTRON_2) = NEGI; 
      GenePhase6_EXPL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = NEGI;    
      }  


    start_reporting("GenePhase6 Matrix calculation: ");  
    for(j=startj;j<endj;j++) {  
      auto int score;    
      auto int temp;     
      for(i=starti;i<endi;i++)   {  
        /* Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed */ 
        should_calc = 0; 
        for(k=0;k<dpenv->len;k++)    {  
          auto DPUnit * u;   
          u = dpenv->dpu[k]; 
          switch(u->type)    {  
            case DPENV_RECT :    
              if( i >= u->starti && j >= u->startj && i <= (u->starti+u->height) && j <= (u->startj+u->length))  
                should_calc = 1;     
              break; 
            case DPENV_DIAG :    
              if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)   
                should_calc = 1;     
              break; 
            }  
          if( should_calc == 1 ) 
            break;   
          }  
        if( should_calc == 0)    {  
          GenePhase6_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
          GenePhase6_EXPL_MATRIX(mat,i,j,INSERT) = NEGI; 
          GenePhase6_EXPL_MATRIX(mat,i,j,DELETE) = NEGI; 
          GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0) = NEGI;   
          GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1) = NEGI;   
          GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-1,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);     
        /* From state INSERT to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = GenePhase6_EXPL_SPECIAL(mat,i-1,j-3,START) + ((mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BEFORE_MATCH_CODING to state MATCH */ 
        temp = GenePhase6_EXPL_SPECIAL(mat,i-1,j-3,BEFORE_MATCH_CODING) + (mat->query->gws->seg[i]->transition[GW_START2MATCH]+mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-2,MATCH) + mat->gp->transition[GP4_DELETE_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-1,MATCH) + mat->gp->transition[GP4_DELETE_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-4,MATCH) + mat->gp->transition[GP4_INSERT_1_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-5,MATCH) + mat->gp->transition[GP4_INSERT_2_BASE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_EXPL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) )   {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* state MATCH is a source for special AFTER_MATCH_CODING */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_MATCH2END]) + (0) ;   
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) )     {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,AFTER_MATCH_CODING) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-3,MATCH) + (mat->query->gws->seg[i]->transition[GW_MATCH2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);   
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-3,INSERT) + (mat->query->gws->seg[i]->transition[GW_INSERT2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-3,DELETE) + (mat->query->gws->seg[i]->transition[GW_DELETE2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state LOOP to state INSERT */ 
        temp = GenePhase6_EXPL_SPECIAL(mat,i-0,j-3,LOOP) + ((mat->query->gws->seg[i]->transition[GW_START2INSERT]+mat->query->gws->seg[i]->insert[CSEQ_GENOMIC_CODON(mat->target,j)])+mat->gp->transition[GP4_LOOP2MODEL]);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-6,INTRON_0) + (mat->query->gws->seg[i]->match[CSEQ_GENOMIC_CODON(mat->target,j)]+CSEQ_GENOMIC_3SS(mat->target,(j-3)));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-5,INTRON_1) + CSEQ_GENOMIC_3SS(mat->target,(j-2));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-4,INTRON_2) + CSEQ_GENOMIC_3SS(mat->target,(j-1));   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-2,INSERT) + mat->gp->transition[GP4_DELETE_1_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INSERT to state INSERT */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->gp->transition[GP4_DELETE_2_BASE];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         score += CSEQ_GENOMIC_CDSPOT(mat->target,j);    
         GenePhase6_EXPL_MATRIX(mat,i,j,INSERT) = score; 


        /* state INSERT is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_INSERT2END]) + (0) ;  
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) )   {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->query->gws->seg[i]->transition[GW_MATCH2DELETE];    
        /* From state INSERT to state DELETE */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-0,INSERT) + mat->query->gws->seg[i]->transition[GW_INSERT2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state DELETE */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->query->gws->seg[i]->transition[GW_DELETE2DELETE];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,DELETE) = score; 


        /* state DELETE is a source for special LOOP */ 
        temp = score + (mat->query->gws->seg[i]->transition[GW_DELETE2END]) + (0) ;  
        if( temp > GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) )   {  
          GenePhase6_EXPL_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* Finished calculating state DELETE */ 


        /* For state INTRON_0 */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-8,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_0);    
        /* From state INSERT to state INTRON_0 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-8,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_0 to state INTRON_0 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INTRON_0) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_0 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_0) = score;   


        /* Finished calculating state INTRON_0 */ 


        /* For state INTRON_1 */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-9,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_1);    
        /* From state INSERT to state INTRON_1 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-9,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_1 to state INTRON_1 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INTRON_1) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_1 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_1) = score;   


        /* Finished calculating state INTRON_1 */ 


        /* For state INTRON_2 */ 
        /* setting first movement to score */ 
        score = GenePhase6_EXPL_MATRIX(mat,i-0,j-10,MATCH) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->intron_2);   
        /* From state INSERT to state INTRON_2 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-10,INSERT) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))+mat->query->phase[i]->insert_intron);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state INTRON_2 to state INTRON_2 */ 
        temp = GenePhase6_EXPL_MATRIX(mat,i-0,j-1,INTRON_2) + (mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON]);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INTRON_2 */ 
        /* Add any movement independant score and put away */ 
         GenePhase6_EXPL_MATRIX(mat,i,j,INTRON_2) = score;   


        /* Finished calculating state INTRON_2 */ 
        }  


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,LOOP); 


      /* Source START is a special source for LOOP */ 
      /* Has restricted position */ 
      if( (j-1) == 0  )  {  
        temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,START) + (0) + (0);   
        if( temp > score )   
          score = temp;  
        }  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source INSERT for state LOOP is not special... already calculated */ 
      /* Source DELETE for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (mat->gp->transition[GP4_LOOP2LOOP]) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,LOOP) = score; 
      /* Finished updating state LOOP */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,END);  


      /* Source LOOP is a special source for END */ 
      /* Has restricted position */ 
      if( j == mat->lenj-1 ) {  
        temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (0) + (0);    
        if( temp > score )   
          score = temp;  
        }  


      /* Source AFTER_MATCH_CODING is a special source for END */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->stop->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,END) = score;  
      /* Finished updating state END */ 




      /* Special state BEFORE_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_CODING);  


      /* Source START is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,START) + (mat->general_model->start->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,BEFORE_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 6,BEFORE_MATCH_INTRON_0) + (CSEQ_GENOMIC_3SS(mat->target,(j-3))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 5,BEFORE_MATCH_INTRON_1) + (CSEQ_GENOMIC_3SS(mat->target,(j-2))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 4,BEFORE_MATCH_INTRON_2) + (CSEQ_GENOMIC_3SS(mat->target,(j-1))) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_CODING) = score;  
      /* Finished updating state BEFORE_MATCH_CODING */ 




      /* Special state BEFORE_MATCH_INTRON_0 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0);    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,j))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_0 is a special source for BEFORE_MATCH_INTRON_0 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_0) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_0) = score;    
      /* Finished updating state BEFORE_MATCH_INTRON_0 */ 




      /* Special state BEFORE_MATCH_INTRON_1 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1);    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 9,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_1 is a special source for BEFORE_MATCH_INTRON_1 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_1) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_1) = score;    
      /* Finished updating state BEFORE_MATCH_INTRON_1 */ 




      /* Special state BEFORE_MATCH_INTRON_2 has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2);    


      /* Source BEFORE_MATCH_CODING is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 10,BEFORE_MATCH_CODING) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+CSEQ_GENOMIC_5SS(mat->target,(j-7)))) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Source BEFORE_MATCH_INTRON_2 is a special source for BEFORE_MATCH_INTRON_2 */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 1,BEFORE_MATCH_INTRON_2) + ((mat->gp->intron[CSEQ_GENOMIC_BASE(mat->target,j)]+mat->gp->transition[GP4_INTRON2INTRON])) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,BEFORE_MATCH_INTRON_2) = score;    
      /* Finished updating state BEFORE_MATCH_INTRON_2 */ 




      /* Special state AFTER_MATCH_CODING has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = GenePhase6_EXPL_SPECIAL(mat,0,j,AFTER_MATCH_CODING);   


      /* Source MATCH for state AFTER_MATCH_CODING is not special... already calculated */ 
      /* Source AFTER_MATCH_CODING is a special source for AFTER_MATCH_CODING */ 
      temp = GenePhase6_EXPL_SPECIAL(mat,0,j - 3,AFTER_MATCH_CODING) + (mat->general_model->general->codon[CSEQ_GENOMIC_CODON(mat->target,j)]) + (0);    
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      GenePhase6_EXPL_SPECIAL(mat,0,j,AFTER_MATCH_CODING) = score;   
      /* Finished updating state AFTER_MATCH_CODING */ 


      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  GenePhase6_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenePhase6 *]
 *
 */
GenePhase6 * GenePhase6_alloc(void) 
{
    GenePhase6 * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenePhase6 *) ckalloc (sizeof(GenePhase6))) == NULL)    {  
      warn("GenePhase6_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->basematrix = NULL;  
    out->shatter = NULL; 
    out->leni = 0;   
    out->lenj = 0;   


    return out;  
}    


/* Function:  free_GenePhase6(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenePhase6 *]
 *
 * Return [UNKN ]  Undocumented return value [GenePhase6 *]
 *
 */
GenePhase6 * free_GenePhase6(GenePhase6 * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenePhase6 obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   
    if( obj->basematrix != NULL) 
      free_BaseMatrix(obj->basematrix);  
    if( obj->shatter != NULL)    
      free_ShatterMatrix(obj->shatter);  
    /* obj->query is linked in */ 
    /* obj->target is linked in */ 
    /* obj->gp is linked in */ 
    /* obj->general_model is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
