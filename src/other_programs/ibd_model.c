#ifdef _cplusplus
extern "C" {
#endif
#include "ibd_model.h"


# line 82 "ibd_model.dy"
void dump_naive_AlnBlock_IBD(char * strain_a,char * strain_b,GenoVarChr * chr,AlnBlock * alb,FILE * ofp)
{
  AlnColumn * bstart;
  AlnColumn * e;
  int j;
  int jj;

  e = bstart = alb->start;
  j = 0;
  jj = 0;

  while( e != NULL ) {
    e = e->next;
    jj++;

    if( strcmp(e->alu[1]->text_label,bstart->alu[1]->text_label) == 0 ) {
      continue;
    }



    /* break */

    fprintf(ofp,"%s\t%ld\t%ld\t%d\t%d\t%s\t%s\t%s\n",
	    chr->chr,
	    chr->loci[j]->var->pos,
	    chr->loci[jj]->var->pos,
	    j,
	    jj,
	    strain_a,
	    strain_b,
	    bstart->alu[1]->text_label);

    if( strcmp(e->alu[0]->text_label,"END") == 0 ) {
      break;
    }
    
    bstart = e;
    j = jj;
  }
  
  


}

# line 128 "ibd_model.dy"
PairMatch * new_PairMatch(GenoVarSet * set,GenoVarChr * chr,char * strain_a,char * strain_b)
{
  PairMatch * out;
  int index_a;
  int index_b;

  int i;

  index_a = individual_index_from_string_GenoVarSet(set,strain_a);
  assert(index_a != -1);

  index_b = individual_index_from_string_GenoVarSet(set,strain_b);
  assert(index_b != -1);

  

  out = PairMatch_alloc();


  out->match = calloc(chr->len,sizeof(char));
  out->len = chr->len;

  for(i=0;i<chr->len;i++) {
    if( chr->loci[i]->ind[index_a] == chr->loci[i]->ind[index_b] ) {
      out->match[i] = 1;
    } else {
      out->match[i] = 0;
    }
  }
  
  return(out);

}


# line 163 "ibd_model.dy"
GenomePairPara * new_GenomePairPara_null(Probability ibd_mismatch,Probability diff_mismatch,Probability ibd_switch,Probability diff_switch)
{
  GenomePairPara * out;


  out = GenomePairPara_alloc();

  out->ibd_mismatch = Probability2Score(1);
  out->ibd_match = Probability2Score(1);
  out->diff_mismatch = Probability2Score(diff_mismatch/ ibd_mismatch);
  out->diff_match = Probability2Score((1.0 - diff_mismatch)/(1.0 - ibd_mismatch));

  out->bswitch_ibd = Probability2Score(ibd_switch);
  out->extend_ibd  = Probability2Score(1.0 - ibd_switch);


  out->bswitch_diff = Probability2Score(diff_switch);
  out->extend_diff  = Probability2Score(1.0 - diff_switch);
  out->len = 1;
  return(out);

}

# line 186 "ibd_model.dy"
GenomePairPara * new_GenomePairPara(Probability ibd_mismatch,Probability diff_mismatch,Probability ibd_switch,Probability diff_switch)
{
  GenomePairPara * out;


  out = GenomePairPara_alloc();

  out->ibd_mismatch = Probability2Score(ibd_mismatch);
  out->ibd_match = Probability2Score(1.0 - ibd_mismatch);
  out->diff_mismatch = Probability2Score(diff_mismatch);
  out->diff_match = Probability2Score(1.0 - diff_mismatch);

  out->bswitch_ibd = Probability2Score(ibd_switch);
  out->extend_ibd  = Probability2Score(1.0 - ibd_switch);


  out->bswitch_diff = Probability2Score(diff_switch);
  out->extend_diff  = Probability2Score(1.0 - diff_switch);
  out->len = 1;
  return(out);

}

/* Function:  hard_link_MousePairMatch_Posterior(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MousePairMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch_Posterior *]
 *
 */
MousePairMatch_Posterior * hard_link_MousePairMatch_Posterior(MousePairMatch_Posterior * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MousePairMatch_Posterior object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MousePairMatch_Posterior_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch_Posterior *]
 *
 */
MousePairMatch_Posterior * MousePairMatch_Posterior_alloc(void) 
{
    MousePairMatch_Posterior * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MousePairMatch_Posterior *) ckalloc (sizeof(MousePairMatch_Posterior))) == NULL)    {  
      warn("MousePairMatch_Posterior_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->forward = NULL; 
    out->backward = NULL;    


    return out;  
}    


/* Function:  free_MousePairMatch_Posterior(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MousePairMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch_Posterior *]
 *
 */
MousePairMatch_Posterior * free_MousePairMatch_Posterior(MousePairMatch_Posterior * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MousePairMatch_Posterior obj. Should be trappable");  
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
    if( obj->forward != NULL)    
      free_MousePairMatch(obj->forward);     
    if( obj->backward != NULL)   
      free_MousePairMatch(obj->backward);    


    ckfree(obj); 
    return NULL; 
}    


# line 229 "ibd_model.c"
/* Function:  hard_link_PairMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PairMatch *]
 *
 */
PairMatch * hard_link_PairMatch(PairMatch * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairMatch object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairMatch *]
 *
 */
PairMatch * PairMatch_alloc(void) 
{
    PairMatch * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairMatch *) ckalloc (sizeof(PairMatch))) == NULL)  {  
      warn("PairMatch_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->len = 0;    
    out->match = NULL;   


    return out;  
}    


/* Function:  free_PairMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PairMatch *]
 *
 */
PairMatch * free_PairMatch(PairMatch * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairMatch obj. Should be trappable"); 
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
    if( obj->match != NULL)  
      ckfree(obj->match);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GenomePairPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomePairPara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePairPara *]
 *
 */
GenomePairPara * hard_link_GenomePairPara(GenomePairPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomePairPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomePairPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomePairPara *]
 *
 */
GenomePairPara * GenomePairPara_alloc(void) 
{
    GenomePairPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomePairPara *) ckalloc (sizeof(GenomePairPara))) == NULL)    {  
      warn("GenomePairPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->ibd_match = 0;  
    out->ibd_mismatch = 0;   
    out->diff_match = 0; 
    out->diff_mismatch = 0;  
    out->extend_ibd = 0; 
    out->bswitch_ibd = 0;    
    out->extend_diff = 0;    
    out->bswitch_diff = 0;   
    out->len = 0;    


    return out;  
}    


/* Function:  free_GenomePairPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomePairPara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePairPara *]
 *
 */
GenomePairPara * free_GenomePairPara(GenomePairPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomePairPara obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    




  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Fri Feb  4 18:45:52 2011           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define IBD 0    
#define DIFF 1   


#define Start 0  
#define End 1    


#define MousePairMatch_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*2)+STATE][i+0]  
#define MousePairMatch_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1] 
#define MousePairMatch_READ_OFF_ERROR -2
    


#define MousePairMatch_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+0)*2)+STATE]    
#define MousePairMatch_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE] 




#define MousePairMatch_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]  
#define MousePairMatch_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)  


/* Function:  PackAln_read_Shatter_MousePairMatch(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_MousePairMatch(MousePairMatch * mat) 
{
    MousePairMatch_access_func_holder holder;    


    holder.access_main    = MousePairMatch_shatter_access_main;  
    holder.access_special = MousePairMatch_shatter_access_special;   
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_MousePairMatch(mat,holder);  
}    


/* Function:  MousePairMatch_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MousePairMatch_shatter_access_main(MousePairMatch * mat,int i,int j,int state) 
{
    return MousePairMatch_SHATTER_MATRIX(mat,i,j,state); 
}    


/* Function:  MousePairMatch_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MousePairMatch_shatter_access_special(MousePairMatch * mat,int i,int j,int state) 
{
    return MousePairMatch_SHATTER_SPECIAL(mat,i,j,state);    
}    


/* Function:  calculate_shatter_MousePairMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MousePairMatch matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [MousePairMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_MousePairMatch(MousePairMatch * mat,DPEnvelope * dpenv) 
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
    int * SIG_0_1;   


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,2,lenj,2);    
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


    start_reporting("MousePairMatch Matrix calculation: ");  
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
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   




        /* For state IBD */ 
        /* setting first movement to score */ 
        score = SIG_0_1[IBD] + mat->query->extend_ibd;   
        /* From state DIFF to state IBD */ 
        temp = SIG_0_1[DIFF] + mat->query->bswitch_diff;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state IBD */ 
          temp = MousePairMatch_SHATTER_SPECIAL(mat,i-0,j-1,Start) + 0;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         SIG_0_0[IBD] = score;   


        /* state IBD is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_SHATTER_SPECIAL(mat,i,j,End) )   {  
            MousePairMatch_SHATTER_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = SIG_0_1[DIFF] + mat->query->extend_diff;     
        /* From state IBD to state DIFF */ 
        temp = SIG_0_1[IBD] + mat->query->bswitch_ibd;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state DIFF */ 
          temp = MousePairMatch_SHATTER_SPECIAL(mat,i-0,j-1,Start) + 0;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         SIG_0_0[DIFF] = score;  


        /* state DIFF is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_SHATTER_SPECIAL(mat,i,j,End) )   {  
            MousePairMatch_SHATTER_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state DIFF */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_MousePairMatch(dbsi,out,query,target)
 *
 * Descrip:    This function makes a database search of MousePairMatch
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePairPara*]
 * Arg:        target [UNKN ] Undocumented argument [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_MousePairMatch(DBSearchImpl * dbsi,Hscore * out,GenomePairPara* query,PairMatch* target ) 
{
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_MousePairMatch. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_MousePairMatch. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 0 )  
      warn("Although you are asking at run-time for database tracing, the MousePairMatch matrix was not compiled with database tracing. No tracing will be made");   
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_MousePairMatch(out,query,target );  
      case DBSearchImpl_Pthreads :   
        warn("This matrix MousePairMatch was not dyc compiled with thread support"); 
        return SEARCH_ERROR; 
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from MousePairMatch",impl_string_DBSearchImpl(dbsi)); 
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  score_only_logsum_MousePairMatch(query,target)
 *
 * Descrip:    This function calculates the score over all paths
 *             This is using a logsum method to sort it all out
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score score_only_logsum_MousePairMatch(GenomePairPara* query,PairMatch* target ) 
{
    int i;   
    int j;   
    int bestscore = 0;   
    int k;   
    MousePairMatch * mat;    


    mat = allocate_MousePairMatch_only(query, target );  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DISASTER!");   
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 0) * 2,2,2)) == NULL)  {  
      warn("Score only matrix for MousePairMatch cannot be allocated, (asking for 1  by %d  cells)",mat->leni*2);    
      mat = free_MousePairMatch(mat);    
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->leni;i++)    {  
        for(k=0;k<2;k++) 
          MousePairMatch_VSMALL_MATRIX(mat,i,j,k) = NEGI;    
        }  
      MousePairMatch_VSMALL_SPECIAL(mat,i,j,Start) = 0;  
      MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) = NEGI; 
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;  
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;    
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state IBD */ 
          temp = MousePairMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + 0;   
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_VSMALL_MATRIX(mat,i,j,IBD) = score;  


        /* state IBD is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) = Probability_logsum(MousePairMatch_VSMALL_SPECIAL(mat,i,j,End),temp);  
          }  


        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;    
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;  
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state DIFF */ 
          temp = MousePairMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + 0;   
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_VSMALL_MATRIX(mat,i,j,DIFF) = score; 


        /* state DIFF is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) = Probability_logsum(MousePairMatch_VSMALL_SPECIAL(mat,i,j,End),temp);  
          }  


        /* Finished calculating state DIFF */ 
        } /* end of for all query positions */ 




      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      } /* end of for all target positions */ 


    mat = free_MousePairMatch(mat);  
    return bestscore;    
}    


/* Function:  forward_logsum_MousePairMatch(query,target,dpri)
 *
 * Descrip:    This function calculates the matrix over all paths
 *             This is using a logsum method to sort it all out
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * forward_logsum_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPRunImpl * dpri) 
{
    MousePairMatch * mat;    
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    assert((mat=allocate_MousePairMatch_only(query, target )) != NULL);  
    mat->basematrix = BaseMatrix_alloc_matrix_specials_score_offset((mat->lenj+1)*2,(mat->leni+0),2,mat->lenj+1);    
    assert(mat->basematrix != NULL);     
    leni = mat->leni;    
    lenj = mat->lenj;    
    mat->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_MousePairMatch(mat);    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("MousePairMatch Matrix calculation: ");  
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;    
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;  
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state IBD */ 
          temp = MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + 0;     
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = score;    


        /* state IBD is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = Probability_logsum(MousePairMatch_EXPL_SPECIAL(mat,i,j,End),temp);  
          }  


        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;  
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;    
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state DIFF */ 
          temp = MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + 0;     
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = score;   


        /* state DIFF is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = Probability_logsum(MousePairMatch_EXPL_SPECIAL(mat,i,j,End),temp);  
          }  


        /* Finished calculating state DIFF */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return mat;  
}    


/* Function:  backward_logsum_MousePairMatch(query,target,dpri)
 *
 * Descrip:    This function calculates the matrix over all paths
 *             This is using a logsum method to sort it all out
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * backward_logsum_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPRunImpl * dpri) 
{
    MousePairMatch * mat;    
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 
    int max_score;   
    int min_score;   


    assert((mat=allocate_MousePairMatch_only(query, target )) != NULL);  
    mat->basematrix = BaseMatrix_alloc_matrix_specials_score_offset((mat->lenj+1)*2,(mat->leni+0),2,mat->lenj+1);    
    assert(mat->basematrix != NULL);     
    leni = mat->leni;    
    lenj = mat->lenj;    
    mat->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_MousePairMatch(mat);    
    tot = leni * lenj;   
    num = 0; 
    /* On the backward view, we need to have 0 at end, and NEGI at start */ 
    for(j= -1;j<lenj;j++)    {  
      MousePairMatch_EXPL_SPECIAL(mat,0,j,Start) = NEGI;     
      MousePairMatch_EXPL_SPECIAL(mat,0,j,End) = 0;  
      }  


    start_reporting("MousePairMatch Matrix calculation: ");  


    for(j=lenj-1;j>=0;j--)   {  
      auto int score;    
      auto int temp;     
      /* We look for underflow or overflow on this j column. */ 
      /* We do this at the start of the loop because in the backwards case at this point we know nothing will _update_ this column */ 
      min_score = SCORE_OVERFLOW;    
      max_score = SCORE_UNDERFLOW;   
      for(i=0;i<mat->leni;i++)   {  
        if( MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) > max_score) 
          max_score =  MousePairMatch_EXPL_MATRIX(mat,i,j,IBD);  
        if( MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) < min_score) 
          min_score =  MousePairMatch_EXPL_MATRIX(mat,i,j,IBD);  
        if( MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) > max_score)    
          max_score =  MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF);     
        if( MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) < min_score)    
          min_score =  MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF);     
        }  
      if( MousePairMatch_EXPL_SPECIAL(mat,i,j,Start) > max_score)    
        max_score =  MousePairMatch_EXPL_SPECIAL(mat,i,j,Start);     
      if( MousePairMatch_EXPL_SPECIAL(mat,i,j,Start) < min_score)    
        min_score =  MousePairMatch_EXPL_SPECIAL(mat,i,j,Start);     
      if( MousePairMatch_EXPL_SPECIAL(mat,i,j,End) > max_score)  
        max_score =  MousePairMatch_EXPL_SPECIAL(mat,i,j,End);   
      if( MousePairMatch_EXPL_SPECIAL(mat,i,j,End) < min_score)  
        min_score =  MousePairMatch_EXPL_SPECIAL(mat,i,j,End);   
      if( max_score > SCORE_OVERFLOW )   {  
        if( min_score < SCORE_UNDERFLOW) 
          fatal("Both overflow and underflow on the same column");   
        mat->basematrix->score_offset[j] = -1;   
        for(i=0;i<mat->leni;i++) {  
          MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) -= SCORE_OVERFLOW;     
          MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) -= SCORE_OVERFLOW;    
          }  
        MousePairMatch_EXPL_SPECIAL(mat,i,j,Start) -= SCORE_OVERFLOW;    
        MousePairMatch_EXPL_SPECIAL(mat,i,j,End) -= SCORE_OVERFLOW;  
        }  
      if( min_score < SCORE_UNDERFLOW )  {  
        mat->basematrix->score_offset[j] = +1;   
        for(i=0;i<mat->leni;i++) {  
          MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) -= SCORE_UNDERFLOW;    
          MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) -= SCORE_UNDERFLOW;   
          }  
        MousePairMatch_EXPL_SPECIAL(mat,i,j,Start) -= SCORE_UNDERFLOW;   
        MousePairMatch_EXPL_SPECIAL(mat,i,j,End) -= SCORE_UNDERFLOW;     
        }  
      for(i=leni-1;i>=0;i--) {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   




        /* special state End recieves from IBD */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = MousePairMatch_EXPL_SPECIAL(mat,i,j,End) + 0 + (0);     
          MousePairMatch_EXPL_MATRIX(mat,i-0,j-0,IBD) = Probability_logsum(MousePairMatch_EXPL_MATRIX(mat,i-0,j-0,IBD),temp);    
          }  
        /* special state End recieves from DIFF */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = MousePairMatch_EXPL_SPECIAL(mat,i,j,End) + 0 + (0);     
          MousePairMatch_EXPL_MATRIX(mat,i-0,j-0,DIFF) = Probability_logsum(MousePairMatch_EXPL_MATRIX(mat,i-0,j-0,DIFF),temp);  
          }  


        /* Doing calculations which end at this state */ 
        /* Reversed from state IBD to state IBD */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) + mat->query->extend_ibd + (TargetMatchIBD(j));   
        MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) = Probability_logsum(MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD),temp);  
        /* Reversed from state IBD to state DIFF */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) + mat->query->bswitch_diff + (TargetMatchIBD(j));     
        MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) = Probability_logsum(MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF),temp);    
        /* Has restricted position */ 
        if( (j-0) == 0  )    {  
          /* Reversed from state IBD to state Start */ 
          temp = MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) + 0 + (TargetMatchIBD(j));  
          MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) = Probability_logsum(MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start),temp);  
          }  


        /* Doing calculations which end at this state */ 
        /* Reversed from state DIFF to state DIFF */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) + mat->query->extend_diff + (TargetMatchDIFF(j));    
        MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) = Probability_logsum(MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF),temp);    
        /* Reversed from state DIFF to state IBD */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) + mat->query->bswitch_ibd + (TargetMatchDIFF(j));    
        MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) = Probability_logsum(MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD),temp);  
        /* Has restricted position */ 
        if( (j-0) == 0  )    {  
          /* Reversed from state DIFF to state Start */ 
          temp = MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) + 0 + (TargetMatchDIFF(j));    
          MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) = Probability_logsum(MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start),temp);  
          }  
        }  
      }  
    stop_reporting();    
    return mat;  
}    


/* Function:  serial_search_MousePairMatch(out,query,target)
 *
 * Descrip:    This function makes a database search of MousePairMatch
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePairPara*]
 * Arg:        target [UNKN ] Undocumented argument [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_MousePairMatch(Hscore * out,GenomePairPara* query,PairMatch* target ) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_MousePairMatch(query, target );   
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("MousePairMatch search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure"); 
        return SEARCH_ERROR; 
        }  
      /* Now: add query/target information to the entry */ 
      ds->score = score;     
      add_Hscore(out,ds);    
      } /* end of if storing datascore */ 
    pop_errormsg_stack();    
    push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);    


    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_MousePairMatch(query,target)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_MousePairMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_MousePairMatch(GenomePairPara* query,PairMatch* target ) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    MousePairMatch * mat;    


    mat = allocate_MousePairMatch_only(query, target );  
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 0) * 2,2,2)) == NULL)  {  
      warn("Score only matrix for MousePairMatch cannot be allocated, (asking for 1  by %d  cells)",mat->leni*2);    
      mat = free_MousePairMatch(mat);    
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->leni;i++)    {  
        for(k=0;k<2;k++) 
          MousePairMatch_VSMALL_MATRIX(mat,i,j,k) = NEGI;    
        }  
      MousePairMatch_VSMALL_SPECIAL(mat,i,j,Start) = 0;  
      MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) = NEGI; 
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;  
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state IBD */ 
          temp = MousePairMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + 0;   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_VSMALL_MATRIX(mat,i,j,IBD) = score;  


        /* state IBD is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) )    {  
            MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) = temp;   
            }  


          }  


        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;    
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_VSMALL_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state DIFF */ 
          temp = MousePairMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + 0;   
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_VSMALL_MATRIX(mat,i,j,DIFF) = score; 


        /* state DIFF is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) )    {  
            MousePairMatch_VSMALL_SPECIAL(mat,i,j,End) = temp;   
            }  


          }  


        /* Finished calculating state DIFF */ 
        } /* end of for all query positions */ 




      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      if( bestscore < MousePairMatch_VSMALL_SPECIAL(mat,0,j,End) )   
        bestscore = MousePairMatch_VSMALL_SPECIAL(mat,0,j,End);  
      } /* end of for all target positions */ 


    mat = free_MousePairMatch(mat);  
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_MousePairMatch(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_MousePairMatch
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    MousePairMatch * mat;    
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->len * target->len;    
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 2 * sizeof(int)) > 1000*dpri->kbyte_size) {  
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
        if( (mat=allocate_Expl_MousePairMatch(query, target ,dpri)) == NULL )    {  
          warn("Unable to allocate large MousePairMatch version");   
          return NULL;   
          }  
        calculate_dpenv_MousePairMatch(mat,dpenv);   
        out =  PackAln_read_Expl_MousePairMatch(mat);    
        }  
      else   {  
        mat = allocate_MousePairMatch_only(query, target );  
        calculate_shatter_MousePairMatch(mat,dpenv);     
        out = PackAln_read_Shatter_MousePairMatch(mat);  
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_MousePairMatch(query, target )) == NULL )    {  
          warn("Unable to allocate small MousePairMatch version");   
          return NULL;   
          }  
        out = PackAln_calculate_Small_MousePairMatch(mat,dpenv);     
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_MousePairMatch(query, target ,dpri)) == NULL )    {  
          warn("Unable to allocate large MousePairMatch version");   
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_MousePairMatch(mat);     
          out =  PackAln_read_Expl_MousePairMatch(mat);  
          }  
        }  
      }  


    mat = free_MousePairMatch(mat);  
    return out;  
}    


/* Function:  allocate_MousePairMatch_only(query,target)
 *
 * Descrip:    This function only allocates the MousePairMatch structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * allocate_MousePairMatch_only(GenomePairPara* query,PairMatch* target ) 
{
    MousePairMatch * out;    


    if((out= MousePairMatch_alloc()) == NULL)    {  
      warn("Allocation of basic MousePairMatch structure failed...");    
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->leni = query->len;  
    out->lenj = target->len;     
    return out;  
}    


/* Function:  allocate_Expl_MousePairMatch(query,target,dpri)
 *
 * Descrip:    This function allocates the MousePairMatch structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_MousePairMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * allocate_Expl_MousePairMatch(GenomePairPara* query,PairMatch* target ,DPRunImpl * dpri) 
{
    MousePairMatch * out;    


    out = allocate_MousePairMatch_only(query, target );  
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*2 && dpri->cache->maxlenj >= (out->leni+0))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*2,(out->leni+0),2,out->lenj+1)) == NULL) {  
        warn("Explicit matrix MousePairMatch cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);    
        free_MousePairMatch(out);    
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_MousePairMatch(out);    
    return out;  
}    


/* Function:  init_MousePairMatch(mat)
 *
 * Descrip:    This function initates MousePairMatch matrix when in explicit mode
 *             Called in /allocate_Expl_MousePairMatch
 *
 *
 * Arg:        mat [UNKN ] MousePairMatch which contains explicit basematrix memory [MousePairMatch *]
 *
 */
void init_MousePairMatch(MousePairMatch * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-0);i<mat->query->len;i++)   {  
      for(j= (-1);j<2;j++)   {  
        MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = NEGI;  
        MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = NEGI; 
        }  
      }  
    for(j= (-1);j<mat->target->len;j++)  {  
      for(i= (-0);i<1;i++)   {  
        MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = NEGI;  
        MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = NEGI; 
        }  
      MousePairMatch_EXPL_SPECIAL(mat,i,j,Start) = 0;    
      MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = NEGI;   
      }  
    return;  
}    


/* Function:  recalculate_PackAln_MousePairMatch(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by MousePairMatch
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 */
void recalculate_PackAln_MousePairMatch(PackAln * pal,MousePairMatch * mat) 
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
        case IBD :   
          if( offi == 0 && offj == 1 && prev->state == IBD ) {  
            pau->score = mat->query->extend_ibd + (TargetMatchIBD(j));   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == DIFF )    {  
            pau->score = mat->query->bswitch_diff + (TargetMatchIBD(j));     
            continue;    
            }  
          if( offj == 1 && prev->state == (Start+2) )    {  
            pau->score = 0 + (TargetMatchIBD(j));    
            continue;    
            }  
          warn("In recaluclating PackAln with state IBD, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case DIFF :  
          if( offi == 0 && offj == 1 && prev->state == DIFF )    {  
            pau->score = mat->query->extend_diff + (TargetMatchDIFF(j));     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == IBD ) {  
            pau->score = mat->query->bswitch_ibd + (TargetMatchDIFF(j));     
            continue;    
            }  
          if( offj == 1 && prev->state == (Start+2) )    {  
            pau->score = 0 + (TargetMatchDIFF(j));   
            continue;    
            }  
          warn("In recaluclating PackAln with state DIFF, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case (Start+2) :     
          warn("In recaluclating PackAln with state Start, got a bad source state. Error!"); 
          break; 
        case (End+2) :   
          if( offj == 0 && prev->state == IBD )  {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == DIFF ) {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state End, got a bad source state. Error!");   
          break; 
        default :    
          warn("In recaluclating PackAln got a bad recipient state. Error!");    
        }  
      prev = pau;    
      }  
    return;  
}    
/* divide and conquor macros are next */ 
#define MousePairMatch_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+0)*2+state])    
#define MousePairMatch_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+0)*2+state])    
#define MousePairMatch_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])   
#define MousePairMatch_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])  
#define MousePairMatch_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+0)*2 + state]) 
#define MousePairMatch_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])  
#define MousePairMatch_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 2) + ((i+0) * 2) + (state)])  
#define MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 16) + ((i+0) * 16) + (state * 8) + shadow+1])  
#define MousePairMatch_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])  
/* Function:  allocate_Small_MousePairMatch(query,target)
 *
 * Descrip:    This function allocates the MousePairMatch structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_MousePairMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePairPara*]
 * Arg:        target [UNKN ] target data structure [PairMatch*]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
#define MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])  
MousePairMatch * allocate_Small_MousePairMatch(GenomePairPara* query,PairMatch* target ) 
{
    MousePairMatch * out;    


    out = allocate_MousePairMatch_only(query, target );  
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 0) * 2,16,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix MousePairMatch cannot be allocated, (asking for 2 by %d main cells)",out->leni+1);   
      free_MousePairMatch(out);  
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_MousePairMatch(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for MousePairMatch structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_MousePairMatch 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_MousePairMatch 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_MousePairMatch(MousePairMatch * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for MousePairMatch due to wrong type of matrix");  
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_MousePairMatch(mat,dpenv);   
    score = start_end_find_end_MousePairMatch(mat,&endj);    
    out->score = score;  
    stopstate = End;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,0);  
    startj = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,1);  
    startstate = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,2);  
    stopi = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,3);   
    stopj = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,4);   
    stopstate = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,5);   
    temp = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,6);    
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);  


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_MousePairMatch(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL); 
    if( stopi == MousePairMatch_READ_OFF_ERROR || stopstate != Start )   {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_MousePairMatch(out,mat);   
      return out;    
      }  


    /* Ok. Put away start start... */ 
    pau = PackAlnUnit_alloc();   
    pau->i = stopi;  
    pau->j = stopj;  
    pau->state = stopstate + 2;  
    add_PackAln(out,pau);    


    log_full_error(REPORT,0,"Alignment recovered");  
    stop_reporting();    
    invert_PackAln(out); 
    recalculate_PackAln_MousePairMatch(out,mat); 
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_MousePairMatch(mat)
 *
 * Descrip:    This function calculates an alignment for MousePairMatch structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_MousePairMatch 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_MousePairMatch
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_MousePairMatch 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_MousePairMatch(MousePairMatch * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_MousePairMatch(mat,NULL);    
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_MousePairMatch(mat);  
    return out;  
}    


/* Function:  AlnRangeSet_from_MousePairMatch(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for MousePairMatch structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_MousePairMatch 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_MousePairMatch
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_MousePairMatch(MousePairMatch * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_MousePairMatch");    
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_MousePairMatch(mat,&jpos);   
    state = End; 


    while( (temp = AlnRange_build_MousePairMatch(mat,jpos,state,&jpos,&state)) != NULL)  
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_MousePairMatch(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_MousePairMatch
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_MousePairMatch(MousePairMatch * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_MousePairMatch");    
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_MousePairMatch(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE) {  
      warn("In AlnRanger_build_MousePairMatch alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj); 
      return NULL;   
      }  
    if( state == Start || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);   
    out->startj = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);   
    out->startstate = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);   
    out->stopi = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);    
    out->stopj = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);    
    out->stopstate = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);    
    out->startscore = MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);   
    out->stopscore = MousePairMatch_DC_SHADOW_SPECIAL(mat,0,jpos,state); 


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_MousePairMatch(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);   
    if( jpos == MousePairMatch_READ_OFF_ERROR)   {  
      warn("In AlnRange_build_MousePairMatch alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);  
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = MousePairMatch_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);   
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
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
boolean read_hidden_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_MousePairMatch(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  


      if( i == MousePairMatch_READ_OFF_ERROR)    {  
        warn("In MousePairMatch hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);    
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
        warn("In MousePairMatch hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);  
        return FALSE;    
        }  
      }  
    warn("In MousePairMatch hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);   
    return FALSE;    
}    


/* Function:  max_hidden_MousePairMatch(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
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
int max_hidden_MousePairMatch(MousePairMatch * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = MousePairMatch_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->len)   {  
      warn("In MousePairMatch matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = MousePairMatch_HIDDEN_MATRIX(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case IBD :     
        /* Not allowing special sources.. skipping Start */ 
        temp = cscore - (mat->query->bswitch_diff) -  (TargetMatchIBD(j));   
        if( temp == MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,DIFF) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DIFF;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,DIFF);    
            }  
          return MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,DIFF);     
          }  
        temp = cscore - (mat->query->extend_ibd) -  (TargetMatchIBD(j)); 
        if( temp == MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,IBD) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = IBD;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,IBD); 
            }  
          return MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,IBD);  
          }  
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DIFF :    
        /* Not allowing special sources.. skipping Start */ 
        temp = cscore - (mat->query->bswitch_ibd) -  (TargetMatchDIFF(j));   
        if( temp == MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,IBD) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = IBD;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,IBD); 
            }  
          return MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,IBD);  
          }  
        temp = cscore - (mat->query->extend_diff) -  (TargetMatchDIFF(j));   
        if( temp == MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,DIFF) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DIFF;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,DIFF);    
            }  
          return MousePairMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,DIFF);     
          }  
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_MousePairMatch(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
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
boolean read_special_strip_MousePairMatch(MousePairMatch * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != Start)   { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 2; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_MousePairMatch(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);  
      if( i == MousePairMatch_READ_OFF_ERROR)    {  
        warn("In special strip read MousePairMatch, got a bad read off error. Sorry!");  
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != Start && j < MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))  {  
      warn("In special strip read MousePairMatch, at special [%d] state [%d] overshot!",j,state);    
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 2;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_MousePairMatch(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
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
int max_special_strip_MousePairMatch(MousePairMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = MousePairMatch_READ_OFF_ERROR;   
    if( isspecial == FALSE ) {  
      warn("In special strip max function for MousePairMatch, got a non special start point. Problem! (bad!)");  
      return (-1);   
      }  


    if( j < 0 || j > mat->target->len)   {  
      warn("In MousePairMatch matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    cscore = MousePairMatch_DC_SHADOW_SPECIAL(mat,i,j,state);    
    switch(state)    { /*switch on special states*/ 
      case Start :   
      case End :     
        /* Source DIFF is not a special */ 
        /* Source IBD is not a special */ 
      default:   
        warn("Major problem (!) - in MousePairMatch special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_MousePairMatch(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
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
int max_matrix_to_special_MousePairMatch(MousePairMatch * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = MousePairMatch_READ_OFF_ERROR;   


    if( j < 0 || j > mat->lenj)  {  
      warn("In MousePairMatch matrix to special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case IBD :     
        temp = cscore - (0) -  (TargetMatchIBD(j));  
        if( temp == MousePairMatch_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,Start) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Start; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - MousePairMatch_DC_SHADOW_SPECIAL(mat,i-0,j-1,Start);   
            }  
          return MousePairMatch_DC_SHADOW_MATRIX(mat,i - 0,j - 1,Start) ;    
          }  
        /* Source DIFF is not a special, should not get here! */ 
        /* Source IBD is not a special, should not get here! */ 
        warn("Major problem (!) - in MousePairMatch matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case DIFF :    
        temp = cscore - (0) -  (TargetMatchDIFF(j));     
        if( temp == MousePairMatch_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,Start) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Start; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - MousePairMatch_DC_SHADOW_SPECIAL(mat,i-0,j-1,Start);   
            }  
          return MousePairMatch_DC_SHADOW_MATRIX(mat,i - 0,j - 1,Start) ;    
          }  
        /* Source IBD is not a special, should not get here! */ 
        /* Source DIFF is not a special, should not get here! */ 
        warn("Major problem (!) - in MousePairMatch matrix to special read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_MousePairMatch(mat,starti,startj,stopi,stopj);   


    MousePairMatch_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;  


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          MousePairMatch_HIDDEN_MATRIX(mat,i,j,IBD) = NEGI;  
          MousePairMatch_HIDDEN_MATRIX(mat,i,j,DIFF) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;  
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_HIDDEN_MATRIX(mat,i,j,IBD) = score;  
        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;    
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_HIDDEN_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_HIDDEN_MATRIX(mat,i,j,DIFF) = score; 
        /* Finished calculating state DIFF */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_MousePairMatch(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_MousePairMatch(MousePairMatch * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-0);i<=stopi;i++) {  
        MousePairMatch_HIDDEN_MATRIX(mat,i,j,IBD) = NEGI;
   
        MousePairMatch_HIDDEN_MATRIX(mat,i,j,DIFF) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_MousePairMatch
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_MousePairMatch to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [MousePairMatch *]
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
boolean full_dc_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_MousePairMatch");  
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);     
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)   {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)    {  
      warn("In divide and conquor for MousePairMatch, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);   
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= MousePairMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);    
    lstartj= MousePairMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);    
    lstate = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);    


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_MousePairMatch(mat,MousePairMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),MousePairMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),MousePairMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_MousePairMatch(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE) {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_MousePairMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
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
boolean do_dc_single_pass_MousePairMatch(MousePairMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_MousePairMatch(mat); 


    MousePairMatch_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;   
    run_up_dc_MousePairMatch(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);   
    push_dc_at_merge_MousePairMatch(mat,starti,stopi,halfj,&halfj,dpenv);    
    follow_on_dc_MousePairMatch(mat,starti,stopi,halfj,stopj,dpenv,perc_done);   
    return TRUE; 
}    


/* Function:  push_dc_at_merge_MousePairMatch(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_MousePairMatch(MousePairMatch * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<1;count++,j++)    {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = NEGI;   
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,0) = (-100);    
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,1) = (-100);    
          MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = NEGI;  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,0) = (-100);   
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state IBD, pushing when j - offj <= mergej */ 
        score = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;   
        if( j - 1 <= mergej) {  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,0) = i-0;   
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,1) = j-1;   
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,2) = IBD;   
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,3) = i; 
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,4) = j; 
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,5) = IBD;   
          }  
        else {  
          for(k=0;k<7;k++)   
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,k) = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,IBD,k);   
          }  


        temp = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,0) = i-0; 
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,1) = j-1; 
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,2) = DIFF;    
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,3) = i;   
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,4) = j;   
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,5) = IBD; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,k) = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DIFF,k);    
            }  
          }  
        /* Add any movement independant score */ 
        score += TargetMatchIBD(j);  
        MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = score;    
        /* Finished with state IBD */ 


        /* For state DIFF, pushing when j - offj <= mergej */ 
        score = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;     
        if( j - 1 <= mergej) {  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,0) = i-0;  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,1) = j-1;  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,2) = DIFF; 
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,3) = i;    
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,4) = j;    
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,5) = DIFF; 
          }  
        else {  
          for(k=0;k<7;k++)   
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,k) = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DIFF,k); 
          }  


        temp = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,0) = i-0;    
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,1) = j-1;    
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,2) = IBD;    
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,3) = i;  
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,4) = j;  
            MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,5) = DIFF;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,k) = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,IBD,k);    
            }  
          }  
        /* Add any movement independant score */ 
        score += TargetMatchDIFF(j);     
        MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = score;   
        /* Finished with state DIFF */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_MousePairMatch(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_MousePairMatch(MousePairMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = NEGI;   
          MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,IBD,k);    
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DIFF,k); 
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = score;   
        for(k=0;k<7;k++) 
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DIFF,k);   
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,IBD,k);  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = score;  
        for(k=0;k<7;k++) 
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DIFF */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_MousePairMatch(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_MousePairMatch(MousePairMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = NEGI;   
          MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;   
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = score;   
        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;     
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = score;  
        /* Finished calculating state DIFF */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_MousePairMatch(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 */
}    
void init_dc_MousePairMatch(MousePairMatch * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->query->len;i++)  {  
        MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = NEGI; 
        MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = NEGI;    
        for(k=0;k<7;k++) {  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,k) = (-1);  
          MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_MousePairMatch(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [MousePairMatch *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_MousePairMatch(MousePairMatch * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = MousePairMatch_DC_SHADOW_SPECIAL(mat,0,mat->target->len-1,End);    
    maxj = mat->target->len-1;   
    for(j= mat->target->len-2 ;j >= 0 ;j--)  {  
      if( MousePairMatch_DC_SHADOW_SPECIAL(mat,0,j,End) > max )  {  
        max = MousePairMatch_DC_SHADOW_SPECIAL(mat,0,j,End); 
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_MousePairMatch(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [MousePairMatch]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_MousePairMatch(MousePairMatch *mat,DPEnvelope * dpenv) 
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
    lenj = mat->target->len; 
    total = leni * lenj; 


    score_pointers = (int *) calloc (1 * (leni + 0) * 2,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 0) * 2 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,IBD) = NEGI;   
          MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,DIFF) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd + (TargetMatchIBD(j));     
        /* assign local shadown pointer */ 
        localsp = &(MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,IBD,0));  
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff +(TargetMatchIBD(j));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DIFF,0));   
          }  
        /* From state Start to state IBD */ 
        temp = MousePairMatch_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,Start) + 0 + (TargetMatchIBD(j));    
        if( temp  > score )  {  
          score = temp;  
          /* This state [Start] is a special for IBD... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= IBD;   
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,IBD) = score;   
        for(k=0;k<7;k++) 
          MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,IBD,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state IBD is a source for special End */ 
        temp = score + (0) + (0) ;   
        if( temp > MousePairMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) )   {  
          MousePairMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,k) = MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,IBD,k);  
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,6) = MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,IBD,6);    
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,3) = i;    
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,4) = j;    
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,5) = IBD;  
          }  




        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff + (TargetMatchDIFF(j));  
        /* assign local shadown pointer */ 
        localsp = &(MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,DIFF,0)); 
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd +(TargetMatchDIFF(j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,IBD,0));    
          }  
        /* From state Start to state DIFF */ 
        temp = MousePairMatch_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,Start) + 0 + (TargetMatchDIFF(j));   
        if( temp  > score )  {  
          score = temp;  
          /* This state [Start] is a special for DIFF... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= DIFF;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         MousePairMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,DIFF) = score;  
        for(k=0;k<7;k++) 
          MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DIFF,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* state DIFF is a source for special End */ 
        temp = score + (0) + (0) ;   
        if( temp > MousePairMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) )   {  
          MousePairMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) = temp;  
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,k) = MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DIFF,k); 
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,6) = MousePairMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DIFF,6);   
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,3) = i;    
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,4) = j;    
          MousePairMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,5) = DIFF; 
          }  




        /* Finished calculating state DIFF */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_MousePairMatch(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 */
void init_start_end_linear_MousePairMatch(MousePairMatch * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->query->len;i++)  {  
        MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,IBD) = NEGI; 
        MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,IBD,0) = (-1);    
        MousePairMatch_DC_SHADOW_MATRIX(mat,i,j,DIFF) = NEGI;    
        MousePairMatch_DC_SHADOW_MATRIX_SP(mat,i,j,DIFF,0) = (-1);   
        }  
      }  


    for(j=(-1);j<mat->target->len;j++)   {  
      MousePairMatch_DC_SHADOW_SPECIAL(mat,0,j,Start) = 0;   
      MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,j,Start,0) = j;  
      MousePairMatch_DC_SHADOW_SPECIAL(mat,0,j,End) = NEGI;  
      MousePairMatch_DC_SHADOW_SPECIAL_SP(mat,0,j,End,0) = (-1); 
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_MousePairMatch(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_MousePairMatch(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_MousePairMatch();    
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "Q","END" };    
/* Function:  AlnConvertSet_MousePairMatch(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "IBD","DIFF","END" };  
AlnConvertSet * AlnConvertSet_MousePairMatch(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = IBD;   
    acu->state2 = IBD;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DIFF;  
    acu->state2 = IBD;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Start + 2; 
    acu->is_from_special = TRUE; 
    acu->state2 = IBD;   
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DIFF;  
    acu->state2 = DIFF;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = IBD;   
    acu->state2 = DIFF;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Start + 2; 
    acu->is_from_special = TRUE; 
    acu->state2 = DIFF;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = IBD;   
    acu->state2 = End + 2;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DIFF;  
    acu->state2 = End + 2;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    return out;  
}    


/* Function:  PackAln_read_Expl_MousePairMatch(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_MousePairMatch(MousePairMatch * mat) 
{
    MousePairMatch_access_func_holder holder;    


    holder.access_main    = MousePairMatch_explicit_access_main; 
    holder.access_special = MousePairMatch_explicit_access_special;  
    return PackAln_read_generic_MousePairMatch(mat,holder);  
}    


/* Function:  MousePairMatch_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MousePairMatch_explicit_access_main(MousePairMatch * mat,int i,int j,int state) 
{
    return MousePairMatch_EXPL_MATRIX(mat,i,j,state);    
}    


/* Function:  MousePairMatch_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MousePairMatch_explicit_access_special(MousePairMatch * mat,int i,int j,int state) 
{
    return MousePairMatch_EXPL_SPECIAL(mat,i,j,state);   
}    


/* Function:  PackAln_read_generic_MousePairMatch(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:          h [UNKN ] Undocumented argument [MousePairMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_MousePairMatch(MousePairMatch * mat,MousePairMatch_access_func_holder h) 
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


    out->score =  find_end_MousePairMatch(mat,&i,&j,&state,&isspecial,h);    


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in MousePairMatch_basic_read, returning a mess.(Sorry!)",out->len);   
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 2;     
    prev=pau;    
    while( state != Start || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_MousePairMatch(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      else   
        max_calc_MousePairMatch(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);   
      if(i == MousePairMatch_READ_OFF_ERROR || j == MousePairMatch_READ_OFF_ERROR || state == MousePairMatch_READ_OFF_ERROR )    {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in MousePairMatch_basic_read, returning partial alignment",out->len);   
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 2;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_MousePairMatch(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [MousePairMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_MousePairMatch(MousePairMatch * mat,int * ri,int * rj,int * state,boolean * isspecial,MousePairMatch_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->target->len-1,End); 
    maxj = mat->target->len-1;   
    for(j= mat->target->len-2 ;j >= 0 ;j--)  {  
      if( (temp =(*h.access_special)(mat,0,j,End)) > max )   {  
        max = temp;  
        maxj = j;    
        }  
      }  


    if( ri != NULL)  
       *ri = 0;  
    if( rj != NULL)  
       *rj = maxj;   
    if( state != NULL)   
       *state = End; 
    if( isspecial != NULL)   
       *isspecial = TRUE;    


    return max;  
}    


/* Function:  MousePairMatch_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void MousePairMatch_debug_show_matrix(MousePairMatch * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->len;i++) {  
      for(j=startj;j<stopj && j < mat->target->len;j++)  {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State IBD %d\n",MousePairMatch_EXPL_MATRIX(mat,i,j,IBD));   
        fprintf(ofp,"State DIFF %d\n",MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF)); 
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_MousePairMatch(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [MousePairMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_MousePairMatch(MousePairMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MousePairMatch_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = MousePairMatch_READ_OFF_ERROR;   


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->len)   {  
      warn("In MousePairMatch matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);  
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case IBD :     
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          temp = cscore - (0) -  (TargetMatchIBD(j));    
          if( temp == (*h.access_special)(mat,i - 0,j - 1,Start) )   {  
            *reti = i - 0;   
            *retj = j - 1;   
            *retstate = Start;   
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,Start);  
              }  
            return (*h.access_main)(mat,i - 0,j - 1,Start);  
            }  
          }  
        temp = cscore - (mat->query->bswitch_diff) -  (TargetMatchIBD(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,DIFF) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DIFF;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,DIFF);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,DIFF);     
          }  
        temp = cscore - (mat->query->extend_ibd) -  (TargetMatchIBD(j)); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,IBD) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = IBD;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,IBD); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,IBD);  
          }  
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case DIFF :    
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          temp = cscore - (0) -  (TargetMatchDIFF(j));   
          if( temp == (*h.access_special)(mat,i - 0,j - 1,Start) )   {  
            *reti = i - 0;   
            *retj = j - 1;   
            *retstate = Start;   
            *retspecial = TRUE;  
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,Start);  
              }  
            return (*h.access_main)(mat,i - 0,j - 1,Start);  
            }  
          }  
        temp = cscore - (mat->query->bswitch_ibd) -  (TargetMatchDIFF(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,IBD) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = IBD;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,IBD); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,IBD);  
          }  
        temp = cscore - (mat->query->extend_diff) -  (TargetMatchDIFF(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,DIFF) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = DIFF;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,DIFF);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,DIFF);     
          }  
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_MousePairMatch(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MousePairMatch *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [MousePairMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_MousePairMatch(MousePairMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MousePairMatch_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = MousePairMatch_READ_OFF_ERROR;   


    if( j < 0 || j > mat->target->len)   {  
      warn("In MousePairMatch matrix special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case Start :   
      case End :     
        /* source DIFF is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,DIFF) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = DIFF;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,DIFF);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,DIFF) ;  
            }  
          } /* end of for i >= 0 */ 
        /* source IBD is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,IBD) )    {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = IBD; 
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,IBD);   
              }  
            return (*h.access_main)(mat,i - 0,j - 0,IBD) ;   
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in MousePairMatch read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_MousePairMatch(mat)
 *
 * Descrip:    This function calculates the MousePairMatch matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_MousePairMatch
 *
 *
 * Arg:        mat [UNKN ] MousePairMatch which contains explicit basematrix memory [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_MousePairMatch(MousePairMatch * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_MousePairMatch, passed a non Explicit matrix type, cannot calculate!"); 
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("MousePairMatch Matrix calculation: ");  
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;    
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state IBD */ 
          temp = MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + 0;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = score;    


        /* state IBD is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_EXPL_SPECIAL(mat,i,j,End) )  {  
            MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = temp;     
            }  


          }  


        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;  
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state DIFF */ 
          temp = MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + 0;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = score;   


        /* state DIFF is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_EXPL_SPECIAL(mat,i,j,End) )  {  
            MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = temp;     
            }  


          }  


        /* Finished calculating state DIFF */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_MousePairMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MousePairMatch matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] MousePairMatch which contains explicit basematrix memory [MousePairMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_MousePairMatch(MousePairMatch * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_MousePairMatch, passed a non Explicit matrix type, cannot calculate!"); 
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


    for(j=startj-1;j<endj;j++)   {  
      for(i=0;i<mat->leni;i++)   {  
        MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = NEGI;  
        MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = NEGI; 
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      MousePairMatch_EXPL_SPECIAL(mat,i,j,Start) = 0;    
      MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = NEGI;   
      }  


    start_reporting("MousePairMatch Matrix calculation: ");  
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
          MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = NEGI;    
          MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state IBD */ 
        /* setting first movement to score */ 
        score = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) + mat->query->extend_ibd;    
        /* From state DIFF to state IBD */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->bswitch_diff;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state IBD */ 
          temp = MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + 0;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for IBD */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchIBD(j);     
         MousePairMatch_EXPL_MATRIX(mat,i,j,IBD) = score;    


        /* state IBD is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_EXPL_SPECIAL(mat,i,j,End) )  {  
            MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = temp;     
            }  


          }  


        /* Finished calculating state IBD */ 


        /* For state DIFF */ 
        /* setting first movement to score */ 
        score = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,DIFF) + mat->query->extend_diff;  
        /* From state IBD to state DIFF */ 
        temp = MousePairMatch_EXPL_MATRIX(mat,i-0,j-1,IBD) + mat->query->bswitch_ibd;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state DIFF */ 
          temp = MousePairMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + 0;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for DIFF */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDIFF(j);    
         MousePairMatch_EXPL_MATRIX(mat,i,j,DIFF) = score;   


        /* state DIFF is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MousePairMatch_EXPL_SPECIAL(mat,i,j,End) )  {  
            MousePairMatch_EXPL_SPECIAL(mat,i,j,End) = temp;     
            }  


          }  


        /* Finished calculating state DIFF */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  MousePairMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * MousePairMatch_alloc(void) 
{
    MousePairMatch * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MousePairMatch *) ckalloc (sizeof(MousePairMatch))) == NULL)    {  
      warn("MousePairMatch_alloc failed ");  
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


/* Function:  free_MousePairMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MousePairMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MousePairMatch *]
 *
 */
MousePairMatch * free_MousePairMatch(MousePairMatch * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MousePairMatch obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
