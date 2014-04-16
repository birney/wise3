#ifdef _cplusplus
extern "C" {
#endif
#include "mousethreestate.h"



# line 113 "mousethreestate.dy"
AncestralVarSet * ancestral_prediction_MouseThreeState(GenoVarSet * gvs,char * mus,char * cast,char * dom,GenomePara * gp,DPRunImpl * dpri)
{
  AncestralVarSet * out;
  AncestralIndividual * ai;
  char strain_state[1024];
  
  int i;
  int j;


  int mus_index;
  int cast_index;
  int dom_index;
 
  int ind;

  assert(gvs != NULL);
  assert(mus != NULL);
  assert(dom != NULL);
  assert(cast != NULL);
  assert(gp != NULL);

  fprintf(stderr,"Into ancestral set\n");

  out= AncestralVarSet_alloc_len(gvs->len);

  for(i=0;i<gvs->ind_len;i++) {

    if( strcmp(gvs->ind[i]->id,mus) == 0 ) {
      ai = AncestralIndividual_alloc();
      ai->name = stringalloc(gvs->ind[i]->id);
      ai->single_letter = 'M';
      strain_state[i] = out->anc_len;
      mus_index = out->anc_len;
      add_anc_AncestralVarSet(out,ai);
      continue;
    }

    if( strcmp(gvs->ind[i]->id,dom) == 0 ) {
      ai = AncestralIndividual_alloc();
      ai->name = stringalloc(gvs->ind[i]->id);
      ai->single_letter = 'D';
      dom_index = out->anc_len;
      strain_state[i] = out->anc_len;
      add_anc_AncestralVarSet(out,ai);
      continue;
    }

    if( strcmp(gvs->ind[i]->id,cast) == 0 ) {
      ai = AncestralIndividual_alloc();
      ai->name = stringalloc(gvs->ind[i]->id);
      ai->single_letter = 'C';
      cast_index = out->anc_len;
      strain_state[i] = out->anc_len;
      add_anc_AncestralVarSet(out,ai);
      continue;
    }

    strain_state[i] = 100;
    add_ind_AncestralVarSet(out,hard_link_Individual(gvs->ind[i]));
  }

  fprintf(stderr,"About to make chromosomes\n");

  create_blank_AncestralChromosomes(out,gvs);
    
  /* now loop over each chromosome, 
     run HMM and put in the designation
  */

  fprintf(stderr,"Entering main loop\n");

  for(j=0;j<gvs->len;j++) {
    ind = 0;
    for(i=0;i<gvs->ind_len;i++) {
      if( strain_state[i] != 100 ) {
	/* it is an ancestor, skip it! */
	continue;
      }
      
      /* now make a packaln for this */
      
      auto SnpMatch * snpm;
      auto AlnBlock * alb;
      auto PackAln * pal;

      info("Going to handle strain %s",gvs->ind[i]->id);      
      
      snpm = new_SnpMatch(gvs->chr[j],gvs,gvs->ind[i]->id,mus,cast,dom);

      assert(snpm != NULL);

      pal = PackAln_bestmemory_MouseSNPMatch(gp,snpm,NULL,dpri);

      alb = convert_PackAln_to_AlnBlock_MouseSNPMatch(pal);

      populate_individual_Ancestral_decoding_ThreeState(alb,out->chr[j],ind,mus_index,cast_index,dom_index);
      
      ind++;
      free_AlnBlock(alb);
      free_PackAln(pal);
      free_SnpMatch(snpm);
    }
  }



  return out;
}


# line 224 "mousethreestate.dy"
void populate_individual_Ancestral_decoding_ThreeState(AlnBlock * alb,AncestralVarChr * avc,int ind_number,short int mus_index,short int cast_index,short int dom_index)
{
  AlnColumn * alc;
  int j;

  
  assert(alb != NULL);
  assert(avc != NULL);
  
  alc = alb->start->next;

  assert(alc != NULL);

  j = 0;
  while( strcmp(alc->alu[0]->text_label,"END") != 0 ) {
    
   	 
    if( strcmp(alc->alu[1]->text_label,"MUSCULUS") == 0 ) {
      avc->loci[j]->ind[ind_number] = mus_index;
    } else if ( strcmp(alc->alu[1]->text_label,"CAST") == 0 ) {
      avc->loci[j]->ind[ind_number] = cast_index;
    } else if ( strcmp(alc->alu[1]->text_label,"DOMESTICUS") == 0 ) {
      avc->loci[j]->ind[ind_number] = dom_index;
    } else {
      avc->loci[j]->ind[ind_number] = 'U';
    }

/*    fprintf(stderr,"Looking at %s - produced %d\n",alc->alu[1]->text_label,avc->loci[j]->ind[ind_number]); */

    j++;
    alc = alc->next;

    if( j > avc->len ) {
      warn("Ran over chromosome lenght - AlnBlock is out of sync from Locus block");
    }
  }

}


# line 264 "mousethreestate.dy"
void create_blank_AncestralChromosomes(AncestralVarSet * avs,GenoVarSet * gvs)
{
  int i;
  int j;
  int k;

  assert(avs != NULL);
  assert(gvs != NULL);
  
  for(i=0;i<gvs->len;i++) {
    auto AncestralVarChr * avc;
    
    avc = AncestralVarChr_alloc_len(gvs->chr[i]->len);
    add_AncestralVarSet(avs,avc);

    avc->chr = stringalloc(gvs->chr[i]->chr);

    for(j=0;j < gvs->chr[i]->len;j++) {
      auto AncestralVarLocus * avl;

      avl = AncestralVarLocus_alloc();
      avl->var = hard_link_VarLocus(gvs->chr[i]->loci[j]->var);
      avl->ind = calloc(avs->ind_len,sizeof(char));
      avl->anc = calloc(avs->anc_len,sizeof(char));

      for(k=0;k<avs->ind_len;k++) {
	avl->ind[k] = -1;
      }
      add_AncestralVarChr(avc,avl);
    }
  }

}


# line 299 "mousethreestate.dy"
GenomePara * new_GenomePara(Probability match,Probability bswitch)
{
  GenomePara * out;

  out = GenomePara_alloc();

  out->match = Probability2Score(match);
  out->mismatch = Probability2Score(1.0 - match);

  out->extend = Probability2Score(1.0 - 3*(bswitch));
  out->bswitch = Probability2Score(bswitch);

  out->len = 1;
  return(out);
}

# line 315 "mousethreestate.dy"
void show_SnpMatchStats(SnpMatch * snpm,FILE * ofp)
{
  int i;
  long int mus  = 0;
  long int dom = 0;
  long int cast = 0;

  assert(snpm != NULL);
  assert(ofp != NULL);

  for(i=0;i<snpm->len;i++) {
    if( snpm->musculus[i] == 1 ) {
      mus++;
    }
    if( snpm->domesticus[i] == 1 ) {
      dom++;
    }
    if( snpm->cast[i] == 1 ) {
      cast++;
    }
  }

  fprintf(ofp,"Musculus %.2f Domesticus %.2f Cast %.2f %ld\n",
	  ((double)mus/(double)snpm->len),
	  ((double)dom/(double)snpm->len),
	  ((double)cast/(double)snpm->len),
	  snpm->len
	  );
}
    
  

# line 347 "mousethreestate.dy"
SnpMatch * new_SnpMatch(GenoVarChr * chr,GenoVarSet * set,char * test_strain,char * musculus,char * cast,char * domesticus)
{
  int strain_index;
  int musculus_index;
  int cast_index;
  int domesticus_index;

  int i;
  int j;

  SnpMatch * out;
  int len;
  
  len = number_of_simple_snp_loci_GenoVarChr(chr);

  out = SnpMatch_alloc();
  out->len = len;
  out->musculus = calloc(len,sizeof(char));
  out->cast = calloc(len,sizeof(char));
  out->domesticus = calloc(len,sizeof(char));
  
  strain_index = individual_index_from_string_GenoVarSet(set,test_strain);
  assert(strain_index != -1);

  musculus_index = individual_index_from_string_GenoVarSet(set,musculus);
  assert(musculus_index != -1);

  cast_index = individual_index_from_string_GenoVarSet(set,cast);
  assert(cast_index != -1);

  domesticus_index = individual_index_from_string_GenoVarSet(set,domesticus);
  assert(domesticus_index != -1);

  /*  fprintf(stdout,"Got indexes %d %d %d %d\n",
	  strain_index,
	  musculus_index,
	  cast_index,
	  domesticus_index);
  */

  j = 0;
  for(i=0;i<chr->len;i++) {
    auto GenoVarLocus * l;
    l = chr->loci[i];
    
    if( l->var->locus_type != SIMPLE_SNP_LOCUS ) {
      /*fprintf(stdout,"Non simple locus at %d with %d!\n",i,l->var->locus_type);*/
      continue;
    }

    if( l->ind[strain_index] == l->ind[musculus_index] ) {
      out->musculus[j] = 1;
    } else {
      out->musculus[j] = 0;
    }


    if( l->ind[strain_index] == l->ind[cast_index] ) {
      out->cast[j] = 1;
    } else {
      out->cast[j] = 0;
    }

    if( l->ind[strain_index] == l->ind[domesticus_index] ) {
      out->domesticus[j] = 1;
    } else {
      out->domesticus[j] = 0;
    }

    /*    fprintf(stdout,"%d,%d Testing %d vs %d - outcome %d %d %d\n",i,j,l->ind[strain_index],l->ind[musculus_index],out->musculus[j],out->cast[j],out->domesticus[j]); */

    j++;
  }



  return out;
}

/* Function:  hard_link_MouseSNPMatch_Posterior(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MouseSNPMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch_Posterior *]
 *
 */
MouseSNPMatch_Posterior * hard_link_MouseSNPMatch_Posterior(MouseSNPMatch_Posterior * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MouseSNPMatch_Posterior object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MouseSNPMatch_Posterior_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch_Posterior *]
 *
 */
MouseSNPMatch_Posterior * MouseSNPMatch_Posterior_alloc(void) 
{
    MouseSNPMatch_Posterior * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MouseSNPMatch_Posterior *) ckalloc (sizeof(MouseSNPMatch_Posterior))) == NULL)  {  
      warn("MouseSNPMatch_Posterior_alloc failed "); 
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


/* Function:  free_MouseSNPMatch_Posterior(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MouseSNPMatch_Posterior *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch_Posterior *]
 *
 */
MouseSNPMatch_Posterior * free_MouseSNPMatch_Posterior(MouseSNPMatch_Posterior * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MouseSNPMatch_Posterior obj. Should be trappable");   
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
      free_MouseSNPMatch(obj->forward);  
    if( obj->backward != NULL)   
      free_MouseSNPMatch(obj->backward);     


    ckfree(obj); 
    return NULL; 
}    


# line 418 "mousethreestate.c"
/* Function:  hard_link_SnpMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SnpMatch *]
 *
 * Return [UNKN ]  Undocumented return value [SnpMatch *]
 *
 */
SnpMatch * hard_link_SnpMatch(SnpMatch * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SnpMatch object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SnpMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SnpMatch *]
 *
 */
SnpMatch * SnpMatch_alloc(void) 
{
    SnpMatch * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SnpMatch *) ckalloc (sizeof(SnpMatch))) == NULL)    {  
      warn("SnpMatch_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->len = 0;    
    out->musculus = NULL;    
    out->cast = NULL;    
    out->domesticus = NULL;  


    return out;  
}    


/* Function:  free_SnpMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SnpMatch *]
 *
 * Return [UNKN ]  Undocumented return value [SnpMatch *]
 *
 */
SnpMatch * free_SnpMatch(SnpMatch * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SnpMatch obj. Should be trappable");  
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
    if( obj->musculus != NULL)   
      ckfree(obj->musculus);     
    if( obj->cast != NULL)   
      ckfree(obj->cast);     
    if( obj->domesticus != NULL) 
      ckfree(obj->domesticus);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GenomePara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenomePara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePara *]
 *
 */
GenomePara * hard_link_GenomePara(GenomePara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenomePara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenomePara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenomePara *]
 *
 */
GenomePara * GenomePara_alloc(void) 
{
    GenomePara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenomePara *) ckalloc (sizeof(GenomePara))) == NULL)    {  
      warn("GenomePara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->match = 0;  
    out->mismatch = 0;   
    out->extend = 0; 
    out->bswitch = 0;    
    out->len = 0;    


    return out;  
}    


/* Function:  free_GenomePara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenomePara *]
 *
 * Return [UNKN ]  Undocumented return value [GenomePara *]
 *
 */
GenomePara * free_GenomePara(GenomePara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenomePara obj. Should be trappable");    
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
  /*            Wed Jan 26 11:49:09 2011           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define Musculus 0   
#define Cast 1   
#define Domesticus 2 


#define Start 0  
#define End 1    


#define MouseSNPMatch_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*3)+STATE][i+0]   
#define MouseSNPMatch_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]  
#define MouseSNPMatch_READ_OFF_ERROR -2
 


#define MouseSNPMatch_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+0)*3)+STATE] 
#define MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE]  




#define MouseSNPMatch_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]   
#define MouseSNPMatch_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)   


/* Function:  PackAln_read_Shatter_MouseSNPMatch(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_MouseSNPMatch(MouseSNPMatch * mat) 
{
    MouseSNPMatch_access_func_holder holder;     


    holder.access_main    = MouseSNPMatch_shatter_access_main;   
    holder.access_special = MouseSNPMatch_shatter_access_special;    
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_MouseSNPMatch(mat,holder);   
}    


/* Function:  MouseSNPMatch_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MouseSNPMatch_shatter_access_main(MouseSNPMatch * mat,int i,int j,int state) 
{
    return MouseSNPMatch_SHATTER_MATRIX(mat,i,j,state);  
}    


/* Function:  MouseSNPMatch_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MouseSNPMatch_shatter_access_special(MouseSNPMatch * mat,int i,int j,int state) 
{
    return MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,state); 
}    


/* Function:  calculate_shatter_MouseSNPMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MouseSNPMatch matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [MouseSNPMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_MouseSNPMatch(MouseSNPMatch * mat,DPEnvelope * dpenv) 
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


    mat->shatter = new_ShatterMatrix(dpenv,3,lenj,2);    
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


    start_reporting("MouseSNPMatch Matrix calculation: ");   
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




        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = SIG_0_1[Musculus] + mat->query->extend;  
        /* From state Cast to state Musculus */ 
        temp = SIG_0_1[Cast] + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = SIG_0_1[Domesticus] + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Musculus */ 
          temp = MouseSNPMatch_SHATTER_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         SIG_0_0[Musculus] = score;  


        /* state Musculus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,End) )    {  
            MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,End) = temp;   
            }  


          }  


        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = SIG_0_1[Musculus] + mat->query->bswitch;     
        /* From state Cast to state Cast */ 
        temp = SIG_0_1[Cast] + mat->query->extend;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Cast */ 
        temp = SIG_0_1[Domesticus] + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Cast */ 
          temp = MouseSNPMatch_SHATTER_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         SIG_0_0[Cast] = score;  


        /* state Cast is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,End) )    {  
            MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,End) = temp;   
            }  


          }  


        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = SIG_0_1[Musculus] + mat->query->bswitch;     
        /* From state Cast to state Domesticus */ 
        temp = SIG_0_1[Cast] + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = SIG_0_1[Domesticus] + mat->query->extend;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Domesticus */ 
          temp = MouseSNPMatch_SHATTER_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;     
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         SIG_0_0[Domesticus] = score;    


        /* state Domesticus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,End) )    {  
            MouseSNPMatch_SHATTER_SPECIAL(mat,i,j,End) = temp;   
            }  


          }  


        /* Finished calculating state Domesticus */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_MouseSNPMatch(dbsi,out,query,target)
 *
 * Descrip:    This function makes a database search of MouseSNPMatch
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePara*]
 * Arg:        target [UNKN ] Undocumented argument [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_MouseSNPMatch(DBSearchImpl * dbsi,Hscore * out,GenomePara* query,SnpMatch* target ) 
{
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_MouseSNPMatch. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_MouseSNPMatch. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 0 )  
      warn("Although you are asking at run-time for database tracing, the MouseSNPMatch matrix was not compiled with database tracing. No tracing will be made");    
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_MouseSNPMatch(out,query,target );   
      case DBSearchImpl_Pthreads :   
        warn("This matrix MouseSNPMatch was not dyc compiled with thread support");  
        return SEARCH_ERROR; 
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from MouseSNPMatch",impl_string_DBSearchImpl(dbsi));  
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  score_only_logsum_MouseSNPMatch(query,target)
 *
 * Descrip:    This function calculates the score over all paths
 *             This is using a logsum method to sort it all out
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Score]
 *
 */
Score score_only_logsum_MouseSNPMatch(GenomePara* query,SnpMatch* target ) 
{
    int i;   
    int j;   
    int bestscore = 0;   
    int k;   
    MouseSNPMatch * mat;     


    mat = allocate_MouseSNPMatch_only(query, target );   
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DISASTER!");   
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 0) * 3,2,2)) == NULL)  {  
      warn("Score only matrix for MouseSNPMatch cannot be allocated, (asking for 1  by %d  cells)",mat->leni*3); 
      mat = free_MouseSNPMatch(mat);     
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          MouseSNPMatch_VSMALL_MATRIX(mat,i,j,k) = NEGI; 
        }  
      MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,Start) = 0;   
      MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = NEGI;  
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;  
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;  
        score = Probability_logsum(score,temp);  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;    
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Musculus */ 
          temp = MouseSNPMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;  
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_VSMALL_MATRIX(mat,i,j,Musculus) = score;  


        /* state Musculus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = Probability_logsum(MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End),temp);    
          }  


        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;     
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;   
        score = Probability_logsum(score,temp);  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;    
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Cast */ 
          temp = MouseSNPMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;  
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_VSMALL_MATRIX(mat,i,j,Cast) = score;  


        /* state Cast is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = Probability_logsum(MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End),temp);    
          }  


        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;     
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;  
        score = Probability_logsum(score,temp);  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;     
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Domesticus */ 
          temp = MouseSNPMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;  
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_VSMALL_MATRIX(mat,i,j,Domesticus) = score;    


        /* state Domesticus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = Probability_logsum(MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End),temp);    
          }  


        /* Finished calculating state Domesticus */ 
        } /* end of for all query positions */ 




      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      } /* end of for all target positions */ 


    mat = free_MouseSNPMatch(mat);   
    return bestscore;    
}    


/* Function:  forward_logsum_MouseSNPMatch(query,target,dpri)
 *
 * Descrip:    This function calculates the matrix over all paths
 *             This is using a logsum method to sort it all out
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * forward_logsum_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPRunImpl * dpri) 
{
    MouseSNPMatch * mat; 
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    assert((mat=allocate_MouseSNPMatch_only(query, target )) != NULL);   
    mat->basematrix = BaseMatrix_alloc_matrix_specials_score_offset((mat->lenj+1)*3,(mat->leni+0),2,mat->lenj+1);    
    assert(mat->basematrix != NULL);     
    leni = mat->leni;    
    lenj = mat->lenj;    
    mat->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_MouseSNPMatch(mat);     
    tot = leni * lenj;   
    num = 0; 


    start_reporting("MouseSNPMatch Matrix calculation: ");   
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;    
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;    
        score = Probability_logsum(score,temp);  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;  
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Musculus */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = score;    


        /* state Musculus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = Probability_logsum(MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End),temp);    
          }  


        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;   
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;     
        score = Probability_logsum(score,temp);  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;  
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Cast */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = score;    


        /* state Cast is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = Probability_logsum(MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End),temp);    
          }  


        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;   
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;    
        score = Probability_logsum(score,temp);  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;   
        score = Probability_logsum(score,temp);  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Domesticus */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          score = Probability_logsum(score,temp);    
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = score;  


        /* state Domesticus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = Probability_logsum(MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End),temp);    
          }  


        /* Finished calculating state Domesticus */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return mat;  
}    


/* Function:  backward_logsum_MouseSNPMatch(query,target,dpri)
 *
 * Descrip:    This function calculates the matrix over all paths
 *             This is using a logsum method to sort it all out
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * backward_logsum_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPRunImpl * dpri) 
{
    MouseSNPMatch * mat; 
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 
    int max_score;   
    int min_score;   


    assert((mat=allocate_MouseSNPMatch_only(query, target )) != NULL);   
    mat->basematrix = BaseMatrix_alloc_matrix_specials_score_offset((mat->lenj+1)*3,(mat->leni+0),2,mat->lenj+1);    
    assert(mat->basematrix != NULL);     
    leni = mat->leni;    
    lenj = mat->lenj;    
    mat->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_MouseSNPMatch(mat);     
    tot = leni * lenj;   
    num = 0; 
    /* On the backward view, we need to have 0 at end, and NEGI at start */ 
    for(j= -1;j<lenj;j++)    {  
      MouseSNPMatch_EXPL_SPECIAL(mat,0,j,Start) = NEGI;  
      MouseSNPMatch_EXPL_SPECIAL(mat,0,j,End) = 0;   
      }  


    start_reporting("MouseSNPMatch Matrix calculation: ");   


    for(j=lenj-1;j>=0;j--)   {  
      auto int score;    
      auto int temp;     
      /* We look for underflow or overflow on this j column. */ 
      /* We do this at the start of the loop because in the backwards case at this point we know nothing will _update_ this column */ 
      min_score = SCORE_OVERFLOW;    
      max_score = SCORE_UNDERFLOW;   
      for(i=0;i<mat->leni;i++)   {  
        if( MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) > max_score) 
          max_score =  MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus);  
        if( MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) < min_score) 
          min_score =  MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus);  
        if( MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) > max_score) 
          max_score =  MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast);  
        if( MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) < min_score) 
          min_score =  MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast);  
        if( MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) > max_score)   
          max_score =  MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus);    
        if( MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) < min_score)   
          min_score =  MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus);    
        }  
      fprintf(stderr,"Before i position %d, max %d, min %d\n",max_score,min_score);
      if( MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start) > max_score) 
        max_score =  MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start);  
      fprintf(stderr,"Post start i position %d, max %d, min %d\n",max_score,min_score);
      if( MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start) < min_score) 
        min_score =  MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start);  
      if( MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) > max_score)   
        max_score =  MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End);    
      fprintf(stderr,"Post end i position %d, max %d, min %d %d\n",max_score,min_score,MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End));
        if( MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) < min_score)   
        min_score =  MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End);   


      fprintf(stderr,"At position %d, max %d, min %d\n",max_score,min_score);
      if( max_score > SCORE_OVERFLOW )   {  
        if( min_score < SCORE_UNDERFLOW) 
          fatal("Both overflow and underflow on the same column");   
        mat->basematrix->score_offset[j] = -1;   
        for(i=0;i<mat->leni;i++) {  
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) -= SCORE_OVERFLOW;     
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) -= SCORE_OVERFLOW;     
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) -= SCORE_OVERFLOW;   
          }  
        MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start) -= SCORE_OVERFLOW;     
        MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) -= SCORE_OVERFLOW;   
        }  
      if( min_score < SCORE_UNDERFLOW )  {  
	fprintf(stderr,"On position %d, underflowing\n");
        mat->basematrix->score_offset[j] = +1;   
        for(i=0;i<mat->leni;i++) {  
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) -= SCORE_UNDERFLOW;    
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) -= SCORE_UNDERFLOW;    
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) -= SCORE_UNDERFLOW;  
          }  
        MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start) -= SCORE_UNDERFLOW;    
        MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) -= SCORE_UNDERFLOW;  
        }  
      for(i=leni-1;i>=0;i--) {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   




        /* special state End recieves from Musculus */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) + 0 + (0);  
          MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-0,Musculus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-0,Musculus),temp);    
          }  
        /* special state End recieves from Cast */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) + 0 + (0);  
          MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-0,Cast) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-0,Cast),temp);    
          }  
        /* special state End recieves from Domesticus */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) + 0 + (0);  
          MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-0,Domesticus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-0,Domesticus),temp);    
          }  


        /* Doing calculations which end at this state */ 
        /* Reversed from state Musculus to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) + mat->query->extend + (TargetMatchMusculus(j));  
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus),temp);  
        /* Reversed from state Musculus to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) + mat->query->bswitch + (TargetMatchMusculus(j));     
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast),temp);  
        /* Reversed from state Musculus to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) + mat->query->bswitch + (TargetMatchMusculus(j));     
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus),temp);  
        /* Has restricted position */ 
        if( (j-0) == 0  )    {  
          /* Reversed from state Musculus to state Start */ 
          temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) + mat->query->bswitch + (TargetMatchMusculus(j));   
          MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) = Probability_logsum(MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start),temp);    
          }  


        /* Doing calculations which end at this state */ 
        /* Reversed from state Cast to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) + mat->query->bswitch + (TargetMatchCast(j));     
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus),temp);  
        /* Reversed from state Cast to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) + mat->query->extend + (TargetMatchCast(j));  
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast),temp);  
        /* Reversed from state Cast to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) + mat->query->bswitch + (TargetMatchCast(j));     
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus),temp);  
        /* Has restricted position */ 
        if( (j-0) == 0  )    {  
          /* Reversed from state Cast to state Start */ 
          temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) + mat->query->bswitch + (TargetMatchCast(j));   
          MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) = Probability_logsum(MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start),temp);    
          }  


        /* Doing calculations which end at this state */ 
        /* Reversed from state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) + mat->query->bswitch + (TargetMatchDomesticus(j));     
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus),temp);  
        /* Reversed from state Domesticus to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) + mat->query->bswitch + (TargetMatchDomesticus(j));     
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast),temp);  
        /* Reversed from state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) + mat->query->extend + (TargetMatchDomesticus(j));  
        MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) = Probability_logsum(MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus),temp);  
        /* Has restricted position */ 
        if( (j-0) == 0  )    {  
          /* Reversed from state Domesticus to state Start */ 
          temp = MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) + mat->query->bswitch + (TargetMatchDomesticus(j));   
          MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) = Probability_logsum(MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start),temp);    
          }  
        }  
      }  
    stop_reporting();    
    return mat;  
}    


/* Function:  serial_search_MouseSNPMatch(out,query,target)
 *
 * Descrip:    This function makes a database search of MouseSNPMatch
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [GenomePara*]
 * Arg:        target [UNKN ] Undocumented argument [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_MouseSNPMatch(Hscore * out,GenomePara* query,SnpMatch* target ) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_MouseSNPMatch(query, target );    
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("MouseSNPMatch search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");  
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


/* Function:  score_only_MouseSNPMatch(query,target)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_MouseSNPMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_MouseSNPMatch(GenomePara* query,SnpMatch* target ) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    MouseSNPMatch * mat;     


    mat = allocate_MouseSNPMatch_only(query, target );   
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 0) * 3,2,2)) == NULL)  {  
      warn("Score only matrix for MouseSNPMatch cannot be allocated, (asking for 1  by %d  cells)",mat->leni*3); 
      mat = free_MouseSNPMatch(mat);     
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          MouseSNPMatch_VSMALL_MATRIX(mat,i,j,k) = NEGI; 
        }  
      MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,Start) = 0;   
      MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = NEGI;  
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;  
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Musculus */ 
          temp = MouseSNPMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_VSMALL_MATRIX(mat,i,j,Musculus) = score;  


        /* state Musculus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) )     {  
            MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = temp;    
            }  


          }  


        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;     
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Cast */ 
          temp = MouseSNPMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_VSMALL_MATRIX(mat,i,j,Cast) = score;  


        /* state Cast is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) )     {  
            MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = temp;    
            }  


          }  


        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;     
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_VSMALL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Domesticus */ 
          temp = MouseSNPMatch_VSMALL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;  
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_VSMALL_MATRIX(mat,i,j,Domesticus) = score;    


        /* state Domesticus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) )     {  
            MouseSNPMatch_VSMALL_SPECIAL(mat,i,j,End) = temp;    
            }  


          }  


        /* Finished calculating state Domesticus */ 
        } /* end of for all query positions */ 




      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      if( bestscore < MouseSNPMatch_VSMALL_SPECIAL(mat,0,j,End) )    
        bestscore = MouseSNPMatch_VSMALL_SPECIAL(mat,0,j,End);   
      } /* end of for all target positions */ 


    mat = free_MouseSNPMatch(mat);   
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_MouseSNPMatch(query,target,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_MouseSNPMatch
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    MouseSNPMatch * mat; 
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->len * target->len;    
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 3 * sizeof(int)) > 1000*dpri->kbyte_size) {  
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
        if( (mat=allocate_Expl_MouseSNPMatch(query, target ,dpri)) == NULL ) {  
          warn("Unable to allocate large MouseSNPMatch version");    
          return NULL;   
          }  
        calculate_dpenv_MouseSNPMatch(mat,dpenv);    
        out =  PackAln_read_Expl_MouseSNPMatch(mat); 
        }  
      else   {  
        mat = allocate_MouseSNPMatch_only(query, target );   
        calculate_shatter_MouseSNPMatch(mat,dpenv);  
        out = PackAln_read_Shatter_MouseSNPMatch(mat);   
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_MouseSNPMatch(query, target )) == NULL ) {  
          warn("Unable to allocate small MouseSNPMatch version");    
          return NULL;   
          }  
        out = PackAln_calculate_Small_MouseSNPMatch(mat,dpenv);  
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_MouseSNPMatch(query, target ,dpri)) == NULL ) {  
          warn("Unable to allocate large MouseSNPMatch version");    
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_MouseSNPMatch(mat);  
          out =  PackAln_read_Expl_MouseSNPMatch(mat);   
          }  
        }  
      }  


    mat = free_MouseSNPMatch(mat);   
    return out;  
}    


/* Function:  allocate_MouseSNPMatch_only(query,target)
 *
 * Descrip:    This function only allocates the MouseSNPMatch structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * allocate_MouseSNPMatch_only(GenomePara* query,SnpMatch* target ) 
{
    MouseSNPMatch * out;     


    if((out= MouseSNPMatch_alloc()) == NULL) {  
      warn("Allocation of basic MouseSNPMatch structure failed..."); 
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->leni = query->len;  
    out->lenj = target->len;     
    return out;  
}    


/* Function:  allocate_Expl_MouseSNPMatch(query,target,dpri)
 *
 * Descrip:    This function allocates the MouseSNPMatch structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_MouseSNPMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * allocate_Expl_MouseSNPMatch(GenomePara* query,SnpMatch* target ,DPRunImpl * dpri) 
{
    MouseSNPMatch * out; 


    out = allocate_MouseSNPMatch_only(query, target );   
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*3 && dpri->cache->maxlenj >= (out->leni+0))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*3,(out->leni+0),2,out->lenj+1)) == NULL) {  
        warn("Explicit matrix MouseSNPMatch cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj); 
        free_MouseSNPMatch(out);     
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_MouseSNPMatch(out);     
    return out;  
}    


/* Function:  init_MouseSNPMatch(mat)
 *
 * Descrip:    This function initates MouseSNPMatch matrix when in explicit mode
 *             Called in /allocate_Expl_MouseSNPMatch
 *
 *
 * Arg:        mat [UNKN ] MouseSNPMatch which contains explicit basematrix memory [MouseSNPMatch *]
 *
 */
void init_MouseSNPMatch(MouseSNPMatch * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-0);i<mat->query->len;i++)   {  
      for(j= (-1);j<2;j++)   {  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = NEGI;  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = NEGI;  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = NEGI;    
        }  
      }  
    for(j= (-1);j<mat->target->len;j++)  {  
      for(i= (-0);i<1;i++)   {  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = NEGI;  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = NEGI;  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = NEGI;    
        }  
      MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start) = 0; 
      MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = NEGI;    
      }  
    return;  
}    


/* Function:  recalculate_PackAln_MouseSNPMatch(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by MouseSNPMatch
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 */
void recalculate_PackAln_MouseSNPMatch(PackAln * pal,MouseSNPMatch * mat) 
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
        case Musculus :  
          if( offi == 0 && offj == 1 && prev->state == Musculus )    {  
            pau->score = mat->query->extend + (TargetMatchMusculus(j));  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == Cast )    {  
            pau->score = mat->query->bswitch + (TargetMatchMusculus(j));     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == Domesticus )  {  
            pau->score = mat->query->bswitch + (TargetMatchMusculus(j));     
            continue;    
            }  
          if( offj == 1 && prev->state == (Start+3) )    {  
            pau->score = mat->query->bswitch + (TargetMatchMusculus(j));     
            continue;    
            }  
          warn("In recaluclating PackAln with state Musculus, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case Cast :  
          if( offi == 0 && offj == 1 && prev->state == Musculus )    {  
            pau->score = mat->query->bswitch + (TargetMatchCast(j));     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == Cast )    {  
            pau->score = mat->query->extend + (TargetMatchCast(j));  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == Domesticus )  {  
            pau->score = mat->query->bswitch + (TargetMatchCast(j));     
            continue;    
            }  
          if( offj == 1 && prev->state == (Start+3) )    {  
            pau->score = mat->query->bswitch + (TargetMatchCast(j));     
            continue;    
            }  
          warn("In recaluclating PackAln with state Cast, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case Domesticus :    
          if( offi == 0 && offj == 1 && prev->state == Musculus )    {  
            pau->score = mat->query->bswitch + (TargetMatchDomesticus(j));   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == Cast )    {  
            pau->score = mat->query->bswitch + (TargetMatchDomesticus(j));   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == Domesticus )  {  
            pau->score = mat->query->extend + (TargetMatchDomesticus(j));    
            continue;    
            }  
          if( offj == 1 && prev->state == (Start+3) )    {  
            pau->score = mat->query->bswitch + (TargetMatchDomesticus(j));   
            continue;    
            }  
          warn("In recaluclating PackAln with state Domesticus, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case (Start+3) :     
          warn("In recaluclating PackAln with state Start, got a bad source state. Error!"); 
          break; 
        case (End+3) :   
          if( offj == 0 && prev->state == Musculus ) {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == Domesticus )   {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offj == 0 && prev->state == Cast ) {  
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
#define MouseSNPMatch_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+0)*3+state]) 
#define MouseSNPMatch_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+0)*3+state]) 
#define MouseSNPMatch_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])    
#define MouseSNPMatch_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])   
#define MouseSNPMatch_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+0)*3 + state])  
#define MouseSNPMatch_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])   
#define MouseSNPMatch_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 3) + ((i+0) * 3) + (state)])   
#define MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 24) + ((i+0) * 24) + (state * 8) + shadow+1])   
#define MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])   
/* Function:  allocate_Small_MouseSNPMatch(query,target)
 *
 * Descrip:    This function allocates the MouseSNPMatch structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_MouseSNPMatch_only
 *
 *
 * Arg:         query [UNKN ] query data structure [GenomePara*]
 * Arg:        target [UNKN ] target data structure [SnpMatch*]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
#define MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])   
MouseSNPMatch * allocate_Small_MouseSNPMatch(GenomePara* query,SnpMatch* target ) 
{
    MouseSNPMatch * out; 


    out = allocate_MouseSNPMatch_only(query, target );   
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 0) * 3,16,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix MouseSNPMatch cannot be allocated, (asking for 2 by %d main cells)",out->leni+1);    
      free_MouseSNPMatch(out);   
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_MouseSNPMatch(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for MouseSNPMatch structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_MouseSNPMatch 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_MouseSNPMatch 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_MouseSNPMatch(MouseSNPMatch * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for MouseSNPMatch due to wrong type of matrix");   
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_MouseSNPMatch(mat,dpenv);    
    score = start_end_find_end_MouseSNPMatch(mat,&endj); 
    out->score = score;  
    stopstate = End;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,0);   
    startj = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,1);   
    startstate = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,2);   
    stopi = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,3);    
    stopj = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,4);    
    stopstate = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,5);    
    temp = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,endj,End,6); 
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);   


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_MouseSNPMatch(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);  
    if( stopi == MouseSNPMatch_READ_OFF_ERROR || stopstate != Start )    {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_MouseSNPMatch(out,mat);    
      return out;    
      }  


    /* Ok. Put away start start... */ 
    pau = PackAlnUnit_alloc();   
    pau->i = stopi;  
    pau->j = stopj;  
    pau->state = stopstate + 3;  
    add_PackAln(out,pau);    


    log_full_error(REPORT,0,"Alignment recovered");  
    stop_reporting();    
    invert_PackAln(out); 
    recalculate_PackAln_MouseSNPMatch(out,mat);  
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_MouseSNPMatch(mat)
 *
 * Descrip:    This function calculates an alignment for MouseSNPMatch structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_MouseSNPMatch 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_MouseSNPMatch
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_MouseSNPMatch 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_MouseSNPMatch(MouseSNPMatch * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_MouseSNPMatch(mat,NULL); 
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_MouseSNPMatch(mat);   
    return out;  
}    


/* Function:  AlnRangeSet_from_MouseSNPMatch(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for MouseSNPMatch structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_MouseSNPMatch 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_MouseSNPMatch
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_MouseSNPMatch(MouseSNPMatch * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_MouseSNPMatch"); 
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_MouseSNPMatch(mat,&jpos);    
    state = End; 


    while( (temp = AlnRange_build_MouseSNPMatch(mat,jpos,state,&jpos,&state)) != NULL)   
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_MouseSNPMatch(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_MouseSNPMatch
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_MouseSNPMatch(MouseSNPMatch * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_MouseSNPMatch"); 
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_MouseSNPMatch(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)  {  
      warn("In AlnRanger_build_MouseSNPMatch alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);  
      return NULL;   
      }  
    if( state == Start || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);    
    out->startj = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);    
    out->startstate = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);    
    out->stopi = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3); 
    out->stopj = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4); 
    out->stopstate = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5); 
    out->startscore = MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);    
    out->stopscore = MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,jpos,state);  


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_MouseSNPMatch(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);    
    if( jpos == MouseSNPMatch_READ_OFF_ERROR)    {  
      warn("In AlnRange_build_MouseSNPMatch alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);   
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);    
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
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
boolean read_hidden_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_MouseSNPMatch(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   


      if( i == MouseSNPMatch_READ_OFF_ERROR) {  
        warn("In MouseSNPMatch hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj); 
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
        warn("In MouseSNPMatch hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);   
        return FALSE;    
        }  
      }  
    warn("In MouseSNPMatch hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);    
    return FALSE;    
}    


/* Function:  max_hidden_MouseSNPMatch(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
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
int max_hidden_MouseSNPMatch(MouseSNPMatch * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = MouseSNPMatch_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->len)   {  
      warn("In MouseSNPMatch matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,state); 
    switch(state)    { /*Switch state */ 
      case Musculus :    
        /* Not allowing special sources.. skipping Start */ 
        temp = cscore - (mat->query->bswitch) -  (TargetMatchMusculus(j));   
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Domesticus) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Domesticus;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Domesticus);   
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Domesticus);    
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchMusculus(j));   
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Cast) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Cast;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Cast); 
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Cast);  
          }  
        temp = cscore - (mat->query->extend) -  (TargetMatchMusculus(j));    
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Musculus) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Musculus;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Musculus); 
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Musculus);  
          }  
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case Cast :    
        /* Not allowing special sources.. skipping Start */ 
        temp = cscore - (mat->query->bswitch) -  (TargetMatchCast(j));   
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Domesticus) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Domesticus;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Domesticus);   
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Domesticus);    
          }  
        temp = cscore - (mat->query->extend) -  (TargetMatchCast(j));    
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Cast) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Cast;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Cast); 
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Cast);  
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchCast(j));   
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Musculus) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Musculus;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Musculus); 
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Musculus);  
          }  
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case Domesticus :  
        /* Not allowing special sources.. skipping Start */ 
        temp = cscore - (mat->query->extend) -  (TargetMatchDomesticus(j));  
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Domesticus) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Domesticus;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Domesticus);   
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Domesticus);    
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchDomesticus(j)); 
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Cast) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Cast;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Cast); 
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Cast);  
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchDomesticus(j)); 
        if( temp == MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Musculus) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Musculus;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Musculus); 
            }  
          return MouseSNPMatch_HIDDEN_MATRIX(mat,i - 0,j - 1,Musculus);  
          }  
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_MouseSNPMatch(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
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
boolean read_special_strip_MouseSNPMatch(MouseSNPMatch * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != Start)    { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 3; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_MouseSNPMatch(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);   
      if( i == MouseSNPMatch_READ_OFF_ERROR) {  
        warn("In special strip read MouseSNPMatch, got a bad read off error. Sorry!");   
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != Start && j < MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))   {  
      warn("In special strip read MouseSNPMatch, at special [%d] state [%d] overshot!",j,state); 
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 3;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_MouseSNPMatch(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
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
int max_special_strip_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = MouseSNPMatch_READ_OFF_ERROR;    
    if( isspecial == FALSE ) {  
      warn("In special strip max function for MouseSNPMatch, got a non special start point. Problem! (bad!)");   
      return (-1);   
      }  


    if( j < 0 || j > mat->target->len)   {  
      warn("In MouseSNPMatch matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1; 
      }  


    cscore = MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case Start :   
      case End :     
        /* Source Cast is not a special */ 
        /* Source Domesticus is not a special */ 
        /* Source Musculus is not a special */ 
      default:   
        warn("Major problem (!) - in MouseSNPMatch special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_MouseSNPMatch(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
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
int max_matrix_to_special_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = MouseSNPMatch_READ_OFF_ERROR;    


    if( j < 0 || j > mat->lenj)  {  
      warn("In MouseSNPMatch matrix to special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case Musculus :    
        temp = cscore - (mat->query->bswitch) -  (TargetMatchMusculus(j));   
        if( temp == MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,Start) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Start; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i-0,j-1,Start);    
            }  
          return MouseSNPMatch_DC_SHADOW_MATRIX(mat,i - 0,j - 1,Start) ;     
          }  
        /* Source Domesticus is not a special, should not get here! */ 
        /* Source Cast is not a special, should not get here! */ 
        /* Source Musculus is not a special, should not get here! */ 
        warn("Major problem (!) - in MouseSNPMatch matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case Cast :    
        temp = cscore - (mat->query->bswitch) -  (TargetMatchCast(j));   
        if( temp == MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,Start) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Start; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i-0,j-1,Start);    
            }  
          return MouseSNPMatch_DC_SHADOW_MATRIX(mat,i - 0,j - 1,Start) ;     
          }  
        /* Source Domesticus is not a special, should not get here! */ 
        /* Source Cast is not a special, should not get here! */ 
        /* Source Musculus is not a special, should not get here! */ 
        warn("Major problem (!) - in MouseSNPMatch matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case Domesticus :  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchDomesticus(j));     
        if( temp == MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,Start) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Start; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - MouseSNPMatch_DC_SHADOW_SPECIAL(mat,i-0,j-1,Start);    
            }  
          return MouseSNPMatch_DC_SHADOW_MATRIX(mat,i - 0,j - 1,Start) ;     
          }  
        /* Source Domesticus is not a special, should not get here! */ 
        /* Source Cast is not a special, should not get here! */ 
        /* Source Musculus is not a special, should not get here! */ 
        warn("Major problem (!) - in MouseSNPMatch matrix to special read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_MouseSNPMatch(mat,starti,startj,stopi,stopj);    


    MouseSNPMatch_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;   


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Musculus) = NEGI;  
          MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Cast) = NEGI;  
          MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Domesticus) = NEGI;    
          continue;  
          } /* end of Is not in envelope */ 


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;  
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Musculus) = score;  
        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;     
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Cast) = score;  
        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;     
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_HIDDEN_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Domesticus) = score;    
        /* Finished calculating state Domesticus */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_MouseSNPMatch(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-0);i<=stopi;i++) {  
        MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Musculus) = NEGI;
   
        MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Cast) = NEGI;
   
        MouseSNPMatch_HIDDEN_MATRIX(mat,i,j,Domesticus) = NEGI;
 
        }  
      }  


    return;  
}    


/* Function:  full_dc_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_MouseSNPMatch
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_MouseSNPMatch to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [MouseSNPMatch *]
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
boolean full_dc_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_MouseSNPMatch");   
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);  
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)    {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE) {  
      warn("In divide and conquor for MouseSNPMatch, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);    
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);     
    lstartj= MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);     
    lstate = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);     


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_MouseSNPMatch(mat,MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_MouseSNPMatch(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)  {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_MouseSNPMatch(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
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
boolean do_dc_single_pass_MouseSNPMatch(MouseSNPMatch * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_MouseSNPMatch(mat);  


    MouseSNPMatch_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;    
    run_up_dc_MouseSNPMatch(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);    
    push_dc_at_merge_MouseSNPMatch(mat,starti,stopi,halfj,&halfj,dpenv);     
    follow_on_dc_MouseSNPMatch(mat,starti,stopi,halfj,stopj,dpenv,perc_done);    
    return TRUE; 
}    


/* Function:  push_dc_at_merge_MouseSNPMatch(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_MouseSNPMatch(MouseSNPMatch * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = NEGI;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,0) = (-100);    
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,1) = (-100);    
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = NEGI;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,0) = (-100);    
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,1) = (-100);    
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = NEGI;     
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,0) = (-100);  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,1) = (-100);  
          continue;  
          } /* end of Is not in envelope */ 


        /* For state Musculus, pushing when j - offj <= mergej */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;   
        if( j - 1 <= mergej) {  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,0) = i-0;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,1) = j-1;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,2) = Musculus;  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,3) = i; 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,4) = j; 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,5) = Musculus;  
          }  
        else {  
          for(k=0;k<7;k++)   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,k);   
          }  


        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,0) = i-0; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,1) = j-1; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,2) = Cast;    
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,3) = i;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,4) = j;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,5) = Musculus;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,k); 
            }  
          }  


        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,0) = i-0; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,1) = j-1; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,2) = Domesticus;  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,3) = i;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,4) = j;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,5) = Musculus;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += TargetMatchMusculus(j);     
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = score;    
        /* Finished with state Musculus */ 


        /* For state Cast, pushing when j - offj <= mergej */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;  
        if( j - 1 <= mergej) {  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,0) = i-0;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,1) = j-1;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,2) = Musculus;  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,3) = i; 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,4) = j; 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,5) = Cast;  
          }  
        else {  
          for(k=0;k<7;k++)   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,k);   
          }  


        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,0) = i-0; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,1) = j-1; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,2) = Cast;    
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,3) = i;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,4) = j;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,5) = Cast;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,k); 
            }  
          }  


        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,0) = i-0; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,1) = j-1; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,2) = Domesticus;  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,3) = i;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,4) = j;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,5) = Cast;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,k);   
            }  
          }  
        /* Add any movement independant score */ 
        score += TargetMatchCast(j);     
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = score;    
        /* Finished with state Cast */ 


        /* For state Domesticus, pushing when j - offj <= mergej */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;  
        if( j - 1 <= mergej) {  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,0) = i-0; 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,1) = j-1; 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,2) = Musculus;    
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,3) = i;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,4) = j;   
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,5) = Domesticus;  
          }  
        else {  
          for(k=0;k<7;k++)   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,k); 
          }  


        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,0) = i-0;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,1) = j-1;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,2) = Cast;  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,3) = i; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,4) = j; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,5) = Domesticus;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,k);   
            }  
          }  


        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,0) = i-0;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,1) = j-1;   
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,2) = Domesticus;    
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,3) = i; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,4) = j; 
            MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,5) = Domesticus;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k) = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += TargetMatchDomesticus(j);   
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = score;  
        /* Finished with state Domesticus */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_MouseSNPMatch(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_MouseSNPMatch(MouseSNPMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = NEGI;   
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = NEGI;   
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,k);    
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,k);  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,k);    
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = score;   
        for(k=0;k<7;k++) 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,k);    
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,k);  
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,k);    
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = score;   
        for(k=0;k<7;k++) 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,k);    
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,k);  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,k);    
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = score; 
        for(k=0;k<7;k++) 
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state Domesticus */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_MouseSNPMatch(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_MouseSNPMatch(MouseSNPMatch * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = NEGI;   
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = NEGI;   
          MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;   
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = score;   
        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;  
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = score;   
        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;  
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_DC_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = score; 
        /* Finished calculating state Domesticus */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_MouseSNPMatch(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 */
}    
void init_dc_MouseSNPMatch(MouseSNPMatch * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->query->len;i++)  {  
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = NEGI; 
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = NEGI; 
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = NEGI;   
        for(k=0;k<7;k++) {  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,k) = (-1);  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,k) = (-1);  
          MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k) = (-1);    
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_MouseSNPMatch(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [MouseSNPMatch *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_MouseSNPMatch(MouseSNPMatch * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,mat->target->len-1,End); 
    maxj = mat->target->len-1;   
    for(j= mat->target->len-2 ;j >= 0 ;j--)  {  
      if( MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,j,End) > max )   {  
        max = MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,j,End);  
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_MouseSNPMatch(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [MouseSNPMatch]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_MouseSNPMatch(MouseSNPMatch *mat,DPEnvelope * dpenv) 
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


    score_pointers = (int *) calloc (1 * (leni + 0) * 3,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 0) * 3 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,Musculus) = NEGI;   
          MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,Cast) = NEGI;   
          MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,Domesticus) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend + (TargetMatchMusculus(j));    
        /* assign local shadown pointer */ 
        localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,0));  
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch +(TargetMatchMusculus(j));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,0));    
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch +(TargetMatchMusculus(j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,0));  
          }  
        /* From state Start to state Musculus */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch + (TargetMatchMusculus(j));  
        if( temp  > score )  {  
          score = temp;  
          /* This state [Start] is a special for Musculus... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= Musculus;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,Musculus) = score;   
        for(k=0;k<7;k++) 
          MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Musculus,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state Musculus is a source for special End */ 
        temp = score + (0) + (0) ;   
        if( temp > MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) )    {  
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,k) = MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Musculus,k);   
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,6) = MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Musculus,6); 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,3) = i; 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,4) = j; 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,5) = Musculus;  
          }  




        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch + (TargetMatchCast(j));   
        /* assign local shadown pointer */ 
        localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,0));  
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend +(TargetMatchCast(j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,0));    
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch +(TargetMatchCast(j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,0));  
          }  
        /* From state Start to state Cast */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch + (TargetMatchCast(j));  
        if( temp  > score )  {  
          score = temp;  
          /* This state [Start] is a special for Cast... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= Cast;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,Cast) = score;   
        for(k=0;k<7;k++) 
          MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Cast,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state Cast is a source for special End */ 
        temp = score + (0) + (0) ;   
        if( temp > MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) )    {  
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,k) = MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Cast,k);   
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,6) = MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Cast,6); 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,3) = i; 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,4) = j; 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,5) = Cast;  
          }  




        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch + (TargetMatchDomesticus(j));     
        /* assign local shadown pointer */ 
        localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Musculus,0));  
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch +(TargetMatchDomesticus(j));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Cast,0));    
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend +(TargetMatchDomesticus(j));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,Domesticus,0));  
          }  
        /* From state Start to state Domesticus */ 
        temp = MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch + (TargetMatchDomesticus(j));    
        if( temp  > score )  {  
          score = temp;  
          /* This state [Start] is a special for Domesticus... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= Domesticus;    
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         MouseSNPMatch_DC_OPT_SHADOW_MATRIX(mat,i,j,Domesticus) = score; 
        for(k=0;k<7;k++) 
          MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* state Domesticus is a source for special End */ 
        temp = score + (0) + (0) ;   
        if( temp > MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) )    {  
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL(mat,i,j,End) = temp;   
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,k) = MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Domesticus,k); 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,6) = MouseSNPMatch_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,Domesticus,6);   
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,3) = i; 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,4) = j; 
          MouseSNPMatch_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,End,5) = Domesticus;    
          }  




        /* Finished calculating state Domesticus */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_MouseSNPMatch(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 */
void init_start_end_linear_MouseSNPMatch(MouseSNPMatch * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-0);i<mat->query->len;i++)  {  
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Musculus) = NEGI; 
        MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Musculus,0) = (-1);    
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Cast) = NEGI; 
        MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Cast,0) = (-1);    
        MouseSNPMatch_DC_SHADOW_MATRIX(mat,i,j,Domesticus) = NEGI;   
        MouseSNPMatch_DC_SHADOW_MATRIX_SP(mat,i,j,Domesticus,0) = (-1);  
        }  
      }  


    for(j=(-1);j<mat->target->len;j++)   {  
      MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,j,Start) = 0;    
      MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,j,Start,0) = j;   
      MouseSNPMatch_DC_SHADOW_SPECIAL(mat,0,j,End) = NEGI;   
      MouseSNPMatch_DC_SHADOW_SPECIAL_SP(mat,0,j,End,0) = (-1);  
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_MouseSNPMatch(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_MouseSNPMatch(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_MouseSNPMatch(); 
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "Q","END" };    
/* Function:  AlnConvertSet_MouseSNPMatch(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "MUSCULUS","CAST","DOMESTICUS","END" };    
AlnConvertSet * AlnConvertSet_MouseSNPMatch(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Musculus;  
    acu->state2 = Musculus;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Cast;  
    acu->state2 = Musculus;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Domesticus;    
    acu->state2 = Musculus;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Start + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = Musculus;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Musculus;  
    acu->state2 = Cast;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Cast;  
    acu->state2 = Cast;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Domesticus;    
    acu->state2 = Cast;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Start + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = Cast;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Musculus;  
    acu->state2 = Domesticus;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Cast;  
    acu->state2 = Domesticus;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Domesticus;    
    acu->state2 = Domesticus;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Start + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = Domesticus;    
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Musculus;  
    acu->state2 = End + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Domesticus;    
    acu->state2 = End + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = Cast;  
    acu->state2 = End + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[3];   
    return out;  
}    


/* Function:  PackAln_read_Expl_MouseSNPMatch(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_MouseSNPMatch(MouseSNPMatch * mat) 
{
    MouseSNPMatch_access_func_holder holder;     


    holder.access_main    = MouseSNPMatch_explicit_access_main;  
    holder.access_special = MouseSNPMatch_explicit_access_special;   
    return PackAln_read_generic_MouseSNPMatch(mat,holder);   
}    


/* Function:  MouseSNPMatch_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MouseSNPMatch_explicit_access_main(MouseSNPMatch * mat,int i,int j,int state) 
{
    return MouseSNPMatch_EXPL_MATRIX(mat,i,j,state); 
}    


/* Function:  MouseSNPMatch_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int MouseSNPMatch_explicit_access_special(MouseSNPMatch * mat,int i,int j,int state) 
{
    return MouseSNPMatch_EXPL_SPECIAL(mat,i,j,state);    
}    


/* Function:  PackAln_read_generic_MouseSNPMatch(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:          h [UNKN ] Undocumented argument [MouseSNPMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_MouseSNPMatch(MouseSNPMatch * mat,MouseSNPMatch_access_func_holder h) 
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


    out->score =  find_end_MouseSNPMatch(mat,&i,&j,&state,&isspecial,h); 


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in MouseSNPMatch_basic_read, returning a mess.(Sorry!)",out->len);    
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 3;     
    prev=pau;    
    while( state != Start || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_MouseSNPMatch(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      else   
        max_calc_MouseSNPMatch(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);    
      if(i == MouseSNPMatch_READ_OFF_ERROR || j == MouseSNPMatch_READ_OFF_ERROR || state == MouseSNPMatch_READ_OFF_ERROR )   {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in MouseSNPMatch_basic_read, returning partial alignment",out->len);    
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 3;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_MouseSNPMatch(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [MouseSNPMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_MouseSNPMatch(MouseSNPMatch * mat,int * ri,int * rj,int * state,boolean * isspecial,MouseSNPMatch_access_func_holder h) 
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


/* Function:  MouseSNPMatch_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void MouseSNPMatch_debug_show_matrix(MouseSNPMatch * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->len;i++) {  
      for(j=startj;j<stopj && j < mat->target->len;j++)  {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State Musculus %d\n",MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus));  
        fprintf(ofp,"State Cast %d\n",MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast));  
        fprintf(ofp,"State Domesticus %d\n",MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus));  
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_MouseSNPMatch(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [MouseSNPMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MouseSNPMatch_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = MouseSNPMatch_READ_OFF_ERROR;    


    if( i < 0 || j < 0 || i > mat->query->len || j > mat->target->len)   {  
      warn("In MouseSNPMatch matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);   
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case Musculus :    
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          temp = cscore - (mat->query->bswitch) -  (TargetMatchMusculus(j)); 
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
        temp = cscore - (mat->query->bswitch) -  (TargetMatchMusculus(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Domesticus) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Domesticus;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Domesticus);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Domesticus);   
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchMusculus(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Cast) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Cast;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Cast);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Cast);     
          }  
        temp = cscore - (mat->query->extend) -  (TargetMatchMusculus(j));    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Musculus) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Musculus;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Musculus);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Musculus);     
          }  
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case Cast :    
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          temp = cscore - (mat->query->bswitch) -  (TargetMatchCast(j)); 
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
        temp = cscore - (mat->query->bswitch) -  (TargetMatchCast(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Domesticus) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Domesticus;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Domesticus);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Domesticus);   
          }  
        temp = cscore - (mat->query->extend) -  (TargetMatchCast(j));    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Cast) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Cast;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Cast);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Cast);     
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchCast(j));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Musculus) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Musculus;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Musculus);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Musculus);     
          }  
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case Domesticus :  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          temp = cscore - (mat->query->bswitch) -  (TargetMatchDomesticus(j));   
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
        temp = cscore - (mat->query->extend) -  (TargetMatchDomesticus(j));  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Domesticus) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Domesticus;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Domesticus);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Domesticus);   
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchDomesticus(j)); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Cast) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Cast;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Cast);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Cast);     
          }  
        temp = cscore - (mat->query->bswitch) -  (TargetMatchDomesticus(j)); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,Musculus) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = Musculus;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,Musculus);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,Musculus);     
          }  
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_MouseSNPMatch(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [MouseSNPMatch *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [MouseSNPMatch_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_MouseSNPMatch(MouseSNPMatch * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,MouseSNPMatch_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = MouseSNPMatch_READ_OFF_ERROR;    


    if( j < 0 || j > mat->target->len)   {  
      warn("In MouseSNPMatch matrix special read off - out of bounds on matrix [j is %d in special]",j); 
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case Start :   
      case End :     
        /* source Cast is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,Cast) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = Cast;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,Cast);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,Cast) ;  
            }  
          } /* end of for i >= 0 */ 
        /* source Domesticus is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,Domesticus) ) {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = Domesticus;  
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,Domesticus);    
              }  
            return (*h.access_main)(mat,i - 0,j - 0,Domesticus) ;    
            }  
          } /* end of for i >= 0 */ 
        /* source Musculus is from main matrix */ 
        for(i= mat->query->len-1;i >= 0 ;i--)    { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,Musculus) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = Musculus;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,Musculus);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,Musculus) ;  
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in MouseSNPMatch read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_MouseSNPMatch(mat)
 *
 * Descrip:    This function calculates the MouseSNPMatch matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_MouseSNPMatch
 *
 *
 * Arg:        mat [UNKN ] MouseSNPMatch which contains explicit basematrix memory [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_MouseSNPMatch(MouseSNPMatch * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_MouseSNPMatch, passed a non Explicit matrix type, cannot calculate!");  
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("MouseSNPMatch Matrix calculation: ");   
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;    
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Musculus */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = score;    


        /* state Musculus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) )   {  
            MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;   
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Cast */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = score;    


        /* state Cast is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) )   {  
            MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;   
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Domesticus */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = score;  


        /* state Domesticus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) )   {  
            MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state Domesticus */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_MouseSNPMatch(mat,dpenv)
 *
 * Descrip:    This function calculates the MouseSNPMatch matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] MouseSNPMatch which contains explicit basematrix memory [MouseSNPMatch *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_MouseSNPMatch(MouseSNPMatch * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_MouseSNPMatch, passed a non Explicit matrix type, cannot calculate!");  
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
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = NEGI;  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = NEGI;  
        MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = NEGI;    
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      MouseSNPMatch_EXPL_SPECIAL(mat,i,j,Start) = 0; 
      MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = NEGI;    
      }  


    start_reporting("MouseSNPMatch Matrix calculation: ");   
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
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = NEGI;    
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = NEGI;    
          MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = NEGI;  
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state Musculus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->extend;    
        /* From state Cast to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Musculus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Musculus */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Musculus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchMusculus(j);    
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Musculus) = score;    


        /* state Musculus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) )   {  
            MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state Musculus */ 


        /* For state Cast */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;   
        /* From state Cast to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->extend;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Cast */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->bswitch;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Cast */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Cast */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchCast(j);    
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Cast) = score;    


        /* state Cast is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) )   {  
            MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state Cast */ 


        /* For state Domesticus */ 
        /* setting first movement to score */ 
        score = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Musculus) + mat->query->bswitch;   
        /* From state Cast to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Cast) + mat->query->bswitch;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state Domesticus to state Domesticus */ 
        temp = MouseSNPMatch_EXPL_MATRIX(mat,i-0,j-1,Domesticus) + mat->query->extend;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* Has restricted position */ 
        if( (j-1) == 0  )    {  
          /* From state Start to state Domesticus */ 
          temp = MouseSNPMatch_EXPL_SPECIAL(mat,i-0,j-1,Start) + mat->query->bswitch;    
          if( temp  > score )    {  
            score = temp;    
            }  
          }  


        /* Ok - finished max calculation for Domesticus */ 
        /* Add any movement independant score and put away */ 
         score += TargetMatchDomesticus(j);  
         MouseSNPMatch_EXPL_MATRIX(mat,i,j,Domesticus) = score;  


        /* state Domesticus is a source for special End */ 
        /* Has restricted position */ 
        if( j == mat->lenj-1 )   {  
          temp = score + (0) + (0) ;     
          if( temp > MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) )   {  
            MouseSNPMatch_EXPL_SPECIAL(mat,i,j,End) = temp;  
            }  


          }  


        /* Finished calculating state Domesticus */ 
        }  


      /* Special state Start has no special to special movements */ 


      /* Special state End has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  MouseSNPMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * MouseSNPMatch_alloc(void) 
{
    MouseSNPMatch * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MouseSNPMatch *) ckalloc (sizeof(MouseSNPMatch))) == NULL)  {  
      warn("MouseSNPMatch_alloc failed ");   
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


/* Function:  free_MouseSNPMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MouseSNPMatch *]
 *
 * Return [UNKN ]  Undocumented return value [MouseSNPMatch *]
 *
 */
MouseSNPMatch * free_MouseSNPMatch(MouseSNPMatch * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MouseSNPMatch obj. Should be trappable"); 
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
