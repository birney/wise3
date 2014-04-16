#ifdef _cplusplus
extern "C" {
#endif
#include "locus_model_estimators.h"


# line 43 "locus_model_estimators.dy"
ModelProbabilitySet * estimate_model_ModelProbabilitySet(LocusGenotypeCount * lgc,double sweep,int points,LocusModelType type,double sel_sweep,int sweep_points)
{
  int i;
  int pop = lgc->len-1;
  int pop_point[MAX_POP];
  double central[MAX_POP];
  double incr;
  ModelProbabilityPoint mpp;
  ModelProbabilityPoint submpp;
  
  int overflow;
  ModelProbabilitySet * mps;

  double sweep_incr;
  int ii;
  int k;

  mps = ModelProbabilitySet_alloc_std();
  mps->keep = 0;
  mps->best = ModelProbabilityPoint_alloc();
  mps->best->likelihood_score = NEGI;

  assert(lgc != NULL);
  assert(lgc->len > 0);
  assert(sweep > 0 && sweep < 1.0);
  assert(points > 2);


  /*  fprintf(stderr,"Looking at locus %s with length %d\n",lgc->bl->locus_id,lgc->len);*/

  for(i=0;i<lgc->len;i++) {
    pop_point[i] = 0;
    central[i] = central_first_frequency_PopulationGenotypeCount(lgc->pgc[i]);
  }

  incr = sweep*2 / points;
  sweep_incr = sel_sweep*2 / sweep_points;

  while( 1 ) {
    /* we test this position first */


    overflow = 0;

    for(i=0;i<lgc->len;i++) {
      mpp.estimate_first_freq[i] = 
	central[i] - sweep + (pop_point[i]*incr);
      if( mpp.estimate_first_freq[i] > 1.0 || mpp.estimate_first_freq[i] < 0.0 ) {
	overflow = 1;
      }
    }

    if( overflow == 0 ) {

      if( type == LocusModel_Normal_Multi ) {
	estimate_normal_multi(lgc,&mpp);
	add_and_update_ModelProbabilitySet(mps,&mpp);
      } else {
	for(k=0;k < sweep_points;k++) {

	  for(ii=0;ii<lgc->len;ii++) {
	    submpp.estimate_first_freq[ii] = mpp.estimate_first_freq[ii];
	  }
	  submpp.estimate_selection_hemi = 0.0 - sel_sweep + (sweep_incr*k);
	  
	  switch( type ) {
	  case LocusModel_SingleHomozygous_First_Selection :
	    estimate_homozygous_first_selection_multi(lgc,&submpp,0);
	    break;
	  case LocusModel_SingleHomozygous_Second_Selection :
	    estimate_homozygous_second_selection_multi(lgc,&submpp);
	    break;
	  case LocusModel_Hetrozygous_Selection :
	    estimate_hetreozygous_selection_multi(lgc,&submpp);
	    break;
	  default :
	    fatal("No model for type %d yet",type);
	    break;
	  }

	  add_and_update_ModelProbabilitySet(mps,&submpp);

	} /* end of for (k) selection sweep */


      } /* end of else is a selection type */
    } else {
      /* overflow */

    }


    /* we end when all populations hit position */
    for(i=0;i<lgc->len;i++) {
      if( pop_point[i] != points ) 
	break;
    }

    if( i >= lgc->len ) {
      break;
    }

    /* if we have reached here - we need to update positions 
     * start at the end add one; if it is at position, update
     * next to the left and set to zero - recurse to the end
     */

    for(i=lgc->len-1;i >= 0;i--) {
      if( pop_point[i] != points ) {
	pop_point[i]++;
	break;
      } else {
	pop_point[i] = 0;
	/* return to for loop */
      }

    }
    
  }

  return mps;
       
}

# line 167 "locus_model_estimators.dy"
void estimate_normal_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp)
{
  int i;
  Score sc= 0;

  assert(lgc != NULL);
  assert(mpp != NULL);


  for(i=0;i < lgc->len;i++) {

    sc += 
      (2*Probability2Score(mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousFirst] +
      (2*Probability2Score(1.0-mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousSecond] +
      (Probability2Score(2*mpp->estimate_first_freq[i]*(1.0-mpp->estimate_first_freq[i])))*lgc->pgc[i]->count[BiGenotypeHetrozygous];
  }
  
  /*  fprintf(stderr,"For model with %f and %f, got %d\n",mpp->estimate_first_freq[0],mpp->estimate_first_freq[1],sc);*/


  mpp->likelihood_score = sc;

  return;

}


# line 194 "locus_model_estimators.dy"
void estimate_homozygous_first_selection_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp,int show)
{
  int i;
  Score sc= 0;
  double sel_prob;
  double rem_prob;

  Score first;
  Score het;
  Score second;

  Score first_p;
  Score second_p;
  Score het_p;

  assert(lgc != NULL);
  assert(mpp != NULL);


  sel_prob = 1.0-mpp->estimate_selection_hemi;
  rem_prob = 1.0+mpp->estimate_selection_hemi;

  for(i=0;i < lgc->len;i++) {

    first = ((Probability2Score(sel_prob)+2*Probability2Score(mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousFirst]);

    second = (Probability2Score(rem_prob)+2*Probability2Score(1.0-mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousSecond];
    
    het =       (Probability2Score(rem_prob)+Probability2Score(2*mpp->estimate_first_freq[i]*(1.0-mpp->estimate_first_freq[i])))*lgc->pgc[i]->count[BiGenotypeHetrozygous];

    first_p = ((Probability2Score(rem_prob)+2*Probability2Score(mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousFirst]);

    second_p = (Probability2Score(sel_prob)+2*Probability2Score(1.0-mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousSecond];
    
    het_p =       (Probability2Score(rem_prob)+Probability2Score(2*mpp->estimate_first_freq[i]*(1.0-mpp->estimate_first_freq[i])))*lgc->pgc[i]->count[BiGenotypeHetrozygous];

    if( show) {
      fprintf(stderr,"In pop %d\t%d %d,\t%d %d,\t%d %d\n",
	      i,
	      first,lgc->pgc[i]->count[BiGenotypeHomozygousFirst],
	      het,lgc->pgc[i]->count[BiGenotypeHetrozygous],
	      second,lgc->pgc[i]->count[BiGenotypeHomozygousSecond]);
      fprintf(stderr,"Prime  %d\t%d %d,\t%d %d,\t%d %d\n",
	      i,
	      first_p,lgc->pgc[i]->count[BiGenotypeHomozygousFirst],
	      het_p,lgc->pgc[i]->count[BiGenotypeHetrozygous],
	      second_p,lgc->pgc[i]->count[BiGenotypeHomozygousSecond]);
      
    }

    sc += first + second + het;
     
  }
  
  /*  fprintf(stderr,"For model with %f and %f, sel, %f got %d\n",mpp->estimate_first_freq[0],mpp->estimate_first_freq[1],sel_prob,sc);*/


  mpp->likelihood_score = sc;

  return;

}


# line 258 "locus_model_estimators.dy"
void estimate_hetreozygous_selection_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp)
{
  int i;
  Score sc= 0;
  double sel_prob;
  double rem_prob;

  Score first;
  Score het;
  Score second;

  assert(lgc != NULL);
  assert(mpp != NULL);


  sel_prob = 1.0-mpp->estimate_selection_hemi;
  rem_prob = 1.0+mpp->estimate_selection_hemi;

  for(i=0;i < lgc->len;i++) {

    first = ((Probability2Score(rem_prob)+2*Probability2Score(mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousFirst]);

    second = (Probability2Score(rem_prob)+2*Probability2Score(1.0-mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousSecond];
    
    het =       (Probability2Score(sel_prob)+Probability2Score(2*mpp->estimate_first_freq[i]*(1.0-mpp->estimate_first_freq[i])))*lgc->pgc[i]->count[BiGenotypeHetrozygous];

    sc += first + second + het;
     
  }
  


  mpp->likelihood_score = sc;

  return;

}


# line 297 "locus_model_estimators.dy"
void estimate_homozygous_second_selection_multi(LocusGenotypeCount * lgc,ModelProbabilityPoint * mpp)
{
  int i;
  Score sc= 0;
  double sel_prob;
  double rem_prob;

  Score first;
  Score second;
  Score het;

  assert(lgc != NULL);
  assert(mpp != NULL);


  sel_prob = 1.0-mpp->estimate_selection_hemi;
  rem_prob = 1.0+mpp->estimate_selection_hemi;

  for(i=0;i < lgc->len;i++) {

    first = ((Probability2Score(rem_prob)+2*Probability2Score(mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousFirst]);

    second = (Probability2Score(sel_prob)+2*Probability2Score(1.0-mpp->estimate_first_freq[i]))*lgc->pgc[i]->count[BiGenotypeHomozygousSecond];
    
    het =       (Probability2Score(rem_prob)+Probability2Score(2*mpp->estimate_first_freq[i]*(1.0-mpp->estimate_first_freq[i])))*lgc->pgc[i]->count[BiGenotypeHetrozygous];

    sc += first + second + het;

  }
  
  /*  fprintf(stderr,"For model with %f and %f, got %d\n",mpp->estimate_first_freq[0],mpp->estimate_first_freq[1],sc);*/


  mpp->likelihood_score = sc;

  return;

}


# line 337 "locus_model_estimators.dy"
void add_and_update_ModelProbabilitySet(ModelProbabilitySet * out,ModelProbabilityPoint * mpp)
{
  int i;

  if( out->best->likelihood_score < mpp->likelihood_score ) {
    out->best->estimate_selection_hemi = mpp->estimate_selection_hemi;
    out->best->likelihood_score = mpp->likelihood_score;
    for(i=0;i<MAX_POP;i++) {
      out->best->estimate_first_freq[i] = mpp->estimate_first_freq[i]; 
    }

  }


}



# line 355 "locus_model_estimators.dy"
ModelProbabilitySet * estimate_normal_set(LocusGenotypeCount * lgc,double sweep, double points)
{
  int i;
  int first_total = 0;
  int defined_total = 0;

  double central;
  double incr;

  ModelProbabilitySet * mps;
  ModelProbabilityPoint * mpp;

  mps = ModelProbabilitySet_alloc_std();

  for(i=0;i<lgc->len;i++) {
    first_total   += (lgc->pgc[i]->count[BiGenotypeHomozygousFirst]*2 + lgc->pgc[i]->count[BiGenotypeHetrozygous]);
    defined_total += lgc->pgc[i]->defined_total*2;
  }

  central = (double) first_total / (double) defined_total;




  incr = sweep/points;

  mps->best = NULL;
  for(i=0;i<points;i++) {
    mpp = estimate_normal(lgc,central-sweep+(incr*i));
    fprintf(stderr,"Estimate for %f is %f\n",(central-sweep+(incr*i)),Score2Bits(mpp->likelihood_score));

    add_ModelProbabilitySet(mps,mpp);
    if( mps->best == NULL || mps->best->likelihood_score < mpp->likelihood_score ) {
      mps->best = mpp;
    }
  }

  for(i=0;i<points;i++) {
    mpp = estimate_normal(lgc,central+sweep-(incr*i));
    add_ModelProbabilitySet(mps,mpp);
    if( mps->best == NULL ||  mps->best->likelihood_score < mpp->likelihood_score ) {
      mps->best = mpp;
    }
  }

  return mps;

}

# line 404 "locus_model_estimators.dy"
ModelProbabilityPoint * estimate_normal(LocusGenotypeCount * lgc,double estimate_first_freq)
{
  int i;
  Score sub_likelihood;

  ModelProbabilityPoint * mpp;

  mpp = ModelProbabilityPoint_alloc();

  mpp->estimate_selection_hemi = 0.0;
  mpp->estimate_first_freq[0] = estimate_first_freq;

  /* likelihood is
   * Prob(Data|Model)
   */
     
  sub_likelihood = 0;
  for(i=0;i<lgc->len;i++) {

    sub_likelihood += 
      (2*Probability2Score(estimate_first_freq))*lgc->pgc[i]->count[BiGenotypeHomozygousFirst] +
      
      (2*Probability2Score(1.0-estimate_first_freq))*lgc->pgc[i]->count[BiGenotypeHomozygousSecond] +
      (Probability2Score(estimate_first_freq*(1.0-estimate_first_freq))*lgc->pgc[i]->count[BiGenotypeHetrozygous]);
  }

  mpp->likelihood_score = sub_likelihood;

  return mpp;
}


# line 407 "locus_model_estimators.c"
/* Function:  hard_link_ModelProbabilityPoint(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModelProbabilityPoint *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilityPoint *]
 *
 */
ModelProbabilityPoint * hard_link_ModelProbabilityPoint(ModelProbabilityPoint * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ModelProbabilityPoint object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ModelProbabilityPoint_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilityPoint *]
 *
 */
ModelProbabilityPoint * ModelProbabilityPoint_alloc(void) 
{
    ModelProbabilityPoint * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ModelProbabilityPoint *) ckalloc (sizeof(ModelProbabilityPoint))) == NULL)  {  
      warn("ModelProbabilityPoint_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->estimate_selection_hemi = 0;    
    /* estimate_first_freq[MAX_POP] is an array: no default possible */ 
    out->likelihood_score = 0;   


    return out;  
}    


/* Function:  free_ModelProbabilityPoint(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModelProbabilityPoint *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilityPoint *]
 *
 */
ModelProbabilityPoint * free_ModelProbabilityPoint(ModelProbabilityPoint * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ModelProbabilityPoint obj. Should be trappable"); 
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


/* Function:  swap_ModelProbabilitySet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ModelProbabilitySet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [ModelProbabilityPoint **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ModelProbabilitySet(ModelProbabilityPoint ** list,int i,int j)  
{
    ModelProbabilityPoint * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ModelProbabilitySet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ModelProbabilitySet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [ModelProbabilityPoint **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ModelProbabilitySet(ModelProbabilityPoint ** list,int left,int right,int (*comp)(ModelProbabilityPoint * ,ModelProbabilityPoint * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ModelProbabilitySet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ModelProbabilitySet (list,++last,i);    
      }  
    swap_ModelProbabilitySet (list,left,last);   
    qsort_ModelProbabilitySet(list,left,last-1,comp);    
    qsort_ModelProbabilitySet(list,last+1,right,comp);   
}    


/* Function:  sort_ModelProbabilitySet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ModelProbabilitySet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ModelProbabilitySet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ModelProbabilitySet(ModelProbabilitySet * obj,int (*comp)(ModelProbabilityPoint *, ModelProbabilityPoint *)) 
{
    qsort_ModelProbabilitySet(obj->mpp,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_ModelProbabilitySet(obj,len)
 *
 * Descrip:    Really an internal function for add_ModelProbabilitySet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ModelProbabilitySet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ModelProbabilitySet(ModelProbabilitySet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ModelProbabilitySet called with no need");    
      return TRUE;   
      }  


    if( (obj->mpp = (ModelProbabilityPoint ** ) ckrealloc (obj->mpp,sizeof(ModelProbabilityPoint *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_ModelProbabilitySet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ModelProbabilitySet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ModelProbabilitySet *]
 * Arg:        add [OWNER] Object to add to the list [ModelProbabilityPoint *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ModelProbabilitySet(ModelProbabilitySet * obj,ModelProbabilityPoint * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ModelProbabilitySet(obj,obj->len + ModelProbabilitySetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->mpp[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ModelProbabilitySet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ModelProbabilitySet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ModelProbabilitySet(ModelProbabilitySet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->mpp[i] != NULL)   {  
        free_ModelProbabilityPoint(obj->mpp[i]); 
        obj->mpp[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ModelProbabilitySet_alloc_std(void)
 *
 * Descrip:    Equivalent to ModelProbabilitySet_alloc_len(ModelProbabilitySetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * ModelProbabilitySet_alloc_std(void) 
{
    return ModelProbabilitySet_alloc_len(ModelProbabilitySetLISTLENGTH); 
}    


/* Function:  ModelProbabilitySet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * ModelProbabilitySet_alloc_len(int len) 
{
    ModelProbabilitySet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ModelProbabilitySet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->mpp = (ModelProbabilityPoint ** ) ckcalloc (len,sizeof(ModelProbabilityPoint *))) == NULL)  {  
      warn("Warning, ckcalloc failed in ModelProbabilitySet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ModelProbabilitySet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ModelProbabilitySet *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * hard_link_ModelProbabilitySet(ModelProbabilitySet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ModelProbabilitySet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ModelProbabilitySet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * ModelProbabilitySet_alloc(void) 
{
    ModelProbabilitySet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ModelProbabilitySet *) ckalloc (sizeof(ModelProbabilitySet))) == NULL)  {  
      warn("ModelProbabilitySet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->keep = 0;   
    out->mpp = NULL; 
    out->len = out->maxlen = 0;  
    out->best = NULL;    
    out->total_likelihood = 0;   


    return out;  
}    


/* Function:  free_ModelProbabilitySet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ModelProbabilitySet *]
 *
 * Return [UNKN ]  Undocumented return value [ModelProbabilitySet *]
 *
 */
ModelProbabilitySet * free_ModelProbabilitySet(ModelProbabilitySet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ModelProbabilitySet obj. Should be trappable");   
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
    /* obj->data is linked in */ 
    if( obj->mpp != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->mpp[i] != NULL) 
          free_ModelProbabilityPoint(obj->mpp[i]);   
        }  
      ckfree(obj->mpp);  
      }  
    if( obj->best != NULL)   
      free_ModelProbabilityPoint(obj->best);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
