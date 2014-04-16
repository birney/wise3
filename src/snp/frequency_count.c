#ifdef _cplusplus
extern "C" {
#endif
#include "frequency_count.h"


# line 32 "frequency_count.dy"
double chisquared_stat_LocusGenotypeCount(LocusGenotypeCount * lgc)
{
  double ret = 0.0;
  int i;

  for(i=0;i<lgc->len;i++) {
    ret += chisquared_stat_PopulationGenotypeCount(lgc->pgc[i]);
  }

  return ret;
}


# line 45 "frequency_count.dy"
double chisquared_stat_PopulationGenotypeCount(PopulationGenotypeCount * pgc)
{
  double first;
  double second;
  double ret = 0.0;
  double diff;
  double expec;

  assert(pgc->defined_total > 5);

  first = central_first_frequency_PopulationGenotypeCount(pgc);

  if( first < 0.0000000000000001 || first > 0.999999999999999 ) {
    warn("Problem in chisquared statistic with central first freq at %f, returning 0.0",first);
    return 0.0;
  }


  second = 1.0 - first;


  expec =  (first*first * (double)pgc->defined_total);

  /*  fprintf(stderr,"first is %.2f, totalis %d, expect is %f\n",first,pgc->defined_total,expec);*/

  assert(expec > 0.0);

  diff = pgc->count[BiGenotypeHomozygousFirst] - expec;
  
  ret += (diff*diff)/expec;

  expec =  (first*second*2 * (double)pgc->defined_total);
  assert(expec > 0.0);

  diff = pgc->count[BiGenotypeHetrozygous] - expec;
  
  ret += (diff*diff)/expec;


  expec =  (second*second * (double)pgc->defined_total);
  assert(expec > 0.0);

  diff = pgc->count[BiGenotypeHetrozygous] - expec;
  
  ret += (diff*diff)/expec;

  return ret;

}

# line 95 "frequency_count.dy"
boolean seen_each_genotype_in_all_populations_LocusGenotypeCount(LocusGenotypeCount * lgc)
{
  int i;


  for(i=0;i<lgc->len;i++) {
    if( lgc->pgc[i]->count[BiGenotypeHomozygousFirst] == 0 ||
	lgc->pgc[i]->count[BiGenotypeHomozygousSecond] == 0 ||
	lgc->pgc[i]->count[BiGenotypeHetrozygous] == 0 ) {
      return FALSE;
    }
  }

  return TRUE;
}


# line 112 "frequency_count.dy"
double central_first_frequency_PopulationGenotypeCount(PopulationGenotypeCount * pgc)
{
  assert(pgc != NULL);
  assert(pgc->defined_total > 0);

  return ((double)(pgc->count[BiGenotypeHomozygousFirst]*2+pgc->count[BiGenotypeHetrozygous]))/
    (2*(double)(pgc->defined_total));
}

# line 121 "frequency_count.dy"
double smallest_minor_allele_LocusGenotypeCount(LocusGenotypeCount * lgc) 
{
  int i;
  double smallest = 0.5;
  double central;


  for(i=0;i<lgc->len;i++) {
    central = central_first_frequency_PopulationGenotypeCount(lgc->pgc[i]);
    if( central > 0.5 ) {
      if( smallest > (1.0 - central)) {
	smallest = (1.0 - central);
      }
    } else {
      if( smallest > central ) {
	smallest = central;
      }
    }
  }

  return smallest;
}

# line 144 "frequency_count.dy"
LocusGenotypeCount * resampled_LocusGenotypeCount(LocusGenotypeCount * lgc)
{
  int i,j;
  LocusGenotypeCount * out;
  PopulationGenotypeCount * pgc;
  double rnd1;
  double rnd2;
  double central;

  assert(lgc != NULL);

  out = LocusGenotypeCount_alloc_len(lgc->len);
  out->bl = lgc->bl;

  for(i=0;i<lgc->len;i++) {


    pgc  = PopulationGenotypeCount_alloc();

    pgc->count[BiGenotypeHomozygousFirst] = 0;
    pgc->count[BiGenotypeHetrozygous] = 0;
    pgc->count[BiGenotypeHomozygousSecond] = 0;

    central = central_first_frequency_PopulationGenotypeCount(lgc->pgc[i]);
    for(j=0;j<lgc->pgc[i]->defined_total;j++) {
      rnd1 = random_0_to_1();
      rnd2 = random_0_to_1();
      if( rnd1 < central && rnd2 < central ) {
	pgc->count[BiGenotypeHomozygousFirst]++;
      } else if ( rnd1 > central && rnd2 > central ) {
	pgc->count[BiGenotypeHomozygousSecond]++;
      } else {
	pgc->count[BiGenotypeHetrozygous]++;
      }
    }
    pgc->defined_total = lgc->pgc[i]->defined_total;
    add_LocusGenotypeCount(out,pgc);
  }

  return out;
}

# line 186 "frequency_count.dy"
LocusGenotypeCountSet * LocusGenotypeCountSet_from_BiGenotypeSet(BiGenotypeSet * bgs)
{
  int i;
  int j;
  int p;
  LocusGenotypeCount * lgc;
  PopulationGenotypeCount* pgc;
  LocusGenotypeCountSet * out;

  out = LocusGenotypeCountSet_alloc_len(bgs->len);

  for(i=0;i<bgs->len;i++) {
    lgc = LocusGenotypeCount_alloc_len(bgs->ps->len);
    lgc->bl = bgs->locus[i]->locus;

    
    for(p=0;p<bgs->ps->len;p++) {

      pgc = PopulationGenotypeCount_alloc();
      pgc->pop = bgs->ps->pop[p];
      add_LocusGenotypeCount(lgc,pgc);

      pgc->count[BiGenotypeHomozygousFirst] = 0;
      pgc->count[BiGenotypeHomozygousSecond] = 0;
      pgc->count[BiGenotypeHetrozygous] = 0;
      pgc->count[BiGenotypeUnknown] = 0;

      pgc->defined_total = 0;

      for(j=0;j<bgs->locus[i]->len;j++) {

	if( bgs->locus[i]->big[j]->person->pop != pgc->pop ) {
	  continue;
	}

	pgc->count[bgs->locus[i]->big[j]->type]++;
	if( bgs->locus[i]->big[j]->type != BiGenotypeUnknown ) {
	  pgc->defined_total++;
	}
	
      }


    }
    add_LocusGenotypeCountSet(out,lgc);

  }

  return out;
}



# line 220 "frequency_count.c"
/* Function:  hard_link_PopulationGenotypeCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PopulationGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationGenotypeCount *]
 *
 */
PopulationGenotypeCount * hard_link_PopulationGenotypeCount(PopulationGenotypeCount * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PopulationGenotypeCount object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PopulationGenotypeCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PopulationGenotypeCount *]
 *
 */
PopulationGenotypeCount * PopulationGenotypeCount_alloc(void) 
{
    PopulationGenotypeCount * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PopulationGenotypeCount *) ckalloc (sizeof(PopulationGenotypeCount))) == NULL)  {  
      warn("PopulationGenotypeCount_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* count[BiGenotypeLength] is an array: no default possible */ 
    out->defined_total = 0;  


    return out;  
}    


/* Function:  free_PopulationGenotypeCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PopulationGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationGenotypeCount *]
 *
 */
PopulationGenotypeCount * free_PopulationGenotypeCount(PopulationGenotypeCount * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PopulationGenotypeCount obj. Should be trappable");   
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
    /* obj->pop is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LocusGenotypeCount(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LocusGenotypeCount
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PopulationGenotypeCount **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LocusGenotypeCount(PopulationGenotypeCount ** list,int i,int j)  
{
    PopulationGenotypeCount * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LocusGenotypeCount(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LocusGenotypeCount which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PopulationGenotypeCount **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LocusGenotypeCount(PopulationGenotypeCount ** list,int left,int right,int (*comp)(PopulationGenotypeCount * ,PopulationGenotypeCount * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LocusGenotypeCount(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LocusGenotypeCount (list,++last,i); 
      }  
    swap_LocusGenotypeCount (list,left,last);    
    qsort_LocusGenotypeCount(list,left,last-1,comp); 
    qsort_LocusGenotypeCount(list,last+1,right,comp);    
}    


/* Function:  sort_LocusGenotypeCount(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LocusGenotypeCount
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LocusGenotypeCount *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LocusGenotypeCount(LocusGenotypeCount * obj,int (*comp)(PopulationGenotypeCount *, PopulationGenotypeCount *)) 
{
    qsort_LocusGenotypeCount(obj->pgc,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_LocusGenotypeCount(obj,len)
 *
 * Descrip:    Really an internal function for add_LocusGenotypeCount
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocusGenotypeCount *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LocusGenotypeCount(LocusGenotypeCount * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LocusGenotypeCount called with no need"); 
      return TRUE;   
      }  


    if( (obj->pgc = (PopulationGenotypeCount ** ) ckrealloc (obj->pgc,sizeof(PopulationGenotypeCount *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_LocusGenotypeCount, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LocusGenotypeCount(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocusGenotypeCount *]
 * Arg:        add [OWNER] Object to add to the list [PopulationGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LocusGenotypeCount(LocusGenotypeCount * obj,PopulationGenotypeCount * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LocusGenotypeCount(obj,obj->len + LocusGenotypeCountLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->pgc[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_LocusGenotypeCount(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LocusGenotypeCount(LocusGenotypeCount * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pgc[i] != NULL)   {  
        free_PopulationGenotypeCount(obj->pgc[i]);   
        obj->pgc[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LocusGenotypeCount_alloc_std(void)
 *
 * Descrip:    Equivalent to LocusGenotypeCount_alloc_len(LocusGenotypeCountLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * LocusGenotypeCount_alloc_std(void) 
{
    return LocusGenotypeCount_alloc_len(LocusGenotypeCountLISTLENGTH);   
}    


/* Function:  LocusGenotypeCount_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * LocusGenotypeCount_alloc_len(int len) 
{
    LocusGenotypeCount * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LocusGenotypeCount_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pgc = (PopulationGenotypeCount ** ) ckcalloc (len,sizeof(PopulationGenotypeCount *))) == NULL)  {  
      warn("Warning, ckcalloc failed in LocusGenotypeCount_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LocusGenotypeCount(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * hard_link_LocusGenotypeCount(LocusGenotypeCount * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocusGenotypeCount object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocusGenotypeCount_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * LocusGenotypeCount_alloc(void) 
{
    LocusGenotypeCount * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocusGenotypeCount *) ckalloc (sizeof(LocusGenotypeCount))) == NULL)    {  
      warn("LocusGenotypeCount_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pgc = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_LocusGenotypeCount(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCount *]
 *
 */
LocusGenotypeCount * free_LocusGenotypeCount(LocusGenotypeCount * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocusGenotypeCount obj. Should be trappable");    
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
    /* obj->bl is linked in */ 
    if( obj->pgc != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pgc[i] != NULL) 
          free_PopulationGenotypeCount(obj->pgc[i]); 
        }  
      ckfree(obj->pgc);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_LocusGenotypeCountSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_LocusGenotypeCountSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [LocusGenotypeCount **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_LocusGenotypeCountSet(LocusGenotypeCount ** list,int i,int j)  
{
    LocusGenotypeCount * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_LocusGenotypeCountSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_LocusGenotypeCountSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [LocusGenotypeCount **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_LocusGenotypeCountSet(LocusGenotypeCount ** list,int left,int right,int (*comp)(LocusGenotypeCount * ,LocusGenotypeCount * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_LocusGenotypeCountSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_LocusGenotypeCountSet (list,++last,i);  
      }  
    swap_LocusGenotypeCountSet (list,left,last); 
    qsort_LocusGenotypeCountSet(list,left,last-1,comp);  
    qsort_LocusGenotypeCountSet(list,last+1,right,comp); 
}    


/* Function:  sort_LocusGenotypeCountSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_LocusGenotypeCountSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [LocusGenotypeCountSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_LocusGenotypeCountSet(LocusGenotypeCountSet * obj,int (*comp)(LocusGenotypeCount *, LocusGenotypeCount *)) 
{
    qsort_LocusGenotypeCountSet(obj->lgc,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_LocusGenotypeCountSet(obj,len)
 *
 * Descrip:    Really an internal function for add_LocusGenotypeCountSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocusGenotypeCountSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_LocusGenotypeCountSet(LocusGenotypeCountSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_LocusGenotypeCountSet called with no need");  
      return TRUE;   
      }  


    if( (obj->lgc = (LocusGenotypeCount ** ) ckrealloc (obj->lgc,sizeof(LocusGenotypeCount *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_LocusGenotypeCountSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_LocusGenotypeCountSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [LocusGenotypeCountSet *]
 * Arg:        add [OWNER] Object to add to the list [LocusGenotypeCount *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_LocusGenotypeCountSet(LocusGenotypeCountSet * obj,LocusGenotypeCount * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_LocusGenotypeCountSet(obj,obj->len + LocusGenotypeCountSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->lgc[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_LocusGenotypeCountSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [LocusGenotypeCountSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_LocusGenotypeCountSet(LocusGenotypeCountSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->lgc[i] != NULL)   {  
        free_LocusGenotypeCount(obj->lgc[i]);    
        obj->lgc[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  LocusGenotypeCountSet_alloc_std(void)
 *
 * Descrip:    Equivalent to LocusGenotypeCountSet_alloc_len(LocusGenotypeCountSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * LocusGenotypeCountSet_alloc_std(void) 
{
    return LocusGenotypeCountSet_alloc_len(LocusGenotypeCountSetLISTLENGTH); 
}    


/* Function:  LocusGenotypeCountSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * LocusGenotypeCountSet_alloc_len(int len) 
{
    LocusGenotypeCountSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = LocusGenotypeCountSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->lgc = (LocusGenotypeCount ** ) ckcalloc (len,sizeof(LocusGenotypeCount *))) == NULL)    {  
      warn("Warning, ckcalloc failed in LocusGenotypeCountSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_LocusGenotypeCountSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [LocusGenotypeCountSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * hard_link_LocusGenotypeCountSet(LocusGenotypeCountSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a LocusGenotypeCountSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  LocusGenotypeCountSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * LocusGenotypeCountSet_alloc(void) 
{
    LocusGenotypeCountSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocusGenotypeCountSet *) ckalloc (sizeof(LocusGenotypeCountSet))) == NULL)  {  
      warn("LocusGenotypeCountSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->lgc = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_LocusGenotypeCountSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocusGenotypeCountSet *]
 *
 * Return [UNKN ]  Undocumented return value [LocusGenotypeCountSet *]
 *
 */
LocusGenotypeCountSet * free_LocusGenotypeCountSet(LocusGenotypeCountSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocusGenotypeCountSet obj. Should be trappable"); 
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
    if( obj->lgc != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->lgc[i] != NULL) 
          free_LocusGenotypeCount(obj->lgc[i]);  
        }  
      ckfree(obj->lgc);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
