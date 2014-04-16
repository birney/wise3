#ifdef _cplusplus
extern "C" {
#endif
#include "locus_framework.h"

# line 36 "locus_framework.dy"
BiGenotypeType type_from_hapmap_string_BiLocus(BiLocus * bi,char * g)
{
  if( g[0] == bi->first_allele_char && g[1] == bi->first_allele_char ) {
    return BiGenotypeHomozygousFirst;
  }

  if( g[0] == bi->first_allele_char && g[1] == bi->second_allele_char ) {
    return BiGenotypeHetrozygous;
  }

  if( g[1] == bi->first_allele_char && g[0] == bi->second_allele_char ) {
    return BiGenotypeHetrozygous;
  }

  if( g[0] == bi->second_allele_char && g[1] == bi->second_allele_char ) {
    return BiGenotypeHomozygousSecond;
  }

  if( g[0] == 'N' && g[1] == 'N' ) {
    return BiGenotypeUnknown;
  }

  return BiGenotypeError;

}

# line 62 "locus_framework.dy"
BiLocusFramework * new_BiLocusFramework(void)
{
  BiLocusFramework * out;

  out = BiLocusFramework_alloc_std();

  out->hash = g_hash_table_new(g_str_hash,g_str_equal);

  return out;
}


# line 74 "locus_framework.dy"
BiLocus * find_or_new_BiLocus_from_BiLocusFramework(BiLocusFramework * bgf,char * ref_id,char first,char second)
{
  BiLocus * out;
  int i;


  if( (out = (BiLocus*) g_hash_table_lookup(bgf->hash,(gconstpointer)ref_id)) != NULL ) {
    return out;
  }

  out = BiLocus_alloc();
  out->locus_id = stringalloc(ref_id);
  out->first_allele_char = first;
  out->second_allele_char = second;

  add_BiLocusFramework(bgf,out);
  
  g_hash_table_insert(bgf->hash,(gpointer)out->locus_id,(void*)out);

  return out;
}



# line 70 "locus_framework.c"
/* Function:  hard_link_BiLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocus *]
 *
 */
BiLocus * hard_link_BiLocus(BiLocus * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BiLocus object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BiLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiLocus *]
 *
 */
BiLocus * BiLocus_alloc(void) 
{
    BiLocus * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BiLocus *) ckalloc (sizeof(BiLocus))) == NULL)  {  
      warn("BiLocus_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->locus_id = NULL;    
    out->first_allele_char = 'u';    
    out->second_allele_char = 'u';   


    return out;  
}    


/* Function:  free_BiLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocus *]
 *
 */
BiLocus * free_BiLocus(BiLocus * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a BiLocus obj. Should be trappable");   
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
    if( obj->locus_id != NULL)   
      ckfree(obj->locus_id);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_BiLocusFramework(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_BiLocusFramework
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [BiLocus **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_BiLocusFramework(BiLocus ** list,int i,int j)  
{
    BiLocus * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_BiLocusFramework(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_BiLocusFramework which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [BiLocus **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_BiLocusFramework(BiLocus ** list,int left,int right,int (*comp)(BiLocus * ,BiLocus * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_BiLocusFramework(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_BiLocusFramework (list,++last,i);   
      }  
    swap_BiLocusFramework (list,left,last);  
    qsort_BiLocusFramework(list,left,last-1,comp);   
    qsort_BiLocusFramework(list,last+1,right,comp);  
}    


/* Function:  sort_BiLocusFramework(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_BiLocusFramework
 *
 *
 * Arg:         obj [UNKN ] Object containing list [BiLocusFramework *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_BiLocusFramework(BiLocusFramework * obj,int (*comp)(BiLocus *, BiLocus *)) 
{
    qsort_BiLocusFramework(obj->locus,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_BiLocusFramework(obj,len)
 *
 * Descrip:    Really an internal function for add_BiLocusFramework
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiLocusFramework *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_BiLocusFramework(BiLocusFramework * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_BiLocusFramework called with no need");   
      return TRUE;   
      }  


    if( (obj->locus = (BiLocus ** ) ckrealloc (obj->locus,sizeof(BiLocus *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_BiLocusFramework, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_BiLocusFramework(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiLocusFramework *]
 * Arg:        add [OWNER] Object to add to the list [BiLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_BiLocusFramework(BiLocusFramework * obj,BiLocus * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_BiLocusFramework(obj,obj->len + BiLocusFrameworkLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->locus[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_BiLocusFramework(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [BiLocusFramework *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_BiLocusFramework(BiLocusFramework * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->locus[i] != NULL) {  
        free_BiLocus(obj->locus[i]); 
        obj->locus[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  BiLocusFramework_alloc_std(void)
 *
 * Descrip:    Equivalent to BiLocusFramework_alloc_len(BiLocusFrameworkLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * BiLocusFramework_alloc_std(void) 
{
    return BiLocusFramework_alloc_len(BiLocusFrameworkLISTLENGTH);   
}    


/* Function:  BiLocusFramework_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * BiLocusFramework_alloc_len(int len) 
{
    BiLocusFramework * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = BiLocusFramework_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->locus = (BiLocus ** ) ckcalloc (len,sizeof(BiLocus *))) == NULL)    {  
      warn("Warning, ckcalloc failed in BiLocusFramework_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_BiLocusFramework(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiLocusFramework *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * hard_link_BiLocusFramework(BiLocusFramework * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BiLocusFramework object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BiLocusFramework_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * BiLocusFramework_alloc(void) 
{
    BiLocusFramework * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BiLocusFramework *) ckalloc (sizeof(BiLocusFramework))) == NULL)    {  
      warn("BiLocusFramework_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->locus = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_BiLocusFramework(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiLocusFramework *]
 *
 * Return [UNKN ]  Undocumented return value [BiLocusFramework *]
 *
 */
BiLocusFramework * free_BiLocusFramework(BiLocusFramework * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a BiLocusFramework obj. Should be trappable");  
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
    if( obj->locus != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->locus[i] != NULL)   
          free_BiLocus(obj->locus[i]);   
        }  
      ckfree(obj->locus);    
      }  
    /* obj->hash is linked in */ 


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
