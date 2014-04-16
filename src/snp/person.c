#ifdef _cplusplus
extern "C" {
#endif
#include "person.h"

# line 28 "person.dy"
boolean is_in_Population(Population * p,char * person_id)
{
  int i;

  for(i=0;i<p->len;i++) {
    if( strcmp(p->person[i]->person_id,person_id) == 0 ) {
      return TRUE;
    }
  }

  return FALSE;
}

# line 41 "person.dy"
Population * find_or_new_Population_in_PopulationSet(PopulationSet * ps,char * pop_name)
{
  int i;
  Population * pe = NULL;
  
  assert(ps != NULL);

  for(i=0;i<ps->len;i++) {
    if( strcmp(ps->pop[i]->pop_name,pop_name) == 0 ) {
      pe = ps->pop[i];
      break;
    }
  }

  if( pe != NULL ) {
    return pe;
  }

  pe = Population_alloc_std();
  pe->pop_name = stringalloc(pop_name);

  add_PopulationSet(ps,pe);
  
  return pe;

}

# line 68 "person.dy"
Person * find_or_new_Person_in_Population(Population * p,char * person_id)
{
  int i;
  Person * pe = NULL;
  
  assert(p != NULL);

  for(i=0;i<p->len;i++) {
    if( strcmp(p->person[i]->person_id,person_id) == 0 ) {
      pe = p->person[i];
      break;
    }
  }

  if( pe != NULL ) {
    return pe;
  }

  pe = Person_alloc();
  pe->person_id = stringalloc(person_id);
  pe->pop = p;

  add_Population(p,pe);
  
  return pe;
}

# line 75 "person.c"
/* Function:  hard_link_Person(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Person *]
 *
 * Return [UNKN ]  Undocumented return value [Person *]
 *
 */
Person * hard_link_Person(Person * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Person object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Person_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Person *]
 *
 */
Person * Person_alloc(void) 
{
    Person * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Person *) ckalloc (sizeof(Person))) == NULL)    {  
      warn("Person_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->person_id = NULL;   


    return out;  
}    


/* Function:  free_Person(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Person *]
 *
 * Return [UNKN ]  Undocumented return value [Person *]
 *
 */
Person * free_Person(Person * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Person obj. Should be trappable");    
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
    if( obj->person_id != NULL)  
      ckfree(obj->person_id);    
    /* obj->pop is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Population(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Population
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Person **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Population(Person ** list,int i,int j)  
{
    Person * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Population(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Population which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Person **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Population(Person ** list,int left,int right,int (*comp)(Person * ,Person * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Population(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Population (list,++last,i); 
      }  
    swap_Population (list,left,last);    
    qsort_Population(list,left,last-1,comp); 
    qsort_Population(list,last+1,right,comp);    
}    


/* Function:  sort_Population(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Population
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Population *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Population(Population * obj,int (*comp)(Person *, Person *)) 
{
    qsort_Population(obj->person,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_Population(obj,len)
 *
 * Descrip:    Really an internal function for add_Population
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Population *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Population(Population * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Population called with no need"); 
      return TRUE;   
      }  


    if( (obj->person = (Person ** ) ckrealloc (obj->person,sizeof(Person *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Population, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Population(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Population *]
 * Arg:        add [OWNER] Object to add to the list [Person *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Population(Population * obj,Person * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Population(obj,obj->len + PopulationLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->person[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_Population(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Population *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Population(Population * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->person[i] != NULL)    {  
        free_Person(obj->person[i]); 
        obj->person[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Population_alloc_std(void)
 *
 * Descrip:    Equivalent to Population_alloc_len(PopulationLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Population_alloc_std(void) 
{
    return Population_alloc_len(PopulationLISTLENGTH);   
}    


/* Function:  Population_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Population_alloc_len(int len) 
{
    Population * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Population_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->person = (Person ** ) ckcalloc (len,sizeof(Person *))) == NULL) {  
      warn("Warning, ckcalloc failed in Population_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Population(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Population *]
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * hard_link_Population(Population * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Population object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Population_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * Population_alloc(void) 
{
    Population * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Population *) ckalloc (sizeof(Population))) == NULL)    {  
      warn("Population_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pop_name = NULL;    
    out->person = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_Population(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Population *]
 *
 * Return [UNKN ]  Undocumented return value [Population *]
 *
 */
Population * free_Population(Population * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Population obj. Should be trappable");    
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
    if( obj->pop_name != NULL)   
      ckfree(obj->pop_name);     
    if( obj->person != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->person[i] != NULL)  
          free_Person(obj->person[i]);   
        }  
      ckfree(obj->person);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_PopulationSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PopulationSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Population **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PopulationSet(Population ** list,int i,int j)  
{
    Population * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PopulationSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PopulationSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Population **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PopulationSet(Population ** list,int left,int right,int (*comp)(Population * ,Population * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PopulationSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PopulationSet (list,++last,i);  
      }  
    swap_PopulationSet (list,left,last); 
    qsort_PopulationSet(list,left,last-1,comp);  
    qsort_PopulationSet(list,last+1,right,comp); 
}    


/* Function:  sort_PopulationSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PopulationSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PopulationSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PopulationSet(PopulationSet * obj,int (*comp)(Population *, Population *)) 
{
    qsort_PopulationSet(obj->pop,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_PopulationSet(obj,len)
 *
 * Descrip:    Really an internal function for add_PopulationSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PopulationSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PopulationSet(PopulationSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PopulationSet called with no need");  
      return TRUE;   
      }  


    if( (obj->pop = (Population ** ) ckrealloc (obj->pop,sizeof(Population *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_PopulationSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PopulationSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PopulationSet *]
 * Arg:        add [OWNER] Object to add to the list [Population *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PopulationSet(PopulationSet * obj,Population * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PopulationSet(obj,obj->len + PopulationSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->pop[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_PopulationSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PopulationSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PopulationSet(PopulationSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pop[i] != NULL)   {  
        free_Population(obj->pop[i]);    
        obj->pop[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PopulationSet_alloc_std(void)
 *
 * Descrip:    Equivalent to PopulationSet_alloc_len(PopulationSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * PopulationSet_alloc_std(void) 
{
    return PopulationSet_alloc_len(PopulationSetLISTLENGTH); 
}    


/* Function:  PopulationSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * PopulationSet_alloc_len(int len) 
{
    PopulationSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PopulationSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pop = (Population ** ) ckcalloc (len,sizeof(Population *))) == NULL)    {  
      warn("Warning, ckcalloc failed in PopulationSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PopulationSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PopulationSet *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * hard_link_PopulationSet(PopulationSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PopulationSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PopulationSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * PopulationSet_alloc(void) 
{
    PopulationSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PopulationSet *) ckalloc (sizeof(PopulationSet))) == NULL)  {  
      warn("PopulationSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pop = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_PopulationSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PopulationSet *]
 *
 * Return [UNKN ]  Undocumented return value [PopulationSet *]
 *
 */
PopulationSet * free_PopulationSet(PopulationSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PopulationSet obj. Should be trappable"); 
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
    if( obj->pop != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pop[i] != NULL) 
          free_Population(obj->pop[i]);  
        }  
      ckfree(obj->pop);  
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
