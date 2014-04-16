#ifdef _cplusplus
extern "C" {
#endif
#include "graphseq.h"


# line 41 "graphseq.dy"
boolean prepare_GraphSeq(GraphSeq * gs)
{
  int i;
  GraphSeqEdgeHolder * th;

  /* check every edge has a sequence */

  for(i=0;i<gs->len;i++) {
    if( gs->comp[i]->seq == NULL ) {
      return FALSE;
    }

    /* if this has a free 5 or 3 end, add it. Add it twice if it has both */

    if( gs->comp[i]->five_len == 0 ) {
      th = GraphSeqEdgeHolder_alloc();
      th->edge = gs->comp[i];
      add_free_GraphSeq(gs,th);
    }

    if( gs->comp[i]->three_len == 0 ) {
      th = GraphSeqEdgeHolder_alloc();
      th->edge = gs->comp[i];
      add_free_GraphSeq(gs,th);
    }

  }

  
  return TRUE;

}

# line 74 "graphseq.dy"
GraphSeq * read_simple_GraphSeq(FILE * ifp)
{
  GraphSeq * out;
  GraphSeqEdge * e;
  char buffer[MAXLINE];
  char tempid[MAXLINE];
  char tempid2[MAXLINE];
  int first;
  int second;
  int templen;

  GraphSeqEdge * one;
  GraphSeqEdge * two;

  GraphSeqEdgeHolder * th;

  Sequence * seq;



  assert(ifp != NULL);

  out = GraphSeq_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '/' && buffer[1] == '/' ) {
      break;
    }

    if( strstartcmp(buffer,"edge") == 0 ) {
      if( sscanf(buffer,"edge %s %d",tempid,&templen) != 2 ) {
	warn("Bad edge line %s",buffer);
      } else {
        /* not testing that the ID is unique */
	e = GraphSeqEdge_alloc_std();
	e->id = stringalloc(tempid);
	e->len = templen;
	
	add_GraphSeq(out,e);
      }
    }

    if( strstartcmp(buffer,"join") == 0 ) {
      if( sscanf(buffer,"join %s %d %s %d",tempid,&first,tempid2,&second) != 4 ) {
	warn("Bad join line %s",buffer);
      } else {
	/* check that both edges exist */
	one = find_GraphSeqEdge(out,tempid);
	two = find_GraphSeqEdge(out,tempid2);

	if( one == NULL || two == NULL ) {
	  warn("Join line %s has unknown edges",buffer);
	  continue; /* back to while */
	}

	if( (first != 5 && first != 3) || (second != 5 && second != 3) ) {
	  warn("Unable to read join line %s, must be 5 or 3",buffer);
	  continue;
	}

	th = GraphSeqEdgeHolder_alloc();
	th->edge = two;
	if( first == 5 ) {
	  add_five_GraphSeqEdge(one,th);
	} else {
	  add_three_GraphSeqEdge(one,th);
	}

	th = GraphSeqEdgeHolder_alloc();
	th->edge = one;
	if( second == 5 ) {
	  add_five_GraphSeqEdge(two,th);
	} else {
	  add_three_GraphSeqEdge(two,th);
	}

      }
    } /* end of join */
  }
  
  while( (seq = read_fasta_Sequence(ifp)) != NULL ) {
    one = find_GraphSeqEdge(out,seq->name);
    if( one == NULL ) {
      warn("No sequence of id %s, not assigning",seq->name);
    } else {
      one->seq = seq;
    }
  }



  return(out);
}

# line 168 "graphseq.dy"
GraphSeqEdge * find_GraphSeqEdge(GraphSeq * gs,char * id) 
{
  int i;


  assert(gs != NULL);

  for(i=0;i<gs->len;i++) {
    if( strcmp(gs->comp[i]->id,id) == 0 ) {
      return(gs->comp[i]);
    }
  }

  return(NULL);
}





# line 156 "graphseq.c"
/* Function:  hard_link_GraphSeqEdgeHolder(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdgeHolder *]
 *
 */
GraphSeqEdgeHolder * hard_link_GraphSeqEdgeHolder(GraphSeqEdgeHolder * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GraphSeqEdgeHolder object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GraphSeqEdgeHolder_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdgeHolder *]
 *
 */
GraphSeqEdgeHolder * GraphSeqEdgeHolder_alloc(void) 
{
    GraphSeqEdgeHolder * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GraphSeqEdgeHolder *) ckalloc (sizeof(GraphSeqEdgeHolder))) == NULL)    {  
      warn("GraphSeqEdgeHolder_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   


    return out;  
}    


/* Function:  free_GraphSeqEdgeHolder(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdgeHolder *]
 *
 */
GraphSeqEdgeHolder * free_GraphSeqEdgeHolder(GraphSeqEdgeHolder * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GraphSeqEdgeHolder obj. Should be trappable");    
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
    /* obj->edge is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_five_GraphSeqEdge(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_five_GraphSeqEdge
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GraphSeqEdgeHolder **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_five_GraphSeqEdge(GraphSeqEdgeHolder ** list,int i,int j)  
{
    GraphSeqEdgeHolder * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_five_GraphSeqEdge(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_five_GraphSeqEdge which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GraphSeqEdgeHolder **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_five_GraphSeqEdge(GraphSeqEdgeHolder ** list,int left,int right,int (*comp)(GraphSeqEdgeHolder * ,GraphSeqEdgeHolder * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_five_GraphSeqEdge(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_five_GraphSeqEdge (list,++last,i);  
      }  
    swap_five_GraphSeqEdge (list,left,last); 
    qsort_five_GraphSeqEdge(list,left,last-1,comp);  
    qsort_five_GraphSeqEdge(list,last+1,right,comp); 
}    


/* Function:  sort_five_GraphSeqEdge(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_five_GraphSeqEdge
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GraphSeqEdge *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_five_GraphSeqEdge(GraphSeqEdge * obj,int (*comp)(GraphSeqEdgeHolder *, GraphSeqEdgeHolder *)) 
{
    qsort_five_GraphSeqEdge(obj->fivep,0,obj->five_len-1,comp);  
    return;  
}    


/* Function:  expand_five_GraphSeqEdge(obj,len)
 *
 * Descrip:    Really an internal function for add_five_GraphSeqEdge
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeqEdge *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_five_GraphSeqEdge(GraphSeqEdge * obj,int len) 
{


    if( obj->five_maxlen > obj->five_len )   {  
      warn("expand_GraphSeqEdgefive_ called with no need");  
      return TRUE;   
      }  


    if( (obj->fivep = (GraphSeqEdgeHolder ** ) ckrealloc (obj->fivep,sizeof(GraphSeqEdgeHolder *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GraphSeqEdge, returning FALSE"); 
      return FALSE;  
      }  
    obj->five_maxlen = len;  
    return TRUE; 
}    


/* Function:  add_five_GraphSeqEdge(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeqEdge *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_five_GraphSeqEdge(GraphSeqEdge * obj,GraphSeqEdgeHolder * add) 
{
    if( obj->five_len >= obj->five_maxlen)   {  
      if( expand_five_GraphSeqEdge(obj,obj->five_len + GraphSeqEdgeLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->fivep[obj->five_len++]=add; 
    return TRUE; 
}    


/* Function:  flush_five_GraphSeqEdge(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_five_GraphSeqEdge(GraphSeqEdge * obj) 
{
    int i;   


    for(i=0;i<obj->five_len;i++) { /*for i over list length*/ 
      if( obj->fivep[i] != NULL) {  
        free_GraphSeqEdgeHolder(obj->fivep[i]);  
        obj->fivep[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->five_len = 0;   
    return i;    
}    


/* Function:  swap_three_GraphSeqEdge(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_three_GraphSeqEdge
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GraphSeqEdgeHolder **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_three_GraphSeqEdge(GraphSeqEdgeHolder ** list,int i,int j)  
{
    GraphSeqEdgeHolder * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_three_GraphSeqEdge(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_three_GraphSeqEdge which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GraphSeqEdgeHolder **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_three_GraphSeqEdge(GraphSeqEdgeHolder ** list,int left,int right,int (*comp)(GraphSeqEdgeHolder * ,GraphSeqEdgeHolder * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_three_GraphSeqEdge(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_three_GraphSeqEdge (list,++last,i); 
      }  
    swap_three_GraphSeqEdge (list,left,last);    
    qsort_three_GraphSeqEdge(list,left,last-1,comp); 
    qsort_three_GraphSeqEdge(list,last+1,right,comp);    
}    


/* Function:  sort_three_GraphSeqEdge(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_three_GraphSeqEdge
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GraphSeqEdge *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_three_GraphSeqEdge(GraphSeqEdge * obj,int (*comp)(GraphSeqEdgeHolder *, GraphSeqEdgeHolder *)) 
{
    qsort_three_GraphSeqEdge(obj->threep,0,obj->three_len-1,comp);   
    return;  
}    


/* Function:  expand_three_GraphSeqEdge(obj,len)
 *
 * Descrip:    Really an internal function for add_three_GraphSeqEdge
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeqEdge *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_three_GraphSeqEdge(GraphSeqEdge * obj,int len) 
{


    if( obj->three_maxlen > obj->three_len )     {  
      warn("expand_GraphSeqEdgethree_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->threep = (GraphSeqEdgeHolder ** ) ckrealloc (obj->threep,sizeof(GraphSeqEdgeHolder *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GraphSeqEdge, returning FALSE"); 
      return FALSE;  
      }  
    obj->three_maxlen = len; 
    return TRUE; 
}    


/* Function:  add_three_GraphSeqEdge(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeqEdge *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_three_GraphSeqEdge(GraphSeqEdge * obj,GraphSeqEdgeHolder * add) 
{
    if( obj->three_len >= obj->three_maxlen) {  
      if( expand_three_GraphSeqEdge(obj,obj->three_len + GraphSeqEdgeLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->threep[obj->three_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_three_GraphSeqEdge(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_three_GraphSeqEdge(GraphSeqEdge * obj) 
{
    int i;   


    for(i=0;i<obj->three_len;i++)    { /*for i over list length*/ 
      if( obj->threep[i] != NULL)    {  
        free_GraphSeqEdgeHolder(obj->threep[i]); 
        obj->threep[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->three_len = 0;  
    return i;    
}    


/* Function:  GraphSeqEdge_alloc_std(void)
 *
 * Descrip:    Equivalent to GraphSeqEdge_alloc_len(GraphSeqEdgeLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * GraphSeqEdge_alloc_std(void) 
{
    return GraphSeqEdge_alloc_len(GraphSeqEdgeLISTLENGTH);   
}    


/* Function:  GraphSeqEdge_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * GraphSeqEdge_alloc_len(int len) 
{
    GraphSeqEdge * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GraphSeqEdge_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->fivep = (GraphSeqEdgeHolder ** ) ckcalloc (len,sizeof(GraphSeqEdgeHolder *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GraphSeqEdge_alloc_len");    
      return NULL;   
      }  
    out->five_len = 0;   
    out->five_maxlen = len;  


    if((out->threep = (GraphSeqEdgeHolder ** ) ckcalloc (len,sizeof(GraphSeqEdgeHolder *))) == NULL) {  
      warn("Warning, ckcalloc failed in GraphSeqEdge_alloc_len");    
      return NULL;   
      }  
    out->three_len = 0;  
    out->three_maxlen = len; 


    return out;  
}    


/* Function:  hard_link_GraphSeqEdge(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * hard_link_GraphSeqEdge(GraphSeqEdge * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GraphSeqEdge object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GraphSeqEdge_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * GraphSeqEdge_alloc(void) 
{
    GraphSeqEdge * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GraphSeqEdge *) ckalloc (sizeof(GraphSeqEdge))) == NULL)    {  
      warn("GraphSeqEdge_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->id = NULL;  
    out->len = 0;    
    out->seq = NULL; 
    out->fivep = NULL;   
    out->five_len = out->five_maxlen = 0;    
    out->threep = NULL;  
    out->three_len = out->three_maxlen = 0;  


    return out;  
}    


/* Function:  free_GraphSeqEdge(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeqEdge *]
 *
 */
GraphSeqEdge * free_GraphSeqEdge(GraphSeqEdge * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GraphSeqEdge obj. Should be trappable");  
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
    if( obj->id != NULL) 
      ckfree(obj->id);   
    if( obj->seq != NULL)    
      free_Sequence(obj->seq);   
    if( obj->fivep != NULL)  {  
      for(i=0;i<obj->five_len;i++)   {  
        if( obj->fivep[i] != NULL)   
          free_GraphSeqEdgeHolder(obj->fivep[i]);    
        }  
      ckfree(obj->fivep);    
      }  
    if( obj->threep != NULL) {  
      for(i=0;i<obj->three_len;i++)  {  
        if( obj->threep[i] != NULL)  
          free_GraphSeqEdgeHolder(obj->threep[i]);   
        }  
      ckfree(obj->threep);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_free_GraphSeq(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_free_GraphSeq
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GraphSeqEdgeHolder **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_free_GraphSeq(GraphSeqEdgeHolder ** list,int i,int j)  
{
    GraphSeqEdgeHolder * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_free_GraphSeq(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_free_GraphSeq which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GraphSeqEdgeHolder **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_free_GraphSeq(GraphSeqEdgeHolder ** list,int left,int right,int (*comp)(GraphSeqEdgeHolder * ,GraphSeqEdgeHolder * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_free_GraphSeq(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_free_GraphSeq (list,++last,i);  
      }  
    swap_free_GraphSeq (list,left,last); 
    qsort_free_GraphSeq(list,left,last-1,comp);  
    qsort_free_GraphSeq(list,last+1,right,comp); 
}    


/* Function:  sort_free_GraphSeq(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_free_GraphSeq
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GraphSeq *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_free_GraphSeq(GraphSeq * obj,int (*comp)(GraphSeqEdgeHolder *, GraphSeqEdgeHolder *)) 
{
    qsort_free_GraphSeq(obj->free_ends,0,obj->free_len-1,comp);  
    return;  
}    


/* Function:  expand_free_GraphSeq(obj,len)
 *
 * Descrip:    Really an internal function for add_free_GraphSeq
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeq *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_free_GraphSeq(GraphSeq * obj,int len) 
{


    if( obj->free_maxlen > obj->free_len )   {  
      warn("expand_GraphSeqfree_ called with no need");  
      return TRUE;   
      }  


    if( (obj->free_ends = (GraphSeqEdgeHolder ** ) ckrealloc (obj->free_ends,sizeof(GraphSeqEdgeHolder *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GraphSeq, returning FALSE"); 
      return FALSE;  
      }  
    obj->free_maxlen = len;  
    return TRUE; 
}    


/* Function:  add_free_GraphSeq(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeq *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdgeHolder *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_free_GraphSeq(GraphSeq * obj,GraphSeqEdgeHolder * add) 
{
    if( obj->free_len >= obj->free_maxlen)   {  
      if( expand_free_GraphSeq(obj,obj->free_len + GraphSeqLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->free_ends[obj->free_len++]=add; 
    return TRUE; 
}    


/* Function:  flush_free_GraphSeq(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_free_GraphSeq(GraphSeq * obj) 
{
    int i;   


    for(i=0;i<obj->free_len;i++) { /*for i over list length*/ 
      if( obj->free_ends[i] != NULL) {  
        free_GraphSeqEdgeHolder(obj->free_ends[i]);  
        obj->free_ends[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->free_len = 0;   
    return i;    
}    


/* Function:  swap_GraphSeq(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GraphSeq
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GraphSeqEdge **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GraphSeq(GraphSeqEdge ** list,int i,int j)  
{
    GraphSeqEdge * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GraphSeq(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GraphSeq which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GraphSeqEdge **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GraphSeq(GraphSeqEdge ** list,int left,int right,int (*comp)(GraphSeqEdge * ,GraphSeqEdge * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GraphSeq(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GraphSeq (list,++last,i);   
      }  
    swap_GraphSeq (list,left,last);  
    qsort_GraphSeq(list,left,last-1,comp);   
    qsort_GraphSeq(list,last+1,right,comp);  
}    


/* Function:  sort_GraphSeq(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GraphSeq
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GraphSeq *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GraphSeq(GraphSeq * obj,int (*comp)(GraphSeqEdge *, GraphSeqEdge *)) 
{
    qsort_GraphSeq(obj->comp,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_GraphSeq(obj,len)
 *
 * Descrip:    Really an internal function for add_GraphSeq
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeq *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GraphSeq(GraphSeq * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GraphSeq called with no need");   
      return TRUE;   
      }  


    if( (obj->comp = (GraphSeqEdge ** ) ckrealloc (obj->comp,sizeof(GraphSeqEdge *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GraphSeq, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GraphSeq(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GraphSeq *]
 * Arg:        add [OWNER] Object to add to the list [GraphSeqEdge *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GraphSeq(GraphSeq * obj,GraphSeqEdge * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GraphSeq(obj,obj->len + GraphSeqLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->comp[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_GraphSeq(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GraphSeq(GraphSeq * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->comp[i] != NULL)  {  
        free_GraphSeqEdge(obj->comp[i]); 
        obj->comp[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GraphSeq_alloc_std(void)
 *
 * Descrip:    Equivalent to GraphSeq_alloc_len(GraphSeqLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * GraphSeq_alloc_std(void) 
{
    return GraphSeq_alloc_len(GraphSeqLISTLENGTH);   
}    


/* Function:  GraphSeq_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * GraphSeq_alloc_len(int len) 
{
    GraphSeq * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GraphSeq_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->free_ends = (GraphSeqEdgeHolder ** ) ckcalloc (len,sizeof(GraphSeqEdgeHolder *))) == NULL)  {  
      warn("Warning, ckcalloc failed in GraphSeq_alloc_len");    
      return NULL;   
      }  
    out->free_len = 0;   
    out->free_maxlen = len;  


    if((out->comp = (GraphSeqEdge ** ) ckcalloc (len,sizeof(GraphSeqEdge *))) == NULL)   {  
      warn("Warning, ckcalloc failed in GraphSeq_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GraphSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * hard_link_GraphSeq(GraphSeq * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GraphSeq object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GraphSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * GraphSeq_alloc(void) 
{
    GraphSeq * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GraphSeq *) ckalloc (sizeof(GraphSeq))) == NULL)    {  
      warn("GraphSeq_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->free_ends = NULL;   
    out->free_len = out->free_maxlen = 0;    
    out->comp = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_GraphSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GraphSeq *]
 *
 * Return [UNKN ]  Undocumented return value [GraphSeq *]
 *
 */
GraphSeq * free_GraphSeq(GraphSeq * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GraphSeq obj. Should be trappable");  
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
    if( obj->free_ends != NULL)  {  
      for(i=0;i<obj->free_len;i++)   {  
        if( obj->free_ends[i] != NULL)   
          free_GraphSeqEdgeHolder(obj->free_ends[i]);    
        }  
      ckfree(obj->free_ends);    
      }  
    if( obj->comp != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->comp[i] != NULL)    
          free_GraphSeqEdge(obj->comp[i]);   
        }  
      ckfree(obj->comp); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
