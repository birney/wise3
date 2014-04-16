#ifdef _cplusplus
extern "C" {
#endif
#include "peak.h"

/* Function:  write_tagged_MergedPeakList(mpl,ofp)
 *
 * Descrip:    writes tagged peak list
 *
 *
 * Arg:        mpl [UNKN ] Undocumented argument [MergedPeakList *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 47 "peak.dy"
void write_tagged_MergedPeakList(MergedPeakList * mpl,FILE * ofp)
{
  int i;
  int j;
  int count = 1;
 
  for(i=0;i<mpl->len;i++) {
    fprintf(ofp,"%s\t%d\t%d\tmerged.%d",mpl->peak[i]->chr,mpl->peak[i]->start,mpl->peak[i]->end,count);
    for(j=0;j<mpl->peak[i]->len;j++) {
      fprintf(ofp,"\t%s=%f",mpl->peak[i]->comp[j]->para->assay_name,mpl->peak[i]->comp[j]->score);
    }
    fprintf(ofp,"\n");
    count++;
  }
      

}


/* Function:  make_MergedPeakList(lop,stderr_logging_level)
 *
 * Descrip:    makes a merged peak list (single linkage clustering)
 *             from a list of peaklists
 *
 *
 * Arg:                         lop [UNKN ] Undocumented argument [ListofPeakList *]
 * Arg:        stderr_logging_level [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
# line 70 "peak.dy"
MergedPeakList * make_MergedPeakList(ListofPeakList * lop,int stderr_logging_level)
{
  MergedPeakList * out;
  MergedPeak * p;

  PeakList * temp;
  int i;
  int j;
  int count  = 0;
 
  assert(lop != NULL);

  temp = PeakList_alloc_std();

  for(i=0;i<lop->len;i++) {
    for(j=0;j<lop->set[i]->len;j++) {
      add_PeakList(temp,hard_link_Peak(lop->set[i]->peak[j]));
    }
  }

  fprintf(stderr,"Sorting\n");
  sort_PeakList_by_chr_start(temp);

  out = MergedPeakList_alloc_std();
  
  /* first peak into merged set */
  p = MergedPeak_alloc_std();
  p->chr = stringalloc(temp->peak[0]->chr);
  p->start = temp->peak[0]->start;
  p->end   = temp->peak[0]->end;
  add_MergedPeak(p,hard_link_Peak(temp->peak[0]));
  add_MergedPeakList(out,p);

  for(i=1;i<temp->len;i++) {
    if( count %1000 == 0  && stderr_logging_level > 0 ) {
    	fprintf(stderr,"Merging %s %d-%d\n",temp->peak[i]->chr,temp->peak[i]->start,temp->peak[i]->end);
    }			
    count++;
    if( strcmp(temp->peak[i]->chr,p->chr) != 0 ||
	temp->peak[i]->start > p->end ) {
      /* new merged peak */
      
        p = MergedPeak_alloc_std();
	p->chr = stringalloc(temp->peak[i]->chr);
	p->start = temp->peak[i]->start;
	p->end   = temp->peak[i]->end;
	add_MergedPeak(p,hard_link_Peak(temp->peak[i]));
	add_MergedPeakList(out,p);
    } else {
      /* merge in */
      if( p->end < temp->peak[i]->end ) {
	p->end = temp->peak[i]->end;
      }
      add_MergedPeak(p,hard_link_Peak(temp->peak[i]));
    }
  }

  free_PeakList(temp);

  return out;
}

/* Function:  sort_PeakList_by_chr_start(pl)
 *
 * Descrip:    sorts a peak list
 *
 *
 * Arg:        pl [UNKN ] Undocumented argument [PeakList *]
 *
 */
# line 135 "peak.dy"
void sort_PeakList_by_chr_start(PeakList * pl)
{

  sort_PeakList(pl,&comp_Peak);

}

/* Function:  comp_Peak(one,two)
 *
 * Descrip:    internal for sorting
 *             !internal
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [Peak *]
 * Arg:        two [UNKN ] Undocumented argument [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
# line 146 "peak.dy"
int comp_Peak(Peak * one,Peak * two ) 
{
  int cmp;

  if( (cmp = strcmp(one->chr,two->chr)) == 0 ) {
    return one->start - two->start;
  } else {
    return cmp;
  }
}


/* Function:  read_bed_PeakList(ifp,assay_name)
 *
 * Descrip:    creates a peak list from a 4 column bed file
 *
 *
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:        assay_name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
# line 161 "peak.dy"
PeakList * read_bed_PeakList(FILE * ifp,char * assay_name)
{
  PeakPara * para;
  PeakList * out;
  Peak * p;
  char line[MAXLINE];

  int start;
  int end;
  char chr[MAXLINE];
  char id[MAXLINE];
  float score;

  assert(assay_name != NULL);
  assert(ifp != NULL);

  out = PeakList_alloc_std();
  
  para = PeakPara_alloc();
  para->assay_name = stringalloc(assay_name);
  
  while( fgets(line,MAXLINE,ifp) ) {
    if( line[0] == '#' ) {
      continue;
    }
    sscanf(line,"%s %d %d %s %f",chr,&start,&end,id,&score);
    
    p = Peak_alloc();
    p->chr = stringalloc(chr);
    p->id = stringalloc(id);
    p->start = start;
    p->end   = end;
    p->score = (double) score;
    p->para = para;

    add_PeakList(out,p);
  }


  return out;
}


/* Function:  read_npf_PeakList(ifp,assay_name)
 *
 * Descrip:    creates a peak list from a 4 column npf file
 *
 *
 * Arg:               ifp [UNKN ] Undocumented argument [FILE *]
 * Arg:        assay_name [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
# line 207 "peak.dy"
PeakList * read_npf_PeakList(FILE * ifp,char * assay_name)
{
  PeakPara * para;
  PeakList * out;
  Peak * p;
  char line[MAXLINE];

  int start;
  int end;
  char chr[MAXLINE];
  char str;
  int phase;
  char id[MAXLINE];
  float score;

  assert(assay_name != NULL);
  assert(ifp != NULL);

  out = PeakList_alloc_std();
  
  para = PeakPara_alloc();
  para->assay_name = stringalloc(assay_name);
  
  while( fgets(line,MAXLINE,ifp) ) {
    if( line[0] == '#' ) {
      continue;
    }
    sscanf(line,"%s %d %d %s %d %c %f",chr,&start,&end,id,&phase,&str,&score);
    
    p = Peak_alloc();
    p->chr = stringalloc(chr);
    p->id = stringalloc(id);
    p->start = start;
    p->end   = end;
    p->score = (double) score;
    p->para = para;

    add_PeakList(out,p);
  }


  return out;
}

/* Function:  write_4col_bed_PeakList(pl,out)
 *
 * Descrip:    writes out a 4-column bed
 *
 *
 * Arg:         pl [UNKN ] Undocumented argument [PeakList *]
 * Arg:        out [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 254 "peak.dy"
void write_4col_bed_PeakList(PeakList * pl,FILE * out)
{
  int i;

  for(i=0;i<pl->len;i++) {
    fprintf(out,"%s\t%d\t%d\t%s\t%f\n",pl->peak[i]->chr,pl->peak[i]->start,pl->peak[i]->end,pl->peak[i]->id,pl->peak[i]->score);
  }

}
# line 266 "peak.c"
/* Function:  hard_link_PeakPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PeakPara *]
 *
 * Return [UNKN ]  Undocumented return value [PeakPara *]
 *
 */
PeakPara * hard_link_PeakPara(PeakPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PeakPara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PeakPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PeakPara *]
 *
 */
PeakPara * PeakPara_alloc(void) 
{
    PeakPara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PeakPara *) ckalloc (sizeof(PeakPara))) == NULL)    {  
      warn("PeakPara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->assay_name = NULL;  


    return out;  
}    


/* Function:  free_PeakPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PeakPara *]
 *
 * Return [UNKN ]  Undocumented return value [PeakPara *]
 *
 */
PeakPara * free_PeakPara(PeakPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PeakPara obj. Should be trappable");  
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
    if( obj->assay_name != NULL) 
      ckfree(obj->assay_name);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_Peak(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [Peak *]
 *
 */
Peak * hard_link_Peak(Peak * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Peak object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Peak_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Peak *]
 *
 */
Peak * Peak_alloc(void) 
{
    Peak * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Peak *) ckalloc (sizeof(Peak))) == NULL)    {  
      warn("Peak_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->chr = NULL; 
    out->start = 0;  
    out->end = 0;    
    out->score = 0;  
    out->id = NULL;  


    return out;  
}    


/* Function:  free_Peak(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [Peak *]
 *
 */
Peak * free_Peak(Peak * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Peak obj. Should be trappable");  
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
    if( obj->chr != NULL)    
      ckfree(obj->chr);  
    if( obj->id != NULL) 
      ckfree(obj->id);   
    /* obj->para is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_PeakList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_PeakList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Peak **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_PeakList(Peak ** list,int i,int j)  
{
    Peak * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_PeakList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_PeakList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Peak **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_PeakList(Peak ** list,int left,int right,int (*comp)(Peak * ,Peak * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_PeakList(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_PeakList (list,++last,i);   
      }  
    swap_PeakList (list,left,last);  
    qsort_PeakList(list,left,last-1,comp);   
    qsort_PeakList(list,last+1,right,comp);  
}    


/* Function:  sort_PeakList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_PeakList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [PeakList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_PeakList(PeakList * obj,int (*comp)(Peak *, Peak *)) 
{
    qsort_PeakList(obj->peak,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_PeakList(obj,len)
 *
 * Descrip:    Really an internal function for add_PeakList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PeakList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_PeakList(PeakList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_PeakList called with no need");   
      return TRUE;   
      }  


    if( (obj->peak = (Peak ** ) ckrealloc (obj->peak,sizeof(Peak *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_PeakList, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_PeakList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [PeakList *]
 * Arg:        add [OWNER] Object to add to the list [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_PeakList(PeakList * obj,Peak * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_PeakList(obj,obj->len + PeakListLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->peak[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_PeakList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_PeakList(PeakList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->peak[i] != NULL)  {  
        free_Peak(obj->peak[i]); 
        obj->peak[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  PeakList_alloc_std(void)
 *
 * Descrip:    Equivalent to PeakList_alloc_len(PeakListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * PeakList_alloc_std(void) 
{
    return PeakList_alloc_len(PeakListLISTLENGTH);   
}    


/* Function:  PeakList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * PeakList_alloc_len(int len) 
{
    PeakList * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = PeakList_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->peak = (Peak ** ) ckcalloc (len,sizeof(Peak *))) == NULL)   {  
      warn("Warning, ckcalloc failed in PeakList_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_PeakList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * hard_link_PeakList(PeakList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PeakList object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PeakList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * PeakList_alloc(void) 
{
    PeakList * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PeakList *) ckalloc (sizeof(PeakList))) == NULL)    {  
      warn("PeakList_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->para = NULL;    
    out->peak = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_PeakList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [PeakList *]
 *
 */
PeakList * free_PeakList(PeakList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PeakList obj. Should be trappable");  
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
    if( obj->para != NULL)   
      free_PeakPara(obj->para);  
    if( obj->peak != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->peak[i] != NULL)    
          free_Peak(obj->peak[i]);   
        }  
      ckfree(obj->peak); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_ListofPeakList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ListofPeakList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [PeakList **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ListofPeakList(PeakList ** list,int i,int j)  
{
    PeakList * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ListofPeakList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ListofPeakList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [PeakList **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ListofPeakList(PeakList ** list,int left,int right,int (*comp)(PeakList * ,PeakList * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ListofPeakList(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ListofPeakList (list,++last,i); 
      }  
    swap_ListofPeakList (list,left,last);    
    qsort_ListofPeakList(list,left,last-1,comp); 
    qsort_ListofPeakList(list,last+1,right,comp);    
}    


/* Function:  sort_ListofPeakList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ListofPeakList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [ListofPeakList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ListofPeakList(ListofPeakList * obj,int (*comp)(PeakList *, PeakList *)) 
{
    qsort_ListofPeakList(obj->set,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_ListofPeakList(obj,len)
 *
 * Descrip:    Really an internal function for add_ListofPeakList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ListofPeakList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ListofPeakList(ListofPeakList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_ListofPeakList called with no need"); 
      return TRUE;   
      }  


    if( (obj->set = (PeakList ** ) ckrealloc (obj->set,sizeof(PeakList *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_ListofPeakList, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ListofPeakList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [ListofPeakList *]
 * Arg:        add [OWNER] Object to add to the list [PeakList *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ListofPeakList(ListofPeakList * obj,PeakList * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_ListofPeakList(obj,obj->len + ListofPeakListLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->set[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ListofPeakList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [ListofPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ListofPeakList(ListofPeakList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->set[i] != NULL)   {  
        free_PeakList(obj->set[i]);  
        obj->set[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  ListofPeakList_alloc_std(void)
 *
 * Descrip:    Equivalent to ListofPeakList_alloc_len(ListofPeakListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * ListofPeakList_alloc_std(void) 
{
    return ListofPeakList_alloc_len(ListofPeakListLISTLENGTH);   
}    


/* Function:  ListofPeakList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * ListofPeakList_alloc_len(int len) 
{
    ListofPeakList * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = ListofPeakList_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->set = (PeakList ** ) ckcalloc (len,sizeof(PeakList *))) == NULL)    {  
      warn("Warning, ckcalloc failed in ListofPeakList_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_ListofPeakList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [ListofPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * hard_link_ListofPeakList(ListofPeakList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a ListofPeakList object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  ListofPeakList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * ListofPeakList_alloc(void) 
{
    ListofPeakList * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(ListofPeakList *) ckalloc (sizeof(ListofPeakList))) == NULL)    {  
      warn("ListofPeakList_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->set = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_ListofPeakList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [ListofPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [ListofPeakList *]
 *
 */
ListofPeakList * free_ListofPeakList(ListofPeakList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a ListofPeakList obj. Should be trappable");    
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
    if( obj->set != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->set[i] != NULL) 
          free_PeakList(obj->set[i]);    
        }  
      ckfree(obj->set);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_MergedPeak(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_MergedPeak
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Peak **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_MergedPeak(Peak ** list,int i,int j)  
{
    Peak * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_MergedPeak(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_MergedPeak which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Peak **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_MergedPeak(Peak ** list,int left,int right,int (*comp)(Peak * ,Peak * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_MergedPeak(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_MergedPeak (list,++last,i); 
      }  
    swap_MergedPeak (list,left,last);    
    qsort_MergedPeak(list,left,last-1,comp); 
    qsort_MergedPeak(list,last+1,right,comp);    
}    


/* Function:  sort_MergedPeak(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_MergedPeak
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MergedPeak *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_MergedPeak(MergedPeak * obj,int (*comp)(Peak *, Peak *)) 
{
    qsort_MergedPeak(obj->comp,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_MergedPeak(obj,len)
 *
 * Descrip:    Really an internal function for add_MergedPeak
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MergedPeak *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MergedPeak(MergedPeak * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_MergedPeak called with no need"); 
      return TRUE;   
      }  


    if( (obj->comp = (Peak ** ) ckrealloc (obj->comp,sizeof(Peak *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MergedPeak, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_MergedPeak(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MergedPeak *]
 * Arg:        add [OWNER] Object to add to the list [Peak *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_MergedPeak(MergedPeak * obj,Peak * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_MergedPeak(obj,obj->len + MergedPeakLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->comp[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_MergedPeak(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_MergedPeak(MergedPeak * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->comp[i] != NULL)  {  
        free_Peak(obj->comp[i]); 
        obj->comp[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  MergedPeak_alloc_std(void)
 *
 * Descrip:    Equivalent to MergedPeak_alloc_len(MergedPeakLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * MergedPeak_alloc_std(void) 
{
    return MergedPeak_alloc_len(MergedPeakLISTLENGTH);   
}    


/* Function:  MergedPeak_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * MergedPeak_alloc_len(int len) 
{
    MergedPeak * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = MergedPeak_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->comp = (Peak ** ) ckcalloc (len,sizeof(Peak *))) == NULL)   {  
      warn("Warning, ckcalloc failed in MergedPeak_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_MergedPeak(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * hard_link_MergedPeak(MergedPeak * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MergedPeak object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MergedPeak_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * MergedPeak_alloc(void) 
{
    MergedPeak * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MergedPeak *) ckalloc (sizeof(MergedPeak))) == NULL)    {  
      warn("MergedPeak_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->chr = NULL; 
    out->start = 0;  
    out->end = 0;    
    out->comp = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_MergedPeak(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeak *]
 *
 */
MergedPeak * free_MergedPeak(MergedPeak * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MergedPeak obj. Should be trappable");    
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
    if( obj->chr != NULL)    
      ckfree(obj->chr);  
    if( obj->comp != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->comp[i] != NULL)    
          free_Peak(obj->comp[i]);   
        }  
      ckfree(obj->comp); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_MergedPeakList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_MergedPeakList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [MergedPeak **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_MergedPeakList(MergedPeak ** list,int i,int j)  
{
    MergedPeak * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_MergedPeakList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_MergedPeakList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [MergedPeak **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_MergedPeakList(MergedPeak ** list,int left,int right,int (*comp)(MergedPeak * ,MergedPeak * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_MergedPeakList(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_MergedPeakList (list,++last,i); 
      }  
    swap_MergedPeakList (list,left,last);    
    qsort_MergedPeakList(list,left,last-1,comp); 
    qsort_MergedPeakList(list,last+1,right,comp);    
}    


/* Function:  sort_MergedPeakList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_MergedPeakList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [MergedPeakList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_MergedPeakList(MergedPeakList * obj,int (*comp)(MergedPeak *, MergedPeak *)) 
{
    qsort_MergedPeakList(obj->peak,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_MergedPeakList(obj,len)
 *
 * Descrip:    Really an internal function for add_MergedPeakList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MergedPeakList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_MergedPeakList(MergedPeakList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_MergedPeakList called with no need"); 
      return TRUE;   
      }  


    if( (obj->peak = (MergedPeak ** ) ckrealloc (obj->peak,sizeof(MergedPeak *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_MergedPeakList, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_MergedPeakList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [MergedPeakList *]
 * Arg:        add [OWNER] Object to add to the list [MergedPeak *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_MergedPeakList(MergedPeakList * obj,MergedPeak * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_MergedPeakList(obj,obj->len + MergedPeakListLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->peak[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_MergedPeakList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [MergedPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_MergedPeakList(MergedPeakList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->peak[i] != NULL)  {  
        free_MergedPeak(obj->peak[i]);   
        obj->peak[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  MergedPeakList_alloc_std(void)
 *
 * Descrip:    Equivalent to MergedPeakList_alloc_len(MergedPeakListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * MergedPeakList_alloc_std(void) 
{
    return MergedPeakList_alloc_len(MergedPeakListLISTLENGTH);   
}    


/* Function:  MergedPeakList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * MergedPeakList_alloc_len(int len) 
{
    MergedPeakList * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = MergedPeakList_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->peak = (MergedPeak ** ) ckcalloc (len,sizeof(MergedPeak *))) == NULL)   {  
      warn("Warning, ckcalloc failed in MergedPeakList_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_MergedPeakList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [MergedPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * hard_link_MergedPeakList(MergedPeakList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a MergedPeakList object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  MergedPeakList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * MergedPeakList_alloc(void) 
{
    MergedPeakList * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(MergedPeakList *) ckalloc (sizeof(MergedPeakList))) == NULL)    {  
      warn("MergedPeakList_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->peak = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_MergedPeakList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [MergedPeakList *]
 *
 * Return [UNKN ]  Undocumented return value [MergedPeakList *]
 *
 */
MergedPeakList * free_MergedPeakList(MergedPeakList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a MergedPeakList obj. Should be trappable");    
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
    if( obj->peak != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->peak[i] != NULL)    
          free_MergedPeak(obj->peak[i]); 
        }  
      ckfree(obj->peak); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
