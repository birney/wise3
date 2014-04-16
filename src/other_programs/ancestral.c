#ifdef _cplusplus
extern "C" {
#endif
#include "ancestral.h"

# line 55 "ancestral.dy"
void dump_AncestralBlockSet(AncestralBlockSet * abs,FILE * ofp)
{
  int i;
  int j;

  assert(abs != NULL);
  assert(ofp != NULL);

  fprintf(ofp,"strain\tchr\tstart\tend\tancestor\tsnp\tdiff\n");
  for(i=0;i<abs->len;i++) {
    for(j=0;j<abs->strain[i]->len;j++) {
      fprintf(ofp,"%s\t%s\t%ld\t%ld\t%s\t%d\t%d\n",
	      abs->strain[i]->ind->id,
	      abs->strain[i]->block[j]->chr,
	      abs->strain[i]->block[j]->start,
	      abs->strain[i]->block[j]->end,
	      abs->strain[i]->block[j]->ai->name,
	      abs->strain[i]->block[j]->snp_count,
	      abs->strain[i]->block[j]->anc_diff
	      );
    }
  }
}


# line 80 "ancestral.dy"
AncestralBlockSet * AncestralBlockSet_from_AncestralVarSet(AncestralVarSet * avs,GenoVarSet * gvs)
{
  int i;
  int j;
  int k;

  int snp_count;
  int anc_diff;

  AncestralBlockSet * out;
  AncestralBlockStrain * strain;
  AncestralBlock * b;


  assert(avs != NULL);
  assert(gvs != NULL);

  out = AncestralBlockSet_alloc_std();

  for(i=0;i<avs->ind_len;i++) {
    strain = AncestralBlockStrain_alloc_std();
    add_AncestralBlockSet(out,strain);

    strain->ind = hard_link_Individual(avs->ind[i]);

    for(j=0;j<avs->len;j++) {
      auto AncestralVarLocus * loci;

      loci = avs->chr[j]->loci[0];
      snp_count = 1;

      for(k=1;k<avs->chr[j]->len;k++) {
	if( loci->ind[i] != avs->chr[j]->loci[k]->ind[i] ) {
	  /* new block from loci to here */
	  b = AncestralBlock_alloc();
	  b->chr = stringalloc(avs->chr[j]->chr);
	  b->start = loci->var->pos;
	  b->end = avs->chr[j]->loci[k-1]->var->pos;
	  b->snp_count = snp_count;
	  b->anc_diff = anc_diff;

	  assert(loci->ind[i] < avs->anc_len);

	  b->ai = hard_link_AncestralIndividual(avs->anc[loci->ind[i]]);
	  add_AncestralBlockStrain(strain,b);


	  loci =  avs->chr[j]->loci[k];
	  snp_count = 0;
	  anc_diff = 0;
	} else {
	  snp_count++;
	}
      }
    }
  }

  return out;
}

# line 140 "ancestral.dy"
void dump_as_pairwise_AncestralVarSet(AncestralVarSet * avs,GenoVarSet * gvs,FILE * ofp)
{
  int i;
  int j;
  
  int a;
  int b;
  

  /* loop over chromosomes first */

  for(i=0;i<avs->len;i++) {
    auto AncestralVarChr * chr;

    chr = avs->chr[i];

    /* loop over pairwise individuals */

    for(a=0;a<avs->ind_len;a++) {
      for(b=a+1;b<avs->ind_len;b++) {
	
	auto AncestralVarLocus * loci = NULL;

	auto int snp_count = 0;
	auto int diff = 0;

	auto int original_a;
	auto int original_b;

	/* the ancestral indexes are not necessarily the same
	   as the original indexes */

	original_a = individual_index_from_string_GenoVarSet(gvs,avs->ind[a]->id);
	original_b = individual_index_from_string_GenoVarSet(gvs,avs->ind[b]->id);
	

	loci = chr->loci[0];
	
	for(j=1;j<chr->len;j++) {
	 
	  if( loci->ind[a] != chr->loci[j]->ind[a] ||
	      loci->ind[b] != chr->loci[j]->ind[b] ) {
	    /* means this block has ended; dump */

	    fprintf(ofp,"%s\t%ld\t%ld\t%s\t%s\t%s\t%s\t%d\t%d\n",
		    chr->chr,
		    loci->var->pos,
		    chr->loci[j]->var->pos,
		    avs->ind[a]->id,
		    avs->ind[b]->id,
		    avs->anc[loci->ind[a]]->name,
		    avs->anc[loci->ind[b]]->name,
		    snp_count,
		    diff
		    );
	    loci = chr->loci[j];
	    snp_count = 0;
	    diff = 0;
	    continue;
	  }

	  /* block continues - look at snp difference */
	  
	  snp_count++;

	  /* we assumme everything matches perfectly */
	  assert(gvs->chr[i]->loci[j]->var == chr->loci[j]->var);
	  
	  if( gvs->chr[i]->loci[j]->ind[original_a] != 
	      gvs->chr[i]->loci[j]->ind[original_b] ) {
	    diff++;
	  } else {
	    /* no difference */
	  }
	} /* end of for loci */
      } /* end of inner loop, ind b */
    } /* end of outer loop, ind a */
  }
	  
}


# line 222 "ancestral.dy"
void write_simple_AncestralVarSet(AncestralVarSet * avs,FILE * ofp)
{
  int i;
  int j;
  int k;


  assert(avs != NULL);
  assert(ofp != NULL);

  fprintf(ofp,"Chr\tPos\tRef\tAlt");
  
  for(i=0;i<avs->ind_len;i++) {
    fprintf(ofp,"\t%s",avs->ind[i]->id);
  }

  fprintf(ofp,"\n");

  for(i=0;i<avs->len;i++) {
    auto AncestralVarChr * chr = avs->chr[i];
    for(j=0;j<chr->len;j++) {
      auto AncestralVarLocus * l = chr->loci[j];

      assert(l != NULL);
      assert(l->var != NULL);

      fprintf(ofp,"%s\t%ld\t%s\t%s",chr->chr,l->var->pos,l->var->ref_allele,l->var->alt_allele);
      
      for(k=0;k<avs->ind_len;k++) {
	auto char index = l->ind[k];
	auto char c;

	if( index == -1 ) {
	  c = 'U';
	} else {
	  c = avs->anc[index]->single_letter;
	}

	fprintf(ofp,"\t%c",c);
      }

      fprintf(ofp,"\n");
    }
  }

}


# line 270 "ancestral.dy"
FrameSet * AncestralVarChr_as_FrameSet(AncestralVarSet * avs,AncestralVarChr * avc)
{
  FrameSet * fs;
  int i;
  int j;
  AncestralVarLocus * l;
  BoxGlyph * box;

  long int start_x;
  double scale_factor = 20000.0;

  Colour * set[256];
  char buffer[128];

  set[0] = std_Colour("red");
  set[1] = std_Colour("blue");
  set[2] = std_Colour("green");



  start_x = avc->loci[0]->var->pos;

  fs = FrameSet_alloc_std();

  sprintf(buffer,"%ld",start_x);
  box = text_BoxGlyph(50,10,buffer,7.0);

  add_FrameSet(fs,box);

  sprintf(buffer,"%ld",avc->loci[avc->len-1]->var->pos);

  box = text_BoxGlyph(50+(avc->loci[avc->len-1]->var->pos - start_x)/scale_factor,10,buffer,7.0);

  add_FrameSet(fs,box);

  for(i=0;i<avs->ind_len;i++) {
    l = avc->loci[0];

    box = text_BoxGlyph(15,(i+2)*10,avs->ind[i]->id,10.0);
    add_FrameSet(fs,box);
    for(j=1;j<avc->len;j++) {
      if( l->ind[i] != avc->loci[j]->ind[i] ) {
	/* make a rectangle */
	box = filled_rect_BoxGlyph(50 + ((l->var->pos - start_x) / scale_factor),
				   (i+2)*10,
				   50 + ((avc->loci[j]->var->pos - start_x) / scale_factor),
				   (i+2.9) * 10,
				   set[l->ind[i]],
				   NULL);
	l = avc->loci[j];
	add_FrameSet(fs,box);
      }
    }
  }

      
  return(fs);
}

# line 284 "ancestral.c"
/* Function:  hard_link_AncestralVarLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarLocus *]
 *
 */
AncestralVarLocus * hard_link_AncestralVarLocus(AncestralVarLocus * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralVarLocus object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralVarLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarLocus *]
 *
 */
AncestralVarLocus * AncestralVarLocus_alloc(void) 
{
    AncestralVarLocus * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralVarLocus *) ckalloc (sizeof(AncestralVarLocus))) == NULL)  {  
      warn("AncestralVarLocus_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->var = NULL; 
    out->ind = NULL; 
    out->anc = NULL; 


    return out;  
}    


/* Function:  free_AncestralVarLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarLocus *]
 *
 */
AncestralVarLocus * free_AncestralVarLocus(AncestralVarLocus * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralVarLocus obj. Should be trappable"); 
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
    if( obj->var != NULL)    
      free_VarLocus(obj->var);   
    if( obj->ind != NULL)    
      ckfree(obj->ind);  
    if( obj->anc != NULL)    
      ckfree(obj->anc);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AncestralVarChr(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AncestralVarChr
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AncestralVarLocus **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AncestralVarChr(AncestralVarLocus ** list,int i,int j)  
{
    AncestralVarLocus * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AncestralVarChr(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AncestralVarChr which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AncestralVarLocus **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AncestralVarChr(AncestralVarLocus ** list,int left,int right,int (*comp)(AncestralVarLocus * ,AncestralVarLocus * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AncestralVarChr(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AncestralVarChr (list,++last,i);    
      }  
    swap_AncestralVarChr (list,left,last);   
    qsort_AncestralVarChr(list,left,last-1,comp);    
    qsort_AncestralVarChr(list,last+1,right,comp);   
}    


/* Function:  sort_AncestralVarChr(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AncestralVarChr
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AncestralVarChr *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AncestralVarChr(AncestralVarChr * obj,int (*comp)(AncestralVarLocus *, AncestralVarLocus *)) 
{
    qsort_AncestralVarChr(obj->loci,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_AncestralVarChr(obj,len)
 *
 * Descrip:    Really an internal function for add_AncestralVarChr
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarChr *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AncestralVarChr(AncestralVarChr * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AncestralVarChr called with no need");    
      return TRUE;   
      }  


    if( (obj->loci = (AncestralVarLocus ** ) ckrealloc (obj->loci,sizeof(AncestralVarLocus *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AncestralVarChr, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AncestralVarChr(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarChr *]
 * Arg:        add [OWNER] Object to add to the list [AncestralVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AncestralVarChr(AncestralVarChr * obj,AncestralVarLocus * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AncestralVarChr(obj,obj->len + AncestralVarChrLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->loci[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_AncestralVarChr(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AncestralVarChr(AncestralVarChr * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->loci[i] != NULL)  {  
        free_AncestralVarLocus(obj->loci[i]);    
        obj->loci[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AncestralVarChr_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralVarChr_alloc_len(AncestralVarChrLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * AncestralVarChr_alloc_std(void) 
{
    return AncestralVarChr_alloc_len(AncestralVarChrLISTLENGTH); 
}    


/* Function:  AncestralVarChr_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * AncestralVarChr_alloc_len(int len) 
{
    AncestralVarChr * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AncestralVarChr_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->loci = (AncestralVarLocus ** ) ckcalloc (len,sizeof(AncestralVarLocus *))) == NULL) {  
      warn("Warning, ckcalloc failed in AncestralVarChr_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AncestralVarChr(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * hard_link_AncestralVarChr(AncestralVarChr * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralVarChr object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralVarChr_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * AncestralVarChr_alloc(void) 
{
    AncestralVarChr * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralVarChr *) ckalloc (sizeof(AncestralVarChr))) == NULL)  {  
      warn("AncestralVarChr_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->chr = NULL; 
    out->loci = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_AncestralVarChr(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarChr *]
 *
 */
AncestralVarChr * free_AncestralVarChr(AncestralVarChr * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralVarChr obj. Should be trappable");   
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
    if( obj->loci != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->loci[i] != NULL)    
          free_AncestralVarLocus(obj->loci[i]);  
        }  
      ckfree(obj->loci); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AncestralIndividual(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralIndividual *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralIndividual *]
 *
 */
AncestralIndividual * hard_link_AncestralIndividual(AncestralIndividual * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralIndividual object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralIndividual_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralIndividual *]
 *
 */
AncestralIndividual * AncestralIndividual_alloc(void) 
{
    AncestralIndividual * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralIndividual *) ckalloc (sizeof(AncestralIndividual))) == NULL)  {  
      warn("AncestralIndividual_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->single_letter = 'U';    


    return out;  
}    


/* Function:  free_AncestralIndividual(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralIndividual *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralIndividual *]
 *
 */
AncestralIndividual * free_AncestralIndividual(AncestralIndividual * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralIndividual obj. Should be trappable");   
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
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AncestralVarSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AncestralVarSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AncestralVarChr **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AncestralVarSet(AncestralVarChr ** list,int i,int j)  
{
    AncestralVarChr * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AncestralVarSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AncestralVarSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AncestralVarChr **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AncestralVarSet(AncestralVarChr ** list,int left,int right,int (*comp)(AncestralVarChr * ,AncestralVarChr * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AncestralVarSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AncestralVarSet (list,++last,i);    
      }  
    swap_AncestralVarSet (list,left,last);   
    qsort_AncestralVarSet(list,left,last-1,comp);    
    qsort_AncestralVarSet(list,last+1,right,comp);   
}    


/* Function:  sort_AncestralVarSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AncestralVarSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AncestralVarSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AncestralVarSet(AncestralVarSet * obj,int (*comp)(AncestralVarChr *, AncestralVarChr *)) 
{
    qsort_AncestralVarSet(obj->chr,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_AncestralVarSet(obj,len)
 *
 * Descrip:    Really an internal function for add_AncestralVarSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AncestralVarSet(AncestralVarSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AncestralVarSet called with no need");    
      return TRUE;   
      }  


    if( (obj->chr = (AncestralVarChr ** ) ckrealloc (obj->chr,sizeof(AncestralVarChr *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AncestralVarSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AncestralVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        add [OWNER] Object to add to the list [AncestralVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AncestralVarSet(AncestralVarSet * obj,AncestralVarChr * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AncestralVarSet(obj,obj->len + AncestralVarSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->chr[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_AncestralVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AncestralVarSet(AncestralVarSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->chr[i] != NULL)   {  
        free_AncestralVarChr(obj->chr[i]);   
        obj->chr[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_ind_AncestralVarSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ind_AncestralVarSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Individual **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ind_AncestralVarSet(Individual ** list,int i,int j)  
{
    Individual * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ind_AncestralVarSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ind_AncestralVarSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Individual **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ind_AncestralVarSet(Individual ** list,int left,int right,int (*comp)(Individual * ,Individual * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ind_AncestralVarSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ind_AncestralVarSet (list,++last,i);    
      }  
    swap_ind_AncestralVarSet (list,left,last);   
    qsort_ind_AncestralVarSet(list,left,last-1,comp);    
    qsort_ind_AncestralVarSet(list,last+1,right,comp);   
}    


/* Function:  sort_ind_AncestralVarSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ind_AncestralVarSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AncestralVarSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ind_AncestralVarSet(AncestralVarSet * obj,int (*comp)(Individual *, Individual *)) 
{
    qsort_ind_AncestralVarSet(obj->ind,0,obj->ind_len-1,comp);   
    return;  
}    


/* Function:  expand_ind_AncestralVarSet(obj,len)
 *
 * Descrip:    Really an internal function for add_ind_AncestralVarSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ind_AncestralVarSet(AncestralVarSet * obj,int len) 
{


    if( obj->ind_maxlen > obj->ind_len )     {  
      warn("expand_AncestralVarSetind_ called with no need");    
      return TRUE;   
      }  


    if( (obj->ind = (Individual ** ) ckrealloc (obj->ind,sizeof(Individual *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AncestralVarSet, returning FALSE");  
      return FALSE;  
      }  
    obj->ind_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ind_AncestralVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        add [OWNER] Object to add to the list [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ind_AncestralVarSet(AncestralVarSet * obj,Individual * add) 
{
    if( obj->ind_len >= obj->ind_maxlen) {  
      if( expand_ind_AncestralVarSet(obj,obj->ind_len + AncestralVarSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->ind[obj->ind_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ind_AncestralVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ind_AncestralVarSet(AncestralVarSet * obj) 
{
    int i;   


    for(i=0;i<obj->ind_len;i++)  { /*for i over list length*/ 
      if( obj->ind[i] != NULL)   {  
        free_Individual(obj->ind[i]);    
        obj->ind[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->ind_len = 0;    
    return i;    
}    


/* Function:  swap_anc_AncestralVarSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_anc_AncestralVarSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AncestralIndividual **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_anc_AncestralVarSet(AncestralIndividual ** list,int i,int j)  
{
    AncestralIndividual * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_anc_AncestralVarSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_anc_AncestralVarSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AncestralIndividual **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_anc_AncestralVarSet(AncestralIndividual ** list,int left,int right,int (*comp)(AncestralIndividual * ,AncestralIndividual * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_anc_AncestralVarSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_anc_AncestralVarSet (list,++last,i);    
      }  
    swap_anc_AncestralVarSet (list,left,last);   
    qsort_anc_AncestralVarSet(list,left,last-1,comp);    
    qsort_anc_AncestralVarSet(list,last+1,right,comp);   
}    


/* Function:  sort_anc_AncestralVarSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_anc_AncestralVarSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AncestralVarSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_anc_AncestralVarSet(AncestralVarSet * obj,int (*comp)(AncestralIndividual *, AncestralIndividual *)) 
{
    qsort_anc_AncestralVarSet(obj->anc,0,obj->anc_len-1,comp);   
    return;  
}    


/* Function:  expand_anc_AncestralVarSet(obj,len)
 *
 * Descrip:    Really an internal function for add_anc_AncestralVarSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_anc_AncestralVarSet(AncestralVarSet * obj,int len) 
{


    if( obj->anc_maxlen > obj->anc_len )     {  
      warn("expand_AncestralVarSetanc_ called with no need");    
      return TRUE;   
      }  


    if( (obj->anc = (AncestralIndividual ** ) ckrealloc (obj->anc,sizeof(AncestralIndividual *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AncestralVarSet, returning FALSE");  
      return FALSE;  
      }  
    obj->anc_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_anc_AncestralVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralVarSet *]
 * Arg:        add [OWNER] Object to add to the list [AncestralIndividual *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_anc_AncestralVarSet(AncestralVarSet * obj,AncestralIndividual * add) 
{
    if( obj->anc_len >= obj->anc_maxlen) {  
      if( expand_anc_AncestralVarSet(obj,obj->anc_len + AncestralVarSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->anc[obj->anc_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_anc_AncestralVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_anc_AncestralVarSet(AncestralVarSet * obj) 
{
    int i;   


    for(i=0;i<obj->anc_len;i++)  { /*for i over list length*/ 
      if( obj->anc[i] != NULL)   {  
        free_AncestralIndividual(obj->anc[i]);   
        obj->anc[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->anc_len = 0;    
    return i;    
}    


/* Function:  AncestralVarSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralVarSet_alloc_len(AncestralVarSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * AncestralVarSet_alloc_std(void) 
{
    return AncestralVarSet_alloc_len(AncestralVarSetLISTLENGTH); 
}    


/* Function:  AncestralVarSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * AncestralVarSet_alloc_len(int len) 
{
    AncestralVarSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AncestralVarSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->chr = (AncestralVarChr ** ) ckcalloc (len,sizeof(AncestralVarChr *))) == NULL)  {  
      warn("Warning, ckcalloc failed in AncestralVarSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->ind = (Individual ** ) ckcalloc (len,sizeof(Individual *))) == NULL)    {  
      warn("Warning, ckcalloc failed in AncestralVarSet_alloc_len"); 
      return NULL;   
      }  
    out->ind_len = 0;    
    out->ind_maxlen = len;   


    if((out->anc = (AncestralIndividual ** ) ckcalloc (len,sizeof(AncestralIndividual *))) == NULL)  {  
      warn("Warning, ckcalloc failed in AncestralVarSet_alloc_len"); 
      return NULL;   
      }  
    out->anc_len = 0;    
    out->anc_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AncestralVarSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * hard_link_AncestralVarSet(AncestralVarSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralVarSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralVarSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * AncestralVarSet_alloc(void) 
{
    AncestralVarSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralVarSet *) ckalloc (sizeof(AncestralVarSet))) == NULL)  {  
      warn("AncestralVarSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->chr = NULL; 
    out->len = out->maxlen = 0;  
    out->ind = NULL; 
    out->ind_len = out->ind_maxlen = 0;  
    out->anc = NULL; 
    out->anc_len = out->anc_maxlen = 0;  


    return out;  
}    


/* Function:  free_AncestralVarSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralVarSet *]
 *
 */
AncestralVarSet * free_AncestralVarSet(AncestralVarSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralVarSet obj. Should be trappable");   
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
    if( obj->chr != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->chr[i] != NULL) 
          free_AncestralVarChr(obj->chr[i]); 
        }  
      ckfree(obj->chr);  
      }  
    if( obj->ind != NULL)    {  
      for(i=0;i<obj->ind_len;i++)    {  
        if( obj->ind[i] != NULL) 
          free_Individual(obj->ind[i]);  
        }  
      ckfree(obj->ind);  
      }  
    if( obj->anc != NULL)    {  
      for(i=0;i<obj->anc_len;i++)    {  
        if( obj->anc[i] != NULL) 
          free_AncestralIndividual(obj->anc[i]); 
        }  
      ckfree(obj->anc);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AncestralBlock(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlock *]
 *
 */
AncestralBlock * hard_link_AncestralBlock(AncestralBlock * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralBlock object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralBlock_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlock *]
 *
 */
AncestralBlock * AncestralBlock_alloc(void) 
{
    AncestralBlock * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralBlock *) ckalloc (sizeof(AncestralBlock))) == NULL)    {  
      warn("AncestralBlock_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->chr = NULL; 
    out->start = 0;  
    out->end = 0;    
    out->snp_count = 0;  
    out->anc_diff = 0;   
    out->ai = NULL;  


    return out;  
}    


/* Function:  free_AncestralBlock(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralBlock *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlock *]
 *
 */
AncestralBlock * free_AncestralBlock(AncestralBlock * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralBlock obj. Should be trappable");    
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
    if( obj->ai != NULL) 
      free_AncestralIndividual(obj->ai);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AncestralBlockStrain(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AncestralBlockStrain
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AncestralBlock **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AncestralBlockStrain(AncestralBlock ** list,int i,int j)  
{
    AncestralBlock * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AncestralBlockStrain(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AncestralBlockStrain which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AncestralBlock **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AncestralBlockStrain(AncestralBlock ** list,int left,int right,int (*comp)(AncestralBlock * ,AncestralBlock * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AncestralBlockStrain(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AncestralBlockStrain (list,++last,i);   
      }  
    swap_AncestralBlockStrain (list,left,last);  
    qsort_AncestralBlockStrain(list,left,last-1,comp);   
    qsort_AncestralBlockStrain(list,last+1,right,comp);  
}    


/* Function:  sort_AncestralBlockStrain(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AncestralBlockStrain
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AncestralBlockStrain *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AncestralBlockStrain(AncestralBlockStrain * obj,int (*comp)(AncestralBlock *, AncestralBlock *)) 
{
    qsort_AncestralBlockStrain(obj->block,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_AncestralBlockStrain(obj,len)
 *
 * Descrip:    Really an internal function for add_AncestralBlockStrain
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralBlockStrain *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AncestralBlockStrain(AncestralBlockStrain * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AncestralBlockStrain called with no need");   
      return TRUE;   
      }  


    if( (obj->block = (AncestralBlock ** ) ckrealloc (obj->block,sizeof(AncestralBlock *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_AncestralBlockStrain, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AncestralBlockStrain(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralBlockStrain *]
 * Arg:        add [OWNER] Object to add to the list [AncestralBlock *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AncestralBlockStrain(AncestralBlockStrain * obj,AncestralBlock * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AncestralBlockStrain(obj,obj->len + AncestralBlockStrainLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->block[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_AncestralBlockStrain(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AncestralBlockStrain(AncestralBlockStrain * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->block[i] != NULL) {  
        free_AncestralBlock(obj->block[i]);  
        obj->block[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AncestralBlockStrain_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralBlockStrain_alloc_len(AncestralBlockStrainLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * AncestralBlockStrain_alloc_std(void) 
{
    return AncestralBlockStrain_alloc_len(AncestralBlockStrainLISTLENGTH);   
}    


/* Function:  AncestralBlockStrain_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * AncestralBlockStrain_alloc_len(int len) 
{
    AncestralBlockStrain * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AncestralBlockStrain_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->block = (AncestralBlock ** ) ckcalloc (len,sizeof(AncestralBlock *))) == NULL)  {  
      warn("Warning, ckcalloc failed in AncestralBlockStrain_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AncestralBlockStrain(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * hard_link_AncestralBlockStrain(AncestralBlockStrain * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralBlockStrain object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralBlockStrain_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * AncestralBlockStrain_alloc(void) 
{
    AncestralBlockStrain * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralBlockStrain *) ckalloc (sizeof(AncestralBlockStrain))) == NULL)    {  
      warn("AncestralBlockStrain_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->block = NULL;   
    out->len = out->maxlen = 0;  
    out->ind = NULL; 


    return out;  
}    


/* Function:  free_AncestralBlockStrain(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockStrain *]
 *
 */
AncestralBlockStrain * free_AncestralBlockStrain(AncestralBlockStrain * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralBlockStrain obj. Should be trappable");  
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
    if( obj->block != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->block[i] != NULL)   
          free_AncestralBlock(obj->block[i]);    
        }  
      ckfree(obj->block);    
      }  
    if( obj->ind != NULL)    
      free_Individual(obj->ind);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_AncestralBlockSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_AncestralBlockSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [AncestralBlockStrain **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_AncestralBlockSet(AncestralBlockStrain ** list,int i,int j)  
{
    AncestralBlockStrain * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_AncestralBlockSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_AncestralBlockSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [AncestralBlockStrain **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_AncestralBlockSet(AncestralBlockStrain ** list,int left,int right,int (*comp)(AncestralBlockStrain * ,AncestralBlockStrain * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_AncestralBlockSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_AncestralBlockSet (list,++last,i);  
      }  
    swap_AncestralBlockSet (list,left,last); 
    qsort_AncestralBlockSet(list,left,last-1,comp);  
    qsort_AncestralBlockSet(list,last+1,right,comp); 
}    


/* Function:  sort_AncestralBlockSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_AncestralBlockSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [AncestralBlockSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_AncestralBlockSet(AncestralBlockSet * obj,int (*comp)(AncestralBlockStrain *, AncestralBlockStrain *)) 
{
    qsort_AncestralBlockSet(obj->strain,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_AncestralBlockSet(obj,len)
 *
 * Descrip:    Really an internal function for add_AncestralBlockSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralBlockSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_AncestralBlockSet(AncestralBlockSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_AncestralBlockSet called with no need");  
      return TRUE;   
      }  


    if( (obj->strain = (AncestralBlockStrain ** ) ckrealloc (obj->strain,sizeof(AncestralBlockStrain *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_AncestralBlockSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_AncestralBlockSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [AncestralBlockSet *]
 * Arg:        add [OWNER] Object to add to the list [AncestralBlockStrain *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_AncestralBlockSet(AncestralBlockSet * obj,AncestralBlockStrain * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_AncestralBlockSet(obj,obj->len + AncestralBlockSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->strain[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_AncestralBlockSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [AncestralBlockSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_AncestralBlockSet(AncestralBlockSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->strain[i] != NULL)    {  
        free_AncestralBlockStrain(obj->strain[i]);   
        obj->strain[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  AncestralBlockSet_alloc_std(void)
 *
 * Descrip:    Equivalent to AncestralBlockSet_alloc_len(AncestralBlockSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * AncestralBlockSet_alloc_std(void) 
{
    return AncestralBlockSet_alloc_len(AncestralBlockSetLISTLENGTH); 
}    


/* Function:  AncestralBlockSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * AncestralBlockSet_alloc_len(int len) 
{
    AncestralBlockSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = AncestralBlockSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->strain = (AncestralBlockStrain ** ) ckcalloc (len,sizeof(AncestralBlockStrain *))) == NULL) {  
      warn("Warning, ckcalloc failed in AncestralBlockSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_AncestralBlockSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AncestralBlockSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * hard_link_AncestralBlockSet(AncestralBlockSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AncestralBlockSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AncestralBlockSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * AncestralBlockSet_alloc(void) 
{
    AncestralBlockSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AncestralBlockSet *) ckalloc (sizeof(AncestralBlockSet))) == NULL)  {  
      warn("AncestralBlockSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->strain = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_AncestralBlockSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AncestralBlockSet *]
 *
 * Return [UNKN ]  Undocumented return value [AncestralBlockSet *]
 *
 */
AncestralBlockSet * free_AncestralBlockSet(AncestralBlockSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AncestralBlockSet obj. Should be trappable"); 
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
    if( obj->strain != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->strain[i] != NULL)  
          free_AncestralBlockStrain(obj->strain[i]); 
        }  
      ckfree(obj->strain);   
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
