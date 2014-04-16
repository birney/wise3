#ifdef _cplusplus
extern "C" {
#endif
#include "snplocus.h"


# line 49 "snplocus.dy"
GenoVarSet * only_simple_snp_loci_GenoVarSet(GenoVarSet * gvs)
{
  int i,j;
  GenoVarSet * out;
  GenoVarChr * chr;

  out = GenoVarSet_alloc_len(gvs->len);
  
  for(i=0;i<gvs->ind_len;i++) {
    add_ind_GenoVarSet(out,hard_link_Individual(gvs->ind[i]));
  }

  for(i=0;i<gvs->len;i++) {
    chr = GenoVarChr_alloc_len(gvs->chr[i]->len);
    chr->chr = stringalloc(gvs->chr[i]->chr);
    for(j=0;j<gvs->chr[i]->len;j++) {
      if( gvs->chr[i]->loci[j]->var->locus_type == SIMPLE_SNP_LOCUS ) {
	add_GenoVarChr(chr,hard_link_GenoVarLocus(gvs->chr[i]->loci[j]));
      }
    }
    add_GenoVarSet(out,chr);
  }

  return(out);
}

# line 75 "snplocus.dy"
int individual_index_from_string_GenoVarSet(GenoVarSet * s, char * strain)
{
  int i;

  assert(s != NULL);
  assert(strain != NULL);
  
  for(i=0;i<s->ind_len;i++) {
    if( strcmp(s->ind[i]->id,strain) == 0 ) {
      return i;
    }
  }


  return -1;
    
}

# line 93 "snplocus.dy"
int number_of_simple_snp_loci_GenoVarChr(GenoVarChr * c)
{
  int i;
  int out = 0;

  assert(c != NULL);

  for(i=0;i<c->len;i++) {
    if( c->loci[i]->var->locus_type == SIMPLE_SNP_LOCUS ) {
      out++;
    }
  }

  return(out);

}

# line 110 "snplocus.dy"
void write_sanger_GenoVarSet(GenoVarSet * s,FILE * ofp)
{
  int k;
  int i;
  int j;

  assert(s != NULL);
  assert(ofp != NULL);

  fprintf(ofp,"#CHROM\tPOS\tREF");
  
  for(j=0;j<s->ind_len;j++) {
    fprintf(ofp,"\t%s",s->ind[j]->id);
  }
  fprintf(ofp,"\n");

  for(k=0;k<s->len;k++) {
    for(i=0;i<s->chr[k]->len;i++) {
      auto GenoVarLocus * l;
      l = s->chr[k]->loci[i];
      if( l->var->locus_type == COMPLEX_SNP_LOCUS ) {
	continue;
      }


      fprintf(ofp,"%s\t%ld\t%s",l->var->chr->chr,l->var->pos,l->var->ref_allele);

      for(j=0;j<s->ind_len;j++) {
	if( l->ind[j] == HOMO_REF ) {
	  fprintf(ofp,"\t..");
	} else if( l->ind[j] == HOMO_ALT ) {
	  fprintf(ofp,"\t%c%c",l->var->alt_allele[0],l->var->alt_allele[0]);
	} else {
	  fprintf(ofp,"\t%c%c",l->var->alt_allele[0],l->var->ref_allele[0]);
	}
      }

      fprintf(ofp,"\n");

    }
  }
  

}


# line 156 "snplocus.dy"
GenoVarSet * read_sanger_genotype_file(FILE * ifp)
{
  GenoVarSet * out;
  GenoVarChr * chr = NULL;
  Individual * ind;

  GenoVarLocus * g;
  VarLocus * l;

  char buffer[MAXLINE];
  char minibuf[256];
  
  char ** base;
  char ** brk; 

  int i;

  
  assert(ifp != NULL);

  fgets(buffer,MAXLINE,ifp);

  if( buffer[0] != '#' ) {
    warn("Bad file - expecting a header # line first");
    return NULL;
  }

  base = brk = breakstring(buffer,spacestr);

  if( strcmp(*brk,"#CHROM") != 0 ) {
    warn("Bad file - expecting #CHROM");
    return NULL;
  }
  brk++;

  if( strcmp(*brk,"POS") != 0 ) {
    warn("Bad file - expecting POS");
    return NULL;
  }
  brk++;

  if( strcmp(*brk,"REF") != 0 ) {
    warn("Bad file - expecting POS");
    return NULL;
  }
  brk++;

  out = GenoVarSet_alloc_std();

  while( *brk != NULL ) {
    ind = Individual_alloc();
    ind->id = stringalloc(*brk);
    
    add_ind_GenoVarSet(out,ind);
    brk++;
  }

  ckfree(base);


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {

    base = brk = breakstring(buffer,spacestr);
  
    if( chr == NULL || strcmp(chr->chr,*brk) != 0 ) {
      chr = GenoVarChr_alloc_std();
      chr->chr = stringalloc(*brk);
      add_GenoVarSet(out,chr);
    }
    
    g = GenoVarLocus_alloc();
    l = VarLocus_alloc();
    l->locus_type = SIMPLE_SNP_LOCUS;
    g->var = l;
    g->var->chr = chr;
    g->ind = calloc(out->ind_len,sizeof(char));

    add_GenoVarChr(chr,g);
    brk++;

    l->pos = strtol(*brk,NULL,0);

    brk++;

    l->ref_allele = stringalloc(*brk);
    
    brk++;
    i = 0;
    for(i=0; *brk != NULL; brk++, i++ ) {
      if( strcmp(*brk,"..") == 0 ) {
	g->ind[i] = HOMO_REF;
      } else {
	auto char alt;

	if( (*brk)[0] == l->ref_allele[0] ) {
	  alt = (*brk)[1];
	} else {
	  alt = (*brk)[0];
	}
	
	if( l->alt_allele != NULL && l->alt_allele[0] != alt ) {
	  /*warn("Triallelic state. Not good. Unclear what's going to happen");*/
	  l->locus_type = COMPLEX_SNP_LOCUS;
	} else {
	  minibuf[0] = alt;
	  minibuf[1] = '\0';

	  l->alt_allele = stringalloc(minibuf);
	}

	if( (*brk)[0] == (*brk)[1] == l->ref_allele[0] ) {
	  g->ind[i] = HOMO_REF;
	} else if( (*brk)[0] == (*brk)[1] ) {
	  g->ind[i] = HOMO_ALT;
	} else {
	  g->ind[i] = SIMPLE_HET;
	}

      }
    }

    if( i != out->ind_len ) {
      warn("For %s %ld mismatched index %d vs %d\n",g->var->chr->chr,g->var->pos,i,out->ind_len);
    }
  }
	
  return(out);

}



# line 250 "snplocus.c"
/* Function:  hard_link_VarLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [VarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [VarLocus *]
 *
 */
VarLocus * hard_link_VarLocus(VarLocus * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a VarLocus object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  VarLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [VarLocus *]
 *
 */
VarLocus * VarLocus_alloc(void) 
{
    VarLocus * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(VarLocus *) ckalloc (sizeof(VarLocus))) == NULL)    {  
      warn("VarLocus_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pos = 0;    
    out->ref_allele = NULL;  
    out->alt_allele = NULL;  
    out->locus_type = SIMPLE_SNP_LOCUS;  


    return out;  
}    


/* Function:  free_VarLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [VarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [VarLocus *]
 *
 */
VarLocus * free_VarLocus(VarLocus * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a VarLocus obj. Should be trappable");  
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
    /* obj->chr is linked in */ 
    if( obj->ref_allele != NULL) 
      ckfree(obj->ref_allele);   
    if( obj->alt_allele != NULL) 
      ckfree(obj->alt_allele);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_GenoVarLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenoVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarLocus *]
 *
 */
GenoVarLocus * hard_link_GenoVarLocus(GenoVarLocus * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenoVarLocus object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenoVarLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarLocus *]
 *
 */
GenoVarLocus * GenoVarLocus_alloc(void) 
{
    GenoVarLocus * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenoVarLocus *) ckalloc (sizeof(GenoVarLocus))) == NULL)    {  
      warn("GenoVarLocus_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->var = NULL; 
    out->ind = NULL; 


    return out;  
}    


/* Function:  free_GenoVarLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenoVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarLocus *]
 *
 */
GenoVarLocus * free_GenoVarLocus(GenoVarLocus * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenoVarLocus obj. Should be trappable");  
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_Individual(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [Individual *]
 *
 */
Individual * hard_link_Individual(Individual * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Individual object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Individual_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Individual *]
 *
 */
Individual * Individual_alloc(void) 
{
    Individual * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Individual *) ckalloc (sizeof(Individual))) == NULL)    {  
      warn("Individual_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->id = NULL;  


    return out;  
}    


/* Function:  free_Individual(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [Individual *]
 *
 */
Individual * free_Individual(Individual * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Individual obj. Should be trappable");    
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


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_GenoVarChr(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenoVarChr
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenoVarLocus **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenoVarChr(GenoVarLocus ** list,int i,int j)  
{
    GenoVarLocus * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenoVarChr(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenoVarChr which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenoVarLocus **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenoVarChr(GenoVarLocus ** list,int left,int right,int (*comp)(GenoVarLocus * ,GenoVarLocus * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenoVarChr(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenoVarChr (list,++last,i); 
      }  
    swap_GenoVarChr (list,left,last);    
    qsort_GenoVarChr(list,left,last-1,comp); 
    qsort_GenoVarChr(list,last+1,right,comp);    
}    


/* Function:  sort_GenoVarChr(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenoVarChr
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenoVarChr *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenoVarChr(GenoVarChr * obj,int (*comp)(GenoVarLocus *, GenoVarLocus *)) 
{
    qsort_GenoVarChr(obj->loci,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_GenoVarChr(obj,len)
 *
 * Descrip:    Really an internal function for add_GenoVarChr
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarChr *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenoVarChr(GenoVarChr * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenoVarChr called with no need"); 
      return TRUE;   
      }  


    if( (obj->loci = (GenoVarLocus ** ) ckrealloc (obj->loci,sizeof(GenoVarLocus *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_GenoVarChr, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenoVarChr(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarChr *]
 * Arg:        add [OWNER] Object to add to the list [GenoVarLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenoVarChr(GenoVarChr * obj,GenoVarLocus * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenoVarChr(obj,obj->len + GenoVarChrLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->loci[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_GenoVarChr(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenoVarChr(GenoVarChr * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->loci[i] != NULL)  {  
        free_GenoVarLocus(obj->loci[i]); 
        obj->loci[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  GenoVarChr_alloc_std(void)
 *
 * Descrip:    Equivalent to GenoVarChr_alloc_len(GenoVarChrLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * GenoVarChr_alloc_std(void) 
{
    return GenoVarChr_alloc_len(GenoVarChrLISTLENGTH);   
}    


/* Function:  GenoVarChr_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * GenoVarChr_alloc_len(int len) 
{
    GenoVarChr * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenoVarChr_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->loci = (GenoVarLocus ** ) ckcalloc (len,sizeof(GenoVarLocus *))) == NULL)   {  
      warn("Warning, ckcalloc failed in GenoVarChr_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenoVarChr(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * hard_link_GenoVarChr(GenoVarChr * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenoVarChr object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenoVarChr_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * GenoVarChr_alloc(void) 
{
    GenoVarChr * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenoVarChr *) ckalloc (sizeof(GenoVarChr))) == NULL)    {  
      warn("GenoVarChr_alloc failed ");  
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


/* Function:  free_GenoVarChr(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarChr *]
 *
 */
GenoVarChr * free_GenoVarChr(GenoVarChr * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenoVarChr obj. Should be trappable");    
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
          free_GenoVarLocus(obj->loci[i]);   
        }  
      ckfree(obj->loci); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_GenoVarSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_GenoVarSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [GenoVarChr **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_GenoVarSet(GenoVarChr ** list,int i,int j)  
{
    GenoVarChr * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_GenoVarSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_GenoVarSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [GenoVarChr **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_GenoVarSet(GenoVarChr ** list,int left,int right,int (*comp)(GenoVarChr * ,GenoVarChr * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_GenoVarSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_GenoVarSet (list,++last,i); 
      }  
    swap_GenoVarSet (list,left,last);    
    qsort_GenoVarSet(list,left,last-1,comp); 
    qsort_GenoVarSet(list,last+1,right,comp);    
}    


/* Function:  sort_GenoVarSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_GenoVarSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenoVarSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_GenoVarSet(GenoVarSet * obj,int (*comp)(GenoVarChr *, GenoVarChr *)) 
{
    qsort_GenoVarSet(obj->chr,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_GenoVarSet(obj,len)
 *
 * Descrip:    Really an internal function for add_GenoVarSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_GenoVarSet(GenoVarSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_GenoVarSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->chr = (GenoVarChr ** ) ckrealloc (obj->chr,sizeof(GenoVarChr *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GenoVarSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_GenoVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarSet *]
 * Arg:        add [OWNER] Object to add to the list [GenoVarChr *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_GenoVarSet(GenoVarSet * obj,GenoVarChr * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_GenoVarSet(obj,obj->len + GenoVarSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->chr[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_GenoVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_GenoVarSet(GenoVarSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->chr[i] != NULL)   {  
        free_GenoVarChr(obj->chr[i]);    
        obj->chr[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_ind_GenoVarSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ind_GenoVarSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Individual **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ind_GenoVarSet(Individual ** list,int i,int j)  
{
    Individual * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ind_GenoVarSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ind_GenoVarSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Individual **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ind_GenoVarSet(Individual ** list,int left,int right,int (*comp)(Individual * ,Individual * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ind_GenoVarSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ind_GenoVarSet (list,++last,i); 
      }  
    swap_ind_GenoVarSet (list,left,last);    
    qsort_ind_GenoVarSet(list,left,last-1,comp); 
    qsort_ind_GenoVarSet(list,last+1,right,comp);    
}    


/* Function:  sort_ind_GenoVarSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ind_GenoVarSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [GenoVarSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ind_GenoVarSet(GenoVarSet * obj,int (*comp)(Individual *, Individual *)) 
{
    qsort_ind_GenoVarSet(obj->ind,0,obj->ind_len-1,comp);    
    return;  
}    


/* Function:  expand_ind_GenoVarSet(obj,len)
 *
 * Descrip:    Really an internal function for add_ind_GenoVarSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ind_GenoVarSet(GenoVarSet * obj,int len) 
{


    if( obj->ind_maxlen > obj->ind_len )     {  
      warn("expand_GenoVarSetind_ called with no need"); 
      return TRUE;   
      }  


    if( (obj->ind = (Individual ** ) ckrealloc (obj->ind,sizeof(Individual *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_GenoVarSet, returning FALSE");   
      return FALSE;  
      }  
    obj->ind_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ind_GenoVarSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [GenoVarSet *]
 * Arg:        add [OWNER] Object to add to the list [Individual *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ind_GenoVarSet(GenoVarSet * obj,Individual * add) 
{
    if( obj->ind_len >= obj->ind_maxlen) {  
      if( expand_ind_GenoVarSet(obj,obj->ind_len + GenoVarSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->ind[obj->ind_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ind_GenoVarSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ind_GenoVarSet(GenoVarSet * obj) 
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


/* Function:  GenoVarSet_alloc_std(void)
 *
 * Descrip:    Equivalent to GenoVarSet_alloc_len(GenoVarSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * GenoVarSet_alloc_std(void) 
{
    return GenoVarSet_alloc_len(GenoVarSetLISTLENGTH);   
}    


/* Function:  GenoVarSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * GenoVarSet_alloc_len(int len) 
{
    GenoVarSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = GenoVarSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->chr = (GenoVarChr ** ) ckcalloc (len,sizeof(GenoVarChr *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenoVarSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->ind = (Individual ** ) ckcalloc (len,sizeof(Individual *))) == NULL)    {  
      warn("Warning, ckcalloc failed in GenoVarSet_alloc_len");  
      return NULL;   
      }  
    out->ind_len = 0;    
    out->ind_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_GenoVarSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * hard_link_GenoVarSet(GenoVarSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a GenoVarSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  GenoVarSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * GenoVarSet_alloc(void) 
{
    GenoVarSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(GenoVarSet *) ckalloc (sizeof(GenoVarSet))) == NULL)    {  
      warn("GenoVarSet_alloc failed ");  
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


    return out;  
}    


/* Function:  free_GenoVarSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GenoVarSet *]
 *
 * Return [UNKN ]  Undocumented return value [GenoVarSet *]
 *
 */
GenoVarSet * free_GenoVarSet(GenoVarSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a GenoVarSet obj. Should be trappable");    
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
          free_GenoVarChr(obj->chr[i]);  
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


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
