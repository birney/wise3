#ifdef _cplusplus
extern "C" {
#endif
#include "genotype.h"

# line 54 "genotype.dy"
BiGenotype * BiGenotype_alloc_slab(BiGenotypeAllocSlab * slab)
{
  if( slab->array_y < GENOTYPE_ALLOC_LINE ) {
    return &(slab->alloc_array[slab->array_pos][slab->array_y++]);
  }

  if( slab->array_pos >= MAX_GENOTYPE_ALLOC ) {
    fatal("Asking for more genotypes than the product of MAX_GENOTYPE_ALLOC and GENOTYPE_ALLOC_LINE. Not feasible");
  }

  slab->array_pos++;
  slab->alloc_array[slab->array_pos] = calloc(GENOTYPE_ALLOC_LINE,sizeof(BiGenotype));
  slab->array_y = 0;

  return &(slab->alloc_array[slab->array_pos][slab->array_y++]);
}


# line 72 "genotype.dy"
BiGenotypeAllocSlab * new_BiGenotypeAllocSlab(void)
{
  BiGenotypeAllocSlab * out;



  out = malloc(sizeof(BiGenotypeAllocSlab));
  out->array_pos = 0;
  out->array_y = 0;
  out->alloc_array[0] = calloc(GENOTYPE_ALLOC_LINE,sizeof(BiGenotype));

  return out;
}

# line 86 "genotype.dy"
BiGenotypeAllocSlab * free_BiGenotypeAllocSlab(BiGenotypeAllocSlab * s)
{
  int i;
  for(i=0;i<= s->array_pos;i++) {
    free(s->alloc_array[i]);
  }

  free(s);
  return NULL;
}

# line 97 "genotype.dy"
BiGenotype * BiGenotype_alloc(void)
{
  BiGenotype * out;

  warn("You should be using a SLAB allocation now!");

  out = malloc(sizeof(BiGenotype));
  
  return out;
}

# line 108 "genotype.dy"
BiGenotype * free_BiGenotype(BiGenotype * bi)
{
  if( bi == NULL ){
    return NULL;
  }
  free(bi);
  return NULL;
}


# line 118 "genotype.dy"
void write_simple_genotype(BiGenotypeSet * bgs,FILE * ofp)
{
  int i;
  int j;

   
  for(i=0;i<bgs->len;i++) {
    fprintf(ofp,"# Locus %s\n",bgs->locus[i]->locus->locus_id);
    for(j=0;j<bgs->locus[i]->len;j++) {
      fprintf(ofp,"%s\t%s\t%s\t",bgs->locus[i]->locus->locus_id,
	      bgs->locus[i]->big[j]->person->person_id,
	      bgs->locus[i]->big[j]->person->pop->pop_name);
      if( bgs->locus[i]->big[j]->type == BiGenotypeHomozygousFirst ) {
	fprintf(ofp,"%c\t%c\n",bgs->locus[i]->locus->first_allele_char,
		bgs->locus[i]->locus->first_allele_char
		);
      } else if ( bgs->locus[i]->big[j]->type == BiGenotypeHomozygousSecond ) {
	fprintf(ofp,"%c\t%c\n",bgs->locus[i]->locus->second_allele_char,
		bgs->locus[i]->locus->second_allele_char
		);
      } else if ( bgs->locus[i]->big[j]->type == BiGenotypeHetrozygous ) {
	fprintf(ofp,"%c\t%c\n",bgs->locus[i]->locus->first_allele_char,
		bgs->locus[i]->locus->second_allele_char
		);
      } else {
	fprintf(ofp,"N\tN\n");
      }
    }
  }

}


# line 151 "genotype.dy"
boolean read_hapmap_genotype_file(BiGenotypeSet * bgs,char * pop_name,FILE * ifp)
{
  BiLocus * bi;
  BiGenotypeLocus * bgl;
  BiGenotype * big;
  Population * pop;
  Person * p;

  int report_number = 2000;
  int count = 0;

  int people_flag[1028];
  Person * people[1028];


  int people_total;

  char buffer[MAX_HAPMAP_LINE];
  char * runner;
  char * allele;
  
  char * old_locus = NULL;

  int i;



  /* first line is description of columns */
  fgets(buffer,MAX_HAPMAP_LINE,ifp);

  runner = strtok(buffer,spacestr);
  i = 1;
  for(;runner != NULL && i < 11;runner = strtok(NULL,spacestr)) {
    i++;
  }

  /* runner should be "QC_code" */
  if( strcmp(runner,"QC_code") != 0 ) {
    warn("HapMap file - first line, column 11 is not QC_code!");
    return FALSE;
  }


  pop = Population_alloc_std();
  pop->pop_name = stringalloc(pop_name);
  add_PopulationSet(bgs->ps,pop);


  /* count people */
  people_total = 0;
  while( (runner = strtok(NULL,spacestr)) != NULL ) {

    if( strstr(runner,".dup") != NULL ) {
      people_flag[people_total] = 0;
    } else {
      people_flag[people_total] = 1;
      people[people_total] = Person_alloc();
      people[people_total]->person_id = stringalloc(runner);
      people[people_total]->pop = pop;
      add_Population(pop,people[people_total]);
    }

    people_total++;
  }

  /* loop over locus lines */

  while( fgets(buffer,MAX_HAPMAP_LINE,ifp) != NULL ) {
    runner = strtok(buffer,spacestr);

    if( old_locus != NULL && strcmp(old_locus,runner) == 0 ) {
      /* only read a locus in once! Duplications for error checking */
      continue;
    }

    if( count%report_number == 0 ) {
      fprintf(stderr,"Handled %d loci, %s\n",count,runner);
    }
    count++;

    allele = strtok(NULL,spacestr);

    if( allele[1] != '/' ) {
      warn("Allele string looks all wrong %s",allele);
    }

    bi  = find_or_new_BiLocus_from_BiLocusFramework(bgs->framework,runner,allele[0],allele[2]);
    /* this way we don't need to track old_locus memory */
    old_locus = bi->locus_id;

    bgl = find_or_new_BiGenotypeLocus_from_BiGenotypeSet(bgs,bi);

    i = 2;

    for(;runner != NULL && i < 11;runner = strtok(NULL,spacestr)) {
      i++;
    }

    /* this should have QC as first two letters */
    if( runner[0] != 'Q' || runner[1] != 'C' ) {
      warn("Column 11 does not have QC ... %s",runner);
    }
    
    
    
    /* now we have the people genotypes */
    
    i = 0;
    for(runner = strtok(NULL,spacestr);runner != NULL; runner = strtok(NULL,spacestr)) {
      if( people_flag[i] == 1 ) {
	big = BiGenotype_alloc_slab(bgs->slab);
	big->type = type_from_hapmap_string_BiLocus(bi,runner);
	big->person = people[i];
	big->phase  = PhaseStatus_Unknown;
	add_BiGenotypeLocus(bgl,big);
      }
      i++;
    }

    
  }
    

  return TRUE;
    
}


# line 279 "genotype.dy"
BiGenotypeLocus * find_or_new_BiGenotypeLocus_from_BiGenotypeSet(BiGenotypeSet * bgs,BiLocus * bi)
{
  BiGenotypeLocus * bgl;
  int i;
  
  for(i=0;i<bgs->len;i++) {
    if( bgs->locus[i]->locus == bi ) {
      return bgs->locus[i];
    }
  }
  
  bgl = BiGenotypeLocus_alloc_std();
  bgl->locus = bi;
  
  add_BiGenotypeSet(bgs,bgl);

  return bgl;
}


# line 299 "genotype.dy"
BiGenotypeSet * read_simple_genotype_file(FILE * ifp)
{
  char buffer[MAXLINE];
  BiGenotypeSet * out;
  BiLocus * bi;
  BiGenotypeLocus * bgl;

  BiGenotype * bg;
  
  Population * pop;
  Person * p;
  


  char * runner;
  char * locus;
  char * person_id;
  char * population_id;
  char * oldlocus = NULL;
  char * allele1;
  char * allele2;

  int no_alleles = 0;

  BiGenotypeType gt;

  char first;
  char second;


  out = new_BiGenotypeSet();


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' ) {
      continue;
    }
    
    runner = strtok(buffer,spacestr);
    locus = runner;
    runner = strtok(NULL,spacestr);
    person_id = runner;
    runner = strtok(NULL,spacestr);
    population_id = runner;
    runner = strtok(NULL,spacestr);
    allele1 = runner;
    runner = strtok(NULL,spacestr);
    allele2 = runner;


    fprintf(stderr,"Got locus %s person %s with old %s\n",locus,person_id,oldlocus);

    if ( allele2 == NULL ) {
      continue;
    }


    /*
     * locus-id
     * person-id
     * allele-1
     * allele-2
     */
    
    fprintf(stderr,"comparing %s vs %s  with %s - %s\n",oldlocus,locus,allele1,allele2);

    if( oldlocus == NULL || strcmp(locus,oldlocus) != 0 ) {
      fprintf(stderr,"new locus, %s\n",locus);
      if( oldlocus != NULL ) {
	fprintf(stderr,"switching locus\n");

	if( no_alleles != 2 ) {
	  warn("Not good. not 2 allele locus but a %d locus with %c %c\n",no_alleles,first,second);
	} else {
	  bi->first_allele_char  = first;
	  bi->second_allele_char = second;
	  /* add this bgl */
	  add_BiGenotypeSet(out,bgl);
	}
      } else {
	fprintf(stderr,"No switch\n");
      }


      /* make new locus data */

       
      bgl = BiGenotypeLocus_alloc_std();
      bi  = BiLocus_alloc();
      bi->locus_id = stringalloc(locus);
      bgl->locus = bi;
      
      add_BiLocusFramework(out->framework,bi);

      oldlocus = bi->locus_id;
      no_alleles = 0;
      first = 'U';
      second = 'U';
       
    }

    
    /* check if we've seen these alleles yet */

    if( no_alleles == 0 ) {

      first = allele1[0];
      fprintf(stderr,"Assigning from 0 alleles %c %s\n",first,allele2);

      if( allele2[0] != first ) {
	fprintf(stderr,"Assigning second allele\n");
	second = allele2[0];
	no_alleles = 2;
      } else {
	no_alleles = 1;
      }
    } else if ( no_alleles == 1 ) {
      if( first != allele1[0] ) {
	second = allele1[0];
	no_alleles = 2;
	if( allele2[0] != first && allele2[0] != second ) {
	  warn("Mismatched alleles %s,%s with %c,%c",allele1,allele2,first,second);
	  continue;
	} else {
	  /* this is fine, allele2 is either first or second allele */
	}
      } else if( first != allele2[0] ) {
	second = allele2[0];
	no_alleles = 2;
      }
    } else if ( (allele1[0] != first && allele1[0] != second) || 
		(allele2[0] != first && allele2[0] != second) ) {
      warn("More alleles than 2 for locus %s, allele %s,%s\n",locus,allele1,allele2);
      continue;
    }


    fprintf(stderr,"Currently got %d alleles, %c %c\n",no_alleles,first,second);

    /* now figure out what genotype we have */
    
    if( allele1[0] == first && allele2[0] == first ) {
      gt = BiGenotypeHomozygousFirst;
    } else if ( allele1[0] == second && allele2[0] == second ) {
      gt = BiGenotypeHomozygousSecond;
    } else {
      /* as we have checked above that it is either first or second, must
	 be a Het */
      gt = BiGenotypeHetrozygous;
    }

    bg = BiGenotype_alloc_slab(out->slab);
    
    /* check we have this population */

    pop = find_or_new_Population_in_PopulationSet(out->ps,population_id);

    bg->person = find_or_new_Person_in_Population(pop,person_id);
    bg->type = gt;

    add_BiGenotypeLocus(bgl,bg);
  }

  if( no_alleles != 2 ) {
    warn("Not good. not 2 allele locus");
  } else {
    bi->first_allele_char  = first;
    bi->second_allele_char = second;
    /* add this bgl */
    add_BiGenotypeSet(out,bgl);
  }

  fprintf(stderr,"Returning with %d loci\n",out->len);

  return out;

}


# line 478 "genotype.dy"
BiGenotypeSet * new_BiGenotypeSet(void)
{
  BiGenotypeSet * out;

  out = BiGenotypeSet_alloc_std();
  out->framework = new_BiLocusFramework();
  out->ps = PopulationSet_alloc_std();
  out->slab = new_BiGenotypeAllocSlab();

  return out;


}

# line 453 "genotype.c"
/* Function:  swap_BiGenotypeLocus(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_BiGenotypeLocus
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [BiGenotype **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_BiGenotypeLocus(BiGenotype ** list,int i,int j)  
{
    BiGenotype * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_BiGenotypeLocus(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_BiGenotypeLocus which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [BiGenotype **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_BiGenotypeLocus(BiGenotype ** list,int left,int right,int (*comp)(BiGenotype * ,BiGenotype * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_BiGenotypeLocus(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_BiGenotypeLocus (list,++last,i);    
      }  
    swap_BiGenotypeLocus (list,left,last);   
    qsort_BiGenotypeLocus(list,left,last-1,comp);    
    qsort_BiGenotypeLocus(list,last+1,right,comp);   
}    


/* Function:  sort_BiGenotypeLocus(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_BiGenotypeLocus
 *
 *
 * Arg:         obj [UNKN ] Object containing list [BiGenotypeLocus *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_BiGenotypeLocus(BiGenotypeLocus * obj,int (*comp)(BiGenotype *, BiGenotype *)) 
{
    qsort_BiGenotypeLocus(obj->big,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_BiGenotypeLocus(obj,len)
 *
 * Descrip:    Really an internal function for add_BiGenotypeLocus
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiGenotypeLocus *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_BiGenotypeLocus(BiGenotypeLocus * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_BiGenotypeLocus called with no need");    
      return TRUE;   
      }  


    if( (obj->big = (BiGenotype ** ) ckrealloc (obj->big,sizeof(BiGenotype *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_BiGenotypeLocus, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_BiGenotypeLocus(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiGenotypeLocus *]
 * Arg:        add [OWNER] Object to add to the list [BiGenotype *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_BiGenotypeLocus(BiGenotypeLocus * obj,BiGenotype * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_BiGenotypeLocus(obj,obj->len + BiGenotypeLocusLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->big[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_BiGenotypeLocus(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_BiGenotypeLocus(BiGenotypeLocus * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->big[i] != NULL)   {  
        free_BiGenotype(obj->big[i]);    
        obj->big[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  BiGenotypeLocus_alloc_std(void)
 *
 * Descrip:    Equivalent to BiGenotypeLocus_alloc_len(BiGenotypeLocusLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * BiGenotypeLocus_alloc_std(void) 
{
    return BiGenotypeLocus_alloc_len(BiGenotypeLocusLISTLENGTH); 
}    


/* Function:  BiGenotypeLocus_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * BiGenotypeLocus_alloc_len(int len) 
{
    BiGenotypeLocus * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = BiGenotypeLocus_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->big = (BiGenotype ** ) ckcalloc (len,sizeof(BiGenotype *))) == NULL)    {  
      warn("Warning, ckcalloc failed in BiGenotypeLocus_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_BiGenotypeLocus(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * hard_link_BiGenotypeLocus(BiGenotypeLocus * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BiGenotypeLocus object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BiGenotypeLocus_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * BiGenotypeLocus_alloc(void) 
{
    BiGenotypeLocus * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BiGenotypeLocus *) ckalloc (sizeof(BiGenotypeLocus))) == NULL)  {  
      warn("BiGenotypeLocus_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->big = NULL; 
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_BiGenotypeLocus(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeLocus *]
 *
 */
BiGenotypeLocus * free_BiGenotypeLocus(BiGenotypeLocus * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a BiGenotypeLocus obj. Should be trappable");   
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
    if( obj->big != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->big[i] != NULL) 
          free_BiGenotype(obj->big[i]);  
        }  
      ckfree(obj->big);  
      }  
    /* obj->locus is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_BiGenotypeSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_BiGenotypeSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [BiGenotypeLocus **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_BiGenotypeSet(BiGenotypeLocus ** list,int i,int j)  
{
    BiGenotypeLocus * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_BiGenotypeSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_BiGenotypeSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [BiGenotypeLocus **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_BiGenotypeSet(BiGenotypeLocus ** list,int left,int right,int (*comp)(BiGenotypeLocus * ,BiGenotypeLocus * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_BiGenotypeSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_BiGenotypeSet (list,++last,i);  
      }  
    swap_BiGenotypeSet (list,left,last); 
    qsort_BiGenotypeSet(list,left,last-1,comp);  
    qsort_BiGenotypeSet(list,last+1,right,comp); 
}    


/* Function:  sort_BiGenotypeSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_BiGenotypeSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [BiGenotypeSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_BiGenotypeSet(BiGenotypeSet * obj,int (*comp)(BiGenotypeLocus *, BiGenotypeLocus *)) 
{
    qsort_BiGenotypeSet(obj->locus,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_BiGenotypeSet(obj,len)
 *
 * Descrip:    Really an internal function for add_BiGenotypeSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiGenotypeSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_BiGenotypeSet(BiGenotypeSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_BiGenotypeSet called with no need");  
      return TRUE;   
      }  


    if( (obj->locus = (BiGenotypeLocus ** ) ckrealloc (obj->locus,sizeof(BiGenotypeLocus *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_BiGenotypeSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_BiGenotypeSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [BiGenotypeSet *]
 * Arg:        add [OWNER] Object to add to the list [BiGenotypeLocus *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_BiGenotypeSet(BiGenotypeSet * obj,BiGenotypeLocus * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_BiGenotypeSet(obj,obj->len + BiGenotypeSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->locus[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_BiGenotypeSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [BiGenotypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_BiGenotypeSet(BiGenotypeSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->locus[i] != NULL) {  
        free_BiGenotypeLocus(obj->locus[i]); 
        obj->locus[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  BiGenotypeSet_alloc_std(void)
 *
 * Descrip:    Equivalent to BiGenotypeSet_alloc_len(BiGenotypeSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * BiGenotypeSet_alloc_std(void) 
{
    return BiGenotypeSet_alloc_len(BiGenotypeSetLISTLENGTH); 
}    


/* Function:  BiGenotypeSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * BiGenotypeSet_alloc_len(int len) 
{
    BiGenotypeSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = BiGenotypeSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->locus = (BiGenotypeLocus ** ) ckcalloc (len,sizeof(BiGenotypeLocus *))) == NULL)    {  
      warn("Warning, ckcalloc failed in BiGenotypeSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_BiGenotypeSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [BiGenotypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * hard_link_BiGenotypeSet(BiGenotypeSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a BiGenotypeSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  BiGenotypeSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * BiGenotypeSet_alloc(void) 
{
    BiGenotypeSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(BiGenotypeSet *) ckalloc (sizeof(BiGenotypeSet))) == NULL)  {  
      warn("BiGenotypeSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->locus = NULL;   
    out->len = out->maxlen = 0;  
    out->framework = NULL;   
    out->ps = NULL;  
    out->slab = NULL;    


    return out;  
}    


/* Function:  free_BiGenotypeSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [BiGenotypeSet *]
 *
 * Return [UNKN ]  Undocumented return value [BiGenotypeSet *]
 *
 */
BiGenotypeSet * free_BiGenotypeSet(BiGenotypeSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a BiGenotypeSet obj. Should be trappable"); 
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
          free_BiGenotypeLocus(obj->locus[i]);   
        }  
      ckfree(obj->locus);    
      }  
    if( obj->framework != NULL)  
      free_BiLocusFramework(obj->framework);     
    if( obj->ps != NULL) 
      free_PopulationSet(obj->ps);   
    if( obj->slab != NULL)   
      free_BiGenotypeAllocSlab(obj->slab);   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
