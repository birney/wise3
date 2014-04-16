#ifdef _cplusplus
extern "C" {
#endif
#include "two_track.h"

static const char * two_track_alphabet = "ACGT";

/* Function:  write_TwoTrackSetStats(st,ofp)
 *
 * Descrip:    writes out stats
 *
 *
 * Arg:         st [UNKN ] Undocumented argument [TwoTrackSetStats *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 60 "two_track.dy"
void write_TwoTrackSetStats(TwoTrackSetStats * st,FILE * ofp)
{
  int i;

  assert(ofp != NULL);
  assert(st != NULL);

  fprintf(ofp,"Total reads       %ld\n",st->reads);
  fprintf(ofp,"Total main bases  %ld Avg: (%.2f)\n",st->main_bases,(double) st->main_bases / (double) st->reads);
  fprintf(ofp,"Total side bases  %ld Avg: (%.2f)\n",st->side_bases,(double) st->side_bases / (double) st->reads);
  fprintf(ofp,"Main avg likelihood %.6f\n",st->main_avg_likelihood);
  fprintf(ofp,"Side avg likelihood %.6f\n",st->side_avg_likelihood);

  fprintf(ofp,"  Called base distribution\n");
  for(i=0;i<4;i++) {
    fprintf(ofp,"   %c %.4f",char_from_base(i),st->called_main_bases[i]/(double)st->main_bases);
  }

  fprintf(ofp,"\n  Most likely main base distribution\n");
  for(i=0;i<4;i++) {
    fprintf(ofp,"   %c %.4f",char_from_base(i),st->most_likely_main_bases[i]/(double)st->main_bases);
  }

  fprintf(ofp,"\n  Most likely side base distribution\n");
  for(i=0;i<4;i++) {
    fprintf(ofp,"   %c %.4f",char_from_base(i),st->most_likely_side_bases[i]/(double)st->main_bases);
  }

}

/* Function:  TwoTrackSetStats_from_TwoTrackSet(set)
 *
 * Descrip:    Makes stats from a set
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
# line 93 "two_track.dy"
TwoTrackSetStats * TwoTrackSetStats_from_TwoTrackSet(TwoTrackSet * set)
{
  int i,j,k;

  TwoTrackSetStats * out;
  long double total_main_lik = 0.0;
  long double total_side_lik = 0.0;


  out = TwoTrackSetStats_alloc();
  out->reads = set->len;
  out->main_bases = 0;
  out->side_bases = 0;

  out->called_main_bases[0] = out->called_main_bases[1] = out->called_main_bases[2] = out->called_main_bases[3] = 0;

  out->most_likely_main_bases[0] = out->most_likely_main_bases[1] = out->most_likely_main_bases[2] = out->most_likely_main_bases[3] = 0;

  out->most_likely_side_bases[0] = out->most_likely_side_bases[1] = out->most_likely_side_bases[2] = out->most_likely_side_bases[3] = 0;


  for(i=0;i<set->len;i++) {
    for(j=0;j<set->read[i]->seq->len;j++) {
      out->called_main_bases[base_from_char(set->read[i]->seq->seq[j])]++;
    }

    for(j=0;j<set->read[i]->len;j++) {
      if( set->read[i]->unit[j]->type == TwoTrackUnit_main ) {
	out->main_bases++;
	out->most_likely_main_bases[most_likely_base_TwoTrackUnit(set->read[i]->unit[j])]++;

	for(k=0;k<4;k++) {
	  total_main_lik += set->read[i]->unit[j]->emission[k];
	}
      } else if( set->read[i]->unit[j]->type == TwoTrackUnit_side ) {
	out->most_likely_side_bases[most_likely_base_TwoTrackUnit(set->read[i]->unit[j])]++;

	out->side_bases++;
	for(k=0;k<4;k++) {
	  total_side_lik += set->read[i]->unit[j]->emission[k];
	}
      } else {
	warn("Not recognisable type");
      }
    }
  }

  out->main_avg_likelihood = (double) total_main_lik / (double) out->main_bases;

  out->side_avg_likelihood = (double) total_side_lik / (double) out->side_bases;

  return out;
}

# line 147 "two_track.dy"
base most_likely_base_TwoTrackUnit(TwoTrackUnit * u)
{
  int i;
  double max;
  int best;
  
  for(i=1, best = 0, max = u->emission[0];i<4;i++) {
    if( u->emission[i] > max ) {
      max = u->emission[i];
      best = i;
    }
  }

  return best;
}

/* Function:  TwoTrackScore_from_TwoTrack(t,emission_ratio)
 *
 * Descrip:    complete conversion
 *
 *
 * Arg:                     t [UNKN ] Undocumented argument [TwoTrack *]
 * Arg:        emission_ratio [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
# line 166 "two_track.dy"
TwoTrackScore * TwoTrackScore_from_TwoTrack(TwoTrack * t,Probability emission_ratio)
{
  TwoTrackScore * out;
  int i;

  assert(t != NULL);
  assert(emission_ratio != 0.0 );

  out = TwoTrackScore_alloc_len(t->len);
  
  for(i=0;i<t->len;i++) {
    add_TwoTrackScore(out,TwoTrackScoreUnit_from_TwoTrackUnit(t->unit[i],emission_ratio));
  }
  
  return out;
}


/* Function:  TwoTrackScoreUnit_from_TwoTrackUnit(u,ratio)
 *
 * Descrip:    conversion
 *
 *
 * Arg:            u [UNKN ] Undocumented argument [TwoTrackUnit *]
 * Arg:        ratio [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
# line 187 "two_track.dy"
TwoTrackScoreUnit * TwoTrackScoreUnit_from_TwoTrackUnit(TwoTrackUnit * u,Probability ratio)
{
  TwoTrackScoreUnit * out;
  int j;

  out = TwoTrackScoreUnit_alloc();
  out->type = u->type;

  for(j=0;j<4;j++) {
    out->emission[j] = Probability2Score(u->emission[j]/ratio);
  }

  return out;
}

/* Function:  write_plain_TwoTrackSet(set,ofp)
 *
 * Descrip:    plain output
 *
 *
 * Arg:        set [UNKN ] Undocumented argument [TwoTrackSet *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 205 "two_track.dy"
void write_plain_TwoTrackSet(TwoTrackSet * set,FILE * ofp)
{
  int j;

  assert(set != NULL);
  assert(ofp != NULL);

  for(j=0;j<set->len;j++) {
    write_plain_TwoTrack(set->read[j],ofp);
  }
}

/* Function:  TwoTrackUnit_type_to_string(type)
 *
 * Descrip:    type string
 *
 *
 * Arg:        type [UNKN ] Undocumented argument [TwoTrackUnit_type]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
# line 220 "two_track.dy"
char * TwoTrackUnit_type_to_string(TwoTrackUnit_type type)
{
  switch(type) {
  case TwoTrackUnit_main :
      return("Main");
    
  case TwoTrackUnit_side :
      return("Side");

  default :
    return("UnknownTwoTypeUnit_type");
  }
}


/* Function:  write_plain_TwoTrack(two,ofp)
 *
 * Descrip:    output in columns - base, and 4 likes
 *
 *
 *
 * Arg:        two [UNKN ] Undocumented argument [TwoTrack *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 239 "two_track.dy"
void write_plain_TwoTrack(TwoTrack * two,FILE * ofp)
{
  int i;


  assert(two != NULL);
  assert(ofp != NULL);

  fprintf(ofp,"Sequence %s\n",two->seq->name);
  for(i=0;i<two->seq->len;i++) {
    fprintf(ofp,"%c\t%s\t%c\t%f\t%c\t%f\t%c\t%f\t%c\t%f\n",two->seq->seq[i],TwoTrackUnit_type_to_string(two->unit[i]->type),
	    char_from_base(0),two->unit[i]->emission[0],
	    char_from_base(1),two->unit[i]->emission[1],
	    char_from_base(2),two->unit[i]->emission[2],
	    char_from_base(3),two->unit[i]->emission[3]);

  }

  fprintf(ofp,"//\n");

}


/* Function:  read_two_track_Tim_style(seq_file,like_file)
 *
 * Descrip:    Reads a two track file with a default reporting of every
 *             1000 sequences read and the entire file
 *
 *
 * Arg:         seq_file [UNKN ] Undocumented argument [char *]
 * Arg:        like_file [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
# line 266 "two_track.dy"
TwoTrackSet * read_two_track_Tim_style(char * seq_file,char * like_file)
{
  return read_two_track_Tim_style_report(seq_file,like_file,1000,0);
}


/* Function:  read_two_track_Tim_style_report(seq_file,like_file,report_freq,truncate_after)
 *
 * Descrip:    Reads a two track file from Solexa style information
 *             from Tim Massinghams lik files
 *
 *             Assummes all positions are main emissions
 *
 *             report freq means the function will report every xx lines
 *
 *
 * Arg:              seq_file [UNKN ] Undocumented argument [char *]
 * Arg:             like_file [UNKN ] Undocumented argument [char *]
 * Arg:           report_freq [UNKN ] Undocumented argument [int]
 * Arg:        truncate_after [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
# line 280 "two_track.dy"
TwoTrackSet * read_two_track_Tim_style_report(char * seq_file,char * like_file,int report_freq, int truncate_after)
{
  FILE * sfp;
  FILE * lfp;

  TwoTrackSet * out;
  TwoTrack * t;
  TwoTrackUnit * u;

  char sbuffer[MAXLINE];
  char lbuffer[40000];

  long int totallen  = 0;
  long int count = 0;

  char nameb[MAXLINE];

  char * strp;
  char * tok;
  char * split[40];

  int i;
  int j;

  double val;

  assert(seq_file != NULL);
  assert(like_file != NULL);

  sfp = openfile(seq_file,"r");
  lfp = openfile(like_file,"r");

  out = TwoTrackSet_alloc_std();

  while( fgets(sbuffer,MAXLINE,sfp) != NULL ) {
    if( truncate_after != 0 && count > truncate_after ) {
      break;
    }

    if( fgets(lbuffer,40000,lfp) == NULL ) {
      warn("Bad scenario got a sequence line, no likelihood, exiting");
      break;
    }

    
    /* parse sbuffer */
    strp = sbuffer;
    for(i=0; (tok = strsep(&strp,spacestr)) != NULL && i < 5;i++) {
      split[i] = tok;
    }

    sprintf(nameb,"Solexa_%s_%s_%s_%s",split[0],split[1],split[2],split[3]);



    t = TwoTrack_alloc_std();
    t->seq = new_Sequence_from_strings(nameb,split[4]);

    totallen += t->seq->len;
    count++;

    if( report_freq != 0 && count % report_freq == 0 ) {
      info("Read %d sequences %s, average lenght %.2f",count,t->seq->name,(double) totallen / count);
    }
    
   

    /* now parse like */
    strp = lbuffer;
    for(i=0; i < 4 && (tok = strsep(&strp,spacestr)) != NULL;i++) {
       /* could check it matched */
    }

    for(i=0;i<t->seq->len;i++) {
      u = TwoTrackUnit_alloc();
      add_TwoTrack(t,u);

      for(j = 0;j< 4;j++) {
	tok = strsep(&strp,spacestr);
	if( tok == NULL ) {
	  warn("Ran out likelihood values on sequence %s, position %d, nucleotide %d. Going to exit ugly",nameb,i,j);
	  /* this is not ideal */
	  return NULL;
	}

	val = strtod(tok,NULL);

	u->emission[base_from_char(two_track_alphabet[j])] = val;
	u->type = TwoTrackUnit_main;
      }
    }

    add_TwoTrackSet(out,t);
  }


  return out;
}


# line 388 "two_track.c"
/* Function:  hard_link_TwoTrackUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackUnit *]
 *
 */
TwoTrackUnit * hard_link_TwoTrackUnit(TwoTrackUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TwoTrackUnit object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TwoTrackUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackUnit *]
 *
 */
TwoTrackUnit * TwoTrackUnit_alloc(void) 
{
    TwoTrackUnit * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TwoTrackUnit *) ckalloc (sizeof(TwoTrackUnit))) == NULL)    {  
      warn("TwoTrackUnit_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    /* emission[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_TwoTrackUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackUnit *]
 *
 */
TwoTrackUnit * free_TwoTrackUnit(TwoTrackUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TwoTrackUnit obj. Should be trappable");  
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


/* Function:  swap_TwoTrack(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TwoTrack
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TwoTrackUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TwoTrack(TwoTrackUnit ** list,int i,int j)  
{
    TwoTrackUnit * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TwoTrack(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TwoTrack which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TwoTrackUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TwoTrack(TwoTrackUnit ** list,int left,int right,int (*comp)(TwoTrackUnit * ,TwoTrackUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TwoTrack(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TwoTrack (list,++last,i);   
      }  
    swap_TwoTrack (list,left,last);  
    qsort_TwoTrack(list,left,last-1,comp);   
    qsort_TwoTrack(list,last+1,right,comp);  
}    


/* Function:  sort_TwoTrack(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TwoTrack
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TwoTrack *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TwoTrack(TwoTrack * obj,int (*comp)(TwoTrackUnit *, TwoTrackUnit *)) 
{
    qsort_TwoTrack(obj->unit,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_TwoTrack(obj,len)
 *
 * Descrip:    Really an internal function for add_TwoTrack
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrack *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TwoTrack(TwoTrack * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TwoTrack called with no need");   
      return TRUE;   
      }  


    if( (obj->unit = (TwoTrackUnit ** ) ckrealloc (obj->unit,sizeof(TwoTrackUnit *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_TwoTrack, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TwoTrack(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrack *]
 * Arg:        add [OWNER] Object to add to the list [TwoTrackUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TwoTrack(TwoTrack * obj,TwoTrackUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TwoTrack(obj,obj->len + TwoTrackLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_TwoTrack(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TwoTrack(TwoTrack * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_TwoTrackUnit(obj->unit[i]); 
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TwoTrack_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoTrack_alloc_len(TwoTrackLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * TwoTrack_alloc_std(void) 
{
    return TwoTrack_alloc_len(TwoTrackLISTLENGTH);   
}    


/* Function:  TwoTrack_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * TwoTrack_alloc_len(int len) 
{
    TwoTrack * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TwoTrack_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (TwoTrackUnit ** ) ckcalloc (len,sizeof(TwoTrackUnit *))) == NULL)   {  
      warn("Warning, ckcalloc failed in TwoTrack_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TwoTrack(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * hard_link_TwoTrack(TwoTrack * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TwoTrack object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TwoTrack_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * TwoTrack_alloc(void) 
{
    TwoTrack * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TwoTrack *) ckalloc (sizeof(TwoTrack))) == NULL)    {  
      warn("TwoTrack_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->seq = NULL; 
    out->unit = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_TwoTrack(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrack *]
 *
 */
TwoTrack * free_TwoTrack(TwoTrack * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TwoTrack obj. Should be trappable");  
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
    if( obj->seq != NULL)    
      free_Sequence(obj->seq);   
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_TwoTrackUnit(obj->unit[i]);   
        }  
      ckfree(obj->unit); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TwoTrackSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TwoTrackSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TwoTrack **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TwoTrackSet(TwoTrack ** list,int i,int j)  
{
    TwoTrack * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TwoTrackSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TwoTrackSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TwoTrack **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TwoTrackSet(TwoTrack ** list,int left,int right,int (*comp)(TwoTrack * ,TwoTrack * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TwoTrackSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TwoTrackSet (list,++last,i);    
      }  
    swap_TwoTrackSet (list,left,last);   
    qsort_TwoTrackSet(list,left,last-1,comp);    
    qsort_TwoTrackSet(list,last+1,right,comp);   
}    


/* Function:  sort_TwoTrackSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TwoTrackSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TwoTrackSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TwoTrackSet(TwoTrackSet * obj,int (*comp)(TwoTrack *, TwoTrack *)) 
{
    qsort_TwoTrackSet(obj->read,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_TwoTrackSet(obj,len)
 *
 * Descrip:    Really an internal function for add_TwoTrackSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrackSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TwoTrackSet(TwoTrackSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TwoTrackSet called with no need");    
      return TRUE;   
      }  


    if( (obj->read = (TwoTrack ** ) ckrealloc (obj->read,sizeof(TwoTrack *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_TwoTrackSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TwoTrackSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrackSet *]
 * Arg:        add [OWNER] Object to add to the list [TwoTrack *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TwoTrackSet(TwoTrackSet * obj,TwoTrack * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TwoTrackSet(obj,obj->len + TwoTrackSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->read[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_TwoTrackSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TwoTrackSet(TwoTrackSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->read[i] != NULL)  {  
        free_TwoTrack(obj->read[i]); 
        obj->read[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TwoTrackSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoTrackSet_alloc_len(TwoTrackSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * TwoTrackSet_alloc_std(void) 
{
    return TwoTrackSet_alloc_len(TwoTrackSetLISTLENGTH); 
}    


/* Function:  TwoTrackSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * TwoTrackSet_alloc_len(int len) 
{
    TwoTrackSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TwoTrackSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->read = (TwoTrack ** ) ckcalloc (len,sizeof(TwoTrack *))) == NULL)   {  
      warn("Warning, ckcalloc failed in TwoTrackSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TwoTrackSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * hard_link_TwoTrackSet(TwoTrackSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TwoTrackSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TwoTrackSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * TwoTrackSet_alloc(void) 
{
    TwoTrackSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TwoTrackSet *) ckalloc (sizeof(TwoTrackSet))) == NULL)  {  
      warn("TwoTrackSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->read = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_TwoTrackSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackSet *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSet *]
 *
 */
TwoTrackSet * free_TwoTrackSet(TwoTrackSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TwoTrackSet obj. Should be trappable");   
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
    if( obj->read != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->read[i] != NULL)    
          free_TwoTrack(obj->read[i]);   
        }  
      ckfree(obj->read); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_TwoTrackScoreUnit(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * hard_link_TwoTrackScoreUnit(TwoTrackScoreUnit * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TwoTrackScoreUnit object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TwoTrackScoreUnit_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * TwoTrackScoreUnit_alloc(void) 
{
    TwoTrackScoreUnit * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TwoTrackScoreUnit *) ckalloc (sizeof(TwoTrackScoreUnit))) == NULL)  {  
      warn("TwoTrackScoreUnit_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    /* emission[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_TwoTrackScoreUnit(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScoreUnit *]
 *
 */
TwoTrackScoreUnit * free_TwoTrackScoreUnit(TwoTrackScoreUnit * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TwoTrackScoreUnit obj. Should be trappable"); 
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


/* Function:  swap_TwoTrackScore(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TwoTrackScore
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TwoTrackScoreUnit **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TwoTrackScore(TwoTrackScoreUnit ** list,int i,int j)  
{
    TwoTrackScoreUnit * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TwoTrackScore(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TwoTrackScore which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TwoTrackScoreUnit **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TwoTrackScore(TwoTrackScoreUnit ** list,int left,int right,int (*comp)(TwoTrackScoreUnit * ,TwoTrackScoreUnit * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TwoTrackScore(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TwoTrackScore (list,++last,i);  
      }  
    swap_TwoTrackScore (list,left,last); 
    qsort_TwoTrackScore(list,left,last-1,comp);  
    qsort_TwoTrackScore(list,last+1,right,comp); 
}    


/* Function:  sort_TwoTrackScore(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TwoTrackScore
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TwoTrackScore *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TwoTrackScore(TwoTrackScore * obj,int (*comp)(TwoTrackScoreUnit *, TwoTrackScoreUnit *)) 
{
    qsort_TwoTrackScore(obj->unit,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_TwoTrackScore(obj,len)
 *
 * Descrip:    Really an internal function for add_TwoTrackScore
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrackScore *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TwoTrackScore(TwoTrackScore * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TwoTrackScore called with no need");  
      return TRUE;   
      }  


    if( (obj->unit = (TwoTrackScoreUnit ** ) ckrealloc (obj->unit,sizeof(TwoTrackScoreUnit *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TwoTrackScore, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TwoTrackScore(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TwoTrackScore *]
 * Arg:        add [OWNER] Object to add to the list [TwoTrackScoreUnit *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TwoTrackScore(TwoTrackScore * obj,TwoTrackScoreUnit * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TwoTrackScore(obj,obj->len + TwoTrackScoreLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->unit[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_TwoTrackScore(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TwoTrackScore *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TwoTrackScore(TwoTrackScore * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->unit[i] != NULL)  {  
        free_TwoTrackScoreUnit(obj->unit[i]);    
        obj->unit[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TwoTrackScore_alloc_std(void)
 *
 * Descrip:    Equivalent to TwoTrackScore_alloc_len(TwoTrackScoreLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * TwoTrackScore_alloc_std(void) 
{
    return TwoTrackScore_alloc_len(TwoTrackScoreLISTLENGTH); 
}    


/* Function:  TwoTrackScore_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * TwoTrackScore_alloc_len(int len) 
{
    TwoTrackScore * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TwoTrackScore_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->unit = (TwoTrackScoreUnit ** ) ckcalloc (len,sizeof(TwoTrackScoreUnit *))) == NULL) {  
      warn("Warning, ckcalloc failed in TwoTrackScore_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TwoTrackScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackScore *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * hard_link_TwoTrackScore(TwoTrackScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TwoTrackScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TwoTrackScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * TwoTrackScore_alloc(void) 
{
    TwoTrackScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TwoTrackScore *) ckalloc (sizeof(TwoTrackScore))) == NULL)  {  
      warn("TwoTrackScore_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->unit = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_TwoTrackScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackScore *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackScore *]
 *
 */
TwoTrackScore * free_TwoTrackScore(TwoTrackScore * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TwoTrackScore obj. Should be trappable"); 
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
    if( obj->unit != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->unit[i] != NULL)    
          free_TwoTrackScoreUnit(obj->unit[i]);  
        }  
      ckfree(obj->unit); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_TwoTrackSetStats(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TwoTrackSetStats *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * hard_link_TwoTrackSetStats(TwoTrackSetStats * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TwoTrackSetStats object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TwoTrackSetStats_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * TwoTrackSetStats_alloc(void) 
{
    TwoTrackSetStats * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TwoTrackSetStats *) ckalloc (sizeof(TwoTrackSetStats))) == NULL)    {  
      warn("TwoTrackSetStats_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->reads = 0;  
    out->main_bases = 0; 
    out->side_bases = 0; 
    out->main_avg_likelihood = 0;    
    out->side_avg_likelihood = 0;    
    /* called_main_bases[5] is an array: no default possible */ 
    /* most_likely_main_bases[5] is an array: no default possible */ 
    /* most_likely_side_bases[5] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_TwoTrackSetStats(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TwoTrackSetStats *]
 *
 * Return [UNKN ]  Undocumented return value [TwoTrackSetStats *]
 *
 */
TwoTrackSetStats * free_TwoTrackSetStats(TwoTrackSetStats * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TwoTrackSetStats obj. Should be trappable");  
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



#ifdef _cplusplus
}
#endif
