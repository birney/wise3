#ifdef _cplusplus
extern "C" {
#endif
#include "signalalign.h"


# line 108 "signalalign.dy"
SignalEventList * implied_event_from_RawSignalSeq_align(Sequence * seq,RawSignalSeq * rss,SignalMap * sm,AlnBlock * alb)
{

  AlnColumn * alc;
  int i = 1;

  char kbuf[64];

  double sigtime_off;
  int signal_start;
  int signal_end;
  long int rawscore_acc;
  double totalsignal;
  int rogue_number;
  int sequence_count;

  double sigtimestart;
  double sigtimelen;
  int t;
  double mean_sig;
  double sig_var;
  double diff;

  int seq_start;
  int seq_end;

  int kmer;
  double pred;

  SignalEventList * out;
  SignalEvent * event;

  
  out = SignalEventList_alloc_std();

  out->name = stringalloc(rss->name);

  assert(seq  != NULL);
  assert(rss  != NULL);
  assert(alb  != NULL);
  assert(sm   != NULL);


  for(alc=alb->start;alc != NULL;alc = alc->next) {
    
    
    if( strcmp(alc->alu[0]->text_label,"END") == 0 ) {
      break;
    }

    i++;
    
    /* this should be a SEQUENCE MATCH - move forward until it is */

    for(;alc != NULL && strcmp(alc->alu[0]->text_label,"SEQUENCE") != 0;alc = alc->next)
      ;

    if( alc == NULL ) {
      break;
    }

    seq_start = alc->alu[0]->end;

    /* work out starting signal time point */
    sigtime_off = 0.0;
    for(t = 0;t < alc->alu[1]->start;t++) {
      sigtime_off += rss->time_len[t];
    }

    signal_start = alc->alu[1]->start;
    sigtimestart = rss->start_time + sigtime_off;

    rawscore_acc = 0;
    totalsignal = 0.0;
    rogue_number = 0;
    sequence_count = 0;

    sigtimelen = 0.0;

    /* now move along, eating up SEQUENCE and ROGUE_SEQUENCE cases */

    
    for(;alc != NULL && (strcmp(alc->alu[0]->text_label,"SEQUENCE") == 0 || strcmp(alc->alu[0]->text_label,"ROGUE_SEQUENCE") == 0);alc = alc->next) {
      rawscore_acc += alc->alu[0]->score[0];
      if( strcmp(alc->alu[0]->text_label,"ROGUE_SEQUENCE") == 0 ) {
	rogue_number++;
      }
      /* we have to step through the span */
      for(t=alc->alu[1]->start;t < alc->alu[1]->end;t++) {
	  sigtimelen += rss->time_len[t];
      }

      totalsignal += rss->signal[alc->alu[1]->end];
      sequence_count++;
    }
      
    seq_end = alc->alu[0]->end;
    
    /* this is now at the end of this sequence */

    signal_end = alc->alu[1]->end;

    /* calculate mean */

    mean_sig = totalsignal / sequence_count;

    /* now go back to calculate variance */

    sig_var = 0.0;
    
    for(t = signal_start;t < signal_end;t++) {
      diff = (rss->signal[t] - mean_sig);

      sig_var += (diff * diff);
    }

    sig_var = sig_var / sequence_count;
    

    if( alc->alu[0]->start > 5 && alc->alu[0]->end +1 < seq->len) {
      strncpy(kbuf,seq->seq+alc->alu[0]->start-4,5);
      kmer = forward_dna_number_from_string(seq->seq+alc->alu[0]->end - sm->kbasis,sm->kbasis);
      kbuf[5] = '\0';
    } else {
      kmer = 1;
      strncpy(kbuf,"??????",5);
    }
	    
    event = SignalEvent_alloc();

    event->base = kbuf[4];
    strncpy(event->kmer,kbuf,5);
    event->mean = mean_sig;
    event->std = sqrt(sig_var);
    event->time_pos = sigtimestart;
    event->time_length = sigtimelen;

    add_SignalEventList(out,event);
  }


  return(out);
}


# line 253 "signalalign.dy"
void show_fit_RawSignalMat(Sequence * seq,RawSignalSeq * rss,RawSignalMatParaScore * rsmp,SignalMap * sm,AlnBlock * alb,FILE * ofp)
{
  AlnColumn * alc;
  int i = 1;

  char kbuf[64];

  double sigtime_off;
  int signal_start;
  int signal_end;
  long int rawscore_acc;
  double totalsignal;
  int rogue_number;
  int sequence_count;

  double sigtimestart;
  double sigtimelen;
  int t;
  double mean_sig;
  double sig_var;
  double diff;

  int seq_start;
  int seq_end;

  int kmer;
  double pred;

  assert(seq  != NULL);
  assert(rss  != NULL);
  assert(rsmp != NULL);
  assert(alb  != NULL);
  assert(ofp  != NULL);

  
  fprintf(ofp,"AlnNum\tSigLabel\tSeqLabel\tRawScore\tBitsScore\tSigStart\tSigEnd\tSigMean\tSigStdev\tSigBase\tSigKmer\tSigTimeStart\tSigTimeLen\tSignalFit\tSeqStart\tSeqEnd\tSeqBase\tSeqKmer\tModelMean\tRogueSignal\n");
  
  for(alc=alb->start;alc != NULL;alc = alc->next) {
    
    
    if( strcmp(alc->alu[0]->text_label,"END") == 0 ) {
      break;
    }

    fprintf(ofp,"%d\t",i);
    i++;
    
    /* this should be a SEQUENCE MATCH - move forward until it is */

    for(;alc != NULL && strcmp(alc->alu[0]->text_label,"SEQUENCE") != 0;alc = alc->next)
      ;

    if( alc == NULL ) {
      break;
    }

    seq_start = alc->alu[0]->end;

    /* work out starting signal time point */
    sigtime_off = 0.0;
    for(t = 0;t < alc->alu[1]->start;t++) {
      sigtime_off += rss->time_len[t];
    }

    signal_start = alc->alu[1]->start;
    sigtimestart = rss->start_time + sigtime_off;

    rawscore_acc = 0;
    totalsignal = 0.0;
    rogue_number = 0;
    sequence_count = 0;

    sigtimelen = 0.0;

    /* now move along, eating up SEQUENCE and ROGUE_SEQUENCE cases */

    
    for(;alc != NULL && (strcmp(alc->alu[0]->text_label,"SEQUENCE") == 0 || strcmp(alc->alu[0]->text_label,"ROGUE_SEQUENCE") == 0);alc = alc->next) {
      rawscore_acc += alc->alu[0]->score[0];
      if( strcmp(alc->alu[0]->text_label,"ROGUE_SEQUENCE") == 0 ) {
	rogue_number++;
      }
      /* we have to step through the span */
      for(t=alc->alu[1]->start;t < alc->alu[1]->end;t++) {
	  sigtimelen += rss->time_len[t];
      }

      totalsignal += rss->signal[alc->alu[1]->end];
      sequence_count++;
    }
      
    seq_end = alc->alu[0]->end;
    
    /* this is now at the end of this sequence */

    signal_end = alc->alu[1]->end;

    /* calculate mean */

    mean_sig = totalsignal / sequence_count;

    /* now go back to calculate variance */

    sig_var = 0.0;
    
    for(t = signal_start;t < signal_end;t++) {
      diff = (rss->signal[t] - mean_sig);

      sig_var += (diff * diff);
    }

    sig_var = sig_var / sequence_count;
    

    if( alc->alu[0]->start > 5 && alc->alu[0]->end +1 < seq->len) {
      strncpy(kbuf,seq->seq+alc->alu[0]->start-4,5);
      kmer = forward_dna_number_from_string(seq->seq+alc->alu[0]->end - sm->kbasis,sm->kbasis);
      kbuf[5] = '\0';
    } else {
      kmer = 1;
      strncpy(kbuf,"??????",5);
    }
	    


    if( kmer >= 1 && kmer < sm->len ) 
      pred = sm->comp[kmer]->mean;


    
    /* ok, now ready to output ! */

    fprintf(ofp,"SIGNAL\tSEQUENCE\t%ld\t%f\t%d\t%d\t%f\t%f\t?\t?\t%f\t%f\t0.0\t%d\t%d\t%c\t%s\t%f\t%d\n",
	    rawscore_acc,
	    Score2Bits(rawscore_acc),
	    signal_start,
	    signal_end,
	    mean_sig,
	    sqrt(sig_var),
	    sigtimestart,
	    sigtimelen,
	    seq_start,
	    seq_end,
	    kbuf[4],
	    kbuf,
	    pred,
	    rogue_number);
	    
	    

  }

  fprintf(ofp,"//\n");

}
# line 408 "signalalign.dy"
void show_help_RawSignalMatParaProb_from_argv(FILE * ofp)
{
  assert(ofp != NULL);

  fprintf(ofp,"Raw Signal alignment parameters\n");
  fprintf(ofp,"  -raw_single [0.05] - Probability of a single signal kmer\n");
  fprintf(ofp,"  -raw_double [0.15] - Probability of a double signal kmer\n");
  fprintf(ofp,"  -raw_triple [0.4]  - Probability of a triple signal kmer\n");
  fprintf(ofp,"  -raw_quad   [0.4]  - Probability of a quad or more signal kmer\n");
  fprintf(ofp,"  -raw_event_ext   [0.5]   - Extension beyond 4 signal kmer\n");
  fprintf(ofp,"  -raw_between     [0.05]  - Probability of a between kmer signal\n");
  fprintf(ofp,"  -raw_between_ext [0.2]   - Extension of inter-kmer signal\n");
  fprintf(ofp,"  -raw_rogue       [0.01]  - Probability of rogue signal inside a kmer\n");
  fprintf(ofp,"  -raw_rogue_ext [0.5]  - Extension beyond 4 signal kmer\n");




}

# line 428 "signalalign.dy"
RawSignalMatParaScore * make_RawSignalMatParaScore(RawSignalMatParaProb * rsmp)
{
  RawSignalMatParaScore * out;

  assert(rsmp != NULL);

  out = RawSignalMatParaScore_alloc();

  out->single_event_signal  = Probability2Score(rsmp->single_event_signal);
  out->double_event_signal  = Probability2Score(rsmp->double_event_signal);
  out->triple_event_signal  = Probability2Score(rsmp->triple_event_signal);
  out->quad_event_signal    = Probability2Score(rsmp->quad_event_signal);
  out->event_ext            = Probability2Score(rsmp->event_ext);
  out->between_event_signal = Probability2Score(rsmp->between_event_signal);
  out->between_event_ext    = Probability2Score(rsmp->between_event_ext);
  out->rogue_signal         = Probability2Score(rsmp->rogue_signal);
  out->rogue_ext            = Probability2Score(rsmp->rogue_ext);
  out->dark_kmer            = Probability2Score(rsmp->dark_kmer);

  return(out);

}

# line 451 "signalalign.dy"
RawSignalMatParaProb * new_RawSignalMatParaProb_from_argv(int * argc,char ** argv)
{
  RawSignalMatParaProb * out;

  out = RawSignalMatParaProb_alloc();

  out->single_event_signal = 0.05;
  out->double_event_signal = 0.15;
  out->triple_event_signal = 0.4;
  out->quad_event_signal   = 0.4;

  out->event_ext = 0.5;

  out->between_event_signal = 0.05;
  out->between_event_ext = 0.2;

  out->rogue_signal = 0.01;
  out->rogue_ext    = 0.05;

  out->dark_kmer    = 0.01;


  strip_out_float_argument(argc,argv,"raw_single",&out->single_event_signal);
  strip_out_float_argument(argc,argv,"raw_double",&out->double_event_signal);
  strip_out_float_argument(argc,argv,"raw_triple",&out->triple_event_signal);
  strip_out_float_argument(argc,argv,"raw_quad",&out->quad_event_signal);
  strip_out_float_argument(argc,argv,"raw_event_ext",&out->event_ext);
  strip_out_float_argument(argc,argv,"raw_between",&out->between_event_signal);
  strip_out_float_argument(argc,argv,"raw_between_ext",&out->between_event_ext);
  strip_out_float_argument(argc,argv,"raw_rogue",&out->rogue_signal);
  strip_out_float_argument(argc,argv,"raw_rogue_ext",&out->rogue_ext);
  strip_out_float_argument(argc,argv,"raw_dark",&out->dark_kmer);


  return(out);
}

# line 391 "signalalign.c"
/* Function:  hard_link_RawSignalMatParaProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RawSignalMatParaProb *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaProb *]
 *
 */
RawSignalMatParaProb * hard_link_RawSignalMatParaProb(RawSignalMatParaProb * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RawSignalMatParaProb object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RawSignalMatParaProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaProb *]
 *
 */
RawSignalMatParaProb * RawSignalMatParaProb_alloc(void) 
{
    RawSignalMatParaProb * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RawSignalMatParaProb *) ckalloc (sizeof(RawSignalMatParaProb))) == NULL)    {  
      warn("RawSignalMatParaProb_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->single_event_signal = 0.0;  
    out->double_event_signal = 0.0;  
    out->triple_event_signal = 0.0;  
    out->quad_event_signal = 0.0;    
    out->event_ext = 0.0;    
    out->between_event_signal = 0.0; 
    out->between_event_ext = 0.0;    
    out->rogue_signal = 0.0; 
    out->rogue_ext = 0.0;    
    out->dark_kmer = 0.0;    


    return out;  
}    


/* Function:  free_RawSignalMatParaProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalMatParaProb *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaProb *]
 *
 */
RawSignalMatParaProb * free_RawSignalMatParaProb(RawSignalMatParaProb * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RawSignalMatParaProb obj. Should be trappable");  
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


/* Function:  hard_link_RawSignalMatParaScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RawSignalMatParaScore *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaScore *]
 *
 */
RawSignalMatParaScore * hard_link_RawSignalMatParaScore(RawSignalMatParaScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RawSignalMatParaScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RawSignalMatParaScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaScore *]
 *
 */
RawSignalMatParaScore * RawSignalMatParaScore_alloc(void) 
{
    RawSignalMatParaScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RawSignalMatParaScore *) ckalloc (sizeof(RawSignalMatParaScore))) == NULL)  {  
      warn("RawSignalMatParaScore_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->single_event_signal = 0;    
    out->double_event_signal = 0;    
    out->triple_event_signal = 0;    
    out->quad_event_signal = 0;  
    out->event_ext = 0;  
    out->between_event_signal = 0;   
    out->between_event_ext = 0;  
    out->rogue_signal = 0;   
    out->rogue_ext = 0;  
    out->dark_kmer = 0;  


    return out;  
}    


/* Function:  free_RawSignalMatParaScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalMatParaScore *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMatParaScore *]
 *
 */
RawSignalMatParaScore * free_RawSignalMatParaScore(RawSignalMatParaScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RawSignalMatParaScore obj. Should be trappable"); 
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




  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Thu May 22 15:42:15 2014           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define MATCH 0  
#define BETWEEN 1    
#define ROGUE 2  


#define START 0  
#define END 1    


#define RawSignalMat_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+4)*3)+STATE][i+1]    
#define RawSignalMat_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+4]   
#define RawSignalMat_READ_OFF_ERROR -6
  


#define RawSignalMat_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+5)%5][((i+1)*3)+STATE]  
#define RawSignalMat_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+5)%5][STATE]   




#define RawSignalMat_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]    
#define RawSignalMat_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)    


/* Function:  PackAln_read_Shatter_RawSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_RawSignalMat(RawSignalMat * mat) 
{
    RawSignalMat_access_func_holder holder;  


    holder.access_main    = RawSignalMat_shatter_access_main;    
    holder.access_special = RawSignalMat_shatter_access_special; 
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_RawSignalMat(mat,holder);    
}    


/* Function:  RawSignalMat_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int RawSignalMat_shatter_access_main(RawSignalMat * mat,int i,int j,int state) 
{
    return RawSignalMat_SHATTER_MATRIX(mat,i,j,state);   
}    


/* Function:  RawSignalMat_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int RawSignalMat_shatter_access_special(RawSignalMat * mat,int i,int j,int state) 
{
    return RawSignalMat_SHATTER_SPECIAL(mat,i,j,state);  
}    


/* Function:  calculate_shatter_RawSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the RawSignalMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [RawSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_RawSignalMat(RawSignalMat * mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int should_calc;     
    int leni;    
    int lenj;    
    int tot; 
    int num; 
    int starti;  
    int startj;  
    int endi;    
    int endj;    


    int * SIG_0_0;   
    int * SIG_0_1;   
    int * SIG_0_2;   
    int * SIG_0_3;   
    int * SIG_0_4;   
    int * SIG_1_1;   
    int * SIG_1_0;   


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,3,lenj,2);    
    prepare_DPEnvelope(dpenv);   
    starti = dpenv->starti;  
    if( starti < 0 ) 
      starti = 0;    
    startj = dpenv->startj;  
    if( startj < 0 ) 
      startj = 0;    
    endi = dpenv->endi;  
    if( endi > mat->leni )   
      endi = mat->leni;  
    endj = dpenv->endj;  
    if( endj > mat->lenj )   
      endj = mat->lenj;  
    tot = (endi-starti) * (endj-startj); 
    num = 0; 


    start_reporting("RawSignalMat Matrix calculation: ");    
    for(j=startj;j<endj;j++) {  
      auto int score;    
      auto int temp;     
      for(i=starti;i<endi;i++)   {  
        /* Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed */ 
        should_calc = 0; 
        for(k=0;k<dpenv->len;k++)    {  
          auto DPUnit * u;   
          u = dpenv->dpu[k]; 
          switch(u->type)    {  
            case DPENV_RECT :    
              if( i >= u->starti && j >= u->startj && i < (u->starti+u->height) && j < (u->startj+u->length))    
                should_calc = 1;     
              break; 
            case DPENV_DIAG :    
              if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)   
                should_calc = 1;     
              break; 
            }  
          if( should_calc == 1 ) 
            break;   
          }  
        if( should_calc == 0)    
          continue;  


        SIG_0_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i,j);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_0_2 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-2);   
        SIG_0_3 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-3);   
        SIG_0_4 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-4);   
        SIG_1_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-1);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = SIG_0_1[BETWEEN] + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);     
        /* From state BETWEEN to state MATCH */ 
        temp = SIG_0_2[BETWEEN] + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = SIG_0_3[BETWEEN] + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = SIG_0_4[BETWEEN] + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = SIG_0_1[MATCH] + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = SIG_0_1[ROGUE] + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = RawSignalMat_SHATTER_SPECIAL(mat,i-1,j-1,START) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > RawSignalMat_SHATTER_SPECIAL(mat,i,j,END) )   {  
          RawSignalMat_SHATTER_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH] + 0;  
        /* From state BETWEEN to state BETWEEN */ 
        temp = SIG_0_1[BETWEEN] + mat->para->between_event_signal;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[BETWEEN] = score;   


        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = SIG_0_1[MATCH] + mat->para->rogue_signal;    
        /* From state ROGUE to state ROGUE */ 
        temp = SIG_0_1[ROGUE] + mat->para->rogue_ext;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[ROGUE] = score; 


        /* Finished calculating state ROGUE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_RawSignalMat(dbsi,out,seq,signal,sm,para)
 *
 * Descrip:    This function makes a database search of RawSignalMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           seq [UNKN ] Undocumented argument [Sequence*]
 * Arg:        signal [UNKN ] Undocumented argument [RawSignalSeq*]
 * Arg:            sm [UNKN ] Undocumented argument [SignalMap*]
 * Arg:          para [UNKN ] Undocumented argument [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_RawSignalMat(DBSearchImpl * dbsi,Hscore * out,Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para) 
{
#ifdef PTHREAD   
    int i;   
    int thr_no;  
    pthread_attr_t pat;  
    struct thread_pool_holder_RawSignalMat * holder;     
#endif   
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_RawSignalMat. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_RawSignalMat. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for RawSignalMat, but it was compiled with a trace level of 1783633488. Not all trace statements can be shown",dbsi->trace_level);   
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_RawSignalMat(out,seq,signal ,sm,para);  
      case DBSearchImpl_Pthreads :   
#ifdef PTHREAD   
        holder = (struct thread_pool_holder_RawSignalMat *) ckalloc(sizeof(struct thread_pool_holder_RawSignalMat)); 
        if( holder == NULL )     {  
          warn("Unable to allocated thread pool datastructure...");  
          return SEARCH_ERROR;   
          }  
        holder->out = out;   
        holder->dbsi = dbsi; 
        holder->seq = seq;   
        holder->signal = signal; 
        holder->sm = sm; 
        holder->para = para; 
        if( pthread_mutex_init(&(holder->input_lock),NULL) != 0 )    
        fatal("Unable to iniated input mutex lock"); 
        if( pthread_mutex_init(&(holder->output_lock),NULL) != 0 )   
        fatal("Unable to iniated output mutex lock");    
        /* Let us rock! */ 
        thr_no = number_of_threads_DBSearchImpl(dbsi);   
        holder->pool = ckcalloc (thr_no,sizeof(pthread_t));  
        if( holder->pool == NULL )   {  
          warn("Unable to allocated thread pools");  
          return SEARCH_ERROR;   
          }  
        /* Build a thread attribute to make sure we get the most out of SMP boxes */ 
        pthread_attr_init(&pat);     
        /* Give thread libraries a hint that threads should be kernel threads */ 
#ifndef __sgi /* SGI can't set system scope ... */   
#ifdef  HAS_PTHREAD_SETSCOPE 
        pthread_attr_setscope(&pat, PTHREAD_SCOPE_SYSTEM);   
#endif /* set scope */   
#endif /* sgi */ 
        /* Give thread libraries a hint that there are num of threads to run */ 
#ifdef HAS_PTHREAD_SETCONCURRENCY    
        pthread_setconcurrency(thr_no+1);    
#endif /* set concurrency */ 
        for(i=0;i<thr_no;i++)    {  
          if( pthread_create(holder->pool+i,&pat,thread_loop_RawSignalMat,(void *)holder) )  
            fatal("Unable to create a thread!"); 
          }  
        /* Now - wait for all the threads to exit */ 
        for(i=0;i<thr_no;i++)    {  
          if( pthread_join(holder->pool[i],NULL) != 0 )  
            fatal("Unable to join a thread!");   
          }  
        /* Deallocate the thread structures */ 
        ckfree(holder->pool);    
        ckfree(holder);  
        return SEARCH_OK;    
#else /* not compiled with threads */    
        warn("You did not specifiy the PTHREAD compile when compiled the C code for RawSignalMat");  
#endif /* finished threads */    
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from RawSignalMat",impl_string_DBSearchImpl(dbsi));   
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  thread_loop_RawSignalMat(ptr)
 *
 * Descrip:    dummy loop code foreach thread for RawSignalMat
 *
 *
 * Arg:        ptr [UNKN ] Undocumented argument [void *]
 *
 * Return [UNKN ]  Undocumented return value [void *]
 *
 */
void * thread_loop_RawSignalMat(void * ptr) 
{
    fatal("dummy thread loop function"); 
}    


/* Function:  serial_search_RawSignalMat(out,seq,signal,sm,para)
 *
 * Descrip:    This function makes a database search of RawSignalMat
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           seq [UNKN ] Undocumented argument [Sequence*]
 * Arg:        signal [UNKN ] Undocumented argument [RawSignalSeq*]
 * Arg:            sm [UNKN ] Undocumented argument [SignalMap*]
 * Arg:          para [UNKN ] Undocumented argument [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_RawSignalMat(Hscore * out,Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_RawSignalMat(seq, signal , sm, para);     
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("RawSignalMat search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");   
        return SEARCH_ERROR; 
        }  
      /* Now: add query/target information to the entry */ 
      ds->score = score;     
      add_Hscore(out,ds);    
      } /* end of if storing datascore */ 
    pop_errormsg_stack();    
    push_errormsg_stack("DB searching: just finished [Query Pos: %d] [Target Pos: %d]",query_pos,target_pos);    


    pop_errormsg_stack();    
    return SEARCH_OK;    
}    


/* Function:  score_only_RawSignalMat(seq,signal,sm,para)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_RawSignalMat_only
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    RawSignalMat * mat;  


    mat = allocate_RawSignalMat_only(seq, signal , sm, para);    
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(5,(mat->leni + 1) * 3,5,2)) == NULL)  {  
      warn("Score only matrix for RawSignalMat cannot be allocated, (asking for 4  by %d  cells)",mat->leni*3);  
      mat = free_RawSignalMat(mat);  
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<6;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          RawSignalMat_VSMALL_MATRIX(mat,i,j,k) = NEGI;  
        }  
      RawSignalMat_VSMALL_SPECIAL(mat,i,j,START) = 0;    
      RawSignalMat_VSMALL_SPECIAL(mat,i,j,END) = NEGI;   
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = RawSignalMat_VSMALL_SPECIAL(mat,i-1,j-1,START) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_VSMALL_MATRIX(mat,i,j,MATCH) = score;  


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > RawSignalMat_VSMALL_SPECIAL(mat,i,j,END) )    {  
          RawSignalMat_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + 0;   
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_VSMALL_MATRIX(mat,i,j,BETWEEN) = score;    


        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;     
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_VSMALL_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_VSMALL_MATRIX(mat,i,j,ROGUE) = score;  


        /* Finished calculating state ROGUE */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      if( bestscore < RawSignalMat_VSMALL_SPECIAL(mat,0,j,END) ) 
        bestscore = RawSignalMat_VSMALL_SPECIAL(mat,0,j,END);    
      } /* end of for all target positions */ 


    mat = free_RawSignalMat(mat);    
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_RawSignalMat(seq,signal,sm,para,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_RawSignalMat
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    RawSignalMat * mat;  
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = seq->len * signal->len;  
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 3 * sizeof(int)) > 1000*dpri->kbyte_size) {  
        strategy = DPIM_Linear;  
        }  
      else   {  
        strategy = DPIM_Explicit;    
        }  
      }  
    else {  
      strategy = dpri->memory;   
      }  


    if( dpenv != NULL )  {  
      if( strategy == DPIM_Explicit) {  
        if( (mat=allocate_Expl_RawSignalMat(seq, signal , sm, para,dpri)) == NULL )  {  
          warn("Unable to allocate large RawSignalMat version"); 
          return NULL;   
          }  
        calculate_dpenv_RawSignalMat(mat,dpenv);     
        out =  PackAln_read_Expl_RawSignalMat(mat);  
        }  
      else   {  
        mat = allocate_RawSignalMat_only(seq, signal , sm, para);    
        calculate_shatter_RawSignalMat(mat,dpenv);   
        out = PackAln_read_Shatter_RawSignalMat(mat);    
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_RawSignalMat(seq, signal , sm, para)) == NULL )  {  
          warn("Unable to allocate small RawSignalMat version"); 
          return NULL;   
          }  
        out = PackAln_calculate_Small_RawSignalMat(mat,dpenv);   
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_RawSignalMat(seq, signal , sm, para,dpri)) == NULL )  {  
          warn("Unable to allocate large RawSignalMat version"); 
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_RawSignalMat(mat);   
          out =  PackAln_read_Expl_RawSignalMat(mat);    
          }  
        }  
      }  


    mat = free_RawSignalMat(mat);    
    return out;  
}    


/* Function:  allocate_RawSignalMat_only(seq,signal,sm,para)
 *
 * Descrip:    This function only allocates the RawSignalMat structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * allocate_RawSignalMat_only(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para) 
{
    RawSignalMat * out;  


    if((out= RawSignalMat_alloc()) == NULL)  {  
      warn("Allocation of basic RawSignalMat structure failed...");  
      return NULL;   
      }  


    out->seq = seq;  
    out->signal = signal;    
    out->sm = sm;    
    out->para = para;    
    out->leni = seq->len;    
    out->lenj = signal->len;     
    return out;  
}    


/* Function:  allocate_Expl_RawSignalMat(seq,signal,sm,para,dpri)
 *
 * Descrip:    This function allocates the RawSignalMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_RawSignalMat_only
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * allocate_Expl_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para,DPRunImpl * dpri) 
{
    RawSignalMat * out;  


    out = allocate_RawSignalMat_only(seq, signal , sm, para);    
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+4)*3 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+4)*3,(out->leni+1),2,out->lenj+4)) == NULL) {  
        warn("Explicit matrix RawSignalMat cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);  
        free_RawSignalMat(out);  
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_RawSignalMat(out);  
    return out;  
}    


/* Function:  init_RawSignalMat(mat)
 *
 * Descrip:    This function initates RawSignalMat matrix when in explicit mode
 *             Called in /allocate_Expl_RawSignalMat
 *
 *
 * Arg:        mat [UNKN ] RawSignalMat which contains explicit basematrix memory [RawSignalMat *]
 *
 */
void init_RawSignalMat(RawSignalMat * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->seq->len;i++) {  
      for(j= (-4);j<5;j++)   {  
        RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN) = NEGI;    
        RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE) = NEGI;  
        }  
      }  
    for(j= (-4);j<mat->signal->len;j++)  {  
      for(i= (-1);i<2;i++)   {  
        RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN) = NEGI;    
        RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE) = NEGI;  
        }  
      RawSignalMat_EXPL_SPECIAL(mat,i,j,START) = 0;  
      RawSignalMat_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  
    return;  
}    


/* Function:  recalculate_PackAln_RawSignalMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by RawSignalMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 */
void recalculate_PackAln_RawSignalMat(PackAln * pal,RawSignalMat * mat) 
{
    int i,j,k,offi,offj; 
    PackAlnUnit * prev;  
    PackAlnUnit * pau;   


    for(k=1,prev=pal->pau[0];k < pal->len;k++,prev=pau)  {  
      pau = pal->pau[k]; 
      i = pau->i;    
      j = pau->j;    
      offi = pau->i - prev->i;   
      offj = pau->j - prev->j;   
      switch(pau->state) {  
        case MATCH :     
          if( offi == 0 && offj == 1 && prev->state == BETWEEN ) {  
            pau->score = (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 2 && prev->state == BETWEEN ) {  
            pau->score = ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 3 && prev->state == BETWEEN ) {  
            pau->score = (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 4 && prev->state == BETWEEN ) {  
            pau->score = ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal) + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext) + (0);   
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == ROGUE )   {  
            pau->score = Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i) + (0);  
            continue;    
            }  
          if( offj == 1 && prev->state == (START+3) )    {  
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case BETWEEN :   
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = 0 + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == BETWEEN ) {  
            pau->score = mat->para->between_event_signal + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state BETWEEN, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case ROGUE :     
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->para->rogue_signal + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == ROGUE )   {  
            pau->score = mat->para->rogue_ext + (0);     
            continue;    
            }  
          warn("In recaluclating PackAln with state ROGUE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case (START+3) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+3) :   
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = 0 + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state END, got a bad source state. Error!");   
          break; 
        default :    
          warn("In recaluclating PackAln got a bad recipient state. Error!");    
        }  
      prev = pau;    
      }  
    return;  
}    
/* divide and conquor macros are next */ 
#define RawSignalMat_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+4)][(i+1)*3+state])  
#define RawSignalMat_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+5)*8) % 40][(i+1)*3+state])  
#define RawSignalMat_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+4)]) 
#define RawSignalMat_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+4)])    
#define RawSignalMat_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+5)*8)+(shadow+1)) % 40)][(i+1)*3 + state])   
#define RawSignalMat_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+4)])    
#define RawSignalMat_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+4)% 4) * (leni+1) * 3) + ((i+1) * 3) + (state)])    
#define RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+4)% 4) * (leni+1) * 24) + ((i+1) * 24) + (state * 8) + shadow+1])    
#define RawSignalMat_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+4)])    
/* Function:  allocate_Small_RawSignalMat(seq,signal,sm,para)
 *
 * Descrip:    This function allocates the RawSignalMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_RawSignalMat_only
 *
 *
 * Arg:           seq [UNKN ] query data structure [Sequence*]
 * Arg:        signal [UNKN ] target data structure [RawSignalSeq*]
 * Arg:            sm [UNKN ] Resource [SignalMap*]
 * Arg:          para [UNKN ] Resource [RawSignalMatParaScore*]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
#define RawSignalMat_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+4)])    
RawSignalMat * allocate_Small_RawSignalMat(Sequence* seq,RawSignalSeq* signal ,SignalMap* sm,RawSignalMatParaScore* para) 
{
    RawSignalMat * out;  


    out = allocate_RawSignalMat_only(seq, signal , sm, para);    
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(40,(out->leni + 1) * 3,16,out->lenj+4);   
    if(out == NULL)  {  
      warn("Small shadow matrix RawSignalMat cannot be allocated, (asking for 5 by %d main cells)",out->leni+2); 
      free_RawSignalMat(out);    
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_RawSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for RawSignalMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_RawSignalMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_RawSignalMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_RawSignalMat(RawSignalMat * mat,DPEnvelope * dpenv) 
{
    int endj;    
    int score;   
    PackAln * out;   
    PackAlnUnit * pau;   
    int starti;  
    int startj;  
    int startstate;  
    int stopi;   
    int stopj;   
    int stopstate;   
    int temp;    
    int donej;  /* This is for reporting, will be passed as a & arg in */ 
    int totalj; /* This also is for reporting, but as is not changed, can be passed by value */ 


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW )    {  
      warn("Could not calculate packaln small for RawSignalMat due to wrong type of matrix");    
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_RawSignalMat(mat,dpenv); 
    score = start_end_find_end_RawSignalMat(mat,&endj);  
    out->score = score;  
    stopstate = END;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,0);    
    startj = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,1);    
    startstate = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,2);    
    stopi = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,3); 
    stopj = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,4); 
    stopstate = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,5); 
    temp = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,6);  
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);    


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_RawSignalMat(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);   
    if( stopi == RawSignalMat_READ_OFF_ERROR || stopstate != START ) {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_RawSignalMat(out,mat); 
      return out;    
      }  


    /* Ok. Put away start start... */ 
    pau = PackAlnUnit_alloc();   
    pau->i = stopi;  
    pau->j = stopj;  
    pau->state = stopstate + 3;  
    add_PackAln(out,pau);    


    log_full_error(REPORT,0,"Alignment recovered");  
    stop_reporting();    
    invert_PackAln(out); 
    recalculate_PackAln_RawSignalMat(out,mat);   
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_RawSignalMat(mat)
 *
 * Descrip:    This function calculates an alignment for RawSignalMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_RawSignalMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_RawSignalMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_RawSignalMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_RawSignalMat(RawSignalMat * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_RawSignalMat(mat,NULL);  
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_RawSignalMat(mat);    
    return out;  
}    


/* Function:  AlnRangeSet_from_RawSignalMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for RawSignalMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_RawSignalMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_RawSignalMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_RawSignalMat(RawSignalMat * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_RawSignalMat");  
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_RawSignalMat(mat,&jpos); 
    state = END; 


    while( (temp = AlnRange_build_RawSignalMat(mat,jpos,state,&jpos,&state)) != NULL)    
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_RawSignalMat(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_RawSignalMat
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_RawSignalMat(RawSignalMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_RawSignalMat");  
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_RawSignalMat(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)   {  
      warn("In AlnRanger_build_RawSignalMat alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);   
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0); 
    out->startj = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1); 
    out->startstate = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2); 
    out->stopi = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);  
    out->stopj = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);  
    out->stopstate = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);  
    out->startscore = RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6); 
    out->stopscore = RawSignalMat_DC_SHADOW_SPECIAL(mat,0,jpos,state);   


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_RawSignalMat(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL); 
    if( jpos == RawSignalMat_READ_OFF_ERROR) {  
      warn("In AlnRange_build_RawSignalMat alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);    
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = RawSignalMat_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate); 
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:               out [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean read_hidden_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
{
    int i;   
    int j;   
    int state;   
    int cellscore;   
    int isspecial;   
    /* We don't need hiddenj here, 'cause matrix access handled by max funcs */ 
    PackAlnUnit * pau;   


    /* stop position is on the path */ 
    i = stopi;   
    j = stopj;   
    state= stopstate;    
    isspecial = FALSE;   


    while( i >= starti && j >= startj)   {  
      /* Put away current i,j,state */ 
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state;   
      add_PackAln(out,pau);  


      max_hidden_RawSignalMat(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    


      if( i == RawSignalMat_READ_OFF_ERROR)  {  
        warn("In RawSignalMat hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);  
        return FALSE;    
        }  


      if( i == starti && j == startj && state == startstate) {  
/* Put away final state (start of this block) */ 
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state; 
        add_PackAln(out,pau);    
          return TRUE;   
        }  
      if( i == starti && j == startj)    {  
        warn("In RawSignalMat hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);    
        return FALSE;    
        }  
      }  
    warn("In RawSignalMat hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state); 
    return FALSE;    
}    


/* Function:  max_hidden_RawSignalMat(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:           hiddenj [UNKN ] Undocumented argument [int]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_hidden_RawSignalMat(RawSignalMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = RawSignalMat_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->seq->len || j > mat->signal->len) {  
      warn("In RawSignalMat matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = RawSignalMat_HIDDEN_MATRIX(mat,i,j,state);  
    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)) -  (0);    
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,ROGUE) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = ROGUE; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,ROGUE); 
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,ROGUE);  
          }  
        temp = cscore - ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext)) -  (0); 
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH); 
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);  
          }  
        temp = cscore - (((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal)) -  (0);   
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 4,BETWEEN) )    {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-4,BETWEEN);   
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 4,BETWEEN);    
          }  
        temp = cscore - ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal)) -  (0);   
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 3,BETWEEN) )    {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-3,BETWEEN);   
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 3,BETWEEN);    
          }  
        temp = cscore - (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal)) -  (0); 
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 2,BETWEEN) )    {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-2,BETWEEN);   
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 2,BETWEEN);    
          }  
        temp = cscore - ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal)) -  (0);   
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,BETWEEN) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,BETWEEN);   
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,BETWEEN);    
          }  
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case BETWEEN :     
        temp = cscore - (mat->para->between_event_signal) -  (0);    
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,BETWEEN) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,BETWEEN);   
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,BETWEEN);    
          }  
        temp = cscore - (0) -  (0);  
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-1,j-0,MATCH); 
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);  
          }  
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case ROGUE :   
        temp = cscore - (mat->para->rogue_ext) -  (0);   
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,ROGUE) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = ROGUE; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,ROGUE); 
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,ROGUE);  
          }  
        temp = cscore - (mat->para->rogue_signal) -  (0);    
        if( temp == RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH); 
            }  
          return RawSignalMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);  
          }  
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_RawSignalMat(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int *]
 * Arg:        startstate [UNKN ] Undocumented argument [int *]
 * Arg:               out [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean read_special_strip_RawSignalMat(RawSignalMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
{
    int i;   
    int j;   
    int state;   
    int cellscore;   
    int isspecial;   
    PackAlnUnit * pau;   


    /* stop position is on the path */ 
    i = stopi;   
    j = stopj;   
    state= stopstate;    
    isspecial = TRUE;    


    /* Loop until state has the same j as its stop in shadow pointers */ 
    /* This will be the state is came out from, OR it has hit !start */ 
    /* We may not want to get the alignment, in which case out will be NULL */ 
    while( j > RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START) { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 3; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_RawSignalMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    
      if( i == RawSignalMat_READ_OFF_ERROR)  {  
        warn("In special strip read RawSignalMat, got a bad read off error. Sorry!");    
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))    {  
      warn("In special strip read RawSignalMat, at special [%d] state [%d] overshot!",j,state);  
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 3;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_RawSignalMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_special_strip_RawSignalMat(RawSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = RawSignalMat_READ_OFF_ERROR; 
    if( isspecial == FALSE ) {  
      warn("In special strip max function for RawSignalMat, got a non special start point. Problem! (bad!)");    
      return (-1);   
      }  


    if( j < 0 || j > mat->signal->len)   {  
      warn("In RawSignalMat matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    cscore = RawSignalMat_DC_SHADOW_SPECIAL(mat,i,j,state);  
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* Source MATCH is not a special */ 
      default:   
        warn("Major problem (!) - in RawSignalMat special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_RawSignalMat(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:            cscore [UNKN ] Undocumented argument [int]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_matrix_to_special_RawSignalMat(RawSignalMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = RawSignalMat_READ_OFF_ERROR; 


    if( j < 0 || j > mat->lenj)  {  
      warn("In RawSignalMat matrix to special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (0);  
        if( temp == RawSignalMat_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - RawSignalMat_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);     
            }  
          return RawSignalMat_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;  
          }  
        /* Source ROGUE is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        /* Source BETWEEN is not a special, should not get here! */ 
        /* Source BETWEEN is not a special, should not get here! */ 
        /* Source BETWEEN is not a special, should not get here! */ 
        /* Source BETWEEN is not a special, should not get here! */ 
        warn("Major problem (!) - in RawSignalMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case BETWEEN :     
        /* Source BETWEEN is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in RawSignalMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case ROGUE :   
        /* Source ROGUE is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in RawSignalMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_RawSignalMat(mat,starti,startj,stopi,stopj);     


    RawSignalMat_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;    


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          RawSignalMat_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;  
          RawSignalMat_HIDDEN_MATRIX(mat,i,j,BETWEEN) = NEGI;    
          RawSignalMat_HIDDEN_MATRIX(mat,i,j,ROGUE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_HIDDEN_MATRIX(mat,i,j,MATCH) = score;  
        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + 0;   
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_HIDDEN_MATRIX(mat,i,j,BETWEEN) = score;    
        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;     
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_HIDDEN_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_HIDDEN_MATRIX(mat,i,j,ROGUE) = score;  
        /* Finished calculating state ROGUE */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_RawSignalMat(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_RawSignalMat(RawSignalMat * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-4);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        RawSignalMat_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
   
        RawSignalMat_HIDDEN_MATRIX(mat,i,j,BETWEEN) = NEGI;
 
        RawSignalMat_HIDDEN_MATRIX(mat,i,j,ROGUE) = NEGI;
   
        }  
      }  


    return;  
}    


/* Function:  full_dc_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_RawSignalMat
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_RawSignalMat to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [RawSignalMat *]
 * Arg:            starti [UNKN ] Start position in i [int]
 * Arg:            startj [UNKN ] Start position in j [int]
 * Arg:        startstate [UNKN ] Start position state number [int]
 * Arg:             stopi [UNKN ] Stop position in i [int]
 * Arg:             stopj [UNKN ] Stop position in j [int]
 * Arg:         stopstate [UNKN ] Stop position state number [int]
 * Arg:               out [UNKN ] PackAln structure to put alignment into [PackAln *]
 * Arg:             donej [UNKN ] pointer to a number with the amount of alignment done [int *]
 * Arg:            totalj [UNKN ] total amount of alignment to do (in j coordinates) [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean full_dc_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_RawSignalMat");    
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 20) {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);   
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE) {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)  {  
      warn("In divide and conquor for RawSignalMat, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj); 
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= RawSignalMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);  
    lstartj= RawSignalMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);  
    lstate = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);  


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_RawSignalMat(mat,RawSignalMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),RawSignalMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),RawSignalMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_RawSignalMat(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_RawSignalMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:         perc_done [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean do_dc_single_pass_RawSignalMat(RawSignalMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_RawSignalMat(mat);   


    RawSignalMat_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0; 
    run_up_dc_RawSignalMat(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);     
    push_dc_at_merge_RawSignalMat(mat,starti,stopi,halfj,&halfj,dpenv);  
    follow_on_dc_RawSignalMat(mat,starti,stopi,halfj,stopj,dpenv,perc_done);     
    return TRUE; 
}    


/* Function:  push_dc_at_merge_RawSignalMat(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_RawSignalMat(RawSignalMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<4;count++,j++)    {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100);    
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100);    
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = NEGI;     
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,0) = (-100);  
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,1) = (-100);  
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = NEGI;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,0) = (-100);    
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,1) = (-100);    
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);   
        if( j - 1 <= mergej) {  
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = BETWEEN;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH; 
          }  
        else {  
          for(k=0;k<7;k++)   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BETWEEN,k); 
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);  
        if( temp > score)    {  
          score = temp;  


          if( j - 2 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-2; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = BETWEEN; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,BETWEEN,k);   
            }  
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);    
        if( temp > score)    {  
          score = temp;  


          if( j - 3 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-3; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = BETWEEN; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,BETWEEN,k);   
            }  
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);    
        if( temp > score)    {  
          score = temp;  


          if( j - 4 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-4; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = BETWEEN; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,BETWEEN,k);   
            }  
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k); 
            }  
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-0; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = ROGUE;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,ROGUE,k); 
            }  
          }  
        /* Add any movement independant score */ 
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;    
        /* Finished with state MATCH */ 


        /* For state BETWEEN, pushing when j - offj <= mergej */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0;    
        if( j - 0 <= mergej) {  
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,0) = i-1; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,1) = j-0; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,2) = MATCH;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,3) = i;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,4) = j;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,5) = BETWEEN; 
          }  
        else {  
          for(k=0;k<7;k++)   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k); 
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,0) = i-0;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,1) = j-1;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,2) = BETWEEN;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,3) = i; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,4) = j; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,5) = BETWEEN;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BETWEEN,k); 
            }  
          }  
        /* Add any movement independant score */ 
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = score;  
        /* Finished with state BETWEEN */ 


        /* For state ROGUE, pushing when j - offj <= mergej */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;  
        if( j - 1 <= mergej) {  
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,0) = i-0;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,1) = j-1;   
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,2) = MATCH; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,3) = i; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,4) = j; 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,5) = ROGUE; 
          }  
        else {  
          for(k=0;k<7;k++)   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);   
          }  


        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,0) = i-0; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,1) = j-1; 
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,2) = ROGUE;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,3) = i;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,4) = j;   
            RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,5) = ROGUE;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,k) = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,ROGUE,k); 
            }  
          }  
        /* Add any movement independant score */ 
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = score;    
        /* Finished with state ROGUE */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_RawSignalMat(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_RawSignalMat(RawSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
{
    int i;   
    int j;   
    int k;   
    int score;   
    int temp;    
    int localshadow[7];  
    long int total;  
    long int num;    


    total = (stopi - starti+1) * (stopj - startj+1); 
    num = 0;     


    for(j=startj;j<=stopj;j++)   { /*for each valid j column*/ 
      for(i=starti;i<=stopi;i++) { /*this is strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = NEGI;     
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BETWEEN,k);  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 2,BETWEEN,k);    
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 3,BETWEEN,k);    
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 4,BETWEEN,k);    
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,ROGUE,k);  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        for(k=0;k<7;k++) 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0;    
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);    
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BETWEEN,k);    
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = score; 
        for(k=0;k<7;k++) 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);    
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,ROGUE,k);  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = score;   
        for(k=0;k<7;k++) 
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state ROGUE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_RawSignalMat(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_RawSignalMat(RawSignalMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    long int total;  
    long int num;    


    total = (stopi - starti+1) * (stopj - startj+1); 
    if( total <= 0 ) 
      total = 1; 
    num = 0;     


    for(j=startj;j<=stopj;j++)   { /*for each valid j column*/ 
      for(i=starti;i<=stopi;i++) { /*this is strip*/ 
        if( j == startj && i == starti)  
          continue;  
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = NEGI;     
          RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);   
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0;    
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = score; 
        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;  
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_DC_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = score;   
        /* Finished calculating state ROGUE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_RawSignalMat(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 */
}    
void init_dc_RawSignalMat(RawSignalMat * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<6;j++) {  
      for(i=(-1);i<mat->seq->len;i++)    {  
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI; 
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = NEGI;   
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = NEGI; 
        for(k=0;k<7;k++) {  
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1);  
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,k) = (-1);    
          RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,k) = (-1);  
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_RawSignalMat(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [RawSignalMat *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_RawSignalMat(RawSignalMat * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = RawSignalMat_DC_SHADOW_SPECIAL(mat,0,mat->signal->len-1,END);  
    maxj = mat->signal->len-1;   
    for(j= mat->signal->len-2 ;j >= 0 ;j--)  {  
      if( RawSignalMat_DC_SHADOW_SPECIAL(mat,0,j,END) > max )    {  
        max = RawSignalMat_DC_SHADOW_SPECIAL(mat,0,j,END);   
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_RawSignalMat(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [RawSignalMat]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_RawSignalMat(RawSignalMat *mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int score;   
    int temp;    
    int leni;    
    int lenj;    
    int localshadow[7];  
    long int total;  
    long int num=0;  
    int * score_pointers;    
    int * shadow_pointers;   
    int * localsp;   
    leni = mat->seq->len;    
    lenj = mat->signal->len; 
    total = leni * lenj; 


    score_pointers = (int *) calloc (4 * (leni + 1) * 3,sizeof(int));    
    shadow_pointers = (int *) calloc (4 * (leni + 1) * 3 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;   
          RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i,j,BETWEEN) = NEGI;     
          RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i,j,ROGUE) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal) + (0);     
        /* assign local shadown pointer */ 
        localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BETWEEN,0));    
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal) +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 2,BETWEEN,0));  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 3,BETWEEN,0));  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 4,BETWEEN,0));  
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0));    
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i) +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,ROGUE,0));    
          }  
        /* From state START to state MATCH */ 
        temp = RawSignalMat_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + 0 + (0);  
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= MATCH; 
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score;   
        for(k=0;k<7;k++) 
          RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > RawSignalMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          RawSignalMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            RawSignalMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);    
          RawSignalMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);  
          RawSignalMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          RawSignalMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          RawSignalMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + 0 + (0);  
        /* assign local shadown pointer */ 
        localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));  
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,BETWEEN,0));  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i,j,BETWEEN) = score; 
        for(k=0;k<7;k++) 
          RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal + (0);    
        /* assign local shadown pointer */ 
        localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0));  
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,ROGUE,0));    
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         RawSignalMat_DC_OPT_SHADOW_MATRIX(mat,i,j,ROGUE) = score;   
        for(k=0;k<7;k++) 
          RawSignalMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,ROGUE,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state ROGUE */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_RawSignalMat(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 */
void init_start_end_linear_RawSignalMat(RawSignalMat * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<6;j++) {  
      for(i=(-1);i<mat->seq->len;i++)    {  
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI; 
        RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1);    
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,BETWEEN) = NEGI;   
        RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,BETWEEN,0) = (-1);  
        RawSignalMat_DC_SHADOW_MATRIX(mat,i,j,ROGUE) = NEGI; 
        RawSignalMat_DC_SHADOW_MATRIX_SP(mat,i,j,ROGUE,0) = (-1);    
        }  
      }  


    for(j=(-4);j<mat->signal->len;j++)   {  
      RawSignalMat_DC_SHADOW_SPECIAL(mat,0,j,START) = 0; 
      RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;    
      RawSignalMat_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;    
      RawSignalMat_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);   
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_RawSignalMat(pal)
 *
 * Descrip:    Converts a path alignment to a label alignment
 *             The label alignment is probably much more useful than the path
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * convert_PackAln_to_AlnBlock_RawSignalMat(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_RawSignalMat();  
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "SEQUENCE","BETWEEN_NO_MOVE","BETWEEN_SEQ","ROGUE_SEQUENCE","END" };    
/* Function:  AlnConvertSet_RawSignalMat(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "SIGNAL_MATCH","BETWEEN_NO_MOVE","BETWEEN_SIGNAL","ROGUE_SIGNAL","END" };  
AlnConvertSet * AlnConvertSet_RawSignalMat(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BETWEEN;   
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BETWEEN;   
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 2;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BETWEEN;   
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 3;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BETWEEN;   
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 4;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = ROGUE; 
    acu->state2 = MATCH;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = BETWEEN;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = BETWEEN;   
    acu->state2 = BETWEEN;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = ROGUE;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = ROGUE; 
    acu->state2 = ROGUE;     
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    return out;  
}    


/* Function:  PackAln_read_Expl_RawSignalMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_RawSignalMat(RawSignalMat * mat) 
{
    RawSignalMat_access_func_holder holder;  


    holder.access_main    = RawSignalMat_explicit_access_main;   
    holder.access_special = RawSignalMat_explicit_access_special;    
    return PackAln_read_generic_RawSignalMat(mat,holder);    
}    


/* Function:  RawSignalMat_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int RawSignalMat_explicit_access_main(RawSignalMat * mat,int i,int j,int state) 
{
    return RawSignalMat_EXPL_MATRIX(mat,i,j,state);  
}    


/* Function:  RawSignalMat_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int RawSignalMat_explicit_access_special(RawSignalMat * mat,int i,int j,int state) 
{
    return RawSignalMat_EXPL_SPECIAL(mat,i,j,state); 
}    


/* Function:  PackAln_read_generic_RawSignalMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:          h [UNKN ] Undocumented argument [RawSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_RawSignalMat(RawSignalMat * mat,RawSignalMat_access_func_holder h) 
{
    register PackAln * out;  
    int i;   
    int j;   
    int state;   
    int cellscore = (-1);    
    boolean isspecial;   
    PackAlnUnit * pau = NULL;    
    PackAlnUnit * prev = NULL;   


    assert(mat);     
    assert(h.access_main);   
    assert(h.access_special);    


    out = PackAln_alloc_std();   
    if( out == NULL )    
      return NULL;   


    out->score =  find_end_RawSignalMat(mat,&i,&j,&state,&isspecial,h);  


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in RawSignalMat_basic_read, returning a mess.(Sorry!)",out->len); 
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 3;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_RawSignalMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      else   
        max_calc_RawSignalMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      if(i == RawSignalMat_READ_OFF_ERROR || j == RawSignalMat_READ_OFF_ERROR || state == RawSignalMat_READ_OFF_ERROR )  {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in RawSignalMat_basic_read, returning partial alignment",out->len); 
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 3;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_RawSignalMat(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [RawSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_RawSignalMat(RawSignalMat * mat,int * ri,int * rj,int * state,boolean * isspecial,RawSignalMat_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->signal->len-1,END); 
    maxj = mat->signal->len-1;   
    for(j= mat->signal->len-2 ;j >= 0 ;j--)  {  
      if( (temp =(*h.access_special)(mat,0,j,END)) > max )   {  
        max = temp;  
        maxj = j;    
        }  
      }  


    if( ri != NULL)  
       *ri = 0;  
    if( rj != NULL)  
       *rj = maxj;   
    if( state != NULL)   
       *state = END; 
    if( isspecial != NULL)   
       *isspecial = TRUE;    


    return max;  
}    


/* Function:  RawSignalMat_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void RawSignalMat_debug_show_matrix(RawSignalMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->seq->len;i++)   {  
      for(j=startj;j<stopj && j < mat->signal->len;j++)  {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH)); 
        fprintf(ofp,"State BETWEEN %d\n",RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN)); 
        fprintf(ofp,"State ROGUE %d\n",RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE)); 
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_RawSignalMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [RawSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_RawSignalMat(RawSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,RawSignalMat_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = RawSignalMat_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->seq->len || j > mat->signal->len) {  
      warn("In RawSignalMat matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_special)(mat,i - 1,j - 1,START) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-1,j-1,START);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,START);    
          }  
        temp = cscore - (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,ROGUE) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = ROGUE; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,ROGUE);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,ROGUE);    
          }  
        temp = cscore - ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH);    
          }  
        temp = cscore - (((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 4,BETWEEN) )  {  
          *reti = i - 0; 
          *retj = j - 4; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-4,BETWEEN); 
            }  
          return (*h.access_main)(mat,i - 0,j - 4,BETWEEN);  
          }  
        temp = cscore - ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 3,BETWEEN) )  {  
          *reti = i - 0; 
          *retj = j - 3; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-3,BETWEEN); 
            }  
          return (*h.access_main)(mat,i - 0,j - 3,BETWEEN);  
          }  
        temp = cscore - (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 0,j - 2,BETWEEN) )  {  
          *reti = i - 0; 
          *retj = j - 2; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-2,BETWEEN); 
            }  
          return (*h.access_main)(mat,i - 0,j - 2,BETWEEN);  
          }  
        temp = cscore - ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal)) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,BETWEEN) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,BETWEEN); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,BETWEEN);  
          }  
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case BETWEEN :     
        temp = cscore - (mat->para->between_event_signal) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,BETWEEN) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = BETWEEN;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,BETWEEN); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,BETWEEN);  
          }  
        temp = cscore - (0) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MATCH);    
          }  
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case ROGUE :   
        temp = cscore - (mat->para->rogue_ext) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,ROGUE) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = ROGUE; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,ROGUE);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,ROGUE);    
          }  
        temp = cscore - (mat->para->rogue_signal) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MATCH);    
          }  
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_RawSignalMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [RawSignalMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [RawSignalMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_RawSignalMat(RawSignalMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,RawSignalMat_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = RawSignalMat_READ_OFF_ERROR; 


    if( j < 0 || j > mat->signal->len)   {  
      warn("In RawSignalMat matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source MATCH is from main matrix */ 
        for(i= mat->seq->len-1;i >= 0 ;i--)  { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,MATCH) )  {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = MATCH;   
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,MATCH);     
              }  
            return (*h.access_main)(mat,i - 0,j - 0,MATCH) ;     
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in RawSignalMat read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_RawSignalMat(mat)
 *
 * Descrip:    This function calculates the RawSignalMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_RawSignalMat
 *
 *
 * Arg:        mat [UNKN ] RawSignalMat which contains explicit basematrix memory [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_RawSignalMat(RawSignalMat * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    long int tot;    
    long int num;    


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_RawSignalMat, passed a non Explicit matrix type, cannot calculate!");   
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("RawSignalMat Matrix calculation: ");    
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);    
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = RawSignalMat_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > RawSignalMat_EXPL_SPECIAL(mat,i,j,END) )  {  
          RawSignalMat_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_EXPL_MATRIX(mat,i-1,j-0,MATCH) + 0;     
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN) = score;  


        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;   
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE) = score;    


        /* Finished calculating state ROGUE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_RawSignalMat(mat,dpenv)
 *
 * Descrip:    This function calculates the RawSignalMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] RawSignalMat which contains explicit basematrix memory [RawSignalMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_RawSignalMat(RawSignalMat * mat,DPEnvelope * dpenv) 
{
    int i;   
    int j;   
    int k;   
    int starti;  
    int startj;  
    int endi;    
    int endj;    
    int tot; 
    int num; 
    int should_calc; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_RawSignalMat, passed a non Explicit matrix type, cannot calculate!");   
      return FALSE;  
      }  


    prepare_DPEnvelope(dpenv);   
    starti = dpenv->starti;  
    if( starti < 0 ) 
      starti = 0;    
    startj = dpenv->startj;  
    if( startj < 0 ) 
      startj = 0;    
    endi = dpenv->endi;  
    if( endi > mat->leni )   
      endi = mat->leni;  
    endj = dpenv->endj;  
    if( endj > mat->lenj )   
      endj = mat->lenj;  
    tot = (endi-starti) * (endj-startj); 
    num = 0; 


    for(j=startj-4;j<endj;j++)   {  
      for(i=1;i<mat->leni;i++)   {  
        RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;  
        RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN) = NEGI;    
        RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE) = NEGI;  
        }  
      }  
    for(j=-4;j<mat->lenj;j++)    {  
      RawSignalMat_EXPL_SPECIAL(mat,i,j,START) = 0;  
      RawSignalMat_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  


    start_reporting("RawSignalMat Matrix calculation: ");    
    for(j=startj;j<endj;j++) {  
      auto int score;    
      auto int temp;     
      for(i=starti;i<endi;i++)   {  
        /* Check if is in envelope - code identical to is_in_DPEnvelope, but aggressively inlined here for speed */ 
        should_calc = 0; 
        for(k=0;k<dpenv->len;k++)    {  
          auto DPUnit * u;   
          u = dpenv->dpu[k]; 
          switch(u->type)    {  
            case DPENV_RECT :    
              if( i >= u->starti && j >= u->startj && i <= (u->starti+u->height) && j <= (u->startj+u->length))  
                should_calc = 1;     
              break; 
            case DPENV_DIAG :    
              if(  abs( (i-j) - (u->starti-u->startj)) <= u->height && i+j >= u->starti+u->startj && i+j+u->length >= u->starti+u->startj)   
                should_calc = 1;     
              break; 
            }  
          if( should_calc == 1 ) 
            break;   
          }  
        if( should_calc == 0)    {  
          RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;    
          RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN) = NEGI;  
          RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE) = NEGI;    
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,BETWEEN) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->single_event_signal);    
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-2,BETWEEN) + ((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+mat->para->double_event_signal);   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-3,BETWEEN) + (((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+mat->para->triple_event_signal);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state BETWEEN to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-4,BETWEEN) + ((((Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-1),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-2),mat->seq,i))+Score_offset_RawSignalMap(mat->sm,mat->signal,(j-3),mat->seq,i))+mat->para->quad_event_signal);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + (Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i)+mat->para->event_ext);     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state ROGUE to state MATCH */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,ROGUE) + Score_offset_RawSignalMap(mat->sm,mat->signal,j,mat->seq,i);    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = RawSignalMat_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_EXPL_MATRIX(mat,i,j,MATCH) = score;    


        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > RawSignalMat_EXPL_SPECIAL(mat,i,j,END) )  {  
          RawSignalMat_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state MATCH */ 


        /* For state BETWEEN */ 
        /* setting first movement to score */ 
        score = RawSignalMat_EXPL_MATRIX(mat,i-1,j-0,MATCH) + 0;     
        /* From state BETWEEN to state BETWEEN */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,BETWEEN) + mat->para->between_event_signal;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for BETWEEN */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_EXPL_MATRIX(mat,i,j,BETWEEN) = score;  


        /* Finished calculating state BETWEEN */ 


        /* For state ROGUE */ 
        /* setting first movement to score */ 
        score = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->para->rogue_signal;   
        /* From state ROGUE to state ROGUE */ 
        temp = RawSignalMat_EXPL_MATRIX(mat,i-0,j-1,ROGUE) + mat->para->rogue_ext;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for ROGUE */ 
        /* Add any movement independant score and put away */ 
         RawSignalMat_EXPL_MATRIX(mat,i,j,ROGUE) = score;    


        /* Finished calculating state ROGUE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  RawSignalMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * RawSignalMat_alloc(void) 
{
    RawSignalMat * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RawSignalMat *) ckalloc (sizeof(RawSignalMat))) == NULL)    {  
      warn("RawSignalMat_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->basematrix = NULL;  
    out->shatter = NULL; 
    out->leni = 0;   
    out->lenj = 0;   


    return out;  
}    


/* Function:  free_RawSignalMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalMat *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalMat *]
 *
 */
RawSignalMat * free_RawSignalMat(RawSignalMat * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RawSignalMat obj. Should be trappable");  
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
    if( obj->basematrix != NULL) 
      free_BaseMatrix(obj->basematrix);  
    if( obj->shatter != NULL)    
      free_ShatterMatrix(obj->shatter);  
    /* obj->seq is linked in */ 
    /* obj->signal is linked in */ 
    /* obj->sm is linked in */ 
    /* obj->para is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
