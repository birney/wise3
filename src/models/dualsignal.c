#ifdef _cplusplus
extern "C" {
#endif
#include "dualsignal.h"



# line 83 "dualsignal.dy"
ComplexSequenceEvalSet * kmer_ComplexSequenceEvalSet(int kmer_size)
{
  ComplexSequenceEvalSet * out;

  out = ComplexSequenceEvalSet_alloc_len(10);

  /* put in bases and codons by default - we may well need them */

  add_ComplexSequenceEvalSet(out,base_number_ComplexSequenceEval());
  add_ComplexSequenceEvalSet(out,codon_number_ComplexSequenceEval());

  /* this is the more complex kmer */

  add_ComplexSequenceEvalSet(out,kmer_number_ComplexSequenceEval(kmer_size));

  
  prepare_ComplexSequenceEvalSet(out);

  return (out);

}


# line 106 "dualsignal.dy"
ComplexSequenceEval * kmer_number_ComplexSequenceEval(int kmer_size)
{
  ComplexSequenceEval * out;
  int * d;

  assert(kmer_size > 0 && kmer_size < 21);

  out = ComplexSequenceEval_alloc();

  out->left_window = kmer_size;
  out->right_window = 0;
  out->left_lookback = kmer_size +1;
  out->outside_score = 0;
  out->eval_func       = kmer_complexseq_eval;

  /* memory leak here like all ComplexSeq Evals. Sigh */
  d = malloc(sizeof(int));
  *d = kmer_size;
  (out->data) = (void *) d;

  return(out);

}

# line 130 "dualsignal.dy"
int kmer_complexseq_eval(int type,void * data,char * seq)
{
  int kmer_size;

  assert(seq != NULL);
  assert(data != NULL);

  kmer_size = * (int *) data;

  return( forward_dna_number_from_string(seq - kmer_size,kmer_size));

}

# line 143 "dualsignal.dy"
void show_alignment_with_fit_SimpleSignalMat(AlnBlock * alb,SignalEventList * sel,Sequence * comp,SignalMap * sm,FILE * ofp)
{
  AlnColumn * alc;
  int i = 1;

  char kbuf[10];
  
  double obs;
  double pred;

  kmer_t kmer;


  assert(alb != NULL);
  assert(sel != NULL);
  assert(comp != NULL);
  assert(ofp != NULL);



  fprintf(ofp,"AlnNum\tSigLabel\tSeqLabel\tRawScore\tBitsScore\tSigStart\tSigEnd\tSigMean\tSigStdev\tSigBase\tSigKmer\tSigTimeStart\tSigTimeLen\tSignalFit\tSeqStart\tSeqEnd\tSeqBase\tSeqKmer\tModelMean\n");
  for(alc=alb->start;alc != NULL;alc = alc->next) {


    if( strcmp(alc->alu[0]->text_label,"END") == 0 ) {
      break;
    }

    fprintf(ofp,"%d\t",i);

    
    kmer = forward_dna_number_from_string(comp->seq+alc->alu[1]->end - sm->kbasis,sm->kbasis);

    obs = sel->event[alc->alu[0]->end]->mean;
    pred = sm->comp[kmer]->mean;


    fprintf(ofp,"%s\t%s\t",alc->alu[0]->text_label,alc->alu[1]->text_label);

    fprintf(ofp,"%d\t%f\t",alc->alu[0]->score[0],Score2Bits(alc->alu[0]->score[0]));
    


    if( alc->alu[1]->start > 5 ) {
      strncpy(kbuf,comp->seq+alc->alu[1]->start-4,5);
    } else {
      strncpy(kbuf,"???????????",5);
    }



    fprintf(ofp,"%d\t%d\t%f\t%f\t%c\t%s\t%f\t%f\t",
	    alc->alu[0]->start,
	    alc->alu[0]->end,
	    sel->event[alc->alu[0]->end]->mean,
	    sel->event[alc->alu[0]->end]->std,
	    sel->event[alc->alu[0]->end]->base,
	    sel->event[alc->alu[0]->end]->kmer,
	    sel->event[alc->alu[0]->end]->time_pos,
	    sel->event[alc->alu[0]->end]->time_length);

    fprintf(ofp,"%f\t",(obs-pred)/sm->comp[kmer]->sd);

    fprintf(ofp,"%d\t%d\t%c\t%.5s\t%f",
	    alc->alu[1]->start,
	    alc->alu[1]->end,
	    comp->seq[alc->alu[1]->start],
	    kbuf,
	    pred);
      
    fprintf(ofp,"\n");
	    

    i++;
  }

}



# line 223 "dualsignal.dy"
SignalMap * reverse_SignalMap(SignalMap * sm)
{
  SignalMap * out;
  SignalComp * comp;
  int i;
  int j;
  
  assert(sm != NULL);

  out = SignalMap_alloc();
  out->comp = (SignalComp **) calloc(sizeof(SignalComp *),out->len);

  out->kbasis = sm->kbasis;
  out->len = sm->len;
  out->emission_length = sm->emission_length;
  out->emission_start  = sm->emission_start;
  out->emission_end    = sm->emission_end;
  out->emission_step   = sm->emission_step;

  for(i=0;i<sm->emission_length;i++) {
    out->em_null[i] = sm->em_null[i];
  }

  for(i=0;i<out->len;i++) {
    comp = SignalComp_alloc();
    out->comp[i] = comp;

    comp->kmer = reverse_complement_dna_number(sm->comp[i]->kmer,out->kbasis);
    reverse_map_dna_number(comp->kmer,out->kbasis,comp->kseq);
    comp->mean = sm->comp[i]->mean;
    comp->sd   = sm->comp[i]->sd;
    for(j=0;j<out->emission_length;j++) {
      comp->em_prob[j]   = sm->comp[i]->em_prob[j];
      comp->em_score[j]  = sm->comp[i]->em_score[j];
      comp->em_weight[j] = sm->comp[i]->em_weight[j];
    }
  }

  return(out);

}





# line 269 "dualsignal.dy"
Score Score_offset_RawSignalMap(SignalMap * sm,RawSignalSeq * raw,int i,Sequence * comp,int j)
{
  kmer_t kmer;
  int step;

  assert(sm != NULL);
  assert(raw != NULL);
  assert(comp != NULL);

  if( j <= sm->kbasis ) {
    return NEGI;
  }

  /*  kmer = forward_dna_number_from_string(comp->seq - sm->kbasis + j,sm->kbasis); */
  kmer = forward_dna_number_from_string(comp->seq - sm->kbasis + j,sm->kbasis);
  
  step = (int)((raw->signal[i] - sm->emission_start) / sm->emission_step);
/*
  fprintf(stdout,"Position %d,%d  signal %f, kmer step %d (start %f) score %d\n",i,j,sseq->signal[j],step,sm->emission_start, sm->comp[kmer]->em_score[step]);
*/



  if( step < 0 || step >= sm->emission_length ) {
    return NEGI;
  } else {
    /* fprintf(stderr,"for %d,%d, returning kmer %ld step %d - emission %d length %d \n",i,j,kmer,step,sm->emission_length,sm->len); */

    return( sm->comp[kmer]->em_score[step]);
  }

}

# line 302 "dualsignal.dy"
Score Score_offset_SignalMap(SignalMap * sm,SignalSeq * sseq,int i,Sequence * comp,int j)
{
  kmer_t kmer;
  int step;

  assert(sm != NULL);
  assert(sseq != NULL);
  assert(comp != NULL);

  if( j <= sm->kbasis ) {
    return NEGI;
  }

  /*  kmer = forward_dna_number_from_string(comp->seq - sm->kbasis + j,sm->kbasis); */
  kmer = forward_dna_number_from_string(comp->seq - sm->kbasis + j,sm->kbasis);
  
  step = (int)((sseq->signal[i] - sm->emission_start) / sm->emission_step);
/*
  fprintf(stdout,"Position %d,%d  signal %f, kmer step %d (start %f) score %d\n",i,j,sseq->signal[j],step,sm->emission_start, sm->comp[kmer]->em_score[step]);
*/



  if( step < 0 || step >= sm->emission_length ) {
    return NEGI;
  } else {
    /* fprintf(stderr,"for %d,%d, returning kmer %ld step %d - emission %d length %d \n",i,j,kmer,step,sm->emission_length,sm->len); */

    return( sm->comp[kmer]->em_score[step]);
  }

}

# line 335 "dualsignal.dy"
boolean flip_coin(Probability p)
{
  long rnd;

  if( p < 0.00000000001 ) {
    return 0;
  }

  rnd = random();

  if( ((double)rnd) / MAXRANDOM < p ) {
    return 1;
  } else {
    return 0;
  }

}

# line 353 "dualsignal.dy"
Probability generate_uniform(void) 
{
  long myrnd;

  myrnd = random();

  return ((double)myrnd / MAXRANDOM);
}

# line 362 "dualsignal.dy"
double draw_from_normal(double mean,double sd)
{
  double range;

  /*** this is not right ***/

  range = sd * 4;
  
  return(mean - (sd *2) + (generate_uniform() * range));
}

# line 373 "dualsignal.dy"
SignalSim * new_SignalSim(double insert,double stay_insert,double skip)
{
  SignalSim * out;

  out = SignalSim_alloc();

  out->skip = skip;
  out->insert = insert;
  out->stay_insert = stay_insert;

  return(out);
}

# line 386 "dualsignal.dy"
RawSignalSeq * simulated_RawSignalSeq_from_Seq(Sequence * seq,SignalMap * sm,SignalSim * sim,int length_of_event)
{
  Sequence * true_call;
  kmer_t kmer;
  int i;
  int j;
  int k;
  RawSignalSeq * out;
  
  char * tempseq;
  double * signal;
  double * signal_var;
  double * time_a;

  /** ERROR - not pleasant **/

  /*  fprintf(stderr,"Doing a hack here... %d",seq->len);*/


  i = 0;
  j = 0;


  signal = calloc(seq->len*2 * length_of_event,sizeof(double));
  time_a = calloc(seq->len*2 * length_of_event,sizeof(double));
  signal_var = calloc(seq->len*2 * length_of_event,sizeof(double));

  while( i < seq->len - sm->kbasis -1  ) {
    if( flip_coin(sim->skip) == 1 ) {
      i++;
      tempseq[j] = tolower(seq->seq[i]);
      j++;
      continue;
    }
    

    kmer = forward_dna_number_from_string(seq->seq+i,sm->kbasis);

    for(k=0;k<length_of_event;k++) {
        signal[j] = draw_from_normal(sm->comp[kmer]->mean,sm->comp[kmer]->sd);
	time_a[j] = 0.02;
	signal_var[j] = 0.0;
	j++;
    }

    i++;
  }


  out = RawSignalSeq_alloc();
  out->signal = signal;
  out->time_len = time_a;
  out->signal_var = signal_var;
  out->name = stringalloc(seq->name);
  out->start_time = 0.0;
  out->len = j;
  return(out);
 
}


# line 447 "dualsignal.dy"
SignalSeq * simulated_SignalSeq_from_Seq(Sequence * seq,SignalMap * sm,SignalSim * sim)
{
  Sequence * true_call;
  kmer_t kmer;
  int i;
  int j;
  SignalSeq * out;
  
  char * tempseq;
  double * signal;

  /** ERROR - not pleasant **/

  /*  fprintf(stderr,"Doing a hack here... %d",seq->len);*/


  i = 0;
  j = 0;

  tempseq = calloc(seq->len*2,sizeof(char));
  signal = calloc(seq->len*2,sizeof(double));

  while( i < seq->len - sm->kbasis -1  ) {
    if( flip_coin(sim->skip) == 1 ) {
      i++;
      tempseq[j] = tolower(seq->seq[i]);
      signal[j] = 0.0;
      j++;
      continue;
    }
    

    kmer = forward_dna_number_from_string(seq->seq+i,sm->kbasis);

    signal[j] = draw_from_normal(sm->comp[kmer]->mean,sm->comp[kmer]->sd);
    tempseq[j] = seq->seq[i];
    
    j++;
    i++;
  }


  out = SignalSeq_alloc();
  out->signal = signal;
  out->called = new_Sequence_from_strings(seq->name,tempseq);
  out->len = j;
  out->true_len = j;
  return(out);
 
}

# line 498 "dualsignal.dy"
void write_RawSignalSeq(RawSignalSeq * rss,FILE *ofp)
{
  int i;
  double time_i;

  assert(rss != NULL);
  assert(ofp != NULL);

  fprintf(ofp,">%s\n",rss->name);
  fprintf(ofp,"start\tlength\tmean\tstdev\n");

  time_i = rss->start_time;

  for(i=0;i<rss->len;i++) {
    fprintf(ofp,"%lf\t%lf\t%lf\t%lf\n",time_i,rss->time_len[i],rss->signal[i],rss->signal_var[i]);
    time_i += rss->time_len[i];
  }

  fprintf(ofp,"//\n");


}


# line 522 "dualsignal.dy"
RawSignalSeq * read_RawSignalSeq(FILE * ifp)
{
  char buffer[MAXLINE];
  RawSignalSeq * out;
  int i;
  char * seq;
  char name[MAXLINE];
  double * signal;	
  double * time_len;
  double * signal_var;
  char * runner;
  int len;
  char base;
  double sig;
  double time_start;
  double time_l;
  double sig_var;

  double time_cum;
  
  double start_time;

  

  assert(ifp != NULL);

  if( fgets(name,MAXLINE,ifp) == NULL ) {
    return NULL;
  }

  if( name[0] != '>' ) {
    warn("Bad SignalSequence, not starting with >");
  } else {
    runner = name+1;
    while( isalnum(*runner) || *runner == '_' ) {
      runner++;
    }

    (*runner) = '\0';
  }
  
  seq = calloc(DUALSIGNAL_INIT_SIZE,sizeof(char));
  signal = calloc(DUALSIGNAL_INIT_SIZE,sizeof(double));
  signal_var = calloc(DUALSIGNAL_INIT_SIZE,sizeof(double));
  time_len = calloc(DUALSIGNAL_INIT_SIZE,sizeof(double));
  len = DUALSIGNAL_INIT_SIZE;
    
  i = 0;

  assert(ifp != NULL);


  while( fgets(buffer,MAXLINE,ifp) != NULL) {

    if( buffer[0] == '/' && buffer[1] == '/' ) {
      break;
    }

    if( strstartcmp(buffer,"start") == 0 ) {
      continue;
    }

    /* this is what we expect start-time time-length mean-signal stddev */

    sscanf(buffer,"%lf %lf %lf %lf\n",&time_start,&time_l,&sig,&sig_var);

    

    /*    fprintf(stderr,"reading %c %f\n",base,sig); */

    signal[i] = sig;
    time_len[i] = time_l;
    signal_var[i] = sig_var;
    if( i == 0 ) {
      start_time = time_start;
      time_cum = start_time;
    } else {
      if( abs(time_cum + time_l - time_start) > 0.001 ) {
	warn("No adjacent times in Raw Signal input. Offsets are going to come out wrong!");
      } else {
	/* nothing */
      }
      
      time_cum += time_l;
    }

    

    i++;
    if( i >= len ) {
      signal = realloc(signal,len*2*sizeof(double));
      time_len = realloc(time_len,len*2*sizeof(double));
      signal_var = realloc(signal_var,len*2*sizeof(double));
      len = len*2;
    }
  }


  out = RawSignalSeq_alloc();
 
  out->name = stringalloc(name+1);

  out->signal = signal;
  out->time_len = time_len;
  out->signal_var = signal_var;
  out->start_time = start_time;

  out->len = i;

  return(out);

}



# line 637 "dualsignal.dy"
SignalSeq * read_SignalSeq(FILE * ifp)
{
  char buffer[MAXLINE];
  SignalSeq * out;
  int i;
  char * seq;
  char name[MAXLINE];
  double * signal;
  char * runner;
  int len;
  char base;
  double sig;

  

  assert(ifp != NULL);

  if( fgets(name,MAXLINE,ifp) == NULL ) {
    return NULL;
  }

  if( name[0] != '>' ) {
    warn("Bad SignalSequence, not starting with >");
  } else {
    runner = name+1;
    while( isalnum(*runner) || *runner == '_' ) {
      runner++;
    }

    (*runner) = '\0';
  }
  
  seq = calloc(1024,sizeof(char));
  signal = calloc(1024,sizeof(double));
  len = 1024;
    
  i = 0;

  assert(ifp != NULL);


  while( fgets(buffer,MAXLINE,ifp) != NULL) {

    if( buffer[0] == '/' && buffer[1] == '/' ) {
      break;
    }

    sscanf(buffer,"%c %lf\n",&base,&sig);

    if( !is_base(base) ) {
      warn("Bad base %c is signal sequence, skipping, but not ideal");
      continue;
    }

    /*    fprintf(stderr,"reading %c %f\n",base,sig); */

    seq[i] = base;
    signal[i] = sig;
    i++;
    if( i >= len ) {
      seq = realloc(seq,len*2*sizeof(char));
      signal = realloc(signal,len*2*sizeof(double));
      len = len*2;
    }
  }


  out = SignalSeq_alloc();
  out->called = Sequence_alloc();
 
  out->called->name = stringalloc(name+1);
  out->called->seq  = seq;
  out->called->len  = i;

  fprintf(stderr,"Got sequence %c%c%c length %d\n",out->called->seq[0],out->called->seq[1],out->called->seq[2],i);

  out->signal = signal;
  out->len = i;
  out->true_len = i;
  return(out);

}

# line 720 "dualsignal.dy"
void write_SignalSeq(SignalSeq * ss,FILE * ofp)
{
  int i;

  assert(ss != NULL);
  assert(ofp != NULL);

  fprintf(ofp,">%s\n",ss->called->name);
  for(i=0;i<ss->len;i++) {
    fprintf(ofp,"%c\t%f\n",ss->called->seq[i],ss->signal[i]);
  }

  fprintf(ofp,"//\n");
 
}

# line 736 "dualsignal.dy"
void write_SignalMap_absolute(SignalMap * sm,FILE * ofp)
{
  int i;
  int j;

  assert(sm != NULL);
  assert(ofp != NULL);

  
  fprintf(ofp,"# Written by DualSignal\n");
  fprintf(ofp,"kbasis %d\n",sm->kbasis);
  fprintf(ofp,"estart %lf\n",sm->emission_start);
  fprintf(ofp,"eend %lf\n",sm->emission_end);
  fprintf(ofp,"estep %lf\n",sm->emission_step);
  fprintf(ofp,"//\n");
  

  for(i=0;i<sm->len;i++) {
    fprintf(ofp,"%s",sm->comp[i]->kseq);
    for(j=0;j<sm->emission_length;j++) {
      fprintf(ofp,"\t%lf",sm->comp[i]->em_weight[j]);
    }
    fprintf(ofp,"\n");
  }

  return;
}

# line 764 "dualsignal.dy"
void dump_SignalMap_Scores(SignalMap * sm, FILE * ofp)
{
  int i;
  int j;

  assert(sm != NULL);
  assert(ofp != NULL);

  
  fprintf(ofp,"# Written by DualSignal\n");
  fprintf(ofp,"kbasis %d\n",sm->kbasis);
  fprintf(ofp,"estart %lf\n",sm->emission_start);
  fprintf(ofp,"eend %lf\n",sm->emission_end);
  fprintf(ofp,"estep %lf\n",sm->emission_step);
  fprintf(ofp,"//\n");
  

  for(i=0;i<sm->len;i++) {
    fprintf(ofp,"%s",sm->comp[i]->kseq);
    for(j=0;j<sm->emission_length;j++) {
      fprintf(ofp,"\t%d",sm->comp[i]->em_score[j]);
    }
    fprintf(ofp,"\n");
  }

  return;  


}

# line 794 "dualsignal.dy"
void write_SignalMap_normal(SignalMap * sm,FILE * ofp)
{
  int i;


  assert(sm != NULL);
  assert(ofp != NULL);

  
  fprintf(ofp,"# DualSignal dump\n");
  fprintf(ofp,"kbasis %d\n",sm->kbasis);

  for(i=0;i<sm->len;i++) {
    fprintf(ofp,"%s %lf %lf\n",sm->comp[i]->kseq,sm->comp[i]->mean,sm->comp[i]->sd);
  }

}

# line 812 "dualsignal.dy"
void prepare_SignalMap(SignalMap * sm,double pseudocount)
{
  int i;
  int j;
  double total;
  double weight_total[MAX_SIGNAL_EMISSION];
  double overall_total = 0.0;
   
  assert(sm != NULL);

  for(i=0;i<sm->emission_length;i++) {
    weight_total[i] = 0;
  }

  for(i=0;i<sm->len;i++) {
    total = 0.0;

    for(j=0;j<sm->emission_length;j++) {
      total += sm->comp[i]->em_weight[j] + pseudocount;
      weight_total[j] += sm->comp[i]->em_weight[j] + pseudocount;
      overall_total += sm->comp[i]->em_weight[j] + pseudocount;
    }

    for(j=0;j<sm->emission_length;j++) {
      sm->comp[i]->em_prob[j] = (sm->comp[i]->em_weight[j] + pseudocount) / total;
    }
    
  }

  for(j=0;j<sm->emission_length;j++) {
    sm->em_null[j] = weight_total[j]/ overall_total;
  }
  
  for(i=0;i<sm->len;i++) {

    for(j=0;j<sm->emission_length;j++) {
      sm->comp[i]->em_score[j] = Probability2Score(sm->comp[i]->em_prob[j] / sm->em_null[j]);
    }
    
  }
  
}

# line 855 "dualsignal.dy"
void convert_normal_to_absolute_SignalMap(SignalMap * sm,double weight)
{
  int i,j;
  double min_em;
  double max_em;
  int emission_length = 400;
  double total;
  
  assert(sm != NULL);

  /* find emission ranges */
  min_em = sm->comp[0]->mean - (sm->comp[0]->sd * 4) ;
  max_em = sm->comp[0]->mean + (sm->comp[0]->sd * 4);
  for(i=1;i<sm->len;i++) {
    if( sm->comp[i]->mean - (sm->comp[i]->sd * 4) < min_em ) {
      min_em = sm->comp[i]->mean - (sm->comp[i]->sd * 4);
    } 
    if( sm->comp[i]->mean + (sm->comp[i]->sd * 4) > max_em ) {
      max_em = sm->comp[i]->mean + (sm->comp[i]->sd * 4);
    } 
  }

  sm->emission_start = min_em;
  sm->emission_end   = max_em;

  sm->emission_length = emission_length;
  sm->emission_step = (max_em - min_em) / emission_length;

  /* now put in weights into each comp */

  for(i=0;i< sm->len;i++) {
    for(j=0;j< emission_length;j++) {
      /* everything gets a 0.0 weight - pseudocounts are added in prepare */
      sm->comp[i]->em_weight[j] = 0.0;
      /*
      fprintf(stdout,"Got %d step %d mean %f sd %f measure %f\n",
	      i,j,sm->comp[i]->mean,sm->comp[i]->sd,sm->emission_start + j*sm->emission_step);
      */


      if( (sm->emission_start + j*sm->emission_step) > (sm->comp[i]->mean - 4*sm->comp[i]->sd) &&
	  (sm->emission_start + j*sm->emission_step) < (sm->comp[i]->mean + 4*sm->comp[i]->sd) ) {
	/* is in distribution */
	 auto double prob_weight;

	/* take j to j+1 probability and multiply by weight */

	 prob_weight = ((crude_normal_probability(sm->comp[i]->mean,sm->comp[i]->sd,sm->emission_start + (j+1)*sm->emission_step)) - (crude_normal_probability(sm->comp[i]->mean,sm->comp[i]->sd,sm->emission_start + (j)*sm->emission_step))) ;
	 /*
	 fprintf(stdout,"For position %d, mean %f sd %f, step %d, measure %f, weight %f\n",
		 i,
		 sm->comp[i]->mean,
		 sm->comp[i]->sd,
		 j,
		 sm->emission_start + j*sm->emission_step,
		 prob_weight);
	 */

	sm->comp[i]->em_weight[j] += weight * prob_weight ;
	  }
    }
  }

  /*  convert weights to probabilities */
    
  convert_weight_to_Probability_SignalMap(sm);

  /*  make an average Probability for the NULL over k */


  for(j=0;j<sm->emission_length;j++) {
    total = 0.0;
    for(i=0;i<sm->len;i++) {
      total += sm->comp[i]->em_prob[j];
    }
    sm->em_null[j] = total / sm->len;
  }


  /* flip to score */

 
  for(i=0;i<sm->len;i++) {
    
    for(j=0;j<sm->emission_length;j++) {
      sm->comp[i]->em_score[j] = Probability2Score(sm->comp[i]->em_prob[j] / sm->em_null[j]);
    }
  }
  
  /* return */    


}


# line 950 "dualsignal.dy"
void convert_weight_to_Probability_SignalMap(SignalMap * sm)
{
  int i;
  int j;
  double total;

  assert(sm != NULL);

  for(i=0;i<sm->len;i++) {
    total = 0.0;

    for(j=0;j<sm->emission_length;j++) {
      total += sm->comp[i]->em_weight[j];
    }

    for(j=0;j<sm->emission_length;j++) {
      sm->comp[i]->em_prob[j] = sm->comp[i]->em_weight[j] / total;
    }

  }

}

# line 973 "dualsignal.dy"
SignalMap * read_SignalMap_absolute(FILE * ifp)
{
  char buffer[MAXLINE];
  int kb  = -1;
  double estart;
  double eend;
  double estep;
  boolean seen_start = 0;
  boolean seen_end   = 0;
  boolean seen_step  = 0;

  char kstring[MAXLINE];
  
  int i;
  int len;

  SignalMap * out = NULL;
  SignalComp * comp;

  char * runner;
  double temp;
  int count;
  
  assert(ifp != NULL);

 

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' ) {
      continue;
    } 
    if( strstartcmp(buffer,"kbasis") == 0 ) {
      sscanf(buffer,"kbasis %d",&kb);
      assert(kb < 8);
      assert(kb > 0);
    }
    if( strstartcmp(buffer,"estart") == 0 ) {
      sscanf(buffer,"estart %lf",&estart);
      seen_start = 1;
    }
    if( strstartcmp(buffer,"eend") == 0 ) {
      sscanf(buffer,"eend %lf",&eend);
      seen_end = 1;
    }
    if( strstartcmp(buffer,"estep") == 0 ) {
      sscanf(buffer,"estep %lf",&estep);
      seen_step = 1;
    }
    if( strstartcmp(buffer,"//") == 0 ) {
      break;
    }
  }


  for(len=1, i = 0;i<kb;i++) {
    len = len * 4;
  }
  
  out = SignalMap_alloc_len(len);
  out->kbasis = kb;
  

  assert(seen_start != 0 && seen_end != 0 && seen_step != 0 && kb != -1 );

  out = SignalMap_alloc_len(len);
  out->kbasis = kb;
  out->emission_start = estart;
  out->emission_end   = eend;
  out->emission_step  = estep;
  out->emission_length = ceil((eend-estart) / estep);

  assert(out->emission_length < MAX_SIGNAL_EMISSION);


  assert(out != NULL);


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    sscanf(buffer,"%s",kstring);
    comp = SignalComp_alloc();
    strcpy(comp->kseq,kstring);
    comp->kmer = forward_dna_number_from_string(comp->kseq,out->kbasis);
    out->comp[comp->kmer] = comp;

    /* move past the first string */

    for(runner = buffer;*runner && !isspace(*runner);runner++)
      ;

    /* now we expect emission length numbers */

    count = 0;
    while( count < out->emission_length ) {

      /* eat white space */
      for(;*runner && isspace(*runner);runner++)
	;
      
      /* slurp up number */
      if( sscanf(runner,"%lf",&temp) == 0 ) {
	warn("Bad number in count, bad line offset?, count %d, line %s",count,buffer);
	return NULL;
      }

      /* go through number */
      for(;*runner && !isspace(*runner);runner++)
	;

      /* eat white space */
      for(;*runner && isspace(*runner);runner++)
	;

      comp->em_weight[count] = temp;
      count++;
    }

  }

  return(out);
}

# line 1094 "dualsignal.dy"
SignalMap * read_SignalMap_tsv(int ksize,FILE * ifp)
{
 char buffer[MAXLINE];
  int kb;
  char kstring[MAXLINE];
  double mean;
  double sd;
  int dummy;

  int len = 0;
  int i;

  SignalMap * out = NULL;
  SignalComp * comp;

  assert(ifp != NULL);
  assert(ksize > 0);


  for(len=1, i = 0;i<ksize;i++) {
    len = len * 4;
  }
  
  out = SignalMap_alloc_len(len);
  for(i=0;i<len;i++) {
    out->comp[i] = NULL;
  }

  out->kbasis = ksize;
  
  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    /* tab delimited, kmer, variant, level_mean level_stdy */
    sscanf(buffer,"%s\t%d\t%lf\t%lf",kstring,&dummy,&mean,&sd);
    if( dummy != 1 ) {
      warn("Bad tsv model string %s",buffer);
      continue;
    }

    comp = SignalComp_alloc();
    strcpy(comp->kseq,kstring);
    comp->kmer = forward_dna_number_from_string(comp->kseq,out->kbasis);
    comp->mean = mean;
    comp->sd = sd;
    /* fprintf(stderr,"Storing for kmer %d, %s %d\n",comp->kmer,comp->kseq,out->kbasis);*/
    out->comp[comp->kmer] = comp;

  }

  

  return out;
}


# line 1148 "dualsignal.dy"
boolean validate_SignalMap_with_warnings(SignalMap * sm)
{
  int i;
  boolean ret = TRUE;

  char retstr[1024];


  assert(sm != NULL);

  retstr[sm->kbasis] = '\0';

  for(i=0;i<sm->len;i++) {
    if( sm->comp[i] == NULL ) {
      reverse_map_dna_number(i,sm->kbasis,retstr);
      warn("Kmer %d (%s) is NULL. SignalMap incomplete",i,retstr);
      ret = FALSE;
    }
  }
 

  return(ret);
}


# line 1173 "dualsignal.dy"
SignalMap * read_SignalMap_normal(FILE * ifp)
{
  char buffer[MAXLINE];
  int kb;
  char kstring[MAXLINE];
  double mean;
  double sd;

  int len = 0;
  int i;

  SignalMap * out = NULL;
  SignalComp * comp;

  assert(ifp != NULL);

  

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' ) {
      continue;
    } 
    if( strstartcmp(buffer,"kbasis") == 0 ) {
      sscanf(buffer,"kbasis %d",&kb);
      assert(kb < 8);
      assert(kb > 0);

      for(len=1, i = 0;i<kb;i++) {
	len = len * 4;
      }

      out = SignalMap_alloc_len(len);
      out->kbasis = kb;

      
      break;
    }
  }


  assert(out != NULL);

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    sscanf(buffer,"%s %lf %lf",kstring,&mean,&sd);
    comp = SignalComp_alloc();
    strcpy(comp->kseq,kstring);
    comp->kmer = forward_dna_number_from_string(comp->kseq,out->kbasis);
    comp->mean = mean;
    comp->sd = sd;
    /* fprintf(stderr,"Storing for kmer %d, %s %d\n",comp->kmer,comp->kseq,out->kbasis);*/
    out->comp[comp->kmer] = comp;
  }


  return(out);
}



# line 1232 "dualsignal.dy"
SignalMap * SignalMap_alloc_len(int len)
{
  SignalMap * out = NULL;

  assert(len > 0 );

  out = SignalMap_alloc();

  out->comp = calloc(len,sizeof(SignalComp*));
  out->len = len;

  return out;
}


# line 1247 "dualsignal.dy"
SignalSeq * SignalSeq_from_SignalEventList(SignalEventList * sel)
{
  int i;
  SignalSeq * out;


  assert(sel != NULL);

  out = SignalSeq_alloc();
  out->signal = calloc(sizeof(double),sel->len);
  out->len = sel->len;
  
  out->called = Sequence_alloc();
  out->called->seq = calloc(sizeof(char),sel->len);

  out->called->name = stringalloc(sel->name);

  
  for(i=0;i<sel->len;i++) {
    out->signal[i] = sel->event[i]->mean;
    out->called->seq[i] = sel->event[i]->base;
  }

  
  return(out);

}

# line 1275 "dualsignal.dy"
void write_SignalEventList(SignalEventList * sel,FILE * ofp)
{
  int i;


  assert(ofp != NULL);
  assert(sel != NULL);

  fprintf(ofp,">%s\n",sel->name);

  for(i=0;i<sel->len;i++) {
    fprintf(ofp,"%c\t%lf\t%lf\t%s\t%lf\t%lf\n",
	    sel->event[i]->base,
	    sel->event[i]->mean,
	    sel->event[i]->std,
	    sel->event[i]->kmer,
	    sel->event[i]->time_pos,
	    sel->event[i]->time_length);
  }

  fprintf(ofp,"//\n");
}


# line 1299 "dualsignal.dy"
SignalEventList * read_SignalEventList(FILE * ifp)
{
  char buffer[MAXLINE];
  SignalEventList * out;
  SignalEvent * event;
  int i;

  char name[MAXLINE];
  
  double * signal;
  char * runner;
  int len;
  char base;
  double sig;
  double std;
  char kmer_str[MAXLINE];
  double time_pos;
  double time_len;

  

  assert(ifp != NULL);

  out = SignalEventList_alloc_std();


  if( fgets(name,MAXLINE,ifp) == NULL ) {
    return NULL;
  }

  if( name[0] != '>' ) {
    warn("Bad SignalSequence, not starting with >");
  } else {
    runner = name+1;
    while( isalnum(*runner) || *runner == '_' ) {
      runner++;
    }

    (*runner) = '\0';
  }
  
  out->name = stringalloc(name+1);
    
  while( fgets(buffer,MAXLINE,ifp) != NULL) {

    if( buffer[0] == '/' && buffer[1] == '/' ) {
      break;
    }

    sscanf(buffer,"%c %lf %lf %s %lf %lf\n",&base,&sig,&std,kmer_str,&time_pos,&time_len);

    if( !is_base(base) ) {
      warn("Bad base %c is signal sequence, skipping, but not ideal");
      continue;
    }

    event = SignalEvent_alloc();

    event->base  = base;
    event->mean = sig;
    event->std     = std;
    strncpy(event->kmer,kmer_str,MAX_EVENT_KMER);
    event->time_pos = time_pos;
    event->time_length = time_len;

    add_SignalEventList(out,event); 

  }


  return(out);

}



# line 1330 "dualsignal.c"
/* Function:  hard_link_SignalEvent(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalEvent *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEvent *]
 *
 */
SignalEvent * hard_link_SignalEvent(SignalEvent * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SignalEvent object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SignalEvent_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalEvent *]
 *
 */
SignalEvent * SignalEvent_alloc(void) 
{
    SignalEvent * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SignalEvent *) ckalloc (sizeof(SignalEvent))) == NULL)  {  
      warn("SignalEvent_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->base = 'u'; 
    out->mean = 0;   
    out->std = 0;    
    /* kmer[MAX_EVENT_KMER] is an array: no default possible */ 
    out->time_pos = 0;   
    out->time_length = 0;    


    return out;  
}    


/* Function:  free_SignalEvent(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalEvent *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEvent *]
 *
 */
SignalEvent * free_SignalEvent(SignalEvent * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SignalEvent obj. Should be trappable");   
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


/* Function:  swap_SignalEventList(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_SignalEventList
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [SignalEvent **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_SignalEventList(SignalEvent ** list,int i,int j)  
{
    SignalEvent * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_SignalEventList(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_SignalEventList which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [SignalEvent **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_SignalEventList(SignalEvent ** list,int left,int right,int (*comp)(SignalEvent * ,SignalEvent * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_SignalEventList(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_SignalEventList (list,++last,i);    
      }  
    swap_SignalEventList (list,left,last);   
    qsort_SignalEventList(list,left,last-1,comp);    
    qsort_SignalEventList(list,last+1,right,comp);   
}    


/* Function:  sort_SignalEventList(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_SignalEventList
 *
 *
 * Arg:         obj [UNKN ] Object containing list [SignalEventList *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_SignalEventList(SignalEventList * obj,int (*comp)(SignalEvent *, SignalEvent *)) 
{
    qsort_SignalEventList(obj->event,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_SignalEventList(obj,len)
 *
 * Descrip:    Really an internal function for add_SignalEventList
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SignalEventList *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_SignalEventList(SignalEventList * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_SignalEventList called with no need");    
      return TRUE;   
      }  


    if( (obj->event = (SignalEvent ** ) ckrealloc (obj->event,sizeof(SignalEvent *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_SignalEventList, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_SignalEventList(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [SignalEventList *]
 * Arg:        add [OWNER] Object to add to the list [SignalEvent *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_SignalEventList(SignalEventList * obj,SignalEvent * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_SignalEventList(obj,obj->len + SignalEventListLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->event[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_SignalEventList(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [SignalEventList *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_SignalEventList(SignalEventList * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->event[i] != NULL) {  
        free_SignalEvent(obj->event[i]); 
        obj->event[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  SignalEventList_alloc_std(void)
 *
 * Descrip:    Equivalent to SignalEventList_alloc_len(SignalEventListLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * SignalEventList_alloc_std(void) 
{
    return SignalEventList_alloc_len(SignalEventListLISTLENGTH); 
}    


/* Function:  SignalEventList_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * SignalEventList_alloc_len(int len) 
{
    SignalEventList * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = SignalEventList_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->event = (SignalEvent ** ) ckcalloc (len,sizeof(SignalEvent *))) == NULL)    {  
      warn("Warning, ckcalloc failed in SignalEventList_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_SignalEventList(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalEventList *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * hard_link_SignalEventList(SignalEventList * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SignalEventList object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SignalEventList_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * SignalEventList_alloc(void) 
{
    SignalEventList * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SignalEventList *) ckalloc (sizeof(SignalEventList))) == NULL)  {  
      warn("SignalEventList_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->event = NULL;   
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_SignalEventList(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalEventList *]
 *
 * Return [UNKN ]  Undocumented return value [SignalEventList *]
 *
 */
SignalEventList * free_SignalEventList(SignalEventList * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SignalEventList obj. Should be trappable");   
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
    if( obj->event != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->event[i] != NULL)   
          free_SignalEvent(obj->event[i]);   
        }  
      ckfree(obj->event);    
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_RawSignalSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [RawSignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalSeq *]
 *
 */
RawSignalSeq * hard_link_RawSignalSeq(RawSignalSeq * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a RawSignalSeq object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  RawSignalSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [RawSignalSeq *]
 *
 */
RawSignalSeq * RawSignalSeq_alloc(void) 
{
    RawSignalSeq * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(RawSignalSeq *) ckalloc (sizeof(RawSignalSeq))) == NULL)    {  
      warn("RawSignalSeq_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->signal = NULL;  
    out->time_len = NULL;    
    out->signal_var = NULL;  
    out->len = 0;    
    out->start_time = 0; 
    out->name = NULL;    


    return out;  
}    


/* Function:  free_RawSignalSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [RawSignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [RawSignalSeq *]
 *
 */
RawSignalSeq * free_RawSignalSeq(RawSignalSeq * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a RawSignalSeq obj. Should be trappable");  
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
    if( obj->signal != NULL) 
      ckfree(obj->signal);   
    if( obj->time_len != NULL)   
      ckfree(obj->time_len);     
    if( obj->signal_var != NULL) 
      ckfree(obj->signal_var);   
    if( obj->name != NULL)   
      ckfree(obj->name);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_SignalSeq(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSeq *]
 *
 */
SignalSeq * hard_link_SignalSeq(SignalSeq * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SignalSeq object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SignalSeq_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalSeq *]
 *
 */
SignalSeq * SignalSeq_alloc(void) 
{
    SignalSeq * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SignalSeq *) ckalloc (sizeof(SignalSeq))) == NULL)  {  
      warn("SignalSeq_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->signal = NULL;  
    out->len = 0;    
    out->true_len = 0;   
    out->called = NULL;  


    return out;  
}    


/* Function:  free_SignalSeq(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalSeq *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSeq *]
 *
 */
SignalSeq * free_SignalSeq(SignalSeq * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SignalSeq obj. Should be trappable"); 
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
    if( obj->signal != NULL) 
      ckfree(obj->signal);   
    if( obj->called != NULL) 
      free_Sequence(obj->called);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_SignalComp(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalComp *]
 *
 * Return [UNKN ]  Undocumented return value [SignalComp *]
 *
 */
SignalComp * hard_link_SignalComp(SignalComp * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SignalComp object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SignalComp_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalComp *]
 *
 */
SignalComp * SignalComp_alloc(void) 
{
    SignalComp * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SignalComp *) ckalloc (sizeof(SignalComp))) == NULL)    {  
      warn("SignalComp_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->kmer = 0;   
    /* kseq[20] is an array: no default possible */ 
    out->mean = 0;   
    out->sd = 0; 
    /* em_prob[MAX_SIGNAL_EMISSION] is an array: no default possible */ 
    /* em_score[MAX_SIGNAL_EMISSION] is an array: no default possible */ 
    /* em_weight[MAX_SIGNAL_EMISSION] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_SignalComp(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalComp *]
 *
 * Return [UNKN ]  Undocumented return value [SignalComp *]
 *
 */
SignalComp * free_SignalComp(SignalComp * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SignalComp obj. Should be trappable");    
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


/* Function:  hard_link_SignalMap(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalMap *]
 *
 * Return [UNKN ]  Undocumented return value [SignalMap *]
 *
 */
SignalMap * hard_link_SignalMap(SignalMap * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SignalMap object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SignalMap_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalMap *]
 *
 */
SignalMap * SignalMap_alloc(void) 
{
    SignalMap * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SignalMap *) ckalloc (sizeof(SignalMap))) == NULL)  {  
      warn("SignalMap_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->kbasis = 0; 
    out->comp = NULL;    
    out->len = 0;    
    out->emission_length = 0;    
    out->emission_start = 0; 
    out->emission_end = 0;   
    out->emission_step = 0;  
    /* em_null[MAX_SIGNAL_EMISSION] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_SignalMap(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalMap *]
 *
 * Return [UNKN ]  Undocumented return value [SignalMap *]
 *
 */
SignalMap * free_SignalMap(SignalMap * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SignalMap obj. Should be trappable"); 
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
    if( obj->comp != NULL)   
      ckfree(obj->comp);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_SignalSim(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SignalSim *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSim *]
 *
 */
SignalSim * hard_link_SignalSim(SignalSim * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SignalSim object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SignalSim_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SignalSim *]
 *
 */
SignalSim * SignalSim_alloc(void) 
{
    SignalSim * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SignalSim *) ckalloc (sizeof(SignalSim))) == NULL)  {  
      warn("SignalSim_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->insert = 0.0;   
    out->stay_insert = 0.0;  
    out->skip = 0.0; 


    return out;  
}    


/* Function:  free_SignalSim(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SignalSim *]
 *
 * Return [UNKN ]  Undocumented return value [SignalSim *]
 *
 */
SignalSim * free_SignalSim(SignalSim * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SignalSim obj. Should be trappable"); 
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
