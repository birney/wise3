#ifdef _cplusplus
extern "C" {
#endif
#include "startendmodel.h"


# line 31 "startendmodel.dy"
void free_Probability(Probability * p)
{
  return(free(p));
}

# line 36 "startendmodel.dy"
void free_Score(Score * s)
{
  return(free(s));
}

# line 41 "startendmodel.dy"
StartEndModelMode * StartEndModelMode_from_argc(int * argc,char ** argv)
{
  StartEndModelMode * out;
  char * style;

  style = strip_out_assigned_argument(argc,argv,"startendstyle");

  out = StartEndModelMode_alloc();
  
  if( strcmp(style,"edge") == 0 ) {
    out->style = stringalloc(style);
  } else if ( strcmp(style,"local") == 0 ) {
    out->style = stringalloc(style);
  } else if ( strcmp(style,"ramp") == 0 ) {
    out->style = stringalloc(style);
  } else {
    warn("Unrecognised start/end style");
    return NULL;
  }
    
  out->i_ramp = 5;
  out->j_ramp = 5;

  strip_out_integer_argument(argc,argv,"se_query_ramp",&out->i_ramp);
  strip_out_integer_argument(argc,argv,"se_query_ramp",&out->j_ramp);

  return(out);

}


# line 72 "startendmodel.dy"
StartEndScore * StartEndScore_from_StartEndModel(StartEndModel * sem)
{
  StartEndScore * out;

  int i;

  assert(sem != NULL);

  out = StartEndScore_alloc();

  out->i_entry = calloc(sem->i_len,sizeof(Score));
  out->j_entry = calloc(sem->j_len,sizeof(Score));

  out->i_len = sem->i_len;
  out->j_len = sem->j_len;

  for(i=0;i<out->i_len;i++) {
    out->i_entry[i] = Probability2Score(sem->i_entry[i]);
  }

  for(i=0;i<out->j_len;i++) {
    out->j_entry[i] = Probability2Score(sem->j_entry[i]);
  }

  return(out);

}




# line 82 "startendmodel.c"
/* Function:  hard_link_StartEndModelMode(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StartEndModelMode *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModelMode *]
 *
 */
StartEndModelMode * hard_link_StartEndModelMode(StartEndModelMode * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StartEndModelMode object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StartEndModelMode_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StartEndModelMode *]
 *
 */
StartEndModelMode * StartEndModelMode_alloc(void) 
{
    StartEndModelMode * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StartEndModelMode *) ckalloc (sizeof(StartEndModelMode))) == NULL)  {  
      warn("StartEndModelMode_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->style = NULL;   
    out->i_ramp = 0; 
    out->j_ramp = 0; 


    return out;  
}    


/* Function:  free_StartEndModelMode(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StartEndModelMode *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModelMode *]
 *
 */
StartEndModelMode * free_StartEndModelMode(StartEndModelMode * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StartEndModelMode obj. Should be trappable"); 
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
    if( obj->style != NULL)  
      ckfree(obj->style);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_StartEndModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StartEndModel *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModel *]
 *
 */
StartEndModel * hard_link_StartEndModel(StartEndModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StartEndModel object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StartEndModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StartEndModel *]
 *
 */
StartEndModel * StartEndModel_alloc(void) 
{
    StartEndModel * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StartEndModel *) ckalloc (sizeof(StartEndModel))) == NULL)  {  
      warn("StartEndModel_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->i_len = 0;  
    out->j_len = 0;  
    out->i_entry = NULL; 
    out->j_entry = NULL; 


    return out;  
}    


/* Function:  free_StartEndModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StartEndModel *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndModel *]
 *
 */
StartEndModel * free_StartEndModel(StartEndModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StartEndModel obj. Should be trappable"); 
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
    if( obj->i_entry != NULL)    
      free_Probability(obj->i_entry);    
    if( obj->j_entry != NULL)    
      free_Probability(obj->j_entry);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_StartEndScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [StartEndScore *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndScore *]
 *
 */
StartEndScore * hard_link_StartEndScore(StartEndScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a StartEndScore object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  StartEndScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [StartEndScore *]
 *
 */
StartEndScore * StartEndScore_alloc(void) 
{
    StartEndScore * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(StartEndScore *) ckalloc (sizeof(StartEndScore))) == NULL)  {  
      warn("StartEndScore_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->i_len = 0;  
    out->j_len = 0;  
    out->i_entry = NULL; 
    out->j_entry = NULL; 


    return out;  
}    


/* Function:  free_StartEndScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [StartEndScore *]
 *
 * Return [UNKN ]  Undocumented return value [StartEndScore *]
 *
 */
StartEndScore * free_StartEndScore(StartEndScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a StartEndScore obj. Should be trappable"); 
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
    if( obj->i_entry != NULL)    
      free_Score(obj->i_entry);  
    if( obj->j_entry != NULL)    
      free_Score(obj->j_entry);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
