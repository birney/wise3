#ifdef _cplusplus
extern "C" {
#endif
#include "motifmatrixdp.h"


# line 6 "motifmatrixdp.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Mon Jun 19 12:42:33 2006           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define SPACER_MATCH 0   
#define MOTIF_MATCH 1    
#define CONS_MATCH 2 


#define START 0  
#define END 1    


#define LocalMotifMatrix_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*3)+STATE][i+1]    
#define LocalMotifMatrix_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]   
#define LocalMotifMatrix_READ_OFF_ERROR -3
  


#define LocalMotifMatrix_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+1)*3)+STATE]  
#define LocalMotifMatrix_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE]   




#define LocalMotifMatrix_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j]    
#define LocalMotifMatrix_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE)    


/* Function:  PackAln_read_Shatter_LocalMotifMatrix(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    LocalMotifMatrix_access_func_holder holder;  


    holder.access_main    = LocalMotifMatrix_shatter_access_main;    
    holder.access_special = LocalMotifMatrix_shatter_access_special; 
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_LocalMotifMatrix(mat,holder);    
}    


/* Function:  LocalMotifMatrix_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalMotifMatrix_shatter_access_main(LocalMotifMatrix * mat,int i,int j,int state) 
{
    return LocalMotifMatrix_SHATTER_MATRIX(mat,i,j,state);   
}    


/* Function:  LocalMotifMatrix_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalMotifMatrix_shatter_access_special(LocalMotifMatrix * mat,int i,int j,int state) 
{
    return LocalMotifMatrix_SHATTER_SPECIAL(mat,i,j,state);  
}    


/* Function:  calculate_shatter_LocalMotifMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalMotifMatrix matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [LocalMotifMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_LocalMotifMatrix(LocalMotifMatrix * mat,DPEnvelope * dpenv) 
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
    int * SIG_1_1;   
    int * SIG_1_0;   
    int * SIG_0_1;   


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


    start_reporting("LocalMotifMatrix Matrix calculation: ");    
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
        SIG_1_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-1);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   




        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_1[SPACER_MATCH] + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);   
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = SIG_1_0[SPACER_MATCH] + mat->ms->spacer_indel;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = SIG_0_1[SPACER_MATCH] + mat->ms->spacer_indel;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = SIG_1_1[MOTIF_MATCH] + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = SIG_1_1[CONS_MATCH] + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_SHATTER_SPECIAL(mat,i-1,j-1,START) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[SPACER_MATCH] = score;  


        /* state SPACER_MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalMotifMatrix_SHATTER_SPECIAL(mat,i,j,END) )   {  
          LocalMotifMatrix_SHATTER_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MOTIF_MATCH] + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);   
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = SIG_1_0[MOTIF_MATCH] + mat->ms->motif_indel;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = SIG_0_1[MOTIF_MATCH] + mat->ms->motif_indel;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = SIG_1_1[SPACER_MATCH] + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[MOTIF_MATCH] = score;   


        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_1[CONS_MATCH] + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);    
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = SIG_1_0[CONS_MATCH] + mat->ms->cons_indel;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = SIG_0_1[CONS_MATCH] + mat->ms->cons_indel;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = SIG_1_1[SPACER_MATCH] + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[CONS_MATCH] = score;    


        /* Finished calculating state CONS_MATCH */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_LocalMotifMatrix(dbsi,out,query,target,motif,ms)
 *
 * Descrip:    This function makes a database search of LocalMotifMatrix
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:          dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         motif [UNKN ] Undocumented argument [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Undocumented argument [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_LocalMotifMatrix(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms) 
{
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_LocalMotifMatrix. Can't process results!");   
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_LocalMotifMatrix. Can't process results!"); 
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for LocalMotifMatrix, but it was compiled with a trace level of -2139062144. Not all trace statements can be shown",dbsi->trace_level);  
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_LocalMotifMatrix(out,query,target ,motif,ms);   
      case DBSearchImpl_Pthreads :   
        warn("This matrix LocalMotifMatrix was not dyc compiled with thread support");   
        return SEARCH_ERROR; 
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from LocalMotifMatrix",impl_string_DBSearchImpl(dbsi));   
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  serial_search_LocalMotifMatrix(out,query,target,motif,ms)
 *
 * Descrip:    This function makes a database search of LocalMotifMatrix
 *             It is a single processor implementation
 *
 *
 * Arg:           out [UNKN ] Undocumented argument [Hscore *]
 * Arg:         query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:        target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         motif [UNKN ] Undocumented argument [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Undocumented argument [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_LocalMotifMatrix(Hscore * out,ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_LocalMotifMatrix(query, target , motif, ms);  
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("LocalMotifMatrix search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");   
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


/* Function:  score_only_LocalMotifMatrix(query,target,motif,ms)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_LocalMotifMatrix_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    LocalMotifMatrix * mat;  


    mat = allocate_LocalMotifMatrix_only(query, target , motif, ms); 
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 1) * 3,2,2)) == NULL)  {  
      warn("Score only matrix for LocalMotifMatrix cannot be allocated, (asking for 1  by %d  cells)",mat->leni*3);  
      mat = free_LocalMotifMatrix(mat);  
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<3;k++) 
          LocalMotifMatrix_VSMALL_MATRIX(mat,i,j,k) = NEGI;  
        }  
      LocalMotifMatrix_VSMALL_SPECIAL(mat,i,j,START) = 0;    
      LocalMotifMatrix_VSMALL_SPECIAL(mat,i,j,END) = NEGI;   
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);    
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_SPECIAL(mat,i-1,j-1,START) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_VSMALL_MATRIX(mat,i,j,SPACER_MATCH) = score;   


        /* state SPACER_MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalMotifMatrix_VSMALL_SPECIAL(mat,i,j,END) )    {  
          LocalMotifMatrix_VSMALL_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);    
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_VSMALL_MATRIX(mat,i,j,MOTIF_MATCH) = score;    


        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);     
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_VSMALL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_VSMALL_MATRIX(mat,i,j,CONS_MATCH) = score; 


        /* Finished calculating state CONS_MATCH */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      if( bestscore < LocalMotifMatrix_VSMALL_SPECIAL(mat,0,j,END) ) 
        bestscore = LocalMotifMatrix_VSMALL_SPECIAL(mat,0,j,END);    
      } /* end of for all target positions */ 


    mat = free_LocalMotifMatrix(mat);    
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_LocalMotifMatrix(query,target,motif,ms,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_LocalMotifMatrix
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    LocalMotifMatrix * mat;  
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->seq->len * target->seq->len;  
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
        if( (mat=allocate_Expl_LocalMotifMatrix(query, target , motif, ms,dpri)) == NULL )   {  
          warn("Unable to allocate large LocalMotifMatrix version"); 
          return NULL;   
          }  
        calculate_dpenv_LocalMotifMatrix(mat,dpenv);     
        out =  PackAln_read_Expl_LocalMotifMatrix(mat);  
        }  
      else   {  
        mat = allocate_LocalMotifMatrix_only(query, target , motif, ms);     
        calculate_shatter_LocalMotifMatrix(mat,dpenv);   
        out = PackAln_read_Shatter_LocalMotifMatrix(mat);    
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_LocalMotifMatrix(query, target , motif, ms)) == NULL )   {  
          warn("Unable to allocate small LocalMotifMatrix version"); 
          return NULL;   
          }  
        out = PackAln_calculate_Small_LocalMotifMatrix(mat,dpenv);   
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_LocalMotifMatrix(query, target , motif, ms,dpri)) == NULL )   {  
          warn("Unable to allocate large LocalMotifMatrix version"); 
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_LocalMotifMatrix(mat);   
          out =  PackAln_read_Expl_LocalMotifMatrix(mat);    
          }  
        }  
      }  


    mat = free_LocalMotifMatrix(mat);    
    return out;  
}    


/* Function:  allocate_LocalMotifMatrix_only(query,target,motif,ms)
 *
 * Descrip:    This function only allocates the LocalMotifMatrix structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * allocate_LocalMotifMatrix_only(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms) 
{
    LocalMotifMatrix * out;  


    if((out= LocalMotifMatrix_alloc()) == NULL)  {  
      warn("Allocation of basic LocalMotifMatrix structure failed...");  
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->motif = motif;  
    out->ms = ms;    
    out->leni = query->seq->len;     
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_LocalMotifMatrix(query,target,motif,ms,dpri)
 *
 * Descrip:    This function allocates the LocalMotifMatrix structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_LocalMotifMatrix_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 * Arg:          dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * allocate_Expl_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms,DPRunImpl * dpri) 
{
    LocalMotifMatrix * out;  


    out = allocate_LocalMotifMatrix_only(query, target , motif, ms); 
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*3 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*3,(out->leni+1),2,out->lenj+1)) == NULL) {  
        warn("Explicit matrix LocalMotifMatrix cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);  
        free_LocalMotifMatrix(out);  
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_LocalMotifMatrix(out);  
    return out;  
}    


/* Function:  init_LocalMotifMatrix(mat)
 *
 * Descrip:    This function initates LocalMotifMatrix matrix when in explicit mode
 *             Called in /allocate_Expl_LocalMotifMatrix
 *
 *
 * Arg:        mat [UNKN ] LocalMotifMatrix which contains explicit basematrix memory [LocalMotifMatrix *]
 *
 */
void init_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->seq->len;i++)  {  
      for(j= (-1);j<2;j++)   {  
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;   
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;    
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH) = NEGI; 
        }  
      }  
    for(j= (-1);j<mat->target->seq->len;j++) {  
      for(i= (-1);i<2;i++)   {  
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;   
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;    
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH) = NEGI; 
        }  
      LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,START) = 0;  
      LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  
    return;  
}    


/* Function:  recalculate_PackAln_LocalMotifMatrix(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by LocalMotifMatrix
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 */
void recalculate_PackAln_LocalMotifMatrix(PackAln * pal,LocalMotifMatrix * mat) 
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
        case SPACER_MATCH :  
          if( offi == 1 && offj == 1 && prev->state == SPACER_MATCH )    {  
            pau->score = (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration) + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == SPACER_MATCH )    {  
            pau->score = mat->ms->spacer_indel + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == SPACER_MATCH )    {  
            pau->score = mat->ms->spacer_indel + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == MOTIF_MATCH ) {  
            pau->score = mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)] + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == CONS_MATCH )  {  
            pau->score = mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)] + (0);   
            continue;    
            }  
          if( offj == 1 && prev->state == (START+3) )    {  
            pau->score = mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)] + (0);   
            continue;    
            }  
          warn("In recaluclating PackAln with state SPACER_MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);   
          break; 
        case MOTIF_MATCH :   
          if( offi == 1 && offj == 1 && prev->state == MOTIF_MATCH ) {  
            pau->score = (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == MOTIF_MATCH ) {  
            pau->score = mat->ms->motif_indel + (0);     
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == MOTIF_MATCH ) {  
            pau->score = mat->ms->motif_indel + (0);     
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == SPACER_MATCH )    {  
            pau->score = (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif) + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state MOTIF_MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case CONS_MATCH :    
          if( offi == 1 && offj == 1 && prev->state == CONS_MATCH )  {  
            pau->score = (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration) + (0);   
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == CONS_MATCH )  {  
            pau->score = mat->ms->cons_indel + (0);  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == CONS_MATCH )  {  
            pau->score = mat->ms->cons_indel + (0);  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == SPACER_MATCH )    {  
            pau->score = (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons) + (0);  
            continue;    
            }  
          warn("In recaluclating PackAln with state CONS_MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case (START+3) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (END+3) :   
          if( offj == 0 && prev->state == SPACER_MATCH ) {  
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
#define LocalMotifMatrix_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+1)*3+state])  
#define LocalMotifMatrix_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+1)*3+state])  
#define LocalMotifMatrix_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)]) 
#define LocalMotifMatrix_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])    
#define LocalMotifMatrix_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+1)*3 + state])   
#define LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])    
#define LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 3) + ((i+1) * 3) + (state)])    
#define LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 24) + ((i+1) * 24) + (state * 8) + shadow+1])    
#define LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)])    
/* Function:  allocate_Small_LocalMotifMatrix(query,target,motif,ms)
 *
 * Descrip:    This function allocates the LocalMotifMatrix structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_LocalMotifMatrix_only
 *
 *
 * Arg:         query [UNKN ] query data structure [ComplexSequence*]
 * Arg:        target [UNKN ] target data structure [ComplexSequence*]
 * Arg:         motif [UNKN ] Resource [MotifConsMatrix*]
 * Arg:            ms [UNKN ] Resource [MotifMatrixScore*]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
#define LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)])    
LocalMotifMatrix * allocate_Small_LocalMotifMatrix(ComplexSequence* query,ComplexSequence* target ,MotifConsMatrix* motif,MotifMatrixScore* ms) 
{
    LocalMotifMatrix * out;  


    out = allocate_LocalMotifMatrix_only(query, target , motif, ms); 
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 1) * 3,16,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix LocalMotifMatrix cannot be allocated, (asking for 2 by %d main cells)",out->leni+2); 
      free_LocalMotifMatrix(out);    
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_LocalMotifMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for LocalMotifMatrix structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_LocalMotifMatrix 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_LocalMotifMatrix 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_LocalMotifMatrix(LocalMotifMatrix * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for LocalMotifMatrix due to wrong type of matrix");    
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_LocalMotifMatrix(mat,dpenv); 
    score = start_end_find_end_LocalMotifMatrix(mat,&endj);  
    out->score = score;  
    stopstate = END;
    
    /* No special to specials: one matrix alignment: simply remove and get */ 
    starti = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,0);    
    startj = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,1);    
    startstate = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,2);    
    stopi = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,3); 
    stopj = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,4); 
    stopstate = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,5); 
    temp = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,endj,END,6);  
    log_full_error(REPORT,0,"[%d,%d][%d,%d] Score %d",starti,startj,stopi,stopj,score);  
    stop_reporting();    
    start_reporting("Recovering alignment: ");   


    /* Figuring how much j we have to align for reporting purposes */ 
    donej = 0;   
    totalj = stopj - startj; 
    full_dc_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv);    


    /* Although we have no specials, need to get start. Better to check than assume */ 


    max_matrix_to_special_LocalMotifMatrix(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);   
    if( stopi == LocalMotifMatrix_READ_OFF_ERROR || stopstate != START ) {  
      warn("Problem in reading off special state system, hit a non start state (or an internal error) in a single alignment mode");  
      invert_PackAln(out);   
      recalculate_PackAln_LocalMotifMatrix(out,mat); 
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
    recalculate_PackAln_LocalMotifMatrix(out,mat);   
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_LocalMotifMatrix(mat)
 *
 * Descrip:    This function calculates an alignment for LocalMotifMatrix structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_LocalMotifMatrix 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_LocalMotifMatrix
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_LocalMotifMatrix 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_LocalMotifMatrix(mat,NULL);  
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_LocalMotifMatrix(mat);    
    return out;  
}    


/* Function:  AlnRangeSet_from_LocalMotifMatrix(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for LocalMotifMatrix structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_LocalMotifMatrix 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_LocalMotifMatrix
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_LocalMotifMatrix");  
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_LocalMotifMatrix(mat,&jpos); 
    state = END; 


    while( (temp = AlnRange_build_LocalMotifMatrix(mat,jpos,state,&jpos,&state)) != NULL)    
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_LocalMotifMatrix(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_LocalMotifMatrix
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_LocalMotifMatrix(LocalMotifMatrix * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_LocalMotifMatrix");  
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_LocalMotifMatrix(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)   {  
      warn("In AlnRanger_build_LocalMotifMatrix alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);   
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0); 
    out->startj = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1); 
    out->startstate = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2); 
    out->stopi = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);  
    out->stopj = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);  
    out->stopstate = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);  
    out->startscore = LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6); 
    out->stopscore = LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,jpos,state);   


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_LocalMotifMatrix(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL); 
    if( jpos == LocalMotifMatrix_READ_OFF_ERROR) {  
      warn("In AlnRange_build_LocalMotifMatrix alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj);    
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate); 
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
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
boolean read_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_LocalMotifMatrix(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    


      if( i == LocalMotifMatrix_READ_OFF_ERROR)  {  
        warn("In LocalMotifMatrix hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);  
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
        warn("In LocalMotifMatrix hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj);    
        return FALSE;    
        }  
      }  
    warn("In LocalMotifMatrix hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state); 
    return FALSE;    
}    


/* Function:  max_hidden_LocalMotifMatrix(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
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
int max_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = LocalMotifMatrix_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->query->seq->len || j > mat->target->seq->len) {  
      warn("In LocalMotifMatrix matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,state);  
    switch(state)    { /*Switch state */ 
      case SPACER_MATCH :    
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]) -  (0); 
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,CONS_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,CONS_MATCH);    
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,CONS_MATCH);     
          }  
        temp = cscore - (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]) -  (0); 
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,MOTIF_MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,MOTIF_MATCH);   
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,MOTIF_MATCH);    
          }  
        temp = cscore - (mat->ms->spacer_indel) -  (0);  
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_MATCH);  
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 0,j - 1,SPACER_MATCH);   
          }  
        temp = cscore - (mat->ms->spacer_indel) -  (0);  
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 0,SPACER_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-0,SPACER_MATCH);  
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 0,SPACER_MATCH);   
          }  
        temp = cscore - ((mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration)) -  (0);  
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,SPACER_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,SPACER_MATCH);  
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,SPACER_MATCH);   
          }  
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case MOTIF_MATCH :     
        temp = cscore - ((MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif)) -  (0);    
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,SPACER_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,SPACER_MATCH);  
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,SPACER_MATCH);   
          }  
        temp = cscore - (mat->ms->motif_indel) -  (0);   
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 0,j - 1,MOTIF_MATCH) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-0,j-1,MOTIF_MATCH);   
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 0,j - 1,MOTIF_MATCH);    
          }  
        temp = cscore - (mat->ms->motif_indel) -  (0);   
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 0,MOTIF_MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-0,MOTIF_MATCH);   
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 0,MOTIF_MATCH);    
          }  
        temp = cscore - ((MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration)) -  (0); 
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,MOTIF_MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,MOTIF_MATCH);   
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,MOTIF_MATCH);    
          }  
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CONS_MATCH :  
        temp = cscore - ((mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons)) -  (0);    
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,SPACER_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,SPACER_MATCH);  
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,SPACER_MATCH);   
          }  
        temp = cscore - (mat->ms->cons_indel) -  (0);    
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 0,j - 1,CONS_MATCH) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-0,j-1,CONS_MATCH);    
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 0,j - 1,CONS_MATCH);     
          }  
        temp = cscore - (mat->ms->cons_indel) -  (0);    
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 0,CONS_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-0,CONS_MATCH);    
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 0,CONS_MATCH);     
          }  
        temp = cscore - ((mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration)) -  (0); 
        if( temp == LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,CONS_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,CONS_MATCH);    
            }  
          return LocalMotifMatrix_HIDDEN_MATRIX(mat,i - 1,j - 1,CONS_MATCH);     
          }  
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_LocalMotifMatrix(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
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
boolean read_special_strip_LocalMotifMatrix(LocalMotifMatrix * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START) { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 3; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_LocalMotifMatrix(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore);    
      if( i == LocalMotifMatrix_READ_OFF_ERROR)  {  
        warn("In special strip read LocalMotifMatrix, got a bad read off error. Sorry!");    
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4))    {  
      warn("In special strip read LocalMotifMatrix, at special [%d] state [%d] overshot!",j,state);  
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


/* Function:  max_special_strip_LocalMotifMatrix(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
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
int max_special_strip_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = LocalMotifMatrix_READ_OFF_ERROR; 
    if( isspecial == FALSE ) {  
      warn("In special strip max function for LocalMotifMatrix, got a non special start point. Problem! (bad!)");    
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In LocalMotifMatrix matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1; 
      }  


    cscore = LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,i,j,state);  
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* Source SPACER_MATCH is not a special */ 
      default:   
        warn("Major problem (!) - in LocalMotifMatrix special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);   
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_LocalMotifMatrix(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
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
int max_matrix_to_special_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = LocalMotifMatrix_READ_OFF_ERROR; 


    if( j < 0 || j > mat->lenj)  {  
      warn("In LocalMotifMatrix matrix to special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case SPACER_MATCH :    
        temp = cscore - (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]) -  (0);     
        if( temp == LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);     
            }  
          return LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;  
          }  
        /* Source CONS_MATCH is not a special, should not get here! */ 
        /* Source MOTIF_MATCH is not a special, should not get here! */ 
        /* Source SPACER_MATCH is not a special, should not get here! */ 
        /* Source SPACER_MATCH is not a special, should not get here! */ 
        /* Source SPACER_MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalMotifMatrix matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case MOTIF_MATCH :     
        /* Source SPACER_MATCH is not a special, should not get here! */ 
        /* Source MOTIF_MATCH is not a special, should not get here! */ 
        /* Source MOTIF_MATCH is not a special, should not get here! */ 
        /* Source MOTIF_MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalMotifMatrix matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      case CONS_MATCH :  
        /* Source SPACER_MATCH is not a special, should not get here! */ 
        /* Source CONS_MATCH is not a special, should not get here! */ 
        /* Source CONS_MATCH is not a special, should not get here! */ 
        /* Source CONS_MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in LocalMotifMatrix matrix to special read off, position %d,%d state %d no source found!",i,j,state);  
        return (-1); 
      default:   
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_LocalMotifMatrix(mat,starti,startj,stopi,stopj);     


    LocalMotifMatrix_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0;    


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;   
          LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;    
          LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,CONS_MATCH) = NEGI;     
          continue;  
          } /* end of Is not in envelope */ 


        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);    
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,SPACER_MATCH) = score;   
        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);    
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,MOTIF_MATCH) = score;    
        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);     
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_HIDDEN_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,CONS_MATCH) = score; 
        /* Finished calculating state CONS_MATCH */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_LocalMotifMatrix(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;
    
        LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;
 
        LocalMotifMatrix_HIDDEN_MATRIX(mat,i,j,CONS_MATCH) = NEGI;
  
        }  
      }  


    return;  
}    


/* Function:  full_dc_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_LocalMotifMatrix
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_LocalMotifMatrix to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [LocalMotifMatrix *]
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
boolean full_dc_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_LocalMotifMatrix");    
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);   
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE) {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)  {  
      warn("In divide and conquor for LocalMotifMatrix, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj); 
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);  
    lstartj= LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);  
    lstate = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);  


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_LocalMotifMatrix(mat,LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_LocalMotifMatrix(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_LocalMotifMatrix(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
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
boolean do_dc_single_pass_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_LocalMotifMatrix(mat);   


    LocalMotifMatrix_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0; 
    run_up_dc_LocalMotifMatrix(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);     
    push_dc_at_merge_LocalMotifMatrix(mat,starti,stopi,halfj,&halfj,dpenv);  
    follow_on_dc_LocalMotifMatrix(mat,starti,stopi,halfj,stopj,dpenv,perc_done);     
    return TRUE; 
}    


/* Function:  push_dc_at_merge_LocalMotifMatrix(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int k;  
    register int count;  
    register int mergej;/* Sources below this j will be stamped by triples */ 
    register int score;  
    register int temp;   


    mergej = startj -1;  
    for(count=0,j=startj;count<1;count++,j++)    {  
      for(i=starti;i<=stopi;i++) {  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;    
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = (-100); 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,1) = (-100); 
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;     
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,0) = (-100);  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,1) = (-100);  
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = NEGI;  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,0) = (-100);   
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,1) = (-100);   
          continue;  
          } /* end of Is not in envelope */ 


        /* For state SPACER_MATCH, pushing when j - offj <= mergej */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);     
        if( j - 1 <= mergej) {  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = i-1;    
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,1) = j-1;    
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,2) = SPACER_MATCH;   
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,3) = i;  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,4) = j;  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,5) = SPACER_MATCH;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,k); 
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = i-1;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,1) = j-0;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,2) = SPACER_MATCH; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,3) = i;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,4) = j;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,5) = SPACER_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SPACER_MATCH,k);   
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = i-0;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,1) = j-1;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,2) = SPACER_MATCH; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,3) = i;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,4) = j;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,5) = SPACER_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_MATCH,k);   
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = i-1;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,1) = j-1;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,2) = MOTIF_MATCH;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,3) = i;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,4) = j;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,5) = SPACER_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MOTIF_MATCH,k);    
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = i-1;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,1) = j-1;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,2) = CONS_MATCH;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,3) = i;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,4) = j;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,5) = SPACER_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,CONS_MATCH,k); 
            }  
          }  
        /* Add any movement independant score */ 
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = score;     
        /* Finished with state SPACER_MATCH */ 


        /* For state MOTIF_MATCH, pushing when j - offj <= mergej */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);     
        if( j - 1 <= mergej) {  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,0) = i-1; 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,1) = j-1; 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,2) = MOTIF_MATCH; 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,3) = i;   
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,4) = j;   
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,5) = MOTIF_MATCH; 
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MOTIF_MATCH,k);   
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;    
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,0) = i-1;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,1) = j-0;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,2) = MOTIF_MATCH;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,3) = i; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,4) = j; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,5) = MOTIF_MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MOTIF_MATCH,k); 
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,0) = i-0;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,1) = j-1;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,2) = MOTIF_MATCH;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,3) = i; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,4) = j; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,5) = MOTIF_MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MOTIF_MATCH,k); 
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,0) = i-1;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,1) = j-1;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,2) = SPACER_MATCH;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,3) = i; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,4) = j; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,5) = MOTIF_MATCH;   
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,k);    
            }  
          }  
        /* Add any movement independant score */ 
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = score;  
        /* Finished with state MOTIF_MATCH */ 


        /* For state CONS_MATCH, pushing when j - offj <= mergej */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);  
        if( j - 1 <= mergej) {  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,0) = i-1;  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,1) = j-1;  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,2) = CONS_MATCH;   
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,3) = i;    
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,4) = j;    
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,5) = CONS_MATCH;   
          }  
        else {  
          for(k=0;k<7;k++)   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,CONS_MATCH,k); 
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;  
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,0) = i-1;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,1) = j-0;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,2) = CONS_MATCH; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,3) = i;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,4) = j;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,5) = CONS_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,CONS_MATCH,k);   
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,0) = i-0;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,1) = j-1;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,2) = CONS_MATCH; 
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,3) = i;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,4) = j;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,5) = CONS_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CONS_MATCH,k);   
            }  
          }  


        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,0) = i-1;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,1) = j-1;    
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,2) = SPACER_MATCH;   
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,3) = i;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,4) = j;  
            LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,5) = CONS_MATCH; 
            }  
          else   {  
            for(k=0;k<7;k++) 
              LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,k); 
            }  
          }  
        /* Add any movement independant score */ 
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = score;   
        /* Finished with state CONS_MATCH */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_LocalMotifMatrix(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;    
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;     
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,k); 
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SPACER_MATCH,k);   
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_MATCH,k);   
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MOTIF_MATCH,k);    
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,CONS_MATCH,k); 
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = score;    
        for(k=0;k<7;k++) 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MOTIF_MATCH,k);  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MOTIF_MATCH,k);    
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MOTIF_MATCH,k);    
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,k);   
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = score; 
        for(k=0;k<7;k++) 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = localshadow[k];  
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);  
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,CONS_MATCH,k);   
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,CONS_MATCH,k); 
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CONS_MATCH,k); 
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,k);   
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = score;  
        for(k=0;k<7;k++) 
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state CONS_MATCH */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_LocalMotifMatrix(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_LocalMotifMatrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;    
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;     
          LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);     
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = score;    
        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);     
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = score; 
        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = score;  
        /* Finished calculating state CONS_MATCH */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_LocalMotifMatrix(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 */
}    
void init_dc_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->query->seq->len;i++) {  
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;  
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;   
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = NEGI;    
        for(k=0;k<7;k++) {  
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = (-1);   
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = (-1);    
          LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = (-1); 
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_LocalMotifMatrix(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [LocalMotifMatrix *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_LocalMotifMatrix(LocalMotifMatrix * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END); 
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,j,END) > max )    {  
        max = LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,j,END);   
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_LocalMotifMatrix(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [LocalMotifMatrix]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_LocalMotifMatrix(LocalMotifMatrix *mat,DPEnvelope * dpenv) 
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
    leni = mat->query->seq->len; 
    lenj = mat->target->seq->len;    
    total = leni * lenj; 


    score_pointers = (int *) calloc (1 * (leni + 1) * 3,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 1) * 3 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;    
          LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;     
          LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,0));   
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,SPACER_MATCH,0)); 
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,SPACER_MATCH,0)); 
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)] +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MOTIF_MATCH,0));  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)] +(0);  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,CONS_MATCH,0));   
          }  
        /* From state START to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)] + (0);     
        if( temp  > score )  {  
          score = temp;  
          /* This state [START] is a special for SPACER_MATCH... push top shadow pointers here */ 
          localshadow[0]= i; 
          localshadow[1]= j; 
          localshadow[2]= SPACER_MATCH;  
          localshadow[3]= (-1);  
          localshadow[4]= (-1);  
          localshadow[5]= (-1);  
          localshadow[6]= score; 
          localsp = localshadow; 
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = score;    
        for(k=0;k<7;k++) 
          LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state SPACER_MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )     {  
          LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,k); 
          LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,6);   
          LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;  
          LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;  
          LocalMotifMatrix_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = SPACER_MATCH;   
          }  




        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration) + (0);   
        /* assign local shadown pointer */ 
        localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MOTIF_MATCH,0));    
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MOTIF_MATCH,0));  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MOTIF_MATCH,0));  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,0)); 
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = score; 
        for(k=0;k<7;k++) 
          LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,k) = localsp[k];  
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration) + (0);    
        /* assign local shadown pointer */ 
        localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,CONS_MATCH,0)); 
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,CONS_MATCH,0));   
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel +(0);     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,CONS_MATCH,0));   
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons) +(0);   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,SPACER_MATCH,0)); 
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         LocalMotifMatrix_DC_OPT_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = score;  
        for(k=0;k<7;k++) 
          LocalMotifMatrix_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state CONS_MATCH */ 


        } /* end of for each i position in strip */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_LocalMotifMatrix(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 */
void init_start_end_linear_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->query->seq->len;i++) {  
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;  
        LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,SPACER_MATCH,0) = (-1); 
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;   
        LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,MOTIF_MATCH,0) = (-1);  
        LocalMotifMatrix_DC_SHADOW_MATRIX(mat,i,j,CONS_MATCH) = NEGI;    
        LocalMotifMatrix_DC_SHADOW_MATRIX_SP(mat,i,j,CONS_MATCH,0) = (-1);   
        }  
      }  


    for(j=(-1);j<mat->target->seq->len;j++)  {  
      LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,j,START) = 0; 
      LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j;    
      LocalMotifMatrix_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI;    
      LocalMotifMatrix_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);   
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_LocalMotifMatrix(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_LocalMotifMatrix(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_LocalMotifMatrix();  
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "SPACER_MATCH","SPACER_INSERT","MOTIF_MATCH","MOTIF_INSERT","CONS_MATCH","CONS_INSERT","END" }; 
/* Function:  AlnConvertSet_LocalMotifMatrix(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "SPACER_MATCH","SPACER_INSERT","MOTIF_MATCH","MOTIF_INSERT","CONS_MATCH","CONS_INSERT","END" };    
AlnConvertSet * AlnConvertSet_LocalMotifMatrix(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_MATCH;  
    acu->state2 = SPACER_MATCH;  
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_MATCH;  
    acu->state2 = SPACER_MATCH;  
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_MATCH;  
    acu->state2 = SPACER_MATCH;  
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MOTIF_MATCH;   
    acu->state2 = SPACER_MATCH;  
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CONS_MATCH;    
    acu->state2 = SPACER_MATCH;  
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 3; 
    acu->is_from_special = TRUE; 
    acu->state2 = SPACER_MATCH;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MOTIF_MATCH;   
    acu->state2 = MOTIF_MATCH;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MOTIF_MATCH;   
    acu->state2 = MOTIF_MATCH;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[3];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MOTIF_MATCH;   
    acu->state2 = MOTIF_MATCH;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_MATCH;  
    acu->state2 = MOTIF_MATCH;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CONS_MATCH;    
    acu->state2 = CONS_MATCH;    
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CONS_MATCH;    
    acu->state2 = CONS_MATCH;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[5];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = CONS_MATCH;    
    acu->state2 = CONS_MATCH;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[5];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_MATCH;  
    acu->state2 = CONS_MATCH;    
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[4];    
    acu->label2 = target_label[4];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = SPACER_MATCH;  
    acu->state2 = END + 3;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[6];    
    acu->label2 = target_label[6];   
    return out;  
}    


/* Function:  PackAln_read_Expl_LocalMotifMatrix(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    LocalMotifMatrix_access_func_holder holder;  


    holder.access_main    = LocalMotifMatrix_explicit_access_main;   
    holder.access_special = LocalMotifMatrix_explicit_access_special;    
    return PackAln_read_generic_LocalMotifMatrix(mat,holder);    
}    


/* Function:  LocalMotifMatrix_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalMotifMatrix_explicit_access_main(LocalMotifMatrix * mat,int i,int j,int state) 
{
    return LocalMotifMatrix_EXPL_MATRIX(mat,i,j,state);  
}    


/* Function:  LocalMotifMatrix_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int LocalMotifMatrix_explicit_access_special(LocalMotifMatrix * mat,int i,int j,int state) 
{
    return LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,state); 
}    


/* Function:  PackAln_read_generic_LocalMotifMatrix(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:          h [UNKN ] Undocumented argument [LocalMotifMatrix_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_LocalMotifMatrix(LocalMotifMatrix * mat,LocalMotifMatrix_access_func_holder h) 
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


    out->score =  find_end_LocalMotifMatrix(mat,&i,&j,&state,&isspecial,h);  


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in LocalMotifMatrix_basic_read, returning a mess.(Sorry!)",out->len); 
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
        max_calc_special_LocalMotifMatrix(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      else   
        max_calc_LocalMotifMatrix(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);     
      if(i == LocalMotifMatrix_READ_OFF_ERROR || j == LocalMotifMatrix_READ_OFF_ERROR || state == LocalMotifMatrix_READ_OFF_ERROR )  {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in LocalMotifMatrix_basic_read, returning partial alignment",out->len); 
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


/* Function:  find_end_LocalMotifMatrix(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [LocalMotifMatrix_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_LocalMotifMatrix(LocalMotifMatrix * mat,int * ri,int * rj,int * state,boolean * isspecial,LocalMotifMatrix_access_func_holder h) 
{
    int j;   
    int max; 
    int maxj;    
    int temp;    


    max = (*h.access_special)(mat,0,mat->target->seq->len-1,END);    
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
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


/* Function:  LocalMotifMatrix_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void LocalMotifMatrix_debug_show_matrix(LocalMotifMatrix * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->seq->len;i++)    {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State SPACER_MATCH %d\n",LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH));   
        fprintf(ofp,"State MOTIF_MATCH %d\n",LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH)); 
        fprintf(ofp,"State CONS_MATCH %d\n",LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH));   
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_LocalMotifMatrix(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [LocalMotifMatrix_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalMotifMatrix_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = LocalMotifMatrix_READ_OFF_ERROR; 


    if( i < 0 || j < 0 || i > mat->query->seq->len || j > mat->target->seq->len) {  
      warn("In LocalMotifMatrix matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state);    
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case SPACER_MATCH :    
        temp = cscore - (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]) -  (0); 
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
        temp = cscore - (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,CONS_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,CONS_MATCH);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,CONS_MATCH);   
          }  
        temp = cscore - (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MOTIF_MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MOTIF_MATCH); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MOTIF_MATCH);  
          }  
        temp = cscore - (mat->ms->spacer_indel) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,SPACER_MATCH) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,SPACER_MATCH);    
            }  
          return (*h.access_main)(mat,i - 0,j - 1,SPACER_MATCH);     
          }  
        temp = cscore - (mat->ms->spacer_indel) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 0,SPACER_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,SPACER_MATCH);    
            }  
          return (*h.access_main)(mat,i - 1,j - 0,SPACER_MATCH);     
          }  
        temp = cscore - ((mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration)) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 1,SPACER_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,SPACER_MATCH);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,SPACER_MATCH);     
          }  
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case MOTIF_MATCH :     
        temp = cscore - ((MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,SPACER_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,SPACER_MATCH);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,SPACER_MATCH);     
          }  
        temp = cscore - (mat->ms->motif_indel) -  (0);   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,MOTIF_MATCH) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,MOTIF_MATCH); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,MOTIF_MATCH);  
          }  
        temp = cscore - (mat->ms->motif_indel) -  (0);   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,MOTIF_MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,MOTIF_MATCH); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,MOTIF_MATCH);  
          }  
        temp = cscore - ((MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MOTIF_MATCH) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MOTIF_MATCH;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MOTIF_MATCH); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MOTIF_MATCH);  
          }  
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      case CONS_MATCH :  
        temp = cscore - ((mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons)) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,SPACER_MATCH) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = SPACER_MATCH;  
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,SPACER_MATCH);    
            }  
          return (*h.access_main)(mat,i - 1,j - 1,SPACER_MATCH);     
          }  
        temp = cscore - (mat->ms->cons_indel) -  (0);    
        if( temp == (*h.access_main)(mat,i - 0,j - 1,CONS_MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,CONS_MATCH);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,CONS_MATCH);   
          }  
        temp = cscore - (mat->ms->cons_indel) -  (0);    
        if( temp == (*h.access_main)(mat,i - 1,j - 0,CONS_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,CONS_MATCH);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,CONS_MATCH);   
          }  
        temp = cscore - ((mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration)) -  (0); 
        if( temp == (*h.access_main)(mat,i - 1,j - 1,CONS_MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = CONS_MATCH;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,CONS_MATCH);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,CONS_MATCH);   
          }  
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      default:   
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found!",i,j,state);    
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_LocalMotifMatrix(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [LocalMotifMatrix *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [LocalMotifMatrix_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_LocalMotifMatrix(LocalMotifMatrix * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,LocalMotifMatrix_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = LocalMotifMatrix_READ_OFF_ERROR; 


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In LocalMotifMatrix matrix special read off - out of bounds on matrix [j is %d in special]",j);  
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case END :     
        /* source SPACER_MATCH is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (0) - (0);     
          if( temp == (*h.access_main)(mat,i - 0,j - 0,SPACER_MATCH) )   {  
            *reti = i - 0;   
            *retj = j - 0;   
            *retstate = SPACER_MATCH;    
            *retspecial = FALSE; 
            if( cellscore != NULL)   {  
              *cellscore = cscore - (*h.access_main)(mat,i-0,j-0,SPACER_MATCH);  
              }  
            return (*h.access_main)(mat,i - 0,j - 0,SPACER_MATCH) ;  
            }  
          } /* end of for i >= 0 */ 
      default:   
        warn("Major problem (!) - in LocalMotifMatrix read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state); 
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_LocalMotifMatrix(mat)
 *
 * Descrip:    This function calculates the LocalMotifMatrix matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_LocalMotifMatrix
 *
 *
 * Arg:        mat [UNKN ] LocalMotifMatrix which contains explicit basematrix memory [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_LocalMotifMatrix(LocalMotifMatrix * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_LocalMotifMatrix, passed a non Explicit matrix type, cannot calculate!");   
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("LocalMotifMatrix Matrix calculation: ");    
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_SPECIAL(mat,i-1,j-1,START) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH) = score; 


        /* state SPACER_MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,END) )  {  
          LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH) = score;  


        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);   
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH) = score;   


        /* Finished calculating state CONS_MATCH */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_LocalMotifMatrix(mat,dpenv)
 *
 * Descrip:    This function calculates the LocalMotifMatrix matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] LocalMotifMatrix which contains explicit basematrix memory [LocalMotifMatrix *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_LocalMotifMatrix(LocalMotifMatrix * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_LocalMotifMatrix, passed a non Explicit matrix type, cannot calculate!");   
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


    for(j=startj-1;j<endj;j++)   {  
      for(i=1;i<mat->leni;i++)   {  
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH) = NEGI;   
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;    
        LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH) = NEGI; 
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,START) = 0;  
      LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,END) = NEGI; 
      }  


    start_reporting("LocalMotifMatrix Matrix calculation: ");    
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
          LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH) = NEGI; 
          LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH) = NEGI;  
          LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH) = NEGI;   
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state SPACER_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_duration);  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-0,SPACER_MATCH) + mat->ms->spacer_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-0,j-1,SPACER_MATCH) + mat->ms->spacer_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,CONS_MATCH) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state SPACER_MATCH */ 
        temp = LocalMotifMatrix_EXPL_SPECIAL(mat,i-1,j-1,START) + mat->ms->comp_spacer->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)];    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for SPACER_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_EXPL_MATRIX(mat,i,j,SPACER_MATCH) = score; 


        /* state SPACER_MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,END) )  {  
          LocalMotifMatrix_EXPL_SPECIAL(mat,i,j,END) = temp;     
          }  




        /* Finished calculating state SPACER_MATCH */ 


        /* For state MOTIF_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,MOTIF_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->motif_duration);  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-0,MOTIF_MATCH) + mat->ms->motif_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MOTIF_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-0,j-1,MOTIF_MATCH) + mat->ms->motif_indel;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state MOTIF_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (MOTIF_MATRIX_IN(mat->motif,i,j,mat->ms->comp_in_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)])+mat->ms->spacer_to_motif);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MOTIF_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_EXPL_MATRIX(mat,i,j,MOTIF_MATCH) = score;  


        /* Finished calculating state MOTIF_MATCH */ 


        /* For state CONS_MATCH */ 
        /* setting first movement to score */ 
        score = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,CONS_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->cons_duration);   
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-0,CONS_MATCH) + mat->ms->cons_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state CONS_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-0,j-1,CONS_MATCH) + mat->ms->cons_indel;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state SPACER_MATCH to state CONS_MATCH */ 
        temp = LocalMotifMatrix_EXPL_MATRIX(mat,i-1,j-1,SPACER_MATCH) + (mat->ms->comp_out_motif->score[CSEQ_DNA_BASE(mat->query,i)][CSEQ_DNA_BASE(mat->target,j)]+mat->ms->spacer_to_cons);     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for CONS_MATCH */ 
        /* Add any movement independant score and put away */ 
         LocalMotifMatrix_EXPL_MATRIX(mat,i,j,CONS_MATCH) = score;   


        /* Finished calculating state CONS_MATCH */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  LocalMotifMatrix_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * LocalMotifMatrix_alloc(void) 
{
    LocalMotifMatrix * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(LocalMotifMatrix *) ckalloc (sizeof(LocalMotifMatrix))) == NULL)    {  
      warn("LocalMotifMatrix_alloc failed ");    
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


/* Function:  free_LocalMotifMatrix(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [LocalMotifMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [LocalMotifMatrix *]
 *
 */
LocalMotifMatrix * free_LocalMotifMatrix(LocalMotifMatrix * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a LocalMotifMatrix obj. Should be trappable");  
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
    /* obj->query is linked in */ 
    /* obj->target is linked in */ 
    /* obj->motif is linked in */ 
    /* obj->ms is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
