#ifdef _cplusplus
extern "C" {
#endif
#include "foldbackmat.h"




# line 8 "foldbackmat.c"


  /*****************   C functions  ****************/
  /*             Written using dynamite            */
  /*            Sat Jun 16 08:57:12 2012           */
  /*            email birney@sanger.ac.uk          */
  /* http://www.sanger.ac.uk/Users/birney/dynamite */
  /*************************************************/


  /* Please report any problems or bugs to         */
  /* Ewan Birney, birney@sanger.ac.uk              */


/* basic set of macros to map states to numbers */ 
#define MATCH 0  
#define STUTTER 1    
#define INSERT 2 
#define DELETE 3 


#define START 0  
#define LOOP 1   
#define END 2    


#define FoldBackMat_EXPL_MATRIX(this_matrix,i,j,STATE) this_matrix->basematrix->matrix[((j+1)*4)+STATE][i+1] 
#define FoldBackMat_EXPL_SPECIAL(matrix,i,j,STATE) matrix->basematrix->specmatrix[STATE][j+1]    
#define FoldBackMat_READ_OFF_ERROR -3
   


#define FoldBackMat_VSMALL_MATRIX(mat,i,j,STATE) mat->basematrix->matrix[(j+2)%2][((i+1)*4)+STATE]   
#define FoldBackMat_VSMALL_SPECIAL(mat,i,j,STATE) mat->basematrix->specmatrix[(j+2)%2][STATE]    




#define FoldBackMat_SHATTER_SPECIAL(matrix,i,j,STATE) matrix->shatter->special[STATE][j] 
#define FoldBackMat_SHATTER_MATRIX(matrix,i,j,STATE)  fetch_cell_value_ShatterMatrix(mat->shatter,i,j,STATE) 


/* Function:  PackAln_read_Shatter_FoldBackMat(mat)
 *
 * Descrip:    Reads off PackAln from shatter matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Shatter_FoldBackMat(FoldBackMat * mat) 
{
    FoldBackMat_access_func_holder holder;   


    holder.access_main    = FoldBackMat_shatter_access_main; 
    holder.access_special = FoldBackMat_shatter_access_special;  
    assert(mat);     
    assert(mat->shatter);    
    return PackAln_read_generic_FoldBackMat(mat,holder); 
}    


/* Function:  FoldBackMat_shatter_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int FoldBackMat_shatter_access_main(FoldBackMat * mat,int i,int j,int state) 
{
    return FoldBackMat_SHATTER_MATRIX(mat,i,j,state);    
}    


/* Function:  FoldBackMat_shatter_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int FoldBackMat_shatter_access_special(FoldBackMat * mat,int i,int j,int state) 
{
    return FoldBackMat_SHATTER_SPECIAL(mat,i,j,state);   
}    


/* Function:  calculate_shatter_FoldBackMat(mat,dpenv)
 *
 * Descrip:    This function calculates the FoldBackMat matrix when in shatter mode
 *
 *
 * Arg:          mat [UNKN ] (null) [FoldBackMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_shatter_FoldBackMat(FoldBackMat * mat,DPEnvelope * dpenv) 
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
    int * SIG_0_1;   
    int * SIG_1_0;   


    leni = mat->leni;    
    lenj = mat->lenj;    


    mat->shatter = new_ShatterMatrix(dpenv,4,lenj,3);    
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


    start_reporting("FoldBackMat Matrix calculation: "); 
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
        SIG_0_1 = fetch_cell_from_ShatterMatrix(mat->shatter,i-0,j-1);   
        SIG_1_0 = fetch_cell_from_ShatterMatrix(mat->shatter,i-1,j-0);   




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = SIG_1_1[MATCH] + 0;  
        /* From state INSERT to state MATCH */ 
        temp = SIG_1_1[INSERT] + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = SIG_1_1[DELETE] + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = FoldBackMat_SHATTER_SPECIAL(mat,i-1,j-1,START) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = SIG_1_1[STUTTER] + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         SIG_0_0[MATCH] = score; 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->lgap) + (0) ;   
        if( temp > FoldBackMat_SHATTER_SPECIAL(mat,i,j,LOOP) )   {  
          FoldBackMat_SHATTER_SPECIAL(mat,i,j,LOOP) = temp;  
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > FoldBackMat_SHATTER_SPECIAL(mat,i,j,END) )    {  
          FoldBackMat_SHATTER_SPECIAL(mat,i,j,END) = temp;   
          }  




        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = SIG_1_1[STUTTER] + 0;    
        /* From state STUTTER to state STUTTER */ 
        temp = SIG_0_1[STUTTER] + mat->stgap;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = SIG_1_0[STUTTER] + mat->stgap;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state STUTTER */ 
        temp = SIG_1_1[MATCH] + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         SIG_0_0[STUTTER] = score;   


        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = SIG_0_1[MATCH] + mat->qgap;  
        /* From state INSERT to state INSERT */ 
        temp = SIG_0_1[INSERT] + mat->qext;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[INSERT] = score;    


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = SIG_1_0[MATCH] + mat->tgap;  
        /* From state DELETE to state DELETE */ 
        temp = SIG_1_0[DELETE] + mat->text;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         SIG_0_0[DELETE] = score;    


        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = FoldBackMat_SHATTER_SPECIAL(mat,0,j,LOOP); 


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = FoldBackMat_SHATTER_SPECIAL(mat,0,j - 1,LOOP) + (mat->lext) + (0);  
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      FoldBackMat_SHATTER_SPECIAL(mat,0,j,LOOP) = score; 
      /* Finished updating state LOOP */ 




      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  search_FoldBackMat(dbsi,out,query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function makes a database search of FoldBackMat
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:           dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:            out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:           comp [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        stutter [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:           qgap [UNKN ] Undocumented argument [Score]
 * Arg:           qext [UNKN ] Undocumented argument [Score]
 * Arg:           tgap [UNKN ] Undocumented argument [Score]
 * Arg:           text [UNKN ] Undocumented argument [Score]
 * Arg:           lgap [UNKN ] Undocumented argument [Score]
 * Arg:           lext [UNKN ] Undocumented argument [Score]
 * Arg:          stgap [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type search_FoldBackMat(DBSearchImpl * dbsi,Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap) 
{
    if( out == NULL )    {  
      warn("Passed in a null Hscore object into search_FoldBackMat. Can't process results!");    
      return SEARCH_ERROR;   
      }  
    if( dbsi == NULL )   {  
      warn("Passed in a null DBSearchImpl object into search_FoldBackMat. Can't process results!");  
      return SEARCH_ERROR;   
      }  
    if( dbsi->trace_level > 5 )  
      warn("Asking for trace level of %d in database search for FoldBackMat, but it was compiled with a trace level of 32. Not all trace statements can be shown",dbsi->trace_level);    
    switch(dbsi->type)   { /*switch on implementation*/ 
      case DBSearchImpl_Serial : 
        return serial_search_FoldBackMat(out,query,target ,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap);    
      case DBSearchImpl_Pthreads :   
        warn("This matrix FoldBackMat was not dyc compiled with thread support");    
        return SEARCH_ERROR; 
      default :  
        warn("database search implementation %s was not provided in the compiled dynamite file from FoldBackMat",impl_string_DBSearchImpl(dbsi));    
        return SEARCH_ERROR; 
      } /* end of switch on implementation */ 


}    


/* Function:  serial_search_FoldBackMat(out,query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function makes a database search of FoldBackMat
 *             It is a single processor implementation
 *
 *
 * Arg:            out [UNKN ] Undocumented argument [Hscore *]
 * Arg:          query [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:         target [UNKN ] Undocumented argument [ComplexSequence*]
 * Arg:           comp [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:        stutter [UNKN ] Undocumented argument [DnaMatrix*]
 * Arg:           qgap [UNKN ] Undocumented argument [Score]
 * Arg:           qext [UNKN ] Undocumented argument [Score]
 * Arg:           tgap [UNKN ] Undocumented argument [Score]
 * Arg:           text [UNKN ] Undocumented argument [Score]
 * Arg:           lgap [UNKN ] Undocumented argument [Score]
 * Arg:           lext [UNKN ] Undocumented argument [Score]
 * Arg:          stgap [UNKN ] Undocumented argument [Score]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type serial_search_FoldBackMat(Hscore * out,ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap) 
{
    int db_status;   
    int score;   
    int query_pos = 0;   
    int target_pos = 0;  
    DataScore * ds;  


    push_errormsg_stack("Before any actual search in db searching"); 


    target_pos = 0;  




    /* No maximum length - allocated on-the-fly */ 
    score = score_only_FoldBackMat(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap);    
    if( should_store_Hscore(out,score) == TRUE )     { /*if storing datascore*/ 
      ds = new_DataScore_from_storage(out);  
      if( ds == NULL )   {  
        warn("FoldBackMat search had a memory error in allocating a new_DataScore (?a leak somewhere - DataScore is a very small datastructure");    
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


/* Function:  score_only_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function just calculates the score for the matrix
 *             I am pretty sure we can do this better, but hey, for the moment...
 *             It calls /allocate_FoldBackMat_only
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int score_only_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap) 
{
    int bestscore = NEGI;    
    int i;   
    int j;   
    int k;   
    FoldBackMat * mat;   


    mat = allocate_FoldBackMat_only(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap);   
    if( mat == NULL )    {  
      warn("Memory allocation error in the db search - unable to communicate to calling function. this spells DIASTER!");    
      return NEGI;   
      }  
    if((mat->basematrix = BaseMatrix_alloc_matrix_and_specials(2,(mat->leni + 1) * 4,2,3)) == NULL)  {  
      warn("Score only matrix for FoldBackMat cannot be allocated, (asking for 1  by %d  cells)",mat->leni*4);   
      mat = free_FoldBackMat(mat);   
      return 0;  
      }  
    mat->basematrix->type = BASEMATRIX_TYPE_VERYSMALL;   


    /* Now, initiate matrix */ 
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->leni;i++)    {  
        for(k=0;k<4;k++) 
          FoldBackMat_VSMALL_MATRIX(mat,i,j,k) = NEGI;   
        }  
      FoldBackMat_VSMALL_SPECIAL(mat,i,j,START) = 0; 
      FoldBackMat_VSMALL_SPECIAL(mat,i,j,LOOP) = NEGI;   
      FoldBackMat_VSMALL_SPECIAL(mat,i,j,END) = NEGI;    
      }  


    /* Ok, lets do-o-o-o-o it */ 


    for(j=0;j<mat->lenj;j++) { /*for all target positions*/ 
      auto int score;    
      auto int temp;     
      for(i=0;i<mat->leni;i++)   { /*for all query positions*/ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + 0;    
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-1,INSERT) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-1,DELETE) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = FoldBackMat_VSMALL_SPECIAL(mat,i-1,j-1,START) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-1,STUTTER) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_VSMALL_MATRIX(mat,i,j,MATCH) = score;   


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->lgap) + (0) ;   
        if( temp > FoldBackMat_VSMALL_SPECIAL(mat,i,j,LOOP) )    {  
          FoldBackMat_VSMALL_SPECIAL(mat,i,j,LOOP) = temp;   
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > FoldBackMat_VSMALL_SPECIAL(mat,i,j,END) )     {  
          FoldBackMat_VSMALL_SPECIAL(mat,i,j,END) = temp;    
          }  




        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-1,STUTTER) + 0;  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-1,MATCH) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_VSMALL_MATRIX(mat,i,j,STUTTER) = score; 


        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_VSMALL_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;    
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_VSMALL_MATRIX(mat,i,j,INSERT) = score;  


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;    
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_VSMALL_MATRIX(mat,i-1,j-0,DELETE) + mat->text;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_VSMALL_MATRIX(mat,i,j,DELETE) = score;  


        /* Finished calculating state DELETE */ 
        } /* end of for all query positions */ 




      /* Special state START has no special to special movements */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = FoldBackMat_VSMALL_SPECIAL(mat,0,j,LOOP);  


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = FoldBackMat_VSMALL_SPECIAL(mat,0,j - 1,LOOP) + (mat->lext) + (0);   
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      FoldBackMat_VSMALL_SPECIAL(mat,0,j,LOOP) = score;  
      /* Finished updating state LOOP */ 




      /* Special state END has no special to special movements */ 
      if( bestscore < FoldBackMat_VSMALL_SPECIAL(mat,0,j,END) )  
        bestscore = FoldBackMat_VSMALL_SPECIAL(mat,0,j,END);     
      } /* end of for all target positions */ 


    mat = free_FoldBackMat(mat);     
    return bestscore;    
}    


/* Function:  PackAln_bestmemory_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_FoldBackMat
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 * Arg:          dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:           dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_bestmemory_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap,DPEnvelope * dpenv,DPRunImpl * dpri) 
{
    long long total; 
    FoldBackMat * mat;   
    PackAln * out;   
    DebugMatrix * de;    
    DPRunImplMemory strategy;    
    assert(dpri);    


    total = query->seq->len * target->seq->len;  
    if( dpri->memory == DPIM_Default )   {  
      if( (total * 4 * sizeof(int)) > 1000*dpri->kbyte_size) {  
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
        if( (mat=allocate_Expl_FoldBackMat(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap,dpri)) == NULL ) {  
          warn("Unable to allocate large FoldBackMat version");  
          return NULL;   
          }  
        calculate_dpenv_FoldBackMat(mat,dpenv);  
        out =  PackAln_read_Expl_FoldBackMat(mat);   
        }  
      else   {  
        mat = allocate_FoldBackMat_only(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap);   
        calculate_shatter_FoldBackMat(mat,dpenv);    
        out = PackAln_read_Shatter_FoldBackMat(mat);     
        }  
      }  
    else {  
      if( strategy == DPIM_Linear )  {  
        /* use small implementation */ 
        if( (mat=allocate_Small_FoldBackMat(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap)) == NULL ) {  
          warn("Unable to allocate small FoldBackMat version");  
          return NULL;   
          }  
        out = PackAln_calculate_Small_FoldBackMat(mat,dpenv);    
        }  
      else   {  
        /* use Large implementation */ 
        if( (mat=allocate_Expl_FoldBackMat(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap,dpri)) == NULL ) {  
          warn("Unable to allocate large FoldBackMat version");  
          return NULL;   
          }  
        if( dpri->debug == TRUE) {  
          fatal("Asked for dydebug, but dynamite file not compiled with -g. Need to recompile dynamite source"); 
          }  
        else {  
          calculate_FoldBackMat(mat);    
          out =  PackAln_read_Expl_FoldBackMat(mat); 
          }  
        }  
      }  


    mat = free_FoldBackMat(mat);     
    return out;  
}    


/* Function:  allocate_FoldBackMat_only(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function only allocates the FoldBackMat structure
 *             checks types where possible and determines leni and lenj
 *             The basematrix area is delt with elsewhere
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * allocate_FoldBackMat_only(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap) 
{
    FoldBackMat * out;   


    if((out= FoldBackMat_alloc()) == NULL)   {  
      warn("Allocation of basic FoldBackMat structure failed...");   
      return NULL;   
      }  


    out->query = query;  
    out->target = target;    
    out->comp = comp;    
    out->stutter = stutter;  
    out->qgap = qgap;    
    out->qext = qext;    
    out->tgap = tgap;    
    out->text = text;    
    out->lgap = lgap;    
    out->lext = lext;    
    out->stgap = stgap;  
    out->leni = query->seq->len;     
    out->lenj = target->seq->len;    
    return out;  
}    


/* Function:  allocate_Expl_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap,dpri)
 *
 * Descrip:    This function allocates the FoldBackMat structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_FoldBackMat_only
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 * Arg:           dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * allocate_Expl_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap,DPRunImpl * dpri) 
{
    FoldBackMat * out;   


    out = allocate_FoldBackMat_only(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap);   
    if( out == NULL )    
      return NULL;   
    if( dpri->should_cache == TRUE ) {  
      if( dpri->cache != NULL )  {  
        if( dpri->cache->maxleni >= (out->lenj+1)*4 && dpri->cache->maxlenj >= (out->leni+1))    
          out->basematrix = hard_link_BaseMatrix(dpri->cache);   
        else 
          dpri->cache = free_BaseMatrix(dpri->cache);    
        }  
      }  
    if( out->basematrix == NULL )    {  
      if( (out->basematrix = BaseMatrix_alloc_matrix_and_specials((out->lenj+1)*4,(out->leni+1),3,out->lenj+1)) == NULL) {  
        warn("Explicit matrix FoldBackMat cannot be allocated, (asking for %d by %d main cells)",out->leni,out->lenj);   
        free_FoldBackMat(out);   
        return NULL; 
        }  
      }  
    if( dpri->should_cache == TRUE && dpri->cache == NULL)   
      dpri->cache = hard_link_BaseMatrix(out->basematrix);   
    out->basematrix->type = BASEMATRIX_TYPE_EXPLICIT;    
    init_FoldBackMat(out);   
    return out;  
}    


/* Function:  init_FoldBackMat(mat)
 *
 * Descrip:    This function initates FoldBackMat matrix when in explicit mode
 *             Called in /allocate_Expl_FoldBackMat
 *
 *
 * Arg:        mat [UNKN ] FoldBackMat which contains explicit basematrix memory [FoldBackMat *]
 *
 */
void init_FoldBackMat(FoldBackMat * mat) 
{
    register int i;  
    register int j;  
    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT)   {  
      warn("Cannot iniate matrix, is not an explicit memory type and you have assummed that");   
      return;    
      }  


    for(i= (-1);i<mat->query->seq->len;i++)  {  
      for(j= (-1);j<2;j++)   {  
        FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;   
        FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER) = NEGI; 
        FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;  
        FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;  
        }  
      }  
    for(j= (-1);j<mat->target->seq->len;j++) {  
      for(i= (-1);i<2;i++)   {  
        FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;   
        FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER) = NEGI; 
        FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;  
        FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;  
        }  
      FoldBackMat_EXPL_SPECIAL(mat,i,j,START) = 0;   
      FoldBackMat_EXPL_SPECIAL(mat,i,j,LOOP) = NEGI; 
      FoldBackMat_EXPL_SPECIAL(mat,i,j,END) = NEGI;  
      }  
    return;  
}    


/* Function:  recalculate_PackAln_FoldBackMat(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by FoldBackMat
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 */
void recalculate_PackAln_FoldBackMat(PackAln * pal,FoldBackMat * mat) 
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
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == INSERT )  {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == DELETE )  {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          if( offj == 1 && prev->state == (START+4) )    {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == STUTTER ) {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          warn("In recaluclating PackAln with state MATCH, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);  
          break; 
        case STUTTER :   
          if( offi == 1 && offj == 1 && prev->state == STUTTER ) {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == STUTTER ) {  
            pau->score = mat->stgap + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));     
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == STUTTER ) {  
            pau->score = mat->stgap + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));     
            continue;    
            }  
          if( offi == 1 && offj == 1 && prev->state == MATCH )   {  
            pau->score = 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
            continue;    
            }  
          warn("In recaluclating PackAln with state STUTTER, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state);    
          break; 
        case INSERT :    
          if( offi == 0 && offj == 1 && prev->state == MATCH )   {  
            pau->score = mat->qgap + (0);    
            continue;    
            }  
          if( offi == 0 && offj == 1 && prev->state == INSERT )  {  
            pau->score = mat->qext + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state INSERT, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case DELETE :    
          if( offi == 1 && offj == 0 && prev->state == MATCH )   {  
            pau->score = mat->tgap + (0);    
            continue;    
            }  
          if( offi == 1 && offj == 0 && prev->state == DELETE )  {  
            pau->score = mat->text + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state DELETE, from [%d,%d,%d], got a bad source state. Error!",offi,offj,prev->state); 
          break; 
        case (START+4) :     
          warn("In recaluclating PackAln with state START, got a bad source state. Error!"); 
          break; 
        case (LOOP+4) :  
          if( offj == 0 && prev->state == MATCH )    {  
            /* i here comes from the previous state ;) - not the real one */ 
            i = prev->i; 
            pau->score = mat->lgap + (0);    
            continue;    
            }  
          if( offj == 1 && prev->state == (LOOP+4) ) {  
            pau->score = mat->lext + (0);    
            continue;    
            }  
          warn("In recaluclating PackAln with state LOOP, got a bad source state. Error!");  
          break; 
        case (END+4) :   
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
#define FoldBackMat_HIDDEN_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[(j-hiddenj+1)][(i+1)*4+state])   
#define FoldBackMat_DC_SHADOW_MATRIX(thismatrix,i,j,state) (thismatrix->basematrix->matrix[((j+2)*8) % 16][(i+1)*4+state])   
#define FoldBackMat_HIDDEN_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state][(j+1)])  
#define FoldBackMat_DC_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)]) 
#define FoldBackMat_DC_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->matrix[((((j+2)*8)+(shadow+1)) % 16)][(i+1)*4 + state])    
#define FoldBackMat_DC_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)]) 
#define FoldBackMat_DC_OPT_SHADOW_MATRIX(thismatrix,i,j,state) (score_pointers[(((j+1)% 1) * (leni+1) * 4) + ((i+1) * 4) + (state)]) 
#define FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(thismatrix,i,j,state,shadow) (shadow_pointers[(((j+1)% 1) * (leni+1) * 32) + ((i+1) * 32) + (state * 8) + shadow+1]) 
#define FoldBackMat_DC_OPT_SHADOW_SPECIAL(thismatrix,i,j,state) (thismatrix->basematrix->specmatrix[state*8][(j+1)]) 
/* Function:  allocate_Small_FoldBackMat(query,target,comp,stutter,qgap,qext,tgap,text,lgap,lext,stgap)
 *
 * Descrip:    This function allocates the FoldBackMat structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_FoldBackMat_only
 *
 *
 * Arg:          query [UNKN ] query data structure [ComplexSequence*]
 * Arg:         target [UNKN ] target data structure [ComplexSequence*]
 * Arg:           comp [UNKN ] Resource [DnaMatrix*]
 * Arg:        stutter [UNKN ] Resource [DnaMatrix*]
 * Arg:           qgap [UNKN ] Resource [Score]
 * Arg:           qext [UNKN ] Resource [Score]
 * Arg:           tgap [UNKN ] Resource [Score]
 * Arg:           text [UNKN ] Resource [Score]
 * Arg:           lgap [UNKN ] Resource [Score]
 * Arg:           lext [UNKN ] Resource [Score]
 * Arg:          stgap [UNKN ] Resource [Score]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
#define FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(thismatrix,i,j,state,shadow) (thismatrix->basematrix->specmatrix[state*8 +shadow+1][(j+1)]) 
FoldBackMat * allocate_Small_FoldBackMat(ComplexSequence* query,ComplexSequence* target ,DnaMatrix* comp,DnaMatrix* stutter,Score qgap,Score qext,Score tgap,Score text,Score lgap,Score lext,Score stgap) 
{
    FoldBackMat * out;   


    out = allocate_FoldBackMat_only(query, target , comp, stutter, qgap, qext, tgap, text, lgap, lext, stgap);   
    if( out == NULL )    
      return NULL;   
    out->basematrix = BaseMatrix_alloc_matrix_and_specials(16,(out->leni + 1) * 4,24,out->lenj+1);   
    if(out == NULL)  {  
      warn("Small shadow matrix FoldBackMat cannot be allocated, (asking for 2 by %d main cells)",out->leni+2);  
      free_FoldBackMat(out);     
      return NULL;   
      }  
    out->basematrix->type = BASEMATRIX_TYPE_SHADOW;  
    return out;  
}    


/* Function:  PackAln_calculate_Small_FoldBackMat(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for FoldBackMat structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_FoldBackMat 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_FoldBackMat 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_calculate_Small_FoldBackMat(FoldBackMat * mat,DPEnvelope * dpenv) 
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
      warn("Could not calculate packaln small for FoldBackMat due to wrong type of matrix"); 
      return NULL;   
      }  


    out = PackAln_alloc_std();   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_FoldBackMat(mat,dpenv);  
    score = start_end_find_end_FoldBackMat(mat,&endj);   
    out->score = score;  
    stopstate = END;
    
    /* Special to specials: have to eat up in strip and then drop back to full_dc for intervening bits */ 
    log_full_error(REPORT,0,"End at %d Score %d",endj,score);    
    stop_reporting();    
    for(;;)  { /*while there are more special bits to recover*/ 
      start_reporting("Special cell aln end   %d:",endj);    
      if( read_special_strip_FoldBackMat(mat,0,endj,stopstate,&endj,&startstate,out) == FALSE )  {  
        warn("Problem in reading off special state system... going to return partial alignment");    
        break;   
        }  
      if( startstate == START || endj <= 0)  {  
        log_full_error(REPORT,0,"Recovered complete alignment"); 
        stop_reporting();    
        break;   
        }  


      log_full_error(REPORT,0,"Finished to %d",endj);    
      stop_reporting();  


      /* Ok... have to eat up another piece of matrix <sigh> */ 
      temp = startstate; 
      starti = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,0);  
      startj = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,1);  
      startstate = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,2);  
      stopi = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,3);   
      stopj = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,4);   
      stopstate = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,5);   


      /* Get out the score of this block. V. important! */ 
      temp = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,endj,temp,6);    
      totalj = stopj - startj;   
      donej  = 0;    
      start_reporting("Main matrix  aln [%d,%d]:",startj,stopj);     
      if(full_dc_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out,&donej,totalj,dpenv) == FALSE)   {  
        warn("In the alignment FoldBackMat [%d,%d][%d,%d], got a problem. Please report bug ... giving you back a partial alignment",starti,startj,stopi,stopj); 
        return out;  
        }  


      /* now have to figure out which special we came from... yikes */ 
      max_matrix_to_special_FoldBackMat(mat,starti,startj,startstate,temp,&stopi,&stopj,&stopstate,&temp,NULL);  
      if( stopi == FoldBackMat_READ_OFF_ERROR)   {  
        warn("In FoldBackMat read off ending at %d ... got a bad matrix to special read off... returning partial alignment",startj); 
        invert_PackAln(out); 
        recalculate_PackAln_FoldBackMat(out,mat);    
        return out;  
        }  
      /* if at start, break, otherwise, back to eat another strip */ 
      if( stopstate == START)    {  
        log_full_error(REPORT,0,"Recovered complete alignment      ");   
        stop_reporting();    
        break;   
        }  
      log_full_error(REPORT,0,"Finished  alignment to %d           ",startj);    
      stop_reporting();  
      endj = stopj;  
      /* stopstate is correct as it is */ 
      } /* end of while there are more special bits to recover */ 
    invert_PackAln(out); 
    recalculate_PackAln_FoldBackMat(out,mat);    
    return out;  


}    


/* Function:  AlnRangeSet_calculate_Small_FoldBackMat(mat)
 *
 * Descrip:    This function calculates an alignment for FoldBackMat structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_FoldBackMat 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_FoldBackMat
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_FoldBackMat 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_calculate_Small_FoldBackMat(FoldBackMat * mat) 
{
    AlnRangeSet * out;   


    start_reporting("Find start end points: ");  
    dc_optimised_start_end_calc_FoldBackMat(mat,NULL);   
    log_full_error(REPORT,0,"Calculated");   


    out = AlnRangeSet_from_FoldBackMat(mat); 
    return out;  
}    


/* Function:  AlnRangeSet_from_FoldBackMat(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for FoldBackMat structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_FoldBackMat 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_FoldBackMat
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * AlnRangeSet_from_FoldBackMat(FoldBackMat * mat) 
{
    AlnRangeSet * out;   
    AlnRange * temp; 
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_FoldBackMat");   
      return NULL;   
      }  


    out = AlnRangeSet_alloc_std();   
    /* Find the end position */ 
    out->score = start_end_find_end_FoldBackMat(mat,&jpos);  
    state = END; 


    while( (temp = AlnRange_build_FoldBackMat(mat,jpos,state,&jpos,&state)) != NULL) 
      add_AlnRangeSet(out,temp); 
    return out;  
}    


/* Function:  AlnRange_build_FoldBackMat(mat,stopj,stopspecstate,startj,startspecstate)
 *
 * Descrip:    This function calculates a single start/end set in linear space
 *             Really a sub-routine for /AlnRangeSet_from_PackAln_FoldBackMat
 *
 *
 * Arg:                   mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:                 stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopspecstate [UNKN ] Undocumented argument [int]
 * Arg:                startj [UNKN ] Undocumented argument [int *]
 * Arg:        startspecstate [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRange *]
 *
 */
AlnRange * AlnRange_build_FoldBackMat(FoldBackMat * mat,int stopj,int stopspecstate,int * startj,int * startspecstate) 
{
    AlnRange * out;  
    int jpos;    
    int state;   


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("Bad error! - non shadow matrix type in AlnRangeSet_from_FoldBackMat");   
      return NULL;   
      }  


    /* Assumme that we have specials (we should!). Read back along the specials till we have the finish point */ 
    if( read_special_strip_FoldBackMat(mat,0,stopj,stopspecstate,&jpos,&state,NULL) == FALSE)    {  
      warn("In AlnRanger_build_FoldBackMat alignment ending at %d, unable to read back specials. Will (evenutally) return a partial range set... BEWARE!",stopj);    
      return NULL;   
      }  
    if( state == START || jpos <= 0) 
      return NULL;   


    out = AlnRange_alloc();  


    out->starti = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,0);  
    out->startj = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,1);  
    out->startstate = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,2);  
    out->stopi = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,3);   
    out->stopj = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,4);   
    out->stopstate = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,5);   
    out->startscore = FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,jpos,state,6);  
    out->stopscore = FoldBackMat_DC_SHADOW_SPECIAL(mat,0,jpos,state);    


    /* Now, we have to figure out where this state came from in the specials */ 
    max_matrix_to_special_FoldBackMat(mat,out->starti,out->startj,out->startstate,out->startscore,&jpos,startj,startspecstate,&state,NULL);  
    if( jpos == FoldBackMat_READ_OFF_ERROR)  {  
      warn("In AlnRange_build_FoldBackMat alignment ending at %d, with aln range between %d-%d in j, unable to find source special, returning this range, but this could get tricky!",stopj,out->startj,out->stopj); 
      return out;    
      }  


    /* Put in the correct score for startstate, from the special */ 
    out->startscore = FoldBackMat_DC_SHADOW_SPECIAL(mat,0,*startj,*startspecstate);  
    /* The correct j coords have been put into startj, startspecstate... so just return out */ 
    return out;  
}    


/* Function:  read_hidden_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
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
boolean read_hidden_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out) 
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


      max_hidden_FoldBackMat(mat,startj,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore); 


      if( i == FoldBackMat_READ_OFF_ERROR)   {  
        warn("In FoldBackMat hidden read off, between %d:%d,%d:%d - at got bad read off. Problem!",starti,startj,stopi,stopj);   
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
        warn("In FoldBackMat hidden read off, between %d:%d,%d:%d - hit start cell, but not in start state. Can't be good!.",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      }  
    warn("In FoldBackMat hidden read off, between %d:%d,%d:%d - gone past start cell (now in %d,%d,%d), can't be good news!.",starti,startj,stopi,stopj,i,j,state);  
    return FALSE;    
}    


/* Function:  max_hidden_FoldBackMat(mat,hiddenj,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
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
int max_hidden_FoldBackMat(FoldBackMat * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = FoldBackMat_READ_OFF_ERROR;  


    if( i < 0 || j < 0 || i > mat->query->seq->len || j > mat->target->seq->len) {  
      warn("In FoldBackMat matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state); 
      return -1; 
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = FoldBackMat_HIDDEN_MATRIX(mat,i,j,state);   
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,STUTTER) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,STUTTER);    
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,STUTTER);     
          }  
        /* Not allowing special sources.. skipping START */ 
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,DELETE); 
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,DELETE);  
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,INSERT); 
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,INSERT);  
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);  
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);   
          }  
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case STUTTER :     
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH);  
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,MATCH);   
          }  
        temp = cscore - (mat->stgap) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));   
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 0,STUTTER) ) {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-0,STUTTER);    
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 0,STUTTER);     
          }  
        temp = cscore - (mat->stgap) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));   
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 0,j - 1,STUTTER) ) {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-0,j-1,STUTTER);    
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 0,j - 1,STUTTER);     
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,STUTTER) ) {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,STUTTER);    
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 1,STUTTER);     
          }  
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->qext) -  (0);  
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-0,j-1,INSERT); 
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 0,j - 1,INSERT);  
          }  
        temp = cscore - (mat->qgap) -  (0);  
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH);  
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 0,j - 1,MATCH);   
          }  
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->text) -  (0);  
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-0,DELETE); 
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 0,DELETE);  
          }  
        temp = cscore - (mat->tgap) -  (0);  
        if( temp == FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-0,MATCH);  
            }  
          return FoldBackMat_HIDDEN_MATRIX(mat,i - 1,j - 0,MATCH);   
          }  
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  read_special_strip_FoldBackMat(mat,stopi,stopj,stopstate,startj,startstate,out)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
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
boolean read_special_strip_FoldBackMat(FoldBackMat * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out) 
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
    while( j > FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4) && state != START)  { /*while more specials to eat up*/ 
      /* Put away current state, if we should */ 
      if(out != NULL)    {  
        pau = PackAlnUnit_alloc();  /* Should deal with memory overflow */ 
        pau->i = i;  
        pau->j = j;  
        pau->state =  state + 4; 
        add_PackAln(out,pau);    
        }  


      max_special_strip_FoldBackMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore); 
      if( i == FoldBackMat_READ_OFF_ERROR)   {  
        warn("In special strip read FoldBackMat, got a bad read off error. Sorry!"); 
        return FALSE;    
        }  
      } /* end of while more specials to eat up */ 


    /* check to see we have not gone too far! */ 
    if( state != START && j < FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,i,j,state,4)) {  
      warn("In special strip read FoldBackMat, at special [%d] state [%d] overshot!",j,state);   
      return FALSE;  
      }  
    /* Put away last state */ 
    if(out != NULL)  {  
      pau = PackAlnUnit_alloc();/* Should deal with memory overflow */ 
      pau->i = i;    
      pau->j = j;    
      pau->state =  state + 4;   
      add_PackAln(out,pau);  
      }  


    /* Put away where we are in startj and startstate */ 
    *startj = j; 
    *startstate = state; 
    return TRUE; 
}    


/* Function:  max_special_strip_FoldBackMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip:    A pretty intense internal function. Deals with read-off only in specials
 *
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
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
int max_special_strip_FoldBackMat(FoldBackMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    int cscore;  


    *reti = (*retj) = (*retstate) = FoldBackMat_READ_OFF_ERROR;  
    if( isspecial == FALSE ) {  
      warn("In special strip max function for FoldBackMat, got a non special start point. Problem! (bad!)"); 
      return (-1);   
      }  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In FoldBackMat matrix special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1; 
      }  


    cscore = FoldBackMat_DC_SHADOW_SPECIAL(mat,i,j,state);   
    switch(state)    { /*switch on special states*/ 
      case START :   
      case LOOP :    
        /* source LOOP is a special */ 
        temp = cscore - (mat->lext) - (0);   
        if( temp == FoldBackMat_DC_SHADOW_SPECIAL(mat,i - 0,j - 1,LOOP) )    {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_DC_SHADOW_SPECIAL(mat,i-0,j-1,LOOP);   
            }  
          return FoldBackMat_DC_SHADOW_MATRIX(mat,i - 0,j - 1,LOOP) ;    
          }  
        /* Source MATCH is not a special */ 
      case END :     
        /* Source MATCH is not a special */ 
      default:   
        warn("Major problem (!) - in FoldBackMat special strip read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);    
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  max_matrix_to_special_FoldBackMat(mat,i,j,state,cscore,reti,retj,retstate,retspecial,cellscore)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
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
int max_matrix_to_special_FoldBackMat(FoldBackMat * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore) 
{
    int temp;    
    *reti = (*retj) = (*retstate) = FoldBackMat_READ_OFF_ERROR;  


    if( j < 0 || j > mat->lenj)  {  
      warn("In FoldBackMat matrix to special read off - out of bounds on matrix [j is %d in special]",j);    
      return -1; 
      }  


    switch(state)    { /*Switch state */ 
      case MATCH :   
        /* Source STUTTER is not a special, should not get here! */ 
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == FoldBackMat_DC_SHADOW_SPECIAL(mat,i - 1,j - 1,START) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = START; 
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - FoldBackMat_DC_SHADOW_SPECIAL(mat,i-1,j-1,START);  
            }  
          return FoldBackMat_DC_SHADOW_MATRIX(mat,i - 1,j - 1,START) ;   
          }  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in FoldBackMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case STUTTER :     
        /* Source MATCH is not a special, should not get here! */ 
        /* Source STUTTER is not a special, should not get here! */ 
        /* Source STUTTER is not a special, should not get here! */ 
        /* Source STUTTER is not a special, should not get here! */ 
        warn("Major problem (!) - in FoldBackMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case INSERT :  
        /* Source INSERT is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in FoldBackMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      case DELETE :  
        /* Source DELETE is not a special, should not get here! */ 
        /* Source MATCH is not a special, should not get here! */ 
        warn("Major problem (!) - in FoldBackMat matrix to special read off, position %d,%d state %d no source found!",i,j,state);   
        return (-1); 
      default:   
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 


}    


/* Function:  calculate_hidden_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:            starti [UNKN ] Undocumented argument [int]
 * Arg:            startj [UNKN ] Undocumented argument [int]
 * Arg:        startstate [UNKN ] Undocumented argument [int]
 * Arg:             stopi [UNKN ] Undocumented argument [int]
 * Arg:             stopj [UNKN ] Undocumented argument [int]
 * Arg:         stopstate [UNKN ] Undocumented argument [int]
 * Arg:             dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void calculate_hidden_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv) 
{
    register int i;  
    register int j;  
    register int score;  
    register int temp;   
    register int hiddenj;    


    hiddenj = startj;    


    init_hidden_FoldBackMat(mat,starti,startj,stopi,stopj);  


    FoldBackMat_HIDDEN_MATRIX(mat,starti,startj,startstate) = 0; 


    for(j=startj;j<=stopj;j++)   {  
      for(i=starti;i<=stopi;i++) {  
        /* Should *not* do very first cell as this is the one set to zero in one state! */ 
        if( i == starti && j == startj ) 
          continue;  
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          FoldBackMat_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;   
          FoldBackMat_HIDDEN_MATRIX(mat,i,j,STUTTER) = NEGI;     
          FoldBackMat_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;  
          FoldBackMat_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;  
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + 0;    
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,INSERT) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,DELETE) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,STUTTER) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_HIDDEN_MATRIX(mat,i,j,MATCH) = score;   
        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,STUTTER) + 0;  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-1,MATCH) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_HIDDEN_MATRIX(mat,i,j,STUTTER) = score; 
        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_HIDDEN_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;    
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_HIDDEN_MATRIX(mat,i,j,INSERT) = score;  
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;    
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_HIDDEN_MATRIX(mat,i-1,j-0,DELETE) + mat->text;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_HIDDEN_MATRIX(mat,i,j,DELETE) = score;  
        /* Finished calculating state DELETE */ 
        }  
      }  


    return;  
}    


/* Function:  init_hidden_FoldBackMat(mat,starti,startj,stopi,stopj)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 *
 */
void init_hidden_FoldBackMat(FoldBackMat * mat,int starti,int startj,int stopi,int stopj) 
{
    register int i;  
    register int j;  
    register int hiddenj;    


    hiddenj = startj;    
    for(j=(startj-1);j<=stopj;j++)   {  
      for(i=(starti-1);i<=stopi;i++) {  
        FoldBackMat_HIDDEN_MATRIX(mat,i,j,MATCH) = NEGI;
    
        FoldBackMat_HIDDEN_MATRIX(mat,i,j,STUTTER) = NEGI;
  
        FoldBackMat_HIDDEN_MATRIX(mat,i,j,INSERT) = NEGI;
   
        FoldBackMat_HIDDEN_MATRIX(mat,i,j,DELETE) = NEGI;
   
        }  
      }  


    return;  
}    


/* Function:  full_dc_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out,donej,totalj,dpenv)
 *
 * Descrip:    The main divide-and-conquor routine. Basically, call /PackAln_calculate_small_FoldBackMat
 *             Not this function, which is pretty hard core. 
 *             Function is given start/end points (in main matrix) for alignment
 *             It does some checks, decides whether start/end in j is small enough for explicit calc
 *               - if yes, calculates it, reads off into PackAln (out), adds the j distance to donej and returns TRUE
 *               - if no,  uses /do_dc_single_pass_FoldBackMat to get mid-point
 *                          saves midpoint, and calls itself to do right portion then left portion
 *             right then left ensures PackAln is added the 'right' way, ie, back-to-front
 *             returns FALSE on any error, with a warning
 *
 *
 * Arg:               mat [UNKN ] Matrix with small memory implementation [FoldBackMat *]
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
boolean full_dc_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv) 
{
    int lstarti; 
    int lstartj; 
    int lstate;  


    if( mat->basematrix->type != BASEMATRIX_TYPE_SHADOW) {  
      warn("*Very* bad error! - non shadow matrix type in full_dc_FoldBackMat"); 
      return FALSE;  
      }  


    if( starti == -1 || startj == -1 || startstate == -1 || stopi == -1 || stopstate == -1)  {  
      warn("In full dc program, passed bad indices, indices passed were %d:%d[%d] to %d:%d[%d]\n",starti,startj,startstate,stopi,stopj,stopstate);   
      return FALSE;  
      }  


    if( stopj - startj < 5)  {  
      log_full_error(REPORT,0,"[%d,%d][%d,%d] Explicit read off",starti,startj,stopi,stopj);/* Build hidden explicit matrix */ 
      calculate_hidden_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv);    
      *donej += (stopj - startj);   /* Now read it off into out */ 
      if( read_hidden_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,out) == FALSE)  {  
        warn("In full dc, at %d:%d,%d:%d got a bad hidden explicit read off... ",starti,startj,stopi,stopj); 
        return FALSE;    
        }  
      return TRUE;   
      }  


/* In actual divide and conquor */ 
    if( do_dc_single_pass_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,(int)(*donej*100)/totalj) == FALSE)   {  
      warn("In divide and conquor for FoldBackMat, at bound %d:%d to %d:%d, unable to calculate midpoint. Problem!",starti,startj,stopi,stopj);  
      return FALSE;  
      }  


/* Ok... now we have to call on each side of the matrix */ 
/* We have to retrieve left hand side positions, as they will be vapped by the time we call LHS */ 
    lstarti= FoldBackMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,0);   
    lstartj= FoldBackMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,1);   
    lstate = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,2);   


/* Call on right hand side: this lets us do the correct read off */ 
    if( full_dc_FoldBackMat(mat,FoldBackMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,3),FoldBackMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,4),FoldBackMat_DC_SHADOW_MATRIX_SP(mat,stopi,stopj,stopstate,5),stopi,stopj,stopstate,out,donej,totalj,dpenv) == FALSE)   {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  
/* Call on left hand side */ 
    if( full_dc_FoldBackMat(mat,starti,startj,startstate,lstarti,lstartj,lstate,out,donej,totalj,dpenv) == FALSE)    {  
/* Warning already issued, simply chained back up to top */ 
      return FALSE;  
      }  


    return TRUE;     
}    


/* Function:  do_dc_single_pass_FoldBackMat(mat,starti,startj,startstate,stopi,stopj,stopstate,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
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
boolean do_dc_single_pass_FoldBackMat(FoldBackMat * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done) 
{
    int halfj;   
    halfj = startj + ((stopj - startj)/2);   


    init_dc_FoldBackMat(mat);    


    FoldBackMat_DC_SHADOW_MATRIX(mat,starti,startj,startstate) = 0;  
    run_up_dc_FoldBackMat(mat,starti,stopi,startj,halfj-1,dpenv,perc_done);  
    push_dc_at_merge_FoldBackMat(mat,starti,stopi,halfj,&halfj,dpenv);   
    follow_on_dc_FoldBackMat(mat,starti,stopi,halfj,stopj,dpenv,perc_done);  
    return TRUE; 
}    


/* Function:  push_dc_at_merge_FoldBackMat(mat,starti,stopi,startj,stopj,dpenv)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int *]
 * Arg:         dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 */
void push_dc_at_merge_FoldBackMat(FoldBackMat * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv) 
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
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-100); 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = (-100); 
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = NEGI;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,0) = (-100);   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,1) = (-100);   
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-100);    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = (-100);    
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-100);    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = (-100);    
          continue;  
          } /* end of Is not in envelope */ 


        /* For state MATCH, pushing when j - offj <= mergej */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;     
        if( j - 1 <= mergej) {  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = MATCH;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;  
          }  
        else {  
          for(k=0;k<7;k++)   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = INSERT;   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);  
            }  
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = DELETE;   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k);  
            }  
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0;    
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = i-1;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,1) = j-1;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,2) = STUTTER;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,3) = i;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,4) = j;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,5) = MATCH;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,STUTTER,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));    
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;     
        /* Finished with state MATCH */ 


        /* For state STUTTER, pushing when j - offj <= mergej */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0;   
        if( j - 1 <= mergej) {  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,0) = i-1;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,1) = j-1;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,2) = STUTTER;  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,3) = i;    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,4) = j;    
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,5) = STUTTER;  
          }  
        else {  
          for(k=0;k<7;k++)   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,STUTTER,k); 
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;   
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,0) = i-0;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,1) = j-1;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,2) = STUTTER;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,3) = i;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,4) = j;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,5) = STUTTER;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,STUTTER,k);   
            }  
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;   
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,0) = i-1;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,1) = j-0;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,2) = STUTTER;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,3) = i;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,4) = j;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,5) = STUTTER;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,STUTTER,k);   
            }  
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,0) = i-1;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,1) = j-1;    
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,2) = MATCH;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,3) = i;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,4) = j;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,5) = STUTTER;    
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
            }  
          }  
        /* Add any movement independant score */ 
        score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));    
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = score;   
        /* Finished with state STUTTER */ 


        /* For state INSERT, pushing when j - offj <= mergej */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;     
        if( j - 1 <= mergej) {  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0;   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1;   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = MATCH; 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i; 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j; 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;    
          }  
        else {  
          for(k=0;k<7;k++)   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k);    
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;     
        if( temp > score)    {  
          score = temp;  


          if( j - 1 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = i-0; 
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,1) = j-1; 
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,2) = INSERT;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,3) = i;   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,4) = j;   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,5) = INSERT;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k); 
            }  
          }  
        /* Add any movement independant score */ 
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;    
        /* Finished with state INSERT */ 


        /* For state DELETE, pushing when j - offj <= mergej */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;     
        if( j - 0 <= mergej) {  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1;   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0;   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = MATCH; 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i; 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j; 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;    
          }  
        else {  
          for(k=0;k<7;k++)   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k);    
          }  


        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->text;     
        if( temp > score)    {  
          score = temp;  


          if( j - 0 <= mergej)   {  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = i-1; 
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,1) = j-0; 
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,2) = DELETE;  
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,3) = i;   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,4) = j;   
            FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,5) = DELETE;  
            }  
          else   {  
            for(k=0;k<7;k++) 
              FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k); 
            }  
          }  
        /* Add any movement independant score */ 
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;    
        /* Finished with state DELETE */ 
        }  
      }  
    /* Put back j into * stop j so that calling function gets it correct */ 
    if( stopj == NULL)   
      warn("Bad news... NULL stopj pointer in push dc function. This means that calling function does not know how many cells I have done!");    
    else 
      *stopj = j;    


    return;  
}    


/* Function:  follow_on_dc_FoldBackMat(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
void follow_on_dc_FoldBackMat(FoldBackMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = NEGI;  
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;   
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]After  mid-j %5d Cells done %d%%%%",perc_done,startj,(num*100)/total);   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k); 
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,k);  
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,k);  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0;    
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,STUTTER,k); 
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;    
        for(k=0;k<7;k++) 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localshadow[k]; 
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0;   
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,STUTTER,k);   
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,STUTTER,k); 
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;   
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,STUTTER,k); 
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,k);   
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = score;  
        for(k=0;k<7;k++) 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = localshadow[k];   
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,k); 
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,k);  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;   
        for(k=0;k<7;k++) 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;     
        /* shift first shadow numbers */ 
        for(k=0;k<7;k++) 
          localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,k); 
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->text;     
        if( temp  > score )  {  
          score = temp;  
          for(k=0;k<7;k++)   
            localshadow[k] = FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,k);  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;   
        for(k=0;k<7;k++) 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localshadow[k];    
        /* Now figure out if any specials need this score */ 
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  run_up_dc_FoldBackMat(mat,starti,stopi,startj,stopj,dpenv,perc_done)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:           starti [UNKN ] Undocumented argument [int]
 * Arg:            stopi [UNKN ] Undocumented argument [int]
 * Arg:           startj [UNKN ] Undocumented argument [int]
 * Arg:            stopj [UNKN ] Undocumented argument [int]
 * Arg:            dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:        perc_done [UNKN ] Undocumented argument [int]
 *
 */
}    
void run_up_dc_FoldBackMat(FoldBackMat * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done) 
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
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = NEGI;  
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;   
          FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num % 1000 == 0 )    
          log_full_error(REPORT,0,"[%d%%%% done]Before mid-j %5d Cells done %d%%%%",perc_done,stopj,(num*100)/total);    


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;     
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0;    
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = score;    
        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0;   
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;   
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = score;  
        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;     
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = score;   
        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;     
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_DC_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->text;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = score;   
        /* Finished calculating state DELETE */ 
        } /* end of this is strip */ 
      } /* end of for each valid j column */ 


/* Function:  init_dc_FoldBackMat(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 */
}    
void init_dc_FoldBackMat(FoldBackMat * mat) 
{
    register int i;  
    register int j;  
    register int k;  


    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->query->seq->len;i++) {  
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = NEGI;    
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI; 
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI; 
        for(k=0;k<7;k++) {  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = (-1);   
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = (-1); 
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = (-1);  
          FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = (-1);  
          }  
        }  
      }  


    return;  
}    


/* Function:  start_end_find_end_FoldBackMat(mat,endj)
 *
 * Descrip:    First function used to find end of the best path in the special state !end
 *
 *
 * Arg:         mat [UNKN ] Matrix in small mode [FoldBackMat *]
 * Arg:        endj [WRITE] position of end in j (meaningless in i) [int *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int start_end_find_end_FoldBackMat(FoldBackMat * mat,int * endj) 
{
    register int j;  
    register int max;    
    register int maxj;   


    max = FoldBackMat_DC_SHADOW_SPECIAL(mat,0,mat->target->seq->len-1,END);  
    maxj = mat->target->seq->len-1;  
    for(j= mat->target->seq->len-2 ;j >= 0 ;j--) {  
      if( FoldBackMat_DC_SHADOW_SPECIAL(mat,0,j,END) > max ) {  
        max = FoldBackMat_DC_SHADOW_SPECIAL(mat,0,j,END);    
        maxj = j;    
        }  
      }  


    if( endj != NULL)    
      *endj = maxj;  


    return max;  
}    


/* Function:  dc_optimised_start_end_calc_FoldBackMat(*mat,dpenv)
 *
 * Descrip:    Calculates special strip, leaving start/end/score points in shadow matrix
 *             Works off specially laid out memory from steve searle
 *
 *
 * Arg:         *mat [UNKN ] Undocumented argument [FoldBackMat]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean dc_optimised_start_end_calc_FoldBackMat(FoldBackMat *mat,DPEnvelope * dpenv) 
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


    score_pointers = (int *) calloc (1 * (leni + 1) * 4,sizeof(int));    
    shadow_pointers = (int *) calloc (1 * (leni + 1) * 4 * 8,sizeof(int));   


    for(j=0;j<lenj;j++)  { /*for each j strip*/ 
      for(i=0;i<leni;i++)    { /*for each i position in strip*/ 
        num++;   
        if( dpenv != NULL && is_in_DPEnvelope(dpenv,i,j) == FALSE )  { /*Is not in envelope*/ 
          FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;    
          FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,STUTTER) = NEGI;  
          FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI;   
          FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI;   
          continue;  
          } /* end of Is not in envelope */ 
        if( num%1000 == 0)   
          log_full_error(REPORT,0,"%6d Cells done [%2d%%%%]",num,num*100/total); 




        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));     
        /* assign local shadown pointer */ 
        localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0));   
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,INSERT) + 0 +(DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,INSERT,0));    
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,DELETE) + 0 +(DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));  
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,DELETE,0));    
          }  
        /* From state START to state MATCH */ 
        temp = FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,i-1,j-1,START) + 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));     
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
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0 +(DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));     
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,STUTTER,0));   
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,MATCH) = score;    
        for(k=0;k<7;k++) 
          FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k) = localsp[k]; 
        /* Now figure out if any specials need this score */ 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->lgap) + (0) ;   
        if( temp > FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) )     {  
          FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,LOOP) = temp;    
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k); 
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,6) = FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);   
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,3) = i;  
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,4) = j;  
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,5) = MATCH;  
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) )  {  
          FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,i,j,END) = temp;     
          /* Have to push only bottem half of system here */ 
          for(k=0;k<3;k++)   
            FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,k) = FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,k);  
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,6) = FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,MATCH,6);    
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,3) = i;   
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,4) = j;   
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,END,5) = MATCH;   
          }  




        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,STUTTER) + 0 + (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));   
        /* assign local shadown pointer */ 
        localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,STUTTER,0)); 
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap +(DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,STUTTER,0));   
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap +(DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,STUTTER,0));   
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-1,MATCH) + 0 +(DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));   
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 1,MATCH,0)); 
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,STUTTER) = score;  
        for(k=0;k<7;k++) 
          FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,STUTTER,k) = localsp[k];   
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap + (0);   
        /* assign local shadown pointer */ 
        localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,MATCH,0));   
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-0,j-1,INSERT) + mat->qext +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 0,j - 1,INSERT,0));    
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,INSERT) = score;   
        for(k=0;k<7;k++) 
          FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,INSERT,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap + (0);   
        /* assign local shadown pointer */ 
        localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,MATCH,0));   
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i-1,j-0,DELETE) + mat->text +(0);    
        if( temp  > score )  {  
          score = temp;  
          /* assign local shadown pointer */ 
          localsp = &(FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i - 1,j - 0,DELETE,0));    
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
        /* Actually, already done inside scores */ 
         FoldBackMat_DC_OPT_SHADOW_MATRIX(mat,i,j,DELETE) = score;   
        for(k=0;k<7;k++) 
          FoldBackMat_DC_OPT_SHADOW_MATRIX_SP(mat,i,j,DELETE,k) = localsp[k];    
        /* Now figure out if any specials need this score */ 


        /* Finished calculating state DELETE */ 


        } /* end of for each i position in strip */ 


      /* Special state START has no special to special movements */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,0,j,LOOP);   


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,0,j - 1,LOOP) + (mat->lext) + (0);    
      if( temp > score ) {  
        score = temp;    
        /* Also got to propagate shadows  */ 
        for(k=0;k<7;k++) 
          FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i,j,LOOP,k) = FoldBackMat_DC_OPT_SHADOW_SPECIAL_SP(mat,i - 0,j - 1,LOOP,k);   
        }  


      /* Put back score... (now updated!) */ 
      FoldBackMat_DC_OPT_SHADOW_SPECIAL(mat,0,j,LOOP) = score;   
      /* Finished updating state LOOP */ 




      /* Special state END has no special to special movements */ 
      } /* end of for each j strip */ 
    free(score_pointers);    
    free(shadow_pointers);   
    return TRUE;     
}    


/* Function:  init_start_end_linear_FoldBackMat(mat)
 *
 * Descrip: No Description
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 */
void init_start_end_linear_FoldBackMat(FoldBackMat * mat) 
{
    register int i;  
    register int j;  
    for(j=0;j<3;j++) {  
      for(i=(-1);i<mat->query->seq->len;i++) {  
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,MATCH) = NEGI;  
        FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,MATCH,0) = (-1); 
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,STUTTER) = NEGI;    
        FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,STUTTER,0) = (-1);   
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,INSERT) = NEGI; 
        FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,INSERT,0) = (-1);    
        FoldBackMat_DC_SHADOW_MATRIX(mat,i,j,DELETE) = NEGI; 
        FoldBackMat_DC_SHADOW_MATRIX_SP(mat,i,j,DELETE,0) = (-1);    
        }  
      }  


    for(j=(-1);j<mat->target->seq->len;j++)  {  
      FoldBackMat_DC_SHADOW_SPECIAL(mat,0,j,START) = 0;  
      FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,j,START,0) = j; 
      FoldBackMat_DC_SHADOW_SPECIAL(mat,0,j,LOOP) = NEGI;    
      FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,j,LOOP,0) = (-1);   
      FoldBackMat_DC_SHADOW_SPECIAL(mat,0,j,END) = NEGI; 
      FoldBackMat_DC_SHADOW_SPECIAL_SP(mat,0,j,END,0) = (-1);    
      }  


    return;  
}    


/* Function:  convert_PackAln_to_AlnBlock_FoldBackMat(pal)
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
AlnBlock * convert_PackAln_to_AlnBlock_FoldBackMat(PackAln * pal) 
{
    AlnConvertSet * acs; 
    AlnBlock * alb;  


    acs = AlnConvertSet_FoldBackMat();   
    alb = AlnBlock_from_PackAln(acs,pal);    
    free_AlnConvertSet(acs); 
    return alb;  
}    


 static char * query_label[] = { "SEQUENCE","INSERT","END" };    
/* Function:  AlnConvertSet_FoldBackMat(void)
 *
 * Descrip: No Description
 *
 *
 * Return [UNKN ]  Undocumented return value [AlnConvertSet *]
 *
 */
 static char * target_label[] = { "SEQUENCE","INSERT","LOOP","END" };    
AlnConvertSet * AlnConvertSet_FoldBackMat(void) 
{
    AlnConvertUnit * acu;    
    AlnConvertSet  * out;    


    out = AlnConvertSet_alloc_std(); 


    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = START + 4; 
    acu->is_from_special = TRUE; 
    acu->state2 = MATCH;     
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = STUTTER;   
    acu->state2 = MATCH;     
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = STUTTER;   
    acu->state2 = STUTTER;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = STUTTER;   
    acu->state2 = STUTTER;   
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = STUTTER;   
    acu->state2 = STUTTER;   
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = STUTTER;   
    acu->offi = 1;   
    acu->offj = 1;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = INSERT;    
    acu->state2 = INSERT;    
    acu->offi = 0;   
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = DELETE;    
    acu->state2 = DELETE;    
    acu->offi = 1;   
    acu->offj = 0;   
    acu->label1 = query_label[0];    
    acu->label2 = target_label[1];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = LOOP + 4;  
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[2];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = LOOP + 4;  
    acu->state2 = LOOP + 4;  
    acu->offi = (-1);    
    acu->offj = 1;   
    acu->label1 = query_label[1];    
    acu->label2 = target_label[0];   
    acu = AlnConvertUnit_alloc();    
    add_AlnConvertSet(out,acu);  
    acu->state1 = MATCH; 
    acu->state2 = END + 4;   
    acu->offi = (-1);    
    acu->offj = 0;   
    acu->label1 = query_label[2];    
    acu->label2 = target_label[3];   
    return out;  
}    


/* Function:  PackAln_read_Expl_FoldBackMat(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_Expl_FoldBackMat(FoldBackMat * mat) 
{
    FoldBackMat_access_func_holder holder;   


    holder.access_main    = FoldBackMat_explicit_access_main;    
    holder.access_special = FoldBackMat_explicit_access_special; 
    return PackAln_read_generic_FoldBackMat(mat,holder); 
}    


/* Function:  FoldBackMat_explicit_access_main(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int FoldBackMat_explicit_access_main(FoldBackMat * mat,int i,int j,int state) 
{
    return FoldBackMat_EXPL_MATRIX(mat,i,j,state);   
}    


/* Function:  FoldBackMat_explicit_access_special(mat,i,j,state)
 *
 * Descrip: No Description
 *
 * Arg:          mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:            i [UNKN ] Undocumented argument [int]
 * Arg:            j [UNKN ] Undocumented argument [int]
 * Arg:        state [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int FoldBackMat_explicit_access_special(FoldBackMat * mat,int i,int j,int state) 
{
    return FoldBackMat_EXPL_SPECIAL(mat,i,j,state);  
}    


/* Function:  PackAln_read_generic_FoldBackMat(mat,h)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:          h [UNKN ] Undocumented argument [FoldBackMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * PackAln_read_generic_FoldBackMat(FoldBackMat * mat,FoldBackMat_access_func_holder h) 
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


    out->score =  find_end_FoldBackMat(mat,&i,&j,&state,&isspecial,h);   


    /* Add final end transition (at the moment we have not got the score! */ 
    if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE )   {  
      warn("Failed the first PackAlnUnit alloc, %d length of Alignment in FoldBackMat_basic_read, returning a mess.(Sorry!)",out->len);  
      return out;    
      }  


    /* Put in positions for end trans. Remember that coordinates in C style */ 
    pau->i = i;  
    pau->j = j;  
    if( isspecial != TRUE)   
      pau->state = state;    
    else pau->state = state + 4;     
    prev=pau;    
    while( state != START || isspecial != TRUE)  { /*while state != START*/ 


      if( isspecial == TRUE )    
        max_calc_special_FoldBackMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);  
      else   
        max_calc_FoldBackMat(mat,i,j,state,isspecial,&i,&j,&state,&isspecial,&cellscore,h);  
      if(i == FoldBackMat_READ_OFF_ERROR || j == FoldBackMat_READ_OFF_ERROR || state == FoldBackMat_READ_OFF_ERROR ) {  
        warn("Problem - hit bad read off system, exiting now");  
        break;   
        }  
      if( (pau= PackAlnUnit_alloc()) == NULL  || add_PackAln(out,pau) == FALSE ) {  
        warn("Failed a PackAlnUnit alloc, %d length of Alignment in FoldBackMat_basic_read, returning partial alignment",out->len);  
        break;   
        }  


      /* Put in positions for block. Remember that coordinates in C style */ 
      pau->i = i;    
      pau->j = j;    
      if( isspecial != TRUE)     
        pau->state = state;  
      else pau->state = state + 4;   
      prev->score = cellscore;   
      prev = pau;    
      } /* end of while state != START */ 


    invert_PackAln(out); 
    return out;  
}    


/* Function:  find_end_FoldBackMat(mat,ri,rj,state,isspecial,h)
 *
 * Descrip: No Description
 *
 * Arg:              mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:               ri [UNKN ] Undocumented argument [int *]
 * Arg:               rj [UNKN ] Undocumented argument [int *]
 * Arg:            state [UNKN ] Undocumented argument [int *]
 * Arg:        isspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:                h [UNKN ] Undocumented argument [FoldBackMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int find_end_FoldBackMat(FoldBackMat * mat,int * ri,int * rj,int * state,boolean * isspecial,FoldBackMat_access_func_holder h) 
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


/* Function:  FoldBackMat_debug_show_matrix(mat,starti,stopi,startj,stopj,ofp)
 *
 * Descrip: No Description
 *
 * Arg:           mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:        starti [UNKN ] Undocumented argument [int]
 * Arg:         stopi [UNKN ] Undocumented argument [int]
 * Arg:        startj [UNKN ] Undocumented argument [int]
 * Arg:         stopj [UNKN ] Undocumented argument [int]
 * Arg:           ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void FoldBackMat_debug_show_matrix(FoldBackMat * mat,int starti,int stopi,int startj,int stopj,FILE * ofp) 
{
    register int i;  
    register int j;  


    for(i=starti;i<stopi && i < mat->query->seq->len;i++)    {  
      for(j=startj;j<stopj && j < mat->target->seq->len;j++) {  
        fprintf(ofp,"Cell [%d - %d]\n",i,j);     
        fprintf(ofp,"State MATCH %d\n",FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH));  
        fprintf(ofp,"State STUTTER %d\n",FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER));  
        fprintf(ofp,"State INSERT %d\n",FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT));    
        fprintf(ofp,"State DELETE %d\n",FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE));    
        fprintf(ofp,"\n\n"); 
        }  
      }  


}    


/* Function:  max_calc_FoldBackMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [FoldBackMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_FoldBackMat(FoldBackMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,FoldBackMat_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = FoldBackMat_READ_OFF_ERROR;  


    if( i < 0 || j < 0 || i > mat->query->seq->len || j > mat->target->seq->len) {  
      warn("In FoldBackMat matrix special read off - out of bounds on matrix [i,j is %d,%d state %d in standard matrix]",i,j,state); 
      return -1;     
      }  


    /* Then you have to select the correct switch statement to figure out the readoff      */ 
    /* Somewhat odd - reverse the order of calculation and return as soon as it is correct */ 
    cscore = (*h.access_main)(mat,i,j,state);    
    switch(state)    { /*Switch state */ 
      case MATCH :   
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,STUTTER) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,STUTTER); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,STUTTER);  
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
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
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,DELETE);   
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,INSERT) )   {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,INSERT);  
            }  
          return (*h.access_main)(mat,i - 1,j - 1,INSERT);   
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH);    
          }  
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case STUTTER :     
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,MATCH) )    {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = MATCH; 
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,MATCH);   
            }  
          return (*h.access_main)(mat,i - 1,j - 1,MATCH);    
          }  
        temp = cscore - (mat->stgap) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));   
        if( temp == (*h.access_main)(mat,i - 1,j - 0,STUTTER) )  {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,STUTTER); 
            }  
          return (*h.access_main)(mat,i - 1,j - 0,STUTTER);  
          }  
        temp = cscore - (mat->stgap) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));   
        if( temp == (*h.access_main)(mat,i - 0,j - 1,STUTTER) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,STUTTER); 
            }  
          return (*h.access_main)(mat,i - 0,j - 1,STUTTER);  
          }  
        temp = cscore - (0) -  (DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j)));    
        if( temp == (*h.access_main)(mat,i - 1,j - 1,STUTTER) )  {  
          *reti = i - 1; 
          *retj = j - 1; 
          *retstate = STUTTER;   
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-1,STUTTER); 
            }  
          return (*h.access_main)(mat,i - 1,j - 1,STUTTER);  
          }  
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case INSERT :  
        temp = cscore - (mat->qext) -  (0);  
        if( temp == (*h.access_main)(mat,i - 0,j - 1,INSERT) )   {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = INSERT;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-0,j-1,INSERT);  
            }  
          return (*h.access_main)(mat,i - 0,j - 1,INSERT);   
          }  
        temp = cscore - (mat->qgap) -  (0);  
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
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      case DELETE :  
        temp = cscore - (mat->text) -  (0);  
        if( temp == (*h.access_main)(mat,i - 1,j - 0,DELETE) )   {  
          *reti = i - 1; 
          *retj = j - 0; 
          *retstate = DELETE;    
          *retspecial = FALSE;   
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_main)(mat,i-1,j-0,DELETE);  
            }  
          return (*h.access_main)(mat,i - 1,j - 0,DELETE);   
          }  
        temp = cscore - (mat->tgap) -  (0);  
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
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      default:   
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found!",i,j,state); 
        return (-1); 
      } /* end of Switch state  */ 
}    


/* Function:  max_calc_special_FoldBackMat(mat,i,j,state,isspecial,reti,retj,retstate,retspecial,cellscore,h)
 *
 * Descrip: No Description
 *
 * Arg:               mat [UNKN ] Undocumented argument [FoldBackMat *]
 * Arg:                 i [UNKN ] Undocumented argument [int]
 * Arg:                 j [UNKN ] Undocumented argument [int]
 * Arg:             state [UNKN ] Undocumented argument [int]
 * Arg:         isspecial [UNKN ] Undocumented argument [boolean]
 * Arg:              reti [UNKN ] Undocumented argument [int *]
 * Arg:              retj [UNKN ] Undocumented argument [int *]
 * Arg:          retstate [UNKN ] Undocumented argument [int *]
 * Arg:        retspecial [UNKN ] Undocumented argument [boolean *]
 * Arg:         cellscore [UNKN ] Undocumented argument [int *]
 * Arg:                 h [UNKN ] Undocumented argument [FoldBackMat_access_func_holder]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int max_calc_special_FoldBackMat(FoldBackMat * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore,FoldBackMat_access_func_holder h) 
{
    register int temp;   
    register int cscore; 


    *reti = (*retj) = (*retstate) = FoldBackMat_READ_OFF_ERROR;  


    if( j < 0 || j > mat->target->seq->len)  {  
      warn("In FoldBackMat matrix special read off - out of bounds on matrix [j is %d in special]",j);   
      return -1;     
      }  


    cscore = (*h.access_special)(mat,i,j,state); 
    switch(state)    { /*switch on special states*/ 
      case START :   
      case LOOP :    
        /* source LOOP is a special */ 
        temp = cscore - (mat->lext) - (0);   
        if( temp == (*h.access_special)(mat,i - 0,j - 1,LOOP) )  {  
          *reti = i - 0; 
          *retj = j - 1; 
          *retstate = LOOP;  
          *retspecial = TRUE;    
          if( cellscore != NULL) {  
            *cellscore = cscore - (*h.access_special)(mat,i-0,j-1,LOOP);     
            }  
          return (*h.access_special)(mat,i - 0,j - 1,LOOP) ;     
          }  
        /* source MATCH is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
          temp = cscore - (mat->lgap) - (0);     
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
      case END :     
        /* source MATCH is from main matrix */ 
        for(i= mat->query->seq->len-1;i >= 0 ;i--)   { /*for i >= 0*/ 
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
        warn("Major problem (!) - in FoldBackMat read off, position %d,%d state %d no source found  dropped into default on source switch!",i,j,state);  
        return (-1); 
      } /* end of switch on special states */ 
}    


/* Function:  calculate_FoldBackMat(mat)
 *
 * Descrip:    This function calculates the FoldBackMat matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_FoldBackMat
 *
 *
 * Arg:        mat [UNKN ] FoldBackMat which contains explicit basematrix memory [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_FoldBackMat(FoldBackMat * mat) 
{
    int i;   
    int j;   
    int leni;    
    int lenj;    
    int tot; 
    int num; 


    if( mat->basematrix->type != BASEMATRIX_TYPE_EXPLICIT )  {  
      warn("in calculate_FoldBackMat, passed a non Explicit matrix type, cannot calculate!");    
      return FALSE;  
      }  


    leni = mat->leni;    
    lenj = mat->lenj;    
    tot = leni * lenj;   
    num = 0; 


    start_reporting("FoldBackMat Matrix calculation: "); 
    for(j=0;j<lenj;j++)  {  
      auto int score;    
      auto int temp;     
      for(i=0;i<leni;i++)    {  
        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,INSERT) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,DELETE) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = FoldBackMat_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,STUTTER) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH) = score; 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->lgap) + (0) ;   
        if( temp > FoldBackMat_EXPL_SPECIAL(mat,i,j,LOOP) )  {  
          FoldBackMat_EXPL_SPECIAL(mat,i,j,LOOP) = temp;     
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > FoldBackMat_EXPL_SPECIAL(mat,i,j,END) )   {  
          FoldBackMat_EXPL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,STUTTER) + 0;    
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER) = score;   


        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;  
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT) = score;    


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;  
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->text;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE) = score;    


        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = FoldBackMat_EXPL_SPECIAL(mat,0,j,LOOP);    


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = FoldBackMat_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (mat->lext) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      FoldBackMat_EXPL_SPECIAL(mat,0,j,LOOP) = score;    
      /* Finished updating state LOOP */ 




      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  calculate_dpenv_FoldBackMat(mat,dpenv)
 *
 * Descrip:    This function calculates the FoldBackMat matrix when in explicit mode, subject to the envelope
 *
 *
 * Arg:          mat [UNKN ] FoldBackMat which contains explicit basematrix memory [FoldBackMat *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean calculate_dpenv_FoldBackMat(FoldBackMat * mat,DPEnvelope * dpenv) 
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
      warn("in calculate_FoldBackMat, passed a non Explicit matrix type, cannot calculate!");    
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
        FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI;   
        FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER) = NEGI; 
        FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;  
        FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;  
        }  
      }  
    for(j=-1;j<mat->lenj;j++)    {  
      FoldBackMat_EXPL_SPECIAL(mat,i,j,START) = 0;   
      FoldBackMat_EXPL_SPECIAL(mat,i,j,LOOP) = NEGI; 
      FoldBackMat_EXPL_SPECIAL(mat,i,j,END) = NEGI;  
      }  


    start_reporting("FoldBackMat Matrix calculation: "); 
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
          FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH) = NEGI; 
          FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER) = NEGI;   
          FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT) = NEGI;    
          FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE) = NEGI;    
          continue;  
          }  


        if( num%1000 == 0 )  
          log_full_error(REPORT,0,"[%7d] Cells %2d%%%%",num,num*100/tot);    
        num++;   


        /* For state MATCH */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + 0;  
        /* From state INSERT to state MATCH */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,INSERT) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state DELETE to state MATCH */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,DELETE) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state START to state MATCH */ 
        temp = FoldBackMat_EXPL_SPECIAL(mat,i-1,j-1,START) + 0;  
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state MATCH */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,STUTTER) + 0;     
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for MATCH */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_EXPL_MATRIX(mat,i,j,MATCH) = score; 


        /* state MATCH is a source for special LOOP */ 
        temp = score + (mat->lgap) + (0) ;   
        if( temp > FoldBackMat_EXPL_SPECIAL(mat,i,j,LOOP) )  {  
          FoldBackMat_EXPL_SPECIAL(mat,i,j,LOOP) = temp;     
          }  




        /* state MATCH is a source for special END */ 
        temp = score + (0) + (0) ;   
        if( temp > FoldBackMat_EXPL_SPECIAL(mat,i,j,END) )   {  
          FoldBackMat_EXPL_SPECIAL(mat,i,j,END) = temp;  
          }  




        /* Finished calculating state MATCH */ 


        /* For state STUTTER */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,STUTTER) + 0;    
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-0,j-1,STUTTER) + mat->stgap;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state STUTTER to state STUTTER */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-0,STUTTER) + mat->stgap;    
        if( temp  > score )  {  
          score = temp;  
          }  
        /* From state MATCH to state STUTTER */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-1,MATCH) + 0;   
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for STUTTER */ 
        /* Add any movement independant score and put away */ 
         score += DnaMatrix_MATCH(mat->comp,CSEQ_DNA_BASE(mat->query,i),CSEQ_DNA_BASE(mat->target,j));   
         FoldBackMat_EXPL_MATRIX(mat,i,j,STUTTER) = score;   


        /* Finished calculating state STUTTER */ 


        /* For state INSERT */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-0,j-1,MATCH) + mat->qgap;  
        /* From state INSERT to state INSERT */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-0,j-1,INSERT) + mat->qext;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for INSERT */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_EXPL_MATRIX(mat,i,j,INSERT) = score;    


        /* Finished calculating state INSERT */ 


        /* For state DELETE */ 
        /* setting first movement to score */ 
        score = FoldBackMat_EXPL_MATRIX(mat,i-1,j-0,MATCH) + mat->tgap;  
        /* From state DELETE to state DELETE */ 
        temp = FoldBackMat_EXPL_MATRIX(mat,i-1,j-0,DELETE) + mat->text;  
        if( temp  > score )  {  
          score = temp;  
          }  


        /* Ok - finished max calculation for DELETE */ 
        /* Add any movement independant score and put away */ 
         FoldBackMat_EXPL_MATRIX(mat,i,j,DELETE) = score;    


        /* Finished calculating state DELETE */ 
        }  


      /* Special state START has no special to special movements */ 


      /* Special state LOOP has special to speical */ 
      /* Set score to current score (remember, state probably updated during main loop */ 
      score = FoldBackMat_EXPL_SPECIAL(mat,0,j,LOOP);    


      /* Source MATCH for state LOOP is not special... already calculated */ 
      /* Source LOOP is a special source for LOOP */ 
      temp = FoldBackMat_EXPL_SPECIAL(mat,0,j - 1,LOOP) + (mat->lext) + (0);     
      if( temp > score ) 
        score = temp;    


      /* Put back score... (now updated!) */ 
      FoldBackMat_EXPL_SPECIAL(mat,0,j,LOOP) = score;    
      /* Finished updating state LOOP */ 




      /* Special state END has no special to special movements */ 
      }  
    stop_reporting();    
    return TRUE;     
}    


/* Function:  FoldBackMat_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * FoldBackMat_alloc(void) 
{
    FoldBackMat * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(FoldBackMat *) ckalloc (sizeof(FoldBackMat))) == NULL)  {  
      warn("FoldBackMat_alloc failed "); 
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


/* Function:  free_FoldBackMat(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [FoldBackMat *]
 *
 * Return [UNKN ]  Undocumented return value [FoldBackMat *]
 *
 */
FoldBackMat * free_FoldBackMat(FoldBackMat * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a FoldBackMat obj. Should be trappable");   
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
    /* obj->comp is linked in */ 
    /* obj->stutter is linked in */ 
    /* obj->qgap is linked in */ 
    /* obj->qext is linked in */ 
    /* obj->tgap is linked in */ 
    /* obj->text is linked in */ 
    /* obj->lgap is linked in */ 
    /* obj->lext is linked in */ 
    /* obj->stgap is linked in */ 


    ckfree(obj); 
    return NULL; 
}    





#ifdef _cplusplus
}
#endif
