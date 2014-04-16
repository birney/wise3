#ifdef _cplusplus
extern "C" {
#endif
#include "dpvector.h"
#include "dynafunc.h"


/* Function:  write_naive_vector(dfp,*dpi,gm)
 *
 * Descrip:    writes out naive vectorised loop
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 19 "dpvector.dy"
boolean write_naive_vector(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm)
{

  write_naive_vector_alloc(dfp,dpi,gm);

  return TRUE;
}



/* Function:  write_naive_vector_calc(dfp,*dpi,gm)
 *
 * Descrip:    writes out an alloc+calc structure for naive vectorisation calc lines
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 32 "dpvector.dy"
boolean write_naive_vector_calc(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm)
{
  FuncInfo * fi;
  int i;
  int j;

  char * arg_str;


  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"create_naive_vector_calc_%s",gm->name);
  add_line_to_Ftext(fi->ft,"creating and restructuring calc lines %s",gm->name);

  arg_str = get_argstr_GenericMatrix(gm);
  add_args_GenericMatrix_FuncInfo(fi,gm);

  start_function_FuncInfo(fi,dfp,"struct naive_vector_calc_%s * create_naive_vector_calc_%s(%s)",gm->name,gm->name,arg_str);
  expr(dfp,"struct naive_vector_calc_%s * out;",gm->name);
  expr(dfp,"int total_depth;");
  expr(dfp,"int j;");
  expr(dfp,"int amino_acid;");

  add_break(dfp);

  expr(dfp,"total_depth = %d * (int)(1+floor((%s->%s / %d))) ;",dpi->vector_depth,gm->query->name,gm->query_len,dpi->vector_depth);
  
  expr(dfp,"out = malloc(sizeof(struct naive_vector_calc_%s));",gm->name);
  for(i=0;i<gm->len;i++) {
    expr(dfp,"out->state_%s_source_%s = calloc(26,sizeof(int *));",gm->state[i]->name,gm->state[i]->source[j]->state_source);
    
  }

 
  
  

}



/* Function:  write_naive_vector_alloc(dfp,*dpi,gm)
 *
 * Descrip:    writes out an allocation structure for naive vectorisation
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 74 "dpvector.dy"
boolean write_naive_vector_alloc(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm)
{
  FuncInfo * fi;
  char * arg_str;


  fi = FuncInfo_named_from_varstr(FI_INTERNAL,"alloc_naive_vector_%s",gm->name);
  add_line_to_Ftext(fi->ft,"allocating for naive vectorisation %s",gm->name);

  arg_str = get_argstr_GenericMatrix(gm);
  add_args_GenericMatrix_FuncInfo(fi,gm);

  start_function_FuncInfo(fi,dfp,"struct naive_vector_%s * alloc_naive_vector_%s(%s)",gm->name,gm->name,arg_str);
  expr(dfp,"struct naive_vector_%s * out;",gm->name);
  expr(dfp,"int total_depth;");
  expr(dfp,"int j;");

  add_break(dfp);
  expr(dfp,"out->shift = calloc(%d,sizeof(int*));",gm->window_j+1);
  expr(dfp,"out->non_shift = calloc(%d,sizeof(int*));",gm->window_j+1);
  add_break(dfp);

  expr(dfp,"total_depth = %d * (int)(1+floor((%s->%s / %d))) ;",dpi->vector_depth,gm->query->name,gm->query_len,dpi->vector_depth);
  expr(dfp,"for(j=0;j<%d;j++)");
  startbrace(dfp);
  expr(dfp,"out->shift[j] = calloc(total_depth*%d,sizeof(int));",gm->len);
  expr(dfp,"out->non_shift[j] = calloc(total_depth*%d,sizeof(int));",gm->len);
  closebrace(dfp);
  
  
  expr(dfp,"return out");
  close_function(dfp);
  add_break(dfp);
}

/* Function:  write_naive_vector_calc_struct(dfp,dpi,gm)
 *
 * Descrip:    writes out a calc structure for calc layouts
 *
 *
 * Arg:        dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        dpi [UNKN ] Undocumented argument [DPImplementation *]
 * Arg:         gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 112 "dpvector.dy"
boolean write_naive_vector_calc_struct(DYNFILE * dfp,DPImplementation * dpi,GenericMatrix * gm)
{
  int i;
  int j;

  start_struct(dfp,"struct naive_vector_calc_%s",gm->name);
  for(i=0;i<gm->len;i++) {
    for(j=0;j<gm->state[i]->len;j++) {
      struct_expr(dfp,"int ** state_%s_source_%s;",gm->state[i]->name,gm->state[i]->source[j]->state_source);
    }
  }

  close_struct(dfp,";");

}


/* Function:  write_naive_vector_struct(dfp,*dpi,gm)
 *
 * Descrip:    writes out vector layout structure for naive vectorisation
 *
 *
 * Arg:         dfp [UNKN ] Undocumented argument [DYNFILE *]
 * Arg:        *dpi [UNKN ] Undocumented argument [DPImplementation]
 * Arg:          gm [UNKN ] Undocumented argument [GenericMatrix *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
# line 132 "dpvector.dy"
boolean write_naive_vector_struct(DYNFILE * dfp,DPImplementation *dpi,GenericMatrix * gm)
{
  int i;
  
  start_struct(dfp,"struct naive_vector_%s",gm->name);
  struct_expr(dfp,"%s %s;",gm->query->element_type,gm->query->name);
  struct_expr(dfp,"%s %s;",gm->target->element_type,gm->target->name);
  
  for(i=0;i<gm->res_len;i++) {
    struct_expr(dfp,"%s %s",gm->resource[i]->element_type,gm->resource[i]->name);
  }

  struct_expr(dfp,"int ** non_shift;");
  struct_expr(dfp,"int ** shift;");
  
  close_struct(dfp,";");

  write_naive_vector_calc_struct(dfp,dpi,gm);
}





# line 187 "dpvector.c"

#ifdef _cplusplus
}
#endif
