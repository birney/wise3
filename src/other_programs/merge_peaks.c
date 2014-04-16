
#include "peak.h"



int main(int argc,char ** argv)
{
  ListofPeakList * lop;
  PeakList * pl;
  MergedPeakList * mpl;
  int i;
  FILE * ifp;
  
  int j;

  lop = ListofPeakList_alloc_std();

  for(i=1;i<argc;i++) {
    
    fprintf(stderr,"Reading %s\n",argv[i]);

    ifp = openfile(argv[i],"r");

    for(j=0;argv[i][j];j++) {
      if( argv[i][j] == '.' ) {
	argv[i][j] = '\0';
	break;
      }
    }

    pl = read_npf_PeakList(ifp,argv[i]);

    add_ListofPeakList(lop,pl);
  }

  mpl = make_MergedPeakList(lop,1);
  
  write_tagged_MergedPeakList(mpl,stdout);

}
