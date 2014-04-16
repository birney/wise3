
#include "peak.h"



int main(int argc,char ** argv)
{
  PeakList * pl;

  pl = read_bed_PeakList(stdin,"test");

  write_4col_bed_PeakList(pl,stdout);

}
