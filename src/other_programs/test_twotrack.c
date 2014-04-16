#include "two_track.h"



int main(int argc,char ** argv) {
  
  TwoTrackSet * set;
  TwoTrackSetStats * stats;

  set = read_two_track_Tim_style_report(argv[1],argv[2],2000,20000);

  assert(set != NULL);

  stats = TwoTrackSetStats_from_TwoTrackSet(set);
  
  write_TwoTrackSetStats(stats,stdout);


  exit(0);
}

