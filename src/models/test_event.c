

#include "dualsignal.h"



int main(int argc,char **argv)
{
  SignalEventList * sel;

  sel = read_SignalEventList(stdin);


  write_SignalEventList(sel,stdout);
    


}
