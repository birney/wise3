

.SUFFIXES :  .dy

CC = cc
CFLAGS =   -g3 -c -pthread  -DUNIX  -I../base/ -I../dynlibsrc `glib-config --cflags`
LFLAGS =  -g -L../base/ -L../dynlibsrc  -lm `glib-config --libs`
AR_OPTIONS = ru

RANLIB_NEEDED = 0

test_peak : test_peak.o peak.o
	$(CC) -o test_peak test_peak.o peak.o $(LFLAGS) -ldyna -lwisebase -lm


merge_peaks : merge_peaks.o peak.o
	$(CC) -o merge_peaks merge_peaks.o peak.o $(LFLAGS) -ldyna -lwisebase -lm

test_twotrack : test_twotrack.o two_track.o 
	$(CC) -o test_twotrack test_twotrack.o two_track.o $(LFLAGS) -ldyna -lwisebase -lm

test_locus : test_locus.o snplocus.o
	$(CC) -o test_locus test_locus.o snplocus.o $(LFLAGS) -ldyna -lwisebase -lm

mouse_hmm : mouse_hmm.o snplocus.o mousethreestate.o ancestral.o plot.o
	$(CC) -o mouse_hmm mouse_hmm.o snplocus.o ancestral.o plot.o mousethreestate.o $(LFLAGS) -ldyna -lwisebase -lm


mouse_anc : mouse_anc.o snplocus.o ancestral.o mousethreestate.o plot.o 
	$(CC) -o mouse_anc mouse_anc.o ancestral.o snplocus.o plot.o mousethreestate.o $(LFLAGS) -ldyna -lwisebase -lm


mouse_pair : mouse_pair.o snplocus.o ancestral.o ibd_model.o plot.o 
	$(CC) -o mouse_pair mouse_pair.o ancestral.o snplocus.o plot.o ibd_model.o $(LFLAGS) -ldyna -lwisebase -lm

.c.o :
	$(CC) $(CFLAGS) $(INCFLAGS) $?

DFLAGS = -l -n Wise2_ -a _api.h -b _api.t -latex -perl -prob

.dy.c : 
	dyc  $(DFLAGS) $<	

.dy.o :
	dyc $(DFLAGS) $<	
	$(CC) $(CFLAGS) $*.c

init:
	dyc $(DFLAGS) *.dy
