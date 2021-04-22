#   SIGMA-THETA Project - Makefile for the binaries of the SIGMA-THETA set :
#	    - 1col2col
#           - X2Y
#           - DriRem
#           - SigmaTheta 
#           - ADev
#           - GCoDev
#           - MDev
#	    - HDev
#	    - PDev
#	    - Aver
#           - uncertainties
#           - RaylConfInt
#           - Asymptote
#	    - Asym2Alpha
#           - AVarDOF
#           - ADUncert
#           - ADGraph
#           - bruiteur
#
#                                                    SIGMA-THETA Project
# 
# Copyright or © or Copr. Université de Franche-Comté, Besançon, France
# Contributor: François Vernotte, UTINAM/OSU THETA (2012/07/17)
# Contact: francois.vernotte@obs-besancon.fr
# 
# This software, SigmaTheta, is a collection of computer programs for
# time and frequency metrology. 
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
# 
# 

BIN = ./
OBJ = obj/
SOURCE = source/
INSTALLDIR = /usr/local/bin/
INSTALL = install
CC = gcc
GIT_VERSION := "$(shell git describe --abbrev=4 --dirty --always --tags)"
GIT_DATE := "$(shell git log | head -3 | grep Date | sed 's/Date:   //')"
CFLAGS = -g -O3 -DST_VERSION=\"$(GIT_VERSION)\" -DST_DATE=\"$(GIT_DATE)\"
# link statically against libfftw3
#FFTW3 = -static -lfftw3
# link dynamically against libfftw3
FFTW3 = -lfftw3

TARGETS = $(BIN)1col2col $(BIN)X2Y $(BIN)DriRem $(BIN)SigmaTheta $(BIN)ADev $(BIN)find_tau_xdev $(BIN)GCoDev $(BIN)MDev $(BIN)HDev $(BIN)PDev $(BIN)Aver $(BIN)uncertainties $(BIN)RaylConfInt $(BIN)Asymptote $(BIN)Asym2Alpha $(BIN)AVarDOF $(BIN)ADUncert $(BIN)ADGraph $(BIN)PSDGraph $(BIN)YkGraph $(BIN)XtGraph $(BIN)DevGraph $(BIN)3CHGraph $(BIN)bruiteur $(BIN)GCUncert $(BIN)3CorneredHat

all: $(TARGETS)

clean:
	rm -f $(TARGETS)
	rm -f $(OBJ)/*.o

$(BIN)1col2col : $(OBJ)1col2col.o $(OBJ)stio_sbr.o  
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)X2Y : $(OBJ)x2y.o $(OBJ)stio_sbr.o  
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)DriRem : $(OBJ)drirem.o $(OBJ)tchebyfit.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)SigmaTheta : $(OBJ)sigma_theta.o $(OBJ)dev_sbr.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)ADev : $(OBJ)adev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)find_tau_xdev : $(OBJ)find_tau_xdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)GCoDev : $(OBJ)gcodev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)MDev : $(OBJ)mdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)HDev : $(OBJ)hdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)PDev : $(OBJ)pdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)Aver : $(OBJ)aver.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -Wl,-Bdynamic -lgsl -lgslcblas -lm

$(BIN)uncertainties : $(OBJ)uncertainties.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)RaylConfInt : $(OBJ)rayleigh.o $(OBJ)rayleigh_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)Asymptote : $(OBJ)asymptote.o $(OBJ)asymptote_sbr.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)Asym2Alpha : $(OBJ)asym2alpha.o $(OBJ)asymptote_sbr.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)AVarDOF : $(OBJ)avardof.o $(OBJ)avardof_sbr.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)ADUncert : $(OBJ)aduncert.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)ADGraph : $(OBJ)adgraph.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)PSDGraph : $(OBJ)psdgraph.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)YkGraph : $(OBJ)ykgraph.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)XtGraph : $(OBJ)xtgraph.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)DevGraph : $(OBJ)devgraph.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)3CHGraph : $(OBJ)3chgraph.o $(OBJ)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)bruiteur : $(OBJ)bruiteur.o $(OBJ)filtre.o $(OBJ)splitmix64.o $(OBJ)xorshift1024star.o $(OBJ)ziggurat.o
	$(CC) $(CFLAGS) -o $@ $^ $(FFTW3) -Wl,-Bdynamic -lm

$(BIN)GCUncert : $(OBJ)gcuncert.o $(OBJ)stio_sbr.o $(OBJ)filtre.o $(OBJ)3ch_sbr.o $(OBJ)splitmix64.o $(OBJ)xorshift1024star.o $(OBJ)ziggurat.o
	$(CC) $(CFLAGS) -o $@ $^ $(FFTW3) -Wl,-Bdynamic -lgsl -lgslcblas -lm

$(BIN)3CorneredHat : $(OBJ)3_cornered_hat.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)stio_sbr.o  $(OBJ)dev_sbr.o $(OBJ)tchebyfit.o $(OBJ)filtre.o $(OBJ)3ch_sbr.o $(OBJ)splitmix64.o $(OBJ)xorshift1024star.o $(OBJ)ziggurat.o
	$(CC) $(CFLAGS) -o $@ $^ $(FFTW3) -Wl,-Bdynamic -lgsl -lgslcblas -lm

$(OBJ)1col2col.o : $(SOURCE)1col2col.c $(SOURCE)sigma_theta.h

$(OBJ)x2y.o : $(SOURCE)x2y.c $(SOURCE)sigma_theta.h

$(OBJ)drirem.o : $(SOURCE)drirem.c $(SOURCE)sigma_theta.h

$(OBJ)sigma_theta.o : $(SOURCE)sigma_theta.c $(SOURCE)sigma_theta.h

$(OBJ)adev.o : $(SOURCE)adev.c $(SOURCE)sigma_theta.h

$(OBJ)gcodev.o : $(SOURCE)gcodev.c $(SOURCE)sigma_theta.h

$(OBJ)mdev.o : $(SOURCE)mdev.c $(SOURCE)sigma_theta.h

$(OBJ)hdev.o : $(SOURCE)hdev.c $(SOURCE)sigma_theta.h

$(OBJ)pdev.o : $(SOURCE)pdev.c $(SOURCE)sigma_theta.h

$(OBJ)aver.o : $(SOURCE)aver.c $(SOURCE)sigma_theta.h

$(OBJ)uncertainties.o : $(SOURCE)uncertainties.c $(SOURCE)sigma_theta.h

$(OBJ)rayleigh.o : $(SOURCE)rayleigh.c $(SOURCE)sigma_theta.h

$(OBJ)asymptote.o : $(SOURCE)asymptote.c $(SOURCE)sigma_theta.h

$(OBJ)asym2alpha.o : $(SOURCE)asym2alpha.c $(SOURCE)sigma_theta.h

$(OBJ)avardof.o : $(SOURCE)avardof.c $(SOURCE)sigma_theta.h

$(OBJ)aduncert.o : $(SOURCE)aduncert.c $(SOURCE)sigma_theta.h

$(OBJ)adgraph.o : $(SOURCE)adgraph.c $(SOURCE)sigma_theta.h

$(OBJ)psdgraph.o : $(SOURCE)psdgraph.c $(SOURCE)sigma_theta.h

$(OBJ)ykgraph.o : $(SOURCE)ykgraph.c $(SOURCE)sigma_theta.h

$(OBJ)xtgraph.o : $(SOURCE)xtgraph.c $(SOURCE)sigma_theta.h

$(OBJ)devgraph.o : $(SOURCE)devgraph.c $(SOURCE)sigma_theta.h

$(OBJ)3chgraph.o : $(SOURCE)3chgraph.c $(SOURCE)sigma_theta.h

$(OBJ)bruiteur.o : $(SOURCE)bruiteur.c $(SOURCE)filtre.h

$(OBJ)gcuncert.o : $(SOURCE)gcuncert.c $(SOURCE)sigma_theta.h $(SOURCE)filtre.h $(SOURCE)3ch_sbr.h

$(OBJ)3_cornered_hat.o : $(SOURCE)3_cornered_hat.c $(SOURCE)sigma_theta.h $(SOURCE)filtre.h $(SOURCE)3ch_sbr.h

$(OBJ)ziggurat.o : $(SOURCE)ziggurat.c $(SOURCE)zigtables.h

install: $(TARGETS)
	$(INSTALL) -c -m 755 $(TARGETS) $(INSTALLDIR)

$(OBJ)%.o: $(SOURCE)%.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: all clean install
