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
INSTALLDIR = /raid0/bin/
BINPREFIX = st_
INSTALL = install
CC = gcc
GIT_VERSION := "$(shell git describe --abbrev=4 --dirty --always --tags)"
GIT_DATE := "$(shell git log | head -3 | grep Date | sed 's/Date:   //')"
CFLAGS = -g -O3 -DST_VERSION=\"$(GIT_VERSION)\" -DST_DATE=\"$(GIT_DATE)\"
# link statically against libfftw3
#FFTW3 = -static -lfftw3
# link dynamically against libfftw3
FFTW3 = -lfftw3

TARGETS = $(BIN)$(BINPREFIX)1col2col $(BIN)$(BINPREFIX)X2Y $(BIN)$(BINPREFIX)DriRem $(BIN)$(BINPREFIX)SigmaTheta $(BIN)$(BINPREFIX)ADev $(BIN)$(BINPREFIX)find_tau_xdev $(BIN)$(BINPREFIX)GCoDev $(BIN)$(BINPREFIX)MDev $(BIN)$(BINPREFIX)HDev $(BIN)$(BINPREFIX)PDev $(BIN)$(BINPREFIX)Aver $(BIN)$(BINPREFIX)uncertainties $(BIN)$(BINPREFIX)RaylConfInt $(BIN)$(BINPREFIX)Asymptote $(BIN)$(BINPREFIX)Asym2Alpha $(BIN)$(BINPREFIX)AVarDOF $(BIN)$(BINPREFIX)ADUncert $(BIN)$(BINPREFIX)ADGraph $(BIN)$(BINPREFIX)PSDGraph $(BIN)$(BINPREFIX)YkGraph $(BIN)$(BINPREFIX)XtGraph $(BIN)$(BINPREFIX)DevGraph $(BIN)$(BINPREFIX)3CHGraph $(BIN)$(BINPREFIX)bruiteur $(BIN)$(BINPREFIX)GCUncert $(BIN)$(BINPREFIX)3CorneredHat

all: $(TARGETS)

clean:
	rm -f $(TARGETS)
	rm -f $(OBJ)/*.o

$(BIN)$(BINPREFIX)1col2col : $(OBJ)$(BINPREFIX)1col2col.o $(OBJ)$(BINPREFIX)stio_sbr.o  
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)X2Y : $(OBJ)$(BINPREFIX)x2y.o $(OBJ)$(BINPREFIX)stio_sbr.o  
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)DriRem : $(OBJ)$(BINPREFIX)drirem.o $(OBJ)$(BINPREFIX)tchebyfit.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)SigmaTheta : $(OBJ)$(BINPREFIX)sigma_theta.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)asymptote_sbr.o $(OBJ)$(BINPREFIX)avardof_sbr.o $(OBJ)$(BINPREFIX)rayleigh_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)ADev : $(OBJ)$(BINPREFIX)adev.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)find_tau_xdev : $(OBJ)$(BINPREFIX)find_tau_xdev.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)GCoDev : $(OBJ)$(BINPREFIX)gcodev.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)MDev : $(OBJ)$(BINPREFIX)mdev.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)HDev : $(OBJ)$(BINPREFIX)hdev.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)PDev : $(OBJ)$(BINPREFIX)pdev.o $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)Aver : $(OBJ)$(BINPREFIX)aver.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -Wl,-Bdynamic -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)uncertainties : $(OBJ)$(BINPREFIX)uncertainties.o $(OBJ)$(BINPREFIX)asymptote_sbr.o $(OBJ)$(BINPREFIX)avardof_sbr.o $(OBJ)$(BINPREFIX)rayleigh_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)RaylConfInt : $(OBJ)$(BINPREFIX)rayleigh.o $(OBJ)$(BINPREFIX)rayleigh_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)Asymptote : $(OBJ)$(BINPREFIX)asymptote.o $(OBJ)$(BINPREFIX)asymptote_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)Asym2Alpha : $(OBJ)$(BINPREFIX)asym2alpha.o $(OBJ)$(BINPREFIX)asymptote_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)AVarDOF : $(OBJ)$(BINPREFIX)avardof.o $(OBJ)$(BINPREFIX)avardof_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)ADUncert : $(OBJ)$(BINPREFIX)aduncert.o $(OBJ)$(BINPREFIX)rayleigh_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)ADGraph : $(OBJ)$(BINPREFIX)adgraph.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)PSDGraph : $(OBJ)$(BINPREFIX)psdgraph.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)YkGraph : $(OBJ)$(BINPREFIX)ykgraph.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)XtGraph : $(OBJ)$(BINPREFIX)xtgraph.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)DevGraph : $(OBJ)$(BINPREFIX)devgraph.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)3CHGraph : $(OBJ)$(BINPREFIX)3chgraph.o $(OBJ)$(BINPREFIX)stio_sbr.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BIN)$(BINPREFIX)bruiteur : $(OBJ)$(BINPREFIX)bruiteur.o $(OBJ)$(BINPREFIX)filtre.o $(OBJ)$(BINPREFIX)splitmix64.o $(OBJ)$(BINPREFIX)xorshift1024star.o $(OBJ)$(BINPREFIX)ziggurat.o
	$(CC) $(CFLAGS) -o $@ $^ $(FFTW3) -Wl,-Bdynamic -lm

$(BIN)$(BINPREFIX)GCUncert : $(OBJ)$(BINPREFIX)gcuncert.o $(OBJ)$(BINPREFIX)stio_sbr.o $(OBJ)$(BINPREFIX)filtre.o $(OBJ)$(BINPREFIX)3ch_sbr.o $(OBJ)$(BINPREFIX)splitmix64.o $(OBJ)$(BINPREFIX)xorshift1024star.o $(OBJ)$(BINPREFIX)ziggurat.o
	$(CC) $(CFLAGS) -o $@ $^ $(FFTW3) -Wl,-Bdynamic -lgsl -lgslcblas -lm

$(BIN)$(BINPREFIX)3CorneredHat : $(OBJ)$(BINPREFIX)3_cornered_hat.o $(OBJ)$(BINPREFIX)asymptote_sbr.o $(OBJ)$(BINPREFIX)avardof_sbr.o $(OBJ)$(BINPREFIX)stio_sbr.o  $(OBJ)$(BINPREFIX)dev_sbr.o $(OBJ)$(BINPREFIX)tchebyfit.o $(OBJ)$(BINPREFIX)filtre.o $(OBJ)$(BINPREFIX)3ch_sbr.o $(OBJ)$(BINPREFIX)splitmix64.o $(OBJ)$(BINPREFIX)xorshift1024star.o $(OBJ)$(BINPREFIX)ziggurat.o
	$(CC) $(CFLAGS) -o $@ $^ $(FFTW3) -Wl,-Bdynamic -lgsl -lgslcblas -lm

$(OBJ)$(BINPREFIX)1col2col.o : $(SOURCE)1col2col.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)x2y.o : $(SOURCE)x2y.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)drirem.o : $(SOURCE)drirem.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)sigma_theta.o : $(SOURCE)sigma_theta.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)adev.o : $(SOURCE)adev.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)gcodev.o : $(SOURCE)gcodev.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)mdev.o : $(SOURCE)mdev.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)hdev.o : $(SOURCE)hdev.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)pdev.o : $(SOURCE)pdev.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)aver.o : $(SOURCE)aver.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)uncertainties.o : $(SOURCE)uncertainties.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)rayleigh.o : $(SOURCE)rayleigh.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)asymptote.o : $(SOURCE)asymptote.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)asym2alpha.o : $(SOURCE)asym2alpha.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)avardof.o : $(SOURCE)avardof.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)aduncert.o : $(SOURCE)aduncert.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)adgraph.o : $(SOURCE)adgraph.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)psdgraph.o : $(SOURCE)psdgraph.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)ykgraph.o : $(SOURCE)ykgraph.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)xtgraph.o : $(SOURCE)xtgraph.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)devgraph.o : $(SOURCE)devgraph.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)3chgraph.o : $(SOURCE)3chgraph.c $(SOURCE)sigma_theta.h

$(OBJ)$(BINPREFIX)bruiteur.o : $(SOURCE)bruiteur.c $(SOURCE)filtre.h

$(OBJ)$(BINPREFIX)gcuncert.o : $(SOURCE)gcuncert.c $(SOURCE)sigma_theta.h $(SOURCE)filtre.h $(SOURCE)3ch_sbr.h

$(OBJ)$(BINPREFIX)3_cornered_hat.o : $(SOURCE)3_cornered_hat.c $(SOURCE)sigma_theta.h $(SOURCE)filtre.h $(SOURCE)3ch_sbr.h

$(OBJ)$(BINPREFIX)ziggurat.o : $(SOURCE)ziggurat.c $(SOURCE)zigtables.h

install: $(TARGETS)
	$(INSTALL) -c -m 755 $(TARGETS) $(INSTALLDIR)

$(OBJ)$(BINPREFIX)%.o: $(SOURCE)%.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: all clean install
