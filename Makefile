#   SIGMA-THETA Project - Makefile for the binaries of the SIGMA-THETA set :
#	    - 1col2col
#           - X2Y
#           - DriRem
#           - SigmaTheta 
#           - ADev
#           - MDev
#	    - HDev
#	    - PDev
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
UB = /usr/bin/

all: $(BIN)1col2col $(BIN)X2Y $(BIN)DriRem $(BIN)SigmaTheta $(BIN)ADev $(BIN)MDev $(BIN)HDev $(BIN)PDev $(BIN)uncertainties $(BIN)RaylConfInt $(BIN)Asymptote $(BIN)Asym2Alpha $(BIN)AVarDOF $(BIN)ADUncert $(BIN)ADGraph $(BIN)PSDGraph $(BIN)YkGraph $(BIN)XtGraph $(BIN)bruiteur

$(BIN)1col2col : $(OBJ)1col2col.o $(OBJ)stio_sbr.o  
	gcc -o $(BIN)1col2col $(OBJ)1col2col.o $(OBJ)stio_sbr.o -lm

$(BIN)X2Y : $(OBJ)x2y.o $(OBJ)stio_sbr.o  
	gcc -o $(BIN)X2Y $(OBJ)x2y.o $(OBJ)stio_sbr.o -lm

$(BIN)DriRem : $(OBJ)drirem.o $(OBJ)tchebyfit.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)DriRem $(OBJ)drirem.o  $(OBJ)tchebyfit.o $(OBJ)stio_sbr.o -lm

$(BIN)SigmaTheta : $(OBJ)sigma_theta.o $(OBJ)dev_sbr.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o
	gcc -o $(BIN)SigmaTheta $(OBJ)sigma_theta.o $(OBJ)dev_sbr.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o -lgsl -lgslcblas -lm

$(BIN)ADev : $(OBJ)adev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	gcc -o $(BIN)ADev $(OBJ)adev.o  $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)MDev : $(OBJ)mdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	gcc -o $(BIN)MDev $(OBJ)mdev.o  $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)HDev : $(OBJ)hdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	gcc -o $(BIN)HDev $(OBJ)hdev.o  $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)PDev : $(OBJ)pdev.o $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o
	gcc -o $(BIN)PDev $(OBJ)pdev.o  $(OBJ)dev_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)uncertainties : $(OBJ)uncertainties.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)uncertainties $(OBJ)uncertainties.o $(OBJ)asymptote_sbr.o $(OBJ)avardof_sbr.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o -lgsl -lgslcblas -lm

$(BIN)RaylConfInt : $(OBJ)rayleigh.o $(OBJ)rayleigh_sbr.o 
	gcc -o $(BIN)RaylConfInt $(OBJ)rayleigh.o $(OBJ)rayleigh_sbr.o -lgsl -lgslcblas -lm

$(BIN)Asymptote : $(OBJ)asymptote.o $(OBJ)asymptote_sbr.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)Asymptote $(OBJ)asymptote.o $(OBJ)asymptote_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)Asym2Alpha : $(OBJ)asym2alpha.o $(OBJ)asymptote_sbr.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)Asym2Alpha $(OBJ)asym2alpha.o $(OBJ)asymptote_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)AVarDOF : $(OBJ)avardof.o $(OBJ)avardof_sbr.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)AVarDOF $(OBJ)avardof.o $(OBJ)avardof_sbr.o $(OBJ)stio_sbr.o -lm

$(BIN)ADUncert : $(OBJ)aduncert.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)ADUncert $(OBJ)aduncert.o $(OBJ)rayleigh_sbr.o $(OBJ)stio_sbr.o -lgsl -lgslcblas -lm

$(BIN)ADGraph : $(OBJ)adgraph.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)ADGraph $(OBJ)adgraph.o $(OBJ)stio_sbr.o -lm

$(BIN)PSDGraph : $(OBJ)psdgraph.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)PSDGraph $(OBJ)psdgraph.o $(OBJ)stio_sbr.o -lgsl -lgslcblas -lm

$(BIN)YkGraph : $(OBJ)ykgraph.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)YkGraph $(OBJ)ykgraph.o $(OBJ)stio_sbr.o -lm

$(BIN)XtGraph : $(OBJ)xtgraph.o $(OBJ)stio_sbr.o 
	gcc -o $(BIN)XtGraph $(OBJ)xtgraph.o $(OBJ)stio_sbr.o -lm

$(BIN)bruiteur : $(OBJ)bruiteur.o $(OBJ)filtre.o
	gcc -o $(BIN)bruiteur $(OBJ)bruiteur.o $(OBJ)filtre.o -lm

$(OBJ)1col2col.o : $(SOURCE)1col2col.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)1col2col.o $(SOURCE)1col2col.c

$(OBJ)x2y.o : $(SOURCE)x2y.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)x2y.o $(SOURCE)x2y.c

$(OBJ)drirem.o : $(SOURCE)drirem.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)drirem.o $(SOURCE)drirem.c

$(OBJ)sigma_theta.o : $(SOURCE)sigma_theta.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)sigma_theta.o $(SOURCE)sigma_theta.c

$(OBJ)adev.o : $(SOURCE)adev.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)adev.o $(SOURCE)adev.c

$(OBJ)mdev.o : $(SOURCE)mdev.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)mdev.o $(SOURCE)mdev.c

$(OBJ)hdev.o : $(SOURCE)hdev.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)hdev.o $(SOURCE)hdev.c

$(OBJ)pdev.o : $(SOURCE)pdev.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)pdev.o $(SOURCE)pdev.c

$(OBJ)uncertainties.o : $(SOURCE)uncertainties.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)uncertainties.o $(SOURCE)uncertainties.c

$(OBJ)tchebyfit.o : $(SOURCE)tchebyfit.c
	gcc -c -o $(OBJ)tchebyfit.o $(SOURCE)tchebyfit.c

$(OBJ)dev_sbr.o : $(SOURCE)dev_sbr.c
	gcc -c -o $(OBJ)dev_sbr.o $(SOURCE)dev_sbr.c

$(OBJ)asymptote_sbr.o : $(SOURCE)asymptote_sbr.c
	gcc -c -o $(OBJ)asymptote_sbr.o $(SOURCE)asymptote_sbr.c

$(OBJ)avardof_sbr.o : $(SOURCE)avardof_sbr.c
	gcc -c -o $(OBJ)avardof_sbr.o $(SOURCE)avardof_sbr.c

$(OBJ)rayleigh_sbr.o : $(SOURCE)rayleigh_sbr.c
	gcc -c -o $(OBJ)rayleigh_sbr.o $(SOURCE)rayleigh_sbr.c

$(OBJ)stio_sbr.o : $(SOURCE)stio_sbr.c
	gcc -c -o $(OBJ)stio_sbr.o $(SOURCE)stio_sbr.c

$(OBJ)rayleigh.o : $(SOURCE)rayleigh.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)rayleigh.o $(SOURCE)rayleigh.c

$(OBJ)asymptote.o : $(SOURCE)asymptote.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)asymptote.o $(SOURCE)asymptote.c

$(OBJ)asym2alpha.o : $(SOURCE)asym2alpha.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)asym2alpha.o $(SOURCE)asym2alpha.c

$(OBJ)avardof.o : $(SOURCE)avardof.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)avardof.o $(SOURCE)avardof.c

$(OBJ)aduncert.o : $(SOURCE)aduncert.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)aduncert.o $(SOURCE)aduncert.c

$(OBJ)adgraph.o : $(SOURCE)adgraph.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)adgraph.o $(SOURCE)adgraph.c

$(OBJ)psdgraph.o : $(SOURCE)psdgraph.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)psdgraph.o $(SOURCE)psdgraph.c

$(OBJ)ykgraph.o : $(SOURCE)ykgraph.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)ykgraph.o $(SOURCE)ykgraph.c

$(OBJ)xtgraph.o : $(SOURCE)xtgraph.c $(SOURCE)sigma_theta.h
	gcc -c -o $(OBJ)xtgraph.o $(SOURCE)xtgraph.c

$(OBJ)bruiteur.o : $(SOURCE)bruiteur.c $(SOURCE)filtre.h
	gcc -c -o $(OBJ)bruiteur.o $(SOURCE)bruiteur.c

$(OBJ)filtre.o : $(SOURCE)filtre.c
	gcc -c -o $(OBJ)filtre.o $(SOURCE)filtre.c

install: ~/.SigmaTheta.conf ~/.randinit2 $(UB)1col2col $(UB)X2Y $(UB)DriRem $(UB)SigmaTheta $(UB)ADev $(UB)MDev $(UB)HDev $(UB)PDev $(UB)uncertainties $(UB)RaylConfInt $(UB)Asymptote $(UB)Asym2Alpha $(UB)AVarDOF $(UB)ADUncert $(UB)ADGraph $(UB)PSDGraph $(UB)YkGraph $(UB)XtGraph $(UB)bruiteur

~/.SigmaTheta.conf: $(BIN).SigmaTheta.conf
	cp $(BIN).SigmaTheta.conf ~/.SigmaTheta.conf

~/.randinit2: $(BIN).randinit2
	cp $(BIN).randinit2 ~/.randinit2
	chown vernotte ~/.randinit2
	chgrp vernotte ~/.randinit2

$(UB)1col2col: $(BIN)1col2col
	cp $(BIN)1col2col $(UB)1col2col

$(UB)X2Y: $(BIN)X2Y
	cp $(BIN)X2Y $(UB)X2Y

$(UB)DriRem: $(BIN)DriRem 
	cp $(BIN)DriRem $(UB)DriRem

$(UB)SigmaTheta: $(BIN)SigmaTheta
	cp $(BIN)SigmaTheta $(UB)SigmaTheta

$(UB)ADev: $(BIN)ADev
	cp $(BIN)ADev $(UB)ADev

$(UB)MDev: $(BIN)MDev
	cp $(BIN)MDev $(UB)MDev

$(UB)HDev: $(BIN)HDev
	cp $(BIN)HDev $(UB)HDev

$(UB)PDev: $(BIN)PDev
	cp $(BIN)PDev $(UB)PDev

$(UB)uncertainties: $(BIN)uncertainties
	cp $(BIN)uncertainties $(UB)uncertainties

$(UB)RaylConfInt: $(BIN)RaylConfInt
	cp $(BIN)RaylConfInt $(UB)RaylConfInt

$(UB)Asymptote: $(BIN)Asymptote
	cp $(BIN)Asymptote $(UB)Asymptote

$(UB)Asym2Alpha: $(BIN)Asym2Alpha
	cp $(BIN)Asym2Alpha $(UB)Asym2Alpha

$(UB)AVarDOF: $(BIN)AVarDOF
	cp $(BIN)AVarDOF $(UB)AVarDOF

$(UB)ADUncert: $(BIN)ADUncert
	cp $(BIN)ADUncert $(UB)ADUncert

$(UB)ADGraph: $(BIN)ADGraph
	cp $(BIN)ADGraph $(UB)ADGraph

$(UB)PSDGraph: $(BIN)PSDGraph
	cp $(BIN)PSDGraph $(UB)PSDGraph

$(UB)YkGraph: $(BIN)YkGraph
	cp $(BIN)YkGraph $(UB)YkGraph

$(UB)XtGraph: $(BIN)XtGraph
	cp $(BIN)XtGraph $(UB)XtGraph

$(UB)bruiteur: $(BIN)bruiteur
	cp $(BIN)bruiteur $(UB)bruiteur
