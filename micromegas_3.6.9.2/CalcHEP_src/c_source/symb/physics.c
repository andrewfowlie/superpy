/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su
*/

#include "syst2.h"
#include "process.h"
#include "physics.h"


vcsect     vcs;


int    subproc_f, subproc_sq;

char  modelmenu[STRSIZ];
int   maxmodel;

int nsub=1, ndiagr=0, n_model=1;
