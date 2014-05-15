/* Fichier include contenant les declarations des fonctions des recipies
 * de la librairie librecipies.a
 * C'est le seul fichier a inclure dans un programme utilisant
 * cette librairie
 *
 * mars 98 Ch Muller Guillaume*/

#ifndef __RECIPIES_H__
#define __RECIPIES_H__


#include <stdio.h>


typedef struct FCOMPLEX {float r,i;} fcomplex;
typedef struct IMMENSE {unsigned long l,r;} immense;
typedef struct GREAT {unsigned short l,c,r;} great;

#include "recipies_ext.h"      /* Fonctions */
#include "nrutil_ext.h"	/* gestion des erreurs et de la memoire */	
#include "complex.h"		

#endif
