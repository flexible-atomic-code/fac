#ifndef _INIT_H_
#define _INIT_H_ 1

#include <stdio.h>
#include <string.h>

#include "global.h"
#include "array.h"
#include "coulomb.h"
#include "config.h"
#include "cfp.h"
#include "angular.h"
#include "recouple.h"
#include "radial.h"
#include "nucleus.h"
#include "dbase.h"
#include "structure.h"
#include "transition.h"
#include "excitation.h"
#include "recombination.h"
#include "ionization.h"

int Info(void);
int InitFac(void);
int ReinitFac(int, int, int, int, int, int, int, int);

#endif
