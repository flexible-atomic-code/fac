/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _INIT_H_
#define _INIT_H_ 1

#include <stdio.h>
#include <string.h>

#include "array.h"
#include "global.h"
#include "coulomb.h"
#include "config.h"
#include "cfp.h"
#include "angular.h"
#include "recouple.h"
#include "radial.h"
#include "nucleus.h"
#include "dbase.h"
#include "structure.h"
#include "mbpt.h"
#include "transition.h"
#include "excitation.h"
#include "recombination.h"
#include "ionization.h"
#include "rmatrix.h"

int Info(void);
int InitFac(void);
int ReinitFac(int, int, int, int, int, int, int, int);
void SetOption(char *s, char *sp, int ip, double dp);
#endif
