#ifndef _DBASE_H_
#define _DBASE_H_

#include <string.h>
#include <stdio.h>

#include "global.h"

#define FILE_NAME_LEN 256 

#define DB_NONE  0x00
#define DB_READ  0x01
#define DB_WRITE 0x02
#define DB_ALL   0xFF

#define DBReadable(mode)  (((mode) & DB_READ)  == DB_READ )
#define DBWritable(mode)  (((mode) & DB_WRITE) == DB_WRITE)

int SetDBase(int i, char *s, char mode);
int ClearDBase(int i);

#if FAC_DEBUG
FILE *debug_log;
#endif

int InitDBase();

#endif
