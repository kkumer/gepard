/*
 *  Flush output from fotran, so that one can easier
 *  observe output in terminal, while saving it to file
 *   */

#include <stdio.h>

void flushout_() 
{
          fflush(NULL);
} 
