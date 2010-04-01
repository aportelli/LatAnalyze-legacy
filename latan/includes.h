#ifndef LATAN_INCLUDES_H_
#define LATAN_INCLUDES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* replacement functions */
#ifndef HAVE_ACOSH
#define acosh gsl_acosh
#endif

/* system includes */
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <gsl/gsl_math.h>

#endif