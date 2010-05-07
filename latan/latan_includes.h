#ifndef LATAN_INCLUDES_H_
#define LATAN_INCLUDES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* replacement functions */
#ifndef HAVE_ACOSH
#ifdef acosh
#undef acosh
#endif
#define acosh gsl_acosh
#endif

/* system includes */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <gsl/gsl_math.h>

#endif