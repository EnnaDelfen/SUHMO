#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _ICECONSTANTS_H_
#define _ICECONSTANTS_H_
#include "BaseNamespaceHeader.H"


#define GROUNDEDMASKVAL (1)
#define groundedmaskval (1)
#define FLOATINGMASKVAL (2)
#define floatingmaskval (2)
#define OPENSEAMASKVAL (4)
#define openseamaskval (4)
#define OPENLANDMASKVAL (8)
#define openlandmaskval (8)
#define OCEANMASKVAL OPENSEAMASKVAL
#define oceanmaskval openseamaskval

#if defined(CH_USE_FLOAT) || defined(CH_LANG_CC)

#define seconds_per_tropical_year     (31556926)
#define SECONDS_PER_TROPICAL_YEAR     (31556926)

#define HUGE_NORM (1.0e10)
#define TINY_NORM (1.0e-10)

#elif defined(CH_USE_DOUBLE) && !defined(CH_LANG_CC)

#define seconds_per_tropical_year     (31556926d0)
#define SECONDS_PER_TROPICAL_YEAR     (31556926d0)

#define HUGE_NORM (1.0d10)
#define TINY_NORM (1.0d-10)

#else

#error One of CH_USE_{FLOAT,DOUBLE} must be defined

#endif

#include "BaseNamespaceFooter.H"
#endif /*CH_ICECONSTANTS_H*/
