#ifndef _DIVERGENCEF_F_H_
#define _DIVERGENCEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDDIVERGENCE 
#define GUARDDIVERGENCE 
// Prototype for Fortran procedure DIVERGENCE ...
//
void FORTRAN_NAME( DIVERGENCE ,divergence )(
      CHFp_CONST_FRA(uEdge)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(idir) );

#define FORT_DIVERGENCE FORTRAN_NAME( inlineDIVERGENCE, inlineDIVERGENCE)
#define FORTNT_DIVERGENCE FORTRAN_NAME( DIVERGENCE, divergence)

inline void FORTRAN_NAME(inlineDIVERGENCE, inlineDIVERGENCE)(
      CHFp_CONST_FRA(uEdge)
      ,CHFp_FRA(div)
      ,CHFp_BOX(gridInt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(idir) )
{
 CH_TIMELEAF("FORT_DIVERGENCE");
 FORTRAN_NAME( DIVERGENCE ,divergence )(
      CHFt_CONST_FRA(uEdge)
      ,CHFt_FRA(div)
      ,CHFt_BOX(gridInt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(idir) );
}
#endif  // GUARDDIVERGENCE 

}

#endif
