#ifndef _AMRUTILF_F_H_
#define _AMRUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUNDIVIDEDGRAD 
#define GUARDUNDIVIDEDGRAD 
// Prototype for Fortran procedure UNDIVIDEDGRAD ...
//
void FORTRAN_NAME( UNDIVIDEDGRAD ,undividedgrad )(
      CHFp_FRA1(dq)
      ,CHFp_CONST_FRA(q)
      ,CHFp_BOX(dInterior)
      ,CHFp_BOX(loedge)
      ,CHFp_BOX(hiedge)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi) );

#define FORT_UNDIVIDEDGRAD FORTRAN_NAME( inlineUNDIVIDEDGRAD, inlineUNDIVIDEDGRAD)
#define FORTNT_UNDIVIDEDGRAD FORTRAN_NAME( UNDIVIDEDGRAD, undividedgrad)

inline void FORTRAN_NAME(inlineUNDIVIDEDGRAD, inlineUNDIVIDEDGRAD)(
      CHFp_FRA1(dq)
      ,CHFp_CONST_FRA(q)
      ,CHFp_BOX(dInterior)
      ,CHFp_BOX(loedge)
      ,CHFp_BOX(hiedge)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_CONST_INT(hashi) )
{
 CH_TIMELEAF("FORT_UNDIVIDEDGRAD");
 FORTRAN_NAME( UNDIVIDEDGRAD ,undividedgrad )(
      CHFt_FRA1(dq)
      ,CHFt_CONST_FRA(q)
      ,CHFt_BOX(dInterior)
      ,CHFt_BOX(loedge)
      ,CHFt_BOX(hiedge)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_CONST_INT(hashi) );
}
#endif  // GUARDUNDIVIDEDGRAD 

}

#endif
