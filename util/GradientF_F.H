#ifndef _GRADIENTF_F_H_
#define _GRADIENTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDNEWMACGRAD 
#define GUARDNEWMACGRAD 
// Prototype for Fortran procedure NEWMACGRAD ...
//
void FORTRAN_NAME( NEWMACGRAD ,newmacgrad )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir)
      ,CHFp_INT(edgeDir) );

#define FORT_NEWMACGRAD FORTRAN_NAME( inlineNEWMACGRAD, inlineNEWMACGRAD)
#define FORTNT_NEWMACGRAD FORTRAN_NAME( NEWMACGRAD, newmacgrad)

inline void FORTRAN_NAME(inlineNEWMACGRAD, inlineNEWMACGRAD)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir)
      ,CHFp_INT(edgeDir) )
{
 CH_TIMELEAF("FORT_NEWMACGRAD");
 FORTRAN_NAME( NEWMACGRAD ,newmacgrad )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(edgeGrid)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(dir)
      ,CHFt_INT(edgeDir) );
}
#endif  // GUARDNEWMACGRAD 

#ifndef GUARDCRSEONESIDEGRAD 
#define GUARDCRSEONESIDEGRAD 
// Prototype for Fortran procedure CRSEONESIDEGRAD ...
//
void FORTRAN_NAME( CRSEONESIDEGRAD ,crseonesidegrad )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_FIA1(mask)
      ,CHFp_BOX(loEdgeBox)
      ,CHFp_BOX(hiEdgeBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(doLo)
      ,CHFp_INT(doHi) );

#define FORT_CRSEONESIDEGRAD FORTRAN_NAME( inlineCRSEONESIDEGRAD, inlineCRSEONESIDEGRAD)
#define FORTNT_CRSEONESIDEGRAD FORTRAN_NAME( CRSEONESIDEGRAD, crseonesidegrad)

inline void FORTRAN_NAME(inlineCRSEONESIDEGRAD, inlineCRSEONESIDEGRAD)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_FIA1(mask)
      ,CHFp_BOX(loEdgeBox)
      ,CHFp_BOX(hiEdgeBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(doLo)
      ,CHFp_INT(doHi) )
{
 CH_TIMELEAF("FORT_CRSEONESIDEGRAD");
 FORTRAN_NAME( CRSEONESIDEGRAD ,crseonesidegrad )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_FIA1(mask)
      ,CHFt_BOX(loEdgeBox)
      ,CHFt_BOX(hiEdgeBox)
      ,CHFt_INT(dir)
      ,CHFt_INT(doLo)
      ,CHFt_INT(doHi) );
}
#endif  // GUARDCRSEONESIDEGRAD 

#ifndef GUARDGRADCC 
#define GUARDGRADCC 
// Prototype for Fortran procedure GRADCC ...
//
void FORTRAN_NAME( GRADCC ,gradcc )(
      CHFp_FRA1(gradPhi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) );

#define FORT_GRADCC FORTRAN_NAME( inlineGRADCC, inlineGRADCC)
#define FORTNT_GRADCC FORTRAN_NAME( GRADCC, gradcc)

inline void FORTRAN_NAME(inlineGRADCC, inlineGRADCC)(
      CHFp_FRA1(gradPhi)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_BOX(gridBox)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_INT(dir) )
{
 CH_TIMELEAF("FORT_GRADCC");
 FORTRAN_NAME( GRADCC ,gradcc )(
      CHFt_FRA1(gradPhi)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_BOX(gridBox)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_INT(dir) );
}
#endif  // GUARDGRADCC 

}

#endif
