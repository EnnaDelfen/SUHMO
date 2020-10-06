#ifndef _EXTRAPBCF_F_H_
#define _EXTRAPBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSIMPLEEXTRAPBC 
#define GUARDSIMPLEEXTRAPBC 
// Prototype for Fortran procedure SIMPLEEXTRAPBC ...
//
void FORTRAN_NAME( SIMPLEEXTRAPBC ,simpleextrapbc )(
      CHFp_FRA(phi)
      ,CHFp_BOX(bcbox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) );

#define FORT_SIMPLEEXTRAPBC FORTRAN_NAME( inlineSIMPLEEXTRAPBC, inlineSIMPLEEXTRAPBC)
#define FORTNT_SIMPLEEXTRAPBC FORTRAN_NAME( SIMPLEEXTRAPBC, simpleextrapbc)

inline void FORTRAN_NAME(inlineSIMPLEEXTRAPBC, inlineSIMPLEEXTRAPBC)(
      CHFp_FRA(phi)
      ,CHFp_BOX(bcbox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) )
{
 CH_TIMELEAF("FORT_SIMPLEEXTRAPBC");
 FORTRAN_NAME( SIMPLEEXTRAPBC ,simpleextrapbc )(
      CHFt_FRA(phi)
      ,CHFt_BOX(bcbox)
      ,CHFt_INT(dir)
      ,CHFt_INT(hiLo) );
}
#endif  // GUARDSIMPLEEXTRAPBC 

#ifndef GUARDOLDSIMPLEREFLECTBC 
#define GUARDOLDSIMPLEREFLECTBC 
// Prototype for Fortran procedure OLDSIMPLEREFLECTBC ...
//
void FORTRAN_NAME( OLDSIMPLEREFLECTBC ,oldsimplereflectbc )(
      CHFp_FRA(phi)
      ,CHFp_BOX(bcbox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) );

#define FORT_OLDSIMPLEREFLECTBC FORTRAN_NAME( inlineOLDSIMPLEREFLECTBC, inlineOLDSIMPLEREFLECTBC)
#define FORTNT_OLDSIMPLEREFLECTBC FORTRAN_NAME( OLDSIMPLEREFLECTBC, oldsimplereflectbc)

inline void FORTRAN_NAME(inlineOLDSIMPLEREFLECTBC, inlineOLDSIMPLEREFLECTBC)(
      CHFp_FRA(phi)
      ,CHFp_BOX(bcbox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) )
{
 CH_TIMELEAF("FORT_OLDSIMPLEREFLECTBC");
 FORTRAN_NAME( OLDSIMPLEREFLECTBC ,oldsimplereflectbc )(
      CHFt_FRA(phi)
      ,CHFt_BOX(bcbox)
      ,CHFt_INT(dir)
      ,CHFt_INT(hiLo) );
}
#endif  // GUARDOLDSIMPLEREFLECTBC 

#ifndef GUARDSIMPLEREFLECTBC 
#define GUARDSIMPLEREFLECTBC 
// Prototype for Fortran procedure SIMPLEREFLECTBC ...
//
void FORTRAN_NAME( SIMPLEREFLECTBC ,simplereflectbc )(
      CHFp_FRA(phi)
      ,CHFp_BOX(bcbox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) );

#define FORT_SIMPLEREFLECTBC FORTRAN_NAME( inlineSIMPLEREFLECTBC, inlineSIMPLEREFLECTBC)
#define FORTNT_SIMPLEREFLECTBC FORTRAN_NAME( SIMPLEREFLECTBC, simplereflectbc)

inline void FORTRAN_NAME(inlineSIMPLEREFLECTBC, inlineSIMPLEREFLECTBC)(
      CHFp_FRA(phi)
      ,CHFp_BOX(bcbox)
      ,CHFp_INT(dir)
      ,CHFp_INT(hiLo) )
{
 CH_TIMELEAF("FORT_SIMPLEREFLECTBC");
 FORTRAN_NAME( SIMPLEREFLECTBC ,simplereflectbc )(
      CHFt_FRA(phi)
      ,CHFt_BOX(bcbox)
      ,CHFt_INT(dir)
      ,CHFt_INT(hiLo) );
}
#endif  // GUARDSIMPLEREFLECTBC 

#ifndef GUARDEXTRAPCORNER2D 
#define GUARDEXTRAPCORNER2D 
// Prototype for Fortran procedure EXTRAPCORNER2D ...
//
void FORTRAN_NAME( EXTRAPCORNER2D ,extrapcorner2d )(
      CHFp_FRA(phi)
      ,CHFp_BOX(box) );

#define FORT_EXTRAPCORNER2D FORTRAN_NAME( inlineEXTRAPCORNER2D, inlineEXTRAPCORNER2D)
#define FORTNT_EXTRAPCORNER2D FORTRAN_NAME( EXTRAPCORNER2D, extrapcorner2d)

inline void FORTRAN_NAME(inlineEXTRAPCORNER2D, inlineEXTRAPCORNER2D)(
      CHFp_FRA(phi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_EXTRAPCORNER2D");
 FORTRAN_NAME( EXTRAPCORNER2D ,extrapcorner2d )(
      CHFt_FRA(phi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDEXTRAPCORNER2D 

}

#endif
