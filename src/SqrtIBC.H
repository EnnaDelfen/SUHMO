#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  SqrtIBC.H
// ============
//
// PhysIBC-derived class for basic (simple) IBC
//

#ifndef _SQRTIBC_H_
#define _SQRTIBC_H_

#include "HydroIBC.H"
#include "suhmo_params.H"

#include "NamespaceHeader.H"


/// Physical/domain initial and boundary conditions
/**

 */
class SqrtIBC : public HydroIBC
{
public:
    /// Constructor
    /**
     */
    SqrtIBC();

    /// Destructor
    /**
     */
    virtual ~SqrtIBC();

    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object.
     */
    virtual void define(const ProblemDomain& a_domain, const Real& a_dx);

    /// Factory method - this object is its own factory
    /**
       Returns a pointer to a new SqrtIBC object with
       m_isDefined = false (i.e., its define() must be called before it
       is used).
     */
    virtual HydroIBC* new_hydroIBC();

    /// Set up initial conditions
    /**
     */
    virtual void initialize(LevelData<FArrayBox>& a_U);

    virtual void initializeData(RealVect& a_dx,
                                suhmo_params params,
                                LevelData<FArrayBox>& a_head,
                                LevelData<FArrayBox>& a_gapHeight,
                                LevelData<FArrayBox>& a_Pw,
                                LevelData<FArrayBox>& a_qw,
                                LevelData<FArrayBox>& a_Re,
                                LevelData<FArrayBox>& a_meltRate,
                                LevelData<FArrayBox>& a_zbed,
                                LevelData<FArrayBox>& a_Pi,
                                LevelData<FArrayBox>& a_iceHeight,
                                LevelData<FArrayBox>& a_bumpHeight,
                                LevelData<FArrayBox>& a_bumpSpacing);

    virtual void initializePi(RealVect& a_dx,
                              suhmo_params Params,     
                              LevelData<FArrayBox>& a_Pi);

    virtual void initializePi(RealVect& a_dx,
                              suhmo_params Params,     
                              LevelData<FArrayBox>& a_head,
                              LevelData<FArrayBox>& a_gapHeight,
                              LevelData<FArrayBox>& a_Pw,
                              LevelData<FArrayBox>& a_zbed,
                              LevelData<FArrayBox>& a_Pi,
                              LevelData<FArrayBox>& a_iceHeight,
                              LevelData<FArrayBox>& a_bumpHeight,
                              LevelData<FArrayBox>& a_bumpSpacing);

    virtual void initializeBed(RealVect& a_dx,
                               suhmo_params Params,     
                               LevelData<FArrayBox>& a_zbed,
                               LevelData<FArrayBox>& a_bumpHeight,
                               LevelData<FArrayBox>& a_bumpSpacing);

protected:

private:
    // Disallowed for all the usual reasons
    void operator=(const SqrtIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    SqrtIBC(const SqrtIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
