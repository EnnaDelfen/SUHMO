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
//  MountainIBC.H
// ============
//
// PhysIBC-derived class for basic (simple) IBC
//

#ifndef _MOUNTAINIBC_H_
#define _MOUNTAINIBC_H_

#include "HydroIBC.H"
#include "suhmo_params.H"

#include "NamespaceHeader.H"

//class VCPoissonParameters
//{
//public:
//
//  void getPoissonParameters (VCPoissonParameters&  a_params);
//
//  IntVect      nCells;
//  int          maxGridSize;
//  int          blockFactor;
//  int          bufferSize;
//  Real         fillRatio;
//  Real         refineThresh;
//  int          domBcType;
//  int          maxLevel;
//  int          numLevels;
//  Vector<int>  refRatio;
//  ProblemDomain coarsestDomain;
//  Real         coarsestDx;
//  RealVect     domainLength;
//  RealVect     probLo;
//  RealVect     probHi;
//  Real         alpha;
//  Real         beta;
//};

/// Physical/domain initial and boundary conditions
/**

 */
class MountainIBC : public HydroIBC
{
public:
    /// Constructor
    /**
     */
    MountainIBC();

    /// Destructor
    /**
     */
    virtual ~MountainIBC();

    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object.
     */
    virtual void define(const ProblemDomain& a_domain, const Real& a_dx);

    /// Factory method - this object is its own factory
    /**
       Returns a pointer to a new MountainIBC object with
       m_isDefined = false (i.e., its define() must be called before it
       is used).
     */
    virtual HydroIBC* new_hydroIBC();

    /* set parameters */ 
    void setParameters(const Real& a_gamma);

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
                                LevelData<FArrayBox>& a_bumpSpacing,
                                LevelData<FArrayBox>& a_levelmagVel);

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

    virtual void resetCovered(RealVect& a_dx,
                              suhmo_params Params,
                              LevelData<FArrayBox>& a_head,
                              LevelData<FArrayBox>& a_iceHeight,
                              LevelData<FArrayBox>& a_Pi);


protected:
    Real m_gamma;

private:
    // Disallowed for all the usual reasons
    void operator=(const MountainIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    MountainIBC(const MountainIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
