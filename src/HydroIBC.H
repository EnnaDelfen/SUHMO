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
//  HydroIBC.H
// ============
//
// PhysIBC-derived class for basic (simple) IBC
//

#ifndef _HYDROIBC_H_
#define _HYDROIBC_H_

#include "PhysIBC.H"
#include "BCFunc.H"
#include "suhmo_params.H"

#include "NamespaceHeader.H"


/// Physical/domain initial and boundary conditions
/**

 */
class HydroIBC : public PhysIBC
{
public:
    /// Constructor
    /**
     */
    HydroIBC();

    /// Destructor
    /**
     */
    virtual ~HydroIBC();

    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object.
     */
    virtual void define(const ProblemDomain& a_domain, const Real& a_dx);

    /// Factory method - this object is its own factory
    /**
       Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
       its define() must be called before it is used).
     */
    virtual PhysIBC* new_physIBC()
    {
        return static_cast<PhysIBC*>(new_hydroIBC());
    }

    /// Factory method - this object is its own factory
    /**
       Returns a pointer to a new HydroIBC object with
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

    virtual void resetCovered(suhmo_params Params,
                              LevelData<FArrayBox>& a_head,
                              LevelData<FArrayBox>& a_Pi);

    virtual void setup_iceMask(LevelData<FArrayBox>& a_Pi,
                               LevelData<FArrayBox>& a_iceMask);

    virtual void setup_iceMask_EC(LevelData<FArrayBox>& a_iceMask,
                                  LevelData<FluxBox>&   a_iceMask_ec);

    /// Set boundary fluxes
    /**
     */
    virtual void primBC(FArrayBox& a_WGdnv,
                        const FArrayBox& a_Wextrap,
                        const FArrayBox& a_W,
                        const int& a_dir,
                        const Side::LoHiSide& a_side,
                        const Real& a_time);

    /// Set boundary slopes
    /**
       The boundary slopes in a_dW are already set to one sided difference
       approximations.  If this function doesn't change them they will be
       used for the slopes at the boundaries.
     */
    virtual void setBdrySlopes(FArrayBox& a_dW, const FArrayBox& a_W, const int& a_dir, const Real& a_time);

    /// Adjust boundary fluxes to account for artificial viscosity
    /**
     */
    virtual void
    artViscBC(FArrayBox& a_F, const FArrayBox& a_U, const FArrayBox& a_divVel, const int& a_dir, const Real& a_time);


protected:

private:
    // Disallowed for all the usual reasons
    void operator=(const HydroIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    HydroIBC(const HydroIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
