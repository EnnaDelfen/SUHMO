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
//  DerivedIBC.H
// ============
//
// PhysIBC-derived class for basic (simple) IBC
//

#ifndef _DERIVEDIBC_H_
#define _DERIVEDIBC_H_

#include "BasicIBC.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**

 */
class DerivedIBC : public BasicIBC
{
public:
    /// Constructor
    /**
     */
    DerivedIBC();

    /// Destructor
    /**
     */
    virtual ~DerivedIBC();

    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object.
     */
    virtual void define(const ProblemDomain& a_domain, const Real& a_dx);

    /// Factory method - this object is its own factory
    /**
       Return a pointer to a new BasicIBC object with
       m_isDefined = false (i.e., its define() must be called before it
       is used).
     */
    virtual BasicIBC* new_thicknessIBC();

#if 0
  /// set boundary phi
  void setBoundaryPhi(Real a_boundaryThickness)
  { m_boundaryPhi = a_boundaryPhi; }
#endif

    /// Set up initial conditions
    /**
     */
    virtual void initialize(LevelData<FArrayBox>& a_U);

    virtual void initializeData(RealVect& a_dx,
                             LevelData<FArrayBox>& a_head,
                             LevelData<FArrayBox>& a_gapHeight,
                             LevelData<FArrayBox>& a_Pw,
                             LevelData<FArrayBox>& a_qw,
                             LevelData<FArrayBox>& a_Re,
                             LevelData<FArrayBox>& a_meltRate,
                             LevelData<FArrayBox>& a_zbed,
                             LevelData<FArrayBox>& a_Pi) override;

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

    /// return boundary condition for Ice velocity solve
    /**
     */
    virtual RefCountedPtr<CompGridVTOBC> velocitySolveBC();

protected:
    // phi at domain boundaries
    Real m_boundaryPhi;

    RealVect m_domainSize;

    // have bc's been set up?
    bool m_isBCsetUp;

    /// set up bc's based on parmparse inputs
    void setupBCs();

    RefCountedPtr<CompGridVTOBC> m_velBCs;

    // Parameters
    Real m_slope;
    Real m_H;
    Real m_gapInit;
   
    // Cst of model
    Real m_rho_ice   = 910.; // kg.m-3
    Real m_rho_water = 1000.; // kg.m-3
    Real m_gravity   = 9.8;  // m.s-1
    Real m_G         = 0.05;   // W.m-2
    Real m_L         = 3.34e-5; // J.KG-1

private:
    // Disallowed for all the usual reasons
    void operator=(const DerivedIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    DerivedIBC(const DerivedIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
