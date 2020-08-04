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
//  BasicIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _BASICIBC_H_
#define _BASICIBC_H_

#include "PhysIBC.H"
#include "PetscCompGridVTO.H"
#include "BCFunc.H"

#include "NamespaceHeader.H"

/// this BCFunction simply wraps a BCFunc
class IceBCFuncWrapper : public CompGridVTOBC
{
public:
    IceBCFuncWrapper() : m_funcptr(NULL)
    {
    }

    IceBCFuncWrapper(BCFunc funcptr) : m_funcptr(funcptr)
    {
    }

    ///
    virtual ~IceBCFuncWrapper()
    {
    }

    /// simply calls through to bcFunc
    virtual void
    operator()(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous)
    {
        CH_assert(m_funcptr != NULL);
        m_funcptr(a_state, a_valid, a_domain, a_dx, a_homogeneous);
    }

    // explictly allow copy and assignment operators
    IceBCFuncWrapper(const IceBCFuncWrapper& src) : m_funcptr(src.m_funcptr)
    {
    }

    virtual IceBCFuncWrapper& operator=(const IceBCFuncWrapper& src)
    {
        m_funcptr = src.m_funcptr;
        return *this;
    }

protected:
    BCFunc m_funcptr;
};

/// Physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies initial and
   boundary conditions. Derived from PhysIBC, adds functions specific
   for this application.
 */
class BasicIBC : public PhysIBC
{
public:
    /// Constructor
    /**
     */
    BasicIBC()
    {
    }

    /// Destructor
    /**
     */
    virtual ~BasicIBC()
    {
    }

    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object. Default is to just call
       PhysIBC::define
     */
    virtual void define(const ProblemDomain& a_domain, const Real& a_dx)
    {
        PhysIBC::define(a_domain, a_dx);
    }

    /// Factory method - this object is its own factory
    /**
       Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
       its define() must be called before it is used).
     */
    virtual PhysIBC* new_physIBC()
    {
        return static_cast<PhysIBC*>(new_thicknessIBC());
    }

    /// same as new_physIBC, except pointer cast to an BasicIBC
    virtual BasicIBC* new_thicknessIBC() = 0;

    /// Set up initial conditions
    /**
     */
    virtual void initialize(LevelData<FArrayBox>& a_U) = 0;

    virtual void initializeData(RealVect& a_dx,
                             LevelData<FArrayBox>& a_Head,
                             LevelData<FArrayBox>& a_GapHeight,
                             LevelData<FArrayBox>& a_zbed, 
                             LevelData<FArrayBox>& a_Pi) = 0;

    /// Set boundary fluxes
    /**
     */
    virtual void primBC(FArrayBox& a_WGdnv,
                        const FArrayBox& a_Wextrap,
                        const FArrayBox& a_W,
                        const int& a_dir,
                        const Side::LoHiSide& a_side,
                        const Real& a_time) = 0;

    /// Set boundary slopes
    /**
       The boundary slopes in a_dW are already set to one sided difference
       approximations.  If this function doesn't change them they will be
       used for the slopes at the boundaries.
     */
    virtual void setBdrySlopes(FArrayBox& a_dW, const FArrayBox& a_W, const int& a_dir, const Real& a_time) = 0;

    /// Adjust boundary fluxes to account for artificial viscosity
    /**
     */
    virtual void artViscBC(FArrayBox& a_F,
                           const FArrayBox& a_U,
                           const FArrayBox& a_divVel,
                           const int& a_dir,
                           const Real& a_time) = 0;

    /// return boundary condition for Ice velocity solve
    /**
     */
    virtual RefCountedPtr<CompGridVTOBC> velocitySolveBC() = 0;

    /// fill ghost cells on velocity
    /** this is a separate function from the BCFunction returned by
        the velocitySolveBC function to take into account the case
        where the boundary condition used in the solve is different
        from the one used when working with the velocity. One example
        of this is the case where we incorporate the normal-stress boundary
        condition into the RHS for the velocity solve; in that case, the
        BCFunction will be a homogeneous BC. However we want to use the
        non-homogeneous form of the BC's when we do things like compute
        derivatives of the velocity at physical domain boundaries
        (say, for example, when computing the strain invariant in order
        to compute the viscosity).  This function fills ghost cells on the
        boundary in the generic case. Default implementation simply calls
        the BCFunction implementation... */
    virtual void
    velocityGhostBC(LevelData<FArrayBox>& a_velocity, const ProblemDomain& a_domain, const RealVect& a_dx, Real a_time)
    {
        RefCountedPtr<CompGridVTOBC> velBC = velocitySolveBC();
        for (DataIterator dit = a_velocity.dataIterator(); dit.ok(); ++dit)
            (*velBC)(a_velocity[dit], a_velocity.getBoxes()[dit], a_domain, a_dx[0], dit(), false);
    }

    /// set AMR grid hierarchy (for the BC's which need this)
    virtual void setGridHierarchy(Vector<DisjointBoxLayout>& a_vectGrids, Vector<ProblemDomain>& a_vectDomain)
    {
        ;
    }

protected:
private:
    // Disallowed for all the usual reasons
    void operator=(const BasicIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    BasicIBC(const BasicIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

/// BCFunction-derived class to handle velocity-solve BC's
/** */
class VelBCFunction : public CompGridVTOBC
{
public:
    VelBCFunction()
    {
        ;
    }

    VelBCFunction(Vector<ProblemDomain>& a_vectDomain)
    {
        ;
    }

    ~VelBCFunction()
    {
        ;
    }

    virtual void define(Vector<ProblemDomain>& a_vectDomain) = 0;

    virtual void operator()(FArrayBox& a_state,
                            const Box& a_valid,
                            const ProblemDomain& a_domain,
                            Real a_dx,
                            bool a_homogeneous) = 0;

protected:
    Vector<ProblemDomain> m_vectDomain;
};

#include "NamespaceFooter.H"
#endif
