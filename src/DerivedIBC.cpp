#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "DerivedIBC.H"
#include "ParmParse.H"
#include "FluxBox.H"
#include "suhmo.H"

#include "NamespaceHeader.H"

void
zeroBCValue(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values)
{
    a_values[0] = 0.0;
}

void
iceNeumannBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous)
{
    if (!a_domain.domainBox().contains(a_state.box()))
    {
        Box valid = a_valid;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
            // don't do anything if periodic
            if (!a_domain.isPeriodic(dir))
            {
                Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
                Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
                if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                    // Real bcVal = 0.0;
                    NeumBC(a_state, valid, a_dx, a_homogeneous, zeroBCValue, dir, Side::Lo);
                }

                if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                    NeumBC(a_state, valid, a_dx, a_homogeneous, zeroBCValue, dir, Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}

void
iceDirichletBC(FArrayBox& a_state, const Box& a_valid, const ProblemDomain& a_domain, Real a_dx, bool a_homogeneous)
{
    if (!a_domain.domainBox().contains(a_state.box()))
    {
        Box valid = a_valid;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
            // don't do anything if periodic
            if (!a_domain.isPeriodic(dir))
            {
                Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
                Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
                if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                    // Real bcVal = 0.0;
                    DiriBC(a_state, valid, a_dx, a_homogeneous, zeroBCValue, dir, Side::Lo);
                }

                if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                    // Real bcVal = 0.0;
                    DiriBC(a_state, valid, a_dx, a_homogeneous, zeroBCValue, dir, Side::Hi);
                }

            } // end if is not periodic in ith direction
        }
    }
}

// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
DerivedIBC::DerivedIBC() : m_boundaryPhi(0.0)
{
    m_isBCsetUp = false;
    m_isDefined = false;
}

DerivedIBC::~DerivedIBC()
{
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
DerivedIBC::define(const ProblemDomain& a_domain, const Real& a_dx)
{
    pout() << "DerivedIBC::define" << endl;
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
BasicIBC*
DerivedIBC::new_thicknessIBC()
{
    DerivedIBC* retval = new DerivedIBC();

    retval->m_boundaryPhi = m_boundaryPhi;

    return static_cast<BasicIBC*>(retval);
}

/** Set up initial conditions 
 */
void
DerivedIBC::initialize(LevelData<FArrayBox>& a_head)
{
    pout() << "DerivedIBC::initialize" << endl;
    // for now, just initialize to a square subregion
    IntVect regionLo = m_domain.domainBox().bigEnd();
    IntVect regionHi = m_domain.domainBox().bigEnd();

    DataIterator dit = a_head.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox& thisHead = a_head[dit];

        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            if ((iv > regionLo) && (iv < regionHi))
            {
                thisHead(iv, 0) = 1.0;
            }
            else
            {
                thisHead(iv, 0) = 2.0;
            }
        } // end loop over cells
    }     // end loop over boxes
}

/** Set up initial conditions 
 */
void
DerivedIBC::initializeData(RealVect& a_dx,
                        LevelData<FArrayBox>& a_head,
                        LevelData<FArrayBox>& a_gapHeight,
                        LevelData<FArrayBox>& a_Pw,
                        LevelData<FArrayBox>& a_qw,
                        LevelData<FArrayBox>& a_Re,
                        LevelData<FArrayBox>& a_meltRate,
                        LevelData<FArrayBox>& a_zbed,
                        LevelData<FArrayBox>& a_Pi)
{
    
    pout() << "DerivedIBC::initializeData" << endl;

    ParmParse ppBC("icbc");
    ppBC.get("H", m_H);
    ppBC.get("slope", m_slope);
    ppBC.get("GapInit", m_gapInit);

    IntVect regionLo = m_domain.domainBox().bigEnd();
    IntVect regionHi = m_domain.domainBox().bigEnd();

    DataIterator dit = a_head.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox& thisHead      = a_head[dit];
        FArrayBox& thisGapHeight = a_gapHeight[dit];
        FArrayBox& thisPw        = a_Pw[dit];
        FArrayBox& thisqw        = a_qw[dit];
        FArrayBox& thisRe        = a_Re[dit];
        FArrayBox& thismeltRate  = a_meltRate[dit];
        FArrayBox& thiszbed      = a_zbed[dit];
        FArrayBox& thispi        = a_Pi[dit];

        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            // gapHeight - initially constant 
            thisGapHeight(iv, 0) = m_gapInit;
            // Initial Head is such that the water pressure is equal to 50% of the ice overburden pressure
            Real P_ice        = hydro_params::m_rho_i * m_gravity * m_H;
            Real P_water_init = P_ice * 0.5;
            Real Fact         = 1./(hydro_params::m_rho_w * m_gravity);
            thiszbed(iv, 0)      = m_slope*(iv[0]+0.5)*a_dx[0];
            thisHead(iv, 0)      = P_water_init * Fact + thiszbed(iv, 0) ;
            thispi(iv, 0)        = P_ice;
            thisPw(iv, 0)        = hydro_params::m_rho_w * m_gravity * (thisHead(iv, 0)  - thiszbed(iv, 0));
            // dummy stuff 
            thisqw(iv, 0)        = 0.0;
            thisRe(iv, 0)        = 0.0;
            thismeltRate(iv, 0)  = m_G / m_L;
        } // end loop over cells
    }     // end loop over boxes

    pout() << "Done with DerivedIBC::initializeData" << endl;
}

/// Set boundary fluxes
/**
 */
void
DerivedIBC::primBC(FArrayBox& a_WGdnv,
                   const FArrayBox& a_Wextrap,
                   const FArrayBox& a_W,
                   const int& a_dir,
                   const Side::LoHiSide& a_side,
                   const Real& a_time)
{
    // do nothing in periodic case
    if (!m_domain.isPeriodic(a_dir))
    {
        int lohisign;
        Box tmp = a_WGdnv.box();

        // Determine which side and thus shifting directions
        lohisign = sign(a_side);
        tmp.shiftHalf(a_dir, lohisign);

        // Is there a domain boundary next to this grid
        if (!m_domain.contains(tmp))
        {
            tmp &= m_domain;

            Box boundaryBox;

            if (a_side == Side::Lo)
            {
                boundaryBox = bdryLo(tmp, a_dir);
            }
            else
            {
                boundaryBox = bdryHi(tmp, a_dir);
            }

            // Set the boundary values
            a_WGdnv.setVal(m_boundaryPhi, boundaryBox, 0, 1);
        }
    }
}

/// Set boundary slopes
/**
   The boundary slopes in a_dW are already set to one sided difference
   approximations.  If this function doesn't change them they will be
   used for the slopes at the boundaries.
*/

void
DerivedIBC::setBdrySlopes(FArrayBox& a_dW, const FArrayBox& a_W, const int& a_dir, const Real& a_time)
{
    // one-sided differences sounds fine with me, so do nothing...
}

/// Adjust boundary fluxes to account for artificial viscosity
/**
 */
void
DerivedIBC::artViscBC(FArrayBox& a_F,
                      const FArrayBox& a_U,
                      const FArrayBox& a_divVel,
                      const int& a_dir,
                      const Real& a_time)
{
    // don't anticipate being here -- if we wind up here, need to
    // give it some thought
    MayDay::Error("DerivedIBC::artViscBC not implemented");
}

/// return boundary condition for Ice velocity solve
/** eventually would like this to be a BCHolder
 */
RefCountedPtr<CompGridVTOBC>
DerivedIBC::velocitySolveBC()
{
    if (!m_isBCsetUp)
    {
        setupBCs();
    }

    return m_velBCs;
}

void
DerivedIBC::setupBCs()
{
    ParmParse ppBC("bc");

    // get boundary conditions
    Vector<int> loBCvect(SpaceDim), hiBCvect(SpaceDim);
    ppBC.getarr("lo_bc", loBCvect, 0, SpaceDim);
    ppBC.getarr("hi_bc", hiBCvect, 0, SpaceDim);

    // this is a placeholder until I can get a BCHolder to work...
    // require all directions to have the same BC for now
    CH_assert(loBCvect[0] == loBCvect[1]);
    CH_assert(hiBCvect[0] == hiBCvect[1]);
    CH_assert(loBCvect[0] == hiBCvect[0]);

    if (loBCvect[0] == 0)
    {
        // BCFuncWrapper* newBCPtr = new BCFuncWrapper(iceDirichletBC);
        m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceDirichletBC));
    }
    else if (loBCvect[0] == 1)
    {
        m_velBCs = RefCountedPtr<CompGridVTOBC>(new IceBCFuncWrapper(iceNeumannBC));
    }
    else
    {
        MayDay::Error("bad BC type");
    }
    m_isBCsetUp = true;
}

#include "NamespaceFooter.H"
