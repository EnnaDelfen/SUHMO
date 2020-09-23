#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "HydroIBC.H"
#include "ParmParse.H"
#include "FluxBox.H"
//#include "AmrHydro.H"

#include "NamespaceHeader.H"

void
zeroBCValue(Real* pos, 
            int* dir, 
            Side::LoHiSide* side, 
            Real* a_values)
{
    a_values[0] = 0.0;
}

// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
HydroIBC::HydroIBC()
{
}

HydroIBC::~HydroIBC()
{
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
HydroIBC::define(const ProblemDomain& a_domain, const Real& a_dx)
{
    pout() << "HydroIBC::define" << endl;
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
HydroIBC*
HydroIBC::new_hydroIBC()
{
    HydroIBC* retval = new HydroIBC();

    return retval;
}

/** Set up initial conditions 
 */
void
HydroIBC::initializeData(RealVect& a_dx,
                           suhmo_params Params,     
                           LevelData<FArrayBox>& a_head,
                           LevelData<FArrayBox>& a_gapHeight,
                           LevelData<FArrayBox>& a_Pw,
                           LevelData<FArrayBox>& a_qw,
                           LevelData<FArrayBox>& a_Re,
                           LevelData<FArrayBox>& a_meltRate,
                           LevelData<FArrayBox>& a_zbed,
                           LevelData<FArrayBox>& a_Pi,
                           LevelData<FArrayBox>& a_iceHeight)
{
    
    pout() << "HydroIBC::initializeData" << endl;

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
        FArrayBox& thisiceHeight = a_iceHeight[dit];

        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            // gapHeight - initially constant 
            thisGapHeight(iv, 0) = Params.m_gapInit;
            // Initial Head is such that the water pressure is equal to 50% of the ice overburden pressure
            thisiceHeight(iv, 0) = Params.m_H;
            Real P_ice        = Params.m_rho_i * Params.m_gravity * Params.m_H;
            Real P_water_init = P_ice * 0.5;
            Real Fact         = 1./(Params.m_rho_w * Params.m_gravity);
            thiszbed(iv, 0)      = Params.m_slope*(iv[0]+0.5)*a_dx[0];
            thisHead(iv, 0)      = P_water_init * Fact + thiszbed(iv, 0) ;
            thispi(iv, 0)        = P_ice; // cst for now ...
            thisPw(iv, 0)        = P_water_init; //Params.m_rho_w * Params.m_gravity * (thisHead(iv, 0)  - thiszbed(iv, 0));
            // dummy stuff 
            thisRe(iv, 0)        = Params.m_ReInit;
            Real denom_q         = 12.0 * Params.m_nu * (1 + Params.m_omega * Params.m_ReInit);
            Real num_q           = - Params.m_gravity * Params.m_gapInit * Params.m_gapInit * Params.m_gapInit * Params.m_slope ;
            thisqw(iv, 0)        = num_q/denom_q; 
            thisqw(iv, 1)        = 0.0;
            thismeltRate(iv, 0)  = (Params.m_G - Params.m_gravity * Params.m_rho_w * thisqw(iv, 0) * Params.m_slope)/ Params.m_L;
        } // end loop over cells
    }     // end loop over boxes

    pout() << "(Done with HydroIBC::initializeData)" << endl;
}


/// Set boundary fluxes
/**
 */
void
HydroIBC::primBC(FArrayBox& a_WGdnv,
                 const FArrayBox& a_Wextrap,
                 const FArrayBox& a_W,
                 const int& a_dir,
                 const Side::LoHiSide& a_side,
                 const Real& a_time)
{
}

/** Set up initial conditions 
 */
void
HydroIBC::initialize(LevelData<FArrayBox>& a_head)
{
}

void
HydroIBC::setBdrySlopes(FArrayBox& a_dW, const FArrayBox& a_W, const int& a_dir, const Real& a_time)
{
    // one-sided differences sounds fine with me, so do nothing...
}

/// Adjust boundary fluxes to account for artificial viscosity
/**
 */
void
HydroIBC::artViscBC(FArrayBox& a_F,
                    const FArrayBox& a_U,
                    const FArrayBox& a_divVel,
                    const int& a_dir,
                    const Real& a_time)
{
    // don't anticipate being here -- if we wind up here, need to
    // give it some thought
    MayDay::Error("HydroIBC::artViscBC not implemented");
}

#include "NamespaceFooter.H"
