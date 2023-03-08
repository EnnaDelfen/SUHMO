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
//#include <random>

#include "NamespaceHeader.H"

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
    pout() << "Initializing a new HydroIBC"<< endl; 
    HydroIBC* retval = new HydroIBC();

    return retval;
}

/** Set up mask 
 */
void
HydroIBC::initializePi(RealVect& a_dx,
                       suhmo_params Params,     
                       LevelData<FArrayBox>& a_Pi)
{
    // Do nothing 
}

void 
HydroIBC::initializePi(RealVect& a_dx,
                       suhmo_params Params,     
                       LevelData<FArrayBox>& a_head,
                       LevelData<FArrayBox>& a_gapHeight,
                       LevelData<FArrayBox>& a_Pw,
                       LevelData<FArrayBox>& a_zbed,
                       LevelData<FArrayBox>& a_Pi,
                       LevelData<FArrayBox>& a_iceHeight,
                       LevelData<FArrayBox>& a_bumpHeight,
                       LevelData<FArrayBox>& a_bumpSpacing)
{
    // Do nothing 
}

void 
HydroIBC::initializeBed(RealVect& a_dx,
                        suhmo_params Params,     
                        LevelData<FArrayBox>& a_zbed,
                        LevelData<FArrayBox>& a_bumpHeight,
                        LevelData<FArrayBox>& a_bumpSpacing)
{

    pout() << "HydroIBC::initializeBed" << endl;

    DataIterator dit = a_zbed.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thiszbed          = a_zbed[dit];
        FArrayBox& thisbumpHeight    = a_bumpHeight[dit];
        FArrayBox& thisbumpSpacing   = a_bumpSpacing[dit];

        BoxIterator bit(thiszbed.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];

            /* bed topography */
            // typical
            thiszbed(iv, 0)      = Params.m_slope*x_loc;
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with HydroIBC::initializeBed)" << endl;
    }
}

void 
HydroIBC::resetCovered(suhmo_params Params,     
                       LevelData<FArrayBox>& a_head,
                       LevelData<FArrayBox>& a_Pi)
{
    // Do nothing 
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
                         LevelData<FArrayBox>& a_iceHeight,
                         LevelData<FArrayBox>& a_bumpHeight,
                         LevelData<FArrayBox>& a_bumpSpacing,
                         LevelData<FArrayBox>& a_levelmagVel)
{
    
    pout() << "HydroIBC::initializeData" << endl;

    // randomization
    //const double mean2 = 0.0;
    //const double stddev2 = 0.1;
    //std::default_random_engine generator;
    //std::normal_distribution<double> dist2(mean2, stddev2);

    DataIterator dit = a_head.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thisHead      = a_head[dit];
        FArrayBox& thisGapHeight = a_gapHeight[dit];
        FArrayBox& thisPw        = a_Pw[dit];
        FArrayBox& thisqw        = a_qw[dit];
        FArrayBox& thisRe        = a_Re[dit];
        FArrayBox& thismeltRate  = a_meltRate[dit];
        FArrayBox& thiszbed      = a_zbed[dit];
        FArrayBox& thispi        = a_Pi[dit];
        FArrayBox& thisiceHeight = a_iceHeight[dit];
        FArrayBox& thisbumpHeight    = a_bumpHeight[dit];
        FArrayBox& thisbumpSpacing   = a_bumpSpacing[dit];


        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];

            /* bed topography */
            // typical
            thiszbed(iv, 0)      = Params.m_slope*x_loc;
            // add randomness
            //thiszbed(iv, 0)      = std::max(Params.m_slope*x_loc + dist2(generator), 0.0);
            /* initial gap height */
            thisGapHeight(iv, 0) = Params.m_gapInit;
            /* Ice height (should be ice only, so surface - (bed + gap)) */
            thisiceHeight(iv, 0) = Params.m_H;
            /* Ice overburden pressure : rho_i * g * H */
            thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * thisiceHeight(iv, 0);
            
            /* option 1: guess Pw, find head */
            // Water press ?? No idea --> Pi/2.0
            thisPw(iv, 0)        = thispi(iv, 0) * 0.5;
            Real Fact            = 1./(Params.m_rho_w * Params.m_gravity);
            thisHead(iv, 0)      = thisPw(iv, 0) * Fact + thiszbed(iv, 0) ;

            /* Bed randomization */
            // typical
            thisbumpHeight(iv, 0)      = Params.m_br;
            thisbumpSpacing(iv, 0)     = Params.m_lr;
            // if randomness
            //thisbumpHeight(iv, 0)      = std::max(thiszbed(iv, 0) - Params.m_slope*x_loc, 0.0); 

            /* dummy stuff -- will be erased right away */
            thisRe(iv, 0)        = std::max(100.0 - 0.1*Params.m_ReInit*x_loc, 10.0);
            Real denom_q         = 12.0 * Params.m_nu * (1 + Params.m_omega * Params.m_ReInit);
            // grad head = slope for now
            Real num_q           = - Params.m_gravity * std::pow(thisGapHeight(iv, 0),3) * Params.m_slope ;
            thisqw(iv, 0)        = num_q/denom_q; 
            thisqw(iv, 1)        = 0.0;
            thismeltRate(iv, 0)  = (Params.m_G/Params.m_L); 
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with HydroIBC::initializeData)" << endl;
    }
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
