#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "SqrtIBC.H"
#include "ParmParse.H"
#include "FluxBox.H"
#include <random>

#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
SqrtIBC::SqrtIBC()
{
}

SqrtIBC::~SqrtIBC()
{
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
SqrtIBC::define(const ProblemDomain& a_domain, const Real& a_dx)
{
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
HydroIBC*
SqrtIBC::new_hydroIBC()
{
    pout() << "Initializing a new SqrtIBC"<< endl; 
    SqrtIBC* retval = new SqrtIBC();

    return static_cast<HydroIBC*>(retval);
}

/** Set up mask 
 */
void
SqrtIBC::initializePi(RealVect& a_dx,
                       suhmo_params Params,     
                       LevelData<FArrayBox>& a_Pi)
{
    // Do nothing 
}

void 
SqrtIBC::initializePi(RealVect& a_dx,
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
SqrtIBC::initializeBed(RealVect& a_dx,
                       suhmo_params Params,     
                       LevelData<FArrayBox>& a_zbed,
                       LevelData<FArrayBox>& a_bumpHeight,
                       LevelData<FArrayBox>& a_bumpSpacing)
{

    pout() << "SqrtIBC::initializeBed" << endl;

    bool AGU_test_case = false;
    //if (AGU_test_case) {
        // randomization
        const double mean    = 0.0;
        const double stddev2 = 0.5;
        std::default_random_engine generator;
        std::default_random_engine generator2;
        std::normal_distribution<double> dist2(mean, stddev2);
    //}

    DataIterator dit = a_zbed.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thiszbed          = a_zbed[dit];
        FArrayBox& thisbumpHeight    = a_bumpHeight[dit];
        FArrayBox& thisbumpSpacing   = a_bumpSpacing[dit];

        BoxIterator bit(thiszbed.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];
            Real y_loc = (iv[1]+0.5)*a_dx[1];

            /* bed topography */
            if (AGU_test_case) {
                // add randomness
                thiszbed(iv, 0)      = std::max(Params.m_slope*x_loc + dist2(generator2) 
                                       + 30.0*std::sin(x_loc/10.0) 
                                       //+ 1000*std::sin(3.0*x_loc) 
                                       //+ 50*std::sin(x_loc+2.0)
                                       //+ 5000*std::sin(x_loc/100.0), 0.0); 
                                       //+ 500.0*std::sin(y_loc/50.0+1.0)
                                       //+ 0.01*(y_loc - 1000) * (y_loc - 1000) , 0.0);
                                       , 0.0);

                //if (x_loc < 400.0) {
                //    thiszbed(iv, 0)      = Params.m_slope*x_loc;
                //} else if  (x_loc < 1000.0) {
                //    thiszbed(iv, 0)      = std::max(thiszbed(iv, 0) , Params.m_slope*x_loc);
                //}
                thiszbed(iv, 0)      = thiszbed(iv, 0) / 10000.0;
            } else {
                // typical
                thiszbed(iv, 0)      = Params.m_slope*x_loc;
            }
            /* Ice height (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            //thisiceHeight(iv, 0) = 6.0 * (std::sqrt(x_loc + Params.m_H) - std::sqrt(Params.m_H)) + 1.0;
            /* Ice overburden pressure : rho_i * g * H */
            //thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * thisiceHeight(iv, 0);
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with SqrtIBC::initializeBed)" << endl;
    }

}

/** Set up initial conditions 
 */
void
SqrtIBC::initializeData(RealVect& a_dx,
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
                         LevelData<FArrayBox>& a_bumpSpacing)
{
    
    pout() << "SqrtIBC::initializeData" << endl;

    bool AGU_test_case = false;
    //if (AGU_test_case) {
        // randomization
        const double mean    = 0.0;
        const double stddev2 = 0.5;
        std::default_random_engine generator;
        std::default_random_engine generator2;
        std::normal_distribution<double> dist2(mean, stddev2);
    //}

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
            Real y_loc = (iv[1]+0.5)*a_dx[1];

            /* bed topography */
            if (AGU_test_case) {
                // add randomness
                thiszbed(iv, 0)      = std::max(Params.m_slope*x_loc + dist2(generator2) 
                                       + 30.0*std::sin(x_loc/10.0) 
                                       //+ 1000*std::sin(3.0*x_loc) 
                                       //+ 50*std::sin(x_loc+2.0)
                                       //+ 5000*std::sin(x_loc/100.0) 
                                       //+ 500.0*std::sin(y_loc/50.0+1.0)
                                       //+ 0.01*(y_loc - 1000) * (y_loc - 1000) 
                                       , 0.0);

                //if (x_loc < 400.0) {
                //    thiszbed(iv, 0)      = Params.m_slope*x_loc;
                //} else if  (x_loc < 1000.0) {
                //    thiszbed(iv, 0)      = std::max(thiszbed(iv, 0) , Params.m_slope*x_loc);
                //}
                thiszbed(iv, 0)      = thiszbed(iv, 0) / 10000.0;
            } else {
                // typical
                thiszbed(iv, 0)      = Params.m_slope*x_loc;
            }
            /* initial gap height */
            thisGapHeight(iv, 0) = Params.m_gapInit;
            /* Ice height (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            thisiceHeight(iv, 0) = 6.0 * (std::sqrt(x_loc + Params.m_H) - std::sqrt(Params.m_H)) + 1.0;
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
            //thisbumpHeight(iv, 0)      = std::max(Params.m_br + dist(generator), 0.0); 

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
        pout() << "(Done with SqrtIBC::initializeData)" << endl;
    }
}


/** Set up initial conditions 
 */
void
SqrtIBC::initialize(LevelData<FArrayBox>& a_head)
{
}

#include "NamespaceFooter.H"
