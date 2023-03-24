#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ValleyIBC.H"
#include "ParmParse.H"
#include "FluxBox.H"
#include <random>

#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
ValleyIBC::ValleyIBC()
{
}

ValleyIBC::~ValleyIBC()
{
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
ValleyIBC::define(const ProblemDomain& a_domain, const Real& a_dx)
{
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
HydroIBC*
ValleyIBC::new_hydroIBC()
{
    pout() << "Initializing a new ValleyIBC"<< endl; 
    ValleyIBC* retval = new ValleyIBC();
    retval->m_gamma = m_gamma;

    return static_cast<HydroIBC*>(retval);
}

/* Set Parameters */
void
ValleyIBC::setParameters(const Real& a_gamma) {
    m_gamma = a_gamma;
}

/** Set up mask 
 */
void
ValleyIBC::initializePi(RealVect& a_dx,
                        suhmo_params Params,     
                        LevelData<FArrayBox>& a_Pi)
{
    
    pout() << "ValleyIBC::initializePi" << endl;

    DataIterator dit = a_Pi.dataIterator();

    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thispi        = a_Pi[dit];

        BoxIterator bit(thispi.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];
            Real y_loc = (iv[1]+0.5)*a_dx[1] - 750.0;

            /* bed topography */
            /* Ice height = ICE from 0 (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            Real IceHeight = 100.0 * std::pow(x_loc + 200.0, 0.25) + x_loc/60.0 - std::pow(2.0e10, 0.25) + 1.0;
            // Weird bed
            Real gamma_b         = 0.05;
            Real H_6000          = 100.0 * std::pow(6000. + 200.0, 0.25) + 6000./60.0 - std::pow(2.0e10, 0.25) + 1.0;
            Real fx              = (H_6000 - 6000.*m_gamma) * x_loc * x_loc / (6000.0*6000.0) + m_gamma*x_loc;
            Real fx_gb           = (H_6000 - 6000.*gamma_b) * x_loc * x_loc / (6000.0*6000.0) + gamma_b*x_loc;
            Real gy              = 0.5e-6 * std::abs(y_loc*y_loc*y_loc); 
            Real hx              = (-4.5 * x_loc / 6000. + 5.0) * (IceHeight - fx) / (IceHeight - fx_gb + 1.0e-16);
            Real Bed             = fx + gy*hx;

            /* Ice overburden pressure : rho_i * g * H WITH H = Ice height - (bed + gap)*/
            thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * std::max(IceHeight - Bed, 0.0);
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with ValleyIBC::initializePi)" << endl;
    }
}

void
ValleyIBC::initializePi(RealVect& a_dx,
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
    
    pout() << "ValleyIBC::initializePi" << endl;

    DataIterator dit = a_Pi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thisHead      = a_head[dit];
        FArrayBox& thisGapHeight = a_gapHeight[dit];
        FArrayBox& thiszbed      = a_zbed[dit];
        FArrayBox& thispi        = a_Pi[dit];
        FArrayBox& thisPw        = a_Pw[dit];
        FArrayBox& thisiceHeight = a_iceHeight[dit];
        FArrayBox& thisbumpHeight    = a_bumpHeight[dit];
        FArrayBox& thisbumpSpacing   = a_bumpSpacing[dit];

        BoxIterator bit(thispi.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];
            Real y_loc = (iv[1]+0.5)*a_dx[1] - 750.0;

            /* bed topography */
            /* Ice height = ICE from 0 (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            thisiceHeight(iv, 0) = 100.0 * std::pow(x_loc + 200.0, 0.25) + x_loc/60.0 - std::pow(2.0e10, 0.25) + 1.0;
            // Weird bed
            Real gamma_b         = 0.05;
            Real H_6000          = 100.0 * std::pow(6000. + 200.0, 0.25) + 6000./60.0 - std::pow(2.0e10, 0.25) + 1.0;
            Real fx              = (H_6000 - 6000.*m_gamma) * x_loc * x_loc / (6000.0*6000.0) + m_gamma*x_loc;
            Real fx_gb           = (H_6000 - 6000.*gamma_b) * x_loc * x_loc / (6000.0*6000.0) + gamma_b*x_loc;
            Real gy              = 0.5e-6 * std::abs(y_loc*y_loc*y_loc); 
            Real hx              = (-4.5 * x_loc / 6000. + 5.0) * (thisiceHeight(iv, 0) - fx) / (thisiceHeight(iv, 0) - fx_gb + 1.0e-16);
            thiszbed(iv, 0)      = fx + gy*hx;
            /* Ice overburden pressure : rho_i * g * H WITH H = Ice height - (bed + gap)*/
            thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * std::max(thisiceHeight(iv, 0) - thiszbed(iv, 0), 0.0);

            /* initial gap height */
            //if (thispi(iv, 0) == 0.0) {
            //    thisGapHeight(iv, 0) = 1.0e-16;
            //} 

        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with ValleyIBC::initializePi)" << endl;
    }
}


void 
ValleyIBC::initializeBed(RealVect& a_dx,
                         suhmo_params Params,     
                         LevelData<FArrayBox>& a_zbed,
                         LevelData<FArrayBox>& a_bumpHeight,
                         LevelData<FArrayBox>& a_bumpSpacing)
{
    // Do nothing 
}


void 
ValleyIBC::resetCovered(suhmo_params Params,     
                        LevelData<FArrayBox>& a_head,
                        LevelData<FArrayBox>& a_Pi)
{
    // Clean neg h in unused covered zones
    DataIterator dit = a_Pi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thispi        = a_Pi[dit];
        FArrayBox& thisHead      = a_head[dit];

        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            if (thisHead(iv, 0) < 0.0) {
                thisHead(iv, 0) = 0.0;
            }
        }
    }
}


/** Set up initial conditions 
 */
void
ValleyIBC::initializeData(RealVect& a_dx,
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
    
    pout() << "ValleyIBC::initializeData" << endl;

    // randomization
    const double mean = 0.0;
    const double stddev  = 0.05;
    const double stddev2 = 0.001;
    std::default_random_engine generator;
    std::default_random_engine generator2;
    std::normal_distribution<double> dist(mean, stddev);
    std::normal_distribution<double> dist2(mean, stddev2);

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

        FArrayBox& thismagVel    = a_levelmagVel[dit]; 


        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];
            Real y_loc = (iv[1]+0.5)*a_dx[1] - 750.0;

            /* bed topography */
            /* Ice height = ICE from 0 (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            thisiceHeight(iv, 0) = 100.0 * std::pow(x_loc + 200.0, 0.25) + x_loc/60.0 - std::pow(2.0e10, 0.25) + 1.0;
            // Weird bed
            Real gamma_b         = 0.05;
            Real H_6000          = 100.0 * std::pow(6000. + 200.0, 0.25) + 6000./60.0 - std::pow(2.0e10, 0.25) + 1.0;
            Real fx              = (H_6000 - 6000.*m_gamma) * x_loc * x_loc / (6000.0*6000.0) + m_gamma*x_loc;
            Real fx_gb           = (H_6000 - 6000.*gamma_b) * x_loc * x_loc / (6000.0*6000.0) + gamma_b*x_loc;
            Real gy              = 0.5e-6 * std::abs(y_loc*y_loc*y_loc); 
            Real hx              = (-4.5 * x_loc / 6000. + 5.0) * (thisiceHeight(iv, 0) - fx) / (thisiceHeight(iv, 0) - fx_gb + 1.0e-16);
            thiszbed(iv, 0)      = fx + gy*hx;
            /* Ice overburden pressure : rho_i * g * H WITH H = Ice height - (bed + gap)*/
            thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * std::max(thisiceHeight(iv, 0) - thiszbed(iv, 0), 0.0);
            //Real smooth_Pi       = Params.m_rho_i * Params.m_gravity * std::max(thisiceHeight(iv, 0) - (fx + gy*hx), 0.0);
            /* initial gap height */
            if (thispi(iv, 0) == 0.0) {
                thisGapHeight(iv, 0) = 1.0e-16;
                //thiszbed(iv, 0) = thisiceHeight(iv, 0);
                //thisiceHeight(iv, 0) = 0.0;
            } else {
                //thiszbed(iv, 0)      = std::max(fx + gy*hx + dist2(generator), 0.0);
                //thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * std::max(thisiceHeight(iv, 0) - thiszbed(iv, 0), 0.0);
                thisGapHeight(iv, 0) = Params.m_gapInit;
                //thisGapHeight(iv, 0) = std::max(Params.m_gapInit + dist2(generator2), 0.0001);

            }
            
            /* option 1: guess Pw, find head */
            // Water press ?? No idea --> Pi/2.0
            thisPw(iv, 0)        = thispi(iv, 0) * 0.5;
            Real Fact            = 1./(Params.m_rho_w * Params.m_gravity);
            thisHead(iv, 0)      = thisPw(iv, 0) * Fact + thiszbed(iv, 0) ;
            //if (thispi(iv, 0) == 0.0) {
            //    thisHead(iv, 0)      = 0.0;
            //}

            /* Bed randomization */
            // typical
            thisbumpHeight(iv, 0)      = Params.m_br;
            thisbumpSpacing(iv, 0)     = Params.m_lr;
            // if randomness
            //thisbumpHeight(iv, 0)      = std::min(std::max(Params.m_br + dist(generator), 0.0), 0.1); 
            //thisbumpHeight(iv, 0)      = std::max(thisbumpHeight(iv, 0)  + 0.01*(thiszbed(iv, 0) - (fx + gy*hx)), 0.0);
            if (thispi(iv, 0) == 0.0) {
                thisbumpHeight(iv, 0)      = 0.1 * Params.m_br;
            }

            /* dummy stuff -- will be erased right away */
            thisRe(iv, 0)        = std::max(100.0 - 0.1*Params.m_ReInit*x_loc, 10.0);
            Real denom_q         = 12.0 * Params.m_nu * (1 + Params.m_omega * Params.m_ReInit);
            // grad head = slope for now
            Real num_q           = - Params.m_gravity * std::pow(thisGapHeight(iv, 0),3) * Params.m_slope ;
            thisqw(iv, 0)        = num_q/denom_q; 
            thisqw(iv, 1)        = 0.0;
            thismeltRate(iv, 0)  = (Params.m_G/Params.m_L); 

            thismagVel(iv, 0)  = std::sqrt(  Params.m_ub[0]*Params.m_ub[0] 
                                           + Params.m_ub[1]*Params.m_ub[1] ); 
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with ValleyIBC::initializeData)" << endl;
    }
}


/** Set up initial conditions 
 */
void
ValleyIBC::initialize(LevelData<FArrayBox>& a_head)
{
}

#include "NamespaceFooter.H"
