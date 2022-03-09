#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "MountainSetupIBC.H"
#include "ParmParse.H"
#include "FluxBox.H"
#include <random>

#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
MountainIBC::MountainIBC()
{
}

MountainIBC::~MountainIBC()
{
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
MountainIBC::define(const ProblemDomain& a_domain, const Real& a_dx)
{
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
HydroIBC*
MountainIBC::new_hydroIBC()
{
    pout() << "Initializing a new MountainIBC"<< endl; 
    MountainIBC* retval = new MountainIBC();
    retval->m_gamma = m_gamma;

    return static_cast<HydroIBC*>(retval);
}

/* Set Parameters */
void
MountainIBC::setParameters(const Real& a_gamma) {
    m_gamma = a_gamma;
}

/** Set up mask 
 */
void
MountainIBC::initializePi(RealVect& a_dx,
                        suhmo_params Params,     
                        LevelData<FArrayBox>& a_Pi)
{
    // Do nothing 
}

void
MountainIBC::initializePi(RealVect& a_dx,
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
    DataIterator dit = a_head.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& thisHead      = a_head[dit];
        FArrayBox& thisGapHeight = a_gapHeight[dit];
        FArrayBox& thispi        = a_Pi[dit];
        FArrayBox& thisiceHeight = a_iceHeight[dit];

        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];
            Real y_loc = (iv[1]+0.5)*a_dx[1];
            /* bed topography */
            Real ax       = 1.5e-3;
            Real by       = -1.5e-3;
            Real cst      = 100.0;
            Real step_1   = ax*x_loc + by*y_loc + cst + 100.0 ;
            /* Ice height = ICE from 0 (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            thisiceHeight(iv, 0) = step_1;
            /* Ice overburden pressure : rho_i * g * H WITH H = Ice height - (bed + gap)*/
            thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * std::max(thisiceHeight(iv, 0), 0.0);
            /* initial gap height */
            if (thispi(iv, 0) == 0.0) {
                thisGapHeight(iv, 0) = 1.0e-16;
            } 
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with MountainIBC::initializePi)" << endl;
    }
}


void 
MountainIBC::initializeBed(RealVect& a_dx,
                         suhmo_params Params,     
                         LevelData<FArrayBox>& a_zbed,
                         LevelData<FArrayBox>& a_bumpHeight,
                         LevelData<FArrayBox>& a_bumpSpacing)
{
    // Do nothing 
}


void 
MountainIBC::resetCovered(suhmo_params Params,     
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
MountainIBC::initializeData(RealVect& a_dx,
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
    
    pout() << "MountainIBC::initializeData" << endl;

    // randomization
    //const double mean2 = 0.0;
    //const double stddev2 = 10;
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
            Real y_loc = (iv[1]+0.5)*a_dx[1];

            /* bed topography */

            // Weird bed
            //STEP 1 -- tilted base 
            Real ax       = 1.5e-3;
            Real by       = -1.5e-3;
            Real cst      = 100.0;
            Real step_1   = std::max(ax*x_loc + by*y_loc + cst, 0.0);
            /* Ice height = ICE from 0 (should be ice only, so surface - (bed + gap)) */
            // parabolic profile
            thisiceHeight(iv, 0) = ax*x_loc + by*y_loc + cst + 100.0 ; //Params.m_H;

            //STEP 2 -- middle finger of the dino
            Real A_max = 2.0*(10000 - (100000-x_loc + y_loc)*(100000 - x_loc + y_loc)/2000000.0);
            Real step_2           = std::max( ( A_max - ( ( x_loc - 50000 + (y_loc - x_loc)/2.0 )*( x_loc - 50000 + (y_loc - x_loc)/2.0 )/30000  
                                                      + ( y_loc - 50000 - (y_loc - x_loc)/2.0 )*( y_loc - 50000 - (y_loc - x_loc)/2.0 )/30000 )) / 1e2, 0.0);

            //STEP 3 -- side finger
            A_max = .5*(10000000 - y_loc*y_loc*y_loc*y_loc/1e12 + 40*y_loc*y_loc/8e3 );
            Real step_3 = std::max( (A_max - ( ( x_loc - 100000 )*( x_loc - 100000 )/40 )) / 4e4, 0.0);


            //STEP 4 -- side finger 
            A_max = .8*(+ x_loc*x_loc*x_loc*x_loc/1e12 - 50*x_loc*x_loc/7e3 );
            //A_max = .8*(-x_loc*x_loc*x_loc*x_loc/1e12 + 50*x_loc*x_loc/7e3 );
            Real step_4 = std::max( (A_max - y_loc *y_loc/50 ) / 5e4, 0.0);

            //STEP 5 -- other side finger
            A_max = .8*(-2.*(x_loc-30000)*(x_loc-30000)*(x_loc-30000)*(x_loc-30000)/1e12 + 70*(x_loc)*(x_loc)/7e3 );
            Real step_5 = std::max( (A_max - (y_loc+2000) *(y_loc+2000)/20 ) / 17e4, 0.0);
            //save
            //A_max = .8*(-x_loc*x_loc*x_loc*x_loc/1e12 + 50*x_loc*x_loc/7e3 );
            //Real step_5 = std::max( (A_max - y_loc *y_loc/40 ) / 3e4, 0.0);

            // Putting all together
            thiszbed(iv, 0)    = step_1 + std::max(std::max(step_3,step_2), std::max(step_4, step_5));
             
            /* Ice overburden pressure : rho_i * g * H WITH H = Ice height - (bed + gap)*/
            thispi(iv, 0)        = Params.m_rho_i * Params.m_gravity * std::max(thisiceHeight(iv, 0), 0.0);
            /* initial gap height */
            if (thispi(iv, 0) == 0.0) {
                thisGapHeight(iv, 0) = 1.0e-16;
            } else {
                thisGapHeight(iv, 0) = Params.m_gapInit;
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
        } // end loop over cells
    }     // end loop over boxes

    if (Params.m_verbosity > 3) {
        pout() << "(Done with MountainIBC::initializeData)" << endl;
    }
}


/** Set up initial conditions 
 */
void
MountainIBC::initialize(LevelData<FArrayBox>& a_head)
{
}

#include "NamespaceFooter.H"
