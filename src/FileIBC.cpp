#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "FileIBC.H"
#include "ParmParse.H"
#include "FluxBox.H"
#include "FillFromReference.H"
#include "ReadLevelData.H"
#include <random>

#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
// set default thickness at domain edge to be zero
FileIBC::FileIBC()
{
}

FileIBC::~FileIBC()
{
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this
   initial and boundary condition object. Just calls base-class define
*/
void
FileIBC::define(const ProblemDomain& a_domain, const Real& a_dx)
{
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
HydroIBC*
FileIBC::new_hydroIBC()
{
    pout() << "Initializing a new FileIBC"<< endl; 
    FileIBC* retval = new FileIBC();
    retval->m_refDx = m_refDx;
    retval->m_refThickness = m_refThickness;
    retval->m_refTopography = m_refTopography;
    
    return static_cast<HydroIBC*>(retval);
}


/* set up file input and names */
void FileIBC::setup(std::string a_geomFile,
                    std::string a_thicknessName,
                    std::string a_topographyName)
{
  /* read level data from external source */
  RefCountedPtr<LevelData<FArrayBox> > levelThck(new LevelData<FArrayBox>());
  RefCountedPtr<LevelData<FArrayBox> > levelTopg(new LevelData<FArrayBox>());
  
  Real dx;
  
  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectData;
  vectData.push_back(levelThck);
  vectData.push_back(levelTopg);
  
  Vector<std::string> names(2);
  names[0] = a_thicknessName;
  names[1] = a_topographyName;
  readLevelData(vectData,dx,a_geomFile,names,1);

  /* now copy into persistent storage */
  m_refDx         = dx*RealVect::Unit;
  m_refThickness  = levelThck;   // ice thickness
  m_refTopography = levelTopg;  // bed topography

}

/** Set up mask 
 */
void
FileIBC::initializePi(RealVect& a_dx,
                      suhmo_params Params,     
                      LevelData<FArrayBox>& a_Pi)
{
    // Do nothing 
}

void
FileIBC::initializePi(RealVect& a_dx,
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
FileIBC::initializeBed(RealVect& a_dx,
                       suhmo_params Params,     
                       LevelData<FArrayBox>& a_zbed,
                       LevelData<FArrayBox>& a_bumpHeight,
                       LevelData<FArrayBox>& a_bumpSpacing)
{
    // Do nothing 
}


/** Set up initial conditions 
 */
void
FileIBC::initializeData(RealVect& a_dx,
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
    
    pout() << "FileIBC::initializeData" << endl;

    bool verbose = true;

    /* initialize bed and ice thickness */
    const LevelData<FArrayBox>& topoRef = *(m_refTopography);
    const LevelData<FArrayBox>& thickRef = *(m_refThickness);    
    FillFromReference(a_zbed, topoRef, a_dx, m_refDx,verbose);
    FillFromReference(a_iceHeight, thickRef, a_dx, m_refDx, verbose);

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

        BoxIterator bit(thisHead.box()); 
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real x_loc = (iv[0]+0.5)*a_dx[0];

            /* bed topography --> done  */
            // add randomness
            //thiszbed(iv, 0)      = std::max(thiszbed(iv, 0) + dist2(generator), 0.0);
            /* Ice height  --> done*/
            /* Ice overburden pressure : rho_i * g * H */
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
        pout() << "(Done with FileIBC::initializeData)" << endl;
    }
}

/** Set up initial conditions 
 */
void
FileIBC::initialize(LevelData<FArrayBox>& a_head)
{
}

#include "NamespaceFooter.H"
