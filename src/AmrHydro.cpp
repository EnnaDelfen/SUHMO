#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>

using std::ifstream;
using std::ios;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::string;

#include "FORT_PROTO.H"
#include "AmrHydro.H"
#include "AmrHydroF_F.H"

#include "Box.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "parstream.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "PiecewiseLinearFillPatch.H"
#include "QuadCFInterp.H"
#include "CoarseAverageFace.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "DivergenceF_F.H"
#include "AMRUtilF_F.H"
#include "CH_HDF5.H"
#include "computeNorm.H" 
#include "MayDay.H"
#include "CONSTANTS.H"
#include "Gradient.H"
#include "ExtrapGhostCells.H"

#include "AMRFASMultiGrid.H"
#include "VCAMRNonLinearPoissonOp.H"

#include "NamespaceHeader.H"

// small parameter defining when times are equal
#define TIME_EPS 1.0e-12
// small parameter defining when a norm is "zero"
#define TINY_NORM 1.0e-8


/* BC FOR HEAD */

std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
bool              GlobalBCRS::s_areBCsParsed= false;
Real              GlobalBCRS::s_xlo_diri= -1000;
Real              GlobalBCRS::s_xhi_diri= -1000;
Real              GlobalBCRS::s_ylo_diri= -1000;
Real              GlobalBCRS::s_yhi_diri= -1000;
Real              GlobalBCRS::s_xlo_neum= -1000;
Real              GlobalBCRS::s_xhi_neum= -1000;
Real              GlobalBCRS::s_ylo_neum= -1000;
Real              GlobalBCRS::s_yhi_neum= -1000;

/* Parse boundary conditions and values in input file -- only once */
void ParseBC() 
{
    ParmParse ppBC("bc");
    ppBC.getarr("lo_bc", GlobalBCRS::s_bcLo, 0, CH_SPACEDIM);
    ppBC.getarr("hi_bc", GlobalBCRS::s_bcHi, 0, CH_SPACEDIM);

    ParmParse ppAmr("AmrHydro");
    std::vector<int> isPerio(CH_SPACEDIM, 0);
    ppAmr.getarr("is_periodic", isPerio, 0, CH_SPACEDIM);

    ParmParse pp;
    if (isPerio[0] == 0) {
        // x lo
        if (GlobalBCRS::s_bcLo[0] == 0) {
            pp.get("x.lo_dirich_val",GlobalBCRS::s_xlo_diri);
        } else if (GlobalBCRS::s_bcLo[0] == 1) {
            pp.get("x.lo_neumann_val",GlobalBCRS::s_xlo_neum);
        }
        // x hi
        if (GlobalBCRS::s_bcHi[0] == 0) {
            pp.get("x.hi_dirich_val",GlobalBCRS::s_xhi_diri);
        } else if (GlobalBCRS::s_bcHi[0] == 1) {
            pp.get("x.hi_neumann_val",GlobalBCRS::s_xhi_neum);
        }
    }
    if (isPerio[1] == 0) {
        // y lo
        if (GlobalBCRS::s_bcLo[1] == 0) {
            pp.get("y.lo_dirich_val",GlobalBCRS::s_ylo_diri);
        } else if (GlobalBCRS::s_bcLo[1] == 1) {
            pp.get("y.lo_neumann_val",GlobalBCRS::s_ylo_neum);
        }
        // y hi
        if (GlobalBCRS::s_bcHi[1] == 0) {
            pp.get("y.hi_dirich_val",GlobalBCRS::s_yhi_diri);
        } else if (GlobalBCRS::s_bcHi[1] == 1) {
            pp.get("y.hi_neumann_val",GlobalBCRS::s_yhi_neum);
        }
    }
    GlobalBCRS::s_areBCsParsed = true;
}


/* Return Neumann boundary val -- based on parsed value */
void NeumannValue(Real* pos,
                  int* dir, 
                  Side::LoHiSide* side, 
                  Real* a_values)
{
    ParmParse pp;
    Real bcVal = 0.0;
    if ( *dir == 0 ) {
       if (*side == Side::Lo) {
          bcVal = GlobalBCRS::s_xlo_neum;
       } else { 
          bcVal = GlobalBCRS::s_xhi_neum;
       }    
    } else if ( *dir == 1 ) {
       if (*side == Side::Lo) {
          bcVal = GlobalBCRS::s_ylo_neum;
       } else { 
          bcVal = GlobalBCRS::s_yhi_neum;
       }    
    }
    a_values[0] = bcVal;
}


Real DirichletValue(int dir, 
                    Side::LoHiSide side)
{
    Real bcVal = 0.0;
    if ( dir == 0 ) {
       if (side == Side::Lo) {
          bcVal = GlobalBCRS::s_xlo_diri;
       } else { 
          bcVal = GlobalBCRS::s_xhi_diri;
       }    
    } else if ( dir == 1 ) {
       if (side == Side::Lo) {
          bcVal = GlobalBCRS::s_ylo_diri;
       } else { 
          bcVal = GlobalBCRS::s_yhi_diri;
       }    
    }
    return bcVal;
}

/* Return Dirichlet boundary val -- based on parsed value */
void DirichletValue(Real* pos,
                    int* dir, 
                    Side::LoHiSide* side, 
                    Real* a_values)
{
    a_values[0] = DirichletValue(*dir, *side);
}


/* Return BC val -- based on parsed conditions and values */
void mixBCValues(FArrayBox& a_state,
                 const Box& a_valid,
                 const ProblemDomain& a_domain,
                 Real a_dx,
                 bool a_homogeneous)
{
  if(!a_domain.domainBox().contains(a_state.box())) {

      if (!GlobalBCRS::s_areBCsParsed) {
          ParseBC();
      }

      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir) {
          // don't do anything if periodic -- should be perio in y dir 1
          if (!a_domain.isPeriodic(dir)) {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              // box of ghost cells is outside of domain bounds ?
              if(!a_domain.domainBox().contains(ghostBoxLo)) {
                  if (GlobalBCRS::s_bcLo[dir] == 0) {
		              DiriBC(a_state,
		                     valid,
		                     a_dx,
		                     a_homogeneous,
		                     DirichletValue,
		                     dir,
		                     Side::Lo, 1);
                  } else if (GlobalBCRS::s_bcLo[dir] == 1) {
		              NeumBC(a_state,
		                     valid,
		                     a_dx,
		                     a_homogeneous,
		                     NeumannValue,
		                     dir,
		                     Side::Lo);
                  }
              }
              // box of ghost cells is outside of domain bounds ?
              if(!a_domain.domainBox().contains(ghostBoxHi)) {
                  if (GlobalBCRS::s_bcHi[dir] == 0) {
		              DiriBC(a_state,
		                     valid,
		                     a_dx,
		                     a_homogeneous,
		                     DirichletValue,
		                     dir,
		                     Side::Hi, 1);
                  } else if (GlobalBCRS::s_bcHi[dir] == 1) {
		              NeumBC(a_state,
		                     valid,
		                     a_dx,
		                     a_homogeneous,
		                     NeumannValue,
		                     dir,
                             Side::Hi);
                  }
              }
          } // end if is not periodic in ith direction
      } // end dir loop
  }
}


/* BC FOR ANY VAR */

/* For any var -- not just Head because based on adj val */
void FixedNeumBCFill(FArrayBox& a_state,
                     const Box& a_valid,
                     const ProblemDomain& a_domain,
                     Real a_dx)
{
  if(!a_domain.domainBox().contains(a_state.box())) {
      for(int dir=0; dir<CH_SPACEDIM; ++dir) {
          // don't do anything if periodic -- should be perio in y dir 1
          if (!a_domain.isPeriodic(dir)) {
              Box ghostBoxLo = adjCellBox(a_valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(a_valid, dir, Side::Hi, 1);

              if ((!a_domain.domainBox().contains(ghostBoxLo)) && (a_state.box().contains(ghostBoxLo)) ) {
                  // Neum
                  Box fromRegion = ghostBoxLo;
                  int isign = sign(Side::Lo);
                  fromRegion.shift(dir, -isign);
                  a_state.copy(a_state, fromRegion, 0, ghostBoxLo, 0, a_state.nComp());
              } // End ghostBoxLo in dir

              if ((!a_domain.domainBox().contains(ghostBoxHi)) && (a_state.box().contains(ghostBoxHi)) ) {
                  // Neum
                  Box fromRegion = ghostBoxHi;
                  int isign = sign(Side::Hi);
                  fromRegion.shift(dir, -isign);
                  a_state.copy(a_state, fromRegion, 0, ghostBoxHi, 0, a_state.nComp());
              } // End ghostBoxHi in dir

          } // end if is not periodic in ith direction
      } // end dir loop
  }
}


/* This function forces the ghost cells on the domain boundaries of a_state to be 0 -- no linear interpolation */ 
void NullBCFill(FArrayBox& a_state,
            const Box& a_valid,
            const ProblemDomain& a_domain,
            Real a_dx)
{
  // If box is outside of domain bounds ?
  if(!a_domain.domainBox().contains(a_state.box())) {
      for(int dir=0; dir<CH_SPACEDIM; ++dir) {
          // don't do anything if periodic -- should be perio in y dir 1
          if (!a_domain.isPeriodic(dir)) {

              Box ghostBoxLo = adjCellBox(a_valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(a_valid, dir, Side::Hi, 1);

              if ((!a_domain.domainBox().contains(ghostBoxLo)) && (a_state.box().contains(ghostBoxLo)) ) {
                  ghostBoxLo &= a_state.box();
                  for (BoxIterator bit(ghostBoxLo); bit.ok(); ++bit) {
                      IntVect ivTo = bit();
                      a_state(ivTo, 0) = 0.0;
                  }
              } // End ghostBoxLo in dir

              if ((!a_domain.domainBox().contains(ghostBoxHi)) && (a_state.box().contains(ghostBoxHi)) ) {
                  ghostBoxHi &= a_state.box();
                  for (BoxIterator bit(ghostBoxHi); bit.ok(); ++bit) {
                      IntVect ivTo = bit();
                      a_state(ivTo, 0) = 0.0;
                  }
              } // End ghostBoxHi in dir

          } // end if is not periodic in ith direction
      } // end dir loop
  }
}



/* AmrHydro class functions */

/* Solve for head with NL piece of operator provided as ext func */
void
AmrHydro::SolveForHead_nl(const Vector<DisjointBoxLayout>&               a_grids,
                          Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                          Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                          Vector<ProblemDomain>& a_domains,
                          Vector<int>& refRatio,
                          Real& coarsestDx,
                          Vector<LevelData<FArrayBox>*>& a_head, 
                          Vector<LevelData<FArrayBox>*>& a_RHS)
{
    Vector<RefCountedPtr<LevelData<FArrayBox> > > B(m_finest_level + 1);
    Vector<RefCountedPtr<LevelData<FArrayBox> > > Pri(m_finest_level + 1);
    Vector<RefCountedPtr<LevelData<FArrayBox> > > zb(m_finest_level + 1);

    IntVect nGhost  = m_num_head_ghost * IntVect::Unit;

    for (int lev = 0; lev <= m_finest_level; lev++) {
        B[lev]   = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, nGhost));
        Pri[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, nGhost));
        zb[lev]  = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, nGhost));

        LevelData<FArrayBox>& levelGapHeight    = *m_gapheight[lev];
        LevelData<FArrayBox>& levelPressi       = *m_overburdenpress[lev];
        LevelData<FArrayBox>& levelBedElevation = *m_bedelevation[lev];
        DataIterator dit = levelGapHeight.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            (*B[lev])[dit].copy(  levelGapHeight[dit], 0, 0, 1);
            (*Pri[lev])[dit].copy(levelPressi[dit], 0, 0, 1);
            (*zb[lev])[dit].copy( levelBedElevation[dit], 0, 0, 1);
        }
    }

    VCAMRNonLinearPoissonOpFactory poissonOpF_head;
    NL_level NLfunctTmp        = &AmrHydro::NonLinear_level;
    waterFlux_level wFfunctTmp = &AmrHydro::WFlx_level;
    poissonOpF_head.define(a_domains[0],
                           a_grids,
                           refRatio,
                           coarsestDx,
                           &mixBCValues,
                           0.0, a_aCoef,
                           - 1.0, a_bCoef,
                          this, NLfunctTmp, wFfunctTmp,
                          B, Pri, zb, m_compute_Bcoeff);

    AMRLevelOpFactory<LevelData<FArrayBox> >& opFactoryPtr = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) poissonOpF_head;

    AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;
    amrSolver = new AMRFASMultiGrid<LevelData<FArrayBox> >();

    // bottom solver ?
    BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
    if (m_verbosity > 3) {
        bottomSolver.m_verbosity = 4;
    } else {
        bottomSolver.m_verbosity = 1;
    }

    int numLevels = m_finest_level + 1;
    amrSolver->define(a_domains[0], opFactoryPtr,
                      &bottomSolver, numLevels);

    int numSmooth = 4;  // number of relax before averaging
    int numBottom = 8;  // num of bottom smoothings
    int numMG     = 1;  // Vcycle selected
    int maxIter   = 100; // max number of v cycles
    Real eps        =  1.0e-8;  // solution tolerance
    Real normThresh =  1.0e-10;  // abs tol
    Real hang       =  0.01;      // next rnorm should be < (1-m_hang)*norm_last 
    if (m_cur_step < 50) { 
        numBottom  = 10;  // num of bottom smoothings
        eps        =  1.0e-10;  // solution tolerance
        hang       =  0.0001;      // next rnorm should be < (1-m_hang)*norm_last 
    }
    amrSolver->setSolverParameters(numSmooth, numSmooth, numBottom,
                                   numMG, maxIter, 
                                   eps, hang, normThresh);
    if (m_cur_step < 50) { 
        amrSolver->m_imin = 20;
    }

    if (m_verbosity > 3) {
        amrSolver->m_verbosity = 4;
    } else {
        amrSolver->m_verbosity = 1;
    } 

    // solve !
    bool zeroInitialGuess = false;
    amrSolver->solve(a_head, a_RHS, m_finest_level, 0, zeroInitialGuess);

    delete amrSolver;
}

AmrHydro::AmrHydro() : m_IBCPtr(NULL)
{
    setParams();
    setDefaults();
}

//  Default for hydro_params
void
AmrHydro::setParams()
{
    m_suhmoParm = new suhmo_params;
    m_suhmoParm->readInputs();
}


void
AmrHydro::setDefaults()
{
    pout() << "AmrHydro::setDefaults()" << endl;

    // set some bogus values as defaults
    m_PrintCustom = false;
    m_post_proc   = false;
    m_is_defined = false;
    m_verbosity = 4;
    m_max_level = -1;
    m_finest_level = -1;
    m_block_factor = -1;
    m_max_box_size = 32;
    m_max_base_grid_size = 32;
    m_fill_ratio = -1;
    m_do_restart = false;
    m_restart_step = -1;

    m_domainSize = -1 * RealVect::Unit;

    // set the rest of these to reasonable defaults
    m_regrid_lbase = 0;
    m_nesting_radius = 1;
    m_tag_defined = false;
    m_n_tag_var = 0;
    m_tags_grow = 0;
    m_tags_grow_dir = IntVect::Zero;

    m_plot_prefix = "plot";
    m_plot_interval = 10000000;
    m_plot_time_interval = 1.0e+12;
    //m_write_gradPhi = false;

    m_check_prefix = "chk";
    m_check_interval = -1;
    m_check_overwrite = true;
    m_check_exit = false;

    m_fixed_dt = 0.0;
    m_num_head_ghost = 1;

    m_stable_dt = 0.0;

    m_offsetTime = 0.0;
}

AmrHydro::~AmrHydro()
{
    if (m_verbosity > 4)
    {
        pout() << "AmrHydro::~AmrHydro()" << endl;
    }

    // clean up memory
    for (int lev = 0; lev < m_head.size(); lev++)
    {
        if (m_head[lev] != NULL)
        {
            delete m_head[lev];
            m_head[lev] = NULL;
        }
        if (m_gradhead[lev] != NULL)
        {
            delete m_gradhead[lev];
            m_gradhead[lev] = NULL;
        }
        if (m_gradhead_ec[lev] != NULL)
        {
            delete m_gradhead_ec[lev];
            m_gradhead_ec[lev] = NULL;
        }
        if (m_old_head[lev] != NULL)
        {
            delete m_old_head[lev];
            m_old_head[lev] = NULL;
        }
        if (m_gapheight[lev] != NULL)
        {
            delete m_gapheight[lev];
            m_gapheight[lev] = NULL;
        }
        if (m_old_gapheight[lev] != NULL)
        {
            delete m_old_gapheight[lev];
            m_old_gapheight[lev] = NULL;
        }
        if (m_Pw[lev] != NULL)
        {
            delete m_Pw[lev];
            m_Pw[lev] = NULL;
        }
        if (m_qw[lev] != NULL)
        {
            delete m_qw[lev];
            m_qw[lev] = NULL;
        }
        if (m_Re[lev] != NULL)
        {
            delete m_Re[lev];
            m_Re[lev] = NULL;
        }
        if (m_meltRate[lev] != NULL)
        {
            delete m_meltRate[lev];
            m_meltRate[lev] = NULL;
        }
        if (m_iceheight[lev] != NULL)
        {
            delete m_iceheight[lev];
            m_iceheight[lev] = NULL;
        }
        if (m_bedelevation[lev] != NULL)
        {
            delete m_bedelevation[lev];
            m_bedelevation[lev] = NULL;
        }
        if (m_overburdenpress[lev] != NULL)
        {
            delete m_overburdenpress[lev];
            m_overburdenpress[lev] = NULL;
        }
        // Bed randomization
        if (m_bumpHeight[lev] != NULL)
        {
            delete m_bumpHeight[lev];
            m_bumpHeight[lev] = NULL;
        }
        if (m_bumpSpacing[lev] != NULL)
        {
            delete m_bumpSpacing[lev];
            m_bumpSpacing[lev] = NULL;
        }
    }

    delete m_suhmoParm;
    m_suhmoParm = NULL;

    if (m_IBCPtr != NULL)
    {
        delete m_IBCPtr;
        m_IBCPtr = NULL;
    }
    
    // that should be it!
}

void
AmrHydro::initialize()
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::initialize" << endl;
    }

    /* PARSING */
    ParmParse ppSolver("solver");
    m_use_FAS        = false;
    m_compute_Bcoeff = false;
    m_use_NL         = false;
    ppSolver.query("use_fas", m_use_FAS); // use FAS scheme for head 
    if (m_use_FAS) {
        ppSolver.query("use_NL", m_use_NL); // use FAS formulation with NL portion 
        if (m_use_NL) {
            ppSolver.query("bcoeff_otf", m_compute_Bcoeff); // compute B coeffs on the fly 
        }
    }
    m_eps_PicardIte = 1.0e-6;
    ppSolver.query("eps_PicardIte", m_eps_PicardIte);  

    m_time     = 0.0;  // start at t = 0
    m_cur_step = 0;    // start at dt = 0
    m_cur_PicardIte = 0;

    ParmParse ppAmr("AmrHydro");
    ppAmr.get("max_level", m_max_level);  // max level
    if (m_max_level > 0) {
        m_refinement_ratios.resize(m_max_level, -1);
        ppAmr.getarr("ref_ratios", m_refinement_ratios, 0, m_max_level);
    } else {
        m_refinement_ratios.resize(1);
        m_refinement_ratios[0] = -1;
    }
    ppAmr.query("PrintCustom", m_PrintCustom); // To plot after each Picard ite -- debug mode
    ppAmr.query("post_proc", m_post_proc); // To print some post-proc analysis-- debug mode
    ppAmr.query("verbosity", m_verbosity);
    ppAmr.query("block_factor", m_block_factor);
    ppAmr.query("regrid_lbase", m_regrid_lbase);  // smaller lev subject to regridding
    ppAmr.query("regrid_interval", m_regrid_interval);
    ppAmr.query("block_factor", m_block_factor);
    ppAmr.query("max_box_size", m_max_box_size);
    m_max_base_grid_size = m_max_box_size;
    ppAmr.query("max_base_grid_size", m_max_base_grid_size);
    Vector<int> ancells(SpaceDim);
    ppAmr.getarr("num_cells", ancells, 0, ancells.size());  // num cells in each dir
    // this one doesn't have a vertical dimension
    Vector<int> domLoIndex(SpaceDim, 0);
    ppAmr.queryarr("domainLoIndex", domLoIndex, 0, SpaceDim);  // allows for domains with lower indices which are not positive
    bool is_periodic[SpaceDim];
    for (int dir = 0; dir < SpaceDim; dir++) is_periodic[dir] = false; // assumption is that domains are not periodic
    Vector<int> is_periodic_int(SpaceDim, 0);
    ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; dir++) {
        is_periodic[dir] = (is_periodic_int[dir] == 1);
    }
    ppAmr.query("cfl", m_cfl);
    m_initial_cfl = m_cfl;  // cfl subject to change
    ppAmr.query("initial_cfl", m_initial_cfl);
    ppAmr.query("max_dt_grow_factor", m_max_dt_grow);   // not sure
    ppAmr.query("fixed_dt", m_fixed_dt);
    ppAmr.query("offsetTime", m_offsetTime);    // not sure
    ppAmr.query("plot_interval", m_plot_interval);
    ppAmr.query("plot_time_interval", m_plot_time_interval);
    ppAmr.query("plot_prefix", m_plot_prefix);
    //ppAmr.query("write_gradPhi", m_write_gradPhi);  // do not bother with outputing grad for now
    ppAmr.query("check_interval", m_check_interval);   // not sure
    ppAmr.query("check_prefix", m_check_prefix);
    ppAmr.query("check_overwrite", m_check_overwrite);
    ppAmr.query("check_exit", m_check_exit);   // no idea
    ppAmr.get("fill_ratio", m_fill_ratio);     // no idea
    ppAmr.query("nestingRadius", m_nesting_radius);
    ppAmr.query("tags_grow", m_tags_grow);
    {
        Vector<int> tgd(SpaceDim, 0);
        ppAmr.queryarr("tags_grow_dir", tgd, 0, tgd.size());
        for (int dir = 0; dir < SpaceDim; dir++) {
            m_tags_grow_dir[dir] = tgd[dir];
        }
    }
    ppAmr.query("n_tag_variables", m_n_tag_var); // how many variables are we using to tag ?
    if (m_n_tag_var != 0) {
        //
        m_tag_var.resize(m_n_tag_var); 
        ppAmr.getarr("tag_variables", m_tag_var, 0, m_n_tag_var);
        //
        m_tag_min.resize(m_n_tag_var,0.0); 
        ppAmr.getarr("tagging_mins", m_tag_min, 0, m_n_tag_var);
        //
        m_tag_cap.resize(m_n_tag_var,0.0); 
        ppAmr.getarr("tagging_caps", m_tag_cap, 0, m_n_tag_var);
        //
        m_tagging_val.resize(m_n_tag_var,0.0); 
        ppAmr.getarr("tagging_values", m_tagging_val, 0, m_n_tag_var);
        //
        m_tag_defined = true;
    } else {
        if (m_max_level > 0) {
            // abort if there isn't a tagging criterion
            MayDay::Error("No Tagging criterion defined");
        }
    }

    // Set up problem domains
    {
        IntVect loVect = IntVect(D_DECL(domLoIndex[0], domLoIndex[1], domLoIndex[3]));
        IntVect hiVect( D_DECL(domLoIndex[0] + ancells[0] - 1, 
                               domLoIndex[1] + ancells[1] - 1, 
                               domLoIndex[2] + ancells[2] - 1));
        ProblemDomain baseDomain(loVect, hiVect);
        // Periodicity
        for (int dir = 0; dir < SpaceDim; dir++) baseDomain.setPeriodic(dir, is_periodic[dir]);

        // Set up vector of domains
        m_amrDomains.resize(m_max_level + 1);
        m_amrDx.resize(m_max_level + 1);
        m_amrDomains[0] = baseDomain;
        m_amrDx[0] = m_domainSize[0] / baseDomain.domainBox().size(0) * RealVect::Unit;
        for (int lev = 1; lev <= m_max_level; lev++)
        {
            m_amrDomains[lev] = refine(m_amrDomains[lev - 1], m_refinement_ratios[lev - 1]);
            m_amrDx[lev] = m_amrDx[lev - 1] / m_refinement_ratios[lev - 1];
        }
    } // leaving problem domain setup scope

    // Tagging via an external file of boxes ? not sure I get this 
    std::string tagSubsetBoxesFile = "";
    m_vectTagSubset.resize(m_max_level);
    ppAmr.query("tagSubsetBoxesFile", tagSubsetBoxesFile);
    if (tagSubsetBoxesFile != "") {
        if (procID() == uniqueProc(SerialTask::compute)) {
            ifstream is(tagSubsetBoxesFile.c_str(), ios::in);
            int lineno = 1;
            if (is.fail()) {
                pout() << "Can't open " << tagSubsetBoxesFile << std::endl;
                MayDay::Error("Cannot open refine boxes file");
            }
            for (int lev = 0; lev < m_max_level; lev++) {
                // allowable tokens to identify levels in tag subset file
                const char level[6] = "level";
                const char domain[7] = "domain";
                char s[6];
                is >> s;
                if (std::string(level) == std::string(s)) {
                    int inlev;
                    is >> inlev;
                    if (inlev != lev)
                    {
                        pout() << "expected ' " << lev << "' at line " << lineno << std::endl;
                        MayDay::Error("bad input file");
                    }
                } else if (std::string(domain) == std::string(s)) {
                    // basic idea here is that we read in domain box
                    // (domains must be ordered from coarse->fine)
                    // until we get to a domain box which matches ours.
                    // This lets us make a single list of subset regions
                    // which we can use for any coarsening/refining of the domain
                    const Box& levelDomainBox = m_amrDomains[lev].domainBox();
                    bool stillLooking = true;
                    while (stillLooking) {
                        Box domainBox;
                        is >> domainBox;
                        if (domainBox == levelDomainBox) {
                            pout() << "Found a domain matching level " << lev << endl;
                            stillLooking = false;
                        } else { // move on until we find our level
                            // read in info for the level we're ignoring
                            // advance to next line
                            while (is.get() != '\n')
                                ;
                            lineno++;
                            int nboxes;
                            is >> nboxes;
                            if (nboxes > 0) {
                                for (int i = 0; i < nboxes; ++i) {
                                    Box box;
                                    is >> box;
                                    while (is.get() != '\n')
                                        ;
                                    lineno++;
                                }
                            }
                            is >> s;
                            if (std::string(domain) != std::string(s)) {
                                pout() << "expected '" << domain << "' at line " << lineno << ", got " << s
                                       << std::endl;
                                MayDay::Error("bad input file");
                            }
                        }
                    }
                } else {
                    pout() << "expected '" << level << "' or '" << domain << "' at line " << lineno << ", got " << s
                           << std::endl;
                    MayDay::Error("bad input file");
                }
                // advance to next line
                while (is.get() != '\n')
                    ;
                lineno++;
                int nboxes;
                is >> nboxes;
                if (nboxes > 0) {
                    for (int i = 0; i < nboxes; ++i) {
                        Box box;
                        is >> box;
                        while (is.get() != '\n')
                            ;
                        lineno++;
                        m_vectTagSubset[lev] |= box;
                        pout() << " level " << lev << " refine box : " << box << std::endl;
                    }
                }
                // advance to next line
                while (is.get() != '\n')
                    ;
                lineno++;

                if (lev > 0) {
                    // add lower level's subset to this subset
                    IntVectSet crseSet(m_vectTagSubset[lev - 1]);
                    if (!crseSet.isEmpty()) {
                        crseSet.refine(m_refinement_ratios[lev - 1]);
                        // crseSet.nestingRegion(m_block_factor,m_amrDomains[lev]);
                        if (m_vectTagSubset[lev].isEmpty()) {
                            m_vectTagSubset[lev] = crseSet;
                        } else {
                            m_vectTagSubset[lev] &= crseSet;
                        }
                    }
                }
            } // end loop over levels

        } // end if serial compute
        for (int lev = 0; lev < m_max_level; lev++) broadcast(m_vectTagSubset[lev], uniqueProc(SerialTask::compute));
    } // end if (tagSubsetBoxesFile != "")

    // check to see if we're using predefined grids
    bool usePredefinedGrids = false;
    std::string gridFile;
    if (ppAmr.contains("grids_file")) {
        usePredefinedGrids = true;
        ppAmr.get("grids_file", gridFile);
    }

    // check to see if we're restarting from a checkpoint file
    if (!ppAmr.contains("restart_file")) {
        // we are NOT restarting from a checkpoint file
        if (m_verbosity > 3) {
            pout() << "\nAmrHydro::initialize - Initializing data containers" << endl;
        }

        //---------------------------------
        // Set AMR hierarchy vector
        //---------------------------------
        // Time-dependent variables    
        m_old_head.resize(m_max_level + 1, NULL);
        m_head.resize(m_max_level + 1, NULL);

        m_gradhead.resize(m_max_level + 1, NULL);
        m_gradhead_ec.resize(m_max_level + 1, NULL);

        m_old_gapheight.resize(m_max_level + 1, NULL);
        m_gapheight.resize(m_max_level + 1, NULL);
        
        m_Pw.resize(m_max_level + 1, NULL);
        m_qw.resize(m_max_level + 1, NULL);
        m_Re.resize(m_max_level + 1, NULL);
        m_meltRate.resize(m_max_level + 1, NULL);
        m_iceheight.resize(m_max_level + 1, NULL);

        // Time constant variables
        m_bedelevation.resize(m_max_level + 1, NULL);
        m_overburdenpress.resize(m_max_level + 1, NULL);

        // Bed randomization
        m_bumpHeight.resize(m_max_level + 1, NULL);
        m_bumpSpacing.resize(m_max_level + 1, NULL);

        //-------------------------------------------------
        // For each level, define a collection of FArray/FluxBox
        //-------------------------------------------------
        for (int lev = 0; lev <= m_max_level; lev++) {
            m_old_head[lev] = new LevelData<FArrayBox>;
            m_head[lev] =  new LevelData<FArrayBox>;
            
            m_gradhead[lev] = new LevelData<FArrayBox>;
            m_gradhead_ec[lev] = new LevelData<FluxBox>;

            m_old_gapheight[lev] = new LevelData<FArrayBox>;
            m_gapheight[lev] = new LevelData<FArrayBox>;
           
            m_Pw[lev] = new LevelData<FArrayBox>;
            m_qw[lev] = new LevelData<FArrayBox>;
            m_Re[lev] = new LevelData<FArrayBox>;
            m_meltRate[lev] = new LevelData<FArrayBox>;
            m_iceheight[lev] = new LevelData<FArrayBox>;

            m_bedelevation[lev] = new LevelData<FArrayBox>;
            m_overburdenpress[lev] = new LevelData<FArrayBox>;

            m_bumpHeight[lev] = new LevelData<FArrayBox>;
            m_bumpSpacing[lev] = new LevelData<FArrayBox>;
        }


        int finest_level = -1;
        if (usePredefinedGrids) {
            setupFixedGrids(gridFile);
        } else {
            // now create  grids
            initGrids(finest_level);
        }

        /* Null the unused levels ... */
        for (int lev = m_finest_level+1; lev <= m_max_level; lev++) {
            m_old_head[lev]    = NULL;
            m_head[lev]        = NULL;
            
            m_gradhead[lev]    = NULL;
            m_gradhead_ec[lev] = NULL;

            m_old_gapheight[lev] = NULL;
            m_gapheight[lev]     = NULL;
           
            m_Pw[lev] = NULL;
            m_qw[lev] = NULL;
            m_Re[lev] = NULL;
            m_meltRate[lev]  = NULL;
            m_iceheight[lev] = NULL;

            m_bedelevation[lev]    = NULL;
            m_overburdenpress[lev] = NULL;

            m_bumpHeight[lev] = NULL;
            m_bumpSpacing[lev] = NULL;
        }

    } else {
        // we're restarting from a checkpoint file
        string restart_file;
        ppAmr.get("restart_file", restart_file);
        m_do_restart = true;
#ifdef CH_USE_HDF5
        restart(restart_file);
#endif
        // once we've set up everything, this lets us over-ride the
        // time and step number in the restart checkpoint file with
        // one specified in the inputs
        // AF dunno if we still need that
        if (ppAmr.contains("restart_time")) {
            Real restart_time;
            ppAmr.get("restart_time", restart_time);
            m_time = restart_time;
        }

        if (ppAmr.contains("restart_step")) {
            int restart_step;
            ppAmr.get("restart_step", restart_step);
            m_cur_step = restart_step;
            m_restart_step = restart_step;
        }
    }

    // set up counter of number of cells
    m_num_cells.resize(m_max_level + 1, 0);
    for (int lev = 0; lev <= m_finest_level; lev++) {
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        if (m_verbosity > 3) {
            // write out initial grids
            pout() << "    Level " << lev << " grids: " << levelGrids << endl;
        }
        LayoutIterator lit = levelGrids.layoutIterator();
        for (lit.begin(); lit.ok(); ++lit) {
            const Box& thisBox = levelGrids.get(lit());
            m_num_cells[lev] += thisBox.numPts();
        }
    }

    // finally, set up covered_level flags
    m_covered_level.resize(m_max_level + 1, 0);
    // note that finest level can't be covered.
    for (int lev = m_finest_level - 1; lev >= 0; lev--) {
        if (m_covered_level[lev + 1] == 1) {
            // if the next finer level is covered, then this one is too.
            m_covered_level[lev] = 1;
        } else {
            // see if the grids finer than this level completely cover it
            IntVectSet fineUncovered(m_amrDomains[lev + 1].domainBox());
            const DisjointBoxLayout& fineGrids = m_amrGrids[lev + 1];
            LayoutIterator lit = fineGrids.layoutIterator();
            for (lit.begin(); lit.ok(); ++lit) {
                const Box& thisBox = fineGrids.get(lit());
                fineUncovered.minus_box(thisBox);
            }
            if (fineUncovered.isEmpty()) {
                m_covered_level[lev] = 1;
            }
        }
    } // end loop over levels to determine covered levels

    if (m_verbosity > 3) {
        pout() << "Done with AmrHydro::initialize\n" << endl;
    }
}

// Init conditions
void
AmrHydro::setIBC(HydroIBC* a_IBC)
{
    m_IBCPtr = a_IBC->new_hydroIBC();
}

/* Main advance function */
void
AmrHydro::run(Real a_max_time, int a_max_step)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::run -- max_time= " << a_max_time << ", max_step = " << a_max_step << endl;
    }

    Real dt;
    // only call computeInitialDt if we're not doing restart
    if (!m_do_restart) {
        dt = computeInitialDt();
    } else {
        dt = computeDt();
    }

    // advance solution until done
    if (!(m_plot_time_interval > TIME_EPS) || m_plot_time_interval > a_max_time) m_plot_time_interval = a_max_time;

    while (a_max_time > m_time && (m_cur_step < a_max_step)) {
        Real next_plot_time = m_plot_time_interval * (1.0 + Real(int((m_time / m_plot_time_interval))));
        if (!(next_plot_time > m_time)) {
            // trap case where machine precision results in (effectively)
            // m_plot_time_interval * (1.0 + Real(int((m_time/m_plot_time_interval)))) == m_time
            next_plot_time += m_plot_time_interval;
        }

        next_plot_time = std::min(next_plot_time, a_max_time);

        while ((next_plot_time > m_time) && (m_cur_step < a_max_step) && (dt > TIME_EPS)) {
            if ((m_cur_step % m_plot_interval == 0) && m_plot_interval > 0) {
#ifdef CH_USE_HDF5
                writePlotFile();
#endif
            }

            if ( (m_cur_step != 0) && (m_cur_step != m_restart_step) && (m_cur_step % m_regrid_interval == 0)) {
                regrid();
            }

            if (m_cur_step != 0) {
                dt = computeDt();
            }

            if (next_plot_time - m_time + TIME_EPS < dt) dt = std::max(2 * TIME_EPS, next_plot_time - m_time);

            if ((m_cur_step % m_check_interval == 0) && (m_check_interval > 0) && (m_cur_step != m_restart_step)) {
#ifdef CH_USE_HDF5
                writeCheckpointFile();
#endif
                if (m_cur_step > 0 && m_check_exit) {
                    if (m_verbosity > 2) {
                        return;
                    }
                }
            }

            /* core */
            timeStepFAS(dt);               // we really want to use this

        } // end of plot_time_interval
#ifdef CH_USE_HDF5
        if (m_plot_interval >= 0) writePlotFile();
#endif
    } // end timestepping loop

    // dump out final plotfile, if appropriate
    if (m_plot_interval >= 0) {
#ifdef CH_USE_HDF5
        writePlotFile();
#endif
    }

    // dump out final checkpoint file, if appropriate
    if (m_check_interval >= 0) {
#ifdef CH_USE_HDF5
        writeCheckpointFile();
#endif
    }

    if (m_verbosity > 2)
    {
        pout() << "AmrHydro::run finished" << endl;
    }
}

/* Needed routines for timeStep */
void AmrHydro::WFlx_level(LevelData<FluxBox>&          a_bcoef, 
                          const LevelData<FArrayBox>&  a_u,
                          const LevelData<FArrayBox>*  a_ucoarse,
                          LevelData<FArrayBox>&        a_B,
                          Real                         a_dx,
                          bool                         a_print_WFX,
                          int                          a_smooth,
                          int                          a_depth)
{
    const DisjointBoxLayout& levelGrids  = a_u.disjointBoxLayout();
    const ProblemDomain&     levelDomain = levelGrids.physDomain();

    // Copy phi
    LevelData<FArrayBox> lcl_u(levelGrids, 1, a_u.ghostVect() ); 
    DataIterator dit = a_u.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        lcl_u[dit].copy(  a_u[dit], 0, 0, 1);
    }

    // Compute CC gradient of phi
    LevelData<FArrayBox> lvlgradH(levelGrids, SpaceDim, a_u.ghostVect() ); 
    LevelData<FArrayBox>* crsePsiPtr = NULL;
    LevelData<FArrayBox>* finePsiPtr = NULL;
    int nRefCrse=-1;
    int nRefFine=-1;
    // CC version -- GC ???
    Gradient::compGradientCC(lvlgradH, lcl_u,
                             crsePsiPtr, finePsiPtr,
                             a_dx, nRefCrse, nRefFine,
                             levelDomain);
    lvlgradH.exchange();
    ExtrapGhostCells( lvlgradH, levelDomain);

    // handle ghost cells on the coarse-fine interface
    if ( a_ucoarse != NULL) {
        const DisjointBoxLayout& levelGridscoarse  = a_ucoarse->getBoxes();
        const ProblemDomain&     levelDomaincoarse = levelGridscoarse.physDomain();

        // Copy phi coarse
        LevelData<FArrayBox> lcl_ucoarse(levelGridscoarse, 1, a_ucoarse->ghostVect() ); 
        DataIterator dit2 = (*a_ucoarse).dataIterator();
        for (dit2.begin(); dit2.ok(); ++dit2) {
            lcl_ucoarse[dit2].copy( (*a_ucoarse)[dit2], 0, 0, 1);
        }

        // Compute CC gradient of phi coarse
        Real a_dxcoarse = a_dx*2;
        LevelData<FArrayBox> lvlgradHcoarse(levelGridscoarse, SpaceDim, a_ucoarse->ghostVect() ); 
        Gradient::compGradientCC(lvlgradHcoarse, lcl_ucoarse,
                                 crsePsiPtr, finePsiPtr,
                                 a_dxcoarse, nRefCrse, nRefFine,
                                 levelDomaincoarse);
        lvlgradHcoarse.exchange();
        ExtrapGhostCells( lvlgradHcoarse, levelDomaincoarse);

        QuadCFInterp qcfi(levelGrids, &levelGridscoarse,
                          a_dx, 2.,  
                          2,  // ncomps
                          levelDomain);
        qcfi.coarseFineInterp(lvlgradH, lvlgradHcoarse);
    }

    // Compute Re
    LevelData<FArrayBox> lvlRe(levelGrids, 1, a_u.ghostVect() ); 
    for (dit.begin(); dit.ok(); ++dit) {

        const Box& region = lvlRe[dit].box();

        FORT_COMPUTERE( CHF_FRA(a_B[dit]),
                        CHF_FRA(lvlgradH[dit]),
                        CHF_BOX(region),
                        CHF_FRA(lvlRe[dit]),
                        CHF_CONST_REAL(m_suhmoParm->m_omega),
                        CHF_CONST_REAL(m_suhmoParm->m_nu) );
    }

    //if (a_print_WFX) {
    //    /* custom plt here -- debug print */
    //        int nStuffToPlot = 4;
    //        Vector<std::string> vectName;
    //        vectName.resize(nStuffToPlot);
    //        vectName[0]="GapHeight";
    //        vectName[1]="Head";
    //        vectName[2]="Reynolds";
    //        vectName[3]="GradH";

    //        Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
    //        stuffToPlot.resize(nStuffToPlot);
    //        for (int zz = 0; zz < nStuffToPlot; zz++) {
    //            stuffToPlot[zz].resize(1, NULL);
    //        }

    //        stuffToPlot[0][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());
    //        stuffToPlot[1][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());
    //        stuffToPlot[2][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());
    //        stuffToPlot[3][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());

    //        LevelData<FArrayBox>& levelGapSTP   = *stuffToPlot[0][0];
    //        LevelData<FArrayBox>& levelHeadSTP  = *stuffToPlot[1][0];
    //        LevelData<FArrayBox>& levelReSTP    = *stuffToPlot[2][0];
    //        LevelData<FArrayBox>& levelGradHSTP = *stuffToPlot[3][0];

    //        DataIterator dit = levelHeadSTP.dataIterator();
    //        for (dit.begin(); dit.ok(); ++dit) {
    //            levelGapSTP[dit].copy(a_B[dit], 0, 0, 1);
    //            levelHeadSTP[dit].copy(lcl_u[dit], 0, 0, 1);
    //            levelReSTP[dit].copy(lvlRe[dit], 0, 0, 1);
    //            levelGradHSTP[dit].copy(lvlgradH[dit], 0, 0, 1);
    //        }
    //        writePltWFX(nStuffToPlot, vectName, stuffToPlot, 
    //                    levelDomain, levelGrids, 
    //                    a_smooth, a_depth);
    //}
    
    // CC -> EC
    LevelData<FluxBox>   lvlB_ec(levelGrids, 1, IntVect::Zero);
    LevelData<FluxBox>   lvlRe_ec(levelGrids, 1, IntVect::Zero);
    CellToEdge(lvlRe, lvlRe_ec);
    CellToEdge(a_B, lvlB_ec);

    // Update Bcoef
    dit = a_bcoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FluxBox& bC        = a_bcoef[dit];
        FluxBox& B_ec      = lvlB_ec[dit];
        FluxBox& Re_ec     = lvlRe_ec[dit];
        // loop over directions
        for (int dir = 0; dir<SpaceDim; dir++) {

            const Box& region = bC[dir].box();

            FORT_COMPUTEBCOEFF( CHF_FRA(B_ec[dir]),
                                CHF_FRA(Re_ec[dir]),
                                CHF_BOX(region),
                                CHF_FRA(bC[dir]),
                                CHF_CONST_REAL(m_suhmoParm->m_omega),
                                CHF_CONST_REAL(m_suhmoParm->m_nu) );
        }
    }

}


void AmrHydro::NonLinear_level(LevelData<FArrayBox>&        a_NL, 
                               LevelData<FArrayBox>&        a_dNL,
                               const LevelData<FArrayBox>&  a_u,
                               LevelData<FArrayBox>&        a_B,
                               LevelData<FArrayBox>&        a_Pi,
                               LevelData<FArrayBox>&        a_zb)
{

  if (m_use_NL) {

      DataIterator levelDit = a_NL.dataIterator();
      const DisjointBoxLayout& dbl = a_NL.disjointBoxLayout();

      for (levelDit.begin(); levelDit.ok(); ++levelDit) {

          const Box& region = dbl[levelDit];

          FORT_COMPUTENONLINEARTERMS( CHF_FRA(a_u[levelDit]),
                                      CHF_FRA(a_B[levelDit]),
                                      CHF_FRA(a_Pi[levelDit]),
                                      CHF_FRA(a_zb[levelDit]),
                                      CHF_BOX(region),
                                      CHF_FRA(a_NL[levelDit]),
                                      CHF_FRA(a_dNL[levelDit]),
                                      CHF_CONST_REAL(m_suhmoParm->m_A),
                                      CHF_CONST_REAL(m_suhmoParm->m_cutOffbr),
                                      CHF_CONST_REAL(m_suhmoParm->m_maxOffbr));

          //    // Correction Colin 03/18
          //    if ( m_suhmoParm->m_br > thisB(iv,0) ) {
          //        //thisNL(iv, 0) = thisNL(iv, 0)  * std::tanh( std::pow(thisB(iv,0),3) / std::pow(m_suhmoParm->m_br,3) );
          //        //thisdNL(iv, 0) = thisdNL(iv, 0) * std::tanh( std::pow(thisB(iv,0),3) / std::pow(m_suhmoParm->m_br,3) );
          //    }
      } // end loop over grids on this level
  }
}


void
AmrHydro::compute_grad_zb_ec(int                 lev,
                             LevelData<FluxBox>& a_levelgradZb_ec)
{
    // Compute grad(zb) -EC only - 
    LevelData<FArrayBox>& levelZb       = *m_bedelevation[lev];
    LevelData<FArrayBox>* crsePsiPtr = NULL;
    LevelData<FArrayBox>* finePsiPtr = NULL;
    int nRefCrse=-1;
    int nRefFine=-1;
    if (lev > 0) {
        crsePsiPtr = m_bedelevation[lev-1];
        nRefCrse = m_refinement_ratios[lev-1];
    }
    if (lev < m_finest_level) {
        finePsiPtr = m_bedelevation[lev+1];  // What does it do with the fine stuff ???
        nRefFine = m_refinement_ratios[lev];
    }
    Real dx = m_amrDx[lev][0];  
    // EC version
    Gradient::compGradientMAC(a_levelgradZb_ec, levelZb,
                              crsePsiPtr, finePsiPtr,
                              dx, nRefCrse, nRefFine,
                              m_amrDomains[lev]);
}

void
AmrHydro::compute_grad_head(int lev)
{
    // Compute grad(h) -EC and CC- 
    LevelData<FArrayBox>& levelcurrentH = *m_head[lev];
    LevelData<FArrayBox>& levelgradH    = *m_gradhead[lev];
    LevelData<FluxBox>&   levelgradH_ec = *m_gradhead_ec[lev];

    LevelData<FArrayBox>* crsePsiPtr = NULL;
    LevelData<FArrayBox>* finePsiPtr = NULL;
    int nRefCrse=-1;
    int nRefFine=-1;
    if (lev > 0) {
        crsePsiPtr = m_head[lev-1];
        nRefCrse = m_refinement_ratios[lev-1];
    }
    if (lev < m_finest_level) {
        finePsiPtr = m_head[lev+1];  // What does it do with the fine stuff ???
        nRefFine = m_refinement_ratios[lev];
    }
    Real dx = m_amrDx[lev][0];  
    // CC version
    Gradient::compGradientCC(levelgradH, levelcurrentH,
                             crsePsiPtr, finePsiPtr,
                             dx, nRefCrse, nRefFine,
                             m_amrDomains[lev]);
    // handle ghost cells on the coarse-fine interface
    if (lev > 0) {
        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                          m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                          2,  // num comps
                          m_amrDomains[lev]);
        qcfi.coarseFineInterp(*m_gradhead[lev], *m_gradhead[lev-1]);
    }
    // Need to fill the ghost cells of gradH -- extrapolate on the no perio boundaries   
    levelgradH.exchange();
    ExtrapGhostCells( levelgradH, m_amrDomains[lev]);

    // EC version
    Gradient::compGradientMAC(levelgradH_ec, levelcurrentH,
                             crsePsiPtr, finePsiPtr,
                             dx, nRefCrse, nRefFine,
                             m_amrDomains[lev]);
}


void
AmrHydro::compute_grad_var(int lev,
                           Vector<LevelData<FArrayBox>*>&  a_var,
                           Vector<LevelData<FArrayBox>*>&  a_gradvar,
                           Vector<LevelData<FluxBox>*>&    a_gradvar_ec)
{
    // Compute grad(h) -EC and CC- 
    LevelData<FArrayBox>& levelvar        = *a_var[lev];
    LevelData<FArrayBox>& levelgradvar    = *a_gradvar[lev];
    LevelData<FluxBox>&   levelgradvar_ec = *a_gradvar_ec[lev];

    LevelData<FArrayBox>* crsePsiPtr = NULL;
    LevelData<FArrayBox>* finePsiPtr = NULL;
    int nRefCrse=-1;
    int nRefFine=-1;
    if (lev > 0) {
        crsePsiPtr = a_var[lev-1];
        nRefCrse = m_refinement_ratios[lev-1];
    }
    if (lev < m_finest_level) {
        finePsiPtr = a_var[lev+1];  // What does it do with the fine stuff ???
        nRefFine = m_refinement_ratios[lev];
    }
    Real dx = m_amrDx[lev][0];  
    // CC version
    Gradient::compGradientCC(levelgradvar, levelvar,
                             crsePsiPtr, finePsiPtr,
                             dx, nRefCrse, nRefFine,
                             m_amrDomains[lev]);
    // handle ghost cells on the coarse-fine interface
    if (lev > 0) {
        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                          m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                          2,  // num comps
                          m_amrDomains[lev]);
        qcfi.coarseFineInterp(*a_gradvar[lev], *a_gradvar[lev-1]);
    }
    // Need to fill the ghost cells of gradH -- extrapolate on the no perio boundaries   
    levelgradvar.exchange();
    ExtrapGhostCells( levelgradvar, m_amrDomains[lev]);

    // EC version
    Gradient::compGradientMAC(levelgradvar_ec, levelvar,
                             crsePsiPtr, finePsiPtr,
                             dx, nRefCrse, nRefFine,
                             m_amrDomains[lev]);
}

void 
AmrHydro::evaluate_Qw_ec(int lev, 
                         LevelData<FluxBox>& levelQw_ec, 
                         LevelData<FluxBox>& levelB_ec,
                         LevelData<FluxBox>& levelRe_ec)
{
    // For Qw
    LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];

    DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
    DataIterator dit                    = levelGrids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        // EC quantities
        FluxBox& Qwater_ec = levelQw_ec[dit];
        FluxBox& currB_ec  = levelB_ec[dit];
        FluxBox& Re_ec     = levelRe_ec[dit];
        FluxBox& gradH_ec  = levelgradH_ec[dit];

        // loop over directions
        for (int dir = 0; dir<SpaceDim; dir++) {

            const Box& region = Qwater_ec[dir].box();

            FORT_COMPUTEQW( CHF_FRA(currB_ec[dir]),
                            CHF_FRA(Re_ec[dir]),
                            CHF_FRA(gradH_ec[dir]),
                            CHF_BOX(region),
                            CHF_FRA(Qwater_ec[dir]),
                            CHF_CONST_REAL(m_suhmoParm->m_omega),
                            CHF_CONST_REAL(m_suhmoParm->m_nu) );
        } 
    }
}

void
AmrHydro::evaluate_Re_quadratic(int lev, bool computeGrad)
{
    // For Re
    LevelData<FArrayBox>& levelRe       = *m_Re[lev]; 
    LevelData<FArrayBox>& levelgradH    = *m_gradhead[lev];
    LevelData<FArrayBox>& levelB        = *m_gapheight[lev];    

    if (computeGrad) {
        // Compute grad(h) -EC and CC- 
        LevelData<FArrayBox>& levelcurrentH = *m_head[lev];
        LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];

        LevelData<FArrayBox>* crsePsiPtr = NULL;
        LevelData<FArrayBox>* finePsiPtr = NULL;
        int nRefCrse=-1;
        int nRefFine=-1;
        if (lev > 0) {
            crsePsiPtr = m_head[lev-1];
            nRefCrse = m_refinement_ratios[lev-1];
        }
        if (lev < m_finest_level) {
            finePsiPtr = m_head[lev+1];  // What does it do with the fine stuff ???
            nRefFine = m_refinement_ratios[lev];
        }
        Real dx = m_amrDx[lev][0];  
        // CC version
        Gradient::compGradientCC(levelgradH, levelcurrentH,
                                 crsePsiPtr, finePsiPtr,
                                 dx, nRefCrse, nRefFine,
                                 m_amrDomains[lev]);
        // handle ghost cells on the coarse-fine interface
        if (lev > 0) {
            QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                              m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                              2,  // num comps
                              m_amrDomains[lev]);
            qcfi.coarseFineInterp(*m_gradhead[lev], *m_gradhead[lev-1]);
        }
        // Need to fill the ghost cells of gradH -- extrapolate on the no perio boundaries   
        levelgradH.exchange();
        ExtrapGhostCells( levelgradH, m_amrDomains[lev]);

        // EC version
        Gradient::compGradientMAC(levelgradH_ec, levelcurrentH,
                                 crsePsiPtr, finePsiPtr,
                                 dx, nRefCrse, nRefFine,
                                 m_amrDomains[lev]);
    }

    // Compute Re at CC
    DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
    DataIterator dit                    = levelGrids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {

        const Box& region = levelRe[dit].box();

        FORT_COMPUTERE( CHF_FRA(levelB[dit]),
                            CHF_FRA(levelgradH[dit]),
                            CHF_BOX(region),
                            CHF_FRA(levelRe[dit]),
                            CHF_CONST_REAL(m_suhmoParm->m_omega),
                            CHF_CONST_REAL(m_suhmoParm->m_nu) );
    }
}


void
AmrHydro::aCoeff_bCoeff(LevelData<FArrayBox>&  levelacoef, 
                        LevelData<FluxBox>&    levelbcoef, 
                        LevelData<FluxBox>&    levelRe, 
                        LevelData<FluxBox>&    levelB)
{
    DataIterator dit = levelbcoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {

        FluxBox& B        = levelB[dit];
        FArrayBox& aC     = levelacoef[dit];
        FluxBox& bC       = levelbcoef[dit];
        FluxBox& Re       = levelRe[dit];

        aC.setVal(0.0);

        // loop over directions
        for (int dir = 0; dir<SpaceDim; dir++) {

            const Box& region = bC[dir].box();

            FORT_COMPUTEBCOEFF( CHF_FRA(B[dir]),
                                CHF_FRA(Re[dir]),
                                CHF_BOX(region),
                                CHF_FRA(bC[dir]),
                                CHF_CONST_REAL(m_suhmoParm->m_omega),
                                CHF_CONST_REAL(m_suhmoParm->m_nu) );
        }
    }
}


void
AmrHydro::dCoeff(LevelData<FluxBox>&    leveldcoef, 
                 Real                   a_dx,
                 Real                   a_dt)
{
    DataIterator dit = leveldcoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {

        FluxBox& bC       = leveldcoef[dit];

        // loop over directions
        for (int dir = 0; dir<SpaceDim; dir++) {

            const Box& region = bC[dir].box();

            FORT_COMPUTEDCOEFF( CHF_BOX(region),
                                CHF_FRA(bC[dir]),
                                CHF_CONST_REAL(a_dx),
                                CHF_CONST_REAL(a_dt));

        }
    }
}


void 
AmrHydro::Calc_moulin_source_term_distributed (LevelData<FArrayBox>& levelMoulinSrc) 
{
   // UNITS FOR m_moulin_flux SHOULD BE M3/S
   // WOrk on coarser level
   int curr_level = 0;

   // keep track of each moulin contrib separately
   LevelData<FArrayBox> levelMoulinSrcTmp(m_amrGrids[curr_level], m_suhmoParm->m_n_moulins, IntVect::Zero);
   std::vector<Real> a_moulinsInteg;    // m.s-1 Sliding velocity
   a_moulinsInteg.resize(m_suhmoParm->m_n_moulins, 0.0);

   // Loop over points and moulins to get influx at each location
   DataIterator dit = levelMoulinSrcTmp.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
       FArrayBox& moulinSrcTmp   = levelMoulinSrcTmp[dit];

       moulinSrcTmp.setVal(0.0);

       BoxIterator bit(moulinSrcTmp.box());
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           // Position of centroid
           Real x_loc = (iv[0]+0.5)*m_amrDx[curr_level][0];
           Real y_loc = (iv[1]+0.5)*m_amrDx[curr_level][1];
           // Loop over all moulins
           for (int m = 0; m<m_suhmoParm->m_n_moulins; m++) {
               // Radius
               //Real rad_loc = (x_loc - m_suhmoParm->m_moulin_position[m*2 + 0])**2 + (y_loc - m_suhmoParm->m_moulin_position[m*2 + 1])**2;
               moulinSrcTmp(iv, m) = 1.0 /( m_suhmoParm->m_sigma[m] * std::sqrt(2.0 * 3.14) )  * 
                                                 std::exp(- 1.0 / (2.0 * m_suhmoParm->m_sigma[m] * m_suhmoParm->m_sigma[m]) * ( 
                                                 std::pow(x_loc - m_suhmoParm->m_moulin_position[m*2 + 0], 2)
                                              +  std::pow(y_loc - m_suhmoParm->m_moulin_position[m*2 + 1], 2) ) );
               if (moulinSrcTmp(iv, m) < 1.0e-09) {
                   moulinSrcTmp(iv, m) = 0.0;
               }
               a_moulinsInteg[m] += moulinSrcTmp(iv,m) * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]; 
           }
       }
   }

#ifdef CH_MPI
   Real recv;
   for (int m = 0; m<m_suhmoParm->m_n_moulins; m++) {
       int result = MPI_Allreduce(&a_moulinsInteg[m], &recv, 1, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       a_moulinsInteg[m] = recv;
   }
#endif
 
   // Rescale each moulins
   std::vector<Real> a_moulinsIntegFinal;    // m.s-1 Sliding velocity
   a_moulinsIntegFinal.resize(m_suhmoParm->m_n_moulins, 0.0);
   for (dit.begin(); dit.ok(); ++dit) {
       FArrayBox& moulinSrc      = levelMoulinSrc[dit];
       FArrayBox& moulinSrcTmp   = levelMoulinSrcTmp[dit];

       moulinSrc.setVal(m_suhmoParm->m_distributed_input);

       BoxIterator bit(moulinSrcTmp.box());
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           for (int m = 0; m<m_suhmoParm->m_n_moulins; m++) {
               moulinSrc(iv,0) += moulinSrcTmp(iv,m) * ( 1.0 - m_suhmoParm->m_runoff * std::sin(2.0 * Pi * m_time / 86400. ) )  
                                  / a_moulinsInteg[m] *  m_suhmoParm->m_moulin_flux[m] ; 
               a_moulinsIntegFinal[m] += moulinSrcTmp(iv,m) * ( 1.0 - m_suhmoParm->m_runoff * std::sin(2.0 * Pi * m_time / 86400. ) )   
                                  / a_moulinsInteg[m] * m_suhmoParm->m_moulin_flux[m] * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]; 
           }
       }
   }

#ifdef CH_MPI
   for (int m = 0; m<m_suhmoParm->m_n_moulins; m++) {
       recv = 0;
       int result = MPI_Allreduce(&a_moulinsIntegFinal[m], &recv, 1, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       a_moulinsIntegFinal[m] = recv;
   }
#endif
   // Check
   Real totalRechar = 0.0;
   for (int m = 0; m<m_suhmoParm->m_n_moulins; m++) {
       totalRechar += a_moulinsIntegFinal[m];
   }
   if (m_verbosity > 3) {
       pout() << "level is: " << curr_level << ", Integral of all moulins is (vol. rate): " << totalRechar << endl; 
   } 

}


void
AmrHydro::CalcRHS_gapHeightFAS(LevelData<FArrayBox>& levelRHS_b, 
                               LevelData<FArrayBox>& levelPi, 
                               LevelData<FArrayBox>& levelPw, 
                               LevelData<FArrayBox>& levelmR, 
                               LevelData<FArrayBox>& levelB,
                               LevelData<FArrayBox>& levelDT,
                               LevelData<FArrayBox>& levelbumpHeight,
                               LevelData<FArrayBox>& levelbumpSpacing,
                               LevelData<FArrayBox>& levelRHS_b_A, 
                               LevelData<FArrayBox>& levelRHS_b_B,
                               LevelData<FArrayBox>& levelRHS_b_C,
                               LevelData<FArrayBox>& levelchanDegree)
{
   DataIterator dit = levelRHS_b.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

       FArrayBox& B       = levelB[dit];
       FArrayBox& RHS     = levelRHS_b[dit];

       FArrayBox& DT      = levelDT[dit];

       FArrayBox& RHS_A   = levelRHS_b_A[dit];
       FArrayBox& RHS_B   = levelRHS_b_B[dit];
       FArrayBox& RHS_C   = levelRHS_b_C[dit];
       FArrayBox& CD      = levelchanDegree[dit];

       FArrayBox& Pressi  = levelPi[dit];
       FArrayBox& Pw      = levelPw[dit];
       FArrayBox& meltR   = levelmR[dit];
       FArrayBox& BH      = levelbumpHeight[dit];
       FArrayBox& BL      = levelbumpSpacing[dit];
       
       // initialize RHS for h
       RHS.setVal(0.0);
       RHS_A.setVal(0.0);
       RHS_B.setVal(0.0);
       RHS_C.setVal(0.0);

       // first term
       RHS.copy(meltR, 0, 0, 1);
       RHS_A.copy(meltR, 0, 0, 1);
       RHS *= 1.0 / m_suhmoParm->m_rho_i;
       RHS_A *= 1.0 / m_suhmoParm->m_rho_i;
       // third term ...
       Real ub_norm = std::sqrt(  m_suhmoParm->m_ub[0]*m_suhmoParm->m_ub[0] 
                                + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]); 
       BoxIterator bit(RHS.box()); // can use gridBox? 
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           if ( B(iv,0) < BH(iv,0)) {
               RHS(iv,0) +=  ub_norm * (BH(iv,0) - B(iv,0)) / BL(iv,0);
               RHS_B(iv,0) = ub_norm * (BH(iv,0) - B(iv,0)) / BL(iv,0);
           }
           // second term ... assume  n = 3 !!
           Real PimPw = (Pressi(iv,0) - Pw(iv,0));
           Real AbsPimPw = std::abs(PimPw);
           if ( m_suhmoParm->m_cutOffbr > B(iv,0) ) {
               RHS(iv,0)   -= m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0) * ( 1.0 - (m_suhmoParm->m_cutOffbr - B(iv,0)) / m_suhmoParm->m_cutOffbr );
               RHS_C(iv,0) =- m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0) * ( 1.0 - (m_suhmoParm->m_cutOffbr - B(iv,0)) / m_suhmoParm->m_cutOffbr );

               // Add a Diffusive term to mdot
               RHS(iv,0)   += m_suhmoParm->m_DiffFactor * DT(iv,0) * B(iv,0) * ( 1.0 - (m_suhmoParm->m_cutOffbr - B(iv,0)) / m_suhmoParm->m_cutOffbr );

           } else if ( m_suhmoParm->m_maxOffbr < B(iv,0) ) {
               RHS(iv,0)   -= m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0) * ( 1.0 - (m_suhmoParm->m_maxOffbr  - B(iv,0)) / m_suhmoParm->m_maxOffbr  );
               RHS_C(iv,0) =- m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0) * ( 1.0 - (m_suhmoParm->m_maxOffbr  - B(iv,0)) / m_suhmoParm->m_maxOffbr  );

               // Add a Diffusive term to mdot
               RHS(iv,0)   += m_suhmoParm->m_DiffFactor * DT(iv,0) * B(iv,0) * ( 1.0 - (m_suhmoParm->m_maxOffbr  - B(iv,0)) / m_suhmoParm->m_maxOffbr  );

           } else {
               RHS(iv,0)   -= m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0); 
               RHS_C(iv,0) =- m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0);

               // Add a Diffusive term to mdot
               RHS(iv,0)   += m_suhmoParm->m_DiffFactor * DT(iv,0) * B(iv,0);
           }
           CD(iv,0) = RHS_A(iv,0) / (RHS_A(iv,0) + RHS_B(iv,0));
       }
   }
}


void
AmrHydro::timeStepFAS(Real a_dt)
{
    CH_TIME("AmrHydro::timestepFAS");

    m_cur_step += 1;

    if (m_verbosity > 2) {
        pout() << "\n\n-- Timestep " << m_cur_step << " Advancing solution from time " << m_time << " ( " << time()
               << ")"
                  " with dt = "
               << a_dt << endl;
    }
    
    /* Sketch of FAS */
    //
    // NOTE: all computation done at CC. When edge qties needed, use CC->Edge
    // I Copy new h and new b into old h and old b
    //     Take care of GC/BC etc.
    // II h calc
    //     Fill perio GC and BC of h and b
    //     Put h into h_lag (GC too) and b into b_lag
    //     Interp CC b to EC
    //     Compute grad(h) and grad(zb)
    //     IV Re/Qw dependency !! TWO options 
    //             One: inner loop with fixed nb of iter
    //             Second: Solve for quad equation to get Re
    //     Update  RHS = f(Qw, grad(zb))
    //     Evaluate mR 
    //     Compute lagged adv term = f((Qw, grad(h))
    //     Compute alpha, aCoeff (0) and beta, bCoeff
    //     Solve for h with FAS scheme
    // III b calc
    //     Form RHS for b
    //     Solve for b using Forward Euler simple scheme
    //  
    /* End comments */

    IntVect HeadGhostVect = m_num_head_ghost * IntVect::Unit;
    Real coarsestDx = m_amrDx[0][0];  

    /* I Copy new h and new b into old h and old b */

    // Also create and initialize tmp vectors
    Vector<LevelData<FArrayBox>*> a_head_lagged;
    Vector<LevelData<FArrayBox>*> a_head_curr;
    Vector<LevelData<FArrayBox>*> a_gapheight_lagged;
    Vector<LevelData<FArrayBox>*> RHS_h;
    Vector<LevelData<FArrayBox>*> moulin_source_term;
    Vector<LevelData<FArrayBox>*> RHS_b;
    Vector<LevelData<FArrayBox>*> a_ReQwIter;
    Vector<LevelData<FArrayBox>*> a_diffusiveTerm;
    Vector<LevelData<FArrayBox>*> a_chanDegree;
    // DEBUG 
    Vector<LevelData<FArrayBox>*> RHS_b_A;
    Vector<LevelData<FArrayBox>*> RHS_b_B;
    Vector<LevelData<FArrayBox>*> RHS_b_C;

    a_head_lagged.resize(m_finest_level + 1, NULL);
    a_head_curr.resize(m_finest_level + 1, NULL);
    a_gapheight_lagged.resize(m_finest_level + 1, NULL);
    RHS_h.resize(m_finest_level + 1, NULL);
    moulin_source_term.resize(m_finest_level + 1, NULL);
    RHS_b.resize(m_finest_level + 1, NULL);
    a_ReQwIter.resize(m_finest_level + 1, NULL);
    a_diffusiveTerm.resize(m_finest_level + 1, NULL);
    a_chanDegree.resize(m_finest_level + 1, NULL);
    RHS_b_A.resize(m_finest_level + 1, NULL);
    RHS_b_B.resize(m_finest_level + 1, NULL);
    RHS_b_C.resize(m_finest_level + 1, NULL);
    // FACE CENTERED STUFF
    Vector<LevelData<FluxBox>*> a_Qw_ec;
    Vector<LevelData<FluxBox>*> a_Re_ec;
    Vector<LevelData<FluxBox>*> a_GapHeight_ec;
    Vector<LevelData<FluxBox>*> a_gradZb_ec;
    Vector<LevelData<FluxBox>*> a_Dcoef;
    a_Qw_ec.resize(m_finest_level + 1, NULL);
    a_Re_ec.resize(m_finest_level + 1, NULL);
    a_GapHeight_ec.resize(m_finest_level + 1, NULL);
    a_gradZb_ec.resize(m_finest_level + 1, NULL);
    a_Dcoef.resize(m_finest_level + 1, NULL);
    // alpha*aCoef(x)*I - beta*Div(bCoef(x)*Grad) -- note for us: alpha = 0 beta = - 1 
    Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(m_finest_level + 1);
    Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(m_finest_level + 1);

    //     Take care of GC/BC etc.
    if (m_verbosity > 3) {
        pout() <<"   ...Copy current into old & take care of ghost cells and BCs "<< endl;
    }

        /* custom plt here -- debug print */
        //if (m_PrintCustom && (m_cur_step > 100045)) {
        //    int nStuffToPlot = 6;
        //    Vector<std::string> vectName;
        //    vectName.resize(nStuffToPlot);
        //    // H/B/Pice/Zb/Re/Hice/BH/BL are stuff that are read in case of restart
        //    vectName[0]="head";
        //    vectName[1]="gapHeight";
        //    vectName[2]="Pwater";
        //    vectName[3]="metingRate";
        //    vectName[4]="Reynolds";
        //    vectName[5]="Zbed";

        //    Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
        //    stuffToPlot.resize(nStuffToPlot);
        //    for (int zz = 0; zz < nStuffToPlot; zz++) {
        //        stuffToPlot[zz].resize(m_max_level + 1, NULL);
        //    }

        //    for (int lev = 0; lev <= m_finest_level; lev++) {
        //        stuffToPlot[0][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        //        stuffToPlot[1][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        //        stuffToPlot[2][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        //        stuffToPlot[3][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        //        stuffToPlot[4][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        //        stuffToPlot[5][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);

        //        LevelData<FArrayBox>& levelHead      = *m_head[lev];
        //        LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];

        //        LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
        //        LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];

        //        LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
        //        LevelData<FArrayBox>& levelPwSTP     = *stuffToPlot[2][lev];

        //        LevelData<FArrayBox>& levelMR        = *m_meltRate[lev];
        //        LevelData<FArrayBox>& levelMRSTP     = *stuffToPlot[3][lev];

        //        LevelData<FArrayBox>& levelRe        = *m_Re[lev];
        //        LevelData<FArrayBox>& levelReSTP     = *stuffToPlot[4][lev];

        //        LevelData<FArrayBox>& levelZb        = *m_bedelevation[lev];
        //        LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[5][lev];

        //        DataIterator dit = levelHead.dataIterator();
        //        for (dit.begin(); dit.ok(); ++dit) {
        //            levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
        //            levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
        //            levelPwSTP[dit].copy(levelPw[dit], 0, 0, 1);
        //            levelMRSTP[dit].copy(levelMR[dit], 0, 0, 1);
        //            levelReSTP[dit].copy(levelRe[dit], 0, 0, 1);
        //            levelZbSTP[dit].copy(levelZb[dit], 0, 0, 1);
        //        }
        //    } // loop on levs
        //    writePltCustom(nStuffToPlot, vectName, stuffToPlot, std::to_string(-10));
        //} // end customPlt


    for (int lev = 0; lev <= m_finest_level; lev++) {
        LevelData<FArrayBox>& oldH       = *m_old_head[lev];
        LevelData<FArrayBox>& currentH   = *m_head[lev];

        LevelData<FArrayBox>& oldB       = *m_old_gapheight[lev];
        LevelData<FArrayBox>& currentB   = *m_gapheight[lev];

        LevelData<FArrayBox>& currentRe  = *m_Re[lev];

        // handle ghost cells on the coarse-fine interface
        // Linear interp better
        if (lev > 0) {
            int nGhost = currentRe.ghostVect()[0];
            PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                m_amrGrids[lev - 1],
                                                1, // ncomps
                                                m_amrDomains[lev - 1],
                                                m_refinement_ratios[lev - 1],
                                                nGhost);
            headFiller.fillInterp(*m_head[lev], *m_head[lev-1], *m_head[lev-1], 0.0, 0, 0, 1);
            headFiller.fillInterp(*m_gapheight[lev], *m_gapheight[lev-1], *m_gapheight[lev-1], 0.0, 0, 0, 1);
            headFiller.fillInterp(*m_Re[lev], *m_Re[lev-1], *m_Re[lev-1], 0.0, 0, 0, 1);
        }

        // fill perio boundaries
        currentH.exchange();
        currentB.exchange();
        currentRe.exchange();


        // Head and b RHS
        RHS_h[lev]                = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        moulin_source_term[lev]   = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        RHS_b[lev]                = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_diffusiveTerm[lev]      = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_chanDegree[lev]         = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // Head and B lagged for iterations
        a_head_lagged[lev]      = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        a_head_curr[lev]        = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        a_gapheight_lagged[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // Re/Qw iterations -- testing
        a_ReQwIter[lev]         = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // DEBUG
        RHS_b_A[lev]            = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        RHS_b_B[lev]            = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        RHS_b_C[lev]            = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        // Stuff for OpLin
        aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero));
        bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero));
        // Face centered stuff
        a_Qw_ec[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_Re_ec[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_GapHeight_ec[lev]  = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_gradZb_ec[lev]     = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_Dcoef[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);

        // Get the valid boxes
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        DataIterator dit                    = currentH.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            // get the validBox
            const Box& validBox = levelGrids.get(dit);

            // Fill BC ghost cells of h and b
            mixBCValues(currentH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0], false);
            FixedNeumBCFill(currentB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

            // Copy curr into old -- copy ghost cells too 
            oldH[dit].copy(currentH[dit], 0, 0, 1);
            oldB[dit].copy(currentB[dit], 0, 0, 1);
        }
        ExtrapGhostCells( currentRe, m_amrDomains[lev]);

    } // there. We should start with consistent b and h, GC BC and all ...


    /* II h calc */

    if (m_verbosity > 3) {
        pout() <<"   ...Solve for h ! "<< endl;
    }
    bool converged_h = false;
    int ite_idx = 0;
    m_cur_PicardIte = 0;
    while (!converged_h) { 
        // Solve for h using lagged (iteration lagged) qtities
        //         Fill perio GC and BC of h and b
        //         Put h into h_lag (GC too) and b into b_lag
        //         Interp CC b to EC
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelcurH      = *m_head[lev];
            LevelData<FArrayBox>& levelnewH_lag  = *a_head_lagged[lev];

            LevelData<FArrayBox>& levelcurB      = *m_gapheight[lev];
            LevelData<FArrayBox>& levelnewB_lag  = *a_gapheight_lagged[lev];
            LevelData<FluxBox>&   levelnewB_ec   = *a_GapHeight_ec[lev];

            LevelData<FArrayBox>& levelcurRe     = *m_Re[lev];

            // Get the valid boxes
            const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

            // deal with c-f interf if any
            if (lev > 0) {
                int nGhost = levelcurH.ghostVect()[0];
                PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                    m_amrGrids[lev - 1],
                                                    1, // ncomps
                                                    m_amrDomains[lev - 1],
                                                    m_refinement_ratios[lev - 1],
                                                    nGhost);
                headFiller.fillInterp(*m_head[lev], *m_head[lev-1], *m_head[lev-1], 0.0, 0, 0, 1);
                headFiller.fillInterp(*m_gapheight[lev], *m_gapheight[lev-1], *m_gapheight[lev-1], 0.0, 0, 0, 1);
                headFiller.fillInterp(*m_Re[lev], *m_Re[lev-1], *m_Re[lev-1], 0.0, 0, 0, 1);
            }

            // fill perio boundaries
            levelcurH.exchange();
            levelcurB.exchange();
            levelcurRe.exchange();

            // Fill BC and put h into h_lag
            DataIterator dit = levelnewH_lag.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                // get the validBox & fill BC ghost cells
                const Box& validBox = levelGrids.get(dit);
                mixBCValues(levelcurH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0], false);
                FixedNeumBCFill(levelcurB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

                levelnewH_lag[dit].copy(levelcurH[dit], 0, 0, 1); // should copy ghost cells too !
                levelnewB_lag[dit].copy(levelcurB[dit], 0, 0, 1); // should copy ghost cells too !
            }
            ExtrapGhostCells( levelcurRe, m_amrDomains[lev]);

            // Interpolate b to edges
            CellToEdge(levelcurB, levelnewB_ec);
        }  // loop on levs -- same thing, we should start with consistent b and h GC/BC and all


        //         Compute grad(h) and grad(zb)
        if (m_verbosity > 3) {
            pout() <<"        Compute grad(h) and grad(Zb) "<< endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            // Compute grad(h) -EC and CC- 
            compute_grad_head(lev);

            // Compute grad(Zb) -EC-
            LevelData<FluxBox>& levelgradZb_ec  = *a_gradZb_ec[lev];
            compute_grad_zb_ec(lev, levelgradZb_ec);

            // Compute dCoeff 
            LevelData<FluxBox>&   levelDcoef  = *a_Dcoef[lev]; 
            Real dx = m_amrDx[lev][0];  
            dCoeff(levelDcoef, dx, a_dt);
        } // end loop levels


        //     IV Re/Qw dependency !! TWO options 
        bool m_Re_Q_loop = false;
        if (m_verbosity > 3) {
            pout() <<"        Re/Qw dependency "<< endl;
        }
        //             One: inner loop with fixed nb of iter
        //               x Update VECTOR Qw = f(Re, grad(h))
        //               x Update Re = f(Qw)
        if (m_Re_Q_loop) {
            int max_ite_Re = 50; 
            Real max_Re_diff = 0.0;
            for (int it = 0; it <= max_ite_Re; it++) {
                if (m_verbosity > 5) {
                    pout() << "           ------------------------" << endl;
                    pout() << "           ite " << it << endl;
                }
            
                max_Re_diff = -10000;

                for (int lev = 0; lev <= m_finest_level; lev++) {

                    LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 
                    LevelData<FArrayBox>& levelRe       = *m_Re[lev]; 
               
                    // EC quantities
                    LevelData<FluxBox>& levelB_ec       = *a_GapHeight_ec[lev];    
                    LevelData<FluxBox>& levelQw_ec      = *a_Qw_ec[lev]; 
                    LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];
                    LevelData<FluxBox>& levelRe_ec      = *a_Re_ec[lev];

                    // tmp holder
                    LevelData<FArrayBox>& levelReQwIter = *a_ReQwIter[lev];

                    DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
                    DataIterator dit                    = levelGrids.dataIterator();

                    CellToEdge(levelRe, levelRe_ec);

                    // Get Qw at EC and interp at CC
                    for (dit.begin(); dit.ok(); ++dit) {
                        // EC quantities
                        FluxBox& currB_ec  = levelB_ec[dit];
                        FluxBox& Qwater_ec = levelQw_ec[dit];
                        FluxBox& gradH_ec  = levelgradH_ec[dit];
                        FluxBox& Re_ec     = levelRe_ec[dit];

                        // loop over directions
                        for (int dir = 0; dir<SpaceDim; dir++) {
                            FArrayBox& Qwater_ecFab = Qwater_ec[dir];
                            FArrayBox& gradH_ecFab  = gradH_ec[dir];
                            FArrayBox& currB_ecFab  = currB_ec[dir];
                            FArrayBox& Re_ecFab     = Re_ec[dir];

                            BoxIterator bitEC(Qwater_ecFab.box()); // can use gridBox? 

                            for (bitEC.begin(); bitEC.ok(); ++bitEC) {
                                IntVect iv = bitEC();
                                Real num_q = - std::pow(currB_ecFab(iv, 0),3) * m_suhmoParm->m_gravity * gradH_ecFab(iv, 0);
                                Real denom_q = 12.0 * m_suhmoParm->m_nu * (1 + m_suhmoParm->m_omega * Re_ecFab(iv, 0));
                                Qwater_ecFab(iv, 0) = num_q/denom_q;
                            }
                        } 
                    }
                    EdgeToCell(levelQw_ec, levelQw); 
                    // handle ghost cells on the coarse-fine interface
                    if (lev > 0) {
                        int nGhost = levelQw.ghostVect()[0];
                        PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                            m_amrGrids[lev - 1],
                                                            2, // ncomps
                                                            m_amrDomains[lev - 1],
                                                            m_refinement_ratios[lev - 1],
                                                            nGhost);
                        headFiller.fillInterp(*m_qw[lev], *m_qw[lev-1], *m_qw[lev-1], 0.0, 0, 0, 1);
                    }
                    // Need to fill the ghost cells -- extrapolate on perio boundaries   
                    levelQw.exchange();
                    ExtrapGhostCells( levelQw, m_amrDomains[lev]);

                    //Get Re at CC
                    for (dit.begin(); dit.ok(); ++dit) {
                        FArrayBox& Qwater  = levelQw[dit];
                        FArrayBox& Re      = levelRe[dit];

                        levelReQwIter[dit].copy(Re, 0, 0, 1);

                        BoxIterator bit(Qwater.box()); // can use gridBox? 
                        for (bit.begin(); bit.ok(); ++bit) {
                            IntVect iv = bit();
                            // Update Re using this new Qw
                            Re(iv, 0) = std::sqrt( Qwater(iv, 0) * Qwater(iv, 0) 
                                                 + Qwater(iv, 1) * Qwater(iv, 1)) / m_suhmoParm->m_nu;
                        }

                        // Get max of diff bet prev ite Re and freshly computed one
                        levelReQwIter[dit].minus(Re, 0, 0, 1);
                        levelReQwIter[dit].abs();
                        max_Re_diff = std::max(max_Re_diff, levelReQwIter[dit].max());
                    }
                    // handle ghost cells on the coarse-fine interface
                    if (lev > 0) {
                        int nGhost = levelRe.ghostVect()[0];
                        PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                            m_amrGrids[lev - 1],
                                                            1, // ncomps
                                                            m_amrDomains[lev - 1],
                                                            m_refinement_ratios[lev - 1],
                                                            nGhost);
                        headFiller.fillInterp(*m_Re[lev], *m_Re[lev-1], *m_Re[lev-1], 0.0, 0, 0, 1);
                    }
                    // Need to fill the ghost cells of Re -- extrapolate on perio boundaries   
                    levelRe.exchange();
                    ExtrapGhostCells( levelRe, m_amrDomains[lev]);

                    CellToEdge(levelRe, levelRe_ec);
                } // end loop on levs

                if (m_verbosity > 5) {
                    pout() << "           ------------------------" << endl;
                }
            } // end Qw/Re ites
        //             Second: Solve for quad equation to get Re
        //               x Re = f(grad(h)) at CC
        //               x Re CC->EC
        //               x Update EC VECTOR Qw = f(Re, grad(h))
        //               x Qw EC->CC
        } else {
            // Solve quadratic pb. Re first, q second
            for (int lev = 0; lev <= m_finest_level; lev++) {
                // For Re
                // CC quantities
                LevelData<FArrayBox>& levelRe       = *m_Re[lev]; 
                // EC quantities
                LevelData<FluxBox>& levelRe_ec      = *a_Re_ec[lev];
                evaluate_Re_quadratic(lev, false);
                // handle ghost cells on the coarse-fine interface
                if (lev > 0) {
                    int nGhost = levelRe.ghostVect()[0];
                    PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                        m_amrGrids[lev - 1],
                                                        1, // ncomps
                                                        m_amrDomains[lev - 1],
                                                        m_refinement_ratios[lev - 1],
                                                        nGhost);
                    headFiller.fillInterp(*m_Re[lev], *m_Re[lev-1], *m_Re[lev-1], 0.0, 0, 0, 1);
                }
                // Need to fill the ghost cells of Re -- extrapolate on perio boundaries   
                levelRe.exchange();
                ExtrapGhostCells( levelRe, m_amrDomains[lev]);
                CellToEdge(levelRe, levelRe_ec);

                // For Qw
                // CC quantities
                LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 
                // EC quantities
                LevelData<FluxBox>& levelQw_ec      = *a_Qw_ec[lev]; 
                LevelData<FluxBox>& levelB_ec       = *a_GapHeight_ec[lev];    
                evaluate_Qw_ec(lev, levelQw_ec, levelB_ec, levelRe_ec);
                EdgeToCell(levelQw_ec, levelQw); 
                // handle ghost cells on the coarse-fine interface
                if (lev > 0) {
                    int nGhost = levelQw.ghostVect()[0];
                    PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                        m_amrGrids[lev - 1],
                                                        2, // ncomps
                                                        m_amrDomains[lev - 1],
                                                        m_refinement_ratios[lev - 1],
                                                        nGhost);
                    headFiller.fillInterp(*m_qw[lev], *m_qw[lev-1], *m_qw[lev-1], 0.0, 0, 0, 1);
                }
                // Need to fill the ghost cells -- extrapolate on boundaries   
                levelQw.exchange();
                ExtrapGhostCells( levelQw, m_amrDomains[lev]);
            } // end loop on levs
        }
        // Just in case the quad resolution gives funky results
        Real minRe = computeMin(m_Re, m_refinement_ratios, Interval(0,0), 0);
        if (minRe < 0.0) {
            pout() <<"        Re min is NEG !! ("<< minRe << ") abort ... "<< endl;
            MayDay::Error("Abort");
        }


        //     Update  RHS = f(Qw, grad(zb)) 
        //     Compute lagged adv term = f((Qw, grad(h))
        if (m_verbosity > 3) {
            pout() <<"        Update RHS(h) and adv term "<< endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            // CC quantities
            LevelData<FArrayBox>& levelB                  = *m_gapheight[lev];    
            LevelData<FArrayBox>& levelRHSh               = *RHS_h[lev];
            LevelData<FArrayBox>& levelmoulin_source_term = *moulin_source_term[lev];
            LevelData<FArrayBox>& levelPi                 = *m_overburdenpress[lev];
            LevelData<FArrayBox>& levelPw                 = *m_Pw[lev];
            LevelData<FArrayBox>& levelZb                 = *m_bedelevation[lev];
            LevelData<FArrayBox>& levelH                  = *m_head[lev];
            LevelData<FArrayBox>& levelBH                 = *m_bumpHeight[lev];
            LevelData<FArrayBox>& levelBL                 = *m_bumpSpacing[lev];

            LevelData<FArrayBox>& levelDterm              = *a_diffusiveTerm[lev];    

            // EC quantities
            LevelData<FluxBox>&   levelQw_ec     = *a_Qw_ec[lev]; 
            LevelData<FluxBox>&   levelgradH_ec  = *m_gradhead_ec[lev];
            LevelData<FluxBox>&   levelgradZb_ec = *a_gradZb_ec[lev];

            LevelData<FluxBox>&   levelDcoef     = *a_Dcoef[lev]; 

            DisjointBoxLayout& levelGrids        = m_amrGrids[lev];
            DataIterator dit                     = levelGrids.dataIterator();

            // tmp arrays for vect operations
            LevelData<FluxBox>   leveltmp_ec(levelGrids, 1, IntVect::Zero);
            LevelData<FArrayBox> leveltmp_cc(levelGrids, 1*SpaceDim, HeadGhostVect);
            LevelData<FluxBox>   leveltmp2_ec(levelGrids, 1, IntVect::Zero);
            LevelData<FArrayBox> leveltmp2_cc(levelGrids, 1*SpaceDim, HeadGhostVect);

           
            // UNITS FOR m_moulin_flux SHOULD BE M3/S
            // UNITS FOR m_distributed_input SHOULD BE M/S
            if (m_suhmoParm->m_n_moulins > 0) {
                if (lev == 0) {
                    Calc_moulin_source_term_distributed(levelmoulin_source_term);
                } else {
                    FineInterp interpolator(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1], m_amrDomains[lev]);
                    interpolator.interpToFine(*moulin_source_term[lev], *moulin_source_term[lev - 1]);
                }
            } else if (m_suhmoParm->m_n_moulins < 0) {
                if (m_suhmoParm->m_time_varying_input) {

                    Real T_K   = -16.0 * std::cos( 2.0 * Pi * m_time / ( 365.*24*60*60.) ) - 5.0 + m_suhmoParm->m_deltaT;

                    LevelData<FArrayBox>& levelIceHeight = *m_iceheight[lev];
                    for (dit.begin(); dit.ok(); ++dit) {

                        const Box& region = levelmoulin_source_term[dit].box();

                        FORT_COMPUTE_TIMEVARYINGRECHARGE( CHF_FRA(levelIceHeight[dit]),
                                                          CHF_BOX(region),
                                                          CHF_FRA(levelmoulin_source_term[dit]),
                                                          CHF_CONST_REAL(T_K),
                                                          CHF_CONST_REAL(m_suhmoParm->m_distributed_input) );
                    }
                } else {
                    for (dit.begin(); dit.ok(); ++dit) {
                        FArrayBox& moulinSrc = levelmoulin_source_term[dit];
                        moulinSrc.setVal(m_suhmoParm->m_distributed_input);
                    }
                }
            } else {
                for (dit.begin(); dit.ok(); ++dit) {
                    FArrayBox& moulinSrc = levelmoulin_source_term[dit];
                    moulinSrc.setVal(0.0);
                }
            }

            for (dit.begin(); dit.ok(); ++dit) {
                // EC quantities
                FluxBox& Qwater_ec = levelQw_ec[dit];
                FluxBox& gradH_ec  = levelgradH_ec[dit];
                FluxBox& tmp_ec    = leveltmp_ec[dit];
                FluxBox& gradZb_ec = levelgradZb_ec[dit];
                FluxBox& tmp2_ec   = leveltmp2_ec[dit];

                // loop over directions
                for (int dir = 0; dir<SpaceDim; dir++) {

                    const Box& region = tmp_ec[dir].box();

                    FORT_COMPUTESCAPROD( CHF_FRA(Qwater_ec[dir]),
                                         CHF_FRA(gradH_ec[dir]),
                                         CHF_FRA(gradZb_ec[dir]),
                                         CHF_BOX(region),
                                         CHF_FRA(tmp_ec[dir]),
                                         CHF_FRA(tmp2_ec[dir]) );

                } // loop over dir
            }
            // Qw gradH
            EdgeToCell(leveltmp_ec,  leveltmp_cc);
            // Qw gradZb
            EdgeToCell(leveltmp2_ec, leveltmp2_cc);

            // Compute diffusive term  here ?
            for (dit.begin(); dit.ok(); ++dit) {

                const Box& region = levelDterm[dit].box();
                const FluxBox& thisDcoef  = levelDcoef[dit];

                FORT_COMPUTEDIFTERM2D( CHF_FRA(levelB[dit]),
                                       CHF_BOX(region),
                                       CHF_FRA(levelDterm[dit]),
                                       CHF_CONST_FRA(thisDcoef[0]),
                                       CHF_CONST_FRA(thisDcoef[1]));
            }

            Real rho_coef = (1.0 /  m_suhmoParm->m_rho_w - 1.0 / m_suhmoParm->m_rho_i);
            Real ub_norm = std::sqrt(  m_suhmoParm->m_ub[0]*m_suhmoParm->m_ub[0] 
                           + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]); // / m_suhmoParm->m_lr;
            for (dit.begin(); dit.ok(); ++dit) {
                // CC
                FArrayBox& B       = levelB[dit];
                FArrayBox& RHSh    = levelRHSh[dit];
                FArrayBox& tmp_cc  = leveltmp_cc[dit];
                FArrayBox& tmp2_cc = leveltmp2_cc[dit];
                FArrayBox& Pressi  = levelPi[dit];
                FArrayBox& Pressw  = levelPw[dit];
                FArrayBox& zb      = levelZb[dit];
                FArrayBox& currH   = levelH[dit];
                FArrayBox& bumpHeight    = levelBH[dit];
                FArrayBox& bumpSpacing   = levelBL[dit];
                FArrayBox& DiffusiveTerm = levelDterm[dit];

                FArrayBox& moulinSrc = levelmoulin_source_term[dit];

                BoxIterator bit(RHSh.box());
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit();
                    // Actual RHS
                    if ( B(iv,0) < bumpHeight(iv,0)) {
                        RHSh(iv,0) = rho_coef / m_suhmoParm->m_L * (  m_suhmoParm->m_G 
                                     + m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                                      (tmp2_cc(iv, 0) + tmp2_cc(iv, 1)) )  
                                     - ub_norm * (bumpHeight(iv,0) - B(iv,0)) / bumpSpacing(iv,0);
                    } else {
                        RHSh(iv,0) = rho_coef / m_suhmoParm->m_L * (  m_suhmoParm->m_G 
                                     + m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                                      (tmp2_cc(iv, 0) + tmp2_cc(iv, 1)) ); 
                    }

                    // Adv term right now in RHS
                    RHSh(iv,0) -= rho_coef / m_suhmoParm->m_L * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * (1.0 
                                  + m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w) * 
                                  (tmp_cc(iv, 0) + tmp_cc(iv, 1));

                    // Add sliding term 
                    // mR -- turn off fric heat and geot heat
                    Pressw(iv,0) = m_suhmoParm->m_gravity * m_suhmoParm->m_rho_w * (currH(iv,0) - zb(iv,0));
                    Real sca_prod = 20. * 20. * m_suhmoParm->m_ub[0] * std::abs(Pressi(iv,0) - Pressw(iv,0)) * m_suhmoParm->m_ub[0];
                    RHSh(iv,0) += sca_prod * rho_coef / m_suhmoParm->m_L ;

                    // Add moulin 
                    RHSh(iv,0) += moulinSrc(iv,0);

                    if ( m_suhmoParm->m_cutOffbr > B(iv,0) ) {
                        // Add a Diffusive term to mdot
                        RHSh(iv,0)   += m_suhmoParm->m_DiffFactor * DiffusiveTerm(iv,0) * B(iv,0) * ( 1.0 - (m_suhmoParm->m_cutOffbr - B(iv,0)) / m_suhmoParm->m_cutOffbr );
                    } else if ( m_suhmoParm->m_maxOffbr < B(iv,0) ) {
                        // Add a Diffusive term to mdot
                        RHSh(iv,0)   += m_suhmoParm->m_DiffFactor * DiffusiveTerm(iv,0) * B(iv,0) * ( 1.0 - (m_suhmoParm->m_maxOffbr  - B(iv,0)) / m_suhmoParm->m_maxOffbr  );
                    } else {
                        // Add a Diffusive term to mdot
                        RHSh(iv,0) += m_suhmoParm->m_DiffFactor * DiffusiveTerm(iv,0) * B(iv,0);
                    }
                }
            }
        }// loop on levs


        //     Compute alpha, aCoeff (0) and beta, bCoeff
        Vector<DisjointBoxLayout> m_amrGrids_curr;
        Vector<ProblemDomain> m_amrDomains_curr;
        m_amrGrids_curr.resize(m_finest_level + 1);
        m_amrDomains_curr.resize(m_finest_level + 1);
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelacoef = *aCoef[lev];
            //EC quantities
            LevelData<FluxBox>& levelbcoef  = *bCoef[lev];
            LevelData<FluxBox>& levelRe_ec  = *a_Re_ec[lev]; 
            LevelData<FluxBox>& levelB_ec   = *a_GapHeight_ec[lev];    

            // Compute aCoeff and bCoeff using updated qtites
            aCoeff_bCoeff(levelacoef, levelbcoef, levelRe_ec, levelB_ec);

            //stuff for solve -- to get SAME number of levs and not max/finest
            m_amrGrids_curr[lev]   = m_amrGrids[lev];
            m_amrDomains_curr[lev] = m_amrDomains[lev];
            LevelData<FArrayBox>& levelcurH      = *m_head[lev];
            LevelData<FArrayBox>& levelcurHlcl   = *a_head_curr[lev];
            DisjointBoxLayout& levelGrids        = m_amrGrids[lev];
            DataIterator dit                     = levelGrids.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                // Copy curr into old -- copy ghost cells too 
                levelcurHlcl[dit].copy(levelcurH[dit], 0, 0, 1);
            }
        } // loop on levs

        //     Solve for h with FAS scheme
        if (m_verbosity > 3) {
            pout() <<"        Poisson solve for h "<< endl;
        }

        SolveForHead_nl(m_amrGrids_curr, aCoef, bCoef,
                        m_amrDomains_curr, m_refinement_ratios, coarsestDx,
                        a_head_curr, RHS_h);

        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelcurH      = *m_head[lev];
            LevelData<FArrayBox>& levelcurHlcl   = *a_head_curr[lev];
            DisjointBoxLayout& levelGrids        = m_amrGrids[lev];
            DataIterator dit                     = levelGrids.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                // Copy curr into old -- copy ghost cells too 
                levelcurH[dit].copy(levelcurHlcl[dit], 0, 0, 1);
            }
        }

        /* Averaging down and fill in ghost cells */
        if (m_verbosity > 3) {
            pout() <<"   ...Average down "<< endl;
        }
        for (int lev = m_finest_level; lev > 0; lev--) {
            CoarseAverage averager(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1]);
            averager.averageToCoarse(*m_head[lev - 1], *m_head[lev]);
        }
        // handle ghost cells on the coarse-fine interface
        for (int lev = 1; lev <= m_finest_level; lev++) {
            QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                              m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                              1,  // num comps
                              m_amrDomains[lev]);
            qcfi.coarseFineInterp(*m_head[lev], *m_head[lev-1]);
        }



        /* CONVERGENCE TESTS with LAGGED quantities */
        Real maxHead = computeMax(m_head, m_refinement_ratios, Interval(0,0), 0);
        Real minHead = computeMin(m_head, m_refinement_ratios, Interval(0,0), 0);
        if (m_verbosity > 3) {
            pout() <<"        Check for convergence of h (max = " << maxHead<<", min = " << minHead << ")"<<endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelnewH_lag  = *a_head_lagged[lev];
            LevelData<FArrayBox>& levelcurrentH  = *m_head[lev];

            DataIterator dit = levelnewH_lag.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                levelnewH_lag[dit].minus(levelcurrentH[dit], 0, 0, 1);
                levelnewH_lag[dit].divide(maxHead);
                levelnewH_lag[dit].abs();
            }
        }
        Real max_resH = computeMax(a_head_lagged, m_refinement_ratios, Interval(0,0), 0);
        if (m_verbosity > 3) {
            pout() <<ite_idx<< "         x(h) "<<max_resH<<endl;
        }

        /* custom plt here -- debug print */
        //if (m_PrintCustom && (m_cur_step > 45)) {
        //    int nStuffToPlot = 12;
        //    Vector<std::string> vectName;
        //    vectName.resize(nStuffToPlot);
        //    vectName[0]="head";
        //    vectName[1]="gapHeight";
        //    vectName[2]="head_residual";
        //    vectName[3]="RHS_head";
        //    vectName[4]="RHS_b";
        //    vectName[5]="Qw_X";
        //    vectName[6]="Qw_Y";
        //    vectName[7]="sourceMoulin";
        //    vectName[8]="Pwater";
        //    vectName[9]="metingRate";
        //    vectName[10]="Reynolds";
        //    vectName[11]="Zbed";

        //    Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
        //    stuffToPlot.resize(nStuffToPlot);
        //    for (int zz = 0; zz < nStuffToPlot; zz++) {
        //        stuffToPlot[zz].resize(m_max_level + 1, NULL);
        //    }

        //    for (int lev = 0; lev <= m_finest_level; lev++) {
        //        stuffToPlot[0][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[1][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[2][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[3][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[4][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[5][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[6][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[7][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[8][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[9][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[10][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        //        stuffToPlot[11][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);

        //        LevelData<FArrayBox>& levelHead      = *m_head[lev];
        //        LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];

        //        LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
        //        LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];

        //        LevelData<FArrayBox>& levelRes       = *a_head_lagged[lev];
        //        LevelData<FArrayBox>& levelResSTP    = *stuffToPlot[2][lev];

        //        LevelData<FArrayBox>& levelRHS       = *RHS_h[lev];
        //        LevelData<FArrayBox>& levelRHSSTP    = *stuffToPlot[3][lev];

        //        LevelData<FArrayBox>& levelRHSB       = *RHS_b[lev];
        //        LevelData<FArrayBox>& levelRHSBSTP    = *stuffToPlot[4][lev];

        //        LevelData<FArrayBox>& levelQw         = *m_qw[lev];
        //        LevelData<FArrayBox>& levelQwXSTP     = *stuffToPlot[5][lev];
        //        LevelData<FArrayBox>& levelQwYSTP     = *stuffToPlot[6][lev];

        //        LevelData<FArrayBox>& levelSourceM    = *moulin_source_term[lev];
        //        LevelData<FArrayBox>& levelSourceMSTP = *stuffToPlot[7][lev];

        //        LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
        //        LevelData<FArrayBox>& levelPwSTP     = *stuffToPlot[8][lev];

        //        LevelData<FArrayBox>& levelMR        = *m_meltRate[lev];
        //        LevelData<FArrayBox>& levelMRSTP     = *stuffToPlot[9][lev];

        //        LevelData<FArrayBox>& levelRe        = *m_Re[lev];
        //        LevelData<FArrayBox>& levelReSTP     = *stuffToPlot[10][lev];

        //        LevelData<FArrayBox>& levelZb        = *m_bedelevation[lev];
        //        LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[11][lev];

        //        DataIterator dit = levelHead.dataIterator();
        //        for (dit.begin(); dit.ok(); ++dit) {
        //            levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
        //            levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
        //            levelResSTP[dit].copy(levelRes[dit], 0, 0, 1);
        //            levelRHSSTP[dit].copy(levelRHS[dit], 0, 0, 1);
        //            levelRHSBSTP[dit].copy(levelRHSB[dit], 0, 0, 1);
        //            levelQwXSTP[dit].copy(levelQw[dit], 0, 0, 1);
        //            levelQwYSTP[dit].copy(levelQw[dit], 1, 0, 1);
        //            levelSourceMSTP[dit].copy(levelSourceM[dit], 0, 0, 1);
        //            levelPwSTP[dit].copy(levelPw[dit], 0, 0, 1);
        //            levelMRSTP[dit].copy(levelMR[dit], 0, 0, 1);
        //            levelReSTP[dit].copy(levelRe[dit], 0, 0, 1);
        //            levelZbSTP[dit].copy(levelZb[dit], 0, 0, 1);
        //        }
        //    } // loop on levs
        //    writePltCustom(nStuffToPlot, vectName, stuffToPlot, std::to_string(ite_idx));
        //} // end customPlt

        if (ite_idx > 500) {
            pout() <<"        does not converge (Picard iterations > 500)."<< endl;
            if (m_verbosity > 0) {
                pout() <<ite_idx<< "         x(h) "<<max_resH<<endl;
            }
            MayDay::Error("Abort");
        } else {
            if (m_cur_step < 50) {
                if (max_resH < 0.01){
                    if (m_verbosity > 0) {
                        pout() <<"        converged( it = "<< ite_idx << ", x(h) = " <<max_resH<< ")."<< endl;
                    }
                    if (m_verbosity > 3) {
                        pout() <<"        Check h (max = "<< maxHead<<", min = " << minHead << ")"<<endl;
                    }
                    converged_h = true;
                }
            } else {
                if (max_resH < m_eps_PicardIte){
                    if (m_verbosity > 0) {
                        pout() <<"        converged( it = "<< ite_idx << ", x(h) = " <<max_resH<< ")."<< endl;
                    }
                    if (m_verbosity > 3) {
                        pout() <<"        Check h (max = "<< maxHead<<", min = " << minHead << ")"<<endl;
                    }
                    converged_h = true;
                }
            }
        }
     
        ite_idx++;
        m_cur_PicardIte++;
        if (m_verbosity > 3) {
            pout() << endl;
        }
    } // end while h solve


   // III b calc

   /* SOLVE FOR B OUTSIDE THE LOOP */
   //     Form RHS for b -- using oldb
   if (m_verbosity > 3) {
       pout() <<"   ...Solve for b ! "<< endl;
       pout() <<"        Update gap height with expl Euler scheme"<< endl;
    }

    //int gh_method = 0; // 0: backward Euler, 1:...
    for (int lev = 0; lev <= m_finest_level; lev++) {
        /* Make sure head GC are properly filled*/
        LevelData<FArrayBox>& levelH   = *m_head[lev];
        // handle ghost cells on the coarse-fine interface
        if (lev > 0) {
            int nGhost = levelH.ghostVect()[0];
            PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                m_amrGrids[lev - 1],
                                                1, // ncomps
                                                m_amrDomains[lev - 1],
                                                m_refinement_ratios[lev - 1],
                                                nGhost);
            headFiller.fillInterp(*m_head[lev], *m_head[lev-1], *m_head[lev-1], 0.0, 0, 0, 1);
        }
        // Need to fill the ghost cells of H -- use BC function 
        levelH.exchange();
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        DataIterator dit = levelGrids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            // get the validBox
            const Box& validBox = levelGrids.get(dit);

            // Fill BC ghost cells of h and b
            mixBCValues(levelH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0], false);
        }

        /* Re evaluate Re with fresh h */
        // CC quantities
        LevelData<FArrayBox>& levelRe       = *m_Re[lev]; 
        // EC quantities
        LevelData<FluxBox>& levelRe_ec      = *a_Re_ec[lev];
        evaluate_Re_quadratic(lev, true);
        // handle ghost cells on the coarse-fine interface
        if (lev > 0) {
            int nGhost = levelRe.ghostVect()[0];
            PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                m_amrGrids[lev - 1],
                                                1, // ncomps
                                                m_amrDomains[lev - 1],
                                                m_refinement_ratios[lev - 1],
                                                nGhost);
            headFiller.fillInterp(*m_Re[lev], *m_Re[lev-1], *m_Re[lev-1], 0.0, 0, 0, 1);
        }
        // Need to fill the ghost cells of Re -- extrapolate on boundaries   
        levelRe.exchange();
        ExtrapGhostCells( levelRe, m_amrDomains[lev]);
        CellToEdge(levelRe, levelRe_ec);

        /* Update Qw */
        // CC quantities
        LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 
        // EC quantities
        LevelData<FluxBox>& levelQw_ec      = *a_Qw_ec[lev]; 
        LevelData<FluxBox>& levelB_ec       = *a_GapHeight_ec[lev];    
        evaluate_Qw_ec(lev, levelQw_ec, levelB_ec, levelRe_ec);
        EdgeToCell(levelQw_ec, levelQw); 
        // handle ghost cells on the coarse-fine interface
        if (lev > 0) {
            int nGhost = levelQw.ghostVect()[0];
            PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                                m_amrGrids[lev - 1],
                                                2, // ncomps
                                                m_amrDomains[lev - 1],
                                                m_refinement_ratios[lev - 1],
                                                nGhost);
            headFiller.fillInterp(*m_qw[lev], *m_qw[lev-1], *m_qw[lev-1], 0.0, 0, 0, 1);
        }
        // Need to fill the ghost cells -- extrapolate on boundaries   
        levelQw.exchange();
        ExtrapGhostCells( levelQw, m_amrDomains[lev]);

        /* Compute grad(Zb) -EC- */ 
        // NOTE does not change right now so useless
        //LevelData<FluxBox>& levelgradZb_ec  = *a_gradZb_ec[lev];
        //compute_grad_zb_ec(lev, levelgradZb_ec);

        /* Reconstruct mR for RHS of B */
        LevelData<FluxBox>&   levelgradH_ec  = *m_gradhead_ec[lev];
        LevelData<FluxBox>&   levelgradZb_ec = *a_gradZb_ec[lev];
        //DisjointBoxLayout&    levelGrids     = m_amrGrids[lev];
        // tmp arrays for vect operations
        LevelData<FluxBox>   leveltmp_ec(levelGrids, 1, IntVect::Zero);
        LevelData<FArrayBox> leveltmp_cc(levelGrids, 1*SpaceDim, HeadGhostVect);
        LevelData<FluxBox>   leveltmp2_ec(levelGrids, 1, IntVect::Zero);
        LevelData<FArrayBox> leveltmp2_cc(levelGrids, 1*SpaceDim, HeadGhostVect);

        //DataIterator dit = levelGrids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            // EC quantities
            FluxBox& Qwater_ec = levelQw_ec[dit];
            FluxBox& gradH_ec  = levelgradH_ec[dit];
            FluxBox& gradZb_ec = levelgradZb_ec[dit];
            FluxBox& tmp_ec    = leveltmp_ec[dit];
            FluxBox& tmp2_ec   = leveltmp2_ec[dit];

            // loop over directions
            for (int dir = 0; dir<SpaceDim; dir++) {

                const Box& region = tmp_ec[dir].box();

                FORT_COMPUTESCAPROD( CHF_FRA(Qwater_ec[dir]),
                                     CHF_FRA(gradH_ec[dir]),
                                     CHF_FRA(gradZb_ec[dir]),
                                     CHF_BOX(region),
                                     CHF_FRA(tmp_ec[dir]),
                                     CHF_FRA(tmp2_ec[dir]) );
            } // loop over dir
        }
        // Qw gradH
        EdgeToCell(leveltmp_ec,  leveltmp_cc);
        // Qw gradZb
        EdgeToCell(leveltmp2_ec, leveltmp2_cc);

        LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
        LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];
        LevelData<FArrayBox>& levelZb    = *m_bedelevation[lev];
        LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
        for (dit.begin(); dit.ok(); ++dit) {
            FArrayBox& newH    = levelH[dit];
            FArrayBox& Pressw  = levelPw[dit];
            FArrayBox& Pressi  = levelPi[dit];
            FArrayBox& zb      = levelZb[dit];
            FArrayBox& mR      = levelmR[dit];
            FArrayBox& tmp_cc  = leveltmp_cc[dit];
            FArrayBox& tmp2_cc = leveltmp2_cc[dit];

            BoxIterator bit(Pressw.box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit(); 
                //Pw
                Pressw(iv,0) = m_suhmoParm->m_gravity * m_suhmoParm->m_rho_w * (newH(iv,0) - zb(iv,0));
                // mR -- turn off fric heat and geot heat
                Real sca_prod = 20. * 20. * m_suhmoParm->m_ub[0] * std::abs(Pressi(iv,0) - Pressw(iv,0)) * m_suhmoParm->m_ub[0];
                mR(iv,0)   = m_suhmoParm->m_G + 
                             sca_prod
                             - m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                             (tmp_cc(iv, 0) + tmp_cc(iv, 1) - tmp2_cc(iv, 0) - tmp2_cc(iv, 1)) -
                             m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * (tmp_cc(iv, 0) + tmp_cc(iv, 1));
                mR(iv,0)   = mR(iv,0) / m_suhmoParm->m_L;
            }
        }

        // 2. a : Get the RHS of gh eqs:
        LevelData<FArrayBox>& levelRHS_b = *RHS_b[lev];
        LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    
        LevelData<FArrayBox>& levelBH    = *m_bumpHeight[lev];    
        LevelData<FArrayBox>& levelBL    = *m_bumpSpacing[lev];    
        LevelData<FArrayBox>& levelDT    = *a_diffusiveTerm[lev];    

        // DEBUG
        LevelData<FArrayBox>& levelRHS_b_A = *RHS_b_A[lev];
        LevelData<FArrayBox>& levelRHS_b_B = *RHS_b_B[lev];
        LevelData<FArrayBox>& levelRHS_b_C = *RHS_b_C[lev];
        LevelData<FArrayBox>& levelchanDegree = *a_chanDegree[lev];

        CalcRHS_gapHeightFAS(levelRHS_b, levelPi, 
                             levelPw, levelmR, 
                             levelB, levelDT, levelBH, levelBL,
                             levelRHS_b_A, levelRHS_b_B, levelRHS_b_C, levelchanDegree); 

        LevelData<FArrayBox>& leveloldB  = *m_old_gapheight[lev];    

        // 2. b : update gap height
        for (dit.begin(); dit.ok(); ++dit) {
 
            FArrayBox& oldB    = leveloldB[dit];
            FArrayBox& newB    = levelB[dit];
            FArrayBox& RHS     = levelRHS_b[dit];

            // DO NOT TOUCH GHOST CELLS
            BoxIterator bit(RHS.box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit(); 
                // does not work
                //newB(iv,0) = std::max(RHS(iv,0) * a_dt + oldB(iv,0), 1.0e-12);
                newB(iv,0) = RHS(iv,0) * a_dt + oldB(iv,0);
            }
        }
    }  // loop on levs


    /* FINAL custom plt here -- debug print */
    if ((m_PrintCustom) && (m_cur_step % m_plot_interval == 0) ) {
        int nStuffToPlot = 10;
        Vector<std::string> vectName;
        vectName.resize(nStuffToPlot);
        vectName[0]="head";
        vectName[1]="gapHeight";
        vectName[2]="RHS_head";
        vectName[3]="RHS_b";
        vectName[4]="RHS_bA";
        vectName[5]="RHS_bB";
        vectName[6]="RHS_bC";
        vectName[7]="RHS_hmoul";
        vectName[8]="Channelization";
        vectName[9]="Qx";

        Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
        stuffToPlot.resize(nStuffToPlot);
        for (int zz = 0; zz < nStuffToPlot; zz++) {
            stuffToPlot[zz].resize(m_max_level + 1, NULL);
        }

        for (int lev = 0; lev <= m_finest_level; lev++) {
            stuffToPlot[0][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[1][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[2][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[3][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[4][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[5][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[6][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[7][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[8][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
            stuffToPlot[9][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);

            LevelData<FArrayBox>& levelHead      = *m_head[lev];
            LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];

            LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
            LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];

            LevelData<FArrayBox>& levelRHS       = *RHS_h[lev];
            LevelData<FArrayBox>& levelRHSSTP    = *stuffToPlot[2][lev];

            LevelData<FArrayBox>& levelRHSB      = *RHS_b[lev];
            LevelData<FArrayBox>& levelRHSBSTP   = *stuffToPlot[3][lev];

            LevelData<FArrayBox>& levelRHSB1     = *RHS_b_A[lev];
            LevelData<FArrayBox>& levelRHSB1STP  = *stuffToPlot[4][lev];

            LevelData<FArrayBox>& levelRHSB2     = *RHS_b_B[lev];
            LevelData<FArrayBox>& levelRHSB2STP  = *stuffToPlot[5][lev];

            LevelData<FArrayBox>& levelRHSB3     = *RHS_b_C[lev];
            LevelData<FArrayBox>& levelRHSB3STP  = *stuffToPlot[6][lev];

            LevelData<FArrayBox>& levelRHSH1     = *moulin_source_term[lev];
            LevelData<FArrayBox>& levelRHSH1STP  = *stuffToPlot[7][lev];

            LevelData<FArrayBox>& levelCD        = *a_chanDegree[lev];
            LevelData<FArrayBox>& levelCDSTP     = *stuffToPlot[8][lev];

            LevelData<FArrayBox>& levelqw        = *m_qw[lev];
            LevelData<FArrayBox>& levelqwSTP     = *stuffToPlot[9][lev];

            DataIterator dit = levelHead.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
                levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
                levelRHSSTP[dit].copy(levelRHS[dit], 0, 0, 1);
                levelRHSBSTP[dit].copy(levelRHSB[dit], 0, 0, 1);
                levelRHSB1STP[dit].copy(levelRHSB1[dit], 0, 0, 1);
                levelRHSB2STP[dit].copy(levelRHSB2[dit], 0, 0, 1);
                levelRHSB3STP[dit].copy(levelRHSB3[dit], 0, 0, 1);
                levelRHSH1STP[dit].copy(levelRHSH1[dit], 0, 0, 1);
                levelCDSTP[dit].copy(levelCD[dit], 0, 0, 1);
                levelqwSTP[dit].copy(levelqw[dit], 0, 0, 1);
            }
        } // loop on levs
        writePltCustom(nStuffToPlot, vectName, stuffToPlot, "");

        for (int lev = 0; lev <= m_finest_level; lev++) {
            delete stuffToPlot[0][lev];
            delete stuffToPlot[1][lev];
            delete stuffToPlot[2][lev];
            delete stuffToPlot[3][lev];
            delete stuffToPlot[4][lev];
            delete stuffToPlot[5][lev];
            delete stuffToPlot[6][lev];
            delete stuffToPlot[7][lev];
            delete stuffToPlot[8][lev];
            delete stuffToPlot[9][lev];
        }
    } // end customPlt



    /* POST PROC -- 1 LEVEL  */
    if (m_post_proc) {
        pout() << "\n\n\n";
        pout() <<"oo POST PROC analysis (domsize is "<< m_amrDomains[0].domainBox().size(0) <<") oo "<< endl;
        
        int DomSize = m_amrDomains[0].domainBox().size(0); 
        int DomSizeY = m_amrDomains[0].domainBox().size(1); 
 
        // DISCHARGE -- Q
        Vector<Real> out_water_flux_x_tot(DomSize , 0.0);
        Vector<Real> out_water_flux_x_distrib(DomSize, 0.0);        // INEFFICIENT 
        Vector<Real> out_water_flux_x_chan(DomSize, 0.0);           // EFFICIENT
        Vector<Real> out_water_flux_x_chan_Bands(3, 0.0);
        Vector<Real> out_water_flux_x_distrib_Bands(3, 0.0);
        int idx_bandMin_lo = 0;
        int idx_bandMin_hi = 0;
        int idx_bandMed_lo = 0;
        int idx_bandMed_hi = 0;
        Real xloc_bandMin_lo = 9999999999;
        Real xloc_bandMin_hi = -1;
        Real xloc_bandMed_lo = 9999999999;
        Real xloc_bandMed_hi = -1;
        // RHS B -- Moulins and mRate
        Vector<Real> out_recharge_tot(DomSize, 0.0);
        Vector<Real> out_recharge(DomSize, 0.0);
        // Pressures
        Vector<Real> avPressure(DomSize, 0.0);
        Vector<Real> out_avgN(4, 0.0);
        int count_avgNLow  = 0;
        int count_avgNMed  = 0;
        int count_avgNHigh = 0;
        int count_avgN     = 0;

        DisjointBoxLayout&    levelGrids = m_amrGrids[0];

        LevelData<FArrayBox>& levelchanDegree = *a_chanDegree[0];
        LevelData<FluxBox>   levelCD_ec(levelGrids, 1, IntVect::Zero);
        levelchanDegree.exchange();
        CellToEdge(levelchanDegree, levelCD_ec);

        LevelData<FArrayBox>& levelMoulinSrc  = *moulin_source_term[0];
        LevelData<FArrayBox>& levelmR         = *m_meltRate[0];
        LevelData<FArrayBox>& levelPressi     = *m_overburdenpress[0];
        LevelData<FArrayBox>& levelPw         = *m_Pw[0];


        LevelData<FluxBox>&   levelQw_ec = *a_Qw_ec[0]; 

        DataIterator dit                 = levelGrids.dataIterator();

        for (dit.begin(); dit.ok(); ++dit) {
            FluxBox&   Qwater_ec    = levelQw_ec[dit];
            FArrayBox& Qwater_ecFab = Qwater_ec[0];
            FluxBox&   CD_ec        = levelCD_ec[dit];
            FArrayBox& CD_ecFab     = CD_ec[0];

            FArrayBox&  MS          = levelMoulinSrc[dit];
            FArrayBox&  MR          = levelmR[dit];
            FArrayBox& Pressi       = levelPressi[dit];
            FArrayBox& Pw           = levelPw[dit];

            const Box& validBox = levelGrids.get(dit);

            BoxIterator bitEC(validBox); 
            for (bitEC.begin(); bitEC.ok(); ++bitEC) {
                IntVect iv = bitEC();
                // DISCHARGE -- Q
                out_water_flux_x_tot[iv[0]]     += Qwater_ecFab(iv, 0) * m_amrDx[0][0];
                out_water_flux_x_chan[iv[0]]    += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * CD_ecFab(iv, 0);
                out_water_flux_x_distrib[iv[0]] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * (1.0 - CD_ecFab(iv, 0));

                // RHS B -- Moulins and mRate
                //CC
                Real xloc = (iv[0]+0.5)*m_amrDx[0][0];    
                // Because of mass conservation eq -basically recharge = RHS of this eq in my head.
                // If I only use the MS, then for run B its sum of A1 and A5 which in this case does NOT equal total discharge !!
                out_recharge[iv[0]]      += MS(iv, 0) * m_amrDx[0][0] * m_amrDx[0][0];  
                out_recharge_tot[iv[0]]  += (MS(iv, 0) + MR(iv, 0)/m_suhmoParm->m_rho_w) * m_amrDx[0][0] * m_amrDx[0][0];

                // Pressures
                avPressure[iv[0]]        += (Pressi(iv, 0) - Pw(iv, 0)); 
                out_avgN[3]              += (Pressi(iv, 0) - Pw(iv, 0)); 
                count_avgN += 1; 
                if ((xloc > 10000) && (xloc < 15000)){
                    out_avgN[0]   += (Pressi(iv, 0) - Pw(iv, 0));
                    count_avgNLow += 1; 
                    xloc_bandMin_lo = std::min(xloc, xloc_bandMin_lo);
                    xloc_bandMin_hi = std::max(xloc, xloc_bandMin_hi);
                    out_water_flux_x_chan_Bands[0] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * CD_ecFab(iv, 0);
                    out_water_flux_x_distrib_Bands[0] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * (1.0 - CD_ecFab(iv, 0));
                } else if ((xloc > 50000) && (xloc < 55000)) {
                    out_avgN[1]   += (Pressi(iv, 0) - Pw(iv, 0));
                    count_avgNMed += 1; 
                    xloc_bandMed_lo = std::min(xloc, xloc_bandMed_lo);
                    xloc_bandMed_hi = std::max(xloc, xloc_bandMed_hi);
                    out_water_flux_x_chan_Bands[1] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * CD_ecFab(iv, 0);
                    out_water_flux_x_distrib_Bands[1] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * (1.0 - CD_ecFab(iv, 0));
                } else if ((xloc > 85000) && (xloc < 90000)) {
                    out_avgN[2]    += (Pressi(iv, 0) - Pw(iv, 0));
                    count_avgNHigh += 1; 
                    out_water_flux_x_chan_Bands[2] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * CD_ecFab(iv, 0);
                    out_water_flux_x_distrib_Bands[2] += Qwater_ecFab(iv, 0) * m_amrDx[0][0] * (1.0 - CD_ecFab(iv, 0));
                }  

                // TEMPORARY
                int time_tmp = (int) m_time;
                if ((iv[0] == 0) && (iv[1] == 2) && (time_tmp % 86400 == 0)) {
                    pout() << "Time(d) moulin = " <<  m_time / 86400. << " " << MS(iv, 0) << endl;
                }
            }
        } 

#ifdef CH_MPI
       // DISCHARGE -- Q
       Vector<Real> recv(DomSize);
       int result = MPI_Allreduce(&out_water_flux_x_tot[0], &recv[0], DomSize, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_water_flux_x_tot = recv;

       Vector<Real> recvChan(DomSize);
       result = MPI_Allreduce(&out_water_flux_x_chan[0], &recvChan[0], DomSize, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_water_flux_x_chan = recvChan;

       Vector<Real> recvDist(DomSize);
       result = MPI_Allreduce(&out_water_flux_x_distrib[0], &recvDist[0], DomSize, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_water_flux_x_distrib = recvDist;

       // RHS B -- Moulins and mRate
       Vector<Real> recvMoulin(DomSize);
       result = MPI_Allreduce(&out_recharge[0], &recvMoulin[0], DomSize, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_recharge = recvMoulin;

       Vector<Real> recvMoulinTot(DomSize);
       result = MPI_Allreduce(&out_recharge_tot[0], &recvMoulinTot[0], DomSize, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_recharge_tot = recvMoulinTot;

       // Pressures
       Vector<Real> recvavP(DomSize);
       result = MPI_Allreduce(&avPressure[0], &recvavP[0], DomSize, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       avPressure = recvavP;

       Vector<Real> recvavN(4);
       result = MPI_Allreduce(&out_avgN[0], &recvavN[0], 4, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_avgN = recvavN;

       Vector<Real> recChanBand(3);
       result = MPI_Allreduce(&out_water_flux_x_chan_Bands[0], &recChanBand[0], 4, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_water_flux_x_chan_Bands= recChanBand;

       Vector<Real> recChanDistr(3);
       result = MPI_Allreduce(&out_water_flux_x_distrib_Bands[0], &recChanDistr[0], 4, MPI_CH_REAL,
                                  MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_SUM");
       }
       out_water_flux_x_distrib_Bands= recChanDistr;

       int recvInt;
       result = MPI_Allreduce(&count_avgNLow, &recvInt, 1, MPI_INT,
                              MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MAX");
       }
       count_avgNLow = recvInt;
       recvInt = 0;
       result = MPI_Allreduce(&count_avgNMed, &recvInt, 1, MPI_INT,
                              MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MAX");
       }
       count_avgNMed = recvInt;
       recvInt = 0;
       result = MPI_Allreduce(&count_avgNHigh, &recvInt, 1, MPI_INT,
                              MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MAX");
       }
       count_avgNHigh = recvInt;
       recvInt = 0;
       result = MPI_Allreduce(&count_avgN, &recvInt, 1, MPI_INT,
                              MPI_SUM, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MAX");
       }
       count_avgN = recvInt;

       Real recvMax;
       result = MPI_Allreduce(&xloc_bandMin_hi, &recvMax, 1, MPI_CH_REAL,
                              MPI_MAX, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MAX");
       }
       xloc_bandMin_hi = recvMax;
       result = MPI_Allreduce(&xloc_bandMed_hi, &recvMax, 1, MPI_CH_REAL,
                              MPI_MAX, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MAX");
       }
       xloc_bandMed_hi = recvMax;

       Real recvMin;
       result = MPI_Allreduce(&xloc_bandMin_lo, &recvMin, 1, MPI_CH_REAL,
                              MPI_MIN, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MIN");
       }
       xloc_bandMin_lo = recvMin;
       result = MPI_Allreduce(&xloc_bandMed_lo, &recvMin, 1, MPI_CH_REAL,
                              MPI_MIN, Chombo_MPI::comm);
       if (result != MPI_SUCCESS)
       {
           MayDay::Error("communication error on MPI_MIN");
       }
       xloc_bandMed_lo = recvMin;
#endif

        for (int xi = (DomSize - 2); xi > -1 ; xi--) {
            out_recharge[xi]  = out_recharge[xi]  + out_recharge[xi+1] ;
            out_recharge_tot[xi]  = out_recharge_tot[xi]  + out_recharge_tot[xi+1] ;
        }
        
        idx_bandMin_lo = (int) xloc_bandMin_lo/m_amrDx[0][0] - 0.5;
        idx_bandMin_hi = (int) xloc_bandMin_hi/m_amrDx[0][0] - 0.5;
        idx_bandMed_lo = (int) xloc_bandMed_lo/m_amrDx[0][0] - 0.5;
        idx_bandMed_hi = (int) xloc_bandMed_hi/m_amrDx[0][0] - 0.5;

        // TEMPORAL POSTPROC
        //int time_tmp = (int) m_time + a_dt;
        //if (time_tmp % 86400 == 0) {
        //    pout() << idx_bandMin_lo << " " << idx_bandMin_hi << " " << idx_bandMed_lo << " "<< idx_bandMed_hi << endl;
        //    pout() << "1-Time  2-disTOT  3-disEFF  4-disINEF  " 
        //           << "5-disEFF_LB_lo  6-disEFF_LB_hi 7-disEFF_LB  "
        //           << "8-disINEF_LB_lo  9-disINEF_LB_hi  10-disINEF_LB"
        //           << " " <<  time_tmp / 86400. << " " << - out_water_flux_x_tot[1]   << " " << -out_water_flux_x_chan[1]   << " " << -out_water_flux_x_distrib[1]
        //           << " " <<  - out_water_flux_x_chan[idx_bandMin_lo]      << " " << - out_water_flux_x_chan[idx_bandMin_hi] 
        //           //<< " " <<  - out_water_flux_x_chan[idx_bandMin_lo] + out_water_flux_x_chan[idx_bandMin_hi]  
        //           << " " << DomSizeY*out_water_flux_x_chan_Bands[0]/count_avgNLow  
        //           << " " <<  - out_water_flux_x_distrib[idx_bandMin_lo]      << " " << - out_water_flux_x_distrib[idx_bandMin_hi] 
        //           //<< " " <<  - out_water_flux_x_distrib[idx_bandMin_lo] + out_water_flux_x_distrib[idx_bandMin_hi] 
        //           << " " << DomSizeY*out_water_flux_x_distrib_Bands[0]/count_avgNLow  
        //           << endl;
        //    pout() << "Time  rechTOT  rechMS  "
        //           << "rechMS_LB_lo    rechMS_LB_hi "
        //           << " " <<  time_tmp / 86400.  << " " <<  out_recharge_tot[0]  << " "  << out_recharge[0] 
        //           << " " <<  out_recharge[idx_bandMin_lo] << " " << out_recharge[idx_bandMin_hi]
        //           << endl;
        //    pout() << "Time  N_LB  N_MB  N_HB  avN  "
        //           << " " << time_tmp / 86400. 
        //           << " " << out_avgN[0]/count_avgNLow << " " << out_avgN[1]/count_avgNMed << " " << out_avgN[2]/count_avgNHigh 
        //           << " " << out_avgN[3] /count_avgN 
        //           << endl;
        //}
 
        // SPATIAL POSTPROC
        for (int xi = 0; xi < DomSize; xi++) {
            Real x_loc = (xi+0.5)*m_amrDx[0][0];    
            pout() << "Xpos avN recharge discharge EFFdischarge INEFdischarge "
                   << " " << x_loc << " " << avPressure[xi]/DomSizeY << " " << out_recharge[xi] 
                   << " " << -out_water_flux_x_tot[xi] << " " << -out_water_flux_x_chan[xi] << " " << -out_water_flux_x_distrib[xi]
                   << endl;
        }
        //pout() <<" -          contrib from CHANNEL (%)   = "<< channel_water_flux_x << " (" << channel_water_flux_x/out_water_flux_x * 100. << ")"<< endl;
        //pout() <<" -          contrib from DISTRIB (%)   = "<< distrib_water_flux_x << " (" << distrib_water_flux_x/out_water_flux_x * 100. << ")"<< endl;
        //pout() << "END\n\n\n";
    }

    /* Averaging down and fill in ghost cells */
    if (m_verbosity > 3) {
        pout() <<"   ...Average down "<< endl;
    }
    for (int lev = m_finest_level; lev > 0; lev--) {
        CoarseAverage averager(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1]);
        averager.averageToCoarse(*m_gapheight[lev - 1], *m_gapheight[lev]);
    }
    // handle ghost cells on the coarse-fine interface
    for (int lev = 1; lev <= m_finest_level; lev++) {
        PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                            m_amrGrids[lev - 1],
                                            1, // ncomps
                                            m_amrDomains[lev - 1],
                                            m_refinement_ratios[lev - 1],
                                            1);
        headFiller.fillInterp(*m_gapheight[lev], *m_gapheight[lev-1], *m_gapheight[lev-1], 0.0, 0, 0, 1);
    }

    /* CONVERGENCE TESTS with LAGGED quantities */
    Real maxGap = computeMax(m_gapheight, m_refinement_ratios, Interval(0,0), 0);
    Real minGap = computeMin(m_gapheight, m_refinement_ratios, Interval(0,0), 0);
    if (m_verbosity > 3) {
        pout() <<"        Check b (max = "<< maxGap<<", min = " << minGap << ")"<<endl;
    }
    for (int lev = 0; lev <= m_finest_level; lev++) {
        LevelData<FArrayBox>& levelnewB_lag  = *a_gapheight_lagged[lev]; // should be the same as the OLD
        LevelData<FArrayBox>& levelcurrentB  = *m_gapheight[lev];

        DataIterator dit = levelnewB_lag.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            levelnewB_lag[dit].minus(levelcurrentB[dit], 0, 0, 1);
            levelnewB_lag[dit].divide(maxGap);
            levelnewB_lag[dit].abs();
        }
    }
    Real max_resB = computeMax(a_gapheight_lagged, m_refinement_ratios, Interval(0,0), 0);
    if (m_verbosity > 3) {
        pout() << "          x(b) "<<max_resB<<endl;
    }

    // finally, update to new time and increment current step
    m_dt = a_dt;
    m_time += a_dt;
    //m_cur_step += 1;

    /* clean up memory */
    pout() <<"   ...clean up memory "<< endl;
    for (int lev = 0; lev <= m_finest_level; lev++) {
            delete a_head_lagged[lev];
            delete a_head_curr[lev];
            delete a_gapheight_lagged[lev];
            delete RHS_h[lev];
            delete moulin_source_term[lev];
            delete RHS_b[lev];
            delete a_ReQwIter[lev];
            delete a_diffusiveTerm[lev];
            delete a_chanDegree[lev];
            delete RHS_b_A[lev];
            delete RHS_b_B[lev];
            delete RHS_b_C[lev];
            delete a_Qw_ec[lev];
            delete a_Re_ec[lev];
            delete a_GapHeight_ec[lev];
            delete a_gradZb_ec[lev];
            delete a_Dcoef[lev];
            //delete aCoef[lev];
            //delete bCoef[lev];
    }

    /* write diagnostic info */
    if (m_verbosity > 0) {
        pout() << "VERBOSE: AmrHydro::timestep " << m_cur_step<< " --     end time = "
               //<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
               << m_time << " ( " << time() << " )"
               //<< " (" << m_time/m_seconds_per_year << " yr)"
               << ", dt = "
               //<< setiosflags(ios::fixed) << setprecision(6) << setw(12)
               << a_dt
               //<< " ( " << a_dt/m_seconds_per_year << " yr )"
               << endl;
    }

    int totalCellsAdvanced = 0;
    for (int lev = 0; lev < m_num_cells.size(); lev++) {
        totalCellsAdvanced += m_num_cells[lev];
    }

    if (m_verbosity > 0) {
        pout() << "Time = " << m_time << " cells advanced = " << totalCellsAdvanced << endl;

        for (int lev = 0; lev < m_num_cells.size(); lev++) {
            pout() << "Time = " << m_time << "  level " << lev << " cells advanced = " << m_num_cells[lev] << endl;
        }
        pout() << endl;
    }
}


// do regridding
void
AmrHydro::regrid()
{
    CH_TIME("AmrHydro::regrid");

    if (m_verbosity > 2) {
        pout() << "AmrHydro::regrid" << endl;
    }

    if (!m_tag_defined) {
        MayDay::Error(" Need a tagging variable ... ");
    }

    // only do any of this if the max level > 0
    if (m_max_level > 0) {
        m_n_regrids++;

        // first generate tags
        Vector<IntVectSet> tagVect(m_max_level);
        tagCells(tagVect);

        // now generate new boxes
        int top_level = min(m_finest_level, m_max_level - 1);
        Vector<Vector<Box> > old_grids(m_finest_level + 1);
        Vector<Vector<Box> > new_grids;

        // this is clunky, but i don't know of a better way to turn
        // a DisjointBoxLayout into a Vector<Box>
        for (int lev = 0; lev <= m_finest_level; lev++) {
            const DisjointBoxLayout& levelDBL = m_amrGrids[lev];
            old_grids[lev].resize(levelDBL.size());
            LayoutIterator lit = levelDBL.layoutIterator();
            int boxIndex = 0;
            for (lit.begin(); lit.ok(); ++lit, ++boxIndex) {
                old_grids[lev][boxIndex] = levelDBL[lit()];
            }
        }

        int new_finest_level;

        BRMeshRefine meshrefine( m_amrDomains[0],  m_refinement_ratios, 
                                 m_fill_ratio,     m_block_factor, 
                                 m_nesting_radius, m_max_box_size);

        new_finest_level = meshrefine.regrid(new_grids,      tagVect, 
                                             m_regrid_lbase, top_level, old_grids);
        if (m_verbosity > 3) {
            pout() << " Old finest level was " << m_finest_level << endl;
            pout() << " New finest level is " << new_finest_level << endl;
        }

        // test to see if grids have changed
        bool gridsSame = true;
        for (int lev = m_regrid_lbase + 1; lev <= new_finest_level; ++lev) {
            int numGridsNew = new_grids[lev].size();
            Vector<int> procIDs(numGridsNew);
            LoadBalance(procIDs, new_grids[lev]);
            const DisjointBoxLayout newDBL(new_grids[lev], procIDs, m_amrDomains[lev]);
            const DisjointBoxLayout oldDBL = m_amrGrids[lev];
            gridsSame &= oldDBL.sameBoxes(newDBL);
        }
        if (gridsSame) {
            if (m_verbosity > 3) {
                pout() << "AmrHydro::regrid -- grids unchanged" << endl;
            }
        } else {
            if (m_verbosity > 3) {
                pout() << "AmrHydro::regrid -- grids changed" << endl;
            }
        }

        // now loop through levels and redefine if necessary
        for (int lev = m_regrid_lbase + 1; lev <= new_finest_level; ++lev) {
            int numGridsNew = new_grids[lev].size();
            Vector<int> procIDs(numGridsNew);
            LoadBalance(procIDs, new_grids[lev]);

            if (m_verbosity > 3) {
                pout() << " Re-defining data on level " << lev << 
                          ". numGrids: " << numGridsNew << endl;
            }

            const DisjointBoxLayout newDBL(new_grids[lev], procIDs, 
                                           m_amrDomains[lev]);

            const DisjointBoxLayout oldDBL = m_amrGrids[lev];

            m_amrGrids[lev] = newDBL;
             
            /* Make sure level exists */
            if (m_head[lev] == NULL) {
                /* NEED */
                m_head[lev]            = new LevelData<FArrayBox>;
                m_gapheight[lev]       = new LevelData<FArrayBox>;
                m_Re[lev]              = new LevelData<FArrayBox>;
                m_bumpHeight[lev]      = new LevelData<FArrayBox>;
                m_bumpSpacing[lev]     = new LevelData<FArrayBox>;
                // CST for now
                m_iceheight[lev]       = new LevelData<FArrayBox>;
                m_bedelevation[lev]    = new LevelData<FArrayBox>;
                m_overburdenpress[lev] = new LevelData<FArrayBox>;
            }

          
            /* NEED VARIABLE */
            // HEAD
            LevelData<FArrayBox>* old_headDataPtr     = m_head[lev];
            LevelData<FArrayBox>* new_headDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_head[0]->nComp(), m_head[0]->ghostVect());
            // old Head will be replaced by copy of head
            m_old_head[lev]        = new LevelData<FArrayBox>(newDBL, m_old_head[0]->nComp(), m_old_head[0]->ghostVect());
            // Gap Height
            LevelData<FArrayBox>* old_gapheightDataPtr    = m_gapheight[lev];
            LevelData<FArrayBox>* new_gapheightDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_gapheight[0]->nComp(), m_gapheight[0]->ghostVect());
            // old gapheight will be replaced by copy of head
            m_old_gapheight[lev]   = new LevelData<FArrayBox>(newDBL, m_old_gapheight[0]->nComp(), m_old_gapheight[0]->ghostVect());
            // Re
            LevelData<FArrayBox>* old_ReDataPtr = m_Re[lev];
            LevelData<FArrayBox>* new_ReDataPtr =
                new LevelData<FArrayBox>(newDBL, m_Re[0]->nComp(), m_Re[0]->ghostVect());
            // Ice Height
            LevelData<FArrayBox>* old_iceheightDataPtr    = m_iceheight[lev];
            LevelData<FArrayBox>* new_iceheightDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_iceheight[0]->nComp(), m_iceheight[0]->ghostVect());
            // Bed elevation
            LevelData<FArrayBox>* old_bedelevationDataPtr    = m_bedelevation[lev];
            LevelData<FArrayBox>* new_bedelevationDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_bedelevation[0]->nComp(), m_bedelevation[0]->ghostVect());
            // Pi
            LevelData<FArrayBox>* old_overburdenpressDataPtr = m_overburdenpress[lev];
            LevelData<FArrayBox>* new_overburdenpressDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_overburdenpress[0]->nComp(), m_overburdenpress[0]->ghostVect());
            // Bed randomization
            LevelData<FArrayBox>* old_bumpHeightDataPtr     = m_bumpHeight[lev];
            LevelData<FArrayBox>* new_bumpHeightDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_bumpHeight[0]->nComp(), m_bumpHeight[0]->ghostVect());
            LevelData<FArrayBox>* old_bumpSpacingDataPtr     = m_bumpSpacing[lev];
            LevelData<FArrayBox>* new_bumpSpacingDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_bumpSpacing[0]->nComp(), m_bumpSpacing[0]->ghostVect());

            // Other vars: grads / Pw / Qw / mR
            if (m_gradhead[lev] != NULL) {
                delete m_gradhead[lev];
                delete m_gradhead_ec[lev];
                delete m_Pw[lev];
                delete m_qw[lev];
                delete m_meltRate[lev];
            }
            m_gradhead[lev]        = new LevelData<FArrayBox>(newDBL, m_gradhead[0]->nComp(), m_gradhead[0]->ghostVect());
            m_gradhead_ec[lev]     = new LevelData<FluxBox>(newDBL, m_gradhead_ec[0]->nComp(), IntVect::Zero);
            m_Pw[lev]              = new LevelData<FArrayBox>(newDBL, m_Pw[0]->nComp(), m_Pw[0]->ghostVect());
            m_qw[lev]              = new LevelData<FArrayBox>(newDBL, m_qw[0]->nComp(), m_qw[0]->ghostVect());
            m_meltRate[lev]        = new LevelData<FArrayBox>(newDBL, m_meltRate[0]->nComp(), m_meltRate[0]->ghostVect());


            /* DEBUG VAR */
            LevelData<FArrayBox>& head  = *new_headDataPtr;
            LevelData<FArrayBox>& gH    = *new_gapheightDataPtr;
            LevelData<FArrayBox>& Re    = *new_ReDataPtr;
            LevelData<FArrayBox>& zB    = *new_bedelevationDataPtr;
            LevelData<FArrayBox>& overP = *new_overburdenpressDataPtr;
            LevelData<FArrayBox>& iceH  = *new_iceheightDataPtr;

            DataIterator dit = newDBL.dataIterator();

            // DEBUG
            if (m_verbosity > 20) {

                pout() << " Checking data beg of process  " << endl;

                for (dit.begin(); dit.ok(); ++dit) {
                    BoxIterator bit(head[dit].box()); 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit(); 
                        pout() << iv << " h: " << head[dit](iv,0) 
                                     << ",b: " << gH[dit](iv,0)           
                                     << ",Re: " << Re[dit](iv,0)           
                                     << ",zB: " << zB[dit](iv,0)           
                                     << ",oP: " << overP[dit](iv,0)           
                                     << ",iceH: " << iceH[dit](iv,0) << endl;
                    }
                }

            } // End verbosity

            // Fill with interpolated data from coarser level
            {
                // may eventually want to do post-regrid smoothing on this!
                FineInterp interpolator(newDBL, 1, m_refinement_ratios[lev - 1], m_amrDomains[lev]);

                //limitSlopes          =  0,
                //noSlopeLimiting      =  1,   XXXXXXXXX
                //PCInterp             =  2,
                //limitTangentialOnly  =  3,
                interpolator.m_boundary_limit_type     = 3;

                /* NEED */
                // HEAD
                interpolator.interpToFine(*new_headDataPtr, *m_head[lev - 1]);
                // Gap Height
                interpolator.interpToFine(*new_gapheightDataPtr, *m_gapheight[lev - 1]);
                // Re
                interpolator.interpToFine(*new_ReDataPtr, *m_Re[lev - 1]);
                // Bed Randomization
                interpolator.interpToFine(*new_bumpHeightDataPtr, *m_bumpHeight[lev - 1]);
                interpolator.interpToFine(*new_bumpSpacingDataPtr, *m_bumpSpacing[lev - 1]);
                /* CST FOR NOW */
                // Ice Height
                interpolator.interpToFine(*new_iceheightDataPtr, *m_iceheight[lev - 1]);
                // Bed elevation
                interpolator.interpToFine(*new_bedelevationDataPtr, *m_bedelevation[lev - 1]);
                // Pi
                interpolator.interpToFine(*new_overburdenpressDataPtr, *m_overburdenpress[lev - 1]);
            }

            // now potentially copy old-grid data on this level into new holder
            if (old_headDataPtr!= NULL) {
                if (oldDBL.isClosed()) {
                    /* NEED */
                    // HEAD
                    old_headDataPtr->copyTo(*new_headDataPtr);
                    // Gap Height
                    old_gapheightDataPtr->copyTo(*new_gapheightDataPtr);
                    // Re
                    old_ReDataPtr->copyTo(*new_ReDataPtr);
                    // Bed Randomization
                    old_bumpHeightDataPtr->copyTo(*new_bumpHeightDataPtr);
                    old_bumpSpacingDataPtr->copyTo(*new_bumpSpacingDataPtr);
                    /* CST FOR NOW */
                    // Ice Height
                    old_iceheightDataPtr->copyTo(*new_iceheightDataPtr);
                    // Bed elevation
                    old_bedelevationDataPtr->copyTo(*new_bedelevationDataPtr);
                    // Pi 
                    old_overburdenpressDataPtr->copyTo(*new_overburdenpressDataPtr);
                }
                // can now delete old data
                delete old_headDataPtr;
                delete old_gapheightDataPtr;
                delete old_ReDataPtr;
                delete old_bumpHeightDataPtr;
                delete old_bumpSpacingDataPtr;
                delete old_iceheightDataPtr;
                delete old_bedelevationDataPtr;
                delete old_overburdenpressDataPtr;
            } 

            // handle ghost cells on the coarse-fine interface
            QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                              m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                              1,  // num comps
                              m_amrDomains[lev]);
            /* NEED */
            // HEAD -- should be useless
            //qcfi.coarseFineInterp(new_headDataPtr, *m_head[lev-1]);           
            // Gap Height -- should be useless
            //qcfi.coarseFineInterp(new_gapheightDataPtr, *m_gapheight[lev-1]);
            // Re -- should be useless
            //qcfi.coarseFineInterp(new_ReDataPtr, *m_Re[lev-1]);
            // Bed Randomization
            qcfi.coarseFineInterp(*new_bumpHeightDataPtr, *m_bumpHeight[lev-1]);
            qcfi.coarseFineInterp(*new_bumpSpacingDataPtr, *m_bumpSpacing[lev-1]);
            /* CST FOR NOW */
            // Ice Height
            qcfi.coarseFineInterp(*new_iceheightDataPtr, *m_iceheight[lev-1]);
            // Bed elevation
            qcfi.coarseFineInterp(*new_bedelevationDataPtr, *m_bedelevation[lev-1]);
            // Pi 
            qcfi.coarseFineInterp(*new_overburdenpressDataPtr, *m_overburdenpress[lev-1]);


            // exchange is necessary to fill periodic ghost cells
            // which aren't filled by the copyTo from oldLevelH or by
            // CF interp OR for when there are several Grids
            /* NEED */
            //new_headDataPtr->exchange();
            //new_gapheightDataPtr->exchange();
            //new_ReDataPtr->exchange();
            new_bumpHeightDataPtr->exchange();
            new_bumpSpacingDataPtr->exchange();
            // CST FOR NOW
            new_iceheightDataPtr->exchange();
            new_bedelevationDataPtr->exchange();
            new_overburdenpressDataPtr->exchange();

            ExtrapGhostCells( *new_bumpHeightDataPtr, m_amrDomains[lev]);
            ExtrapGhostCells( *new_bumpSpacingDataPtr, m_amrDomains[lev]);
            ExtrapGhostCells( *new_iceheightDataPtr, m_amrDomains[lev]);
            ExtrapGhostCells( *new_bedelevationDataPtr, m_amrDomains[lev]);
            ExtrapGhostCells( *new_overburdenpressDataPtr, m_amrDomains[lev]);


            // now place new holders into multilevel arrays
            /* NEED */
            m_head[lev]            = new_headDataPtr;
            m_gapheight[lev]       = new_gapheightDataPtr;
            m_Re[lev]              = new_ReDataPtr;
            m_bumpHeight[lev]      = new_bumpHeightDataPtr;
            m_bumpSpacing[lev]     = new_bumpSpacingDataPtr;
            // CST FOR NOW
            m_iceheight[lev]       = new_iceheightDataPtr;
            m_bedelevation[lev]    = new_bedelevationDataPtr;
            m_overburdenpress[lev] = new_overburdenpressDataPtr;

            // can now delete holders
            delete new_headDataPtr;
            delete new_gapheightDataPtr;
            delete new_ReDataPtr;
            delete new_iceheightDataPtr;
            delete new_bedelevationDataPtr;
            delete new_overburdenpressDataPtr;
            delete new_bumpHeightDataPtr;
            delete new_bumpSpacingDataPtr;


            if (m_verbosity > 20) {

                pout() << " Checking data after Exchanges " << endl;

                for (dit.begin(); dit.ok(); ++dit) {
                    BoxIterator bit(head[dit].box()); 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit(); 
                        pout() << iv << " h: " << head[dit](iv,0) 
                                     << ",b: " << gH[dit](iv,0)           
                                     << ",Re: " << Re[dit](iv,0)           
                                     << ",zB: " << zB[dit](iv,0)           
                                     << ",oP: " << overP[dit](iv,0)           
                                     << ",iceH: " << iceH[dit](iv,0) << endl;
                    }
                }

            } // End verbosity

        } // end loop over currently defined levels

        // now ensure that any remaining levels are null pointers
        // (in case of de-refinement)
        for (int lev = new_finest_level + 1; lev < m_old_head.size(); lev++) {
            if (m_old_head[lev] != NULL) {
                delete m_old_head[lev];
                m_old_head[lev] = NULL;
            }
            if (m_head[lev] != NULL) {
                delete m_head[lev];
                m_head[lev] = NULL;
            }
            if (m_gradhead[lev] != NULL) {
                delete m_gradhead[lev];
                m_gradhead[lev] = NULL;
            }
            if (m_gradhead_ec[lev] != NULL) {
                delete m_gradhead_ec[lev];
                m_gradhead_ec[lev] = NULL;
            }
            if (m_old_gapheight[lev] != NULL) {
                delete m_old_gapheight[lev];
                m_old_gapheight[lev] = NULL;
            }
            if (m_gapheight[lev] != NULL) {
                delete m_gapheight[lev];
                m_gapheight[lev] = NULL;
            }
            if (m_Pw[lev] != NULL) {
                delete m_Pw[lev];
                m_Pw[lev] = NULL;
            }
            if (m_qw[lev] != NULL) {
                delete m_qw[lev];
                m_qw[lev] = NULL;
            }
            if (m_Re[lev] != NULL) {
                delete m_Re[lev];
                m_Re[lev] = NULL;
            }
            if (m_meltRate[lev] != NULL) {
                delete m_meltRate[lev];
                m_meltRate[lev] = NULL;
            }
            if (m_iceheight[lev] != NULL) {
                delete m_iceheight[lev];
                m_iceheight[lev] = NULL;
            }
            if (m_bedelevation[lev] != NULL) {
                delete m_bedelevation[lev];
                m_bedelevation[lev] = NULL;
            }
            if (m_overburdenpress[lev] != NULL) {
                delete m_overburdenpress[lev];
                m_overburdenpress[lev] = NULL;
            }
            if (m_bumpHeight[lev] != NULL) {
                delete m_bumpHeight[lev];
                m_bumpHeight[lev] = NULL;
            }
            if (m_bumpSpacing[lev] != NULL) {
                delete m_bumpSpacing[lev];
                m_bumpSpacing[lev] = NULL;
            }

            DisjointBoxLayout emptyDBL;
            m_amrGrids[lev] = emptyDBL;
        }

        m_finest_level = new_finest_level;

        // set up counter of number of cells
        for (int lev = 0; lev <= m_max_level; lev++) {
            m_num_cells[lev] = 0;
            if (lev <= m_finest_level) {
                const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
                LayoutIterator lit = levelGrids.layoutIterator();
                for (lit.begin(); lit.ok(); ++lit) {
                    const Box& thisBox = levelGrids.get(lit());
                    m_num_cells[lev] += thisBox.numPts();
                }
            }
        }

        // finally, set up covered_level flags
        m_covered_level.resize(m_max_level + 1, 0);
        // note that finest level can't be covered.
        for (int lev = m_finest_level - 1; lev >= 0; lev--) {
            // if the next finer level is covered, then this one is too.
            if (m_covered_level[lev + 1] == 1) {
                m_covered_level[lev] = 1;
            } else {
                // see if the grids finer than this level completely cover it
                IntVectSet fineUncovered(m_amrDomains[lev + 1].domainBox());
                const DisjointBoxLayout& fineGrids = m_amrGrids[lev + 1];

                LayoutIterator lit = fineGrids.layoutIterator();
                for (lit.begin(); lit.ok(); ++lit) {
                    const Box& thisBox = fineGrids.get(lit());
                    fineUncovered.minus_box(thisBox);
                }

                if (fineUncovered.isEmpty()) {
                    m_covered_level[lev] = 1;
                }
            }
        } // end loop over levels to determine covered levels

    } // end if max level > 0 in the first place
}

void
AmrHydro::tagCells(Vector<IntVectSet>& a_tags)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::tagCells" << endl;
    }

    for (int m = 0; m < m_n_tag_var; m++) {
        int top_level = a_tags.size();
        top_level = min(m_tag_cap[m], min(top_level - 1, m_finest_level));
        int min_level = max(m_tag_min[m], 0);
        // loop over levels
        for (int lev = min_level; lev <= top_level; lev++) {
            IntVectSet& levelTags = a_tags[lev];
            tagCellsLevel(levelTags, lev, m);
            IntVectSet& tagSubset = m_vectTagSubset[lev];
            if (tagSubset.numPts() > 0) {
                levelTags &= tagSubset;
            }
        }
    }

    // throw away any coarse level tags outside m_tag_subset
    // if (m_verbosity > 3)
    //   {
    //     pout() << "AmrHydro::tagCells, subset II" << endl;
    //   }
    // if (m_tag_subset.numPts() > 0)
    //   {
    //     IntVectSet tag_subset = m_tag_subset;
    //     a_tags[0] &= tag_subset;
    //     for (int lev = 1; lev <= top_level; lev++)
    // 	{
    // 	  tag_subset.refine(m_refinement_ratios[lev-1]);
    // 	  a_tags[lev] &= tag_subset;
    // 	}

    //   }
}

void
AmrHydro::tagCellsLevel(IntVectSet& a_tags, int a_level, int a_tag)
{
    if (m_verbosity > 4) {
        pout() << "AmrHydro::tagCellsLevel " << a_level << endl;
    }

    // first stab -- don't do BC's; just do one-sided
    // stencils at box edges (hopefully good enough),
    // since doing BC's properly is somewhat expensive.
        if (m_tag_var[a_tag] == "meltingRate") { 
            DataIterator dit = m_meltRate[a_level]->dataIterator();
            LevelData<FArrayBox>& levelPhi = *m_meltRate[a_level];
            const DisjointBoxLayout& levelGrids = m_amrGrids[a_level];

            // need to ensure that ghost cells are set properly
            levelPhi.exchange(levelPhi.interval());

            IntVectSet local_tags;
            for (dit.begin(); dit.ok(); ++dit) {
                FArrayBox& phi   = levelPhi[dit];
                // now tag cells based on values
                BoxIterator bit(levelGrids[dit()]);
                for (bit.begin(); bit.ok(); ++bit) {
                    const IntVect& iv = bit();
                    if (fabs(phi(iv, 0)) > m_tagging_val[a_tag]) local_tags |= iv;
                } // end loop over cells
            }         // end loop over grids

            // now buffer tags
            local_tags.grow(m_tags_grow);
            for (int dir = 0; dir < SpaceDim; dir++) {
                if (m_tags_grow_dir[dir] > m_tags_grow) local_tags.grow(dir, std::max(0, m_tags_grow_dir[dir] - m_tags_grow));
            }
            local_tags &= m_amrDomains[a_level];
            if (a_tags.isEmpty()) {
                a_tags = local_tags;
            } else {
                a_tags |= local_tags;
            }
        } else if (m_tag_var[a_tag] == "someother") {
            /* to do */
        } else {
            MayDay::Error(" Wrong tagging value ... ");
        }

}

void
AmrHydro::tagCellsInit(Vector<IntVectSet>& a_tags)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::tagCellsInit" << endl;
    }

    if (!m_tag_defined) {
        MayDay::Error(" Need a tagging variable ... ");
    }

    tagCells(a_tags);
    m_vectTags = a_tags;
}

void
AmrHydro::initGrids(int a_finest_level) {

    if (m_verbosity > 2) {
        pout() << "AmrHydro::initGrids" << endl;
    }

    m_finest_level = 0;
    // first create base level
    Vector<Box> baseBoxes;
    domainSplit(m_amrDomains[0], baseBoxes, m_max_base_grid_size, m_block_factor);

    Vector<int> procAssign(baseBoxes.size());
    LoadBalance(procAssign, baseBoxes);

    DisjointBoxLayout baseGrids(baseBoxes, procAssign, m_amrDomains[0]);

    if (m_verbosity > 2) {
        long long numCells0 = baseGrids.numCells();
        pout() << "    Level 0: " << numCells0 << " cells: " << baseGrids << endl;
    }

    // Set collection of boxes hierarchy, init level 0 to base grid
    m_amrGrids.resize(m_max_level + 1);
    m_amrGrids[0] = baseGrids;

    // levelSetup now create the data container for each box in level 0
    levelSetup(0, baseGrids);

    // initialize base level data
    initData(m_head);

    int numLevels = 1;
    bool moreLevels = (m_max_level > 0);

    int baseLevel = 0;

    BRMeshRefine meshrefine;
    if (moreLevels) {
        meshrefine.define(
            m_amrDomains[0], m_refinement_ratios, m_fill_ratio, m_block_factor, m_nesting_radius, m_max_box_size);
    }

    Vector<IntVectSet> tagVect(m_max_level);

    Vector<Vector<Box> > oldBoxes(1);
    Vector<Vector<Box> > newBoxes;
    oldBoxes[0] = baseBoxes;
    newBoxes = oldBoxes;
    int new_finest_level = 0;

    while (moreLevels) {
        // default is moreLevels = false
        // (only repeat loop in the case where a new level is generated
        // which is still coarser than maxLevel)
        moreLevels = false;
        tagCellsInit(tagVect);

        // two possibilities -- need to generate grids
        // level-by-level, or we are refining all the
        // way up for the initial time.  check to
        // see which it is by seeing if the finest-level
        // tags are empty
        if (tagVect[m_max_level - 1].isEmpty()) {
            int top_level = m_finest_level;
            int old_top_level = top_level;
            new_finest_level = meshrefine.regrid(newBoxes, tagVect, baseLevel, top_level, oldBoxes);

            if (new_finest_level > top_level) top_level++;
            oldBoxes = newBoxes;

            // now see if we need another pass through grid generation
            if ((top_level < m_max_level) && (top_level > old_top_level)) {
                moreLevels = true;
            }
        } else {
            // for now, define old_grids as just domains
            oldBoxes.resize(m_max_level + 1);
            for (int lev = 1; lev <= m_max_level; lev++) {
                oldBoxes[lev].push_back(m_amrDomains[lev].domainBox());
            }

            int top_level = m_max_level - 1;
            new_finest_level = meshrefine.regrid(newBoxes, tagVect, baseLevel, top_level, oldBoxes);
        }

        numLevels = Min(new_finest_level, m_max_level) + 1;
        if (m_verbosity > 3) {
            pout() << "numLevels " << numLevels << endl;
        }

        // now loop through levels and define
        for (int lev = baseLevel + 1; lev <= new_finest_level; ++lev) {
            int numGridsNew = newBoxes[lev].size();
            Vector<int> procIDs(numGridsNew);
            LoadBalance(procIDs, newBoxes[lev]);
            const DisjointBoxLayout newDBL(newBoxes[lev], procIDs, m_amrDomains[lev]);
            m_amrGrids[lev] = newDBL;

            if (m_verbosity > 2) {
                long long levelNumCells = newDBL.numCells();
                pout() << "   Level " << lev << ": " << levelNumCells << " cells: " << m_amrGrids[lev] << endl;
            }

            levelSetup(lev, m_amrGrids[lev]);

        } // end loop over levels

        m_finest_level = new_finest_level;

        // finally, initialize data on final hierarchy
        // only do this if we've created new levels
        if (m_finest_level > 0) {
            initData(m_head);
        }
    } // end while more levels to do

    for (int lev = m_finest_level+1; lev <= m_max_level; lev++) {
        DisjointBoxLayout emptyDBL;
        m_amrGrids[lev] = emptyDBL;
    }
}


void
AmrHydro::setupFixedGrids(const std::string& a_gridFile)
{
    if (m_verbosity > 2) {
        pout() << "AmrHydro::setupFixedGrids" << endl;
    }
    Vector<Vector<Box> > gridvect;

    if (procID() == uniqueProc(SerialTask::compute)) {
        gridvect.push_back(Vector<Box>(1, m_amrDomains[0].domainBox()));
        // read in predefined grids
        ifstream is(a_gridFile.c_str(), ios::in);
        if (is.fail()) {
            MayDay::Error("Cannot open grids file");
        }

        // format of file:
        //   number of levels, then for each level (starting with level 1):
        //   number of grids on level, list of boxes
        int inNumLevels;
        is >> inNumLevels;
        CH_assert(inNumLevels <= m_max_level + 1);
        if (m_verbosity > 3) {
            pout() << "numLevels = " << inNumLevels << endl;
        }
        while (is.get() != '\n')
            ;
        gridvect.resize(inNumLevels);

        // check to see if coarsest level needs to be broken up
        domainSplit(m_amrDomains[0], gridvect[0], m_max_base_grid_size, m_block_factor);

        // some printing
        if (m_verbosity > 2) {
            pout() << "level 0: ";
            for (int n = 0; n < gridvect[0].size(); n++) {
                pout() << gridvect[0][n] << endl;
            }
        }

        // now loop over levels, starting with level 1
        int numGrids = 0;
        for (int lev = 1; lev < inNumLevels; lev++) {
            is >> numGrids;
            if (m_verbosity >= 3) {
                pout() << "level " << lev << " numGrids = " << numGrids << endl;
                pout() << "Grids: ";
            }
            while (is.get() != '\n')
                ;
            gridvect[lev].resize(numGrids);

            for (int i = 0; i < numGrids; i++) {
                Box bx;
                is >> bx;
                while (is.get() != '\n')
                    ;
                // quick check on box size
                Box bxRef(bx);
                if (bxRef.longside() > m_max_box_size) {
                    pout() << "Grid " << bx << " too large" << endl;
                    MayDay::Error();
                }
                if (m_verbosity >= 3) {
                    pout() << bx << endl;
                }
                gridvect[lev][i] = bx;
            } // end loop over boxes on this level
        }     // end loop over levels
    }         // end if serial proc

    // broadcast results
    broadcast(gridvect, uniqueProc(SerialTask::compute));

    m_amrGrids.resize(m_max_level + 1);
    RealVect dx = m_amrDx[0] * RealVect::Unit;
    for (int lev = 0; lev < gridvect.size(); lev++) {
        int numGridsLev = gridvect[lev].size();
        Vector<int> procIDs(numGridsLev);
        LoadBalance(procIDs, gridvect[lev]);
        const DisjointBoxLayout newDBL(gridvect[lev], procIDs, m_amrDomains[lev]);
        m_amrGrids[lev] = newDBL;
        // build storage for this level
        levelSetup(lev, m_amrGrids[lev]);
        if (lev < gridvect.size() - 1) {
            dx /= m_refinement_ratios[lev];
        }
    }

    // finally set finest level and initialize data on hierarchy
    m_finest_level = gridvect.size() - 1;

    initData(m_head);

    for (int lev = m_finest_level+1; lev <= m_max_level; lev++) {
        DisjointBoxLayout emptyDBL;
        m_amrGrids[lev] = emptyDBL;
    }
}

void
AmrHydro::levelSetup(int a_level, const DisjointBoxLayout& a_grids)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::levelSetup - level " << a_level << endl;
    }
    int nPhiComp = 1;
    IntVect ghostVect = IntVect::Unit;
    IntVect HeadGhostVect = m_num_head_ghost * IntVect::Unit;

    m_old_head[a_level]->define(a_grids, nPhiComp, HeadGhostVect);
    m_old_gapheight[a_level]->define(a_grids, nPhiComp, HeadGhostVect);

    m_head[a_level]         = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_gradhead[a_level]     = new LevelData<FArrayBox>(a_grids, SpaceDim*nPhiComp, HeadGhostVect);
    m_gradhead_ec[a_level]  = new LevelData<FluxBox>(a_grids, nPhiComp, IntVect::Zero);
    m_gapheight[a_level]    = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_Pw[a_level]           = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_qw[a_level]           = new LevelData<FArrayBox>(a_grids, SpaceDim*nPhiComp, HeadGhostVect);
    m_Re[a_level]           = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_meltRate[a_level]     = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_iceheight[a_level]    = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_bedelevation[a_level]    = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_overburdenpress[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_bumpHeight[a_level]      = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
    m_bumpSpacing[a_level]      = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
}

void
AmrHydro::initData(Vector<LevelData<FArrayBox>*>& a_head)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::initData" << endl;
    }

    for (int lev = 0; lev <= m_finest_level; lev++) {
        if (m_verbosity > 3) {
            pout()<< "  .. init lev " << lev << endl;
        }
        LevelData<FArrayBox>& levelHead      = *m_head[lev];
        LevelData<FArrayBox>& levelGapHeight = *m_gapheight[lev];
        LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
        LevelData<FArrayBox>& levelqw        = *m_qw[lev];
        LevelData<FArrayBox>& levelRe        = *m_Re[lev];
        LevelData<FArrayBox>& levelmR        = *m_meltRate[lev];
        LevelData<FArrayBox>& levelIceHeight = *m_iceheight[lev];
        LevelData<FArrayBox>& levelzBed      = *m_bedelevation[lev];
        LevelData<FArrayBox>& levelPi        = *m_overburdenpress[lev];
        LevelData<FArrayBox>& levelbumpHeight   = *m_bumpHeight[lev];
        LevelData<FArrayBox>& levelbumpSpacing  = *m_bumpSpacing[lev];
        //if (lev > 0) {
        //    // fill the ghost cells of a_vectCoordSys[lev]->getH();
        //    // m_head
        //    LevelData<FArrayBox>& coarseHead = *m_head[lev - 1];
        //    int nGhost = levelHead.ghostVect()[0];
        //    PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
        //                                       m_amrGrids[lev - 1],
        //                                       levelHead.nComp(),
        //                                       m_amrDomains[lev - 1],
        //                                       m_refinement_ratios[lev - 1],
        //                                       nGhost);
        //    headFiller.fillInterp(levelHead, coarseHead, coarseHead, 0.0, 0, 0, 1);

        //    // m_gapheight
        //    LevelData<FArrayBox>& coarseGapHeight = *m_gapheight[lev - 1];
        //    //int nGhost = levelGapHeight.ghostVect()[0]; Same number of ghost cells ?
        //    PiecewiseLinearFillPatch gapHeightFiller(m_amrGrids[lev],
        //                                       m_amrGrids[lev - 1],
        //                                       levelGapHeight.nComp(),
        //                                       m_amrDomains[lev - 1],
        //                                       m_refinement_ratios[lev - 1],
        //                                       nGhost);
        //    gapHeightFiller.fillInterp(levelGapHeight, coarseGapHeight, coarseGapHeight, 0.0, 0, 0, 1);
        //}

        RealVect levelDx = m_amrDx[lev] * RealVect::Unit;
        m_IBCPtr->define(m_amrDomains[lev], levelDx[0]);
        // int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:0;

        m_IBCPtr->initializeData(levelDx, 
                                 *m_suhmoParm,       
                                 levelHead, levelGapHeight, 
                                 levelPw, levelqw,
                                 levelRe, levelmR,
                                 levelzBed, levelPi,
                                 levelIceHeight,
                                 levelbumpHeight, levelbumpSpacing); 

        // initialize old h and b to be the current value
        levelHead.copyTo(*m_old_head[lev]);
        levelGapHeight.copyTo(*m_old_gapheight[lev]);
    }

    //writePlotFile();

    // may be necessary to average down here
    //for (int lev = m_finest_level; lev > 0; lev--)
    //{
    //    CoarseAverage avgDown(m_amrGrids[lev], m_head[lev]->nComp(), m_refinement_ratios[lev - 1]);
    //    avgDown.averageToCoarse(*m_head[lev - 1], *m_head[lev]);
    //}

}



void
AmrHydro::BCDataRegrid(int a_level,
                       LevelData<FArrayBox>& a_levzBed,
                       LevelData<FArrayBox>& a_levPi,
                       LevelData<FArrayBox>& a_levIceHeight)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::BCDataRegrid" << endl;
        pout()<<  "  .. lev " << a_level << endl;
    }


    RealVect levelDx = m_amrDx[a_level] * RealVect::Unit;
    m_IBCPtr->BCData(levelDx, m_amrGrids[a_level],
                      m_amrDomains[a_level],  
                      *m_suhmoParm,       
                      a_levzBed, a_levPi,
                      a_levIceHeight);
}


// compute timestep
Real
AmrHydro::computeDt()
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::computeDt" << endl;
    }

    if (m_fixed_dt > TINY_NORM) return m_fixed_dt;

    Real dt = 1.0e50;
    for (int lev = 0; lev <= m_finest_level; lev++) {
        Real dtLev = dt;
        // AF todo
        dt = min(dt, dtLev);
    }

#ifdef CH_MPI
    Real tmp = 1.;
    int result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
    if (result != MPI_SUCCESS)
    {
        MayDay::Error("communication error on norm");
    }
    dt = tmp;
#endif

    if (m_cur_step == 0)
    {
        dt *= m_initial_cfl;
    }
    else
    {
        dt *= m_cfl;
    }

    // also check to see if max grow rate applies
    // (m_dt > 0 test screens out initial time, when we set m_dt to a negative
    // number by default)
    // Use the value stored in m_stable_dt in case dt was altered to hit a plot interval
    // m_max_dt_grow < 0 implies that we don't enforce this.
    if ((m_max_dt_grow > 0) && (dt > m_max_dt_grow * m_stable_dt) && (m_stable_dt > 0))
        dt = m_max_dt_grow * m_stable_dt;

    m_stable_dt = dt;

    return dt; // min(dt,2.0);
}

Real
AmrHydro::computeInitialDt()
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::computeInitialDt" << endl;
    }

    // for now, just call computeDt;
    Real dt = computeDt();
    return dt;
}

#ifdef CH_USE_HDF5
void
AmrHydro::writePltWFX(int numPlotComps, 
                      Vector<std::string>& vectName,
                      Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot,
                      const ProblemDomain& pbDom,
                      const DisjointBoxLayout& dbl,
                      int nbSmooth, int nbDepth)
{
    if (m_verbosity > 4) {
        pout() << "AmrHydro::writePltWFX, AMRFASMGit/ depth = " << nbSmooth << " " << nbDepth <<endl;
    }

    Box domain = pbDom.domainBox();
    Real dt = 1.;
    int numLevels = 1;

    LevelData<FArrayBox>& example   = *stuffToPlot[0][0];

    // Use plot data container for all vars
    Vector<LevelData<FArrayBox>*> plotData(1, NULL);
    IntVect ghostVect(example.ghostVect());

    plotData[0] = new LevelData<FArrayBox>(dbl, numPlotComps, ghostVect);
    LevelData<FArrayBox>& plotDataLev      = *plotData[0];

    DataIterator dit = dbl.dataIterator();
    for (int kk = 0; kk < numPlotComps; kk++) {
        LevelData<FArrayBox>& stuffToPlotLev  = *stuffToPlot[kk][0];
        for (dit.begin(); dit.ok(); ++dit) {
            FArrayBox& thisPlotData      = plotDataLev[dit];
            FArrayBox& thisstuffToPlot   = stuffToPlotLev[dit];
            thisPlotData.copy(thisstuffToPlot, 0, kk, 1);
        }  // end loop over boxes on this level
    } // end loop over vars to add to pltData

    // generate plotfile name
    char iter_str[100];
    sprintf(iter_str, "%s%06d_PI%03d_WFX-AMRFASMG%01d_DEPTH%01d", m_plot_prefix.c_str(), m_cur_step, m_cur_PicardIte, nbSmooth, nbDepth);

    string filename(iter_str);
    
    //filename.append(namePlot);
    filename.append(".hdf5");

    WriteAMRHierarchyHDF5(
        filename, m_amrGrids, plotData, vectName, domain, m_amrDx[0][0], dt, time(), m_refinement_ratios, numLevels);

    // need to delete plotData
    for (int lev = 0; lev < numLevels; lev++)
    {
        if (plotData[lev] != NULL)
        {
            delete plotData[lev];
            plotData[lev] = NULL;
        }
    }
}

// write custom debug hdf5 plotfile to the standard location
void
AmrHydro::writePltCustom(int numPlotComps, 
                         Vector<std::string>& vectName,
                         Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot,
                         string namePlot)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::writePltCustom" << endl;
    }

    Box domain = m_amrDomains[0].domainBox();
    Real dt = 1.;
    int numLevels = m_finest_level + 1;

    // Use plot data container for all vars
    Vector<LevelData<FArrayBox>*> plotData(m_head.size(), NULL);
    IntVect ghostVect(IntVect::Unit);
    for (int lev = 0; lev < numLevels; lev++) {
        plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], numPlotComps, ghostVect);

        LevelData<FArrayBox>& plotDataLev      = *plotData[lev];

        DataIterator dit = m_amrGrids[lev].dataIterator();
        for (int kk = 0; kk < numPlotComps; kk++) {
            LevelData<FArrayBox>& stuffToPlotLev  = *stuffToPlot[kk][lev];
            for (dit.begin(); dit.ok(); ++dit) {
                FArrayBox& thisPlotData      = plotDataLev[dit];
                FArrayBox& thisstuffToPlot   = stuffToPlotLev[dit];
                thisPlotData.copy(thisstuffToPlot, 0, kk, 1);
            }  // end loop over boxes on this level
        } // end loop over vars to add to pltData

    } // end loop over levels for computing plot data

    // generate plotfile name
    char iter_str[100];
    sprintf(iter_str, "%sCustom%06d", m_plot_prefix.c_str(), m_cur_step);

    string filename(iter_str);
    
    filename.append(namePlot);

    if (SpaceDim == 1)
    {
        filename.append(".hdf5");
    }
    else if (SpaceDim == 2)
    {
        filename.append(".hdf5");
    }
    else if (SpaceDim == 3)
    {
        filename.append(".hdf5");
    }

    WriteAMRHierarchyHDF5(
        filename, m_amrGrids, plotData, vectName, domain, m_amrDx[0][0], dt, time(), m_refinement_ratios, numLevels);

    // need to delete plotData
    for (int lev = 0; lev < numLevels; lev++)
    {
        if (plotData[lev] != NULL)
        {
            delete plotData[lev];
            plotData[lev] = NULL;
        }
    }
}


// write hdf5 plotfile to the standard location
void
AmrHydro::writePlotFile()
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::writePlotFile" << endl;
    }

    // plot comps: head + gapHeight + bedelevation + overburdenPress
    int numPlotComps = 12;

    // generate data names
    string headName("head");
    string gapHeightName("gapHeight");
    string zbedName("bedelevation");
    string piName("overburdenPress");
    string pwName("Pw"); 
    string qwName("Qw_x"); 
    string qwNameY("Qw_y"); 
    string ReName("Re"); 
    string meltRateName("meltRate"); 
    string xGradName("GradHead_x");
    string yGradName("GradHead_y");
    string BumpHeight("bumpHeight");

    Vector<string> vectName(numPlotComps);
    vectName[0] = headName;
    vectName[1] = gapHeightName;
    vectName[2] = zbedName;
    vectName[3] = piName;
    vectName[4] = pwName;
    vectName[5] = qwName;
    vectName[6] = qwNameY;
    vectName[7] = ReName;
    vectName[8] = meltRateName;
    vectName[9] = xGradName;
    vectName[10] = yGradName;
    vectName[11] = BumpHeight;

    Box domain = m_amrDomains[0].domainBox();
    Real dt = 1.;
    int numLevels = m_finest_level + 1;
    // compute plot data
    Vector<LevelData<FArrayBox>*> plotData(m_head.size(), NULL);

    // ghost vect makes things simpler
    IntVect ghostVect(IntVect::Unit);

    for (int lev = 0; lev < numLevels; lev++) {
        // first allocate storage
        plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], numPlotComps, ghostVect);
    }

    for (int lev = 0; lev < numLevels; lev++) {
        // now copy new-time solution into plotData
        LevelData<FArrayBox>& plotDataLev = *plotData[lev];

        const LevelData<FArrayBox>& levelHead      = *m_head[lev];
        const LevelData<FArrayBox>& levelGradHead  = *m_gradhead[lev];
        const LevelData<FArrayBox>& levelgapHeight = *m_gapheight[lev];
        const LevelData<FArrayBox>& levelzbed      = *m_bedelevation[lev];
        const LevelData<FArrayBox>& levelpi        = *m_overburdenpress[lev];
        const LevelData<FArrayBox>& levelpw        = *m_Pw[lev];
        const LevelData<FArrayBox>& levelqw        = *m_qw[lev];
        const LevelData<FArrayBox>& levelRe        = *m_Re[lev];
        const LevelData<FArrayBox>& levelmR        = *m_meltRate[lev];
        const LevelData<FArrayBox>& levelBH        = *m_bumpHeight[lev];
        //LevelData<FArrayBox> levelGradPhi;
        //if (m_write_gradPhi)
        //{
        //    levelGradPhi.define(m_amrGrids[lev], SpaceDim, ghostVect);
        //    // compute cc gradient here
        //    DataIterator dit = levelGradPhi.dataIterator();
        //    for (dit.begin(); dit.ok(); ++dit)
        //    {
        //        const FArrayBox& thisPhi = levelHead[dit];
        //        FArrayBox& thisGrad = levelGradPhi[dit];
        //        // just loop over interiors
        //        BoxIterator bit(m_amrGrids[lev][dit]);
        //        for (bit.begin(); bit.ok(); ++bit)
        //        {
        //            IntVect iv = bit();
        //            for (int dir = 0; dir < SpaceDim; dir++)
        //            {
        //                IntVect plusIV = iv + BASISV(dir);
        //                IntVect minusIV = iv - BASISV(dir);
        //                thisGrad(iv, dir) = 0.5 * (thisPhi(plusIV, 0) - thisPhi(minusIV, 0)) / m_amrDx[lev][dir];
        //            }
        //        } // end loop over cells
        //    }     // end loop over grids
        //}         // end if write_gradPhi

        DataIterator dit = m_amrGrids[lev].dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            //const Box& gridBox        = m_amrGrids[lev][dit];
            FArrayBox& thisPlotData   = plotDataLev[dit];
            int comp = 0;
            const FArrayBox& thisHead       = levelHead[dit];
            const FArrayBox& thisGradHead   = levelGradHead[dit];
            const FArrayBox& thisGapHeight  = levelgapHeight[dit];
            const FArrayBox& thiszbed       = levelzbed[dit];
            const FArrayBox& thisPi         = levelpi[dit];
            const FArrayBox& thisPw         = levelpw[dit];
            const FArrayBox& thisqw         = levelqw[dit];
            const FArrayBox& thisRe         = levelRe[dit];
            const FArrayBox& thismR         = levelmR[dit];
            const FArrayBox& thisBH         = levelBH[dit];

            thisPlotData.copy(thisHead, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisGapHeight, 0, comp, 1);
            comp++;
            thisPlotData.copy(thiszbed, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisPi, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisPw, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisqw, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisqw, 1, comp, 1);
            comp++;
            thisPlotData.copy(thisRe, 0, comp, 1);
            comp++;
            thisPlotData.copy(thismR, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisGradHead, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisGradHead, 1, comp, 1);
            comp++;
            thisPlotData.copy(thisBH, 0, comp, 1);

        } // end loop over boxes on this level

        // this is just so that visit surface plots look right
        // fill coarse-fine ghost-cell values with interpolated data
        //if (lev > 0)
        //{
        //    PiecewiseLinearFillPatch interpolator(m_amrGrids[lev],
        //                                          m_amrGrids[lev - 1],
        //                                          numPlotComps,
        //                                          m_amrDomains[lev - 1],
        //                                          m_refinement_ratios[lev - 1],
        //                                          ghostVect[0]);

        //    // no interpolation in time
        //    Real time_interp_coeff = 0.0;
        //    interpolator.fillInterp(
        //        *plotData[lev], *plotData[lev - 1], *plotData[lev - 1], time_interp_coeff, 0, 0, numPlotComps);
        //}
        // just in case...
        // plotData[lev]->exchange();
    } // end loop over levels for computing plot data

    // generate plotfile name
    char iter_str[100];
    sprintf(iter_str, "%s%06d.", m_plot_prefix.c_str(), m_cur_step);

    string filename(iter_str);

    if (SpaceDim == 1) {
        filename.append("1d.hdf5");
    }
    else if (SpaceDim == 2) {
        filename.append("2d.hdf5");
    }
    else if (SpaceDim == 3) {
        filename.append("3d.hdf5");
    }

    WriteAMRHierarchyHDF5(
        filename, m_amrGrids, plotData, vectName, domain, m_amrDx[0][0], dt, time(), m_refinement_ratios, numLevels);

    // need to delete plotData
    for (int lev = 0; lev < numLevels; lev++) {
        if (plotData[lev] != NULL) {
            delete plotData[lev];
            plotData[lev] = NULL;
        }
    }
}

// write checkpoint file out for later restarting
void
AmrHydro::writeCheckpointFile() const
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::writeCheckpointfile" << endl;
    }

    CH_TIME("AmrHydro::writeCheckpointFile");

    char* iter_str;
    if (m_check_overwrite) {
        // overwrite the same checkpoint file, rather than re-writing them
        std::string fs("%s.%dd.hdf5");
        iter_str = new char[m_check_prefix.size() + fs.size() + 16];
        sprintf(iter_str, "%s.%dd.hdf5", m_check_prefix.c_str(), SpaceDim);
    } else {
        // or hang on to them
        std::string fs;
        fs.assign("%s%06d.%dd.hdf5");
        iter_str = new char[m_check_prefix.size() + fs.size() + 16];
        sprintf(iter_str, "%s%06d.%dd.hdf5", m_check_prefix.c_str(), m_cur_step, SpaceDim);
    }

    if (m_verbosity > 3) {
        pout() << "checkpoint file name = " << iter_str << endl;
    }

    // Equivalent of AmrIce::writeCheckpointFile(const string& a_file) in Bisicles

    HDF5Handle handle(iter_str, HDF5Handle::CREATE);

    // write amr data -- only dump out things which are essential
    // to restarting the computation (i.e. max_level, finest_level,
    // time, refinement ratios, etc.).  Other parameters (regrid
    // intervals, block-factor, etc can be changed by the inputs
    // file of the new run.
    // At the moment, the maximum level is not allowed to change,
    // although in principle, there is no real reason why it couldn't
    //
    HDF5HeaderData header;
    header.m_int["max_level"] = m_max_level;
    header.m_int["finest_level"] = m_finest_level;
    header.m_int["current_step"] = m_cur_step;
    header.m_real["time"] = m_time;
    header.m_real["dt"] = m_dt;
    header.m_int["num_comps"] = 8; // H/B/Pice/Zb/Re/Hice/BH/BL
    header.m_real["cfl"] = m_cfl;

    // periodicity info
    D_TERM(if (m_amrDomains[0].isPeriodic(0)) 
               header.m_int["is_periodic_0"] = 1; 
           else 
               header.m_int["is_periodic_0"] = 0; ,

           if (m_amrDomains[0].isPeriodic(1)) 
               header.m_int["is_periodic_1"] = 1;
           else 
               header.m_int["is_periodic_1"] = 0; ,

           if (m_amrDomains[0].isPeriodic(2)) 
               header.m_int["is_periodic_2"] = 1;
           else 
               header.m_int["is_periodic_2"] = 0;
           );

    // set up component names
    string headName("head");
    string gapHeightName("gapHeight");
    string piName("overburdenPress");
    string zbedName("bedelevation");
    string ReName("Re"); 
    string IceHeightName("iceHeight"); 
    string bumpHeightName("bumpHeight"); 
    string bumpSpacingName("bumpSpacing"); 

    int nComp = 0;
    char compStr[30];
    // HEAD
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = headName;
    nComp++;
    // GapHeight
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = gapHeightName;
    nComp++;
    // Pi
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = piName;
    nComp++;
    // Zb
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = zbedName;
    nComp++;
    // Re
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = ReName;
    nComp++;
    // Hice
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = IceHeightName;
    nComp++;
    // Bump Height (BH)
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = bumpHeightName;
    nComp++;
    // Bump Spacing (BL)
    sprintf(compStr, "component_%04d", nComp);
    header.m_string[compStr] = bumpSpacingName;

    header.writeToFile(handle);

    // now loop over levels and write out each level's data
    // note that we loop over all allowed levels, even if they
    // are not defined at the moment.
    for (int lev = 0; lev <= m_max_level; lev++) {
        // set up the level string
        char levelStr[20];
        sprintf(levelStr, "%d", lev);
        const std::string label = std::string("level_") + levelStr;

        handle.setGroup(label);

        // set up the header info
        HDF5HeaderData levelHeader;
        if (lev < m_max_level) {
            levelHeader.m_int["ref_ratio"] = m_refinement_ratios[lev];
        }
        levelHeader.m_real["dx"] = m_amrDx[lev][0];
        levelHeader.m_box["prob_domain"] = m_amrDomains[lev].domainBox();

        levelHeader.writeToFile(handle);

        // now write the data for this level
        // only try to write data if level is defined.
        if (lev <= m_finest_level) {
            write(handle, m_amrGrids[lev]);
            const IntVect ghost = IntVect::Unit * 1;
            write(handle, *m_head[lev], "headData", m_head[lev]->ghostVect());
            write(handle, *m_gapheight[lev], "gapHeightData", m_gapheight[lev]->ghostVect());
            write(handle, *m_overburdenpress[lev], "overburdenPressData", m_overburdenpress[lev]->ghostVect());
            write(handle, *m_bedelevation[lev], "bedelevationData", m_bedelevation[lev]->ghostVect());
            write(handle, *m_Re[lev], "ReData", m_Re[lev]->ghostVect());
            write(handle, *m_iceheight[lev], "iceHeightData", m_iceheight[lev]->ghostVect());
            write(handle, *m_bumpHeight[lev], "bumpHeightData", m_bumpHeight[lev]->ghostVect());
            write(handle, *m_bumpSpacing[lev], "bumpSpacingData", m_bumpSpacing[lev]->ghostVect());
        } // end loop over levels
    }

    handle.close();
}

// read checkpoint file for restart
void
AmrHydro::readCheckpointFile(HDF5Handle& a_handle)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::readCheckpointFile" << endl;
    }

    CH_TIME("AmrHydro::readCheckpointFile");

    HDF5HeaderData& header = m_headerData;
    header.readFromFile(a_handle);
    ParmParse ppAmr("amr");
    
    if (m_verbosity >= 3) {
        pout() << "hdf5 header data: " << endl;
        pout() << header << endl;
    }

    // read max level
    if (header.m_int.find("max_level") == header.m_int.end()) {
        MayDay::Error("checkpoint file does not contain max_level");
    }
    // we can change max level upon restart
    int max_level_check = header.m_int["max_level"];
    if (max_level_check != m_max_level) {
        if (m_verbosity > 0) {
            pout() << "Restart file has a different max level than inputs file" << endl;
            pout() << "     ** max level from inputs file = "   << m_max_level << endl;
            pout() << "     ** max level in checkpoint file = " << max_level_check << endl;
            pout() << "Using max level from inputs file" << endl;
        }
    }
    // read finest level
    if (header.m_int.find("finest_level") == header.m_int.end()) {
        MayDay::Error("checkpoint file does not contain finest_level");
    }
    m_finest_level = header.m_int["finest_level"];
    if (m_finest_level > m_max_level) {
        MayDay::Error("finest level in restart file > inputs file max allowable level!");
    }

    // read current step
    if (header.m_int.find("current_step") == header.m_int.end()){
        MayDay::Error("checkpoint file does not contain current_step");
    }
    m_cur_step = header.m_int["current_step"];
    m_restart_step = m_cur_step;
    // optionally, over-ride the step number in the restart checkpoint file with one specified in the inputs
    if (ppAmr.contains("restart_step") ) {
        int restart_step;
        ppAmr.get("restart_step", restart_step);
        m_cur_step = restart_step;
        m_restart_step = restart_step;
    }
    if (m_verbosity >= 3) {
        pout() << "restart step ? " << m_cur_step << endl;
    }

    // read time
    if (header.m_real.find("time") == header.m_real.end()) {
        MayDay::Error("checkpoint file does not contain time");
    }
    m_time = header.m_real["time"];
    //optionally, over-ride the time in the restart checkpoint file with one specified in the inputs 
    if (ppAmr.contains("restart_time") ) { 
        bool set_time = true;
        ppAmr.query("restart_set_time",set_time); // set amr.restart_set_time = false to prevent time reset
        if (set_time){
            Real restart_time;
            ppAmr.get("restart_time", restart_time);
            m_time = restart_time;
        }
    }

    // read timestep
    if (header.m_real.find("dt") == header.m_real.end()) {
        MayDay::Error("checkpoint file does not contain dt");
    }
    m_dt = header.m_real["dt"];

    // read num comps
    if (header.m_int.find("num_comps") == header.m_int.end()) {
        MayDay::Error("checkpoint file does not contain num_comps");
    }

    // read cfl
    if (header.m_real.find("cfl") == header.m_real.end()) {
        MayDay::Error("checkpoint file does not contain cfl");
    }
    Real check_cfl = header.m_real["cfl"];
    ParmParse ppCheck("AMRDriver");
    if (ppCheck.contains("cfl")) {
        // check for consistency and warn if different
        if (check_cfl != m_cfl) {
            if (m_verbosity > 0) {
                pout() << "CFL in checkpoint file different from inputs file" << endl;
                pout() << "     cfl in inputs file = " << m_cfl << endl;
                pout() << "     cfl in checkpoint file = " << check_cfl << endl;
                pout() << "Using cfl from inputs file" << endl;
            }
        } // end if cfl numbers differ
    } else { //cfl not present in inputs file
        m_cfl = check_cfl;
    }

    // read periodicity info
    // Get the periodicity info -- this is more complicated than it really
    // needs to be in order to preserve backward compatibility
    bool isPeriodic[SpaceDim];
    D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end())) 
               isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
           else 
               isPeriodic[0] = false; ,

           if (!(header.m_int.find("is_periodic_1") == header.m_int.end())) 
               isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
           else 
               isPeriodic[1] = false; ,

           if (!(header.m_int.find("is_periodic_2") == header.m_int.end())) 
               isPeriodic[2] = (header.m_int["is_periodic_2"] == 1);
           else 
               isPeriodic[2] = false;);

    // now resize stuff
    m_amrDomains.resize(m_max_level + 1);
    m_amrGrids.resize(m_max_level + 1);
    m_amrDx.resize(m_max_level + 1);
    m_old_head.resize(m_max_level + 1, NULL);
    m_head.resize(m_max_level + 1, NULL);
    m_gradhead.resize(m_max_level + 1, NULL);
    m_gradhead_ec.resize(m_max_level + 1, NULL);
    m_old_gapheight.resize(m_max_level + 1, NULL);
    m_gapheight.resize(m_max_level + 1, NULL);
    m_Pw.resize(m_max_level + 1, NULL);
    m_qw.resize(m_max_level + 1, NULL);
    m_Re.resize(m_max_level + 1, NULL);
    m_meltRate.resize(m_max_level + 1, NULL);
    m_iceheight.resize(m_max_level + 1, NULL);
    m_bedelevation.resize(m_max_level + 1, NULL);
    m_overburdenpress.resize(m_max_level + 1, NULL);
    m_bumpHeight.resize(m_max_level + 1, NULL);
    m_bumpSpacing.resize(m_max_level + 1, NULL);

    // now read in level-by-level data -- go to Max lev of checkfile
    for (int lev = 0; lev <= max_level_check; lev++) {
        // set up the level string
        char levelStr[20];
        sprintf(levelStr, "%d", lev);
        const std::string label = std::string("level_") + levelStr;

        a_handle.setGroup(label);

        // read header info
        HDF5HeaderData levheader;
        levheader.readFromFile(a_handle);

        if (m_verbosity >= 3) {
            pout() << "level " << lev << " header data" << endl;
            pout() << levheader << endl;
        }

        // Get the refinement ratio
        if (lev < max_level_check) {
            int checkRefRatio;
            if (levheader.m_int.find("ref_ratio") == levheader.m_int.end()) {
                MayDay::Error("checkpoint file does not contain ref_ratio");
            }
            checkRefRatio = levheader.m_int["ref_ratio"];

            // check for consistency
            if (checkRefRatio != m_refinement_ratios[lev]){
                MayDay::Error("inputs file and checkpoint file ref ratios inconsistent");
            }
        }

        // read dx
        if (levheader.m_real.find("dx") == levheader.m_real.end()) {
            MayDay::Error("checkpoint file does not contain dx");
        }
        // TODO check why that Abs doesn't work
        //if (lev <= max_level_check) {
        //    if ( Abs(m_amrDx[lev] - levheader.m_real["dx"]) > TINY_NORM ) {
        //        MayDay::Error("restart file dx != input file dx");
        //    }
        //}
        m_amrDx[lev] = RealVect::Unit * (levheader.m_real["dx"]);

        // read problem domain box
        if (levheader.m_box.find("prob_domain") == levheader.m_box.end()) {
            MayDay::Error("checkpoint file does not contain prob_domain");
        }
        Box domainBox = levheader.m_box["prob_domain"];

        m_amrDomains[lev] = ProblemDomain(domainBox, isPeriodic);

        // the rest is only applicable if this level is defined
        if (lev <=  m_finest_level && lev <= max_level_check) {
            // read grids
            Vector<Box> grids;
            const int grid_status = read(a_handle, grids);
            if (grid_status != 0) {
                MayDay::Error("checkpoint file does not contain a Vector<Box>");
            }
            // do load balancing
            int numGrids = grids.size();
            Vector<int> procIDs(numGrids);
            LoadBalance(procIDs, grids);
            DisjointBoxLayout levelDBL(grids, procIDs, m_amrDomains[lev]);
            m_amrGrids[lev] = levelDBL;

            // allocate this level's storage
            IntVect nGhost  = m_num_head_ghost * IntVect::Unit;
            int nPhiComp = 1;
            m_old_head[lev]      = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_head[lev]          = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_gradhead[lev]      = new LevelData<FArrayBox>(levelDBL, SpaceDim*nPhiComp, nGhost);
            m_gradhead_ec[lev]   = new LevelData<FluxBox>(levelDBL, nPhiComp, IntVect::Zero);
            m_old_gapheight[lev] = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_gapheight[lev]     = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_Pw[lev]            = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_qw[lev]            = new LevelData<FArrayBox>(levelDBL, SpaceDim*nPhiComp, nGhost);
            m_Re[lev]            = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_meltRate[lev]      = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_iceheight[lev]     = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_bedelevation[lev]  = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_overburdenpress[lev]   = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);
            m_bumpHeight[lev]        = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost); 
            m_bumpSpacing[lev]       = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost); 

            // read this level's data
            /* HEAD */
            LevelData<FArrayBox>& old_head = *m_old_head[lev];
            LevelData<FArrayBox> tmpHead;
            tmpHead.define(old_head);
            int dataStatus = read<FArrayBox>(a_handle, tmpHead, "headData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain head data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                old_head[dit].copy(tmpHead[dit]);
                (*m_head[lev])[dit].copy(tmpHead[dit]);
            }
            /* B */
            LevelData<FArrayBox> tmpB;
            tmpB.define(old_head);
            dataStatus = read<FArrayBox>(a_handle, tmpB, "gapHeightData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain gap height data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_old_gapheight[lev])[dit].copy(tmpB[dit]);
                (*m_gapheight[lev])[dit].copy(tmpB[dit]);
            }
            /* Re */
            LevelData<FArrayBox> tmpRe;
            tmpRe.define(old_head);
            dataStatus = read<FArrayBox>(a_handle, tmpRe, "ReData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain Re data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_Re[lev])[dit].copy(tmpRe[dit]);
            }
            /* Ice Height */
            LevelData<FArrayBox> tmpIceHeight;
            tmpIceHeight.define(old_head);
            dataStatus = read<FArrayBox>(a_handle, tmpIceHeight, "iceHeightData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain H data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_iceheight[lev])[dit].copy(tmpIceHeight[dit]);
            }
            /* Ice Pressure */
            LevelData<FArrayBox> tmpPi;
            tmpPi.define(old_head);
            dataStatus = read<FArrayBox>(a_handle, tmpPi, "overburdenPressData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain Pi data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_overburdenpress[lev])[dit].copy(tmpPi[dit]);
            }
            /* Bed Topo */
            LevelData<FArrayBox> tmpZb;
            tmpZb.define(old_head);
            dataStatus = read<FArrayBox>(a_handle, tmpZb, "bedelevationData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain ZB data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_bedelevation[lev])[dit].copy(tmpZb[dit]);
            }
            /* Bump characteristic length scales */
            LevelData<FArrayBox> tmpBH;
            tmpBH.define(old_head);
            dataStatus = read<FArrayBox>(a_handle, tmpBH, "bumpHeightData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain Bump Height data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_bumpHeight[lev])[dit].copy(tmpBH[dit]);
            }
            dataStatus = read<FArrayBox>(a_handle, tmpBH, "bumpSpacingData", levelDBL);
            if (dataStatus != 0) {
                MayDay::Error("checkpoint file does not contain Bump Spacing data");
            }
            for (DataIterator dit(levelDBL); dit.ok(); ++dit) {
                (*m_bumpSpacing[lev])[dit].copy(tmpBH[dit]);
            }

        } // end if this level is defined
    }     // end loop over levels

}

// set up for restart
void
AmrHydro::restart(string& a_restart_file)
{
    if (m_verbosity > 2) {
        pout() << "AmrHydro::restart" << endl;
    }

    HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
    // first read in data from checkpoint file
    readCheckpointFile(handle);
    handle.close();
}

#endif

#include "NamespaceFooter.H"
