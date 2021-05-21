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

#include "AmrHydro.H"

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

void
ParseBC() 
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

void DirichletValue(Real* pos,
                    int* dir, 
                    Side::LoHiSide* side, 
                    Real* a_values)
{
    a_values[0] = DirichletValue(*dir, *side);
}

// could be same as below
void NeumBCForHcoarse(FArrayBox& a_state,
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


void FixedNeumBCFill(FArrayBox& a_state,
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


void BCFill(FArrayBox& a_state,
            const Box& a_valid,
            const ProblemDomain& a_domain,
            Real a_dx)
{
  // If box is outside of domain bounds ?
  if(!a_domain.domainBox().contains(a_state.box())) {

      if (!GlobalBCRS::s_areBCsParsed) {
          ParseBC();
      }

      for(int dir=0; dir<CH_SPACEDIM; ++dir) {
          // don't do anything if periodic -- should be perio in y dir 1
          if (!a_domain.isPeriodic(dir)) {
              Box ghostBoxLo = adjCellBox(a_valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(a_valid, dir, Side::Hi, 1);

              if ((!a_domain.domainBox().contains(ghostBoxLo)) && (a_state.box().contains(ghostBoxLo)) ) {
                  if (GlobalBCRS::s_bcLo[dir] == 0) {
                      // Diri
                      ghostBoxLo &= a_state.box();
                      int isign = sign(Side::Lo);
                      for (BoxIterator bit(ghostBoxLo); bit.ok(); ++bit) {
                          IntVect ivTo = bit();
                          IntVect ivClose = ivTo -   isign*BASISV(dir);
                          IntVect ivFar   = ivTo - 2*isign*BASISV(dir);
                          //pout() << "iv ivclose ivFar" << ivTo << " " << ivClose << " " << ivFar <<endl;
                          Real nearVal = a_state(ivClose, 0);
                          Real farVal  = a_state(ivFar,   0);
                          Real inhomogVal = DirichletValue(dir, Side::Lo);
                          // linear or quad
                          //Real ghostVal =  2.0 * inhomogVal - nearVal;
                          Real ghostVal =  (8.0 / 3.0) * inhomogVal + (1.0 / 3.0) * farVal - 2.0 * nearVal;
                          a_state(ivTo, 0) = ghostVal;
                      }
                  } else if (GlobalBCRS::s_bcLo[dir] == 1) {
                      // Neum
                      Box fromRegion = ghostBoxLo;
                      int isign = sign(Side::Lo);
                      fromRegion.shift(dir, -isign);
                      a_state.copy(a_state, fromRegion, 0, ghostBoxLo, 0, a_state.nComp());
                  } // End BC options
              } // End ghostBoxLo in dir

              if ((!a_domain.domainBox().contains(ghostBoxHi)) && (a_state.box().contains(ghostBoxHi)) ) {
                  if (GlobalBCRS::s_bcHi[dir] == 0) {
                      // Diri
                      ghostBoxHi &= a_state.box();
                      int isign = sign(Side::Hi);
                      for (BoxIterator bit(ghostBoxHi); bit.ok(); ++bit) {
                          IntVect ivTo = bit();
                          IntVect ivClose = ivTo -   isign*BASISV(dir);
                          IntVect ivFar   = ivTo - 2*isign*BASISV(dir);
                          Real nearVal = a_state(ivClose, 0);
                          Real farVal  = a_state(ivFar,   0);
                          Real inhomogVal = DirichletValue(dir, Side::Hi);
                          // linear or quad
                          //Real ghostVal =  2.0 * inhomogVal - nearVal;
                          Real ghostVal =  (8.0 / 3.0) * inhomogVal + (1.0 / 3.0) * farVal - 2.0 * nearVal;
                          a_state(ivTo, 0) = ghostVal;
                      }
                  } else if (GlobalBCRS::s_bcHi[dir] == 1) {
                      // Neum
                      Box fromRegion = ghostBoxHi;
                      int isign = sign(Side::Hi);
                      fromRegion.shift(dir, -isign);
                      a_state.copy(a_state, fromRegion, 0, ghostBoxHi, 0, a_state.nComp());
                  } // End BC options
              } // End ghostBoxHi in dir

          } // end if is not periodic in ith direction
      } // end dir loop
  }
}

void 
mixBCValues(FArrayBox& a_state,
            const Box& a_valid,
            const ProblemDomain& a_domain,
            Real a_dx,
            bool a_homogeneous)
{
  // If box is outside of domain bounds ?
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
    Real eps        =  1.0e-10;  // solution tolerance
    Real hang       =  0.01;      // next rnorm should be < (1-m_hang)*norm_last 
    Real normThresh =  1.0e-16;  // abs tol
    amrSolver->setSolverParameters(numSmooth, numSmooth, numBottom,
                                   numMG, maxIter, 
                                   eps, hang, normThresh);

    if (m_verbosity > 3) {
        amrSolver->m_verbosity = 4;
    } else {
        amrSolver->m_verbosity = 1;
    } 

    // solve !
    bool zeroInitialGuess = false;
    amrSolver->solve(a_head, a_RHS, m_finest_level, 0, zeroInitialGuess);
}

/* AmrHydro class functions */
void
AmrHydro::SolveForHead(const Vector<DisjointBoxLayout>&               a_grids,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                      Vector<ProblemDomain>& a_domains,
                      Vector<int>& refRatio,
                      Real& coarsestDx,
                      Vector<LevelData<FArrayBox>*>& a_head, 
                      Vector<LevelData<FArrayBox>*>& a_RHS)
{
    VCAMRPoissonOp2Factory* poissonOpF_head = new VCAMRPoissonOp2Factory;

    //BCHolder bc(ConstDiriNeumBC(IntVect::Unit, RealVect::Unit,  IntVect::Zero, RealVect::Zero));
    poissonOpF_head->define(a_domains[0],
                      a_grids,
                      refRatio,
                      coarsestDx,
                      &mixBCValues, //bc
                      0.0,
                      a_aCoef,
                      - 1.0,
                      a_bCoef);

    RefCountedPtr< AMRLevelOpFactory<LevelData<FArrayBox> > > opFactoryPtr(poissonOpF_head);

    MultilevelLinearOp<FArrayBox> poissonOp;
    // options ?
    poissonOp.m_num_mg_iterations = 3;
    poissonOp.m_num_mg_smooth = 4;
    poissonOp.m_preCondSolverDepth = -1;
    poissonOp.define(a_grids, refRatio, a_domains, m_amrDx, opFactoryPtr, 0);
    // bottom solver ?
    BiCGStabSolver<Vector<LevelData<FArrayBox>* > > solver;
    bool homogeneousBC = false;  
    solver.define(&poissonOp, homogeneousBC); 
    solver.m_normType = 0;
    if (m_verbosity > 3) {
        solver.m_verbosity = 4;
    } else {
        solver.m_verbosity = 1;
    }
    //solver.m_eps = 1.0e-6;
    solver.m_reps = 1.0e-10;
    solver.setConvergenceMetrics(1.0, 1.0e-10);
    solver.m_imax = 100;
    //
    solver.solve(a_head, a_RHS);
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
    m_is_defined = false;
    m_verbosity = 4;
    m_max_level = -1;
    m_finest_level = -1;
    m_tag_cap = 100;
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
    m_tag_var = "none";
    m_tagging_val = 1.0;
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
    ppAmr.query("verbosity", m_verbosity);
    ppAmr.query("block_factor", m_block_factor);
    ppAmr.query("regrid_lbase", m_regrid_lbase);  // smaller lev subject to regridding
    ppAmr.query("regrid_interval", m_regrid_interval);
    ppAmr.query("tagCap", m_tag_cap);               // no idea yet
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
    ppAmr.query("tag_variable", m_tag_var);
    // if we set this to be true, require that we also provide the threshold
    if (m_tag_var.compare("none") != 0) {
        ppAmr.query("tagging_val", m_tagging_val);
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
    // AF that part is a little fuzzy
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
    //MayDay::Error("STOP");
}

// set BC for head ?
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
            if ((m_use_FAS) && (m_use_NL)) {
                timeStepFAS(dt);
            } else {
                timeStep(dt);
            }

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
        FArrayBox& GradH   = lvlgradH[dit];
        FArrayBox& B       = a_B[dit];
        FArrayBox& Re      = lvlRe[dit];

        BoxIterator bit(Re.box()); // can use gridBox? 
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            Real sqrt_gradH_cc = std::sqrt(GradH(iv, 0) * GradH(iv, 0) + GradH(iv, 1) * GradH(iv, 1));
            Real discr = 1.0 + 4.0 * m_suhmoParm->m_omega * (
                        std::pow(B(iv, 0), 3) * m_suhmoParm->m_gravity * sqrt_gradH_cc) / (
                        12.0 * m_suhmoParm->m_nu * m_suhmoParm->m_nu);  
            Re(iv, 0) = (- 1.0 + std::sqrt(discr)) / (2.0 * m_suhmoParm->m_omega) ; 
        }
    }
    if (a_print_WFX) {
        /* custom plt here -- debug print */
            int nStuffToPlot = 4;
            Vector<std::string> vectName;
            vectName.resize(nStuffToPlot);
            vectName[0]="GapHeight";
            vectName[1]="Head";
            vectName[2]="Reynolds";
            vectName[3]="GradH";

            Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
            stuffToPlot.resize(nStuffToPlot);
            for (int zz = 0; zz < nStuffToPlot; zz++) {
                stuffToPlot[zz].resize(1, NULL);
            }

            stuffToPlot[0][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());
            stuffToPlot[1][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());
            stuffToPlot[2][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());
            stuffToPlot[3][0]  = new LevelData<FArrayBox>(levelGrids, 1, a_u.ghostVect());

            LevelData<FArrayBox>& levelGapSTP   = *stuffToPlot[0][0];
            LevelData<FArrayBox>& levelHeadSTP  = *stuffToPlot[1][0];
            LevelData<FArrayBox>& levelReSTP    = *stuffToPlot[2][0];
            LevelData<FArrayBox>& levelGradHSTP = *stuffToPlot[3][0];

            DataIterator dit = levelHeadSTP.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                levelGapSTP[dit].copy(a_B[dit], 0, 0, 1);
                levelHeadSTP[dit].copy(lcl_u[dit], 0, 0, 1);
                levelReSTP[dit].copy(lvlRe[dit], 0, 0, 1);
                levelGradHSTP[dit].copy(lvlgradH[dit], 0, 0, 1);
            }
            writePltWFX(nStuffToPlot, vectName, stuffToPlot, 
                        levelDomain, levelGrids, 
                        a_smooth, a_depth);
    }
    
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
            FArrayBox& bCFab = bC[dir];
            FArrayBox& BFab  = B_ec[dir];
            FArrayBox& ReFab = Re_ec[dir];

            // initialize 
            bCFab.setVal(0.0);

            BoxIterator bit(bCFab.box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit();
                // Update b coeff
                Real num_q = - std::pow(BFab(iv, 0),3) * m_suhmoParm->m_gravity;
                Real denom_q = 12.0 * m_suhmoParm->m_nu * (1 + m_suhmoParm->m_omega * ReFab(iv, 0));
                bCFab(iv, 0) = num_q/denom_q;
                if (bCFab(iv, 0) >= 0.0) {
                    pout() << " --(dx= "<< a_dx<<") cell dir Bc-EC " << iv << " " << dir << " " << bCFab(iv, 0) << "\n";
                    //MayDay::Error("Abort");
                }
            }
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
      for (levelDit.begin(); levelDit.ok(); ++levelDit) {

          FArrayBox& thisNL          = a_NL[levelDit];
          FArrayBox& thisdNL         = a_dNL[levelDit];
          const FArrayBox& thisU     = a_u[levelDit];
          FArrayBox& thisB           = a_B[levelDit];
          FArrayBox& thisPi          = a_Pi[levelDit];
          FArrayBox& thiszb          = a_zb[levelDit];

          BoxIterator bit(thisNL.box());
          for (bit.begin(); bit.ok(); ++bit) {
              IntVect iv = bit();
              thisNL(iv, 0)  = - m_suhmoParm->m_A * thisB(iv,0) * 
                               std::pow( (thisPi(iv,0) - m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                               (thisU(iv,0) -thiszb(iv,0))), 3);
              thisdNL(iv, 0) = 3.0 * m_suhmoParm->m_A * thisB(iv,0) * m_suhmoParm->m_rho_w *m_suhmoParm->m_gravity *
                                 ( std::pow( (thisPi(iv,0) - m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * 
                                 (thisU(iv,0) - thiszb(iv,0))), 2) );
          }
      } // end loop over grids on this level
  }
}

void AmrHydro::NonLinear(Vector<LevelData<FArrayBox>* >        a_NL,
                         Vector<LevelData<FArrayBox>* >        a_dNL, 
                         const Vector<LevelData<FArrayBox>* >  a_u,
                         int a_finestLevel)
{

  for (int lev=0; lev<=a_finestLevel; lev++) {

    LevelData<FArrayBox>& levelNL   = *(a_NL[lev]);
    LevelData<FArrayBox>& leveldNL  = *(a_dNL[lev]);
    LevelData<FArrayBox>& levelU    = *(a_u[lev]);

    LevelData<FArrayBox>& levelB     = *m_gapheight[lev];
    LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];
    LevelData<FArrayBox>& levelZb    = *m_bedelevation[lev];

    DataIterator levelDit = levelNL.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit) {

        FArrayBox& thisNL    = levelNL[levelDit];
        FArrayBox& thisdNL   = leveldNL[levelDit];
        FArrayBox& thisU     = levelU[levelDit];

        FArrayBox& B         = levelB[levelDit];
        FArrayBox& Pressi    = levelPi[levelDit];
        FArrayBox& Zb        = levelZb[levelDit];

        BoxIterator bit(thisNL.box());
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            pout() << "Press ? " << iv << " " << Pressi(iv,0) << "\n";
            thisNL(iv, 0)  = m_suhmoParm->m_A * B(iv,0) * 
                             std::pow( (Pressi(iv,0) - m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                             (thisU(iv,0) - Zb(iv,0))), 3);
            thisdNL(iv, 0) = - 3.0 * m_suhmoParm->m_A * B(iv,0) * m_suhmoParm->m_rho_w *m_suhmoParm->m_gravity *
                             ( std::pow( (Pressi(iv,0) - m_suhmoParm->m_rho_w * m_suhmoParm->m_G * 
                             (thisU(iv,0) - Zb(iv,0))), 2));
        }
    } // end loop over grids on this level
  } // end loop over levels
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

            FArrayBox& bCFab = bC[dir];
            FArrayBox& BFab  = B[dir];
            FArrayBox& ReFab = Re[dir];

            // initialize 
            bCFab.setVal(0.0);

            BoxIterator bit(bCFab.box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit();
                // Update b coeff
                Real num_q = - std::pow(BFab(iv, 0),3) * m_suhmoParm->m_gravity;
                Real denom_q = 12.0 * m_suhmoParm->m_nu * (1 + m_suhmoParm->m_omega * ReFab(iv, 0));
                bCFab(iv, 0) = num_q/denom_q;
            }
        }
    }
}


void
AmrHydro::aCoeff_bCoeff_CC(LevelData<FArrayBox>&  levelacoef, 
                           LevelData<FArrayBox>&  levelbcoef_cc, 
                           LevelData<FArrayBox>&  levelRe, 
                           LevelData<FArrayBox>&  levelB)
{
    DataIterator dit = levelbcoef_cc.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {

        FArrayBox& B      = levelB[dit];
        FArrayBox& aC     = levelacoef[dit];
        FArrayBox& bC_cc  = levelbcoef_cc[dit];
        FArrayBox& Re     = levelRe[dit];

        // initialize 
        aC.setVal(0.0);
        bC_cc.setVal(0.0);

        BoxIterator bit(bC_cc.box()); 
        for (bit.begin(); bit.ok(); ++bit) {
            IntVect iv = bit();
            // Update b coeff
            Real num_q = - std::pow(B(iv, 0),3) * m_suhmoParm->m_gravity;
            Real denom_q = 12.0 * m_suhmoParm->m_nu * (1 + m_suhmoParm->m_omega * Re(iv, 0));
            bC_cc(iv, 0) = num_q/denom_q;
        }
    }
}

void 
AmrHydro::Calc_moulin_source_term (LevelData<FArrayBox>& levelMoulinSrc) 
{
   // This is a huge hack for placing Moulin source term. 
   // Always make sure you have enough cells to distribute mdot ! 4 at lev 0, should be in same bit
   int curr_level = 0;
   DataIterator dit = levelMoulinSrc.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {
       FArrayBox& moulinSrc = levelMoulinSrc[dit];

       moulinSrc.setVal(0.0);

       BoxIterator bit(moulinSrc.box());
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           // low face
           Real x_lo = (iv[0])*m_amrDx[curr_level][0];
           Real y_lo = (iv[1])*m_amrDx[curr_level][1];
           // hi face
           Real x_hi = (iv[0]+1)*m_amrDx[curr_level][0];
           Real y_hi = (iv[1]+1)*m_amrDx[curr_level][1];

           // Loop over all moulins
           for (int m = 0; m<m_suhmoParm->m_n_moulins; m++) {
               if ( (m_suhmoParm->m_moulin_position[m*2 + 0] <= x_hi) 
                     && (m_suhmoParm->m_moulin_position[m*2 + 0] > x_lo) ) {
                   if ( (m_suhmoParm->m_moulin_position[m*2 + 1] <= y_hi) 
                         && (m_suhmoParm->m_moulin_position[m*2 + 1] > y_lo) ) {
                        moulinSrc(iv,0) = m_suhmoParm->m_moulin_flux[m] / 
                                 ( m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        // Y line
                        //     o
                        //     o
                        // - - X - - 
                        //     o
                        //     o
                        //int coord_TL[3];
                        //coord_TL[0] = iv[0]; 
                        //coord_TL[1] = iv[1] + 1; 
                        //coord_TL[2] = 0;
                        //IntVect iv_TL(coord_TL);
                        //moulinSrc(iv_TL,0) = m_suhmoParm->m_moulin_flux / 
                        //             (5.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //coord_TL[1] = iv[1] - 1; 
                        //IntVect iv_TL2(coord_TL);
                        //moulinSrc(iv_TL2,0) = m_suhmoParm->m_moulin_flux / 
                        //             (5.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //coord_TL[1] = iv[1] - 2; 
                        //IntVect iv_TL3(coord_TL);
                        //moulinSrc(iv_TL3,0) = m_suhmoParm->m_moulin_flux / 
                        //             (9.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //coord_TL[1] = iv[1] + 2; 
                        //IntVect iv_TL4(coord_TL);
                        //moulinSrc(iv_TL4,0) = m_suhmoParm->m_moulin_flux / 
                        //             (9.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        // X line
                        //     I
                        //     I
                        // o o X o o 
                        //     I
                        //     I
                        //int coord_TR[3];
                        //coord_TR[0] = iv[0] + 1; 
                        //coord_TR[1] = iv[1]; 
                        //coord_TR[2] = 0;
                        //IntVect iv_TR(coord_TR);
                        //moulinSrc(iv_TR,0) = m_suhmoParm->m_moulin_flux / 
                        //             (3.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //coord_TR[0] = iv[0] - 1; 
                        //IntVect iv_TR2(coord_TR);
                        //moulinSrc(iv_TR2,0) = m_suhmoParm->m_moulin_flux / 
                        //             (3.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //coord_TR[0] = iv[0] - 2; 
                        //IntVect iv_TR3(coord_TR);
                        //moulinSrc(iv_TR3,0) = m_suhmoParm->m_moulin_flux / 
                        //             (9.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //coord_TR[0] = iv[0] + 2; 
                        //IntVect iv_TR4(coord_TR);
                        //moulinSrc(iv_TR4,0) = m_suhmoParm->m_moulin_flux / 
                        //             (9.0 * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]);  // m3s-1 / m2
                        //pout() << "level is: " << curr_level << ", IV: " << iv << endl; 
                   }
               }
           }
       }
   }
   for (dit.begin(); dit.ok(); ++dit) {
       FArrayBox& moulinSrc = levelMoulinSrc[dit];
       BoxIterator bit(moulinSrc.box());
       Real integ_moulin = 0.0;
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           integ_moulin += moulinSrc(iv,0) * m_amrDx[curr_level][0] * m_amrDx[curr_level][1]; 
       }
       if (m_verbosity > 3) {
           pout() << "level is: " << curr_level << ", Integral of moulin vol. rate: " << integ_moulin << endl; 
       } 
   }
}

void
AmrHydro::CalcRHS_head(int curr_level, 
                       LevelData<FArrayBox>& levelRHS_h, 
                       LevelData<FArrayBox>& levelPi, 
                       LevelData<FArrayBox>& levelPw, 
                       LevelData<FArrayBox>& levelmR, 
                       LevelData<FArrayBox>& levelB,
                       LevelData<FArrayBox>& levelMoulinSrc)
{
   DataIterator dit = levelRHS_h.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

       FArrayBox& B       = levelB[dit];
       FArrayBox& RHS     = levelRHS_h[dit];

       FArrayBox& Pressi    = levelPi[dit];
       FArrayBox& Pw        = levelPw[dit];
       FArrayBox& meltR     = levelmR[dit];
       FArrayBox& moulinSrc = levelMoulinSrc[dit];
       
       // initialize RHS for h
       RHS.setVal(0.0);
       // first term
       RHS.copy(meltR, 0, 0, 1);
       RHS *= (1.0 /  m_suhmoParm->m_rho_w - 1.0 / m_suhmoParm->m_rho_i);
       // third term ...
       Real ub_norm = std::sqrt(  m_suhmoParm->m_ub[0]*m_suhmoParm->m_ub[0] 
                                + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]) / m_suhmoParm->m_lr;
       BoxIterator bit(RHS.box()); // can use gridBox? 
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           
           if ( B(iv,0) < m_suhmoParm->m_br) {
               RHS(iv,0) -= ub_norm * (m_suhmoParm->m_br - B(iv,0));
           }
           // second term ... assume  n = 3 !!
           Real PimPw = (Pressi(iv,0) - Pw(iv,0));
           Real AbsPimPw = std::abs(PimPw);
           RHS(iv,0) += m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0);
           // Add moulin 
           RHS(iv,0) += moulinSrc(iv,0);
       }
   }
}

void
AmrHydro::CalcRHS_gapHeight_semiExpl(LevelData<FArrayBox>& levelRHS_b, 
                                     LevelData<FArrayBox>& levelDENOM_b, 
                                     LevelData<FArrayBox>& levelPi, 
                                     LevelData<FArrayBox>& levelPw, 
                                     LevelData<FArrayBox>& levelmR, 
                                     LevelData<FArrayBox>& levelB,
                                     Real a_dt,
                                     LevelData<FArrayBox>& levelRHS_b_A)
{

   DataIterator dit = levelRHS_b.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

       FArrayBox& B       = levelB[dit];
       FArrayBox& RHS     = levelRHS_b[dit];
       FArrayBox& DENOM   = levelDENOM_b[dit];

       FArrayBox& RHS_A   = levelRHS_b_A[dit];

       FArrayBox& Pressi  = levelPi[dit];
       FArrayBox& Pw      = levelPw[dit];
       FArrayBox& meltR   = levelmR[dit];

       // initialize RHS for h
       RHS.setVal(0.0);
       DENOM.setVal(0.0);
       RHS_A.setVal(0.0);

       // first term
       RHS.copy(meltR, 0, 0, 1);
       RHS *= 1.0 / m_suhmoParm->m_rho_i;

       // third term ...
       Real ub_norm = std::sqrt(  m_suhmoParm->m_ub[0]*m_suhmoParm->m_ub[0] 
                                + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]) / m_suhmoParm->m_lr;
       BoxIterator bit(RHS.box()); // can use gridBox? 
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           if ( B(iv,0) < m_suhmoParm->m_br) {
               RHS(iv,0) += ub_norm * (m_suhmoParm->m_br - B(iv,0));
           }
           RHS(iv,0) = B(iv,0) + RHS(iv,0) * a_dt;
           // second term ... assume  n = 3 !!
           Real PimPw = (Pressi(iv,0) - Pw(iv,0));
           Real AbsPimPw = std::abs(PimPw);
           DENOM(iv,0) = 1.0 + a_dt * m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw;
           RHS_A(iv,0) = 1.0 + a_dt * m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw;
       }
   }
}  


void
AmrHydro::CalcRHS_gapHeightFAS(LevelData<FArrayBox>& levelRHS_b, 
                               LevelData<FArrayBox>& levelPi, 
                               LevelData<FArrayBox>& levelPw, 
                               LevelData<FArrayBox>& levelmR, 
                               LevelData<FArrayBox>& levelB)
{
   DataIterator dit = levelRHS_b.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

       FArrayBox& B       = levelB[dit];
       FArrayBox& RHS     = levelRHS_b[dit];

       FArrayBox& Pressi  = levelPi[dit];
       FArrayBox& Pw      = levelPw[dit];
       FArrayBox& meltR   = levelmR[dit];
       
       // initialize RHS for h
       RHS.setVal(0.0);

       // first term
       RHS.copy(meltR, 0, 0, 1);
       RHS *= 1.0 / m_suhmoParm->m_rho_i;
       // third term ...
       Real ub_norm = std::sqrt(  m_suhmoParm->m_ub[0]*m_suhmoParm->m_ub[0] 
                                + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]) / m_suhmoParm->m_lr;
       BoxIterator bit(RHS.box()); // can use gridBox? 
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           if ( B(iv,0) < m_suhmoParm->m_br) {
               RHS(iv,0) += ub_norm * (m_suhmoParm->m_br - B(iv,0));
           }
           // second term ... assume  n = 3 !!
           Real PimPw = (Pressi(iv,0) - Pw(iv,0));
           Real AbsPimPw = std::abs(PimPw);
           RHS(iv,0) -= m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0);
       }
   }
}


void
AmrHydro::CalcRHS_gapHeight(LevelData<FArrayBox>& levelRHS_b, 
                            LevelData<FArrayBox>& levelPi, 
                            LevelData<FArrayBox>& levelPw, 
                            LevelData<FArrayBox>& levelmR, 
                            LevelData<FArrayBox>& levelB,
                            LevelData<FArrayBox>& levelRHS_b_A,
                            LevelData<FArrayBox>& levelRHS_b_B,
                            LevelData<FArrayBox>& levelRHS_b_C)
{
   DataIterator dit = levelRHS_b.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

       FArrayBox& B       = levelB[dit];
       FArrayBox& RHS     = levelRHS_b[dit];

       FArrayBox& RHS_A   = levelRHS_b_A[dit];
       FArrayBox& RHS_B   = levelRHS_b_B[dit];
       FArrayBox& RHS_C   = levelRHS_b_C[dit];

       FArrayBox& Pressi  = levelPi[dit];
       FArrayBox& Pw      = levelPw[dit];
       FArrayBox& meltR   = levelmR[dit];
       
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
                                + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]) / m_suhmoParm->m_lr;
       BoxIterator bit(RHS.box()); // can use gridBox? 
       for (bit.begin(); bit.ok(); ++bit) {
           IntVect iv = bit();
           if ( B(iv,0) < m_suhmoParm->m_br) {
               RHS(iv,0) += ub_norm * (m_suhmoParm->m_br - B(iv,0));
               RHS_B(iv,0) = ub_norm * (m_suhmoParm->m_br - B(iv,0));
           }
           // second term ... assume  n = 3 !!
           Real PimPw = (Pressi(iv,0) - Pw(iv,0));
           Real AbsPimPw = std::abs(PimPw);
           RHS(iv,0) -= m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0);
           RHS_C(iv,0) = - m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0);
       }
   }
}

/* core of advance routine */
void
AmrHydro::timeStep(Real a_dt)
{
    CH_TIME("AmrHydro::timestep");

    m_cur_step += 1;

    if (m_verbosity > 2) {
        pout() << "\n\n-- Timestep " << m_cur_step << " Advancing solution from time " << m_time << " ( " << time()
               << ")"
                  " with dt = "
               << a_dt << endl;
    }
    
    /* Sketch of loops: Fig 1 of Sommers */
    //
    // NOTE: all computation done at CC. When edge qties needed, use CC->Edge
    // I Copy new h and new b into old h and old b
    //     Take care of GC/BC etc.
    // II LOOP: h and b calc
    //     Fill perio GC and BC of h and b
    //     Put h into h_lag (GC too) and b into b_lag
    //     Interp CC b to EC
    //     Update water pressure Pw=f(h)
    //     Compute grad(h) and grad(Pw)
    //     IV INNER LOOP: Re/Qw !! (fixed nb of iter for now)
    //             Update VECTOR Qw = f(Re, grad(h))
    //             Update Re = f(Qw)
    //     Update melting rate = f(Qw, grad(h), grad(Pw))
    //     Compute aCoeff and bCoeff, CC RHS for h and solve for h again
    //     Form RHS for b
    //     Solve for b using Forward Euler simple scheme
    //     Check convergence using h and h_lag AND b and b_lag
    //  
    /* End comments */

    IntVect HeadGhostVect = m_num_head_ghost * IntVect::Unit;
    Real coarsestDx = m_amrDx[0][0];  

    /* I Copy new h and new b into old h and old b */

    // Also create and initialize tmp vectors
    Vector<LevelData<FArrayBox>*> a_head_lagged;
    Vector<LevelData<FArrayBox>*> a_gapheight_lagged;
    Vector<LevelData<FArrayBox>*> RHS_h;
    Vector<LevelData<FArrayBox>*> a_head_curr;
    Vector<LevelData<FArrayBox>*> moulin_source_term;
    Vector<LevelData<FArrayBox>*> RHS_b;
    Vector<LevelData<FArrayBox>*> DENOM_b;
    Vector<LevelData<FArrayBox>*> a_ReQwIter;
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
    DENOM_b.resize(m_finest_level + 1, NULL);
    a_ReQwIter.resize(m_finest_level + 1, NULL);
    // DEBUG
    RHS_b_A.resize(m_finest_level + 1, NULL);
    RHS_b_B.resize(m_finest_level + 1, NULL);
    RHS_b_C.resize(m_finest_level + 1, NULL);
    // FACE CENTERED STUFF
    Vector<LevelData<FluxBox>*> a_Qw_ec;
    Vector<LevelData<FluxBox>*> a_Re_ec;
    Vector<LevelData<FluxBox>*> a_GapHeight_ec;
    Vector<LevelData<FluxBox>*> a_gradPw_ec;
    a_Qw_ec.resize(m_finest_level + 1, NULL);
    a_Re_ec.resize(m_finest_level + 1, NULL);
    a_GapHeight_ec.resize(m_finest_level + 1, NULL);
    a_gradPw_ec.resize(m_finest_level + 1, NULL);
    // alpha*aCoef(x)*I - beta*Div(bCoef(x)*Grad) -- note for us: alpha = 0 beta = - 1 
    Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(m_finest_level + 1);
    Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(m_finest_level + 1);

    //     Take care of GC/BC etc.
    if (m_verbosity > 3) {
        pout() <<"   ...Copy current into old & take care of ghost cells and BCs "<< endl;
    }

    for (int lev = 0; lev <= m_finest_level; lev++) {
        LevelData<FArrayBox>& oldH       = *m_old_head[lev];
        LevelData<FArrayBox>& currentH   = *m_head[lev];

        LevelData<FArrayBox>& oldB       = *m_old_gapheight[lev];
        LevelData<FArrayBox>& currentB   = *m_gapheight[lev];

        LevelData<FArrayBox>& currentRe  = *m_Re[lev];

        // handle ghost cells on the coarse-fine interface
        if (lev > 0) {
            QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                              m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                              1,  // ncomps
                              m_amrDomains[lev]);
            qcfi.coarseFineInterp(*m_head[lev], *m_head[lev-1]);
            qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
            qcfi.coarseFineInterp(*m_Re[lev], *m_Re[lev-1]);
        }

        // fill perio boundaries
        currentH.exchange();
        currentB.exchange();
        currentRe.exchange();

        // Head and b RHS
        RHS_h[lev]                = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        moulin_source_term[lev]   = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        RHS_b[lev]                = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        DENOM_b[lev]              = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
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
        a_gradPw_ec[lev]     = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);

        // Get the valid boxes
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        DataIterator dit                    = oldH.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            // get the validBox
            const Box& validBox = levelGrids.get(dit);

            // Fill BC ghost cells of h and b
            BCFill(currentH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);
            FixedNeumBCFill(currentB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

            // Copy curr into old -- copy ghost cells too 
            oldH[dit].copy(currentH[dit], 0, 0, 1);
            oldB[dit].copy(currentB[dit], 0, 0, 1);
        }
        ExtrapGhostCells( currentRe, m_amrDomains[lev]);

        // Keep DEBUG for now
        if (lev > 0 && (m_cur_step == 406)) {
           pout() << " Checking data in new grid after regrid in operations " << endl;
           for (dit.begin(); dit.ok(); ++dit) {
               BoxIterator bit(currentH[dit].box()); 
               for (bit.begin(); bit.ok(); ++bit) {
                   IntVect iv = bit(); 
                   pout() << iv << " h: " << currentH[dit](iv,0) 
                                << ",b: " << currentB[dit](iv,0) << endl;
               }
           }
        }
    } // there. We should start with consistent b and h, GC BC and all ...


    /* II LOOP: h and b calc */
    if (m_verbosity > 3) {
        pout() <<"   ...Solve for h (update b too) ! "<< endl;
    }
    bool converged_h = false;
    int ite_idx = 0;
    m_cur_PicardIte = 0;
    while (!converged_h) { 
        if (m_verbosity > 3) {
            pout() <<"   ------------------------------------- "<< endl;
            pout() <<"     Iteration "<< ite_idx << endl;
            pout() <<"   ------------------------------------- "<< endl;
        }
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
                QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                  m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                  1,  // ncomps
                                  m_amrDomains[lev]);
                qcfi.coarseFineInterp(*m_head[lev], *m_head[lev-1]);
                qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
                qcfi.coarseFineInterp(*m_Re[lev], *m_Re[lev-1]);
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
                BCFill(levelcurH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);
                FixedNeumBCFill(levelcurB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

                levelnewH_lag[dit].copy(levelcurH[dit], 0, 0, 1); // should copy ghost cells too !
                levelnewB_lag[dit].copy(levelcurB[dit], 0, 0, 1); // should copy ghost cells too !
            }
            ExtrapGhostCells( levelcurRe, m_amrDomains[lev]);

            // Interpolate b to edges
            CellToEdge(levelcurB, levelnewB_ec);
        }  // loop on levs -- same thing, we should start with consistent b and h GC/BC and all


        //         Update water pressure Pw=f(h)
        if (m_verbosity > 3) {
            pout() <<"        Update water pressure "<< endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelcurrentH = *m_head[lev];

            LevelData<FArrayBox>& levelPw       = *m_Pw[lev];
            LevelData<FArrayBox>& levelzBed     = *m_bedelevation[lev];
           
            DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
            DataIterator dit                    = levelGrids.dataIterator();
            // Update water pressure Pw=f(h)
            for (dit.begin(); dit.ok(); ++dit) {
                FArrayBox& Pw       = levelPw[dit];
                FArrayBox& currH    = levelcurrentH[dit];
                FArrayBox& zbed     = levelzBed[dit];

                BoxIterator bit(Pw.box());
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit();
                    Pw(iv,0) = (currH(iv,0) - zbed(iv,0)) * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity;
                }
            }
        } // end loop on levs

        //         Compute grad(h) and grad(Pw)
        if (m_verbosity > 3) {
            pout() <<"        Compute grad(h) and grad(Pw) "<< endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelcurrentH = *m_head[lev];
            LevelData<FArrayBox>& levelgradH    = *m_gradhead[lev];

            LevelData<FArrayBox>& levelPw       = *m_Pw[lev];
           
            // EC quantities
            LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];
            LevelData<FluxBox>& levelgradPw_ec  = *a_gradPw_ec[lev];

            DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
            DataIterator dit                    = levelGrids.dataIterator();

            // Compute grad(h) -EC and CC- 
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
            // levelgradH is just for debug purposes
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

            // Compute grad(Pw) -EC-
            crsePsiPtr = NULL;
            finePsiPtr = NULL;
            nRefCrse=-1;
            nRefFine=-1;
            if (lev > 0) {
                crsePsiPtr = m_Pw[lev-1];
                nRefCrse = m_refinement_ratios[lev-1];
            }
            if (lev < m_finest_level) {
                finePsiPtr = m_Pw[lev+1];  // What does it do with the fine stuff ???
                nRefFine = m_refinement_ratios[lev];
            }
            // EC version
            Gradient::compGradientMAC(levelgradPw_ec, levelPw,
                                     crsePsiPtr, finePsiPtr,
                                     dx, nRefCrse, nRefFine,
                                     m_amrDomains[lev]);
        } // end loop levels


        //     IV INNER LOOP: Re/Qw !! (fixed nb of iter for now)
        //             Update VECTOR Qw = f(Re, grad(h))
        //             Update Re = f(Qw)
        if (m_verbosity > 3) {
            pout() <<"        Re/Qw dependency "<< endl;
        }
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
                    QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                      m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                      2,  // num comps
                                      m_amrDomains[lev]);
                    qcfi.coarseFineInterp(*m_qw[lev], *m_qw[lev-1]);
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
                    QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                      m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                      1,  // num comps
                                      m_amrDomains[lev]);
                    qcfi.coarseFineInterp(*m_Re[lev], *m_Re[lev-1]);
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


        //         Update melting rate = f(Qw, grad(h), grad(Pw))
        if (m_verbosity > 3) {
            pout() <<"        Update melting rate "<< endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {

            LevelData<FArrayBox>& levelmR       = *m_meltRate[lev];

            // EC quantities
            LevelData<FluxBox>& levelQw_ec      = *a_Qw_ec[lev]; 
            LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];
            LevelData<FluxBox>& levelgradPw_ec  = *a_gradPw_ec[lev];

            DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
            DataIterator dit                    = levelGrids.dataIterator();

            // tmp arrays for meltingRate computation
            LevelData<FluxBox>   leveltmp_ec(levelGrids, 1, IntVect::Zero);
            LevelData<FArrayBox> leveltmp_cc(levelGrids, 1*SpaceDim, HeadGhostVect);

            for (dit.begin(); dit.ok(); ++dit) {
                // EC quantities
                FluxBox& Qwater_ec = levelQw_ec[dit];
                FluxBox& gradH_ec  = levelgradH_ec[dit];
                FluxBox& gradPw_ec = levelgradPw_ec[dit];
                FluxBox& tmp_ec    = leveltmp_ec[dit];
                // loop over directions
                for (int dir = 0; dir<SpaceDim; dir++) {
                    FArrayBox& Qwater_ecFab = Qwater_ec[dir];
                    FArrayBox& gradH_ecFab  = gradH_ec[dir];
                    FArrayBox& gradPw_ecFab = gradPw_ec[dir];
                    FArrayBox& tmp_ecFab    = tmp_ec[dir];

                    BoxIterator bitEC(Qwater_ecFab.box()); // can use gridBox? 

                    for (bitEC.begin(); bitEC.ok(); ++bitEC) {
                        IntVect iv = bitEC();
                        tmp_ecFab(iv, 0)  = - m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * (
                                              Qwater_ecFab(iv, 0) * gradH_ecFab(iv, 0) );
                        tmp_ecFab(iv, 0)  -=  m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * (
                                              Qwater_ecFab(iv, 0) * gradPw_ecFab(iv, 0) );
                    }
                } // loop over dir
            }
            EdgeToCell(leveltmp_ec, leveltmp_cc);
            // handle ghost cells on the coarse-fine interface
            //if (lev > 0) {
            //    QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
            //                      m_amrDx[lev][0], m_refinement_ratios[lev-1],  
            //                      2,  // num comps
            //                      m_amrDomains[lev]);
            //    qcfi.coarseFineInterp(leveltmp_ec[lev], leveltmp_ec[lev-1]);
            //}
            //leveltmp_cc.exchange();
            //ExtrapGhostCells( leveltmp_cc, m_amrDomains[lev]);

            for (dit.begin(); dit.ok(); ++dit) {
                // CC
                FArrayBox& meltR   = levelmR[dit];
                FArrayBox& tmp_cc  = leveltmp_cc[dit];
                BoxIterator bit(meltR.box());
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit();
                    //meltR(iv, 0) += term in ub and stress  <-- TODO
                    meltR(iv, 0)  = (m_suhmoParm->m_G + tmp_cc(iv, 0) + tmp_cc(iv, 1)) / m_suhmoParm->m_L;
                }
            }
            // handle ghost cells on the coarse-fine interface
            if (lev > 0) {
                QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                  m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                  1,  // num comps
                                  m_amrDomains[lev]);
                qcfi.coarseFineInterp(*m_meltRate[lev], *m_meltRate[lev-1]);
            }
            // Need to fill the ghost cells of melting rate -- extrapolate on perio boundaries   
            levelmR.exchange();
            ExtrapGhostCells( levelmR, m_amrDomains[lev]);
        }// loop on levs


        //         Compute aCoeff and bCoeff, RHS and solve for h again
        Vector<DisjointBoxLayout> m_amrGrids_curr;
        Vector<ProblemDomain> m_amrDomains_curr;
        m_amrGrids_curr.resize(m_finest_level + 1);
        m_amrDomains_curr.resize(m_finest_level + 1);
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    

            LevelData<FArrayBox>& levelacoef = *aCoef[lev];
            LevelData<FArrayBox>& levelRHS_h = *RHS_h[lev];
            LevelData<FArrayBox>& levelmoulin_source_term = *moulin_source_term[lev];

            LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
            LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
            LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];

            //EC quantities
            LevelData<FluxBox>& levelbcoef  = *bCoef[lev];
            LevelData<FluxBox>& levelRe_ec  = *a_Re_ec[lev]; 
            LevelData<FluxBox>& levelB_ec   = *a_GapHeight_ec[lev];    

            // Compute aCoeff and bCoeff using updated qtites
            aCoeff_bCoeff(levelacoef, levelbcoef, levelRe_ec, levelB_ec);
            // Calc moulin source term -- only level 0, then interpolate on coarser levs
            if (lev == 0) {
                Calc_moulin_source_term(levelmoulin_source_term);
            } else {
                FineInterp interpolator(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1], m_amrDomains[lev]);
                interpolator.interpToFine(*moulin_source_term[lev], *moulin_source_term[lev - 1]);
            }
            // Form RHS for h using updated qtites
            CalcRHS_head(lev, levelRHS_h, levelPi, 
                         levelPw, levelmR, 
                         levelB, levelmoulin_source_term);

            m_amrGrids_curr[lev]   = m_amrGrids[lev];
            m_amrDomains_curr[lev] = m_amrDomains[lev];
            a_head_curr[lev]       = m_head[lev];
        } // loop on levs

        // Solve for h using updated qtites
        if (m_verbosity > 3) {
            pout() <<"        Poisson solve for h "<< endl;
        }

        if (m_use_FAS) {
            SolveForHead_nl(m_amrGrids_curr, aCoef, bCoef,
                            m_amrDomains_curr, m_refinement_ratios, coarsestDx,
                            a_head_curr, RHS_h);
        } else {
            SolveForHead(m_amrGrids_curr, aCoef, bCoef,
                     m_amrDomains_curr, m_refinement_ratios, coarsestDx,
                     a_head_curr, RHS_h);
        }

        for (int lev = 0; lev <= m_finest_level; lev++) {
            m_head[lev] = a_head_curr[lev];
        }

        /* TRY TO SOLVE FOR B IN THE LOOP */
        //     Form RHS for b -- using b of current Picard iteration
        //if (m_verbosity > 3) {
        //    pout() <<"   ...Solve for b ! "<< endl;
        //} 
        //int gh_method = 0; // 0: backward Euler, 1:...
        //for (int lev = 0; lev <= m_finest_level; lev++) {
        //    LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    

        //    LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
        //    LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
        //    LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];

        //    LevelData<FArrayBox>& levelRHS_b = *RHS_b[lev];
        //    LevelData<FArrayBox>& levelDENOM_b = *DENOM_b[lev];

        //    // DEBUG
        //    LevelData<FArrayBox>& levelRHS_b_A = *RHS_b_A[lev];
        //    LevelData<FArrayBox>& levelRHS_b_B = *RHS_b_B[lev];
        //    LevelData<FArrayBox>& levelRHS_b_C = *RHS_b_C[lev];

        //    // 2. a : Get the RHS of gh eqs:
        //    CalcRHS_gapHeight(levelRHS_b, levelPi, 
        //                      levelPw, levelmR, 
        //                      levelB, 
        //                      levelRHS_b_A, levelRHS_b_B, levelRHS_b_C);
        //    //CalcRHS_gapHeight_semiExpl(levelRHS_b, levelDENOM_b, 
        //    //                           levelPi, levelPw, levelmR, 
        //    //                           levelB, a_dt, 
        //    //                           levelRHS_b_A);    
        //}  // loop on levs

        //     Solve for b using Forward Euler simple scheme -- use OLD b here
        //if (m_verbosity > 3) {
        //    pout() <<"        Update gap height with expl Euler scheme"<< endl;
        //}
        //for (int lev = m_finest_level; lev >= 0; lev--)
        //for (int lev = 0; lev <= m_finest_level; lev++) {
        //    LevelData<FArrayBox>& leveloldB  = *m_old_gapheight[lev];    
        //    LevelData<FArrayBox>& levelnewB  = *m_gapheight[lev];    

        //    LevelData<FArrayBox>& levelRHS_b   = *RHS_b[lev];
        //    LevelData<FArrayBox>& levelDENOM_b = *DENOM_b[lev];

        //    // 2. b : update gap height
        //    DisjointBoxLayout& levelGrids    = m_amrGrids[lev];
        //    DataIterator dit = levelGrids.dataIterator();
        //    for (dit.begin(); dit.ok(); ++dit) {
 
        //        FArrayBox& oldB    = leveloldB[dit];
        //        FArrayBox& newB    = levelnewB[dit];
        //        FArrayBox& RHS     = levelRHS_b[dit];
        //        FArrayBox& DENOM   = levelDENOM_b[dit];

        //        // DO NOT TOUCH GHOST CELLS
        //        BoxIterator bit(RHS.box()); 
        //        for (bit.begin(); bit.ok(); ++bit) {
        //            IntVect iv = bit(); 
        //            newB(iv,0) = RHS(iv,0) * a_dt + oldB(iv,0);
        //            //newB(iv,0) = RHS(iv,0) / DENOM(iv,0);
        //        }
        //    }
        //}  // loop on levs
        

        /* Averaging down and fill in ghost cells */
        if (m_verbosity > 3) {
            pout() <<"   ...Average down "<< endl;
        }
        for (int lev = m_finest_level; lev > 0; lev--) {
            CoarseAverage averager(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1]);
            averager.averageToCoarse(*m_head[lev - 1], *m_head[lev]);
            //averager.averageToCoarse(*m_gapheight[lev - 1], *m_gapheight[lev]);
        }
        // handle ghost cells on the coarse-fine interface
        for (int lev = 1; lev <= m_finest_level; lev++) {
            QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                              m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                              1,  // num comps
                              m_amrDomains[lev]);
            qcfi.coarseFineInterp(*m_head[lev], *m_head[lev-1]);
            //qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
        }

        /* CONVERGENCE TESTS with LAGGED quantities */
        Real maxHead = computeMax(m_head, m_refinement_ratios, Interval(0,0), 0);
        Real maxb    = computeMax(m_gapheight, m_refinement_ratios, Interval(0,0), 0);
        if (m_verbosity > 3) {
            pout() <<"        Check for convergence of h,b (max are "<< maxHead<<" "<< maxb <<")"<<endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            LevelData<FArrayBox>& levelnewH_lag  = *a_head_lagged[lev];
            LevelData<FArrayBox>& levelcurrentH  = *m_head[lev];
            LevelData<FArrayBox>& levelnewB_lag  = *a_gapheight_lagged[lev];
            LevelData<FArrayBox>& levelcurrentB  = *m_gapheight[lev];

            DataIterator dit = levelnewH_lag.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                levelnewH_lag[dit].minus(levelcurrentH[dit], 0, 0, 1);
                levelnewH_lag[dit].divide(maxHead);
                levelnewH_lag[dit].abs();

                levelnewB_lag[dit].minus(levelcurrentB[dit], 0, 0, 1);
                levelnewB_lag[dit].divide(maxb);
                levelnewB_lag[dit].abs();
            }
        }
        Real max_resH = computeMax(a_head_lagged, m_refinement_ratios, Interval(0,0), 0);
        Real max_resB = computeMax(a_gapheight_lagged, m_refinement_ratios, Interval(0,0), 0);
        if (m_verbosity > 3) {
            pout() <<ite_idx<< "         x(h,b) "<<max_resH<<" "<<max_resB<<endl;
        }

        if (ite_idx > 500) {
            pout() <<"        does not converge (Picard iterations > 500)."<< endl;
            if (m_verbosity > 0) {
                pout() <<ite_idx<< "         x(h,b) "<<max_resH<<" "<<max_resB<<endl;
            }
            MayDay::Error("Abort");
        } else {
            //if ((max_resH < 1.0e-6) && (max_resB < 1.0e-6)) {
            if (max_resH < 1.0e-6){
                if (m_verbosity > 0) {
                    pout() <<"        converged( it = "<< ite_idx << ", x(h,b) = " <<max_resH<<" "<<max_resB<< ")."<< endl;
                }
                converged_h = true;
            }
        }
     
        /* custom plt here -- debug print */
        if (m_PrintCustom && (m_cur_step == 51)) {
            int nStuffToPlot = 17;
            Vector<std::string> vectName;
            vectName.resize(nStuffToPlot);
            vectName[0]="head";
            vectName[1]="gapHeight";
            vectName[2]="gapHeightOld";
            vectName[3]="head_residual";
            vectName[4]="gapHeight_residual";
            vectName[5]="RHS_head";
            vectName[6]="RHS_b";
            vectName[7]="RHS_b_A";
            vectName[8]="RHS_b_B";
            vectName[9]="RHS_b_C";
            vectName[10]="Qw_X";
            vectName[11]="Qw_Y";
            vectName[12]="sourceMoulin";
            vectName[13]="Pwater";
            vectName[14]="metingRate";
            vectName[15]="Reynolds";
            vectName[16]="Zbed";

            Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
            stuffToPlot.resize(nStuffToPlot);
            for (int zz = 0; zz < nStuffToPlot; zz++) {
                stuffToPlot[zz].resize(m_max_level + 1, NULL);
            }

            for (int lev = 0; lev <= m_finest_level; lev++) {
                stuffToPlot[0][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[1][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[2][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[3][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[4][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[5][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[6][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[7][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[8][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[9][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[10][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[11][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[12][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[13][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, m_num_head_ghost * IntVect::Unit);
                stuffToPlot[14][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_num_head_ghost * IntVect::Unit );
                stuffToPlot[15][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_num_head_ghost * IntVect::Unit );
                stuffToPlot[16][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1,m_num_head_ghost * IntVect::Unit );

                LevelData<FArrayBox>& levelHead      = *m_head[lev];
                LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];

                LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
                LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];

                LevelData<FArrayBox>& levelOldGap    = *m_old_gapheight[lev]; 
                LevelData<FArrayBox>& levelOldGapSTP = *stuffToPlot[2][lev];

                LevelData<FArrayBox>& levelRes       = *a_head_lagged[lev];
                LevelData<FArrayBox>& levelResSTP    = *stuffToPlot[3][lev];

                LevelData<FArrayBox>& levelResB      = *a_gapheight_lagged[lev];
                LevelData<FArrayBox>& levelResBSTP   = *stuffToPlot[4][lev];

                LevelData<FArrayBox>& levelRHS       = *RHS_h[lev];
                LevelData<FArrayBox>& levelRHSSTP    = *stuffToPlot[5][lev];

                LevelData<FArrayBox>& levelRHSB       = *RHS_b[lev];
                LevelData<FArrayBox>& levelRHSBSTP    = *stuffToPlot[6][lev];
                LevelData<FArrayBox>& levelRHSBA      = *RHS_b_A[lev];
                LevelData<FArrayBox>& levelRHSBASTP   = *stuffToPlot[7][lev];
                LevelData<FArrayBox>& levelRHSBB      = *RHS_b_B[lev];
                LevelData<FArrayBox>& levelRHSBBSTP   = *stuffToPlot[8][lev];
                LevelData<FArrayBox>& levelRHSBC      = *RHS_b_C[lev];
                LevelData<FArrayBox>& levelRHSBCSTP   = *stuffToPlot[9][lev];

                LevelData<FArrayBox>& levelQw         = *m_qw[lev];
                LevelData<FArrayBox>& levelQwXSTP     = *stuffToPlot[10][lev];
                LevelData<FArrayBox>& levelQwYSTP     = *stuffToPlot[11][lev];

                LevelData<FArrayBox>& levelSourceM    = *moulin_source_term[lev];
                LevelData<FArrayBox>& levelSourceMSTP = *stuffToPlot[12][lev];

                LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
                LevelData<FArrayBox>& levelPwSTP     = *stuffToPlot[13][lev];

                LevelData<FArrayBox>& levelMR        = *m_meltRate[lev];
                LevelData<FArrayBox>& levelMRSTP     = *stuffToPlot[14][lev];

                LevelData<FArrayBox>& levelRe        = *m_Re[lev];
                LevelData<FArrayBox>& levelReSTP     = *stuffToPlot[15][lev];

                LevelData<FArrayBox>& levelZb        = *m_bedelevation[lev];
                LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[16][lev];

                DataIterator dit = levelHead.dataIterator();
                for (dit.begin(); dit.ok(); ++dit) {
                    levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
                    levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
                    levelOldGapSTP[dit].copy(levelOldGap[dit], 0, 0, 1);
                    levelResSTP[dit].copy(levelRes[dit], 0, 0, 1);
                    levelResBSTP[dit].copy(levelResB[dit], 0, 0, 1);
                    levelRHSSTP[dit].copy(levelRHS[dit], 0, 0, 1);
                    levelRHSBSTP[dit].copy(levelRHSB[dit], 0, 0, 1);
                    levelRHSBASTP[dit].copy(levelRHSBA[dit], 0, 0, 1);
                    levelRHSBBSTP[dit].copy(levelRHSBB[dit], 0, 0, 1);
                    levelRHSBCSTP[dit].copy(levelRHSBC[dit], 0, 0, 1);
                    levelQwXSTP[dit].copy(levelQw[dit], 0, 0, 1);
                    levelQwYSTP[dit].copy(levelQw[dit], 1, 0, 1);
                    levelSourceMSTP[dit].copy(levelSourceM[dit], 0, 0, 1);
                    levelPwSTP[dit].copy(levelPw[dit], 0, 0, 1);
                    levelMRSTP[dit].copy(levelMR[dit], 0, 0, 1);
                    levelReSTP[dit].copy(levelRe[dit], 0, 0, 1);
                    levelZbSTP[dit].copy(levelZb[dit], 0, 0, 1);
                }
            } // loop on levs
            writePltCustom(nStuffToPlot, vectName, stuffToPlot, std::to_string(ite_idx));
        } // end customPlt

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
        LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    

        LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
        LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
        LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];

        LevelData<FArrayBox>& levelRHS_b = *RHS_b[lev];

        // DEBUG
        LevelData<FArrayBox>& levelRHS_b_A = *RHS_b_A[lev];
        LevelData<FArrayBox>& levelRHS_b_B = *RHS_b_B[lev];
        LevelData<FArrayBox>& levelRHS_b_C = *RHS_b_C[lev];

        // 2. a : Get the RHS of gh eqs:
        CalcRHS_gapHeight(levelRHS_b, levelPi, 
                          levelPw, levelmR, 
                          levelB, 
                          levelRHS_b_A, levelRHS_b_B, levelRHS_b_C);

        LevelData<FArrayBox>& leveloldB  = *m_old_gapheight[lev];    

        // 2. b : update gap height
        DisjointBoxLayout& levelGrids    = m_amrGrids[lev];
        DataIterator dit = levelGrids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
 
            FArrayBox& oldB    = leveloldB[dit];
            FArrayBox& newB    = levelB[dit];
            FArrayBox& RHS     = levelRHS_b[dit];

            // DO NOT TOUCH GHOST CELLS
            BoxIterator bit(RHS.box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit(); 
                newB(iv,0) = RHS(iv,0) * a_dt + oldB(iv,0);
            }
        }
    }  // loop on levs

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
        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                          m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                          1,  // num comps
                          m_amrDomains[lev]);
        qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
    }

    /* Temporal probe ... very ugly */
    //for (int lev = m_finest_level; lev >= 0; lev--) {
    //    LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 

    //    DisjointBoxLayout& levelGrids    = m_amrGrids[lev];
    //    DataIterator dit = levelGrids.dataIterator();
    //    
    //    for (dit.begin(); dit.ok(); ++dit) {
    //        BoxIterator bit(levelQw[dit].box()); 
    //        for (bit.begin(); bit.ok(); ++bit) {
    //            IntVect iv = bit(); 
    //            if (iv == IntVect::Zero) {
    //                pout() << iv << " ** TimeStep " << m_cur_step << " Time " << m_time 
    //                       <<" Temporal Qw " << levelQw[dit](iv,0) << " **"<< endl;
    //            }
    //        }
    //    }
    //}  // loop on levs

    // finally, update to new time and increment current step
    m_dt = a_dt;
    m_time += a_dt;

    // write diagnostic info
    if (m_verbosity > 0) {
        pout() << "VERBOSE: AmrHydro::timestep " << m_cur_step - 1 << " --     end time = "
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
    // DEBUG 
    a_head_lagged.resize(m_finest_level + 1, NULL);
    a_head_curr.resize(m_finest_level + 1, NULL);
    a_gapheight_lagged.resize(m_finest_level + 1, NULL);
    RHS_h.resize(m_finest_level + 1, NULL);
    moulin_source_term.resize(m_finest_level + 1, NULL);
    RHS_b.resize(m_finest_level + 1, NULL);
    a_ReQwIter.resize(m_finest_level + 1, NULL);
    // FACE CENTERED STUFF
    Vector<LevelData<FluxBox>*> a_Qw_ec;
    Vector<LevelData<FluxBox>*> a_Re_ec;
    Vector<LevelData<FluxBox>*> a_GapHeight_ec;
    Vector<LevelData<FluxBox>*> a_gradZb_ec;
    a_Qw_ec.resize(m_finest_level + 1, NULL);
    a_Re_ec.resize(m_finest_level + 1, NULL);
    a_GapHeight_ec.resize(m_finest_level + 1, NULL);
    a_gradZb_ec.resize(m_finest_level + 1, NULL);
    // alpha*aCoef(x)*I - beta*Div(bCoef(x)*Grad) -- note for us: alpha = 0 beta = - 1 
    Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(m_finest_level + 1);
    Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(m_finest_level + 1);

    //     Take care of GC/BC etc.
    if (m_verbosity > 3) {
        pout() <<"   ...Copy current into old & take care of ghost cells and BCs "<< endl;
    }

    for (int lev = 0; lev <= m_finest_level; lev++) {
        LevelData<FArrayBox>& oldH       = *m_old_head[lev];
        LevelData<FArrayBox>& currentH   = *m_head[lev];

        LevelData<FArrayBox>& oldB       = *m_old_gapheight[lev];
        LevelData<FArrayBox>& currentB   = *m_gapheight[lev];

        LevelData<FArrayBox>& currentRe  = *m_Re[lev];

        // handle ghost cells on the coarse-fine interface
        if (lev > 0) {
            QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                              m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                              1,  // ncomps
                              m_amrDomains[lev]);
            qcfi.coarseFineInterp(*m_head[lev], *m_head[lev-1]);
            qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
            qcfi.coarseFineInterp(*m_Re[lev], *m_Re[lev-1]);
        }

        // fill perio boundaries
        currentH.exchange();
        currentB.exchange();
        currentRe.exchange();

        // Head and b RHS
        RHS_h[lev]                = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        moulin_source_term[lev]   = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        RHS_b[lev]                = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        // Head and B lagged for iterations
        a_head_lagged[lev]      = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        a_head_curr[lev]        = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        a_gapheight_lagged[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // Re/Qw iterations -- testing
        a_ReQwIter[lev]         = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // Stuff for OpLin
        aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero));
        bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero));
        // Face centered stuff
        a_Qw_ec[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_Re_ec[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_GapHeight_ec[lev]  = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_gradZb_ec[lev]     = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);

        // Get the valid boxes
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        DataIterator dit                    = oldH.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            // get the validBox
            const Box& validBox = levelGrids.get(dit);

            // Fill BC ghost cells of h and b
            BCFill(currentH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);
            FixedNeumBCFill(currentB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

            // Copy curr into old -- copy ghost cells too 
            oldH[dit].copy(currentH[dit], 0, 0, 1);
            oldB[dit].copy(currentB[dit], 0, 0, 1);
        }
        ExtrapGhostCells( currentRe, m_amrDomains[lev]);

    } // there. We should start with consistent b and h, GC BC and all ...

        /* custom plt here -- debug print */
        if (m_PrintCustom && (m_cur_step == 51)) {
            int nStuffToPlot = 12;
            Vector<std::string> vectName;
            vectName.resize(nStuffToPlot);
            vectName[0]="head";
            vectName[1]="gapHeight";
            vectName[2]="head_residual";
            vectName[3]="RHS_head";
            vectName[4]="RHS_b";
            vectName[5]="Qw_X";
            vectName[6]="Qw_Y";
            vectName[7]="sourceMoulin";
            vectName[8]="Pwater";
            vectName[9]="metingRate";
            vectName[10]="Reynolds";
            vectName[11]="Zbed";

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
                stuffToPlot[10][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[11][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);

                LevelData<FArrayBox>& levelHead      = *m_head[lev];
                LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];

                LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
                LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];

                LevelData<FArrayBox>& levelRes       = *a_head_lagged[lev];
                LevelData<FArrayBox>& levelResSTP    = *stuffToPlot[2][lev];

                LevelData<FArrayBox>& levelRHS       = *RHS_h[lev];
                LevelData<FArrayBox>& levelRHSSTP    = *stuffToPlot[3][lev];

                LevelData<FArrayBox>& levelRHSB       = *RHS_b[lev];
                LevelData<FArrayBox>& levelRHSBSTP    = *stuffToPlot[4][lev];

                LevelData<FArrayBox>& levelQw         = *m_qw[lev];
                LevelData<FArrayBox>& levelQwXSTP     = *stuffToPlot[5][lev];
                LevelData<FArrayBox>& levelQwYSTP     = *stuffToPlot[6][lev];

                LevelData<FArrayBox>& levelSourceM    = *moulin_source_term[lev];
                LevelData<FArrayBox>& levelSourceMSTP = *stuffToPlot[7][lev];

                LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
                LevelData<FArrayBox>& levelPwSTP     = *stuffToPlot[8][lev];

                LevelData<FArrayBox>& levelMR        = *m_meltRate[lev];
                LevelData<FArrayBox>& levelMRSTP     = *stuffToPlot[9][lev];

                LevelData<FArrayBox>& levelRe        = *m_Re[lev];
                LevelData<FArrayBox>& levelReSTP     = *stuffToPlot[10][lev];

                LevelData<FArrayBox>& levelZb        = *m_bedelevation[lev];
                LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[11][lev];

                DataIterator dit = levelHead.dataIterator();
                for (dit.begin(); dit.ok(); ++dit) {
                    levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
                    levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
                    levelResSTP[dit].copy(levelRes[dit], 0, 0, 1);
                    levelRHSSTP[dit].copy(levelRHS[dit], 0, 0, 1);
                    levelRHSBSTP[dit].copy(levelRHSB[dit], 0, 0, 1);
                    levelQwXSTP[dit].copy(levelQw[dit], 0, 0, 1);
                    levelQwYSTP[dit].copy(levelQw[dit], 1, 0, 1);
                    levelSourceMSTP[dit].copy(levelSourceM[dit], 0, 0, 1);
                    levelPwSTP[dit].copy(levelPw[dit], 0, 0, 1);
                    levelMRSTP[dit].copy(levelMR[dit], 0, 0, 1);
                    levelReSTP[dit].copy(levelRe[dit], 0, 0, 1);
                    levelZbSTP[dit].copy(levelZb[dit], 0, 0, 1);
                }
            } // loop on levs
            writePltCustom(nStuffToPlot, vectName, stuffToPlot, std::to_string(-1));
        } // end customPlt


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
                QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                  m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                  1,  // ncomps
                                  m_amrDomains[lev]);
                qcfi.coarseFineInterp(*m_head[lev], *m_head[lev-1]);
                qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
                qcfi.coarseFineInterp(*m_Re[lev], *m_Re[lev-1]);
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
                BCFill(levelcurH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);
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
            LevelData<FArrayBox>& levelcurrentH = *m_head[lev];
            LevelData<FArrayBox>& levelgradH    = *m_gradhead[lev];

            LevelData<FArrayBox>& levelZb       = *m_bedelevation[lev];
           
            // EC quantities
            LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];
            LevelData<FluxBox>& levelgradZb_ec  = *a_gradZb_ec[lev];

            DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
            DataIterator dit                    = levelGrids.dataIterator();

            // Compute grad(h) -EC and CC- 
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

            // Compute grad(Zb) -EC-
            crsePsiPtr = NULL;
            finePsiPtr = NULL;
            nRefCrse=-1;
            nRefFine=-1;
            if (lev > 0) {
                crsePsiPtr = m_bedelevation[lev-1];
                nRefCrse = m_refinement_ratios[lev-1];
            }
            if (lev < m_finest_level) {
                finePsiPtr = m_bedelevation[lev+1];  // What does it do with the fine stuff ???
                nRefFine = m_refinement_ratios[lev];
            }
            // EC version
            Gradient::compGradientMAC(levelgradZb_ec, levelZb,
                                     crsePsiPtr, finePsiPtr,
                                     dx, nRefCrse, nRefFine,
                                     m_amrDomains[lev]);
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
                        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                          m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                          2,  // num comps
                                          m_amrDomains[lev]);
                        qcfi.coarseFineInterp(*m_qw[lev], *m_qw[lev-1]);
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
                        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                          m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                          1,  // num comps
                                          m_amrDomains[lev]);
                        qcfi.coarseFineInterp(*m_Re[lev], *m_Re[lev-1]);
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
                LevelData<FArrayBox>& levelRe       = *m_Re[lev]; 
                LevelData<FArrayBox>& levelgradH    = *m_gradhead[lev];
                LevelData<FArrayBox>& levelB        = *m_gapheight[lev];    
                // For Qw
                LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 
               
                DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
                DataIterator dit                    = levelGrids.dataIterator();

                // Compute Re at CC
                for (dit.begin(); dit.ok(); ++dit) {
                    FArrayBox& Re      = levelRe[dit];
                    FArrayBox& GradHcc = levelgradH[dit];
                    FArrayBox& B       = levelB[dit];

                    BoxIterator bit(Re.box()); // can use gridBox? 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit();
                        // Update Re using GradH at CC
                        Real sqrt_gradH_cc = std::sqrt(GradHcc(iv, 0) * GradHcc(iv, 0) + GradHcc(iv, 1) * GradHcc(iv, 1));
                        Real discr = 1.0 + 4.0 * m_suhmoParm->m_omega * (
                                     std::pow(B(iv, 0), 3) * m_suhmoParm->m_gravity * sqrt_gradH_cc) / (
                                     12.0 * m_suhmoParm->m_nu * m_suhmoParm->m_nu);  
                        Re(iv, 0) = (- 1.0 + std::sqrt(discr)) / (2.0 * m_suhmoParm->m_omega) ; 
                    }
                }

                // EC quantities
                LevelData<FluxBox>& levelRe_ec      = *a_Re_ec[lev];
                LevelData<FluxBox>& levelQw_ec      = *a_Qw_ec[lev]; 
                LevelData<FluxBox>& levelB_ec       = *a_GapHeight_ec[lev];    
                LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];

                CellToEdge(levelRe, levelRe_ec);
 
                for (dit.begin(); dit.ok(); ++dit) {
                    // EC quantities
                    FluxBox& Re_ec     = levelRe_ec[dit];
                    FluxBox& Qwater_ec = levelQw_ec[dit];
                    FluxBox& currB_ec  = levelB_ec[dit];
                    FluxBox& gradH_ec  = levelgradH_ec[dit];

                    // loop over directions
                    for (int dir = 0; dir<SpaceDim; dir++) {
                        FArrayBox& Re_ecFab     = Re_ec[dir];
                        FArrayBox& Qwater_ecFab = Qwater_ec[dir];
                        FArrayBox& currB_ecFab  = currB_ec[dir];
                        FArrayBox& gradH_ecFab  = gradH_ec[dir];

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
                    QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                                      m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                                      2,  // num comps
                                      m_amrDomains[lev]);
                    qcfi.coarseFineInterp(*m_qw[lev], *m_qw[lev-1]);
                }
                // Need to fill the ghost cells -- extrapolate on boundaries   
                levelQw.exchange();
                ExtrapGhostCells( levelQw, m_amrDomains[lev]);

            } // end loop on levs
        }
        // Just in case the quad resolution gives funky results
        Real minRe = computeMin(m_Re, m_refinement_ratios, Interval(0,0), 0);
        if (minRe < 0.0) {
            pout() <<"        Re min is NEG !! abort ... "<< endl;
            MayDay::Error("Abort");
        }


        //     Update  RHS = f(Qw, grad(zb)) 
        //     compute mR 
        //     Compute lagged adv term = f((Qw, grad(h))
        if (m_verbosity > 3) {
            pout() <<"        Update RHS(h) and adv term "<< endl;
        }
        for (int lev = 0; lev <= m_finest_level; lev++) {
            // CC quantities
            LevelData<FArrayBox>& levelB         = *m_gapheight[lev];    
            LevelData<FArrayBox>& levelRHSh      = *RHS_h[lev];
            LevelData<FArrayBox>& levelmR        = *m_meltRate[lev];
            LevelData<FArrayBox>& levelmoulin_source_term = *moulin_source_term[lev];

            // EC quantities
            LevelData<FluxBox>&   levelQw_ec     = *a_Qw_ec[lev]; 
            LevelData<FluxBox>&   levelgradH_ec  = *m_gradhead_ec[lev];
            LevelData<FluxBox>&   levelgradZb_ec = *a_gradZb_ec[lev];

            DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
            DataIterator dit                    = levelGrids.dataIterator();

            // tmp arrays for vect operations
            LevelData<FluxBox>   leveltmp_ec(levelGrids, 1, IntVect::Zero);
            LevelData<FArrayBox> leveltmp_cc(levelGrids, 1*SpaceDim, HeadGhostVect);
            LevelData<FluxBox>   leveltmp2_ec(levelGrids, 1, IntVect::Zero);
            LevelData<FArrayBox> leveltmp2_cc(levelGrids, 1*SpaceDim, HeadGhostVect);

            if (lev == 0) {
                Calc_moulin_source_term(levelmoulin_source_term);
            } else {
                FineInterp interpolator(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1], m_amrDomains[lev]);
                interpolator.interpToFine(*moulin_source_term[lev], *moulin_source_term[lev - 1]);
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
                    FArrayBox& Qwater_ecFab = Qwater_ec[dir];
                    FArrayBox& gradH_ecFab  = gradH_ec[dir];
                    FArrayBox& tmp_ecFab    = tmp_ec[dir];
                    FArrayBox& gradZb_ecFab = gradZb_ec[dir];
                    FArrayBox& tmp2_ecFab   = tmp2_ec[dir];

                    BoxIterator bitEC(Qwater_ecFab.box()); 

                    for (bitEC.begin(); bitEC.ok(); ++bitEC) {
                        IntVect iv = bitEC();
                        tmp_ecFab(iv, 0)  = Qwater_ecFab(iv, 0) * gradH_ecFab(iv, 0);
                        tmp2_ecFab(iv, 0) = Qwater_ecFab(iv, 0) * gradZb_ecFab(iv, 0) ;
                    }
                } // loop over dir
            }
            // Qw gradH
            EdgeToCell(leveltmp_ec,  leveltmp_cc);
            // Qw gradZb
            EdgeToCell(leveltmp2_ec, leveltmp2_cc);

            Real rho_coef = (1.0 /  m_suhmoParm->m_rho_w - 1.0 / m_suhmoParm->m_rho_i);
            Real ub_norm = std::sqrt(  m_suhmoParm->m_ub[0]*m_suhmoParm->m_ub[0] 
                           + m_suhmoParm->m_ub[1]*m_suhmoParm->m_ub[1]) / m_suhmoParm->m_lr;
            for (dit.begin(); dit.ok(); ++dit) {
                // CC
                FArrayBox& B       = levelB[dit];
                FArrayBox& RHSh    = levelRHSh[dit];
                FArrayBox& mR      = levelmR[dit];
                FArrayBox& tmp_cc  = leveltmp_cc[dit];
                FArrayBox& tmp2_cc = leveltmp2_cc[dit];

                FArrayBox& moulinSrc = levelmoulin_source_term[dit];

                BoxIterator bit(RHSh.box());
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit();
                    // Actual RHS
                    if ( B(iv,0) < m_suhmoParm->m_br) {
                        RHSh(iv,0) = rho_coef / m_suhmoParm->m_L * (  m_suhmoParm->m_G 
                                     + m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                                      (tmp2_cc(iv, 0) + tmp2_cc(iv, 1)) )  
                                     - ub_norm * (m_suhmoParm->m_br - B(iv,0));
                    } else {
                        RHSh(iv,0) = rho_coef / m_suhmoParm->m_L * (  m_suhmoParm->m_G 
                                     + m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                                      (tmp2_cc(iv, 0) + tmp2_cc(iv, 1)) ); 
                    }
                    // Need to reconstruct mR
                    mR(iv,0)   = m_suhmoParm->m_G - 
                                 m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity *
                                 (tmp_cc(iv, 0) + tmp_cc(iv, 1) - tmp2_cc(iv, 0) - tmp2_cc(iv, 1)) -
                                 m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * (tmp_cc(iv, 0) + tmp_cc(iv, 1));
                    mR(iv,0)   = mR(iv,0) / m_suhmoParm->m_L;

                    // Adv term right now in RHS
                    RHSh(iv,0) -= rho_coef / m_suhmoParm->m_L * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * (1.0 
                                  + m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w) * 
                                  (tmp_cc(iv, 0) + tmp_cc(iv, 1));

                    // Add moulin 
                    RHSh(iv,0) += moulinSrc(iv,0);
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

            m_amrGrids_curr[lev]   = m_amrGrids[lev];
            m_amrDomains_curr[lev] = m_amrDomains[lev];
            a_head_curr[lev]       = m_head[lev];
        } // loop on levs

        //     Solve for h with FAS scheme
        if (m_verbosity > 3) {
            pout() <<"        Poisson solve for h "<< endl;
        }

        SolveForHead_nl(m_amrGrids_curr, aCoef, bCoef,
                        m_amrDomains_curr, m_refinement_ratios, coarsestDx,
                        a_head_curr, RHS_h);

        for (int lev = 0; lev <= m_finest_level; lev++) {
            m_head[lev] = a_head_curr[lev];
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
        if (m_verbosity > 3) {
            pout() <<"        Check for convergence of h (max is "<< maxHead<<" )"<<endl;
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
        if (m_PrintCustom && (m_cur_step == 51)) {
            int nStuffToPlot = 12;
            Vector<std::string> vectName;
            vectName.resize(nStuffToPlot);
            vectName[0]="head";
            vectName[1]="gapHeight";
            vectName[2]="head_residual";
            vectName[3]="RHS_head";
            vectName[4]="RHS_b";
            vectName[5]="Qw_X";
            vectName[6]="Qw_Y";
            vectName[7]="sourceMoulin";
            vectName[8]="Pwater";
            vectName[9]="metingRate";
            vectName[10]="Reynolds";
            vectName[11]="Zbed";

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
                stuffToPlot[10][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
                stuffToPlot[11][lev]  = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);

                LevelData<FArrayBox>& levelHead      = *m_head[lev];
                LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];

                LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
                LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];

                LevelData<FArrayBox>& levelRes       = *a_head_lagged[lev];
                LevelData<FArrayBox>& levelResSTP    = *stuffToPlot[2][lev];

                LevelData<FArrayBox>& levelRHS       = *RHS_h[lev];
                LevelData<FArrayBox>& levelRHSSTP    = *stuffToPlot[3][lev];

                LevelData<FArrayBox>& levelRHSB       = *RHS_b[lev];
                LevelData<FArrayBox>& levelRHSBSTP    = *stuffToPlot[4][lev];

                LevelData<FArrayBox>& levelQw         = *m_qw[lev];
                LevelData<FArrayBox>& levelQwXSTP     = *stuffToPlot[5][lev];
                LevelData<FArrayBox>& levelQwYSTP     = *stuffToPlot[6][lev];

                LevelData<FArrayBox>& levelSourceM    = *moulin_source_term[lev];
                LevelData<FArrayBox>& levelSourceMSTP = *stuffToPlot[7][lev];

                LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
                LevelData<FArrayBox>& levelPwSTP     = *stuffToPlot[8][lev];

                LevelData<FArrayBox>& levelMR        = *m_meltRate[lev];
                LevelData<FArrayBox>& levelMRSTP     = *stuffToPlot[9][lev];

                LevelData<FArrayBox>& levelRe        = *m_Re[lev];
                LevelData<FArrayBox>& levelReSTP     = *stuffToPlot[10][lev];

                LevelData<FArrayBox>& levelZb        = *m_bedelevation[lev];
                LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[11][lev];

                DataIterator dit = levelHead.dataIterator();
                for (dit.begin(); dit.ok(); ++dit) {
                    levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
                    levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
                    levelResSTP[dit].copy(levelRes[dit], 0, 0, 1);
                    levelRHSSTP[dit].copy(levelRHS[dit], 0, 0, 1);
                    levelRHSBSTP[dit].copy(levelRHSB[dit], 0, 0, 1);
                    levelQwXSTP[dit].copy(levelQw[dit], 0, 0, 1);
                    levelQwYSTP[dit].copy(levelQw[dit], 1, 0, 1);
                    levelSourceMSTP[dit].copy(levelSourceM[dit], 0, 0, 1);
                    levelPwSTP[dit].copy(levelPw[dit], 0, 0, 1);
                    levelMRSTP[dit].copy(levelMR[dit], 0, 0, 1);
                    levelReSTP[dit].copy(levelRe[dit], 0, 0, 1);
                    levelZbSTP[dit].copy(levelZb[dit], 0, 0, 1);
                }
            } // loop on levs
            writePltCustom(nStuffToPlot, vectName, stuffToPlot, std::to_string(ite_idx));
        } // end customPlt

        if (ite_idx > 500) {
            pout() <<"        does not converge (Picard iterations > 500)."<< endl;
            if (m_verbosity > 0) {
                pout() <<ite_idx<< "         x(h) "<<max_resH<<endl;
            }
            MayDay::Error("Abort");
        } else {
            if (max_resH < 1.0e-6){
                if (m_verbosity > 0) {
                    pout() <<"        converged( it = "<< ite_idx << ", x(h) = " <<max_resH<< ")."<< endl;
                }
                converged_h = true;
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
        LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    
        LevelData<FArrayBox>& levelH     = *m_head[lev];

        LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
        LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
        LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];
        LevelData<FArrayBox>& levelZb    = *m_bedelevation[lev];

        LevelData<FArrayBox>& levelRHS_b = *RHS_b[lev];

        DisjointBoxLayout& levelGrids    = m_amrGrids[lev];
        DataIterator dit = levelGrids.dataIterator();
        // Actually compute Pw for output 
        for (dit.begin(); dit.ok(); ++dit) {
            FArrayBox& newH    = levelH[dit];
            FArrayBox& Pressw  = levelPw[dit];
            FArrayBox& zb      = levelZb[dit];

            BoxIterator bit(Pressw.box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit(); 
                Pressw(iv,0) = m_suhmoParm->m_gravity * m_suhmoParm->m_rho_w * (newH(iv,0) - zb(iv,0));
            }
        }

        // 2. a : Get the RHS of gh eqs:
        CalcRHS_gapHeightFAS(levelRHS_b, levelPi, 
                             levelPw, levelmR, 
                             levelB); 

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
                newB(iv,0) = RHS(iv,0) * a_dt + oldB(iv,0);
            }
        }
    }  // loop on levs

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
        QuadCFInterp qcfi(m_amrGrids[lev], &m_amrGrids[lev-1],
                          m_amrDx[lev][0], m_refinement_ratios[lev-1],  
                          1,  // num comps
                          m_amrDomains[lev]);
        qcfi.coarseFineInterp(*m_gapheight[lev], *m_gapheight[lev-1]);
    }

    /* CONVERGENCE TESTS with LAGGED quantities */
    Real maxGap = computeMax(m_gapheight, m_refinement_ratios, Interval(0,0), 0);
    if (m_verbosity > 3) {
        pout() <<"        Check for convergence of b (max is "<< maxGap<<" )"<<endl;
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

    // write diagnostic info
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
             
            if (m_old_head[lev] == NULL) {
                m_old_head[lev]        = new LevelData<FArrayBox>;
                m_head[lev]            = new LevelData<FArrayBox>;
                m_old_gapheight[lev]   = new LevelData<FArrayBox>;
                m_gapheight[lev]       = new LevelData<FArrayBox>;
                m_Re[lev]              = new LevelData<FArrayBox>;
                m_iceheight[lev]       = new LevelData<FArrayBox>;
                m_bedelevation[lev]    = new LevelData<FArrayBox>;
                m_overburdenpress[lev] = new LevelData<FArrayBox>;
                
                m_gradhead[lev]      = new LevelData<FArrayBox>;
                m_gradhead_ec[lev]   = new LevelData<FluxBox>;
                m_Pw[lev]            = new LevelData<FArrayBox>;
                m_qw[lev]            = new LevelData<FArrayBox>;
                m_meltRate[lev]      = new LevelData<FArrayBox>;
            }

            // HEAD
            LevelData<FArrayBox>* old_oldheadDataPtr  = m_old_head[lev];
            LevelData<FArrayBox>* old_headDataPtr     = m_head[lev];
            LevelData<FArrayBox>* new_oldheadDataPtr =
                new LevelData<FArrayBox>(newDBL, m_old_head[0]->nComp(), m_old_head[0]->ghostVect());
            LevelData<FArrayBox>* new_headDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_head[0]->nComp(), m_head[0]->ghostVect());
            // Gap Height
            LevelData<FArrayBox>* old_oldgapheightDataPtr = m_old_gapheight[lev];
            LevelData<FArrayBox>* old_gapheightDataPtr    = m_gapheight[lev];
            LevelData<FArrayBox>* new_oldgapheightDataPtr =
                new LevelData<FArrayBox>(newDBL, m_old_gapheight[0]->nComp(), m_old_gapheight[0]->ghostVect());
            LevelData<FArrayBox>* new_gapheightDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_gapheight[0]->nComp(), m_gapheight[0]->ghostVect());
            // Re
            LevelData<FArrayBox>* old_ReDataPtr = m_Re[lev];
            LevelData<FArrayBox>* new_ReDataPtr =
                new LevelData<FArrayBox>(newDBL, m_Re[0]->nComp(), m_Re[0]->ghostVect());
            // Ice Height
            LevelData<FArrayBox>* old_iceheightDataPtr = m_iceheight[lev];
            LevelData<FArrayBox>* new_iceheightDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_iceheight[0]->nComp(), m_iceheight[0]->ghostVect());
            // Bed elevation
            LevelData<FArrayBox>* old_bedelevationDataPtr = m_bedelevation[lev];
            LevelData<FArrayBox>* new_bedelevationDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_bedelevation[0]->nComp(), m_bedelevation[0]->ghostVect());
            // Pi
            LevelData<FArrayBox>* old_overburdenpressDataPtr = m_overburdenpress[lev];
            LevelData<FArrayBox>* new_overburdenpressDataPtr    =
                new LevelData<FArrayBox>(newDBL, m_overburdenpress[0]->nComp(), m_overburdenpress[0]->ghostVect());

            // Other vars: grad / Pw / Qw / mR
            LevelData<FArrayBox>* old_gradheadDataPtr = m_gradhead[lev];
            LevelData<FArrayBox>* new_gradheadDataPtr = 
                new LevelData<FArrayBox>(newDBL, m_gradhead[0]->nComp(), m_gradhead[0]->ghostVect());
            //LevelData<FluxBox>* old_gradhead_ecDataPtr = m_gradhead_ec[lev];
            LevelData<FluxBox>* new_gradhead_ecDataPtr =
                new LevelData<FluxBox>(newDBL, m_gradhead_ec[0]->nComp(), IntVect::Zero);
            LevelData<FArrayBox>* old_PwDataPtr = m_Pw[lev];
            LevelData<FArrayBox>* new_PwDataPtr =
                new LevelData<FArrayBox>(newDBL, m_Pw[0]->nComp(), m_Pw[0]->ghostVect());
            LevelData<FArrayBox>* old_qwDataPtr = m_qw[lev];
            LevelData<FArrayBox>* new_qwDataPtr = 
                new LevelData<FArrayBox>(newDBL, m_qw[0]->nComp(), m_qw[0]->ghostVect());
            LevelData<FArrayBox>* old_meltRateDataPtr = m_meltRate[lev];
            LevelData<FArrayBox>* new_meltRateDataPtr =
                new LevelData<FArrayBox>(newDBL, m_meltRate[0]->nComp(), m_meltRate[0]->ghostVect());

            //call initData to take care of BC and initialize the levels... 
            initDataRegrid(lev,
                           *new_headDataPtr,
                           *new_gapheightDataPtr,
                           *new_PwDataPtr,
                           *new_qwDataPtr,
                           *new_ReDataPtr,
                           *new_meltRateDataPtr,
                           *new_bedelevationDataPtr,
                           *new_overburdenpressDataPtr,
                           *new_iceheightDataPtr);
            //BCDataRegrid(lev, *new_bedelevationDataPtr, *new_overburdenpressDataPtr,*new_iceheightDataPtr);

            LevelData<FArrayBox>& head  = *new_headDataPtr;
            LevelData<FArrayBox>& gH    = *new_gapheightDataPtr;
            LevelData<FArrayBox>& Pw    = *new_PwDataPtr;
            LevelData<FArrayBox>& qw    = *new_qwDataPtr;
            LevelData<FArrayBox>& Re    = *new_ReDataPtr;
            LevelData<FArrayBox>& melt  = *new_meltRateDataPtr;
            LevelData<FArrayBox>& zB    = *new_bedelevationDataPtr;
            LevelData<FArrayBox>& overP = *new_overburdenpressDataPtr;
            LevelData<FArrayBox>& iceH  = *new_iceheightDataPtr;

            DataIterator dit = newDBL.dataIterator();

            // DEBUG
            if (m_verbosity > 20) {

                pout() << " Checking data after initDataRegrid " << endl;

                for (dit.begin(); dit.ok(); ++dit) {
                    BoxIterator bit(head[dit].box()); 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit(); 
                        pout() << iv << " h: " << head[dit](iv,0) 
                                     << ",b: " << gH[dit](iv,0)           
                                     << ",Pw: " << Pw[dit](iv,0)           
                                     << ",qw: " << qw[dit](iv,0)           
                                     << ",Re: " << Re[dit](iv,0)           
                                     << ",mdot: " << melt[dit](iv,0)           
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
                FineInterp interpolatorGrad(newDBL, SpaceDim, m_refinement_ratios[lev - 1], m_amrDomains[lev]);

                interpolator.m_boundary_limit_type     = 1;
                interpolatorGrad.m_boundary_limit_type = 1;

                // HEAD
                interpolator.interpToFine(*new_oldheadDataPtr, *m_old_head[lev - 1]);
                interpolator.interpToFine(*new_headDataPtr, *m_head[lev - 1]);
                // Gap Height
                interpolator.interpToFine(*new_oldgapheightDataPtr, *m_old_gapheight[lev - 1]);
                interpolator.interpToFine(*new_gapheightDataPtr, *m_gapheight[lev - 1]);
                // Re
                interpolator.interpToFine(*new_ReDataPtr, *m_Re[lev - 1]);
                // Ice Height
                //interpolator.interpToFine(*new_iceheightDataPtr, *m_iceheight[lev - 1]);
                // Bed elevation
                //interpolator.interpToFine(*new_bedelevationDataPtr, *m_bedelevation[lev - 1]);
                // Pi
                //interpolator.interpToFine(*new_overburdenpressDataPtr, *m_overburdenpress[lev - 1]);
                // Other vars: grad / Pw / Qw / mR
                interpolatorGrad.interpToFine(*new_gradheadDataPtr, *m_gradhead[lev - 1]);
                interpolator.interpToFine(*new_PwDataPtr, *m_Pw[lev - 1]);
                interpolatorGrad.interpToFine(*new_qwDataPtr, *m_qw[lev - 1]);
                interpolator.interpToFine(*new_meltRateDataPtr, *m_meltRate[lev - 1]); 
            }

            // DEBUG
            if (m_verbosity > 20) {

                pout() << " Checking data after interpToFine " << endl;

                for (dit.begin(); dit.ok(); ++dit) {
                    BoxIterator bit(head[dit].box()); 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit(); 
                        pout() << iv << " h: " << head[dit](iv,0) 
                                     << ",b: " << gH[dit](iv,0)           
                                     << ",Pw: " << Pw[dit](iv,0)           
                                     << ",qw: " << qw[dit](iv,0)           
                                     << ",Re: " << Re[dit](iv,0)           
                                     << ",mdot: " << melt[dit](iv,0)           
                                     << ",zB: " << zB[dit](iv,0)           
                                     << ",oP: " << overP[dit](iv,0)           
                                     << ",iceH: " << iceH[dit](iv,0) << endl;
                    }
                }
            } // End verbosity

            // now potentially copy old-grid data on this level into new holder
            if (old_oldheadDataPtr!= NULL) {
                if (oldDBL.isClosed()) {
                    // HEAD
                    old_oldheadDataPtr->copyTo(*new_oldheadDataPtr);
                    old_headDataPtr->copyTo(*new_headDataPtr);
                    // Gap Height
                    old_oldgapheightDataPtr->copyTo(*new_oldgapheightDataPtr);
                    old_gapheightDataPtr->copyTo(*new_gapheightDataPtr);
                    // Re
                    old_ReDataPtr->copyTo(*new_ReDataPtr);
                    // Ice Height
                    //old_iceheightDataPtr->copyTo(*new_iceheightDataPtr);
                    // Bed elevation
                    //old_bedelevationDataPtr->copyTo(*new_bedelevationDataPtr);
                    // Pi
                    //old_overburdenpressDataPtr->copyTo(*new_overburdenpressDataPtr);
                    // 
                    // Other vars: grad / Pw / Qw / mR
                    old_gradheadDataPtr->copyTo(*new_gradheadDataPtr);  
                    old_PwDataPtr->copyTo(*new_PwDataPtr); 
                    old_qwDataPtr->copyTo(*new_qwDataPtr);
                    old_meltRateDataPtr->copyTo(*new_meltRateDataPtr);
                }
                // can now delete old data
                delete old_oldheadDataPtr;
                delete old_headDataPtr;
                delete old_oldgapheightDataPtr;
                delete old_gapheightDataPtr;
                delete old_ReDataPtr;
                delete old_iceheightDataPtr;
                delete old_bedelevationDataPtr;
                delete old_overburdenpressDataPtr;
                //
                delete old_gradheadDataPtr;
                delete old_PwDataPtr;
                delete old_qwDataPtr;
                delete old_meltRateDataPtr;
            } 

            if (m_verbosity > 20) {

                pout() << " Checking data after copyTo Old->New " << endl;

                for (dit.begin(); dit.ok(); ++dit) {
                    BoxIterator bit(head[dit].box()); 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit(); 
                        pout() << iv << " h: " << head[dit](iv,0) 
                                     << ",b: " << gH[dit](iv,0)           
                                     << ",Pw: " << Pw[dit](iv,0)           
                                     << ",qw: " << qw[dit](iv,0)           
                                     << ",Re: " << Re[dit](iv,0)           
                                     << ",mdot: " << melt[dit](iv,0)           
                                     << ",zB: " << zB[dit](iv,0)           
                                     << ",oP: " << overP[dit](iv,0)           
                                     << ",iceH: " << iceH[dit](iv,0) << endl;
                    }
                }

            } // End verbosity

            // exchange is necessary to fill periodic ghost cells
            // which aren't filled by the copyTo from oldLevelH or by
            // CF interp OR for when there are several Grids
            new_oldheadDataPtr->exchange();
            new_headDataPtr->exchange();
            new_oldgapheightDataPtr->exchange();
            new_gapheightDataPtr->exchange();
            new_ReDataPtr->exchange();
            new_iceheightDataPtr->exchange();
            new_bedelevationDataPtr->exchange();
            new_overburdenpressDataPtr->exchange();
            //
            new_gradheadDataPtr->exchange();
            new_PwDataPtr->exchange();
            new_qwDataPtr->exchange();
            new_meltRateDataPtr->exchange();

            // now place new holders into multilevel arrays
            m_old_head[lev]      = new_oldheadDataPtr;
            m_head[lev]          = new_headDataPtr;
            m_old_gapheight[lev] = new_oldgapheightDataPtr;
            m_gapheight[lev]     = new_gapheightDataPtr;
            m_Re[lev]            = new_ReDataPtr;
            m_iceheight[lev]     = new_iceheightDataPtr;
            m_bedelevation[lev]  = new_bedelevationDataPtr;
            m_overburdenpress[lev] = new_overburdenpressDataPtr;
            // 
            m_gradhead[lev]     = new_gradheadDataPtr;     // this one aint filled
            m_gradhead_ec[lev]  = new_gradhead_ecDataPtr;  // this one aint filled
            m_Pw[lev]       = new_PwDataPtr;
            m_qw[lev]       = new_qwDataPtr;
            m_meltRate[lev] = new_meltRateDataPtr;


            if (m_verbosity > 20) {

                pout() << " Checking data after Exchanges " << endl;

                for (dit.begin(); dit.ok(); ++dit) {
                    BoxIterator bit(head[dit].box()); 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit(); 
                        pout() << iv << " h: " << head[dit](iv,0) 
                                     << ",b: " << gH[dit](iv,0)           
                                     << ",Pw: " << Pw[dit](iv,0)           
                                     << ",qw: " << qw[dit](iv,0)           
                                     << ",Re: " << Re[dit](iv,0)           
                                     << ",mdot: " << melt[dit](iv,0)           
                                     << ",zB: " << zB[dit](iv,0)           
                                     << ",oP: " << overP[dit](iv,0)           
                                     << ",iceH: " << iceH[dit](iv,0) << endl;
                    }
                }

            } // End verbosity

        } // end loop over currently defined levels

        // now ensure that any remaining levels are null pointers
        // (in case of de-refinement)
        for (int lev = new_finest_level + 1; lev < m_old_head.size(); lev++)
        {
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

    int top_level = a_tags.size();
    top_level = min(m_tag_cap, min(top_level - 1, m_finest_level));
    // loop over levels
    for (int lev = 0; lev <= top_level; lev++) {
        IntVectSet& levelTags = a_tags[lev];
        tagCellsLevel(levelTags, lev);
        IntVectSet& tagSubset = m_vectTagSubset[lev];
        if (tagSubset.numPts() > 0) {
            levelTags &= tagSubset;
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
AmrHydro::tagCellsLevel(IntVectSet& a_tags, int a_level)
{
    if (m_verbosity > 4) {
        pout() << "AmrHydro::tagCellsLevel " << a_level << endl;
    }

    // first stab -- don't do BC's; just do one-sided
    // stencils at box edges (hopefully good enough),
    // since doing BC's properly is somewhat expensive.
    
    if (m_tag_var == "meltingRate") { 
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
                if (fabs(phi(iv, 0)) > m_tagging_val) local_tags |= iv;
            } // end loop over cells
        }         // end loop over grids

        // now buffer tags
        local_tags.grow(m_tags_grow);
        for (int dir = 0; dir < SpaceDim; dir++) {
            if (m_tags_grow_dir[dir] > m_tags_grow) local_tags.grow(dir, std::max(0, m_tags_grow_dir[dir] - m_tags_grow));
        }
        local_tags &= m_amrDomains[a_level];
        a_tags = local_tags;
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
            if ((top_level < m_max_level) && (top_level > old_top_level) && (new_finest_level <= m_tag_cap)) {
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
                                 levelIceHeight); //, levelgradzBed);

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
AmrHydro::initDataRegrid(int a_level,
                         LevelData<FArrayBox>& a_levhead,
                         LevelData<FArrayBox>& a_levgapheight,
                         LevelData<FArrayBox>& a_levPw,
                         LevelData<FArrayBox>& a_levqw,
                         LevelData<FArrayBox>& a_levRe,
                         LevelData<FArrayBox>& a_levmR,
                         LevelData<FArrayBox>& a_levzBed,
                         LevelData<FArrayBox>& a_levPi,
                         LevelData<FArrayBox>& a_levIceHeight)
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::initDataRegrid" << endl;
        pout()<<  "  .. lev " << a_level << endl;

    }


    RealVect levelDx = m_amrDx[a_level] * RealVect::Unit;
    m_IBCPtr->initializeData(levelDx, 
                             *m_suhmoParm,       
                             a_levhead, a_levgapheight, 
                             a_levPw, a_levqw,
                             a_levRe, a_levmR,
                             a_levzBed, a_levPi,
                             a_levIceHeight);
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

/// write hdf5 plotfile to the standard location
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
    sprintf(iter_str, "%s%06d_custom", m_plot_prefix.c_str(), m_cur_step);

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


/// write hdf5 plotfile to the standard location
void
AmrHydro::writePlotFile()
{
    if (m_verbosity > 3) {
        pout() << "AmrHydro::writePlotFile" << endl;
    }

    // plot comps: head + gapHeight + bedelevation + overburdenPress
    int numPlotComps = 11;

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

/// write checkpoint file out for later restarting
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
    header.m_int["num_comps"] = 6; // H/B/Pice/Zb/Re/Hice
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
        } // end loop over levels
    }

    handle.close();
}

/// read checkpoint file for restart
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
            m_overburdenpress[lev] = new LevelData<FArrayBox>(levelDBL, nPhiComp, nGhost);

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

        } // end if this level is defined
    }     // end loop over levels

}

/// set up for restart
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
