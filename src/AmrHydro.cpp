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
#include "CoarseAverageFace.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "AMRPoissonOpF_F.H"
#include "DivergenceF_F.H"
#include "AMRUtilF_F.H"
#include "CH_HDF5.H"
#include "computeNorm.H" 
#include "MayDay.H"
#include "CONSTANTS.H"
#include "Gradient.H"
#include "ExtrapGhostCells.H"

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
                      //BoxIterator bit(ghostBoxLo);
                      //for (bit.begin(); bit.ok(); ++bit) {
                      //    IntVect iv = bit();
                      //    a_state(iv, 0) = DirichletValue(dir, Side::Lo);
                      //}
                      ghostBoxLo &= a_state.box();
                      int isign = sign(Side::Lo);
                      for (BoxIterator bit(ghostBoxLo); bit.ok(); ++bit) {
                          IntVect ivTo = bit();
                          IntVect ivClose = ivTo -   isign*BASISV(dir);
                          //IntVect ivFar   = ivTo - 2*isign*BASISV(dir);
                          Real nearVal = a_state(ivClose, 0);
                          //Real farVal  = a_state(ivFar,   0);
                          Real inhomogVal = DirichletValue(dir, Side::Lo);
                          // linear or quad
                          Real ghostVal =  2.0 * inhomogVal - nearVal;
                          //Real ghostVal =  (8.0 / 3.0) * inhomogVal + (1.0 / 3.0) * farVal - 2.0 * nearVal;
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
                      //BoxIterator bit(ghostBoxHi);
                      //for (bit.begin(); bit.ok(); ++bit) {
                      //    IntVect iv = bit();
                      //    a_state(iv, 0) = DirichletValue(dir, Side::Hi);
                      //}
                      ghostBoxHi &= a_state.box();
                      int isign = sign(Side::Hi);
                      for (BoxIterator bit(ghostBoxHi); bit.ok(); ++bit) {
                          IntVect ivTo = bit();
                          IntVect ivClose = ivTo -   isign*BASISV(dir);
                          //IntVect ivFar   = ivTo - 2*isign*BASISV(dir);
                          Real nearVal = a_state(ivClose, 0);
                          //Real farVal  = a_state(ivFar,   0);
                          Real inhomogVal = DirichletValue(dir, Side::Hi);
                          // linear or quad
                          Real ghostVal =  2.0 * inhomogVal - nearVal;
                          //Real ghostVal =  (8.0 / 3.0) * inhomogVal + (1.0 / 3.0) * farVal - 2.0 * nearVal;
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
		                     Side::Lo);
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
		                     Side::Hi);
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

void
AmrHydro::SolveForHead(
                      const Vector<DisjointBoxLayout>&               a_grids,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                      ProblemDomain& coarsestDomain,
                      Vector<int>& refRatio,
                      Real& coarsestDx,
                      Vector<LevelData<FArrayBox>*>& a_head, 
                      Vector<LevelData<FArrayBox>*>& a_RHS)
{
    VCAMRPoissonOp2Factory* poissonOpF_head = new VCAMRPoissonOp2Factory;

    //BCHolder bc(ConstDiriNeumBC(IntVect::Unit, RealVect::Unit,  IntVect::Zero, RealVect::Zero));
    poissonOpF_head->define(coarsestDomain,
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
    poissonOp.define(m_amrGrids, m_refinement_ratios, m_amrDomains, m_amrDx, opFactoryPtr, 0);
    // bottom solver ?
    BiCGStabSolver<Vector<LevelData<FArrayBox>* > > solver;
    bool homogeneousBC = false;  
    solver.define(&poissonOp, homogeneousBC); 
    solver.m_normType = 0;
    solver.m_verbosity = 4;
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
    m_tagging_val = 0.1;
    m_tagOnGradPhi = true;
    m_tagging_val = 1.0;
    m_tags_grow = 0;
    m_tags_grow_dir = IntVect::Zero;

    m_plot_prefix = "plot";
    m_plot_interval = 10000000;
    m_plot_time_interval = 1.0e+12;
    m_write_gradPhi = true;

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
        if (m_gradZb[lev] != NULL)
        {
            delete m_gradZb[lev];
            m_gradZb[lev] = NULL;
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
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::initialize" << endl;
    }

    // time stepping
    m_time = 0.0;
    m_cur_step = 0;

    ParmParse ppAmr("AmrHydro");
    // num cells in each dir
    Vector<int> ancells(SpaceDim);
    // allows for domains with lower indices which are not positive
    Vector<int> domLoIndex(SpaceDim, 0);
    // assumption is that domains are not periodic
    bool is_periodic[SpaceDim];
    for (int dir = 0; dir < SpaceDim; dir++) is_periodic[dir] = false;
    Vector<int> is_periodic_int(SpaceDim, 0);

    // max level
    ppAmr.get("maxLevel", m_max_level);
    if (m_max_level > 0)
    {
        m_refinement_ratios.resize(m_max_level, -1);
        ppAmr.getarr("ref_ratios", m_refinement_ratios, 0, m_max_level);
    }
    else
    {
        m_refinement_ratios.resize(1);
        m_refinement_ratios[0] = -1;
    }


    ppAmr.query("PrintCustom", m_PrintCustom);

    ppAmr.query("verbosity", m_verbosity);

    ppAmr.query("block_factor", m_block_factor);
    ppAmr.query("regrid_lbase", m_regrid_lbase);
    ppAmr.query("regrid_interval", m_regrid_interval);
    ppAmr.query("tagCap", m_tag_cap);
    ppAmr.query("block_factor", m_block_factor);
    ppAmr.query("max_box_size", m_max_box_size);
    m_max_base_grid_size = m_max_box_size;
    ppAmr.query("max_base_grid_size", m_max_base_grid_size);

    ppAmr.getarr("num_cells", ancells, 0, ancells.size());

    // this one doesn't have a vertical dimension
    ppAmr.queryarr("domainLoIndex", domLoIndex, 0, SpaceDim);

    ppAmr.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        is_periodic[dir] = (is_periodic_int[dir] == 1);
    }

    ppAmr.query("cfl", m_cfl);
    m_initial_cfl = m_cfl;
    ppAmr.query("initial_cfl", m_initial_cfl);

    ppAmr.query("max_dt_grow_factor", m_max_dt_grow);

    ppAmr.query("fixed_dt", m_fixed_dt);
    ppAmr.query("offsetTime", m_offsetTime);

    ppAmr.query("plot_interval", m_plot_interval);
    ppAmr.query("plot_time_interval", m_plot_time_interval);
    ppAmr.query("plot_prefix", m_plot_prefix);

    ppAmr.query("write_gradPhi", m_write_gradPhi);

    ppAmr.query("check_interval", m_check_interval);
    ppAmr.query("check_prefix", m_check_prefix);
    ppAmr.query("check_overwrite", m_check_overwrite);
    ppAmr.query("check_exit", m_check_exit);

    ppAmr.get("fill_ratio", m_fill_ratio);

    ppAmr.query("nestingRadius", m_nesting_radius);

    ppAmr.query("tags_grow", m_tags_grow);
    {
        Vector<int> tgd(SpaceDim, 0);
        ppAmr.queryarr("tags_grow_dir", tgd, 0, tgd.size());
        for (int dir = 0; dir < SpaceDim; dir++)
        {
            m_tags_grow_dir[dir] = tgd[dir];
        }
    }
    bool isThereATaggingCriterion = false;
    ppAmr.query("tag_on_phi", m_tagOnGradPhi);
    isThereATaggingCriterion |= m_tagOnGradPhi;
    // if we set this to be true, require that we also provide the threshold
    if (m_tagOnGradPhi)
    {
        ppAmr.query("tagging_val", m_tagging_val);
    }
    // abort if there isn't a tagging criterion
    if (m_max_level > 0 && !isThereATaggingCriterion)
    {
        MayDay::Error("No Tagging criterion defined");
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

    // Tagging crap, see later
    std::string tagSubsetBoxesFile = "";
    m_vectTagSubset.resize(m_max_level);
    ppAmr.query("tagSubsetBoxesFile", tagSubsetBoxesFile);
    if (tagSubsetBoxesFile != "")
    {
        if (procID() == uniqueProc(SerialTask::compute))
        {
            ifstream is(tagSubsetBoxesFile.c_str(), ios::in);
            int lineno = 1;
            if (is.fail())
            {
                pout() << "Can't open " << tagSubsetBoxesFile << std::endl;
                MayDay::Error("Cannot open refine boxes file");
            }
            for (int lev = 0; lev < m_max_level; lev++)
            {
                // allowable tokens to identify levels in tag subset file
                const char level[6] = "level";
                const char domain[7] = "domain";
                char s[6];
                is >> s;
                if (std::string(level) == std::string(s))
                {
                    int inlev;
                    is >> inlev;
                    if (inlev != lev)
                    {
                        pout() << "expected ' " << lev << "' at line " << lineno << std::endl;
                        MayDay::Error("bad input file");
                    }
                }
                else if (std::string(domain) == std::string(s))
                {
                    // basic idea here is that we read in domain box
                    // (domains must be ordered from coarse->fine)
                    // until we get to a domain box which matches ours.
                    // This lets us make a single list of subset regions
                    // which we can use for any coarsening/refining of the domain
                    const Box& levelDomainBox = m_amrDomains[lev].domainBox();
                    bool stillLooking = true;
                    while (stillLooking)
                    {
                        Box domainBox;
                        is >> domainBox;
                        if (domainBox == levelDomainBox)
                        {
                            pout() << "Found a domain matching level " << lev << endl;
                            stillLooking = false;
                        }
                        else // move on until we find our level
                        {
                            // read in info for the level we're ignoring
                            // advance to next line
                            while (is.get() != '\n')
                                ;
                            lineno++;
                            int nboxes;
                            is >> nboxes;
                            if (nboxes > 0)
                            {
                                for (int i = 0; i < nboxes; ++i)
                                {
                                    Box box;
                                    is >> box;
                                    while (is.get() != '\n')
                                        ;
                                    lineno++;
                                }
                            }
                            is >> s;
                            if (std::string(domain) != std::string(s))
                            {
                                pout() << "expected '" << domain << "' at line " << lineno << ", got " << s
                                       << std::endl;
                                MayDay::Error("bad input file");
                            }
                        }
                    }
                }
                else
                {
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
                if (nboxes > 0)
                {
                    for (int i = 0; i < nboxes; ++i)
                    {
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

                if (lev > 0)
                {
                    // add lower level's subset to this subset
                    IntVectSet crseSet(m_vectTagSubset[lev - 1]);
                    if (!crseSet.isEmpty())
                    {
                        crseSet.refine(m_refinement_ratios[lev - 1]);
                        // crseSet.nestingRegion(m_block_factor,m_amrDomains[lev]);
                        if (m_vectTagSubset[lev].isEmpty())
                        {
                            m_vectTagSubset[lev] = crseSet;
                        }
                        else
                        {
                            m_vectTagSubset[lev] &= crseSet;
                        }
                    }
                }

            } // end loop over levels

        } // end if serial compute
        for (int lev = 0; lev < m_max_level; lev++) broadcast(m_vectTagSubset[lev], uniqueProc(SerialTask::compute));
    }

    // check to see if we're using predefined grids
    bool usePredefinedGrids = false;
    std::string gridFile;
    if (ppAmr.contains("gridsFile")) {
        usePredefinedGrids = true;
        ppAmr.get("gridsFile", gridFile);
    }

    // check to see if we're restarting from a checkpoint file
    if (!ppAmr.contains("restart_file")) {
        // we are restarting from a checkpoint file
        if (m_verbosity > 3) {
            pout() << "\nAmrHydro::initialize - Initializing data containers" << endl;
        }

        //---------------------------------
        // Set AMR hierarchy vector
        //---------------------------------
        // AF: add vars here
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

        m_gradZb.resize(m_max_level + 1, NULL);

        //-------------------------------------------------
        // For each level, define a collection of FArrayBox
        //-------------------------------------------------
        for (int lev = 0; lev < m_head.size(); lev++) {
            m_old_head[lev] = new LevelData<FArrayBox>;
            m_head[lev] = new LevelData<FArrayBox>;
            
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

            m_gradZb[lev] = new LevelData<FArrayBox>;
        }

        int finest_level = -1;
        if (usePredefinedGrids) {
            setupFixedGrids(gridFile);
        } else {
            // now create  grids
            initGrids(finest_level);
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
        if (m_verbosity > 4)
        {
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

    if (m_verbosity > 3)
    {
        pout() << "Done with AmrHydro::initialize\n" << endl;
    }
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
            // dump plotfile before regridding
            if ((m_cur_step % m_plot_interval == 0) && m_plot_interval > 0) {
#ifdef CH_USE_HDF5
                writePlotFile();
#endif
            }

            if ((m_cur_step != 0) && (m_cur_step % m_regrid_interval == 0)) {
                regrid();
            }

            if (m_cur_step != 0) {
                // compute dt after regridding in case number of levels has changed
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
            timeStep(dt);

        } // end of plot_time_interval
#ifdef CH_USE_HDF5
        if (m_plot_interval >= 0) writePlotFile();
#endif
    } // end timestepping loop

    // dump out final plotfile, if appropriate
    if (m_plot_interval >= 0)
    {
#ifdef CH_USE_HDF5
        writePlotFile();
#endif
    }

    // dump out final checkpoint file, if appropriate
    if (m_check_interval >= 0)
    {
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
AmrHydro::CalcRHS_head(LevelData<FArrayBox>& levelRHS_h, 
                       LevelData<FArrayBox>& levelPi, 
                       LevelData<FArrayBox>& levelPw, 
                       LevelData<FArrayBox>& levelmR, 
                       LevelData<FArrayBox>& levelB)
{
   DataIterator dit = levelRHS_h.dataIterator();
   for (dit.begin(); dit.ok(); ++dit) {

       FArrayBox& B       = levelB[dit];
       FArrayBox& RHS     = levelRHS_h[dit];

       FArrayBox& Pressi  = levelPi[dit];
       FArrayBox& Pw      = levelPw[dit];
       FArrayBox& meltR   = levelmR[dit];
       
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
           //Real x_loc = (iv[0]+0.5)*m_amrDx[0][0];
           //Real y_loc = (iv[1]+0.5)*m_amrDx[0][1];
           if ( B(iv,0) < m_suhmoParm->m_br) {
               RHS(iv,0) -= ub_norm * (m_suhmoParm->m_br - B(iv,0));
           }
           // second term ... assume  n = 3 !!
           Real PimPw = (Pressi(iv,0) - Pw(iv,0));
           Real AbsPimPw = std::abs(PimPw);
           RHS(iv,0) += m_suhmoParm->m_A * std::pow(AbsPimPw, 2) * PimPw * B(iv,0);
           if ( (iv[0] == 32) && (iv[1] == 32) ) {
               //pout() << "Moulin location "<< iv << endl;
               RHS(iv,0) += 4.0 / (m_amrDx[0][0] * m_amrDx[0][1]);  // m3s-1 / m2
           }
       }
   }
}

void
AmrHydro::CalcRHS_gapHeight(LevelData<FArrayBox>& levelRHS_b, 
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

/* core of advance routine */
void
AmrHydro::timeStep(Real a_dt)
{
    CH_TIME("AmrHydro::timestep");

    if (m_verbosity >= 2)
    {
        pout() << "\n\n\n-- Timestep " << m_cur_step << " Advancing solution from time " << m_time << " ( " << time()
               << ")"
                  " with dt = "
               << a_dt << endl;
    }
    
    /* Sketch of loops: Fig 1 of Sommers */
    //
    // NOTE: all computation done at CC. When edge qties needed, use CC->Edge
    // I Copy new h and new b into old h and old b
    // II BIG OUTER LOOP: h and b ...
    //     III SMALL OUTER LOOP: h calc
    //         Fill GC of h
    //         Put h into h_lag (GC too)
    //         Update water pressure Pw=f(h)
    //         Compute grad(h) and grad(Pw)
    //         IV INNER LOOP: Re/Qw !! (fixed nb of iter for now)
    //             Update VECTOR Qw = f(Re, grad(h))
    //             Update Re = f(Qw)
    //         Update melting rate = f(Qw, grad(h), grad(Pw))
    //         Compute aCoeff and bCoeff_cc, RHS and solve for h again
    //         Check convergence using h and h_lag
    //          
    //     Form RHS for b
    //     Solve for b using Forward Euler simple scheme
    //  
    /* End comments */

    IntVect HeadGhostVect = m_num_head_ghost * IntVect::Unit;
    Real coarsestDx = m_amrDx[0][0];  

    /* I Copy new h and new b into old h and old b */

    // Also create and initialize tmp vectors
    Vector<LevelData<FArrayBox>*> a_head_lagged;
    Vector<LevelData<FArrayBox>*> a_gapheight_lagged;
    Vector<LevelData<FArrayBox>*> RHS_h;
    Vector<LevelData<FArrayBox>*> RHS_b;
    Vector<LevelData<FArrayBox>*> a_ReQwIter;
    a_head_lagged.resize(m_max_level + 1, NULL);
    a_gapheight_lagged.resize(m_max_level + 1, NULL);
    RHS_h.resize(m_max_level + 1, NULL);
    RHS_b.resize(m_max_level + 1, NULL);
    a_ReQwIter.resize(m_max_level + 1, NULL);
    // FACE CENTERED STUFF
    Vector<LevelData<FluxBox>*> a_Qw_ec;
    Vector<LevelData<FluxBox>*> a_Re_ec;
    Vector<LevelData<FluxBox>*> a_GapHeight_ec;
    a_Qw_ec.resize(m_max_level + 1, NULL);
    a_Re_ec.resize(m_max_level + 1, NULL);
    a_GapHeight_ec.resize(m_max_level + 1, NULL);
    // alpha*aCoef(x)*I - beta*Div(bCoef(x)*Grad) -- note for us: alpha = 0 beta = - 1 
    Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(m_max_level + 1);
    Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(m_max_level + 1);

    pout() <<"   ...Copy current into old & take care of ghost cells and BCs "<< endl;
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        LevelData<FArrayBox>& oldH       = *m_old_head[lev];
        LevelData<FArrayBox>& currentH   = *m_head[lev];

        LevelData<FArrayBox>& oldB       = *m_old_gapheight[lev];
        LevelData<FArrayBox>& currentB   = *m_gapheight[lev];

        // fill internal ghost cells
        currentH.exchange();
        currentB.exchange();

        // Head and b RHS
        RHS_h[lev]         = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        RHS_b[lev]         = new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero);
        // Head and B lagged for iterations
        a_head_lagged[lev]      = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        a_gapheight_lagged[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // Re/Qw iterations -- testing
        a_ReQwIter[lev]    = new LevelData<FArrayBox>(m_amrGrids[lev], 1, HeadGhostVect);
        // Stuff for OpLin
        aCoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_amrGrids[lev], 1, IntVect::Zero));
        bCoef[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero));
        // Face centered stuff
        a_Qw_ec[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_Re_ec[lev]         = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);
        a_GapHeight_ec[lev]  = new LevelData<FluxBox>(m_amrGrids[lev], 1, IntVect::Zero);

        // Get the valid boxes
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        DataIterator dit                    = oldH.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            // get the validBox
            const Box& validBox = levelGrids.get(dit);

            // Fill BC ghost cells of h and b
            BCFill(currentH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);
            FixedNeumBCFill(currentB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

            // Copy curr into old -- copy ghost cells too 
            oldH[dit].copy(currentH[dit], 0, 0, 1);
            oldB[dit].copy(currentB[dit], 0, 0, 1);
        }

    }

    /* II BIG OUTER LOOP: h and b ... */

    /*     III SMALL OUTER LOOP: h calc */
    pout() <<"   ...Solve for h (update b too) ! "<< endl;
    bool converged_h = false;
    int ite_idx = 0;
    while (!converged_h)
    { 
        pout() <<"   ------------------------------------- "<< endl;
        pout() <<"     Iteration "<< ite_idx << endl;
        pout() <<"   ------------------------------------- "<< endl;
        // Solve for h using lagged (iteration lagged) qtities
        //         Fill GC of h and b
        //         Put h into h_lag (GC too) and b into b_lag
        for (int lev = 0; lev <= m_finest_level; lev++)
        {
            LevelData<FArrayBox>& levelcurH      = *m_head[lev];
            LevelData<FArrayBox>& levelnewH_lag  = *a_head_lagged[lev];

            LevelData<FArrayBox>& levelcurB      = *m_gapheight[lev];
            LevelData<FArrayBox>& levelnewB_lag  = *a_gapheight_lagged[lev];
            LevelData<FluxBox>& levelnewB_ec     = *a_GapHeight_ec[lev];

            // Get the valid boxes
            const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

            // fill internal ghost cells
            levelcurH.exchange();
            levelcurB.exchange();

            // Put h into h_lag
            DataIterator dit = levelnewH_lag.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                // get the validBox & fill BC ghost cells
                const Box& validBox = levelGrids.get(dit);
                BCFill(levelcurH[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);
                FixedNeumBCFill(levelcurB[dit], validBox, m_amrDomains[lev], m_amrDx[lev][0]);

                levelnewH_lag[dit].copy(levelcurH[dit], 0, 0, 1); // should copy ghost cells too !
                levelnewB_lag[dit].copy(levelcurB[dit], 0, 0, 1); // should copy ghost cells too !
            }

            // Interpolate b to edges
            CellToEdge(levelcurB, levelnewB_ec);
        }  // loop on levs

        //         Update water pressure Pw=f(h)
        //         Compute grad(h) and grad(Pw)
        //         IV INNER LOOP: Re/Qw !!
        //             Update VECTOR Qw = f(Re, grad(h))
        //             Update Re = f(Qw)
        //         Update melting rate = f(Qw, grad(h), grad(Pw))
        for (int lev = 0; lev <= m_finest_level; lev++)
        {
            LevelData<FArrayBox>& levelB        = *m_gapheight[lev];    

            LevelData<FArrayBox>& levelcurrentH = *m_head[lev];
            LevelData<FArrayBox>& levelgradH    = *m_gradhead[lev];

            LevelData<FArrayBox>& levelmR       = *m_meltRate[lev];
            LevelData<FArrayBox>& levelPw       = *m_Pw[lev];
            LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 
            LevelData<FArrayBox>& levelRe       = *m_Re[lev]; 
            LevelData<FArrayBox>& levelzBed     = *m_bedelevation[lev];

           
            // EC quantities
            LevelData<FluxBox>& levelB_ec       = *a_GapHeight_ec[lev];    
            LevelData<FluxBox>& levelQw_ec      = *a_Qw_ec[lev]; 
            LevelData<FluxBox>& levelgradH_ec   = *m_gradhead_ec[lev];
            LevelData<FluxBox>& levelRe_ec      = *a_Re_ec[lev];

            // tmp holder
            LevelData<FArrayBox>& levelReQwIter = *a_ReQwIter[lev];

            DisjointBoxLayout& levelGrids       = m_amrGrids[lev];
            DataIterator dit                    = levelGrids.dataIterator();

            // Update water pressure Pw=f(h)
            pout() <<"        Update water pressure "<< endl;
            for (dit.begin(); dit.ok(); ++dit) {
                FArrayBox& Pw       = levelPw[dit];
                FArrayBox& currH    = levelcurrentH[dit];
                FArrayBox& zbed     = levelzBed[dit];

                BoxIterator bit(Pw.box());
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit();
                    Pw(iv,0) = (currH(iv,0) - zbed(iv,0)) * m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity;
                }
                // Other solution
                //Pw.copy(currH, 0, 0, 1);
                //Pw.minus(zbed, 0, 0, 1);
                //Pw *= m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity;;
            }
        
            // Compute grad(h) -EC and CC- and grad(Pw) -EC-
            pout() <<"        Compute grad(h) and grad(Pw) "<< endl;
            // tmp holder for grad(Pw)
            LevelData<FluxBox>   levelgradPw_ec(levelGrids, 1, IntVect::Zero);

            LevelData<FArrayBox>* crsePsiPtr = NULL;
            LevelData<FArrayBox>* finePsiPtr = NULL;
            int nRefCrse=-1;
            int nRefFine=-1;
            // one level only for now ...
            //if (lev > 0) {
            //    ...
            //}
            //if (lev < m_level) {
            //    ...
            //}
            Real dx = m_amrDx[lev][0];  
            Gradient::compGradientCC(levelgradH, levelcurrentH,
                                     crsePsiPtr, finePsiPtr,
                                     dx, nRefCrse, nRefFine,
                                     m_amrDomains[lev]);
            // Need to fill the ghost cells of gradH -- extrapolate on the no perio boundaries   
            levelgradH.exchange();
            ExtrapGhostCells( levelgradH, m_amrDomains[0]);

            // EC version
            Gradient::compGradientMAC(levelgradH_ec, levelcurrentH,
                                     crsePsiPtr, finePsiPtr,
                                     dx, nRefCrse, nRefFine,
                                     m_amrDomains[lev]);

            // EC version
            Gradient::compGradientMAC(levelgradPw_ec, levelPw,
                                     crsePsiPtr, finePsiPtr,
                                     dx, nRefCrse, nRefFine,
                                     m_amrDomains[lev]);

            //         IV INNER LOOP: Re/Qw !! --> only reev Re now
            //             Update VECTOR Qw = f(Re, grad(h))
            //             Update Re = f(Qw)
            pout() <<"        Re/Qw dependency "<< endl;
            int max_ite_Re = 10; 
            Real max_Re_diff = 0.0;
            for (int it = 0; it <= max_ite_Re; it++) {
                //pout() << "           ------------------------" << endl;
                //pout() << "           ite " << it << endl;
           
                max_Re_diff = -10000;

                CellToEdge(levelRe, levelRe_ec);

                // Get Qw at EC
                for (dit.begin(); dit.ok(); ++dit) {

                    FArrayBox& Qwater  = levelQw[dit];
                    FArrayBox& Re      = levelRe[dit];

                    levelReQwIter[dit].copy(Re, 0, 0, 1);
                    //pout() << "BEG Max Re "<< levelRe[dit].max() << endl;

                    // EC quantities
                    FluxBox& Qwater_ec = levelQw_ec[dit];
                    FluxBox& gradH_ec  = levelgradH_ec[dit];
                    FluxBox& currB_ec  = levelB_ec[dit];
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
                    EdgeToCell(Qwater_ec, Qwater); 

                    //Get Re at CC
                    BoxIterator bit(Qwater.box()); // can use gridBox? 
                    for (bit.begin(); bit.ok(); ++bit) {
                        IntVect iv = bit();
                        // Update Re using this new Qw
                        Re(iv, 0) = std::sqrt( Qwater(iv, 0) * Qwater(iv, 0) 
                                             + Qwater(iv, 1) * Qwater(iv, 1)) / m_suhmoParm->m_nu;
                    }
                    levelReQwIter[dit].minus(Re, 0, 0, 1);
                    levelReQwIter[dit].abs();
                    max_Re_diff = std::max(max_Re_diff, levelReQwIter[dit].max());
                    //pout() << "END Max levelReQwIter, Re "<< levelReQwIter[dit].max() << " " << levelRe[dit].max() << endl;
                }
                // Need to fill the ghost cells of Re -- extrapolate on the no perio boundaries   
                levelRe.exchange();
                ExtrapGhostCells( levelRe, m_amrDomains[0]);

                //pout() << "           ------------------------" << endl;

            } // end Qw/Re ites

            //         Update melting rate = f(Qw, grad(h), grad(Pw))
            pout() <<"        Update melting rate "<< endl;
            for (dit.begin(); dit.ok(); ++dit) {
                FArrayBox& gradH   = levelgradH[dit];

                FArrayBox& meltR   = levelmR[dit];
                FArrayBox& Qwater  = levelQw[dit];

                FArrayBox& gradPw  = levelgradPw[dit];

                BoxIterator bit(meltR.box());
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit();
                    meltR(iv, 0)  = m_suhmoParm->m_G; // / m_suhmoParm->m_L;
                    //meltR(iv, 0) += term in ub and stress  <-- TODO
                    meltR(iv, 0) -= m_suhmoParm->m_rho_w * m_suhmoParm->m_gravity * (
                                    Qwater(iv, 0) * gradH(iv, 0) + 
                                    Qwater(iv, 1) * gradH(iv, 1) ); 
                    meltR(iv, 0) -=  m_suhmoParm->m_ct * m_suhmoParm->m_cw * m_suhmoParm->m_rho_w * (
                                    Qwater(iv, 0) * gradPw(iv, 0) + 
                                    Qwater(iv, 1) * gradPw(iv, 1) );
                    meltR(iv, 0) = meltR(iv, 0) / m_suhmoParm->m_L;
                }
            }

        }// loop on levs


        //         Compute aCoeff and bCoeff_cc, RHS and solve for h again
        for (int lev = 0; lev <= m_finest_level; lev++)
        {
            LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    

            LevelData<FArrayBox>& levelacoef = *aCoef[lev];
            LevelData<FluxBox>&   levelbcoef = *bCoef[lev];
            LevelData<FArrayBox>& levelRHS_h = *RHS_h[lev];

            LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
            LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
            LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];
            LevelData<FArrayBox>& levelRe    = *m_Re[lev]; 

            LevelData<FArrayBox> levelbcoef_cc(m_amrGrids[lev], 1, IntVect::Unit);

            // Compute aCoeff and bCoeff_cc using updated qtites
            aCoeff_bCoeff_CC(levelacoef, levelbcoef_cc, levelRe, levelB);
            // bCoeff_cc -> bCoeff via CC->Edge
            CellToEdge(levelbcoef_cc, levelbcoef);
            // Form RHS for h using updated qtites
            CalcRHS_head(levelRHS_h, levelPi, 
                         levelPw, levelmR, 
                         levelB);
        } // loop on levs

        // Solve for h using updated qtites
        pout() <<"        Poisson solve for h "<< endl;
        SolveForHead(m_amrGrids, aCoef, bCoef,
                     m_amrDomains[0], m_refinement_ratios, coarsestDx,
                     m_head, RHS_h);


        /* TRY TO SOLVE FOR B IN THE LOOP */
        //     Form RHS for b -- using b of current Picard iteration
        pout() <<"   ...Solve for b ! "<< endl;
        //int gh_method = 0; // 0: backward Euler, 1:...
        for (int lev = m_finest_level; lev >= 0; lev--)
        {
            LevelData<FArrayBox>& levelB     = *m_gapheight[lev];    

            LevelData<FArrayBox>& levelmR    = *m_meltRate[lev];
            LevelData<FArrayBox>& levelPw    = *m_Pw[lev];
            LevelData<FArrayBox>& levelPi    = *m_overburdenpress[lev];

            LevelData<FArrayBox>& levelRHS_b = *RHS_b[lev];

            // 2. a : Get the RHS of gh eqs:
            CalcRHS_gapHeight(levelRHS_b, levelPi, 
                              levelPw, levelmR, 
                              levelB);
        }  // loop on levs

        //     Solve for b using Forward Euler simple scheme -- use OLD b here
        pout() <<"        Update gap height with expl Euler scheme"<< endl;
        for (int lev = m_finest_level; lev >= 0; lev--)
        {
            LevelData<FArrayBox>& leveloldB  = *m_old_gapheight[lev];    
            LevelData<FArrayBox>& levelnewB  = *m_gapheight[lev];    

            LevelData<FArrayBox>& levelRHS_b = *RHS_b[lev];

            // 2. b : update gap height
            DisjointBoxLayout& levelGrids    = m_amrGrids[lev];
            DataIterator dit = levelGrids.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
 
                FArrayBox& oldB    = leveloldB[dit];
                FArrayBox& newB    = levelnewB[dit];
                FArrayBox& RHS     = levelRHS_b[dit];

                // DO NOT TOUCH GHOST CELLS
                BoxIterator bit(RHS.box()); 
                for (bit.begin(); bit.ok(); ++bit) {
                    IntVect iv = bit(); 
                    newB(iv,0) = RHS(iv,0) * a_dt + oldB(iv,0);
                    //pout() << "old B " << oldB(iv,0) <<  "new B " << newB(iv,0) << endl;
                }
            }
        }  // loop on levs


        /* CONVERGENCE TESTS with LAGGED quantities */
        Real maxHead = computeMax(m_head, m_refinement_ratios, Interval(0,0), 0);
        Real maxb    = computeMax(m_gapheight, m_refinement_ratios, Interval(0,0), 0);
        pout() <<"        Check for convergence of h,b (max are "<< maxHead<<" "<< maxb <<")"<<endl;
        for (int lev = 0; lev <= m_finest_level; lev++)
        {
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
        //Real max_res = computeNorm(a_head_lagged, m_refinement_ratios , coarsestDx, Interval(0,0), 0, 0);
        Real max_resH = computeMax(a_head_lagged, m_refinement_ratios, Interval(0,0), 0);
        Real max_resB = computeMax(a_gapheight_lagged, m_refinement_ratios, Interval(0,0), 0);
        pout() <<ite_idx<< "         x(h,b) "<<max_resH<<" "<<max_resB<<endl;

        //if (ite_idx > 1) {
            if (ite_idx > 500) {
                pout() <<"        does not converge."<< endl;
                MayDay::Error("Abort");
            } else {
                if ((max_resH < 1.0e-6) && (max_resB < 1.0e-6)) {
                    pout() <<"        converged."<< endl;
                    converged_h = true;
                }
            }
        //}
     
        // custom plt here
        // debug print
        if (m_PrintCustom) {
            Vector<std::string> vectName;
            vectName.resize(8);
            vectName[0]="head";
            vectName[1]="gapHeight";
            vectName[2]="bedelevation";
            vectName[3]="head_residual";
            vectName[4]="gapHeight_residual";
            vectName[5]="RHS_head";
            vectName[6]="RHS_b";
            vectName[7]="gradH";
            //vectName[8]="gradPw";
            Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
            stuffToPlot.resize(8);
            stuffToPlot[0].resize(m_max_level + 1, NULL);
            stuffToPlot[1].resize(m_max_level + 1, NULL);
            stuffToPlot[2].resize(m_max_level + 1, NULL);
            stuffToPlot[3].resize(m_max_level + 1, NULL);
            stuffToPlot[4].resize(m_max_level + 1, NULL);
            stuffToPlot[5].resize(m_max_level + 1, NULL);
            stuffToPlot[6].resize(m_max_level + 1, NULL);
            stuffToPlot[7].resize(m_max_level + 1, NULL);
            //stuffToPlot[8].resize(m_max_level + 1, NULL);
            stuffToPlot[0][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
            stuffToPlot[1][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
            stuffToPlot[2][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
            stuffToPlot[3][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
            stuffToPlot[4][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
            stuffToPlot[5][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, IntVect::Zero);
            stuffToPlot[6][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, IntVect::Zero);
            stuffToPlot[7][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
            //stuffToPlot[8][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, IntVect::Zero);
            for (int lev = 0; lev <= m_finest_level; lev++)
            {
                LevelData<FArrayBox>& levelHead      = *m_head[lev];
                LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];
                LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
                LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];
                LevelData<FArrayBox>& levelzBed      = *m_bedelevation[lev];
                LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[2][lev];

                LevelData<FArrayBox>& levelRes       = *a_head_lagged[lev];
                LevelData<FArrayBox>& levelResSTP    = *stuffToPlot[3][lev];
                LevelData<FArrayBox>& levelResB      = *a_gapheight_lagged[lev];
                LevelData<FArrayBox>& levelResBSTP   = *stuffToPlot[4][lev];

                LevelData<FArrayBox>& levelRHS       = *RHS_h[lev];
                LevelData<FArrayBox>& levelRHSSTP    = *stuffToPlot[5][lev];
                LevelData<FArrayBox>& levelRHSB      = *RHS_b[lev];
                LevelData<FArrayBox>& levelRHSBSTP   = *stuffToPlot[6][lev];

                LevelData<FArrayBox>& levelGradH     = *m_gradhead[lev];
                LevelData<FArrayBox>& levelgHSTP     = *stuffToPlot[7][lev];
                //LevelData<FArrayBox>& levelGradPw    = *RHS_b[lev];
                //LevelData<FArrayBox>& levelgPwSTP    = *stuffToPlot[8][lev];

                DataIterator dit = levelHead.dataIterator();
                for (dit.begin(); dit.ok(); ++dit) {

                    levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
                    levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
                    levelZbSTP[dit].copy(levelzBed[dit], 0, 0, 1);

                    levelResSTP[dit].copy(levelRes[dit], 0, 0, 1);
                    levelResBSTP[dit].copy(levelResB[dit], 0, 0, 1);

                    levelRHSSTP[dit].copy(levelRHS[dit], 0, 0, 1);
                    levelRHSBSTP[dit].copy(levelRHSB[dit], 0, 0, 1);

                    levelgHSTP[dit].copy(levelGradH[dit], 0, 0, 1);
                }
            } // loop on levs
            writePltCustom(8, vectName, stuffToPlot, std::to_string(ite_idx));
        }

        ite_idx++;
        pout() << endl;
    } // end while h solve


     // Temporal probe ...
    for (int lev = m_finest_level; lev >= 0; lev--)
    {
        LevelData<FArrayBox>& levelQw       = *m_qw[lev]; 

        DisjointBoxLayout& levelGrids    = m_amrGrids[lev];
        DataIterator dit = levelGrids.dataIterator();
        
        for (dit.begin(); dit.ok(); ++dit) {
            BoxIterator bit(levelQw[dit].box()); 
            for (bit.begin(); bit.ok(); ++bit) {
                IntVect iv = bit(); 
                if (iv == IntVect::Zero) {
                    pout() << iv << " ** TimeStep " << m_cur_step << " Time " << m_time 
                           <<" Temporal Qw " << levelQw[dit](iv,0) << " **"<< endl;
                }
            }
        }
    }  // loop on levs

    // allocate face-centered storage
    //Vector<LevelData<FluxBox>*> vectFluxes(m_head.size(), NULL);
    //for (int lev = m_finest_level; lev >= 0; lev--)
    //{
    //    const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

    //    IntVect ghostVect = IntVect::Unit;

    //    // if we're doing AMR, we'll need to average these fluxes
    //    // down to coarser levels. As things stand now,
    //    // CoarseAverageFace requires that the coarse LevelData<FluxBox>
    //    // have a ghost cell.
    //    vectFluxes[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, ghostVect);

    //    LevelData<FArrayBox>& levelOldHead = *m_old_head[lev];

    //    // ensure that ghost cells for head are filled in
    //    if (lev > 0)
    //    {
    //        int nGhost = levelOldHead.ghostVect()[0];
    //        PiecewiseLinearFillPatch headFiller(
    //            levelGrids, m_amrGrids[lev - 1], 1, m_amrDomains[lev - 1], m_refinement_ratios[lev - 1], nGhost);

    //        // since we're not subcycling, don't need to interpolate in time
    //        Real time_interp_coeff = 0.0;
    //        headFiller.fillInterp(levelOldHead, *m_old_head[lev - 1], *m_old_head[lev - 1], time_interp_coeff, 0, 0, 1);
    //    }
    //    // these are probably unnecessary...
    //    levelOldHead.exchange();
    //}

    // compute flux
    //computeFluxes(vectFluxes, a_dt);

    // compute div(F) and update solution
    //updateHead(m_head, m_old_head, vectFluxes, a_dt);

    // clean up temp storage
    //for (int lev = 0; lev <= m_finest_level; lev++)
    //{
    //    if (vectFluxes[lev] != NULL)
    //    {
    //        delete vectFluxes[lev];
    //        vectFluxes[lev] = NULL;
    //    }
    //}

    // finally, update to new time and increment current step
    m_dt = a_dt;
    m_time += a_dt;
    m_cur_step += 1;

// write diagnostic info, like sum of ice
    if (m_verbosity > 0) {
        pout() << "VERBOSE: AmrHydro::timestep " << m_cur_step << " --     end time = "
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
    }

    // debug print
    //if (m_PrintCustom) {
    //    Vector<std::string> vectName;
    //    vectName.resize(3);
    //    vectName[0]="head";
    //    vectName[1]="gapHeight";
    //    vectName[2]="bedelevation";
    //    Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot;
    //    stuffToPlot.resize(3);
    //    stuffToPlot[0].resize(m_max_level + 1, NULL);
    //    stuffToPlot[1].resize(m_max_level + 1, NULL);
    //    stuffToPlot[2].resize(m_max_level + 1, NULL);
    //    stuffToPlot[0][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
    //    stuffToPlot[1][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
    //    stuffToPlot[2][0]  = new LevelData<FArrayBox>(m_amrGrids[0], 1, m_num_head_ghost * IntVect::Unit);
    //    for (int lev = 0; lev <= m_finest_level; lev++)
    //    {
    //        LevelData<FArrayBox>& levelHead      = *m_head[lev];
    //        LevelData<FArrayBox>& levelHeadSTP   = *stuffToPlot[0][lev];
    //        LevelData<FArrayBox>& levelGap       = *m_gapheight[lev];    
    //        LevelData<FArrayBox>& levelGapSTP    = *stuffToPlot[1][lev];
    //        LevelData<FArrayBox>& levelzBed      = *m_bedelevation[lev];
    //        LevelData<FArrayBox>& levelZbSTP     = *stuffToPlot[2][lev];
    //        // Put h into h_lag
    //        DataIterator dit = levelHead.dataIterator();
    //        for (dit.begin(); dit.ok(); ++dit) {
    //            levelHeadSTP[dit].copy(levelHead[dit], 0, 0, 1);
    //            levelGapSTP[dit].copy(levelGap[dit], 0, 0, 1);
    //            levelZbSTP[dit].copy(levelzBed[dit], 0, 0, 1);
    //        }
    //    } // loop on levs
    //    writePltCustom(3, vectName, stuffToPlot, "final_");
    //}


}

// compute half-time face-centered thickness using unsplit PPM
void
AmrHydro::computeFluxes(Vector<LevelData<FluxBox>*>& a_Flux, Real a_dt)
{
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        LevelData<FluxBox>& levelFlux = *a_Flux[lev];

        DataIterator dit = m_amrGrids[lev].dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            // this is just an example, after all..,
            levelFlux[dit].setVal(1.0);

        } // end loop over grids

    } // end loop over levels for computing fluxes

    // coarse average new fluxes to covered regions
    for (int lev = m_finest_level; lev > 0; lev--)
    {
        CoarseAverageFace faceAverager(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1]);
        faceAverager.averageToCoarse(*a_Flux[lev - 1], *a_Flux[lev]);
    }
}

/* update head */
void
AmrHydro::updateHead(Vector<LevelData<FArrayBox>*>& a_new_head,
                     const Vector<LevelData<FArrayBox>*>& a_old_head,
                     const Vector<LevelData<FluxBox>*>& a_vectFluxes,
                     Real a_dt)
{
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        LevelData<FluxBox>& levelFlux = *a_vectFluxes[lev];
        LevelData<FArrayBox>& levelNewH = *a_new_head[lev];
        LevelData<FArrayBox>& levelOldH = *a_old_head[lev];

        const RealVect& dx = m_amrDx[lev];

        DataIterator dit = levelGrids.dataIterator();

        for (dit.begin(); dit.ok(); ++dit)
        {
            const Box& gridBox = levelGrids[dit];
            FArrayBox& newH    = levelNewH[dit];
            FArrayBox& oldH    = levelOldH[dit];
            FluxBox& thisFlux  = levelFlux[dit];
            newH.setVal(0.0);

            // loop over directions and increment with div(F)
            for (int dir = 0; dir < SpaceDim; dir++)
            {
                // use the divergence from
                // Chombo/example/fourthOrderMappedGrids/util/DivergenceF.ChF
                FORT_DIVERGENCE(CHF_CONST_FRA(thisFlux[dir]),
                                CHF_FRA(newH),
                                CHF_BOX(gridBox),
                                CHF_CONST_REAL(dx[dir]),
                                CHF_INT(dir));
            }

            newH *= -1 * a_dt;
            newH.plus(oldH, 0, 0, 1);

        } // end loop over grids
    }     // end loop over levels

    // average down thickness to coarser levels and fill in ghost cells
    // before calling recomputeGeometry.
    int headGhost = a_new_head[0]->ghostVect()[0];
    // average from the finest level down
    for (int lev = m_finest_level; lev > 0; lev--)
    {
        CoarseAverage averager(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1]);
        averager.averageToCoarse(*a_new_head[lev - 1], *a_new_head[lev]);
    }

    // now pass back over and do PiecewiseLinearFillPatch
    for (int lev = 1; lev <= m_finest_level; lev++)
    {
        PiecewiseLinearFillPatch filler(
            m_amrGrids[lev], m_amrGrids[lev - 1], 1, m_amrDomains[lev - 1], m_refinement_ratios[lev - 1], headGhost);

        Real interp_coef = 0;
        filler.fillInterp(*a_new_head[lev], *a_new_head[lev - 1], *a_new_head[lev - 1], interp_coef, 0, 0, 1);
    }

    // average from the finest level down
    for (int lev = m_finest_level; lev > 0; lev--)
    {
        CoarseAverage averager(m_amrGrids[lev], 1, m_refinement_ratios[lev - 1]);
        averager.averageToCoarse(*a_new_head[lev - 1], *a_new_head[lev]);
    }

    for (int lev = 1; lev <= m_finest_level; lev++)
    {
        PiecewiseLinearFillPatch filler(
            m_amrGrids[lev], m_amrGrids[lev - 1], 1, m_amrDomains[lev - 1], m_refinement_ratios[lev - 1], headGhost);

        Real interp_coef = 0;
        filler.fillInterp(*a_new_head[lev], *a_new_head[lev - 1], *a_new_head[lev - 1], interp_coef, 0, 0, 1);
    }
}

// do regridding
void
AmrHydro::regrid()
{
    CH_TIME("AmrHydro::regrid");

    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::regrid" << endl;
    }

    // only do any of this if the max level > 0
    if (m_max_level > 0)
    {
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
        for (int lev = 0; lev <= m_finest_level; lev++)
        {
            const DisjointBoxLayout& levelDBL = m_amrGrids[lev];
            old_grids[lev].resize(levelDBL.size());
            LayoutIterator lit = levelDBL.layoutIterator();
            int boxIndex = 0;
            for (lit.begin(); lit.ok(); ++lit, ++boxIndex)
            {
                old_grids[lev][boxIndex] = levelDBL[lit()];
            }
        }

        int new_finest_level;

        BRMeshRefine meshrefine(
            m_amrDomains[0], m_refinement_ratios, m_fill_ratio, m_block_factor, m_nesting_radius, m_max_box_size);

        new_finest_level = meshrefine.regrid(new_grids, tagVect, m_regrid_lbase, top_level, old_grids);

        // test to see if grids have changed
        bool gridsSame = true;
        for (int lev = m_regrid_lbase + 1; lev <= new_finest_level; ++lev)
        {
            int numGridsNew = new_grids[lev].size();
            Vector<int> procIDs(numGridsNew);
            LoadBalance(procIDs, new_grids[lev]);
            const DisjointBoxLayout newDBL(new_grids[lev], procIDs, m_amrDomains[lev]);
            const DisjointBoxLayout oldDBL = m_amrGrids[lev];
            gridsSame &= oldDBL.sameBoxes(newDBL);
        }
        if (gridsSame)
        {
            if (m_verbosity > 3)
            {
                pout() << "AmrHydro::regrid -- grids unchanged" << endl;
            }
        }

        // now loop through levels and redefine if necessary
        for (int lev = m_regrid_lbase + 1; lev <= new_finest_level; ++lev)
        {
            int numGridsNew = new_grids[lev].size();
            Vector<int> procIDs(numGridsNew);
            LoadBalance(procIDs, new_grids[lev]);

            const DisjointBoxLayout newDBL(new_grids[lev], procIDs, m_amrDomains[lev]);

            const DisjointBoxLayout oldDBL = m_amrGrids[lev];

            m_amrGrids[lev] = newDBL;

            // build new storage
            LevelData<FArrayBox>* old_oldDataPtr = m_old_head[lev];
            LevelData<FArrayBox>* old_headDataPtr = m_head[lev];

            LevelData<FArrayBox>* new_oldDataPtr =
                new LevelData<FArrayBox>(newDBL, m_old_head[0]->nComp(), m_old_head[0]->ghostVect());

            LevelData<FArrayBox>* new_headDataPtr =
                new LevelData<FArrayBox>(newDBL, m_head[0]->nComp(), m_head[0]->ghostVect());

            // first fill with interpolated data from coarser level

            {
                // may eventually want to do post-regrid smoothing on this!
                FineInterp interpolator(newDBL, 1, m_refinement_ratios[lev - 1], m_amrDomains[lev]);

                interpolator.interpToFine(*new_oldDataPtr, *m_old_head[lev - 1]);
                interpolator.interpToFine(*new_headDataPtr, *m_head[lev - 1]);
            }

            // now copy old-grid data on this level into new holder
            if (old_oldDataPtr != NULL)
            {
                if (oldDBL.isClosed())
                {
                    old_oldDataPtr->copyTo(*new_oldDataPtr);
                    old_headDataPtr->copyTo(*new_headDataPtr);
                }
                // can now delete old data
                delete old_oldDataPtr;
                delete old_headDataPtr;
            }

            // exchange is necessary to fill periodic ghost cells
            // which aren't filled by the copyTo from oldLevelH
            new_oldDataPtr->exchange();
            new_headDataPtr->exchange();

            // now place new holders into multilevel arrays
            m_old_head[lev] = new_oldDataPtr;
            m_head[lev] = new_headDataPtr;

        } // end loop over currently defined levels

        // now ensure that any remaining levels are null pointers
        // (in case of de-refinement)
        for (int lev = new_finest_level + 1; lev < m_old_head.size(); lev++)
        {
            if (m_old_head[lev] != NULL)
            {
                delete m_old_head[lev];
                m_old_head[lev] = NULL;
            }

            if (m_head[lev] != NULL)
            {
                delete m_head[lev];
                m_head[lev] = NULL;
            }

            DisjointBoxLayout emptyDBL;
            m_amrGrids[lev] = emptyDBL;
        }

        m_finest_level = new_finest_level;

        // set up counter of number of cells
        for (int lev = 0; lev <= m_max_level; lev++)
        {
            m_num_cells[lev] = 0;
            if (lev <= m_finest_level)
            {
                const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
                LayoutIterator lit = levelGrids.layoutIterator();
                for (lit.begin(); lit.ok(); ++lit)
                {
                    const Box& thisBox = levelGrids.get(lit());
                    m_num_cells[lev] += thisBox.numPts();
                }
            }
        }

        // finally, set up covered_level flags
        m_covered_level.resize(m_max_level + 1, 0);
        // note that finest level can't be covered.
        for (int lev = m_finest_level - 1; lev >= 0; lev--)
        {
            // if the next finer level is covered, then this one is too.
            if (m_covered_level[lev + 1] == 1)
            {
                m_covered_level[lev] = 1;
            }
            else
            {
                // see if the grids finer than this level completely cover it
                IntVectSet fineUncovered(m_amrDomains[lev + 1].domainBox());
                const DisjointBoxLayout& fineGrids = m_amrGrids[lev + 1];

                LayoutIterator lit = fineGrids.layoutIterator();
                for (lit.begin(); lit.ok(); ++lit)
                {
                    const Box& thisBox = fineGrids.get(lit());
                    fineUncovered.minus_box(thisBox);
                }

                if (fineUncovered.isEmpty())
                {
                    m_covered_level[lev] = 1;
                }
            }
        } // end loop over levels to determine covered levels

    } // end if max level > 0 in the first place
}

void
AmrHydro::tagCells(Vector<IntVectSet>& a_tags)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::tagCells" << endl;
    }

    int top_level = a_tags.size();
    top_level = min(m_tag_cap, min(top_level - 1, m_finest_level));
    // loop over levels
    for (int lev = 0; lev <= top_level; lev++)
    {
        IntVectSet& levelTags = a_tags[lev];
        tagCellsLevel(levelTags, lev);
        IntVectSet& tagSubset = m_vectTagSubset[lev];
        if (tagSubset.numPts() > 0)
        {
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
    if (m_verbosity > 4)
    {
        pout() << "AmrHydro::tagCellsLevel " << a_level << endl;
    }

    // base tags on undivided gradient of phi
    // first stab -- don't do BC's; just do one-sided
    // stencils at box edges (hopefully good enough),
    // since doing BC's properly is somewhat expensive.

    DataIterator dit = m_head[a_level]->dataIterator();

    LevelData<FArrayBox>& levelPhi = *m_head[a_level];

    const DisjointBoxLayout& levelGrids = m_amrGrids[a_level];

    // need to ensure that ghost cells are set properly
    levelPhi.exchange(levelPhi.interval());

    IntVectSet local_tags;
    if (m_tagOnGradPhi)
    {
        for (dit.begin(); dit.ok(); ++dit)
        {
            // note that we only need one component here
            // because the fortran subroutine stores the max(abs(grad))
            // over all components into the 0th position
            FArrayBox gradPhi(levelGrids[dit()], 1);

            for (int dir = 0; dir < SpaceDim; dir++)
            {
                const Box b = levelGrids[dit()];
                const Box bcenter = b & grow(m_amrDomains[a_level], -BASISV(dir));
                const Box blo = b & adjCellLo(bcenter, dir);
                const Box bhi = b & adjCellHi(bcenter, dir);
                const int haslo = !blo.isEmpty();
                const int hashi = !bhi.isEmpty();
                FORT_UNDIVIDEDGRAD(CHF_FRA1(gradPhi, 0),
                                   CHF_CONST_FRA(levelPhi[dit()]),
                                   CHF_BOX(bcenter),
                                   CHF_BOX(blo),
                                   CHF_BOX(bhi),
                                   CHF_CONST_INT(dir),
                                   CHF_CONST_INT(haslo),
                                   CHF_CONST_INT(hashi));

                // now tag cells based on values
                BoxIterator bit(levelGrids[dit()]);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    const IntVect& iv = bit();
                    if (fabs(gradPhi(iv, 0)) > m_tagging_val) local_tags |= iv;
                } // end loop over cells
            }     // end loop over directions
        }         // end loop over grids
    }             // end if tag on grad vel

    // tag on laplacian(velocity)
    //if (m_tagOnLapPhi)
    //{
    //    for (dit.begin(); dit.ok(); ++dit)
    //    {
    //        FArrayBox lapPhi(levelGrids[dit()], levelPhi.nComp());
    //        lapPhi.setVal(0.0);
    //        Real alpha = 0;
    //        Real beta = 1.0;

    //        // use undivided laplacian (set dx = 1)
    //        Real bogusDx = 1.0;
    //        Box lapBox = levelPhi[dit].box();
    //        lapBox.grow(-2);
    //        lapBox &= levelGrids[dit];
    //        // assumes that ghost cells boundary conditions are properly set
    //        FORT_OPERATORLAP(CHF_FRA(lapPhi),
    //                         CHF_FRA(levelPhi[dit]),
    //                         CHF_BOX(lapBox),
    //                         CHF_CONST_REAL(bogusDx),
    //                         CHF_CONST_REAL(alpha),
    //                         CHF_CONST_REAL(beta));

    //        // now tag cells based on values
    //        BoxIterator bit(lapBox);

    //        for (bit.begin(); bit.ok(); ++bit)
    //        {
    //            const IntVect& iv = bit();
    //            for (int comp = 0; comp < lapPhi.nComp(); comp++)
    //            {
    //                if (fabs(lapPhi(iv, comp)) > m_laplacian_tagging_val)
    //                {
    //                    local_tags |= iv;
    //                }
    //            }

    //        } // end loop over cells
    //    }     // end loop over grids
    //}         // end if tag on laplacian(phi)

    // now buffer tags

    local_tags.grow(m_tags_grow);
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        if (m_tags_grow_dir[dir] > m_tags_grow) local_tags.grow(dir, std::max(0, m_tags_grow_dir[dir] - m_tags_grow));
    }
    local_tags &= m_amrDomains[a_level];

    a_tags = local_tags;
}

void
AmrHydro::tagCellsInit(Vector<IntVectSet>& a_tags)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::tagCellsInit" << endl;
    }

    // default is to just call regular tagging
    tagCells(a_tags);
    m_vectTags = a_tags;
}

void
AmrHydro::initGrids(int a_finest_level)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::initGrids" << endl;
    }

    m_finest_level = 0;
    // first create base level
    Vector<Box> baseBoxes;
    domainSplit(m_amrDomains[0], baseBoxes, m_max_base_grid_size, m_block_factor);

    Vector<int> procAssign(baseBoxes.size());
    LoadBalance(procAssign, baseBoxes);

    DisjointBoxLayout baseGrids(baseBoxes, procAssign, m_amrDomains[0]);

    if (m_verbosity > 3)
    {
        long long numCells0 = baseGrids.numCells();
        pout() << "    Level 0: " << numCells0 << " cells: " << baseGrids << endl;
    }

    // Set collection of boxes hierarchy, init level 0 to base grid
    m_amrGrids.resize(m_max_level + 1);
    m_amrGrids[0] = baseGrids;

    // levelSetup now create the data container for each box in level 0
    levelSetup(0, baseGrids);

    // AF: really necessary ?
    LevelData<FArrayBox>& baseLevelVel = *m_head[0];
    DataIterator baseDit = baseGrids.dataIterator();
    for (baseDit.begin(); baseDit.ok(); ++baseDit)
    {
        // initial guess at base-level velocity is zero
        baseLevelVel[baseDit].setVal(0.0);
    }

    // initialize base level data
    initData(m_head);

    int numLevels = 1;
    bool moreLevels = (m_max_level > 0);

    int baseLevel = 0;

    BRMeshRefine meshrefine;
    if (moreLevels)
    {
        meshrefine.define(
            m_amrDomains[0], m_refinement_ratios, m_fill_ratio, m_block_factor, m_nesting_radius, m_max_box_size);
    }

    Vector<IntVectSet> tagVect(m_max_level);

    Vector<Vector<Box> > oldBoxes(1);
    Vector<Vector<Box> > newBoxes;
    oldBoxes[0] = baseBoxes;
    newBoxes = oldBoxes;
    int new_finest_level = 0;

    while (moreLevels)
    {
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
        if (tagVect[m_max_level - 1].isEmpty())
        {
            int top_level = m_finest_level;
            int old_top_level = top_level;
            new_finest_level = meshrefine.regrid(newBoxes, tagVect, baseLevel, top_level, oldBoxes);

            if (new_finest_level > top_level) top_level++;
            oldBoxes = newBoxes;

            // now see if we need another pass through grid generation
            if ((top_level < m_max_level) && (top_level > old_top_level) && (new_finest_level <= m_tag_cap))
            {
                moreLevels = true;
            }
        }
        else
        {
            // for now, define old_grids as just domains
            oldBoxes.resize(m_max_level + 1);
            for (int lev = 1; lev <= m_max_level; lev++)
            {
                oldBoxes[lev].push_back(m_amrDomains[lev].domainBox());
            }

            int top_level = m_max_level - 1;
            new_finest_level = meshrefine.regrid(newBoxes, tagVect, baseLevel, top_level, oldBoxes);
        }

        numLevels = Min(new_finest_level, m_max_level) + 1;

        // now loop through levels and define
        for (int lev = baseLevel + 1; lev <= new_finest_level; ++lev)
        {
            int numGridsNew = newBoxes[lev].size();
            Vector<int> procIDs(numGridsNew);
            LoadBalance(procIDs, newBoxes[lev]);
            const DisjointBoxLayout newDBL(newBoxes[lev], procIDs, m_amrDomains[lev]);
            m_amrGrids[lev] = newDBL;

            if (m_verbosity > 2)
            {
                long long levelNumCells = newDBL.numCells();
                pout() << "   Level " << lev << ": " << levelNumCells << " cells: " << m_amrGrids[lev] << endl;
            }

            levelSetup(lev, m_amrGrids[lev]);

        } // end loop over levels

        m_finest_level = new_finest_level;

        // finally, initialize data on final hierarchy
        // only do this if we've created new levels
        if (m_finest_level > 0)
        {
            defineSolver();

            initData(m_head);
        }
    } // end while more levels to do
}

void
AmrHydro::setupFixedGrids(const std::string& a_gridFile)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::setupFixedGrids" << endl;
    }
    Vector<Vector<Box> > gridvect;

    if (procID() == uniqueProc(SerialTask::compute))
    {
        gridvect.push_back(Vector<Box>(1, m_amrDomains[0].domainBox()));

        // read in predefined grids
        ifstream is(a_gridFile.c_str(), ios::in);

        if (is.fail())
        {
            MayDay::Error("Cannot open grids file");
        }

        // format of file:
        //   number of levels, then for each level (starting with level 1):
        //   number of grids on level, list of boxes
        int inNumLevels;
        is >> inNumLevels;

        CH_assert(inNumLevels <= m_max_level + 1);

        if (m_verbosity > 3)
        {
            pout() << "numLevels = " << inNumLevels << endl;
        }

        while (is.get() != '\n')
            ;

        gridvect.resize(inNumLevels);

        // check to see if coarsest level needs to be broken up
        domainSplit(m_amrDomains[0], gridvect[0], m_max_base_grid_size, m_block_factor);

        if (m_verbosity >= 3)
        {
            pout() << "level 0: ";
            for (int n = 0; n < gridvect[0].size(); n++)
            {
                pout() << gridvect[0][n] << endl;
            }
        }

        // now loop over levels, starting with level 1
        int numGrids = 0;
        for (int lev = 1; lev < inNumLevels; lev++)
        {
            is >> numGrids;

            if (m_verbosity >= 3)
            {
                pout() << "level " << lev << " numGrids = " << numGrids << endl;
                pout() << "Grids: ";
            }

            while (is.get() != '\n')
                ;

            gridvect[lev].resize(numGrids);

            for (int i = 0; i < numGrids; i++)
            {
                Box bx;
                is >> bx;

                while (is.get() != '\n')
                    ;

                // quick check on box size
                Box bxRef(bx);

                if (bxRef.longside() > m_max_box_size)
                {
                    pout() << "Grid " << bx << " too large" << endl;
                    MayDay::Error();
                }

                if (m_verbosity >= 3)
                {
                    pout() << bx << endl;
                }

                gridvect[lev][i] = bx;
            } // end loop over boxes on this level
        }     // end loop over levels
    }         // end if serial proc

    // broadcast results
    broadcast(gridvect, uniqueProc(SerialTask::compute));

    // now create disjointBoxLayouts and allocate grids

    m_amrGrids.resize(m_max_level + 1);

    // probably eventually want to do this differently
    RealVect dx = m_amrDx[0] * RealVect::Unit;

    for (int lev = 0; lev < gridvect.size(); lev++)
    {
        int numGridsLev = gridvect[lev].size();
        Vector<int> procIDs(numGridsLev);
        LoadBalance(procIDs, gridvect[lev]);
        const DisjointBoxLayout newDBL(gridvect[lev], procIDs, m_amrDomains[lev]);

        m_amrGrids[lev] = newDBL;

        // build storage for this level

        levelSetup(lev, m_amrGrids[lev]);
        if (lev < gridvect.size() - 1)
        {
            dx /= m_refinement_ratios[lev];
        }
    }

    // finally set finest level and initialize data on hierarchy
    m_finest_level = gridvect.size() - 1;

    initData(m_head);
}

void
AmrHydro::levelSetup(int a_level, const DisjointBoxLayout& a_grids)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::levelSetup - level " << endl;
    }
    int nPhiComp = 1;
    IntVect ghostVect = IntVect::Unit;
    IntVect HeadGhostVect = m_num_head_ghost * IntVect::Unit;
    m_old_head[a_level]->define(a_grids, nPhiComp, HeadGhostVect);
    m_old_gapheight[a_level]->define(a_grids, nPhiComp, HeadGhostVect);

    // AF: setup struct for level a_level
    if (a_level == 0 || m_head[a_level] == NULL)
    {
        // First pass or first level
        m_head[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_gradhead[a_level] = new LevelData<FArrayBox>(a_grids, SpaceDim*nPhiComp, HeadGhostVect);
        m_gradhead_ec[a_level] = new LevelData<FluxBox>(a_grids, nPhiComp, IntVect::Zero);
        m_gapheight[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_Pw[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_qw[a_level] = new LevelData<FArrayBox>(a_grids, SpaceDim*nPhiComp, HeadGhostVect);
        m_Re[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_meltRate[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_iceheight[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_bedelevation[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_overburdenpress[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_gradZb[a_level] = new LevelData<FArrayBox>(a_grids, SpaceDim*nPhiComp, HeadGhostVect);
    }
    else
    {
        // AF: not here yet but at some point to do
        // do phi a bit differently in order to use previously
        // computed velocity field as an initial guess
        {
            LevelData<FArrayBox>* newHeadPtr = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
            LevelData<FArrayBox>* newGapHeightPtr = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);

            // first do interp from coarser level
            FineInterp HeadInterp(a_grids, nPhiComp, m_refinement_ratios[a_level - 1], m_amrDomains[a_level]);
            FineInterp GapHeightInterp(a_grids, nPhiComp, m_refinement_ratios[a_level - 1], m_amrDomains[a_level]);

            HeadInterp.interpToFine(*newHeadPtr, *m_head[a_level - 1]);
            GapHeightInterp.interpToFine(*newGapHeightPtr, *m_gapheight[a_level - 1]);

            // can only copy from existing level if we're not on the
            // newly created level
            // if (a_level != new_finest_level)
            if (m_head[a_level]->isDefined())
            {
                m_head[a_level]->copyTo(*newHeadPtr);
                m_gapheight[a_level]->copyTo(*newGapHeightPtr);
            }

            // finally, do an exchange (this may wind up being unnecessary)
            newHeadPtr->exchange();
            newGapHeightPtr->exchange();

            delete (m_head[a_level]);
            delete (m_gapheight[a_level]);
            m_head[a_level] = newHeadPtr;
            m_gapheight[a_level] = newGapHeightPtr;
        }
    } // end interpolate/copy new velocity

    // probably eventually want to do this differently
    RealVect dx = m_amrDx[a_level] * RealVect::Unit;
}

void
AmrHydro::initData(Vector<LevelData<FArrayBox>*>& a_head)

{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::initData" << endl;
    }

    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        LevelData<FArrayBox>& levelHead      = *m_head[lev];
        LevelData<FArrayBox>& levelGapHeight = *m_gapheight[lev];
        LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
        LevelData<FArrayBox>& levelqw        = *m_qw[lev];
        LevelData<FArrayBox>& levelRe        = *m_Re[lev];
        LevelData<FArrayBox>& levelmR        = *m_meltRate[lev];
        LevelData<FArrayBox>& levelIceHeight = *m_iceheight[lev];
        LevelData<FArrayBox>& levelzBed      = *m_bedelevation[lev];
        LevelData<FArrayBox>& levelPi        = *m_overburdenpress[lev];
        LevelData<FArrayBox>& levelgradzBed  = *m_gradZb[lev];
        if (lev > 0)
        {
            // AF: !!! deal with multiple levels later !!!
            // fill the ghost cells of a_vectCoordSys[lev]->getH();
            // m_head
            LevelData<FArrayBox>& coarseHead = *m_head[lev - 1];
            int nGhost = levelHead.ghostVect()[0];
            PiecewiseLinearFillPatch headFiller(m_amrGrids[lev],
                                               m_amrGrids[lev - 1],
                                               levelHead.nComp(),
                                               m_amrDomains[lev - 1],
                                               m_refinement_ratios[lev - 1],
                                               nGhost);
            headFiller.fillInterp(levelHead, coarseHead, coarseHead, 0.0, 0, 0, 1);

            // m_gapheight
            LevelData<FArrayBox>& coarseGapHeight = *m_gapheight[lev - 1];
            //int nGhost = levelGapHeight.ghostVect()[0]; Same number of ghost cells ?
            PiecewiseLinearFillPatch gapHeightFiller(m_amrGrids[lev],
                                               m_amrGrids[lev - 1],
                                               levelGapHeight.nComp(),
                                               m_amrDomains[lev - 1],
                                               m_refinement_ratios[lev - 1],
                                               nGhost);
            gapHeightFiller.fillInterp(levelGapHeight, coarseGapHeight, coarseGapHeight, 0.0, 0, 0, 1);
        }

        RealVect levelDx = m_amrDx[lev] * RealVect::Unit;
        m_IBCPtr->define(m_amrDomains[lev], levelDx[0]);
        // int refRatio = (lev > 0)?m_refinement_ratios[lev-1]:0;

        m_IBCPtr->initializeData(levelDx, 
                                 *m_suhmoParm,       
                                 levelHead, levelGapHeight, 
                                 levelPw, levelqw,
                                 levelRe, levelmR,
                                 levelzBed, levelPi,
                                 levelIceHeight, levelgradzBed);

        // initialize oldPhi to be the current value
        levelHead.copyTo(*m_old_head[lev]);
        levelGapHeight.copyTo(*m_old_gapheight[lev]);
    }

    // may be necessary to average down here
    for (int lev = m_finest_level; lev > 0; lev--)
    {
        CoarseAverage avgDown(m_amrGrids[lev], m_head[lev]->nComp(), m_refinement_ratios[lev - 1]);
        avgDown.averageToCoarse(*m_head[lev - 1], *m_head[lev]);
    }

//#define writePlotsImmediately
#ifdef writePlotsImmediately
    if (m_plot_interval >= 0)
    {
#ifdef CH_USE_HDF5
        writePlotFile();
#endif
    }
#endif
}

void
AmrHydro::defineSolver()
{
    // just a stub for when we might actually need it...
}

// compute timestep
Real
AmrHydro::computeDt()
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::computeDt" << endl;
    }

    if (m_fixed_dt > TINY_NORM) return m_fixed_dt;

    Real dt = 1.0e50;
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        Real dtLev = dt;
        //// pretend phi is velocity here --> totally wrong for problem at hand ...
        //const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        //const LevelData<FArrayBox>& levelVel = *m_head[lev];
        //DataIterator levelDit = levelVel.dataIterator();
        //for (levelDit.reset(); levelDit.ok(); ++levelDit)
        //{
        //    int p = 0;
        //    const Box& gridBox = levelGrids[levelDit];
        //    Real maxVel = 1.0 + levelVel[levelDit].norm(gridBox, p, 0, 1);
        //    maxVel = max(maxVel,1.0);    
        //    Real localDt = m_amrDx[lev][0] / maxVel;
        //    dtLev = min(dtLev, localDt);
        //}

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
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::computeInitialDt" << endl;
    }

    // for now, just call computeDt;
    Real dt = computeDt();
    return dt;
}

#ifdef CH_USE_HDF5
/// write hdf5 plotfile to the standard location
void
AmrHydro::writePltCustom(int numPlotComps, 
                         Vector<std::string>& vectName,
                         Vector<Vector<LevelData<FArrayBox>*>> stuffToPlot,
                         string namePlot)
{
    if (m_verbosity > 3)
    {
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
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::writePlotFile" << endl;
    }

    // plot comps: head + gapHeight + bedelevation + overburdenPress
    int numPlotComps = 9;
    //int numPlotComps = 3;

    // add in grad(head) if desired
    //if (m_write_gradPhi)
    //{
    //    numPlotComps += SpaceDim;
    //}

    // generate data names

    string headName("head");
    string gapHeightName("gapHeight");
    string zbedName("bedelevation");
    string piName("overburdenPress");
    string pwName("Pw"); 
    string qwName("Qw_x"); 
    string ReName("Re"); 
    string meltRateName("meltRate"); 
    string xGradName("GradHead_x");

    Vector<string> vectName(numPlotComps);
    // int dThicknessComp;

    vectName[0] = headName;
    vectName[1] = gapHeightName;
    vectName[2] = zbedName;
    vectName[3] = piName;
    vectName[4] = pwName;
    vectName[5] = qwName;
    vectName[6] = ReName;
    vectName[7] = meltRateName;
    vectName[8] = xGradName;
    //if (m_write_gradPhi)
    //{
    //    vectName[numPlotComps] = xGradName;
    //    if (SpaceDim > 1) vectName[numPlotComps + 1] = yGradName;
    //    if (SpaceDim > 2) vectName[numPlotComps + 2] = zGradName;
    //}

    Box domain = m_amrDomains[0].domainBox();
    Real dt = 1.;
    int numLevels = m_finest_level + 1;
    // compute plot data
    Vector<LevelData<FArrayBox>*> plotData(m_head.size(), NULL);

    // ghost vect makes things simpler
    IntVect ghostVect(IntVect::Unit);

    for (int lev = 0; lev < numLevels; lev++)
    {
        // first allocate storage
        plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev], numPlotComps, ghostVect);
    }

    for (int lev = 0; lev < numLevels; lev++)
    {
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
        for (dit.begin(); dit.ok(); ++dit)
        {
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
            thisPlotData.copy(thisRe, 0, comp, 1);
            comp++;
            thisPlotData.copy(thismR, 0, comp, 1);
            comp++;
            thisPlotData.copy(thisGradHead, 0, comp, 1);
            comp++;
            // now copy for grad(head)
            //if (m_write_gradPhi)
            //{
            //    const FArrayBox& thisGradPhi = levelGradPhi[dit];
            //    thisPlotData.copy(thisGradPhi, 0, comp, SpaceDim);
            //    comp += SpaceDim;

            //} // end if we are writing grad(head)

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

    if (SpaceDim == 1)
    {
        filename.append("1d.hdf5");
    }
    else if (SpaceDim == 2)
    {
        filename.append("2d.hdf5");
    }
    else if (SpaceDim == 3)
    {
        filename.append("3d.hdf5");
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

/// write checkpoint file out for later restarting
void
AmrHydro::writeCheckpointFile() const
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::writeCheckpointfile" << endl;
    }

    CH_TIME("AmrHydro::writeCheckpointFile");

    string headName("head");
    Vector<string> vectName(1);
    for (int comp = 0; comp < 1; comp++)
    {
        char idx[4];
        sprintf(idx, "%d", comp);
        vectName[comp] = headName + string(idx);
    }
    Box domain = m_amrDomains[0].domainBox();
    // int numLevels = m_finest_level +1;

    // generate checkpointfile name
    char(iter_str[100]);

    if (m_check_overwrite)
    {
        // overwrite the same checkpoint file, rather than re-writing them
        sprintf(iter_str, "%s.%dd.hdf5", m_check_prefix.c_str(), SpaceDim);
    }
    else
    {
        // or hang on to them, if you are a bit sentimental. It's better than keeping
        // every core dump you generate.
        sprintf(iter_str, "%s%06d.%dd.hdf5", m_check_prefix.c_str(), m_cur_step, SpaceDim);
    }

    if (m_verbosity > 3)
    {
        pout() << "checkpoint file name = " << iter_str << endl;
    }

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
    header.m_int["num_comps"] = 1;
    // at the moment, save cfl, but it can be changed by the inputs
    // file if desired.
    header.m_real["cfl"] = m_cfl;

    // periodicity info
    D_TERM(if (m_amrDomains[0].isPeriodic(0)) header.m_int["is_periodic_0"] = 1; else header.m_int["is_periodic_0"] = 0;
           ,

           if (m_amrDomains[0].isPeriodic(1)) header.m_int["is_periodic_1"] = 1;
           else header.m_int["is_periodic_1"] = 0;
           ,

           if (m_amrDomains[0].isPeriodic(2)) header.m_int["is_periodic_2"] = 1;
           else header.m_int["is_periodic_2"] = 0;);

    // set up component names
    char compStr[30];
    // string thicknessName("thickness");
    string compName;
    int nComp = 0;
    for (int comp = 0; comp < 1; comp++)
    {
        // first generate component name
        char idx[5];
        sprintf(idx, "%04d", comp);
        compName = headName + string(idx);
        sprintf(compStr, "component_%04d", comp);
        header.m_string[compStr] = compName;
    }
    nComp++;

    header.writeToFile(handle);

    // now loop over levels and write out each level's data
    // note that we loop over all allowed levels, even if they
    // are not defined at the moment.
    for (int lev = 0; lev <= m_max_level; lev++)
    {
        // set up the level string
        char levelStr[20];
        sprintf(levelStr, "%d", lev);
        const std::string label = std::string("level_") + levelStr;

        handle.setGroup(label);

        // set up the header info
        HDF5HeaderData levelHeader;
        if (lev < m_max_level)
        {
            levelHeader.m_int["ref_ratio"] = m_refinement_ratios[lev];
        }
        levelHeader.m_real["dx"] = m_amrDx[lev][0];
        levelHeader.m_box["prob_domain"] = m_amrDomains[lev].domainBox();

        levelHeader.writeToFile(handle);

        // now write the data for this level
        // only try to write data if level is defined.
        if (lev <= m_finest_level)
        {
            write(handle, m_amrGrids[lev]);

            const IntVect ghost = IntVect::Unit * 2;

            write(handle, *m_head[lev], "headData", m_head[lev]->ghostVect());

        } // end loop over levels
    }

    handle.close();
}

/// read checkpoint file for restart
void
AmrHydro::readCheckpointFile(HDF5Handle& a_handle)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::readCheckpointFile" << endl;
    }

    CH_TIME("AmrHydro::readCheckpointFile");

    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if (m_verbosity >= 3)
    {
        pout() << "hdf5 header data: " << endl;
        pout() << header << endl;
    }

    // read max level
    if (header.m_int.find("max_level") == header.m_int.end())
    {
        MayDay::Error("checkpoint file does not contain max_level");
    }
    // we can change max level upon restart
    int max_level_check = header.m_int["max_level"];
    if (max_level_check != m_max_level)
    {
        if (m_verbosity > 0)
        {
            pout() << "Restart file has a different max level than inputs file" << endl;
            pout() << "     max level from inputs file = " << m_max_level << endl;
            pout() << "     max level in checkpoint file = " << max_level_check << endl;
            pout() << "Using max level from inputs file" << endl;
        }
    }
    // read finest level
    if (header.m_int.find("finest_level") == header.m_int.end())
    {
        MayDay::Error("checkpoint file does not contain finest_level");
    }

    m_finest_level = header.m_int["finest_level"];
    if (m_finest_level > m_max_level)
    {
        MayDay::Error("finest level in restart file > max allowable level!");
    }

    // read current step
    if (header.m_int.find("current_step") == header.m_int.end())
    {
        MayDay::Error("checkpoint file does not contain current_step");
    }

    m_cur_step = header.m_int["current_step"];
    m_restart_step = m_cur_step;

    // read time
    if (header.m_real.find("time") == header.m_real.end())
    {
        MayDay::Error("checkpoint file does not contain time");
    }

    m_time = header.m_real["time"];

    // read timestep
    if (header.m_real.find("dt") == header.m_real.end())
    {
        MayDay::Error("checkpoint file does not contain dt");
    }

    m_dt = header.m_real["dt"];

    // read num comps
    if (header.m_int.find("num_comps") == header.m_int.end())
    {
        MayDay::Error("checkpoint file does not contain num_comps");
    }

    // read cfl
    if (header.m_real.find("cfl") == header.m_real.end())
    {
        MayDay::Error("checkpoint file does not contain cfl");
    }

    Real check_cfl = header.m_real["cfl"];
    ParmParse ppCheck("AMRDriver");

    if (ppCheck.contains("cfl"))
    {
        // check for consistency and warn if different
        if (check_cfl != m_cfl)
        {
            if (m_verbosity > 0)
            {
                pout() << "CFL in checkpoint file different from inputs file" << endl;
                pout() << "     cfl in inputs file = " << m_cfl << endl;
                pout() << "     cfl in checkpoint file = " << check_cfl << endl;
                pout() << "Using cfl from inputs file" << endl;
            }
        } // end if cfl numbers differ
    }     // end if cfl present in inputs file
    else
    {
        m_cfl = check_cfl;
    }

    // read periodicity info
    // Get the periodicity info -- this is more complicated than it really
    // needs to be in order to preserve backward compatibility
    bool isPeriodic[SpaceDim];
    D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end())) isPeriodic[0] =
               (header.m_int["is_periodic_0"] == 1);
           else isPeriodic[0] = false;
           ,

           if (!(header.m_int.find("is_periodic_1") == header.m_int.end())) isPeriodic[1] =
               (header.m_int["is_periodic_1"] == 1);
           else isPeriodic[1] = false;
           ,

           if (!(header.m_int.find("is_periodic_2") == header.m_int.end())) isPeriodic[2] =
               (header.m_int["is_periodic_2"] == 1);
           else isPeriodic[2] = false;);

    // now resize stuff
    m_amrDomains.resize(m_max_level + 1);
    m_amrGrids.resize(m_max_level + 1);
    m_amrDx.resize(m_max_level + 1);
    m_old_head.resize(m_max_level + 1, NULL);
    m_head.resize(m_max_level + 1, NULL);

    // now read in level-by-level data
    for (int lev = 0; lev <= max_level_check; lev++)
    {
        // set up the level string
        char levelStr[20];
        sprintf(levelStr, "%d", lev);
        const std::string label = std::string("level_") + levelStr;

        a_handle.setGroup(label);

        // read header info
        HDF5HeaderData header;
        header.readFromFile(a_handle);

        if (m_verbosity >= 3)
        {
            pout() << "level " << lev << " header data" << endl;
            pout() << header << endl;
        }

        // Get the refinement ratio
        if (lev < max_level_check)
        {
            int checkRefRatio;
            if (header.m_int.find("ref_ratio") == header.m_int.end())
            {
                MayDay::Error("checkpoint file does not contain ref_ratio");
            }
            checkRefRatio = header.m_int["ref_ratio"];

            // check for consistency
            if (checkRefRatio != m_refinement_ratios[lev])
            {
                MayDay::Error("inputs file and checkpoint file ref ratios inconsistent");
            }
        }

        // read dx
        if (header.m_real.find("dx") == header.m_real.end())
        {
            MayDay::Error("checkpoint file does not contain dx");
        }

        m_amrDx[lev] = RealVect::Unit * (header.m_real["dx"]);

        // read problem domain box
        if (header.m_box.find("prob_domain") == header.m_box.end())
        {
            MayDay::Error("checkpoint file does not contain prob_domain");
        }
        Box domainBox = header.m_box["prob_domain"];

        m_amrDomains[lev] = ProblemDomain(domainBox, isPeriodic);

        // the rest is only applicable if this level is defined
        if (lev <= m_finest_level)
        {
            // read grids
            Vector<Box> grids;
            const int grid_status = read(a_handle, grids);
            if (grid_status != 0)
            {
                MayDay::Error("checkpoint file does not contain a Vector<Box>");
            }
            // do load balancing
            int numGrids = grids.size();
            Vector<int> procIDs(numGrids);
            LoadBalance(procIDs, grids);
            DisjointBoxLayout levelDBL(grids, procIDs, m_amrDomains[lev]);
            m_amrGrids[lev] = levelDBL;

            // allocate this level's storage
            IntVect nGhost = m_num_head_ghost * IntVect::Unit;
            m_old_head[lev] = new LevelData<FArrayBox>(levelDBL, 1, nGhost);
            m_head[lev]     = new LevelData<FArrayBox>(levelDBL, SpaceDim, nGhost);

            // read this level's data

            LevelData<FArrayBox>& old_head = *m_old_head[lev];
            LevelData<FArrayBox> tmpHead;
            tmpHead.define(old_head);

            int dataStatus = read<FArrayBox>(a_handle, tmpHead, "headData", levelDBL);
            for (DataIterator dit(levelDBL); dit.ok(); ++dit)
            {
                old_head[dit].copy(tmpHead[dit]);
            }

            if (dataStatus != 0)
            {
                MayDay::Error("checkpoint file does not contain head data");
            }

        } // end if this level is defined
    }     // end loop over levels

    // do we need to close the handle?

    // defineSolver();
}

/// set up for restart
void
AmrHydro::restart(string& a_restart_file)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrHydro::restart" << endl;
    }

    HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
    // first read in data from checkpoint file
    readCheckpointFile(handle);
    handle.close();
    // don't think I need to do anything else, do I?
}

#endif

#include "NamespaceFooter.H"
