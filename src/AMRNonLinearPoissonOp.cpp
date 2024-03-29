#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FORT_PROTO.H"
#include "BoxIterator.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "FineInterp.H"
#include "CoarseAverage.H"
#include "CH_OpenMP.H"
#include "AMRMultiGrid.H"
#include "Misc.H"
#include "computeSum.H"

#include "AMRNonLinearPoissonOp.H"
#include "AMRNonLinearPoissonOpF_F.H"
#include "CCProjectorF_F.H"
#include "MACProjectorF_F.H"

#include "NamespaceHeader.H"

int AMRNonLinearPoissonOp::s_exchangeMode = 1; // 1: no overlap (default); 0: ...
int AMRNonLinearPoissonOp::s_relaxMode = 1; // 1: GSRB (not implemented); 4: Jacobi(not implemented), 6: GS (not red black)
int AMRNonLinearPoissonOp::s_maxCoarse = 2;

// ---------------------------------------------------------
static void
amrpgetMultiColors(Vector<IntVect>& a_colors)
{

#if CH_SPACEDIM==2
  a_colors.resize(4);
  a_colors[0] = IntVect::Zero;             // (0,0)
  a_colors[1] = IntVect::Unit;             // (1,1)
  a_colors[2] = IntVect::Zero + BASISV(1); // (0,1)
  a_colors[3] = IntVect::Zero + BASISV(0); // (1,0)
#elif CH_SPACEDIM==3
  a_colors.resize(8);
  a_colors[0] = IntVect::Zero;                         // (0,0,0)
  a_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1); // (1,1,0)
  a_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2); // (0,1,1)
  a_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2); // (1,0,1)
  a_colors[4] = IntVect::Zero + BASISV(1);             // (0,1,0)
  a_colors[5] = IntVect::Zero + BASISV(0);             // (1,0,0)
  a_colors[6] = IntVect::Zero + BASISV(2);             // (0,0,1)
  a_colors[7] = IntVect::Unit;                         // (1,1,1)
#endif
}

// ---------------------------------------------------------

AMRNonLinearPoissonOp::AMRNonLinearPoissonOp()
    :m_verbosity(3),
     m_print(false),
     m_use_FAS(true)
{
}

AMRNonLinearPoissonOp::~AMRNonLinearPoissonOp()
{
}

/** full define function for AMRLevelOp with both coarser and finer levels */
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                                   const DisjointBoxLayout& a_gridsFiner,
                                   const DisjointBoxLayout& a_gridsCoarser,
                                   const RealVect&          a_dxLevel,
                                   int                      a_refRatio,
                                   int                      a_refRatioFiner,
                                   const ProblemDomain&     a_domain,
                                   BCHolder                 a_bc,
                                   const Copier&            a_exchange,
                                   const CFRegion&          a_cfregion,
                                   const int                a_nComp)

{
  CH_TIME("AMRNonLinearPoissonOp::define1");

  amrpgetMultiColors(m_colors);

  this->define(a_grids, a_gridsCoarser, a_dxLevel, a_refRatio, a_domain, a_bc,
               a_exchange, a_cfregion, a_nComp);
  m_refToFiner = a_refRatioFiner;

  if (a_gridsFiner.isClosed())
    {
      ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
      m_levfluxreg.define(a_gridsFiner,
                          a_grids,
                          fineDomain,
                          m_refToFiner,
                          a_nComp /* ncomp*/);
    }
}

// ---------------------------------------------------------
/** full define function for AMRLevelOp<LevelData<FArrayBox> > with finer levels, but no coarser */
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                                   const DisjointBoxLayout& a_gridsFiner,
                                   const RealVect&          a_dxLevel,
                                   int                      a_refRatio, // dummy arg
                                   int                      a_refRatioFiner,
                                   const ProblemDomain&     a_domain,
                                   BCHolder                 a_bc,
                                   const Copier&            a_exchange,
                                   const CFRegion&          a_cfregion,
                                   const int                a_nComp)
{
  CH_TIME("AMRNonLinearPoissonOp::define2");

  CH_assert(a_refRatio == 1);

  amrpgetMultiColors(m_colors);

  // calls the MG version of define
  this->define(a_grids, a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion);
  m_refToFiner = a_refRatioFiner;

  ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
  m_levfluxreg.define(a_gridsFiner,
                      a_grids,
                      fineDomain,
                      m_refToFiner,
                      a_nComp /* ncomp*/);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                                   const DisjointBoxLayout& a_coarse,
                                   const RealVect&          a_dxLevel,
                                   int                      a_refRatio,
                                   const ProblemDomain&     a_domain,
                                   BCHolder                 a_bc,
                                   const Copier&            a_exchange,
                                   const CFRegion&          a_cfregion,
                                   int a_numComp)
{
  CH_TIME("AMRNonLinearPoissonOp::define3");

  amrpgetMultiColors(m_colors);

  this->define(a_grids, a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion);
  m_refToCoarser = a_refRatio;

  m_dxCrse      = a_refRatio*a_dxLevel[0];
  m_dxCrse_vect = a_refRatio*a_dxLevel;
  m_refToFiner = 1;

  m_interpWithCoarser.define(a_grids, &a_coarse, a_dxLevel,
                             m_refToCoarser, a_numComp, m_domain);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                                   const RealVect&          a_dx,
                                   const ProblemDomain&     a_domain,
                                   BCHolder                 a_bc,
                                   const Copier&            a_exchange,
                                   const CFRegion&          a_cfregion)
{
  CH_TIME("AMRNonLinearPoissonOp::define4");

  amrpgetMultiColors(m_colors);

  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx           = a_dx[0];
  m_dx_vect      = a_dx;
  m_dxCrse       = 2*a_dx[0];
  m_dxCrse_vect  = 2*a_dx;

  // redefined in AMRLevelOp<LevelData<FArrayBox> >::define virtual function.
  m_refToCoarser = 2;
  m_refToFiner   = 2;

  // these get set again after define is called
  m_alpha = 0.0;
  m_beta  = 1.0;

  m_use_FAS = true;
  m_exchangeCopier = a_exchange;
  // m_exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  // m_exchangeCopier.trimEdges(a_grids, IntVect::Unit);

  m_cfregion = a_cfregion;

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                                   const RealVect&          a_dx,
                                   const ProblemDomain&     a_domain,
                                   BCHolder                 a_bc)
{
  CH_TIME("AMRNonLinearPoissonOp::define5");

  amrpgetMultiColors(m_colors);

  // 
  Copier copier;
  copier.define(a_grids, a_grids, IntVect::Unit, true);

  CFRegion cfregion(a_grids, a_domain);

  this->define(a_grids, a_dx, a_domain, a_bc, copier, cfregion);
}

void AMRNonLinearPoissonOp::define(const DisjointBoxLayout& a_grids,
                                   const DisjointBoxLayout* a_baseBAPtr,
                                   RealVect                 a_dx,
                                   int                      a_refRatio,
                                   const ProblemDomain&     a_domain,
                                   BCHolder                 a_bc)
{
  Copier exchangeCopier;
  exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
  exchangeCopier.trimEdges(a_grids, IntVect::Unit);

  CFRegion cfregion;
  cfregion.define(a_grids, a_domain);

  if (a_baseBAPtr != NULL)
    {
      define(a_grids, *a_baseBAPtr, a_dx, a_refRatio, a_domain, a_bc, exchangeCopier, cfregion);
    }
  else
    {
      define(a_grids, a_dx, a_domain, a_bc, exchangeCopier, cfregion);
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::residual(LevelData<FArrayBox>&       a_lhs,
                                     const LevelData<FArrayBox>& a_phi,
                                     const LevelData<FArrayBox>& a_rhs,
                                     bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::residual");

  if (a_homogeneous && (!m_use_FAS)) {
      homogeneousCFInterp((LevelData<FArrayBox>&)a_phi);
      residualI(a_lhs,a_phi,a_rhs,a_homogeneous);
  } else {
      residualI(a_lhs,a_phi,a_rhs,false);
  }
}

void AMRNonLinearPoissonOp::residualNF(LevelData<FArrayBox>&                a_lhs,
                                       LevelData<FArrayBox>&                a_phi,
                                       const LevelData<FArrayBox>*          a_phiCoarse,
                                       const LevelData<FArrayBox>&          a_rhs,
                                       bool                                 a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::residualNF");
  
  if (a_homogeneous && (!m_use_FAS)) {
      homogeneousCFInterp((LevelData<FArrayBox>&)a_phi);
  } else {
      if (a_phiCoarse != NULL) {
          m_interpWithCoarser.coarseFineInterp(a_phi, *a_phiCoarse);
      }
 }

  residualI(a_lhs,a_phi,a_rhs,a_homogeneous);
}


// ---------------------------------------------------------
void AMRNonLinearPoissonOp::UpdateOperator(const LevelData<FArrayBox>&   a_phi, 
                                          const LevelData<FArrayBox>*   a_phicoarsePtr, 
                                          int                           a_depth,
                                          int                           a_AMRFASMGiter,
                                          bool                          a_homogeneous)
{
    pout() << "AMRNonLinearPoissonOp::UpdateOperator does nothing \n";
} 

void AMRNonLinearPoissonOp::AverageOperator(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                                            int a_depth) 
{
    pout() <<"No AverageOperator in AMRNonLinearPoissonOp \n";
}

void AMRNonLinearPoissonOp::residualI(LevelData<FArrayBox>&       a_lhs,
                                      const LevelData<FArrayBox>& a_phi,
                                      const LevelData<FArrayBox>& a_rhs,
                                      bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::residualI");
  if (m_verbosity > 3) {
      pout() << "AMRNonLinearPoissonOp::residualI\n";
  }

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  if (s_exchangeMode == 0)
    phi.exchange(phi.interval(), m_exchangeCopier);
  else if (s_exchangeMode == 1)
    phi.exchangeNoOverlap(m_exchangeCopier);
  else
    MayDay::Abort("exchangeMode");

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                        *m_B, *m_iceMask, *m_Pi, *m_zb);

  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("AMRNonLinearPoissonOp::BCs");

    for (dit.begin(); dit.ok(); ++dit)
      {
        if (!m_use_FAS) {
            m_bc(phi[dit], dbl[dit], m_domain, m_dx_vect, a_homogeneous);
        } else {
            m_bc(phi[dit], dbl[dit], m_domain, m_dx_vect, false);
        }
      }
  }

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit];
      FORT_OPERATORLAPRESNL(CHF_FRA(a_lhs[dit]),
                          CHF_CONST_FRA(phi[dit]),
                          CHF_CONST_FRA(a_rhs[dit]),
                          CHF_BOX(region),
                          CHF_CONST_REALVECT(m_dx_vect),
                          CHF_CONST_REAL(m_alpha),
                          CHF_CONST_REAL(m_beta),
                          CHF_CONST_FRA(a_nlfunc[dit]));
    }
}

// ---------------------------------------------------------
/**************************/
// Preconditionior is not used for FAS solve, but retained so
// regular multigrid can still be used (which should be identical to FAS
// multigrid for linear problems)
void AMRNonLinearPoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                                    const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::preCond");

  // diagonal term of this operator is (alpha - 4 * beta/h/h) in 2D,
  // (alpha - 6 * beta/h/h) in 3D,
  // so inverse of this is our initial multiplier

  CH_assert(a_phi.nComp() == a_rhs.nComp());

  Real sum_b = ( D_TERM(  2.0 * m_beta / (m_dx_vect[0]*m_dx_vect[0]),
                        + 2.0 * m_beta / (m_dx_vect[1]*m_dx_vect[1]),
                        + 2.0 * m_beta / (m_dx_vect[2]*m_dx_vect[2]) )
               );
  Real mult = 1.0 / (m_alpha - sum_b);

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  int nbox = dit.size();

#pragma omp parallel for
    for(int ibox=0; ibox<nbox; ibox++) {
        a_phi[dit[ibox]].copy(a_rhs[dit[ibox]]);
        a_phi[dit[ibox]] *= mult;
    }
 //end pragma
  int dummyDepth = 0;
  m_print = false;
  relax(a_phi, a_rhs, 2, dummyDepth, dummyDepth);
}

void AMRNonLinearPoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                                    const LevelData<FArrayBox>& a_res,
                                    const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::preCond");

  // diagonal term of this operator is (alpha - 4 * beta/h/h) in 2D,
  // (alpha - 6 * beta/h/h) in 3D,
  // so inverse of this is our initial multiplier

  CH_assert(a_phi.nComp() == a_rhs.nComp());

  Real sum_b = ( D_TERM(  2.0 * m_beta / (m_dx_vect[0]*m_dx_vect[0]),
                        + 2.0 * m_beta / (m_dx_vect[1]*m_dx_vect[1]),
                        + 2.0 * m_beta / (m_dx_vect[2]*m_dx_vect[2]) )
               );
  Real mult = 1.0 / (m_alpha - sum_b);

  incr(a_phi, a_res, mult);

  int dummyDepth = 0;
  m_print = false;
  relax(a_phi, a_rhs, 2, dummyDepth, dummyDepth);
}

void AMRNonLinearPoissonOp::applyOpMg(LevelData<FArrayBox>& a_lhs,
                                      LevelData<FArrayBox>& a_phi,
                                      LevelData<FArrayBox>* a_phiCoarse,
                                      bool a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOpMg");

  // Do CF stuff if we have a coarser level that's not just a single grid cell
   if (a_phiCoarse != NULL)
   {
     const ProblemDomain& probDomain = a_phiCoarse->disjointBoxLayout().physDomain();
     const Box& domBox = probDomain.domainBox();
     //    IntVect hi = domBox.b
     if (domBox.bigEnd() != domBox.smallEnd())
     {
       m_interpWithCoarser.coarseFineInterp(a_phi, *a_phiCoarse);
     }
   }

   applyOpI(a_lhs, a_phi, a_homogeneous);

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::applyOp(LevelData<FArrayBox>&       a_lhs,
                                    const LevelData<FArrayBox>& a_phi,
                                    bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOp");

  if (a_homogeneous && (!m_use_FAS)) {
      homogeneousCFInterp((LevelData<FArrayBox>&)a_phi);
      applyOpI(a_lhs,a_phi,a_homogeneous);
  } else {
      applyOpI(a_lhs,a_phi,false);
  }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::applyOpI(LevelData<FArrayBox>&       a_lhs,
                                     const LevelData<FArrayBox>& a_phi,
                                     bool                        a_homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOpI");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  RealVect dx = m_dx_vect;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  for (dit.begin(); dit.ok(); ++dit) {
      if (!m_use_FAS) {
          m_bc(phi[dit], dbl[dit()], m_domain, dx, a_homogeneous);
      } else {
          m_bc(phi[dit], dbl[dit()], m_domain, dx, false);
     }
  }

  phi.exchange(phi.interval(), m_exchangeCopier);

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                        *m_B, *m_iceMask, *m_Pi, *m_zb);

  for (dit.begin(); dit.ok(); ++dit) {
        const Box& region = dbl[dit];
        FORT_OPERATORLAPNL(CHF_FRA(a_lhs[dit]),
                           CHF_CONST_FRA(phi[dit]),
                           CHF_BOX(region),
                           CHF_CONST_REALVECT(m_dx_vect),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_REAL(m_beta),
                           CHF_CONST_FRA(a_nlfunc[dit]));
  }
}

void AMRNonLinearPoissonOp::applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                              const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("AMRNonLinearPoissonOp::applyOpNoBoundary");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                        *m_B, *m_iceMask, *m_Pi, *m_zb);

  DataIterator dit = phi.dataIterator();
  int nbox=dit.size();
  phi.exchange(phi.interval(), m_exchangeCopier);

#pragma omp parallel 
  {
#pragma omp for 
    for (int ibox = 0; ibox < nbox; ibox++)
    {
      const Box& region = dbl[dit[ibox]];
      FORT_OPERATORLAPNL(CHF_FRA(a_lhs[dit[ibox]]),
                         CHF_CONST_FRA(phi[dit[ibox]]),
                         CHF_BOX(region),
                         CHF_CONST_REALVECT(m_dx_vect),
                         CHF_CONST_REAL(m_alpha),
                         CHF_CONST_REAL(m_beta),
                         CHF_CONST_FRA(a_nlfunc[dit[ibox]]));
    }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::create(LevelData<FArrayBox>&       a_lhs,
                                   const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::create");

  m_levelOps.create(a_lhs, a_rhs);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                            const LevelData<FArrayBox>& a_rhs,
                                            const int &                 a_refRat)
{
  CH_TIME("AMRNonLinearPoissonOp::createCoarsened");

  int ncomp = a_rhs.nComp();
  IntVect ghostVect = a_rhs.ghostVect();

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();
  CH_assert(dbl.coarsenable(a_refRat));

  DisjointBoxLayout dblCoarsenedFine;

  if (a_refRat == 2)
    {
      if (m_coarsenedMGrids.size() == 0)
        coarsen(m_coarsenedMGrids, dbl, 2);
      dblCoarsenedFine = m_coarsenedMGrids;
    }
  else
    {
      coarsen(dblCoarsenedFine, dbl, a_refRat);
    }

  a_lhs.define(dblCoarsenedFine, ncomp, ghostVect );
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::assign(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::assign");

  m_levelOps.assign(a_lhs, a_rhs);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::assignLocal(LevelData<FArrayBox>&        a_lhs,
                               const  LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::assignLocal");

  for (DataIterator dit= a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      a_lhs[dit].copy(a_rhs[dit]);
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::buildCopier(Copier&                      a_copier,
                               const  LevelData<FArrayBox>& a_lhs,
                               const  LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::buildCopier");

  const DisjointBoxLayout& dbl=a_lhs.disjointBoxLayout();

  a_copier.define(a_rhs.disjointBoxLayout(), dbl, IntVect::Zero);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::assignCopier( LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_rhs,
                                 const Copier&               a_copier)
{
  CH_TIME("AMRNonLinearPoissonOp::assignCopier");
  a_rhs.copyTo(a_rhs.interval(), a_lhs, a_lhs.interval(), a_copier);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::zeroCovered( LevelData<FArrayBox>& a_lhs,
                                LevelData<FArrayBox>& a_rhs,
                                const Copier&         a_copier)
{
  CH_TIME("AMRNonLinearPoissonOp::zeroCovered");

  m_levelOps.copyToZero(a_lhs, a_copier);
}

// ---------------------------------------------------------
Real AMRNonLinearPoissonOp::dotProduct(const LevelData<FArrayBox>& a_1,
                                       const LevelData<FArrayBox>& a_2)
{
  CH_TIME("AMRNonLinearPoissonOp::dotProduct");

  return m_levelOps.dotProduct(a_1, a_2);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::mDotProduct(const LevelData<FArrayBox>& a_1,
                                        const int a_sz,
                                        const LevelData<FArrayBox> a_2[],
                                        Real a_mdots[])
{
  CH_TIME("AMRNonLinearPoissonOp::mDotProduct");

  m_levelOps.mDotProduct(a_1, a_sz, a_2, a_mdots);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::incr( LevelData<FArrayBox>&       a_lhs,
                                  const LevelData<FArrayBox>& a_x,
                                  Real                        a_scale)
{
  CH_TIME("AMRNonLinearPoissonOp::incr");

  m_levelOps.incr(a_lhs, a_x, a_scale);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::axby( LevelData<FArrayBox>&       a_lhs,
                                  const LevelData<FArrayBox>& a_x,
                                  const LevelData<FArrayBox>& a_y,
                                  Real                        a_a,
                                  Real                        a_b)
{
  CH_TIME("AMRNonLinearPoissonOp::axby");

  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::scale(LevelData<FArrayBox>& a_lhs,
                                  const Real&           a_scale)
{
  CH_TIME("AMRNonLinearPoissonOp::scale");

  m_levelOps.scale(a_lhs, a_scale);
}

// ---------------------------------------------------------
Real AMRNonLinearPoissonOp::norm(const LevelData<FArrayBox>& a_x,
                                 int                         a_ord)
{
  CH_TIME("AMRNonLinearPoissonOp::norm");

  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}

// ---------------------------------------------------------
Real AMRNonLinearPoissonOp::localMaxNorm(const LevelData<FArrayBox>& a_x)
{
  CH_TIME("AMRNonLinearPoissonOp::localMaxNorm");

  Real localMax = 0;
  int nComp=a_x.nComp();
  for (DataIterator dit=a_x.dataIterator(); dit.ok(); ++dit)
    {
      localMax = Max(localMax, a_x[dit].norm(a_x.box(dit()), 0, 0, nComp));
    }
  return localMax;
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::setToZero(LevelData<FArrayBox>& a_lhs)
{
  CH_TIME("AMRNonLinearPoissonOp::setToZero");

  m_levelOps.setToZero(a_lhs);
}

void AMRNonLinearPoissonOp::relaxNF(LevelData<FArrayBox>&       a_e,
                                    const LevelData<FArrayBox>* a_eCoarse,
                                    const LevelData<FArrayBox>& a_residual,
                                    int                         a_iterations,
                                    int                         a_AMRFASMGiter,
                                    int                         a_depth,
                                    bool                        a_print)
{
  if (a_eCoarse != NULL) {
    m_interpWithCoarser.coarseFineInterp(a_e, *a_eCoarse);
  }

  m_print = a_print;
  relax(a_e, a_residual, a_iterations, a_AMRFASMGiter, a_depth);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::relax(LevelData<FArrayBox>&       a_e,
                                  const LevelData<FArrayBox>& a_residual,
                                  int                         a_iterations,
                                  int                         a_AMRFASMGiter,
                                  int                         a_depth)
{
  CH_TIME("AMRNonLinearPoissonOp::relax");
  if (m_verbosity > 3) {
      pout() << "AMRNonLinearPoissonOp::relax depth = "<< a_depth <<" \n"; 
  } 

  for (int i = 0; i < a_iterations; i++)
    {
      if (m_verbosity > 3) {
          pout() <<" levelGSRB "<< i <<"\n";
      }
      switch (s_relaxMode)
        {
        case 0:
          looseGSRB(a_e, a_residual);
          break;
        case 1:
          levelGSRB(a_e, a_residual, i, a_AMRFASMGiter, a_depth);
          break;
        case 2:
          overlapGSRB(a_e, a_residual);
          break;
        case 3:
          levelGSRBLazy(a_e, a_residual);
          break;
        case 4:
          levelJacobi(a_e, a_residual);
          break;
        case 5:
          levelMultiColor(a_e, a_residual);
          break;
        case 6:
          levelGS(a_e, a_residual);
          break;
        default:
          MayDay::Abort("unrecognized relaxation mode");
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::createCoarser(LevelData<FArrayBox>&       a_coarse,
                                          const LevelData<FArrayBox>& a_fine,
                                          bool                        a_ghosted)
{
  CH_TIME("AMRNonLinearPoissonOp::createCoarser");

  // CH_assert(!a_ghosted);
  IntVect ghost = a_fine.ghostVect();

  CH_assert(a_fine.disjointBoxLayout().coarsenable(2));
  if (m_coarsenedMGrids.size() == 0)
    coarsen(m_coarsenedMGrids, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
  a_coarse.define(m_coarsenedMGrids, a_fine.nComp(), ghost);
}


void AMRNonLinearPoissonOp::restrictR(LevelData<FArrayBox>& a_phiCoarse,
                                      const LevelData<FArrayBox>& a_phiFine)
{
  //    a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox&       phiFine = a_phiFine[dit];
    FArrayBox&       phiCoarse = a_phiCoarse[dit];

    Box region = dblFine.get(dit());
    const IntVect& iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    phiCoarse.setVal(0.0);

    // m_dx not used here
    FORT_RESTRICTNL(CHF_FRA_SHIFT(phiCoarse, civ),
                    CHF_CONST_FRA_SHIFT(phiFine, iv),
                    CHF_BOX_SHIFT(region, iv),
                    CHF_CONST_REAL(m_dx));
  }
}


// ---------------------------------------------------------
void AMRNonLinearPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                             LevelData<FArrayBox>&       a_phiFine,
                                             const LevelData<FArrayBox>& a_rhsFine)
{
  // default impl
  restrictResidual(a_resCoarse, a_phiFine, nullptr, a_rhsFine, true);
}


void AMRNonLinearPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                             LevelData<FArrayBox>&       a_phiFine,
                                             const LevelData<FArrayBox>* a_phiCoarse,
                                             const LevelData<FArrayBox>& a_rhsFine,
                                             bool homogeneous)
{
  CH_TIME("AMRNonLinearPoissonOp::restrictResidual");

  if (a_phiCoarse != nullptr) {
      m_interpWithCoarser.coarseFineInterp(a_phiFine, *a_phiCoarse);
  }

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit) {
      FArrayBox &phi = a_phiFine[dit];
      m_bc(phi, dblFine[dit()], m_domain, m_dx_vect, homogeneous);
  }

  a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  LevelData<FArrayBox>  a_nlfunc(dblFine, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dblFine, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phiFine,
                                        *m_B, *m_iceMask, *m_Pi, *m_zb);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit) {
      FArrayBox&       phi = a_phiFine[dit];
      const FArrayBox& rhs = a_rhsFine[dit];
      FArrayBox&       res = a_resCoarse[dit];

      FArrayBox&       nlfunc = a_nlfunc[dit];
        
      Box region = dblFine[dit];
      const IntVect& iv = region.smallEnd();
      IntVect civ = coarsen(iv, 2);
        
      res.setVal(0.0);
        
      FORT_RESTRICTRESNL(CHF_FRA_SHIFT(res, civ),
                         CHF_CONST_FRA_SHIFT(phi, iv),
                         CHF_CONST_FRA_SHIFT(rhs, iv),
                         CHF_CONST_REAL(m_alpha),
                         CHF_CONST_REAL(m_beta),
                         CHF_CONST_FRA_SHIFT(nlfunc, iv),
                         CHF_BOX_SHIFT(region, iv),
                         CHF_CONST_REALVECT(m_dx_vect));
  }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                             const LevelData<FArrayBox>& a_correctCoarse)
{
  CH_TIME("AMRNonLinearPoissonOp::prolongIncrement");
  if (m_verbosity > 3) {
      pout() << "AMRNonLinearPoissonOp::prolongIncrement\n";
  }
  
  DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
  int mgref = 2; //this is a multigrid func
  DataIterator dit = a_phiThisLevel.dataIterator();
  int nbox=dit.size();
  
#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_phiThisLevel[dit[ibox]];
        const FArrayBox& coarse = a_correctCoarse[dit[ibox]];
        Box region = dbl[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ=coarsen(iv, 2);
        
        FORT_PROLONGNL(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(mgref));
      }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRResidual(LevelData<FArrayBox>&              a_residual,
                                        const LevelData<FArrayBox>&        a_phiFine,
                                        const LevelData<FArrayBox>&        a_phi,
                                        const LevelData<FArrayBox>&        a_phiCoarse,
                                        const LevelData<FArrayBox>&        a_rhs,
                                        bool                               a_homogeneousPhysBC,
                                        AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRResidual");

  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
              a_homogeneousPhysBC, a_finerOp);

  axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRResidualNC(LevelData<FArrayBox>&              a_residual,
                                          const LevelData<FArrayBox>&        a_phiFine,
                                          const LevelData<FArrayBox>&        a_phi,
                                          const LevelData<FArrayBox>&        a_rhs,
                                          bool                               a_homogeneousPhysBC,
                                          AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRResidualNC");

   AMROperatorNC(a_residual, a_phiFine, a_phi,
                 a_homogeneousPhysBC, a_finerOp);

   axby(a_residual, a_residual, a_rhs, -1.0, 1.0);

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                                          const LevelData<FArrayBox>& a_phi,
                                          const LevelData<FArrayBox>& a_phiCoarse,
                                          const LevelData<FArrayBox>& a_rhs,
                                          bool                        a_homogeneousPhysBC)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRResidualNF");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined()) {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  }

  //apply boundary conditions
  this->residualI(a_residual, a_phi, a_rhs, a_homogeneousPhysBC);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMROperator(LevelData<FArrayBox>&              a_LofPhi,
                                        const LevelData<FArrayBox>&        a_phiFine,
                                        const LevelData<FArrayBox>&              a_phi,
                                        const LevelData<FArrayBox>&        a_phiCoarse,
                                        bool                               a_homogeneousPhysBC,
                                        AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMROperator");

  CH_assert(a_phi.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined()) {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  }

  // apply physical boundary conditions in applyOpI
  applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);

  if (a_phiFine.isDefined())
  {
    CH_assert(a_finerOp != NULL);
    reflux(a_phiFine, a_phi, a_LofPhi, a_finerOp);
  }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMROperatorNC(LevelData<FArrayBox>&              a_LofPhi,
                                          const LevelData<FArrayBox>&        a_phiFine,
                                          const LevelData<FArrayBox>&        a_phi,
                                          bool                               a_homogeneousPhysBC,
                                          AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIME("AMRNonLinearPoissonOp::AMROperatorNC");

  CH_assert(a_phi.isDefined());

  // apply physical boundary conditions in applyOpI
  applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);

  if (a_phiFine.isDefined())
  {
    CH_assert(a_finerOp != NULL);
    reflux(a_phiFine, a_phi, a_LofPhi, a_finerOp);
  }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                                          const LevelData<FArrayBox>& a_phi,
                                          const LevelData<FArrayBox>& a_phiCoarse,
                                          bool                        a_homogeneousPhysBC)
{
  CH_TIME("AMRNonLinearPoissonOp::AMROperatorNF");

  CH_assert(a_phi.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  if (a_phiCoarse.isDefined()) {
      m_interpWithCoarser.coarseFineInterp(phi, a_phiCoarse);
  }

  // apply physical boundary conditions in applyOpI
  this->applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                                        const LevelData<FArrayBox>& a_residual,
                                        const LevelData<FArrayBox>& a_correction,
                                        const LevelData<FArrayBox>& a_coarseCorrection,
                                        bool a_skip_res )
{
  CH_TIME("AMRNonLinearPoissonOp::AMRRestrict");

  LevelData<FArrayBox> r;
  create(r, a_residual);

  AMRRestrictS(a_resCoarse, a_residual, a_correction, a_coarseCorrection, r, a_skip_res );
}

// hacks to do BCs right
// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRRestrictS(LevelData<FArrayBox>&       a_resCoarse, // output
                                         const LevelData<FArrayBox>& a_residual,  // input
                                         const LevelData<FArrayBox>& a_correction, // for residual
                                         const LevelData<FArrayBox>& a_coarseCorrection, // C-F interp.
                                         LevelData<FArrayBox>&       a_scratch,   // temp buffer
                                         bool a_skip_res // flag skipping computing residual, used for lhs restrictions (FAS)
                                ) 
{
  CH_TIME("AMRNonLinearPoissonOp::AMRRestrictS");

  if ( !a_skip_res )
    {
      AMRResidualNF( a_scratch, a_correction, a_coarseCorrection, a_residual, false );
    }
  else 
    {
      // just copy data (phi in this case, even if its called residual)
      assignLocal( a_scratch, a_residual );
    }

  DisjointBoxLayout dblCoar = a_resCoarse.disjointBoxLayout();

#pragma omp parallel 
  {
    DataIterator dit = a_residual.dataIterator();
    int nbox=dit.size();
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& coarse = a_resCoarse[dit[ibox]];
        const FArrayBox& fine = a_scratch[dit[ibox]];
        const Box& b = dblCoar[dit[ibox]];
        Box refbox(IntVect::Zero,
                   (m_refToCoarser-1)*IntVect::Unit);
        FORT_AVERAGE( CHF_FRA(coarse),
                      CHF_CONST_FRA(fine),
                      CHF_BOX(b),
                      CHF_CONST_INT(m_refToCoarser),
                      CHF_BOX(refbox)
                      );
      }
  }//end pragma
}

// ---------------------------------------------------------
/** a_correction += I[2h->h](a_coarseCorrection) */
void AMRNonLinearPoissonOp::AMRProlong(LevelData<FArrayBox>&       a_correction,
                                       const LevelData<FArrayBox>& a_coarseCorrection)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRProlong");

  DisjointBoxLayout c;
  coarsen(c, a_correction.disjointBoxLayout(), m_refToCoarser);

  LevelData<FArrayBox> eCoar(c, a_correction.nComp(), a_coarseCorrection.ghostVect());
  a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  DataIterator dit = a_correction.dataIterator();
  int nbox=dit.size();

#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_correction[dit[ibox]];
        const FArrayBox& coarse = eCoar[dit[ibox]];
        
        Box region = dbl[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ = coarsen(iv, m_refToCoarser);
        
        FORT_PROLONGNL(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(m_refToCoarser));
      }
  }// end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRProlongS(LevelData<FArrayBox>&       a_correction,
                                        const LevelData<FArrayBox>& a_coarseCorrection,
                                        LevelData<FArrayBox>&       a_temp,
                                        const Copier&               a_copier)
{
  CH_TIME("AMRNonLinearPoissonOp::AMRProlongS");

  a_coarseCorrection.copyTo(a_temp.interval(), a_temp, a_temp.interval(), a_copier);

  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  DataIterator dit = a_correction.dataIterator();
  int nbox=dit.size();
  
#pragma omp parallel 
  {
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_correction[dit[ibox]];
        const FArrayBox& coarse = a_temp[dit[ibox]];
        
        Box region = dbl[dit[ibox]];
        const IntVect& iv =  region.smallEnd();
        IntVect civ= coarsen(iv, m_refToCoarser);
        
        FORT_PROLONGNL(CHF_FRA_SHIFT(phi, iv),
                     CHF_CONST_FRA_SHIFT(coarse, civ),
                     CHF_BOX_SHIFT(region, iv),
                     CHF_CONST_INT(m_refToCoarser));
      }
  }//end pragma
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRProlongS_2(LevelData<FArrayBox>&       a_correction,
                                         const LevelData<FArrayBox>& a_coarseCorrection,
                                         LevelData<FArrayBox>&       a_temp,
                                         const Copier&               a_copier,
                                         const Copier&         a_cornerCopier,
                                         const AMRLevelOp<LevelData<FArrayBox> >* a_crsOp )
{
  CH_TIME("AMRNonLinearPoissonOp::AMRProlongS_2");
  
  DisjointBoxLayout dbl = a_correction.disjointBoxLayout();
  DisjointBoxLayout cdbl = a_temp.disjointBoxLayout();
  AMRNonLinearPoissonOp* coarserAMRPOp = (AMRNonLinearPoissonOp*) a_crsOp;
  
  a_coarseCorrection.copyTo( a_temp.interval(), a_temp, a_temp.interval(), a_copier );

  // I think we should be using a coarse data iterator for applying the coarse BC?
  for (DataIterator cdit = a_temp.dataIterator(); cdit.ok(); ++cdit)
  {
    FArrayBox& coarse = a_temp[cdit];
    if (!m_use_FAS) {
        coarserAMRPOp->m_bc( coarse, cdbl[cdit], coarserAMRPOp->m_domain, coarserAMRPOp->m_dx_vect, true );
    } else {
        coarserAMRPOp->m_bc( coarse, cdbl[cdit], coarserAMRPOp->m_domain, coarserAMRPOp->m_dx_vect, false );
    }
  }

  // The corner copier passed in as an argument doesn't always work, whilst this one seems better
  CornerCopier cornerCopy(cdbl, cdbl, coarserAMRPOp->m_domain, a_temp.ghostVect(), true);

  a_temp.exchange(cornerCopy);
  //a_temp.exchange( a_temp.interval(), a_cornerCopier ); -- needed for AMR
#pragma omp parallel
  {
    DataIterator dit = a_correction.dataIterator();
    int nbox = dit.size();
#pragma omp for 
    for(int ibox = 0; ibox < nbox; ibox++)
      {
        FArrayBox& phi =  a_correction[dit[ibox]];
        FArrayBox& coarse = a_temp[dit[ibox]];
        
        Box region = dbl[dit[ibox]];
        const IntVect& iv = region.smallEnd();
        IntVect civ = coarsen(iv, m_refToCoarser);
        
#if 0
        FORT_PROLONGNL( CHF_FRA_SHIFT(phi, iv),
                      CHF_CONST_FRA_SHIFT(coarse, civ),
                      CHF_BOX_SHIFT(region, iv),
                      CHF_CONST_INT(m_refToCoarser));
#else
        FORT_PROLONG_2_NL( CHF_FRA_SHIFT(phi, iv),
                        CHF_CONST_FRA_SHIFT(coarse, civ),
                        CHF_BOX_SHIFT(region, iv),
                        CHF_CONST_INT(m_refToCoarser) 
                        );
#endif
      }
  }//end pragma

  // debug
  //write(&a_temp,"z_src.hdf5"); 
  //write(&a_correction,"z_prol2.hdf5"); exit(12);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                                              const LevelData<FArrayBox>& a_correction,
                                              const LevelData<FArrayBox>& a_coarseCorrection)
{
//   LevelData<FArrayBox> r;
//   this->create(r, a_residual);
//   this->AMRResidualNF(r, a_correction, a_coarseCorrection, a_residual, true);
//   this->assign(a_residual, r);
  this->AMRResidualNF(a_residual, a_correction, a_coarseCorrection, a_residual, false);
}

// ---------------------------------------------------------
///compute norm over all cells on coarse not covered by finer
Real AMRNonLinearPoissonOp::AMRNorm(const LevelData<FArrayBox>& a_coarResid,
                                    const LevelData<FArrayBox>& a_fineResid,
                                    const int&                  a_refRat,
                                    const int&                  a_ord)

{
  CH_TIME("AMRNonLinearPoissonOp::AMRNorm");

  //create temp and zero out under finer grids
  LevelData<FArrayBox> coarTemp;
  m_levelOps.create(coarTemp, a_coarResid);
  m_levelOps.assign(coarTemp, a_coarResid);

  if (a_fineResid.isDefined())
  {
    const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
    const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();

    int ncomp = coarTemp.nComp();

    for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
      {
        FArrayBox& coarTempFAB = coarTemp[dit];
        LayoutIterator litFine = fineGrids.layoutIterator();

        for (litFine.reset(); litFine.ok(); ++litFine)
          {
            Box overlayBox = coarTempFAB.box();
            Box coarsenedGrid = coarsen(fineGrids[litFine], a_refRat);

            overlayBox &= coarsenedGrid;

            if (!overlayBox.isEmpty())
              {
                coarTempFAB.setVal(0.0, overlayBox, 0, ncomp);
              }
          }
      }
  }

  // return norm of temp
  return norm(coarTemp, a_ord);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::setAlphaAndBeta(const Real& a_alpha,
                                   const Real& a_beta)
{
  m_alpha = a_alpha; 
  m_beta  = a_beta;  
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::setBC(const BCHolder& a_bc)
{
  m_bc = a_bc;
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                                   const LevelData<FArrayBox>&        a_phi,
                                   LevelData<FArrayBox>&              a_residual,
                                   AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIMERS("AMRNonLinearPoissonOp::reflux");
  CH_TIMER("AMRNonLinearPoissonOp::reflux::incrementCoarse", t2);

  m_levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  CH_START(t2);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab = a_phi[dit];

      if (m_levfluxreg.hasCF(dit()))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox coarflux;
              getFlux(coarflux, coarfab, idir);

              Real scale = 1.0;
              for (int i=0; i<SpaceDim; ++i) { 
                  if (idir != i) {
                      scale *= m_dx_vect[i];
                  }
              }
              m_levfluxreg.incrementCoarse(coarflux, scale, dit(),
                                           interv, interv, idir);
            }
        }
    }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = ( LevelData<FArrayBox>&)a_phiFine;

  AMRNonLinearPoissonOp* finerAMRPOp = (AMRNonLinearPoissonOp*) a_finerOp;
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(phiFineRef, a_phi);
  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // phiFineRef.exchange(a_phiFine.interval());
  IntVect phiGhost = phiFineRef.ghostVect();
  int ncomps = a_phiFine.nComp();

  CH_TIMER("AMRNonLinearPoissonOp::reflux::incrementFine", t3);
  CH_START(t3);

  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const Box& gridbox = dblFine[ditf];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              if (m_levfluxreg.hasCF(ditf(), sit()))
                {
                  Side::LoHiSide hiorlo = sit();
                  Box fluxBox = bdryBox(gridbox,idir,hiorlo,1);

                  FArrayBox fineflux(fluxBox,ncomps);
                  getFlux(fineflux, phifFab, fluxBox, idir, m_refToFiner);

                  Real scale = 1.0;
                  for (int i=0; i<SpaceDim; ++i) { 
                      if (idir != i) {
                          scale *= m_dx_vect[i];
                      }
                  }
                  m_levfluxreg.incrementFine(fineflux, scale, ditf(),
                                             interv, interv, idir, hiorlo);
                }
            }
        }
    }

  CH_STOP(t3);

  Real scale = 1.0;
  for (int i=0; i<SpaceDim; ++i) {
      scale *= m_dx_vect[i];
  }
  scale = 1.0/scale;
  m_levfluxreg.reflux(a_residual, scale);
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::write(const LevelData<FArrayBox>* a_data,
                         const char*                 a_filename)
{
#ifdef CH_USE_HDF5
  writeLevelname(a_data, a_filename);
#else
  MayDay::Warning("AMRNonLinearPoissonOp::write unimplemented since CH_USE_HDF5 undefined");
#endif
}

/***/
// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelGSRB( LevelData<FArrayBox>&       a_phi,
                                       const LevelData<FArrayBox>& a_rhs,
                                       int                         a_ite,
                                       int                         a_AMRFASMGiter,
                                       int                         a_depth )
{
  CH_TIME("AMRNonLinearPoissonOp::levelGSRB");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                        *m_B, *m_iceMask, *m_Pi, *m_zb);

  bool a_homo = false;

  DataIterator dit = a_phi.dataIterator();
  int nbox=dit.size();
  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {
      CH_TIME("AMRNonLinearPoissonOp::levelGSRB::Compute");

      // fill in intersection of ghostcells and a_phi's boxes
      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::homogeneousCFInterp");
        // For FAS, we don't want to do this?
        if (!m_use_FAS) {
            homogeneousCFInterp(a_phi);
            a_homo = true;
        }
      }

      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::exchange");
        if (s_exchangeMode == 0)
          a_phi.exchange( a_phi.interval(), m_exchangeCopier );
        else if (s_exchangeMode == 1)
          a_phi.exchangeNoOverlap(m_exchangeCopier);
        else
          MayDay::Abort("exchangeMode");
      }
#pragma omp parallel
      {
#pragma omp for
        for (int ibox=0; ibox < nbox; ibox++)
          {
            const Box& region = dbl[dit[ibox]];
            FArrayBox& phiFab = a_phi[dit[ibox]];
            
            // if you're not homogeneous and not FAS you should have done something to end up homogeneous
            m_bc( phiFab, region, m_domain, m_dx_vect, a_homo );
            
            if (m_alpha == 0.0 && m_beta == 1.0 )
              {
                FORT_GSRBLAPLACIANFUNCNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REALVECT(m_dx_vect),
                                   CHF_CONST_INT(whichPass),
                                   CHF_CONST_FRA(a_nlfunc[dit[ibox]]),
                                   CHF_CONST_FRA(a_nlDfunc[dit[ibox]]));
              }
            else
              {
                FORT_GSRBHELMHOLTZFUNCNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REALVECT(m_dx_vect),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(m_beta),
                                   CHF_CONST_FRA(a_nlfunc[dit[ibox]]),
                                   CHF_CONST_FRA(a_nlDfunc[dit[ibox]]),
                                   CHF_CONST_INT(whichPass));
              }
          } // end loop through grids
      }//end pragma
    } // end loop through red-black
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelMultiColor(LevelData<FArrayBox>&       a_phi,
                                            const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Error("AMRNonLinearPoissonOp::levelMultiColor not implemented");

}


// ---------------------------------------------------------
void AMRNonLinearPoissonOp::looseGSRB(LevelData<FArrayBox>&       a_phi,
                                      const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::looseGSRB");

  MayDay::Error("AMRNonLinearPoissonOp::looseGSRB not implemented");

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::overlapGSRB(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::overlapGSRB");
  MayDay::Error("AMRNonLinearPoissonOp::overlapGSRB not implemented");
}

/***/
// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelGSRBLazy( LevelData<FArrayBox>&       a_phi,
                                           const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::levelGSRBLazy");
  MayDay::Error("AMRNonLinearPoissonOp::levelGSRBLazy not implemented");

}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                                        const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRNonLinearPoissonOp::levelJacobi");
  MayDay::Error("AMRNonLinearPoissonOp::levelJacobi not implemented");
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::levelGS( LevelData<FArrayBox>&       a_phi,
                                     const LevelData<FArrayBox>& a_rhs )
{

  CH_TIME("AMRNonLinearPoissonOp::levelGS");

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  const DisjointBoxLayout& dbl = a_rhs.disjointBoxLayout();

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                        *m_B, *m_iceMask, *m_Pi, *m_zb);

  DataIterator dit = a_phi.dataIterator();
  int nbox=dit.size();
  // do first red, then black passes
  int whichPass = 0;

      //CH_TIME("AMRNonLinearPoissonOp::levelGSRB::Compute");

      // fill in intersection of ghostcells and a_phi's boxes
      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::homogeneousCFInterp");
//        homogeneousCFInterp(a_phi);
      }

      {
        CH_TIME("AMRNonLinearPoissonOp::levelGSRB::exchange");
        if (s_exchangeMode == 0)
          a_phi.exchange( a_phi.interval(), m_exchangeCopier );
        else if (s_exchangeMode == 1)
          a_phi.exchangeNoOverlap(m_exchangeCopier);
        else
          MayDay::Abort("exchangeMode");
      }
#pragma omp parallel
      {
#pragma omp for
        for (int ibox=0; ibox < nbox; ibox++)
          {
            const Box& region = dbl[dit[ibox]];
            FArrayBox& phiFab = a_phi[dit[ibox]];

            m_bc( phiFab, region, m_domain, m_dx_vect, true );

            if (m_alpha == 0.0 && m_beta == 1.0 )
              {
                FORT_GSLAPLACIANFUNCNL(CHF_FRA(phiFab),
                                   CHF_CONST_FRA(a_rhs[dit[ibox]]),
                                   CHF_BOX(region),
                                   CHF_CONST_REALVECT(m_dx_vect),
                                   CHF_CONST_INT(whichPass),
                                   CHF_CONST_FRA(a_nlfunc[dit[ibox]]),
                                   CHF_CONST_FRA(a_nlDfunc[dit[ibox]]));

              }
            else
              {
              MayDay::Error("AMRNonLinearPoissonOp: FORT_GSHELMHOLTZNL not implemented");
              }
          } // end loop through grids
      }//end pragma

}




// ---------------------------------------------------------
void AMRNonLinearPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif)
{
  CH_TIME("AMRNonLinearPoissonOp::CFInterp");

  CH_assert( a_phif.ghostVect() >= IntVect::Unit);

  DataIterator dit = a_phif.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const DataIndex& datInd = dit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              homogeneousCFInterp(a_phif,datInd,idir,sit());
            }
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                                const DataIndex&      a_datInd,
                                                int                   a_idir,
                                                Side::LoHiSide        a_hiorlo)
{
  // CH_TIME("AMRNonLinearPoissonOp::homogeneousCFInterp");

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  //  CH_assert (m_ncomp == a_phif.nComp());

  const CFIVS* cfivs_ptr = NULL;
  if (a_hiorlo == Side::Lo)
    cfivs_ptr = &(m_cfregion.loCFIVS(a_datInd, a_idir));
  else
    cfivs_ptr = &(m_cfregion.hiCFIVS(a_datInd, a_idir));

  if (cfivs_ptr->isPacked())
    {
      int ihiorlo = sign(a_hiorlo);
      FArrayBox& phiFab = a_phif[a_datInd];
      if (phiFab.box().size(a_idir) == 3)
        {
          FORTNT_INTERPHOMOLINEAR(CHF_FRA(phiFab),
                                  CHF_BOX(cfivs_ptr->packedBox()),
                                  CHF_CONST_REAL(m_dx_vect[a_idir]),
                                  CHF_CONST_REAL(m_dxCrse_vect[a_idir]),
                                  CHF_CONST_INT(a_idir),
                                  CHF_CONST_INT(ihiorlo));
        }
      else
        {
          FORTNT_INTERPHOMO(CHF_FRA(phiFab),
                            CHF_BOX(cfivs_ptr->packedBox()),
                            CHF_CONST_REAL(m_dx_vect[a_idir]),
                            CHF_CONST_REAL(m_dxCrse_vect[a_idir]),
                            CHF_CONST_INT(a_idir),
                            CHF_CONST_INT(ihiorlo));
        }
    }
  else
    {
      const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();
      if (!interp_ivs.isEmpty())
        {
          // Assuming homogenous, interpolate on fine ivs
          interpOnIVSHomo(a_phif, a_datInd, a_idir,
                          a_hiorlo, interp_ivs);
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::singleBoxCFInterp(FArrayBox& a_phi)
{
  CH_TIME("AMRNonLinearPoissonOp::singleBoxCFInterp");

  Box region = a_phi.box();
  region.grow(-1);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next())
        {
          Box edge = adjCellBox(region, idir, sit(), 1);
          int ihiorlo = sign(sit());
          FORT_INTERPHOMO(CHF_FRA(a_phi),
                      CHF_BOX(edge),
                      CHF_CONST_REAL(m_dx_vect[idir]),
                      CHF_CONST_REAL(m_dxCrse_vect[idir]),
                      CHF_CONST_INT(idir),
                      CHF_CONST_INT(ihiorlo));
        }
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                                            const DataIndex&      a_datInd,
                                            const int             a_idir,
                                            const Side::LoHiSide  a_hiorlo,
                                            const IntVectSet&     a_interpIVS)
{
  // CH_TIME("AMRNonLinearPoissonOp::interpOnIVSHomo");

  CH_assert((a_idir >= 0) && (a_idir  < SpaceDim));
  CH_assert((a_hiorlo == Side::Lo )||(a_hiorlo == Side::Hi ));
  CH_assert( a_phif.ghostVect() >= IntVect::Unit);
  IVSIterator fine_ivsit(a_interpIVS);
  FArrayBox& a_phi = a_phif[a_datInd];
  int ihilo = sign(a_hiorlo);

  if (a_phi.box().size(a_idir)==3) // we are in a 1-wide  box
    {
      IntVect iv;
      Real pa;
      Real factor = 1-2*m_dx_vect[a_idir]/(m_dx_vect[a_idir]+m_dxCrse_vect[a_idir]);
      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          iv = fine_ivsit();
          iv[a_idir]-=ihilo;
          // linear interpolation
          for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
            {
              pa = a_phi(iv, ivar);
              a_phi(fine_ivsit(), ivar) = factor*pa;
            }
        }
    }
  else if (false) // the old expensive way of computing this  bvs
    {

      // much of these scalar values can be precomputed and stored if
      // we ever need to speed-up this function (ndk)
      Real x1 = m_dx_vect[a_idir];
      Real x2 = 0.5*(3. * x1 + m_dxCrse_vect[a_idir]);
      Real denom = 1.0-((x1+x2)/x1);
      Real idenom = 1/(denom); // divide is more expensive usually
      Real x = 2.*x1;
      Real xsquared = x*x;

      Real m1 = 1/(x1*x1);
      Real m2 = 1/(x1*(x1-x2));

      Real q1 = 1/(x1-x2);
      Real q2 = x1+x2;

      Real pa,pb,a,b;
      IntVect ivf;
      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          ivf = fine_ivsit();
          // quadratic interpolation
          for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
            {
              ivf[a_idir]-=2*ihilo;
              pa = a_phi(ivf, ivar);
              ivf[a_idir]+=ihilo;
              pb = a_phi(ivf, ivar);

              a = ((pb-pa)*m1 - (pb)*m2)*idenom;
              b = (pb)*q1 - a*q2;

              ivf[a_idir]+=ihilo;
              a_phi(fine_ivsit(), ivar) = a*xsquared + b*x + pa;

            } //end loop over components
        } //end loop over fine intvects
    }
  else // symbolic reduced version of CF quadratic stencil
    {
      Real pa, pb;
      Real c1 = 2*(m_dxCrse_vect[a_idir]-m_dx_vect[a_idir])/(m_dxCrse_vect[a_idir]+  m_dx_vect[a_idir]); //first inside point
      Real c2 =  -(m_dxCrse_vect[a_idir]-m_dx_vect[a_idir])/(m_dxCrse_vect[a_idir]+3*m_dx_vect[a_idir]); // next point inward
      IntVect ivf;
      for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit)
        {
          ivf = fine_ivsit();
          // quadratic interpolation
          for (int ivar = 0; ivar < a_phif.nComp(); ivar++)
            {
              ivf[a_idir]-=2*ihilo;
              pa = a_phi(ivf, ivar);
              ivf[a_idir]+=ihilo;
              pb = a_phi(ivf, ivar);

              ivf[a_idir]+=ihilo;
              a_phi(fine_ivsit(), ivar) = c1*pb + c2*pa;

            } //end loop over components
        } //end loop over fine intvects
    }
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::getFlux(FArrayBox&       a_flux,
                                    const FArrayBox& a_data,
                                    const Box&       a_edgebox,
                                    int              a_dir,
                                    int              a_ref) const
{
  // In this version of getFlux, the edgebox is passed in, and the flux array
  // is already defined.

  CH_TIME("AMRNonLinearPoissonOp::getFlux1");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!a_edgebox.isEmpty());

  Real scale = m_beta * a_ref / m_dx_vect[a_dir];

  FORT_NEWGETFLUXNL(CHF_FRA(a_flux),
                  CHF_CONST_FRA(a_data),
                  CHF_BOX(a_edgebox),
                  CHF_CONST_REAL(scale),
                  CHF_CONST_INT(a_dir));
}

// ---------------------------------------------------------
void AMRNonLinearPoissonOp::getFlux(FArrayBox&       a_flux,
                                    const FArrayBox& a_data,
                                    int              a_dir,
                                    int              a_ref) const
{
  CH_TIME("AMRNonLinearPoissonOp::getFlux2");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());

  Box edgebox = surroundingNodes(a_data.box(),a_dir);
  edgebox.grow(a_dir, -1);

  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!edgebox.isEmpty());

  a_flux.define(edgebox, a_data.nComp());

  //FArrayBox fflux(edgebox, a_data.nComp());
  Real scale = m_beta * a_ref / m_dx_vect[a_dir];

  FORT_NEWGETFLUXNL(CHF_FRA(a_flux),
                  CHF_CONST_FRA(a_data),
                  CHF_BOX(edgebox),
                  CHF_CONST_REAL(scale),
                  CHF_CONST_INT(a_dir));
}

// Factory

// ---------------------------------------------------------
//  AMR Factory define function
void AMRNonLinearPoissonOpFactory::define(const ProblemDomain&             a_coarseDomain,
                                          const Vector<DisjointBoxLayout>& a_grids,
                                          const Vector<int>&               a_refRatios,
                                          const RealVect&                  a_coarsedx,
                                          BCHolder                         a_bc,
                                          AmrHydro*                        a_amrHydro,
                                          NL_level                         a_nllevel,
                                          Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_B,
                                          Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_Pi,
                                          Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_zb,
                                          Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_IM,
                                          Real                             a_alpha,
                                          Real                             a_beta)
{
  CH_TIME("AMRNonLinearPoissonOpFactory::define");

  m_boxes = a_grids;

  m_refRatios = a_refRatios;

  m_bc = a_bc;

  m_dx.resize(a_grids.size());
  m_dx[0] = a_coarsedx;
  //D_TERM(m_dx[0][0] = a_coarsedx[0];, m_dx[0][1] = a_coarsedx[1];, m_dx[0][2] = a_coarsedx[2];)

  m_domains.resize(a_grids.size());
  m_domains[0] = a_coarseDomain;

  m_exchangeCopiers.resize(a_grids.size());
  m_exchangeCopiers[0].exchangeDefine(a_grids[0], IntVect::Unit);
  m_exchangeCopiers[0].trimEdges(a_grids[0], IntVect::Unit);

  m_cfregion.resize(a_grids.size());
  m_cfregion[0].define(a_grids[0], m_domains[0]);

  for (int i = 1; i < a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1] / m_refRatios[i-1];
      //D_TERM(m_dx[i][0] = m_dx[i-1][0] / m_refRatios[i-1];, 
      //       m_dx[i][1] = m_dx[i-1][1] / m_refRatios[i-1];,
      //       m_dx[i][2] = m_dx[i-1][2] / m_refRatios[i-1];)

      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);

      if (a_grids[i].isClosed())
        {
          m_exchangeCopiers[i].exchangeDefine(a_grids[i], IntVect::Unit);
          m_exchangeCopiers[i].trimEdges(a_grids[i], IntVect::Unit);
          m_cfregion[i].define(a_grids[i], m_domains[i]);
        }
    }

  m_alpha = a_alpha;
  m_beta = a_beta;

  m_amrHydro = a_amrHydro;
  m_nllevel  = a_nllevel;

  m_verbosity = 3;
  m_print     = false;

  m_use_FAS = true; // May need to make this a define arg ... but I'm only working with FAS usually

  m_B  = a_B;  // Gap Height
  m_Pi = a_Pi; // Overb Press  
  m_zb = a_zb; // Bed Elevation
  m_iceMask = a_IM; // Ice Mask
}

// ---------------------------------------------------------
// MultiGrid define function -- totally not going to work with NL func
void AMRNonLinearPoissonOpFactory::define(const ProblemDomain&     a_domain,
                                          const DisjointBoxLayout& a_grid,
                                          const RealVect&          a_dx,
                                          BCHolder                 a_bc,
                                          int                      a_maxDepth,
                                          Real                     a_alpha,
                                          Real                     a_beta)
{
  AmrHydro* amrHydro;
  NL_level  nllevel; 

  Vector<DisjointBoxLayout> grids(1);
  grids[0] = a_grid;
  Vector<int> refRatio(1, 2);

  Vector<RefCountedPtr<LevelData<FArrayBox> > >  B(grids.size());  // Gap Height
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  Pri(grids.size());// Overb Press  
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  zb(grids.size()); // Bed Elevation
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  iceMask(grids.size()); // Ice Mask
  for (int i = 0; i < grids.size(); ++i) {
      B[i]   = RefCountedPtr<LevelData<FArrayBox> >(
                   new LevelData<FArrayBox>(grids[i], 1, IntVect::Unit));
      Pri[i] = RefCountedPtr<LevelData<FArrayBox> >(
                   new LevelData<FArrayBox>(grids[i], 1, IntVect::Unit));
      zb[i]  = RefCountedPtr<LevelData<FArrayBox> >(
                   new LevelData<FArrayBox>(grids[i], 1, IntVect::Unit));
      iceMask[i]  = RefCountedPtr<LevelData<FArrayBox> >(
                       new LevelData<FArrayBox>(grids[i], 1, IntVect::Unit));

      for (DataIterator dit = B[i]->dataIterator(); dit.ok(); ++dit) {
          (*B[i])[dit()].setVal(1.0);
          (*Pri[i])[dit()].setVal(1.0);
          (*zb[i])[dit()].setVal(1.0);
          (*iceMask[i])[dit()].setVal(1.0);
      }
  }

  define(a_domain, grids, refRatio, a_dx, a_bc, 
         amrHydro, nllevel, B, Pri, zb, iceMask, a_alpha, a_beta);
}

// ---------------------------------------------------------
MGLevelOp<LevelData<FArrayBox> >* AMRNonLinearPoissonOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                        int                  a_depth,
                                                                        bool                 a_homoOnly)
{
  CH_TIME("AMRNonLinearPoissonOpFactory::MGnewOp");

  RealVect dxCrse = -IntVect::Unit;

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++) {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) {
          break;
      }
  }

  CH_assert(ref != m_domains.size()); // didn't find domain

  if (ref > 0)
    {
      dxCrse = m_dx[ref-1];
    }

  ProblemDomain domain(m_domains[ref]);
  RealVect dx = m_dx[ref];
  int coarsening = 1;

  for (int i = 0; i < a_depth; i++)
    {
      coarsening *= 2;
      domain.coarsen(2);
    }

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*AMRNonLinearPoissonOp::s_maxCoarse))
  {
    return NULL;
  }

  dx *= coarsening;

  DisjointBoxLayout layout;
  coarsen_dbl(layout, m_boxes[ref], coarsening);

  Copier ex = m_exchangeCopiers[ref];
  CFRegion cfregion = m_cfregion[ref];

  if (coarsening > 1) {
      ex.coarsen(coarsening);
      cfregion.coarsen(coarsening);
  }

  AMRNonLinearPoissonOp* newOp = new AMRNonLinearPoissonOp;
  newOp->define(layout, dx, domain, m_bc, ex, cfregion);

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_verbosity = m_verbosity;
  newOp->m_print     = m_print;

  newOp->m_amrHydro  = m_amrHydro; 
  newOp->m_nllevel   = m_nllevel;

  if (a_depth == 0) {
      // Problem SPECIFIC
      newOp->m_B  = m_B[ref];  // Gap Height
      newOp->m_Pi = m_Pi[ref]; // Overb Press  
      newOp->m_zb = m_zb[ref]; // Bed Elevation
      newOp->m_iceMask = m_iceMask[ref]; // Ice Mask
  } else {
      // Problem SPECIFIC
      RefCountedPtr<LevelData<FArrayBox> > B( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FArrayBox> > Pri( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FArrayBox> > zb( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FArrayBox> > iceMask( new LevelData<FArrayBox> );
      B->define(layout,   m_B[ref]->nComp(),  m_B[ref]->ghostVect());
      Pri->define(layout, m_Pi[ref]->nComp(), m_Pi[ref]->ghostVect());
      zb->define(layout,  m_zb[ref]->nComp(), m_zb[ref]->ghostVect());
      iceMask->define(layout, m_iceMask[ref]->nComp(), m_iceMask[ref]->ghostVect());

      // average coefficients to coarser level
      // for now, do this with a CoarseAverage --
      // may want to switch to harmonic averaging at some point
      CoarseAverage averagerB(m_B[ref]->getBoxes(),
                             layout, B->nComp(), coarsening);
      CoarseAverage averagerPi(m_Pi[ref]->getBoxes(),
                             layout, Pri->nComp(), coarsening);
      CoarseAverage averagerZb(m_zb[ref]->getBoxes(),
                             layout, zb->nComp(), coarsening);
      CoarseAverage averagerIM(m_iceMask[ref]->getBoxes(),
                             layout, iceMask->nComp(), coarsening);

      averagerB.averageToCoarse(*B,   *(m_B[ref]));
      averagerPi.averageToCoarse(*Pri, *(m_Pi[ref]));
      averagerZb.averageToCoarse(*zb,  *(m_zb[ref]));
      averagerIM.averageToCoarse(*iceMask,  *(m_iceMask[ref]));

      newOp->m_B  = B;   // Gap Height
      newOp->m_Pi = Pri; // Overb Press  
      newOp->m_zb = zb;  // Bed Elevation
      newOp->m_iceMask = iceMask;  // Bed Elevation
    }

  newOp->m_dxCrse      = dxCrse[0];
  newOp->m_dxCrse_vect = dxCrse;

  newOp->m_use_FAS = m_use_FAS;

  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;
}

// ---------------------------------------------------------
AMRLevelOp<LevelData<FArrayBox> >* AMRNonLinearPoissonOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("AMRNonLinearPoissonOpFactory::AMRnewOp");

  AMRNonLinearPoissonOp* newOp = new AMRNonLinearPoissonOp;
  RealVect dxCrse = -IntVect::Unit;

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++) {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) {
          break;
      }
  }

  if (ref == 0) {
      // coarsest AMR level
      if ((m_domains.size() == 1) || (!m_boxes[1].isClosed()) ) {
          // no finer level
          newOp->define(m_boxes[0], m_dx[0],
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
      } else {
          // finer level exists but no coarser
          int dummyRat = 1;  // argument so compiler can find right function
          int refToFiner = m_refRatios[0]; // actual refinement ratio
          newOp->define(m_boxes[0], m_boxes[1], m_dx[0],
                        dummyRat, refToFiner,
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
      }
  } else if ((ref ==  m_domains.size()-1) || (!m_boxes[ref+1].isClosed())) {
      dxCrse = m_dx[ref-1];

      // finest AMR level
      newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[ref], m_cfregion[ref]);
  } else if ( ref == m_domains.size()) {
      MayDay::Abort("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");

  } else {
      dxCrse = m_dx[ref-1];

      // intermediate AMR level, full define
      newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                    m_refRatios[ref-1], m_refRatios[ref],
                    a_indexSpace, m_bc,
                    m_exchangeCopiers[ref], m_cfregion[ref]);
  }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_verbosity = m_verbosity;
  newOp->m_print     = m_print;

  newOp->m_amrHydro  = m_amrHydro; 
  newOp->m_nllevel   = m_nllevel;

  newOp->m_dxCrse      = dxCrse[0];
  newOp->m_dxCrse_vect = dxCrse;

  newOp->m_use_FAS = m_use_FAS;

  // Problem SPECIFIC
  newOp->m_B  = m_B[ref];  // Gap Height
  newOp->m_Pi = m_Pi[ref]; // Overb Press  
  newOp->m_zb = m_zb[ref]; // Bed Elevation
  newOp->m_iceMask = m_iceMask[ref]; // Ice Mask

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

//-----------------------------------------------------------------------
void
AMRNonLinearPoissonOp::finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                                              int a_coarseningFactor)
{
  const AMRNonLinearPoissonOp& op =
    dynamic_cast<const AMRNonLinearPoissonOp&>(a_operator);

  // Problem SPECIFIC
  LevelData<FArrayBox>& BCoar     = *m_B;  // Gap Height
  LevelData<FArrayBox>& PiCoar    = *m_Pi; // Overb Press
  LevelData<FArrayBox>& zbCoar    = *m_zb; // Bed Elevation
  LevelData<FArrayBox>& IMCoar    = *m_iceMask; // Bed Elevation

  // Problem SPECIFIC
  const LevelData<FArrayBox>& BFine  = *(op.m_B);
  const LevelData<FArrayBox>& PiFine = *(op.m_Pi);
  const LevelData<FArrayBox>& zbFine = *(op.m_zb);
  const LevelData<FArrayBox>& IMFine = *(op.m_iceMask);

  if (a_coarseningFactor != 1) {
      // B
      CoarseAverage cellAverageB(BFine.disjointBoxLayout(),
                                BCoar.disjointBoxLayout(),
                                1, a_coarseningFactor);
      for (DataIterator dit = BCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        BCoar[dit()].setVal(0.);
      cellAverageB.averageToCoarse(BCoar, BFine);
      // Pi
      CoarseAverage cellAveragePi(PiFine.disjointBoxLayout(),
                                PiCoar.disjointBoxLayout(),
                                1, a_coarseningFactor);
      for (DataIterator dit = PiCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        PiCoar[dit()].setVal(0.);
      cellAveragePi.averageToCoarse(PiCoar, PiFine);
      // zb
      CoarseAverage cellAverageZb(zbFine.disjointBoxLayout(),
                                zbCoar.disjointBoxLayout(),
                                1, a_coarseningFactor);
      for (DataIterator dit = zbCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        zbCoar[dit()].setVal(0.);
      cellAverageZb.averageToCoarse(zbCoar, zbFine);
      // IM
      CoarseAverage cellAverageIM(IMFine.disjointBoxLayout(),
                                  IMCoar.disjointBoxLayout(),
                                  1, a_coarseningFactor);
      for (DataIterator dit = IMCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        IMCoar[dit()].setVal(0.);
      cellAverageIM.averageToCoarse(IMCoar, IMFine);
  }

  // Handle inter-box ghost cells.
  BCoar.exchange();
  PiCoar.exchange();
  zbCoar.exchange();
  IMCoar.exchange();

  // Notify any observers of this change.
  notifyObserversOfChange();
}

//-----------------------------------------------------------------------

// ---------------------------------------------------------
int AMRNonLinearPoissonOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_domains.size(); ilev++) {
      if (m_domains[ilev].domainBox() == a_domain.domainBox()) {
          retval = m_refRatios[ilev];
          found = true;
      }
  }

  if (!found)
    {
      MayDay::Abort("Domain not found in AMR hierarchy");
    }

  return retval;
}

#include "NamespaceFooter.H"
