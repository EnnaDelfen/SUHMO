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
#include "CoarseAverageFace.H"
#include "AMRMultiGrid.H"
#include "Misc.H"

#include "AMRNonLinearPoissonOpF_F.H"

#include "VCAMRNonLinearPoissonOp.H"
#include "VCAMRNonLinearPoissonOpF_F.H"
#include "DebugOut.H"

#include "ExtrapGhostCells.H"

#include "NamespaceHeader.H"

// Func to call after before the start of each VCycle on finer level. 
// to follow with averaging on subsequent levels
void VCAMRNonLinearPoissonOp::UpdateOperator(const LevelData<FArrayBox>&   a_phi,
                                             const LevelData<FArrayBox>*   a_phicoarsePtr, 
                                             int                           a_depth,
                                             int                           a_AMRFASMGiter,
                                             bool                          a_homogeneous) 
{
  CH_TIME("VCAMRNonLinearPoissonOp::UpdateOperator");

  if(a_homogeneous) {
      MayDay::Abort("VCAMRNonLinearPoissonOp::UpdateOperator homogeneous");
  }

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  phi.exchange(phi.interval(), m_exchangeCopier);

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  DataIterator dit = phi.dataIterator(); 
  for (dit.begin(); dit.ok(); ++dit) {
      m_bc(phi[dit], dbl[dit()],m_domain, m_dx, a_homogeneous);
  }

   //test for updating the Re
   MEMBER_FUNC_PTR(*m_amrHydro, m_waterFluxlevel)(*m_bCoef, a_phi, a_phicoarsePtr,
                                                  *m_B, *m_Pi, m_dx, 
                                                   false, a_AMRFASMGiter, a_depth);

  // Recompute the relaxation coef after updating m_bCoef
  m_lambdaNeedsResetting = true;
  resetLambda();

}

void VCAMRNonLinearPoissonOp::AverageOperator(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                                              int a_depth)
{
  CH_TIME("VCAMRNonLinearPoissonOp::AverageOperator");

  int coarsening = 1;
  for (int i = 0; i < a_depth; i++) {
      coarsening *= 2;
  }

  const VCAMRNonLinearPoissonOp& op = dynamic_cast<const VCAMRNonLinearPoissonOp&>(a_operator);

  // Perform multigrid coarsening on the operator data.
  LevelData<FluxBox>&   bcoefCoar       = *m_bCoef;
  const LevelData<FluxBox>&   bcoefFine = *(op.m_bCoef);

  if (coarsening != 1) {
      // bCoef
      CoarseAverageFace faceAverage(bcoefFine.disjointBoxLayout(),
                                    1, coarsening);
      faceAverage.averageToCoarse(bcoefCoar, bcoefFine);
  }

  // Handle inter-box ghost cells.
  bcoefCoar.exchange();

  // Recompute the relaxation coef after updating m_bCoef
  m_lambdaNeedsResetting = true;
  resetLambda();
}


void VCAMRNonLinearPoissonOp::residualI(LevelData<FArrayBox>&            a_lhs,
                                        const LevelData<FArrayBox>&      a_phi,
                                        const LevelData<FArrayBox>&      a_rhs,
                                        bool                             a_homogeneous)
{
  CH_TIME("VCAMRNonLinearPoissonOp::residualI");
  if (m_verbosity > 3) {
      pout() << "VCAMRNonLinearPoissonOp::residualI with homgeneous"<< a_homogeneous <<"\n";
  }
  if (a_homogeneous) {
      MayDay::Abort("VCAMRNonLinearPoissonOp::residualI homogeneous");
  }

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();
  {
    CH_TIME("VCAMRNonLinearPoissonOp::residualIBC");

    for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, m_dx, a_homogeneous);
    }
  }

  phi.exchange(phi.interval(), m_exchangeCopier);

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                          *m_B, *m_Pi, *m_zb);

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit()];
      const FluxBox& thisBCoef = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
      FORT_VCNLCOMPUTERES1D
#elif CH_SPACEDIM == 2
      FORT_VCNLCOMPUTERES2D
#elif CH_SPACEDIM == 3
      FORT_VCNLCOMPUTERES3D
#else
      This_will_not_compile!
#endif
                         (CHF_FRA(a_lhs[dit]),
                          CHF_CONST_FRA(phi[dit]),
                          CHF_CONST_FRA(a_rhs[dit]),
                          CHF_CONST_REAL(m_alpha),
                          CHF_CONST_FRA((*m_aCoef)[dit]),
                          CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                          CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                          CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                          CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                          This_will_not_compile!
#endif
                          CHF_CONST_FRA(a_nlfunc[dit]),
                          CHF_BOX(region),
                          CHF_CONST_REAL(m_dx));
    } // end loop over boxes
}

/**************************/
// this preconditioner first initializes phihat to (IA)phihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
// Preconditionior is not used for FAS solve !!
void VCAMRNonLinearPoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                                      const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("VCAMRNonLinearPoissonOp::preCond");

  // diagonal term of this operator in:
  //
  //       alpha * a(i)
  //     + beta  * sum_over_dir (b(i-1/2*e_dir) + b(i+1/2*e_dir)) / (dx*dx)
  //
  // The inverse of this is our initial multiplier.

  int ncomp = a_phi.nComp();

  CH_assert(m_lambda.isDefined());
  CH_assert(a_rhs.nComp()    == ncomp);
  CH_assert(m_bCoef->nComp() == ncomp);

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // also need to average and sum face-centered bCoefs to cell-centers
      Box gridBox = a_rhs[dit].box();

      // approximate inverse
      a_phi[dit].copy(a_rhs[dit]);
      a_phi[dit].mult(m_lambda[dit], gridBox, 0, 0, ncomp);
    }
  int dummyDepth = 0;
  m_print = false;
  relax(a_phi, a_rhs, 2, dummyDepth, dummyDepth);
}

void VCAMRNonLinearPoissonOp::applyOpMg(LevelData<FArrayBox>& a_lhs,
                                        LevelData<FArrayBox>& a_phi,
                                        LevelData<FArrayBox>* a_phiCoarse,
                                        bool a_homogeneous)
{
  CH_TIME("VCAMRNonLinearPoissonOp::applyOpMg");

  // Do CF stuff if we have a coarser level that's not just a single grid cell
   if (a_phiCoarse != NULL) {
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

void VCAMRNonLinearPoissonOp::applyOpI(LevelData<FArrayBox>&       a_lhs,
                                       const LevelData<FArrayBox>& a_phi,
                                       bool                        a_homogeneous )
{
  CH_TIME("VCAMRNonLinearPoissonOp::applyOpI");
  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
  Real dx = m_dx;
  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      m_bc(phi[dit], dbl[dit()],m_domain, dx, a_homogeneous);
    }

  applyOpNoBoundary(a_lhs, a_phi);
}

void VCAMRNonLinearPoissonOp::applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                                const LevelData<FArrayBox>& a_phi)
{
  CH_TIME("VCAMRNonLinearPoissonOp::applyOpNoBoundary");
  if (m_verbosity > 3) {
      pout() <<"VCAMRNonLinearPoissonOp::applyOpNoBoundary \n";
  }

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = phi.dataIterator();

  phi.exchange(phi.interval(), m_exchangeCopier);

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                          *m_B, *m_Pi, *m_zb);

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit()];
      const FluxBox& thisBCoef = (*m_bCoef)[dit];
#if CH_SPACEDIM == 1
      FORT_VCNLCOMPUTEOP1D
#elif CH_SPACEDIM == 2
      FORT_VCNLCOMPUTEOP2D
#elif CH_SPACEDIM == 3
      FORT_VCNLCOMPUTEOP3D
#else
      This_will_not_compile!
#endif
                        (CHF_FRA(a_lhs[dit]),
                         CHF_CONST_FRA(phi[dit]),
                         CHF_CONST_REAL(m_alpha),
                         CHF_CONST_FRA((*m_aCoef)[dit]),
                         CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                         CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                         CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                         CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                         This_will_not_compile!
#endif
                         CHF_CONST_FRA(a_nlfunc[dit]),
                         CHF_BOX(region),
                         CHF_CONST_REAL(m_dx));
    } // end loop over boxes
}

void VCAMRNonLinearPoissonOp::restrictR(LevelData<FArrayBox>& a_phiCoarse,
                                        const LevelData<FArrayBox>& a_phiFine)
{
  if (m_verbosity > 3) {
      pout() << "VCAMRNonLinearPoissonOp::restrictR\n"; 
  }

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
  {
    const FArrayBox&       phiFine = a_phiFine[dit];
    FArrayBox&             phiCoarse = a_phiCoarse[dit];

    Box region = dblFine.get(dit());
    const IntVect& iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    phiCoarse.setVal(0.0);

    FORT_RESTRICTVCNL(CHF_FRA_SHIFT(phiCoarse, civ),
                  CHF_CONST_FRA_SHIFT(phiFine, iv),
                  CHF_BOX_SHIFT(region, iv),
                  CHF_CONST_REAL(m_dx));
  }
}


void VCAMRNonLinearPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                             LevelData<FArrayBox>&       a_phiFine,
                                             const LevelData<FArrayBox>& a_rhsFine) 
{
  // default impl
  restrictResidual(a_resCoarse, a_phiFine, nullptr, a_rhsFine, true);
}


void VCAMRNonLinearPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                               LevelData<FArrayBox>&       a_phiFine,
                                               const LevelData<FArrayBox>* a_phiCoarse,
                                               const LevelData<FArrayBox>& a_rhsFine,
                                               bool homogeneous)
{
  CH_TIME("VCAMRNonLinearPoissonOp::restrictResidual");
  if(homogeneous) {
      MayDay::Abort("VCAMRNonLinearPoissonOp::restrictResidual homogeneous");
  }

  if (a_phiCoarse != nullptr) {
      m_interpWithCoarser.coarseFineInterp(a_phiFine, *a_phiCoarse);
  }

  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit) {
      FArrayBox& phi = a_phiFine[dit];
      m_bc(phi, dblFine[dit()], m_domain, m_dx, homogeneous);
  }

  a_phiFine.exchange(a_phiFine.interval(), m_exchangeCopier);

  LevelData<FArrayBox>  a_nlfunc(dblFine, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dblFine, 1, IntVect::Zero);
  MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phiFine,
                                          *m_B, *m_Pi, *m_zb);

  for (DataIterator dit = a_phiFine.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox&       phi = a_phiFine[dit];
      const FArrayBox& rhs = a_rhsFine[dit];
      FArrayBox&       res = a_resCoarse[dit];

      FArrayBox&       nlfunc = a_nlfunc[dit];

      const FArrayBox& thisACoef = (*m_aCoef)[dit];
      const FluxBox&   thisBCoef = (*m_bCoef)[dit];

      Box region = dblFine.get(dit());
      const IntVect& iv = region.smallEnd();
      IntVect civ = coarsen(iv, 2);

      res.setVal(0.0);

#if CH_SPACEDIM == 1
      FORT_RESTRICTRESVCNL1D
#elif CH_SPACEDIM == 2
      FORT_RESTRICTRESVCNL2D
#elif CH_SPACEDIM == 3
      FORT_RESTRICTRESVCNL3D
#else
      This_will_not_compile!
#endif
                          (CHF_FRA_SHIFT(res, civ),
                           CHF_CONST_FRA_SHIFT(phi, iv),
                           CHF_CONST_FRA_SHIFT(rhs, iv),
                           CHF_CONST_REAL(m_alpha),
                           CHF_CONST_FRA_SHIFT(thisACoef, iv),
                           CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                           CHF_CONST_FRA_SHIFT(thisBCoef[0], iv),
#endif
#if CH_SPACEDIM >= 2
                           CHF_CONST_FRA_SHIFT(thisBCoef[1], iv),
#endif
#if CH_SPACEDIM >= 3
                           CHF_CONST_FRA_SHIFT(thisBCoef[2], iv),
#endif
#if CH_SPACEDIM >= 4
                           This_will_not_compile!
#endif
                           CHF_CONST_FRA_SHIFT(nlfunc, iv),
                           CHF_BOX_SHIFT(region, iv),
                           CHF_CONST_REAL(m_dx));
    }
}

void VCAMRNonLinearPoissonOp::setAlphaAndBeta(const Real& a_alpha,
                                              const Real& a_beta)
{
  m_alpha  = a_alpha;
  m_beta   = a_beta;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void VCAMRNonLinearPoissonOp::computeCoeffsOTF(bool a_update_operator)
{
  m_update_operator = a_update_operator;
}


void VCAMRNonLinearPoissonOp::setSpecificParams(const RefCountedPtr<LevelData<FArrayBox> >& a_B,
                                                const RefCountedPtr<LevelData<FArrayBox> >& a_Pi,
                                                const RefCountedPtr<LevelData<FArrayBox> >& a_zb)
{
  // Problem SPECIFIC
  m_B  = a_B;  // Gap Height
  m_Pi = a_Pi; // Overb Press
  m_zb = a_zb; // Bed Elevation
}

void VCAMRNonLinearPoissonOp::setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                                       const RefCountedPtr<LevelData<FluxBox  > >& a_bCoef,
                                       const Real&                                 a_alpha,
                                       const Real&                                 a_beta)
{
  m_alpha = a_alpha;
  m_beta  = a_beta;

  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void VCAMRNonLinearPoissonOp::resetLambda()
{
  if (m_lambdaNeedsResetting) {
      Real scale = 1.0 / (m_dx*m_dx);

      // Compute it box by box, point by point
      for (DataIterator dit = m_lambda.dataIterator(); dit.ok(); ++dit) {
          FArrayBox&       lambdaFab = m_lambda[dit];
          const FArrayBox& aCoefFab  = (*m_aCoef)[dit];
          const FluxBox&   bCoefFab  = (*m_bCoef)[dit];
          const Box& curBox          = lambdaFab.box();

          // Compute the diagonal term
          lambdaFab.copy(aCoefFab);
          lambdaFab.mult(m_alpha);

          for (int dir = 0; dir < SpaceDim; dir++) {
              FORT_SUMFACESNL(CHF_FRA(lambdaFab),
                  CHF_CONST_REAL(m_beta),
                  CHF_CONST_FRA(bCoefFab[dir]),
                  CHF_BOX(curBox),
                  CHF_CONST_INT(dir),
                  CHF_CONST_REAL(scale));
          }
      }

      // Lambda is reset.
      m_lambdaNeedsResetting = false;
  }
}

// Compute the reciprocal of the diagonal entry of the operator matrix
void VCAMRNonLinearPoissonOp::computeLambda()
{
  CH_TIME("VCAMRNonLinearPoissonOp::computeLambda");

  CH_assert(!m_lambda.isDefined());

  // Define lambda
  m_lambda.define(m_aCoef->disjointBoxLayout(),m_aCoef->nComp());

  resetLambda();
}

//
// VCAMRNonLinearPoissonOp::reflux()
//   There are currently the new version (first) and the old version (second)
//   in this file.  Brian asked to preserve the old version in this way for
//   now. - TJL (12/10/2007)
//
void VCAMRNonLinearPoissonOp::reflux(const LevelData<FArrayBox>&        a_phiFine,
                                     const LevelData<FArrayBox>&        a_phi,
                                     LevelData<FArrayBox>&              a_residual,
                                     AMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
  CH_TIMERS("VCAMRNonLinearPoissonOp::reflux");
  CH_TIMER("VCAMRNonLinearPoissonOp::reflux::incrementCoarse", t2);
  CH_TIMER("VCAMRNonLinearPoissonOp::reflux::incrementFine", t3);

  m_levfluxreg.setToZero();
  Interval interv(0,a_phi.nComp()-1);

  CH_START(t2);

  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit) {
      const FArrayBox& coarfab   = a_phi[dit];
      const FluxBox&   coarBCoef = (*m_bCoef)[dit];
      const Box&       gridBox   = a_phi.getBoxes()[dit];

      if (m_levfluxreg.hasCF(dit())) {
          for (int idir = 0; idir < SpaceDim; idir++) {
                FArrayBox coarflux;
                Box faceBox = surroundingNodes(gridBox, idir);

                getFlux(coarflux, coarfab, coarBCoef, faceBox, idir);

                Real scale = 1.0;
                m_levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                             interv, interv, idir);
          }
      }
  }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = ( LevelData<FArrayBox>&)a_phiFine;

  VCAMRNonLinearPoissonOp* finerAMRPOp = (VCAMRNonLinearPoissonOp*) a_finerOp;
  QuadCFInterp& quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(phiFineRef, a_phi);
  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // phiFineRef.exchange(a_phiFine.interval());
  IntVect phiGhost = phiFineRef.ghostVect();
  int ncomps = a_phiFine.nComp();

  CH_START(t3);

  DataIterator ditf = a_phiFine.dataIterator();
  const DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf) {
      const FArrayBox& phifFab   = a_phiFine[ditf];
      const FluxBox&   fineBCoef = (*(finerAMRPOp->m_bCoef))[ditf];
      const Box&       gridbox   = dblFine.get(ditf());

      for (int idir = 0; idir < SpaceDim; idir++) {
          //int normalGhost = phiGhost[idir];
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next()) {
              if (m_levfluxreg.hasCF(ditf(), sit())) {
                  Side::LoHiSide hiorlo = sit();
                  Box fluxBox = bdryBox(gridbox,idir,hiorlo,1);

                  FArrayBox fineflux(fluxBox,ncomps);
                  getFlux(fineflux, phifFab, fineBCoef, fluxBox, idir,
                          m_refToFiner);

                  Real scale = 1.0;
                  m_levfluxreg.incrementFine(fineflux, scale, ditf(),
                                             interv, interv, idir, hiorlo);
              }
          }
      }
  }

  CH_STOP(t3);

  Real scale = 1.0/m_dx;
  m_levfluxreg.reflux(a_residual, scale);
}

void VCAMRNonLinearPoissonOp::levelGSRB(LevelData<FArrayBox>&       a_phi,
                                        const LevelData<FArrayBox>& a_rhs,
                                        int                         a_ite,
                                        int                         a_AMRFASMGiter, 
                                        int                         a_depth)
{
  CH_TIME("VCAMRNonLinearPoissonOp::levelGSRB");
    if (m_print) {
        pout() << "VCAMRNonLinearPoissonOp::levelGSRB, smooth, depth, AMRFASMGit = " << a_ite << " " << a_depth << " " << a_AMRFASMGiter << endl;
    }

  CH_assert(a_phi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_phi.ghostVect() >= IntVect::Unit);
  CH_assert(a_phi.nComp() == a_rhs.nComp());

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();

  LevelData<FArrayBox>  a_nlfunc(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox>  a_nlDfunc(dbl, 1, IntVect::Zero);

  DataIterator dit = a_phi.dataIterator();
  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++)
    {

      // fill in intersection of ghostcells and a_phi's boxes
      {
        CH_TIME("VCAMRNonLinearPoissonOp::levelGSRB::homogeneousCFInterp");
        // For FAS, we don't want to do this?
        //        homogeneousCFInterp(a_phi);
      }

      {
        CH_TIME("VCAMRNonLinearPoissonOp::levelGSRB::exchange");
        a_phi.exchange(a_phi.interval(), m_exchangeCopier);
      }

      {
        CH_TIME("VCAMRNonLinearPoissonOp::levelGSRB::BCs");
        // now step through grids...
        for (dit.begin(); dit.ok(); ++dit)
          {
            // invoke physical BC's where necessary
            m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, false);
          }
      }

      MEMBER_FUNC_PTR(*m_amrHydro, m_nllevel)(a_nlfunc, a_nlDfunc, a_phi,
                                              *m_B, *m_Pi, *m_zb);

      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& region = dbl.get(dit());
          const FluxBox& thisBCoef  = (*m_bCoef)[dit];

#if CH_SPACEDIM == 1
          FORT_GSRBHELMHOLTZVCNL1D
#elif CH_SPACEDIM == 2
          FORT_GSRBHELMHOLTZVCNL2D
#elif CH_SPACEDIM == 3
          FORT_GSRBHELMHOLTZVCNL3D
#else
          This_will_not_compile!
#endif
                                (CHF_FRA(a_phi[dit]),
                                 CHF_CONST_FRA(a_rhs[dit]),
                                 CHF_BOX(region),
                                 CHF_CONST_REAL(m_dx),
                                 CHF_CONST_REAL(m_alpha),
                                 CHF_CONST_FRA((*m_aCoef)[dit]),
                                 CHF_CONST_REAL(m_beta),
#if CH_SPACEDIM >= 1
                                 CHF_CONST_FRA(thisBCoef[0]),
#endif
#if CH_SPACEDIM >= 2
                                 CHF_CONST_FRA(thisBCoef[1]),
#endif
#if CH_SPACEDIM >= 3
                                 CHF_CONST_FRA(thisBCoef[2]),
#endif
#if CH_SPACEDIM >= 4
                                 This_will_not_compile!
#endif
                                 CHF_CONST_FRA(a_nlfunc[dit]),
                                 CHF_CONST_FRA(a_nlDfunc[dit]),
                                 CHF_CONST_FRA(m_lambda[dit]),
                                 CHF_CONST_INT(whichPass));
        } // end loop through grids
    } // end loop through red-black


    // Uodate GC for purpose of recomputing bCoeff
    // not sure we still need this
    a_phi.exchange(a_phi.interval(), m_exchangeCopier);
    if (m_print) {
        //pout() << "VCAMRNonLinearPoissonOp::levelGSRB, smooth, depth, AMRFASMGit = " << a_ite << " " << a_depth << " " << a_AMRFASMGiter << endl;
        MEMBER_FUNC_PTR(*m_amrHydro, m_print_data)(a_phi, a_ite, a_depth, a_AMRFASMGiter);
    }

    for (dit.begin(); dit.ok(); ++dit) {
        m_bc(a_phi[dit], dbl[dit()], m_domain, m_dx, true);
    }
}

void VCAMRNonLinearPoissonOp::levelMultiColor(LevelData<FArrayBox>&       a_phi,
                                     const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Abort("VCAMRNonLinearPoissonOp::levelMultiColor - Not implemented");
}

void VCAMRNonLinearPoissonOp::looseGSRB(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Abort("VCAMRNonLinearPoissonOp::looseGSRB - Not implemented");
}

void VCAMRNonLinearPoissonOp::overlapGSRB(LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Abort("VCAMRNonLinearPoissonOp::overlapGSRB - Not implemented");
}

void VCAMRNonLinearPoissonOp::levelGSRBLazy(LevelData<FArrayBox>&       a_phi,
                                   const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Abort("VCAMRNonLinearPoissonOp::levelGSRBLazy - Not implemented");
}

void VCAMRNonLinearPoissonOp::levelJacobi(LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>& a_rhs)
{
  MayDay::Abort("VCAMRNonLinearPoissonOp::levelJacobi - Not implemented");
}

void VCAMRNonLinearPoissonOp::getFlux(FArrayBox&            a_flux,
                                       const FArrayBox&     a_data,
                                       const FluxBox&       a_bCoef,
                                       const Box&           a_facebox,
                                       int                  a_dir,
                                       int                  a_ref) const
{
  CH_TIME("VCAMRNonLinearPoissonOp::getFlux");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  CH_assert(!a_facebox.isEmpty());

  // probably the simplest way to test centering
  // a_box needs to be face-centered in the a_dir
  Box faceTestBox(IntVect::Zero, IntVect::Unit);
  faceTestBox.surroundingNodes(a_dir);
  CH_assert(a_facebox.type() == faceTestBox.type());

  const FArrayBox& bCoefDir = a_bCoef[a_dir];

  // reality check for bCoef
  CH_assert(bCoefDir.box().contains(a_facebox));

  a_flux.resize(a_facebox, a_data.nComp());
  BoxIterator bit(a_facebox);

  Real scale = m_beta * a_ref / m_dx;

  for ( bit.begin(); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      IntVect shiftiv = BASISV(a_dir);
      IntVect ivlo = iv - shiftiv;
      IntVect ivhi = iv;

      CH_assert(a_data.box().contains(ivlo));
      CH_assert(a_data.box().contains(ivhi));

      for (int ivar = 0; ivar < a_data.nComp(); ivar++)
        {
          Real phihi = a_data(ivhi,ivar);
          Real philo = a_data(ivlo,ivar);
          Real gradphi = (phihi - philo ) * scale;

          a_flux(iv,ivar) = -bCoefDir(iv, ivar) * gradphi;
        }
    }
}

//-----------------------------------------------------------------------
void
VCAMRNonLinearPoissonOp::setTime(Real a_time)
{
  // Jot down the time.
  m_time = a_time;

  // Interpolate the b coefficient data if necessary / possible. If
  // the B coefficient depends upon the solution, the operator is nonlinear
  // and the integrator must decide how to treat it.
  if (!m_bCoefInterpolator.isNull() &&
      !m_bCoefInterpolator->dependsUponSolution())
    m_bCoefInterpolator->interpolate(*m_bCoef, a_time);

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;

  // Set the time on the boundary holder.
  m_bc.setTime(a_time);

  // Notify our observers that the time has been set.
  // FIXME: Must implement response of multigrid operators!
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

// Factory
VCAMRNonLinearPoissonOpFactory::VCAMRNonLinearPoissonOpFactory()
{
  setDefaultValues();
}

//-----------------------------------------------------------------------
//  AMR Factory define function
void VCAMRNonLinearPoissonOpFactory::define(const ProblemDomain&         a_coarseDomain,
                                            const Vector<DisjointBoxLayout>&      a_grids,
                                            const Vector<int>&                    a_refRatios,
                                            const Real&                           a_coarsedx,
                                            BCHolder                              a_bc,
                                            const Real&                           a_alpha,
                                            Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                                            const Real&                           a_beta,
                                            Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                                            AmrHydro* a_amrHydro, NL_level a_nllevel,
                                            waterFlux_level a_wFlvl,
                                            printData a_PDfunc,
                                            Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_B,
                                            Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_Pi,
                                            Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_zb,
                                            bool a_update_operator)
{
  CH_TIME("VCAMRNonLinearPoissonOpFactory::define");

  setDefaultValues();

  m_boxes = a_grids;

  m_refRatios = a_refRatios;

  m_bc = a_bc;

  m_dx.resize(a_grids.size());
  m_dx[0] = a_coarsedx;

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

      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);

      m_exchangeCopiers[i].exchangeDefine(a_grids[i], IntVect::Unit);
      m_exchangeCopiers[i].trimEdges(a_grids[i], IntVect::Unit);

      m_cfregion[i].define(a_grids[i], m_domains[i]);
    }

  m_alpha = a_alpha;
  m_aCoef = a_aCoef;

  m_beta  = a_beta;
  m_bCoef = a_bCoef;

  m_amrHydro   = a_amrHydro;
  m_nllevel    = a_nllevel;
  m_waterFluxlevel = a_wFlvl;
  m_print_data     = a_PDfunc;

  m_verbosity = 3;
  m_update_operator = a_update_operator;
  m_print = false;

  m_B  = a_B;  // Gap Height
  m_Pi = a_Pi; // Overb Press  
  m_zb = a_zb; // Bed Elevation

}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// AMR Factory define function, with coefficient data allocated automagically
// for operators.
void
VCAMRNonLinearPoissonOpFactory::define(const ProblemDomain&   a_coarseDomain,
                                       const Vector<DisjointBoxLayout>& a_grids,
                                       const Vector<int>&               a_refRatios,
                                       const Real&                      a_coarsedx,
                                       BCHolder                         a_bc,
                                       const IntVect&                   a_ghostVect)
{
  // This just allocates coefficient data, sets alpha = beta = 1, and calls
  // the other define() method.
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  aCoef(a_grids.size());
  Vector<RefCountedPtr<LevelData<FluxBox> > >    bCoef(a_grids.size());
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  B(a_grids.size());  // Gap Height
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  Pri(a_grids.size());// Overb Press  
  Vector<RefCountedPtr<LevelData<FArrayBox> > >  zb(a_grids.size()); // Bed Elevation
  AmrHydro* amrHydro;
  NL_level  nllevel; 
  waterFlux_level wFlevel;
  printData PDfunc;
  for (int i = 0; i < a_grids.size(); ++i)
  {
    aCoef[i] = RefCountedPtr<LevelData<FArrayBox> >(
                 new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));
    bCoef[i] = RefCountedPtr<LevelData<FluxBox> >(
                 new LevelData<FluxBox>(a_grids[i], 1, a_ghostVect));

    B[i]   = RefCountedPtr<LevelData<FArrayBox> >(
                 new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));
    Pri[i] = RefCountedPtr<LevelData<FArrayBox> >(
                 new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));
    zb[i]  = RefCountedPtr<LevelData<FArrayBox> >(
                 new LevelData<FArrayBox>(a_grids[i], 1, a_ghostVect));

    // Initialize the a and b coefficients to 1 for starters.
    for (DataIterator dit = aCoef[i]->dataIterator(); dit.ok(); ++dit)
    {
      (*aCoef[i])[dit()].setVal(1.0);
      for (int idir = 0; idir < SpaceDim; ++idir)
        (*bCoef[i])[dit()][idir].setVal(1.0);

      (*B[i])[dit()].setVal(1.0);
      (*Pri[i])[dit()].setVal(1.0);
      (*zb[i])[dit()].setVal(1.0);
    }
  }
  // these choices are weird and not in accordance with the default Lapl
  Real alpha = 1.0, beta = 1.0;
  define(a_coarseDomain, a_grids, a_refRatios, a_coarsedx, a_bc,
         alpha, aCoef, beta, bCoef, amrHydro, nllevel, wFlevel, PDfunc,
         B, Pri, zb, true);
}

//-----------------------------------------------------------------------
MGLevelOp<LevelData<FArrayBox> >* VCAMRNonLinearPoissonOpFactory::MGnewOp(const ProblemDomain& a_indexSpace,
                                                                          int                  a_depth,
                                                                          bool                 a_homoOnly)
{
  CH_TIME("VCAMRNonLinearPoissonOpFactory::MGnewOp");
  if (m_verbosity > 3) {
      pout() << "In VCAMRNonLinearPoissonOpFactory::MGnewOp ! depth "<< a_depth <<"\n"; 
  }

  Real dxCrse = -1.0;

  int ref;
  for (ref = 0; ref < m_domains.size(); ref++) {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) {
          break;
      }
  }

  CH_assert(ref !=  m_domains.size()); // didn't find domain

  if (ref > 0) {
      dxCrse = m_dx[ref-1];
  }

  ProblemDomain domain(m_domains[ref]);
  Real dx = m_dx[ref];
  int coarsening = 1;

  for (int i = 0; i < a_depth; i++) {
      coarsening *= 2;
      domain.coarsen(2);
  }

  if (m_verbosity > 3) {
      pout() << "VCAMRNonLinearPoissonOpFactory::MGnewOp ( depth, ref_ratio = " << a_depth << " " << coarsening << ")\n";
  }

  if (coarsening > 1 && !m_boxes[ref].coarsenable(coarsening*VCAMRNonLinearPoissonOp::s_maxCoarse)) {
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

  VCAMRNonLinearPoissonOp* newOp = new VCAMRNonLinearPoissonOp;
  newOp->define(layout, dx, domain, m_bc, ex, cfregion);

  newOp->m_alpha       = m_alpha;
  newOp->m_beta        = m_beta;

  newOp->m_update_operator = m_update_operator;
  newOp->m_verbosity      = m_verbosity;
  newOp->m_print          = m_print;

  newOp->m_amrHydro    = m_amrHydro; 
  newOp->m_nllevel     = m_nllevel;
  newOp->m_waterFluxlevel = m_waterFluxlevel;
  newOp->m_print_data     = m_print_data;

  if (a_depth == 0) {
      // don't need to coarsen anything for this
      newOp->m_aCoef = m_aCoef[ref];
      newOp->m_bCoef = m_bCoef[ref];
      // Problem SPECIFIC
      newOp->m_B  = m_B[ref];  // Gap Height
      newOp->m_Pi = m_Pi[ref]; // Overb Press  
      newOp->m_zb = m_zb[ref]; // Bed Elevation
  } else {
      // need to coarsen coefficients
      RefCountedPtr<LevelData<FArrayBox> > aCoef( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FluxBox> > bCoef( new LevelData<FluxBox> );
      aCoef->define(layout, m_aCoef[ref]->nComp(), m_aCoef[ref]->ghostVect());
      bCoef->define(layout, m_bCoef[ref]->nComp(), m_bCoef[ref]->ghostVect());

      // Problem SPECIFIC
      RefCountedPtr<LevelData<FArrayBox> > B( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FArrayBox> > Pri( new LevelData<FArrayBox> );
      RefCountedPtr<LevelData<FArrayBox> > zb( new LevelData<FArrayBox> );
      B->define(layout,   m_B[ref]->nComp(),  m_B[ref]->ghostVect());
      Pri->define(layout, m_Pi[ref]->nComp(), m_Pi[ref]->ghostVect());
      zb->define(layout,  m_zb[ref]->nComp(), m_zb[ref]->ghostVect());

      // average coefficients to coarser level
      // for now, do this with a CoarseAverage --
      // may want to switch to harmonic averaging at some point
      CoarseAverage averager(m_aCoef[ref]->getBoxes(),
                             layout, aCoef->nComp(), coarsening);
      CoarseAverageFace faceAverager(m_bCoef[ref]->getBoxes(),
                                     bCoef->nComp(), coarsening);

      CoarseAverage averagerB(m_B[ref]->getBoxes(),
                             layout, B->nComp(), coarsening, m_B[ref]->ghostVect());
      CoarseAverage averagerPi(m_Pi[ref]->getBoxes(),
                             layout, Pri->nComp(), coarsening, m_Pi[ref]->ghostVect());
      CoarseAverage averagerZb(m_zb[ref]->getBoxes(),
                             layout, zb->nComp(), coarsening, m_zb[ref]->ghostVect());

      if (m_coefficient_average_type == CoarseAverage::arithmetic)
        {
          averager.averageToCoarse(*aCoef, *(m_aCoef[ref]));
          faceAverager.averageToCoarse(*bCoef, *(m_bCoef[ref]));

          averagerB.averageToCoarse(*B,   *(m_B[ref]));
          averagerPi.averageToCoarse(*Pri, *(m_Pi[ref]));
          averagerZb.averageToCoarse(*zb,  *(m_zb[ref]));

          // freakin boundaries and ghost cells
          // exchange ? Right now no bec we added the last arg to CoarseAverage def and it seems to do the perio fine
          const DisjointBoxLayout& B_dbl = B->getBoxes();
          const ProblemDomain&     B_dom = B_dbl.physDomain();
          DataIterator dit               = B->dataIterator();
          for (dit.begin(); dit.ok(); ++dit) {
              const Box& validBox = B_dbl.get(dit);
              NeumBCForB((*B)[dit], validBox, B_dom, dx);
          }
        }
      else if (m_coefficient_average_type == CoarseAverage::harmonic)
        {
          averager.averageToCoarseHarmonic(*aCoef, *(m_aCoef[ref]));
          faceAverager.averageToCoarseHarmonic(*bCoef, *(m_bCoef[ref]));

          averagerB.averageToCoarseHarmonic(*B,   *(m_B[ref]));
          averagerPi.averageToCoarseHarmonic(*Pri, *(m_Pi[ref]));
          averagerZb.averageToCoarseHarmonic(*zb,  *(m_zb[ref]));
        }
      else
        {
          MayDay::Abort("VCAMRNonLinearPoissonOpFactory::MGNewOp -- bad averagetype");
        }

      newOp->m_aCoef = aCoef;
      newOp->m_bCoef = bCoef;

      // Problem SPECIFIC
      newOp->m_B  = B;   // Gap Height
      newOp->m_Pi = Pri; // Overb Press  
      newOp->m_zb = zb;  // Bed Elevation
  }

  newOp->computeLambda();

  newOp->m_dxCrse = dxCrse;

  return (MGLevelOp<LevelData<FArrayBox> >*)newOp;
}

AMRLevelOp<LevelData<FArrayBox> >* VCAMRNonLinearPoissonOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
  CH_TIME("VCAMRNonLinearPoissonOpFactory::AMRnewOp");
  if (m_verbosity > 3) {
      pout() << "In VCAMRNonLinearPoissonOpFactory::AMRnewOp ! \n"; 
  }

  VCAMRNonLinearPoissonOp* newOp = new VCAMRNonLinearPoissonOp;
  Real dxCrse = -1.0;

  // Need to know how many components we're solving for
  int nComp = 1;

  // Find number of comps from m_bCoef, which should be defined on at least one level
  for (int lev = 0; lev <m_bCoef.size(); lev++) {
    if (m_bCoef[lev] != NULL) {
        nComp = m_bCoef[lev]->nComp();
        break;
    }

    // Should never reach this point
    CH_assert("VCAMRNonLinearPoissonOpFactory::AMRnewOp - m_bCoef not defined");
  }

  int ref;

  for (ref = 0; ref< m_domains.size(); ref++) {
      if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) {
          break;
      }
  }

  if (ref == 0) {
      // coarsest AMR level
      if (m_domains.size() == 1 || !m_boxes[1].isClosed())  {
          // no finer level
          newOp->define(m_boxes[0], m_dx[0],
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0]);
      } else {
          // finer level exists but no coarser
          int dummyRat = 1;  // argument so compiler can find right function
          int refToFiner = m_refRatios[0]; // actual refinement ratio
          newOp->define(m_boxes[0],  m_boxes[1], m_dx[0],
                        dummyRat, refToFiner,
                        a_indexSpace, m_bc,
                        m_exchangeCopiers[0], m_cfregion[0],
                        nComp);
      }
  } else if ((ref ==  m_domains.size()-1) || (!m_boxes[ref+1].isClosed())) {
      dxCrse = m_dx[ref-1];

      // finest AMR level
      newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                  m_refRatios[ref-1],
                  a_indexSpace, m_bc,
                  m_exchangeCopiers[ref], m_cfregion[ref],
                  nComp);
  } else if ( ref == m_domains.size()) {
      MayDay::Abort("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");

  } else {
      dxCrse = m_dx[ref-1];

      // intermediate AMR level, full define
      newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                  m_refRatios[ref-1], m_refRatios[ref],
                  a_indexSpace, m_bc,
                  m_exchangeCopiers[ref], m_cfregion[ref],
                  nComp);
  }

  newOp->m_alpha = m_alpha;
  newOp->m_beta  = m_beta;

  newOp->m_update_operator = m_update_operator;
  newOp->m_verbosity      = m_verbosity;
  newOp->m_print          = m_print;

  newOp->m_aCoef = m_aCoef[ref];
  newOp->m_bCoef = m_bCoef[ref];

  newOp->m_amrHydro       = m_amrHydro;
  newOp->m_nllevel        = m_nllevel;
  newOp->m_waterFluxlevel = m_waterFluxlevel;
  newOp->m_print_data     = m_print_data;

  // Problem SPECIFIC
  newOp->m_B  = m_B[ref];  // Gap Height
  newOp->m_Pi = m_Pi[ref]; // Overb Press  
  newOp->m_zb = m_zb[ref]; // Bed Elevation

  if (newOp->m_aCoef != NULL) {
      newOp->computeLambda();
  }

  newOp->m_dxCrse = dxCrse;

  return (AMRLevelOp<LevelData<FArrayBox> >*)newOp;
}

int VCAMRNonLinearPoissonOpFactory::refToFiner(const ProblemDomain& a_domain) const
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

//-----------------------------------------------------------------------
void VCAMRNonLinearPoissonOpFactory::NeumBCForB(FArrayBox& a_state,
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


void VCAMRNonLinearPoissonOpFactory::setDefaultValues()
{
  // Default to Laplacian operator
  m_alpha = 0.0;
  m_beta = -1.0;

  m_coefficient_average_type = CoarseAverage::arithmetic;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
VCAMRNonLinearPoissonOp::finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                                              int a_coarseningFactor)
{
  const VCAMRNonLinearPoissonOp& op =
  dynamic_cast<const VCAMRNonLinearPoissonOp&>(a_operator);

  // Perform multigrid coarsening on the operator data.
  LevelData<FArrayBox>& acoefCoar = *m_aCoef;
  LevelData<FluxBox>&   bcoefCoar = *m_bCoef;
  // Problem SPECIFIC
  LevelData<FArrayBox>& BCoar     = *m_B;  // Gap Height
  LevelData<FArrayBox>& PiCoar    = *m_Pi; // Overb Press
  LevelData<FArrayBox>& zbCoar    = *m_zb; // Bed Elevation

  const LevelData<FArrayBox>& acoefFine = *(op.m_aCoef);
  const LevelData<FluxBox>&   bcoefFine = *(op.m_bCoef);
  // Problem SPECIFIC
  const LevelData<FArrayBox>& BFine  = *(op.m_B);
  const LevelData<FArrayBox>& PiFine = *(op.m_Pi);
  const LevelData<FArrayBox>& zbFine = *(op.m_zb);

  if (a_coarseningFactor != 1) {
      // aCoef
      CoarseAverage cellAverage(acoefFine.disjointBoxLayout(),
                                acoefCoar.disjointBoxLayout(),
                                1, a_coarseningFactor);
      for (DataIterator dit = acoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        acoefCoar[dit()].setVal(0.);
      cellAverage.averageToCoarse(acoefCoar, acoefFine);
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
      // bCoef
      CoarseAverageFace faceAverage(bcoefFine.disjointBoxLayout(),
                                    1, a_coarseningFactor);
      for (DataIterator dit = bcoefCoar.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
        bcoefCoar[dit()].setVal(0.);
      faceAverage.averageToCoarse(bcoefCoar, bcoefFine);
  }

  // Handle inter-box ghost cells.
  acoefCoar.exchange();
  bcoefCoar.exchange();

  // Problem SPECIFIC
  BCoar.exchange();
  //ExtrapGhostCells( BCoar, BCoar.disjointBoxLayout().physDomain());
  PiCoar.exchange();
  zbCoar.exchange();

  // Mark the relaxation coefficient dirty.
  m_lambdaNeedsResetting = true;

  // Notify any observers of this change.
  notifyObserversOfChange();
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
