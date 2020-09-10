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
#include "AmrHydro.H"

#include "NamespaceHeader.H"

std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
bool              GlobalBCRS::s_areBCsParsed= false;

void
zeroBCValue(Real* pos, 
            int* dir, 
            Side::LoHiSide* side, 
            Real* a_values)
{
    a_values[0] = 0.0;
}

void ParseNeumannValue(Real* pos,
            int* dir, 
            Side::LoHiSide* side, 
            Real* a_values)
{
    ParmParse pp;
    Real bcVal;
    pp.get("Neumann_bc_value",bcVal);
    a_values[0]=bcVal;
}

void ParseDirichletValue(Real* pos,
            int* dir, 
            Side::LoHiSide* side, 
            Real* a_values)
{
    ParmParse pp;
    Real bcVal;
    pp.get("Dirichlet_bc_value",bcVal);
    a_values[0]=bcVal;
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
          ParmParse ppBC("bc");
          ppBC.getarr("lo_bc", GlobalBCRS::s_bcLo, 0, SpaceDim);
          ppBC.getarr("hi_bc", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
      }

      Box valid = a_valid;
      for(int dir=0; dir<CH_SPACEDIM; ++dir) {
          // don't do anything if periodic -- should be perio in y dir 1
          if (!a_domain.isPeriodic(dir)) {
              Box ghostBoxLo = adjCellBox(valid, dir, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, dir, Side::Hi, 1);
              // box of ghost cells is outside of domain bounds ?
              if(!a_domain.domainBox().contains(ghostBoxLo)) {
                  if (GlobalBCRS::s_bcLo[i] == 0) {
                      pout() << "const diri bcs lo for direction " << i << endl;
		              DiriBC(a_state,
		                 valid,
		                 a_dx,
		                 a_homogeneous,
		                 ParseDirichletValue,
		                 dir,
		                 Side::Lo);
                  } else if (GlobalBCRS::s_bcLo[i] == 1) {
                      pout() << "const neum bcs lo for direction " << i << endl;
		              NeumBC(a_state,
		                 valid,
		                 a_dx,
		                 a_homogeneous,
		                 ParseNeumannValue,
		                 dir,
		                 Side::Lo);
                  }
              }
              // box of ghost cells is outside of domain bounds ?
              if(!a_domain.domainBox().contains(ghostBoxHi)) {
                  if (GlobalBCRS::s_bcHi[i] == 0) {
                      pout() << "const diri bcs hi for direction " << i << endl;
		              DiriBC(a_state,
		                 valid,
		                 a_dx,
		                 a_homogeneous,
		                 ParseDirichletValue,
		                 dir,
		                 Side::Hi);
                  } else if (GlobalBCRS::s_bcLo[i] == 1) {
                      pout() << "const neum bcs lo for direction " << i << endl;
		              NeumBC(a_state,
		                 valid,
		                 a_dx,
		                 a_homogeneous,
		                 ParseNeumannValue,
		                 dir,
                         Side::Hi);
                  }
              }
          } // end if is not periodic in ith direction
      } // end dir loop
  }
}


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
    pout() << "HydroIBC::define" << endl;
    PhysIBC::define(a_domain, a_dx);
}

/// Factory method - this object is its own factory
/**
   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
   its define() must be called before it is used).
*/
HydroIBC*
HydroIBC::new_headIBC()
{
    HydroIBC* retval = new HydroIBC();

    return retval;
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
                           LevelData<FArrayBox>& a_iceHeight)
{
    
    pout() << "HydroIBC::initializeData" << endl;

    IntVect regionLo = m_domain.domainBox().bigEnd();
    IntVect regionHi = m_domain.domainBox().bigEnd();

    DataIterator dit = a_head.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox& thisHead      = a_head[dit];
        FArrayBox& thisGapHeight = a_gapHeight[dit];
        FArrayBox& thisPw        = a_Pw[dit];
        FArrayBox& thisqw        = a_qw[dit];
        FArrayBox& thisRe        = a_Re[dit];
        FArrayBox& thismeltRate  = a_meltRate[dit];
        FArrayBox& thiszbed      = a_zbed[dit];
        FArrayBox& thispi        = a_Pi[dit];
        FArrayBox& thisiceHeight = a_iceHeight[dit];

        BoxIterator bit(thisHead.box()); // Default .box() have ghostcells ?
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            // gapHeight - initially constant 
            thisGapHeight(iv, 0) = Params.m_gapInit;
            // Initial Head is such that the water pressure is equal to 50% of the ice overburden pressure
            thisiceHeight(iv, 0) = Params.m_H;
            Real P_ice        = Params.m_rho_i * Params.m_gravity * Params.m_H;
            Real P_water_init = P_ice * 0.5;
            Real Fact         = 1./(Params.m_rho_w * Params.m_gravity);
            thiszbed(iv, 0)      = Params.m_slope*(iv[0]+0.5)*a_dx[0];
            thisHead(iv, 0)      = P_water_init * Fact + thiszbed(iv, 0) ;
            thispi(iv, 0)        = P_ice; // cst for now ...
            thisPw(iv, 0)        = Params.m_rho_w * Params.m_gravity * (thisHead(iv, 0)  - thiszbed(iv, 0));
            // dummy stuff 
            thisRe(iv, 0)        = Params.m_ReInit;
            thisqw(iv, 0)        = 0.0;
            thisqw(iv, 1)        = 0.0;
            thismeltRate(iv, 0)  = Params.m_G / Params.m_L;
        } // end loop over cells
    }     // end loop over boxes

    m_poissonOpF_head = new VCAMRPoissonOp2Factory;

    pout() << "(Done with HydroIBC::initializeData)" << endl;
}


AMRLevelOpFactory<LevelData<FArrayBox> >* 
HydroIBC::defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&               a_grids,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                      ProblemDomain& coarsestDomain,
                      Vector<int>& refRatio,
                      Real& coarsestDx)
{

    m_poissonOpF_head->define(coarsestDomain,
                      a_grids,
                      refRatio,
                      coarsestDx,
                      &mixBCValues,
                      0.0,
                      a_aCoef,
                      -1.0,
                      a_bCoef);

    return (AMRLevelOpFactory<LevelData<FArrayBox> >*) m_poissonOpF_head;
}

#include "NamespaceFooter.H"
