#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef AMR_HYDRO_H
#define AMR_HYDRO_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "HydroIBC.H"
#include "suhmo_params.H"
#include "MultilevelLinearOp.H"
#include "BiCGStabSolver.H"
#include "CellToEdge.H"

#include "NamespaceHeader.H"

class GlobalBCRS  
{
public: 
    static std::vector<int> s_bcLo, s_bcHi; 
    static bool             s_areBCsParsed;
    static Real             s_xlo_diri;
    static Real             s_xhi_diri;
    static Real             s_ylo_diri;
    static Real             s_yhi_diri;
    static Real             s_xlo_neum;
    static Real             s_xhi_neum;
    static Real             s_ylo_neum;
    static Real             s_yhi_neum;
};

/// framework class to manage non-subcycled AMR model
/**
 */
class AmrHydro
{
public:
    /** At the moment, there is only one constructor, which defines itself
        based on an inputs file.
    */
    AmrHydro();

    /* destructor */
    virtual ~AmrHydro();

    /* set default values before initialization */
    void setDefaults();

    /* set default values for hydro_params */
    void setParams();

    /* initializes new object based on data in ParmParse database */
    void initialize();

    /* set domain size (in meters) */
    void setDomainSize(const RealVect& a_domainSize)
    {
        m_domainSize = a_domainSize;
    }

    /* set Init and Boundary Conditions */
    void setIBC(HydroIBC* a_IBC);

    /* advance the solution until either max_time or max_step are reached */
    void run(Real a_max_time, int a_max_step);

    /* compute one timestep */
    void timeStep(Real a_dt);

    /* compute fluxes */
    void computeFluxes(Vector<LevelData<FluxBox>*>& a_vectFluxes, Real a_dt);

    void updateHead(Vector<LevelData<FArrayBox>*>& a_new_head,
                   const Vector<LevelData<FArrayBox>*>& a_old_head,
                   const Vector<LevelData<FluxBox>*>& a_fluxes,
                   Real a_dt);

    void aCoeff_bCoeff_CC(LevelData<FArrayBox>&  levelacoef, 
                           LevelData<FArrayBox>&  levelbcoef_cc, 
                           LevelData<FArrayBox>&  levelRe, 
                           LevelData<FArrayBox>&  levelB);

    void CalcRHS_head(LevelData<FArrayBox>& levelRHS_h, 
                       LevelData<FArrayBox>& levelPi, 
                       LevelData<FArrayBox>& levelPw, 
                       LevelData<FArrayBox>& levelmR, 
                       LevelData<FArrayBox>& levelB);

    void CalcRHS_gapHeight(LevelData<FArrayBox>& levelRHS_b, 
                       LevelData<FArrayBox>& levelPi, 
                       LevelData<FArrayBox>& levelPw, 
                       LevelData<FArrayBox>& levelmR, 
                       LevelData<FArrayBox>& levelB);

    void SolveForHead(
                      const Vector<DisjointBoxLayout>&               a_grids,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef,
                      ProblemDomain& coarsestDomain,
                      Vector<int>& refRatio,
                      Real& coarsestDx,
                      Vector<LevelData<FArrayBox>*>& a_head, 
                      Vector<LevelData<FArrayBox>*>& a_RHS);

    /* regridding */
    void regrid();


#ifdef CH_USE_HDF5

    void writePltCustom(int numPlotComps,
                        Vector<std::string>& vectName,
                        Vector<Vector<LevelData<FArrayBox>*>> stp,
                        string namePlot);

    /** writes plotfile to the standard location; filename is
        <plot_root>.<step>.<DIM>d.hdf5
     */
    void writePlotFile();

    /** writes checkpoint data to the standard location; filename is
        <check_root>.<step>.<DIM>d.hdf5
     */
    void writeCheckpointFile() const;

    /** read checkpoint data from file pointed to by a_handle
     */
    void readCheckpointFile(HDF5Handle& a_handle);

    /* set up for restart */
    void restart(string& a_restart_file);

#endif

    /* set current solution time */
    void setTime(Real a_time)
    {
        m_time = a_time;
    }

    Real time() const
    {
        return m_time + m_offsetTime;
    };

    Real dt() const
    {
        return m_dt;
    }

    int finestLevel() const
    {
        return m_finest_level;
    }

    /* is this object defined and initialized? */
    bool isDefined() const
    {
        return m_is_defined;
    }

    const RealVect dx(int a_level) const
    {
        return m_amrDx[a_level];
    };

    const DisjointBoxLayout& grids(int a_level) const
    {
        return m_amrGrids[a_level];
    }

    const Vector<int>& refRatios() const
    {
        return m_refinement_ratios;
    }

    const Vector<RealVect>& amrDx() const
    {
        return m_amrDx;
    }
    
    const Vector<ProblemDomain>& amrDomains() const
    {
        return m_amrDomains;
    }

    /* Access functions */
    // m_gapheight
    const LevelData<FArrayBox>* gapHeight(int a_level) const
    {
        return m_gapheight[a_level];
    }
    Vector<LevelData<FArrayBox>*>& amrGapHeight()
    {
        return m_gapheight;
    }

    // m_head
    const LevelData<FArrayBox>* head(int a_level) const
    {
        return m_head[a_level];
    }
    Vector<LevelData<FArrayBox>*>& amrHead()
    {
        return m_head;
    }

    // m_bedelevation
    const LevelData<FArrayBox>* zBed(int a_level) const
    {
        return m_bedelevation[a_level];
    }
    Vector<LevelData<FArrayBox>*>& amrZBed()
    {
        return m_bedelevation;
    }

    /* set physical parameters (useful when calling from another driver) */
    //void setPhysicalConstants(Real a_gravity);

    /* accessor for AMR grids */
    const Vector<DisjointBoxLayout>& amrGrids() const
    {
        return m_amrGrids;
    }

protected:
    /// compute tags for regridding
    void tagCells(Vector<IntVectSet>& a_tags);

    /// compute tags for the level a_level
    void tagCellsLevel(IntVectSet& a_tags, int a_level);

    /// compute tags at initial time
    void tagCellsInit(Vector<IntVectSet>& a_tags);

    /// initialize grids at initial time
    void initGrids(int a_finest_level);

    /// set up grids from grids file
    void setupFixedGrids(const std::string& a_gridFile);

    /// set up storage, etc on an AMR level
    void levelSetup(int a_level, const DisjointBoxLayout& a_grids);

    /// initialize data on hierarchy
    void initData(Vector<LevelData<FArrayBox>*>& a_initialVelocity);

    /// define elliptic solver for computing velocity field
    virtual void defineSolver();

    /// compute timestep
    Real computeDt();

    /// compute timestep at initial time
    Real computeInitialDt();

    /// IBC
    HydroIBC* m_IBCPtr;

    /// max number of levels
    int m_max_level;

    /// current finest level
    int m_finest_level;

    /// regridding lbase -- generally (and by default) this will be 0
    int m_regrid_lbase;

    /// maximum level at which refinement tags are computed
    int m_tag_cap;

    /// blocking factor
    int m_block_factor;

    /// only allow tagging within a subset of the coarse domain
    IntVectSet m_tag_subset;

    Vector<IntVectSet> m_vectTagSubset;

    /// grid efficiency
    Real m_fill_ratio;

    /// proper nesting radius
    int m_nesting_radius;

    /// max box size
    int m_max_box_size;

    /** this allows a different max box size for the level 0 domain,
        which can be useful since the number of boxes on level 0 is often
        driven by load balancing concerns, while max_box_size on other levels
        is generally driven by issues of cache size, etc) */
    int m_max_base_grid_size;

    /* regrid interval */
    int m_regrid_interval;
    int m_n_regrids;

    /* stores tags from previous regrids */
    Vector<IntVectSet> m_vectTags;

    /// tag on grad(phi)?
    bool m_tagOnGradPhi;

    /// tagging value (undivided gradient(phi) threshold for regridding)
    Real m_tagging_val;

    /// tag on Laplacian(phi)
    bool m_tagOnLapPhi;

    /// tagging value for Laplacian(phi)
    Real m_laplacian_tagging_val;

    /// amount to buffer tags used in regridding
    // grow tags in all directions by m_tags_grow
    int m_tags_grow;
    // grow tags in all direction dir by max(0,m_tags_grow_dir[dir] -  m_tags_grow)
    IntVect m_tags_grow_dir;

    /// refinement ratios
    Vector<int> m_refinement_ratios;

    /// cell spacing at each level
    Vector<RealVect> m_amrDx;

    /// problem domains at each level
    Vector<ProblemDomain> m_amrDomains;

    // problem domain size
    RealVect m_domainSize;

    /// current grids
    Vector<DisjointBoxLayout> m_amrGrids;

    /// keeps track of which levels are completely covered
    Vector<int> m_covered_level;

    /// book-keeping; keeps track of number of cells per level
    Vector<int> m_num_cells;

    /// current time
    Real m_time;
    Real m_offsetTime;

    /// most recent timestep
    Real m_dt;

    /// most recent stable timestep (> m_dt if we took a step just to hit a plot interval)
    Real m_stable_dt;

    /// timestep scaling
    Real m_cfl;

    /// cfl number for initial timestep (useful if initial data needs small cfl)
    Real m_initial_cfl;

    /// maximum amount cfl number may grow in one timestep
    Real m_max_dt_grow;

    // debug 
    bool m_PrintCustom;

    /* current step */
    int m_cur_step;

    /** Model outputs include spatio temporal fields of Hydraulic head,  effective  Pressure,  
        subglacial  gap  height  (the  effective  geometry  representative  of  an  entire  element),  
        depth-integrated  water  flux,  and degree of channelization */
    /* current old-time data */
    /* */
    Vector<LevelData<FArrayBox>*> m_old_head;               // Hydraulic head
    Vector<LevelData<FArrayBox>*> m_old_gapheight;          // Subglacial height
    /* */
    Vector<LevelData<FArrayBox>*> m_old_Pw;                 // Water pressure
    Vector<LevelData<FArrayBox>*> m_old_qw;                 // Water flux
    Vector<LevelData<FArrayBox>*> m_old_Re;                 // Reynolds number
    Vector<LevelData<FArrayBox>*> m_old_meltRate;           // Melting rate

    /* new-time data */
    /* */
    Vector<LevelData<FArrayBox>*> m_head;                   // Hydraulic head
    Vector<LevelData<FArrayBox>*> m_gapheight;              // Subglacial height

    Vector<LevelData<FArrayBox>*> m_gradhead;               // Grad of Hydraulic head

    /* */
    Vector<LevelData<FArrayBox>*> m_Pw;                     // Water pressure
    Vector<LevelData<FArrayBox>*> m_qw;                     // Water flux
    Vector<LevelData<FArrayBox>*> m_Re;                     // Reynolds number
    Vector<LevelData<FArrayBox>*> m_meltRate;               // Melting rate
    Vector<LevelData<FArrayBox>*> m_iceheight;              // Ice slab height

    /* Time-constant parameters */
    Vector<LevelData<FArrayBox>*> m_bedelevation;           //
    Vector<LevelData<FArrayBox>*> m_overburdenpress;        // rhoi*g*H


    /* Solver */
    //VCAMRPoissonOp2Factory* m_poisson_h;

    /* verbosity level */
    int m_verbosity;

    /* is this object initialized? */
    bool m_is_defined;

    ///
    bool m_do_restart;

    /// if starting from a restart, timestep of restart
    int m_restart_step;

    ///
    string m_plot_prefix;

    ///
    string m_check_prefix;

    bool m_check_overwrite; // overwrite check point files
    bool m_check_exit;      // exit immediately after writing checkpoint?

    ///
    int m_plot_interval;
    /// write a plot file every m_plot_time_interval years (rather than steps)
    Real m_plot_time_interval;

    /// if true, include grad(phi) in plotfiles
    bool m_write_gradPhi;

    ///
    int m_check_interval;

    // used to store a fixed time step > 0 , if desired
    Real m_fixed_dt;

    /* num ghost cells */
    int m_num_head_ghost;

    /* SUHMO Params */
    suhmo_params* m_suhmoParm;

private:
    template <class T>
    void levelAllocate(LevelData<T>** a_ptrPtr, const DisjointBoxLayout& a_grids, int a_nComp, const IntVect& a_ghost)
    {
        CH_assert(a_ptrPtr != NULL);
        if (*a_ptrPtr != NULL) delete *a_ptrPtr;
        *a_ptrPtr = new LevelData<T>(a_grids, a_nComp, a_ghost);
    }
};

#include "NamespaceFooter.H"

#endif
