#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef AMR_DRIVER_H
#define AMR_DRIVER_H

#include "LevelData.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "PatchGodunov.H"
#include "BasicIBC.H"
#include "AdvectPhysics.H"

#include "NamespaceHeader.H"

/// framework class to manage non-subcycled AMR model
/**
 */
class AmrDriver
{
public:
    /// Default constructor
    /** At the moment, there is only one constructor, which defines itself
        based on an inputs file.  (This may change, of course)
    */
    AmrDriver();

    /// destructor
    /** destructor
     */
    virtual ~AmrDriver();

    /// set default values before initialization
    void setDefaults();

    /// initializes object based on inputs data passed in through parmParse
    /** initializes new object based on data in ParmParse database.  This
        may change to a more explicit initialization if it's thought to be
        a good idea.
    */
    void initialize();

    /// set domain size (in meters)
    void setDomainSize(const RealVect& a_domainSize)
    {
        m_domainSize = a_domainSize;
    }

    /// set IBC
    void setIBC(BasicIBC* a_IBC);

    /// advance solution until either max_time or max_step are reached
    /**
     */
    void run(Real a_max_time, int a_max_step);

    /// compute one timestep
    void timeStep(Real a_dt);

    /// compute fluxes
    void computeFluxes(Vector<LevelData<FluxBox>*>& a_vectFluxes, Real a_dt);

    void updatePhi(Vector<LevelData<FArrayBox>*>& a_new_phi,
                   const Vector<LevelData<FArrayBox>*>& a_old_phi,
                   const Vector<LevelData<FluxBox>*>& a_fluxes,
                   Real a_dt);

    /// do regridding
    void regrid();

    /// is this object defined and initialized?
    /**
     */
    bool isDefined() const;

#ifdef CH_USE_HDF5

    /// write hdf5 plotfile to the standard location
    /** writes plotfile to the standard location; filename is
        <plot_root>.<step>.<DIM>d.hdf5
     */
    void writePlotFile();

    /// write checkpoint file out for later restarting
    /** writes checkpoint data to the standard location; filename is
        <check_root>.<step>.<DIM>d.hdf5
     */
    void writeCheckpointFile() const;

    /// read checkpoint file for restart
    /** read checkpoint data from file pointed to by a_handle
     */
    void readCheckpointFile(HDF5Handle& a_handle);

    /// set up for restart
    void restart(string& a_restart_file);

#endif

    /// set current solution time
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

    const LevelData<FArrayBox>* phi(int a_level) const
    {
        return m_phi[a_level];
    };

    Vector<LevelData<FArrayBox>*>& amrPhi()
    {
        return m_phi;
    }

    /// set physical parameters (useful when calling from another driver)
    void setPhysicalConstants(Real a_gravity);

    // accessor for AMR grids
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

#if 0
  /// set up BC's based on parmparse inputs
  void setupBCs();
#endif

    /// IBC
    BasicIBC* m_IBCPtr;

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

    /// max box size for level 0 (defaults to max_box_size)
    /** this allows a different max box size for the level 0 domain,
        which can be useful since the number of boxes on level 0 is often
        driven by load balancing concerns, while max_box_size on other levels
        is generally driven by issues of cache size, etc) */
    int m_max_base_grid_size;

    /// regrid interval
    int m_regrid_interval;
    int m_n_regrids;

    /// stores tags from previous regrids
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

    // current step
    int m_cur_step;

    /// current old-time data
    Vector<LevelData<FArrayBox>*> m_old_phi;

    /// new-time data
    Vector<LevelData<FArrayBox>*> m_phi;

    /// how verbose should we be?
    int m_verbosity;

    /// is this object initialized?
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

    // physical constant -- gravity
    Real m_gravity;

    /// num ghost cells
    int m_num_phi_ghost;

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
