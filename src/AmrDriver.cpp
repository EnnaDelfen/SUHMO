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

#include "AmrDriver.H"

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
#include "MayDay.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"

// small parameter defining when times are equal
#define TIME_EPS 1.0e-12
// small parameter defining when a norm is "zero"
#define TINY_NORM 1.0e-8

AmrDriver::AmrDriver() : m_IBCPtr(NULL)
{
    setDefaults();
}

void
AmrDriver::setDefaults()
{
    // set some bogus values as defaults
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

AmrDriver::~AmrDriver()
{
    if (m_verbosity > 4)
    {
        pout() << "AmrDriver::~AmrDriver()" << endl;
    }

    // clean up memory
    for (int lev = 0; lev < m_head.size(); lev++)
    {
        if (m_head[lev] != NULL)
        {
            delete m_head[lev];
            m_head[lev] = NULL;
        }
        if (m_gapheight[lev] != NULL)
        {
            delete m_gapheight[lev];
            m_gapheight[lev] = NULL;
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
    }


    if (m_IBCPtr != NULL)
    {
        delete m_IBCPtr;
        m_IBCPtr = NULL;
    }

    // that should be it!
}

void
AmrDriver::initialize()
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::initialize" << endl;
    }

    // set time to be 0 -- do this now in case it needs to be changed later
    m_time = 0.0;
    m_cur_step = 0;

    // first, read in info from parmParse file
    ParmParse ppAmr("AmrDriver");
    Vector<int> ancells(SpaceDim);
    // allows for domains with lower indices which are not positive
    Vector<int> domLoIndex(SpaceDim, 0);
    // slc : SpaceDim == 2 implies poor-mans multidim mode, in which we still
    // care about the number of vertical layers.
    Vector<Real> domsize(SpaceDim);

    // assumption is that domains are not periodic
    bool is_periodic[SpaceDim];
    for (int dir = 0; dir < SpaceDim; dir++) is_periodic[dir] = false;
    Vector<int> is_periodic_int(SpaceDim, 0);

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

    ppAmr.query("verbosity", m_verbosity);
    ppAmr.query("block_factor", m_block_factor);
    ppAmr.query("regrid_lbase", m_regrid_lbase);
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

    // now set up problem domains
    {
        IntVect loVect = IntVect(D_DECL(domLoIndex[0], domLoIndex[1], domLoIndex[3]));
        IntVect hiVect(
            D_DECL(domLoIndex[0] + ancells[0] - 1, domLoIndex[1] + ancells[1] - 1, domLoIndex[2] + ancells[2] - 1));

        ProblemDomain baseDomain(loVect, hiVect);
        // now set periodicity
        for (int dir = 0; dir < SpaceDim; dir++) baseDomain.setPeriodic(dir, is_periodic[dir]);

        // now set up vector of domains
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
    if (ppAmr.contains("gridsFile"))
    {
        usePredefinedGrids = true;
        ppAmr.get("gridsFile", gridFile);
    }

    // check to see if we're restarting from a checkpoint file
    if (!ppAmr.contains("restart_file"))
    {
        // if we're not restarting
        if (m_verbosity > 3)
        {
            pout() << "AmrDriver::initialize - Initializing data containers" << endl;
        }

        //---------------------------------
        // Set AMR hierarchy vector
        //---------------------------------
        // AF: add vars here
        // Time-dependent variables    
        m_old_head.resize(m_max_level + 1, NULL);
        m_head.resize(m_max_level + 1, NULL);

        m_old_gapheight.resize(m_max_level + 1, NULL);
        m_gapheight.resize(m_max_level + 1, NULL);
        
        m_Pw.resize(m_max_level + 1, NULL);
        m_qw.resize(m_max_level + 1, NULL);
        m_Re.resize(m_max_level + 1, NULL);
        m_meltRate.resize(m_max_level + 1, NULL);
        
        // Time constant variables
        m_bedelevation.resize(m_max_level + 1, NULL);
        m_overburdenpress.resize(m_max_level + 1, NULL);

        //-------------------------------------------------
        // For each level, define a collection of FArrayBox
        //-------------------------------------------------
        for (int lev = 0; lev < m_head.size(); lev++)
        {
            m_old_head[lev] = new LevelData<FArrayBox>;
            m_head[lev] = new LevelData<FArrayBox>;

            m_old_gapheight[lev] = new LevelData<FArrayBox>;
            m_gapheight[lev] = new LevelData<FArrayBox>;
           
            m_Pw[lev] = new LevelData<FArrayBox>;
            m_qw[lev] = new LevelData<FArrayBox>;
            m_Re[lev] = new LevelData<FArrayBox>;
            m_meltRate[lev] = new LevelData<FArrayBox>;

            m_bedelevation[lev] = new LevelData<FArrayBox>;
            m_overburdenpress[lev] = new LevelData<FArrayBox>;
        }

        int finest_level = -1;
        if (usePredefinedGrids)
        {
            setupFixedGrids(gridFile);
        }
        else
        {
            // now create  grids
            initGrids(finest_level);
        }

        // that should be it
    }
    else
    {
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
        if (ppAmr.contains("restart_time"))
        {
            Real restart_time;
            ppAmr.get("restart_time", restart_time);
            m_time = restart_time;
        }

        if (ppAmr.contains("restart_step"))
        {
            int restart_step;
            ppAmr.get("restart_step", restart_step);
            m_cur_step = restart_step;
            m_restart_step = restart_step;
        }
    }

    // set up counter of number of cells
    m_num_cells.resize(m_max_level + 1, 0);
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        if (m_verbosity > 4)
        {
            // write out initial grids
            pout() << "Level " << lev << " grids: " << levelGrids << endl;
        }
        LayoutIterator lit = levelGrids.layoutIterator();
        for (lit.begin(); lit.ok(); ++lit)
        {
            const Box& thisBox = levelGrids.get(lit());
            m_num_cells[lev] += thisBox.numPts();
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

    if (m_verbosity > 3)
    {
        pout() << "Done with AmrDriver::initialize" << endl;
    }
}

// set BC for thickness advection
void
AmrDriver::setIBC(BasicIBC* a_IBC)
{
    m_IBCPtr = a_IBC->new_thicknessIBC();
}

/* Main advance function */
void
AmrDriver::run(Real a_max_time, int a_max_step)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::run -- max_time= " << a_max_time << ", max_step = " << a_max_step << endl;
    }

    Real dt;
    // only call computeInitialDt if we're not doing restart
    if (!m_do_restart)
    {
        dt = computeInitialDt();
    }
    else
    {
        dt = computeDt();
    }

    // advance solution until done
    if (!(m_plot_time_interval > TIME_EPS) || m_plot_time_interval > a_max_time) m_plot_time_interval = a_max_time;

    while (a_max_time > m_time && (m_cur_step < a_max_step))
    {
        Real next_plot_time = m_plot_time_interval * (1.0 + Real(int((m_time / m_plot_time_interval))));
        if (!(next_plot_time > m_time))
        {
            // trap case where machine precision results in (effectively)
            // m_plot_time_interval * (1.0 + Real(int((m_time/m_plot_time_interval)))) == m_time
            next_plot_time += m_plot_time_interval;
        }

        next_plot_time = std::min(next_plot_time, a_max_time);

        while ((next_plot_time > m_time) && (m_cur_step < a_max_step) && (dt > TIME_EPS))
        {
            // dump plotfile before regridding
            if ((m_cur_step % m_plot_interval == 0) && m_plot_interval > 0)
            {
#ifdef CH_USE_HDF5
                writePlotFile();
#endif
            }

            if ((m_cur_step != 0) && (m_cur_step % m_regrid_interval == 0))
            {
                regrid();
            }

            if (m_cur_step != 0)
            {
                // compute dt after regridding in case number of levels has changed
                dt = computeDt();
            }

            if (next_plot_time - m_time + TIME_EPS < dt) dt = std::max(2 * TIME_EPS, next_plot_time - m_time);

            if ((m_cur_step % m_check_interval == 0) && (m_check_interval > 0) && (m_cur_step != m_restart_step))
            {
#ifdef CH_USE_HDF5
                writeCheckpointFile();
#endif
                if (m_cur_step > 0 && m_check_exit)
                {
                    if (m_verbosity > 2)
                    {
                        pout() << "AmrDriver::exit on checkpoint" << endl;
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
        pout() << "AmrDriver::run finished" << endl;
    }
}

/* core of advance routine */
void
AmrDriver::timeStep(Real a_dt)
{
    CH_TIME("AmrDriver::timestep");

    if (m_verbosity >= 2)
    {
        pout() << "Timestep " << m_cur_step << " Advancing solution from time " << m_time << " ( " << time()
               << ")"
                  " with dt = "
               << a_dt << endl;
    }

    // first copy head and gap into old
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        LevelData<FArrayBox>& oldH     = *m_old_head[lev];
        LevelData<FArrayBox>& currentH = *m_head[lev];

        LevelData<FArrayBox>& oldB     = *m_old_gapheight[lev];
        LevelData<FArrayBox>& currentB = *m_gapheight[lev];

        // this way we avoid communication and maintain ghost cells...
        DataIterator dit = oldH.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
        {
            oldH[dit].copy(currentH[dit], 0, 0, 1);
            oldB[dit].copy(currentB[dit], 0, 0, 1);
        }
    }

    // allocate face-centered storage
    Vector<LevelData<FluxBox>*> vectFluxes(m_head.size(), NULL);
    for (int lev = m_finest_level; lev >= 0; lev--)
    {
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];

        IntVect ghostVect = IntVect::Unit;

        // if we're doing AMR, we'll need to average these fluxes
        // down to coarser levels. As things stand now,
        // CoarseAverageFace requires that the coarse LevelData<FluxBox>
        // have a ghost cell.
        vectFluxes[lev] = new LevelData<FluxBox>(m_amrGrids[lev], 1, ghostVect);

        LevelData<FArrayBox>& levelOldHead = *m_old_head[lev];

        // ensure that ghost cells for head are filled in
        if (lev > 0)
        {
            int nGhost = levelOldHead.ghostVect()[0];
            PiecewiseLinearFillPatch headFiller(
                levelGrids, m_amrGrids[lev - 1], 1, m_amrDomains[lev - 1], m_refinement_ratios[lev - 1], nGhost);

            // since we're not subcycling, don't need to interpolate in time
            Real time_interp_coeff = 0.0;
            headFiller.fillInterp(levelOldHead, *m_old_head[lev - 1], *m_old_head[lev - 1], time_interp_coeff, 0, 0, 1);
        }
        // these are probably unnecessary...
        levelOldHead.exchange();
    }

    // compute flux
    computeFluxes(vectFluxes, a_dt);

    // compute div(F) and update solution
    updateHead(m_head, m_old_head, vectFluxes, a_dt);

    // clean up temp storage
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        if (vectFluxes[lev] != NULL)
        {
            delete vectFluxes[lev];
            vectFluxes[lev] = NULL;
        }
    }

    // finally, update to new time and increment current step
    m_dt = a_dt;
    m_time += a_dt;
    m_cur_step += 1;

// write diagnostic info, like sum of ice
#if 0
  Real sumPhi = computeTotalPhi();
  Real diffSum = sumPhi - m_lastSumPhi;
  Real totalDiffSum = sumPhi - m_initialSumPhi;
  
  if (m_verbosity > 0) 
    {
      pout() << "Step " << m_cur_step << ", time = " << m_time << " ( " << time() << " ) " 
             << ": sum(phi) = " << sumIce 
             << " (" << diffSum
             << " " << totalDiffSum
             << ")" << endl;
    }
  
  m_lastSumIce = sumIce;
#endif

    if (m_verbosity > 0)
    {
        pout() << "AmrDriver::timestep " << m_cur_step << " --     end time = "
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
    for (int lev = 0; lev < m_num_cells.size(); lev++)
    {
        totalCellsAdvanced += m_num_cells[lev];
    }

    if (m_verbosity > 0)
    {
        pout() << "Time = " << m_time << " cells advanced = " << totalCellsAdvanced << endl;

        for (int lev = 0; lev < m_num_cells.size(); lev++)
        {
            pout() << "Time = " << m_time << "  level " << lev << " cells advanced = " << m_num_cells[lev] << endl;
        }
    }
}

// compute half-time face-centered thickness using unsplit PPM
void
AmrDriver::computeFluxes(Vector<LevelData<FluxBox>*>& a_Flux, Real a_dt)
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
AmrDriver::updateHead(Vector<LevelData<FArrayBox>*>& a_new_head,
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

#if 0
  //re-fill ghost cells ouside the domain
  for (int lev=0; lev <= finestTimestepLevel()  ; ++lev)
    {
      RealVect levelDx = m_amrDx[lev]*RealVect::Unit;
      m_thicknessIBCPtr->setGeometryBCs(*m_vect_coordSys[lev],
                        q                m_amrDomains[lev],levelDx, m_time, m_dt);

    }
#endif

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
AmrDriver::regrid()
{
    CH_TIME("AmrDriver::regrid");

    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::regrid" << endl;
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
                pout() << "AmrDriver::regrid -- grids unchanged" << endl;
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

/// set physical parameters (useful when calling from another driver)
void
AmrDriver::setPhysicalConstants(Real a_gravity)
{
    m_gravity = a_gravity;
}

void
AmrDriver::tagCells(Vector<IntVectSet>& a_tags)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::tagCells" << endl;
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
    //     pout() << "AmrDriver::tagCells, subset II" << endl;
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
AmrDriver::tagCellsLevel(IntVectSet& a_tags, int a_level)
{
    if (m_verbosity > 4)
    {
        pout() << "AmrDriver::tagCellsLevel " << a_level << endl;
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
AmrDriver::tagCellsInit(Vector<IntVectSet>& a_tags)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::tagCellsInit" << endl;
    }

    // default is to just call regular tagging
    tagCells(a_tags);
    m_vectTags = a_tags;
}

void
AmrDriver::initGrids(int a_finest_level)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::initGrids" << endl;
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
        pout() << "Level 0: " << numCells0 << " cells: " << baseGrids << endl;
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
AmrDriver::setupFixedGrids(const std::string& a_gridFile)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::setupFixedGrids" << endl;
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
AmrDriver::levelSetup(int a_level, const DisjointBoxLayout& a_grids)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::levelSetup - level " << endl;
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
        m_gapheight[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_bedelevation[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
        m_overburdenpress[a_level] = new LevelData<FArrayBox>(a_grids, nPhiComp, HeadGhostVect);
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
AmrDriver::initData(Vector<LevelData<FArrayBox>*>& a_head)

{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::initData" << endl;
    }

    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        LevelData<FArrayBox>& levelHead      = *m_head[lev];
        LevelData<FArrayBox>& levelGapHeight = *m_gapheight[lev];
        LevelData<FArrayBox>& levelPw        = *m_Pw[lev];
        LevelData<FArrayBox>& levelqw        = *m_qw[lev];
        LevelData<FArrayBox>& levelRe        = *m_Re[lev];
        LevelData<FArrayBox>& levelmR        = *m_meltRate[lev];
        LevelData<FArrayBox>& levelzBed      = *m_bedelevation[lev];
        LevelData<FArrayBox>& levelPi        = *m_overburdenpress[lev];
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
                                 levelHead, levelGapHeight, 
                                 levelPw, levelqw,
                                 levelRe, levelmR,
                                 levelzBed, levelPi);

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
AmrDriver::defineSolver()
{
    // just a stub for when we might actually need it...
}

// compute timestep
Real
AmrDriver::computeDt()
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::computeDt" << endl;
    }

    if (m_fixed_dt > TINY_NORM) return m_fixed_dt;

    Real dt = 1.0e50;
    for (int lev = 0; lev <= m_finest_level; lev++)
    {
        Real dtLev = dt;
        const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
        // pretend phi is velocity here
        const LevelData<FArrayBox>& levelVel = *m_head[lev];
        DataIterator levelDit = levelVel.dataIterator();
        for (levelDit.reset(); levelDit.ok(); ++levelDit)
        {
            int p = 0;
            const Box& gridBox = levelGrids[levelDit];
            Real maxVel = 1.0 + levelVel[levelDit].norm(gridBox, p, 0, 1);
            maxVel = max(maxVel,1.0);    
            Real localDt = m_amrDx[lev][0] / maxVel;
            dtLev = min(dtLev, localDt);
        }

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
AmrDriver::computeInitialDt()
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::computeInitialDt" << endl;
    }

    // for now, just call computeDt;
    Real dt = computeDt();
    return dt;
}

#ifdef CH_USE_HDF5

/// write hdf5 plotfile to the standard location
void
AmrDriver::writePlotFile()
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::writePlotFile" << endl;
    }

    // plot comps: head + gapHeight + bedelevation + overburdenPress
    int numPlotComps = 8;

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
    string qwName("Qw"); 
    string ReName("Re"); 
    string meltRateName("meltRate"); 
    string xGradName("xGradPhi");
    string yGradName("yGradPhi");
    string zGradName("zGradPhi");

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
            const Box& gridBox        = m_amrGrids[lev][dit];
            FArrayBox& thisPlotData   = plotDataLev[dit];
            int comp = 0;
            const FArrayBox& thisHead       = levelHead[dit];
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
        if (lev > 0)
        {
            PiecewiseLinearFillPatch interpolator(m_amrGrids[lev],
                                                  m_amrGrids[lev - 1],
                                                  numPlotComps,
                                                  m_amrDomains[lev - 1],
                                                  m_refinement_ratios[lev - 1],
                                                  ghostVect[0]);

            // no interpolation in time
            Real time_interp_coeff = 0.0;
            interpolator.fillInterp(
                *plotData[lev], *plotData[lev - 1], *plotData[lev - 1], time_interp_coeff, 0, 0, numPlotComps);
        }
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
AmrDriver::writeCheckpointFile() const
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::writeCheckpointfile" << endl;
    }

    CH_TIME("AmrDriver::writeCheckpointFile");

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
AmrDriver::readCheckpointFile(HDF5Handle& a_handle)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::readCheckpointFile" << endl;
    }

    CH_TIME("AmrDriver::readCheckpointFile");

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

            m_head[lev] = new LevelData<FArrayBox>(levelDBL, SpaceDim, nGhost);

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
AmrDriver::restart(string& a_restart_file)
{
    if (m_verbosity > 3)
    {
        pout() << "AmrDriver::restart" << endl;
    }

    HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
    // first read in data from checkpoint file
    readCheckpointFile(handle);
    handle.close();
    // don't think I need to do anything else, do I?
}

#endif

#include "NamespaceFooter.H"
