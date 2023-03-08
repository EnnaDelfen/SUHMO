#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//
//  FileIBC.H
// ============
//
// PhysIBC-derived class for geometries read from an hdf5 file
//

#ifndef _FILEIBC_H_
#define _FILEIBC_H_

#include <string>
#include "HydroIBC.H"
#include "suhmo_params.H"

#include "NamespaceHeader.H"

/// Physical/domain initial and boundary conditions
/**

 */
class FileIBC : public HydroIBC
{
public:
    /// Constructor
    /**
     */
    FileIBC();

    /// Destructor
    /**
     */
    virtual ~FileIBC();

    /// Define the object
    /**
       Set the problem domain index space and the grid spacing for this
       initial and boundary condition object.
     */
    virtual void define(const ProblemDomain& a_domain, const Real& a_dx);

    /// Factory method - this object is its own factory
    /**
       Returns a pointer to a new FileIBC object with
       m_isDefined = false (i.e., its define() must be called before it
       is used).
     */
    virtual HydroIBC* new_hydroIBC();

    /* set up file input and names */
    void setup(std::string a_geomFile,
               std::string a_thicknessName,
               std::string a_topographyName,
               std::string a_velXFile,
               std::string a_velYFile);
  
    /// Set up initial conditions
    /**
     */
    virtual void initialize(LevelData<FArrayBox>& a_U);

    virtual void initializeData(RealVect& a_dx,
                                suhmo_params params,
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
                                LevelData<FArrayBox>& a_bumpSpacing,
                                LevelData<FArrayBox>& a_levelmagVel);

    virtual void initializePi(RealVect& a_dx,
                              suhmo_params Params,     
                              LevelData<FArrayBox>& a_Pi);

    virtual void initializePi(RealVect& a_dx,
                              suhmo_params Params,     
                              LevelData<FArrayBox>& a_head,
                              LevelData<FArrayBox>& a_gapHeight,
                              LevelData<FArrayBox>& a_Pw,
                              LevelData<FArrayBox>& a_zbed,
                              LevelData<FArrayBox>& a_Pi,
                              LevelData<FArrayBox>& a_iceHeight,
                              LevelData<FArrayBox>& a_bumpHeight,
                              LevelData<FArrayBox>& a_bumpSpacing);

    virtual void initializeBed(RealVect& a_dx,
                               suhmo_params Params,     
                               LevelData<FArrayBox>& a_zbed,
                               LevelData<FArrayBox>& a_bumpHeight,
                               LevelData<FArrayBox>& a_bumpSpacing);

    virtual void resetCovered(suhmo_params Params,
                              LevelData<FArrayBox>& a_head,
                              LevelData<FArrayBox>& a_Pi);


protected:
    /* reference ice thickness */
    RefCountedPtr<LevelData<FArrayBox> > m_refThickness;

    /* reference bedrock topography */
    RefCountedPtr<LevelData<FArrayBox> > m_refTopography;

    /* reference base (or not) velocity*/
    RefCountedPtr<LevelData<FArrayBox> > m_velX;
    RefCountedPtr<LevelData<FArrayBox> > m_velY;

    /* cell spacing of reference data - ice thick and topg */
    RealVect m_refDx;
    /* cell spacing of reference data - velo*/
    RealVect m_refDxVel;

private:
    // Disallowed for all the usual reasons
    void operator=(const FileIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }

    // Disallowed for all the usual reasons
    FileIBC(const FileIBC& a_input)
    {
        MayDay::Error("invalid operator");
    }
};

#include "NamespaceFooter.H"
#endif
