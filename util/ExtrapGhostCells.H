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
//  ExtrapGhostCells.H
// ============
//
// Simple function to fill domain ghost cells with extrapolated values
//

#ifndef _EXTRAPGHOSTCELLS_H_
#define _EXTRAPGHOSTCELLS_H_

#include "ProblemDomain.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

/* Higher level */
void 
ExtrapGhostCellsEC(LevelData<FluxBox>& a_phi,
                   const ProblemDomain& a_domain);

void 
ExtrapGhostCells(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain);

void 
CopyGhostCellsEC(LevelData<FluxBox>& a_phi,
                 const ProblemDomain& a_domain);

void 
CopyGhostCells(LevelData<FArrayBox>& a_phi,
               const ProblemDomain& a_domain);

void 
NullGhostCells(LevelData<FArrayBox>& a_phi,
               const ProblemDomain& a_domain);



void 
ExtrapGhostCellsEC(FluxBox& a_phi,
                   const ProblemDomain& a_domain,
                   const IntVect& a_ghostVect);

void 
ExtrapGhostCells(FArrayBox& a_phi,
                 const ProblemDomain& a_domain,
                 const IntVect& a_ghostVect,
                 int dir_edges);

void 
CopyGhostCells(FArrayBox& a_phi,
               const ProblemDomain& a_domain,
               const IntVect& a_ghostVect,
               int dir_edges);

void 
NullGhostCells(FArrayBox& a_phi,
               const ProblemDomain& a_domain,
               const IntVect& a_ghostVect,
               int dir_edges);

#include "NamespaceFooter.H"

#endif
