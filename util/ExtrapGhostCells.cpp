#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ExtrapGhostCells.H"
#include "ExtrapBCF_F.H"

#include "NamespaceHeader.H"

void 
ExtrapGhostCellsEC(LevelData<FluxBox>& a_phi,
                   const ProblemDomain& a_domain)
{
    DataIterator dit = a_phi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FluxBox& phi = a_phi[dit];
        for (int dir=0; dir<SpaceDim; dir++) { 
            FArrayBox& dirPhi = phi[dir];
            ExtrapGhostCells(dirPhi, a_domain,
                             a_phi.ghostVect(), dir);
        }
    }
}

void 
CopyGhostCellsEC(LevelData<FluxBox>& a_phi,
                 const ProblemDomain& a_domain)
{
    DataIterator dit = a_phi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        FluxBox& phi = a_phi[dit];
        for (int dir=0; dir<SpaceDim; dir++) { 
            FArrayBox& dirPhi = phi[dir];
            CopyGhostCells(dirPhi, a_domain,
                           a_phi.ghostVect(), dir);
        }
    }
}

void 
ExtrapGhostCells(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain)
{
    DataIterator dit = a_phi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
          ExtrapGhostCells(a_phi[dit], a_domain, 
                           a_phi.ghostVect(), -1);
    }
}

void 
CopyGhostCells(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain)
{
    DataIterator dit = a_phi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
          CopyGhostCells(a_phi[dit], a_domain, 
                         a_phi.ghostVect(), -1);
    }
}

void 
NullGhostCells(LevelData<FArrayBox>& a_phi,
                 const ProblemDomain& a_domain)
{
    DataIterator dit = a_phi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
          NullGhostCells(a_phi[dit], a_domain, 
                         a_phi.ghostVect(), -1);
    }
}



void 
ExtrapGhostCellsEC(FluxBox& a_phi,
                   const ProblemDomain& a_domain,
                   const IntVect& a_ghostVect)
{
    for (int dir=0; dir<SpaceDim; dir++) { 
        FArrayBox& dirPhi = a_phi[dir];
        ExtrapGhostCells(dirPhi, a_domain,
                         a_ghostVect, dir);
    }
}

void 
ExtrapGhostCells(FArrayBox& a_phi,
                 const ProblemDomain& a_domain,
                 const IntVect& a_ghostVect,
                 int dir_edges)
{

  Box domainBox = a_domain.domainBox();
  if (dir_edges > -1) {
      IntVect conv_type = IntVect::Zero; 
      conv_type[dir_edges] = 1;
      const IntVect conv_type_cst = conv_type;
      domainBox.convert(conv_type_cst);
  }

  for (int dir=0; dir<SpaceDim; dir++) {

      if ((!a_domain.isPeriodic(dir)) || (domainBox.type() != IntVect::Zero)) {

          int rad = a_ghostVect[dir];

          // lo-side 
          int hiLo = 0;
          // slc :: we need to do the low side one strip at 
          // at time so that strip n is filled before we
          // extrapolate from it into strip n -1.
          Box  stripLo;
          if ( (domainBox.type() != IntVect::Zero) && (dir==dir_edges) )  {
              stripLo = bdryLo(domainBox,  
                               dir, 1);

              IntVect shiftVect = IntVect::Zero; 
              shiftVect[dir] = -1;
              stripLo.shift(shiftVect);
          } else { 
              stripLo = adjCellLo(domainBox,  
                                  dir, 1);
          }

          // do this to try to catch corner cells
          stripLo.grow(rad);
          stripLo.grow(dir,-rad);

          for (int i =0; i < rad; ++i) {
              Box ghostBoxLo = stripLo;
              ghostBoxLo &= a_phi.box();

              if (!ghostBoxLo.isEmpty())
              {
                  FORT_SIMPLEEXTRAPBC(CHF_FRA(a_phi),
                                      CHF_BOX(ghostBoxLo),
                                      CHF_INT(dir), 
                                      CHF_INT(hiLo));
              }
              stripLo.shift(dir,-1);
          }
          
          // hi-side
          hiLo = 1;
          Box ghostBoxHi;
          if ( (domainBox.type() != IntVect::Zero) && (dir==dir_edges) )  {
              ghostBoxHi = bdryHi(domainBox,  
                                  dir, rad);

              IntVect shiftVect = IntVect::Zero; 
              shiftVect[dir] = 1;
              ghostBoxHi.shift(shiftVect);
          } else { 
              ghostBoxHi = adjCellHi(domainBox, 
                                     dir, rad);
          }
          // do this to try to catch corner cells
          ghostBoxHi.grow(1);
          ghostBoxHi.grow(dir,-1);
          ghostBoxHi &= a_phi.box();

          if(!ghostBoxHi.isEmpty())
          {
              FORT_SIMPLEEXTRAPBC(CHF_FRA(a_phi),
                                  CHF_BOX(ghostBoxHi),
                                  CHF_INT(dir), 
                                  CHF_INT(hiLo));
          }
          
      }
  } // end loop over directions
}

void 
CopyGhostCells(FArrayBox& a_phi,
               const ProblemDomain& a_domain,
               const IntVect& a_ghostVect,
               int dir_edges)
{

  Box domainBox = a_domain.domainBox();
  if (dir_edges > -1) {
      IntVect conv_type = IntVect::Zero; 
      conv_type[dir_edges] = 1;
      const IntVect conv_type_cst = conv_type;
      domainBox.convert(conv_type_cst);
  }

  for (int dir=0; dir<SpaceDim; dir++) {

      if ((!a_domain.isPeriodic(dir)) || (domainBox.type() != IntVect::Zero)) {

          int rad = a_ghostVect[dir];

          // lo-side 
          int hiLo = 0;
          // slc :: we need to do the low side one strip at 
          // at time so that strip n is filled before we
          // extrapolate from it into strip n -1.
          Box  stripLo;
          if ( (domainBox.type() != IntVect::Zero) && (dir==dir_edges) )  {
              stripLo = bdryLo(domainBox,  
                               dir, 1);

              IntVect shiftVect = IntVect::Zero; 
              shiftVect[dir] = -1;
              stripLo.shift(shiftVect);
          } else { 
              stripLo = adjCellLo(domainBox,  
                                  dir, 1);
          }

          // do this to try to catch corner cells
          stripLo.grow(rad);
          stripLo.grow(dir,-rad);

          for (int i =0; i < rad; ++i) {
              Box ghostBoxLo = stripLo;
              ghostBoxLo &= a_phi.box();

              if (!ghostBoxLo.isEmpty())
              {
                  FORT_SIMPLECOPYBC(CHF_FRA(a_phi),
                                      CHF_BOX(ghostBoxLo),
                                      CHF_INT(dir), 
                                      CHF_INT(hiLo));
              }
              stripLo.shift(dir,-1);
          }
          
          // hi-side
          hiLo = 1;
          Box ghostBoxHi;
          if ( (domainBox.type() != IntVect::Zero) && (dir==dir_edges) )  {
              ghostBoxHi = bdryHi(domainBox,  
                                  dir, rad);

              IntVect shiftVect = IntVect::Zero; 
              shiftVect[dir] = 1;
              ghostBoxHi.shift(shiftVect);
          } else { 
              ghostBoxHi = adjCellHi(domainBox, 
                                     dir, rad);
          }
          // do this to try to catch corner cells
          ghostBoxHi.grow(1);
          ghostBoxHi.grow(dir,-1);

          // AF: why no loop on rad in hi side ?
          ghostBoxHi &= a_phi.box();

          if(!ghostBoxHi.isEmpty())
          {
              FORT_SIMPLECOPYBC(CHF_FRA(a_phi),
                                CHF_BOX(ghostBoxHi),
                                CHF_INT(dir), 
                                CHF_INT(hiLo));
          }
          
      }
  } // end loop over directions
}

void 
NullGhostCells(FArrayBox& a_phi,
               const ProblemDomain& a_domain,
               const IntVect& a_ghostVect,
               int dir_edges)
{

  Box domainBox = a_domain.domainBox();
  if (dir_edges > -1) {
      IntVect conv_type = IntVect::Zero; 
      conv_type[dir_edges] = 1;
      const IntVect conv_type_cst = conv_type;
      domainBox.convert(conv_type_cst);
  }

  for (int dir=0; dir<SpaceDim; dir++) {

      if ((!a_domain.isPeriodic(dir)) || (domainBox.type() != IntVect::Zero)) {

          int rad = a_ghostVect[dir];

          // lo-side 
          int hiLo = 0;
          // slc :: we need to do the low side one strip at 
          // at time so that strip n is filled before we
          // extrapolate from it into strip n -1.
          Box  stripLo;
          if ( (domainBox.type() != IntVect::Zero) && (dir==dir_edges) )  {
              stripLo = bdryLo(domainBox,  
                               dir, 1);

              IntVect shiftVect = IntVect::Zero; 
              shiftVect[dir] = -1;
              stripLo.shift(shiftVect);
          } else { 
              stripLo = adjCellLo(domainBox,  
                                  dir, 1);
          }

          // do this to try to catch corner cells
          stripLo.grow(rad);
          stripLo.grow(dir,-rad);

          for (int i =0; i < rad; ++i) {
              Box ghostBoxLo = stripLo;
              ghostBoxLo &= a_phi.box();

              if (!ghostBoxLo.isEmpty())
              {
                  FORT_NULLBC(CHF_FRA(a_phi),
                              CHF_BOX(ghostBoxLo),
                              CHF_INT(dir), 
                              CHF_INT(hiLo));
              }
              stripLo.shift(dir,-1);
          }
          
          // hi-side
          hiLo = 1;
          Box ghostBoxHi;
          if ( (domainBox.type() != IntVect::Zero) && (dir==dir_edges) )  {
              ghostBoxHi = bdryHi(domainBox,  
                                  dir, rad);

              IntVect shiftVect = IntVect::Zero; 
              shiftVect[dir] = 1;
              ghostBoxHi.shift(shiftVect);
          } else { 
              ghostBoxHi = adjCellHi(domainBox, 
                                     dir, rad);
          }
          // do this to try to catch corner cells
          ghostBoxHi.grow(1);
          ghostBoxHi.grow(dir,-1);

          // AF: why no loop on rad in hi side ?
          ghostBoxHi &= a_phi.box();

          if(!ghostBoxHi.isEmpty())
          {
              FORT_NULLBC(CHF_FRA(a_phi),
                          CHF_BOX(ghostBoxHi),
                          CHF_INT(dir), 
                          CHF_INT(hiLo));
          }
          
      }
  } // end loop over directions
}

#include "NamespaceFooter.H"
