#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef SUHMO_PARM_H
#define SUHMO_PARM_H

#include "LevelData.H"
#include "FArrayBox.H"

#include "NamespaceHeader.H"

class suhmo_params { 
public:
    suhmo_params() { setDefaults(); };

    virtual ~suhmo_params() {};

    void setDefaults();

    void readInputs();

    // known params
    Real m_rho_i;       // kg.m-3
    Real m_rho_w;       // kg.m-3
    Real m_gravity;     // m.s-1
    // constants of problem
    Real m_G;             // Geothermal flux W.m-2
    Real m_L;             // Latent heat J.kg-1
    Real m_H;             // ice height m -- cst ?? 
    Real m_nu;            // Kinematic viscosity of water m2.s-1
    Real m_ct;            // Change of pressure melting point with temperature K.Pa-1
    Real m_cw;            // Heat capacity of water J.kg-1.K-1
    Real m_omega;         // Parameter controlling nonlinear transition between lam and turbu flow
    // BC
    std::vector<Real> m_ub;    // m.s-1 Sliding velocity
    Real m_br;            // height between bed bumps m
    Real m_lr;            // spacing between bed bumps m
    Real m_A;             // Flow law ??
    // Initialization 
    Real m_slope;         // bed slope
    Real m_gapInit;       // subglacial init gap height 
    Real m_ReInit;        // Initial Reynolds number 
    // Moulin
    std::vector<Real> m_moulin_position;    // m.s-1 Sliding velocity
    Real m_moulin_flux;    // m3.s-1 
    
};

#include "NamespaceFooter.H"
#endif
