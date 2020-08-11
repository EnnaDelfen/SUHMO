#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef SUHMO_H
#define SUHMO_H

namespace hydro_params { 
    // known params
    double m_rho_i;       // kg.m-3
    double m_rho_w;       // kg.m-3
    double m_gravity;     // m.s-1
    // constants of problem
    Real m_G;             // Geothermal flux W.m-2
    Real m_L;             // Latent heat J.kg-1
    Real m_H;             // ice height m -- cst ?? 
    Real m_nu;            // Kinematic viscosity of water m2.s-1
    Real m_ct;            // Change of pressure melting point with temperature K.Pa-1
    Real m_cw;            // Heat capacity of water J.kg-1.K-1
    Real m_omega;         // Parameter controlling nonlinear transition between lam and turbu flow
    // BC
    Vector<Real> m_ub;            // m.s-1 Sliding velocity
    Real m_br;            // height between bed bumps m
    Real m_lr;            // spacing between bed bumps m
    Real m_A;             // Flow law ??
    // Initialization 
    Real m_slope;         // bed slope
    Real m_gapInit;       // subglacial init gap height 
    Real m_ReInit;        // Initial Reynolds number 
}
#endif
