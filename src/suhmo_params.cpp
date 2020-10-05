#include <suhmo_params.H>
#include "ParmParse.H"

#include "NamespaceHeader.H"

void suhmo_params::setDefaults()
{
    m_rho_i = 910.0;
    m_rho_w = 1000.0;
    m_gravity = 9.8;
    // constants of problem
    m_G = 0.0;
    m_L = 0.0;   
    m_H = 0.0;  
    m_nu = 0.0;
    m_ct = 0.0;
    m_cw = 0.0;
    m_omega = 0.0;
    // BC
    //m_ub ??
    m_br = 0.0 ;
    m_lr = 0.0;
    m_A = 0.0;
    // Initialization 
    m_slope = 0.0;
    m_gapInit = 0.0;
    m_ReInit = 0.0;
}

void suhmo_params::readInputs()
{
    pout() << "suhmo_params::readInputs()" << endl;

    // Cst of problem
    m_rho_i = 910.0;
    m_rho_w = 1000.0;
    m_gravity = 9.8;
    // Problem dependent variables
    ParmParse ppParams("suhmo");
    ppParams.get("GeoFlux", m_G);
    ppParams.get("LatHeat", m_L);
    ppParams.get("IceHeight", m_H);
    ppParams.get("WaterViscosity", m_nu);
    ppParams.get("ct", m_ct);
    ppParams.get("cw", m_cw);
    ppParams.get("turbulentParam", m_omega);
    m_ub.resize(SpaceDim,0.0);
    ppParams.getarr("SlidingVelocity", m_ub, 0, SpaceDim);
    ppParams.get("br", m_br);
    ppParams.get("lr", m_lr);
    ppParams.get("A", m_A);
    ppParams.get("slope", m_slope);
    ppParams.get("GapInit", m_gapInit);
    ppParams.get("ReInit", m_ReInit);
}

#include "NamespaceFooter.H"
