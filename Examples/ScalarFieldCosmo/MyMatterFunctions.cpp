/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"
#include "CouplingAndPotential.hpp"



Real ScalarField::V(const Real phi, const Real X) const
{ return m_matter_params.coupling->V(phi, X); }

Real ScalarField::dV_dphi(const Real phi, const Real X) const
{ return m_matter_params.coupling->dV_dphi(phi, X); }

Real ScalarField::G2(const Real phi, const Real X) const
{ return m_matter_params.coupling->G2(phi, X); }

Real ScalarField::dG2_dphi(const Real phi, const Real X) const
{ return m_matter_params.coupling->dG2_dphi(phi, X); }

Real ScalarField::dG2_dX(const Real phi, const Real X) const
{ return m_matter_params.coupling->dG2_dX(phi, X); }

Real ScalarField::d2G2_dXX(const Real phi, const Real X) const
{ return m_matter_params.coupling->d2G2_dXX(phi, X); }

Real ScalarField::d2G2_dXphi(const Real phi, const Real X) const
{ return m_matter_params.coupling->d2G2_dXphi(phi, X); }

Real ScalarField::G3(const Real phi, const Real X) const
{ return m_matter_params.coupling->G3(phi, X); }

Real ScalarField::dG3_dphi(const Real phi, const Real X) const
{ return m_matter_params.coupling->dG3_dphi(phi, X); }

Real ScalarField::dG3_dX(const Real phi, const Real X) const
{ return m_matter_params.coupling->dG3_dX(phi, X); }

Real ScalarField::d2G3_dXX(const Real phi, const Real X) const
{ return m_matter_params.coupling->d2G3_dXX(phi, X); }

Real ScalarField::d2G3_dXphi(const Real phi, const Real X) const
{ return m_matter_params.coupling->d2G3_dXphi(phi, X); }

Real ScalarField::d2G3_dphiphi(const Real phi, const Real X) const
{ return m_matter_params.coupling->d2G3_dphiphi(phi, X); }

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real L = domainLength[0];
    Real dphi_value = m_matter_params.dphi / 3. *
                      (sin(2 * M_PI * loc[0] / L) + sin(2 * M_PI * loc[1] / L) +
                       sin(2 * M_PI * loc[2] / L));
    return m_matter_params.phi_0 + dphi_value;
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real L = domainLength[0];
    Real dpi_value = m_matter_params.dpi / 3. *
                     (sin(2 * M_PI * loc[0] / L) + sin(2 * M_PI * loc[1] / L) +
                      sin(2 * M_PI * loc[2] / L));
    return m_matter_params.pi_0 + dpi_value;
}
