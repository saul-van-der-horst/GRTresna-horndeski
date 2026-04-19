/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"
#include "DerivativeOperators.hpp"
#include "EMTensor.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "IntVect.H"
#include "LevelData.H"
#include "MultigridVariables.hpp"
#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"

void ScalarField::initialise_matter_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                         const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // Iterate over the box and set non zero comps
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            multigrid_vars_box(iv, c_phi_0) = my_phi_function(loc);
            multigrid_vars_box(iv, c_Pi_0) = my_Pi_function(loc);
        }
    }
}

// template <class data_t>
emtensor_t ScalarField::compute_emtensor(const IntVect a_iv,
                                         const RealVect &a_dx,
                                         FArrayBox &a_multigrid_vars_box) const
{
    emtensor_t out;

    DerivativeOperators derivs(a_dx);
    RealVect loc;
    Grids::get_loc(loc, a_iv, a_dx, center);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh = psi_and_Aij_functions->compute_bowenyork_psi(loc);
    Real psi_0 = psi_reg + psi_bh;
    Real Pi_0 = a_multigrid_vars_box(a_iv, c_Pi_0);
    Real phi_0 = a_multigrid_vars_box(a_iv, c_phi_0);
    Real K_0 = a_multigrid_vars_box(a_iv, c_K_0);
    Tensor<2, Real, SpaceDim> h;
    FOR2(i, j){h[i][j] = a_multigrid_vars_box(a_iv, c_h + index_ij(i, j));}
    Tensor<2, Real, SpaceDim> h_UU = compute_inverse(h);
    Tensor<2, Tensor<1, Real, SpaceDim>, SpaceDim> d1_h;
    derivs.get_d1(d1_h, a_iv, a_multigrid_vars_box, c_h);
    auto chris = compute_christoffel(d1_h, h_UU);
    Real chi = pow(psi_0, -4.0);
    
    Tensor<1, Real, SpaceDim> d1_phi;
    derivs.get_d1(d1_phi, a_iv, a_multigrid_vars_box, c_phi_0);
    Real d1_phi_squared = 0;
    FOR1(i) { d1_phi_squared += d1_phi[i] * d1_phi[i]; }
    Tensor<2, Real, SpaceDim> d2_phi;
    derivs.get_d2(d2_phi, a_iv, a_multigrid_vars_box, c_phi_0);
    Real Xplus = 0.5 * (Pi_0 * Pi_0 + pow(psi_0, -4.0) * d1_phi_squared);
    Real X     = 0.5 * (Pi_0 * Pi_0 - pow(psi_0, -4.0) * d1_phi_squared);
    Tensor<1, Real, SpaceDim> d1_psi;
    derivs.get_d1(d1_psi, a_iv, a_multigrid_vars_box, c_psi_reg);
    Tensor<1, Real, SpaceDim> d1_chi;
    FOR1(i) { d1_chi[i] = -4.0 * pow(psi_0, -5.0) * d1_psi[i]; }
    Real dphi_dot_dchi = 0.;
    FOR1(i) { dphi_dot_dchi += d1_phi[i] * d1_chi[i];}
    
    Tensor<2, Real, SpaceDim> covd2phi;
    FOR2(i, j)
    {   covd2phi[i][j] = d2_phi[i][j];}
    Tensor<2, Real, SpaceDim> tau_ij;
    FOR2(i, j)
    {tau_ij[i][j] =
             K_0 * Pi_0 * delta[i][j] / 3. +
            0.5 * (delta[i][j] * dphi_dot_dchi + d1_phi[i] * d1_chi[j] +
                   d1_phi[j] * d1_chi[i] +
                   chi * (covd2phi[i][j] + covd2phi[j][i]));
    }
    Real tau= compute_trace(quantities.tau_ij, h_UU)
    Tensor<2, Real, SpaceDim> tau_i;
    For2(i,j)
    {tau_i[i]=K_0*d1_phi[i] / 3. + d1_Pi[i]+d1_phi[j]*tau_ij[i][j]}
    Real tau_i_dot_dphi = 0.;
    FOR(i, j)
    {
        tau_i_dot_dphi +=
            delta[i][j] * d1_phi[i] * tau_i[j];
    }
    Real tau_ij_dot_dphi2 = 0.;
    FOR(i, j)
    {
        tau_ij_dot_dphi2 +=
            delta[i][j] * d1.phi[i] * tau_ij_dot_dphi[j];
    }
    CouplingAndPotential coupling_and_potential(m_p.coupling_and_potential_params);
    Real V        = coupling_and_potential.V(phi_0, X);
    Real g2       = coupling_and_potential.G2(phi_0, X);
    Real dg2_dX   = coupling_and_potential.dG2_dX(phi_0, X);
    Real dg3_dX   = coupling_and_potential.dG3_dX(phi_0, X);
    Real dg3_dphi = coupling_and_potential.dG3_dphi(phi_0, X);

    out.rho =
        dg3_dX * (tau * Pi_0 * Pi_0 -
                                   tau_ij_dot_dphi2 * chi) +
              dg3_dphi * 2. * Xplus +
              dg2_dX * Pi_0 * Pi_0 - g2 + Xplus +
              V;
    FOR1(i) { out.Si[i] = dg3_dX *
                        (-tau * Pi_0 * d1_phi[i] -
                         Pi_0 * Pi_0 * tau_i[i] +
                         d1_phi[i] * tau_i_dot_dphi * chi +
                         Pi_0 * tau_ij_dot_dphi[i]) -
                    Pi_0 * d1_phi[i] *
                        (1. + dg2_dX + 2. * dg3_dphi); }

    return out;
}
