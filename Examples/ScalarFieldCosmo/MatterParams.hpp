/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"

namespace MatterParams
{

struct params_t
{
    Real phi_0;
    Real dphi;
    Real pi_0;
    Real dpi;
    Real scalar_mass;
    std::string model_name;
    Real g2;
    Real g3;
    Real rbs_g3;
    Real rbs_Lambda;
    Real nu;
    Real Mpl;
    Real usr_v0;
    Real usr_y1;
    Real y2;
    Real ga3_scalar_mass;
    Real mu;
    Real exph_lambda;
    Real v;
    Real y1;
    Real eta;
    Real dbin_lambda1;
    Real dbin_Lambda;
    Real f;
    Real v0;
    Real b;
    Real gamma;
    std::unique_ptr<ICouplingAndPotential> coupling;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    // Always required
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi",  matter_params.dphi);
    pp.get("pi_0",  matter_params.pi_0);
    pp.get("dpi",   matter_params.dpi);
    pp.get("model_name", matter_params.model_name);

    if (matter_params.model_name == "kgb-default")
    {
        CouplingAndPotential<KGBDefault>::params_t m_params;

        pp.get("scalar_mass", m_params.scalar_mass);
        pp.get("g2",          m_params.g2);
        pp.get("g3",          m_params.g3);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-default", m_params);
    }
    else if (matter_params.model_name == "kgb-rbs")
    {
        CouplingAndPotential<KGBRunning_braiding_starobinsky>::params_t m_params;

        pp.get("rbs_g3",     m_params.rbs_g3);
        pp.get("rbs_Lambda", m_params.rbs_Lambda);
        pp.get("nu",  m_params.nu);
        pp.get("Mpl",    m_params.Mpl);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-rbs", m_params);
    }
    else if (matter_params.model_name == "kgb-usr")
    {
        CouplingAndPotential<KGBUltra_slow_roll>::params_t m_params;

        pp.get("usr_v0", m_params.usr_v0);
        pp.get("usr_y1", m_params.usr_y1);
        pp.get("y2", m_params.y2);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-usr", m_params);
    }
    else if (matter_params.model_name == "kgb-ga3")
    {
        CouplingAndPotential<KGBCubic_galileon>::params_t m_params;

        pp.get("ga3_scalar_mass", m_params.ga3_scalar_mass);
        pp.get("mu",          m_params.mu);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-ga3", m_params);
    }
    else if (matter_params.model_name == "kgb-exph")
    {
        CouplingAndPotential<KGBExponential_hilltop>::params_t m_params;

        pp.get("exph_lambda", m_params.exph_lambda);
        pp.get("v",      m_params.v);
        pp.get("y1",     m_params.y1);
        pp.get("eta",   m_params.eta);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-exph", m_params);
    }
    else if (matter_params.model_name == "kgb-dbin")
    {
        CouplingAndPotential<KGBDBI_natural>::params_t m_params;

        pp.get("dbin_lambda1", m_params.dbin_lambda1);
        pp.get("dbin_Lambda", m_params.dbin_Lambda);
        pp.get("f",      m_params.f);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-dbin", m_params);
    }
    else if (matter_params.model_name == "kgb-dbipl")
    {
        CouplingAndPotential<KGBDBI_power_law>::params_t m_params;

        pp.get("v0",    m_params.v0);
        pp.get("b",     m_params.b);
        pp.get("Mpl",   m_params.Mpl);
        pp.get("gamma", m_params.gamma);

        matter_params.coupling =
            makeCouplingAndPotential("kgb-dbipl", m_params);
    }
    else
    {
        throw std::invalid_argument("Unknown model_name: '"
                                    + matter_params.model_name + "'");
    }
}    


}; // namespace MatterParams

#endif
