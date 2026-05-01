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
    Real usr_y2;
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
    Real Mpl;
    Real gamma;
    std::unique_ptr<ICouplingAndPotential> coupling;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi", matter_params.dphi);
    pp.get("pi_0", matter_params.pi_0);
    pp.get("dpi", matter_params.dpi);
    pp.get("Mpl",    matter_params.Mpl);
    pp.get("model_name", matter_params.model_name);
    pp.get("scalar_mass", matter_params.scalar_mass);
    pp.get("g2", matter_params.g2);
    pp.get("g3", matter_params.g3);
    pp.get("rbs_g3",     matter_params.rbs_g3);
    pp.get("rbs_Lambda", matter_params.rbs_Lambda);
    pp.get("eta",  matter_params.eta);
    pp.get("usr_v0", matter_params.usr_v0);
    pp.get("usr_y1", matter_params.usr_y1);
    pp.get("y2", matter_params.y2);
    pp.get("mu", matter_params.mu);
    pp.get("exph_lambda", matter_params.exph_lambda);
    pp.get("v",      matter_params.v);
    pp.get("y1",     matter_params.y1);
    pp.get("eta",   matter_params.eta);
   
    pp.get("dbin_lambda1", matter_params.dbin_lambda1);
    pp.get("dbin_Lambda", matter_params.dbin_Lambda);
    pp.get("f",      matter_params.f);
    
    pp.get("v0",    matter_params.v0);
    pp.get("b",     matter_params.b);
    pp.get("gamma", matter_params.gamma);

    // Build m_params for the chosen model then create coupling object
    CouplingAndPotential<KGBDefault>::params_t m_params;

    if (matter_params.model_name == "kgb-default")
    {
        m_params.scalar_mass = matter_params.scalar_mass;
        m_params.g2          = matter_params.g2;
        m_params.g3          = matter_params.g3;
    }
    else if (matter_params.model_name == "kgb-rbs")
    {
        m_params.g3     = matter_params.rbs_g3;
        m_params.Lambda = matter_params.rbs_Lambda;
        m_params.alpha  = matter_params.nu;
        m_params.Mpl    = matter_params.Mpl;
    }
    else if (matter_params.model_name == "kgb-usr")
    {
        m_params.v0 = matter_params.usr_v0;
        m_params.y1 = matter_params.usr_y1;
        m_params.y2 = matter_params.y2;
    }
    else if (matter_params.model_name == "kgb-ga3")
    {
        m_params.scalar_mass = matter_params.ga3_scalar_mass;
        m_params.mu          = matter_params.mu;
    }
    else if (matter_params.model_name == "kgb-exph")
    {
        m_params.lambda = matter_params.exph_lambda;
        m_params.v      = matter_params.v;
        m_params.y1     = matter_params.exph_y1;
        m_params.beta   = matter_params.eta;
    }
    else if (matter_params.model_name == "kgb-dbin")
    {
        m_params.lambda = matter_params.dbin_lambda1;
        m_params.Lambda = matter_params.dbin_Lambda;
        m_params.f      = matter_params.f;
    }
    else if (matter_params.model_name == "kgb-dbipl")
    {
        m_params.v0    = matter_params.v0;
        m_params.b     = matter_params.b;
        m_params.Mpl   = matter_params.Mpl;
        m_params.gamma = matter_params.gamma;
    }
    else
    {
        throw std::invalid_argument("Unknown model_name: '"
                                    + matter_params.model_name + "'");
    }

    matter_params.coupling = makeCouplingAndPotential(
        matter_params.model_name, m_params);

    
}

}; // namespace MatterParams

#endif
