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
    Real kgb_usr_y2;
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
    pp.get("scalar_mass", matter_params.scalar_mass);
}

}; // namespace MatterParams

#endif
