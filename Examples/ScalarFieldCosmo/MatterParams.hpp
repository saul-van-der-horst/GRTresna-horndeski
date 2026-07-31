/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "CouplingAndPotential.hpp"

namespace MatterParams
{

struct params_t
{
    Real phi_0;
    Real dphi;
    Real pi_0;
    Real dpi;
    CouplingAndPotential coupling;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi", matter_params.dphi);
    pp.get("pi_0", matter_params.pi_0);
    pp.get("dpi", matter_params.dpi);
    CouplingAndPotential::params_t coupling_params;
    pp.get("scalar_mass", coupling_params.scalar_mass);
    pp.get("g2", coupling_params.g2);
    pp.get("g3", coupling_params.g3);
    matter_params.coupling = CouplingAndPotential(coupling_params);
}

}; // namespace MatterParams

#endif
