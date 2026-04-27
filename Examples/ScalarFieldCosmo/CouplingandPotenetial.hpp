#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include <string>
#include <stdexcept>
#include <memory>


class CouplingAndPotential
{
  public:
    struct params_t
    {
        mutable double scalar_mass = 0.0;
        double g2 = 0.0;
        double g3 = 0.0;
    };
    params_t m_params;

    CouplingAndPotential(params_t a_params) : m_params(a_params) {}
    virtual ~CouplingAndPotential() = default;

  protected:
    virtual double V_impl           (double phi, double X) const { return 0.; }
    virtual double dV_dphi_impl     (double phi, double X) const { return 0.; }
    virtual double G2_impl          (double phi, double X) const { return 0.; }
    virtual double dG2_dphi_impl    (double phi, double X) const { return 0.; }
    virtual double dG2_dX_impl      (double phi, double X) const { return 0.; }
    virtual double d2G2_dXX_impl    (double phi, double X) const { return 0.; }
    virtual double d2G2_dXphi_impl  (double phi, double X) const { return 0.; }
    virtual double G3_impl          (double phi, double X) const { return 0.; }
    virtual double dG3_dphi_impl    (double phi, double X) const { return 0.; }
    virtual double dG3_dX_impl      (double phi, double X) const { return 0.; }
    virtual double d2G3_dXX_impl    (double phi, double X) const { return 0.; }
    virtual double d2G3_dXphi_impl  (double phi, double X) const { return 0.; }
    virtual double d2G3_dphiphi_impl(double phi, double X) const { return 0.; }

  public:
    template <class data_t>
    ALWAYS_INLINE data_t V           (const data_t phi, const data_t X) const { return V_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t dV_dphi     (const data_t phi, const data_t X) const { return dV_dphi_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t G2          (const data_t phi, const data_t X) const { return G2_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi    (const data_t phi, const data_t X) const { return dG2_dphi_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX      (const data_t phi, const data_t X) const { return dG2_dX_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX    (const data_t phi, const data_t X) const { return d2G2_dXX_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXphi  (const data_t phi, const data_t X) const { return d2G2_dXphi_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t G3          (const data_t phi, const data_t X) const { return G3_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t dG3_dphi    (const data_t phi, const data_t X) const { return dG3_dphi_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX      (const data_t phi, const data_t X) const { return dG3_dX_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX    (const data_t phi, const data_t X) const { return d2G3_dXX_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXphi  (const data_t phi, const data_t X) const { return d2G3_dXphi_impl(phi, X); }
    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dphiphi(const data_t phi, const data_t X) const { return d2G3_dphiphi_impl(phi, X); }
};


class Canonical : public CouplingAndPotential
{
  public:
    Canonical(params_t p) : CouplingAndPotential(p) {}
};

// -- Canonical leaf classes --
class CanonicalModel1 : public Canonical
{
  public:
    CanonicalModel1(params_t p) : Canonical(p) {}
};

class CanonicalModel2 : public Canonical
{
  public:
    CanonicalModel2(params_t p) : Canonical(p) {}
};


class KEssence : public CouplingAndPotential
{
  public:
    KEssence(params_t p) : CouplingAndPotential(p) {}
};


class KEssenceModel1 : public KEssence
{
  public:
    KEssenceModel1(params_t p) : KEssence(p) {}
};

class KEssenceModel2 : public KEssence
{
  public:
    KEssenceModel2(params_t p) : KEssence(p) {}
};


class KGB : public CouplingAndPotential
{
  public:
    KGB(params_t p) : CouplingAndPotential(p) {}
};

// -- KGB leaf classes --
class KGBModel1 : public KGB
{
  public:
    KGBModel1(params_t p) : KGB(p) {}
};

class KGBModel2 : public KGB
{
  public:
    KGBModel2(params_t p) : KGB(p) {}
};


inline std::unique_ptr<CouplingAndPotential>
makeCouplingAndPotential(const std::string &model_name,
                         CouplingAndPotential::params_t params)
{
    if      (model_name == "canonical-1")  return std::make_unique<CanonicalModel1>(params);
    else if (model_name == "canonical-2")  return std::make_unique<CanonicalModel2>(params);
    else if (model_name == "k-essence-1")  return std::make_unique<KEssenceModel1>(params);
    else if (model_name == "k-essence-2")  return std::make_unique<KEssenceModel2>(params);
    else if (model_name == "kgb-1")        return std::make_unique<KGBModel1>(params);
    else if (model_name == "kgb-2")        return std::make_unique<KGBModel2>(params);
    else
        throw std::invalid_argument("Unknown model: '" + model_name +
                                    "'. Options: canonical-1, canonical-2, "
                                    "k-essence-1, k-essence-2, kgb-1, kgb-2");
}

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
