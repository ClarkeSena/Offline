#ifndef STMMC_LaBrIntrinsicDecays_hh
#define STMMC_LaBrIntrinsicDecays_hh
//
// Decay kinematics for the intrinsic radioactivity of a LaBr3 crystal, shared by
// La138Source_module.cc and LaBrIntrinsicSource_module.cc.
//
// Each decay is emitted as its PRODUCTS (gammas, electrons, alphas) from one vertex at one
// time, so Geant4 transport produces the photopeak/escape/sum structure and coincident
// deposits sum naturally. No G4 radioactive decay is used.
//
//   La-138 (T1/2 1.03e11 yr), two branches:
//     EC   -> Ba-138*: 1435.8 keV gamma + full atomic relaxation of the Ba vacancy
//             (K: X-ray + remaining cascade energy, or all-Auger; L; M). A fully absorbed
//             K-capture event therefore sums to 1435.8 + 37.4 = 1473.2 keV, L capture to
//             1441.8 keV.
//     beta-> Ce-138*: 788.7 keV gamma + beta electron (Q-788.7 = 263.7 keV endpoint), with a
//             2nd-forbidden-unique shape factor and Coulomb (Fermi) correction.
//
//   Ac-227 chain (contaminant; secular equilibrium by default):
//     Ac-227 (b-, 1.4% a) -> Th-227 (a) -> Ra-223 (a) -> Rn-219 (a) -> Po-215 (a)
//       -> Pb-211 (b-) -> Bi-211 (a) -> Tl-207 (b-) -> Pb-207
//     Each alpha line feeds a daughter level; the level energy is emitted as ONE gamma at
//     the same vertex (cascade and internal conversion are not resolved). This conserves the
//     deposited energy when everything is absorbed and approximates the escape probability.
//     The small Fr-223 (1.38%) and Bi-211 b- -> Po-211 (0.28%) branches are omitted.
//     The Rn-219 -> Po-215 (1.78 ms) and later delays are all >> one micropulse, so each
//     member is generated as an independent Poisson process.
//
// Alphas deposit their full energy in Geant4; the scintillation light quenching (alpha/beta
// ratio 0.459 at 6 MeV for the Mu2e crystal, bench-measured; ~0.3 in the literature) is
// applied downstream in the digitizer and truth modules (LaBrTruthHelpers.hh, prolog.fcl).
//
// Decay data sources (checked 2026-09-24):
//   [ENSDF]  Evaluated Nuclear Structure Data File, retrieved through the IAEA LiveChart of
//            Nuclides API (nds.iaea.org/relnsd). A=138 J. Chen (2017 cut-off); Ac-227, Th-227
//            ICTP-2014 group; Ra-223 E. Browne; Rn-219 B. Singh et al.; Po-215 J.K. Tuli et
//            al.; Pb-211 A. Sonzogni et al.; Bi-211 E.A. McCutchan et al.; Tl-207 F.G. Kondev,
//            S. Lalkovski. Q values from the same ground-state tables.
//   [DDEP]   X. Mougeot, M.M. Be, LNHB Table de Radionucleides, La-138 (2016): capture
//            probabilities P_K, P_L, P_M, Ba fluorescence yield and K X-ray relative intensities.
//   [XDB]    X-ray Data Booklet (LBNL), Table 1-1: Ba electron binding energies.
//   [PT]     Sandler et al., PRC 100 (2019) 014308: La-138 beta- Q = 1052.42(41) keV.
// Alpha lines with intensity >= 0.1% are listed; weaker lines are absorbed by renormalisation.
// Original author: Sena Clarke
//

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"

namespace mu2e {
  namespace LaBrIntrinsic {

    constexpr double kElectronMass  = 0.51099895;    // MeV
    constexpr double kAlphaMass     = 3727.3794;     // MeV
    constexpr double kFineStructure = 1.0 / 137.035999;
    constexpr double kAvogadro      = 6.02214076e23;
    constexpr double kSecPerYear    = 3.1557600e7;   // Julian year

    // ---------------------------------------------------------------------------------------
    // Geometry and activity helpers
    // ---------------------------------------------------------------------------------------

    // Solid cylinder volume [cm^3] from radius and half-length in mm.
    inline double cylinderVolume_cm3(double radius_mm, double halfLength_mm) {
      return M_PI * std::pow(radius_mm / 10.0, 2) * (2.0 * halfLength_mm / 10.0);
    }

    // La-138 specific activity [Bq/cm^3] from crystal composition.
    //   laPerMolecule < 1 accounts for Ce (and Sr) substituting on La sites.
    inline double la138SpecificActivity(double density_gcm3, double molarMass_gmol,
                                        double laPerMolecule, double la138Abundance,
                                        double halfLife_yr) {
      const double nLa138PerCm3 =
        density_gcm3 / molarMass_gmol * kAvogadro * laPerMolecule * la138Abundance;
      return std::log(2.0) / (halfLife_yr * kSecPerYear) * nLa138PerCm3;
    }

    // Uniform vertex in a solid cylinder with its axis along z.
    inline CLHEP::Hep3Vector sampleCylinder(CLHEP::RandFlat& flat, const CLHEP::Hep3Vector& centre,
                                            double radius, double halfLength) {
      const double r = radius * std::sqrt(flat.fire());   // uniform per unit area
      const double phi = 2.0 * M_PI * flat.fire();
      const double dz = (2.0 * flat.fire() - 1.0) * halfLength;
      return centre + CLHEP::Hep3Vector(r * std::cos(phi), r * std::sin(phi), dz);
    }

    inline CLHEP::Hep3Vector isotropicDir(CLHEP::RandFlat& flat) {
      const double cosTheta = 2.0 * flat.fire() - 1.0;
      const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
      const double phi = 2.0 * M_PI * flat.fire();
      return CLHEP::Hep3Vector(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    }

    // ---------------------------------------------------------------------------------------
    // Particle emission
    // ---------------------------------------------------------------------------------------

    // Everything one decay needs to write its products.
    struct Emitter {
      CLHEP::RandFlat& flat;
      GenParticleCollection& out;
      GenId genId;
      CLHEP::Hep3Vector pos;
      double time;

      void particle(PDGCode::type pdg, double mass, double KE) {
        if (KE <= 0.0) return;
        const double E = KE + mass;
        const double p = std::sqrt(std::max(0.0, E * E - mass * mass));
        const CLHEP::Hep3Vector dir = isotropicDir(flat);
        out.push_back(GenParticle(pdg, genId, pos, CLHEP::HepLorentzVector(p * dir, E), time));
      }
      void gamma(double E)     { particle(PDGCode::gamma, 0.0, E); }
      void electron(double KE) { particle(PDGCode::e_minus, kElectronMass, KE); }
      void alpha(double KE)    { particle(PDGCode::alpha, kAlphaMass, KE); }
    };

    // ---------------------------------------------------------------------------------------
    // Beta spectrum: N(W) = p W q^2 F(Z,W) S(p,q), sampled by inverse CDF from a table.
    //   W = total electron energy in m_e units, p = momentum, q = W0 - W (neutrino momentum).
    //   F: relativistic point-charge Fermi function evaluated at the nuclear radius R,
    //        F = 2(1+g) (2pR)^(2g-2) exp(pi eta) |Gamma(g + i eta)|^2 / Gamma(2g+1)^2,
    //        g = sqrt(1 - (alpha Z)^2), eta = alpha Z W / p, R = 1.2 A^(1/3) fm (m_e c units).
    //      Atomic screening and the finite-size (L0) and radiative corrections are not
    //      included; screening matters below ~75 keV for La-138 [Quarati 2012].
    //      With aDaughter = 0 the non-relativistic form 2 pi eta / (1 - exp(-2 pi eta)) is used.
    //   S: unique shape factors with Coulomb functions lambda_k = F_{k-1} / F_0:
    //      1 (allowed), q^2 + lambda_2 p^2 (1st-forbidden unique),
    //      q^4 + 10/3 lambda_2 p^2 q^2 + lambda_3 p^4 (2nd-forbidden unique, La-138).
    //      F_{k-1} is the generalized Fermi function (Behrens-Janecke form, k = 1 gives F).
    //      lambda_k = 1 when aDaughter = 0.
    //   Check against ENSDF mean beta energies: La-138 97.9 keV (ENSDF 97.5(37)), Pb-211 g.s.
    //   472.6 (471.3), Tl-207 g.s. 494.4 (492.5), Ac-227 g.s. 11.4 (11.37).
    // ---------------------------------------------------------------------------------------

    // log Gamma(z) for complex z with Re(z) > 0.5 (Lanczos, g = 7, 9 terms).
    inline std::complex<double> logGamma(std::complex<double> z) {
      static const double c[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                                  771.32342877765313, -176.61502916214059, 12.507343278686905,
                                  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
      z -= 1.0;
      std::complex<double> x = c[0];
      for (int i = 1; i < 9; ++i) x += c[i] / (z + double(i));
      const std::complex<double> t = z + 7.5;
      return 0.5 * std::log(2.0 * M_PI) + (z + 0.5) * std::log(t) - t + std::log(x);
    }

    // Generalized Fermi function F_{k-1}(Z, W); k = 1 is the ordinary Fermi function.
    //   F_{k-1} = (k+g_k)/(2k) [k (2k-1)!!]^2 4^k (2pR)^(2(g_k-k)) exp(pi eta)
    //             |Gamma(g_k + i eta)|^2 / Gamma(1 + 2 g_k)^2,   g_k = sqrt(k^2 - (alpha Z)^2)
    inline double fermiFunction(int Z, int A, double W, int k = 1) {
      if (Z == 0) return 1.0;
      const double p = std::sqrt(W * W - 1.0);
      const double eta = Z * kFineStructure * W / p;
      if (A <= 0) {
        const double x = 2.0 * M_PI * eta;
        return x / (1.0 - std::exp(-x));
      }
      const double gk = std::sqrt(double(k * k) - std::pow(Z * kFineStructure, 2));
      const double R = 1.2 * std::cbrt(double(A)) / 386.15926796;   // fm -> hbar/(m_e c)
      double dfact = 1.0;                                           // (2k-1)!!
      for (int j = 2 * k - 1; j > 0; j -= 2) dfact *= j;
      const double lnF = std::log((k + gk) / (2.0 * k)) + 2.0 * std::log(k * dfact)
                       + k * std::log(4.0) + 2.0 * (gk - k) * std::log(2.0 * p * R)
                       + M_PI * eta + 2.0 * logGamma({gk, eta}).real()
                       - 2.0 * std::lgamma(1.0 + 2.0 * gk);
      return std::exp(lnF);
    }
    enum class BetaShape { allowed, firstForbiddenUnique, secondForbiddenUnique };

    class BetaSpectrum {
    public:
      BetaSpectrum() = default;
      BetaSpectrum(double endpoint_MeV, int zDaughter, BetaShape shape, int aDaughter = 0,
                   int nBins = 2000)
        : _endpoint(endpoint_MeV), _dT(endpoint_MeV / nBins), _cdf(nBins + 1, 0.0) {
        if (endpoint_MeV <= 0.0 || nBins < 10)
          throw cet::exception("RANGE") << "BetaSpectrum: endpoint must be > 0\n";
        const double W0 = 1.0 + endpoint_MeV / kElectronMass;
        for (int i = 0; i < nBins; ++i) {
          const double T = (i + 0.5) * _dT;
          const double W = 1.0 + T / kElectronMass;
          const double p = std::sqrt(W * W - 1.0);
          const double q = W0 - W;
          const double F = fermiFunction(zDaughter, aDaughter, W);
          const bool coulomb = (zDaughter != 0 && aDaughter > 0);
          double S = 1.0;
          if (shape == BetaShape::firstForbiddenUnique) {
            const double l2 = coulomb ? fermiFunction(zDaughter, aDaughter, W, 2) / F : 1.0;
            S = q * q + l2 * p * p;
          } else if (shape == BetaShape::secondForbiddenUnique) {
            const double l2 = coulomb ? fermiFunction(zDaughter, aDaughter, W, 2) / F : 1.0;
            const double l3 = coulomb ? fermiFunction(zDaughter, aDaughter, W, 3) / F : 1.0;
            S = std::pow(q, 4) + (10.0 / 3.0) * l2 * p * p * q * q + l3 * std::pow(p, 4);
          }
          _cdf[i + 1] = _cdf[i] + p * W * q * q * F * S;
        }
      }

      // Kinetic energy [MeV]
      double sample(CLHEP::RandFlat& flat) const {
        const double x = flat.fire() * _cdf.back();
        const auto it = std::upper_bound(_cdf.begin(), _cdf.end(), x);
        const size_t i = std::min<size_t>(std::max<long>(it - _cdf.begin() - 1, 0), _cdf.size() - 2);
        const double width = _cdf[i + 1] - _cdf[i];
        const double frac = (width > 0.0) ? (x - _cdf[i]) / width : 0.5;
        return (i + frac) * _dT;
      }

      double endpoint() const { return _endpoint; }

    private:
      double _endpoint = 0.0;
      double _dT = 0.0;
      std::vector<double> _cdf;
    };

    // ---------------------------------------------------------------------------------------
    // La-138
    // ---------------------------------------------------------------------------------------
    struct La138Data {
      // Branching [ENSDF]: EC 65.5(4)%, beta- 34.5(4)%. (DDEP: 65.2(6) / 34.8(6).)
      double ecFraction = 0.655;

      double gammaEC   = 1.435795;   // Ba-138 2+ -> 0+ [MeV] [ENSDF]
      double gammaBeta = 0.788742;   // Ce-138 2+ -> 0+ [MeV] [ENSDF]
      double betaEndpoint = 0.26371; // Q(beta-) 1052.46 keV [ENSDF/PT] - level 788.744 keV [MeV]
      int    zCe = 58;
      int    aCe = 138;

      // Ba atomic relaxation after capture. The EC is 2nd-forbidden unique, which enhances
      // L capture relative to an allowed transition [DDEP]: P_K 0.637(5), P_L 0.275(3),
      // P_M 0.0880(11) (remainder M+ here).
      double pK = 0.637;
      double pL = 0.275;
      double omegaK = 0.900;          // Ba K-shell fluorescence yield, 0.900(4) [DDEP]
      double bindK  = 0.037441;       // Ba K binding [MeV] [XDB]
      double bindL1 = 0.005989;       // Ba L1 binding [MeV] [XDB] (L2 is 5.624; L1 used)
      double bindM1 = 0.001293;       // Ba M1 binding [MeV] [XDB]
      // K X-rays [DDEP] relative: Ka2 31.817 (54.28), Ka1 32.194 (100), Kb1' ~36.35 (29.41),
      // Kb2' ~37.28 (7.41).
      double kAlpha = 0.032061;       // Ba K-alpha, intensity weighted [MeV]
      double kBeta  = 0.036540;       // Ba K-beta, intensity weighted [MeV]
      double kAlphaFraction = 0.807;  // K-alpha share of K X-rays: 154.28 / 191.10
    };

    class La138Decay {
    public:
      explicit La138Decay(const La138Data& d = La138Data(),
                          BetaShape shape = BetaShape::secondForbiddenUnique)
        : _d(d), _beta(d.betaEndpoint, d.zCe, shape, d.aCe) {}

      void emit(Emitter& e) const {
        if (e.flat.fire() < _d.ecFraction) {
          e.gamma(_d.gammaEC);
          // Vacancy relaxation. Whatever is not carried off by an X-ray is released as
          // Auger/Coster-Kronig electrons with ~um range, modelled as one local electron.
          const double u = e.flat.fire();
          if (u < _d.pK) {
            if (e.flat.fire() < _d.omegaK) {
              const double xray = (e.flat.fire() < _d.kAlphaFraction) ? _d.kAlpha : _d.kBeta;
              e.gamma(xray);
              e.electron(_d.bindK - xray);
            } else {
              e.electron(_d.bindK);
            }
          } else if (u < _d.pK + _d.pL) {
            e.electron(_d.bindL1);
          } else {
            e.electron(_d.bindM1);
          }
        } else {
          e.gamma(_d.gammaBeta);
          e.electron(_beta.sample(e.flat));
        }
      }

      const La138Data& data() const { return _d; }

    private:
      La138Data _d;
      BetaSpectrum _beta;
    };

    // ---------------------------------------------------------------------------------------
    // Ac-227 chain
    // ---------------------------------------------------------------------------------------
    struct AlphaLine {
      double energy;     // alpha kinetic energy [MeV]
      double intensity;  // per decay
      double level;      // daughter level fed [MeV], emitted as one gamma
    };

    struct BetaBranch {
      double intensity;             // per decay
      double level;                 // daughter level fed [MeV]
      std::vector<double> gammas;   // de-excitation gammas [MeV]
    };

    // One radioactive species: alpha lines and/or beta branches, intensities per decay of
    // THIS nuclide. Each decay picks exactly one line/branch. The tabulated intensities are
    // renormalised to sum to 1, so unlisted weak lines (and the omitted 0.28% Bi-211 b-
    // branch) are absorbed proportionally rather than producing empty decays.
    class Nuclide {
    public:
      Nuclide(const char* name, double chainFraction, std::vector<AlphaLine> alphas,
              double betaQ, int zDaughter, int aDaughter, std::vector<BetaBranch> betas)
        : _name(name), _chainFraction(chainFraction), _alphas(std::move(alphas)),
          _betas(std::move(betas)) {
        for (const auto& a : _alphas) _total += a.intensity;
        for (const auto& b : _betas) {
          _total += b.intensity;
          _betaSpectra.emplace_back(betaQ - b.level, zDaughter, BetaShape::allowed, aDaughter);
        }
        if (_total <= 0.0)
          throw cet::exception("RANGE") << "Nuclide " << name << " has no decay lines\n";
      }

      void emit(Emitter& e) const {
        double u = e.flat.fire() * _total;
        for (const auto& a : _alphas) {
          if (u < a.intensity) {
            e.alpha(a.energy);
            e.gamma(a.level);
            return;
          }
          u -= a.intensity;
        }
        for (size_t i = 0; i < _betas.size(); ++i) {
          if (u < _betas[i].intensity) {
            e.electron(_betaSpectra[i].sample(e.flat));
            for (double g : _betas[i].gammas) e.gamma(g);
            return;
          }
          u -= _betas[i].intensity;
        }
      }

      const char* name() const { return _name; }
      // Equilibrium activity of this member relative to the Ac-227 chain head.
      double chainFraction() const { return _chainFraction; }

    private:
      const char* _name;
      double _chainFraction;
      double _total = 0.0;
      std::vector<AlphaLine> _alphas;
      std::vector<BetaBranch> _betas;
      std::vector<BetaSpectrum> _betaSpectra;
    };

    // The chain, in decay order. Alpha energies, intensities and daughter levels, beta
    // branches and Q values are ENSDF (see the source list at the top of this file).
    inline std::vector<Nuclide> makeAc227Chain() {
      std::vector<Nuclide> chain;

      // Ac-227 (21.772 yr): b- 98.62% (Q 44.76 keV, first-forbidden non-unique), a 1.38%.
      // Beta branches per decay: g.s. 54, 9.3 keV level 35, 24.5 keV level 10, 37.9 keV level
      // 0.3 (per 100 beta decays) x 0.9862. The two main alpha lines are scaled to the full
      // 1.38% alpha branch.
      chain.emplace_back("Ac227", 1.0,
        std::vector<AlphaLine>{{4.95326, 0.00754, 0.0}, {4.9407, 0.00626, 0.01289}},
        0.04476, 90, 227, std::vector<BetaBranch>{
          {0.5363, 0.0, {}},
          {0.3476, 0.0093, {0.0093}},
          {0.0993, 0.0245, {0.0245}},
          {0.0030, 0.0379, {0.0379}}});

      // Th-227 (18.697 d), fed by the 98.62% b- branch of Ac-227.
      chain.emplace_back("Th227", 0.9862,
        std::vector<AlphaLine>{
          {6.03801, 0.242, 0.0},       {6.0088, 0.029, 0.029858},  {5.97772, 0.235, 0.061424},
          {5.9597, 0.030, 0.079708},   {5.916, 0.00775, 0.123793}, {5.9099, 0.00174, 0.130141},
          {5.8666, 0.0242, 0.174569},  {5.8075, 0.0127, 0.234858}, {5.7955, 0.00311, 0.24739},
          {5.7623, 0.00228, 0.280182}, {5.75687, 0.204, 0.286087}, {5.7132, 0.0489, 0.329856},
          {5.7088, 0.083, 0.334370},   {5.7008, 0.0363, 0.342590}, {5.693, 0.015, 0.35053},
          {5.668, 0.0206, 0.376296},   {5.6133, 0.00216, 0.43224}, {5.6006, 0.0017, 0.445071},
          {5.5859, 0.00176, 0.45993}},
        0.0, 0, 0, std::vector<BetaBranch>{});

      // Ra-223 (11.43 d). Both the Th-227 and Fr-223 branches feed it: fraction 1.
      chain.emplace_back("Ra223", 1.0,
        std::vector<AlphaLine>{
          {5.8713, 0.010, 0.0},        {5.8575, 0.0031, 0.014400}, {5.747, 0.089, 0.126726},
          {5.71623, 0.512, 0.158631},  {5.60673, 0.250, 0.269475}, {5.5398, 0.089, 0.338276},
          {5.528, 0.0014, 0.342792},   {5.5016, 0.010, 0.37619},   {5.4336, 0.018, 0.445037},
          {5.424, 0.0039, 0.44682},    {5.3656, 0.0013, 0.51442},  {5.2873, 0.002, 0.59408},
          {5.2828, 0.00135, 0.598721}},
        0.0, 0, 0, std::vector<BetaBranch>{});

      // Rn-219 (3.96 s). Trapped in the crystal lattice.
      chain.emplace_back("Rn219", 1.0,
        std::vector<AlphaLine>{
          {6.8191, 0.794, 0.0}, {6.5526, 0.129, 0.271228}, {6.530, 0.0011, 0.29356},
          {6.425, 0.075, 0.401812}},
        0.0, 0, 0, std::vector<BetaBranch>{});

      // Po-215 (1.781 ms).
      chain.emplace_back("Po215", 1.0,
        std::vector<AlphaLine>{{7.3861, 1.0, 0.0}},
        0.0, 0, 0, std::vector<BetaBranch>{});

      // Pb-211 (36.1 min): b- Q 1366.1 keV to Bi-211 (first-forbidden non-unique).
      //   g.s. 91.32%, 404.87 keV level 1.63%, 831.96 keV level 6.28%, 1109.49 keV level 0.83%.
      //   The 832 keV level decays directly (832.01 keV, 3.62% incl. conversion) or via
      //   427.09 + 404.85 keV (2.07%): split 0.636 / 0.364. The 1109 keV level is emitted as
      //   one gamma (its cascade is not resolved).
      chain.emplace_back("Pb211", 1.0,
        std::vector<AlphaLine>{},
        1.3661, 83, 211, std::vector<BetaBranch>{
          {0.9132, 0.0, {}},
          {0.0163, 0.404866, {0.404853}},
          {0.0399, 0.831960, {0.83201}},
          {0.0229, 0.831960, {0.427088, 0.404853}},
          {0.0083, 1.109485, {1.109485}}});

      // Bi-211 (2.14 min): a 99.724%. The 0.276% b- branch to Po-211 is omitted.
      chain.emplace_back("Bi211", 1.0,
        std::vector<AlphaLine>{{6.6229, 0.83539, 0.0}, {6.2782, 0.16185, 0.35106}},
        0.0, 0, 0, std::vector<BetaBranch>{});

      // Tl-207 (4.77 min): b- Q 1417.5 keV to Pb-207; 0.271% via the 897.76 keV level.
      chain.emplace_back("Tl207", 0.99724,
        std::vector<AlphaLine>{},
        1.4175, 82, 207, std::vector<BetaBranch>{
          {0.99729, 0.0, {}},
          {0.00271, 0.89776, {0.89776}}});

      return chain;
    }

  } // namespace LaBrIntrinsic
} // namespace mu2e
#endif
