// A configurable primary generator for the modular STM pipeline (LaBrPipeline*.fcl).
// It replaces the separate PhotonGun / PoissonPhotonGun (probably unnecessary) set-ups with three independent
// choices, each set in its own fhicl table:
//
//   mode :
//     "gun"      one particle type (pdgId) at one kinetic energy E (based off of Pawel's photon gun)
//     "spectrum" one particle type, kinetic energy drawn from `spectrum`: (essentially teh photongun but with different energies drawn from a spectrum)
//                  type "lines"     : lineEnergies + lineIntensities (relative)
//                  type "histogram" : binEdges (N+1) + binContents (N), flat within a bin
//                  type "rootHist"  : a TH1 read from fileName/histName (x axis in MeV)
//     "dataset"  particles are resampled from a TTree written by an upstream simulation
//                (STMParticleDump writes one from any VD StepPointMCs). Each entry carries
//                pdgId, position, time, and momentum. Neutrons, photons, electrons etc come
//                with their real upstream spectra, directions, and correlations. Optional
//                pdgIds / EMin / EMax filters and an optional `transform` that maps the
//                upstream plane onto the detector aperture (same idea as
//                ShiftVirtualDetectorStepPointMCs: translate + radial scale, |p| preserved).
//                sampleUnit "event" keeps all particles of one upstream event together.
//
//   emission : how many per art event 
//     poisson : true  -> N = Poisson(mean) ; mean = nPerEvent or rateHz * windowTime (depending on which input you choose)
//     poisson : false -> N = nPerEvent exactly (should just be the normal particle gun which sends photons exactly at the same point in events)
//     timing  : "uniform" (U[0, windowTime); for a dataset upstream event, its first particle
//               is placed uniformly and the rest keep their relative times), "fixed" (= time),
//               or "dataset" (upstream time folded into the window with fmod; dataset only)
//
//   beam (gun/spectrum only) : start point, optional uniform disk radius, direction and an
//     optional uniform cone around it.
//
// Every particle is tagged GenId::particleGun, so the LaBr truth modules classify it as the
// external source (STMMCPulse::photonGun), separate from LaBrIntrinsicSource's La-138/Ac-227.
// Original author: Sena Clarke

// stdlib includes
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

namespace mu2e {
  namespace {
    // A file given as a path that exists (absolute, or relative to the run directory) is
    // used as-is; otherwise it is looked up on MU2E_SEARCH_PATH like other Offline data.
    std::string resolvePath(const std::string& name) {
      if (std::ifstream(name).good())
        return name;
      ConfigFileLookupPolicy findFile;
      return findFile(name);
    }
  }

  class STMParticleSource : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    struct SpectrumConfig {
      fhicl::Atom<std::string> type{Name("type"), Comment("\"lines\", \"histogram\" or \"rootHist\"")};
      fhicl::OptionalSequence<double> lineEnergies{Name("lineEnergies"), Comment("lines: kinetic energies [MeV]")};
      fhicl::OptionalSequence<double> lineIntensities{Name("lineIntensities"), Comment("lines: relative intensities (same length). Default all 1.")};
      fhicl::OptionalSequence<double> binEdges{Name("binEdges"), Comment("histogram: N+1 bin edges [MeV]")};
      fhicl::OptionalSequence<double> binContents{Name("binContents"), Comment("histogram: N relative bin contents (counts per bin, not per MeV)")};
      fhicl::OptionalAtom<std::string> fileName{Name("fileName"), Comment("rootHist: ROOT file (a path, or a name on MU2E_SEARCH_PATH)")};
      fhicl::OptionalAtom<std::string> histName{Name("histName"), Comment("rootHist: TH1 name inside fileName, x axis in MeV")};
    };

    struct BeamConfig {
      fhicl::Atom<double> x{Name("x"), Comment("Start x [mm]")};
      fhicl::Atom<double> y{Name("y"), Comment("Start y [mm]")};
      fhicl::Atom<double> z{Name("z"), Comment("Start z [mm]")};
      fhicl::Atom<double> radius{Name("radius"), Comment("Uniform disk radius around (x,y), transverse to z [mm]. 0 = point."), 0.0};
      fhicl::Sequence<double, 3> direction{Name("direction"), Comment("Beam direction (normalised internally). Default +z."), std::array<double, 3>{0.0, 0.0, 1.0}};
      fhicl::Atom<double> coneHalfAngle{Name("coneHalfAngle"), Comment("Uniform-in-solid-angle cone half-angle around direction [deg]. 0 = pencil, 180 = isotropic."), 0.0};
    };

    struct TransformConfig {
      fhicl::Sequence<double, 3> inputCentre{Name("inputCentre"), Comment("Centre of the upstream plane, e.g. VD101 [mm]")};
      fhicl::Sequence<double, 3> outputCentre{Name("outputCentre"), Comment("Where that centre is moved to, e.g. the LaBr aperture [mm]; z is also the new plane")};
      fhicl::Atom<double> inputRadius{Name("inputRadius"), Comment("Radius of the upstream plane [mm]; particles outside it are dropped. <=0 keeps all."), 0.0};
      fhicl::Atom<double> radialScale{Name("radialScale"), Comment("Transverse position and momentum scale (e.g. 3.989/200 for VD101 -> SSC aperture); |p| is preserved."), 1.0};
    };

    struct DatasetConfig {
      fhicl::Atom<std::string> fileName{Name("fileName"), Comment("ROOT file with the particle TTree (a path, or a name on MU2E_SEARCH_PATH)")};
      fhicl::Atom<std::string> treeName{Name("treeName"), Comment("TTree name, default the STMParticleDump one"), "DumpParticles/particles"};
      fhicl::Atom<std::string> sampleUnit{Name("sampleUnit"), Comment("\"particle\": resample single entries; \"event\": resample whole upstream events (needs the eventId branch)"), "particle"};
      fhicl::OptionalSequence<int> pdgIds{Name("pdgIds"), Comment("Keep only these PDG codes, e.g. [2112] for neutrons. Default: all.")};
      fhicl::Atom<double> EMin{Name("EMin"), Comment("Minimum kinetic energy kept [MeV]"), 0.0};
      fhicl::Atom<double> EMax{Name("EMax"), Comment("Maximum kinetic energy kept [MeV], <=0 = no limit"), 0.0};
      fhicl::OptionalTable<TransformConfig> transform{Name("transform"), Comment("Map the upstream plane onto the aperture")};
    };

    struct EmissionConfig {
      fhicl::Atom<bool> poisson{Name("poisson"), Comment("true: N = Poisson(mean) per event; false: exactly nPerEvent")};
      fhicl::OptionalAtom<double> nPerEvent{Name("nPerEvent"), Comment("Poisson mean, or the fixed count (rounded) when poisson is false")};
      fhicl::OptionalAtom<double> rateHz{Name("rateHz"), Comment("Poisson only: rate [Hz]; mean = rateHz * windowTime. Overrides nPerEvent.")};
      fhicl::Atom<double> windowTime{Name("windowTime"), Comment("Event window [ns]; match the digitizer micropulse (1695)"), 1695.0};
      fhicl::Atom<std::string> timing{Name("timing"), Comment("\"uniform\", \"fixed\" or \"dataset\""), "uniform"};
      fhicl::Atom<double> time{Name("time"), Comment("Emission time for timing \"fixed\" [ns]"), 0.0};
    };

    struct Config {
      fhicl::Atom<std::string> mode{Name("mode"), Comment("\"gun\", \"spectrum\" or \"dataset\"")};
      fhicl::Atom<int> pdgId{Name("pdgId"), Comment("gun/spectrum particle type. Default 22 (photon)."), 22};
      fhicl::OptionalAtom<double> E{Name("E"), Comment("gun: kinetic energy [MeV]")};
      fhicl::OptionalTable<SpectrumConfig> spectrum{Name("spectrum"), Comment("spectrum: energy distribution")};
      fhicl::OptionalTable<BeamConfig> beam{Name("beam"), Comment("gun/spectrum: start position and direction")};
      fhicl::OptionalTable<DatasetConfig> dataset{Name("dataset"), Comment("dataset: upstream particle list")};
      fhicl::Table<EmissionConfig> emission{Name("emission"), Comment("Multiplicity and timing per event")};
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity"), 0};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit STMParticleSource(const Parameters& conf);
    void produce(art::Event& event) override;
    void endJob() override;

  private:
    struct Particle {
      int pdgId;
      CLHEP::Hep3Vector pos;
      CLHEP::Hep3Vector mom; // MeV/c
      double time;           // ns
    };

    void buildSpectrum(const SpectrumConfig& sc);
    void loadDataset(const DatasetConfig& dc);
    double sampleEnergy();
    CLHEP::Hep3Vector sampleDirection();
    CLHEP::Hep3Vector samplePosition();
    double sampleTime(double datasetTime);
    void emit(GenParticleCollection& out, int pdgId, const CLHEP::Hep3Vector& pos,
              const CLHEP::Hep3Vector& mom, double t);

    enum class Mode { gun, spectrum, dataset } _mode;
    enum class Timing { uniform, fixed, dataset } _timing;
    int _pdgId = 22;
    int _verbosity = 0;

    // Multiplicity
    bool _poisson = true;
    double _mean = 0.0;
    int _nFixed = 0;
    double _window = 1695.0;
    double _fixedTime = 0.0;

    // Energy: mono (_E) or a tabulated CDF over _cdf, either discrete lines or flat bins
    double _E = 0.0;
    bool _discrete = true;
    std::vector<double> _values; // line energies, or bin edges (N+1)
    std::vector<double> _cdf;    // cumulative weights, normalised to 1 (N entries)

    // Beam
    CLHEP::Hep3Vector _start, _axis, _u, _v;
    double _radius = 0.0, _cosMax = 1.0;

    // Dataset: particles, grouped into upstream events (_groups[i] = [begin, end))
    std::vector<Particle> _particles;
    std::vector<std::pair<size_t, size_t>> _groups;

    GlobalConstantsHandle<ParticleDataList> _pdt;
    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandFlat _flat;
    CLHEP::RandPoissonQ _randPoisson;
    unsigned long _nEvents = 0, _nEmitted = 0;
  };

  STMParticleSource::STMParticleSource(const Parameters& conf) :
    art::EDProducer(conf),
    _pdgId(conf().pdgId()),
    _verbosity(conf().verbosityLevel()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _flat(_engine),
    _randPoisson(_engine) {
      produces<GenParticleCollection>();

      // ---- mode ----
      const std::string mode = conf().mode();
      if (mode == "gun") {
        _mode = Mode::gun;
        if (!conf().E(_E) || _E <= 0.0)
          throw cet::exception("CONFIG") << "mode \"gun\" needs E > 0 [MeV]\n";
      } else if (mode == "spectrum") {
        _mode = Mode::spectrum;
        SpectrumConfig sc;
        if (!conf().spectrum(sc))
          throw cet::exception("CONFIG") << "mode \"spectrum\" needs a spectrum table\n";
        buildSpectrum(sc);
      } else if (mode == "dataset") {
        _mode = Mode::dataset;
        DatasetConfig dc;
        if (!conf().dataset(dc))
          throw cet::exception("CONFIG") << "mode \"dataset\" needs a dataset table\n";
        loadDataset(dc);
      } else
        throw cet::exception("CONFIG") << "unknown mode \"" << mode << "\" (gun, spectrum, dataset)\n";

      // ---- beam (gun/spectrum) ----
      if (_mode != Mode::dataset) {
        BeamConfig bc;
        if (!conf().beam(bc))
          throw cet::exception("CONFIG") << "modes gun/spectrum need a beam table\n";
        _start.set(bc.x(), bc.y(), bc.z());
        _radius = bc.radius();
        const auto d = bc.direction();
        _axis = CLHEP::Hep3Vector(d[0], d[1], d[2]);
        if (_axis.mag() <= 0.0)
          throw cet::exception("CONFIG") << "beam.direction must be non-zero\n";
        _axis = _axis.unit();
        _u = _axis.orthogonal().unit();
        _v = _axis.cross(_u);
        const double half = bc.coneHalfAngle();
        if (half < 0.0 || half > 180.0)
          throw cet::exception("CONFIG") << "beam.coneHalfAngle must be in [0, 180] deg\n";
        _cosMax = std::cos(half * M_PI / 180.0);
      }

      // ---- emission ----
      const auto& ec = conf().emission();
      _poisson = ec.poisson();
      _window = ec.windowTime();
      _fixedTime = ec.time();
      if (_window <= 0.0)
        throw cet::exception("CONFIG") << "emission.windowTime must be > 0\n";
      double n = 0.0, rate = 0.0;
      const bool haveN = ec.nPerEvent(n);
      if (ec.rateHz(rate)) {
        if (!_poisson)
          throw cet::exception("CONFIG") << "emission.rateHz needs poisson : true\n";
        n = rate * _window * 1e-9;
      } else if (!haveN)
        throw cet::exception("CONFIG") << "emission needs nPerEvent or rateHz\n";
      if (n < 0.0)
        throw cet::exception("CONFIG") << "emission rate must be >= 0\n";
      _mean = n;
      _nFixed = static_cast<int>(std::lround(n));

      const std::string timing = ec.timing();
      if (timing == "uniform") _timing = Timing::uniform;
      else if (timing == "fixed") _timing = Timing::fixed;
      else if (timing == "dataset") {
        if (_mode != Mode::dataset)
          throw cet::exception("CONFIG") << "timing \"dataset\" only works with mode \"dataset\"\n";
        _timing = Timing::dataset;
      } else
        throw cet::exception("CONFIG") << "unknown timing \"" << timing << "\" (uniform, fixed, dataset)\n";

      if (_verbosity > 0)
        mf::LogInfo("STMParticleSource") << "mode " << mode << ", "
          << (_poisson ? "Poisson mean " : "fixed ") << (_poisson ? _mean : _nFixed)
          << " per " << _window << " ns event, timing " << timing
          << (_mode == Mode::dataset ? ", " + std::to_string(_groups.size()) + " dataset units" : "");
    };

  void STMParticleSource::buildSpectrum(const SpectrumConfig& sc) {
    std::vector<double> weights;
    const std::string type = sc.type();
    if (type == "lines") {
      _discrete = true;
      if (!sc.lineEnergies(_values) || _values.empty())
        throw cet::exception("CONFIG") << "spectrum \"lines\" needs lineEnergies\n";
      if (!sc.lineIntensities(weights))
        weights.assign(_values.size(), 1.0);
      if (weights.size() != _values.size())
        throw cet::exception("CONFIG") << "lineIntensities must match lineEnergies in length\n";
    } else if (type == "histogram") {
      _discrete = false;
      if (!sc.binEdges(_values) || !sc.binContents(weights) || _values.size() != weights.size() + 1)
        throw cet::exception("CONFIG") << "spectrum \"histogram\" needs N+1 binEdges and N binContents\n";
    } else if (type == "rootHist") {
      _discrete = false;
      std::string fileName, histName;
      if (!sc.fileName(fileName) || !sc.histName(histName))
        throw cet::exception("CONFIG") << "spectrum \"rootHist\" needs fileName and histName\n";
      std::unique_ptr<TFile> f(TFile::Open(resolvePath(fileName).c_str()));
      TH1* h = f ? dynamic_cast<TH1*>(f->Get(histName.c_str())) : nullptr;
      if (!h)
        throw cet::exception("CONFIG") << "cannot read TH1 " << histName << " from " << fileName << "\n";
      const auto* ax = h->GetXaxis();
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
        _values.push_back(ax->GetBinLowEdge(i));
        weights.push_back(h->GetBinContent(i));
      }
      _values.push_back(ax->GetBinUpEdge(h->GetNbinsX()));
    } else
      throw cet::exception("CONFIG") << "unknown spectrum type \"" << type << "\" (lines, histogram, rootHist)\n";

    for (size_t i = 0; i + 1 < _values.size(); ++i)
      if (!_discrete && _values[i + 1] <= _values[i])
        throw cet::exception("CONFIG") << "spectrum binEdges must increase\n";
    double acc = 0.0;
    for (double w : weights) {
      if (w < 0.0)
        throw cet::exception("CONFIG") << "spectrum weights must be >= 0\n";
      acc += w;
      _cdf.push_back(acc);
    }
    if (acc <= 0.0)
      throw cet::exception("CONFIG") << "spectrum total weight must be > 0\n";
    for (double& c : _cdf)
      c /= acc;
    if (!_discrete && _values.front() < 0.0)
      throw cet::exception("CONFIG") << "spectrum energies must be >= 0\n";
  };

  void STMParticleSource::loadDataset(const DatasetConfig& dc) {
    const std::string path = resolvePath(dc.fileName());
    std::unique_ptr<TFile> f(TFile::Open(path.c_str()));
    TTree* t = f ? dynamic_cast<TTree*>(f->Get(dc.treeName().c_str())) : nullptr;
    if (!t)
      throw cet::exception("CONFIG") << "cannot read TTree " << dc.treeName() << " from " << path << "\n";

    int pdgId = 0;
    Long64_t eventId = 0;
    double x = 0, y = 0, z = 0, px = 0, py = 0, pz = 0, time = 0;
    t->SetBranchAddress("pdgId", &pdgId);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("z", &z);
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("time", &time);
    const bool byEvent = dc.sampleUnit() == "event";
    if (!byEvent && dc.sampleUnit() != "particle")
      throw cet::exception("CONFIG") << "dataset.sampleUnit must be \"particle\" or \"event\"\n";
    if (byEvent) {
      if (!t->GetBranch("eventId"))
        throw cet::exception("CONFIG") << "dataset.sampleUnit \"event\" needs an eventId branch\n";
      t->SetBranchAddress("eventId", &eventId);
    }

    std::vector<int> keep;
    dc.pdgIds(keep);
    const double EMin = dc.EMin(), EMax = dc.EMax();

    TransformConfig tc;
    const bool doTransform = dc.transform(tc);
    CLHEP::Hep3Vector inC, outC;
    double inR = 0.0, scale = 1.0;
    if (doTransform) {
      const auto a = tc.inputCentre(), b = tc.outputCentre();
      inC.set(a[0], a[1], a[2]);
      outC.set(b[0], b[1], b[2]);
      inR = tc.inputRadius();
      scale = tc.radialScale();
      if (scale <= 0.0)
        throw cet::exception("CONFIG") << "transform.radialScale must be > 0\n";
    }

    Long64_t lastEvent = -1;
    unsigned long nRead = t->GetEntries(), nDropped = 0, nUnknown = 0;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!keep.empty() && std::find(keep.begin(), keep.end(), pdgId) == keep.end()) { ++nDropped; continue; }
      CLHEP::Hep3Vector pos(x, y, z), mom(px, py, pz);
      double m = 0.0;
      try {
        m = _pdt->particle(pdgId).mass();
      } catch (const cet::exception&) { // not in the particle table (e.g. a recoil nucleus)
        ++nUnknown;
        continue;
      }
      const double Ek = std::sqrt(mom.mag2() + m * m) - m;
      if (Ek < EMin || (EMax > 0.0 && Ek > EMax)) { ++nDropped; continue; }
      if (doTransform) {
        CLHEP::Hep3Vector d = pos - inC;
        d.setZ(0.0);
        if (inR > 0.0 && d.perp() > inR) { ++nDropped; continue; }
        pos = outC + d * scale;
        // Scale the transverse momentum, keep |p| and the sign of pz (focuses the upstream
        // plane's angular spread onto the smaller aperture, as ShiftVirtualDetectorStepPointMCs).
        const double p = mom.mag();
        const double pxN = mom.x() * scale, pyN = mom.y() * scale;
        const double pzN = std::copysign(std::sqrt(std::max(0.0, p * p - pxN * pxN - pyN * pyN)), mom.z());
        mom.set(pxN, pyN, pzN);
      }
      if (!byEvent || eventId != lastEvent || _groups.empty())
        _groups.emplace_back(_particles.size(), _particles.size());
      _particles.push_back(Particle{pdgId, pos, mom, time});
      _groups.back().second = _particles.size();
      lastEvent = eventId;
    }
    if (_particles.empty())
      throw cet::exception("CONFIG") << "no dataset particles pass the selection in " << path << "\n";
    mf::LogInfo("STMParticleSource") << "dataset " << path << ": " << nRead << " entries, "
      << _particles.size() << " kept (" << nDropped << " filtered, " << nUnknown
      << " with a PDG code missing from the particle table), " << _groups.size()
      << (byEvent ? " upstream events" : " sampling units");
  };

  double STMParticleSource::sampleEnergy() {
    if (_mode == Mode::gun)
      return _E;
    const double u = _flat.fire();
    const size_t i = std::min<size_t>(std::lower_bound(_cdf.begin(), _cdf.end(), u) - _cdf.begin(), _cdf.size() - 1);
    if (_discrete)
      return _values[i];
    return _values[i] + _flat.fire() * (_values[i + 1] - _values[i]);
  };

  CLHEP::Hep3Vector STMParticleSource::sampleDirection() {
    if (_cosMax >= 1.0)
      return _axis;
    const double cosT = 1.0 - _flat.fire() * (1.0 - _cosMax);
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const double phi = _flat.fire(2.0 * M_PI);
    return (_axis * cosT + (_u * std::cos(phi) + _v * std::sin(phi)) * sinT).unit();
  };

  CLHEP::Hep3Vector STMParticleSource::samplePosition() {
    if (_radius <= 0.0)
      return _start;
    const double r = _radius * std::sqrt(_flat.fire());
    const double phi = _flat.fire(2.0 * M_PI);
    return _start + CLHEP::Hep3Vector(r * std::cos(phi), r * std::sin(phi), 0.0);
  };

  double STMParticleSource::sampleTime(double datasetTime) {
    switch (_timing) {
    case Timing::fixed: return _fixedTime;
    case Timing::dataset: {
      const double t = std::fmod(datasetTime, _window);
      return t < 0.0 ? t + _window : t;
    }
    default: return _flat.fire(_window);
    }
  };

  void STMParticleSource::emit(GenParticleCollection& out, int pdgId, const CLHEP::Hep3Vector& pos,
                               const CLHEP::Hep3Vector& mom, double t) {
    const double m = _pdt->particle(pdgId).mass();
    const CLHEP::HepLorentzVector p4(mom, std::sqrt(mom.mag2() + m * m));
    out.push_back(GenParticle(static_cast<PDGCode::type>(pdgId), GenId::particleGun, pos, p4, t));
    ++_nEmitted;
  };

  void STMParticleSource::produce(art::Event& event) {
    auto out = std::make_unique<GenParticleCollection>();
    ++_nEvents;
    const long n = _poisson ? (_mean > 0.0 ? _randPoisson.fire(_mean) : 0) : _nFixed;
    for (long i = 0; i < n; ++i) {
      if (_mode == Mode::dataset) {
        // One sampling unit = one particle, or one whole upstream event. Particles of the
        // same event share a random time shift, measured from the earliest of them, so
        // their relative timing is kept.
        const auto& g = _groups[static_cast<size_t>(_flat.fire() * _groups.size()) % _groups.size()];
        const double shift = _flat.fire(_window);
        double t0 = _particles[g.first].time;
        for (size_t k = g.first; k < g.second; ++k)
          t0 = std::min(t0, _particles[k].time);
        for (size_t k = g.first; k < g.second; ++k) {
          const Particle& p = _particles[k];
          // uniform: the unit starts at a random time and keeps its internal time structure,
          // so late particles can land in the following micropulse (the digitizer buffers
          // microspillBufferLengthCount windows ahead).
          const double t = _timing == Timing::uniform ? shift + (p.time - t0)
                                                      : sampleTime(p.time);
          emit(*out, p.pdgId, p.pos, p.mom, t);
        }
      } else {
        const double Ek = sampleEnergy();
        const double m = _pdt->particle(_pdgId).mass();
        const double p = std::sqrt(Ek * (Ek + 2.0 * m));
        emit(*out, _pdgId, samplePosition(), sampleDirection() * p, sampleTime(0.0));
      }
    }
    event.put(std::move(out));
  };

  void STMParticleSource::endJob() {
    mf::LogInfo("STMParticleSource") << "emitted " << _nEmitted << " particles in " << _nEvents
      << " events (" << (_nEvents ? double(_nEmitted) / _nEvents : 0.0) << " per event)";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::STMParticleSource)
