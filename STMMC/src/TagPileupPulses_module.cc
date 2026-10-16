// Stage 2 of the STM pileup-truth chain.
//
// Reads the per-micropulse STMMCPulse collections from LaBrTruthPulses, buffers
// nMerge micropulses onto one absolute timeline (mirroring how
// ConcatenateDigitizedWaveforms merges nMerge waveforms into one), and tags each
// pulse that piles up with a neighbour it cannot be resolved from.
//
// Pileup window (energy-dependent, analytic). A pulse of energy E has peak
// amplitude peakADCPerKeV*E. Its scintillation tail decays as exp(-dt/tau2), so
// it stays above the resolving threshold for
//     w_after(E) = tau2 * ln( peakADCPerKeV * E / thresholdADC )
// (floored at minWindow, and 0 if the pulse never reaches threshold). A later
// pulse landing within w_after of an earlier one sits on its still-elevated,
// sloped tail and cannot be cleanly separated; a pulse landing within the small
// fixed rise-side window wBefore of a later pulse is caught on its rising edge.
// Because w_after grows only logarithmically with E, the reach past a micropulse
// boundary is at most a few hundred ns even for very large pulses, so buffering a
// batch of micropulses captures essentially all cross-boundary pileup.
//
// Working on the sorted, merged list makes cross-micropulse pileup automatic: the
// last pulse of one event and the first of the next are simply adjacent entries.
// Only pileup straddling a batch boundary (every nMerge events) is missed, the
// same limitation the waveform concatenation already has.
//
// Original author: Sena Clarke

// stdlib includes
#include <algorithm>
#include <cmath>
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/STMMCPulse.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class TagPileupPulses : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> STMMCPulsesTag{Name("STMMCPulsesTag"), Comment("STMMCPulseCollection from LaBrTruthPulses")};
      fhicl::Atom<int> nMerge{Name("nMerge"), Comment("Number of micropulses to buffer onto one timeline before tagging. Should match the waveform concatenation's nMerge. Use the total event count for a single batch."), 1};
      fhicl::Atom<double> pulseTau2{Name("pulseTau2"), Comment("Pulse decay time constant [ns], sets how fast the tail falls to threshold. Default 28.0007 (LaBr)."), 28.0007};
      fhicl::Atom<double> peakADCPerKeV{Name("peakADCPerKeV"), Comment("Peak pulse amplitude per unit energy [ADC/keV] = 1/peakKeVPerADC. Default 1.7544 (=1/0.57)."), 1.7543859649};
      fhicl::Atom<double> thresholdADC{Name("thresholdADC"), Comment("Resolving threshold [ADC]; a pulse tail below this no longer shadows. Should match the ZS/pulse-finder threshold. Default 5."), 5.0};
      fhicl::Atom<double> wBefore{Name("wBefore"), Comment("Fixed rise-side pileup window [ns]: a pulse this close before another is caught on its rising edge. Default 20."), 20.0};
      fhicl::Atom<double> minWindow{Name("minWindow"), Comment("Minimum forward window [ns] applied even to sub-threshold pulses, so coincident small pulses still count. Default 20."), 20.0};
      fhicl::OptionalAtom<bool> makeTTree{Name("makeTTree"), Comment("Write a per-pulse diagnostic TTree")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit TagPileupPulses(const Parameters& conf);
    void produce(art::Event& event) override;
    void endJob() override;
  private:
    void tagAndEmit(art::Event& event);
    void tagBuffer();   // sort the buffer and set the pileup tags
    void fillTree();    // write the buffer to the TTree and count pulses
    double wAfter(double energy_keV) const;

    art::ProductToken<STMMCPulseCollection> STMMCPulsesToken;
    int nMerge = 1;
    double pulseTau2 = 28.0007;
    double peakADCPerKeV = 1.7543859649;
    double thresholdADC = 5.0;
    double wBefore = 20.0;
    double minWindow = 20.0;
    bool makeTTree = false;

    // Buffer of pulses accumulated across the current batch of micropulses.
    std::vector<STMMCPulse> _buffer;
    uint _inputEvents = 0;

    // Diagnostics
    uint64_t n_pulses = 0, n_pileup = 0;

    // TTree
    TTree* ttree = nullptr;
    double t_time = 0, t_energy = 0, t_dtNearest = 0, t_window = 0;
    int t_source = 0, t_nOverlap = 0, t_isPileup = 0;
    uint t_eventId = 0, t_primaryId = 0;
  };

  TagPileupPulses::TagPileupPulses(const Parameters& conf) :
    art::EDProducer(conf),
    STMMCPulsesToken(consumes<STMMCPulseCollection>(conf().STMMCPulsesTag())),
    nMerge(conf().nMerge()),
    pulseTau2(conf().pulseTau2()),
    peakADCPerKeV(conf().peakADCPerKeV()),
    thresholdADC(conf().thresholdADC()),
    wBefore(conf().wBefore()),
    minWindow(conf().minWindow()) {
      produces<STMMCPulseCollection>();
      if (nMerge < 1)
        throw cet::exception("RANGE") << "nMerge must be >= 1\n";
      if (pulseTau2 <= 0.0 || thresholdADC <= 0.0)
        throw cet::exception("RANGE") << "pulseTau2 and thresholdADC must be > 0\n";
      makeTTree = conf().makeTTree() ? *(conf().makeTTree()) : false;
      if (makeTTree) {
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>("pulses", "STM truth pulses with pileup tag");
        ttree->Branch("time", &t_time, "time/D");            // absolute time [ns]
        ttree->Branch("energy", &t_energy, "energy/D");      // deposited energy [keV]
        ttree->Branch("source", &t_source, "source/I");      // 0 external source, 1 La-138, 2 Ac-227
        ttree->Branch("eventId", &t_eventId, "eventId/i");
        ttree->Branch("primaryId", &t_primaryId, "primaryId/i");
        ttree->Branch("isPileup", &t_isPileup, "isPileup/I");
        ttree->Branch("nOverlap", &t_nOverlap, "nOverlap/I");
        ttree->Branch("dtNearest", &t_dtNearest, "dtNearest/D"); // [ns] to nearest overlap
        ttree->Branch("window", &t_window, "window/D");          // w_after(E) [ns]
      }
    }

  double TagPileupPulses::wAfter(double energy_keV) const {
    // Time for the exp(-dt/tau2) tail of a pulse of this energy to fall back to
    // the resolving threshold. Grows only logarithmically with energy.
    double peakADC = peakADCPerKeV * energy_keV;
    if (peakADC <= thresholdADC)
      return minWindow;
    return std::max(minWindow, pulseTau2 * std::log(peakADC / thresholdADC));
  }

  void TagPileupPulses::produce(art::Event& event) {
    ++_inputEvents;
    auto const& pulses = event.getProduct(STMMCPulsesToken);
    _buffer.insert(_buffer.end(), pulses.begin(), pulses.end());

    // Emit on the batch boundary, otherwise emit an empty collection (same
    // cadence as ConcatenateDigitizedWaveforms).
    if ((_inputEvents % nMerge) == 0) {
      tagAndEmit(event);
    } else {
      event.put(std::make_unique<STMMCPulseCollection>());
    }
  }

  void TagPileupPulses::tagAndEmit(art::Event& event) {
    tagBuffer();
    fillTree();
    auto out = std::make_unique<STMMCPulseCollection>(_buffer);
    _buffer.clear();
    event.put(std::move(out));
  }

  void TagPileupPulses::tagBuffer() {
    std::sort(_buffer.begin(), _buffer.end(),
              [](const STMMCPulse& a, const STMMCPulse& b) { return a.time() < b.time(); });

    // Record each pulse's own forward window and reset tags.
    const size_t N = _buffer.size();
    std::vector<double> reach(N);
    for (size_t i = 0; i < N; ++i) {
      double w = wAfter(_buffer[i].energy());
      _buffer[i].setWindow(w);
      _buffer[i].setPileup(false);
      _buffer[i].setNOverlap(0);
      _buffer[i].setDtNearest(0.0);
      reach[i] = std::max(w, wBefore); // furthest a later pulse can be and still overlap i
    }

    // Sorted list: for each i only scan forward until the gap exceeds its reach.
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = i + 1; j < N; ++j) {
        double dt = _buffer[j].time() - _buffer[i].time();
        if (dt >= reach[i])
          break; // sorted: no later j can be closer
        // Overlap if j lands on i's decaying tail (w_after of i) or i lands on
        // j's rising edge (fixed wBefore).
        if (dt < _buffer[i].window() || dt < wBefore) {
          for (size_t k : {i, j}) {
            _buffer[k].setPileup(true);
            _buffer[k].setNOverlap(_buffer[k].nOverlap() + 1);
            double prev = _buffer[k].dtNearest();
            if (prev == 0.0 || dt < prev)
              _buffer[k].setDtNearest(dt);
          }
        }
      }
    }

  }

  void TagPileupPulses::fillTree() {
    n_pulses += _buffer.size();
    for (const auto& p : _buffer) {
      if (p.isPileup())
        ++n_pileup;
      if (makeTTree) {
        t_time = p.time();
        t_energy = p.energy();
        t_source = p.source();
        t_eventId = p.eventId();
        t_primaryId = p.primaryId();
        t_isPileup = p.isPileup() ? 1 : 0;
        t_nOverlap = p.nOverlap();
        t_dtNearest = p.dtNearest();
        t_window = p.window();
        ttree->Fill();
      }
    }
  }

  void TagPileupPulses::endJob() {
    // A job that ends mid-batch (events not a multiple of nMerge) still has pulses in the
    // buffer. There is no event left to put them in, but they are tagged and written to the
    // TTree so the tree always covers the whole run. Pileup across the unseen next batch
    // boundary cannot be checked, exactly as at any batch boundary.
    if (!_buffer.empty()) {
      mf::LogInfo("TagPileupPulses")
        << "end of job mid-batch: tagging and writing " << _buffer.size()
        << " buffered pulses to the TTree (not to the art event)";
      tagBuffer();
      fillTree();
      _buffer.clear();
    }
    double frac = (n_pulses > 0) ? static_cast<double>(n_pileup) / n_pulses : 0.0;
    mf::LogInfo("TagPileupPulses")
      << "pulses: " << n_pulses << ", tagged pileup: " << n_pileup
      << " (" << 100.0 * frac << "%)";
  }
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TagPileupPulses)
