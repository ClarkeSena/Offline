#ifndef MCDataProducts_STMMCPulse_hh
#define MCDataProducts_STMMCPulse_hh
//
// MC-truth record of a single STM scintillator/semiconductor pulse: the total
// energy one decay (all its primaries) deposits in the crystal, its arrival time on the
// concatenated-waveform time axis, and whether it piles up with a neighbouring
// pulse (i.e. arrives close enough that the offline pulse finder cannot resolve
// the two as separate peaks). Produced by the STM pileup-truth chain
// (LaBrTruthPulses -> TagPileupPulses).
//
// Original author: Sena Clarke
//

#include <vector>
#include <Rtypes.h>

namespace mu2e {

  class STMMCPulse {
  public:
    // Source classification of the originating primary. photonGun covers every source
    // EXTERNAL to the crystal (photon gun, calibration source, beam); La138 and Ac227 are
    // the crystal's intrinsic activity.
    enum Source { unknown = -1, photonGun = 0, La138 = 1, Ac227 = 2 };

    STMMCPulse() = default;
    STMMCPulse(double time, double energy, unsigned primaryId, int source, unsigned eventId) :
      _time(time), _energy(energy), _primaryId(primaryId), _source(source), _eventId(eventId) {}

    double   time()      const { return _time; }      // absolute time [ns] on the concatenated grid
    double   energy()    const { return _energy; }    // total ionizing energy deposited [keV]
    unsigned primaryId() const { return _primaryId; } // SimParticle id of the (lowest-id) primary of the decay
    int      source()    const { return _source; }    // STMMCPulse::Source
    unsigned eventId()   const { return _eventId; }   // micropulse (art event) the pulse originated in

    bool   isPileup()  const { return _isPileup; }    // overlaps an unresolvable neighbour
    int    nOverlap()  const { return _nOverlap; }    // number of neighbours it overlaps
    double dtNearest() const { return _dtNearest; }   // |dt| to nearest overlapping neighbour [ns]
    double window()    const { return _window; }      // forward shadow window w_after(E) used [ns]

    void setPileup(bool v)      { _isPileup = v; }
    void setNOverlap(int n)     { _nOverlap = n; }
    void setDtNearest(double d) { _dtNearest = d; }
    void setWindow(double w)    { _window = w; }

  private:
    double   _time      = 0.0;
    double   _energy    = 0.0;
    unsigned _primaryId = 0;
    int      _source    = unknown;
    unsigned _eventId   = 0;
    bool     _isPileup  = false;
    int      _nOverlap  = 0;
    double   _dtNearest = 0.0;
    double   _window    = 0.0;
  };

  typedef std::vector<STMMCPulse> STMMCPulseCollection;
}
#endif
