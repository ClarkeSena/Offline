// Concatenates several GenParticleCollections into one, so multiple independent
// generators (e.g. PhotonGun + La138Source) can feed a single g4run, whose
// primaryTag points at this merged collection.
// Original author: Sena Clarke

// stdlib includes
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Sequence.h"

// Offline includes
#include "Offline/MCDataProducts/inc/GenParticle.hh"

namespace mu2e {
  class GenParticleMerger : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Sequence<art::InputTag> inputTags{Name("inputTags"), Comment("GenParticleCollections to concatenate")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit GenParticleMerger(const Parameters& conf);
    void produce(art::Event& event) override;
  private:
    std::vector<art::ProductToken<GenParticleCollection>> tokens;
  };

  GenParticleMerger::GenParticleMerger(const Parameters& conf):
    art::EDProducer(conf) {
      for (const auto& tag : conf().inputTags())
        tokens.push_back(consumes<GenParticleCollection>(tag));
      produces<GenParticleCollection>();
    };

  void GenParticleMerger::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    for (const auto& token : tokens) {
      auto const& coll = event.getProduct(token);
      output->insert(output->end(), coll.begin(), coll.end());
    }
    event.put(std::move(output));
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenParticleMerger)
