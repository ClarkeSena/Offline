// Produces an empty StepPointMCCollection. Used to fill unused input tags
// when running the HPGe digitizer with a single-source simulation.
// Original author: Sena Clarke

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {
  class EmptyStepPointMCs : public art::EDProducer {
  public:
    struct Config {};
    using Parameters = art::EDProducer::Table<Config>;
    explicit EmptyStepPointMCs(const Parameters& conf)
      : art::EDProducer{conf} {
      produces<StepPointMCCollection>();
    }
  private:
    void produce(art::Event& event) override {
      event.put(std::make_unique<StepPointMCCollection>());
    }
  };
}

DEFINE_ART_MODULE(mu2e::EmptyStepPointMCs)
