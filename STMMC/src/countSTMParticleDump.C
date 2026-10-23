// Particle-type breakdown of an STMParticleDump output (DumpSTMParticles.fcl).
// For each volume id in the dump (virtual detector id, or STMDet crystal volume) prints the number
// of entries per PDG code, the number of distinct input events, and the median / max kinetic energy.
// With firstCrossingOnly : true (the dump default) one entry = one track crossing that volume.
// Usage: root -l -b -q 'countSTMParticleDump.C("nts.dump_vd101.root")'
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
  std::string pdgName(int pdg) {
    static const std::map<int, std::string> names = {
      {22, "gamma"}, {11, "e-"}, {-11, "e+"}, {13, "mu-"}, {-13, "mu+"}, {2112, "n"},
      {2212, "p"}, {211, "pi+"}, {-211, "pi-"}, {111, "pi0"}, {12, "nu_e"}, {-12, "anti-nu_e"},
      {14, "nu_mu"}, {-14, "anti-nu_mu"}, {1000010020, "deuteron"}, {1000010030, "triton"},
      {1000020040, "alpha"}, {1000020030, "He3"}};
    auto it = names.find(pdg);
    return it != names.end() ? it->second : (pdg > 1000000000 ? "ion" : "other");
  }
}

void countSTMParticleDump(const char* filename) {
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) { std::cout << "Cannot open " << filename << std::endl; return; }
  TTree* t = (TTree*)f->Get("DumpParticles/particles");
  TTree* info = (TTree*)f->Get("DumpParticles/info");
  if (!t) { std::cout << "No DumpParticles/particles tree in " << filename << std::endl; return; }

  Long64_t nEvents = 0;
  if (info) { info->SetBranchAddress("nEvents", &nEvents); info->GetEntry(0); }

  Int_t vd = 0, pdg = 0;
  Long64_t eventId = 0;
  Double_t E = 0;
  t->SetBranchAddress("virtualdetectorId", &vd);
  t->SetBranchAddress("pdgId", &pdg);
  t->SetBranchAddress("eventId", &eventId);
  t->SetBranchAddress("E", &E);

  std::map<int, std::map<int, std::vector<double>>> energies; // vd -> pdg -> kinetic E
  std::map<int, std::set<Long64_t>> eventsPerVd;
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    energies[vd][pdg].push_back(E);
    eventsPerVd[vd].insert(eventId);
  }

  std::cout << "File " << filename << ": " << nEvents << " input events, " << t->GetEntries()
            << " entries" << std::endl;
  for (auto& [v, byPdg] : energies) {
    long total = 0;
    for (auto& [p, es] : byPdg) total += es.size();
    std::cout << "\n=== volume " << v << ": " << total << " entries in " << eventsPerVd[v].size()
              << " events" << std::endl;
    std::cout << std::left << std::setw(12) << "pdgId" << std::setw(11) << "name" << std::right
              << std::setw(9) << "count" << std::setw(9) << "frac"
              << std::setw(16) << "median E [MeV]" << std::setw(14) << "max E [MeV]" << std::endl;
    std::vector<std::pair<int, std::vector<double>*>> rows;
    for (auto& [p, es] : byPdg) rows.push_back({p, &es});
    std::sort(rows.begin(), rows.end(),
              [](const auto& a, const auto& b) { return a.second->size() > b.second->size(); });
    for (auto& [p, es] : rows) {
      std::sort(es->begin(), es->end());
      std::cout << std::left << std::setw(12) << p << std::setw(11) << pdgName(p) << std::right
                << std::setw(9) << es->size() << std::setw(9) << std::fixed << std::setprecision(3)
                << double(es->size()) / total << std::setw(16) << std::setprecision(4)
                << (*es)[es->size() / 2] << std::setw(14) << es->back() << std::endl;
    }
  }
}
