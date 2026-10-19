// Dumps STM detector geometry from the GeometryService.
// Prints crystal centre positions in Mu2e coordinates for both HPGe and LaBr,
// along with all crystal dimensions and the mu2eOriginInWorld transform.
// Original author: Sena Clarke

#include <cmath>
#include <iomanip>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "CLHEP/Units/PhysicalConstants.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/STMGeom/inc/STM.hh"

namespace mu2e {
  class STMGeometryDump : public art::EDAnalyzer {
  public:
    struct Config {};
    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit STMGeometryDump(const Parameters& conf)
      : art::EDAnalyzer{conf} {}
  private:
    void analyze(art::Event const& event) override;
    bool printed_ = false;
  };

  void STMGeometryDump::analyze(art::Event const&) {
    if (printed_) return;
    printed_ = true;

    GeomHandle<WorldG4> world;
    GeomHandle<STM> stm;

    auto const& mu2eOrigin = world->mu2eOriginInWorld();

    auto const* hpge = stm->getHPGeDetectorPtr();
    auto const* labr = stm->getLaBrDetectorPtr();
    auto const* ssc  = stm->getSTM_SSCPtr();
    auto const* front = stm->getFrontShieldingPtr();

    int w = 45;
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "  STM Geometry Dump\n";
    std::cout << "================================================================\n";

    std::cout << "\n--- Coordinate Transform ---\n";
    std::cout << std::left << std::setw(w) << "  mu2eOriginInWorld.x [mm]"     << mu2eOrigin.x() << "\n";
    std::cout << std::left << std::setw(w) << "  mu2eOriginInWorld.y [mm]"     << mu2eOrigin.y() << "\n";
    std::cout << std::left << std::setw(w) << "  mu2eOriginInWorld.z [mm]"     << mu2eOrigin.z() << "\n";

    std::cout << "\n--- SSC / Shielding ---\n";
    std::cout << std::left << std::setw(w) << "  SSC originInMu2e [mm]"        << ssc->originInMu2e() << "\n";
    std::cout << std::left << std::setw(w) << "  SSC offset_Spot [mm]"         << ssc->offset_Spot() << "\n";
    std::cout << std::left << std::setw(w) << "  SSC FrontToWall [mm]"         << ssc->FrontToWall() << "\n";
    std::cout << std::left << std::setw(w) << "  SSC ZGap [mm]"               << ssc->ZGap() << "\n";
    std::cout << std::left << std::setw(w) << "  Front_Thickness [mm]"         << front->Front_Thickness() << "\n";
    double Front_T = front->Front_Thickness() + ssc->ZGap();
    std::cout << std::left << std::setw(w) << "  Front_T (Thickness+ZGap) [mm]" << Front_T << "\n";

    std::cout << "\n================================================================\n";
    std::cout << "  HPGe Detector\n";
    std::cout << "================================================================\n";
    std::cout << std::left << std::setw(w) << "  originInMu2e.x [mm]"          << hpge->originInMu2e().x() << "\n";
    std::cout << std::left << std::setw(w) << "  originInMu2e.y [mm]"          << hpge->originInMu2e().y() << "\n";
    std::cout << std::left << std::setw(w) << "  originInMu2e.z [mm]"          << hpge->originInMu2e().z() << "\n";
    std::cout << std::left << std::setw(w) << "  CrystalR [mm]"               << hpge->CrystalR() << "\n";
    std::cout << std::left << std::setw(w) << "  CrystalL [mm]"               << hpge->CrystalL() << "\n";
    std::cout << std::left << std::setw(w) << "  HoleR [mm]"                  << hpge->HoleR() << "\n";
    std::cout << std::left << std::setw(w) << "  HoleL [mm]"                  << hpge->HoleL() << "\n";
    std::cout << std::left << std::setw(w) << "  Z_HPGe [mm]"                 << hpge->Z_HPGe() << "\n";
    std::cout << std::left << std::setw(w) << "  offset_HPGe [mm]"            << hpge->offset_HPGe() << "\n";
    std::cout << std::left << std::setw(w) << "  EndcapR [mm]"                << hpge->EndcapR() << "\n";
    std::cout << std::left << std::setw(w) << "  EndcapL [mm]"                << hpge->EndcapL() << "\n";
    std::cout << std::left << std::setw(w) << "  WindowD [mm]"                << hpge->WindowD() << "\n";
    std::cout << std::left << std::setw(w) << "  EndcapD [mm]"                << hpge->EndcapD() << "\n";
    std::cout << std::left << std::setw(w) << "  AirD [mm]"                   << hpge->AirD() << "\n";
    std::cout << std::left << std::setw(w) << "  Capsule_Windowthick [mm]"    << hpge->Capsule_Windowthick() << "\n";
    std::cout << std::left << std::setw(w) << "  Capsule_Wallthick [mm]"      << hpge->Capsule_Wallthick() << "\n";
    std::cout << std::left << std::setw(w) << "  Capsule_Endthick [mm]"       << hpge->Capsule_Endthick() << "\n";
    std::cout << std::left << std::setw(w) << "  Capsule_Walllength [mm]"     << hpge->Capsule_Walllength() << "\n";

    double hpge_internal = hpge->WindowD() + hpge->AirD() + hpge->Capsule_Windowthick() + hpge->CrystalL() / 2.0;
    double hpge_proj = hpge_internal * std::sqrt(2.0) / 2.0;
    double hpge_crystal_x = hpge->originInMu2e().x() - hpge_proj;
    double hpge_crystal_y = 0.0;
    double hpge_crystal_z = hpge->originInMu2e().z() + Front_T + hpge_proj;

    std::cout << "\n  --- Derived Crystal Centre ---\n";
    std::cout << std::left << std::setw(w) << "  internal depth [mm]"          << hpge_internal << "\n";
    std::cout << std::left << std::setw(w) << "  internal * sqrt(2)/2 [mm]"    << hpge_proj << "\n";
    std::cout << std::left << std::setw(w) << "  crystalCentreX [mm]"          << hpge_crystal_x << "\n";
    std::cout << std::left << std::setw(w) << "  crystalCentreY [mm]"          << hpge_crystal_y << "\n";
    std::cout << std::left << std::setw(w) << "  crystalCentreZ [mm]"          << hpge_crystal_z << "\n";
    std::cout << std::left << std::setw(w) << "  rotation"                     << "45 deg about Y" << "\n";

    std::cout << "\n================================================================\n";
    std::cout << "  LaBr Detector\n";
    std::cout << "================================================================\n";
    std::cout << std::left << std::setw(w) << "  originInMu2e.x [mm]"          << labr->originInMu2e().x() << "\n";
    std::cout << std::left << std::setw(w) << "  originInMu2e.y [mm]"          << labr->originInMu2e().y() << "\n";
    std::cout << std::left << std::setw(w) << "  originInMu2e.z [mm]"          << labr->originInMu2e().z() << "\n";
    std::cout << std::left << std::setw(w) << "  CrystalR [mm]"               << labr->CrystalR() << "\n";
    std::cout << std::left << std::setw(w) << "  CrystalL [mm]"               << labr->CrystalL() << "\n";
    std::cout << std::left << std::setw(w) << "  Z_LaBr [mm]"                 << labr->Z_LaBr() << "\n";
    std::cout << std::left << std::setw(w) << "  offset_LaBr [mm]"            << labr->offset_LaBr() << "\n";
    std::cout << std::left << std::setw(w) << "  EndcapR [mm]"                << labr->EndcapR() << "\n";
    std::cout << std::left << std::setw(w) << "  EndcapL [mm]"                << labr->EndcapL() << "\n";
    std::cout << std::left << std::setw(w) << "  WindowD [mm]"                << labr->WindowD() << "\n";
    std::cout << std::left << std::setw(w) << "  EndcapD [mm]"                << labr->EndcapD() << "\n";
    std::cout << std::left << std::setw(w) << "  AirD [mm]"                   << labr->AirD() << "\n";

    double labr_crystal_x = labr->originInMu2e().x();
    double labr_crystal_y = 0.0;
    double labr_crystal_z = labr->originInMu2e().z() + Front_T + labr->WindowD() + labr->AirD() + labr->CrystalL() / 2.0;

    std::cout << "\n  --- Derived Crystal Centre ---\n";
    std::cout << std::left << std::setw(w) << "  crystalCentreX [mm]"          << labr_crystal_x << "\n";
    std::cout << std::left << std::setw(w) << "  crystalCentreY [mm]"          << labr_crystal_y << "\n";
    std::cout << std::left << std::setw(w) << "  crystalCentreZ [mm]"          << labr_crystal_z << "\n";
    std::cout << std::left << std::setw(w) << "  rotation"                     << "None (identity)" << "\n";

    std::cout << "\n================================================================\n";
    std::cout << "  Comparison\n";
    std::cout << "================================================================\n";
    std::cout << std::left << std::setw(35) << "" << std::setw(15) << "HPGe" << std::setw(15) << "LaBr" << "\n";
    std::cout << std::left << std::setw(35) << "  crystalCentreX [mm]" << std::setw(15) << hpge_crystal_x << std::setw(15) << labr_crystal_x << "\n";
    std::cout << std::left << std::setw(35) << "  crystalCentreY [mm]" << std::setw(15) << hpge_crystal_y << std::setw(15) << labr_crystal_y << "\n";
    std::cout << std::left << std::setw(35) << "  crystalCentreZ [mm]" << std::setw(15) << hpge_crystal_z << std::setw(15) << labr_crystal_z << "\n";
    std::cout << std::left << std::setw(35) << "  crystalR [mm]" << std::setw(15) << hpge->CrystalR() << std::setw(15) << labr->CrystalR() << "\n";
    std::cout << std::left << std::setw(35) << "  crystalL [mm]" << std::setw(15) << hpge->CrystalL() << std::setw(15) << labr->CrystalL() << "\n";
    std::cout << "\n";
  }
}

DEFINE_ART_MODULE(mu2e::STMGeometryDump)
