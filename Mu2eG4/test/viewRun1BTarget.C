// 3D view of the Run 1B hanging target in the full Mu2e geometry.
// Shows only the target parts (Run1BTarget*) plus nearby context: the TS5
// downstream end wall, Collimator 5, the OPA and the DS2 rails.
// Also prints each target part's position in the G4 world frame (mm).
//
// Make the GDML first (mu2e -c Offline/Mu2eG4/fcl/gdmldump_run1_b_hangingTarget.fcl), then
//   root -l 'Offline/Mu2eG4/test/viewRun1BTarget.C("mu2e_run1_b_hangingTarget_closed.gdml")'
// Pass false as the second argument to hide the context volumes.
// Needs OpenGL over XQuartz (see Offline/STMMC/src/viewSTMClean.C for the setup).

#include "TEnv.h"
#include "TGeoManager.h"
#include "TGeoIterator.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TString.h"

#include <iostream>

void viewRun1BTarget(const char* gdmlFile = "mu2e_run1_b_hangingTarget_closed.gdml",
                     bool showContext = true) {
  gEnv->SetValue("Viewer3D.DefaultDrawOption", "ogl");
  TGeoManager* geom = TGeoManager::Import(gdmlFile);
  if (!geom) { std::cerr << "Could not import " << gdmlFile << std::endl; return; }

  struct Piece { const char* name; int color; int transparency; };
  const Piece target[] = {
    {"Run1BTargetPlate",       kMagenta + 1, 0},
    {"Run1BTargetTrackTop",    kGray + 2,    0},
    {"Run1BTargetTrackBottom", kGray + 2,    0},
    {"Run1BTargetTrolleyS",    kOrange,      0},
    {"Run1BTargetTrolleyN",    kOrange,      0},
    {"Run1BTargetHangerS",     kCyan + 1,    0},
    {"Run1BTargetHangerN",     kCyan + 1,    0},
  };
  const Piece context[] = {
    {"TS5DownstreamEndwall",  kRed - 7,   70},
    {"TS5DownstreamEndwall2", kRed - 7,   70},
    {"Coll51",                kGreen - 6, 70},
    {"protonabs3",            kBlue - 7,  80},
    {"NorthRailDS2",          kGray,      60},
    {"SouthRailDS2",          kGray,      60},
  };

  TIter next(geom->GetListOfVolumes());
  TGeoVolume* v = nullptr;
  while ((v = (TGeoVolume*)next())) v->SetVisibility(kFALSE);

  auto show = [&](const Piece& p) {
    TGeoVolume* k = geom->GetVolume(p.name);
    if (!k) { std::cerr << "WARNING: volume '" << p.name << "' not found.\n"; return; }
    k->SetVisibility(kTRUE);
    k->SetLineColor(p.color);
    k->SetFillColor(p.color);
    k->SetTransparency(p.transparency);
  };
  for (const Piece& p : target) show(p);
  if (showContext) for (const Piece& p : context) show(p);

  // Print global positions of the target parts. TGeo works in cm, so x10 for mm.
  TGeoIterator it(geom->GetTopVolume());
  TGeoNode* node = nullptr;
  while ((node = it.Next())) {
    TString name = node->GetVolume()->GetName();
    if (!name.BeginsWith("Run1BTarget")) continue;
    const double* t = it.GetCurrentMatrix()->GetTranslation();
    std::cout << name << "  centre in world, mm (" << 10 * t[0] << ", " << 10 * t[1] << ", "
              << 10 * t[2] << ")"
              << "  material " << node->GetVolume()->GetMaterial()->GetName() << std::endl;
  }

  geom->SetVisLevel(20);
  geom->GetTopVolume()->Draw("ogl");
}
