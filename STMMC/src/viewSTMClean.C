// Clean interactive 3D view of the STM detectors + the photon beam, for screenshots.
// Loads the real Mu2e geometry (mu2e_common.gdml, written by g4run with debug.writeGDML:true)
// and shows ONLY the SSC collimator and the LaBr + HPGe crystals, plus a yellow line for the
// photon beam, hiding everything else so the pieces are framed on their own.
// Uses ROOT's TGeo OpenGL viewer (same approach as the working stmDownstreamEnvelope macro).
//
// The beam endpoints are taken FROM THE GEOMETRY (SSC downstream face -> LaBr back face, at
// LaBr's transverse position), so it is correctly placed and in the same units as the drawn
// volumes - no mm/cm guessing. That is the path the PhotonGun fires along toward LaBr.
//
// Needs working OpenGL over XQuartz. If "glXCreateContext failed":
//   Mac:  defaults write org.xquartz.X11 enable_iglx -bool true   (then fully quit + relaunch XQuartz)
//   gpvm: export LIBGL_ALWAYS_INDIRECT=1                          (then `glxgears` should spin)
//
// Run from the directory that contains the .gdml:
//   cd /exp/mu2e/app/users/sclarke/Simulation
//   root -l 'Offline/STMMC/src/viewSTMClean.C("mu2e_common.gdml")'
// Then rotate/zoom in the GL window and screenshot.
// Original author: Sena Clarke

#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TGeoIterator.h"
#include "TPolyLine3D.h"
#include "TEnv.h"
#include "TList.h"
#include "TString.h"

#include <cstring>
#include <iostream>

namespace {
  // Global translation + bounding-box half-sizes of the first placement of a volume.
  bool findVolume(TGeoManager* geom, const char* name, double pos[3], double half[3]) {
    TGeoIterator it(geom->GetTopVolume());
    TGeoNode* node = nullptr;
    while ((node = it.Next())) {
      if (!strcmp(node->GetVolume()->GetName(), name)) {
        const double* t = it.GetCurrentMatrix()->GetTranslation();
        pos[0] = t[0]; pos[1] = t[1]; pos[2] = t[2];
        TGeoBBox* bb = (TGeoBBox*)node->GetVolume()->GetShape();
        half[0] = bb->GetDX(); half[1] = bb->GetDY(); half[2] = bb->GetDZ();
        return true;
      }
    }
    std::cerr << "WARNING: volume '" << name << "' not found.\n";
    return false;
  }
}

void viewSTMClean(const char* gdmlFile = "mu2e_common.gdml") {
  gEnv->SetValue("Viewer3D.DefaultDrawOption", "ogl");

  TGeoManager* geom = TGeoManager::Import(gdmlFile);
  if (!geom) { std::cerr << "Could not import " << gdmlFile << std::endl; return; }

  // STM logical-volume names (confirmed present in this GDML).
  struct Piece { const char* name; int color; int transparency; };
  const Piece pieces[] = {
    {"TungstenSSCPV", kGray + 1, 30},   // SSC, slightly transparent so the beam shows through
    {"fLaBrPV",       kAzure + 1, 0},   // LaBr crystal
    {"fHPGePV",       kOrange + 7, 0},  // HPGe crystal
  };

  // Hide everything, then re-enable just the three pieces.
  TIter next(geom->GetListOfVolumes());
  TGeoVolume* v = nullptr;
  while ((v = (TGeoVolume*)next())) v->SetVisibility(kFALSE);
  for (const Piece& p : pieces) {
    TGeoVolume* k = geom->GetVolume(p.name);
    if (!k) { std::cerr << "WARNING: volume '" << p.name << "' not found.\n"; continue; }
    k->SetVisibility(kTRUE);
    k->SetLineColor(p.color);
    k->SetFillColor(p.color);
    k->SetTransparency(p.transparency);
  }

  // Build the photon beam from geometry positions: SSC downstream face -> LaBr back face,
  // travelling along z at LaBr's transverse (x,y). Same units as the geometry, so it lines up.
  double sPos[3], sHalf[3], lPos[3], lHalf[3];
  bool haveSSC  = findVolume(geom, "TungstenSSCPV", sPos, sHalf);
  bool haveLaBr = findVolume(geom, "fLaBrPV", lPos, lHalf);
  if (haveSSC && haveLaBr) {
    double zStart = sPos[2] + sHalf[2];   // SSC downstream face
    double zEnd   = lPos[2] + lHalf[2];   // LaBr back face
    TPolyLine3D* beam = new TPolyLine3D(2);
    beam->SetPoint(0, lPos[0], lPos[1], zStart);
    beam->SetPoint(1, lPos[0], lPos[1], zEnd);
    beam->SetLineColor(kYellow);
    beam->SetLineWidth(5);
    beam->Draw();   // added to the pad before the GL viewer is created below
  }

  // Draw the geometry last so the GL viewer builds its scene from the full pad (incl. beam).
  geom->SetVisLevel(10);
  TGeoVolume* env = geom->GetVolume("stmDownstreamEnvelope");
  if (env) {
    env->SetVisibility(kFALSE);
    env->SetVisDaughters(kTRUE);
    env->Draw("ogl");
  } else {
    std::cerr << "WARNING: 'stmDownstreamEnvelope' not found; drawing the top volume instead.\n";
    geom->GetTopVolume()->Draw("ogl");
  }
}
