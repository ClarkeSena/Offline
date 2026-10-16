// Plot where Geant4 deposits energy in the LaBr crystal, in the transverse (x-y) plane,
// overlaid on the crystal's footprint. Each bin is weighted by the ionizing energy deposited
// (keV), NOT by step count: an unweighted step-count map is dominated by the photon gun's
// pencil-beam track, whose zero-deposit transport steps all pile onto the central axis and
// swamp the physically distributed deposition. Weighting by deposit removes that artifact and
// shows the true spatial pattern (Compton recoils + photoabsorption from the gun, and the
// volume-distributed La-138 self-activity).
//
// Positions and deposits come from the "steps" TTree written by the LaBrTruthSpectrum analyzer
// ("TruthLaBr") in FullSimToHPGeAndLaBrWaveforms.fcl. That analyzer already applies the
// LaBr-side x cut and the crystal-cylinder bounds and keeps only steps with ionizingEdep > 0,
// so no further selection is needed here. Branches: x,y,z [mm], edep [keV], source (0=photon
// gun, 1=La-138).
//
// The LaBr crystal is a solid cylinder (G4Tubs, identity rotation) of radius 38.1 mm and
// length 76.2 mm, centred at (x,y,z) = (-3863.36, 0, 40891.50). With identity rotation the
// cylinder axis lies along global z, so its x-y cross-section is a disk of radius 38.1 mm
// centred at (-3863.36, 0).
//
// Output: LaBrTruthPositions_xy.pdf           (combined deposit-weighted map)
//         LaBrTruthPositions_xy_split.pdf     (photon gun vs La-138 deposit-weighted panels)
//         LaBrTruthPositions_xy_overlay.pdf   (both sources, one view, colour-coded scatter)
//         LaBrTruthPositions_singlePhotons.pdf (two gun photons' full deposit chains, x-z)
// Usage: root -l -b -q 'plotLaBrTruthPositions.C("nts.FullSimHPGeLaBr.root")'
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TBox.h"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

namespace {
  // LaBr crystal geometry (global Mu2e coords, mm).
  const double crystalCentreX = -3863.36;
  const double crystalCentreY = 0.0;
  const double crystalR       = 38.1;

  // Draw the crystal x-y cross-section (circle + centre marker) on the current pad. The
  // colour is configurable so it stays visible against either a COLZ palette (red) or the
  // two-colour source overlay (black, to avoid clashing with the point colours).
  void drawCrystal(int color = kRed) {
    TEllipse* circle = new TEllipse(crystalCentreX, crystalCentreY, crystalR, crystalR);
    circle->SetFillStyle(0);
    circle->SetLineColor(color);
    circle->SetLineWidth(3);
    circle->Draw();
    TMarker* centre = new TMarker(crystalCentreX, crystalCentreY, 5); // 5 = cross
    centre->SetMarkerColor(color);
    centre->SetMarkerSize(2);
    centre->Draw();
  }

  // LaBr crystal longitudinal extent (z, mm). Identity rotation => axis along global z.
  const double crystalCentreZ = 40891.50;
  const double crystalL       = 76.2;

  // One deposit step, kept for per-photon grouping after the single tree pass.
  struct Step { double x, y, z, edep; int source; unsigned eventId, primaryId; };

  // Draw all deposits of a single photon in the x-z (longitudinal) plane on the current pad.
  // z is the beam/crystal axis (horizontal), x is transverse (vertical). Marker area scales
  // with deposited energy; deposits are joined in z-order with a faint line to suggest the
  // photon's forward path (the tree has no time ordering, so this is approximate). The
  // crystal's x-z cross-section is drawn as a rectangle outline.
  void drawPhoton(std::vector<Step>& steps, const char* tag) {
    std::sort(steps.begin(), steps.end(),
              [](const Step& a, const Step& b) { return a.z < b.z; });

    double eTot = 0, eMax = 0;
    for (const Step& s : steps) { eTot += s.edep; if (s.edep > eMax) eMax = s.edep; }

    const double zpad = 8.0, xpad = 8.0;
    const double zlo = crystalCentreZ - crystalL / 2, zhi = crystalCentreZ + crystalL / 2;
    TH2D* fr = new TH2D(Form("fr_%s", tag),
                        Form("%s;z [mm] (crystal axis);x [mm]", tag),
                        10, zlo - zpad, zhi + zpad,
                        10, crystalCentreX - crystalR - xpad, crystalCentreX + crystalR + xpad);
    fr->Draw();

    // Crystal x-z cross-section outline.
    TBox* box = new TBox(zlo, crystalCentreX - crystalR, zhi, crystalCentreX + crystalR);
    box->SetFillStyle(0);
    box->SetLineColor(kGray + 2);
    box->SetLineWidth(2);
    box->Draw("l");

    // Approximate path line in z-order.
    TGraph* path = new TGraph();
    for (const Step& s : steps) path->SetPoint(path->GetN(), s.z, s.x);
    path->SetLineColor(kGray + 1);
    path->SetLineStyle(2);
    path->Draw("L SAME");

    // One marker per deposit, area ~ edep.
    for (const Step& s : steps) {
      TMarker* m = new TMarker(s.z, s.x, 20);
      m->SetMarkerColor(kAzure + 1);
      m->SetMarkerSize(0.6 + 2.4 * (eMax > 0 ? s.edep / eMax : 0.0));
      m->Draw();
    }

    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(0.035);
    tl.DrawLatex(0.13, 0.84,
                 Form("%zu deposits, #SigmaE = %.0f keV", steps.size(), eTot));
  }

  // Allocate a deposit-weighted x-y histogram framed a little wider than the crystal.
  TH2D* makeMap(const char* name, const char* title) {
    const double pad = 10.0;
    TH2D* h = new TH2D(name, Form("%s;x [mm];y [mm];deposited energy [keV] / bin", title),
                       120, crystalCentreX - crystalR - pad, crystalCentreX + crystalR + pad,
                       120, crystalCentreY - crystalR - pad, crystalCentreY + crystalR + pad);
    return h;
  }
}

void plotLaBrTruthPositions(const char* filename = "nts.FullSimHPGeLaBr.root") {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return;
  }

  TTree* t = (TTree*)f->Get("TruthLaBr/steps");
  if (!t) {
    std::cerr << "Cannot find TruthLaBr/steps (is LaBrTruthSpectrum in the job?)" << std::endl;
    return;
  }

  double x = 0, y = 0, z = 0, edep = 0;
  int source = 0;
  unsigned eventId = 0, primaryId = 0;
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);
  t->SetBranchAddress("edep", &edep);
  t->SetBranchAddress("source", &source);
  // primaryId/eventId may be absent in older files; grouping by photon is then skipped.
  bool hasPrimaryId = t->GetBranch("primaryId") && t->GetBranch("eventId");
  if (hasPrimaryId) {
    t->SetBranchAddress("eventId", &eventId);
    t->SetBranchAddress("primaryId", &primaryId);
  }

  TH2D* hAll = makeMap("hAll", "LaBr Geant4 energy deposits (all sources)");
  TH2D* hGun = makeMap("hGun", "Photon-gun deposits");
  TH2D* hLa  = makeMap("hLa",  "La-138 deposits");

  // Per-source scatter points for the two-colour overlay (one point per kept step).
  TGraph* gGun = new TGraph();
  TGraph* gLa  = new TGraph();

  // Group deposits by individual primary photon: key = (eventId, primaryId).
  std::map<std::pair<unsigned, unsigned>, std::vector<Step>> byPhoton;

  Long64_t n = t->GetEntries();
  double eGun = 0, eLa = 0;
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    hAll->Fill(x, y, edep);
    if (source == 0) { hGun->Fill(x, y, edep); eGun += edep; gGun->SetPoint(gGun->GetN(), x, y); }
    else             { hLa->Fill(x, y, edep);  eLa  += edep; gLa->SetPoint(gLa->GetN(), x, y); }
    if (hasPrimaryId)
      byPhoton[{eventId, primaryId}].push_back({x, y, z, edep, source, eventId, primaryId});
  }
  std::cout << "Steps: " << n
            << " | photon-gun deposit: " << eGun / 1e3 << " MeV"
            << ", La-138 deposit: " << eLa / 1e3 << " MeV" << std::endl;

  // --- Combined deposit-weighted map ---------------------------------------------------
  TCanvas* c = new TCanvas("c", "LaBr deposit map", 800, 750);
  c->SetRightMargin(0.15);
  hAll->Draw("COLZ");
  drawCrystal();

  TLatex tl;
  tl.SetNDC();
  tl.SetTextSize(0.028);
  tl.DrawLatex(0.15, 0.93,
               Form("crystal centre (x,y) = (%.2f, %.1f) mm", crystalCentreX, crystalCentreY));
  c->SaveAs("LaBrTruthPositions_xy.pdf");

  // --- Source split: photon gun vs La-138 ----------------------------------------------
  TCanvas* cs = new TCanvas("cs", "LaBr deposit map by source", 1400, 650);
  cs->Divide(2, 1);
  cs->cd(1)->SetRightMargin(0.15);
  hGun->Draw("COLZ");
  drawCrystal();
  cs->cd(2)->SetRightMargin(0.15);
  hLa->Draw("COLZ");
  drawCrystal();
  cs->SaveAs("LaBrTruthPositions_xy_split.pdf");

  // --- Two-colour overlay: photon gun vs La-138 in one x-y view ------------------------
  // A scatter (not deposit-weighted): every kept step is one point, coloured by source, so
  // you can see at a glance which deposits come from the gun beam vs the La-138 activity.
  const double pad = 10.0;
  TH2D* frame = new TH2D("frame", "LaBr Geant4 truth deposits by source;x [mm];y [mm]",
                         10, crystalCentreX - crystalR - pad, crystalCentreX + crystalR + pad,
                         10, crystalCentreY - crystalR - pad, crystalCentreY + crystalR + pad);
  gGun->SetMarkerColor(kAzure + 1);
  gGun->SetMarkerStyle(20);
  gGun->SetMarkerSize(0.4);
  gLa->SetMarkerColor(kOrange + 7);
  gLa->SetMarkerStyle(20);
  gLa->SetMarkerSize(0.4);

  TCanvas* co = new TCanvas("co", "LaBr deposits by source", 800, 750);
  frame->Draw();
  gGun->Draw("P SAME");
  gLa->Draw("P SAME");
  drawCrystal(kBlack);

  TLegend* leg = new TLegend(0.15, 0.80, 0.55, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gGun, Form("photon gun (%lld steps)", (Long64_t)gGun->GetN()), "p");
  leg->AddEntry(gLa,  Form("La-138 (%lld steps)", (Long64_t)gLa->GetN()), "p");
  leg->Draw();
  co->SaveAs("LaBrTruthPositions_xy_overlay.pdf");

  // --- Single-photon interaction histories ---------------------------------------------
  // Pick the two photon-gun primaries that deposit in the most distinct steps (the richest
  // interaction chains) and draw each one's deposits on its own pad, x-z view. This shows a
  // single photon Compton-scattering and finally photoabsorbing, leaving several spatially
  // separated deposits - exactly the multi-deposit behaviour the steps tree now resolves.
  if (hasPrimaryId) {
    std::vector<std::pair<unsigned, unsigned>> gunPhotons;
    for (auto& kv : byPhoton) {
      bool isGun = !kv.second.empty() && kv.second.front().source == 0;
      if (isGun && kv.second.size() >= 2) gunPhotons.push_back(kv.first);
    }
    std::sort(gunPhotons.begin(), gunPhotons.end(),
              [&](const std::pair<unsigned, unsigned>& a, const std::pair<unsigned, unsigned>& b) {
                return byPhoton[a].size() > byPhoton[b].size();
              });

    if (gunPhotons.size() < 2) {
      std::cout << "Only " << gunPhotons.size()
                << " multi-deposit gun photon(s) found; skipping per-photon plot." << std::endl;
    } else {
      TCanvas* cp = new TCanvas("cp", "Single-photon interaction histories", 1400, 650);
      cp->Divide(2, 1);
      for (int p = 0; p < 2; ++p) {
        auto key = gunPhotons[p];
        cp->cd(p + 1);
        drawPhoton(byPhoton[key],
                   Form("Gun photon: event %u, primary %u", key.first, key.second));
      }
      cp->SaveAs("LaBrTruthPositions_singlePhotons.pdf");
      std::cout << "Wrote LaBrTruthPositions_singlePhotons.pdf (2 gun photons, "
                << byPhoton[gunPhotons[0]].size() << " and "
                << byPhoton[gunPhotons[1]].size() << " deposits)" << std::endl;
    }
  } else {
    std::cout << "steps tree has no primaryId branch (old file); "
              << "rebuild + rerun the sim for per-photon plots." << std::endl;
  }

  std::cout << "Wrote LaBrTruthPositions_xy.pdf, LaBrTruthPositions_xy_split.pdf, "
            << "and LaBrTruthPositions_xy_overlay.pdf" << std::endl;
}
