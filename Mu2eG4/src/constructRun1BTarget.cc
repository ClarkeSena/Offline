//
// Free function to construct the simplified Run 1B movable ("hanging") target.
//
// Six parts, translated from New_Target/Run1B_target_simplified_*.gdml:
//   Run1BTargetPlate          540 x 540 x 17.5 mm plate, R100 corners
//   Run1BTargetTrackTop       10135A71 track, lips facing down
//   Run1BTargetTrackBottom    10135A71 track, lips facing up
//   Run1BTargetTrolleyS/N     10135A81 trolley: box body + two polycone wheels
//   Run1BTargetHangerS/N      hanger plates, unions of boxes
//
// Local frame: origin on the TS5 axis at the downstream face of the TS5
// downstream end wall, z downstream, y up (Mu2e axes). All parts are placed
// directly in DS2Vacuum. Parameters: run1bTarget.* (Run1BHangingTarget_v01.txt).
//

#include <array>
#include <iostream>
#include <string>
#include <vector>

// art / messagefacility
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"

// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/BeamlineGeom/inc/StraightSection.hh"
#include "Offline/BeamlineGeom/inc/TransportSolenoid.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eG4/inc/constructRun1BTarget.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"

// G4
#include "Geant4/G4Box.hh"
#include "Geant4/G4Colour.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4UnionSolid.hh"

namespace mu2e {

  namespace {

    struct NestFlags {
      bool visible;
      bool solid;
      bool forceAuxEdgeVisible;
      bool placePV;
      bool doSurfaceCheck;
    };

    // Place a solid that was built here (union) into the parent.
    VolumeInfo placeSolid(std::string const& name,
                          G4VSolid* solid,
                          G4Material* material,
                          G4RotationMatrix const* rot,
                          CLHEP::Hep3Vector const& posInMu2e,
                          VolumeInfo const& parent,
                          G4Colour const& colour,
                          NestFlags const& f) {
      VolumeInfo info(name, posInMu2e - parent.centerInMu2e(), parent.centerInWorld);
      info.solid = solid;
      finishNesting(info, material, rot, info.centerInParent, parent.logical, 0,
                    f.visible, colour, f.solid, f.forceAuxEdgeVisible, f.placePV,
                    f.doSurfaceCheck);
      return info;
    }

    // Union of boxes; box 0 is the reference, later boxes are offset from its
    // centre with their x offset multiplied by parkDir.
    G4VSolid* makeBoxUnion(std::string const& name,
                           std::vector<double> const& halfDims,
                           std::vector<double> const& offsets,
                           double parkDir) {
      std::size_t const nBoxes = halfDims.size() / 3;
      G4VSolid* result = new G4Box(name + "_b0", halfDims[0] * CLHEP::mm,
                                   halfDims[1] * CLHEP::mm, halfDims[2] * CLHEP::mm);
      for (std::size_t i = 1; i < nBoxes; ++i) {
        G4Box* box = new G4Box(name + "_b" + std::to_string(i), halfDims[3 * i] * CLHEP::mm,
                               halfDims[3 * i + 1] * CLHEP::mm, halfDims[3 * i + 2] * CLHEP::mm);
        G4ThreeVector const off(parkDir * offsets[3 * i] * CLHEP::mm,
                                offsets[3 * i + 1] * CLHEP::mm, offsets[3 * i + 2] * CLHEP::mm);
        result = new G4UnionSolid(name + "_u" + std::to_string(i), result, box, nullptr, off);
      }
      return result;
    }

  } // namespace

  void constructRun1BTarget(SimpleConfig const& config) {

    int const verbosity = config.getInt("run1bTarget.verbosityLevel", 0);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry(config, "run1bTarget", "run1bTarget");
    NestFlags const flags{geomOptions->isVisible("run1bTarget"),
                          geomOptions->isSolid("run1bTarget"),
                          geomOptions->forceAuxEdgeVisible("run1bTarget"),
                          geomOptions->placePV("run1bTarget"),
                          geomOptions->doSurfaceCheck("run1bTarget")};

    Mu2eG4Helper* helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));
    AntiLeakRegistry& reg = helper->antiLeakRegistry();
    VolumeInfo const& parent = helper->locateVolInfo("DS2Vacuum");

    // ---- Frame origin: TS5 axis at the downstream face of the TS5 end wall
    GeomHandle<Beamline> beamline;
    StraightSection const* ts5 = beamline->getTS().getTSCryo<StraightSection>(
      TransportSolenoid::TSRegion::TS5, TransportSolenoid::TSRadialPart::IN);
    CLHEP::Hep3Vector const origin(
      ts5->getGlobal().x() + config.getDouble("run1bTarget.xOffset", 0.) * CLHEP::mm,
      0.,
      ts5->getGlobal().z() + ts5->getHalfLength() +
        config.getDouble("run1bTarget.zOffset", 0.) * CLHEP::mm);

    // The parts must lie inside DS2Vacuum, which starts at the TS5 end.
    GeomHandle<DetectorSolenoid> ds;
    double const ds2FrontZ = ds->vac_zLocDs23Split() - 2. * ds->vac_halfLengthDs2();
    if (origin.z() < ds2FrontZ - 1.e-6) {
      throw cet::exception("GEOM") << __func__ << ": Run 1B target frame origin z = " << origin.z()
                                   << " is upstream of the DS2Vacuum front, z = " << ds2FrontZ
                                   << "\n";
    }

    double const travel = config.getDouble("run1bTarget.travel") * CLHEP::mm;
    double const parkDir = config.getDouble("run1bTarget.parkDir");
    if (parkDir != 1. && parkDir != -1.) {
      throw cet::exception("GEOM") << __func__ << ": run1bTarget.parkDir must be +1 or -1, got "
                                   << parkDir << "\n";
    }

    // Position in Mu2e of a moving part (plate, trolleys, hangers) and of a fixed one (tracks)
    auto movingPos = [&](double xCAD, double y, double z) {
      return origin + CLHEP::Hep3Vector(parkDir * (xCAD * CLHEP::mm + travel), y * CLHEP::mm,
                                        z * CLHEP::mm);
    };
    auto fixedPos = [&](double xCAD, double y, double z) {
      return origin + CLHEP::Hep3Vector(parkDir * xCAD * CLHEP::mm, y * CLHEP::mm, z * CLHEP::mm);
    };

    G4Material* plateMat = findMaterialOrThrow(config.getString("run1bTarget.plate.materialName"));
    G4Material* trackMat = findMaterialOrThrow(config.getString("run1bTarget.track.materialName"));
    G4Material* trolleyMat =
      findMaterialOrThrow(config.getString("run1bTarget.trolley.materialName"));
    G4Material* hangerMat =
      findMaterialOrThrow(config.getString("run1bTarget.hanger.materialName"));

    std::vector<VolumeInfo> built;

    // ---- Plate
    {
      int const n = config.getInt("run1bTarget.plate.nVertices");
      std::vector<double> u, v;
      config.getVectorDouble("run1bTarget.plate.outlineU", u, n);
      config.getVectorDouble("run1bTarget.plate.outlineV", v, n);
      CLHEP::Hep3Vector const pos = movingPos(config.getDouble("run1bTarget.plate.xCAD"),
                                              config.getDouble("run1bTarget.plate.y"),
                                              config.getDouble("run1bTarget.plate.z"));
      built.push_back(nestExtrudedSolid(
        "Run1BTargetPlate", config.getDouble("run1bTarget.plate.halfThickness") * CLHEP::mm, u, v,
        plateMat, nullptr, pos - parent.centerInMu2e(), parent, 0, flags.visible,
        G4Colour::Magenta(), flags.solid, flags.forceAuxEdgeVisible, flags.placePV,
        flags.doSurfaceCheck));
    }

    // ---- Tracks: profile (U across, V = absolute y) extruded along local z,
    // rotated 90 deg about y so the length runs along x. The profile is
    // symmetric in U, so the sign of the rotation does not matter.
    {
      CLHEP::HepRotation* trackRot = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
      trackRot->rotateY(90. * CLHEP::degree);
      for (std::string const tag : {"trackTop", "trackBottom"}) {
        std::string const key = "run1bTarget." + tag;
        int const n = config.getInt(key + ".nVertices");
        std::vector<double> u, v;
        config.getVectorDouble(key + ".outlineU", u, n);
        config.getVectorDouble(key + ".outlineV", v, n);
        CLHEP::Hep3Vector const pos =
          fixedPos(config.getDouble(key + ".xCAD"), 0., config.getDouble(key + ".z"));
        std::string const name =
          (tag == "trackTop") ? "Run1BTargetTrackTop" : "Run1BTargetTrackBottom";
        built.push_back(nestExtrudedSolid(
          name, config.getDouble(key + ".halfLength") * CLHEP::mm, u, v, trackMat, trackRot,
          pos - parent.centerInMu2e(), parent, 0, flags.visible, G4Colour::Grey(), flags.solid,
          flags.forceAuxEdgeVisible, flags.placePV, flags.doSurfaceCheck));
      }
    }

    // ---- Trolleys
    {
      std::vector<double> body, xCAD, wz, wr;
      config.getVectorDouble("run1bTarget.trolley.bodyHalfDims", body, 3);
      config.getVectorDouble("run1bTarget.trolley.xCAD", xCAD, 2);
      int const nPlanes = config.getInt("run1bTarget.trolley.wheel.nPlanes");
      config.getVectorDouble("run1bTarget.trolley.wheel.z", wz, nPlanes);
      config.getVectorDouble("run1bTarget.trolley.wheel.rOut", wr, nPlanes);
      for (auto& z : wz) z *= CLHEP::mm;
      for (auto& r : wr) r *= CLHEP::mm;
      std::vector<double> const wrIn(nPlanes, 0.);
      double const wx = parkDir * config.getDouble("run1bTarget.trolley.wheel.xOffset") * CLHEP::mm;
      double const wy = config.getDouble("run1bTarget.trolley.wheel.yOffset") * CLHEP::mm;

      std::array<std::string, 2> const side{"S", "N"};
      for (std::size_t i = 0; i < 2; ++i) {
        std::string const name = "Run1BTargetTrolley" + side[i];
        G4Box* bodySolid = new G4Box(name + "_body", body[0] * CLHEP::mm, body[1] * CLHEP::mm,
                                     body[2] * CLHEP::mm);
        G4Polycone* wheel = new G4Polycone(name + "_wheel", 0., CLHEP::twopi, nPlanes, wz.data(),
                                           wrIn.data(), wr.data());
        G4VSolid* u1 =
          new G4UnionSolid(name + "_u1", bodySolid, wheel, nullptr, G4ThreeVector(-wx, wy, 0.));
        G4VSolid* trolley =
          new G4UnionSolid(name + "_u2", u1, wheel, nullptr, G4ThreeVector(wx, wy, 0.));
        CLHEP::Hep3Vector const pos = movingPos(xCAD[i], config.getDouble("run1bTarget.trolley.y"),
                                                config.getDouble("run1bTarget.trolley.z"));
        built.push_back(placeSolid(name, trolley, trolleyMat, nullptr, pos, parent,
                                   G4Colour::Yellow(), flags));
      }
    }

    // ---- Hanger plates
    for (std::string const side : {"S", "N"}) {
      std::string const key = "run1bTarget.hanger" + side;
      int const nBoxes = config.getInt(key + ".nBoxes");
      std::vector<double> halfDims, offsets;
      config.getVectorDouble(key + ".boxHalfDims", halfDims, 3 * nBoxes);
      config.getVectorDouble(key + ".boxOffsets", offsets, 3 * nBoxes);
      std::string const name = "Run1BTargetHanger" + side;
      G4VSolid* hanger = makeBoxUnion(name, halfDims, offsets, parkDir);
      CLHEP::Hep3Vector const pos = movingPos(config.getDouble(key + ".xCAD"),
                                              config.getDouble(key + ".y"),
                                              config.getDouble(key + ".z"));
      built.push_back(
        placeSolid(name, hanger, hangerMat, nullptr, pos, parent, G4Colour::Cyan(), flags));
    }

    if (verbosity > 0) {
      std::cout << __func__ << " Run 1B target frame origin in Mu2e: " << origin
                << ", travel = " << travel << " mm, parkDir = " << parkDir
                << ", DS2Vacuum front z = " << ds2FrontZ << std::endl;
      for (auto const& vi : built) {
        std::cout << __func__ << "   " << vi.name << " centre in Mu2e " << vi.centerInMu2e()
                  << ", material " << vi.logical->GetMaterial()->GetName()
                  << ", mass " << vi.logical->GetMass() / CLHEP::kg << " kg" << std::endl;
      }
    }
  }

} // namespace mu2e
