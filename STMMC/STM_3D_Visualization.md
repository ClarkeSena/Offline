# STM 3D Visualization

How to make 3D figures of the STM geometry (SSC collimator + LaBr/HPGe crystals) and the
photon beam, for screenshots and talks/thesis figures.

## TL;DR (the working recipe)

Render the real geometry from the exported GDML in **ROOT's TGeo OpenGL viewer**, displayed
over XQuartz. Geant4's own visualization does **not** work in this build (see below), so do
not use the `g4run.visualization` / `.mac` route.

```bash
# --- one-time on your Mac: enable indirect GLX so OpenGL works over XQuartz ---
defaults write org.xquartz.X11 enable_iglx -bool true
defaults read  org.xquartz.X11 enable_iglx        # must print 1
killall Xquartz                                    # fully quit
open -a XQuartz                                    # relaunch BEFORE ssh

# --- connect with X forwarding ---
ssh -Y sclarke@mu2egpvm01.fnal.gov

# --- on the gpvm, in that shell ---
export LIBGL_ALWAYS_INDIRECT=1                     # needed every session (add to ~/.bashrc)
glxgears                                            # sanity check: gears must spin, then Ctrl-C
cd /exp/mu2e/app/users/sclarke/Simulation
muse setup
root -l 'Offline/STMMC/src/viewSTMClean.C("mu2e_common.gdml")'
```

A GL window opens with the SSC (gray), LaBr (blue), HPGe (orange), and a yellow photon-beam
line. Rotate/zoom with the mouse and screenshot with `Cmd-Shift-4`.

## Files

- `Offline/STMMC/src/viewSTMClean.C` — the working ROOT macro. Imports `mu2e_common.gdml`,
  hides all volumes, re-shows only `TungstenSSCPV` (SSC), `fLaBrPV` (LaBr), `fHPGePV` (HPGe),
  draws them via the TGeo OpenGL viewer, and overlays the photon beam as a `TPolyLine3D`.
- `mu2e_common.gdml` — the exported geometry (in the Simulation working dir). Produced by a
  g4run job with `physics.producers.g4run.debug.writeGDML: true` and a `GDMLFileName`.

The macro draws `stmDownstreamEnvelope` as the (invisible) mother so the three visible
daughters keep their real relative positions; `SetVisLevel(10)` makes the recursion reach the
crystals.

## The photon beam overlay

The beam endpoints are read **from the geometry itself**: SSC downstream face
(`sPos.z + sHalf.z`) to LaBr back face (`lPos.z + lHalf.z`), at LaBr's transverse `(x, y)`.
This keeps the line in the same units/frame as the drawn volumes (no mm/cm conversion to get
wrong) and matches the path the `PhotonGun` fires along toward LaBr.

It is the **nominal gun direction**, a single straight pencil line (the gun fires with
`px = py = 0`), not simulated tracks. To show where photons actually interact, the next step is
to overlay the real energy deposits inside LaBr as 3D points from the `TruthLaBr/steps` tree in
`nts.FullSimHPGeLaBr.root` (branches `x, y, z, edep, source`); those coordinates are in mm and
must be scaled into the geometry's units.

## Key gotcha: this build has no Geant4 visualization

The Muse build `al9-prof-e29-p092` is compiled **without** any `G4VIS_USE_OPENGL*` flag. In
`Mu2eG4/src/Mu2eG4_module.cc`, the entire vis path — creating `G4VisManager` **and** executing
`visualization.initMacro` — is inside `#if (defined G4VIS_USE_OPENGLX || ...)`. So in this
build:

- any `g4run.visualization.initMacro` (`.mac` file) is silently ignored;
- no window opens, `/vis/list` prints nothing, no VRML/file export is produced;
- the per-event `"Enter a character to go to the next event"` prompt still appears (it is NOT
  behind the guard) but does nothing visual — feed it with `yes n | mu2e ...` to auto-advance.

Do not chase OGLSX / OGLSQt / VRML2FILE / HepRepFile in this build; they all route through the
disabled code. Use the ROOT TGeo + GDML route above instead.

## XQuartz / OpenGL troubleshooting

Symptom seen: `No matching fbConfigs or visuals found` / `glx: failed to create drisw screen` /
`glXCreateContext failed`. Cause: indirect GLX (IGLX) is disabled on XQuartz by default, and
libGL was trying the local software driver instead of sending GL to XQuartz.

Fix (both halves are required):
1. **Mac:** `defaults write org.xquartz.X11 enable_iglx -bool true`, then **fully quit**
   (`killall Xquartz`) and relaunch XQuartz, then reconnect `ssh -Y` *after* it is running.
2. **gpvm:** `export LIBGL_ALWAYS_INDIRECT=1` in the same shell you run ROOT in.
3. Test with `glxgears` — it must spin before ROOT GL will work.

If it still fails, update XQuartz to >= 2.8.5 and redo step 1. (XQuartz lives at
`/Applications/Utilities/XQuartz.app`.)

## Approaches that were tried and rejected

- **Geant4 interactive vis (OGLSQt/OGLSX)** and **Geant4 file drivers (VRML2FILE)** — dead;
  the build has no G4 vis (see above).
- **ROOT TGeo CPU ray-tracer (`Raytrace()`) to a headless PNG** — produced a blank image; the
  ray-traced pixels do not survive the batch `Update`/`SaveAs` repaint. Not reliable here.
- **ROOT `TPolyLine3D` wireframe to a headless PNG** — works without OpenGL, but only gives
  wireframe outlines, not shaded solids. Kept only as a no-display fallback idea; superseded
  once XQuartz OpenGL was fixed.

## Fallbacks if OpenGL ever breaks again

- ParaView (installed on the Mac) is a native 3D viewer needing no XQuartz/GLX, if the geometry
  can be exported to a format it reads.
- Headless `TPolyLine3D` wireframe PNG (no GL/display) for a quick non-interactive figure.
