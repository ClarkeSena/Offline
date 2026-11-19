#ifndef Mu2eG4_constructRun1BTarget_hh
#define Mu2eG4_constructRun1BTarget_hh
//
// Free function to construct the simplified Run 1B movable ("hanging")
// target: Al plate, tracks, trolleys and hanger plates on the downstream
// face of the TS5 downstream end wall. Parameters: run1bTarget.* in
// Mu2eG4/geom/Run1BHangingTarget_v01.txt; enabled by hasRun1BTarget.
//

namespace mu2e {

  class SimpleConfig;

  void constructRun1BTarget(SimpleConfig const& config);

}

#endif /* Mu2eG4_constructRun1BTarget_hh */
