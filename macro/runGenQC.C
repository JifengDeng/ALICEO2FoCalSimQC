#pragma cling load("./O2FoCalGenQC.so")

#include <iostream>
#include "TSystem.h"
// Header lives under src/; path is relative to this macro file.
#include "../src/O2FoCalGenQC.h"

// Argument order matches runQC.sh: kine, hits, output [, geometry].
void runGenQC(const char* kineFile = "input/o2sim_Kine.root",
              const char* hitsFile = "input/o2sim_HitsFOC.root",
              const char* outputFile = "output/Results.root",
              const char* geoFile = "../FOCAL/simulation/geometryFiles/geometry_Sheets.txt") {
  std::cout << "Initializing O2FoCalGenQC..." << std::endl;
  Long64_t init_time = static_cast<Long64_t>(gSystem->Now());
  std::cout << "Start time (ms): " << init_time << std::endl;

  O2FoCalGenQC* genQC = new O2FoCalGenQC(outputFile);

  // Geometry is required to decode col/row/layer/segment from hit positions.
  genQC->LoadGeometry(geoFile);

  // Read kinematics, compute per track, and fill histograms in one pass.
  genQC->ReadKine(kineFile);

  // Read FOCAL hits, fill event/hit/geo/match histograms.
  if (!genQC->ReadFOCHits(hitsFile)) {
    std::cout << "O2FoCalGenQC: ReadFOCHits failed" << std::endl;
  }

  genQC->WriteOutput();

  Long64_t end_time = static_cast<Long64_t>(gSystem->Now());
  std::cout << "End time (ms): " << end_time << std::endl;
  std::cout << "Time taken: " << end_time - init_time << " ms" << std::endl;

  delete genQC;

  std::cout << "O2FoCalGenQC executed successfully" << std::endl;
}
