#ifndef O2_FOCAL_GEN_QC_H
#define O2_FOCAL_GEN_QC_H

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"

#include <vector>

class TTree;

namespace o2
{
template <class T>
class MCTrackT;
using MCTrack = MCTrackT<float>;
}

namespace o2::focal
{
class Geometry;
}

class O2FoCalGenQC {
public:
  O2FoCalGenQC(const TString& outputFilename);
  ~O2FoCalGenQC();

  // Load FOCAL geometry (geometry_Sheets.txt). Called automatically by ReadFOCHits if needed.
  void LoadGeometry(const TString& geometryFile = "");

  // One pass over the MCTrack branch. Kinematics are computed and filled per track.
  void ReadKine(const TString& inputFilename);

  // Read FOCAL hits, decode geometry, match to cached kinematics, fill QC histograms.
  bool ReadFOCHits(const TString& inputFilename);

  void WriteOutput();

private:
  struct Kinematics {
    double pt = 0.;
    double eta = 0.;
    double y = 0.;
    double phi = 0.;
    double theta = 0.;
    double p = 0.;
    double energy = 0.;
    bool etaDefined = false;
    bool yDefined = false;
  };

  // Compact per-track cache used to match hits to MC particles.
  struct TrackKine {
    int pdg = 0;
    int process = 0;
    int mother = -1;
    bool isPrimary = false;
    Kinematics kine;
  };

  void BookKineHistograms(Long64_t nEvents);
  void BookHitHistograms(Long64_t nEvents);
  Kinematics CalculateKinematics(const o2::MCTrack& trk) const;
  void FillTrack(const o2::MCTrack& trk);
  TString ResolveGeometryFile(const TString& geometryFile) const;
  const char* FindHitBranchName(TTree* tree) const;
  void ResetHistPointers();

  TString fOutputFilename;
  TFile* fOutputFile = nullptr;
  o2::focal::Geometry* fGeometry = nullptr; //! owned by Geometry singleton
  int fNLayers = 0;
  int fNSegments = 0;

  std::vector<std::vector<TrackKine>> fTracksByEvent;

  // --- Kine (original MCTrack quantities + derived kinematics) ---
  TH1D* fHPDG = nullptr;
  TH1D* fHPt = nullptr;
  TH1D* fHEta = nullptr;
  TH1D* fHY = nullptr;
  TH1D* fHPhi = nullptr;
  TH1D* fHTheta = nullptr;
  TH1D* fHP = nullptr;
  TH1D* fHEnergy = nullptr;
  TH1D* fHProcess = nullptr;
  TH1D* fHPx = nullptr;
  TH1D* fHPy = nullptr;
  TH1D* fHPz = nullptr;
  TH1D* fHVx = nullptr;
  TH1D* fHVy = nullptr;
  TH1D* fHVz = nullptr;
  TH1D* fHVt = nullptr;
  TH1D* fHIsPrimary = nullptr;

  // --- Event-level hit QC ---
  TH1D* fHEventETot = nullptr;
  TH1D* fHEventEECAL = nullptr;
  TH1D* fHEventEHCAL = nullptr;
  TH1D* fHEventNHits = nullptr;
  TH1D* fHEventNHitPrimaries = nullptr;

  // --- Hit-level stored quantities ---
  TH1D* fHHitTrackId = nullptr;
  TH1D* fHHitX = nullptr;
  TH1D* fHHitY = nullptr;
  TH1D* fHHitZ = nullptr;
  TH1D* fHHitELoss = nullptr;
  TH1D* fHHitTime = nullptr;
  TH1D* fHHitSubsystem = nullptr;
  TH1D* fHHitPrimary = nullptr;
  TH1D* fHHitInitialEnergy = nullptr;
  TH2D* fHHitXY = nullptr;
  TH3D* fHHitXYZ = nullptr;

  // --- Geometry ---
  TH1D* fHGeoCol = nullptr;
  TH1D* fHGeoRow = nullptr;
  TH1D* fHGeoLayer = nullptr;
  TH1D* fHGeoSegment = nullptr;
  TH1D* fHGeoInDetector = nullptr;
  TH2D* fHGeoLayerMult = nullptr;
  std::vector<TH2D*> fHLayerXY;

  // --- Matching / acceptance / response ---
  TH1D* fHMatchPDG = nullptr;
  TH1D* fHMatchPt = nullptr;
  TH1D* fHMatchEta = nullptr;
  TH1D* fHMatchY = nullptr;
  TH1D* fHMatchPhi = nullptr;
  TH1D* fHMatchTheta = nullptr;
  TH1D* fHMatchP = nullptr;
  TH1D* fHMatchEnergy = nullptr;
  TH1D* fHMatchProcess = nullptr;
  TH1D* fHMatchIsPrimary = nullptr;
  TH1D* fHMatchResponse = nullptr;
  TH1D* fHMatchResponsePrimary = nullptr;
  TH2D* fHMatchEdepVsE = nullptr;
  TH1D* fHAccGenEta = nullptr;
  TH1D* fHAccHitEta = nullptr;
  TH1D* fHAccGenE = nullptr;
  TH1D* fHAccHitE = nullptr;
  TH1D* fHAccGenPt = nullptr;
  TH1D* fHAccHitPt = nullptr;
};

#endif // O2_FOCAL_GEN_QC_H
