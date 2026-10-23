#include "O2FoCalGenQC.h"

#include "FOCALBase/Geometry.h"
#include "FOCALBase/Hit.h"
#include "SimulationDataFormat/MCTrack.h"

#include "TBranch.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TInterpreter.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace
{
constexpr double kELossToGeV = 1e-9; // hit energy loss is stored in eV

void AddEnvInclude(const char* envName)
{
  const char* val = gSystem->Getenv(envName);
  if (val && val[0] != '\0') {
    gSystem->AddIncludePath(Form("-I%s/include", val));
  }
}

bool HasCollectionProxy(const char* className)
{
  TClass* cl = TClass::GetClass(className);
  return cl && cl->IsLoaded() && cl->GetCollectionProxy();
}

bool EnsureCollectionProxy(const char* className, const char* headers)
{
  if (HasCollectionProxy(className)) {
    return true;
  }
  std::cout << "O2FoCalGenQC: generating dictionary for " << className << std::endl;
  gInterpreter->GenerateDictionary(className, headers);
  gROOT->GetListOfClasses()->FindObject(className);
  if (!HasCollectionProxy(className)) {
    TString alt(className);
    alt.ReplaceAll("std::", "");
    return HasCollectionProxy(alt.Data());
  }
  return true;
}

bool EnsureHitVectorIO()
{
  AddEnvInclude("O2_ROOT");
  AddEnvInclude("BOOST_ROOT");
  AddEnvInclude("FAIRROOT_ROOT");
  if (const char* o2 = gSystem->Getenv("O2_ROOT")) {
    if (o2[0] != '\0') {
      gSystem->AddIncludePath(Form("-I%s/include/GPU", o2));
    }
  }

  if (gSystem->Load("libO2FOCALBase") < 0) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: cannot load libO2FOCALBase" << std::endl;
    return false;
  }

  if (!EnsureCollectionProxy("std::vector<o2::focal::Hit>", "FOCALBase/Hit.h;vector")) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: no CollectionProxy for vector<o2::focal::Hit>" << std::endl;
    return false;
  }
  return true;
}
}

O2FoCalGenQC::O2FoCalGenQC(const TString& outputFilename)
  : fOutputFilename(outputFilename)
{
  fOutputFile = new TFile(fOutputFilename, "RECREATE");
}

O2FoCalGenQC::~O2FoCalGenQC()
{
  if (fOutputFile) {
    if (fOutputFile->IsOpen()) {
      fOutputFile->Close();
    }
    delete fOutputFile;
    fOutputFile = nullptr;
  }
}

void O2FoCalGenQC::BookKineHistograms(Long64_t nEvents)
{
  if (fHPt) {
    return;
  }
  if (!fOutputFile || !fOutputFile->IsOpen()) {
    std::cerr << "O2FoCalGenQC::BookKineHistograms ERROR: output file is not open" << std::endl;
    return;
  }

  fOutputFile->cd();
  TDirectory* dir = fOutputFile->mkdir("kine");
  dir->cd();

  (void)nEvents;
  fHPDG = new TH1D("hPDG", "PDG code;PDG;tracks", 10001, -5000.5, 5000.5);
  fHPt = new TH1D("hPt", "Transverse momentum;p_{T} (GeV/c);tracks", 400, 0., 20.);
  fHEta = new TH1D("hEta", "Pseudorapidity;#eta;tracks", 280, -6., 8.);
  fHY = new TH1D("hY", "Rapidity;y;tracks", 280, -6., 8.);
  fHPhi = new TH1D("hPhi", "Azimuth;#phi (rad);tracks", 360, 0., TMath::TwoPi());
  fHTheta = new TH1D("hTheta", "Polar angle;#theta (rad);tracks", 180, 0., TMath::Pi());
  fHP = new TH1D("hP", "Momentum;p (GeV/c);tracks", 400, 0., 200.);
  fHEnergy = new TH1D("hEnergy", "Energy;E (GeV);tracks", 400, 0., 200.);
  fHProcess = new TH1D("hProcess", "Production process;TMCProcess;tracks", 64, -0.5, 63.5);

  fHPx = new TH1D("hPx", "p_{x};p_{x} (GeV/c);tracks", 400, -20., 20.);
  fHPy = new TH1D("hPy", "p_{y};p_{y} (GeV/c);tracks", 400, -20., 20.);
  fHPz = new TH1D("hPz", "p_{z};p_{z} (GeV/c);tracks", 400, -200., 200.);
  fHVx = new TH1D("hVx", "Vertex x;v_{x} (cm);tracks", 200, -20., 20.);
  fHVy = new TH1D("hVy", "Vertex y;v_{y} (cm);tracks", 200, -20., 20.);
  fHVz = new TH1D("hVz", "Vertex z;v_{z} (cm);tracks", 400, -50., 50.);
  fHVt = new TH1D("hVt", "Vertex time;t;tracks", 200, 0., 100.);
  fHIsPrimary = new TH1D("hIsPrimary", "Is primary;isPrimary;tracks", 2, -0.5, 1.5);

  fOutputFile->cd();
  TDirectory* dirMatch = fOutputFile->mkdir("match");
  dirMatch->cd();
  fHAccGenEta = new TH1D("hAccGenEta", "Generated primaries;#eta;primaries", 280, -6., 8.);
  fHAccGenE = new TH1D("hAccGenE", "Generated primaries;E (GeV);primaries", 400, 0., 200.);
  fHAccGenPt = new TH1D("hAccGenPt", "Generated primaries;p_{T} (GeV/c);primaries", 400, 0., 20.);

  fOutputFile->cd();
}

void O2FoCalGenQC::BookHitHistograms(Long64_t nEvents)
{
  if (fHHitELoss) {
    return;
  }
  if (!fOutputFile || !fOutputFile->IsOpen()) {
    std::cerr << "O2FoCalGenQC::BookHitHistograms ERROR: output file is not open" << std::endl;
    return;
  }

  double xMin = -60., xMax = 60.;
  double yMin = -60., yMax = 60.;
  double zMin = 700., zMax = 900.;
  int nCol = 2000;
  int nRow = 2000;
  fNLayers = fNLayers > 0 ? fNLayers : 54;
  fNSegments = fNSegments > 0 ? fNSegments : 7;

  if (fGeometry) {
    const double xHalf = 0.5 * fGeometry->getFOCALSizeX() + fGeometry->getMiddleTowerOffset() + 5.;
    const double yHalf = 0.5 * fGeometry->getFOCALSizeY() + 5.;
    const double z0 = fGeometry->getFOCALZ0();
    const double zHalf = 0.5 * fGeometry->getFOCALSizeZ() + 5.;
    xMin = -xHalf;
    xMax = xHalf;
    yMin = -yHalf;
    yMax = yHalf;
    zMin = z0 - zHalf;
    zMax = z0 + zHalf;
    fNLayers = std::max(1, fGeometry->getNumberOfLayers());
    fNSegments = std::max(1, fGeometry->getVirtualNSegments());
    int nColMax = nCol;
    int nRowMax = nRow;
    for (int iseg = 0; iseg < fGeometry->getVirtualNSegments(); ++iseg) {
      auto [ok, nc, nr] = fGeometry->getVirtualNColRow(iseg);
      if (ok) {
        nColMax = std::max(nColMax, nc + 2);
        nRowMax = std::max(nRowMax, nr + 2);
      }
    }
    nCol = nColMax;
    nRow = nRowMax;
  }

  fOutputFile->cd();
  TDirectory* dirEvent = fOutputFile->mkdir("event");
  dirEvent->cd();
  fHEventETot = new TH1D("hEventETot", "Total deposited energy per event;E_{tot} (GeV);events", 400, 0., 200.);
  fHEventEECAL = new TH1D("hEventEECAL", "ECAL deposited energy per event;E_{ECAL} (GeV);events", 400, 0., 200.);
  fHEventEHCAL = new TH1D("hEventEHCAL", "HCAL deposited energy per event;E_{HCAL} (GeV);events", 400, 0., 200.);
  fHEventNHits = new TH1D("hEventNHits", "Hits per event;N_{hit};events", 400, 0., 40000.);
  fHEventNHitPrimaries = new TH1D("hEventNHitPrimaries", "Primaries with hits per event;N_{primary}^{hit};events", 200, -0.5, 199.5);

  fOutputFile->cd();
  TDirectory* dirHit = fOutputFile->mkdir("hit");
  dirHit->cd();
  fHHitTrackId = new TH1D("hHitTrackId", "Hit track id;track id;hits", 20000, -0.5, 19999.5);
  fHHitX = new TH1D("hHitX", "Hit x;x (cm);hits", 400, xMin, xMax);
  fHHitY = new TH1D("hHitY", "Hit y;y (cm);hits", 400, yMin, yMax);
  fHHitZ = new TH1D("hHitZ", "Hit z;z (cm);hits", 400, zMin, zMax);
  fHHitELoss = new TH1D("hHitELoss", "Hit energy loss;E_{loss} (GeV);hits", 400, 0., 0.1);
  fHHitTime = new TH1D("hHitTime", "Hit time;t (ns);hits", 200, 0., 100.);
  fHHitSubsystem = new TH1D("hHitSubsystem", "Hit subsystem (0=EPADS,1=EPIXELS,2=HCAL,3=UNKNOWN);subsystem;hits", 4, -0.5, 3.5);
  fHHitPrimary = new TH1D("hHitPrimary", "Hit primary id;primary;hits", 20000, -1.5, 19998.5);
  fHHitInitialEnergy = new TH1D("hHitInitialEnergy", "Parent energy at FOCAL entrance;E_{initial} (GeV);hits", 400, 0., 200.);
  fHHitXY = new TH2D("hHitXY", "Hit x vs y;x (cm);y (cm)", 200, xMin, xMax, 200, yMin, yMax);
  fHHitXYZ = new TH3D("hHitXYZ", "Hit position;x (cm);y (cm);z (cm)", 60, xMin, xMax, 60, yMin, yMax, 60, zMin, zMax);

  fOutputFile->cd();
  TDirectory* dirGeo = fOutputFile->mkdir("geo");
  dirGeo->cd();
  fHGeoCol = new TH1D("hCol", "Virtual column;col;hits", nCol, -0.5, nCol - 0.5);
  fHGeoRow = new TH1D("hRow", "Virtual row;row;hits", nRow, -0.5, nRow - 0.5);
  fHGeoLayer = new TH1D("hLayer", "Layer;layer;hits", fNLayers + 1, -0.5, fNLayers + 0.5);
  fHGeoSegment = new TH1D("hSegment", "Virtual segment;segment;hits", fNSegments + 1, -0.5, fNSegments + 0.5);
  fHGeoInDetector = new TH1D("hInDetector", "Hit inside virtual geometry;in detector;hits", 2, -0.5, 1.5);
  fHGeoLayerMult = new TH2D("hLayerMult", "Hits per layer per event;layer;N_{hit}(layer)",
                            fNLayers, -0.5, fNLayers - 0.5, 200, 0., 5000.);
  fHLayerXY.resize(static_cast<std::size_t>(fNLayers), nullptr);
  for (int il = 0; il < fNLayers; ++il) {
    fHLayerXY[static_cast<std::size_t>(il)] = new TH2D(Form("hLayerXY_%02d", il),
                                                       Form("Hit x vs y, layer %d;x (cm);y (cm)", il),
                                                       200, xMin, xMax, 200, yMin, yMax);
  }

  fOutputFile->cd();
  TDirectory* dirMatch = fOutputFile->GetDirectory("match");
  if (!dirMatch) {
    dirMatch = fOutputFile->mkdir("match");
  }
  dirMatch->cd();
  fHMatchPDG = new TH1D("hHitPartPDG", "PDG of tracks producing hits;PDG;tracks", 10001, -5000.5, 5000.5);
  fHMatchPt = new TH1D("hHitPartPt", "p_{T} of tracks producing hits;p_{T} (GeV/c);tracks", 400, 0., 20.);
  fHMatchEta = new TH1D("hHitPartEta", "#eta of tracks producing hits;#eta;tracks", 280, -6., 8.);
  fHMatchY = new TH1D("hHitPartY", "y of tracks producing hits;y;tracks", 280, -6., 8.);
  fHMatchPhi = new TH1D("hHitPartPhi", "#phi of tracks producing hits;#phi (rad);tracks", 360, 0., TMath::TwoPi());
  fHMatchTheta = new TH1D("hHitPartTheta", "#theta of tracks producing hits;#theta (rad);tracks", 180, 0., TMath::Pi());
  fHMatchP = new TH1D("hHitPartP", "p of tracks producing hits;p (GeV/c);tracks", 400, 0., 200.);
  fHMatchEnergy = new TH1D("hHitPartEnergy", "E of tracks producing hits;E (GeV);tracks", 400, 0., 200.);
  fHMatchProcess = new TH1D("hHitPartProcess", "Process of tracks producing hits;TMCProcess;tracks", 64, -0.5, 63.5);
  fHMatchIsPrimary = new TH1D("hHitPartIsPrimary", "Is primary (tracks producing hits);isPrimary;tracks", 2, -0.5, 1.5);
  fHMatchResponse = new TH1D("hResponse", "E_{dep}/E_{true} per hit-producing track;E_{dep}/E_{true};tracks", 200, 0., 2.);
  fHMatchResponsePrimary = new TH1D("hResponsePrimary", "E_{dep}/E_{true} per hit-producing primary;E_{dep}/E_{true};primaries", 200, 0., 2.);
  fHMatchEdepVsE = new TH2D("hEdepVsE", "E_{dep} vs E_{true};E_{true} (GeV);E_{dep} (GeV)", 200, 0., 200., 200, 0., 50.);
  fHAccHitEta = new TH1D("hAccHitEta", "Primaries with hits;#eta;primaries", 280, -6., 8.);
  fHAccHitE = new TH1D("hAccHitE", "Primaries with hits;E (GeV);primaries", 400, 0., 200.);
  fHAccHitPt = new TH1D("hAccHitPt", "Primaries with hits;p_{T} (GeV/c);primaries", 400, 0., 20.);

  (void)nEvents;
  fOutputFile->cd();
}

O2FoCalGenQC::Kinematics O2FoCalGenQC::CalculateKinematics(const o2::MCTrack& trk) const
{
  Kinematics kine;
  const double px = trk.Px();
  const double py = trk.Py();
  const double pz = trk.Pz();
  const double energy = trk.GetEnergy();

  const double pt = std::sqrt(px * px + py * py);
  const double p = std::sqrt(pt * pt + pz * pz);

  kine.pt = pt;
  kine.p = p;
  kine.energy = energy;
  kine.phi = TMath::Pi() + TMath::ATan2(-py, -px);
  kine.theta = (p == 0. || pz == 0.) ? TMath::PiOver2() : TMath::ACos(pz / p);
  if (kine.phi >= TMath::TwoPi()) {
    kine.phi = std::nextafter(TMath::TwoPi(), 0.);
  }
  if (kine.theta >= TMath::Pi()) {
    kine.theta = std::nextafter(TMath::Pi(), 0.);
  }

  if (p != TMath::Abs(pz)) {
    kine.eta = 0.5 * std::log((p + pz) / (p - pz));
    kine.etaDefined = true;
  }
  if (energy != TMath::Abs(pz)) {
    kine.y = 0.5 * std::log((energy + pz) / (energy - pz));
    kine.yDefined = true;
  }
  return kine;
}

void O2FoCalGenQC::FillTrack(const o2::MCTrack& trk)
{
  if (!fHPt) {
    std::cerr << "O2FoCalGenQC::FillTrack ERROR: histograms are not booked" << std::endl;
    return;
  }

  const Kinematics kine = CalculateKinematics(trk);
  fHPDG->Fill(trk.GetPdgCode());
  fHPt->Fill(kine.pt);
  if (kine.etaDefined) {
    fHEta->Fill(kine.eta);
  }
  if (kine.yDefined) {
    fHY->Fill(kine.y);
  }
  fHPhi->Fill(kine.phi);
  fHTheta->Fill(kine.theta);
  fHP->Fill(kine.p);
  fHEnergy->Fill(kine.energy);
  fHProcess->Fill(trk.getProcess());
  fHPx->Fill(trk.Px());
  fHPy->Fill(trk.Py());
  fHPz->Fill(trk.Pz());
  fHVx->Fill(trk.Vx());
  fHVy->Fill(trk.Vy());
  fHVz->Fill(trk.Vz());
  fHVt->Fill(trk.T());
  fHIsPrimary->Fill(trk.isPrimary() ? 1. : 0.);

  if (trk.isPrimary()) {
    if (kine.etaDefined) {
      fHAccGenEta->Fill(kine.eta);
    }
    fHAccGenE->Fill(kine.energy);
    fHAccGenPt->Fill(kine.pt);
  }
}

TString O2FoCalGenQC::ResolveGeometryFile(const TString& geometryFile) const
{
  std::vector<TString> candidates;
  if (geometryFile.Length()) {
    candidates.push_back(geometryFile);
  }
  candidates.emplace_back("geometry_Sheets.txt");
  candidates.emplace_back("../FOCAL/simulation/geometryFiles/geometry_Sheets.txt");
  candidates.emplace_back("FOCALcode/O2version/FOCAL/simulation/geometryFiles/geometry_Sheets.txt");

  for (auto path : candidates) {
    gSystem->ExpandPathName(path);
    if (!gSystem->AccessPathName(path.Data())) {
      return path;
    }
  }
  return "";
}

void O2FoCalGenQC::LoadGeometry(const TString& geometryFile)
{
  if (fGeometry) {
    return;
  }

  const TString resolved = ResolveGeometryFile(geometryFile);
  if (resolved.Length() == 0) {
    std::cerr << "O2FoCalGenQC::LoadGeometry ERROR: geometry file not found"
              << (geometryFile.Length() ? Form(" (tried %s)", geometryFile.Data()) : "") << std::endl;
    return;
  }

  fGeometry = o2::focal::Geometry::getInstance(resolved.Data());
  if (!fGeometry) {
    std::cerr << "O2FoCalGenQC::LoadGeometry ERROR: Geometry::getInstance failed for " << resolved.Data() << std::endl;
    return;
  }

  fNLayers = fGeometry->getNumberOfLayers();
  fNSegments = fGeometry->getVirtualNSegments();
  std::cout << "O2FoCalGenQC::LoadGeometry: " << resolved.Data()
            << "  nLayers=" << fNLayers
            << "  nSegments=" << fNSegments
            << "  sizeX=" << fGeometry->getFOCALSizeX()
            << "  sizeY=" << fGeometry->getFOCALSizeY()
            << "  sizeZ=" << fGeometry->getFOCALSizeZ()
            << "  z0=" << fGeometry->getFOCALZ0() << std::endl;
}

const char* O2FoCalGenQC::FindHitBranchName(TTree* tree) const
{
  if (!tree) {
    return nullptr;
  }
  if (tree->GetBranch("FOCHit")) {
    return "FOCHit";
  }
  TObjArray* branches = tree->GetListOfBranches();
  if (!branches) {
    return nullptr;
  }
  for (int i = 0; i < branches->GetEntries(); ++i) {
    auto* br = static_cast<TBranch*>(branches->At(i));
    if (!br) {
      continue;
    }
    const TString cls = br->GetClassName();
    const TString nam = br->GetName();
    if (cls.Contains("focal::Hit") || nam.Contains("Hit")) {
      return br->GetName();
    }
  }
  return nullptr;
}

void O2FoCalGenQC::ReadKine(const TString& inputFilename)
{
  if (!fOutputFile || !fOutputFile->IsOpen()) {
    std::cerr << "O2FoCalGenQC::ReadKine ERROR: output file is not open" << std::endl;
    return;
  }

  TString path = inputFilename;
  gSystem->ExpandPathName(path);
  if (gSystem->AccessPathName(path.Data())) {
    std::cerr << "O2FoCalGenQC::ReadKine ERROR: file not found: " << path.Data() << std::endl;
    return;
  }

  TFile* inFile = TFile::Open(path.Data(), "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "O2FoCalGenQC::ReadKine ERROR: cannot open " << path.Data() << std::endl;
    delete inFile;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(inFile->Get("o2sim"));
  if (!tree) {
    std::cerr << "O2FoCalGenQC::ReadKine ERROR: tree 'o2sim' not found in " << path.Data() << std::endl;
    inFile->Close();
    delete inFile;
    return;
  }

  if (!tree->GetBranch("MCTrack")) {
    std::cerr << "O2FoCalGenQC::ReadKine ERROR: branch 'MCTrack' not found" << std::endl;
    inFile->Close();
    delete inFile;
    return;
  }

  std::vector<o2::MCTrack>* tracks = nullptr;
  if (tree->SetBranchAddress("MCTrack", &tracks) < 0) {
    std::cerr << "O2FoCalGenQC::ReadKine ERROR: SetBranchAddress failed for 'MCTrack'" << std::endl;
    inFile->Close();
    delete inFile;
    return;
  }

  const Long64_t nEvents = tree->GetEntries();
  BookKineHistograms(nEvents);
  fTracksByEvent.assign(static_cast<std::size_t>(nEvents), {});

  long long nTracks = 0;
  for (Long64_t iev = 0; iev < nEvents; ++iev) {
    tree->GetEntry(iev);
    if (!tracks) {
      continue;
    }
    const int nInEvent = static_cast<int>(tracks->size());
    auto& cache = fTracksByEvent[static_cast<std::size_t>(iev)];
    cache.resize(static_cast<std::size_t>(nInEvent));
    for (int itrk = 0; itrk < nInEvent; ++itrk) {
      const auto& trk = (*tracks)[itrk];
      FillTrack(trk);
      TrackKine stored;
      stored.pdg = trk.GetPdgCode();
      stored.process = trk.getProcess();
      stored.mother = trk.getMotherTrackId();
      stored.isPrimary = trk.isPrimary();
      stored.kine = CalculateKinematics(trk);
      cache[static_cast<std::size_t>(itrk)] = stored;
      ++nTracks;
    }
  }

  std::cout << "O2FoCalGenQC::ReadKine: " << path.Data()
            << "  nEvents=" << nEvents << "  nTracks=" << nTracks << std::endl;

  tree->ResetBranchAddresses();
  delete tracks;
  inFile->Close();
  delete inFile;
}

bool O2FoCalGenQC::ReadFOCHits(const TString& inputFilename)
{
  if (!fOutputFile || !fOutputFile->IsOpen()) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: output file is not open" << std::endl;
    return false;
  }

  if (!EnsureHitVectorIO()) {
    return false;
  }

  if (!fGeometry) {
    LoadGeometry();
  }

  TString path = inputFilename;
  gSystem->ExpandPathName(path);
  if (gSystem->AccessPathName(path.Data())) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: file not found: " << path.Data() << std::endl;
    return false;
  }

  TFile* inFile = TFile::Open(path.Data(), "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: cannot open " << path.Data() << std::endl;
    delete inFile;
    return false;
  }

  TTree* tree = dynamic_cast<TTree*>(inFile->Get("o2sim"));
  if (!tree) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: tree 'o2sim' not found in " << path.Data() << std::endl;
    inFile->Close();
    delete inFile;
    return false;
  }

  const char* hitBranch = FindHitBranchName(tree);
  if (!hitBranch) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: FOCAL hit branch not found" << std::endl;
    tree->Print();
    inFile->Close();
    delete inFile;
    return false;
  }

  TBranch* br = tree->GetBranch(hitBranch);
  std::cout << "O2FoCalGenQC::ReadFOCHits: branch=" << hitBranch
            << " class=" << (br ? br->GetClassName() : "<none>") << std::endl;

  std::vector<o2::focal::Hit>* hits = nullptr;
  if (tree->SetBranchAddress(hitBranch, &hits) < 0) {
    std::cout << "O2FoCalGenQC::ReadFOCHits ERROR: SetBranchAddress failed for '" << hitBranch << "'" << std::endl;
    inFile->Close();
    delete inFile;
    return false;
  }

  const Long64_t nEvents = tree->GetEntries();
  BookHitHistograms(nEvents);

  long long nHitTotal = 0;
  long long nInside = 0;
  for (Long64_t iev = 0; iev < nEvents; ++iev) {
    tree->GetEntry(iev);
    if (!hits) {
      continue;
    }

    double eTot = 0.;
    double eECAL = 0.;
    double eHCAL = 0.;
    std::unordered_set<int> hitPrimaries;
    std::unordered_set<int> hitTracks;
    std::unordered_map<int, double> edepByTrack;
    std::unordered_map<int, double> edepByPrimary;
    std::vector<int> nPerLayer(static_cast<std::size_t>(std::max(fNLayers, 1)), 0);

    for (const auto& hit : *hits) {
      const double elossEV = static_cast<double>(hit.GetEnergyLoss());
      const double elossGeV = elossEV * kELossToGeV;
      const double x = hit.GetX();
      const double y = hit.GetY();
      const double z = hit.GetZ();
      const int trackId = hit.GetTrackID();
      const int primary = hit.getPrimary();
      const auto subsystem = hit.getSubsystem();

      fHHitTrackId->Fill(trackId);
      fHHitX->Fill(x);
      fHHitY->Fill(y);
      fHHitZ->Fill(z);
      fHHitELoss->Fill(elossGeV);
      fHHitTime->Fill(hit.GetTime());
      fHHitSubsystem->Fill(static_cast<int>(subsystem));
      fHHitPrimary->Fill(primary);
      fHHitInitialEnergy->Fill(hit.getInitialEnergy());
      fHHitXY->Fill(x, y);
      fHHitXYZ->Fill(x, y, z);

      eTot += elossGeV;
      if (subsystem == o2::focal::Hit::Subsystem_t::HCAL) {
        eHCAL += elossGeV;
      } else if (subsystem == o2::focal::Hit::Subsystem_t::EPADS || subsystem == o2::focal::Hit::Subsystem_t::EPIXELS) {
        eECAL += elossGeV;
      }

      hitTracks.insert(trackId);
      edepByTrack[trackId] += elossGeV;
      if (primary >= 0) {
        hitPrimaries.insert(primary);
        edepByPrimary[primary] += elossGeV;
      }

      bool inside = false;
      int col = -1, row = -1, layer = -1, segment = -1;
      if (fGeometry) {
        auto info = fGeometry->getVirtualInfo(x, y, z);
        inside = std::get<0>(info);
        col = std::get<1>(info);
        row = std::get<2>(info);
        layer = std::get<3>(info);
        segment = std::get<4>(info);
      }
      fHGeoInDetector->Fill(inside ? 1. : 0.);
      if (inside) {
        ++nInside;
        fHGeoCol->Fill(col);
        fHGeoRow->Fill(row);
        fHGeoLayer->Fill(layer);
        fHGeoSegment->Fill(segment);
        if (layer >= 0 && layer < fNLayers) {
          nPerLayer[static_cast<std::size_t>(layer)]++;
          if (fHLayerXY[static_cast<std::size_t>(layer)]) {
            fHLayerXY[static_cast<std::size_t>(layer)]->Fill(x, y);
          }
        }
      }
    }

    fHEventETot->Fill(eTot);
    fHEventEECAL->Fill(eECAL);
    fHEventEHCAL->Fill(eHCAL);
    fHEventNHits->Fill(static_cast<double>(hits->size()));
    fHEventNHitPrimaries->Fill(static_cast<double>(hitPrimaries.size()));
    for (int il = 0; il < fNLayers; ++il) {
      fHGeoLayerMult->Fill(il, nPerLayer[static_cast<std::size_t>(il)]);
    }
    nHitTotal += static_cast<long long>(hits->size());

    if (iev < static_cast<Long64_t>(fTracksByEvent.size())) {
      const auto& cache = fTracksByEvent[static_cast<std::size_t>(iev)];
      for (int trackId : hitTracks) {
        if (trackId < 0 || trackId >= static_cast<int>(cache.size())) {
          continue;
        }
        const auto& stored = cache[static_cast<std::size_t>(trackId)];
        fHMatchPDG->Fill(stored.pdg);
        fHMatchPt->Fill(stored.kine.pt);
        if (stored.kine.etaDefined) {
          fHMatchEta->Fill(stored.kine.eta);
        }
        if (stored.kine.yDefined) {
          fHMatchY->Fill(stored.kine.y);
        }
        fHMatchPhi->Fill(stored.kine.phi);
        fHMatchTheta->Fill(stored.kine.theta);
        fHMatchP->Fill(stored.kine.p);
        fHMatchEnergy->Fill(stored.kine.energy);
        fHMatchProcess->Fill(stored.process);
        fHMatchIsPrimary->Fill(stored.isPrimary ? 1. : 0.);
        const double edep = edepByTrack[trackId];
        fHMatchEdepVsE->Fill(stored.kine.energy, edep);
        if (stored.kine.energy > 0.) {
          fHMatchResponse->Fill(edep / stored.kine.energy);
        }
      }
      for (int primary : hitPrimaries) {
        if (primary < 0 || primary >= static_cast<int>(cache.size())) {
          continue;
        }
        const auto& stored = cache[static_cast<std::size_t>(primary)];
        if (!stored.isPrimary) {
          continue;
        }
        if (stored.kine.etaDefined) {
          fHAccHitEta->Fill(stored.kine.eta);
        }
        fHAccHitE->Fill(stored.kine.energy);
        fHAccHitPt->Fill(stored.kine.pt);
        const double edep = edepByPrimary[primary];
        if (stored.kine.energy > 0.) {
          fHMatchResponsePrimary->Fill(edep / stored.kine.energy);
        }
      }
    }
  }

  std::cout << "O2FoCalGenQC::ReadFOCHits: " << path.Data()
            << "  branch=" << hitBranch
            << "  nEvents=" << nEvents
            << "  nHits=" << nHitTotal
            << "  nInsideGeometry=" << nInside << std::endl;
  if (fTracksByEvent.empty()) {
    std::cout << "O2FoCalGenQC::ReadFOCHits: kinematics cache is empty; call ReadKine() first for matching/acceptance" << std::endl;
  }

  tree->ResetBranchAddresses();
  delete hits;
  inFile->Close();
  delete inFile;
  return true;
}

void O2FoCalGenQC::ResetHistPointers()
{
  fHPDG = nullptr;
  fHPt = nullptr;
  fHEta = nullptr;
  fHY = nullptr;
  fHPhi = nullptr;
  fHTheta = nullptr;
  fHP = nullptr;
  fHEnergy = nullptr;
  fHProcess = nullptr;
  fHPx = nullptr;
  fHPy = nullptr;
  fHPz = nullptr;
  fHVx = nullptr;
  fHVy = nullptr;
  fHVz = nullptr;
  fHVt = nullptr;
  fHIsPrimary = nullptr;
  fHEventETot = nullptr;
  fHEventEECAL = nullptr;
  fHEventEHCAL = nullptr;
  fHEventNHits = nullptr;
  fHEventNHitPrimaries = nullptr;
  fHHitTrackId = nullptr;
  fHHitX = nullptr;
  fHHitY = nullptr;
  fHHitZ = nullptr;
  fHHitELoss = nullptr;
  fHHitTime = nullptr;
  fHHitSubsystem = nullptr;
  fHHitPrimary = nullptr;
  fHHitInitialEnergy = nullptr;
  fHHitXY = nullptr;
  fHHitXYZ = nullptr;
  fHGeoCol = nullptr;
  fHGeoRow = nullptr;
  fHGeoLayer = nullptr;
  fHGeoSegment = nullptr;
  fHGeoInDetector = nullptr;
  fHGeoLayerMult = nullptr;
  fHLayerXY.clear();
  fHMatchPDG = nullptr;
  fHMatchPt = nullptr;
  fHMatchEta = nullptr;
  fHMatchY = nullptr;
  fHMatchPhi = nullptr;
  fHMatchTheta = nullptr;
  fHMatchP = nullptr;
  fHMatchEnergy = nullptr;
  fHMatchProcess = nullptr;
  fHMatchIsPrimary = nullptr;
  fHMatchResponse = nullptr;
  fHMatchResponsePrimary = nullptr;
  fHMatchEdepVsE = nullptr;
  fHAccGenEta = nullptr;
  fHAccHitEta = nullptr;
  fHAccGenE = nullptr;
  fHAccHitE = nullptr;
  fHAccGenPt = nullptr;
  fHAccHitPt = nullptr;
}

void O2FoCalGenQC::WriteOutput()
{
  if (!fOutputFile || !fOutputFile->IsOpen()) {
    std::cerr << "O2FoCalGenQC::WriteOutput ERROR: output file is not open" << std::endl;
    return;
  }

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();

  std::cout << "O2FoCalGenQC::WriteOutput: wrote " << fOutputFilename.Data() << std::endl;

  delete fOutputFile;
  fOutputFile = nullptr;
  ResetHistPointers();
  fTracksByEvent.clear();
}
