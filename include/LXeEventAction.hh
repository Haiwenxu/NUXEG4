//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file NUXE/include/LXeEven tAction.hh
/// \brief Definition of the LXeEventAction class
//

#ifndef LXeEventAction_h
#define LXeEventAction_h 1

#include "LXeEventMessenger.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"

class G4Event;
class LXeDetectorConstruction;

class LXeEventAction : public G4UserEventAction
{
 public:
  LXeEventAction(const LXeDetectorConstruction*);
  ~LXeEventAction();

 public:
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void SetEventVerbose(G4int v) { fVerbose = v; }

  void SetSiPMThreshold(G4int t) { fSiPMThreshold = t; }

  void SetForceDrawPhotons(G4bool b) { fForcedrawphotons = b; }
  void SetForceDrawNoPhotons(G4bool b) { fForcenophotons = b; }

  void IncPhotonCount_Ceren() { ++fPhotonCount_Ceren; }
  void IncEDep(G4double dep) { fTotE += dep; }
  void IncAbsorption() { ++fAbsorptionCount; }
  void IncBoundaryAbsorption() { ++fBoundaryAbsorptionCount; }
  void IncHitCount(G4int i = 1) { fHitCount += i; }

  void SetEWeightPos(const G4ThreeVector& p) { fEWeightPos = p; }
  // void SetReconPos(const G4ThreeVector& p) { fReconPos = p; }
  void SetConvPos(const G4ThreeVector& p)
  {
    fConvPos    = p;
    fConvPosSet = true;
  }
  void SetPosMax(const G4ThreeVector& p, G4double edep)
  {
    fPosMax  = p;
    fEdepMax = edep;
  }

  G4int GetPhotonCount_Ceren() const { return fPhotonCount_Ceren; }
  G4int GetHitCount() const { return fHitCount; }
  G4double GetEDep() const { return fTotE; }
  G4int GetAbsorptionCount() const { return fAbsorptionCount; }
  G4int GetBoundaryAbsorptionCount() const { return fBoundaryAbsorptionCount; }

  G4ThreeVector GetEWeightPos() { return fEWeightPos; }
  // G4ThreeVector GetReconPos() { return fReconPos; }
  G4ThreeVector GetConvPos() { return fConvPos; }
  G4ThreeVector GetPosMax() { return fPosMax; }
  G4double GetEDepMax() { return fEdepMax; }
  G4double IsConvPosSet() { return fConvPosSet; }

  // Gets the total photon count produced

  void IncSiPMSAboveThreshold() { ++fSiPMsAboveThreshold; }
  G4int GetSiPMSAboveThreshold() { return fSiPMsAboveThreshold; }


  void push_back_TOF( G4double tof) { Event_Hits_TOF.push_back(tof); }
  std::vector<G4double>& Get_TOF_vector() { return Event_Hits_TOF; }

  void push_back_Hit_X( G4double hitX ) { Event_Hits_X.push_back(hitX); }
  std::vector<G4double>& Get_Event_Hits_Xvector() { return Event_Hits_X; }

  void push_back_Hit_Y( G4double hitY ) { Event_Hits_Y.push_back(hitY); }
  std::vector<G4double>& Get_Event_Hits_Yvector() { return Event_Hits_Y; }

  void push_back_Hit_Z( G4double hitZ ) { Event_Hits_Z.push_back(hitZ); }
  std::vector<G4double>& Get_Event_Hits_Zvector() { return Event_Hits_Z; }

  void push_back_Hit_sipmN( G4int n ) { sipmN.push_back(n); }
  std::vector<G4int>& Get_Event_Hits_sipmN() { return sipmN; }










  // std::vector<G4double>& GetHitTime() { return HitTime; }

 private:
  LXeEventMessenger* fEventMessenger;
  const LXeDetectorConstruction* fDetector;

  G4int fSiPMCollID;

  G4int fVerbose;

  G4int fSiPMThreshold;

  G4bool fForcedrawphotons;
  G4bool fForcenophotons;

  G4int fHitCount;
  G4int fPhotonCount_Ceren;
  G4int fAbsorptionCount;
  G4int fBoundaryAbsorptionCount;

  G4double fTotE;

  // These only have meaning if totE > 0
  // If totE = 0 then these won't be set by EndOfEventAction
  G4ThreeVector fEWeightPos;
  // G4ThreeVector fReconPos;  // Also relies on hitCount>0
  G4ThreeVector fConvPos;   // true (initial) converstion position
  G4bool fConvPosSet;
  G4ThreeVector fPosMax;
  G4double fEdepMax;

  G4int fSiPMsAboveThreshold;

  std::vector<G4double> Event_Hits_TOF;

  std::vector<G4double> Event_Hits_X;
  std::vector<G4double> Event_Hits_Y;
  std::vector<G4double> Event_Hits_Z;


  std::vector<G4int> sipmN;



};

#endif
