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
/// \file NUXE/src/LXeEventAction.cc
/// \brief Implementation of the LXeEventAction class
//
//
#include "LXeEventAction.hh"

#include "LXeDetectorConstruction.hh"
#include "LXeHistoManager.hh"
#include "LXeSiPMHit.hh"
#include "LXeRun.hh"
#include "LXeTrajectory.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::LXeEventAction(const LXeDetectorConstruction* det)
  : fDetector(det)
  , fSiPMCollID(-1)
  , fVerbose(1)
  , fSiPMThreshold(1)
  , fForcedrawphotons(false)
  , fForcenophotons(false)
{
  fEventMessenger = new LXeEventMessenger(this);

  fHitCount                = 0;
  fPhotonCount_Ceren       = 0;
  fAbsorptionCount         = 0;
  fBoundaryAbsorptionCount = 0;
  fTotE                    = 0.0;

  fConvPosSet = false;
  fEdepMax    = 0.0;

  fSiPMsAboveThreshold = 0;

  // HitTime                  = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::~LXeEventAction() {
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::BeginOfEventAction(const G4Event*)
{
  fHitCount                = 0;
  fPhotonCount_Ceren       = 0;
  fAbsorptionCount         = 0;
  fBoundaryAbsorptionCount = 0;
  fTotE                    = 0.0;
  // Event_Hits_TOF           = ;


  fConvPosSet = false;
  fEdepMax    = 0.0;

  

  fSiPMsAboveThreshold = 0;

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fSiPMCollID < 0)
    fSiPMCollID = SDman->GetCollectionID("sipmHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4TrajectoryContainer* trajectoryContainer =
    anEvent->GetTrajectoryContainer();

  G4int n_trajectories = 0;
  if(trajectoryContainer)
    n_trajectories = trajectoryContainer->entries();
  // G4cout << "------------------------ number of trajectories of this event: " << n_trajectories << G4endl;
  // G4cout << "------------------------ trajectory vector of this event: " << trajectoryContainer->GetVector() << G4endl;

  // extract the trajectories and draw them
  if(G4VVisManager::GetConcreteInstance())
  {
    for(G4int i = 0; i < n_trajectories; ++i)
    {
      LXeTrajectory* trj =
        (LXeTrajectory*) ((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName() == "opticalphoton")
      {
        trj->SetForceDrawTrajectory(fForcedrawphotons);
        trj->SetForceNoDrawTrajectory(fForcenophotons);
      }
      trj->DrawTrajectory();
    }
  }

  LXeSiPMHitsCollection* sipmHC     = nullptr;
  G4HCofThisEvent* hitsCE         = anEvent->GetHCofThisEvent();


  // 0:   Primary Energy
  G4AnalysisManager::Instance()->FillNtupleDColumn(0,anEvent->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy());
  // 1-3: Primary Position(x,y,z)
  G4AnalysisManager::Instance()->FillNtupleDColumn(1,anEvent->GetPrimaryVertex()->GetPosition().x());
  G4AnalysisManager::Instance()->FillNtupleDColumn(2,anEvent->GetPrimaryVertex()->GetPosition().y());
  G4AnalysisManager::Instance()->FillNtupleDColumn(3,anEvent->GetPrimaryVertex()->GetPosition().z());
  // 4:   Secondaries(Hits) Time
  // G4AnalysisManager::Instance()->FillNtupleDColumn(4,);
  //5-7:  Hits Position(x,y,z)
  // G4AnalysisManager::Instance()->FillNtupleDColumn(5,);
  // G4AnalysisManager::Instance()->FillNtupleDColumn(6,);
  // G4AnalysisManager::Instance()->FillNtupleDColumn(7,);
  

  



  // Get the hit collections
  if(hitsCE)
  {
    if(fSiPMCollID >= 0)
    {
      sipmHC = (LXeSiPMHitsCollection*) (hitsCE->GetHC(fSiPMCollID));
    }
  }



  if(sipmHC)
  {
    std::vector<G4ThreeVector> reconPos;
    size_t N_Hits = sipmHC->entries();
    // std::cout << "=============== SiPM# " << SiPMs << " ================" << std::endl;

    // Gather info from all SiPMs
    for(size_t i = 0; i < N_Hits; ++i)
    {
      fHitCount += (*sipmHC)[i]->GetPhotonCount();
      reconPos.push_back((*sipmHC)[i]->GetSiPMPos());
      // reconPos += (*sipmHC)[i]->GetSiPMPos() * (*sipmHC)[i]->GetPhotonCount();

      if((*sipmHC)[i]->GetPhotonCount() >= fSiPMThreshold)
      {
        G4AnalysisManager::Instance()->FillH1(8, (*sipmHC)[i]->GetSiPMNumber());
        push_back_Hit_sipmN((*sipmHC)[i]->GetSiPMNumber());
        push_back_Hit_X((*sipmHC)[i]->GetSiPMPos().x());
        push_back_Hit_Y((*sipmHC)[i]->GetSiPMPos().y());
        push_back_Hit_Z((*sipmHC)[i]->GetSiPMPos().z());


        // G4AnalysisManager::Instance()->FillH3(0,(*sipmHC)[i]->GetSiPMPos().x(),
        //                                   (*sipmHC)[i]->GetSiPMPos().y(),
        //                                   (*sipmHC)[i]->GetSiPMPos().z());

        ++fSiPMsAboveThreshold;
      }
      else
      {  // wasn't above the threshold, turn it back off
        (*sipmHC)[i]->SetDrawit(false);
      }
    }

    G4AnalysisManager::Instance()->FillH1(1, fHitCount);
    G4AnalysisManager::Instance()->FillH1(2, fSiPMsAboveThreshold);

    if(fHitCount > 0)
    {  // don't bother unless there were hits
      // reconPos /= fHitCount;

      // if(fVerbose > 0)
      // {
      //   // G4cout << "\tReconstructed position of hits in LXe : " << reconPos / mm
      //   //        << G4endl;
      // }
      // for(size_t i = 0; i < reconPos.size(); ++i){
      //   G4AnalysisManager::Instance()->FillH3(0, reconPos[i].x(), reconPos[i].y(), reconPos[i].z());
      // }
      // fReconPos = reconPos;
    }
    sipmHC->DrawAllHits();
  }

  G4AnalysisManager::Instance()->AddNtupleRow();


  G4AnalysisManager::Instance()->FillH1(4, fPhotonCount_Ceren);
  G4AnalysisManager::Instance()->FillH1(5, fAbsorptionCount);
  G4AnalysisManager::Instance()->FillH1(6, fBoundaryAbsorptionCount);

  if(fVerbose > 0)
  {
    // G4cout << "*************************** New Event ***************************" << G4endl;

    // End of event output. later to be controlled by a verbose level
    G4cout << "\tNumber of photons that hit SiPMs in this event : " << fHitCount
           << G4endl;
    G4cout << "\tNumber of SiPMs above threshold(" << fSiPMThreshold
           << ") : " << fSiPMsAboveThreshold << G4endl;
    G4cout << "\tNumber of photons produced by cerenkov in this event : "
           << fPhotonCount_Ceren << G4endl;
    G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
           << fAbsorptionCount << G4endl;
    G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
           << "this event : " << fBoundaryAbsorptionCount << G4endl;
    G4cout << "Unaccounted for photons in this event : "
           << (fPhotonCount_Ceren - fAbsorptionCount -
               fHitCount - fBoundaryAbsorptionCount)
           << G4endl;
  }

  // update the run statistics
  LXeRun* run = static_cast<LXeRun*>(
    G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->IncHitCount(fHitCount);
  run->IncPhotonCount_Ceren(fPhotonCount_Ceren);
  run->IncEDep(fTotE);
  run->IncAbsorption(fAbsorptionCount);
  run->IncBoundaryAbsorption(fBoundaryAbsorptionCount);
  run->IncHitsAboveThreshold(fSiPMsAboveThreshold);

  // If we have set the flag to save 'special' events, save here
  if(fPhotonCount_Ceren < fDetector->GetSaveThreshold())
  {
    G4RunManager::GetRunManager()->rndmSaveThisEvent();
  }



  if (Event_Hits_TOF.size()){
    Event_Hits_TOF.clear();
  }

  if (Event_Hits_X.size()){
    Event_Hits_X.clear();
  }
  if (Event_Hits_Y.size()){
    Event_Hits_Y.clear();
  }
  if (Event_Hits_Z.size()){
    Event_Hits_Z.clear();
  }
  if(sipmN.size()){
    sipmN.clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
