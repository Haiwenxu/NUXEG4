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
/// \file NUXE/src/LXeSiPMSD.cc
/// \brief Implementation of the LXeSiPMSD class
//
//
#include "LXeSiPMSD.hh"

#include "LXeDetectorConstruction.hh"
#include "LXeSiPMHit.hh"
#include "LXeUserTrackInformation.hh"

#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSiPMSD::LXeSiPMSD(G4String name)
  : G4VSensitiveDetector(name)
  , fSiPMHitCollection(nullptr)
  , fSiPMPositionsX(nullptr)
  , fSiPMPositionsY(nullptr)
  , fSiPMPositionsZ(nullptr)
  , fHitCID(-1)
{
  collectionName.insert("sipmHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSiPMSD::~LXeSiPMSD()
{
  delete fSiPMPositionsX;
  delete fSiPMPositionsY;
  delete fSiPMPositionsZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSiPMSD::SetSipmPositions(const std::vector<G4ThreeVector>& positions)
{
  for(size_t i = 0; i < positions.size(); ++i)
  {
    if(fSiPMPositionsX)
      fSiPMPositionsX->push_back(positions[i].x());
    if(fSiPMPositionsY)
      fSiPMPositionsY->push_back(positions[i].y());
    if(fSiPMPositionsZ)
      fSiPMPositionsZ->push_back(positions[i].z());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSiPMSD::Initialize(G4HCofThisEvent* hitsCE)
{
  fSiPMHitCollection =
    new LXeSiPMHitsCollection(SensitiveDetectorName, collectionName[0]);

  if(fHitCID < 0)
  {
    fHitCID = G4SDManager::GetSDMpointer()->GetCollectionID(fSiPMHitCollection);
  }
  hitsCE->AddHitsCollection(fHitCID, fSiPMHitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool LXeSiPMSD::ProcessHits(G4Step*, G4TouchableHistory*) { return false; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Generates a hit and uses the postStepPoint's mother volume replica number
// PostStepPoint because the hit is generated manually when the photon is
// absorbed by the photocathode

G4bool LXeSiPMSD::ProcessHits_boundary(const G4Step* aStep, G4TouchableHistory*)
{
  // need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() !=
     G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  // User replica number 1 since photocathode is a daughter volume
  // to the sipm which was replicated
  G4int sipmNumber =
    aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  G4VPhysicalVolume* physVol =
    aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);


  // Find the correct hit collection
  size_t n       = fSiPMHitCollection->entries();
  LXeSiPMHit* hit = nullptr;





  //// Avoid different hits on the same sensitive detector(SiPM in this case):
  // for(size_t i = 0; i < n; ++i)
  // {
  //   if((*fSiPMHitCollection)[i]->GetSiPMNumber() == sipmNumber)
  //   {
  //     hit = (*fSiPMHitCollection)[i];
  //     break;
  //   }
  // }
  // if(hit == nullptr)
  // {                         // this SiPM wasn't previously hit in this event
  //   hit = new LXeSiPMHit();  // so create new hit
  //   hit->SetSiPMNumber(sipmNumber);
  //   hit->SetSiPMPhysVol(physVol);
  //   fSiPMHitCollection->insert(hit);
  //   hit->SetSiPMPos((*fSiPMPositionsX)[sipmNumber], (*fSiPMPositionsY)[sipmNumber],
  //                  (*fSiPMPositionsZ)[sipmNumber]);
  // }

  // Create a hit with every single detected optical photons, duplicated SD allowed:
  hit = new LXeSiPMHit();
  hit->SetSiPMNumber(sipmNumber);
  hit->SetSiPMPhysVol(physVol);
  fSiPMHitCollection->insert(hit);
  hit->SetSiPMPos((*fSiPMPositionsX)[sipmNumber], (*fSiPMPositionsY)[sipmNumber],
                (*fSiPMPositionsZ)[sipmNumber]);


  hit->IncPhotonCount();  // increment hit for the selected sipm
  hit->SetDrawit(true);

  return true;
}
