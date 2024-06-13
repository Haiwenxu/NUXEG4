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
/// \file NUXE/include/LXeSiPMSD.hh
/// \brief Definition of the LXeSiPMSD class
//
//
#ifndef LXeSiPMSD_h
#define LXeSiPMSD_h 1

#include "LXeSiPMHit.hh"

#include "G4VSensitiveDetector.hh"

#include <vector>

class G4DataVector;
class G4HCofThisEvent;
class G4Step;

class LXeSiPMSD : public G4VSensitiveDetector
{
 public:
  LXeSiPMSD(G4String name);
  ~LXeSiPMSD();

  void Initialize(G4HCofThisEvent*) override;
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;

  // A version of processHits active on boundary
  G4bool ProcessHits_boundary(const G4Step*, G4TouchableHistory*);

  // Initialize the arrays to store SiPM possitions
  inline void InitSiPMs()
  {
    if(fSiPMPositionsX)
      delete fSiPMPositionsX;
    if(fSiPMPositionsY)
      delete fSiPMPositionsY;
    if(fSiPMPositionsZ)
      delete fSiPMPositionsZ;
    fSiPMPositionsX = new G4DataVector();
    fSiPMPositionsY = new G4DataVector();
    fSiPMPositionsZ = new G4DataVector();
  }

  // Store a SiPM position
  void SetSipmPositions(const std::vector<G4ThreeVector>& positions);

 private:
  LXeSiPMHitsCollection* fSiPMHitCollection;

  G4DataVector* fSiPMPositionsX;
  G4DataVector* fSiPMPositionsY;
  G4DataVector* fSiPMPositionsZ;

  G4int fHitCID;
};

#endif
