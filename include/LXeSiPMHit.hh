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
/// \file optical/LXe/include/LXeSiPMHit.hh
/// \brief Definition of the LXeSiPMHit class
//
//
#ifndef LXeSiPMHit_h
#define LXeSiPMHit_h 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4VPhysicalVolume.hh"

class LXeSiPMHit : public G4VHit
{
 public:
  LXeSiPMHit();
  LXeSiPMHit(const LXeSiPMHit& right);
  ~LXeSiPMHit();

  const LXeSiPMHit& operator=(const LXeSiPMHit& right);
  G4bool operator==(const LXeSiPMHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);

  virtual void Draw();
  virtual void Print();

  inline void SetDrawit(G4bool b) { fDrawit = b; }
  inline G4bool GetDrawit() { return fDrawit; }

  inline void IncPhotonCount() { ++fPhotons; }
  inline G4int GetPhotonCount() { return fPhotons; }

  inline void SetSiPMNumber(G4int n) { fSipmNumber = n; }
  inline G4int GetSiPMNumber() { return fSipmNumber; }

  inline void SetSiPMPhysVol(G4VPhysicalVolume* physVol)
  {
    this->fPhysVol = physVol;
  }
  inline G4VPhysicalVolume* GetSiPMPhysVol() { return fPhysVol; }

  inline void SetSiPMPos(G4double x, G4double y, G4double z)
  {
    fPos = G4ThreeVector(x, y, z);
  }

  inline G4ThreeVector GetSiPMPos() { return fPos; }

 private:
  G4int fSipmNumber;
  G4int fPhotons;
  G4ThreeVector fPos;
  G4VPhysicalVolume* fPhysVol;
  G4bool fDrawit;
};

typedef G4THitsCollection<LXeSiPMHit> LXeSiPMHitsCollection;

extern G4ThreadLocal G4Allocator<LXeSiPMHit>* LXeSiPMHitAllocator;

inline void* LXeSiPMHit::operator new(size_t)
{
  if(!LXeSiPMHitAllocator)
    LXeSiPMHitAllocator = new G4Allocator<LXeSiPMHit>;
  return (void*) LXeSiPMHitAllocator->MallocSingle();
}

inline void LXeSiPMHit::operator delete(void* aHit)
{
  LXeSiPMHitAllocator->FreeSingle((LXeSiPMHit*) aHit);
}

#endif
