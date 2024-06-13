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
/// \file NUXE/src/LXePrimaryGeneratorAction.cc
/// \brief Implementation of the LXePrimaryGeneratorAction class
//
//
#include "LXePrimaryGeneratorAction.hh"
#include "LXeHistoManager.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::LXePrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun     = new G4ParticleGun(n_particle);
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // G4String particleName;

  //solid angle
  //
  // G4double alphaMin =  0*rad;      //alpha in [0,pi]
  // G4double alphaMax = pi*rad;
  // G4double fCosAlphaMin = std::cos(alphaMin);
  // G4double fCosAlphaMax = std::cos(alphaMax);
  
  // G4double fPsiMin = 0*rad;       //psi in [0, 2*pi]
  // G4double fPsiMax = twopi*rad;

  // G4double cosAlpha = fCosAlphaMin-G4UniformRand()*(fCosAlphaMin-fCosAlphaMax);
  // G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  // G4double psi = fPsiMin + G4UniformRand()*(fPsiMax - fPsiMin);

  // G4double ux = sinAlpha*std::cos(psi),
  //          uy = sinAlpha*std::sin(psi),
  //          uz = cosAlpha;


  // const G4double r = 20.0 * (G4UniformRand()) * mm;
  // const G4double zmax = 4.5 * cm;   
  
  // //vertex 1 uniform on cylinder
  // G4double alpha = twopi*(G4UniformRand());  //alpha uniform in (0, 2*pi)
  // // G4double alpha = twopi*G4UniformRand();  //alpha uniform in (0, 2*pi)
  // G4double mux = std::cos(alpha);
  // G4double muy = std::sin(alpha);
  // G4double z = zmax*(2*G4UniformRand() - 1);  //z uniform in (-zmax, +zmax)

  // fParticleGun->SetParticleDefinition(
  //   particleTable->FindParticle(particleName = "opticalphoton"));
  // // Default energy,position,momentum
  // m_dKineticEnergy = 7.0 * eV;
  // // m_dKineticEnergy = 511 * keV;
  // fParticleGun->SetParticleEnergy(m_dKineticEnergy);




  // fParticleGun->SetParticlePosition(G4ThreeVector(r*mux, r*muy, z));
  // // fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  // // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(2*G4UniformRand()-1, 2*G4UniformRand()-1, 2*G4UniformRand()-1));
  // fParticleGun->SetParticlePolarization(G4ThreeVector(G4UniformRand(),G4UniformRand(), G4UniformRand()));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  fParticleGun->SetParticleDefinition(
    particleTable->FindParticle(particleName = "opticalphoton"));
  // Default energy,position,momentum
  m_dKineticEnergy = 7.0 * eV;
  // m_dKineticEnergy = 511 * keV;
  fParticleGun->SetParticleEnergy(m_dKineticEnergy);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::~LXePrimaryGeneratorAction() { delete fParticleGun; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//NOTE: This function denotes the beginning of each event.
void LXePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  CLHEP::HepRandom::setTheSeed((unsigned)clock());



  const G4double r = 36.5 * (G4UniformRand()) * mm;
  const G4double rmin = 50.0 * mm;
  const G4double zmax = 4.2 * cm;   
  
  //vertex 1 uniform on cylinder
  G4double alpha = twopi*(G4UniformRand());  //alpha uniform in (0, 2*pi)
  // G4double alpha = twopi*G4UniformRand();  //alpha uniform in (0, 2*pi)
  G4double mux = std::cos(alpha);
  G4double muy = std::sin(alpha);
  G4double z = zmax*(2*G4UniformRand() - 1);  //z uniform in (-zmax, +zmax)





  fParticleGun->SetParticlePosition(G4ThreeVector(r*mux, r*muy, z));
  // fParticleGun->SetParticlePosition(G4ThreeVector(r*mux+rmin, r*muy+rmin, z));



  G4ThreeVector Pder = G4ThreeVector(2*G4UniformRand()-1, 2*G4UniformRand()-1, 2*G4UniformRand()-1);
  // G4ThreeVector Pder = G4ThreeVector(-(r*mux+rmin), -(r*muy+rmin), -z);
  fParticleGun->SetParticleMomentumDirection(Pder);
  fParticleGun->SetParticlePolarization(Pder.orthogonal());
  fParticleGun->GeneratePrimaryVertex(anEvent);

}
