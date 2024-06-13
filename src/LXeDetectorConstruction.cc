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
/// \file NUXE/src/LXeDetectorConstruction.cc
/// \brief Implementation of the LXeDetectorConstruction class
//
//
#include "LXeDetectorConstruction.hh"

#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"
#include "LXeSiPMSD.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4OpticalSurface.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::LXeDetectorConstruction()
  : fLXe_mt(nullptr)
  // , fMPTPStyrene(nullptr)
{
  fExperimentalHall_box  = nullptr;
  fExperimentalHall_log  = nullptr;
  fExperimentalHall_phys = nullptr;

  fLXe = fAl = fAir = fVacuum = fGlass = nullptr;
  fPstyrene = fPMMA = fPethylene1 = fPethylene2 = nullptr;

  fN = fO = fC = fH = fF = nullptr;

  fSaveThreshold = 0;
  SetDefaults();

  DefineMaterials();
  fDetectorMessenger = new LXeDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::~LXeDetectorConstruction()
{
  if(fMainVolume)
  {
    delete fMainVolume;
  }
  delete fLXe_mt;
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::DefineMaterials()
{
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyPMMA = 1;
  G4int nC_PMMA  = 3 + 2 * polyPMMA;
  G4int nH_PMMA  = 6 + 2 * polyPMMA;

  G4int polyeth = 1;
  G4int nC_eth  = 2 * polyeth;
  G4int nH_eth  = 4 * polyeth;

  //PTFE
  G4int ptfe = 1;
  G4int nC_ptfe = 2 * ptfe;
  G4int nF_ptfe = 4 * ptfe;

  //***Elements
  fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
  fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
  fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
  fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);
  fF = new G4Element("F", "F", z = 9., a = 19.00 * g / mole);

  //***Materials
  // Liquid Xenon
  fLXe = new G4Material("LXe", z = 54., a = 131.29 * g / mole,
                        density = 3.020 * g / cm3);
  // Aluminum
  fAl = new G4Material("Al", z = 13., a = 26.98 * g / mole,
                       density = 2.7 * g / cm3);
  // Vacuum
  fVacuum = new G4Material("Vacuum", z = 1., a = 1.01 * g / mole,
                           density = universe_mean_density, kStateGas,
                           0.1 * kelvin, 1.e-19 * pascal);
  // Air
  fAir = new G4Material("Air", density = 1.29 * mg / cm3, 2);
  fAir->AddElement(fN, 70 * perCent);
  fAir->AddElement(fO, 30 * perCent);
  // Glass
  fGlass = new G4Material("Glass", density = 1.032 * g / cm3, 2);
  fGlass->AddElement(fC, 91.533 * perCent);
  fGlass->AddElement(fH, 8.467 * perCent);


  // PTFE
  fPTFE = new G4Material("PTFE", density = 2200. * kg / m3, 2);
  fPTFE->AddElement(fC, nC_ptfe);
  fPTFE->AddElement(fF, nF_ptfe);

  //***Material properties tables

  std::vector<G4double> lxe_Energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };

  std::vector<G4double> lxe_SCINT = { 0.1, 1.0, 0.1 };
  std::vector<G4double> lxe_RIND  = { 1.59, 1.57, 1.54 };
  std::vector<G4double> lxe_ABSL  = { 35. * cm, 35. * cm, 35. * cm };
  fLXe_mt = new G4MaterialPropertiesTable();
  fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT1", lxe_Energy, lxe_SCINT);
  fLXe_mt->AddProperty("SCINTILLATIONCOMPONENT2", lxe_Energy, lxe_SCINT);
  fLXe_mt->AddProperty("RINDEX", lxe_Energy, lxe_RIND);
  fLXe_mt->AddProperty("ABSLENGTH", lxe_Energy, lxe_ABSL);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
  fLXe_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
  fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
  fLXe_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
  fLXe->SetMaterialPropertiesTable(fLXe_mt);

  // Set the Birks Constant for the LXe scintillator
  fLXe->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  // PTFE optical properties
  std::vector<G4double> PTFE_energy = { 7.0 * eV, 7.07 * eV, 7.14 * eV };
  std::vector<G4double> PTFE_AbsLength = { 0.1 * cm, 0.1 * cm, 0.1 * cm, };
  // std::vector<G4double> PTFE_reflectivity = { 0.9, 0.9, 0.9 };
  std::vector<G4double> PTFE_refractivity = { 1.51, 1.54, 1.58 };
  // std::vector<G4double> PTFE_SpecularLobe = {0.25, 0.25, 0.25 };
  // std::vector<G4double> PTFE_SpecularSpike = { 0.01, 0.01, 0.01 };
  // std::vector<G4double> PTFE_DiffuseLobe = { 0.55, 0.55, 0.55 };
  // std::vector<G4double> PTFE_BackScatter = { 0.01, 0.01, 0.01 };
  // std::vector<G4double> PTFE_efficiency = { 0.0, 0.0, 0.0 };
  G4MaterialPropertiesTable* PTFE_mt = new G4MaterialPropertiesTable();
  PTFE_mt->AddProperty("ABSLENGTH", PTFE_energy, PTFE_AbsLength);
  printf("Set PTFE optical properties.\n");
  PTFE_mt->AddProperty("RINDEX", PTFE_energy, PTFE_refractivity);
  // PTFE_mt->AddProperty("REFLECTIVITY", PTFE_energy, PTFE_reflectivity);
  // PTFE_mt->AddProperty("SPECULARLOBECONSTANT", PTFE_energy, PTFE_SpecularLobe);
  // PTFE_mt->AddProperty("SPECULARSPIKECONSTANT", PTFE_energy, PTFE_SpecularSpike);
  // PTFE_mt->AddProperty("BACKSCATTERCONSTANT", PTFE_energy, PTFE_BackScatter);
  // PTFE_mt->AddProperty("EFFICIENCY", PTFE_energy, PTFE_efficiency);
  fPTFE->SetMaterialPropertiesTable(PTFE_mt);



  std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm, 420. * cm };
  std::vector<G4double> glass_refractivity = { 1.60, 1.605, 1.61 };
  G4MaterialPropertiesTable* glass_mt   = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH", lxe_Energy, glass_AbsLength);
  glass_mt->AddProperty("RINDEX", lxe_Energy,glass_refractivity);
  fGlass->SetMaterialPropertiesTable(glass_mt);

  G4MaterialPropertiesTable* vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", "Air");
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);  // Give air the same rindex

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::Construct()
{
  // The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = 1.1 * m;
  G4double expHall_y = 1.1 * m;
  G4double expHall_z = 1.1 * m;

  // Create experimental hall
  fExperimentalHall_box =
    new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
  fExperimentalHall_log =
    new G4LogicalVolume(fExperimentalHall_box, fVacuum, "expHall_log", 0, 0, 0);
  fExperimentalHall_phys = new G4PVPlacement(
    0, G4ThreeVector(), fExperimentalHall_log, "expHall", 0, false, 0);

  fExperimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Place the main volume
  if(fMainVolumeOn)
  {
    fMainVolume = new LXeMainVolume(0, G4ThreeVector(), fExperimentalHall_log,
                                    false, 0, this);
  }

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::ConstructSDandField()
{
  if(!fMainVolume)
    return;

  // SiPM SD

  LXeSiPMSD* sipm = fSipm_SD.Get();
  if(!sipm)
  {
    // Created here so it exists as SiPMs are being placed
    G4cout << "Construction /LXeDet/sipmSD" << G4endl;
    LXeSiPMSD* sipm_SD = new LXeSiPMSD("/LXeDet/sipmSD");
    fSipm_SD.Put(sipm_SD);

    sipm_SD->InitSiPMs();
    sipm_SD->SetSipmPositions(fMainVolume->GetSipmPositions());
  }
  else
  {
    sipm->InitSiPMs();
    sipm->SetSipmPositions(fMainVolume->GetSipmPositions());
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fSipm_SD.Get());
  // sensitive detector is not actually on the photocathode.
  // processHits gets done manually by the stepping action.
  // It is used to detect when photons hit and get absorbed & detected at the
  // boundary to the photocathode (which doesn't get done by attaching it to a
  // logical volume.
  // It does however need to be attached to something or else it doesn't get
  // reset at the begining of events

  SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fSipm_SD.Get());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDefaults()
{
  // Resets to default values
  fPolygon_edge_num = 8;
  fPolygon_edge_len = 105.*sin(3.1415926/fPolygon_edge_num) * mm;
  fPolygon_height = 100./2. * mm;
  fHousing_r = (2 + fPolygon_edge_len / (2*sin(3.1415926/fPolygon_edge_num))) * mm;
  fHousing_z = 120./2. * mm;
  
  fAnode_radius = 9e-3 * mm;
  fGate_radius = 5e-2 * mm;
  fAnode_height = 100./2. * mm;
  fGate_height = 100./2. * mm;
  
  fRefl     = 1.0;

  fMainVolumeOn = true;
  fMainVolume   = nullptr;

  G4UImanager::GetUIpointer()->ApplyCommand(
    "/LXe/detector/scintYieldFactor 1.");

  if(fLXe_mt)
    fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", 12000. / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingReflectivity(G4double r)
{
  fRefl = r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainVolumeOn(G4bool b)
{
  fMainVolumeOn = b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainScintYield(G4double y)
{
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD", y / MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetSaveThreshold(G4int save)
{
  // Sets the save threshold for the random number seed. If the number of
  // photons generated in an event is lower than this, then save the seed for
  // this event in a file called run###evt###.rndm

  fSaveThreshold = save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
