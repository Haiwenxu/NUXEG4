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
/// \file optical/LXe/src/LXeDetectorMessenger.cc
/// \brief Implementation of the LXeDetectorMessenger class
//
//
#include "LXeDetectorMessenger.hh"

#include "LXeDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Scintillation.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorMessenger::LXeDetectorMessenger(LXeDetectorConstruction* detector)
  : fLXeDetector(detector)
{
  // Setup a command directory for detector controls with guidance
  fDetectorDir = new G4UIdirectory("/LXe/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fVolumesDir = new G4UIdirectory("/LXe/detector/volumes/");
  fVolumesDir->SetGuidance("Enable/disable volumes");

  fReflectivityCmd = new G4UIcmdWithADouble("/LXe/detector/reflectivity", this);
  fReflectivityCmd->SetGuidance("Set the reflectivity of the housing.");
  fReflectivityCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fReflectivityCmd->SetToBeBroadcasted(false);


  fLxeCmd = new G4UIcmdWithABool("/LXe/detector/volumes/lxe", this);
  fLxeCmd->SetGuidance("Enable/Disable the main detector volume.");
  fLxeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fLxeCmd->SetToBeBroadcasted(false);

  fMainScintYield =
    new G4UIcmdWithADouble("/LXe/detector/MainScintYield", this);
  fMainScintYield->SetGuidance("Set scinitillation yield of main volume.");
  fMainScintYield->SetGuidance("Specified in photons/MeV");
  fMainScintYield->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMainScintYield->SetToBeBroadcasted(false);


  fSaveThresholdCmd = new G4UIcmdWithAnInteger("/LXe/saveThreshold", this);
  fSaveThresholdCmd->SetGuidance(
    "Set the photon count threshold for saving the random number seed");
  fSaveThresholdCmd->SetParameterName("photons", true);
  fSaveThresholdCmd->SetDefaultValue(4500);
  fSaveThresholdCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDefaultsCmd = new G4UIcommand("/LXe/detector/defaults", this);
  fDefaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  fDefaultsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDefaultsCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorMessenger::~LXeDetectorMessenger()
{
  delete fLxeCmd;
  delete fReflectivityCmd;
  delete fMainScintYield;
  delete fSaveThresholdCmd;
  delete fDefaultsCmd;
  delete fDetectorDir;
  delete fVolumesDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

  if(command == fReflectivityCmd)
  {
    fLXeDetector->SetHousingReflectivity(
      fReflectivityCmd->GetNewDoubleValue(newValue));
  }
  // else if(command == fWlsCmd)
  // {
  //   fLXeDetector->SetWLSSlabOn(fWlsCmd->GetNewBoolValue(newValue));
  // }
  else if(command == fLxeCmd)
  {
    fLXeDetector->SetMainVolumeOn(fLxeCmd->GetNewBoolValue(newValue));
  }
  else if(command == fMainScintYield)
  {
    fLXeDetector->SetMainScintYield(
      fMainScintYield->GetNewDoubleValue(newValue));
  }
  else if(command == fSaveThresholdCmd)
  {
    fLXeDetector->SetSaveThreshold(fSaveThresholdCmd->GetNewIntValue(newValue));
  }
  else if(command == fDefaultsCmd)
  {
    fLXeDetector->SetDefaults();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}
