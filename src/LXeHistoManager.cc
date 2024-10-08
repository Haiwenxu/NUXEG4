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
/// \file NUXE/src/LXeHistoManager.cc
/// \brief Implementation of the LXeHistoManager class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LXeHistoManager.hh"


#include "G4RunManager.hh"
#include "LXeEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::LXeHistoManager()
  : fFileName("../Data/lxe")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeHistoManager::~LXeHistoManager() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeHistoManager::Book()
{

  // LXeEventAction* eventAction = static_cast<LXeEventAction*>(
  //       G4RunManager::GetRunManager()->GetUserEventAction());

  LXeEventAction* eventAction = const_cast<LXeEventAction*>(
    static_cast<const LXeEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction())
);



  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  // enable inactivation of histograms
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName(fFileName);


  // Define histogram indices, titles

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins   = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;
  G4double Tmin = 0.;
  G4double Tmax = 300.;

  // 0
  analysisManager->CreateH1("0", "dummy", nbins, vmin, vmax);
  // 1
  analysisManager->CreateH1("hits per event", "hits per event", nbins, vmin,
                            vmax);
  // 2
  analysisManager->CreateH1("hits above threshold",
                            "hits per event above threshold", nbins, vmin,
                            vmax);
  // 3
  analysisManager->CreateH1("scintillation", "scintillation photons per event",
                            nbins, vmin, vmax);
  // 4
  analysisManager->CreateH1("Cerenkov", "Cerenkov photons per event", nbins,
                            vmin, vmax);
  // 5
  analysisManager->CreateH1("absorbed", "absorbed photons per event", nbins,
                            vmin, vmax);
  // 6
  analysisManager->CreateH1("boundary absorbed",
                            "photons absorbed at boundary per event", nbins,
                            vmin, vmax);
  // 7
  analysisManager->CreateH1(
    "E dep", "energy deposition in scintillator per event", nbins, vmin, vmax);

  // 8
  analysisManager->CreateH1(
    "SiPM Hits", "SiPM# Hit Pattern", nbins, vmin, vmax);
    
  // 9
  analysisManager->CreateH1(
    "Event Time", "Time(ns)", nbins, Tmin, Tmax);

  G4int nXbins = 200;
  G4int nYbins = 200;
  G4int nZbins = 200;
  G4double Xmin = -50;
  G4double Ymin = -50;
  G4double Zmin = -60;
  G4double Xmax = 50;
  G4double Ymax = 50;
  G4double Zmax = 60;
  const G4String XunitName = "mm";
  const G4String YunitName = "mm";
  const G4String ZunitName = "mm";
  const G4String XfunctionName = "none";
  const G4String YfunctionName = "none";
  const G4String ZfunctionName = "none";
  const G4String XbinSchemeName = "linear";
  const G4String YbinSchemeName = "linear";
  const G4String ZbinSchemeName = "linear";
  // 9
  analysisManager->CreateH3("Position Reconstruction", "PosReconstruction",
                            nXbins, Xmin, Xmax,
                            nYbins, Ymin, Ymax,
                            nZbins, Zmin, Zmax,
                            XunitName, YunitName, ZunitName,
                            XfunctionName, YfunctionName, ZfunctionName,
                            XbinSchemeName, YbinSchemeName, ZbinSchemeName
  );



  analysisManager->CreateNtuple("NUXE", "Event_Hit");
  analysisManager->CreateNtupleDColumn("Event_Energy"); // column Id = 0
  analysisManager->CreateNtupleDColumn("Event_X"); // column Id = 1
  analysisManager->CreateNtupleDColumn("Event_Y"); // column Id = 2
  analysisManager->CreateNtupleDColumn("Event_Z"); // column Id = 3
  analysisManager->CreateNtupleDColumn("Hit_Time",eventAction->Get_TOF_vector()); // column Id = 4
  analysisManager->CreateNtupleDColumn("Hit_X",eventAction->Get_Event_Hits_Xvector()); // column Id = 5
  analysisManager->CreateNtupleDColumn("Hit_Y",eventAction->Get_Event_Hits_Yvector()); // column Id = 6
  analysisManager->CreateNtupleDColumn("Hit_Z",eventAction->Get_Event_Hits_Zvector()); // column Id = 7
  analysisManager->CreateNtupleIColumn("SiPM#",eventAction->Get_Event_Hits_sipmN());    // column Id = 8
  
  
  // analysisManager->SetNtupleMerging(true);
  analysisManager->FinishNtuple(); 

  analysisManager->SetNtupleFileName(0, "../Data/LXe_ntuple");


  // Create all histograms as activated
  for(G4int i = 0; i < analysisManager->GetNofH1s(); ++i)
  {
    analysisManager->SetH1Activation(i, true);
  }

  analysisManager->SetH3Activation(0, true);
}
