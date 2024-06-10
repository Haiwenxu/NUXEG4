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
/// \file optical/LXe/include/LXeDetectorConstruction.hh
/// \brief Definition of the LXeDetectorConstruction class
//
//
#ifndef LXeDetectorConstruction_h
#define LXeDetectorConstruction_h 1

#include "LXeDetectorMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MaterialTable.hh"
#include "G4Cache.hh"
#include "G4SystemOfUnits.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "globals.hh"
#include <vector>

class LXeMainVolume;
class LXeSiPMSD;
class LXeScintSD;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;
class G4SubtractionSolid;

class LXeDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  LXeDetectorConstruction();
  ~LXeDetectorConstruction();

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

  // Functions to modify the geometry
  void SetDefaults();
  void SetSaveThreshold(G4int);

  // Get values
  G4int GetSaveThreshold() const { return fSaveThreshold; };

	G4double GetHousingR() const { return fHousing_r; }
	G4double GetHousingZ() const { return fHousing_z; }
	G4double GetPolygonEdgeNum() const { return fPolygon_edge_num; }
	G4double GetPolygonEdgeLen() const { return fPolygon_edge_len; }
	G4double GetPolygonHeight() const { return fPolygon_height; }
	G4double GetAnodeRadius() const { return fAnode_radius; }
	G4double GetGateRadius() const { return fGate_radius; }
	G4double GetAnodeHeight() const { return fAnode_height; }
	G4double GetGateHeight() const { return fGate_height; }
  

  void SetHousingReflectivity(G4double);
  G4double GetHousingReflectivity() const { return fRefl; }

  void SetMainVolumeOn(G4bool b);
  G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

  void SetMainScintYield(G4double);




 private:
  void DefineMaterials();

  LXeDetectorMessenger* fDetectorMessenger;

  G4Box* fExperimentalHall_box;
  G4LogicalVolume* fExperimentalHall_log;
  G4VPhysicalVolume* fExperimentalHall_phys;

  // Materials & Elements
  G4Material* fLXe;
  G4Material* fAl;
  G4Element* fN;
  G4Element* fO;
  G4Material* fAir;
  G4Material* fVacuum;
  G4Element* fC;
  G4Element* fH;
  G4Material* fGlass;
  G4Material* fPstyrene;
  G4Material* fPMMA;
  G4Material* fPethylene1;
  G4Material* fPethylene2;

  G4Element* fF;
  G4Material* fPTFE;

  // Geometry
  // G4double fScint_x;
  // G4double fScint_y;
  // G4double fScint_z;

  G4double fHousing_r;
  G4double fHousing_z;

  G4double fPolygon_edge_num;
  G4double fPolygon_edge_len;
  G4double fPolygon_height;
  G4double fAnode_radius;
  G4double fGate_radius;
  G4double fAnode_height;
  G4double fGate_height;


  // G4double Cathode_width;
  // G4double Cathode_height;
  // G4double Cathode_thickness;
  // G4double Cathode_R;
  // G4double Cathode_VEdge;
  // G4double Cathode_HEdge;  
  
  // G4double Cathode_inner_cut_width;
  // G4double Cathode_inner_cut_thickness;
  // G4double Cathode_inner_cut_height;

  // G4double Cathode_HGrid_length;
  // G4double Cathode_HGrid_thickness;

  // G4double Cathode_VGrid_length;
  // G4double Cathode_VGrid_thickness;

  // G4double Cathode_VGrid_pitch;
  // G4double Cathode_HGrid_pitch;

  // G4int N_VGrid;
  // G4int N_HGrid;





  G4int fSaveThreshold;
  G4double fRefl;
  G4bool fMainVolumeOn;




  LXeMainVolume* fMainVolume;

  G4MaterialPropertiesTable* fLXe_mt;

  // Sensitive Detectors
  G4Cache<LXeScintSD*> fScint_SD;
  G4Cache<LXeSiPMSD*> fSipm_SD;
};

#endif
