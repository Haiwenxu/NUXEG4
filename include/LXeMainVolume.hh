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
/// \file optical/LXe/include/LXeMainVolume.hh
/// \brief Definition of the LXeMainVolume class
//
#ifndef LXeMainVolume_h
#define LXeMainVolume_h 1

#include "LXeDetectorConstruction.hh"
#include "G4PVPlacement.hh"
class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;
class G4ExtrudedSolid;
class G4SubtractionSolid;
class G4IntersectionSolid;
class G4UnionSolid;
class LXeMainVolume : public G4PVPlacement
{
 public:
  LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                LXeDetectorConstruction* c);

  G4LogicalVolume* GetLogPhotoCath() { return fPhotocath_log; }
  G4LogicalVolume* GetLogScint() { return fScint_log; }

  std::vector<G4ThreeVector> GetSipmPositions() { return fSipmPositions; }

 private:
  void VisAttributes();
  void SurfaceProperties();

  void PlaceSiPMs(G4LogicalVolume* sipm_Log, G4RotationMatrix* rot, G4double& a,
                 G4double& b, G4double da, G4double db, G4double amin,
                 G4double bmin, G4int na, G4int nb, G4double& x, G4double& y,
                 G4double& z, G4int& k);

  void CopyValues();

  LXeDetectorConstruction* fConstructor;

  G4double fRefl;
  G4double fHousing_r;
  G4double fHousing_z;
  G4int fPolygon_edge_num;
  G4double fPolygon_edge_len;
  G4double fPolygon_height;
  G4double fAnode_radius;
  G4double fGate_radius;
  G4double fAnode_height;
  G4double fGate_height;


  // Basic Volumes
  //
  G4Tubs* fHousing_cylinder;
  G4Box* fScint_box;
  G4Box* fHousing_box;
  G4ExtrudedSolid* fScint_polygon;
  G4Box* fSiPM;
  G4Box* fPhotocath;
  G4Tubs* fAnode_wire;
  G4Tubs* fGate_wire;
  G4Box* fCathode_box;
  G4Box* fCathode_cut;
  G4SubtractionSolid* fCathode_frame;
  G4Box* fCathode_HGrid;
  G4Box* fCathode_VGrid;
  G4UnionSolid* fCathode_Mesh;

  // Logical volumes
  //
  G4LogicalVolume* fScint_log;
  G4LogicalVolume* fHousing_log;
  G4LogicalVolume* fSiPM_log;
  G4LogicalVolume* fPhotocath_log;
  G4LogicalVolume* fAnode_log;
  G4LogicalVolume* fGate_log;
  G4LogicalVolume* fCathode_log;
  G4LogicalVolume* fCathodeMesh_log;


  G4VPhysicalVolume* fLXePV;
  std::vector<G4VPhysicalVolume*> fSipmPVs;
  G4VPhysicalVolume* fAnodePV;
  std::vector<G4VPhysicalVolume*> fGatePVs;
  std::vector<G4VPhysicalVolume*> fCathodePVs;
  std::vector<G4VPhysicalVolume*> fCathodeMeshPVs;


  // Sensitive Detectors positions
  std::vector<G4ThreeVector> fSipmPositions;

};

#endif
