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
/// \file NUXE/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "LXeMainVolume.hh"

#include "globals.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4ExtrudedSolid.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMainVolume::LXeMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                             G4LogicalVolume* pMotherLogical, G4bool pMany,
                             G4int pCopyNo, LXeDetectorConstruction* c)
  : G4PVPlacement(pRot, tlate,
                  // Temp logical volume must be created here
                  new G4LogicalVolume(new G4Box("temp", 1, 1, 1),
                                      G4Material::GetMaterial("Vacuum"), "temp",
                                      0, 0, 0),
                  "housing", pMotherLogical, pMany, pCopyNo)
  , fConstructor(c)
{
    // Get all default detector config values;
  CopyValues();


  G4TwoVector off1(0, 0), off2(0, 0);
  G4double scale1 = 1, scale2 = 1;
  std::vector<G4TwoVector> polygon(fPolygon_edge_num);

  //*** parameters of polygon LXe volume;
  G4double pi = 2 * acos(0.0);
  G4double alpha = pi / fPolygon_edge_num;             //central angle of the polygon
  G4double beta = (pi - alpha) / 2.;                   //half interial angle
  G4double R = fPolygon_edge_len / (2. * sin(alpha));  //Circumradius
  G4double r = fPolygon_edge_len / (2. * tan(alpha));  //Inradius


  //*** Create the vertices of the polygon
  for(G4double i = 0.0; i < fPolygon_edge_num; i++){
      G4double theta = 2 * (i / fPolygon_edge_num) * pi;
      polygon[i].set(R * sin(theta), R * cos(theta));
  }
  std::cout << fPolygon_edge_num << "Polygon Vertices:" << std::endl;
  for(int i = 0; i < fPolygon_edge_num; i ++){
      std::cout << "********" << polygon[i] << std::endl;
  }

  //*** Create the Physical and Logical volume of LXe and Housing PTFE
  fHousing_cylinder = 
      new G4Tubs(
          "housing_tub", 0., 
          fHousing_r, fHousing_z,
          0., 360. * deg);
  fScint_polygon = 
      new G4ExtrudedSolid(
          "scint_polygon",
          polygon,
          fPolygon_height,
          off1, scale1,
          off2, scale2);
  fScint_log = 
      new G4LogicalVolume(
          fScint_polygon, 
          G4Material::GetMaterial("LXe"),
          "scint_log", 0, 0, 0);
  fHousing_log = 
      new G4LogicalVolume(
          fHousing_cylinder, 
          G4Material::GetMaterial("PTFE"), 
          "housing_log", 0, 0, 0);
  //*** placement of polygon LXe inside cylindrical PTFE
  fLXePV = new G4PVPlacement(0, G4ThreeVector(), fScint_log, "scintillator",
                    fHousing_log, false, 0);


  //*** Anode and gate wires:
  fAnode_wire =
      new G4Tubs(
          "anode_tub", 0., 
          fAnode_radius, fAnode_height,
          0., 360. * deg);
  fGate_wire =
      new G4Tubs(
          "gate_tub", 0., 
          fGate_radius, fGate_height,
          0., 360. * deg);
  fAnode_log =
      new G4LogicalVolume(
          fAnode_wire, 
          G4Material::GetMaterial("Al"),
          "anode_log", 0, 0, 0);
  fGate_log =
      new G4LogicalVolume(
          fGate_wire, 
          G4Material::GetMaterial("Al"),
          "gate_log", 0, 0, 0);
  //*** placement of anode and gate in LXe
  fAnodePV = new G4PVPlacement(0, G4ThreeVector(), fAnode_log, "Anode",
                    fScint_log, false, 0);

  G4int Ngate = 8;
  G4double Rgate = 3. * mm;
  std::vector<G4ThreeVector> gatePos[Ngate];

//   for(G4int i = 0; i < Ngate; ++i){
//     gatePos[i] = G4ThreeVector(sin(twopi * i / Ngate), cos(twopi * i / Nagte), 0.)
//   }

  for(G4int i = 0; i < Ngate; ++i){
    fGatePVs.emplace_back(new G4PVPlacement(0, G4ThreeVector(sin(CLHEP::twopi * i / Ngate)*Rgate, cos(CLHEP::twopi * i / Ngate)*Rgate, 0.),
                     fGate_log, "Gate",
                    fScint_log, false, i));
  }


  G4double Cathode_width = 30.5/2. * mm;
  G4double Cathode_height = 100./2. * mm;
  G4double Cathode_thickness = 0.2/2. * mm;
  G4double Cathode_R = 37. * mm;
  G4double Cathode_VEdge = 1. * mm;
  G4double Cathode_HEdge = 4.5 * mm;  
  
  G4double Cathode_inner_cut_width = Cathode_width - Cathode_VEdge;
  G4double Cathode_inner_cut_thickness = Cathode_thickness * 2.;
  G4double Cathode_inner_cut_height = Cathode_height- Cathode_HEdge;

  G4double Cathode_HGrid_length = Cathode_inner_cut_width;
  G4double Cathode_HGrid_thickness = Cathode_thickness;

  G4double Cathode_VGrid_length = Cathode_inner_cut_height;
  G4double Cathode_VGrid_thickness = Cathode_thickness;

  G4double Cathode_VGrid_pitch = 3.9 * mm;
  G4double Cathode_HGrid_pitch = 3.9* mm;

//Change the value of the following two lines to 1 or 2 for Visulization.
// a bigger value may take Forever for the UI to show up
  G4int N_VGrid = 1;
  G4int N_HGrid = 1;

  
  fCathode_box =
      new G4Box(
          "housing_tub", Cathode_width, Cathode_thickness, Cathode_height);

  fCathode_cut =
      new G4Box(
          "cathode_cut", Cathode_inner_cut_width, Cathode_inner_cut_thickness, Cathode_inner_cut_height);

  fCathode_HGrid =
      new G4Box(
          "cathode_HGrid", Cathode_HGrid_length, Cathode_HGrid_thickness,Cathode_HGrid_thickness );
  fCathode_VGrid =
      new G4Box(
          "cathode_HGrid", Cathode_VGrid_length, Cathode_VGrid_thickness,Cathode_VGrid_thickness );          
  
  
  fCathode_frame =
      new G4SubtractionSolid(
        "cathode_frame", fCathode_box, fCathode_cut);
  

  fCathode_Mesh = 
      new G4UnionSolid("cathode_mesh",fCathode_frame,fCathode_HGrid);

  for(G4int i = 0; i < N_HGrid; ++i){
    G4double d = Cathode_HGrid_pitch - Cathode_inner_cut_height + (Cathode_HGrid_pitch+Cathode_VGrid_thickness)*i;
    fCathode_Mesh = 
        new G4UnionSolid("cathode_mesh",fCathode_Mesh,fCathode_HGrid,
                        0,
                        G4ThreeVector(0.,0.,d) );
  }
  



  G4RotationMatrix* Rot_VGrid = new G4RotationMatrix();
  Rot_VGrid->rotateY(90. * deg);

  for(G4int i = 0; i < N_VGrid; ++i){
    G4double d = Cathode_VGrid_pitch - Cathode_inner_cut_width + (Cathode_VGrid_pitch+Cathode_VGrid_thickness)*i;
    fCathode_Mesh = 
        new G4UnionSolid("cathode_mesh",fCathode_Mesh,fCathode_VGrid,
                        Rot_VGrid,
                        G4ThreeVector(d,0.,0.) );
  }

  fCathode_log = 
      new G4LogicalVolume(
          fCathode_Mesh, 
          G4Material::GetMaterial("Al"),
          "cathodeFrame_log", 0, 0, 0);

  




  G4int NCathode = 8;
  std::vector<G4RotationMatrix*> Cathode_rot(NCathode);

  for(G4double i = 0.; i < NCathode; ++i){
    Cathode_rot[i] = new G4RotationMatrix();
    Cathode_rot[i]->rotateZ((CLHEP::twopi*(1./16.) + (CLHEP::twopi*(i/8.))) * rad);
    fCathodePVs.emplace_back(new G4PVPlacement(Cathode_rot[i],
                        G4ThreeVector(sin((CLHEP::twopi*(1./16.) + (CLHEP::twopi*(i/8.))))*Cathode_R,
                                      cos((CLHEP::twopi*(1./16.) + (CLHEP::twopi*(i/8.))))*Cathode_R, 0.),
                        fCathode_log, "Cathode", fScint_log, false, 0));
  }


 



  
  //**** Build SiPMs
  //Hamamatsu S13371
  //https://hamamatsu.su/files/uploads/pdf/3_mppc/s13370_vuv4-mppc_b_(1).pdf
  G4double l_sipm = 15 * mm;              //size; 1.5X1.5 [mm]
  G4double height_sipm = 6.5 * mm;        //thichness: (pin)4[mm] + (ceramic frame)2.5[mm] = 6.5[mm]
  G4double l_sipm_off = 1.8 * mm;         //gap between photocathode and edge: 0.2[mm] + (15[mm] - 2*5.9[mm])/2
  G4double height_photocath = 1.3 * mm;   //photocathode height: 1.3 [mm]
  G4double l_photocath = 5.9 * mm;
  G4double cell_pitch = 0.5 * mm;        
  G4double photocath_gap = 1.2 * mm;      //gap between photocathode and SiPM surface: 2.5[mm] - 1.3[mm] = 1.2 [mm]


  // Create the Physical and Logical Volume of SiPM/photocathode box
  fSiPM = 
      new G4Box("SiPM_box", l_sipm / 2.,
                    l_sipm / 2., height_sipm / 2.);
  fPhotocath = 
      new G4Box("photocath_box",
                l_photocath / 2.
                , l_photocath / 2.
                , height_photocath);
  fSiPM_log = 
      new G4LogicalVolume(fSiPM
                          , G4Material::GetMaterial("Glass"), "sipm_log");
  fPhotocath_log =
      new G4LogicalVolume(fPhotocath
                          , G4Material::GetMaterial("Al"), "photocath_log");

  

  // Placement of the photocathode boxes inside of the SiPM box
  new G4PVPlacement(0, G4ThreeVector(-(cell_pitch / 2.) - l_photocath / 2.
                                    , (cell_pitch / 2.) + l_photocath / 2.
                                    , - height_sipm / 2. + photocath_gap + height_photocath / 2. )
                    , fPhotocath_log, "photocath", fSiPM_log, false, 0);

  new G4PVPlacement(0, G4ThreeVector((cell_pitch / 2.) + l_photocath / 2.
                                    , (cell_pitch / 2.) + l_photocath / 2.
                                    , - height_sipm / 2. + photocath_gap + height_photocath / 2. )
                    , fPhotocath_log, "photocath", fSiPM_log, false, 1);

  new G4PVPlacement(0, G4ThreeVector(-(cell_pitch / 2.) - l_photocath / 2.
                                    , -(cell_pitch / 2.) - l_photocath / 2.
                                    , - height_sipm / 2. + photocath_gap + height_photocath / 2. )
                    , fPhotocath_log, "photocath", fSiPM_log, false, 2);

  new G4PVPlacement(0, G4ThreeVector((cell_pitch / 2.) + l_photocath / 2.
                                    , -(cell_pitch / 2.) - l_photocath / 2.
                                    , - height_sipm / 2. + photocath_gap + height_photocath / 2. )
                    , fPhotocath_log, "photocath", fSiPM_log, false, 3);
  
  

  // Create variables to arrange SiPMs:    
  G4int rowN = 6;                         // number of SiPMs in a row on PCB; 
  G4int colN = 2;                         // number of SiPMs in a column on PCB;
  G4double row_gap = 0.1 * mm;            // gap between SiPMs;
  G4double row_off = l_sipm + row_gap;    // Z offsite for SiPM translation;
  G4int sipm_count = 0;                   // count number of SiPM copies;

  // displacements are relative to the origin of its mother volume
  std::vector<G4double> dx(colN);         // XAxis displacement of SiPMs column on eace PCB;
  std::vector<G4double> dy(colN);         // YAxis displacement of SiPMs column on eace PCB;
  G4double dz;                            // ZAxis displacement of SiPMs row on each PCB;


  //**** Create a small polygon where SiPMs will be placed onto its edges;
  std::vector<G4TwoVector> m_polygon(fPolygon_edge_num);
  G4double r_diff = height_sipm + 1. * mm;                // gap between PTFE and the back of SiPMs;
  G4double m_r = r - r_diff;                              // inradius of the small polygon
  G4double m_R = m_r / cos(alpha);                        // circumradius
  for(G4double i = 0.0; i < fPolygon_edge_num; i++){
      G4double theta = 2 * (i / fPolygon_edge_num) * pi; 
      m_polygon[i].set(m_R * sin(theta), m_R * cos(theta));
  }
//   std::cout << fPolygon_edge_num << "m_Polygon Vertices:" << std::endl;
//   for(int i = 0; i < fPolygon_edge_num; i ++){
//       std::cout << "********" << m_polygon[i] << std::endl;
//   }

  std::vector<G4TwoVector> x_y(colN * fPolygon_edge_num);     // (X,Y) coordinates of SiPMs' positions on XY_plane
  G4int next_edgeN = 1;                                       // edge number counts
  G4int this_edgeN = 0;                                       // edge number counts
  G4int xy_count = 0;                                         // translation counts on XY_plane;
  G4double x1, x2, y1, y2;                                    // two temporary (X,Y) coordinates update to the above vector;
  G4double gap = 0.1 * mm;                                    // gap between two column SiPMs on the PCB;
  G4double d = (l_sipm / 2.) + gap;                           // distance from the center of small polygon edge to the center of SiPM
  
  //Finding the (X,Y) coordinates of SiPMs for placements
  for(G4int i = 0; i < fPolygon_edge_num; i++){               
      G4double slope = 0.;
      G4double C = 0.;

      if(next_edgeN >= fPolygon_edge_num){
          next_edgeN = 0;
      }
        
      slope = (m_polygon[next_edgeN][1] - m_polygon[this_edgeN][1])
            / (m_polygon[next_edgeN][0] - m_polygon[this_edgeN][0]);

    //   std::cout << "This slope is: " << slope << std::endl;

      C = m_polygon[this_edgeN][1] - slope * m_polygon[this_edgeN][0];

    //   std::cout << "This C is: " << C << std::endl; 

      for(G4int j = 0; j < colN; j++){
          x1 = (-2.*C*slope + sqrt(4.*slope*slope*C*C - 4.*(slope*slope+1.)*(C*C-d*d-m_r*m_r)))/(2.*(slope*slope+1.));
          y1 = slope * x1 + C;
          x2 = (-2.*C*slope - sqrt(4.*slope*slope*C*C - 4.*(slope*slope+1.)*(C*C-d*d-m_r*m_r)))/(2.*(slope*slope+1.));
          y2 = slope * x2 + C;

          x_y[xy_count].set(x1, y1);
          x_y[xy_count + 1].set(x2, y2);
      }
      xy_count += 2;
      next_edgeN += 1;
      this_edgeN += 1;
  }
//   std::cout << "************ XY Translation::" << std::endl;
//   for(int i = 0; i < colN*fPolygon_edge_num; i+=2){
//       std::cout << "************ " << x_y[i][0] << "; " << x_y[i][1] << std::endl;
//       std::cout << "************ " << x_y[i+1][0] << "; " << x_y[i+1][1] << std::endl;
    //  }


  //*** Placement of SiPMs inside of LXe and PTFE housing;
  std::vector<G4RotationMatrix*> rot(fPolygon_edge_num);
  G4int mxy_count = 0;
  for(G4double i = 0.; i < fPolygon_edge_num; i++){
      for(G4int j = 0; j < colN; j++){
          rot[i] = new G4RotationMatrix();
          rot[i]->rotateY(90. * degree);
          rot[i]->rotateX((2. * (i+3.) - 1.) * alpha * rad);
          dx[j] = x_y[mxy_count+j][0];
          dy[j] = x_y[mxy_count+j][1];
          for(G4int k = 0; k < rowN; k++){
              dz = ((fPolygon_height / 1.) - (l_sipm/2.) - (fPolygon_height*2.-(rowN*l_sipm + (rowN-1)*row_gap))/2. - k * row_off);
             fSipmPVs.emplace_back(new G4PVPlacement(rot[i], G4ThreeVector(dx[j], dy[j], dz)
                                , fSiPM_log, "SiPM", fScint_log, false, sipm_count));
              sipm_count += 1;
              fSipmPositions.push_back(G4ThreeVector(dx[j], dy[j], dz));
          }
      }
      mxy_count += 2;
  }

  
  

  VisAttributes();

  SurfaceProperties();

  SetLogicalVolume(fHousing_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues()
{
  fHousing_r        = fConstructor->GetHousingR();
  fHousing_z        = fConstructor->GetHousingZ();
  fPolygon_edge_num = fConstructor->GetPolygonEdgeNum();
  fPolygon_edge_len = fConstructor->GetPolygonEdgeLen();
  fPolygon_height   = fConstructor->GetPolygonHeight();
  fRefl             = fConstructor->GetHousingReflectivity();
  fAnode_radius     = fConstructor->GetAnodeRadius();
  fGate_radius      = fConstructor->GetGateRadius();
  fAnode_height     = fConstructor->GetAnodeHeight();
  fGate_height      = fConstructor->GetGateHeight();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes()
{
//   G4VisAttributes* housing_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
//   fHousing_log->SetVisAttributes(housing_va);

  // G4VisAttributes* sphere_va = new G4VisAttributes();
  // sphere_va->SetForceSolid(true);
  // fSphere_log->SetVisAttributes(sphere_va);

G4VisAttributes* housing_va = new G4VisAttributes(G4Color(0.43, 0.43, 0.43, 0.7));
fHousing_log->SetVisAttributes(housing_va);

G4VisAttributes* scint_va = new G4VisAttributes(G4Color(0.05, 0.05, 0.9, 0.4));
// scint_va->SetForceSolid(true);
fScint_log->SetVisAttributes(scint_va);

G4VisAttributes* sipm_va = new G4VisAttributes(G4Color(0.9, 0.9, 0.9, 0.4));
// sipm_va->SetForceSolid(true);
fSiPM_log->SetVisAttributes(sipm_va);

G4VisAttributes* photocath_va = new G4VisAttributes(G4Color(0.1,0.8,0.1));
// photocath_va->SetForceSolid(true);
fPhotocath_log->SetVisAttributes(photocath_va);

G4VisAttributes* anode_va = new G4VisAttributes(G4Color(1.,0.,0.));
anode_va->SetForceSolid(true);
fAnode_log->SetVisAttributes(anode_va);

G4VisAttributes* gate_va = new G4VisAttributes(G4Color(0.5,0.,0.5));
gate_va->SetForceSolid(true);
fGate_log->SetVisAttributes(gate_va);

G4VisAttributes* cathode_va = new G4VisAttributes(G4Color(0.5,0.5,0.));
cathode_va->SetForceSolid(true);
fCathode_log->SetVisAttributes(cathode_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::SurfaceProperties()
{
std::vector<G4double> ephoton = { 7.0 * eV, 7.14 * eV };

//PTFE
std::vector<G4double> ptfe_reflectivity = { 0.99, 0.99 };
std::vector<G4double> ptfe_efficiency = { 0.0, 0.0 };
G4MaterialPropertiesTable* PTFE_pt = new G4MaterialPropertiesTable();
PTFE_pt->AddProperty("REFLECTIVITY", ephoton, ptfe_reflectivity);
PTFE_pt->AddProperty("EFFICIENCY", ephoton, ptfe_efficiency);
G4OpticalSurface* opPTFEHousingSurface = 
    // new G4OpticalSurface("PTFESurface", unified, polishedteflonair, dielectric_metal);
    new G4OpticalSurface("PTFESurface");
opPTFEHousingSurface->SetType(dielectric_metal);
opPTFEHousingSurface->SetModel(unified);
opPTFEHousingSurface->SetFinish(polished);
opPTFEHousingSurface->SetMaterialPropertiesTable(PTFE_pt);
new G4LogicalBorderSurface("LXePTFE",fLXePV,this,opPTFEHousingSurface);
// new G4LogicalSkinSurface("PTFE_surf", fHousing_log, opPTFEHousingSurface);


//Quartz Window
std::vector<G4double> Qz_reflectivity = { 0.0001, 0.0001 };
// std::vector<G4double> Qz_efficiency = { 0.99, 0.99 };
std::vector<G4double> Qz_transmission = { 1., 1., };
G4MaterialPropertiesTable* Qz_pt = new G4MaterialPropertiesTable();
Qz_pt->AddProperty("REFLECTIVITY", ephoton, Qz_reflectivity);
Qz_pt->AddProperty("TRANSMITTANCE", ephoton, Qz_transmission);
// Qz_pt->AddProperty("EFFICIENCY", ephoton, Qz_efficiency);
G4OpticalSurface* opQuartzSurface =
    new G4OpticalSurface("QuartzSurface");
opQuartzSurface->SetType(dielectric_dielectric);
opQuartzSurface->SetModel(unified);
opQuartzSurface->SetFinish(polished);
opQuartzSurface->SetMaterialPropertiesTable(Qz_pt);
// new G4LogicalSkinSurface("LXeQuartz",fSiPM_log,opQuartzSurface);
for(G4int i = 0; i < fSipmPVs.size(); i++){
    new G4LogicalBorderSurface("LXeQuartz",fLXePV,fSipmPVs[i],opQuartzSurface);
}

//Photocathode
std::vector<G4double> photocath_eff = { 1., 1. };
std::vector<G4double> photocath_ReR = { 1.92, 1.92 };
std::vector<G4double> photocath_ImR = { 1.69, 1.69 };
G4MaterialPropertiesTable* photocath_pt = new G4MaterialPropertiesTable();
photocath_pt->AddProperty("EFFICIENCY", ephoton, photocath_eff);
photocath_pt->AddProperty("REALRINDEX", ephoton, photocath_ReR);
photocath_pt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR);
G4OpticalSurface* photocath_opsurf = new G4OpticalSurface(
    "photocath_opsurf", unified, polished, dielectric_metal
);
photocath_opsurf->SetMaterialPropertiesTable(photocath_pt);
new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);

//Anode
std::vector<G4double> Anode_reflectivity = { .9, .9 };
std::vector<G4double> Anode_efficiency = { 0., 0. };
G4MaterialPropertiesTable* Anode_pt = new G4MaterialPropertiesTable();
Anode_pt->AddProperty("REFLECTIVITY", ephoton, Anode_reflectivity);
Anode_pt->AddProperty("EFFICIENCY", ephoton, Anode_efficiency);
G4OpticalSurface* opAnodeSurface =
    new G4OpticalSurface("AnodeSurface");
opAnodeSurface->SetType(dielectric_metal);
opAnodeSurface->SetModel(unified);
opAnodeSurface->SetFinish(polished);
opAnodeSurface->SetMaterialPropertiesTable(Anode_pt);
new G4LogicalBorderSurface("LXeAnode",fLXePV,fAnodePV,opAnodeSurface);

//Gate
std::vector<G4double> Gate_reflectivity = { .7, .7 };
std::vector<G4double> Gate_efficiency = { 0., 0. };
G4MaterialPropertiesTable* Gate_pt = new G4MaterialPropertiesTable();
Gate_pt->AddProperty("REFLECTIVITY", ephoton, Gate_reflectivity);
Gate_pt->AddProperty("EFFICIENCY", ephoton, Gate_efficiency);
G4OpticalSurface* opGateSurface =
    new G4OpticalSurface("GateSurface");
opGateSurface->SetType(dielectric_metal);
opGateSurface->SetModel(unified);
opGateSurface->SetFinish(polished);
opGateSurface->SetMaterialPropertiesTable(Gate_pt);
for(G4int i = 0; i < fGatePVs.size(); ++i){
    new G4LogicalBorderSurface("LXeGate",fLXePV,fGatePVs[i],opGateSurface);
}


//Cathode
std::vector<G4double> Cathode_reflectivity = { 0.7, 0.7 };
std::vector<G4double> Cathode_efficiency = { 0., 0. };
G4MaterialPropertiesTable* Cathode_pt = new G4MaterialPropertiesTable();
Cathode_pt->AddProperty("REFLECTIVITY", ephoton, Cathode_reflectivity);
Cathode_pt->AddProperty("EFFICIENCY", ephoton, Cathode_efficiency);
G4OpticalSurface* opCathodeSurface =
    new G4OpticalSurface("CathodeSurface");
opCathodeSurface->SetType(dielectric_metal);
opCathodeSurface->SetModel(unified);
opCathodeSurface->SetFinish(polished);
opCathodeSurface->SetMaterialPropertiesTable(Cathode_pt);
for(G4int i = 0; i < fCathodePVs.size(); ++i){
  new G4LogicalBorderSurface("LXeCathode",fLXePV,fCathodePVs[i],opCathodeSurface);
}
}
