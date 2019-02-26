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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"

#include "B1DetectorMessenger.hh"

#include "G4UserLimits.hh" //George
#include "G4StepStatus.hh" //George
#include "G4StepPoint.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fScoringVolume(0),
  fScatXPos(0),
  fScatYPos(0),
  fScatPolarR(0),
  fScatPolarPhi(0),
  fScatRotX(0),
  fScatRotY(0),
  fScatRotZ(0),
  fScatRad(0),
  fScatHeight(0),
  fDetXPos(0),
  fDetYPos(0),
  fDetPolarR(0),
  fDetPolarPhi(0),
  fDetRotX(0),
  fDetRotY(0),
  fDetRotZ(0),
  fDetRad(0),
  fDetHeight(0),
  fDetectorMessenger(0)
{
//Material Definition  of Lanthanum(III) Bromide - by Ben
G4Isotope* iso79_Br = new G4Isotope("79_Br", 35, 79, 78.918*g/mole);
G4Isotope* iso81_Br = new G4Isotope("81_Br", 35, 81, 80.916*g/mole);
G4Isotope* iso138_La = new G4Isotope("138_La", 57, 138, 137.907*g/mole);
G4Isotope* iso139_La = new G4Isotope("139_La", 57, 139, 138.906*g/mole);

G4Element* Br = new G4Element("Bromide", "Br", 2);
Br->AddIsotope(iso79_Br, 50.69*perCent);
Br->AddIsotope(iso81_Br, 49.31*perCent);
G4Element* La = new G4Element("Lanthanum", "La", 2);
La->AddIsotope(iso138_La, 0.08881*perCent);
La->AddIsotope(iso139_La, 99.9119*perCent);

G4Material* LaBr = 
new G4Material("Lanthanum_Bromide", 5.06*g/cm3, 2, kStateSolid);
LaBr->AddElement( Br, 0.6036373016896684);
LaBr->AddElement( La, 0.39636269831033155);

//More Materials - George
//G4double a;  // atomic mass
//G4double z;  // atomic number
//G4double density,ncomponents,fractionmass,nel, symbol;

//G4Element* elGe = new G4Element("Ge", 32, 72.64*g/mole);
//G4Element* elBi = new G4Element("Bi", 83, 208.9*g/mole);
//G4Element* O = new G4Element("O", 8, 16.00*g/mole);
//G4Material* BGO = new G4Material("BGO", 7.13*g/cm3, 3);
//BGO->AddElement(elBi, 4);//no. el.
//BGO->AddElement(elGe, 3);
//BGO->AddElement(O, 12);

//G4Element* Pb = new G4Element("Pb", 82, 207.20*g/mole); //different format
//G4Element* W  = new G4Element("W" , 74, 183.84*g/mole);
//G4Material* PWOD = new G4Material("PbWO" , 8.300*g/cm3, 3);
//PWOD->AddElement(Pb, fractionmass=0.455);
//PWOD->AddElement(W , fractionmass=0.404);
//PWOD->AddElement(O , fractionmass=0.141);
 
fDetectorMessenger = new B1DetectorMessenger(this);
 
fScatXPos = 0;
fScatYPos = 0;
fScatPolarR = 0;
fScatPolarPhi = 9000; // 9000 is a default value to get the code to ignore this if we use cartesian
fScatRotX = 90*deg;
fScatRotY = 0;
fScatRotZ = 0;
fScatRad = 7*mm;
fScatHeight = 21.5*mm;
fDetXPos = -12.15*cm;
fDetYPos = 26.205*cm;
fDetPolarR = 0;
fDetPolarPhi = 9000;
fDetRotX = 0;
fDetRotY = 30*deg;
fDetRotZ = 0;
fDetRad = 1.905*cm;
fDetHeight = 1.905*cm;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 100*cm, env_sizeZ = 200*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =  new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
  
  G4VisAttributes* logicWorldVis = new G4VisAttributes();
  logicWorldVis->SetVisibility(false);
  logicWorld->SetVisAttributes(logicWorldVis);
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //
  G4RotationMatrix* rotEnv = new G4RotationMatrix();
  //rotEnv->rotateX(90*deg);
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
  
  G4VisAttributes* logicEnvVis = new G4VisAttributes();
  logicEnvVis->SetVisibility(false);
  logicEnv->SetVisAttributes(logicEnvVis);

               
  new G4PVPlacement(rotEnv,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4RotationMatrix* rot1 = new G4RotationMatrix();
  rot1->rotateX(fScatRotX);
  rot1->rotateY(fScatRotY);
  rot1->rotateZ(fScatRotZ);
  
  
  //
  //     
  // Shape 1
  //  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4ThreeVector* pos1 = new G4ThreeVector(fScatXPos, 0*cm, fScatYPos);
  if(fScatPolarR!=0)
    {
      pos1->setMag(fScatPolarR);
    }
  if(fScatPolarPhi!=9000)
    {
      pos1->setTheta(fScatPolarPhi);
    }
        
  // Tube section shape       
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = fScatRad;
  G4double shape1_hz = fScatHeight;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Tubs* solidShape1 =    
    new G4Tubs("Shape1", shape1_rmina, shape1_rmaxa, shape1_hz, shape1_phimin, shape1_phimax);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Scatterer");           //its name
               
  new G4PVPlacement(rot1,                    //rotation
                    G4ThreeVector(pos1->x(), pos1->y(), pos1->z()),  //at position
                    logicShape1,             //its logical volume
                    "Scatterer",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // Body (George) 
  //

  G4bool wantBody = false; 
  if(wantBody==true){
	G4Material* bodyMat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
	G4ThreeVector* bodyPos = new G4ThreeVector(0.*cm,0.*cm,0.*cm); //needs to be around source
	G4double bodyHeight = 10*cm;
	G4double bodyRadius = 18*cm; //Average according to some website - get proper source
	G4Tubs* bodyShape =    
    new G4Tubs("Body", 0*cm, bodyRadius/2, bodyHeight/2, 0*deg, 360*deg);
	G4LogicalVolume* bodyVolume = new G4LogicalVolume(bodyShape,bodyMat,"Body"); 
	new G4PVPlacement(rot1,                    //rotation
                    G4ThreeVector(bodyPos->x(), bodyPos->y(), bodyPos->z()),  //at position
                    bodyVolume,             //its logical volume
                    "Body",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
	};

  //     
  // Shape 2
  //
  G4ThreeVector* pos2 = new G4ThreeVector(fDetXPos, 0*cm, fDetYPos);
  if(fDetPolarR!=0)
    {
      pos2->setMag(fDetPolarR);
    }
  if(fDetPolarPhi!=9000)
    {
      pos2->setTheta(fDetPolarPhi);
    }
  G4Material* shape2_mat = nist->FindOrBuildMaterial("Lanthanum_Bromide");
  G4RotationMatrix* rot2 = new G4RotationMatrix();
  rot2->rotateX(fDetRotX);
  rot2->rotateY(fDetRotY);
  rot2->rotateZ(fDetRotZ);

  G4double shape2_rmina =  0.*cm, shape2_rmaxa = fDetRad;
  G4double shape2_hz = fDetHeight;
  G4double shape2_phimin = 0.*deg, shape2_phimax = 360.*deg;
  G4Tubs* solidShape2 =    new G4Tubs("Shape2",
    shape2_rmina, shape2_rmaxa, shape2_hz,
    shape2_phimin, shape2_phimax);
                
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Absorber");           //its name
               
  new G4PVPlacement(rot2,                       //no rotation
                    G4ThreeVector(pos2->x(), pos2->y(), pos2->z()),        		//at position
                    logicShape2,             //its logical volume
                    "Absorber",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  
  //Varying step length depending on the logical volume 
  G4double stepLength = 0.0005*mm;
  G4UserLimits* maxStep = new G4UserLimits(stepLength); 
  logicShape1->SetUserLimits(maxStep);
  logicShape2->SetUserLimits(maxStep);
 

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// By Jack
void B1DetectorConstruction::SetScatXPos(G4double val)
{
  fScatXPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatYPos(G4double val)
{
  fScatYPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatPolarR(G4double val)
{
  fScatPolarR = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatPolarPhi(G4double val)
{
  fScatPolarPhi = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatRotX(G4double val)
{
  fScatRotX = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatRotY(G4double val)
{
  fScatRotY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatRotZ(G4double val)
{
  fScatRotZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatRad(G4double val)
{
  fScatRad = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScatHeight(G4double val)
{
  fScatHeight = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetXPos(G4double val)
{
  fDetXPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetYPos(G4double val)
{
  fDetYPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetPolarR(G4double val)
{
  fDetPolarR = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetPolarPhi(G4double val)
{
  fDetPolarPhi = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetRotX(G4double val)
{
  fDetRotX = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetRotY(G4double val)
{
  fDetRotY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetRotZ(G4double val)
{
  fDetRotZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetRad(G4double val)
{
  fDetRad = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDetHeight(G4double val)
{
  fDetHeight = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
