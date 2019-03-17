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

#include <cmath> // Douglas

#include <assert.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fScoringVolume(0),
  fScatXPos(0),
  fScatYPos(0),
  fScatZPos(0),
  fScatPolarR(0),
  fScatPolarPhi(0),
  fScatPolarTheta(0),
  fScatRotX(0),
  fScatRotY(0),
  fScatRotZ(0),
  fScatRad(0),
  fScatHeight(0),
  fScat2XPos(0),
  fScat2YPos(0),
  fScat2ZPos(0),
  fScat2PolarR(0),
  fScat2PolarPhi(0),
  fScat2PolarTheta(0),
  fScat2RotX(0),
  fScat2RotY(0),
  fScat2RotZ(0),
  fScat2Rad(0),
  fScat2Height(0),
  fScat2Bool(0),
  fDetXPos(0),
  fDetYPos(0),
  fDetZPos(0),
  fDetPolarR(0),
  fDetPolarPhi(0),
  fDetPolarTheta(0),
  fDetRotX(0),
  fDetRotY(0),
  fDetRotZ(0),
  fDetRad(0),
  fDetHeight(0),
  fDet2XPos(0),
  fDet2YPos(0),
  fDet2ZPos(0),
  fDet2PolarR(0),
  fDet2PolarPhi(0),
  fDet2PolarTheta(0),
  fDet2RotX(0),
  fDet2RotY(0),
  fDet2RotZ(0),
  fDet2Rad(0),
  fDet2Height(0),
  fDet2Bool(0),
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

//G4Material* LaBr =  THIS IS THE OLD METHOD, NOT SURE OF IT'S RELIABILITY
//new G4Material("Lanthanum_Bromide", 5.06*g/cm3, 2, kStateSolid);
//LaBr->AddElement( Br, 0.6036373016896684);
//LaBr->AddElement( La, 0.39636269831033155);

//More Materials - George
 G4NistManager* nist = G4NistManager::Instance();
 G4Material* CeF3 = nist->FindOrBuildMaterial("G4_CERIUM_FLUORIDE");
 G4Material* CI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
 G4Material* CaF2 = nist->FindOrBuildMaterial("G4_CALCIUM_FLUORIDE");
 G4Material* BaF2 = nist->FindOrBuildMaterial("G4_BARIUM_FLUORIDE");
 
 //George Materials 
//Lanthanum Bromide 
G4Material* LaBr = new G4Material("Lanthanum_Bromide", 5.06*g/cm3, 2, kStateSolid);
LaBr->AddElement(nist->FindOrBuildElement(57),1); //Lanthanum
LaBr->AddElement(nist->FindOrBuildElement(35),3); //Bromine

 //eRes of 10.5% at 662keV  
 G4Material* LSO = new G4Material("LSO", 7.4*g/cm3, 3, kStateSolid);
 LSO->AddElement(nist->FindOrBuildElement(71),2); // Lutetium
 LSO->AddElement(nist->FindOrBuildElement(14),1); //Sulfur
 LSO->AddElement(nist->FindOrBuildElement(8),5); //Oxygen

 //eRes of 6.6% at 662keV for 20×20mm cylindrical crystal - Resource = CdWO4 Crystal in Gamma-Ray Spectrometry 
 G4Material* CdWO4 = new G4Material("CdWO4", 7.9*g/cm3, 3, kStateSolid);
  CdWO4->AddElement(nist->FindOrBuildElement(48),1);  //Cadium
  CdWO4->AddElement(nist->FindOrBuildElement(78),1);  //Tungsten
  CdWO4->AddElement(nist->FindOrBuildElement(8),4);   //Oxygen
 
 //eRes of  13.35% at 662keV for 76.2mm*76.2mm cylindrical cystal - Resource = Comparative Study of Large NaI(Tl) and BGO Scintillators
 G4Material* BGO = new G4Material("BGO", 7.13*g/cm3, 3, kStateSolid);
  BGO->AddElement(nist->FindOrBuildElement(83),4); //Bismuth
  BGO->AddElement(nist->FindOrBuildElement(32),3); //Germanium
  BGO->AddElement(nist->FindOrBuildElement(8),12); //Oxygen

//Potential scatterer? Give it a quick test. 2.9% at 662 keV 
G4Material* Cs2NaLaBr3I3 = new G4Material("theBeast", 4.00*g/cm3, 5, kStateSolid);
  Cs2NaLaBr3I3->AddElement(nist->FindOrBuildElement(55),2); //Caesium
  Cs2NaLaBr3I3->AddElement(nist->FindOrBuildElement(11),1); //Germanium
  Cs2NaLaBr3I3->AddElement(nist->FindOrBuildElement(57),1); //Oxygen
  Cs2NaLaBr3I3->AddElement(nist->FindOrBuildElement(35),3); //Bismuth
  Cs2NaLaBr3I3->AddElement(nist->FindOrBuildElement(53),3); //Germanium

fDetectorMessenger = new B1DetectorMessenger(this);
 
fScatXPos = 0;
fScatYPos = 0;
fScatZPos = 0;
fScatPolarR = 0;
fScatPolarPhi = 9000; // 9000 is a default value to get the code to ignore this if we use cartesian
fScatPolarTheta = 9000;
fScatRotX = 90*deg;
fScatRotY = 0;
fScatRotZ = 0;
fScatRad = 7*mm;
fScatHeight = 21.5*mm;
SetScatMat("G4_SODIUM_IODIDE");

fScat2XPos = 10*cm;
fScat2YPos = 0;
fScat2ZPos = 0;
fScat2PolarR = 0;
fScat2PolarPhi = 9000; // 9000 is a default value to get the code to ignore this if we use cartesian
fScat2PolarTheta = 9000;
fScat2RotX = 90*deg;
fScat2RotY = 0;
fScat2RotZ = 0;
fScat2Rad = 7*mm;
fScat2Height = 21.5*mm;
SetScat2Mat("G4_SODIUM_IODIDE");
fScat2Bool = false;
fScat2Type = "Scatterer";
fScat2CopyNo = 1;

fDetXPos = -12.15*cm;
fDetYPos = 0;
fDetZPos = 26.205*cm;
fDetPolarR = 0;
fDetPolarPhi = 9000;
fDetPolarTheta = 9000;
fDetRotX = 0;
fDetRotY = 0;
fDetRotZ = 0;
fDetRad = 1.905*cm;
fDetHeight = 1.905*cm;
SetDetMat("Lanthanum_Bromide");

fDet2XPos = 12.15*cm;
fDet2YPos = 0;
fDet2ZPos = 26.205*cm;
fDet2PolarR = 0;
fDet2PolarPhi = 9000;
fDet2PolarTheta = 9000;
fDet2RotX = 0;
fDet2RotY = -30*deg;
fDet2RotZ = 0;
fDet2Rad = 1.905*cm;
fDet2Height = 1.905*cm;
SetDet2Mat("Lanthanum_Bromide");
fDet2Bool = false;
fDet2Type = "Absorber";
fDet2CopyNo = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{
  delete fDetectorMessenger;
}

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

  G4RotationMatrix* rot3 = new G4RotationMatrix();
  rot3->rotateX(fScat2RotX);
  rot3->rotateY(fScat2RotY);
  rot3->rotateZ(fScat2RotZ);
  

//These two booleans are for George for Material testing and investigating the Pixelated detector. They are both definitely commented out in the github version.

  //G4bool wantEverything = false;  // for PixelatedDetector testing
  //if(wantEverything==true){

  //G4bool wantScatterer = false;    //for MaterialTesting (Just the one absorber) 
  //if(wantScatterer==true){


  //
  //     
  // Scatterer
  //
  G4ThreeVector* pos1 = new G4ThreeVector(fScatXPos, fScatYPos, fScatZPos);
  if(fScatPolarR!=0)
    {
      pos1->setMag(fScatPolarR);
    }
  if(fScatPolarPhi!=9000)
    {
      pos1->setTheta(fScatPolarPhi);
    }
  if(fScatPolarTheta!=9000)
    {// theta of scatterer by Doug
      pos1->setPhi(fScatPolarTheta);
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

  //}; //End of if want scatterer GEORGE


  //
  //     
  // Scatterer 2 - By Jack
  //
  if (fScat2Bool == true)
    {
      G4ThreeVector* pos3 = new G4ThreeVector(fScat2XPos, fScat2YPos, fScat2ZPos);
      if(fScat2PolarR!=0)
	{
	  pos3->setMag(fScat2PolarR);
	}
      if(fScat2PolarPhi!=9000)
	{
	  pos3->setTheta(fScat2PolarPhi);
	}
      if(fScat2PolarTheta!=9000)
	{
	  pos3->setPhi(fScat2PolarTheta);
	}
        
      // Tube section shape       
      G4double shape3_rmina =  0.*cm, shape3_rmaxa = fScat2Rad;
      G4double shape3_hz = fScat2Height;
      G4double shape3_phimin = 0.*deg, shape3_phimax = 360.*deg;
      G4Tubs* solidShape3 =    
	new G4Tubs("Shape3", shape3_rmina, shape3_rmaxa, shape3_hz, shape3_phimin, shape3_phimax);
                      
      G4LogicalVolume* logicShape3 =                         
	new G4LogicalVolume(solidShape3,         //its solid
			    shape3_mat,          //its material
			    fScat2Type);           //its name
               
      new G4PVPlacement(rot3,                    //rotation
			G4ThreeVector(pos3->x(), pos3->y(), pos3->z()),  //at position
			logicShape3,             //its logical volume
			fScat2Type,                //its name
			logicEnv,                //its mother  volume
			false,                   //no boolean operation
			fScat2CopyNo,            //copy number
			checkOverlaps);          //overlaps checking

    }// End of if scat2Bool
 
  //     
  // Absorber
  //
  G4ThreeVector* pos2 = new G4ThreeVector(fDetXPos, fDetYPos, fDetZPos);
  if(fDetPolarR!=0)
    {
      pos2->setMag(fDetPolarR);
    }
  if(fDetPolarPhi!=9000)
    {
      pos2->setTheta(fDetPolarPhi);
    }
  if(fDetPolarTheta!=9000)
    {
      pos2->setPhi(fDetPolarTheta);
    }
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
  


  //     
  // Absorber 2 - by Jack
  //
  if (fDet2Bool == true)
    {
      G4ThreeVector* pos4 = new G4ThreeVector(fDet2XPos, fDet2YPos, fDet2ZPos);
      if(fDet2PolarR!=0)
	{
	  pos4->setMag(fDet2PolarR);
	}
      if(fDet2PolarPhi!=9000)
	{
	  pos4->setTheta(fDet2PolarPhi);
	}
      if(fDet2PolarTheta!=9000)
	{
	  pos4->setPhi(fDet2PolarTheta);
	} 
      G4RotationMatrix* rot4 = new G4RotationMatrix();
      rot4->rotateX(fDet2RotX);
      rot4->rotateY(fDet2RotY);
      rot4->rotateZ(fDet2RotZ);

      G4double shape4_rmina =  0.*cm, shape4_rmaxa = fDet2Rad;
      G4double shape4_hz = fDet2Height;
      G4double shape4_phimin = 0.*deg, shape4_phimax = 360.*deg;
      G4Tubs* solidShape4 =    new G4Tubs("Shape4",
					  shape4_rmina, shape4_rmaxa, shape4_hz,
					  shape4_phimin, shape4_phimax);
                
      G4LogicalVolume* logicShape4 =                         
	new G4LogicalVolume(solidShape4,         //its solid
			    shape4_mat,          //its material
			    fDet2Type);          //its name
               
      new G4PVPlacement(rot4,                       //no rotation
			G4ThreeVector(pos4->x(), pos4->y(), pos4->z()),        		//at position
			logicShape4,             //its logical volume
			fDet2Type,               //its name
			logicEnv,                //its mother  volume
			false,                   //no boolean operation
			fDet2CopyNo,             //copy number
			checkOverlaps);          //overlaps checking
    

  } // ends activation of absorber 2

  //};// ends wantEverything


  // George's toy (The Pixelated Detectors - made out of 'ideal' materials) 
  G4bool pixelatedOn = false; //turns the pixelated detector on or off
  if(pixelatedOn==true){

  G4double pixelWidth = 0.55*cm; 
  G4double pixelHeight = 1.26*cm;
  G4double pixelDepth = 1*cm;
  G4Box* solidCrystal = new G4Box("Crystal", pixelWidth/2, pixelHeight/2, pixelDepth/2); 	
  G4Material* pixelatedScat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4Material* pixelatedAbsorb = nist->FindOrBuildMaterial("CdWO4");
  G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidCrystal,pixelatedScat,"Scatterer");
  G4LogicalVolume* logicCrystal2 = new G4LogicalVolume(solidCrystal,pixelatedAbsorb,"Absorber");

  int HoriNofCrystals = 8; 
  int VertNofCrystals = 4; 
  for (int j=0;j<VertNofCrystals;j++) { 
  for (int i=0;i<HoriNofCrystals;i++) { 
  G4double x2 = ((i-3.5)*0.56)*cm; 
  G4double y2 = ((j-1.5)*1.27)*cm; 
  G4ThreeVector centreOfPixel = G4ThreeVector(x2,y2,0*cm); 
   new G4PVPlacement(0, centreOfPixel,logicCrystal, //George
                "Scatterer",logicEnv, 
                false,i+8*j,checkOverlaps);};};

  G4double angle = 38*deg; //Desired angle between plain of scatterer and absorber
  G4double separationOfArray = 20; //Desired straight line separation between plains
  G4ThreeVector centreOfArray = G4ThreeVector((-separationOfArray*sin(angle))*cm,0*cm,(separationOfArray*cos(angle))*cm); 
  G4RotationMatrix* rotPixelAbsorb = new G4RotationMatrix();
  rotPixelAbsorb->rotateY(angle); //
  for (int j=0;j<VertNofCrystals;j++) { //George 
  for (int i=0;i<HoriNofCrystals;i++) { //George
  G4double x2 = ((i-3.5)*(0.56)*cos(angle))*cm; 
  G4double y2 = ((j-1.5)*(1.27))*cm;
  G4double z2 = ((i-3.5)*(0.56)*sin(angle));  
  //G4double zSepCentres = 4;
  G4ThreeVector relcentreOfPixel2 = G4ThreeVector(x2,y2,(z2)*cm); 
  G4ThreeVector centreOfPixel2 = relcentreOfPixel2 + centreOfArray;
  new G4PVPlacement(rotPixelAbsorb, centreOfPixel2,logicCrystal2, 
                "Absorber",logicEnv, 
                false,i+8*j,checkOverlaps);};};
 
}//ends the turn Pixelated detector off statement

//Varying step length depending on the logical volume 
  G4double maxStep = 0.01*mm; //0.01 = an acceptable speed but quite slow
  G4UserLimits* stepLimit = new G4UserLimits(); 
  stepLimit->SetMaxAllowedStep(maxStep);
  logicShape1->SetUserLimits(stepLimit);
  logicShape2->SetUserLimits(stepLimit);

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

void B1DetectorConstruction::SetScatZPos(G4double val)
{
  fScatZPos = val;
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

void B1DetectorConstruction::SetScatPolarTheta(G4double val)
{
  fScatPolarTheta = val;
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

void B1DetectorConstruction::SetScatMat(G4String val)
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* temp_mat = nist->FindOrBuildMaterial(val);
  if(temp_mat == NULL){
    std::cout << val << " does not exist in NIST" << std::endl;
    assert(temp_mat != NULL);
  }
  else {
    std::cout << "shape1_mat set to " << val << std::endl;
    shape1_mat = temp_mat;
    std::cout << "shape1_mat->GetName() = " << shape1_mat->GetName() << std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2XPos(G4double val)
{
  fScat2XPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2YPos(G4double val)
{
  fScat2YPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2ZPos(G4double val)
{
  fScat2ZPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2PolarR(G4double val)
{
  fScat2PolarR = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2PolarPhi(G4double val)
{
  fScat2PolarPhi = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2PolarTheta(G4double val)
{
  fScat2PolarTheta = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2RotX(G4double val)
{
  fScat2RotX = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2RotY(G4double val)
{
  fScat2RotY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2RotZ(G4double val)
{
  fScat2RotZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2Rad(G4double val)
{
  fScat2Rad = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2Height(G4double val)
{
  fScat2Height = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2Mat(G4String val)
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* temp_mat = nist->FindOrBuildMaterial(val);
  if(temp_mat == NULL){
    std::cout << val << " does not exist in NIST" << std::endl;
    assert(temp_mat != NULL);
  }
  else {
    std::cout << "shape3_mat set to " << val << std::endl;
    shape3_mat = temp_mat;
    std::cout << "shape3_mat->GetName() = " << shape3_mat->GetName() << std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2Bool(G4bool val)
{
  fScat2Bool = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetScat2Switch(G4bool val)
{
  if (val == false)
    {
      fScat2Type = "Absorber";
      fScat2CopyNo = 2;
    }
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

void B1DetectorConstruction::SetDetZPos(G4double val)
{
  fDetZPos = val;
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

void B1DetectorConstruction::SetDetPolarTheta(G4double val)
{
  fDetPolarTheta = val;
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

void B1DetectorConstruction::SetDetMat(G4String val)
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* temp_mat = nist->FindOrBuildMaterial(val);
  if(temp_mat == NULL){
    std::cout << val << " does not exist in NIST" << std::endl;
    assert(temp_mat != NULL);
  }
  else {
    std::cout << "shape2_mat set to " << val << std::endl;
    shape2_mat = temp_mat;
    std::cout << "shape2_mat->GetName() = " << shape2_mat->GetName() << std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2XPos(G4double val)
{
  fDet2XPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2YPos(G4double val)
{
  fDet2YPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2ZPos(G4double val)
{
  fDet2ZPos = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2PolarR(G4double val)
{
  fDet2PolarR = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2PolarPhi(G4double val)
{
  fDet2PolarPhi = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2PolarTheta(G4double val)
{
  fDet2PolarTheta = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2RotX(G4double val)
{
  fDet2RotX = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2RotY(G4double val)
{
  fDet2RotY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2RotZ(G4double val)
{
  fDet2RotZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2Rad(G4double val)
{
  fDet2Rad = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2Height(G4double val)
{
  fDet2Height = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2Mat(G4String val)
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* temp_mat = nist->FindOrBuildMaterial(val);
  if(temp_mat == NULL){
    std::cout << val << " does not exist in NIST" << std::endl;
    assert(temp_mat != NULL);
  }
  else {
    std::cout << " set to " << val << std::endl;
    shape4_mat = temp_mat;
    std::cout << "shape4_mat->GetName() = " << shape4_mat->GetName() << std::endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2Bool(G4bool val)
{
  fDet2Bool = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::SetDet2Switch(G4bool val)
{
  if (val == false)
    {
      fDet2Type = "Scatterer";
      fDet2CopyNo = 2;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
