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
// BY JACK
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B1DetectorMessenger.hh"

#include "B1DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorMessenger::B1DetectorMessenger(B1DetectorConstruction * det)
:fDetector(det),
 fB1Dir(0),
 fScatDir(0),
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
 fScatMat(0),
 fScat2Dir(0),
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
 fScat2Mat(0),
 fScat2Bool(0),
 fScat2Switch(0),
 fDetDir(0),
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
 fDetMat(0),
 fDet2Dir(0),
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
 fDet2Mat(0),
 fDet2Bool(0),
 fDet2Switch(0),
 fUpdateCmd(0)
{
  // Sets directory for these commands
  fB1Dir = new G4UIdirectory("/B1/");
  fB1Dir->SetGuidance("UI commands for the Compton camera");
  
  // Sets directory for scatterer commands
  fScatDir = new G4UIdirectory("/B1/scat/");
  fScatDir->SetGuidance("scatterer construction commands");
  
  // Set scatterer position in Cartesian coordinates
  fScatXPos = new G4UIcmdWithADoubleAndUnit("/B1/scat/setX", this);
  fScatXPos->SetGuidance("Set X position of the scatterer (zero is in line with source)");
  fScatXPos->SetParameterName("ScatXPos", false);  
  fScatXPos->SetUnitCategory("Length");

  fScatYPos = new G4UIcmdWithADoubleAndUnit("/B1/scat/setY", this);
  fScatYPos->SetGuidance("Set Y position of the scatterer");
  fScatYPos->SetParameterName("ScatYPos", false);  
  fScatYPos->SetUnitCategory("Length");

  fScatZPos = new G4UIcmdWithADoubleAndUnit("/B1/scat/setZ", this);
  fScatZPos->SetGuidance("Set Z position of the scatterer");
  fScatZPos->SetParameterName("ScatZPos", false);  
  fScatZPos->SetUnitCategory("Length");
  
  // Set scatterer position in polar coordinates
  fScatPolarR = new G4UIcmdWithADoubleAndUnit("/B1/scat/setPolarR", this);
  fScatPolarR->SetGuidance("Set R coordinate of the scatterer (polar coordinates)");
  fScatPolarR->SetParameterName("ScatPolarR", false);  
  fScatPolarR->SetUnitCategory("Length");

  fScatPolarPhi = new G4UIcmdWithADoubleAndUnit("/B1/scat/setPolarPhi", this);
  fScatPolarPhi->SetGuidance("Set phi coordinate of the scatterer (polar coordinates)");
  fScatPolarPhi->SetParameterName("ScatPolarPhi", false);  
  fScatPolarPhi->SetUnitCategory("Angle");

  fScatPolarTheta = new G4UIcmdWithADoubleAndUnit("/B1/scat/setPolarTheta", this);
  fScatPolarTheta->SetGuidance("Set theta coordinate of the scatterer (polar coordinates)");
  fScatPolarTheta->SetParameterName("ScatPolarTheta", false);  
  fScatPolarTheta->SetUnitCategory("Angle");

  // Set scatterer rotation
  fScatRotX = new G4UIcmdWithADoubleAndUnit("/B1/scat/setRotX", this);
  fScatRotX->SetGuidance("Set X rotation of scatterer");
  fScatRotX->SetParameterName("ScatRotX", false);  
  fScatRotX->SetUnitCategory("Angle");

  fScatRotY = new G4UIcmdWithADoubleAndUnit("/B1/scat/setRotY", this);
  fScatRotY->SetGuidance("Set Y rotation of scatterer");
  fScatRotY->SetParameterName("ScatRotY", false);  
  fScatRotY->SetUnitCategory("Angle");

  fScatRotZ = new G4UIcmdWithADoubleAndUnit("/B1/scat/setRotZ", this);
  fScatRotZ->SetGuidance("Set Z rotation of scatterer");
  fScatRotZ->SetParameterName("ScatRotZ", false);  
  fScatRotZ->SetUnitCategory("Angle");

  // Set radius and height of scatterer
  fScatRad = new G4UIcmdWithADoubleAndUnit("/B1/scat/setRad", this);
  fScatRad->SetGuidance("Set radius of the scatterer");
  fScatRad->SetParameterName("ScatRad", false);
  fScatRad->SetUnitCategory("Length");

  fScatHeight = new G4UIcmdWithADoubleAndUnit("/B1/scat/setHeight", this);
  fScatHeight->SetGuidance("Set height of the scatterer");
  fScatHeight->SetParameterName("ScatHeight", false);
  fScatHeight->SetUnitCategory("Length");

  // Set material of scatterer
  fScatMat = new G4UIcmdWithAString("/B1/scat/setMat", this);
  fScatMat->SetGuidance("Set material of scatterer");
  fScatMat->SetParameterName("ScatMat", false);

  //-----------------------------------------------------
  // Sets directory for scatterer 2 commands
  fScat2Dir = new G4UIdirectory("/B1/scat2/");
  fScat2Dir->SetGuidance("scatterer 2 construction commands");
  
  // Set scatterer 2 position in Cartesian coordinates
  fScat2XPos = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setX", this);
  fScat2XPos->SetGuidance("Set X position of the scatterer (zero is in line with source)");
  fScat2XPos->SetParameterName("Scat2XPos", false);  
  fScat2XPos->SetUnitCategory("Length");

  fScat2YPos = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setY", this);
  fScat2YPos->SetGuidance("Set Y position of the scatterer");
  fScat2YPos->SetParameterName("Scat2YPos", false);  
  fScat2YPos->SetUnitCategory("Length");

  fScat2ZPos = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setZ", this);
  fScat2ZPos->SetGuidance("Set Z position of the scatterer");
  fScat2ZPos->SetParameterName("Scat2ZPos", false);  
  fScat2ZPos->SetUnitCategory("Length");
  
  // Set scatterer 2 position in polar coordinates
  fScat2PolarR = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setPolarR", this);
  fScat2PolarR->SetGuidance("Set R coordinate of the scatterer (polar coordinates)");
  fScat2PolarR->SetParameterName("Scat2PolarR", false);  
  fScat2PolarR->SetUnitCategory("Length");

  fScat2PolarPhi = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setPolarPhi", this);
  fScat2PolarPhi->SetGuidance("Set phi coordinate of the scatterer (polar coordinates)");
  fScat2PolarPhi->SetParameterName("Scat2PolarPhi", false);  
  fScat2PolarPhi->SetUnitCategory("Angle");

  fScat2PolarTheta = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setPolarTheta", this);
  fScat2PolarTheta->SetGuidance("Set theta coordinate of the scatterer (polar coordinates)");
  fScat2PolarTheta->SetParameterName("Scat2PolarTheta", false);  
  fScat2PolarTheta->SetUnitCategory("Angle");

  // Set scatterer 2 rotation
  fScat2RotX = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setRotX", this);
  fScat2RotX->SetGuidance("Set X rotation of scatterer");
  fScat2RotX->SetParameterName("Scat2RotX", false);  
  fScat2RotX->SetUnitCategory("Angle");

  fScat2RotY = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setRotY", this);
  fScat2RotY->SetGuidance("Set Y rotation of scatterer");
  fScat2RotY->SetParameterName("Scat2RotY", false);  
  fScat2RotY->SetUnitCategory("Angle");

  fScat2RotZ = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setRotZ", this);
  fScat2RotZ->SetGuidance("Set Z rotation of scatterer");
  fScat2RotZ->SetParameterName("Scat2RotZ", false);  
  fScat2RotZ->SetUnitCategory("Angle");

  // Set radius and height of scatterer 2
  fScat2Rad = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setRad", this);
  fScat2Rad->SetGuidance("Set radius of the scatterer");
  fScat2Rad->SetParameterName("Scat2Rad", false);
  fScat2Rad->SetUnitCategory("Length");

  fScat2Height = new G4UIcmdWithADoubleAndUnit("/B1/scat2/setHeight", this);
  fScat2Height->SetGuidance("Set height of the scatterer");
  fScat2Height->SetParameterName("Scat2Height", false);
  fScat2Height->SetUnitCategory("Length");

  // Set material of scatterer 2
  fScat2Mat = new G4UIcmdWithAString("/B1/scat2/setMat", this);
  fScat2Mat->SetGuidance("Set material of scatterer 2");
  fScat2Mat->SetParameterName("Scat2Mat", false);

  // Set whether scatterer 2 is active
  fScat2Bool = new G4UIcmdWithABool("/B1/scat2/setBool", this);
  fScat2Bool->SetGuidance("Set whether scatterer 2 is active or not");
  fScat2Bool->SetParameterName("Scat2Bool", false);

  // Switch scatterer 2 between being a scatterer and an absorber
  fScat2Switch = new G4UIcmdWithABool("/B1/scat2/setSwitch", this);
  fScat2Switch->SetGuidance("Set whether scatterer 2 is a scatterer or an absorber - false switches it to be an absorber");
  fScat2Switch->SetParameterName("Scat2Switch", false);
  
  //-----------------------------------------------------

  // Set directory for detector commands
  fDetDir = new G4UIdirectory("/B1/det/");
  fDetDir->SetGuidance("detector construction commands");

  // Set detector position in Cartesian coordinates
  fDetXPos = new G4UIcmdWithADoubleAndUnit("/B1/det/setX", this);
  fDetXPos->SetGuidance("Set X position of the detector (zero is in line with source)");
  fDetXPos->SetParameterName("DetXPos", false);
  fDetXPos->SetUnitCategory("Length");

  fDetYPos = new G4UIcmdWithADoubleAndUnit("/B1/det/setY", this);
  fDetYPos->SetGuidance("Set Y position of the detector");
  fDetYPos->SetParameterName("DetYPos", false);
  fDetYPos->SetUnitCategory("Length");

  fDetZPos = new G4UIcmdWithADoubleAndUnit("/B1/det/setZ", this);
  fDetZPos->SetGuidance("Set Z position of the detector");
  fDetZPos->SetParameterName("DetZPos", false);
  fDetZPos->SetUnitCategory("Length");

  // Set detector position in polar coordinates
  fDetPolarR = new G4UIcmdWithADoubleAndUnit("/B1/det/setPolarR", this);
  fDetPolarR->SetGuidance("Set R coordinate of the detector (polar coordinates)");
  fDetPolarR->SetParameterName("DetPolarR", false);  
  fDetPolarR->SetUnitCategory("Length");

  fDetPolarPhi = new G4UIcmdWithADoubleAndUnit("/B1/det/setPolarPhi", this);
  fDetPolarPhi->SetGuidance("Set phi coordinate of the detector (polar coordinates)");
  fDetPolarPhi->SetParameterName("DetPolarPhi", false);  
  fDetPolarPhi->SetUnitCategory("Angle");

  fDetPolarTheta = new G4UIcmdWithADoubleAndUnit("/B1/det/setPolarTheta", this);
  fDetPolarTheta->SetGuidance("Set theta coordinate of the detector (polar coordinates)");
  fDetPolarTheta->SetParameterName("DetPolarTheta", false);  
  fDetPolarTheta->SetUnitCategory("Angle");

  // Set detector rotation
  fDetRotX = new G4UIcmdWithADoubleAndUnit("/B1/det/setRotX", this);
  fDetRotX->SetGuidance("Set X rotation of detector");
  fDetRotX->SetParameterName("DetRotX", false);  
  fDetRotX->SetUnitCategory("Angle");

  fDetRotY = new G4UIcmdWithADoubleAndUnit("/B1/det/setRotY", this);
  fDetRotY->SetGuidance("Set Y rotation of detector");
  fDetRotY->SetParameterName("DetRotY", false);  
  fDetRotY->SetUnitCategory("Angle");

  fDetRotZ = new G4UIcmdWithADoubleAndUnit("/B1/det/setRotZ", this);
  fDetRotZ->SetGuidance("Set Z rotation of detector");
  fDetRotZ->SetParameterName("DetRotZ", false);  
  fDetRotZ->SetUnitCategory("Angle");

  // Set radius and height of detector
  fDetRad = new G4UIcmdWithADoubleAndUnit("/B1/det/setRad", this);
  fDetRad->SetGuidance("Set radius of the detector");
  fDetRad->SetParameterName("DetRad", false);
  fDetRad->SetUnitCategory("Length");

  fDetHeight = new G4UIcmdWithADoubleAndUnit("/B1/det/setHeight", this);
  fDetHeight->SetGuidance("Set height of the detector");
  fDetHeight->SetParameterName("DetHeight", false);
  fDetHeight->SetUnitCategory("Length");

  // Set material of detector
  fDetMat = new G4UIcmdWithAString("/B1/det/setMat", this);
  fDetMat->SetGuidance("Set material of detector");
  fDetMat->SetParameterName("DetMat", false);

  //------------------------------------------------------------

  // Set directory for detector 2 commands
  fDet2Dir = new G4UIdirectory("/B1/det2/");
  fDet2Dir->SetGuidance("detector 2 construction commands");

  // Set detector 2 position in Cartesian coordinates
  fDet2XPos = new G4UIcmdWithADoubleAndUnit("/B1/det2/setX", this);
  fDet2XPos->SetGuidance("Set X position of the detector (zero is in line with source)");
  fDet2XPos->SetParameterName("Det2XPos", false);
  fDet2XPos->SetUnitCategory("Length");

  fDet2YPos = new G4UIcmdWithADoubleAndUnit("/B1/det2/setY", this);
  fDet2YPos->SetGuidance("Set Y position of the detector");
  fDet2YPos->SetParameterName("Det2YPos", false);
  fDet2YPos->SetUnitCategory("Length");

  fDet2ZPos = new G4UIcmdWithADoubleAndUnit("/B1/det2/setZ", this);
  fDet2ZPos->SetGuidance("Set Z position of the detector");
  fDet2ZPos->SetParameterName("Det2ZPos", false);
  fDet2ZPos->SetUnitCategory("Length");

  // Set detector 2 position in polar coordinates
  fDet2PolarR = new G4UIcmdWithADoubleAndUnit("/B1/det2/setPolarR", this);
  fDet2PolarR->SetGuidance("Set R coordinate of the detector (polar coordinates)");
  fDet2PolarR->SetParameterName("Det2PolarR", false);  
  fDet2PolarR->SetUnitCategory("Length");

  fDet2PolarPhi = new G4UIcmdWithADoubleAndUnit("/B1/det2/setPolarPhi", this);
  fDet2PolarPhi->SetGuidance("Set phi coordinate of the detector (polar coordinates)");
  fDet2PolarPhi->SetParameterName("Det2PolarPhi", false);  
  fDet2PolarPhi->SetUnitCategory("Angle");

  fDet2PolarTheta = new G4UIcmdWithADoubleAndUnit("/B1/det2/setPolarTheta", this);
  fDet2PolarTheta->SetGuidance("Set theta coordinate of the detector (polar coordinates)");
  fDet2PolarTheta->SetParameterName("Det2PolarTheta", false);  
  fDet2PolarTheta->SetUnitCategory("Angle");

  // Set detector 2 rotation
  fDet2RotX = new G4UIcmdWithADoubleAndUnit("/B1/det2/setRotX", this);
  fDet2RotX->SetGuidance("Set X rotation of detector");
  fDet2RotX->SetParameterName("Det2RotX", false);  
  fDet2RotX->SetUnitCategory("Angle");

  fDet2RotY = new G4UIcmdWithADoubleAndUnit("/B1/det2/setRotY", this);
  fDet2RotY->SetGuidance("Set Y rotation of detector");
  fDet2RotY->SetParameterName("Det2RotY", false);  
  fDet2RotY->SetUnitCategory("Angle");

  fDet2RotZ = new G4UIcmdWithADoubleAndUnit("/B1/det2/setRotZ", this);
  fDet2RotZ->SetGuidance("Set Z rotation of detector");
  fDet2RotZ->SetParameterName("Det2RotZ", false);  
  fDet2RotZ->SetUnitCategory("Angle");

  // Set radius and height of detector 2
  fDet2Rad = new G4UIcmdWithADoubleAndUnit("/B1/det2/setRad", this);
  fDet2Rad->SetGuidance("Set radius of the detector");
  fDet2Rad->SetParameterName("Det2Rad", false);
  fDet2Rad->SetUnitCategory("Length");

  fDet2Height = new G4UIcmdWithADoubleAndUnit("/B1/det2/setHeight", this);
  fDet2Height->SetGuidance("Set height of the detector");
  fDet2Height->SetParameterName("Det2Height", false);
  fDet2Height->SetUnitCategory("Length");

  // Set material of detector 2
  fDet2Mat = new G4UIcmdWithAString("/B1/det2/setMat", this);
  fDet2Mat->SetGuidance("Set material of detector");
  fDet2Mat->SetParameterName("Det2Mat", false);

  // Set whether detector 2 is active
  fDet2Bool = new G4UIcmdWithABool("/B1/det2/setBool", this);
  fDet2Bool->SetGuidance("Set whether detector 2 is active or not");
  fDet2Bool->SetParameterName("Det2Bool", false);

  // Switch absorber 2 between being an absorber and a scatterer
  fDet2Switch = new G4UIcmdWithABool("/B1/det2/setSwitch", this);
  fDet2Switch->SetGuidance("Set whether absorber 2 is a scatterer or an absorber - false switches it to be a scatterer");
  fDet2Switch->SetParameterName("Det2Switch", false);
  
  // Update geometry to apply changes made using the messenger
  fUpdateCmd = new G4UIcmdWithoutParameter("/B1/update",this);
  fUpdateCmd->SetGuidance("Update geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorMessenger::~B1DetectorMessenger()
{
  delete fScatXPos;
  delete fScatYPos;
  delete fScatZPos;
  delete fScatPolarR;
  delete fScatPolarPhi;
  delete fScatPolarTheta;
  delete fScatRotX;
  delete fScatRotY;
  delete fScatRotZ;
  delete fScatRad;
  delete fScatHeight;
  delete fScatMat;

  delete fScat2XPos;
  delete fScat2YPos;
  delete fScat2ZPos;
  delete fScat2PolarR;
  delete fScat2PolarPhi;
  delete fScat2PolarTheta;
  delete fScat2RotX;
  delete fScat2RotY;
  delete fScat2RotZ;
  delete fScat2Rad;
  delete fScat2Height;
  delete fScat2Mat;
  delete fScat2Bool;
  delete fScat2Switch;
  
  delete fDetXPos;
  delete fDetYPos;
  delete fDetZPos;
  delete fDetPolarR;
  delete fDetPolarPhi;
  delete fDetPolarTheta;
  delete fDetRotX;
  delete fDetRotY;
  delete fDetRotZ;
  delete fDetRad;
  delete fDetHeight;
  delete fDetMat;
  
  delete fDet2XPos;
  delete fDet2YPos;
  delete fDet2ZPos;
  delete fDet2PolarR;
  delete fDet2PolarPhi;
  delete fDet2PolarTheta;
  delete fDet2RotX;
  delete fDet2RotY;
  delete fDet2RotZ;
  delete fDet2Rad;
  delete fDet2Height;
  delete fDet2Mat;
  delete fDet2Bool;
  delete fDet2Switch;
  
  delete fUpdateCmd;
  delete fScatDir;  
  delete fB1Dir;
  delete fDetDir;
  delete fDet2Dir;
  delete fScat2Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  // Sets relevant values based on which command is called
  if ( command == fScatXPos )
   {fDetector->SetScatXPos(fScatXPos->GetNewDoubleValue(newValue));}

  if ( command == fScatYPos )
   {fDetector->SetScatYPos(fScatYPos->GetNewDoubleValue(newValue));}

  if ( command == fScatZPos )
   {fDetector->SetScatZPos(fScatZPos->GetNewDoubleValue(newValue));}

  if ( command == fScatPolarR )
   {fDetector->SetScatPolarR(fScatPolarR->GetNewDoubleValue(newValue));}

  if ( command == fScatPolarPhi )
   {fDetector->SetScatPolarPhi(fScatPolarPhi->GetNewDoubleValue(newValue));}

  if ( command == fScatPolarTheta )
   {fDetector->SetScatPolarTheta(fScatPolarTheta->GetNewDoubleValue(newValue));}

  if ( command == fScatRotX )
   {fDetector->SetScatRotX(fScatRotX->GetNewDoubleValue(newValue));}

  if ( command == fScatRotY )
   {fDetector->SetScatRotY(fScatRotY->GetNewDoubleValue(newValue));}

  if ( command == fScatRotZ )
   {fDetector->SetScatRotZ(fScatRotZ->GetNewDoubleValue(newValue));}

  if ( command == fScatRad )
   {fDetector->SetScatRad(fScatRad->GetNewDoubleValue(newValue));}

  if ( command == fScatMat )
   {fDetector->SetScatMat(newValue);}

  if ( command == fScat2XPos )
   {fDetector->SetScat2XPos(fScat2XPos->GetNewDoubleValue(newValue));}

  if ( command == fScat2YPos )
   {fDetector->SetScat2YPos(fScat2YPos->GetNewDoubleValue(newValue));}

  if ( command == fScat2ZPos )
   {fDetector->SetScat2ZPos(fScat2ZPos->GetNewDoubleValue(newValue));}

  if ( command == fScat2PolarR )
   {fDetector->SetScat2PolarR(fScat2PolarR->GetNewDoubleValue(newValue));}

  if ( command == fScat2PolarPhi )
   {fDetector->SetScat2PolarPhi(fScat2PolarPhi->GetNewDoubleValue(newValue));}

  if ( command == fScat2PolarTheta )
   {fDetector->SetScat2PolarTheta(fScat2PolarTheta->GetNewDoubleValue(newValue));}

  if ( command == fScat2RotX )
   {fDetector->SetScat2RotX(fScat2RotX->GetNewDoubleValue(newValue));}

  if ( command == fScat2RotY )
   {fDetector->SetScat2RotY(fScat2RotY->GetNewDoubleValue(newValue));}

  if ( command == fScat2RotZ )
   {fDetector->SetScat2RotZ(fScat2RotZ->GetNewDoubleValue(newValue));}

  if ( command == fScat2Rad )
   {fDetector->SetScat2Rad(fScat2Rad->GetNewDoubleValue(newValue));}

  if ( command == fScat2Height )
   {fDetector->SetScat2Height(fScat2Height->GetNewDoubleValue(newValue));}

  if ( command == fScat2Mat )
   {fDetector->SetScat2Mat(newValue);}

  if ( command == fScat2Bool )
   {fDetector->SetScat2Bool(fScat2Bool->GetNewBoolValue(newValue));}

  if ( command == fScat2Switch )
   {fDetector->SetScat2Switch(fScat2Switch->GetNewBoolValue(newValue));}

  if ( command == fDetXPos )
   {fDetector->SetDetXPos(fDetXPos->GetNewDoubleValue(newValue));}

  if ( command == fDetYPos )
   {fDetector->SetDetYPos(fDetYPos->GetNewDoubleValue(newValue));}

  if ( command == fDetZPos )
   {fDetector->SetDetZPos(fDetZPos->GetNewDoubleValue(newValue));}

  if ( command == fDetPolarR )
   {fDetector->SetDetPolarR(fDetPolarR->GetNewDoubleValue(newValue));}

  if ( command == fDetPolarPhi )
   {fDetector->SetDetPolarPhi(fDetPolarPhi->GetNewDoubleValue(newValue));}

  if ( command == fDetPolarTheta )
   {fDetector->SetDetPolarTheta(fDetPolarTheta->GetNewDoubleValue(newValue));}

  if ( command == fDetRotX )
   {fDetector->SetDetRotX(fDetRotX->GetNewDoubleValue(newValue));}

  if ( command == fDetRotY )
   {fDetector->SetDetRotY(fDetRotY->GetNewDoubleValue(newValue));}

  if ( command == fDetRotZ )
   {fDetector->SetDetRotZ(fDetRotZ->GetNewDoubleValue(newValue));}

  if ( command == fDetRad )
   {fDetector->SetDetRad(fDetRad->GetNewDoubleValue(newValue));}

  if ( command == fDetHeight )
   {fDetector->SetDetHeight(fDetHeight->GetNewDoubleValue(newValue));}

  if ( command == fDetMat )
   {fDetector->SetDetMat(newValue);}

  if ( command == fDet2XPos )
   {fDetector->SetDet2XPos(fDet2XPos->GetNewDoubleValue(newValue));}

  if ( command == fDet2YPos )
   {fDetector->SetDet2YPos(fDet2YPos->GetNewDoubleValue(newValue));}

  if ( command == fDet2ZPos )
   {fDetector->SetDet2ZPos(fDet2ZPos->GetNewDoubleValue(newValue));}

  if ( command == fDet2PolarR )
   {fDetector->SetDet2PolarR(fDet2PolarR->GetNewDoubleValue(newValue));}

  if ( command == fDet2PolarPhi )
   {fDetector->SetDet2PolarPhi(fDet2PolarPhi->GetNewDoubleValue(newValue));}

  if ( command == fDet2PolarTheta )
   {fDetector->SetDet2PolarTheta(fDet2PolarTheta->GetNewDoubleValue(newValue));}

  if ( command == fDet2RotX )
   {fDetector->SetDet2RotX(fDet2RotX->GetNewDoubleValue(newValue));}

  if ( command == fDet2RotY )
   {fDetector->SetDet2RotY(fDet2RotY->GetNewDoubleValue(newValue));}

  if ( command == fDet2RotZ )
   {fDetector->SetDet2RotZ(fDet2RotZ->GetNewDoubleValue(newValue));}

  if ( command == fDet2Rad )
   {fDetector->SetDet2Rad(fDet2Rad->GetNewDoubleValue(newValue));}

  if ( command == fDet2Height )
   {fDetector->SetDet2Height(fDet2Height->GetNewDoubleValue(newValue));}

  if ( command == fDet2Mat )
   {fDetector->SetDet2Mat(newValue);}

  if ( command == fDet2Bool )
   {fDetector->SetDet2Bool(fDet2Bool->GetNewBoolValue(newValue));}

  if ( command == fDet2Switch )
   {fDetector->SetDet2Switch(fDet2Switch->GetNewBoolValue(newValue));}
   
  if  ( command == fUpdateCmd )
   {fDetector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
