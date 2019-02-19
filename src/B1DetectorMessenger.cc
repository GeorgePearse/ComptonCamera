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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorMessenger::B1DetectorMessenger(B1DetectorConstruction * det)
:fDetector(det),
 fB1Dir(0), fScatDir(0), fScatXPos(0), fScatYPos(0), fScatRad(0), fDetDir(0), fDetXPos(0), fDetYPos(0), fUpdateCmd(0)
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
  
  // Set scatterer position in polar coordinates
  fScatPolarR = new G4UIcmdWithADoubleAndUnit("/B1/scat/setPolarR", this);
  fScatPolarR->SetGuidance("Set R coordinate of the scatterer (polar coordinates)");
  fScatPolarR->SetParameterName("ScatPolarR", false);  
  fScatPolarR->SetUnitCategory("Length");

  fScatPolarPhi = new G4UIcmdWithADoubleAndUnit("/B1/scat/setPolarPhi", this);
  fScatPolarPhi->SetGuidance("Set phi coordinate of the scatterer (polar coordinates)");
  fScatPolarPhi->SetParameterName("ScatPolarPhi", false);  
  fScatPolarPhi->SetUnitCategory("Angle");

  // Set radius and height of scatterer
  fScatRad = new G4UIcmdWithADoubleAndUnit("/B1/scat/setRad", this);
  fScatRad->SetGuidance("Set radius of the scatterer");
  fScatRad->SetParameterName("ScatRad", false);
  fScatRad->SetUnitCategory("Length");

  fScatHeight = new G4UIcmdWithADoubleAndUnit("/B1/scat/setHeight", this);
  fScatHeight->SetGuidance("Set height of the scatterer");
  fScatHeight->SetParameterName("ScatHeight", false);
  fScatHeight->SetUnitCategory("Length");
  
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

  // Set detector position in polar coordinates
  fDetPolarR = new G4UIcmdWithADoubleAndUnit("/B1/det/setPolarR", this);
  fDetPolarR->SetGuidance("Set R coordinate of the detector (polar coordinates)");
  fDetPolarR->SetParameterName("DetPolarR", false);  
  fDetPolarR->SetUnitCategory("Length");

  fDetPolarPhi = new G4UIcmdWithADoubleAndUnit("/B1/det/setPolarPhi", this);
  fDetPolarPhi->SetGuidance("Set phi coordinate of the detector (polar coordinates)");
  fDetPolarPhi->SetParameterName("DetPolarPhi", false);  
  fDetPolarPhi->SetUnitCategory("Angle");

  // Set radius and height of detector
  fDetRad = new G4UIcmdWithADoubleAndUnit("/B1/det/setRad", this);
  fDetRad->SetGuidance("Set radius of the detector");
  fDetRad->SetParameterName("DetRad", false);
  fDetRad->SetUnitCategory("Length");

  fDetHeight = new G4UIcmdWithADoubleAndUnit("/B1/det/setHeight", this);
  fDetHeight->SetGuidance("Set height of the detector");
  fDetHeight->SetParameterName("DetHeight", false);
  fDetHeight->SetUnitCategory("Length");
  
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
  delete fScatPolarR;
  delete fScatPolarPhi;
  delete fScatRad;
  delete fScatHeight;
  delete fDetXPos;
  delete fDetYPos;
  delete fDetPolarR;
  delete fDetPolarPhi;
  delete fDetRad;
  delete fDetHeight;
  delete fUpdateCmd;
  delete fScatDir;  
  delete fB1Dir;
  delete fDetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  // Sets relevant values based on which command is called
  if ( command == fScatXPos )
   {fDetector->SetScatXPos(fScatXPos->GetNewDoubleValue(newValue));}

  if ( command == fScatYPos )
   {fDetector->SetScatYPos(fScatYPos->GetNewDoubleValue(newValue));}

  if ( command == fScatPolarR )
   {fDetector->SetScatPolarR(fScatPolarR->GetNewDoubleValue(newValue));}

  if ( command == fScatPolarPhi )
   {fDetector->SetScatPolarPhi(fScatPolarPhi->GetNewDoubleValue(newValue));}

  if ( command == fScatRad )
   {fDetector->SetScatRad(fScatRad->GetNewDoubleValue(newValue));}

  if ( command == fScatHeight )
   {fDetector->SetScatHeight(fScatHeight->GetNewDoubleValue(newValue));}

  if ( command == fDetXPos )
   {fDetector->SetDetXPos(fDetXPos->GetNewDoubleValue(newValue));}

  if ( command == fDetYPos )
   {fDetector->SetDetYPos(fDetYPos->GetNewDoubleValue(newValue));}

  if ( command == fDetPolarR )
   {fDetector->SetDetPolarR(fDetPolarR->GetNewDoubleValue(newValue));}

  if ( command == fDetPolarPhi )
   {fDetector->SetDetPolarPhi(fDetPolarPhi->GetNewDoubleValue(newValue));}

  if ( command == fDetRad )
   {fDetector->SetDetRad(fDetRad->GetNewDoubleValue(newValue));}

  if ( command == fDetHeight )
   {fDetector->SetDetHeight(fDetHeight->GetNewDoubleValue(newValue));}
   
  if  ( command == fUpdateCmd )
   {fDetector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
