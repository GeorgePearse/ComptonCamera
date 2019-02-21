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
/// \file B1PrimaryGeneratorMessenger.cc
/// \brief Implementation of the B1PrimaryGeneratorMessenger class
//
//
// BY JACK
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B1PrimaryGeneratorMessenger.hh"

#include "B1PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorMessenger::B1PrimaryGeneratorMessenger(
                                                  B1PrimaryGeneratorAction* Gun)
:fAction(Gun), fGunDir(0), fXPos(0), fYPos(0)
{
  fGunDir = new G4UIdirectory("/B1/source");
  fGunDir->SetGuidance("source control");
  
  fXPos = new G4UIcmdWithADoubleAndUnit("/B1/source/setXPos",this);
  fXPos->SetGuidance("set X pos of the point source");
  fXPos->SetParameterName("xPos",false);
  fXPos->SetUnitCategory("Length");

  fYPos = new G4UIcmdWithADoubleAndUnit("/B1/source/setYPos",this);
  fYPos->SetGuidance("set Y pos of the point source");
  fYPos->SetParameterName("YPos",false);
  fYPos->SetUnitCategory("Length");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorMessenger::~B1PrimaryGeneratorMessenger()
{
  delete fXPos;
  delete fYPos;
  delete fGunDir;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{    
  if (command == fXPos)
   { fAction->SetXPos(fXPos->GetNewDoubleValue(newValue));}

  if (command == fYPos)
   { fAction->SetYPos(fYPos->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

