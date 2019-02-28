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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{ // track runtime of event by summing delta times in each step - by Jack
  G4double deltaTime = step->GetDeltaTime();
  fEventAction->TotalTime(deltaTime);
  // get physical and logical volume of the current step
  G4VPhysicalVolume* volumePhys 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume();
  G4LogicalVolume* volume = volumePhys->GetLogicalVolume();

  // Code to fix Segmentation Faults
  G4String procName = "";
  G4StepPoint* preStep = step->GetPreStepPoint();
  if(preStep!= nullptr)
  {
    const G4VProcess* proc = preStep->GetProcessDefinedStep();
    if(proc != nullptr)
    {
      procName = proc->GetProcessName();
    }
  }


  // dose in body George
  if (volume->GetName() == "Body")
 	{G4double edepStep = step->GetTotalEnergyDeposit();
      	fEventAction->AddEdepBody(edepStep);} //Have a look at AddEdepDetector. 

 // End of dose in body George


  // check if we are in scoring volume
  if (volume->GetName() != "Scatterer" && volume->GetName() != "Absorber") return;

  // collect energy deposited in step - originally by Jack, generalised with copy number by Douglas
  // scatterer energy
  // get copy number if multiple scatter detectors
  if (volume->GetName() == "Scatterer")
      
    { G4double stepLength = step->GetStepLength(); //George checking step length
      //std::cout<<"stepLength="<<stepLength/mm<< "\n"; //George checking step length
      G4double edepStep = step->GetTotalEnergyDeposit();
      int copyNo = volumePhys->GetCopyNo();
      fEventAction->AddEdepScatterer(edepStep, copyNo);
	if (procName == "compt")
	{
		G4double timeScatterer = step->GetTrack()->GetGlobalTime();
		G4ThreeVector Pos = step->GetPreStepPoint()->GetPosition();
		fEventAction->TimeScatterer(timeScatterer, copyNo);
		fEventAction->Vector(Pos);
		fEventAction->Count();
	}
	// Finding details of processes that cause energy deposition that aren't Compton to solve the 0 scatter coincidence problem - by Jack
	else
	{
	  if (edepStep!=0)
	    {
	      fEventAction->ZeroScatterInfo(edepStep, procName, step->GetPreStepPoint()->GetPosition());
	    }
	}
    }

  // detector energy
  if (volume->GetName() == "Absorber")
    {
      G4double edepStep = step->GetTotalEnergyDeposit();
      int copyNo = volumePhys->GetCopyNo();
      G4double timeDetector = step->GetTrack()->GetGlobalTime();
      fEventAction->AddEdepDetector(edepStep, copyNo);
      fEventAction->TimeDetector(timeDetector, copyNo);
	if (procName == "phot")
	{
		G4ThreeVector Pos2 = step->GetPreStepPoint()->GetPosition();
		fEventAction->Vector2(Pos2);	
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

