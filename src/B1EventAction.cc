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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdepScatterer(0.),
  fEdepDetector(0.),
  fRunTime(0.)
{
fFirstWrite = true;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdepScatterer = 0.;
  fEdepDetector = 0.;
  fBeginTime = fRunTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event*)
{   
  if(fEdepScatterer != 0 && fEdepDetector != 0)
    {
      // Text file writer for Scatterer
      std::ofstream myfile;
      // Special condition for first write to create file
      if(fFirstWrite)
	{
		myfile.open("scatterdata.txt");
	}
	else
	{
		myfile.open ("scatterdata.txt", std::ios::app);
	}
      	if (myfile.is_open())
      	{
          myfile << (fTimeScatterer + fBeginTime)/1000  << " "
	  << fEdepScatterer*1000 << "\n";
	  myfile.close();
        }
	else std::cerr << "Unable to open scatter file" << std::endl;


       // Text file writer for Absorber
       std::ofstream myfile2;
       // Special condition for first write to create file
       if(fFirstWrite)
	{
		myfile2.open("absorbdata.txt");
	}
	else
	{
		myfile2.open ("absorbdata.txt", std::ios::app);
	}
      	if (myfile2.is_open())
      	{
          myfile2 << (fTimeDetector + fBeginTime)/1000 << " "
	  << fEdepDetector*1000 << "\n";
	  myfile2.close();
        }
	else std::cerr << "Unable to open absorb file" << std::endl;

      std::cout << "EndOfEvent fEdepScatterer = " << G4BestUnit(fEdepScatterer, "Energy") <<" at time " << G4BestUnit(fTimeScatterer + fBeginTime, "Time") << std::endl;
      std::cout << "EndOfEvent fEdepDetector = " << G4BestUnit(fEdepDetector, "Energy") << " at time " << G4BestUnit(fTimeDetector + fBeginTime, "Time") << std::endl;

      fFirstWrite = false;
    }
}
// Can use incremented copy number for when have many scatterers and detectors
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
