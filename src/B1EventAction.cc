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
#include <G4SystemOfUnits.hh>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdepScatterer(0.),
  fEdepDetector(0.),
  fRunTime(0.)
{
fFirstWrite = true;
fPeakBroaden = true;
fFirstWritePosCount = true;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Originally by Jack, generalised with copy number by Douglas
void B1EventAction::AddEdepScatterer(G4double edep, int copyNo)
 
{
  {fEdepScatterer += edep;
  fScatCopyNo = std::to_string(copyNo);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Originally by Jack, generalised with copy number by Douglas
void B1EventAction::TimeScatterer(G4double timeScatterer, int copyNo)
{
  fTimeScatterer = timeScatterer;
  fScatCopyNo = std::to_string(copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Originally by Jack, generalised with copy number by Douglas
void B1EventAction::AddEdepDetector(G4double edep, int copyNo)
{
  fEdepDetector += edep;
  fAbsorbCopyNo = std::to_string(copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Originally by Jack, generalised with copy number by Douglas
void B1EventAction::TimeDetector(G4double timeDetector, int copyNo)
{
  fTimeDetector = timeDetector;
  fAbsorbCopyNo = std::to_string(copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// By Douglas
void B1EventAction::PeakBroad(double g, double c, bool scatter = true)
{
if(scatter == true)
    {
    std::cout << " dirac scat peak = " << fEdepScatterer/keV << std::endl;
    double Sigma = std::exp(c)*std::pow(fEdepScatterer*1000,(1-g))/2.35482;
    fEdepScatterer = G4RandGauss::shoot(fEdepScatterer*1000, Sigma)/1000;
    std::cout << " broad scat peak = " << fEdepScatterer << std::endl;
    }
else

    {
    std::cout << " dirac absorb peak = " << fEdepDetector/keV << std::endl;
    double Sigma = std::exp(c)*std::pow(fEdepDetector*1000,1-g)/2.35482;
    fEdepDetector = G4RandGauss::shoot(fEdepDetector*1000, Sigma)/1000;
    std::cout << " broad absorb peak = " << fEdepDetector << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdepScatterer = 0.;
  fEdepDetector = 0.;
  N = 0.;
  fBeginTime = fRunTime;
  posList.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B1EventAction::EndOfEventAction(const G4Event*)
{ // By Douglas
  if(fEdepScatterer != 0 && fEdepDetector != 0)
    { 
      if(fPeakBroaden == true)
      {
	B1EventAction::PeakBroad(0.4209, 0.1962, true);
	B1EventAction::PeakBroad(0.3974, 0.04931, false);
      }
      
      // Text file writer for Scatterer
      std::ofstream myfile;
      scatName = "scatter" + fScatCopyNo + "data.txt";
      absorbName = "absorb" + fAbsorbCopyNo + "data.txt";
      // Special condition for first write to create file
      if(fFirstWrite)
	{
		myfile.open(scatName);
	}
	else
	{
		myfile.open (scatName, std::ios::app);
	}
      	if (myfile.is_open())
      	{
          myfile << (fTimeScatterer + fBeginTime)/1000  << " "
	  << fEdepScatterer/keV << "\n";
	  myfile.close();
        }
	else std::cerr << "Unable to open scatter file" << std::endl;


       // Text file writer for Absorber
       std::ofstream myfile2;
       // Special condition for first write to create file
       if(fFirstWrite)
	{
		myfile2.open(absorbName);
	}
	else
	{
		myfile2.open (absorbName, std::ios::app);
	}
      	if (myfile2.is_open())
      	{
          myfile2 << (fTimeDetector + fBeginTime)/1000 << " "
	  << fEdepDetector/keV << "\n";
	  myfile2.close();
        }
	else std::cerr << "Unable to open absorb file" << std::endl;

      std::cout << "EndOfEvent fEdepScatterer = " << G4BestUnit(fEdepScatterer, "Energy") <<" at time " << G4BestUnit(fTimeScatterer + fBeginTime, "Time") << std::endl;
      std::cout << "EndOfEvent fEdepDetector = " << G4BestUnit(fEdepDetector, "Energy") << " at time " << G4BestUnit(fTimeDetector + fBeginTime, "Time") << std::endl;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  // Text file writer for Scatterer Position and Count - by Ben
  std::ofstream myfile3;
  // Special condition for first write to create file
  if(fFirstWritePosCount)
	{
          myfile3.open("Scat_PosCount.txt");
	}
  else
	{
	  myfile3.open("Scat_PosCount.txt", std::ios::app);
	}
  if (myfile3.is_open())
      	{
	myfile3 << "New Event" << "\n";
	for(unsigned int i=0; i<posList.size(); i++)
	{
	  myfile3 << posList[i] << "\n";
	}
	myfile3.close();
	}
   else std::cerr << "Unable to open Scat_PosCount file" << std::endl;
  fFirstWrite = false;
  fFirstWritePosCount = false;
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.x.....cc
