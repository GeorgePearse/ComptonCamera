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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "time.h" 
#include "iostream"

#include <iomanip>      // put_time
#include <ctime>        // time_t
#include <chrono>       // system_clock
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  count(0.),
  numberUseful(0.),
  numberUseless(0.),
  OneScatter(0.), // need to have in header
  MoreScatter(0.),
  OneScatterEscape(0.), // need to have in header
  MoreScatterEscape(0.),
  photonScattererCount(0.)

{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
  //int counter; //may not be needed
  //int numberUseful;
  //int numberUseless; 
  //G4bool ffirstWrite3 = true;

  // Set initial counter of photons in scatterer and absorber to zero
  photonScattererCount = 0;
  photonAbsorberCount = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

std::ofstream Efficiency;
// using std::chrono::system_clock;
// std::time_t tt = system_clock::to_time_t (system_clock::now());
// struct std::tm * ptm = std::localtime(&tt);
G4double Scent = (numberUseless / count)*100;
G4double Lcent = (numberUseful / count)*100;

  Efficiency.open(fOutput + "Efficiency.txt", std::ios_base::app);
  if (Efficiency.is_open()){
	//File Prints: Time/Date Count Useless Useless(%) Useful Useful(%)
	//put_time(ptm,"%c ")
  	Efficiency<<count<<" "<<numberUseless<<" "<<Scent<<" "<<numberUseful<<" "<<Lcent<<"\n";
  	}
  else Efficiency << "Unable to open file\n";


// George, the other version doesn't work without coincidence -can't analyse one detector
 std::ofstream totalComptons;
 totalComptons.open("totalComptons.txt", std::ios_base::app);
  if (totalComptons.is_open()){
  	totalComptons<<count<<" "<<OneScatter<<" "<<MoreScatter<<" "<<OneScatterEscape<<" "<<MoreScatterEscape<<"\n";
  	}
  else totalComptons << "Unable to open file\n";





  // Save file containing the number of photons in the scatterer and absorber - by Jack
  std::ofstream photonsInVolumeFile;
  photonsInVolumeFile.open(fOutput + "totalPhotonsInUsefulVolumes.txt");
  if(photonsInVolumeFile.is_open())
    {
      photonsInVolumeFile << photonScattererCount << " " << photonAbsorberCount << std::endl;
    }

	//G4ofstreamDestinationBase(Efficiency.txt,true)
    	//	{
        //	  G4cerrbuf.SetDestination(count);
    	//	}

	//G4UImanager::GetUIpointer();
        //G4iosInitialization();
        //G4CoutToFile myout("MyOut_out.log")
        //G4CoutToFile("Efficiency.txt", true);

  std::cout<<"numberUseful="<<numberUseful<< "\n"; //George 
  std::cout<<"numberUseless="<<numberUseless<< "\n"; //George 
  //std::cout<<"count="<<count<< "\n"; //George you were using count somewhere?
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(edep,"Dose") << " rms = " << G4BestUnit(rms,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

