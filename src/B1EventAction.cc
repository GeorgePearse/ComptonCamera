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

#include "G4GenericMessenger.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "G4Track.hh"
#include "G4Step.hh"
#include "B1SteppingAction.hh"

#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdepScatterer(0.),
  fEdepDetector(0.),
  fEdepBody(0.),
  fRunTime(0.),
  fMessenger(0)
{
fFirstWrite = true;
fPeakBroaden = true;   // should be set to true by default.
fFirstWritePosCount = true;
fFirstWritePosCount2 = true;
coincidence = true; //should be set to true unless a material test is being carried out 
fFirstWrite2 = true;
fFirstWriteTotal = true;
fFirstWriteTotal2 = true;
fFirstWriteJackEnergies = true;
fPhotonMomentum = false;
fMomentumWrite = false;
fOutput = "";
counter = 0; 
// Event action generic messenger, used to set output folder - by Jack
 fMessenger = new G4GenericMessenger(this, "/B1/eventAction/", "EventAction control");
 auto& outputCommand = fMessenger->DeclareMethod("setOutput", &B1EventAction::SetOutput, "sets output folder");

 // Makes sure all decimal places are written to files - by Jack
 std::cout.precision(15);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{
  delete fMessenger;
}

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

// Originally by Jack, generalised with copy number by Douglas
void B1EventAction::AddEdepDetector(G4double edep, int copyNo)
{
  fEdepDetector += edep;
  fAbsorbCopyNo = std::to_string(copyNo);
}
//Written by George
void B1EventAction::AddEdepBody(G4double edep)
{fEdepBody += edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Originally by Jack, generalised with copy number by Douglas
void B1EventAction::TimeDetector(G4double timeDetector, int copyNo)
{
  fTimeDetector = timeDetector;
  fAbsorbCopyNo = std::to_string(copyNo);
}


//WrittenByGeorge
void B1EventAction::DeltaMomentum(G4ThreeVector preMomentum, G4ThreeVector postMomentum)
{
fpreMomentum = preMomentum;
fpostMomentum = postMomentum;
}

void B1EventAction::DeltaComptonEnergy(G4double deltaComptonEnergy)
{
fdeltaComptonEnergy = deltaComptonEnergy;
}

void B1EventAction::EnergyExit(G4double energyExit)
{
fenergyExit = energyExit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// By Douglas
void B1EventAction::PeakBroad(double g, double c, bool scatter = true)
{
if(scatter == true)
    {
    //std::cout << " dirac scat peak = " << fEdepScatterer/keV << std::endl;
    double ScatSigma = std::exp(c)*std::pow(fEdepScatterer*1000,(1-g))/2.35482;
    //double ScatFWHM = std::exp(c)*std::pow(fEdepScatterer*1000,(1-g));
    //std::cout << "scatterer FWHM = " << ScatFWHM << std::endl;
    fEdepScatterer = (G4RandGauss::shoot(fEdepScatterer*1000, ScatSigma))/1000;
    //std::cout << " broad scat peak = " << fEdepScatterer << std::endl;
    }
else

    {
    //std::cout << " dirac absorb peak = " << fEdepDetector/keV << std::endl;
    double DetSigma = std::exp(c)*std::pow(fEdepDetector*1000,1-g)/2.35482;
    //double DetFWHM = std::exp(c)*std::pow(fEdepDetector*1000,1-g);
    //std::cout << "absorber FWHM = " << DetFWHM << std::endl;
    fEdepDetector = (G4RandGauss::shoot(fEdepDetector*1000, DetSigma))/1000;
    //std::cout << " broad absorb peak = " << fEdepDetector << std::endl;
    
    //George basic peak broadening 
    //Has no effect provided Sigma is written in "fEdepDetector =" and not Sigma2
    double resBGO = 0.1335; //Second one
    double resLSO = 0.105; //First one 
    double resCdWO4 = 0.066; //Third one
    double restheBeast = 0.029; //Fourth one
    double Sigma2 = resCdWO4*(std::pow(fEdepDetector*1000*662, 0.5))/2.35;
    //endGeorge REMEMBER TO CHANGE BACK TO SIGMA BELOW !!!
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::SetOutput(std::string folderName)
{
  fOutput = folderName;
  system(("mkdir " + fOutput).c_str());
  // Send output folder to run action for the file writers there
  fRunAction->OutputFolder(folderName);
}

void B1EventAction::ZeroScatterInfo(G4double edep, G4String procName, G4ThreeVector pos)
{
  posListNotCompt.push_back(pos);
  procListNotCompt.push_back(procName);
  edepListNotCompt.push_back(edep);
}

// Scatterer and absorber positions - by Ben, generalised using copyNo by Jack
void B1EventAction::Vector(G4ThreeVector Pos, int copyNo)
{
  posList.push_back(Pos);
  fScatCopyNo = std::to_string(copyNo);
}

void B1EventAction::Vector2(G4ThreeVector Pos, int copyNo)
{
  posList2.push_back(Pos);
  fAbsorbCopyNo = std::to_string(copyNo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// By Douglas
void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdepScatterer = 0.;
  fEdepDetector = 0.;
  fEdepBody = 0.;
  fdeltaComptonEnergy = 0.; 
  exitBool = false; // Should be set to false and should follow same rules as N and M 
  N = 0.;
  M = 0.;
   // ComptonScatters for one detector George

  fBeginTime = fRunTime;
  posList.clear();
  posList2.clear();
  posListNotCompt.clear();
  procListNotCompt.clear();
  edepListNotCompt.clear();
  scatterEnergies.clear();
  if (counter%50000 == 0)
  {
  std::cout << " total event counter = " << counter << std::endl;
  }
  counter += 1;
  fRunAction->Count(); //scared this may double count something??
  // Set photon counts in absorber/scatterer to zero at start of events
  photonScattererCount = 0;
  photonAbsorberCount = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B1EventAction::EndOfEventAction(const G4Event*)
{ 

// condition to print all scatter events into a file coincident and non-coincident.
// By Douglas

if (coincidence == false)
   {
     if(fEdepScatterer != 0)
      {
        if(fPeakBroaden == true)
       {
	B1EventAction::PeakBroad(0.5254, 0.7222, true);
       }
      // Text file writer for total Scatterer
      std::ofstream myfiletotal;
      totalscatName = fOutput + "totalscatter" + fScatCopyNo + "data.txt";
      // Special condition for first write to create file
      G4bool switch1 = true;
      if(fFirstWriteTotal && switch1 == true)
	{
		myfiletotal.open(totalscatName);
	}
	else
	{
		myfiletotal.open (totalscatName, std::ios::app);
	}
      	if (myfiletotal.is_open() && switch1 == true)
      	{
          myfiletotal << (fTimeScatterer + fBeginTime)/1000  << " " << fEdepScatterer/keV << "\n";
	  myfiletotal.close();
        }
	else std::cerr << "Unable to open scatter file" << std::endl;
        fFirstWriteTotal = false;
        }

      if(fEdepDetector != 0)
       {
        if(fPeakBroaden == true)
        {
          B1EventAction::PeakBroad(0.3871, -0.5296, false);
         }
        // Text file writer for Absorber
       std::ofstream myfiletotal2;
       totalabsorbName = fOutput + "totalabsorb" + fAbsorbCopyNo + "data.txt";
       // Special condition for first write to create file
       G4bool switch2 = true;
       if(fFirstWriteTotal2 && switch2 == true)
	{
		myfiletotal2.open(totalabsorbName);
	}
	else
	{
		myfiletotal2.open (totalabsorbName, std::ios::app);
	}
      	if (myfiletotal2.is_open() && switch2 == false)
      	{
          myfiletotal2 << (fTimeDetector + fBeginTime)/1000 << " "
	  << fEdepDetector/keV << "\n";
	  myfiletotal2.close();
        }
	else std::cerr << "Unable to open absorb file" << std::endl;
        fFirstWriteTotal2 = false;
        }
        
       if(fPhotonMomentum == true && M==1)
         {
           std::ofstream myfilemomentum;
           momentumchange = fOutput + "momentumchange" + fAbsorbCopyNo + "data.txt";
           if(fMomentumWrite)
             {
              myfilemomentum.open(momentumchange);
              }
           else
              {
               myfilemomentum.open(momentumchange, std::ios::app);
               }
            if(myfilemomentum.is_open())
              {
		if(fdeltaComptonEnergy2/keV > 300 and fdeltaComptonEnergy2/keV < 301)
		{ 
               myfilemomentum << fpreMomentum << " " << fpostMomentum << " " << fdeltaComptonEnergy/keV << "\n";    
		}
               myfilemomentum.close();
               }
            else std::cerr << "Unable to open momentum file" << std::endl;
            fMomentumWrite = false;
          }
   }
//if(M==1){std::cout<<fdeltaMomentum<<" "<<fdeltaComptonEnergy<<"\n";}
//By George - Analysis of effect on height and radius of detector
if(coincidence==false && exitBool == true && fTimeScatterer<fTimeDetector){
if(M==1){fRunAction->Count1ScatterEscape();}
if(M>1){fRunAction->CountMoreScatterEscape();}
}
//if(coincidence==false && fTimeScatterer<fTimeDetector){
//if(M==1){fRunAction->Count1Scatter();};
//if(M>1){fRunAction->CountMoreScatter();};


// By Douglas
  if(fEdepScatterer != 0 && fEdepDetector != 0 && fTimeScatterer < fTimeDetector)
    {
      if(N==1 && fTimeScatterer<fTimeDetector) // By George 
	{
	  fRunAction->CountUseful();
	}
      else
	{
	  fRunAction->CountUseless();
	} // End of by George
      if(fPeakBroaden == true)
      {
	//Sodium Iodide (14*43mm) gradient and intercept of Ln(E) Ln(R) plot
	B1EventAction::PeakBroad(0.5254, 0.7222, true); 
	//Lanthanum Bromide gradient and intercept of Ln(E) Ln(R) plot
	B1EventAction::PeakBroad(0.3871, -0.5296, false); 
      }
      
      // Text file writer for Scatterer
      std::ofstream myfile;
      scatName = fOutput + "scatter" + fScatCopyNo + "data.txt";
      absorbName = fOutput + "absorb" + fAbsorbCopyNo + "data.txt";
      // Special condition for first write to create file
      G4bool switch3 = true;
      if(fFirstWrite && switch3 == true)
	{
		myfile.open(scatName);
	}
	else
	{
		myfile.open (scatName, std::ios::app);
	}
      	if (myfile.is_open() && switch3 == true)
      	{
          myfile << (fTimeScatterer + fBeginTime)/1000  << " "
	  << fEdepScatterer/keV << "\n";
	  myfile.close();
        }
	else std::cerr << "Unable to open scatter file" << std::endl;
       


       // Text file writer for Absorber
       std::ofstream myfile2;
       // Special condition for first write to create file
       G4bool switch4 = true;
       if(fFirstWrite && switch4 == true)
	{
		myfile2.open(absorbName);
	}
	else
	{
		myfile2.open (absorbName, std::ios::app);
	}
      	if (myfile2.is_open() && switch4 == true)
      	{
          myfile2 << (fTimeDetector + fBeginTime)/1000 << " "
	  << fEdepDetector/keV << "\n";
	  myfile2.close();
        }
	else std::cerr << "Unable to open absorb file" << std::endl;
        
       fFirstWrite = false;

      //std::cout << "EndOfEvent fEdepScatterer = " << G4BestUnit(fEdepScatterer, "Energy") <<" at time " << G4BestUnit(fTimeScatterer + fBeginTime, "Time") << std::endl;
      //std::cout << "EndOfEvent fEdepDetector = " << G4BestUnit(fEdepDetector, "Energy") << " at time " << G4BestUnit(fTimeDetector + fBeginTime, "Time") << std::endl;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  // Text file writer for Scatterer Position and Count - by Ben, generalised by Jack w/ copyNo
  // Special condition for first write to create file
	std::ofstream myfile3;
  	if(fFirstWritePosCount)
		{ 
        	  myfile3.open(fOutput + "Scat_PosCount" + fScatCopyNo + ".txt");
		}
  	else
		{
		  myfile3.open(fOutput + "Scat_PosCount" + fScatCopyNo + ".txt", std::ios::app);
		}
  	if (myfile3.is_open())
      		{
		myfile3 << "N" << "\n";
		for(unsigned int i=0; i<posList.size(); i++)
		{
		  myfile3 << posList[i] << "\n";
		}
		myfile3.close();
		}
  	else std::cerr << "Unable to open Scat_PosCount file" << std::endl;
  	fFirstWritePosCount = false;
  G4bool posCountSwitch = true;
  if (posList2.size() > 0 && posCountSwitch == true)
	{
	std::ofstream myfile4;
 	if(fFirstWritePosCount2)
		{
        	  myfile4.open(fOutput + "Abs_PosCount" + fAbsorbCopyNo + ".txt");
		}
  	else
		{
		  myfile4.open(fOutput + "Abs_PosCount" + fAbsorbCopyNo + ".txt", std::ios::app);
		}
  	if (myfile4.is_open())
      		{
		myfile4 << "N" << "\n";
		for(unsigned int j=0; j<posList2.size(); j++)
		{
		  //myfile4 << posList2[j] << " " << procList2[j] << "\n";
		  myfile4 << posList2[j] << " " << "\n";
		}
		myfile4.close();
		}
  	else std::cerr << "Unable to open Abs_PosCount file" << std::endl;
  	fFirstWritePosCount2 = false;
	}
// DEBUG STUFF BECAUSE DOUGLAS CODE RUNNING SLOW
  G4bool debug = false;
  if (procListNotCompt.size() > 0 && debug==true)
    {
      std::ofstream myfileJack;
      if(fFirstWriteNotCompt)
	{
	  myfileJack.open(fOutput + "scatPosProcNameNoCompt.txt");
	}
      else
	{
	  myfileJack.open(fOutput + "scatPosProcNameNoCompt.txt", std::ios::app); 
	}
      if (myfileJack.is_open())
	{
	  myfileJack << "New Event\n";
	  for(unsigned int i=0; i<posListNotCompt.size(); i++)
	    {
	      myfileJack << edepListNotCompt[i] << posListNotCompt[i] << " " << procListNotCompt[i] << std::endl;
	    }
	  myfileJack.close();
      }
    fFirstWriteNotCompt = false;
    }
  // File writer for double Compton scatter investigation
  if(N==2)
    {
      std::ofstream myfileJackEnergies;
      if(fFirstWriteJackEnergies)
	{
	  myfileJackEnergies.open(fOutput + "doubleScatterEnergies.txt");
	}
      else
	{
	  myfileJackEnergies.open(fOutput + "doubleScatterEnergies.txt", std::ios::app);
	}
      if(myfileJackEnergies.is_open())
	{
	  if(scatterEnergies.size() == 2)
	    {
	      myfileJackEnergies << scatterEnergies[0] << " " << scatterEnergies[1] << " " << fEdepScatterer << " " << fEdepDetector << std::endl;
	    }
	  else
	    {
	      // Writes this message if the method I've come up with doesn't work so that I know for certain
	      myfileJackEnergies << "BROKEN METHOD :(" << std::endl;
	    }
	   myfileJackEnergies.close();
	}
      else std::cerr << "Unable to open doubleScatterEnergies file" << std::endl;
      fFirstWriteJackEnergies = false;
    }
  
//Energy deposited in body by George
//if(fEdepBody!=0){
//G4bool patientBody = false;
//if(patientBody == true){
 //std::ofstream myfile5;
//if(fFirstWrite2)
//	{
//		myfile5.open("energyBody.txt");
//	}
//	else
//	{
//		myfile5.open ("energyBody.txt", std::ios::app);
//	}
  //    	if (myfile5.is_open())
    //  	{
     //     myfile5 << fEdepBody/keV << "\n";
//	  myfile5.close();
 //       }
//	else std::cerr << "Unable to open energyBody file" << std::endl;
//	fFirstWrite2 = false;};//}

  }
  // Summing photons that entered scatterer/absorber in run action - by Jack
  if(photonScattererCount>0)
    {
      fRunAction->PhotonScattererCount();
      if(photonAbsorberCount>0)
	{
	  fRunAction->PhotonAbsorberCount();
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.x.....cc
