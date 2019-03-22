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
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class B1RunAction;
class G4GenericMessenger;

/// Event action class
///

class B1EventAction : public G4UserEventAction
{
  public:
    B1EventAction(B1RunAction* runAction);
    virtual ~B1EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    virtual void AddEdepScatterer(G4double edep, int copyNo);
    virtual void AddEdepDetector(G4double edep, int copyNo);
    virtual void AddEdepBody(G4double edep);
    virtual void TimeScatterer(G4double timeScatterer, int copyNo);
    virtual void TimeDetector(G4double timeDetector, int copyNo);
    virtual void PeakBroad(double g, double c, bool scatter);
    virtual void SetOutput(std::string folderName);
    virtual void ZeroScatterInfo(G4double edep, G4String procName, G4ThreeVector pos);
    virtual void Vector(G4ThreeVector Pos, int copyNo);
    virtual void Vector2(G4ThreeVector Pos, int copyNo);
    virtual void DeltaMomentum(G4ThreeVector deltaMomentum);	
    virtual void DeltaComptonEnergy(G4double deltaComptonEnergy);
    virtual void EnergyExit(G4double energyExit);

    void TotalTime(G4double deltaTime){fRunTime += deltaTime;};
    void Proc2(G4String procName){procList2.push_back(procName);};
    void Count(){N += 1;};
    void totalComptons(){M += 1;};
    void exit(){exitBool = true;}; 
    void PhotonScatterer(){photonScattererCount+=1;};
    void PhotonAbsorber(){photonAbsorberCount+=1;};
    //void PhotonDeltaMomentum(deltaMomentum){fPhotonDeltaMomentum=deltaMomentum};

    
  
  private:
    B1RunAction* fRunAction;
    G4double     fEdepScatterer;
    G4double     fEdepDetector;
    G4double	 fEdepBody;
    G4double fTimeScatterer;
    G4double fTimeDetector;
    G4double fRunTime;
    G4double fBeginTime;
    G4bool fFirstWrite;
    G4bool fPeakBroaden;
    G4bool fFirstWritePosCount;
    G4bool fFirstWritePosCount2;
    G4bool coincidence;
    G4bool fFirstWrite2;
    G4bool fFirstWriteNotCompt;
    G4bool fFirstWriteTotal;
    G4bool fFirstWriteTotal2;
    G4bool fPhotonMomentum;
    G4bool fMomentumWrite;


    bool exitBool;
    G4double fenergyExit;

    std::string fScatCopyNo;
    std::string fAbsorbCopyNo;
    std::string absorbName;
    std::string scatName;
    std::string totalscatName;
    std::string totalabsorbName;
    std::string momentumchange;
    int N;
    int M;

    G4ThreeVector fdeltaMomentum;
    G4double fdeltaComptonEnergy;

    int counter;
    int photonScattererCount;
    int photonAbsorberCount;
    std::vector<G4ThreeVector> posList;
    std::vector<G4ThreeVector> posList2;
    std::vector<G4ThreeVector> posListNotCompt;
    std::vector<G4String> procListNotCompt;
    std::vector<G4String> procList2;
    std::vector<G4double> edepListNotCompt;
    G4GenericMessenger* fMessenger;
    std::string fOutput;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
