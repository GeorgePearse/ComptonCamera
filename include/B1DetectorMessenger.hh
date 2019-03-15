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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef B1DetectorMessenger_h
#define B1DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class B1DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B1DetectorMessenger: public G4UImessenger
{
  public:
    B1DetectorMessenger(B1DetectorConstruction* );
   ~B1DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    B1DetectorConstruction*    fDetector;
    
    G4UIdirectory*             fB1Dir;
    G4UIdirectory*             fScatDir;
    G4UIdirectory*             fScat2Dir;
    G4UIdirectory*             fDetDir;
    G4UIdirectory*             fDet2Dir;

    G4UIcmdWithADoubleAndUnit* fScatXPos;
    G4UIcmdWithADoubleAndUnit* fScatYPos;
    G4UIcmdWithADoubleAndUnit* fScatZPos;
    G4UIcmdWithADoubleAndUnit* fScatPolarR;
    G4UIcmdWithADoubleAndUnit* fScatPolarPhi;
    G4UIcmdWithADoubleAndUnit* fScatPolarTheta;
    G4UIcmdWithADoubleAndUnit* fScatRotX;
    G4UIcmdWithADoubleAndUnit* fScatRotY;
    G4UIcmdWithADoubleAndUnit* fScatRotZ;
    G4UIcmdWithADoubleAndUnit* fScatRad;
    G4UIcmdWithADoubleAndUnit* fScatHeight;

    G4UIcmdWithADoubleAndUnit* fScat2XPos;
    G4UIcmdWithADoubleAndUnit* fScat2YPos;
    G4UIcmdWithADoubleAndUnit* fScat2ZPos;
    G4UIcmdWithADoubleAndUnit* fScat2PolarR;
    G4UIcmdWithADoubleAndUnit* fScat2PolarPhi;
    G4UIcmdWithADoubleAndUnit* fScat2PolarTheta;
    G4UIcmdWithADoubleAndUnit* fScat2RotX;
    G4UIcmdWithADoubleAndUnit* fScat2RotY;
    G4UIcmdWithADoubleAndUnit* fScat2RotZ;
    G4UIcmdWithADoubleAndUnit* fScat2Rad;
    G4UIcmdWithADoubleAndUnit* fScat2Height;
    G4UIcmdWithABool*          fScat2Bool;
    G4UIcmdWithABool*          fScat2Switch;
  
    G4UIcmdWithADoubleAndUnit* fDetXPos;
    G4UIcmdWithADoubleAndUnit* fDetYPos;
    G4UIcmdWithADoubleAndUnit* fDetZPos;
    G4UIcmdWithADoubleAndUnit* fDetPolarR;
    G4UIcmdWithADoubleAndUnit* fDetPolarPhi;
    G4UIcmdWithADoubleAndUnit* fDetPolarTheta;
    G4UIcmdWithADoubleAndUnit* fDetRotX;
    G4UIcmdWithADoubleAndUnit* fDetRotY;
    G4UIcmdWithADoubleAndUnit* fDetRotZ;
    G4UIcmdWithADoubleAndUnit* fDetRad;
    G4UIcmdWithADoubleAndUnit* fDetHeight;

    G4UIcmdWithADoubleAndUnit* fDet2XPos;
    G4UIcmdWithADoubleAndUnit* fDet2YPos;
    G4UIcmdWithADoubleAndUnit* fDet2ZPos;
    G4UIcmdWithADoubleAndUnit* fDet2PolarR;
    G4UIcmdWithADoubleAndUnit* fDet2PolarPhi;
    G4UIcmdWithADoubleAndUnit* fDet2PolarTheta;
    G4UIcmdWithADoubleAndUnit* fDet2RotX;
    G4UIcmdWithADoubleAndUnit* fDet2RotY;
    G4UIcmdWithADoubleAndUnit* fDet2RotZ;
    G4UIcmdWithADoubleAndUnit* fDet2Rad;
    G4UIcmdWithADoubleAndUnit* fDet2Height;
    G4UIcmdWithABool*          fDet2Bool;
    G4UIcmdWithABool*          fDet2Switch;
  
    G4UIcmdWithoutParameter*   fUpdateCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

