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
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class B1DetectorMessenger;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual void SetScatXPos(G4double);
    virtual void SetScatYPos(G4double);
    virtual void SetScatZPos(G4double);
    virtual void SetScatPolarR(G4double);
    virtual void SetScatPolarPhi(G4double);
    virtual void SetScatPolarTheta(G4double);
    virtual void SetScatRotX(G4double);
    virtual void SetScatRotY(G4double);
    virtual void SetScatRotZ(G4double);
    virtual void SetScatRad(G4double);
    virtual void SetScatHeight(G4double);
    virtual void SetScatMat(G4String);
  
    virtual void SetScat2XPos(G4double);
    virtual void SetScat2YPos(G4double);
    virtual void SetScat2ZPos(G4double);
    virtual void SetScat2PolarR(G4double);
    virtual void SetScat2PolarPhi(G4double);
    virtual void SetScat2PolarTheta(G4double);
    virtual void SetScat2RotX(G4double);
    virtual void SetScat2RotY(G4double);
    virtual void SetScat2RotZ(G4double);
    virtual void SetScat2Rad(G4double);
    virtual void SetScat2Height(G4double);
    virtual void SetScat2Mat(G4String);
    virtual void SetScat2Bool(G4bool);
    virtual void SetScat2Switch(G4bool);
  
    virtual void SetDetXPos(G4double);
    virtual void SetDetYPos(G4double);
    virtual void SetDetZPos(G4double);
    virtual void SetDetPolarR(G4double);
    virtual void SetDetPolarPhi(G4double);
    virtual void SetDetPolarTheta(G4double);
    virtual void SetDetRotX(G4double);
    virtual void SetDetRotY(G4double);
    virtual void SetDetRotZ(G4double);
    virtual void SetDetRad(G4double);
    virtual void SetDetHeight(G4double);
    virtual void SetDetMat(G4String);
  
    virtual void SetDet2XPos(G4double);
    virtual void SetDet2YPos(G4double);
    virtual void SetDet2ZPos(G4double);
    virtual void SetDet2PolarR(G4double);
    virtual void SetDet2PolarPhi(G4double);
    virtual void SetDet2PolarTheta(G4double);
    virtual void SetDet2RotX(G4double);
    virtual void SetDet2RotY(G4double);
    virtual void SetDet2RotZ(G4double);
    virtual void SetDet2Rad(G4double);
    virtual void SetDet2Height(G4double);
    virtual void SetDet2Mat(G4String);
    virtual void SetDet2Bool(G4bool);
    virtual void SetDet2Switch(G4bool);
  
    virtual G4VPhysicalVolume* Construct();
    virtual G4VPhysicalVolume* ConstructVolumes();
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; } //for Dose

    virtual void UpdateGeometry();

  private:
  G4double fScatXPos;
  G4double fScatYPos;
  G4double fScatZPos;
  G4double fScatPolarR;
  G4double fScatPolarPhi;
  G4double fScatPolarTheta;
  G4double fScatRotX;
  G4double fScatRotY;
  G4double fScatRotZ;
  G4double fScatRad;
  G4double fScatHeight;
  
  G4double fScat2XPos;
  G4double fScat2YPos;
  G4double fScat2ZPos;
  G4double fScat2PolarR;
  G4double fScat2PolarPhi;
  G4double fScat2PolarTheta;
  G4double fScat2RotX;
  G4double fScat2RotY;
  G4double fScat2RotZ;
  G4double fScat2Rad;
  G4double fScat2Height;
  G4bool   fScat2Bool;
  G4String fScat2Type;
  G4int    fScat2CopyNo;
  
  G4double fDetXPos;
  G4double fDetYPos;
  G4double fDetZPos;
  G4double fDetPolarR;
  G4double fDetPolarPhi;
  G4double fDetPolarTheta;
  G4double fDetRotX;
  G4double fDetRotY;
  G4double fDetRotZ;
  G4double fDetRad;
  G4double fDetHeight;
  
  G4double fDet2XPos;
  G4double fDet2YPos;
  G4double fDet2ZPos;
  G4double fDet2PolarR;
  G4double fDet2PolarPhi;
  G4double fDet2PolarTheta;
  G4double fDet2RotX;
  G4double fDet2RotY;
  G4double fDet2RotZ;
  G4double fDet2Rad;
  G4double fDet2Height;
  G4bool   fDet2Bool;
  G4String fDet2Type;
  G4int    fDet2CopyNo;

  G4Material* shape1_mat;
  G4Material* shape2_mat;
  G4Material* shape3_mat;
  G4Material* shape4_mat;
  
  B1DetectorMessenger* fDetectorMessenger;
  G4bool fAluminium;

  protected:
    G4LogicalVolume*  fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

