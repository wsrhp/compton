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
/// \file field/field01/include/F01FieldSetup.hh
/// \brief Definition of the F01FieldSetup class
//
//
// $Id: F01FieldSetup.hh 77115 2013-11-21 15:06:37Z gcosmo $
//
//
//  A class for control of the Magnetic Field of the detector.
//  The field is assumed to be uniform.
//
// Should this be a:
//    i) messenger
//   ii) user class that creates the field       ?
//  iii) simply a derived class of Uniform field ?  <== I have chosen this now.
//   iv) a field manager that creates/updates field    (Prefered?)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F01FieldSetup_h
#define F01FieldSetup_h 1

#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "globals.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class F01FieldMessenger;

class F01FieldSetup: public G4MagneticField
//class F01FieldSetup
{
public:
  F01FieldSetup(G4ThreeVector);  //  The value of the field
  F01FieldSetup();               //  A zero field

  virtual ~F01FieldSetup();

  void SetStepperType( G4int i ) { fStepperType = i; }

  void SetStepper();

  void SetMinStep(G4double sss) { fMinStep = sss; }

  void InitialiseAll();    //  Set parameters and call method below
  void CreateStepperAndChordFinder();

  void SetFieldValue(G4ThreeVector fieldVector);
  void SetFieldValue(G4double      fieldValue);
  G4ThreeVector GetConstantFieldValue();

  G4FieldManager* GetLocalFieldManager() { return fLocalFieldManager;}
  void GetFieldValue(const G4double Point[3], G4double* Bfield) const;

  ///  IsInBoundingBox() returns true if the point is within the
  ///  global bounding box - global coordinates.
  bool IsInBoundingBox(const G4double point[4]) const
  {
    if(point[2] < -75 || point[2] > 80) return false;
    if(point[0] < -20 || point[0] > 20) return false;
    if(point[1] < -105 || point[1] > 105) return false;
    return true;
  }


protected:

  // Find the global Field Manager

  G4FieldManager*         GetGlobalFieldManager();

  G4FieldManager*         fFieldManager;
  G4PropagatorInField*    fFieldPropagator;
  G4FieldManager*         fLocalFieldManager;
  G4ChordFinder*          fChordFinder;
  G4ChordFinder*          fLocalChordFinder;
  G4Mag_UsualEqRhs*       fEquation;
  G4Mag_UsualEqRhs*       fLocalEquation;
  G4MagneticField*        fMagneticField;
  G4MagneticField*        fLocalMagneticField;

  G4MagIntegratorStepper* fStepper;
  G4MagIntegratorStepper* fLocalStepper;
  G4int                   fStepperType;

  G4double                fMinStep;
 
  F01FieldMessenger*      fFieldMessenger;
private:
  static const int imesh=8311; //这里定义磁场矢量场坐标和数值，以及数据点数量
  float x[imesh];
  float y[imesh];
  float z[imesh];
  float Bx[imesh];
  float By[imesh];
  float Bz[imesh];

};

#endif
