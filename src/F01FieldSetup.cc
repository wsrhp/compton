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
/// \file field/field01/src/F01FieldSetup.cc
/// \brief Implementation of the F01FieldSetup class
//
//
// $Id: F01FieldSetup.cc 77115 2013-11-21 15:06:37Z gcosmo $
//
//   User Field setup class implementation.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F01FieldSetup.hh"
#include "F01FieldMessenger.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include<iomanip>
#include <cstring>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//  Constructors:

F01FieldSetup::F01FieldSetup()
 : fFieldManager(0),
   fLocalFieldManager(0),
   fChordFinder(0),
   fLocalChordFinder(0),
   fEquation(0),
   fLocalEquation(0),
   fMagneticField(0),
   fLocalMagneticField(0),
   fStepper(0),
   fLocalStepper(0),
   fStepperType(0),
   fMinStep(0.),
   fFieldMessenger(0)
{
    ifstream bmap;
    string line;
    const char *charline;
    int idx=0;

    bmap.open("Bfield.dat");
    while(bmap) {
        getline(bmap,line);
        charline=line.c_str();
        sscanf(charline,"%f\t%f\t%f\t%f\t%f\t%f\n",
               &z[idx],&y[idx],&x[idx],&Bz[idx],&By[idx],&Bx[idx]);
        //更改数据坐标轴为G4中的坐标轴，因为涉及到后面的数据处理，所以不宜更改G4坐标轴
        //static ofstream ofile("MagneticField.txt",ios_base::out);
        //ofile<< "mag field check i = " << idx << " x = " << x[idx] << " Bx = " << Bx[idx] <<"\r\n"<<endl;
        //G4cout << "mag field check i = " << idx << " x = " << x[idx] << " Bx = " << Bx[idx] <<G4endl;
        idx++;
    }
    InitialiseAll();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::InitialiseAll()
{
    fMinStep     = 0.1*mm; // minimal step of 1 mm is default

    fStepperType = 4;      // ClassicalRK4 is default stepper

    //fLocalMagneticField = new G4UniformMagField(G4ThreeVector(0.3*tesla,0.0,0.0));

    //fEquation = new G4Mag_UsualEqRhs(fMagneticField);
    //  Construct equ. of motion of particles including spin through e.m. fields
    fLocalEquation = new G4Mag_UsualEqRhs(this);

    //  Get transportation, field, and propagator managers
    G4TransportationManager* transportManager = G4TransportationManager::GetTransportationManager();

    fLocalFieldManager = GetGlobalFieldManager();
    //fLocalFieldManager = new G4FieldManager();

    fFieldPropagator = transportManager->GetPropagatorInField();

    fFieldMessenger = new F01FieldMessenger(this);

    CreateStepperAndChordFinder();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F01FieldSetup::~F01FieldSetup()
{
  delete fMagneticField;
  delete fChordFinder;
  delete fStepper;
  delete fFieldMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::CreateStepperAndChordFinder()
{
  // Update field


  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl;

  //fFieldManager->SetDetectorField(fMagneticField );
  fLocalFieldManager->SetDetectorField(this);

  //if (fChordFinder) delete fChordFinder;
  if (fLocalChordFinder) delete fLocalChordFinder;

  //fChordFinder = new G4ChordFinder( fMagneticField, fMinStep,fStepper );
  fLocalChordFinder = new G4ChordFinder(fLocalMagneticField,
                                        fMinStep,fLocalStepper);

  //fFieldManager->SetChordFinder( fChordFinder );
   fLocalFieldManager->SetChordFinder(fLocalChordFinder);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::SetStepper()
{
// Set stepper according to the stepper type

  if (fStepper) delete fStepper;

  switch ( fStepperType )
  {
  case 0:
    fStepper = new G4ExplicitEuler( fEquation );
    fLocalStepper = new G4ExplicitEuler( fLocalEquation );
    G4cout<<"G4ExplicitEuler is called"<<G4endl;
    break;
  case 1:
    fStepper = new G4ImplicitEuler( fEquation );
    fLocalStepper = new G4ImplicitEuler( fLocalEquation );
    G4cout<<"G4ImplicitEuler is called"<<G4endl;
    break;
  case 2:
    fStepper = new G4SimpleRunge( fEquation );
    fLocalStepper = new G4SimpleRunge( fLocalEquation );
    G4cout<<"G4SimpleRunge is called"<<G4endl;
    break;
  case 3:
    fStepper = new G4SimpleHeum( fEquation );
    fLocalStepper = new G4SimpleHeum( fLocalEquation );
    G4cout<<"G4SimpleHeum is called"<<G4endl;
    break;
  case 4:
    fStepper = new G4ClassicalRK4( fEquation );
    fLocalStepper = new G4ClassicalRK4( fLocalEquation );
    G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
    break;
  case 5:
    fStepper = new G4HelixExplicitEuler( fEquation );
    fLocalStepper = new G4HelixExplicitEuler( fLocalEquation );
    G4cout<<"G4HelixExplicitEuler is called"<<G4endl;
    break;
  case 6:
    fStepper = new G4HelixImplicitEuler( fEquation );
    fLocalStepper = new G4HelixImplicitEuler( fLocalEquation );
    G4cout<<"G4HelixImplicitEuler is called"<<G4endl;
    break;
  case 7:
    fStepper = new G4HelixSimpleRunge( fEquation );
    fLocalStepper = new G4HelixSimpleRunge( fLocalEquation );
    G4cout<<"G4HelixSimpleRunge is called"<<G4endl;
    break;
  case 8:
    fStepper = new G4CashKarpRKF45( fEquation );
    fLocalStepper = new G4CashKarpRKF45( fLocalEquation );
    G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
    break;
  case 9:
    fStepper = new G4RKG3_Stepper( fEquation );
    fLocalStepper = new G4RKG3_Stepper( fLocalEquation );
    G4cout<<"G4RKG3_Stepper is called"<<G4endl;
    break;
  default: fStepper = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::SetFieldValue(G4double fieldStrength)
{
  // Set the value of the Global Field to fieldValue along Z

#ifdef G4VERBOSE
  G4cout << "Setting Field strength to "
         << fieldStrength / gauss  << " Gauss."; // << G4endl;
#endif

  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
  this->SetFieldValue( fieldSetVec );

#ifdef G4VERBOSE
  G4double fieldValue[6],  position[4];
  position[0] = position[1] = position[2] = position[3] = 0.0;
  if ( fieldStrength != 0.0 ) {
    fMagneticField->GetFieldValue( position, fieldValue);
    G4ThreeVector fieldVec(fieldValue[0], fieldValue[1], fieldValue[2]);
    // G4cout << " fMagneticField is now " << fMagneticField
    G4cout << " Magnetic field vector is "
           << fieldVec / gauss << " G " << G4endl;
  } else {
    if ( fMagneticField == 0 )
      G4cout << " Magnetic field pointer is null." << G4endl;
    else
      G4Exception("F01FieldSetup::SetFieldValue(double)",
                  "IncorrectForZeroField",
                  FatalException,
                  "fMagneticField ptr should be set to 0 for no field.");
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  // Set the value of the Global Field
  if (fMagneticField) delete fMagneticField;
 
  if (fieldVector != G4ThreeVector(0.,0.,0.))
  {
    fMagneticField = new G4UniformMagField(fieldVector);
  }
  else
  {
    // If the new field's value is Zero, signal it as below
    // so that it is not used for propagation.
    fMagneticField = 0;
  }

  // Set this as the field of the global Field Manager
  GetGlobalFieldManager()->SetDetectorField(fMagneticField);

  // Now notify equation of new field
  fEquation->SetFieldObj( fMagneticField );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F01FieldSetup::GetFieldValue(const G4double Point[3],G4double* Bfield) const
{
    Bfield[0] = Bfield[1] = Bfield[2] = Bfield[3] = Bfield[4] = Bfield[5] = 0.0;
    G4float dmesh= 7.0 ;

    // protect against Geant4 bug that calls us with point[] NaN.
    if(Point[0] != Point[0]) return;
     if(IsInBoundingBox(Point)){
        //static ofstream ofile("MagneticField.txt",ios_base::out);
        //ofile<< "Point[2] = " <<Point[2]<< "Point[0] = " <<Point[0]<< "Point[1] = " <<Point[1]<<"\r\n"<<endl;
   /* Bfield[0] = 0.3*tesla;
    Bfield[1] = 0.;
    Bfield[2] = 0.;

   if(fabs(Point[0]-x[60])<dmesh&&fabs(Point[1]-y[60])<dmesh&&fabs(Point[2]-z[60])<dmesh){
        static ofstream ofile("MagneticField.txt",ios_base::out);
        ofile<< "Point[2] = " <<Point[2]<< "Point[0] = " <<Point[0]<< "Point[1] = " <<Point[1]<<"\r\n"<<endl;
    }*/

       for(G4int i=0; i<8311;i++){
            if(fabs(Point[0]-x[i])<dmesh&&fabs(Point[1]-y[i])<dmesh&&fabs(Point[2]-z[i])<dmesh){
                //static ofstream ofile("MagneticField.txt",ios_base::out);
                //ofile<<"i="<<i<< "  Point[2] = " <<Point[2]<< "Point[0] = " <<Point[0]<< "Point[1] = " <<Point[1]<<"\r\n"<<endl;
                //ofile<<i<<"\t"<<endl;
                Bfield[0]=-Bx[i]*tesla;
                Bfield[1]=By[i]*tesla;
                Bfield[2]=Bz[i]*tesla;

                break;
            }
        }
     }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4FieldManager* F01FieldSetup::GetGlobalFieldManager()
{
  //  Utility method

  return G4TransportationManager::GetTransportationManager()
           ->GetFieldManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
