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
/// \file electromagnetic/TestEm14/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:G4VUserPrimaryGeneratorAction(),gps(0),fDetector(det)
{
  InitializeGPS();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete gps;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  gps->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::InitializeGPS()
{
   gps = new G4GeneralParticleSource();
    // setup details easier via UI commands see gps.mac
    //以下所有的GPS信息均可以通过gps.mac文件改变的。只是初始化这样而已。
    // particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pion = particleTable->FindParticle("e-");
  gps->GetCurrentSource()->SetParticleDefinition(pion);

  // set energy distribution
  G4SPSEneDistribution *eneDist = gps->GetCurrentSource()->GetEneDist() ;
  eneDist->SetEnergyDisType("Mono"); // or gauss
  eneDist->SetMonoEnergy(300.0*MeV);  //默认10MeV

  // set position distribution
  /*
  G4double halfSize = 0.5*(fDetector->GetSize());
  G4double z0 =  halfSize;
  G4double x0 = 0.0*mm;
  G4double y0 = 0.0*mm;
  G4SPSPosDistribution *posDist = gps->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Beam");  // or Point,Plane,Volume,Beam
  posDist->SetCentreCoords(G4ThreeVector(x0,y0,z0));
  posDist->SetBeamSigmaInX(0.1*mm);
  posDist->SetBeamSigmaInY(0.1*mm);
*/

  //G4double z0 = -180.0*mm;
  G4double z0 = -100.0*mm;
  G4double x0 = 0.0*um;
  G4double y0 = 0.0*um;
  G4SPSPosDistribution *posDist = gps->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Point");  // or Point,Plane,Volume,Beam
  posDist->SetCentreCoords(G4ThreeVector(x0,y0,z0));
/*  posDist->SetPosDisShape("Circle");
  posDist->SetRadius(1.0*um);
  posDist->SetBeamSigmaInR(0.0*mm);*/


  // set angular distribution
  G4SPSAngDistribution *angDist = gps->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection( G4ThreeVector(0., 0., 1.) );
  //angDist->SetAngDistType("usr");
  /*angDist->SetAngDistType("beam2d");
  angDist->SetBeamSigmaInAngX(0.1*mrad);
  angDist->SetBeamSigmaInAngY(0.1*mrad);
  angDist->DefineAngRefAxes("angref1",G4ThreeVector(1.,0.,0.));
  */
}
