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
/// \file electromagnetic/TestEm14/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// 
// $Id: PhysicsList.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
 
#include "PhysListEmStandard.hh"
#include "PhysListEmLivermore.hh"
#include "PhysListEmPenelope.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
: G4VModularPhysicsList(),
  fEmPhysicsList(0),
  fMessenger(0)
{
  //G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");
  G4LossTableManager::Instance();
  
  fCurrentDefaultCut   = 1.0*um;
  fCutForGamma         = fCurrentDefaultCut;
  fCutForElectron      = fCurrentDefaultCut;
  fCutForPositron      = fCurrentDefaultCut;

  //G4Cerenkov::SetMaxNumPhotonsPerStep(20);
  //G4Cerenkov::SetMaxBetaChangePerStep(10.0);
  //G4Cerenkov::SetTrackSecondariesFirst(true);

  //G4Scintillation::SetScintillationYieldFactor(1.);
  //G4Scintillation::SetTrackSecondariesFirst(true);

  fMessenger = new PhysicsListMessenger(this);
  fVerboseLebel=1;
  SetVerboseLevel(fVerboseLebel);

  // EM physics
  fEmName = G4String("standard");
  fEmPhysicsList = new PhysListEmStandard(fEmName);
  
  //add new units for cross sections
  // 
  new G4UnitDefinition( "mm2/g", "mm2/g","Surface/Mass", mm2/g);
  new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Mesons
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4EtaPrime.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"

// Baryons
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

// Nuclei
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
// pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
// gamma
  G4Gamma::GammaDefinition();
  
// optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

// leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  

// mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

// barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

// ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmProcessOptions.hh"

void PhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();

  // Electromagnetic physics list
  //
  fEmPhysicsList->ConstructProcess();
  
  // Em options
  //
  G4EmProcessOptions emOptions;
  emOptions.SetIntegral(false);

  // 构建光学过程
  ConstructOp();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4Threading.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"

void PhysicsList::ConstructOp()
{
  G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");
  G4Scintillation* scintillationProcess = new G4Scintillation("Scintillation");
  G4OpAbsorption* absorptionProcess = new G4OpAbsorption();
  G4OpRayleigh* rayleighScatteringProcess = new G4OpRayleigh();
  G4OpMieHG* mieHGScatteringProcess = new G4OpMieHG();
  G4OpBoundaryProcess* boundaryProcess = new G4OpBoundaryProcess();

  cerenkovProcess->SetVerboseLevel(fVerboseLebel);
  scintillationProcess->SetVerboseLevel(fVerboseLebel);
  absorptionProcess->SetVerboseLevel(fVerboseLebel);
  rayleighScatteringProcess->SetVerboseLevel(fVerboseLebel);
  mieHGScatteringProcess->SetVerboseLevel(fVerboseLebel);
  boundaryProcess->SetVerboseLevel(fVerboseLebel);

  // Use Birks Correction in the Scintillation process
  if(!G4Threading::IsWorkerThread())
  {
    G4EmSaturation* emSaturation =
              G4LossTableManager::Instance()->EmSaturation();
    //G4Scintillation::AddSaturation(emSaturation);
  }

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (cerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(cerenkovProcess);
      pmanager->SetProcessOrdering(cerenkovProcess,idxPostStep);
    }
    if (scintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(scintillationProcess);
      pmanager->SetProcessOrderingToLast(scintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(scintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(absorptionProcess);
      pmanager->AddDiscreteProcess(rayleighScatteringProcess);
      pmanager->AddDiscreteProcess(mieHGScatteringProcess);
      pmanager->AddDiscreteProcess(boundaryProcess);
    }

  }
}
void PhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
  //G4Cerenkov::SetMaxNumPhotonsPerStep(MaxNumber);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  
  if (name == fEmName) return;

  if (name == "standard") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandard(name);
   
  } else if (name == "livermore") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmLivermore(name);
    
  } else if (name == "penelope") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmPenelope(name);
                    
  } else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

void PhysicsList::SetCuts()
{
 // fixe lower limit for cut
 G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
  
 // set cut values for gamma at first and for e- second and next for e+,
 // because some processes for e+/e- need cut values for gamma
 SetCutValue(fCutForGamma, "gamma");
 SetCutValue(fCutForElectron, "e-");
 SetCutValue(fCutForPositron, "e+");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
