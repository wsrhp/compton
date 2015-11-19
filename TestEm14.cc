
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4MTHepRandom.hh"

#include "G4Version.hh"
#include "G4VisExecutive.hh"
#if  G4VERSION_NUMBER>=930
#include "G4UIExecutive.hh"
#else
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIQt.hh"
#endif


#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingVerbose.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
int main(int argc,char** argv) {
 
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  long seeds[2];
  time_t systime = time(NULL)*1000;
  seeds[0] = (long) (systime*1000);
  seeds[1] = (long) (systime*G4UniformRand()*1000);
  //seeds[1] = (long) (systime);
  G4int index=(int) (systime*G4UniformRand()); //设置为0-215之间的数，大于215会对215取模
  CLHEP::HepRandom::setTheSeeds(seeds,index);

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  DetectorConstruction* det;
  PrimaryGeneratorAction* prim;
  runManager->SetUserInitialization(det = new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserAction(prim = new PrimaryGeneratorAction(det));
      
  // set user action classes
  RunAction* run;  
  runManager->SetUserAction(run = new RunAction(det,prim)); 
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new SteppingAction);
   

  // Initialize G4 kernel
    runManager->Initialize();  //不运行初始化就可以设置物理模型了

  //Initilize the visualization manager
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {  // batch mode
      //command line contains name of the macro to execute
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
  }
  else {           // interactive mode : define UI session

#if  G4VERSION_NUMBER>=930
      //New since G4 9.3: UI executive setups up
      //correct UI depending on env variables
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      //If UI has graphics execute special macro: opens OpenGL Qt driver
      //if (ui->IsGUI())
        //  UImanager->ApplyCommand("/control/execute visQt.mac");
      //else
          UImanager->ApplyCommand("/control/execute vis.mac");
#else
      //Older versions of G4: UI selected by user
#ifdef G4UI_USE_QT
    G4UIsession * ui = new G4UIterminal(new G4UIQT);
//  #ifdef G4UI_USE_TCSH
//      G4UIsession * ui = new G4UIterminal(new G4UItcsh);
  #else
      G4UIsession * ui = new G4UIterminal();
  #endif
      UImanager->ApplyCommand("/control/execute vis.mac");
#endif
      ui->SessionStart();
      delete ui;
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete runManager;

  return 0;
}
