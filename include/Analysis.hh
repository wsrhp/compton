#ifndef ANALYSIS_HH
#define ANALYSIS_HH 1

#include "globals.hh"
#include "G4Material.hh"

#include "G4Step.hh"
#include "G4Run.hh"
#include "G4THitsMap.hh"
#include <time.h>
#ifdef G4ANALYSIS_USE_ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#endif
class G4Material;
class G4Run;
class G4Event;
class G4ParticleDefinition;
class TFile;
class TTree;
class G4Timer;
class Analysis
{
public:
    //好方法，通过定义一个静态的analysis对象，让整个程序的分析都在这个对象的调用下实现！
    //妈妈再也不用担心局部变量的问题啦！
    static Analysis* GetInstance();
    ~Analysis();
    Analysis();
    void StepAction(const G4Step* astep);
    void BeginOfEvent();
    void EndOfEvent(const G4Event* anEvent);
    void BeginOfRun(const G4Run* aRun);
    void EndOfRun(const G4Run* aRun);
    void ReportRunTime(G4int npercent);

private:
    //注意！！！！将分析类的构造函数定义为私有，就不能再在类外构造对象了
    // 只能在类内部构造对象。这就是singleton的设计模式
//    Analysis();
    static Analysis* singleton;

    G4double particledata[4];
//    particledata[0]标示粒子能量；单位
//    particledata[1,2,3]标示粒子出射方向

    time_t timestart;
    struct tm tmstart;
    char timebuf[80];

    G4long numOfEventnow;
    G4long numOfEventtotal;
    G4int percentOfEvent;//目前运行了多少百分比的event，以%为单位

    // ROOT objects
    G4String  RootfileName;
    TFile*        myROOTfile;
//    TTree*     spectrumtree;
    TTree*     pelectree;
    TTree*     pgamatree;
    TTree*     selectree;
    TTree*     sgamatree;
    TTree*     spositree;
    G4bool    dataout[4];
};

#endif /* ANALYSIS_HH */
