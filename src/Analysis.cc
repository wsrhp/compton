//// Make this appear first!说明中要求把G4Timer.hh放在头文件列表中的第一个位置
#include "G4Timer.hh"

#include "Analysis.hh"
#include <stdlib.h>
#include <iostream>
#include<iomanip>
#include <fstream>

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "G4HCofThisEvent.hh"
#include "G4ProductionCuts.hh"
#include "G4VProcess.hh"
#include "G4StepPoint.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#ifdef G4ANALYSIS_USE_ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#endif

using namespace std;

Analysis* Analysis::singleton = 0;
//控制Analysis起作用的一个静态变量
//Analysis是一个单例类

Analysis* Analysis::GetInstance()
{
    //GetInstance被其他函数调用时就实际是调用一个Analysis的对象analysis.
    if ( singleton == 0 )
    {
        static Analysis analysis;
        singleton = &analysis;
    }
    return singleton;
}

Analysis::Analysis() 
{
    //一定注意啦！要对各个成员数据进行初始化！！！！
    //TotalEnergy = 0;
    myROOTfile = 0;
     pelectree =0;
     pgamatree=0;
    selectree=0;
    sgamatree=0;
    RootfileName="spectrum.root";
    percentOfEvent=1;//必须将其初始化为1
    numOfEventnow=0;
    numOfEventtotal=0;
    for(G4int i=0;i<4;i++){
        particledata[i]=0.;
    }
    dataout[0] = false;
    dataout[1] = true;
    dataout[2] = true;
    dataout[3] = true;  //需要修改这里的布尔值，确定每次保存gamma还是e-或e+

}

Analysis::~Analysis() 
{}
void Analysis::BeginOfEvent()
{}

void Analysis::BeginOfRun(const G4Run* aRun )
{   timestart = time(NULL); //获取run开始时间
   //struct tm* Tmstart = localtime(&timestart);
   //tmstart = *Tmstart;
	tmstart = *(localtime(&timestart));
    strftime(timebuf,80,"%Y-%m-%d %H:%M:%S",&tmstart);  //将时间以2013-08-16 21:59:23的形式输出
    G4cout<<"Starting Run at   " <<timebuf<<G4endl;
    G4cout<<"Starting Run: "<<aRun->GetRunID()<<G4endl;
    G4cout<<"The Number of Event to be processed is :  "<<aRun->GetNumberOfEventToBeProcessed()<<G4endl;
    numOfEventtotal=aRun->GetNumberOfEventToBeProcessed();
    if(numOfEventtotal == 0) { return; }

    //创建保存数据的ROOT文件
#ifdef G4ANALYSIS_USE_ROOT
    G4String rootfilediscripution="ROOT file of beam hit target";
    myROOTfile = new TFile(RootfileName,"RECREATE",rootfilediscripution);
    if(myROOTfile)
    {
        G4cout << "ROOT file "<<RootfileName<<" is created " << G4endl;
    }
    else
    {
        G4Exception("Analysis::PrepareNewRun()", "  ", FatalException, "ROOT file **.root has not been created!");
    }
    //    spectrumtree = new  TTree("spectree","spectrum tree ");
    //    spectrumtree->Branch("particledata",particledata,"particledata[5]/D");
    if (dataout[0]){  pgamatree = new  TTree("pgamatree","Primary Gamma Tree ");
        pgamatree->Branch("particledata",particledata,"particledata[4]/D");}
    if (dataout[1]){ pelectree = new  TTree("pelectree","Primary Electron Tree ");
        pelectree->Branch("particledata",particledata,"particledata[4]/D");}    
    if (dataout[2]){sgamatree = new  TTree("sgamatree","Secendary Gamma Tree ");
        sgamatree->Branch("particledata",particledata,"particledata[4]/D");}
    if (dataout[3]){ selectree = new  TTree("selectree","Secendary Electron Tree ");
        selectree->Branch("particledata",particledata,"particledata[4]/D");}

#endif
}

void Analysis::StepAction(const G4Step* astep)
{
    // save information to ROOT
#ifdef G4ANALYSIS_USE_ROOT
    G4Track* track =  astep->GetTrack();
    if(track->GetTrackStatus()!=fAlive) { return; }
    //需要判断粒子是否已经死掉，如果死掉就不会有当前的各种状态了

    // get volume of the current step
       const G4StepPoint* pre = astep->GetPreStepPoint();
       const G4StepPoint* post = astep->GetPostStepPoint();

       G4String volPre = pre->GetPhysicalVolume()->GetName();
       G4String volPost = post->GetPhysicalVolume()->GetName();
     if(volPost=="Ctube"){
         track->SetTrackStatus(fKillTrackAndSecondaries);
         //track->SetTrackStatus(fKillTrackAndSecondaries);
     }

    G4String  particlename=  track->GetParticleDefinition()->GetParticleName();//得到粒子名称

    particledata[0] = post->GetKineticEnergy()/(1.0*keV);//得到粒子动能
    //G4ThreeVector direction = post->GetMomentumDirection();//得到粒子动量方向
    G4ThreeVector position = post->GetPosition();
   // particledata[1] = direction.getX();
   // particledata[2] = direction.getY();
   // particledata[3] = direction.getZ();


    particledata[1] = position.getX();
    particledata[2] = position.getY();
    particledata[3] = position.getZ();


    if (track->GetParentID() == 0){
        if(particlename=="gamma" && dataout[0]){   pgamatree->Fill(); }
        if(particlename=="e-" && dataout[1]&& (volPre=="Mbox" && volPost=="IPDetector")){
            //TotalEnergy += particledata[0];
            static ofstream ofile("pe-position.txt",ios_base::out);
            ofile<<particledata[1]<<"\t"<<particledata[2]<<"\t"<<particledata[3]<<"\t"<<endl;
            pelectree->Fill();}
    }
    else
    {
        //G4String  process = track->GetCreatorProcess()->GetProcessName();//得到此轨道产生的物理过程
        //ofstream ofile("process.txt",ios_base::out);
        //ofile<<process<<"\t"<<endl;
        if(particlename=="e+" && dataout[2]&& (volPre=="Mbox" && volPost=="IPDetector1")){
            static ofstream ofile("e+position.txt",ios_base::out);
            ofile<<particledata[1]<<"\t"<<particledata[2]<<"\t"<<particledata[3]<<"\t"<<endl;
            sgamatree->Fill();}
        if(particlename=="e-"&& dataout[3]&& (volPre=="Mbox" && volPost=="IPDetector")){
            static ofstream ofile("e-position.txt",ios_base::out);
            ofile<<particledata[1]<<"\t"<<particledata[2]<<"\t"<<particledata[3]<<"\t"<<endl;
            selectree->Fill();}
    }
#endif
}

void Analysis::EndOfEvent(const G4Event* /*anEvent*/)
{
    numOfEventnow+=1;
    //输出运行时间
    if (numOfEventnow==numOfEventtotal*percentOfEvent/100)
    {
        if(percentOfEvent==1)
        {
            ReportRunTime(percentOfEvent);//第一次输出运行时间
            percentOfEvent=10;
        }
        else
        {
            ReportRunTime(percentOfEvent);//输出运行时间及剩余时间
#ifdef G4ANALYSIS_USE_ROOT
            G4cout << "ROOT: files writing  " <<percentOfEvent<<" %"<< G4endl;
            myROOTfile->Write();
#endif
            percentOfEvent+=10;   //这个不能改成20，会出现只运行90%事例数的bug
            //另外曾出现数目超过21474836后显示进度失效，因为需将numOfEventnow和numOfEventtotal设为G4long
        }
    }

}

void Analysis::EndOfRun(const G4Run* /*aRun*/)
{
    // Writing and closing the ROOT file
#ifdef G4ANALYSIS_USE_ROOT
    G4cout << "ROOT: files writed." << G4endl;
    //    myROOTfile->Write();
    G4cout << "ROOT: files closed." << G4endl;
    myROOTfile->Close();
    delete myROOTfile;
#endif
    //ofstream ofile("TotalEnergy.txt",ios_base::out);
    //ofile<<"Total Energy:"<<TotalEnergy<<"keV"<<endl;

    time_t timeend;
    struct tm* tmend;
    timeend = time(NULL); //获取run开始时间
    tmend=localtime(&timeend);
    strftime(timebuf,80,"%Y-%m-%d %H:%M:%S",tmend);  //将时间以2013-08-16 21:59:23的形式输出
    G4cout<<" Run end at   " <<timebuf<<G4endl;
}

void Analysis::ReportRunTime(G4int percent)
{
    //利用C++中的时间函数得到程序运行至某percent of events 的时间和剩余时间
    if(percent<0 ||percent>100)
    {
        G4cout<<"Error in Analysis::ReportRunTime(): illegal input"<<G4endl;
        return;
    }
    G4cout<<"Geant4 has processed "<<percent<<" % events."<<G4endl;
    G4double timecost,timeneed;
    G4double npercent=(G4double) percent;
    G4int dayneed,hourneed,minneed,secneed;

    time_t timenow;
    struct tm *tmnow;
    timenow = time(NULL); //获取当前日历时间
    tmnow=localtime(&timenow);
    strftime(timebuf,80,"%Y-%m-%d %H:%M:%S",&tmstart);  //将时间以2013-08-16 21:59:23的形式输出
    G4cout<<"Starting Run at   " <<timebuf<<G4endl;
    strftime(timebuf,80,"%Y-%m-%d %H:%M:%S",tmnow);  //将时间以2013-08-16 21:59:23的形式输出
    G4cout<<"Time Now is  " <<timebuf<<G4endl;
    //求出时间差。
    timecost=difftime(timenow,timestart);
    timeneed=(100-npercent)/npercent*timecost;
    dayneed=(G4int)(timeneed/(24*3600));
    hourneed=(G4int)((timeneed-dayneed*(24*3600))/3600);
    minneed=(G4int)((timeneed-dayneed*(24*3600)-hourneed*3600)/60);
    secneed=(G4int)(timeneed-dayneed*(24*3600)-hourneed*3600-minneed*60);
    G4cout<<"Simulation need more:  "<<dayneed<<" day "<<hourneed<<" hour "
         <<minneed<<" min "<<secneed<<" sec "<<G4endl;
    //预计结束时间
    int ttemp=0;
    ttemp= tmnow->tm_sec+secneed;
    if(ttemp>60){
        minneed+=1;
        tmnow->tm_sec=ttemp-60;
    }
    else{
        tmnow->tm_sec+=secneed;
    }
    ttemp=  tmnow->tm_min+minneed;
    if(ttemp>60){
        hourneed+=1;
        tmnow->tm_min=ttemp-60;
    }
    else{
        tmnow->tm_min+=minneed;
    }
    ttemp=   tmnow->tm_hour+hourneed;
    if(ttemp>24){
        dayneed+=1;
        tmnow->tm_hour=ttemp-60;
    }
    else{
        tmnow->tm_hour+=hourneed;
    }
    tmnow->tm_mday+=dayneed;
    strftime(timebuf,80,"%Y-%m-%d %H:%M:%S",tmnow);
    G4cout<<"Simulation should end at "<<timebuf<<G4endl;
}
