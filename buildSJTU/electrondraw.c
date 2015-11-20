#include <stdlib.h>
#include <iostream>
#include<iomanip>
#include <fstream>
#ifndef __CINT__
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#endif
#include <math.h>

#include <stdio.h>




//注意，主函数名称和文件名称要一样
int electrondraw(void)
{
/*	double a;
    scanf("%lf",&a);//double类型必须要用%lf
    printf("%2.1lf",a);
    char strtemp[10];
    sprintf(strtemp,"%2.1lf",a);//"%2.1f"用于输出格式说明，把数字转化为字符
    string en=strtemp;*/
	///////////////////////////////////////////////////////////读入ROOT文件/////////////////////////////////////////////////////////////////////
    const double PI = 3.1415926535897932384626433832795;
    string filename="spectrum.root"; 
    string filename1="e-spectrum.txt";                                                       //数据文件名称
    string treename="pelectree";                                                            //root文件中tree的名称，需要和Geant4中的一致
    string branchname = "particledata";                                                 //root文件中branch的名称，需要和Geant4中的一致
     gROOT->Reset();
    cout<<endl;
    cout << "\nStart\n" <<endl;
    TFile *myfile = new TFile(&filename[0]);    //新建对象，指向数据文件
if (!myfile)    {
        cout <<"no input file!" <<endl;
        return;
    }
    cout<<"Open file: "<<filename<<endl;
    TTree *mytree = (TTree *)myfile->Get(&treename[0]);//新建tree,用于处理数据
    if (!mytree)    {
        cout <<"no  tree named:  "<<treename<<endl;
        return;
    }
    cout << "input tree: "<<treename <<endl;
    ///////////////////////////////////////////////////////////能谱统计/////////////////////////////////////////////////////////////////////////////
    //-----------------------------------------参数设置-----------------------------------------------
    double Ebins=6000;   //统计能谱的直方图中bin的个数
    double Elow=0;      //能谱的低能截止,单位keV(对于光子，要有一定的低能截止，因为低能的光子会非常多，而对于正负电子，可以从0开始统计)
    double Eup=600000; //能谱的高能能截止,单位keV
     double eneconv=1000;//1000->MeV,1->keV将输入的能量转换单位

    
    //----------------------------------------直方图---------------------------------------------------
    double data[4]={0};//保存从root文件中取出的数据。
    mytree->SetBranchAddress(&branchname[0],data);//"&"取地址内容，这里给branch里的数据设置了新的指针变量
    //create two histograms
    double EEup=Eup/eneconv;
    double EElow=Elow/eneconv;
    TH1D * h1     = new TH1D("E","Energy Spectrum of source",Ebins,EElow,EEup);             //统计能谱
    Long64_t nentries = mytree->GetEntries();   //注意：因为事例entry很多，所以需要Long64_t
    for (Long64_t i=0;i<nentries;i++) {
        mytree->GetEntry(i);
        //将处理数据的tree指向数据文件中的第i个事例，由于上面设置了的tree中branch的指针为变量data，所以处理data变量就处理了数据。
        h1->Fill(data[0]);     
        }                 //Fill(数据)==统计该数据的个数直方图
        
  //-----------------------------------------画图----------------------------------------------------------
    gStyle->SetPalette(1);        // use precomputed color palette 1
    gStyle->SetOptStat("ner");
    
    TCanvas *c5 = new TCanvas();
    c5->cd();
    h1->Draw("C");  //energy
    h1->SetTitle("Energy spectrum");   //4PI空间出射的次级粒子的能谱分布
    h1->GetXaxis()->SetTitle("Energy(MeV)");
    h1->GetYaxis()->SetTitle("Intensity)");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    gPad->Update();
    
  //----------------------------------------直方图数据输出-------------------------------------------
  cout<<"write data file to disk"<<endl;
 // ofstream ofile1("electron_energy.txt",ios_base::out);//输出数据文件
   ofstream ofile1(&filename1[0],ios_base::out);
   for(i=1;i<(Ebins+1);i++)    {
        ofile1<<h1->GetXaxis()->GetBinCenter(i)<<"\t"<< h1->GetBinContent(i)<<"\t"<<endl;
    }                              //读取数据保存为文本形式
    ofile1.close();
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
     cout<<"\nEnd!\n"<<endl;
    return 0;
}

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
