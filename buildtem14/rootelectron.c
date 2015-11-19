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
#include "TH2.h"sss
#include "TCanvas.h"
#include "TTree.h"
#endif
#include <math.h>



//注意，主函数名称和文件名称要一样
int rootelectron(void)
{
    const double PI = 3.1415926535897932384626433832795;
    ///////////////////////%%%%%%输入参数%%%%%%%%%%%/////////////////////////////////////////////////////////////
    string filename="spectrum.root";                                                        //数据文件名称
    //root文件中tree的名称，需要和Geant4中的一致
    string treename="selectree";
    string branchname = "particledata";                                                 //root文件中branch的名称，需要和Geant4中的一致
    //    string treename[4]={"pelectree","pgamatree","selectree","sgamatree"};
    double Ebins=200;   //统计能谱的直方图中bin的个数，
    double Elow=-75.0;      //能谱的低能截止,单位keV
    double Eup=95.5555673; //能谱的高能能截止,单位keV
    double cbins[3] = {180,180,180};//统计角分布的三个直方图中bin的个数【x,y,z】
    double chlow[3] = {0,-90,0};//统计角分布的三个直方图的左边界【x,y,z】
    double chup[3] = {180,90,90};//统计角分布的三个直方图的左边界【x,y,z】
    double eneconv=1;//1000->MeV,1->keV将输入的能量转换单位

    ///////////////////////%%%%%%读入ROOT文件和Tree%%%%%%%%%%%/////////////////////////////////////////////////////////////
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

    //////////////////////////////%%设置新tree，用于数据处理%%%%%%/////////////////////////////////////%%%%%%%%%%%%%%%

    double data[4]={0};//保存从root文件中取出的数据。
    mytree->SetBranchAddress(&branchname[0],data);//设置新建的tree中branch的指向。
    //create two histograms
    double EEup=Eup/eneconv;
    double EElow=Elow/eneconv;
    double sintheta=0;
    TH1D * h1     = new TH1D("E","Energy Spectrum of source",Ebins,EElow,EEup);             //统计能谱
   // TH2D *h2Ex = new TH2D("h2Ex","h2 Energy vs x",cbins[0],chlow[0],chup[0],Ebins,EElow,EEup);//统计角分布，出射方向与X方向的夹角
   // TH2D *h2Ey = new TH2D("h2Ey","h2 Energy vs y",cbins[1],chlow[1],chup[1],Ebins,EElow,EEup);//统计角分布，出射方向与Y方向的夹角
    TH2D *h2Ez = new TH2D("h2Ez","h2 Energy vs z",cbins[2],chlow[2],chup[2],Ebins,EElow,EEup);//统计角分布，出射方向与Z方向的夹角（横坐标名，纵坐标名，横坐标bins，横坐标低截止，横坐标高截止，纵坐标bins，纵坐标低截止，纵坐标高截止）

    //read all entries and fill the histograms
    Long64_t nentries = mytree->GetEntries();   //注意：因为事例entry很多，所以需要Long64_t
    for (Long64_t i=0;i<nentries;i++) {
        mytree->GetEntry(i);
        //将处理数据的tree指向数据文件中的第i个事例，由于上面设置了的tree中branch的指针为变量data，所以处理data变量就处理了数据。
        h1->Fill(data[3]);                      //统计能谱
      //  sintheta= -sqrt(1-data[3]*data[3]);//data[3]=cos(theta)
       // h2Ex->Fill(acos(data[1]/sintheta)/PI*180,data[0]);
       // h2Ey->Fill(asin(data[2]/sintheta)/PI*180,data[0]);
        h2Ez->Fill(acos(data[1])/PI*180,data[0]);
    }

    ///////////////////////%%%%%%%%%Profile Histogram%%%%%%%/////////////////////////////////////
    TH1D* az;
  //  ax=h2Ex->ProjectionX();   //将二维散点图投影到1维就得到了一维的角分布，
  //  ay=h2Ey->ProjectionX();   //将二维散点图投影到1维就得到了一维的角分布，
    az=h2Ez->ProjectionX();   //将二维散点图投影到1维就得到了一维的角分布，

    ///////////////////////%%%%%%%%%Plot Histogram%%%%%%%/////////////////////////////////////
    //gROOT->SetStyle("Plain");   // uncomment to set a different style
    gStyle->SetPalette(1);        // use precomputed color palette 1
    gStyle->SetOptStat("ner");

 /*   TCanvas *c1 = new TCanvas();
    c1->cd();
    h2Ex->Draw("C");  //energy
    h2Ex->SetTitle("Energy Spectrum");     //散点图
    h2Ex->GetXaxis()->SetTitle("alpha(x)");
    h2Ex->GetYaxis()->SetTitle("Energy (keV)");
    h2Ex->GetXaxis()->CenterTitle();
    h2Ex->GetYaxis()->CenterTitle();
    gPad->Update();
    //    eneh->GetXaxis()->SetRangeUser(0, 1200);
    //    eneh->GetYaxis()->SetRangeUser(1e-8,1e-1);
    //        gPad->SetLogx();
    //            gPad->SetLogy();

    TCanvas *c2 = new TCanvas();
    c2->cd();
    h2Ey->Draw("C");  //energy
    h2Ey->SetTitle("Energy Spectrum");
    h2Ey->GetXaxis()->SetTitle("beta(y)");
    h2Ey->GetYaxis()->SetTitle("Energy (MeV)");
    h2Ey->GetXaxis()->CenterTitle();
    h2Ey->GetYaxis()->CenterTitle();
    gPad->Update();            */
   /*
    TCanvas *c3 = new TCanvas();           
    c3->cd();
    h2Ez->Draw("C");  //energy
    h2Ez->SetTitle("Energy Spectrum");
    h2Ez->GetXaxis()->SetTitle("gamma(z)");
    h2Ez->GetYaxis()->SetTitle("Energy (keV)");
    h2Ez->GetXaxis()->CenterTitle();
    h2Ez->GetYaxis()->CenterTitle();
    gPad->Update();

    TCanvas *c4 = new TCanvas();
    c4->cd();
    az->Draw("C");  //energy
    az->SetTitle("Augler  Distribution");//出射的次级粒子数目与出射角gamma的关系
    az->GetXaxis()->SetTitle("gamma(z)");
    az->GetYaxis()->SetTitle("Intensity)");
    az->GetXaxis()->CenterTitle();
    az->GetYaxis()->CenterTitle();
    gPad->Update();
*/
    TCanvas *c5 = new TCanvas();
    c5->cd();
    h1->Draw("C");  //Particle Distribution
    h1->SetTitle("Particle Distribution");   //IP板上电子沿长边的粒子数分布
    h1->GetXaxis()->SetTitle("X(mm)");
    h1->GetYaxis()->SetTitle("Intensity");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    gPad->Update();
    //    TGraphPolar * grP1 = new TGraphPolar(1000,r,theta);
    //    grP1->SetLineColor(2);
    //    grP1->Draw("AOL");
    //    gPad->Update();

    ///////////////////////%%%%%%%%%Write  Out Put File%%%%%%/////////////////////////////////////
    cout<<"write data file to disk"<<endl;
    // ofstream ofile(&nameofoutfile[0],ios_base::binary);//输出二进制文件
    ofstream ofile1("gamma-e-.txt",ios_base::out);//输出数据文件
    ofstream ofile2("angulargamma.txt",ios_base::out);
    //    第一列为能量，keV
    //    第二列为计数，个
 //  ofstream ofile2("Augler  Distribution(x).txt",ios_base::out);//输出数据文件
    //    第一列为出射角与X轴的夹角，keV
    //    第二列为计数，个
 //   ofstream ofile3("Augler  Distribution(y).txt",ios_base::out);//输出数据文件
 //   ofstream ofile4("Augler  Distribution(z).txt",ios_base::out);//输出数据文件
    //输出数据
    for(i=1;i<(Ebins+1);i++)    {
        double xxx=h1->GetXaxis()->GetBinCenter(i);
        double xxx1=(xxx+75.0)/0.8527778365 ;
        ofile1<<xxx1<<"\t"<< (h1->GetBinContent(i))<<"\t"<<endl;


    }                              //读取数据保存为文本形式
    
    for(j=1;j<(cbins[2]+1);j++)    {
		for(k=1;k<(Ebins+1);k++)   {
			ofile2<<h2Ez->GetXaxis()->GetBinCenter(j)<<"\t"<<h2Ez->GetYaxis()->GetBinCenter(k)<<"\t"<< h2Ez->GetBinContent(j,k)<<"\t"<<endl;
			}
	}

 /*   for(i=1;i<(cbins[0]+1);i++)    {
        ofile2<<ax->GetXaxis()->GetBinCenter(i)<<"\t"<< ax->GetBinContent(i)<<"\t"<<endl;
    }
    for(i=1;i<(cbins[1]+1);i++)    {
        ofile3<<ay->GetXaxis()->GetBinCenter(i)<<"\t"<< ay->GetBinContent(i)<<"\t"<<endl;
    }
    for(i=1;i<(cbins[2]+1);i++)    {
        ofile3<<az->GetXaxis()->GetBinCenter(i)<<"\t"<< az->GetBinContent(i)<<"\t"<<endl;
    }
    ofile1.close();//读写文件后需要关闭文件
 //   ofile2.close();
  //  ofile3.close();
 //   ofile4.close();      */

    cout<<"\nEnd!\n"<<endl;
    return 0;
}

