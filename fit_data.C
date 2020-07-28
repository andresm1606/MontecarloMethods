#include <iostream>
#include "RooDataSet.h"
#include "TMath.h"
#include "TChain.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLegend.h"
#include <vector>

using namespace std;

// AJUSTE DE FUNCIÓN DE LANDAU A HISTOGRAMA DE DATOS BASE 
void fit_data(){

  TCanvas *canvas = new TCanvas("canvas","",50,50,800,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(0);
  
  TChain * ch = new TChain("SystemTree","");

  ch->Add("AndresMunozAcevedo.root");

  TTree* t = (TTree*) ch;

  Int_t N = t->GetEntries();

  cout << "Entries : " << N << endl;

  Float_t  x;
  Int_t nbins = 200;
  t->SetBranchAddress("x",&x);

  TH1D* hist = new TH1D("Histograma","Histograma",nbins,0,2500);

  for(Int_t i =0; i < N; i++){

    t->GetEntry(i);
    
    hist->Fill(x);
 
  }

  TF1 *f = new TF1("fit","[0]*TMath::Landau(x,[1],[2])",0,2500);
  //f->FixParameter(0,N);
  f->SetParNames("Normalization","Mu","c");
  f->SetParameter(1,250);
  f->SetParameter(2,10);
  //f->SetParLimits(0,280,310);
  //f->SetParLimits(1,1,5);
  hist->Fit("fit");
  
  

  // canvas2->SetFillColor(kRed);
  //hist->SetTitle("Fit from a Landau distribution");
  f->SetLineColor(kRed);
  
  hist->GetXaxis()->SetTitle("x");
  hist->GetYaxis()->SetTitle("events");
  hist->Draw();
  f->Draw("same");
  
  Double_t mu = f->GetParameter(1);
  Double_t c = f->GetParameter(2);
  
  TLegend *leg = new TLegend(0.17,0.73,0.92,0.88);
  //TLegend *leg = new TLegend(0.58,0.4,0.85,0.6);// sin parametros
  leg->SetTextSize(0.04);
  //leg->SetHeader("CMS preliminary");                                                                                                 
  leg->SetFillColor(2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("fit","Fit","l");
  leg->AddEntry(hist,"Histogram filled with the original data","f");
  leg->AddEntry("",Form("#mu = %1.1f",mu),"");
  leg->AddEntry("",Form("c = %1.1f",c),"");
  
  leg->Draw();
  canvas->Update();
		   

  cout << "Parameters:" << endl;
  cout << "mu = " << mu << endl;
  cout << "c = " << c << endl;



}
