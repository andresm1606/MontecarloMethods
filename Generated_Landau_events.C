#include <iostream>
#include "TMath.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TLatex.h"
#include <vector>

using namespace std;


vector <Double_t> Metropolis(Int_t N,Double_t xmin,Double_t xmax, Double_t mu, Double_t c,TH1D* sample){

  //Implementación de algoritmo de Metropolis
  vector <Double_t>  M;
  TRandom *g = new TRandom();
  Double_t x0 = (xmax-xmin)*(g->Uniform(1.0)) + xmin;
  
  M.push_back(x0);
  
  while(M.size() <  N){	

    //if (TMath::Landau(M.back(),mu,c,kTRUE) == 0) continue;
    Double_t xprime = (xmax-xmin)*(g->Uniform(1.0)) + xmin;
    Double_t d_S = -TMath::Log(TMath::Landau(xprime,mu,c,kTRUE)/TMath::Landau(M.back(),mu,c,kTRUE));
		
    if(d_S < 0) M.push_back(xprime);
  
    else if(d_S > 0){

      Double_t r = g->Uniform(1.0);
      Double_t p = TMath::Landau(xprime,mu,c,kTRUE)/TMath::Landau(M.back(),mu,c,kTRUE);
      
      if( r < p) M.push_back(xprime);
       
      
      else if(r > p) continue;
    }
  }

  for(Int_t l =0; l < N; l++){
       
    sample->Fill(M.at(l));

  }
  
  return M;
}

void Generated_Landau_events()
{

  TChain * ch = new TChain("SystemTree","");

  ch->Add("AndresMunozAcevedo.root"); // Datos base de archivo .root

  TTree* t = (TTree*) ch;

  Int_t N = t->GetEntries();

  cout << "Entries : " << N << endl;

  Float_t  x;
  Int_t nbins = 200;
  t->SetBranchAddress("x",&x);

  TH1D* hist = new TH1D("Histograma","Histograma",nbins,0,2500);  // Histograma con datos base 

  for(Int_t i =0; i < N; i++){

    t->GetEntry(i);
    
    hist->Fill(x);
 
  }

   // algoritmo de metropolis y ajuste chi2

  vector <Double_t> mu; // "Media"
  vector <Double_t> c ; // "Desviación estandar"

  for(int i = 0; i<=20;i++){
    // vectores con valores donde se espera encontrar el menor chi2 
    mu.push_back(232.1 + 0.5*i);  // Se usó este rango sabiendo los valores del código anterior para mu y c 
    c.push_back(5.42 + 0.5*i);
   
  }
 
  Double_t xmin = 0;
  Double_t xmax = 2500; // Rango donde los datos visualmente se aprecian bien 
  Int_t n = mu.size();

  vector <Double_t> bsample;
  vector <Double_t> chi2 ;
  vector <Double_t> minimo ;
  
  TH1D* sample = new TH1D("events","events",nbins,0,2500);
  TH1D* best_sample = new TH1D("Generated events","Generated events",nbins,0,2500);
  // Iteraciones para implementar el método de chi2
  for(int i=0; i < n; i++ ){

    if (i%2==0) cout << 10*i/2 << "%-"<< flush;
    if(i == n -1) cout << endl;
    for(int j=0; j < n; j++){ 
 
      
      vector <Double_t> M = Metropolis(N,xmin,xmax,mu.at(i),c.at(j),sample);
      
      Double_t suma = 0;
      
      for(int k = 0; k < nbins; k++){

	if (hist->GetBinContent(k) == 0) continue;
	suma += (hist->GetBinContent(k) - sample->GetBinContent(k))*(hist->GetBinContent(k) - sample->GetBinContent(k))/(hist->GetBinContent(k));
	
      }
      
      chi2.push_back(suma);
     
      if (suma <= *min_element(chi2.begin(), chi2.end())){
	bsample.clear();
	bsample = M;
	minimo.clear();
	minimo.push_back(mu.at(i));
	minimo.push_back(c.at(j));
	  
	
      }

      sample->Reset("ICES");
      M.clear();
    }

  }
  
  
  cout << "Valores óptimos :" << endl;

  cout << "Mu = " << minimo.at(0) << endl;

  cout << "c = " << minimo.at(1) << endl;

  cout << "chi2 = " << (*min_element(chi2.begin(), chi2.end())/nbins) << endl;
 
  
   
  for(Int_t p =0; p < N; p++){
	  
    best_sample->Fill(bsample.at(p));
  }

  
  TCanvas *canvas = new TCanvas("canvas","canvas",50,50,800,600);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(0);

  //hist->SetTitle("Events from the tree");
  hist->GetXaxis()->SetTitle("x");
  hist->GetYaxis()->SetTitle("events");
  //hist->SetLineColor(4);
  hist->Draw();
  canvas->Update();
  // canvas2->SetFillColor(kRed);
  // best_sample->SetTitle("Generated Events with Metropolis");
  //best_sample->GetXaxis()->SetTitle("x");
  best_sample->SetLineColor(kRed);
  best_sample->Draw("same");

  TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);

  legpar->SetTextSize(0.035); //text size in pixels                                 
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);

  legpar->AddEntry(hist,"Original Data","f");
  legpar->AddEntry(best_sample,"Metropolis algorithm","f");
  legpar->AddEntry("",Form("#mu = %1.1f",minimo.at(0)),"");
  legpar->AddEntry("",Form("c = %1.1f",minimo.at(1)),"");
  legpar->AddEntry("",Form("#chi^{2}/nbins = %1.1f",(*min_element(chi2.begin(), chi2.end())/nbins)),"");
  legpar->Draw();
}
