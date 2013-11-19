#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TH1.h"
#include <TMath.h>

using namespace std;
void Plot_pool() {
  gROOT->ProcessLine(".L ./tdrstyle.C");
  setTDRStyle();


//TFile *fMuonError_r12 = new TFile("./muonHistosR12.root");
//TH1F *Poolr12 = (TH1F*)fMuonError_r12->Get("demo/DeltaX_over_err");
//TH1F *Pool = (TH1F*)gDirectory->Get(demo/DeltaX);
//Poolr12->SetFillColor(1);
//Poolr12->SetFillStyle(3001);
TFile *fMuonError_12 = new TFile("./muonHistosOldErrors.root");
TH1F *Pool12 = (TH1F*)fMuonError_12->Get("demo/DeltaX_over_err");
Pool12->SetFillColor(2);
Pool12->SetFillStyle(3001);
TFile *fMuonError_cs = new TFile("./muonHistosNewErrors.root");
TH1F *Poolcs = (TH1F*)fMuonError_cs->Get("demo/DeltaX_over_err");
Poolcs->SetFillColor(3);
Poolcs->SetFillStyle(3001);


TCanvas *Canv = new TCanvas("Canv"," ");
//Poolcs->Draw("HIST ");
Pool12->Draw("HIST ");
//Poolr12->Draw("HIST ");
Poolcs->Draw("HIST same");

TLegend* legend= new TLegend( 0.67, 0.6, 0.87, 0.92);
legend->SetTextSize(0.025);
legend->SetFillColor(0);
//legend->AddEntry(Poolr12, "(local_error)^2 = 1/Sqrt(12)","f");
legend->AddEntry(Pool12, "(local_error)^2 = 1/12","f");
legend->AddEntry(Poolcs, "(local_error)^2 = (cluster_size)^2/12)","f");

legend->Draw("same");

}
