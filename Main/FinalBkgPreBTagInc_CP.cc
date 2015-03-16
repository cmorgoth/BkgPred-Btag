/*
Final Background Prediction, June 2014
Includes bug fixed for the kfactors
*/
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TH2F.h"
#include "DM_1DRatio.hh"
#include "DM_2DRatio.hh"
#include "DM_Base.hh"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"
#include <vector>

const int r2Bins = 5;
const float BaseDM::RSQ_BinArr[] = {0.5, 0.6, 0.725, 0.85,  2.50};
const float BaseDM::MR_BinArr[] = {200., 300., 400., 600.,  3500.};

//MR Categories
const int r2B[4] = {11, 6, 6, 4};
float c1B[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0, 1.2};
float c2B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
float c3B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
float c4B[] = {0.50, 0.60, 0.70, .950, 1.2};


const int r2B_tt[4] = {4, 4, 4, 4};
float c1B_tt[] = {0.50, 0.6, 0.75, .9, 1.2};
float c2B_tt[] = {0.50, 0.575, 0.65, 0.8, 1.2};
float c3B_tt[] = {0.50, 0.59, 0.7, 0.85, 1.2};
float c4B_tt[] = {0.50, 0.59, 0.7, 0.85, 1.2};

std::vector<float*> v;
std::vector<float*> v_tt;


void set_plot_style(){
  const int NRGBs = 5;
  const int NCont = 255;
  
  double stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  double red[NRGBs]   = { 0.50, 0.70, 1.00, 1.00, 1.00 };
  double green[NRGBs] = { 0.50, 0.70, 1.00, 0.70, 0.50 };
  double blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.70, 0.50 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

int main(){
  //gROOT->Reset();

  double tt_k = 1.776;
  double z_k = 1.1973;
  double w_k = 1.2401;
  
  set_plot_style();
  TCanvas* cc = new TCanvas("cc", "cc", 640, 640);

  //Pushing Binning
  v.push_back(c1B);
  v.push_back(c2B);
  v.push_back(c3B);
  v.push_back(c4B);
  //Pushing Binning
  v_tt.push_back(c1B_tt);
  v_tt.push_back(c2B_tt);
  v_tt.push_back(c3B_tt);
  v_tt.push_back(c4B_tt);

  double tt_k2factor = 1.0;
  
  ///////////////////////////////////////////////////////
  ///////////////Creating MR categories binned in R2/////
  //////////////////////////////////////////////////////
  
  //2mu-1b
  //TFile* F = new TFile("FinalROOTFiles/OneBtag_Loose_May_2014_NewTriger_AN_4B_Mcut_Z.root");
  
  //TFile* F = new TFile("FinalROOTFiles/OneBtag_Loose_May_2014_NewTriger_AN_NoBtagCorr_MZcut.root");
  //TFile* F = new TFile("FinalROOTFiles/OneLBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE.root");
  TFile* F = new TFile("FinalROOTFiles/OneLBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE_UNCORR_BTAG.root");
  
  TH1F* dy_2mu;
  TH1F* z_2mu;
  TH1F* w_2mu;
  TH1F* tt_2mu;
  TH1F* data_2mu;
  
  //Getting 2mu-1LooseB Histos
  TString name1 = "DY_INC_1D_2mu_Box";
  dy_2mu = (TH1F*)F->Get(name1);
  dy_2mu->Scale(z_k);
  name1 = "Z_INC_1D_2mu_Box";
  z_2mu = (TH1F*)F->Get(name1);
  z_2mu->Scale(z_k);
  name1 = "W_INC_1D_2mu_Box";
  w_2mu = (TH1F*)F->Get(name1);
  w_2mu->Scale(w_k);
  //name1 = "TT_L_INC_1D_2mu_Box";
  name1 = "TT_INC_1D_2mu_Box";
  tt_2mu = (TH1F*)F->Get(name1);
  tt_2mu->Scale(tt_k);
  name1 = "Data_INC_1D_2mu_Box";
  data_2mu = (TH1F*)F->Get(name1);
  
  /////////////////////////
  ////////dy-2mu pred//////
  ////////////////////////
  TH1F* pred_2mu = new TH1F(*data_2mu);
  pred_2mu->Add(tt_2mu, -1.0);
  pred_2mu->Add(w_2mu, -1.0);
  pred_2mu->Add(z_2mu, -1.0);


  //TFile* Ftt = new TFile("FinalROOTFiles/OneTBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE.root");
  
  TFile* Ftt = new TFile("FinalROOTFiles/OneTBtag_ptReweighting_Inclusive.root");
  
  TH1F* dy_2mu_tt;
  TH1F* z_2mu_tt;
  TH1F* w_2mu_tt;
  TH1F* tt_2mu_tt;
  TH1F* data_2mu_tt;
  //Getting 2mu-1TB Histos
  name1 = "DY_INC_1D_2mu_Box";
  dy_2mu_tt = (TH1F*)Ftt->Get(name1);
  dy_2mu_tt->Scale(z_k);
  
  name1 = "Z_INC_1D_2mu_Box";
  z_2mu_tt = (TH1F*)Ftt->Get(name1);
  z_2mu_tt->Scale(z_k);

  name1 = "W_INC_1D_2mu_Box";
  w_2mu_tt = (TH1F*)Ftt->Get(name1);
  w_2mu_tt->Scale(w_k);

  //name1 = "TT_L_INC_1D_2mu_Box";
  name1 = "TT_INC_1D_2mu_Box";
  tt_2mu_tt = (TH1F*)Ftt->Get(name1);
  tt_2mu_tt->Scale(tt_k);
  
  name1 = "Data_INC_1D_2mu_Box";
  data_2mu_tt = (TH1F*)Ftt->Get(name1);
  
  TH1F* pred_2mu_tt = new TH1F(*data_2mu_tt);
  pred_2mu_tt->Add(dy_2mu_tt, -1.0);
  pred_2mu_tt->Add(w_2mu_tt, -1.0);
  pred_2mu_tt->Add(z_2mu_tt, -1.0);

  TH1F* dy;
  TH1F* z;
  TH1F* w;
  TH1F* tt;
  TH1F* data;
  
  std::cout << "debug 1.0" << std::endl;
  //Getting 0,1mu-2TightB Histos
    
  //TFile* F1 = new TFile("FinalROOTFiles/OneTBtag_FixBtag_PFNoPu_Sep2014.root");
  //TFile* F1 = new TFile("FinalROOTFiles/OneTBtag_ptReweighting.root");
  

  //TFile* F1 = new TFile("FinalROOTFiles/TwoTBtag_FixBtag_PFNoPu_Sep2014.root");
  TFile* F1 = new TFile("FinalROOTFiles/TwoTBtag_ptReweighting.root");

  
  
  name1 = "DY_INC_1D_0mu_Box";
  dy = (TH1F*)F1->Get(name1);
  dy->Scale(z_k);
  
  name1 = "Z_INC_1D_0mu_Box";
  z = (TH1F*)F1->Get(name1);
  z->Scale(z_k);
  
  name1 = "W_INC_1D_0mu_Box";
  w = (TH1F*)F1->Get(name1);
  w->Scale(w_k);
  
  //name1 = "TT_L_INC_1D_0mu_Box";
  name1 = "TT_INC_1D_0mu_Box";
  tt = (TH1F*)F1->Get(name1);
  tt->Scale(tt_k);
  
  name1 = "Data_INC_1D_0mu_Box";
  data = (TH1F*)F1->Get(name1);
  
  std::cout << "debug 1.1" << std::endl;
  //////////////////////
  /////TT 0MU PRED//////
  //////////////////////
  TH1F* P_TT_0Mu = new TH1F(*pred_2mu_tt);
  TH1F* TT_R_0Mu = new TH1F(*tt);
  TT_R_0Mu->Divide(tt_2mu_tt);
  P_TT_0Mu->Multiply(TT_R_0Mu);
  
  std::cout << "debug 1.2" << std::endl;
  //////////////////////////////
  /////Z 0MU PRED from 2mu//////
  /////////////////////////////
  TH1F* P_Z_0Mu = new TH1F(*pred_2mu);
  TH1F* Z_R_0Mu = new TH1F(*z);
  Z_R_0Mu->Divide(dy_2mu);
  P_Z_0Mu->Multiply(Z_R_0Mu);

  std::cout << "debug 1.3" << std::endl;
  //////////////////////////////
  /////DY 0MU PRED from 2mu//////
  /////////////////////////////
  TH1F* P_DY_0Mu = new TH1F(*pred_2mu);
  TH1F* DY_R_0Mu = new TH1F(*dy);
  DY_R_0Mu->Divide(dy_2mu);
  P_DY_0Mu->Multiply(DY_R_0Mu);

  std::cout << "debug 1.4" << std::endl;
  //////////////////////////////
  /////W 0MU PRED from 2mu//////
  /////////////////////////////
  TH1F* P_W_0Mu = new TH1F(*pred_2mu);
  TH1F* W_R_0Mu = new TH1F(*w);
  W_R_0Mu->Divide(dy_2mu);
  P_W_0Mu->Multiply(W_R_0Mu);

  //Total Prediction
  TH1F* Pred_0Mu = new TH1F(*P_TT_0Mu);//Creating Prediction to add V+J (from 2mu/z(ll))
  Pred_0Mu->Add(P_W_0Mu, 1.0);
  Pred_0Mu->Add(P_Z_0Mu, 1.0);
  Pred_0Mu->Add(P_DY_0Mu, 1.0);

  //Total MC
  TH1F* TotalMC = new TH1F(*dy);
  TotalMC->Add(tt, 1.0);
  TotalMC->Add(z, 1.0);
  TotalMC->Add(w, 1.0);

  //Getting systemtics
  TFile* f_sys = new TFile("Closure_CP.root");
  TH1F* sys = (TH1F*)f_sys->Get("sys");

  std::cout << "debug 1.5" << std::endl;
  //Include Systematics
  TH1F* Pred_0Mu_SYS = new TH1F(*Pred_0Mu);
  double err;
  for(int j = 1; j <= sys->GetNbinsX(); j++){
    err = sqrt( pow(sys->GetBinContent(j)*Pred_0Mu->GetBinContent(j), 2) + pow(Pred_0Mu->GetBinError(j),2) );
    Pred_0Mu_SYS->SetBinError(j,err);
    
    //Individual SYS
    err = sqrt( pow(sys->GetBinContent(j)*P_TT_0Mu->GetBinContent(j), 2) + pow(P_TT_0Mu->GetBinError(j),2) );
    P_TT_0Mu->SetBinError(j,err);
    err = sqrt( pow(sys->GetBinContent(j)*P_W_0Mu->GetBinContent(j), 2) + pow(P_W_0Mu->GetBinError(j),2) );
    P_W_0Mu->SetBinError(j,err);
    err = sqrt( pow(sys->GetBinContent(j)*P_Z_0Mu->GetBinContent(j), 2) + pow(P_Z_0Mu->GetBinError(j),2) );
    P_Z_0Mu->SetBinError(j,err);
    err = sqrt( pow(sys->GetBinContent(j)*P_DY_0Mu->GetBinContent(j), 2) + pow(P_DY_0Mu->GetBinError(j),2) );
    P_DY_0Mu->SetBinError(j,err);
  }
  
  std::cout << "debug 1.6" << std::endl;
  TFile* fout = new TFile("test_CP_2TB_ptReweighting.root", "RECREATE");
  std::cout << "debug 1.7" << std::endl;
  RatioPlotsBandV2( data, Pred_0Mu, "Data  0#mu", "BKg Pred 0#mu", "PredPlotsFinal/Bkg_0mu2TbEXC_CP", "RSQ", r2B_tt[0], v_tt.at(0),1);
  std::cout << "debug 1.7.1" << std::endl;
  RatioPlotsBandV2( data, Pred_0Mu_SYS, "Data  0#mu", "BKg Pred 0#mu", "PredPlotsFinal/Bkg_0mu2TbEXC_SYS_CP", "RSQ", r2B_tt[0], v_tt.at(0),1);
  RatioPlotsBandV2( data, TotalMC, "Data  0#mu", "Total MC 0#mu", "PredPlotsFinal/MC_0mu2TbEXC_CP", "RSQ", r2B_tt[0], v_tt.at(0),1);
  std::cout << "debug 1.8" << std::endl;
  TH1F* tt0mu = new TH1F(*tt);
  tt0mu->Scale(1.0/tt0mu->Integral());
  TH1F* tt2mu_tt = new TH1F(*tt_2mu_tt);
  tt2mu_tt->Scale(1.0/tt2mu_tt->Integral());
  RatioPlots(tt0mu, tt2mu_tt, "tt_0mu", "tt_2mu", "PredPlotsFinal/test_CP", "RSQ");
  
  TString s = "data0mu";
  data->Write(s);
  s = "Pred0mu";
  Pred_0Mu->Write(s);
  
  Pred_0Mu_SYS->Write("Pred0mu_SYS");
  TotalMC->Write("TotalMC");
  tt->Write("tt");
  dy->Write("dy");
  z->Write("z");
  w->Write("w");
  P_W_0Mu->Write("pw");
  P_DY_0Mu->Write("pdy");
  P_Z_0Mu->Write("pz");
  P_TT_0Mu->Write("ptt");
  
  s = "TT_Ratio";
  TT_R_0Mu->Write(s);
  s = "Z0mu_2muDYRatio";
  Z_R_0Mu->Write(s);
  s = "DY0mu_2muDYRatio";
  DY_R_0Mu->Write(s);
  s = "W0mu_2muDYRatio";
  W_R_0Mu->Write(s);

  fout->Close();
  return 0;
  
}
