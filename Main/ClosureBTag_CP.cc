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


const int r2B_tt[4] = {4, 5, 5, 5};

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

  //Z(ll)-2mu-1Lb
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
  ////////Z(ll)-2mu pred//////
  ////////////////////////
  TH1F* pred_2mu = new TH1F(*data_2mu);
  pred_2mu->Add(tt_2mu, -1.0);
  pred_2mu->Add(w_2mu, -1.0);
  pred_2mu->Add(z_2mu, -1.0);
  
  //TFile* Ftt = new TFile("FinalROOTFiles/OneBtag_Tight_May_2014_NewTriger_AN_4B.root");
  //TFile* Ftt = new TFile("FinalROOTFiles/OneBtag_Tight_May_2014_NewTriger_AN_NoBtagCorr.root");
  TFile* Ftt = new TFile("FinalROOTFiles/OneTBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE.root");
  ///////////////////////////////////////////////////////
  ///////////////Creating MR categories binned in R2/////
  //////////////////////////////////////////////////////	

  TH1F* dy;
  TH1F* z;
  TH1F* w;
  TH1F* tt;
  TH1F* data;

  TH1F* dy_2mu_tt;
  TH1F* z_2mu_tt;
  TH1F* w_2mu_tt;
  TH1F* tt_2mu_tt;
  TH1F* data_2mu_tt;
  ///////////////////////////////
  ////////////2mu-1Tb////////////
  ////////////////////////////////
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
  
  //////////////////////////
  //////////1mu-Files///////
  //////////////////////////
  //TFile* F1 = new TFile("FinalROOTFiles/TwoBtag_Tight_May_2014_NewTriger_AN_4B_McutV3.root");
  //TFile* F1 = new TFile("FinalROOTFiles/OneBtag_Tight_May_2014_NewTriger_AN_4B.root");
  
  //TFile* F1 = new TFile("FinalROOTFiles/TwoBtag_Tight_May_2014_NewTriger_AN_NoBtagCorr.root");

  TFile* F1 = new TFile("FinalROOTFiles/OneTBtag_FixBtag_PFNoPu_Sep2014_INCLUSIVE.root");
  
  /////////////////////////////////
  ///////////1mu-1Tb///////////////
  /////////////////////////////////
  name1 = "DY_INC_1D_1mu_Box";
  dy = (TH1F*)F1->Get(name1);
  dy->Scale(z_k);
  name1 = "Z_INC_1D_1mu_Box";
  z = (TH1F*)F1->Get(name1);
  z->Scale(z_k);
  name1 = "W_INC_1D_1mu_Box";
  w = (TH1F*)F1->Get(name1);
  w->Scale(w_k);
  //name1 = "TT_L_INC_1D_1mu_Box";
  name1 = "TT_INC_1D_1mu_Box";
  tt = (TH1F*)F1->Get(name1);
  tt->Scale(tt_k);
  name1 = "Data_INC_1D_1mu_Box";
  data = (TH1F*)F1->Get(name1);
  
  //1Mu Closure Test
  TH1F* MC_1MU = new TH1F(*tt);
  MC_1MU->Add(dy, 1.0);
  MC_1MU->Add(z, 1.0);
  MC_1MU->Add(w, 1.0);
  
  TH1F* MC_2MU = new TH1F(*tt_2mu_tt);
  MC_2MU->Add(dy_2mu_tt, 1.0);
  MC_2MU->Add(z_2mu_tt, 1.0);
  MC_2MU->Add(w_2mu_tt, 1.0);
  
  //Extra contribution from W
  TH1F* W_extra_1mu = new TH1F(*w);
  W_extra_1mu->Multiply(pred_2mu);
  W_extra_1mu->Divide(dy_2mu);
  
  //Extra contribution from Z(ll)
  TH1F* DY_extra_1mu = new TH1F(*dy);
  DY_extra_1mu->Multiply(pred_2mu);
  DY_extra_1mu->Divide(dy_2mu);
  
  //Extra contribution from Z(nunu)
  TH1F* Z_extra_1mu = new TH1F(*z);
  Z_extra_1mu->Multiply(pred_2mu);
  Z_extra_1mu->Divide(dy_2mu);
  
  /////////////////////////////
  //////////TT PRED 1MU////////
  /////////////////////////////
  TH1F* P_TT_1Mu = new TH1F(*pred_2mu_tt);
  TH1F* TT_R_1Mu = new TH1F(*tt);  
  TT_R_1Mu->Divide(tt_2mu_tt);
  P_TT_1Mu->Multiply(TT_R_1Mu);
  
  //////////////////////////////
  //////Total 1Mu///////////////
  //////////////////////////////
  TH1F* Pred_1Mu = new TH1F(*P_TT_1Mu);
  Pred_1Mu->Add(W_extra_1mu);
  Pred_1Mu->Add(DY_extra_1mu); 
  Pred_1Mu->Add(Z_extra_1mu);
  
  ///////////////////////////////
  //////Total MC 1Mu-1Tb//////////
  //////////////////////////////
  TH1F* MC_1Mu1Tb = new TH1F(*dy);
  MC_1Mu1Tb->Scale(1.1973);
  MC_1Mu1Tb->Add(tt, 1.776);
  MC_1Mu1Tb->Add(z, 1.1973); 
  MC_1Mu1Tb->Add(w, 1.2401);
  
  ///////////////////////////////
  //////Total MC 2Mu-TT//////////
  //////////////////////////////
  TH1F* MC_2MuTT = new TH1F(*dy_2mu_tt);
  MC_2MuTT->Add(tt_2mu_tt, 1.0);
  MC_2MuTT->Add(z_2mu_tt, 1.0); 
  MC_2MuTT->Add(w_2mu_tt, 1.0);

  ///////////////////////////////
  //////Total MC 2Mu-Z//////////
  //////////////////////////////
  TH1F* MC_2MuZ = new TH1F(*tt_2mu);
  MC_2MuZ->Add(dy_2mu, 1.0);
  MC_2MuZ->Add(z_2mu, 1.0); 
  MC_2MuZ->Add(w_2mu, 1.0);
  
  //Creating systematics histo
  TH1F* sys = new TH1F("sys", "sys", r2B_tt[0], v_tt.at(0));
  for(int j = 1; j <= r2B_tt[0]; j++){
    double s_c = fabs(Pred_1Mu->GetBinContent(j)-data->GetBinContent(j))/Pred_1Mu->GetBinContent(j);
    sys->SetBinContent(j, s_c);
  }

  TH1F* Pred_1Mu_SYS = new TH1F(*Pred_1Mu);
  for(int j = 1; j <= r2B_tt[0]; j++){
    double err = sqrt( pow(sys->GetBinContent(j)*Pred_1Mu->GetBinContent(j), 2) + pow(Pred_1Mu->GetBinError(j),2) );
    Pred_1Mu_SYS->SetBinError(j,err);
  }
  
  TFile* fout = new TFile("FinalROOTFiles/Closure_CP.root", "RECREATE");
  
  RatioPlotsBandV2( data, Pred_1Mu, "Data  1#mu", "BKg Pred 1#mu", "PredPlotsFinal/Closure_CP_1mu1Tb_Sep", "RSQ", r2B_tt[0], v_tt.at(0),1);
  RatioPlotsBandV2( data, Pred_1Mu_SYS, "Data  1#mu", "BKg Pred 1#mu", "PredPlotsFinal/Closure_CP_1mu1Tb_SYS_Sep", "RSQ", r2B_tt[0], v_tt.at(0),1);
  RatioPlotsBandV2( data, MC_1Mu1Tb, "Data  1#mu", "MC 1#mu", "PredPlotsFinal/MC_CP_1mu1Tb_Sep", "RSQ", r2B_tt[0], v_tt.at(0),1);
   RatioPlotsBandV2( data_2mu_tt, MC_2MuTT, "Data  2#mu", "BKg Pred 2#mu", "PredPlotsFinal/MC_CP_2Mu1TbTT_Sep", "RSQ", r2B_tt[0], v_tt.at(0),1);
   RatioPlotsBandV2( data_2mu, MC_2MuZ, "Data  2#mu", "BKg Pred 2#mu", "PredPlotsFinal/MC_CP_2Mu1LbZ_Sep", "RSQ", r2B_tt[0], v_tt.at(0),1);
  
   
  W_extra_1mu->Write("w_extra_1mu");
  DY_extra_1mu->Write("dy_extra_1mu");
  Z_extra_1mu->Write("z_extra_1mu");
  P_TT_1Mu->Write("tt_extra_1mu");
  Pred_1Mu_SYS->Write("Pred_1Mu_SYS");
  Pred_1Mu->Write("Pred_1Mu");
  w->Write("w");
  dy->Write("dy");
  z->Write("z");
  tt->Write("tt");
  data->Write("data1mu");
  MC_1Mu1Tb->Write("totalmc1mu");
  
  
  sys->Write();
  fout->Close();
  return 0;

}
