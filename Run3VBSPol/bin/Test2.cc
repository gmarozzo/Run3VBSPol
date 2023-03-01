#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Auxiliary.cpp"
#include "Histodef.cpp"

// include user defined histograms and auxiliary macros
using namespace std;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024

// function to calculate the weight for each event
// the weight is calculated as the product of luminosity and cross section of the process times the genWeight,
// LATER TO BE divided by the number of generated events OF ALL FILES OF THE DATASET(S)

double getWeight(double luminosity, double crossSection, Float_t genWeight, double SumWeights)
{
  return (luminosity * crossSection * genWeight); // / SumWeights;
}


void Test2(string inputFile, string ofile, double crossSection = -1, double IntLuminosity = 59.827879506, bool Signal = false)
{

  if (crossSection < 0. || IntLuminosity < 0.)
    {
      std::cout << "WARNING: crossection " << crossSection << " and Integrated luminosity " << IntLuminosity << endl;
    }

  cout<<"Call completed!"<<endl;

  TFile *fin = TFile::Open(inputFile.c_str());
  TTree *trun = static_cast<TTree *>(fin->Get("Runs"));
  Long64_t genEventCount;
  Double_t genEventSumw;
  trun->SetBranchStatus("*", 0);
  trun->SetBranchStatus("genEventSumw", 1);
  trun->SetBranchStatus("genEventCount", 1);
  trun->SetBranchAddress("genEventSumw", &genEventSumw);
  trun->SetBranchAddress("genEventCount", &genEventCount);


  trun->GetEntry(0);

  TTree *tin = static_cast<TTree *>(fin->Get("Events"));

  // Set all branches to 0
  tin->SetBranchStatus("*", 0);
  // get the pt
  Float_t Muon_pt[MAX_ARRAY_SIZE],  Muon_mass[MAX_ARRAY_SIZE];
  Float_t Muon_eta[MAX_ARRAY_SIZE],  Muon_phi[MAX_ARRAY_SIZE];
  UInt_t nMuon;

  tin->SetBranchStatus("Muon_pt", 1);
  tin->SetBranchAddress("Muon_pt", &Muon_pt);
  tin->SetBranchStatus("nMuon", 1);
  tin->SetBranchAddress("nMuon", &nMuon);  
  tin->SetBranchStatus("Muon_eta", 1);
  tin->SetBranchAddress("Muon_eta", &Muon_eta);
  tin->SetBranchStatus("Muon_phi", 1);
  tin->SetBranchAddress("Muon_phi", &Muon_phi);
  tin->SetBranchStatus("Muon_mass", 1);
  tin->SetBranchAddress("Muon_mass", &Muon_mass);

  // get gen quantities

  Float_t GenPart_pt[GEN_MAX_ARRAY_SIZE];

  tin->SetBranchStatus("GenPart_pt",1);
  tin->SetBranchAddress("GenPart_pt",&GenPart_pt);

  // collect the trigger information
  Bool_t HLT_IsoMu24;
  tin->SetBranchStatus("HLT_IsoMu24", 1);
  tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

  // collect the triggger Ids
  Int_t Muon_charge[MAX_ARRAY_SIZE],Muon_nTrackerLayers[MAX_ARRAY_SIZE], Muon_genPartIdx[MAX_ARRAY_SIZE];
  Bool_t  Muon_triggerIdLoose[MAX_ARRAY_SIZE], Muon_tightId[MAX_ARRAY_SIZE];
  Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE];
  tin->SetBranchStatus("Muon_tightId", 1);
  tin->SetBranchStatus("Muon_charge", 1);
  tin->SetBranchStatus("Muon_triggerIdLoose", 1);
  tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
  tin->SetBranchStatus("Muon_nTrackerLayers", 1);
  tin->SetBranchStatus("Muon_genPartIdx", 1);
  tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
  tin->SetBranchAddress("Muon_charge", &Muon_charge);
  tin->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose);
  tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
  tin->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);
  tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);

  // gen weight
  Float_t genWeight;
  tin->SetBranchStatus("genWeight", 1);
  tin->SetBranchAddress("genWeight", &genWeight);

  int non_matching_muon = 0, non_matching_electron = 0;
  int n_dropped = 0;
  int trigger_dropped = 0;
  UInt_t nEv = tin->GetEntries();
  unsigned int n_events = nEv;
  TLorentzVector *Muon1_p4 = new TLorentzVector();
  TLorentzVector *Muon2_p4 = new TLorentzVector();
  float Weight;
      
  // save the histograms in a new File

  TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
  // create a new tree for the output
  TTree *tout = new TTree("tout", "tout");
  TTree *trun_out = new TTree("Run_out", "Run_out");
  Float_t invMass, muon1_eta, muon1_pt, muon2_pt, muon2_eta, dphi;

  tout->Branch("dphi", &dphi);
  tout->Branch("invMass", &invMass);
  tout->Branch("muon2_eta", &muon2_eta);
  tout->Branch("muon2_pt", &muon2_pt);
  tout->Branch("muon1_eta", &muon1_eta);
  tout->Branch("muon1_pt", &muon1_pt);
  tout->Branch("Weight", &Weight);

  trun_out->Branch("genEventSumw", &genEventSumw);
  trun_out->Branch("IntLumi", &IntLuminosity);
  trun_out->Branch("xs", &crossSection);
  trun_out->Branch("nEv", &n_events);

  trun_out->Fill(); // we already called trun->GetEntry(0);


  bool From2Taus=false, FromTau=false;
  TFile *foutT = new TFile(("Tau"+ofile).c_str(), "RECREATE");
  TTree *toutT = new TTree("toutT", "toutT");
  TTree *trun_outT = new TTree("Run_outT", "Run_outT");

  TRandom3* RndGen= new TRandom3();
  fout->cd();
    #pragma omp parallel for
  for (UInt_t i = 0; i <nEv; i++)
    {
      tin->GetEntry(i);
      if (i % 100000 == 0)
	std::cout << "Processing entry " << i << " of " << nEv << endl;
      // apply triggers

      if (!(HLT_IsoMu24)){
	trigger_dropped++;
	continue;
      };

      bool gotmuplus=false,gotmuminus=false;
      int mu1idx=-1, mu2idx=-1;
      for (UInt_t j = 0; j < nMuon; j++){
	if (abs(Muon_eta[j])<2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15){
	  int NMCparticle=Muon_genPartIdx[j];
	  if(Muon_pt[j]>27.|| ((gotmuplus||gotmuminus) && Muon_pt[j]>25.)){
	    if (!gotmuplus && Muon_charge[j]==1){
	      Muon1_p4->SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
	      gotmuplus=true;
	      mu1idx=j;
	    }
	    if (!gotmuminus && Muon_charge[j]==-1){
	      Muon2_p4->SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
	      gotmuminus=true;
	      mu2idx=j;
	    }
	  }
	}
      }
 
      if(!(gotmuplus && gotmuminus)) {continue;}
      if(Muon1_p4->DeltaR(*Muon2_p4)<0.4) {continue;}


      Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);

      vector<bool> tagged; 
      double t_weight=1.;
      
      Weight*=t_weight; 
           

      dphi=Muon1_p4->DeltaPhi(*Muon2_p4);
       
      muon1_pt = Muon1_p4->Pt();
      muon1_eta = Muon1_p4->Eta();
      muon2_pt = Muon2_p4->Pt();
      muon2_eta = Muon2_p4->Eta();

      h_Muon1_pt->Fill(muon1_pt,Weight);
      h_Muon1_eta->Fill(muon1_eta,Weight);
      h_Muon2_pt->Fill(muon2_pt,Weight);
      h_Muon2_eta->Fill(muon2_eta,Weight);
      h_acopla_mumu->Fill(M_PI-dphi,Weight);

      invMass = (*(Muon1_p4) + *(Muon2_p4)).M();
      h_Muon_Muon_invariant_mass->Fill(invMass,Weight);
      
      tout->Fill();
    }


  std::cout << "NeV = " << nEv << endl;
  std::cout << "trigger dropped = " << trigger_dropped << endl;
  std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

  std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / nEv) << endl;
  int Rem_trigger=nEv-trigger_dropped; //remember the cross trigger in Data
  std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
  std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;


  fout->cd();
  tout->Write();
  trun_out->Write();
  // Write the histograms to the file
  h_Muon1_eta->Write();
  h_Muon1_pt->Write();
  h_Muon2_eta->Write();
  h_Muon2_pt->Write();
    
  h_Muon_Muon_invariant_mass->Write();
  h_acopla_mumu->Write();
  h_NJets->Write();

  fout->Close();
}

int main(int argc, char **argv)
{

  string inputFile = argv[1];
  string outputFile = argv[2];
  double crossSection = atof(argv[3]);
  double IntLuminosity = atof(argv[4]);
  string boolstr = argv[5];
  bool Signal = (boolstr == "true");

  h_Muon1_pt->Sumw2();
  h_Muon1_eta->Sumw2();
  h_Muon2_pt->Sumw2();
  h_Muon2_eta->Sumw2();
  h_Muon_Muon_invariant_mass->Sumw2();    
  h_leading_lepton_pt->Sumw2();
  h_NJets->Sumw2();
  h_acopla_mumu->Sumw2();

  Test2(inputFile, outputFile, crossSection, IntLuminosity, Signal);

  return 0;
}
