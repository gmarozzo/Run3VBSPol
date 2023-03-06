#include "Run3VBSPol/Run3VBSPol/plugins/TreePlanter.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "Math/LorentzVector.h"

#include "TTree.h"
#include "TVectorD.h"
#include "TLorentzVector.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>

#include <typeinfo>


using namespace boost::assign;

using std::cout;
using std::endl;

class Particle;
class Jet;
class Lepton;
class Electron;


TreePlanter::TreePlanter(const edm::ParameterSet &config)
  : theMuonToken     (consumes<pat::MuonCollection>                (config.getParameter<edm::InputTag>("muons"    )))
  , theElectronToken (consumes<pat::ElectronCollection>            (config.getParameter<edm::InputTag>("electrons")))
  , theJetToken      (consumes<std::vector<pat::Jet> >             (config.getParameter<edm::InputTag>("jets"     )))
  , theHLTToken      (consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))){
}


void TreePlanter::beginJob(){

  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("ntp1","ntp1");

  muon_pt_ = new std::vector<float>;
  muon_eta_ = new std::vector<float>;
  muon_recSF_ = new std::vector<float>;
  muon_idSF_ = new std::vector<float>;

  electron_pt_ = new std::vector<float>;
  electron_sceta_ = new std::vector<float>;
  electron_SF_ = new std::vector<float>;

  njets_ = new std::vector<int>;
  mll_ = new std::vector<float>;

  hlt_= new std::vector<string>;

  tree_->Branch("muon_pt",&muon_pt_);
  tree_->Branch("muon_eta",&muon_eta_);
  tree_->Branch("muon_recSF",&muon_recSF_);
  tree_->Branch("muon_idSF",&muon_idSF_);
  tree_->Branch("electron_pt",&electron_pt_);
  tree_->Branch("electron_sceta",&electron_sceta_);
  tree_->Branch("electron_SF",&electron_SF_);
  tree_->Branch("njets",&njets_);
  tree_->Branch("mll",&mll_);
  tree_->Branch("hlt",&hlt_);
}


void TreePlanter::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
  
}


void TreePlanter::endRun(const edm::Run& run, const edm::EventSetup& setup){
  
}

void TreePlanter::endJob(){

  delete muon_pt_;
  delete muon_eta_;
  delete muon_recSF_;
  delete muon_idSF_;
  delete electron_pt_;
  delete electron_sceta_;
  delete electron_SF_;
  delete njets_;
  delete mll_;

}


void TreePlanter::initTree(){
  
}



void TreePlanter::analyze(const edm::Event& event, const edm::EventSetup& setup){

  initTree();

  // Load a bunch of objects from the event
  edm::Handle<pat::MuonCollection>       muons          ; event.getByToken(theMuonToken    ,     muons);
  edm::Handle<pat::ElectronCollection>   electrons      ; event.getByToken(theElectronToken, electrons);
  edm::Handle<std::vector<pat::Jet> >    jets           ; event.getByToken(theJetToken     ,      jets);


  unsigned int muonSize= muons->size();
  unsigned int electronSize = electrons->size();
 
  float muonrecptlimits[]={35.0,50.0,75.0,100.0,500.0};
  float muonrecetalimits[]={0.0,1.2,2.4};
  float muonidptlimits[]={35.0,40.0,45.0,55.0,70.0,100.0,500.0};
  float muonidetalimits[]={0.0,0.9,1.2,2.1,2.4};
  
  float electronidptlimits[]={35.0,40.0,45.0,55.0,70.0,100.0,500.0};
  float electronidetalimits[]={0.0,0.5,1.0,1.444,1.566,2.0,2.4};

  float muonrecSFmap[4][2]={{0.999,0.999},{0.993,0.999},{1.005,1.003},{1.008,1.007}};
  float muonidSFmap[6][4]={{0.982,0.968,0.996,0.999},{0.985,0.977,0.994,0.993},{0.990,0.979,0.994,0.992},{0.989,0.978,0.991,0.989},{1.000,0.987,0.989,1.010},{1.005,1.009,1.007,0.996}};
  float electronidSFmap[6][6]={{0.961,0.942,0.960,1.000,0.991,0.974},{0.967,0.948,0.957,1.000,0.991,0.947},{0.970,0.955,0.964,1.000,0.986,0.939},{0.958,0.940,0.966,1.000,0.983,0.909},{0.959,0.990,0.986,1.000,1.007,0.911},{1.017,0.983,1.027,1.000,0.971,0.921}};


  for (unsigned int imuon=0; imuon<muonSize; imuon++){
    pat::Muon muontofill = (*muons)[imuon];
    muon_pt_->push_back(muontofill.p4().Pt());
    muon_eta_->push_back(muontofill.p4().Eta());
    for(int i=0; i<4;i++){
      if(muontofill.p4().Pt()<muonrecptlimits[i+1]){
        for(int j=0; j<2;j++){
          if(muontofill.p4().Eta()<muonrecetalimits[i+1]){
            muon_recSF_->push_back(muonrecSFmap[i][j]);
            break;}
	}
	break;
      }	   
    }
    
    for(int i=0; i<6;i++){
      if(muontofill.p4().Pt()<muonidptlimits[i+1]){
        for(int j=0; j<4;j++){
          if(muontofill.p4().Eta()<muonidetalimits[i+1]){
            muon_idSF_->push_back(muonidSFmap[i][j]);
            break;}
	}
        break;
      }
    }
  }
  

  for (unsigned int ielectron=0; ielectron<electronSize; ielectron++){
    pat::Electron electrontofill = (*electrons)[ielectron];
    electron_pt_->push_back(electrontofill.p4().Pt());
    electron_sceta_->push_back(electrontofill.p4().Eta());
    for(int i=0; i<6;i++){
      if(electrontofill.p4().Pt()<electronidptlimits[i+1]){
        for(int j=0; j<6;j++){
          if(electrontofill.p4().Eta()<electronidetalimits[i+1]){
            electron_SF_->push_back(electronidSFmap[i][j]);
            break;}
	}
	break;
      }
    }
  }
  
  pat::Muon leadingmuon;
  pat::Electron leadingelectron;

  if(muonSize>0){leadingmuon = (*muons)[0];}
  if(electronSize>0){leadingelectron = (*electrons)[0];}

  unsigned int jetSize = jets->size();

  for(unsigned int k=0; k<jets->size();k++){
    float dR1 = reco::deltaR((*jets)[k],leadingmuon);
    float dR2 = reco::deltaR((*jets)[k],leadingelectron);
    if(dR1 < 0.4  || dR2 < 0.4){
      jetSize--;
    }    
  }

  njets_->push_back(jetSize);

  if(muonSize*electronSize==1 && jetSize>1 ){
    if(leadingmuon.charge()*leadingelectron.charge() == -1){
      ROOT::Math::LorentzVector p4 = leadingmuon.p4()+leadingelectron.p4();
      float mll = p4.M();
      if(mll>20) {mll_->push_back(mll);}
    }   
  }

  tree_->Fill();
  (*muon_pt_).clear();
  (*muon_eta_).clear();
  (*muon_recSF_).clear();
  (*muon_idSF_).clear();
  (*electron_pt_).clear();
  (*electron_sceta_).clear();
  (*electron_SF_).clear();
  (*njets_).clear();
  (*mll_).clear();
}


#include "FWCore/Framework/interface/MakerMacros.h"
// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(TreePlanter);
