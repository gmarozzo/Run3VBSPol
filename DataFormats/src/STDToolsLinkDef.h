#include "Run3VBSPol/DataFormats/interface/Particle.h"
#include "Run3VBSPol/DataFormats/interface/Lepton.h"
#include "Run3VBSPol/DataFormats/interface/Jet.h"
#include "Run3VBSPol/DataFormats/interface/Electron.h"
#include "Run3VBSPol/DataFormats/interface/Photon.h"
#include "Run3VBSPol/DataFormats/interface/Boson.h"
#include "Run3VBSPol/DataFormats/interface/DiBoson.h"
#include "Run3VBSPol/DataFormats/interface/GenEventWeights.h"
#include "Run3VBSPol/DataFormats/interface/MELA.h"
#include "Run3VBSPol/DataFormats/interface/RegionsCounter.h"
#include "Run3VBSPol/DataFormats/interface/RegionTypes.h"
#include "Run3VBSPol/DataFormats/interface/Proton.h"
#include "Run3VBSPol/DataFormats/interface/ProtonPair.h"

#pragma link C++ class  phys::GenEventWeights+;
#pragma link C++ class  phys::MELA+;
#pragma link C++ class  phys::RegionTypes+;
#pragma link C++ class  phys::RegionsCounter+;
#pragma link C++ class  std::map<phys::RegionTypes,Int_t>+;
#pragma link C++ class  phys::Jet::JetScores+;

#pragma link C++ class  phys::Proton+;
#pragma link C++ class  phys::ProtonPair+;
#pragma link C++ class  phys::Particle+;
#pragma link C++ class  phys::Lepton+;
#pragma link C++ class  phys::Jet+;
#pragma link C++ class  phys::Electron+;
#pragma link C++ class  phys::Photon+;
#pragma link C++ class  phys::Boson<phys::Particle>+;
#pragma link C++ class  phys::Boson<phys::Lepton>+;
#pragma link C++ class  phys::Boson<phys::Electron>+;
#pragma link C++ class  phys::Boson<phys::Jet>+;
#pragma link C++ class  phys::Boson<phys::Photon>+;
#pragma link C++ class  phys::DiBoson<phys::Particle, phys::Particle >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Electron >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Electron >+;

#pragma link C++ class  std::vector<phys::Proton>;
#pragma link C++ class  std::vector<phys::ProtonPair>;
#pragma link C++ class  std::vector<phys::Particle>;
#pragma link C++ class  std::vector<phys::Lepton>;
#pragma link C++ class  std::vector<phys::Jet>;
#pragma link C++ class  std::vector<phys::Electron>;
#pragma link C++ class  std::vector<phys::Photon>;
#pragma link C++ class  std::vector<phys::Boson<phys::Particle> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Lepton> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Electron> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Jet> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Photon> >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Particle, phys::Particle > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Electron > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Electron > >;
#pragma link C++ class  std::pair<phys::Boson<phys::Lepton>, phys::Lepton>+;
#pragma link C++ class  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> >;


#endif