#include "Run3VBSPol/DataFormats/src/Particle.C"
#include "Run3VBSPol/DataFormats/src/Lepton.C"
#include "Run3VBSPol/DataFormats/src/Jet.C"
#include "Run3VBSPol/DataFormats/src/Electron.C"


#ifdef __CINT__

#pragma link C++ class  phys::Jet::JetScores+;
#pragma link C++ class  phys::Jet::DeepFlavourScores+;

#pragma link C++ class  phys::Particle+;
#pragma link C++ class  phys::Lepton+;
#pragma link C++ class  phys::Jet+;
#pragma link C++ class  phys::Electron+;

#pragma link C++ class  std::vector<phys::Particle>;
#pragma link C++ class  std::vector<phys::Lepton>;
#pragma link C++ class  std::vector<phys::Jet>;
#pragma link C++ class  std::vector<phys::Electron>;

#endif
