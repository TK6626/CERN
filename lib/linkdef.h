#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"

// rootcling -f lib/dict.cxx -I$(root-config --incdir) lib/linkdef.h
// g++ -shared -fPIC lib/dict.cxx `root-config --cflags --libs` -o lib/libdict.so

#ifdef __CLING__
#pragma link C++ class ROOT::VecOps::RVec<TLorentzVector>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector>+;
#endif

