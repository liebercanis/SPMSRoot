#include "TSpms.hxx"
ClassImp(TSpms)

    TSpms::TSpms() : TNamed("TSpms", "TSpms")
{
  fnvalues = 1;
  // 1500;
  fnflat = 1;
  clear();
  cout << "TSpms instance " << endl;
}