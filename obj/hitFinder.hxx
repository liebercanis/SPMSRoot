////////////////////////////////////////////////////////
//  M.Gold June 2022
// class to make hits from vector data
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex> //includes std::pair, std::make_pair
#include <valarray>
//
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm> // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"
//
#include "TBWave.hxx"
#include "TSpms.hxx"
#include "TBEvent.hxx"

typedef std::vector<std::pair<unsigned, unsigned>> peakType;
typedef std::vector<std::pair<unsigned, unsigned>>::iterator peakTypeIter;
typedef std::map<Double_t, TDetHit, std::less<Double_t>> hitMap;
typedef std::map<Double_t, TDetHit, std::less<Double_t>>::iterator hitMapIter;


class hitFinder
{
public:
  enum
  {
    UPCROSS,
    DOWNCROSS,
    DOUBLEUPCROSS,
    DOUBLEDOWNCROSS
  };
  enum
  {
    NDET = 34
  };

  Double_t qnorm = 1.0;
  double vsign[NDET];
  TFile *fout;
  TTree *ftree;
  TBEvent *bevent;
  TSpms *fspms;
  TString tag;
  int nsamples;
  double rmsCut; 
  TDirectory *fftDir;
  TDirectory *evDir;
  hitFinder(TFile *theFile, TSpms *tspms, TTree *btree, TBEvent *beventInstance, TString theTag);
  virtual ~hitFinder() { ; }
  void plotWave(int idet, Long64_t jentry);
  void plotEvent( unsigned idet,Long64_t ievent);

  void findHits(int idet, Long64_t jentry);
  void differentiate(unsigned diffStep);
  vector<TString> detNames;
  std::vector<double> rdigi;
  std::vector<double> digi;
  std::vector<double> ddigi;
  std::vector<double> hdigi;
  std::vector<double> fdigi; // filtered
  hitMap detHits;
  peakType peakList;
  std::vector<Int_t> peakKind;
  double timeUnit;
  double microSec;
  void fevent( Long64_t ievent, vector<double> rdigi);
  void derivativePeaks(Int_t idet, Double_t rms);
  hitMap makeHits(int idet, Double_t &triggerTime, Double_t &firstCharge);
  void trimPeaks(int idet, std::vector<Double_t> v);
  void extendPeaks(int idet, std::vector<Double_t> v);
  bool getTransforms();
  // hist
  TH1D *hPeakCount;
  TH1I *hHitLength;
  TH1I *hPeakNWidth;

  /// The fft class to take the fourier transform.
  TVirtualFFT *fFFT;
  /// The fft class to take the inverse fourier transform.
  TVirtualFFT *fInverseFFT;

  std::vector<std::complex<double>> FFT(Int_t idet, std::vector<double> rdigi, bool first = true);
  std::vector<Double_t> inverseFFT(Int_t idet, std::vector<std::complex<double>> VectorComplex, std::vector<double> rdigi);
  TGraph *gTransform[NDET];
  bool gotTransforms;

  TH1D *hFFT[NDET];
  TH1D *hInvFFT[NDET];
  TH1D *hFFTFilt[NDET];
  TH1D *hEvWave[NDET];
  TH1D *hEvHitWave[NDET];
  TH1D *hEvDerWave[NDET];
  TH1D *hEvFiltWave[NDET];
  TH1D *hCutHigh;
  TH1D *hCutLow;
};
