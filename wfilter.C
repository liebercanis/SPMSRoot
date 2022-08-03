//  fit for deconvolution function
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TKey.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace TMath;

enum
{
  NDET = 34
};
TH1D *hWave[NDET];
TH1D *hFFT[NDET];
TF1 *fitSig[NDET];
TF1 *fitNoise[NDET];
TFile *fout;
int nfreqencies = 0;
double fnotch = double(2 * nfreqencies) / 291.0;
double sigNotch = 17.5; // 17.5; // was 50

static Double_t notch(Double_t A, Double_t x)
{
  double mean = fnotch; // nfreqencies/fnotch;
  double sig = sigNotch;

  double f = A * Gaus(x, mean, sig);
  // if(x<mean) f=A*1.0E4;
  return f;
}

static Double_t fpmtSig(Double_t *xx, Double_t *par)
{
  double x = xx[0];
  double A = 10.44;
  double a = -0.002841;
  return Exp(A + a * x) + par[0] + par[1] * x + par[2] * x * x;
}

bool getHistos()
{

  TFile *fin = new TFile("20220413-1000.root", "readonly");
  if (!fin)
    return false;
  if (!fin->IsOpen())
    return false;
  TDirectory *fftDir;
  fin->GetObject("fftDir", fftDir);
  printf(" geting histos from file %s dir %s \n", fin->GetName(), fftDir->GetName());
  // fftDir->ls();

  TList *list = fftDir->GetListOfKeys();
  if (!list)
  {
    printf("<E> No keys found in file\n");
    return false;
  }
  TIter next(list);
  TKey *key;
  TObject *obj;
  printf(" keys found in fftDir\n");
  int ihist = 0;
  while ((key = (TKey *)next()))
  {
    if(ihist==NDET)
      break;
    obj = key->ReadObj();
    if (!obj->InheritsFrom("TH1D"))
      continue;
    TString hname(obj->GetName());
    if (hname.Contains("FFTD") && !(hname.Contains("Inv")))
    {
      hFFT[ihist] = (TH1D *)obj;
      //cout << hname << endl;
      ++ihist;
    }
  }

  nfreqencies = hFFT[2]->GetNbinsX();
  printf("got %i %i \n", ihist, nfreqencies);

  return true;
}

void fitHist(int i)
{
  double x0 = 120;
  if (i == 1)
    x0 = 1200;
  fitSig[i] = new TF1(Form("sig%i", i), "Exp([0]+x*[1])", 0, x0);
  /*
     Minimizer is Linear
     Chi2                      =  7.74208e+08
     NDf                       =         4091
     p0                        =  1.76659e+08   +/-   1204.68
     p1                        =      -525305   +/-   6.94655
     p2                        =      699.556   +/-   0.0147496
     p3                        =    -0.519127   +/-   1.51237e-05
     p4                        =    0.0002276   +/-   8.35068e-09
     p5                        = -5.83902e-08   +/-   2.5422e-12
     p6                        =   8.0835e-12   +/-   4.01737e-16
     p7                        =  -4.6552e-16   +/-   2.57313e-20

     p0                        =  1.79086e+08   +/-   1373.93
     p1                        =      -516391   +/-   9.19407
     p2                        =      563.383   +/-   0.0181763
     p3                        =      -0.2173   +/-   1.07035e-05

*/
  if (i == 1)
  {
    double fitHigh[2] = {4.11436e+06, -779.138};
    // double fitLow[4]  = { 1.63757e+08  ,-396577  ,338.182   ,-0.0971809 };
    double fitLow[8] = {1.76659e+08, -525305, 699.556, -0.519127, 0.0002276, -5.83902e-08, 8.0835e-12, -4.6552e-16};
    // double fitLPar[4] = {1.79086e+08 ,  -516391 ,  563.383 , -0.2173 };
    fitSig[i] = new TF1(Form("sig%i", i), "pol7", 0, nfreqencies);
    fitSig[i]->SetParameters(fitLow);
    fitNoise[i] = new TF1(Form("noise%i", i), "pol1", 0, nfreqencies);
    fitNoise[i]->SetParameters(fitHigh);
  }
  else
  {
    fitSig[i] = new TF1(Form("sig%i", i), "Exp([0]+x*[1])", 0, nfreqencies);
    fitSig[i]->SetParNames("const", "slope");
    fitNoise[i] = new TF1(Form("noise%i", i), "[0]+x*[1]", 0, nfreqencies);
    fitSig[i]->SetParameters(10.44, -0.002841);
    fitNoise[i]->SetParNames("const", "slope");
    fitNoise[i]->SetParameters(429.9, -.069);
  }
  fitSig[i]->SetLineColor(kBlue);
  fitSig[i]->GetXaxis()->SetTitle(" frequency ");
  fitSig[i]->GetYaxis()->SetTitle(" power ");
  fitNoise[i]->SetLineColor(kBlack);
  fitNoise[i]->GetXaxis()->SetTitle(" frequency ");
  fitNoise[i]->GetYaxis()->SetTitle(" power ");

  TString hname;
  hname.Form("FFTLowDet%i", i);
  TH1D *histLow = (TH1D *)hFFT[i]->Clone(hname);
  hname.Form("FFTHighDet%i", i);
  TH1D *histHigh = (TH1D *)hFFT[i]->Clone(hname);

  if (i == 1)
    histLow->Fit(fitSig[i], "0", "", 0, nfreqencies);
  else
    histLow->Fit(fitSig[i], "0", "", 0, x0);
  histHigh->Fit(fitNoise[i], "0", "", x0, nfreqencies);

  TString canname;

  canname.Form("FitToDet%i", i);
  gStyle->SetOptFit(1111);
  gStyle->SetOptLogy();
  TCanvas *can = new TCanvas(canname, canname);
  can->Divide(1, 2);
  can->cd(1);
  histLow->Draw("");
  fitSig[i]->Draw("same");
  can->cd(2);
  histHigh->Draw("");
  fitNoise[i]->Draw("same");

  can->Print(".pdf");

  fout->Add(fitSig[i]);
  fout->Add(fitNoise[i]);
}

TGraph *makeWiener(int i)
{
  std::vector<double> yval;
  std::vector<double> xval;
  for (int j = 0; j < nfreqencies; ++j)
  {
    double x = double(j);
    xval.push_back(x);
    double yn, ys;
    yn = fitNoise[i]->Eval(x);
    double fsig = fitSig[i]->Eval(x);
    ys = fsig - yn;
    // boost for sipms.
    if (ys < 0)
      ys = 0;
    double w = ys / (ys + yn);
    if (w < 0)
      w = 0.0;
    if (w > 1)
      w = 1.0;
    yval.push_back(w);
    printf(" bin %i val %E \n",j,ys/(ys+yn));
  }

  TGraph *gr = new TGraph(xval.size(), &xval[0], &yval[0]);
  gr->SetName(Form("WienerTransDet%i", i));
  gr->SetTitle(Form("WienerTransDet%i", i));
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(0.2);
  gr->GetHistogram()->GetXaxis()->SetTitle(" frequency ");

  TString canName;
  canName.Form("WienerTransDet%i", i);
  TCanvas *wcan = new TCanvas(canName, canName);
  wcan->Divide(1, 2);
  wcan->cd(1);
  fitSig[i]->Draw("");
  fitNoise[i]->Draw("sames");
  wcan->cd(2);
  gr->Draw("");
  wcan->Print(".pdf");

  return gr;
}

void wfilter()
{
  fout = new TFile("WienerTransforms.root", "recreate");
  if (!getHistos())
  {
    printf("aint got hist det 2\n");
    return;
  }

  printf("got hist det 2\n");

  /*
  TF1* fpmt = new TF1(Form("sig%i",1),fpmtSig,0,50,3);

  fpmt->SetParameters(10.E8,-525305,-0.519127 );

  TString hname;
  hname.Form("FFTLowDet%i",1);
  TH1D *histLow = (TH1D*) hFFT[1]->Clone(hname);
  histLow->Fit(fpmt,"R");


  TString canname;
  canname.Form("FitToDet%i",1);
  gStyle->SetOptFit(1111);
  gStyle->SetOptLogy();
  TCanvas *can = new TCanvas(canname,canname);
  histLow->Draw("");
  fpmt->Draw("same");

  return;
  */

  for (int i = 2; i < 3; ++i)
  {
    fout->Add(hWave[i]);
    fout->Add(hFFT[i]);
    fitHist(i);
  }

  for (int i = 2; i < 3; ++i)
  {
    TGraph *gr = makeWiener(i);
    fout->Add(gr);
  }

  fout->ls();
  fout->Write();
}
