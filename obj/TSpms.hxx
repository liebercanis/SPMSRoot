/**
** MG, May 2022
**/
#ifndef TSPMS_DEFINED
#define TSPMS_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <TChain.h>
#include <TLeaf.h>
#include <vector>

using namespace std;

// class to store info for the event

class TSpms : public TNamed
{
public:
  TSpms();
  // TSpms::~TSpms(){}
  TTree *fChain;
  Int_t fCurrent; //! current Tree number in a TChain
  // get this from file
  int fnvalues, fnflat;

  // data elements
  UShort_t baseline;
  UInt_t channel;
  UShort_t energy;
  Int_t ievt;
  Int_t numtraces;
  UInt_t packet_id;
  Double_t timestamp;
  UInt_t cumulative_length;
  vector<Short_t> flattened_data;
  Double_t dt;
  Double_t t0;
  // UShort_t  values[fnvalues];
  std::vector<UShort_t> values;
  UShort_t wf_max;
  Float_t wf_std;

  // List of branches
  TBranch *b_baseline; //!
  TBranch *b_channel;  //!

  TBranch *b_energy;    //!
  TBranch *b_ievt;      //!
  TBranch *b_numtraces; //!
  TBranch *b_packet_id; //!
  TBranch *b_timestamp; //!

  // List of branches
  TBranch *b_cumulative_length; //!
  TBranch *b_flattened_data;    //!
  // List of branches
  TBranch *b_dt;     //!
  TBranch *b_t0;     //!
  TBranch *b_values; //!
  TBranch *b_wf_max; //!
  TBranch *b_wf_std; //!

  int getNValues() { return fnvalues; }
  int getNFlat() { return fnflat; }

  void setArrays(int nv, int nf)
  {
    fnvalues = nv;
    fnflat = nf;
    printf(" setArrays fnvalues =%i fnflat =%i \n", fnvalues, fnflat);
  }
  void clear()
  {
    baseline = 0;
    channel = -1;

    energy = 0;
    ievt = 0;
    numtraces = 0;
    packet_id = 0;
    timestamp = 0;
    wf_max = 0;
    wf_std = 0;
    cumulative_length = 0;
    dt = 0;
    t0 = 0;
    values.clear();
    values.resize(fnvalues);
    flattened_data.clear();
    flattened_data.resize(fnflat);
    // cout << " TSpms clear vnvalues = " << values.size() << " fnflat = " << flattened_data.size() << endl;
  }

  void init(TTree *tree)
  {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree)
      return;
    fChain = tree;
    fCurrent = -1;
    // fChain->SetMakeClass(1);

    fChain->SetBranchAddress("baseline", &baseline, &b_baseline);
    fChain->SetBranchAddress("channel", &channel, &b_channel);
    fChain->SetBranchAddress("energy", &energy, &b_energy);
    fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
    fChain->SetBranchAddress("numtraces", &numtraces, &b_numtraces);
    fChain->SetBranchAddress("packet_id", &packet_id, &b_packet_id);
    fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
    fChain->SetBranchAddress("wf_max", &wf_max, &b_wf_max);
    fChain->SetBranchAddress("wf_std", &wf_std, &b_wf_std);

    fChain->SetBranchAddress("cumulative_length", &cumulative_length, &b_cumulative_length);
    fChain->SetBranchAddress("flattened_data", &flattened_data[0], &b_flattened_data);

    fChain->SetBranchAddress("dt", &dt, &b_dt);
    fChain->SetBranchAddress("t0", &t0, &b_t0);
    // fChain->FindLeaf("values")->SetAddress(values);
    fChain->SetBranchAddress("values", &values[0], &b_values);
  }

  void create(TTree *tree)
  {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree)
      return;
    fChain = tree;
    fCurrent = -1;
    // fChain->SetMakeClass(1);

    fChain->Branch("baseline", &baseline);
    fChain->Branch("channel", &channel);
    fChain->Branch("energy", &energy);
    fChain->Branch("ievt", &ievt);
    fChain->Branch("numtraces", &numtraces);
    fChain->Branch("packet_id", &packet_id);
    fChain->Branch("timestamp", &timestamp);
    fChain->Branch("wf_max", &wf_max);
    fChain->Branch("wf_std", &wf_std);

    fChain->Branch("cumulative_length", &cumulative_length);
    fChain->Branch("flattened_data", &flattened_data[0]);

    fChain->Branch("dt", &dt);
    fChain->Branch("t0", &t0);
    // fChain->FindLeaf("values")->SetAddress(values);
    fChain->Branch("values", &values[0]);
  }

  ClassDef(TSpms, 1)
};
#endif
