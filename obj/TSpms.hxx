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
  enum
  {
    NVALUES = 1500,
    NFLAT = 34
  };

  // data elements
  UShort_t baseline;
  UInt_t channel;
  Double_t deadtime;
  Int_t dr_maxticks;
  Int_t dr_start_pps;
  Int_t dr_start_ticks;
  Int_t dr_stop_pps;
  Int_t dr_stop_ticks;
  UShort_t energy;
  Int_t ievt;
  Int_t numtraces;
  UInt_t packet_id;
  Double_t runtime;
  Double_t timestamp;
  Int_t to_abs_mu_usec;
  Int_t to_dt_mu_usec;
  Long_t to_master_sec;
  Long_t to_mu_sec;
  Int_t to_mu_usec;
  Long_t to_start_sec;
  Int_t to_start_usec;
  Int_t ts_maxticks;
  Int_t ts_pps;
  Int_t ts_ticks;
  UShort_t wf_max;
  Float_t wf_std;

  UInt_t cumulative_length;
  vector<Short_t> flattened_data;

  Double_t dt;
  Double_t t0;
  // UShort_t  values[NVALUES];
  std::vector<UShort_t> values;

  // List of branches
  TBranch *b_baseline;          //!
  TBranch *b_channel;           //!
  TBranch *b_deadtime;          //!
  TBranch *b_dr_maxticks;       //!
  TBranch *b_dr_start_pps;      //!
  TBranch *b_dr_start_ticks;    //!
  TBranch *b_dr_stop_pps;       //!
  TBranch *b_dr_stop_ticks;     //!
  TBranch *b_energy;            //!
  TBranch *b_ievt;              //!
  TBranch *b_numtraces;         //!
  TBranch *b_packet_id;         //!
  TBranch *b_runtime;           //!
  TBranch *b_timestamp;         //!
  TBranch *b_to_abs_mu_usec;    //!
  TBranch *b_to_dt_mu_usec;     //!
  TBranch *b_to_master_sec;     //!
  TBranch *b_to_mu_sec;         //!
  TBranch *b_to_mu_usec;        //!
  TBranch *b_to_start_sec;      //!
  TBranch *b_to_start_usec;     //!
  TBranch *b_ts_maxticks;       //!
  TBranch *b_ts_pps;            //!
  TBranch *b_ts_ticks;          //!
  TBranch *b_wf_max;            //!
  TBranch *b_wf_std;            //!
                                // List of branches
  TBranch *b_cumulative_length; //!
  TBranch *b_flattened_data;    //!
  // List of branches
  TBranch *b_dt;     //!
  TBranch *b_t0;     //!
  TBranch *b_values; //!

  void clear()
  {
    baseline = 0;
    channel = -1;
    deadtime = 0;
    dr_maxticks = 0;
    dr_start_pps = 0;
    dr_start_ticks = 0;
    dr_stop_pps = 0;
    dr_stop_ticks = 0;
    energy = 0;
    ievt = 0;
    numtraces = 0;
    packet_id = 0;
    runtime = 0;
    timestamp = 0;
    to_abs_mu_usec = 0;
    to_dt_mu_usec = 0;
    to_master_sec = 0;
    to_mu_sec = 0;
    to_mu_usec = 0;
    to_start_sec = 0;
    to_start_usec = 0;
    ts_maxticks = 0;
    ts_pps = 0;
    ts_ticks = 0;
    wf_max = 0;
    wf_std = 0;
    cumulative_length = 0;
    dt = 0;
    t0 = 0;
    // for(int i=0; i<NVALUES ; ++i) values[i]=0;
    values.clear();
    values.resize(NVALUES);
    flattened_data.clear();
    flattened_data.resize(NFLAT);
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
    fChain->SetBranchAddress("deadtime", &deadtime, &b_deadtime);
    fChain->SetBranchAddress("dr_maxticks", &dr_maxticks, &b_dr_maxticks);
    fChain->SetBranchAddress("dr_start_pps", &dr_start_pps, &b_dr_start_pps);
    fChain->SetBranchAddress("dr_start_ticks", &dr_start_ticks, &b_dr_start_ticks);
    fChain->SetBranchAddress("dr_stop_pps", &dr_stop_pps, &b_dr_stop_pps);
    fChain->SetBranchAddress("dr_stop_ticks", &dr_stop_ticks, &b_dr_stop_ticks);
    fChain->SetBranchAddress("energy", &energy, &b_energy);
    fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
    fChain->SetBranchAddress("numtraces", &numtraces, &b_numtraces);
    fChain->SetBranchAddress("packet_id", &packet_id, &b_packet_id);
    fChain->SetBranchAddress("runtime", &runtime, &b_runtime);
    fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
    fChain->SetBranchAddress("to_abs_mu_usec", &to_abs_mu_usec, &b_to_abs_mu_usec);
    fChain->SetBranchAddress("to_dt_mu_usec", &to_dt_mu_usec, &b_to_dt_mu_usec);
    fChain->SetBranchAddress("to_master_sec", &to_master_sec, &b_to_master_sec);
    fChain->SetBranchAddress("to_mu_sec", &to_mu_sec, &b_to_mu_sec);
    fChain->SetBranchAddress("to_mu_usec", &to_mu_usec, &b_to_mu_usec);
    fChain->SetBranchAddress("to_start_sec", &to_start_sec, &b_to_start_sec);
    fChain->SetBranchAddress("to_start_usec", &to_start_usec, &b_to_start_usec);
    fChain->SetBranchAddress("ts_maxticks", &ts_maxticks, &b_ts_maxticks);
    fChain->SetBranchAddress("ts_pps", &ts_pps, &b_ts_pps);
    fChain->SetBranchAddress("ts_ticks", &ts_ticks, &b_ts_ticks);
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
    fChain->Branch("deadtime", &deadtime);
    fChain->Branch("dr_maxticks", &dr_maxticks);
    fChain->Branch("dr_start_pps", &dr_start_pps);
    fChain->Branch("dr_start_ticks", &dr_start_ticks);
    fChain->Branch("dr_stop_pps", &dr_stop_pps);
    fChain->Branch("dr_stop_ticks", &dr_stop_ticks);
    fChain->Branch("energy", &energy);
    fChain->Branch("ievt", &ievt);
    fChain->Branch("numtraces", &numtraces);
    fChain->Branch("packet_id", &packet_id);
    fChain->Branch("runtime", &runtime);
    fChain->Branch("timestamp", &timestamp);
    fChain->Branch("to_abs_mu_usec", &to_abs_mu_usec);
    fChain->Branch("to_dt_mu_usec", &to_dt_mu_usec);
    fChain->Branch("to_master_sec", &to_master_sec);
    fChain->Branch("to_mu_sec", &to_mu_sec);
    fChain->Branch("to_mu_usec", &to_mu_usec);
    fChain->Branch("to_start_sec", &to_start_sec);
    fChain->Branch("to_start_usec", &to_start_usec);
    fChain->Branch("ts_maxticks", &ts_maxticks);
    fChain->Branch("ts_pps", &ts_pps);
    fChain->Branch("ts_ticks", &ts_ticks);
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
