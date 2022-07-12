/**
** MG, July 2020
**/
#ifndef TBRUN_DEFINED
#define TBRUN_DEFINED
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <TNamed.h>
#include <TTree.h>
#include "TBEvent.hxx"
#include "TDet.hxx"

using namespace std;

// class to store info for the run

class TBRun : public TNamed
{
public:
  TBRun(TString runName = "run0");
  ~TBRun();

  void clear();
  TTree *btree;
  TBEvent *bevent;


  //void dumpEvent(ofstream &dumpFile);

  void fill()
  {
    btree->Fill();
  }

  
  void print()
  {
    bevent->print();
  }

  ClassDef(TBRun, 1)
};
#endif
