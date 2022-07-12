#include "TBEvent.hxx"
ClassImp(TBEvent)

TBEvent::TBEvent(TString runName ): TNamed(runName,runName)
{
  clear();
}


//TBEvent::~TBEvent(){}

