//  MGold, UNM
//  July 2022

#include <iostream>
#include <string>
#include <vector>

#include "H5Cpp.h"

#include <TString.h>
#include <TFile.h>
#include "TSpms.hxx"
#include "TBEvent.hxx"
#include "hitFinder.hxx"

using namespace std;
using namespace H5;

DataSpace dataSpace;
DataSpace *mem_out;
UShort_t *buff;
hsize_t *extdims;
hsize_t *offset;
hsize_t *cnt;
hsize_t ndims;
vector<DataSet *> dataSets;
vector<TString> dataNames;
vector<Group *> groups;
vector<TString> groupNames;
TTree *ftree;
TSpms *tspms;
TBEvent *bevent;
hitFinder *finder;

// this routine not used but here for reference
TString
getType(DataType type, int nout)
{
  TString stype("");
  if (type == PredType::NATIVE_SHORT)
  {
    stype = TString("S");
  }
  if (type == PredType::NATIVE_USHORT)
  {
    stype = TString("s");
  }
  if (type == PredType::NATIVE_INT)
  {
    if (type.getSize() == 16)
      stype = TString("S");
    else
      stype = TString("I");
  }
  if (type == PredType::NATIVE_UINT)
  {
    if (type.getSize() == 16)
      stype = TString("s");
    else
      stype = TString("i");
  }
  if (type == PredType::NATIVE_LONG)
  {
    stype = TString("G");
  }
  if (type == PredType::NATIVE_ULONG)
  {
    stype = string("g");
  }
  if (type == PredType::NATIVE_FLOAT)
  {
    stype = string("F");
  }
  if (type == PredType::NATIVE_DOUBLE)
  {
    stype = string("D");
  }
  if (type == PredType::NATIVE_CHAR)
  {
    stype = string("C");
  }
  if (type == PredType::NATIVE_SCHAR)
  {
    stype = string("c");
  }
  return stype;
}

void getDataNames(H5File *hfile)
{
  Group *main_group = new Group(hfile->openGroup("/spms"));
  for (UInt_t i = 0; i < main_group->getNumObjs(); i++)
  {
    if (H5G_GROUP == main_group->getObjTypeByIdx(i))
    {
      Group *grp = new Group(main_group->openGroup(main_group->getObjnameByIdx(i)));

      groupNames.push_back(TString(main_group->getObjnameByIdx(i)));
      groups.push_back(grp);
    }
  }
  // get other groups
  main_group = new Group(hfile->openGroup("/spms/raw"));
  for (UInt_t i = 0; i < main_group->getNumObjs(); i++)
  {
    if (H5G_GROUP == main_group->getObjTypeByIdx(i))
    {
      Group *grp = new Group(main_group->openGroup(main_group->getObjnameByIdx(i)));
      groupNames.push_back(TString(main_group->getObjnameByIdx(i)));
      groups.push_back(grp);
    }
  }

  for (unsigned i = 0; i < groups.size(); ++i)
  {
    cout << groupNames[i] << "  index  " << i << " " << groups[i]->getObjnameByIdx(0) << " #obj " << groups[i]->getNumObjs();
    for (UInt_t ig = 0; ig < groups[i]->getNumObjs(); ++ig)
    {
      if (H5G_DATASET == groups[i]->getObjTypeByIdx(ig))
      {
        // subGroupIndex.push_back(i);
        //  form data name
        TString theDataName = TString(groups[i]->getObjnameByIdx(ig));
        cout << "..." << ig << " " << theDataName;
      }
    }
    cout << endl;
  }
}

void getDataset(int ig, int index)
{
  Group *grp = groups[ig];
  DataSet *dset = NULL;

  if (!(H5G_DATASET == grp->getObjTypeByIdx(index)))
    return;

  cout << "\t getDataset " << groupNames[ig] << " " << grp->getObjnameByIdx(index);
  dset = new DataSet(grp->openDataSet(grp->getObjnameByIdx(index)));
  dataSpace = dset->getSpace();
  ndims = dataSpace.getSimpleExtentNdims();
  extdims = new hsize_t[ndims];
  ndims = dataSpace.getSimpleExtentDims(extdims, NULL);
  for (int ir = 0; ir < ndims; ++ir)
    printf("  ndims %llu [%i] = %llu ", ndims, ir, extdims[ir]);
  cout << "  ; ";
  offset = new hsize_t[ndims];
  cnt = new hsize_t[ndims];
  cnt[0] = 1;
  offset[0] = 0;
  if (ndims > 1)
  {
    offset[1] = 0;
    cnt[1] = extdims[1];
  }
  dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);
  if (!dataSpace.selectValid())
  {
    printf("ERROR!!\n");
    return;
  }
  cout << "  " << ig << " " << groupNames[ig] << " index " << index << " " << grp->getObjnameByIdx(index) << endl;
  dataSets.push_back(dset);
  dataNames.push_back(grp->getObjnameByIdx(index));
  return;
}

void getEntry(int iset, long long entry)
{
  DataSet *dataSet = dataSets[iset];
  if (!dataSet)
  {
    printf("NULL dataSet\n");
    return;
  }
  char *ptr = NULL;
  ptr = ftree->GetBranch(dataNames[iset])->GetAddress();
  if (!ptr)
  {
    printf("NULL PTR\n");
    return;
  }
  DataSpace dataSpace = dataSet->getSpace();
  H5D_space_status_t hstatus;
  dataSet->getSpaceStatus(hstatus);
  // cout << " dataset stat " << int(hstatus) << " ";
  int ndims = dataSpace.getSimpleExtentNdims();
  extdims = new hsize_t[ndims];
  ndims = dataSpace.getSimpleExtentDims(extdims, NULL);
  offset = new hsize_t[ndims];
  cnt = new hsize_t[ndims];
  cnt[0] = 1;
  offset[0] = 0;
  if (ndims > 1)
  {
    offset[1] = 0;
    cnt[1] = extdims[1];
  }

  // copy before setting for this event
  dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);
  mem_out = new DataSpace();
  mem_out->copy(dataSpace);

  // select this event
  offset[0] = entry;
  dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);

  hsize_t start, end;
  const hssize_t off = entry;
  // cout << " is simple  " << dataSpace.isSimple() << " ";// it is simple
  // for (int ir = 0; ir < ndims; ++ir)
  //   printf("  count [%i] = %llu  off [%i] = %llu ; ", ir, cnt[ir], ir, offset[ir]);
  dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);
  // dataSpace.getSelectBounds(&start, &end);
  // cout << "entry " << entry << " off " << off << " start " << start << " end " << end;
  if (!dataSpace.selectValid())
  {
    printf("ERROR on entry %lli !!\n", entry);
    return;
  }

  // mem_out->getSelectBounds(&start, &end);
  int nout = cnt[0];
  if (ndims > 1)
    nout = cnt[1];
  const DataType type = dataSet->getDataType();
  // cout << " mem start " << start << " end " << end << " type " << int(dataSet->getDataType().getSize()) << "  ";
  dataSet->read(ptr, dataSet->getDataType(), *mem_out, dataSpace);
  // cout << " entry " << entry << "  " <<  dataNames[iset] << " ==  " << ftree->GetBranch(dataNames[iset])->FindLeaf(dataNames[iset])->GetValue()  << endl;
}

void event(Long64_t ientry)
{
  vector<double> digi;
  for (int i = 0; i < tspms->values.size(); ++i)
    digi.push_back(double(tspms->values[i] - tspms->baseline));
  finder->fevent(ientry, digi);
  if (ientry / 10* 10 == ientry)
    printf(" event idet %i  event %lld nhits %d   \n ", tspms->channel, ientry, bevent->nhits);
  if (ientry / 100 * 100 == ientry)
  {
    finder->plotWave(tspms->channel, ientry);
    // finder->plotEvent(tspms->channel,ientry);
  }
}

int main(int argc, char *argv[])
{

  TString tag = TString("20220413-1000");
  cout
      << "reading  " << tag.Data() << endl;

  TString inFileName = tag + TString(".lh5");
  TString outFileName = tag + TString(".root");

  // Open input HDF5 file
  H5File *h5 = new H5File(inFileName, H5F_ACC_RDONLY);

  getDataNames(h5);

  TFile *fout = new TFile(outFileName, "recreate");
  ftree = new TTree("spms", " spms ");
  tspms = new TSpms();
  tspms->create(ftree);
  ftree->GetListOfBranches()->ls();
  // hit analysis
  TTree *btree = new TTree("HitTree", " hit data ");
  bevent = new TBEvent(tag);
  btree->Branch("bev", &bevent);
  btree->GetListOfBranches()->ls();

  for (unsigned i = 0; i < groups.size(); ++i)
  {
    cout << " group " << groupNames[i] << endl;
    for (UInt_t ig = 0; ig < groups[i]->getNumObjs(); ++ig)
      getDataset(i, ig);
  }

  finder = new hitFinder(fout, tspms, btree, bevent, tag);

  // get number of entries in file
  hsize_t entries;
  int ndims = dataSets[0]->getSpace().getSimpleExtentDims(&entries, NULL);

  cout << " \n \t looping over  " << entries << endl;

  for (Long64_t ientry = 0; ientry < entries; ++ientry)
  {
    tspms->clear();
    for (int iset = 0; iset < dataSets.size(); ++iset)
      getEntry(iset, ientry);
    ftree->Fill();
    event(ientry);
    if (ientry / 100 * 100 == ientry)
      cout << ".... " << ientry << " ftree " << ftree->GetEntries() << endl;
  }

  h5->close();
  fout->Write();
  fout->Close();
  cout << " done " << endl;
}
