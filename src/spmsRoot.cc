//  MGold, UNM
//  July 2022
#include <stdio.h>
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
size_t maxEvents;
int nvalue = 1;
int nflat = 48;

unsigned short xval[3404];

bool goodName(TString name)
{
  TObjArray *array = ftree->GetListOfBranches();
  for (int i = 0; i < array->GetEntries(); ++i)
  {
    if (name == TString(array->At(i)->GetName()))
      return true;
  }
  return false;
}

TString getType(DataType type)
{
  TString stype("");
  if (type == PredType::NATIVE_SHORT)
  {
    stype = TString("S");
  }
  if (type == PredType::NATIVE_USHORT)
  {
    stype = TString("sNATIVE_USHORT");
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

int getNFlattened(int iset)
{
  int nf = 0;
  DataSet *dataSet = dataSets[iset];
  if (!dataSet)
  {
    printf("NULL dataSet\n");
    return nf;
  }
  if (dataNames[iset] != TString("flattened_data"))
    return nf;
  //
  cout << iset << "  " << dataNames[iset];
  DataSpace dataSpace = dataSet->getSpace();
  H5D_space_status_t hstatus;
  dataSet->getSpaceStatus(hstatus);
  int ndims = dataSpace.getSimpleExtentNdims();
  extdims = new hsize_t[ndims];
  ndims = dataSpace.getSimpleExtentDims(extdims);
  cout << " ndims " << ndims;
  for (int ir = 0; ir < ndims; ++ir)
    printf("  extdims[%i]=%i ", ir, extdims[ir]);

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
  for (int entry = 0; entry < extdims[0]; ++entry)
  {
    offset[0] = entry;
    dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);

    hsize_t start, end;
    const hssize_t off = entry;
    dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);
    dataSpace.getSelectBounds(&start, &end);
    if (!dataSpace.selectValid())
    {
      printf("ERROR on entry %lli !!\n", entry);
      return nf;
    }

    mem_out->getSelectBounds(&start, &end);
    int nout = cnt[0];
    if (ndims > 1)
      nout = cnt[1];
    const DataType type = dataSet->getDataType();
    int xval;
    dataSet->read(&xval, dataSet->getDataType(), *mem_out, dataSpace);
    cout << entry << ") " << nf << " ";
    if (xval != entry)
      break;
    nf = xval;
  }
  nf += 1;
  cout << " returning nf = " << nf << endl;
  return nf;
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
        if (goodName(theDataName))
          cout << "..." << ig << " " << theDataName;
      }
    }
  }
}

void getDataset(int ig, int index)
{
  Group *grp = groups[ig];
  DataSet *dset = NULL;
  TString theDataName = TString(grp->getObjnameByIdx(index));
  if (!goodName(theDataName))
  {
    cout << " skip " << theDataName << endl;
    return;
  }

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
    nvalue = extdims[1];
  }
  dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);
  if (!dataSpace.selectValid())
  {
    printf("ERROR!!\n");
    return;
  }
  cout << "  " << ig << " "
       << groupNames[ig] << " index " << index << " " << grp->getObjnameByIdx(index) << "ndims= " << ndims;
  for (int ir = 0; ir < ndims; ++ir)
    printf("  extdims[%i] = %i  ", ir, extdims[ir]);
  cout << endl;
  dataSets.push_back(dset);
  dataNames.push_back(grp->getObjnameByIdx(index));
  return;
}

int_fast64_t showDataSets()
{
  int nf = 0;
  for (unsigned iset = 0; iset < dataSets.size(); ++iset)
  {
    DataSet *dataSet = dataSets[iset];
    if (!dataSet)
    {
      printf("NULL dataSet\n");
      return nf;
    }
    char *ptr = NULL;
    ptr = ftree->GetBranch(dataNames[iset])->GetAddress();
    if (!ptr)
    {
      printf("NULL PTR\n");
      return nf;
    }
    if (dataNames[iset] == TString("flattened_data"))
      nf = iset;
    cout << iset << "  " << dataNames[iset] << " Id = " << int(dataSet->getId());
    cout << " type class " << int(dataSet->getTypeClass());
    cout << " type size " << int(dataSet->getDataType().getSize());
    cout << " type  " << getType(dataSet->getDataType()) << endl;
  }
  printf("\n\t total sets %lu nf = %i \n", dataSets.size(), nf);
  return nf;
}

void getEntry(int iset, long long entry)
{
  DataSet *dataSet = dataSets[iset];
  if (!dataSet)
  {
    printf("NULL dataSet\n");
    return;
  }
  if (dataNames[iset] == TString("flattened_data"))
    return;
  char *ptr = NULL;
  ptr = ftree->GetBranch(dataNames[iset])->GetAddress();
  if (!ptr)
  {
    printf("NULL PTR\n");
    return;
  }
  // cout << iset << "  " << dataNames[iset];
  DataSpace dataSpace = dataSet->getSpace();
  H5D_space_status_t hstatus;
  dataSet->getSpaceStatus(hstatus);
  int ndims = dataSpace.getSimpleExtentNdims();
  extdims = new hsize_t[ndims];
  ndims = dataSpace.getSimpleExtentDims(extdims);
  /*
  cout << " ndims " << ndims;
  for (int ir = 0; ir < ndims; ++ir)
    printf("  extdims[%i]=%i ", ir, extdims[ir]);
  */
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
  /*cout << " is simple  " << dataSpace.isSimple() << " "; // it is simple
  for (int ir = 0; ir < ndims; ++ir)
    printf("  count [%i] = %llu  off [%i] = %llu ; ", ir, cnt[ir], ir, offset[ir]);
  */
  dataSpace.selectHyperslab(H5S_SELECT_SET, cnt, offset);
  // dataSpace.getSelectBounds(&start, &end);
  // cout << "entry " << entry << " off " << off << " start " << start << " end " << end;
  if (!dataSpace.selectValid())
  {
    printf("ERROR on entry %lli !!\n", entry);
    return;
  }

  mem_out->getSelectBounds(&start, &end);
  int nout = cnt[0];
  if (ndims > 1)
    nout = cnt[1];
  const DataType type = dataSet->getDataType();
  // cout << " mem start " << start << " end " << end << " type " << int(dataSet->getDataType().getSize()) << "  ";
  dataSet->read(ptr, dataSet->getDataType(), *mem_out, dataSpace);
  // if(dataNames[iset]!=TString("values")) dataSet->read(ptr, dataSet->getDataType(), *mem_out, dataSpace);
  // else  dataSet->read(&xval, dataSet->getDataType(), *mem_out, dataSpace);
  // cout << " entry " << entry << " ==  " << ftree->GetBranch(dataNames[iset])->FindLeaf(dataNames[iset])->GetValue() << endl;

  // for (unsigned iv = 0; iv < tspms->values.size(); ++iv )
  //  tspms->values[iv] = xval[iv];
  /*
  if (dataNames[iset] == TString("values"))
  {
    cout << " values: " << tspms->values[0] << " " << tspms->values[tspms->values.size() - 1] << endl;
  }
  */
}

void event(Long64_t ientry)
{
  vector<double> digi;
  for (int i = 0; i < tspms->values.size(); ++i)
    digi.push_back(double(tspms->values[i] - tspms->baseline));
  finder->fevent(ientry, digi);
  // if (ientry / 100 * 100 == ientry) printf(" event idet %i  event %lld nhits %d   \n ", tspms->channel, ientry, bevent->nhits);
  if (ientry < 10)
  {
    finder->plotWave(tspms->channel, ientry);
    finder->plotEvent(tspms->channel, ientry);
  }
  // printf(" DONE event idet %i  event %lld nhits %d   \n ", tspms->channel, ientry, bevent->nhits);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <file> "
              << " maxEvents " << std::endl;
    return 1;
  }

  TString tag = TString(argv[1]);
  maxEvents = 0;

  if (argc > 2)
    maxEvents = atoi(argv[2]);
  cout
      << "reading  " << tag.Data() << " maxEvents " << maxEvents << endl;

  TString inFileName = TString("spms/") + tag + TString(".lh5");
  TString outFileName = tag + TString(".root");
  cout << inFileName << " " << outFileName << endl;

  // Open input HDF5 file first make sure it exits
  FILE *file;
  if (!(file = fopen(inFileName.Data(), "r")))
  {
    cout << " no file " << inFileName << endl;
    return 0;
  }
  fclose(file);
  H5File *h5 = new H5File(inFileName, H5F_ACC_RDONLY);

  // open output file and make data structure
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

  getDataNames(h5);

  for (unsigned i = 0; i < groups.size(); ++i)
  {
    cout << " group " << groupNames[i] << endl;
    for (UInt_t ig = 0; ig < groups[i]->getNumObjs(); ++ig)
      getDataset(i, ig);
  }

  nflat = getNFlattened(showDataSets());
  // set vector array sizes
  tspms->setArrays(nvalue, nflat);
  tspms->clear();
  tspms->init(ftree);

  finder = new hitFinder(fout, tspms, btree, bevent, tag);

  // get number of entries in file
  hsize_t entries;
  int ndims = dataSets[0]->getSpace().getSimpleExtentDims(&entries, NULL);

  cout << " \n \t file has total of  " << entries;
  ;
  if (maxEvents > 0)
    entries = maxEvents;
  cout << " ... looping over  " << entries << " dataSets.size() " << dataSets.size() << endl;

  for (Long64_t ientry = 0; ientry < entries; ++ientry)
  {
    if (ientry / 100 * 100 == ientry)
      cout << ".... " << ientry << " ftree " << ftree->GetEntries() << endl;
    tspms->clear();
    for (int iset = 0; iset < dataSets.size(); ++iset)
      getEntry(iset, ientry);

    // cout << " " << tspms->values[0] << " " << tspms->values[tspms->values.size() - 1] << endl;
    ftree->Fill();
    event(ientry);
  }

  h5->close();
  cout << " done " << ftree->GetEntries() << endl;
  fout->Write();
  fout->Close();
}
