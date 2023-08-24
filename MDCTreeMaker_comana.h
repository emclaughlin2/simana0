#ifndef MDCTREEMAKER_H
#define MDCTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"

class PHCompositeNode;

class MDCTreeMaker : public SubsysReco
{
 public:

  MDCTreeMaker(const std::string &name = "MDCTreeMaker");

  virtual ~MDCTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


 private:

  TFile *_f;
  TTree *_tree;
  std::string _foutname;

  //int truthjet_n;
  int sectorem;
  int sectorih;
  int sectoroh;
  /*
  float truthjet_pt[1000];
  float truthjet_et[1000];
  float truthjet_ph[1000];
  int reco_jet_n;
  float reco_jet_pt[1000];
  float reco_jet_et[1000];
  float reco_jet_ph[1000];
  int truthpar_n;
  float truthpar_px[100000];
  float truthpar_py[100000];
  float truthpar_pz[100000];
  */
  int truthpar_n1;
  /*
  float truthpar_pt[100000];
  float truthpar_et[100000];
  float truthpar_ph[100000];
  int truthpar_id[100000];
  */
  float truthpar_pt1[100000];
  float truthpar_et1[100000];
  float truthpar_ph1[100000];
  int truthpar_id1[100000];
  /*
  int truthpar_j[100000];
  float emfrac[10000];
  int truthpar_em[100000];
  int truthpar_em1[100000];
  `*/
  float emcalen[100000];
  float ihcalen[100000];
  float ohcalen[100000];
  float emcalet[100000];
  float ihcalet[100000];
  float ohcalet[100000];
  float emcalph[100000];
  float ihcalph[100000];
  float ohcalph[100000];
  float prob[100000];
  float chi2[100000];
  float _cluster_E[10000];
  float _cluster_phi[10000];
  float _cluster_eta[10000];
  int _cluster_ntower[10000];
  float cluster_Ecore[10000];
  float clustowE[10000];
  float clustowet[10000];
  float clustowph[10000];
  int nclus;
  int npart;
  int ncoll;
  float bimp;
  int bestclus;
  int bcnt;
  int bcmt;
  float bctet[100];
  float bctph[100];
  float bcten[100];
  float qsum;
  int npmt;
  float mbdsumE;
  float sumoh;
  float sumih;
  float sumem;
};

#endif // MDCTREEMAKER