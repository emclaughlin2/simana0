#ifndef MDCTREEMAKER_H
#define MDCTREEMAKER_H

#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include <calotrigger/TriggerAnalyzer.h>

class PHCompositeNode;
class CentralityInfo;
class MDCTreeMaker : public SubsysReco
{
 public:

  MDCTreeMaker(const std::string &name = "MDCTreeMaker", const std::string &filename = "output.root", const int dataormc = 0, const int debug = 0, const int correct = 1, const std::string &emcal_node = "TOWERINFO_CALIB_CEMC", const std::string &ihcal_node = "TOWERINFO_CALIB_HCALIN", const std::string &ohcal_node = "TOWERINFO_CALIB_HCALOUT");

  virtual ~MDCTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void set_useMBD(bool flag) {
    use_mbd = flag;
  }

  void set_useEMCal(bool flag) {
    use_emcal = flag; 
  }

  void set_useHCal(bool flag) {
    use_hcal = flag;
  }

 private:
  bool use_emcal = true; 
  bool use_hcal = true; 
  bool use_mbd = true; 
  int _evtct;
  int _correct;
  int _debug;
  int _dataormc;
  TFile *_f;
  TTree *_tree;
  std::string _foutname;
  std::string _emcal_node;
  std::string _ihcal_node;
  std::string _ohcal_node;
  TriggerAnalyzer *triggeranalyzer;
  int sectorem;
  int sectorih;
  int sectoroh;
  int sectormb;
  float emcalen[24576];
  float emcalsyst1[24576];
  float emcalsyst2[24576];
  float emcalsyst3[24576];
  float emcalsyst3d[24576];
  float emcalsyst4[24576];
  float ihcalsyst1[1536];
  float ihcalsyst2[1536];
  float ihcalsyst3[1536];
  float ihcalsyst4[1536];
  float ohcalsyst1[1536];
  float ohcalsyst2[1536];
  float ohcalsyst3[1536];
  float ohcalsyst4[1536];
  float ihcalen[1536];
  float ohcalen[1536];
  float emcalt[24576];
  float ihcalt[1536];
  float ohcalt[1536];
  int emcaladc[24576];
  int ihcaladc[1536];
  int ohcaladc[1536];
  int emcalzsadc[24576];
  int ihcalzsadc[1536];
  int ohcalzsadc[1536];
  float emcalzs[24576];
  float ihcalzs[1536];
  float ohcalzs[1536];
  float emcalpos[24576][3];
  float ihcalpos[1536][3];
  float ohcalpos[1536][3];
  int emcaletabin[24576];
  int ihcaletabin[1536];
  int ohcaletabin[1536];
  int emcalphibin[24576];
  int ihcalphibin[1536];
  int ohcalphibin[1536];
  int emcalstatus[24576];
  int ihcalstatus[1536];
  int ohcalstatus[1536];
  float mbenrgy[256];
  float mbtime[256];
  float mbd_total_q;
  int centbin;
  bool isMinBias;
  int npart;
  int ncoll;
  float bimp;
  float track_vtx[3];
  float svtx_vtx[3];
  float mbd_vtx[3];
  float truth_vtx[3];
  float emetacor[24576];
  float ihetacor[1536];
  float ohetacor[1536];
  float emchi2[24576];
  float ihchi2[1536];
  float ohchi2[1536];
  bool emishot[24576];
  bool ihishot[1536];
  bool ohishot[1536];
  int truthpar_n;
  int truthpar_nh;
  float truthpar_pz[100000];
  float truthpar_pt[100000];
  float truthpar_e[100000];
  float truthpar_eta[100000];
  float truthpar_phi[100000];
  int truthpar_pid[100000];
  float truthparh_e[100000];
  float truthparh_pt[100000];
  float truthparh_pz[100000];
  float truthparh_eta[100000];
  float truthparh_phi[100000];
  int truthparh_id[100000];
  int hotmap[3][96][256] = {0};
  int goodmap[3][96][256] = {0};
  std::vector<int> baryons{2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,
      3224,3214,3114,3322,3312,3324,3314,3334,4122,4222,4212,4112,4224,4214,
      4114,4232,4312,4324,4314,4332,4334,4412,4422,4414,4424,4432,4434,4444,
      5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314,5324,5332,
      5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,
      5514,5524,5532,5534,5542,5544,5554};
};

#endif // MDCTREEMAKER
