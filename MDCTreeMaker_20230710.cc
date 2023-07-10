#include "MDCTreeMaker.h"
#include <ffaobjects/EventHeaderv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <HepMC/GenEvent.h>

#include <jetbackground/TowerBackgroundv1.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <centrality/CentralityInfov1.h>

using namespace std;

//____________________________________________________________________________..
MDCTreeMaker::MDCTreeMaker(const std::string &name):
  SubsysReco(name)
{
  _foutname = name;  
}

//____________________________________________________________________________..
MDCTreeMaker::~MDCTreeMaker()
{

}

//____________________________________________________________________________..
int MDCTreeMaker::Init(PHCompositeNode *topNode)
{
  
  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree = new TTree("ttree","a persevering date tree");
  /*
  _tree->Branch("truthjet_n",&truthjet_n,"truthjet_n/I");
  _tree->Branch("truthjet_pt",truthjet_pt,"truthjet_pt[truthjet_n]/F");
  _tree->Branch("truthjet_et",truthjet_et,"truthjet_et[truthjet_n]/F");
  _tree->Branch("truthjet_ph",truthjet_ph,"truthjet_ph[truthjet_n]/F");
  _tree->Branch("reco_jet_n",&reco_jet_n,"reco_jet_n/I");
  _tree->Branch("reco_jet_pt",reco_jet_pt,"reco_jet_pt[reco_jet_n]/F");
  _tree->Branch("reco_jet_et",reco_jet_et,"reco_jet_et[reco_jet_n]/F");
  _tree->Branch("reco_jet_ph",reco_jet_ph,"reco_jet_ph[reco_jet_n]/F");
  _tree->Branch("truthpar_n",&truthpar_n,"truthpar_n/I");
  _tree->Branch("truthpar_pt",truthpar_pt,"truthpar_pt[truthpar_n]/F");
  _tree->Branch("truthpar_et",truthpar_et,"truthpar_et[truthpar_n]/F");
  _tree->Branch("truthpar_ph",truthpar_ph,"truthpar_ph[truthpar_n]/F");
  _tree->Branch("truthpar_id",truthpar_id,"truthpar_id[truthpar_n]/I");
  */
  _tree->Branch("truthpar_n",&truthpar_n1,"truthpar_n1/I");
  _tree->Branch("truthpar_pt",truthpar_pt1,"truthpar_pt1[truthpar_n1]/F");
  _tree->Branch("truthpar_et",truthpar_et1,"truthpar_et1[truthpar_n1]/F");
  _tree->Branch("truthpar_ph",truthpar_ph1,"truthpar_ph1[truthpar_n1]/F");
  _tree->Branch("truthpar_id",truthpar_id1,"truthpar_id1[truthpar_n1]/I");
  /*
  _tree->Branch("truthpar_j",truthpar_j,"truthpar_j[truthpar_n]/I");
  _tree->Branch("emfrac",emfrac,"emfrac[truthjet_n]/F");
  _tree->Branch("truthpar_em",truthpar_em,"truthpar_em[truthpar_n]/I");
  _tree->Branch("truthpar_em1",truthpar_em1,"truthpar_em1[truthpar_n1]/I");
  */
  _tree->Branch("sectorem",&sectorem,"sectorem/I");
  _tree->Branch("sectorih",&sectorih,"sectorih/I");
  _tree->Branch("sectoroh",&sectoroh,"sectoroh/I");
  _tree->Branch("emcalen",emcalen,"emcalen[sectorem]/F");
  _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F");
  _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
  _tree->Branch("emcalet",emcalet,"emcalet[sectorem]/F");
  _tree->Branch("ihcalet",ihcalet,"ihcalet[sectorih]/F");
  _tree->Branch("ohcalet",ohcalet,"ohcalet[sectoroh]/F");
  _tree->Branch("emcalph",emcalph,"emcalph[sectorem]/F");
  _tree->Branch("ihcalph",ihcalph,"ihcalph[sectorih]/F");
  _tree->Branch("ohcalph",ohcalph,"ihcalph[sectoroh]/F");
  _tree->Branch("nclus",&nclus,"nclus/I");
  _tree->Branch("prob",prob,"prob[nclus]/F");
  _tree->Branch("chi2",chi2,"chi2[nclus]/F");
  _tree->Branch("cluster_E",_cluster_E,"cluster_E[nclus]/F");
  _tree->Branch("cluster_eta",_cluster_eta,"cluster_eta[nclus]/F");
  _tree->Branch("cluster_phi",_cluster_phi,"cluster_phi[nclus]/F");
  _tree->Branch("cluster_ntower",_cluster_ntower,"cluster_ntower[nclus]/I");
  _tree->Branch("cluster_Ecore",cluster_Ecore,"cluster_Ecore[nclus]/F");
  //_tree->Branch("clustowE",clustowE,"clustowE[cluster_ntower[nclus]]/F");
  //_tree->Branch("clustowet",clustowet,"clustowet[cluster_ntower[nclus]]/F");
  //_tree->Branch("clustowph",clustowph,"clustowph[cluster_ntower[nclus]]/F");
  _tree->Branch("npart",&npart,"npart/I");
  _tree->Branch("ncoll",&ncoll,"ncoll/I");
  _tree->Branch("bimp",&bimp,"bimp/F");
  _tree->Branch("bestclus",&bestclus,"bestclus/I");
  _tree->Branch("bcnt",&bcnt,"bcnt/I");
  _tree->Branch("bcmt",&bcmt,"bcmt/I");
  _tree->Branch("bctet",bctet,"bctet[bcnt]/F");
  _tree->Branch("bctph",bctph,"bctph[bcnt]/F");
  _tree->Branch("bcten",bcten,"bcten[bcnt]/F");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
/*
void climbup(int depth, HepMC::GenVertex::particle_iterator part)
{
  ++depth;
  if(!(*part)->production_vertex())
    {
      for(int i=0; i<depth; i++)
	{
	  cout << "   ";
	}
      cout << "pid=" << (*part)->pdg_id() << "barcd=" << (*part)->barcode() << " top level particle" << endl;
    }
  else
    {
      for(HepMC::GenVertex::particle_iterator newpart = (*part)->production_vertex()->particles_begin(HepMC::parents);
	  newpart != (*part)->production_vertex()->particles_end(HepMC::parents); ++newpart)
	{
	  for(int i=0; i<depth; i++)
	    {
	      //cout << "   ";
	    }
	  //cout << "pid=" << (*newpart)->pdg_id() << " barcd=" << (*newpart)->barcode() <<" production position x y z = "
	  //     << (*part)->production_vertex()->position().x() << " " << (*part)->production_vertex()->position().y() <<
	  //  " " << (*part)->production_vertex()->position().z() << endl;
	  //climbup(depth, newpart);
	}
    }
}
*/
//____________________________________________________________________________..
int MDCTreeMaker::process_event(PHCompositeNode *topNode)
{

  // here is where the code goes 
//  truthjet_n = 0;
//  reco_jet_n = 0;
//  truthpar_n = 0;
  truthpar_n1 = 0;
  sectorem = 0;
  sectorih = 0;
  sectoroh = 0;
  nclus = 0;
  for(int i=0; i<sizeof(emcalen)/sizeof(emcalen[0]); i++)
    {
      emcalen[i] = 0;
      emcalet[i] = 0;
      emcalph[i] = 0;
      ihcalen[i] = 0;
      ihcalet[i] = 0;
      ihcalph[i] = 0;
      ohcalen[i] = 0;
      ohcalet[i] = 0;
      ohcalph[i] = 0;
  //    truthpar_px[i] = 0;
  //    truthpar_py[i] = 0;
  //    truthpar_pz[i] = 0;
  //    truthpar_pt[i] = 0;
  //    truthpar_et[i] = 0;
  //    truthpar_ph[i] = 0;
//      truthpar_id[i] = 0;
      truthpar_pt1[i] = 0;
      truthpar_et1[i] = 0;
      truthpar_ph1[i] = 0;
      truthpar_id1[i] = 0;
  //    truthpar_j[i] = 0;
  //    truthpar_em[i] = 0;
  //    truthpar_em1[i] = 0;
      prob[i] = 0;
      chi2[i] = 0;
    }
/*  
  for(int i=0; i<100; i++)
    {
      emfrac[i] = 0;
    }
  */
  {
  RawTowerContainer *towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if(towersEM && geomEM)
    {
      RawTowerContainer::ConstRange begin_end = towersEM->getTowers();
      for(RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
	  float eta = tower_geom->get_eta();
	  //cout << "it: " << sector << " eta: " << eta << " eng: " << tower->get_energy() << endl;
	  emcalen[sectorem] = tower->get_energy();
	  emcalet[sectorem] = eta;
	  emcalph[sectorem] = tower_geom->get_phi();
	  sectorem++;
	}
    }
  
  if(towersIH && geomIH)
    {
      RawTowerContainer::ConstRange begin_end = towersIH->getTowers();
      for(RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
	  float eta = tower_geom->get_eta();
	  ihcalen[sectorih] = tower->get_energy();
	  ihcalet[sectorih] = eta;
	  ihcalph[sectorih] = tower_geom->get_phi();
	  sectorih++;
	}
    }
  
  if(towersOH && geomOH)
    {
      RawTowerContainer::ConstRange begin_end = towersOH->getTowers();
      for(RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawTower *tower = rtiter->second;

	  RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());

	  ohcalen[sectoroh] = tower->get_energy();
	  ohcalet[sectoroh] = tower_geom->get_eta();
	  ohcalph[sectoroh] = tower_geom->get_phi();
	  sectoroh++;
	}
    }
  
  }
  //vector<int> pibarcodes;
  
  /*
  {
    
    JetMap* jets = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r04");

    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");

    for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter) {

      Jet* this_jet = iter->second;
      
      float jet_pt = this_jet->get_pt();
      float jet_et = this_jet->get_eta();
      float jet_ph = this_jet->get_phi();

      truthjet_pt[ truthjet_n ] = jet_pt;
      truthjet_et[ truthjet_n ] = jet_et;
      truthjet_ph[ truthjet_n ] = jet_ph;
      //cout << "BEGIN NEW JET:" << endl;
      
      for (Jet::ConstIter comp = this_jet->begin_comp(); comp !=  this_jet->end_comp(); ++comp) {
        PHG4Particle* g4particle = truthinfo->GetParticle( (*comp).second  );
	
	if(g4particle->get_pid() == 111)
	  {
	    pibarcodes.push_back(g4particle->get_barcode());
	  }
	cout << "Jet particle: pid=" << g4particle->get_pid() << " Embedid= " << truthinfo->isEmbeded(g4particle->get_track_id())
	     << " pT=" << sqrt(g4particle->get_px()*g4particle->get_px() + g4particle->get_py()*g4particle->get_py())
	     << " trkid=" << g4particle->get_track_id() << " vtxid=" << g4particle->get_vtx_id() << " parid=" << g4particle->get_parent_id()
	     << " prmid=" << g4particle->get_primary_id() << " barcd=" << g4particle->get_barcode() << endl;
	
        TLorentzVector t;
        t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());
        float truth_pt = t.Pt ();
        float truth_eta = t.Eta ();
        float truth_phi = t.Phi ();
        int truth_pid = g4particle->get_pid ();

	if(truthinfo->isEmbeded(g4particle->get_track_id()) == 1)
	  {
	    truthpar_em[truthpar_n] = truthinfo->isEmbeded(g4particle->get_track_id());
	    truthpar_pt[truthpar_n] = truth_pt;
	    truthpar_et[truthpar_n] = truth_eta;
	    truthpar_ph[truthpar_n] = truth_phi;
	    truthpar_id[truthpar_n] = truth_pid;
	    truthpar_j[truthpar_n] = truthjet_n;
	    truthpar_n++;
	    if(truth_pid == 11 || truth_pid == 22 || truth_pid == 111)
	      {
		emfrac[truthjet_n] += truth_pt/jet_pt;
	      }
	  }
      }
      truthjet_n++;
    }
  }
  
  {
    JetMap* jets = findNode::getClass<JetMap>(topNode,"AntiKt_Tower_r04");
    
    for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter) {
      
      Jet* this_jet = iter->second;
      
      float jet_pt = this_jet->get_pt();
      float jet_et = this_jet->get_eta();
      float jet_ph = this_jet->get_phi();
      
      
      
      reco_jet_pt[ reco_jet_n ] = jet_pt;
      
      reco_jet_et[ reco_jet_n ] = jet_et;
      
      reco_jet_ph[ reco_jet_n ] = jet_ph;
      
      reco_jet_n++;
      
    }
  }
  */
  {    
    //PHHepMCGenEventMap *genEventMap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    //PHHepMCGenEvent *genEvent = genEventMap->get(1);
    //HepMC::GenEvent *theEvent = genEvent->getEvent();
    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange ();    
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
      {
	PHG4Particle* g4particle = iter->second;
	
	if(truthinfo->isEmbeded(g4particle->get_track_id()) < 1) continue;
	/*
	if(g4particle->get_pid() != 111) continue;
	cout << "Geant4 particle: pid=" << g4particle->get_pid() << " Embedid=" <<  truthinfo->isEmbeded(g4particle->get_track_id())
	     << " pT=" << sqrt(g4particle->get_px()*g4particle->get_px() + g4particle->get_py()*g4particle->get_py()) <<
	  " trkid=" << g4particle->get_track_id() << " vtxid=" << g4particle->get_vtx_id() << " parid=" << g4particle->get_parent_id()
	     << " prmid=" << g4particle->get_primary_id() << " barcd=" << g4particle->get_barcode() << " energy="
	     << g4particle->get_e() << endl;
	*/
	TLorentzVector t;
	t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());
	float truth_pt = t.Pt ();
	float truth_eta = t.Eta ();
	float truth_phi = t.Phi ();
	int truth_pid = g4particle->get_pid (); // particle species     
	// take only particles from primary Pythia event                                                            
//	truthpar_em1[truthpar_n1] = truthinfo->isEmbeded(g4particle->get_track_id());
	truthpar_pt1[truthpar_n1] = truth_pt;
	truthpar_et1[truthpar_n1] = truth_eta;
	truthpar_ph1[truthpar_n1] = truth_phi;
	truthpar_id1[truthpar_n1] = truth_pid;
	truthpar_n1++;
	/*
	if(g4particle->get_pid() != 111) continue;
	if(truthinfo->isEmbeded(g4particle->get_track_id()) != 1) continue;
	counter++;
	for(HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p!= theEvent->particles_end(); ++p)
	  {
	    assert(*p);
	    if((*p)->pdg_id() != 111) continue;
	    if(!(*p)->production_vertex()) continue;
	    if(((*p)->barcode() == g4particle->get_barcode()))
	      {
		if((*p)->production_vertex())
		  {
		    
		    //cout << "Production vertex position x: " << (*p)->production_vertex()->position().x() << endl;
		    //cout << "End vertex pointer: " << (*p)->end_vertex() << endl;
		    //cout << "Production vertex position y: " << (*p)->production_vertex()->position().y() << endl;
		    //cout << "Production vertex position z: " << (*p)->production_vertex()->position().z() << endl;
		    
		  }
		//cout << "Pythia particle: pid=" << (*p)->pdg_id() << " barcd=" << (*p)->barcode() << " pT=" << sqrt((*p)->momentum().px() * (*p)->momentum().px() + (*p)->momentum().py() * (*p)->momentum().py()) << " energy=" << (*p)->momentum().e() << endl;
		for(HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()->particles_begin(HepMC::parents); mother != (*p)->production_vertex()->particles_end(HepMC::parents); ++mother)
		  {
		    //cout << "   mother pid: " << (*mother)->pdg_id() << " barcd: " << (*mother)->barcode() << " production position x y z = " << (*mother)->production_vertex()->position().x() << " " << (*mother)->production_vertex()->position().y() << " " << (*mother)->production_vertex()->position().z() << endl;
		    climbup(1, mother);
		  }
	      }
	  }
	*/
      }
  }
  
  {
    float maxE = 0;
    RawClusterContainer *clusterEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    RawTowerContainer *towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
    RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(clusterEM && geomEM && towersEM)
      {
	RawClusterContainer::Range begin_end = clusterEM->getClusters();
	for(RawClusterContainer::Iterator iter = begin_end.first; iter != begin_end.second; ++iter)
	  {
	    float theta = 3.14159 / 2.0 - atan2( iter->second->get_z() , iter->second->get_r() );
	    float eta = -1 * log( tan( theta / 2.0 ) );
	    RawCluster *clust = iter->second;
	    prob[nclus] = clust->get_prob();
	    chi2[nclus] = clust->get_chi2();
	    _cluster_eta[nclus] = eta;
	    _cluster_phi[nclus] = iter->second->get_phi();
	    _cluster_E[ nclus ] = iter->second->get_energy();
	    _cluster_ntower[ nclus ] =  iter->second->getNTowers();
	    cluster_Ecore[nclus] = iter->second->get_ecore();
	    if(_cluster_E[nclus] > maxE)
	      {
		bestclus = nclus;
		bcnt = _cluster_ntower[nclus];
		maxE = _cluster_E[nclus];
		int cit = 0;
		int tmaxE = 0;
		RawCluster::TowerConstRange ctowers = clust->get_towers();
		RawCluster::TowerConstIterator toweriter;
		for(toweriter = ctowers.first; toweriter != ctowers.second; toweriter++)
		  {
		    RawTower *tower = towersEM->getTower(toweriter->first);
		    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());

		    bcten[cit] = tower->get_energy();
		    bctet[cit] = tower_geom->get_eta();
		    bctph[cit] = tower_geom->get_phi();
		    if(bcten[cit] > tmaxE) 
		      {
			bcmt = cit;
			tmaxE = bcten[cit];
		      }
		    ++cit;
		    if(cit >= 100)
		      {
			cout << "Warning: cluster found with > 100 towers, index " << nclus << endl;
			break;
		      }
		  }
	      }
	    ++nclus;
	  }
      }
  }
  EventHeaderv1 *event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader" );
  if ( event_header ) 
    {
      npart = event_header->get_intval("npart");
      ncoll = event_header->get_intval("ncoll");    
      bimp = event_header->get_floatval("bimp");    
      std::cout << " npart / ncoll / bimp = " << npart << " / " << ncoll << " / " << bimp << std::endl;
    }


  _tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
  
}

//int MDCTreeMaker::CreateNode(PHCompositeNode *topNode)
//{

//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int MDCTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
  std::cout << "MDCTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
  {
  std::cout << "MDCTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
  std::cout << "MDCTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }

  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MDCTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
 std::cout << "MDCTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void MDCTreeMaker::Print(const std::string &what) const
{
  std::cout << "MDCTreeMaker::Print(const std::string &what) const Printing info for " << what << std::endl;
}