// -*- C++ -*-
//
// Package:    GEMLocal
// Class:      GEMLocal
// 
/**\class GEMLocal GEMLocal.cc GEMAnalysis/GEMLocal/plugins/GEMLocal.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu, 24 Oct 2013 11:19:56 GMT
// $Id$
//
//


// system include files
#include <memory>

// root include files
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h" //?
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"

#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h>
#include "Geometry/CommonTopologies/interface/StripTopology.h"
//#include "DataFormats/Math/interface/deltaR.h"

#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/DetLayers/interface/MuRodBarrelLayer.h"
#include "RecoMuon/DetLayers/interface/MuDetRod.h"
#include "RecoMuon/DetLayers/interface/MuRingForwardDoubleLayer.h"
#include "RecoMuon/DetLayers/interface/MuDetRing.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class GEMLocal : public edm::EDAnalyzer {
   public:
      explicit GEMLocal(const edm::ParameterSet&);
      ~GEMLocal();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
      virtual void endJob() ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      std::map<std::string,TH1F*> histContainer_;
      std::map<std::string,TH2F*> histContainer2D_; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GEMLocal::GEMLocal(const edm::ParameterSet& iConfig):
histContainer_()
{
   //now do what ever initialization is needed

}


GEMLocal::~GEMLocal()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GEMLocal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle<GEMRecHitCollection> GEMrecHits;
  iEvent.getByLabel("gemRecHits","",GEMrecHits);

  edm::Handle<edm::PSimHitContainer> GEMsimHits;
  iEvent.getByLabel(edm::InputTag("g4SimHits","MuonGEMHits"), GEMsimHits);

  edm::Handle<GEMDigiCollection> GEMdigis;
  iEvent.getByLabel(edm::InputTag("simMuonGEMDigis"), GEMdigis);

  edm::ESHandle<GEMGeometry> GEMgeo;
  iSetup.get<MuonGeometryRecord>().get(GEMgeo);

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

histContainer_["N_eventi"]->Fill(1);
std::cout<<" ************SIM Hits: "<<std::endl;
for (edm::PSimHitContainer::const_iterator itHit = GEMsimHits->begin(); itHit != GEMsimHits->end(); ++itHit)
  {
	histContainer_["N_SimHit"]->Fill(1);

	//DetId id = DetId(itHit->detUnitId());
	GEMDetId idGem = GEMDetId(itHit->detUnitId());
	LocalPoint lp = itHit->entryPoint();
	double strip_sim = GEMgeo->etaPartition(idGem)->strip(lp);
	
	int roll_sim = idGem.roll();
	histContainer2D_["Occupancy_SIM_partitionVSstrip"]->Fill(strip_sim, roll_sim);
	if(roll_sim==1){histContainer_["Occupncy_SIM_roll1"]->Fill(strip_sim);}
	if(roll_sim==2){histContainer_["Occupncy_SIM_roll2"]->Fill(strip_sim);}
	if(roll_sim==3){histContainer_["Occupncy_SIM_roll3"]->Fill(strip_sim);}
	if(roll_sim==4){histContainer_["Occupncy_SIM_roll4"]->Fill(strip_sim);}
	if(roll_sim==5){histContainer_["Occupncy_SIM_roll5"]->Fill(strip_sim);}
	if(roll_sim==6){histContainer_["Occupncy_SIM_roll6"]->Fill(strip_sim);}
	if(roll_sim==7){histContainer_["Occupncy_SIM_roll7"]->Fill(strip_sim);}
	if(roll_sim==8){histContainer_["Occupncy_SIM_roll8"]->Fill(strip_sim);}
	
	int PID = std::abs(itHit->particleType());
	double x_sim = itHit->localPosition().x();
	double y_sim = itHit->localPosition().y();
	//const GEMDetId id(itHit->detUnitId());
	std::cout<<idGem<<" PiD "<<PID<<std::endl;
	//std::cout<<itHit->detUnitId()<<std::endl;
		
	std::cout<<" sim strip "<<strip_sim<<" SimHit "<<" x: "<<x_sim<<" y "<<y_sim<<std::endl;
	
  }


 int Ndigis = 0;
//////////
std::cout<<" ************DIGI : "<<std::endl;
  GEMDigiCollection::DigiRangeIterator DigiDetUnit;
  for (DigiDetUnit = GEMdigis->begin(); DigiDetUnit != GEMdigis->end(); ++DigiDetUnit)
  {	

//        const GEMEtaPartition * roll = GEMgeo->etaPartition(id);
        int Nroll = (Short_t) (*DigiDetUnit).first.roll();
        std::cout<<(*DigiDetUnit).first<<""<<Nroll<<std::endl;

        GEMDigiCollection::const_iterator digiItr;
        //loop over digis of given roll
        for (digiItr = (*DigiDetUnit).second.first; digiItr != (*DigiDetUnit).second.second; ++digiItr)
        {

                double strip = (Short_t) digiItr->strip();
		histContainer_["N_strips"]->Fill(1);
		histContainer2D_["Occupancy_DIGI_partitionVSstrip"]->Fill(strip, Nroll);
                double bx = (Short_t) digiItr->bx();
                //double size = (Short_t) digiItr->size();
                //std::cout<<strip<<" "<<bx<<" "<<size<<std::endl;
                std::cout<<" strip "<<strip<<" bx "<<bx<<std::endl;
        	double strip_successiva = (Short_t) (digiItr+1)->strip();
                //std::cout<<" strip_successiva "<<strip_successiva<<std::endl;
	if (strip_successiva != (strip+1)){        
	Ndigis++;
        histContainer_["N_Digis"]->Fill(1);
	if (bx==0){histContainer_["N_Digis_bx0"]->Fill(1);}
	if(bx!=0){histContainer_["N_Digis_noise"]->Fill(1);}
	}
		GEMDetId id = (*DigiDetUnit).first;
	        const GEMEtaPartition* roll = GEMgeo->etaPartition(id);
		LocalPoint lp = roll->centreOfStrip(digiItr->strip());
    		double x_digi = (double) lp.x();
                double y_digi = (double) lp.y();
                std::cout<<" Digi "<<" x: "<<x_digi<<" y "<<y_digi<<std::endl;

/*		for (GEMRecHitCollection::const_iterator recHit = GEMrecHits->begin(); recHit != GEMrecHits->end(); ++recHit)
  		{
        	std::cout<<(*recHit).gemId()<<std::endl;
    			float x_reco = recHit->localPosition().x();
    			float y_reco = recHit->localPosition().y();
			double clusterSize = recHit->clusterSize();
                        std::cout<<" RecHit "<<" x: "<<x_reco<<" y "<<y_reco<<" cluster size "<<clusterSize<<std::endl;
		}
*/        }// end digi per roll

  }

/////////
std::cout<<" ************REC Hits: "<<std::endl;

int NrecHit_x_roll1=0, NrecHit_x_roll2=0, NrecHit_x_roll3=0, NrecHit_x_roll4=0, NrecHit_x_roll5=0, NrecHit_x_roll6=0, NrecHit_x_roll7=0, NrecHit_x_roll8=0;
int NrecHit_x_roll1_bx0=0, NrecHit_x_roll2_bx0=0, NrecHit_x_roll3_bx0=0, NrecHit_x_roll4_bx0=0, NrecHit_x_roll5_bx0=0, NrecHit_x_roll6_bx0=0, NrecHit_x_roll7_bx0=0, NrecHit_x_roll8_bx0=0;
int NrecHit_x_roll1_noise=0, NrecHit_x_roll2_noise=0, NrecHit_x_roll3_noise=0, NrecHit_x_roll4_noise=0, NrecHit_x_roll5_noise=0, NrecHit_x_roll6_noise=0, NrecHit_x_roll7_noise=0, NrecHit_x_roll8_noise=0;

for (GEMRecHitCollection::const_iterator recHit = GEMrecHits->begin(); recHit != GEMrecHits->end(); ++recHit)
 {
        std::cout<<(*recHit).gemId()<<std::endl;
	
	double x_reco = recHit->localPosition().x();
        double err_x_reco = recHit->localPositionError().xx();
        double y_reco = recHit->localPosition().y();
        double clusterSize = recHit->clusterSize();
	double bx_reco = (Short_t) recHit->BunchX();
	int recHit_roll = (*recHit).gemId().roll();
	double strip_rec = recHit->firstClusterStrip();
		
	 histContainer_["N_RecHit"]->Fill(1);
	 histContainer_["bx"]->Fill(bx_reco);
	 if (bx_reco==0){histContainer_["N_RecHit_bx0"]->Fill(1);}
	 if (bx_reco!=0){
	 //histContainer_["N_RecHit_noise"]->Fill(1);
	 histContainer_["bx_noise_bx!=0"]->Fill(bx_reco);
	 		/*if(recHit_roll==1){histContainer_["N_RecHit_noise1"]->Fill(1);}
			if(recHit_roll==2){histContainer_["N_RecHit_noise2"]->Fill(1);}
			if(recHit_roll==3){histContainer_["N_RecHit_noise3"]->Fill(1);}
			if(recHit_roll==4){histContainer_["N_RecHit_noise4"]->Fill(1);}
			if(recHit_roll==5){histContainer_["N_RecHit_noise5"]->Fill(1);}
			if(recHit_roll==6){histContainer_["N_RecHit_noise6"]->Fill(1);}
			if(recHit_roll==7){histContainer_["N_RecHit_noise7"]->Fill(1);}
			if(recHit_roll==8){histContainer_["N_RecHit_noise8"]->Fill(1);}*/
	 		}
			
	/*if(GEMsimHits->size()==0 && Ndigis!=0){
	histContainer_["CS_NOISE"]->Fill(clusterSize);
	histContainer_["N_RecHit_noise"]->Fill(1);
	histContainer_["bx_noise"]->Fill(bx_reco);
	 		if(recHit_roll==1){histContainer_["N_RecHit_noise1"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise1"]->GetEntries());}
			if(recHit_roll==2){histContainer_["N_RecHit_noise2"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise2"]->GetEntries());}
			if(recHit_roll==3){histContainer_["N_RecHit_noise3"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise3"]->GetEntries());}
			if(recHit_roll==4){histContainer_["N_RecHit_noise4"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise4"]->GetEntries());}
			if(recHit_roll==5){histContainer_["N_RecHit_noise5"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise5"]->GetEntries());}
			if(recHit_roll==6){histContainer_["N_RecHit_noise6"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise6"]->GetEntries());}
			if(recHit_roll==7){histContainer_["N_RecHit_noise7"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise7"]->GetEntries());}
			if(recHit_roll==8){histContainer_["N_RecHit_noise8"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise8"]->GetEntries());}
	 		}*/
        for (int count =0; count<clusterSize; ++count){       
	histContainer2D_["Occupancy_RECO_partitionVSstrip"]->Fill(strip_rec+count, recHit_roll);
	if(recHit_roll==1){histContainer_["Occupncy_RECO_roll1"]->Fill(strip_rec+count);}
	if(recHit_roll==2){histContainer_["Occupncy_RECO_roll2"]->Fill(strip_rec+count);}
	if(recHit_roll==3){histContainer_["Occupncy_RECO_roll3"]->Fill(strip_rec+count);}
	if(recHit_roll==4){histContainer_["Occupncy_RECO_roll4"]->Fill(strip_rec+count);}
	if(recHit_roll==5){histContainer_["Occupncy_RECO_roll5"]->Fill(strip_rec+count);}
	if(recHit_roll==6){histContainer_["Occupncy_RECO_roll6"]->Fill(strip_rec+count);}
	if(recHit_roll==7){histContainer_["Occupncy_RECO_roll7"]->Fill(strip_rec+count);}
	if(recHit_roll==8){histContainer_["Occupncy_RECO_roll8"]->Fill(strip_rec+count);}}
	
	std::cout<<" strip rec "<<strip_rec<<" RecHit "<<" x: "<<x_reco<<" x_err "<<err_x_reco<<" y "<<y_reco<<" cluster size "<<clusterSize<<" bx "<<bx_reco<<std::endl;
 
	if(recHit_roll==1){NrecHit_x_roll1 ++;if (bx_reco==0){NrecHit_x_roll1_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll1_noise ++;}}
	if(recHit_roll==2){NrecHit_x_roll2 ++;if (bx_reco==0){NrecHit_x_roll2_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll2_noise ++;}}
	if(recHit_roll==3){NrecHit_x_roll3 ++;if (bx_reco==0){NrecHit_x_roll3_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll3_noise ++;}}
	if(recHit_roll==4){NrecHit_x_roll4 ++;if (bx_reco==0){NrecHit_x_roll4_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll4_noise ++;}}
	if(recHit_roll==5){NrecHit_x_roll5 ++;if (bx_reco==0){NrecHit_x_roll5_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll5_noise ++;}}
	if(recHit_roll==6){NrecHit_x_roll6 ++;if (bx_reco==0){NrecHit_x_roll6_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll6_noise ++;}}
	if(recHit_roll==7){NrecHit_x_roll7 ++;if (bx_reco==0){NrecHit_x_roll7_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll7_noise ++;}}
	if(recHit_roll==8){NrecHit_x_roll8 ++;if (bx_reco==0){NrecHit_x_roll8_bx0 ++;} if (GEMsimHits->size()==0 && Ndigis!=0){NrecHit_x_roll8_noise ++;}} 	
/*	if(recHit_roll==1){NrecHit_x_roll1 ++;if (bx_reco==0){NrecHit_x_roll1_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll1_noise ++;}}
	if(recHit_roll==2){NrecHit_x_roll2 ++;if (bx_reco==0){NrecHit_x_roll2_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll2_noise ++;}}
	if(recHit_roll==3){NrecHit_x_roll3 ++;if (bx_reco==0){NrecHit_x_roll3_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll3_noise ++;}}
	if(recHit_roll==4){NrecHit_x_roll4 ++;if (bx_reco==0){NrecHit_x_roll4_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll4_noise ++;}}
	if(recHit_roll==5){NrecHit_x_roll5 ++;if (bx_reco==0){NrecHit_x_roll5_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll5_noise ++;}}
	if(recHit_roll==6){NrecHit_x_roll6 ++;if (bx_reco==0){NrecHit_x_roll6_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll6_noise ++;}}
	if(recHit_roll==7){NrecHit_x_roll7 ++;if (bx_reco==0){NrecHit_x_roll7_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll7_noise ++;}}
	if(recHit_roll==8){NrecHit_x_roll8 ++;if (bx_reco==0){NrecHit_x_roll8_bx0 ++;} if (bx_reco!=0){NrecHit_x_roll8_noise ++;}}*/
 }
 std::cout<<" NrecHit_x_roll1 "<<NrecHit_x_roll1<< " NrecHit_x_roll2 "<<NrecHit_x_roll2<<" NrecHit_x_roll3 "<<NrecHit_x_roll3<<" NrecHit_x_roll4 "<<NrecHit_x_roll4<<" NrecHit_x_roll5 "<<NrecHit_x_roll5<<" NrecHit_x_roll6 "<<NrecHit_x_roll6<<" NrecHit_x_roll7 "<<NrecHit_x_roll7<<" NrecHit_x_roll8 "<<NrecHit_x_roll8<<std::endl;
	if(NrecHit_x_roll1>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll1);
		if(NrecHit_x_roll1_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll1_bx0);
	 	if(NrecHit_x_roll1_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll1_noise);
		histContainer_["N_Cluster_roll1"]->Fill(NrecHit_x_roll1);}
	if(NrecHit_x_roll2>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll2);
		if(NrecHit_x_roll2_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll2_bx0);
	 	if(NrecHit_x_roll2_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll2_noise);
		histContainer_["N_Cluster_roll2"]->Fill(NrecHit_x_roll2);}
	if(NrecHit_x_roll3>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll3);
		if(NrecHit_x_roll3_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll3_bx0);
	 	if(NrecHit_x_roll3_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll3_noise);
		histContainer_["N_Cluster_roll3"]->Fill(NrecHit_x_roll3);}
	if(NrecHit_x_roll4>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll4);
		if(NrecHit_x_roll4_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll4_bx0);
	 	if(NrecHit_x_roll4_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll4_noise);
		histContainer_["N_Cluster_roll4"]->Fill(NrecHit_x_roll4);}
	if(NrecHit_x_roll5>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll5);
		if(NrecHit_x_roll5_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll5_bx0);
	 	if(NrecHit_x_roll5_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll5_noise);
		histContainer_["N_Cluster_roll5"]->Fill(NrecHit_x_roll5);}
	if(NrecHit_x_roll6>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll6);
		if(NrecHit_x_roll6_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll6_bx0);
	 	if(NrecHit_x_roll6_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll6_noise);
		histContainer_["N_Cluster_roll6"]->Fill(NrecHit_x_roll6);}
	if(NrecHit_x_roll7>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll7);
		if(NrecHit_x_roll7_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll7_bx0);
	 	if(NrecHit_x_roll7_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll7_noise);
		histContainer_["N_Cluster_roll7"]->Fill(NrecHit_x_roll7);}
	if(NrecHit_x_roll8>0){
		histContainer_["N_Cluster"]->Fill(NrecHit_x_roll8);
		if(NrecHit_x_roll8_bx0>0)histContainer_["N_Cluster_bx0"]->Fill(NrecHit_x_roll8_bx0);
	 	if(NrecHit_x_roll8_noise>0)histContainer_["N_Cluster_noise"]->Fill(NrecHit_x_roll8_noise);
		histContainer_["N_Cluster_roll8"]->Fill(NrecHit_x_roll8);}

  std::cout<<"n of SimHits: "<<GEMsimHits->size()<<std::endl;
  std::cout<<"n of Digis: "<<Ndigis<<std::endl;
  std::cout<<"n of RecHits: "<<GEMrecHits->size()<<std::endl;
 


// Reco-Sim

for (GEMRecHitCollection::const_iterator recHit = GEMrecHits->begin(); recHit != GEMrecHits->end(); ++recHit)
  {
// local reco
    	double x_reco = recHit->localPosition().x();
    	double err_x_reco = recHit->localPositionError().xx();
	double bx_reco = (Short_t) recHit->BunchX();
 	////   float y_reco = recHit->localPosition().y();
    	unsigned int RecodetId = (unsigned int) (*recHit).gemId();
    	//gem_recHit_.bx = recHit->BunchX();
    	int clusterSize = recHit->clusterSize();
	//double strip_rec = recHit->firstClusterStrip();

    	//GEMDetId id((*recHit).gemId());
 	////std::cout<<" x: "<<x_reco<<" +- "<<err_x_reco<<std::endl;
	 ////std::cout<<" y: "<<y_reco<<std::endl;

// reco Global
	GEMDetId rollId = (GEMDetId)(*recHit).gemId();
	LocalPoint recHitPos=recHit->localPosition();
	const GEMEtaPartition* rollasociated = GEMgeo->etaPartition(rollId);
	const BoundPlane & GEMSurface = rollasociated->surface(); 
	GlobalPoint GEMGlobalPoint = GEMSurface.toGlobal(recHitPos);
	//double X_reco_global = GEMGlobalPoint.x();
	//double Y_reco_global = GEMGlobalPoint.y();
	double recPhi = GEMGlobalPoint.phi();
	double recEta = GEMGlobalPoint.eta();
//std::cout<<" global x "<<x_global<<std::endl; 

	int recHit_roll = rollId.roll();
	histContainer_["CS"]->Fill(clusterSize);
	if (bx_reco==0){histContainer_["CS_bx0"]->Fill(clusterSize);}
	if(bx_reco!=0){histContainer_["CS_noise"]->Fill(clusterSize);}
 	
	if(recHit_roll==1){histContainer_["CS_roll1"]->Fill(clusterSize);}
	if(recHit_roll==2){histContainer_["CS_roll2"]->Fill(clusterSize);}
	if(recHit_roll==3){histContainer_["CS_roll3"]->Fill(clusterSize);}
	if(recHit_roll==4){histContainer_["CS_roll4"]->Fill(clusterSize);}
	if(recHit_roll==5){histContainer_["CS_roll5"]->Fill(clusterSize);}
	if(recHit_roll==6){histContainer_["CS_roll6"]->Fill(clusterSize);}
	if(recHit_roll==7){histContainer_["CS_roll7"]->Fill(clusterSize);}
	if(recHit_roll==8){histContainer_["CS_roll8"]->Fill(clusterSize);}

// solo bx 0
	if (bx_reco!=0) continue;	
  	for (edm::PSimHitContainer::const_iterator itHit = GEMsimHits->begin(); itHit != GEMsimHits->end(); ++itHit)
  	{
// sim local
	

		//GEMDetId idGem = GEMDetId(itHit->detUnitId());
		DetId id = DetId(itHit->detUnitId());
		//LocalPoint lp = itHit->entryPoint();
		//int strip_sim = GEMgeo->etaPartition(idGem)->strip(lp);
		//std::cout<<" sim strip "<<strip_sim<<std::endl;
		
    		double x_sim = itHit->localPosition().x();
    		//float err_x_sim = itHit->localPositionError().xx();
    		double y_sim = itHit->localPosition().y();
		
// sim Global

		GlobalPoint pointSimHit = theTrackingGeometry->idToDetUnit(id)->toGlobal(itHit->localPosition());
		//float X_sim_global = pointSimHit.x();
		//float Y_sim_global = pointSimHit.y();
		float simPhi = pointSimHit.phi();
		float simEta = pointSimHit.eta();
  		//std::cout<<"sim X "<<simX<<std::endl;
		double dPhi_sim_reco = simPhi-recPhi;
		double dR_sim_reco = deltaR(simEta, simPhi, recEta, recPhi);
  		//std::cout<<"dR "<<dR_sim_reco<<std::endl;
		
  		////std::cout<<"sim detId "<<itHit->detUnitId()<<std::endl;
  		////std::cout<<"rec detId "<<detId<<std::endl;

		
	int sim_x_roll = 0;
	for (edm::PSimHitContainer::const_iterator itHit = GEMsimHits->begin(); itHit != GEMsimHits->end(); ++itHit){
	if (itHit->detUnitId() == RecodetId){
	//if (std::abs(itHit->particleType())==13){
	sim_x_roll ++;
	//}
	}
	}

			
  		if (itHit->detUnitId() == RecodetId){// sim rec hit same rawID
			histContainer_["dR_sim_reco_allP"]->Fill(dR_sim_reco);
			histContainer_["dPhi_sim_reco_allP"]->Fill(dPhi_sim_reco);
			histContainer2D_["dPhiVScs_allP"]->Fill(clusterSize, dPhi_sim_reco);

			if(recHit_roll==1 && NrecHit_x_roll1>0){
				histContainer2D_["dPhiVSNc1"]->Fill(NrecHit_x_roll1, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll1, dPhi_sim_reco);}
        		if(recHit_roll==2 && NrecHit_x_roll2>0){
				histContainer2D_["dPhiVSNc2"]->Fill(NrecHit_x_roll2, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll2, dPhi_sim_reco);}
        		if(recHit_roll==3 && NrecHit_x_roll3>0){
				histContainer2D_["dPhiVSNc3"]->Fill(NrecHit_x_roll3, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll3, dPhi_sim_reco);}
        		if(recHit_roll==4 && NrecHit_x_roll4>0){
				histContainer2D_["dPhiVSNc4"]->Fill(NrecHit_x_roll4, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll4, dPhi_sim_reco);}
        		if(recHit_roll==5 && NrecHit_x_roll5>0){
				histContainer2D_["dPhiVSNc5"]->Fill(NrecHit_x_roll5, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll5, dPhi_sim_reco);}
        		if(recHit_roll==6 && NrecHit_x_roll6>0){
				histContainer2D_["dPhiVSNc6"]->Fill(NrecHit_x_roll6, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll6, dPhi_sim_reco);}
        		if(recHit_roll==7 && NrecHit_x_roll7>0){
				histContainer2D_["dPhiVSNc7"]->Fill(NrecHit_x_roll7, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll7, dPhi_sim_reco);}
        		if(recHit_roll==8 && NrecHit_x_roll8>0){
				histContainer2D_["dPhiVSNc8"]->Fill(NrecHit_x_roll8, dPhi_sim_reco);
				histContainer2D_["dPhiVSNc"]->Fill(NrecHit_x_roll8, dPhi_sim_reco);}

		if (std::abs(itHit->particleType())==13){ // se sim hit è muone
			histContainer_["dR_sim_reco_Mu"]->Fill(dR_sim_reco);
			histContainer_["dPhi_sim_reco_Mu"]->Fill(dPhi_sim_reco);
			histContainer2D_["dPhiVScs_Mu"]->Fill(clusterSize, dPhi_sim_reco);

		if (sim_x_roll ==1){

			histContainer_["dPhi_sim_reco_Mu_1simHxroll"]->Fill(dPhi_sim_reco);
			histContainer2D_["dPhiVScs_Mu_1simHxroll"]->Fill(clusterSize, dPhi_sim_reco);

				if(recHit->clusterSize()== 1){histContainer_["dPhi_sim_reco_Mu_1simHxroll_cl1"]->Fill(dPhi_sim_reco);}
		//if(dR_sim_reco > 0.03) continue;// delta R sim rec hits
	std::cout<<"++++++++++++++++"<<std::endl;
    			double dX = x_sim - x_reco;
    			double pull = dX/std::sqrt(err_x_reco);
	std::cout<<" x_sim: "<<x_sim<<" x_reco: "<<x_reco<<" +- "<<err_x_reco<<std::endl;
	std::cout<<" phi_sim: "<<simPhi<<" phi_reco: "<<recPhi<<std::endl;
			 //std::cout<<" x_sim: "<<x_sim<<std::endl;
 			/*std::cout<<" y_reco: "<<y_reco<<std::endl;
			 std::cout<<" y_sim: "<<y_sim<<std::endl;
 			std::cout<<" dX: "<<dX<<std::endl;
			 std::cout<<" pull: "<<pull<<std::endl;*/
		
			 
			 histContainer_["DeltaX"]->Fill(dX);
			 histContainer_["RecoLocalX"]->Fill(x_reco);
			 histContainer_["RecoLocalErrorX"]->Fill(err_x_reco);
			 histContainer_["DeltaX_over_err"]->Fill(pull);
		       	histContainer2D_["hGEM_simVSrec"]->Fill(x_sim, x_reco);
		       	histContainer2D_["hGEM_sim_yVSdx"]->Fill(dX, y_sim);
			
		}// 1 mu per roll
     		}//end muon
     		}//end same det id

   	}// end sim hit
  }// end rec hit
  

//int N_sim=GEMsimHits->size();
//if(N_sim==0 && Ndigis!=0){
for (GEMRecHitCollection::const_iterator recHit = GEMrecHits->begin(); recHit != GEMrecHits->end(); ++recHit)
  {	bool match=false;
  	unsigned int RecodetId = (unsigned int) (*recHit).gemId();
	GEMDetId rollId = (GEMDetId)(*recHit).gemId();
	int recHit_roll = rollId.roll();
  	int clusterSize = recHit->clusterSize();
	//double strip_rec = recHit->firstClusterStrip();
	double bx_reco = (Short_t) recHit->BunchX();

  	for (edm::PSimHitContainer::const_iterator itHit = GEMsimHits->begin(); itHit != GEMsimHits->end(); ++itHit)
  	{
		/*GEMDetId idGem = GEMDetId(itHit->detUnitId());
		LocalPoint lp = itHit->entryPoint();
		int strip_sim = GEMgeo->etaPartition(idGem)->strip(lp);*/
	if (itHit->detUnitId() == RecodetId){
		//if((strip_sim >= (strip_rec -1)-5)&&(strip_sim <= ((strip_rec -1)+(clusterSize-1)+5))){
			match=true;
		//}
	}
	}// loom sim
		
	if (match==false){
		histContainer_["CS_NOISE"]->Fill(clusterSize);
		histContainer_["N_RecHit_noise"]->Fill(1);
		histContainer_["bx_noise"]->Fill(bx_reco);
	 		if(recHit_roll==1){histContainer_["N_RecHit_noise1"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise1"]->GetEntries());}
			if(recHit_roll==2){histContainer_["N_RecHit_noise2"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise2"]->GetEntries());}
			if(recHit_roll==3){histContainer_["N_RecHit_noise3"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise3"]->GetEntries());}
			if(recHit_roll==4){histContainer_["N_RecHit_noise4"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise4"]->GetEntries());}
			if(recHit_roll==5){histContainer_["N_RecHit_noise5"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise5"]->GetEntries());}
			if(recHit_roll==6){histContainer_["N_RecHit_noise6"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise6"]->GetEntries());}
			if(recHit_roll==7){histContainer_["N_RecHit_noise7"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise7"]->GetEntries());}
			if(recHit_roll==8){histContainer_["N_RecHit_noise8"]->Fill(1);
			histContainer_["Noise_rate"]->SetBinContent(recHit_roll , histContainer_["N_RecHit_noise8"]->GetEntries());}
	 		}
   } //loop rec
 //}// if
}// analyze




/*#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif*/



// ------------ method called once each job just before starting event loop  ------------
void 
GEMLocal::beginJob()
{ edm::Service<TFileService> fs;
 histContainer_["N_eventi"]=fs->make<TH1F>("N_eventi", "count",    3,   -0.5, 2.5);
 histContainer_["N_SimHit"]=fs->make<TH1F>("N_SimHit", "count",    3,   -0.5, 2.5);
 histContainer_["N_strips"]=fs->make<TH1F>("N_strips", "count strips",    3,   -0.5, 2.5);
 histContainer_["N_Digis"]=fs->make<TH1F>("N_Digi", "count digis",    3,   -0.5, 2.5);
 histContainer_["N_Digis_bx0"]=fs->make<TH1F>("N_Digi_bx0", "count digis bx=0",    3,   -0.5, 2.5);
 histContainer_["N_Digis_noise"]=fs->make<TH1F>("N_Digi_noise", "count digis bx!=0",    3,   -0.5, 2.5);
 histContainer_["N_RecHit"]=fs->make<TH1F>("N_RecHit", "count recHits",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_bx0"]=fs->make<TH1F>("N_RecHit_bx0", "count recHits bx=0",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise"]=fs->make<TH1F>("N_RecHit_noise", "count recHits bx!=0",    3,   -0.5, 2.5);
 
 histContainer_["Noise_rate"]=fs->make<TH1F>("Noise_rate", "count recHits noise x roll",    9,   0., 9.);
 histContainer_["N_RecHit_noise1"]=fs->make<TH1F>("N_RecHit_noise1", "count recHits bx!=0 roll 1",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise2"]=fs->make<TH1F>("N_RecHit_noise2", "count recHits bx!=0 roll 2",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise3"]=fs->make<TH1F>("N_RecHit_noise3", "count recHits bx!=0 roll 3",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise4"]=fs->make<TH1F>("N_RecHit_noise4", "count recHits bx!=0 roll 4",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise5"]=fs->make<TH1F>("N_RecHit_noise5", "count recHits bx!=0 roll 5",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise6"]=fs->make<TH1F>("N_RecHit_noise6", "count recHits bx!=0 roll 6",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise7"]=fs->make<TH1F>("N_RecHit_noise7", "count recHits bx!=0 roll 7",    3,   -0.5, 2.5);
 histContainer_["N_RecHit_noise8"]=fs->make<TH1F>("N_RecHit_noise8", "count recHits bx!=0 roll 8",    3,   -0.5, 2.5);
 
 histContainer2D_["Occupancy_SIM_partitionVSstrip"] = fs->make<TH2F>("Occupancy_SIM_partitionVSstrip", "Partitions vs strips", 400, 0., 400., 16, 0., 8.);
 histContainer2D_["Occupancy_DIGI_partitionVSstrip"] = fs->make<TH2F>("Occupancy_DIGI_partitionVSstrip", "Partitions vs strips", 400, 0., 400., 16, 0., 8.);
 histContainer2D_["Occupancy_RECO_partitionVSstrip"] = fs->make<TH2F>("Occupancy_RECO_partitionVSstrip", "Partitions vs strips", 400, 0., 400., 16, 0., 8.);
 
 histContainer_["Occupncy_SIM_roll1"]=fs->make<TH1F>("Occupncy_SIM_roll1", "fired strip roll 1",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll2"]=fs->make<TH1F>("Occupncy_SIM_roll2", "fired strip roll 2",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll3"]=fs->make<TH1F>("Occupncy_SIM_roll3", "fired strip roll 3",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll4"]=fs->make<TH1F>("Occupncy_SIM_roll4", "fired strip roll 4",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll5"]=fs->make<TH1F>("Occupncy_SIM_roll5", "fired strip roll 5",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll6"]=fs->make<TH1F>("Occupncy_SIM_roll6", "fired strip roll 6",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll7"]=fs->make<TH1F>("Occupncy_SIM_roll7", "fired strip roll 7",    400,   0, 400);
 histContainer_["Occupncy_SIM_roll8"]=fs->make<TH1F>("Occupncy_SIM_roll8", "fired strip roll 8",    400,   0, 400);
 
 histContainer_["Occupncy_RECO_roll1"]=fs->make<TH1F>("Occupncy_RECO_roll1", "reco fired strip roll 1",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll2"]=fs->make<TH1F>("Occupncy_RECO_roll2", "reco fired strip roll 2",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll3"]=fs->make<TH1F>("Occupncy_RECO_roll3", "reco fired strip roll 3",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll4"]=fs->make<TH1F>("Occupncy_RECO_roll4", "reco fired strip roll 4",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll5"]=fs->make<TH1F>("Occupncy_RECO_roll5", "reco fired strip roll 5",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll6"]=fs->make<TH1F>("Occupncy_RECO_roll6", "reco fired strip roll 6",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll7"]=fs->make<TH1F>("Occupncy_RECO_roll7", "reco fired strip roll 7",    400,   0, 400);
 histContainer_["Occupncy_RECO_roll8"]=fs->make<TH1F>("Occupncy_RECO_roll8", "reco fired strip roll 8",    400,   0, 400);
 
 histContainer_["bx"]=fs->make<TH1F>("bx", "bx",    12,   -6, 6);
 histContainer_["bx_noise"]=fs->make<TH1F>("bx_noise", "bx_noise",    12,   -6, 6);
 histContainer_["bx_noise_bx!=0"]=fs->make<TH1F>("bx_noise_bx!=0", "bx_noise bx!=0",    12,   -6, 6);
 
histContainer_["CS"]=fs->make<TH1F>("CS", "cluster size",    11,   -0.5, 10.5);
histContainer_["CS_NOISE"]=fs->make<TH1F>("CS_NOISE", "cluster size NOISE",    11,   -0.5, 10.5);
histContainer_["CS_bx0"]=fs->make<TH1F>("CS_bx0", "cluster size bx=0",    11,   -0.5, 10.5);
histContainer_["CS_noise"]=fs->make<TH1F>("CS_noise", "cluster size bx!=0",    11,   -0.5, 10.5);
histContainer_["CS_roll1"]=fs->make<TH1F>("CS_roll1", "cluster size roll 1",    11,   -0.5, 10.5);
histContainer_["CS_roll2"]=fs->make<TH1F>("CS_roll2", "cluster size roll 2",    11,   -0.5, 10.5);
histContainer_["CS_roll3"]=fs->make<TH1F>("CS_roll3", "cluster size roll 3",    11,   -0.5, 10.5);
histContainer_["CS_roll4"]=fs->make<TH1F>("CS_roll4", "cluster size roll 4",    11,   -0.5, 10.5);
histContainer_["CS_roll5"]=fs->make<TH1F>("CS_roll5", "cluster size roll 5",    11,   -0.5, 10.5);
histContainer_["CS_roll6"]=fs->make<TH1F>("CS_roll6", "cluster size roll 6",    11,   -0.5, 10.5);
histContainer_["CS_roll7"]=fs->make<TH1F>("CS_roll7", "cluster size roll 7",    11,   -0.5, 10.5);
histContainer_["CS_roll8"]=fs->make<TH1F>("CS_roll8", "cluster size roll 8",    11,   -0.5, 10.5);	 

histContainer_["N_Cluster"]=fs->make<TH1F>("N_Cluster", " N cluster",    8,   -0.5, 7.5);
histContainer_["N_Cluster_bx0"]=fs->make<TH1F>("N_Cluster_bx0", " N cluster bx=0",    8,   -0.5, 7.5);
histContainer_["N_Cluster_noise"]=fs->make<TH1F>("N_Cluster_noise", " N cluster bx!=0",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll1"]=fs->make<TH1F>("N_Cluster_roll1", " N cluster x roll 1",    8,   -0.5, 7.5);	 
histContainer_["N_Cluster_roll2"]=fs->make<TH1F>("N_Cluster_roll2", " N cluster x roll 2",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll3"]=fs->make<TH1F>("N_Cluster_roll3", " N cluster x roll 3",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll4"]=fs->make<TH1F>("N_Cluster_roll4", " N cluster x roll 4",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll5"]=fs->make<TH1F>("N_Cluster_roll5", " N cluster x roll 5",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll6"]=fs->make<TH1F>("N_Cluster_roll6", " N cluster x roll 6",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll7"]=fs->make<TH1F>("N_Cluster_roll7", " N cluster x roll 7",    8,   -0.5, 7.5);
histContainer_["N_Cluster_roll8"]=fs->make<TH1F>("N_Cluster_roll8", " N cluster x roll 8",    8,   -0.5, 7.5);	 

histContainer2D_["dPhiVSNc"]=fs->make<TH2F>("dPhiVSNc", " N cluster x roll  vs dPhi",    12,   0., 6., 500,   -0.001, 0.001);
	 
histContainer2D_["dPhiVSNc1"]=fs->make<TH2F>("dPhiVSNc1", " N cluster x roll 1 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);	 
histContainer2D_["dPhiVSNc2"]=fs->make<TH2F>("dPhiVSNc2", " N cluster x roll 2 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);
histContainer2D_["dPhiVSNc3"]=fs->make<TH2F>("dPhiVSNc3", " N cluster x roll 3 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);
histContainer2D_["dPhiVSNc4"]=fs->make<TH2F>("dPhiVSNc4", " N cluster x roll 4 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);
histContainer2D_["dPhiVSNc5"]=fs->make<TH2F>("dPhiVSNc5", " N cluster x roll 5 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);
histContainer2D_["dPhiVSNc6"]=fs->make<TH2F>("dPhiVSNc6", " N cluster x roll 6 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);
histContainer2D_["dPhiVSNc7"]=fs->make<TH2F>("dPhiVSNc7", " N cluster x roll 7 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);
histContainer2D_["dPhiVSNc8"]=fs->make<TH2F>("dPhiVSNc8", " N cluster x roll 8 vs dPhi",    6,   -0.5, 5.5, 500,   -0.001, 0.001);		 
	 
histContainer_["dR_sim_reco_allP"]=fs->make<TH1F>("dR_sim_reco_allP", "dR_sim_reco allParticles (ele,mu)",    120,   0., 0.12);
histContainer_["dR_sim_reco_Mu"]=fs->make<TH1F>("dR_sim_reco_Mu", "dR_sim_reco Muon only ",    120,   0., 0.12);

histContainer_["dPhi_sim_reco_allP"]=fs->make<TH1F>("dPhi_sim_reco_allP", "dPhi_sim_reco allParticles (ele,mu) ",    500,   -0.001, 0.001);
histContainer_["dPhi_sim_reco_Mu"]=fs->make<TH1F>("dPhi_sim_reco_Mu", "dPhi_sim_reco Muon only ",    500,   -0.001, 0.001);
histContainer_["dPhi_sim_reco_Mu_1simHxroll"]=fs->make<TH1F>("dPhi_sim_reco_Mu_1simHxroll", "dPhi_sim_reco Muon only 1 simHit per roll",    500,   -0.001, 0.001);
histContainer_["dPhi_sim_reco_Mu_1simHxroll_cl1"]=fs->make<TH1F>("dPhi_sim_reco_Mu_1simHxroll_cl1", "dPhi_sim_reco Muon only 1 simHit per roll clusterSize=1 ",    500,   -0.001, 0.001);
 
 histContainer_["RecoLocalX"]=fs->make<TH1F>("RecoLocalX", "x_rec",    400,   -20., 20.);
 histContainer_["RecoLocalErrorX"]=fs->make<TH1F>("RecoLocalErrorX", "xmc - x_rec",    200,   0., 0.5);
 histContainer_["DeltaX"]=fs->make<TH1F>("DeltaX", "xmc - x_rec",    50,   -2., 2.);
 histContainer_["DeltaX_over_err"]=fs->make<TH1F>("DeltaX_over_err", "xmc - x_rec / sigma_err",    500,   -15., 15.);

 histContainer2D_["hGEM_simVSrec"] = fs->make<TH2F>("hGEM_simVSrec", "x_sim vs x_rec", 400, -20, 20, 400, -20, 20);
 histContainer2D_["hGEM_sim_yVSdx"] = fs->make<TH2F>("hGEM_sim_yVSdx", "dX vs y_sim", 200, -1, 1, 500, -10, 10);
 histContainer2D_["dPhiVScs_Mu_1simHxroll"] = fs->make<TH2F>("dPhiVScs_Mu_1simHxroll", "dPhi vs  recHit cluster size mu 1 simHit x roll", 12, 0., 6., 500, -0.001, 0.001);
 histContainer2D_["dPhiVScs_Mu"] = fs->make<TH2F>("dPhiVScs_Mu", "dPhi vs  recHit cluster size mu", 12, 0., 6., 500, -0.001, 0.001);
 histContainer2D_["dPhiVScs_allP"] = fs->make<TH2F>("dPhiVScs_allP", "dPhi vs  recHit cluster size ele + mu", 12, 0., 6., 500, -0.001, 0.001);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GEMLocal::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
GEMLocal::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GEMLocal::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GEMLocal::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GEMLocal::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GEMLocal::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMLocal);
