#define SNIPER_VERSION_2 1
#include "BiPo212_reader.h"
#include "JUNO_PMTs.h"
#include "Identifier/IDService.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "EvtNavigator/NavBuffer.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperLog.h"
#include "Event/SimHeader.h"
#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdLpmtCalibEvt.h"
#include "Event/CdTriggerHeader.h"
#include "Event/CdTriggerEvt.h"
#include "Event/WpCalibHeader.h"
#include "Event/WpCalibEvt.h"
#include "Event/WpTriggerHeader.h"
#include "Event/WpTriggerEvt.h"
#include "Event/CdVertexRecHeader.h"
#include "Event/CdVertexRecEvt.h"
#include "RootWriter/RootWriter.h"
#include <numeric>

#include "TH1F.h"
#include "TTree.h"

int nBins = 200;

DECLARE_ALGORITHM(BiPo212_reader);

BiPo212_reader::BiPo212_reader(const std::string& name) 
	: AlgBase(name),
	  m_iEvt(0),
	  m_buf(0)
{
}

bool BiPo212_reader::initialize() {

    LogDebug << "initializing" << std::endl;
    auto toptask = getRoot();

    idServ = IDService::getIdServ();
    idServ->init();
	
    SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");

    if ( navBuf.invalid() ) {
        LogError << "cannot get the NavBuffer @ /Event" << std::endl;
        return false;
    }
	
    m_buf = navBuf.data();

    SniperPtr<RootWriter> rw(getParent(), "RootWriter");
    if (rw.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enable it in your job option file."
                 << std::endl;
         return false;
    }

    //wp_events_seen = 0;    
    events = rw->bookTree(*m_par,"tree/CdEvents","Events Tree");
    //events->Branch("EvtID",&cdEvtID,"EvtID/I");
    events->Branch("TimeStamp",&timestamp);
    //events->Branch("PMTID",&PMTID);
    //events->Branch("Charge",&charge);
    events -> Branch("ChargeCenter",&ChargeCenter);
    events->Branch("Time",&time);
    events->Branch("NPE",&total_npe,"npe/F");
    events->Branch("TriggerType",&trigger_type);
    //events->Branch("Recox",&CdRecox);
    //events->Branch("Recoy",&CdRecoy);
    //events->Branch("Recoz",&CdRecoz);
    //events->Branch("Recopx",&CdRecopx);
    //events->Branch("Recopy",&CdRecopy);
    //events->Branch("Recopz",&CdRecopz);
    //events->Branch("RecoEnergy",&CdRecoenergy);
    //events->Branch("RecoPE",&CdRecoPESum);
    //events->Branch("Recot0",&CdRecot0);
    //events->Branch("RecoEnergyQuality",&CdRecoEnergyQuality);
    //events->Branch("RecoPositionQuality",&CdRecoPositionQuality);
/*
    events->Branch("x_CM",&x_CM,"x_CM/F");
    events->Branch("y_CM",&y_CM,"y_CM/F");
    events->Branch("z_CM",&z_CM,"z_CM/F");

    wpEvents = rw->bookTree(*m_par,"tree/WpEvents","Water pool Events Tree");
    wpEvents->Branch("EvtID",&wpEvtID,"EvtID/I");
    wpEvents->Branch("TimeStamp",&wptimestamp);
    wpEvents->Branch("PMTID",&wpPMTID);
    wpEvents->Branch("Charge",&wpcharge);
    wpEvents->Branch("Time",&wptime);
    wpEvents->Branch("NPE",&wptotal_npe,"npe/F");
    wpEvents->Branch("TriggerType",&wptrigger_type);
*/
    kernel = new double[nBins];
    PMTs_Pos.SetCdPmts("/cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/Jlatest/data/Detector/Identifier/pmt_CDLPMT_latest.csv");
    return true;
}

bool BiPo212_reader::execute() {

	LogDebug << "=====================================" << std::endl;
	LogDebug << "executing: " << m_iEvt++
             << std::endl;

	JM::CdLpmtCalibEvt* calibevent = 0;
	JM::CdTriggerEvt* triggerevent = 0;
	JM::CdVertexRecEvt* recoevent = 0;

	auto nav = m_buf->curEvt();
//	const auto& paths = nav->getPath();

	auto calibheader = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	if (calibheader) calibevent = (JM::CdLpmtCalibEvt*)calibheader->event();

	auto recoheader = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
	if (recoheader) recoevent = (JM::CdVertexRecEvt*)calibheader->event();
	

	if (!calibevent || !recoevent) {
		LogInfo << "No CalibEvt or RecEvt found, skipping..." << std::endl;
		return true;
	}
	if (calibevent) { //&& recoevent) {

		auto triggerheader = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
		if (triggerheader) triggerevent = (JM::CdTriggerEvt*) triggerheader -> event();

		if (!triggerevent) {
			LogInfo << "No CdTriggerEvt found for this event " << std::endl;
			trigger_type = "None";
		} else {
			const auto& triggerelement = triggerevent -> triggerType();
			if (triggerelement.size() != 1) LogInfo << "------------- Strange trigger element! -------------" << std::endl;
			else trigger_type = triggerelement[0];
		}

	
        	charge.clear();
        	time.clear();
        	PMTID.clear();

		float AccX=0.0,AccY=0.0,AccZ=0.0;
		int charge_in_range = 0;
	
	        for (const auto& element : calibevent->calibPMTCol()) {
		
        	    //PMTID.push_back(element->pmtId());
			auto Times = element -> time();
			int i = 0;
			for (auto pmtChannel : element -> charge() ) {
				charge.push_back(pmtChannel);
				int PmtNo = idServ->id2CopyNo(Identifier(element->pmtId()));
				PMTID.push_back(PmtNo);
				if (Times[i] > 200 && Times[i] < 450) {
					AccX += pmtChannel*PMTs_Pos.GetX(PmtNo);
					AccY += pmtChannel*PMTs_Pos.GetY(PmtNo);
					AccZ += pmtChannel*PMTs_Pos.GetZ(PmtNo);
					charge_in_range += pmtChannel;	
				}
				i++;
			}
			for (auto pmtChannel : element -> time() ) {
				time.push_back(pmtChannel);
			}

        	}

        	total_npe = std::accumulate(charge.begin(),charge.end(),0.0);

		if (total_npe < 6600 && total_npe > 2600) {
			

		}

		/*ChargeCenter = std::make_tuple(AccX/charge_in_range,AccY/charge_in_range,AccZ/charge_in_range);

        	timestamp = nav->TimeStamp();
		const auto& recovertices = recoevent -> vertices();

                if (recovertices.size() != 1) LogInfo << "--------------- Multiple verteces!!!!!!!! -----------------" << endl;
                int ii=0;
                if (recovertices.size() == 2) ii=1;
                CdRecox = recovertices[ii] -> x();
                CdRecoy = recovertices[ii] -> y();
                CdRecoz = recovertices[ii] -> z();
                //CdRecopx = recovertices[ii] -> px();
                //CdRecopy = recovertices[ii] -> py();
                //CdRecopz = recovertices[ii] -> pz();
                //CdRecoenergy = recovertices[ii] -> energy();
                //CdRecoPESum = recovertices[ii] -> peSum();
                //CdRecot0 = recovertices[ii] -> t0();
                //CdRecoPositionQuality = recovertices[ii] -> positionQuality();
                //CdRecoEnergyQuality = recovertices[ii] -> energyQuality();
		events -> Fill();
		//if (total_npe < 6600 && total_npe > 2600) events -> Fill();
		//if ( total_npe > 2600 && total_npe < 3000 ) events -> Fill();
		//else return true;
		*/
    		}
    return true;

}

bool BiPo212_reader::finalize() {

    return true;
    
}
