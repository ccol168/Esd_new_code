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
#include "Event/OecHeader.h"
#include "Event/OecEvt.h"
#include <numeric>
#include <TSpectrum.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include <TGraph.h>

#include "TH1F.h"
#include "TTree.h"
#include "TParameter.h"

int BinsNumber = 200;

std::vector<int> MakeHistogram (std::vector <float> inVec) {

	double MinBinTime = 0.;
	double MaxBinTime = 1200.;
	double BinSize = (MaxBinTime - MinBinTime) / BinsNumber;
	std::vector <int> Histogram(BinsNumber,0);

	for (auto element : inVec) {
		int BinIndex = (element - MinBinTime) / BinSize ;
		Histogram[BinIndex]++;
	}

	return Histogram;
}	

std::vector <double> GetKernel (const char* filename) {
  TFile* file = new TFile(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error opening file" << std::endl;
    return {};
  }
  TH1F* hist = nullptr;
  file->GetObject("histo", hist);

  if (!hist) {
    std::cout << "Histogram named histo not found" << std::endl;
    file->Close();
    return {};
  }

  std::vector<double> entries;
  int nBins = hist->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
      entries.push_back(hist->GetBinContent(i));
  }

  file->Close();
  return entries;

}

double distance (float x, float y, float z, float x1, float y1, float z1) {
	return sqrt(pow(x-x1,2)+pow(y-y1,2)+pow(z-z1,2));
}

double calculate_ToF (float x_int, float y_int, float z_int, float x_PMT, float y_PMT, float z_PMT) {

    float n_water = 1.29;
    float n_scint = 1.54;
    float c = 299.792 ; // mm/ns

    return distance(x_int,y_int,z_int,x_PMT,y_PMT,z_PMT)/c * n_scint;
}

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
    events->Branch("EvtID",&cdEvtID,"EvtID/I");
    events->Branch("TimeStamp",&timestamp);
    //events->Branch("PMTID",&PMTID);
    events->Branch("NPeaks",&NPeaks);
    events->Branch("PeakPositions",&peak_positions);
    events->Branch("Time",&time);
    events->Branch("CorrTime",&corr_time);
    events->Branch("NPE",&total_npe,"npe/F");
    events->Branch("TriggerType",&trigger_type);
    events->Branch("Recox",&CdRecox);
    events->Branch("Recoy",&CdRecoy);
    events->Branch("Recoz",&CdRecoz);
    events->Branch("DeconvolutedSignal",&DeconvolutedSignal);
    events->Branch("PEBi",&PEBi);
    events->Branch("PEPo",&PEPo);
	events->Branch("TimeSinceLastMuon",&TimeSinceLastMuon);
    //events->Branch("Recopx",&CdRecopx);
    //events->Branch("Recopy",&CdRecopy);
    //events->Branch("Recopz",&CdRecopz);
    events->Branch("RecoEnergy",&CdRecoenergy);
    //events->Branch("RecoPE",&CdRecoPESum);
    //events->Branch("Recot0",&CdRecot0);
    //events->Branch("RecoEnergyQuality",&CdRecoEnergyQuality);
    //events->Branch("RecoPositionQuality",&CdRecoPositionQuality);

	// Create summary tree for muon count
	summaryTree = rw->bookTree(*m_par, "tree/summary", "Summary Tree");
	summaryTree->Branch("nMuonsTotal", &nMuonsTotal, "nMuonsTotal/I");
	summaryTree->Branch("runLength", &runLength, "runLength/D");
	minEventTimestamp.Set(0, true, 0, false);
	maxEventTimestamp.Set(0, true, 0, false);
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
    std::vector <double> KernelVector = GetKernel("/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/BiPo212_analyzer/Kernel.root");
    kernel = new double[200];

    for (int i=0; i<200; i++) {
	    kernel[i] = KernelVector[i];
    }

    //spectrum = new TSpectrum();
    CorrTimesHistogramArr = new double[200];

    PMTs_Pos.SetCdPmts("/cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/Jlatest/data/Detector/Identifier/pmt_CDLPMT_latest.csv");
    last_muon_timestamp.Set(0,true,0,false);

    return true;
}

bool BiPo212_reader::execute() {

	LogInfo << "=====================================" << std::endl;
	LogInfo << "executing: " << m_iEvt++ << std::endl;

	JM::CdLpmtCalibEvt* calibevent = 0;
	JM::CdTriggerEvt* triggerevent = 0;
	//JM::CdVertexRecEvt* recoevent = 0;
	JM::OecEvt* oecevent = 0;

	auto nav = m_buf->curEvt();

	if (m_iEvt == 1) {
		LogDebug << "FirstTimeStamp = " << (nav->TimeStamp()).AsString() << endl;
	}

	auto oecheader = JM::getHeaderObject<JM::OecHeader>(nav);
	if (oecheader) oecevent = (JM::OecEvt*)oecheader->event("JM::OecEvt");

	auto calibheader = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
	if (calibheader) calibevent = (JM::CdLpmtCalibEvt*)calibheader->event();

	//auto recoheader = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
	//if (recoheader) recoevent = (JM::CdVertexRecEvt*)calibheader->event();
	timestamp = nav->TimeStamp();
	// Track min and max event timestamps robustly
	double ts = timestamp.GetSec() + timestamp.GetNanoSec() * 1e-9;
	double minTs = minEventTimestamp.GetSec() + minEventTimestamp.GetNanoSec() * 1e-9;
	double maxTs = maxEventTimestamp.GetSec() + maxEventTimestamp.GetNanoSec() * 1e-9;
	if (minEventTimestamp.GetSec() == 0 && minEventTimestamp.GetNanoSec() == 0) {
		minEventTimestamp = timestamp;
		maxEventTimestamp = timestamp;
	} else {
		if (ts < minTs) minEventTimestamp = timestamp;
		if (ts > maxTs) maxEventTimestamp = timestamp;
	}

	if (!calibevent || !oecevent) {
		LogInfo << "No CalibEvt or OecEvt found, skipping..." << std::endl;
		return true;
	}
	if (calibevent && oecevent) { //&& recoevent) {

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

		CdRecox = oecevent -> getVertexX();
                CdRecoy = oecevent -> getVertexY();
                CdRecoz =  oecevent -> getVertexZ();
		CdRecoenergy = oecevent -> getTotalCharge();
	
        	charge.clear();
        	time.clear();
        	PMTID.clear();
			corr_time.clear();

	        for (const auto& element : calibevent->calibPMTCol()) {
		
			for (auto pmtChannel : element -> charge() ) {
				charge.push_back(pmtChannel);
				int PmtNo = idServ->id2CopyNo(Identifier(element->pmtId()));
				PMTID.push_back(PmtNo);
			}
			for (auto pmtChannel : element -> time() ) {
				time.push_back(pmtChannel);
			}

        	}

			total_npe = std::accumulate(charge.begin(),charge.end(),0.0);

			if (total_npe > 30000) {
				last_muon_timestamp = timestamp;
				//LogDebug << "Found muon at " << timestamp.GetSec() << " " << timestamp.GetNanoSec() << endl;
				nMuons++; // Increment muon counter
				return true;
			}

			// Robust check for last_muon_timestamp
			int SecondDifference = 0;
			if (last_muon_timestamp.GetSec() == 0 && last_muon_timestamp.GetNanoSec() == 0) {
				TimeSinceLastMuon = -1;
			} else {
				SecondDifference = timestamp.GetSec() - last_muon_timestamp.GetSec();
				TimeSinceLastMuon = (SecondDifference*1.e9 + (timestamp.GetNanoSec() - last_muon_timestamp.GetNanoSec()))/1.e9;
			}

			if ( (SecondDifference*1000000000 + (timestamp.GetNanoSec() - last_muon_timestamp.GetNanoSec())) < 200000 ) {
				//LogDebug << "Last muon timestamp " << last_muon_timestamp.GetSec() << " " << last_muon_timestamp.GetNanoSec() << endl;  
				//LogDebug << "Last muon timestamp " << timestamp.GetSec() << " " << timestamp.GetNanoSec() << endl;  
				//LogDebug << "Event vetoed due to distance " << SecondDifference*1000000000 + (timestamp.GetNanoSec() - last_muon_timestamp.GetNanoSec()) << " from muon" << endl; 
				return true;
			}

			if (total_npe > 25000 || total_npe < 1100) return true;

			if (time.size() != PMTID.size()) {
				LogInfo << "ERROR: the time and PMTID vectors do not have the same length" << endl;
				return true;
			}
			if (time.empty()) {
				LogInfo << "ERROR: time vector is empty" << endl;
				return true;
			}

			for (int i=0; i<time.size() ; i++) {
				corr_time.push_back(time[i] - calculate_ToF(CdRecox,CdRecoy,CdRecoz,PMTs_Pos.GetX(PMTID[i]),PMTs_Pos.GetY(PMTID[i]),PMTs_Pos.GetZ(PMTID[i])));
			}

			if (!corr_time.empty()) {
				float MinTime = *std::min_element(corr_time.begin(),corr_time.end());
				for (auto& element : corr_time) {
					element -= MinTime;
				}
			}

		peak_positions.clear();

		std::vector <int> CorrTimesHistogram = MakeHistogram(corr_time);
		for (int i=0; i<200; i++) {
			CorrTimesHistogramArr[i] = CorrTimesHistogram[i];
			//LogDebug << i << "\t" << CorrTimesHistogramArr[i] << std::endl;
		}

		double* paddedWaveform = new double[240];
		double* paddedKernel = new double[240];

		for (int i=0; i<240; i++) {
			if (i<40) {
				paddedWaveform[i] = 0;
				paddedKernel[i] = kernel[i];
			} else if (i < 200) {
				paddedWaveform[i] = CorrTimesHistogramArr[i-40];
				paddedKernel[i] = kernel[i];
			} else {
			      	paddedWaveform[i] = CorrTimesHistogramArr[i-40];
				paddedKernel[i] = 0;
			}
		}

		spectrum.Clear();
		spectrum.Deconvolution(paddedWaveform,paddedKernel,240,1000,1,1);

		for (int i=0; i<200; i++) {
			CorrTimesHistogramArr[i] = paddedWaveform[i+40];	
		}

		DeconvolutedSignal.clear();
		for (int i=0; i<200; i++) {
			DeconvolutedSignal.push_back(CorrTimesHistogramArr[i]);
			//LogInfo << i << "\t" << CorrTimesHistogramArr[i] << std::endl;
		}

		delete [] paddedWaveform;
		delete [] paddedKernel;

		std::vector <int> Buffer;
  		std::vector <std::pair<int,int>> PeaksFeatures;
  		for (int i=1; i<200; i++) {
    			if(CorrTimesHistogramArr[i-1] < 5 && CorrTimesHistogramArr[i] >= 5) {
      				Buffer.clear();
      				while(true) {
        				Buffer.push_back(CorrTimesHistogramArr[i]);
        				if (CorrTimesHistogramArr[i] < 5) break;
        				i++;
      				}
				if (Buffer.size() > 20) {
				       	return true;
				}
      				auto MaxItr = std::max_element(Buffer.begin(),Buffer.end());
      				int MaxIndex = std::distance(Buffer.begin(),MaxItr);
      				PeaksFeatures.emplace_back(*MaxItr,i-Buffer.size()+MaxIndex+1);
    			}
 		}

		if (PeaksFeatures.size() == 0) {
		       	return true; 
		}

		std::vector <std::pair <int,int>> SelectedPeakPositions;
  		std::pair <int,int> Provv;
  		bool FirstCycleFlag = true;
  		double sum, mean;

  		auto AbsMax = std::max_element(PeaksFeatures.begin(),PeaksFeatures.end(),[](const std::pair <int,int>&a , const std::pair <int,int>& b){return a.first < b.first;});
  		std::pair<int,int> ProvvMax = {AbsMax->first,AbsMax->second};
  		PeaksFeatures.erase(AbsMax);

  		while (PeaksFeatures.size() > 0) {
    			auto maxIterator = std::max_element(PeaksFeatures.begin(),PeaksFeatures.end(),[](const std::pair <int,int>&a , const std::pair <int,int>& b){return a.first < b.first;});
    			Provv = {maxIterator->first,maxIterator->second};
    			PeaksFeatures.erase(maxIterator);
    			sum = std::accumulate(PeaksFeatures.begin(),PeaksFeatures.end(),0.0,[](double acc, const std::pair <int,int>&a) {return acc + a.first;});
    			if (PeaksFeatures.size() != 0) mean = sum/PeaksFeatures.size();
    			else mean = 0.;
    			if (FirstCycleFlag) {
      				if (ProvvMax.first > mean*7 && ProvvMax.first > 100) {
        				SelectedPeakPositions.emplace_back(ProvvMax.first,ProvvMax.second);
					//LogDebug << "Found a peak at " << ProvvMax.second << std::endl;  
        				FirstCycleFlag = false;
      				} else break;
    			}
    			if (Provv.first > mean*7 && Provv.first > 100) {
      				SelectedPeakPositions.emplace_back(Provv.first,Provv.second);
				//LogDebug << "Found a peak at " << Provv.second << std::endl;  
    			} else break;
  		}

		std::sort(SelectedPeakPositions.begin(),SelectedPeakPositions.end(),[](const std::pair <int,int>&a, const std::pair <int,int>&b){return a.second < b.second;});

  		for (int i=1; i<SelectedPeakPositions.size();i++) {
    			if (SelectedPeakPositions[i].second-SelectedPeakPositions[i-1].second < 33) {
      				if(SelectedPeakPositions[i].first > SelectedPeakPositions[i-1].first) {
        				SelectedPeakPositions.erase(SelectedPeakPositions.begin() + (i-1));
      				} else {
        				SelectedPeakPositions.erase(SelectedPeakPositions.begin() + i);
      				}
				i--;
    			}
  		}

		for(int i=0; i<SelectedPeakPositions.size(); i++) {
			peak_positions.push_back(SelectedPeakPositions[i].second);
		}

		NPeaks = peak_positions.size();
		
		if (NPeaks ==2) {
			//float MidPoint = (PeakPositions[1]-PeakPositions[0])/2.
			PEBi = 0;
			PEPo = 0;
			
			for (int i = 0; i<charge.size();i++) {

				if (corr_time[i] < peak_positions[0]*6 + 6*20 && corr_time[i] > peak_positions[0]*6 - 6*10) {
					PEBi += charge[i];
				}

				if (corr_time[i] < peak_positions[1]*6 + 6*20 && corr_time[i] > peak_positions[1]*6 - 6*10) {
					PEPo += charge[i];
				}
			}

			events->Fill();
			cdEvtID++;
		}
		
    	}

	return true;

}

bool BiPo212_reader::finalize() {

    if (kernel) {
        delete[] kernel;
        kernel = nullptr;
    }
    if (CorrTimesHistogramArr) {
        delete[] CorrTimesHistogramArr;
        CorrTimesHistogramArr = nullptr;
    }
	// Save muon count to output ROOT file using RootWriter
	// Fill summary tree with final muon count and run length
	nMuonsTotal = nMuons;
	// Calculate run length using min and max timestamps
	double minSec = minEventTimestamp.GetSec() + minEventTimestamp.GetNanoSec() * 1e-9;
	double maxSec = maxEventTimestamp.GetSec() + maxEventTimestamp.GetNanoSec() * 1e-9;
	if (minSec > 0 && maxSec > 0 && maxSec >= minSec) {
		runLength = maxSec - minSec;
	} else {
		runLength = 0.0;
	}
	if (summaryTree) summaryTree->Fill();
    return true;
    
}
