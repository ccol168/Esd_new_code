#ifndef RUN_READER_H
#define RUN_READER_H

#include "SniperKernel/AlgBase.h"
#include "EvtNavigator/NavBuffer.h"
#include "Identifier/IDService.h"
#include <fstream>
#include <map>
#include <TSpectrum.h>
#include "JUNO_PMTs.h"

class TH1F;
class TTree;

class BiPo212_reader : public AlgBase //Change the name to something more descriptive
{
    public :
        BiPo212_reader() : BiPo212_reader("BiPo212_reader") {}
        BiPo212_reader(const std::string& name); //Constructor, must have same name as the class

    // Following functions are needed by SNiPER, so they are mandatory
        bool initialize();
        bool execute();
        bool finalize();

    private :

        int m_iEvt, NPeaks; // To count the loops
        JM::NavBuffer* m_buf; // Our buffer with the events
	int cdEvtID = 0;        
        // Define variables that are globally used
        IDService* idServ;
	TTimeStamp timestamp, last_muon_timestamp;
        float total_npe, my_total_npe;
        std::vector<int> PMTID, peak_positions, DeconvolutedSignal; 
        std::vector<float> charge ,time, corr_time;
	//float x_CM, y_CM, z_CM;
        TString trigger_type; //, wptrigger_type;
	float CdRecox, CdRecoy, CdRecoz, CdRecoenergy, PEBi, PEPo; //, CdRecopx, CdRecopy, CdRecopz, CdRecoenergy,CdRecoPESum,CdRecot0, CdRecoPositionQuality, CdRecoEnergyQuality;
	JUNO_PMTs PMTs_Pos;
	double* CorrTimesHistogramArr;
	float Interface_level;

        TTree* events ;
	double* kernel;
	TSpectrum spectrum;
	//TTree* wpEvents;


};

#endif

