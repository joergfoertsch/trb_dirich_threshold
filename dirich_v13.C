#include "trbnet.h"
// #include "dirich_sim.C"
#include "stdint.h"
#include "unistd.h"
#include <iostream>
#include <vector>
#include "TGraph.h"
#include <chrono>
#include "TThread.h"
#include "TMutex.h"
#include <array>
#include <iomanip>
#include <random>
#include <math.h>

//**************************************
//dirich handling routines
//**************************************
// Channel Numbers:
// 0-31: TDC input channels (same index as for threshold setting)

const size_t BUFFER_SIZE4mb = 4194304;	/* 4MByte */
static uint32_t buffer4mb[4194304];


static std::mt19937_64 rnd;

#ifndef NCH
	const int NRCHANNELS = 32;			 //Nr of TDC channels in dirich
	// const int NRCHANNELS = 4;			 //Nr of TDC channels in dirich
	#define NCH
#endif

#ifndef THC
	const int OFFTHRESH = 1000; //Value to switch off channel
	const int THRESHDELAY = 100000; //Delay [mus] for thresh change to succeed
	#define THC
#endif

TMutex TRBAccessMutex;
TMutex GetRateMutex;

// bool find_min_wo_zero(uint16_t i, uint16_t j) { return (i!=0 && i<j); }

class dirich 
{

private:
// public:
	uint16_t gBoardAddress; //Board gBoardAddress
	uint64_t gBoardUID;	//UID of Board

	std::array<uint16_t,NRCHANNELS> fbaseline;
	std::array<uint16_t,NRCHANNELS> fbaseline_old;
	std::array<uint16_t,NRCHANNELS> fnoisewidth;
	std::array<uint16_t,NRCHANNELS> fnoisewidth_old;
	std::array<double,NRCHANNELS> fthresholdmV;

	std::array<uint16_t,NRCHANNELS> fsim_baseline;
	std::array<uint16_t,NRCHANNELS> fsim_noise_sigma;
	std::array<uint16_t,NRCHANNELS> fsim_noise_width;
	std::array<uint16_t,NRCHANNELS> fsim_current_threshold;
	std::array<uint16_t,NRCHANNELS> fsim_singlephotonpeakposition;
	std::array<TH1*,NRCHANNELS> fsim_distro;

	int ReadSingleThreshold(uint8_t channel, uint16_t &thrvalue);
	int ReadThresholds(std::array<uint16_t,NRCHANNELS>& thrarray);
	int WriteSingleThreshold(uint8_t channel, uint16_t thrvalue, bool check);
	int WriteThresholds(std::array<uint16_t,NRCHANNELS> thrarray, bool check);
	int WriteThresholds(uint16_t thrvalue, bool check);
	int ReadSingleScaler(uint8_t channel, uint32_t &scalervalue, std::chrono::high_resolution_clock::time_point& access_time);
	int ReadScalers(uint32_t* scalervalues, std::chrono::high_resolution_clock::time_point& access_time);
	// int GetRates(uint32_t* ratevalues, double delay=1);

	bool fsimulate = true;
	bool fjans_readout = false;
	int gdirichver = 3;

public:
	// constructor and destructor
	dirich();
	dirich(uint16_t gBoardAddress);
	virtual ~dirich();


	// converter functions
	// static double Thr_DtomV(uint16_t value) {return (double)value *2500. / 65536; }//ugly but still better than using unions
	static double Thr_DtomV(uint32_t value) {return (double)value *2500. / 65536; }
	static uint32_t Thr_mVtoD(double value) {return value /2500. *65536+.5; }


	// setter functions
	//threshold
	void SetSingleThresholdmV(uint8_t channel, double thrinmV);
	void SetThresholdsmV(std::array<double,NRCHANNELS> thrarrayinmV);
	void SetThresholdsmV(double thrinmV);
	//baseline
	void SetSingleBaseline(uint8_t channel, uint16_t baseline) {if(channel<NRCHANNELS) fbaseline[channel] = baseline; else std::cerr << "Channel: " << channel << " not specified" << std::endl;}
	void SetSingleBaseline_old(uint8_t channel, uint16_t baseline) {if(channel<NRCHANNELS) fbaseline_old[channel] = baseline; else std::cerr << "Channel: " << channel << " not specified" << std::endl;}
	// void SetBaselines(uint16_t* baseline) {for(uint8_t channel=0;channel<NRCHANNELS;++channel){fbaseline[channel] = baseline[channel];}}
	//noisewidth
	void SetSingleNoisewidth(uint8_t channel, uint16_t noisewidth) {if(channel<NRCHANNELS) fnoisewidth[channel] = noisewidth; else std::cerr << "Channel: " << channel << " not specified" << std::endl;}
	void SetSingleNoisewidth_old(uint8_t channel, uint16_t noisewidth) {if(channel<NRCHANNELS) fnoisewidth_old[channel] = noisewidth; else std::cerr << "Channel: " << channel << " not specified" << std::endl;}
	// void Setnoisewidths(uint16_t* noisewidth) {for(uint8_t channel=0;channel<NRCHANNELS;++channel){fnoisewidth[channel] = noisewidth[channel];}}


	// getter functions

	uint16_t GetBoardAddress() {return gBoardAddress;} //board address
	uint64_t GetBoardUID() {return gBoardUID;} //board address

	double	GetSingleRate(double delay=1, uint8_t channel=0); //rates from scaler
	double* GetRates(double delay=1); //rates from scaler

	// uint16_t ReadSingleThreshold(uint8_t channel){if(channel>=NRCHANNELS){std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0;} uint16_t thrarray[32]={0}; if(ReadThresholds(thrarray)!=-1) return thrarray[channel]; else{std::cerr << "Threshold could not be read" << std::endl; return 0;}};//read threshold from !dirich!
	// std::array<uint16_t,NRCHANNELS> ReadThresholds() {uint16_t* thrarray = (uint16_t*) calloc(NRCHANNELS, sizeof(uint16_t*)); if(ReadThresholds(thrarray)!=-1) return thrarray; else{uint16_t* errthrarray = (uint16_t*) calloc(NRCHANNELS, sizeof(uint16_t*)); return errthrarray;}};

	double GetSingleThresholdmV(uint8_t channel) {if(channel<NRCHANNELS) return fthresholdmV[channel]; else{std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0.;}}//threshold
	std::array<double,NRCHANNELS> GetThresholdsmV() {return fthresholdmV;}

	uint16_t GetSingleBaseline(uint8_t channel) {if(channel<NRCHANNELS) return fbaseline[channel]; else{std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0;}}//baseline
	std::array<uint16_t,NRCHANNELS> GetBaselines() {return fbaseline;}

	uint16_t GetSingleBaseline_old(uint8_t channel) {if(channel<NRCHANNELS) return fbaseline_old[channel]; else{std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0;}}//baseline
	std::array<uint16_t,NRCHANNELS> GetBaselines_old() {return fbaseline_old;}

	uint16_t GetSingleNoisewidth(uint8_t channel) {if(channel<NRCHANNELS) return fnoisewidth[channel]; else{std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0;}}//noisewidth
	std::array<uint16_t,NRCHANNELS> GetNoisewidths() {return fnoisewidth;}

	uint16_t GetSingleNoisewidth_old(uint8_t channel) {if(channel<NRCHANNELS) return fnoisewidth_old[channel]; else{std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0;}}//noisewidth
	std::array<uint16_t,NRCHANNELS> GetNoisewidths_old() {return fnoisewidth_old;}

	// threshold functions
	void DoBaselineScan	( );
	void DoBaselineScan	( uint32_t SearchedNoise, uint16_t MaxStepSize, double MeasureTime, int NrPasses);
	void AnalyzeBaseline				( );
	void AnalyzeBaseline				( uint32_t NoiseThreshold);
	void DoThreshScan	( );
	void DoThreshScan	( uint8_t FirstChannel, uint8_t LastChannel, std::array<uint16_t,NRCHANNELS> FromThr, std::array<uint16_t,NRCHANNELS> ToThr, double MeasureTime, uint16_t StepSize, int NrPasses, int clear_graph);
	void DoFineThreshScan ( );
	void DoThreshSearch	( );
	void DoThreshSearch	( double Perc, bool SPP_SPV /*1==SPP*/, double MeasureTime, int16_t StepSize, int NrPasses);
	void DoThreshScanOverBase ( );
	void DoThreshScanOverBase ( uint8_t FirstChannel, uint8_t LastChannel, std::array<double,NRCHANNELS> ToThrmV, double MeasureTime, double StepSize, int NrPasses);
	void MakeGraphsOverBase ( );
	void MakeGraphsOverBase ( uint16_t ToThr);
	void MakeDiffGraphsOverBase ( );
	// void MakeDiffGraphsOverBase ( uint16_t FromThr, uint16_t ToThr);
	void FindMinThreshScanOverBase(double gThreshold_finding_method);

	bool IsSim (){return fsimulate;}
	bool IsJansReadout (){return fjans_readout;}
	int WhichDirichVersion (){return gdirichver;}

	// threshold visualization items
	std::array<TGraph*,NRCHANNELS> gRateGraphs;
	std::array<TGraph*,NRCHANNELS> gRateGraphsOverBase;
	std::array<TGraph*,NRCHANNELS> gDiffRateGraphsOverBase;
	// void ClearGraphs();


	// settings for threshold measurement
	double 	 gMeasureTime; 	 //Time [s] to determine rate
	int 	 gLowerEdge;		 //start of thresholdscan
	int 	 gUpperEdge;		 //end of thresholdscan
	int 	 gStepsize; 	 //end of thresholdscan
	int 	 gNrPasses; 	 //Nr of passes to scan all channels

	double	 gMeasureTime_over; //Time [s] to determine rate
	int gLowerEdge_over;	//end of thresholdscan
	int gUpperEdge_over;	//end of thresholdscan
	int gStepsize_over;	 //end of thresholdscan
	int gNrPasses_over;	 //Nr of passes to scan all channels

	double	 gThreshold_finding_method; //Method to find perfect threshold: 
																			//0: searches for the minimum in the differentiated spectrum or for the minimal gradient 
																			//0<value<5: tries to find peak and sigma of the single photon distribution and sets the threshold to value*sigma
																			//5<value<100: tries to find the single photon peak and sets the threshold to value% of the spp-position
	int gdirich_reporting_level = 0;

};


dirich::dirich()
{
	dirich(0);
}

dirich::dirich(uint16_t BoardAddress)
	{

	int ret=0;
	for(int tries=0;tries<100;++tries){
		TRBAccessMutex.Lock();
		ret=trb_read_uid(BoardAddress, buffer4mb, BUFFER_SIZE4mb);
		TRBAccessMutex.UnLock();
		if(ret!=4) continue;
		if(buffer4mb[0]==0 || buffer4mb[1]==0 || buffer4mb[3] ==0) continue;
		else break;
	}
	if(ret<=0){
		std::cerr << "No DiRICH found with Address:" << std::hex << BoardAddress << "\nNot adding DiRICH" << std::endl;
		return;
	}
	else if(ret!=4){
		std::cerr << "Too many DiRICH found with Address:" << std::hex << BoardAddress << "(Amount: "<< ret/4 << ")\nNot adding DiRICH" << std::endl;
		return;
	}
	else{
		gBoardAddress = buffer4mb[3];
		uint64_t temp_store = buffer4mb[0];
		uint64_t temp_store2 = buffer4mb[1];
		gBoardUID = temp_store << 32 | temp_store2;
	}
	// std::cout << ret << " " << std::hex << buffer4mb[0] << " " << buffer4mb[1] << " " << gBoardUID << " " << gBoardAddress << std::endl;

	uint32_t cmd = 0x0 | 0xff << 24;
	uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0x10001};

	std::array<uint32_t,18> ret_c;
	int ret1=0, ret2=0;
	for(int failed=0;failed<100;++failed){
		TRBAccessMutex.Lock();
		ret1=trb_register_write_mem(gBoardAddress,0xd400,0,c,18);
		TRBAccessMutex.UnLock();
		// if(ret==-1) continue;
		usleep(100);

		TRBAccessMutex.Lock();
		ret2=trb_register_read(gBoardAddress,0xd412,ret_c.data(),2);
		TRBAccessMutex.UnLock();	
		// std::cout << std::hex << "0x" << gBoardAddress << " 0x" << gBoardUID << " 0x" << (ret_c.at(1) & 0xff00) << std::endl;
		if(ret1!=-1 && ret2==2 && (ret_c.at(1) & 0xff00) == 0x100) break;
	}
	// std::cout << ret << std::endl;
	if(ret1==-1 || ret2!=2 || (ret_c.at(1) & 0xff00) != 0x100){
		std::cerr << "No DiRICH2 Threshold FPGA of newest version detected (Version:0x" << std::hex << (ret_c.at(1) & 0xff00) << ")\nNot adding DiRICH" << std::dec << std::endl;
		return;
	}
	else{
		gdirichver = 3;
	}

	// gBoardAddress=BoardAddress;
	gdirich_reporting_level=0;
	gMeasureTime=.3;
	gLowerEdge=28000;
	gUpperEdge=32000;
	gStepsize=75;
	// gStepsize=25;
	gNrPasses=1;
	gMeasureTime_over=30.;
	gLowerEdge_over=0;
	gUpperEdge_over=600;
	gStepsize_over=10;
	gNrPasses_over=1;
	gThreshold_finding_method=0;
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		gRateGraphs[ichannel]=new TGraph();
		gRateGraphs[ichannel]->SetTitle(Form("Rate graph of dirich 0x%x's channel %i;Threshold;Rate",gBoardAddress,ichannel));
		gRateGraphs[ichannel]->SetName(Form("Rate graph of dirich 0x%x's channel %i",gBoardAddress,ichannel));

		gRateGraphsOverBase[ichannel]=new TGraph();
		gRateGraphsOverBase[ichannel]->SetTitle(Form("Rate graph over baseline of dirich 0x%x's channel %i;Threshold in mV;Rate",gBoardAddress,ichannel));
		gRateGraphsOverBase[ichannel]->SetName(Form("Rate graph over baseline of dirich 0x%x's channel %i",gBoardAddress,ichannel));

		gDiffRateGraphsOverBase[ichannel]=new TGraph();
		gDiffRateGraphsOverBase[ichannel]->SetTitle(Form("Differentiated rate graph over baseline of dirich 0x%x's channel %i;Threshold in mV;Differentiated rate",gBoardAddress,ichannel));
		gDiffRateGraphsOverBase[ichannel]->SetName(Form("Differentiated rate graph over baseline of dirich 0x%x's channel %i",gBoardAddress,ichannel));

		fbaseline[ichannel]=0;
		fbaseline_old[ichannel]=0;
		fnoisewidth[ichannel]=0;
		fnoisewidth_old[ichannel]=0;
		fthresholdmV[ichannel]=0.;
	}

	fsimulate=false;

	if(fsimulate==true){
		std::cout << "dirich " << std::hex << BoardAddress << std::dec << " not found -> simulating" << std::endl;
		rnd.seed(gBoardAddress);

		for(auto& fsim_baseline_it : fsim_baseline){
		std::uniform_int_distribution<int> dist( 29600 , 30400 );
			fsim_baseline_it = dist(rnd);
			std::cout << "fsim_baseline_it " << fsim_baseline_it << std::endl;
		}

		for(auto& fsim_noise_sigma_it : fsim_noise_sigma){
		std::uniform_int_distribution<int> dist( 50 , 800 );
			fsim_noise_sigma_it = dist(rnd);
			std::cout << "fsim_noise_sigma_it " << fsim_noise_sigma_it << std::endl;
		}

		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			std::uniform_int_distribution<int> dist( fsim_baseline.at(ichannel)+2500 , fsim_baseline.at(ichannel)+8000 );
			fsim_singlephotonpeakposition.at(ichannel) = dist(rnd);
			std::cout << "fsim_singlephotonpeakposition.at(ichannel) " << fsim_singlephotonpeakposition.at(ichannel) << std::endl;
		}

		for(auto& fsim_current_threshold_it : fsim_current_threshold){
			fsim_current_threshold_it=0;
		}

		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			int x_pos=65536;
			for(double sum=0;sum<50000;x_pos--){
				sum+=(1.0E12/(fsim_noise_sigma.at(ichannel)*sqrt(2*3.14159265359))*exp(-0.5*((x_pos-fsim_baseline.at(ichannel))/fsim_noise_sigma.at(ichannel))*((x_pos-fsim_baseline.at(ichannel))/fsim_noise_sigma.at(ichannel))));
			}
			fsim_noise_width.at(ichannel)=2*(x_pos-fsim_baseline.at(ichannel));
			std::cout << "fsim_noise_width.at(ichannel) " << fsim_noise_width.at(ichannel) << std::endl;
		}
	}

	else{
		fjans_readout=false;
		// if( BoardAddress == 0x1235 ||	BoardAddress == 0x1231 || BoardAddress == 0x1217 || BoardAddress == 0x1233 || BoardAddress == 0x1238 || BoardAddress == 0x1214 || BoardAddress == 0x1237){
		//	 fjans_readout=true;
		// }
		// else{
		//	 fjans_readout=false;
		// }
		// int ret=0;
		// //check if TDC is implemented!
		// uint32_t temp_tdc_setting[2];
		// TRBAccessMutex.Lock();
		// ret=trb_register_read(BoardAddress, 0xc802, temp_tdc_setting, 2); //get TDC-settings
		// TRBAccessMutex.UnLock();
		// std::cout << "jans readout check phase 1: ret=" << std::dec << ret << " temp_tdc_setting[0]=" << std::hex << temp_tdc_setting[0] << " temp_tdc_setting[2]=" << temp_tdc_setting[1] << std::endl; 
		// if(ret!=2 || temp_tdc_setting[0]!=BoardAddress){
		//	 std::cout << "jans readout check phase 1: is true" << std::endl;
		//	 fjans_readout=true;
		// }

		// if(fjans_readout!=true){ 
			 // TRBAccessMutex.Lock();
		//	 ret=trb_register_write(BoardAddress, 0xc802, temp_tdc_setting[1]); //reset TDC-settings
			 // TRBAccessMutex.UnLock();
		//	 std::cout << "jans readout check phase 2: ret=" << std::dec << ret << " temp_tdc_setting[0]=" << std::hex << temp_tdc_setting[0] << " temp_tdc_setting[2]=" << temp_tdc_setting[1] << std::endl; 
		//	 if(ret==-1){
		//		 std::cout << "jans readout check phase 2: is true" << std::endl;
		//		 fjans_readout=true;
		//	 }
		// }

	//	 // if(fjans_readout!=true){ 
	//	 // ret=-1;
	//	 // uint32_t scaler1[NRCHANNELS];
	//	 // // std::cout << "GetSingleRate" << scaler1 << " " << scaler2 << " " << ratevalues << std::endl;
	//	 // GetRateMutex.Lock();
	//	 // // auto start = std::chrono::high_resolution_clock::now();
	//	 // std::chrono::high_resolution_clock::time_point start1;
	//	 // for(int iterator=0;iterator<100;++iterator){
	//	 // if(ret>=0) break;
	//	 // ret=ReadScalers(scaler1,start1);
	//	 // // std::cout << iterator << "\t";
	//	 // if(iterator==99)fjans_readout=true;
	//	 // }
	//	 // GetRateMutex.UnLock();
	//	 // }
	//	 if(fjans_readout==true){ 
	//		 std::cout << "using jans scaler implementation for dirich "	<< std::hex << BoardAddress << std::dec << std::endl;
	//	 }

	}
	// std::cout << "Initialization of DiRICH " << std::hex << gBoardAddress << " complete" << std::endl;
	return;
}

dirich::~dirich()
{
	for (int i=0; i<NRCHANNELS; i++) {
		if(gRateGraphs[i]) delete gRateGraphs[i];
		if(gDiffRateGraphsOverBase[i]) delete gDiffRateGraphsOverBase[i];
	}
}

int dirich::ReadSingleThreshold(uint8_t channel, uint16_t& thrvalue)
{
	if (channel>NRCHANNELS-1)
		return -1;
	int ret;
	int reg;
	uint8_t real_channel;
	if(gdirichver==0){
		reg=0xa000+31-channel; //old firwmare
	}
	if(gdirichver==1){
		reg=0xa000+channel; //new firwmare 
	}
	if(gdirichver==2){
		real_channel = channel%16+16;
	}
	if(gdirichver==3){
		real_channel = channel%16;
	}
	if(gdirichver<=1){
		uint32_t buffer[2];
		TRBAccessMutex.Lock();
		ret=trb_register_read(gBoardAddress,reg, buffer, 2);
		TRBAccessMutex.UnLock();

		if((gBoardAddress != buffer[0]) || (ret != 2)) return -1;

		thrvalue=(buffer[1] & 0xffff); 
		return 0;
	}	
	else{
		uint32_t cmd = 0x0 << 20 | real_channel << 24 | thrvalue << 0;
		uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(uint32_t)channel/16+1,0x10001}; //evtl. sind auch mehrere Kanäle auf einmal lesbar.
		for(int failed=0;failed<100;++failed){
			TRBAccessMutex.Lock();
			ret=trb_register_write_mem(gBoardAddress,0xd400,0,c,18);
			TRBAccessMutex.UnLock();
			if(ret!=-1) break;
			usleep(1000);
		}
		uint32_t ret_c[2];
		TRBAccessMutex.Lock();
		ret=trb_register_read(gBoardAddress,0xd412,ret_c,2);
		TRBAccessMutex.UnLock();
		if((gBoardAddress != ret_c[0]) || (ret != 2)) return -1;
		thrvalue=(ret_c[1] & 0xffff);
		if(gdirich_reporting_level>3)
			std::cout << std::hex << gBoardAddress << " " << std::dec << (int)channel << " " << (int)real_channel << std::hex << " " << ret_c[0] << " " << thrvalue << std::endl;
		return 0;
	}
}

int dirich::ReadThresholds(std::array<uint16_t,NRCHANNELS>& thrarray){
	int ret=0;
	for(uint8_t ichannel=0;ichannel<NRCHANNELS;++ichannel){
		ret=ReadSingleThreshold(ichannel,thrarray.at(ichannel));
		// std::cout << std::hex << thrarray.at(ichannel) << std::endl;
		if(ret<0) return ret;
	}
	return ret;
}

int dirich::WriteSingleThreshold(uint8_t channel, uint16_t thrvalue, bool check) 
{
	if (channel>NRCHANNELS-1)
		return -1;
	if(fsimulate)
		fsim_current_threshold.at(channel) = thrvalue;
	else{
		int ret=0;

		if(gdirichver==1){
	//	 int reg=0xa000+31-channel; old firwmare
			int reg=0xa000+channel; //new firwmare 
			TRBAccessMutex.Lock();
			ret=trb_register_write(gBoardAddress, reg, (uint32_t)thrvalue);
			TRBAccessMutex.UnLock();
			if(ret==-1) return ret;
		}
		else{
			for(int failed=0;failed<100;++failed){
				uint8_t real_channel = channel%16+16*abs(gdirichver-3);
				// uint8_t real_channel = channel%16+16;
				uint32_t cmd = 0x8 << 20 | real_channel << 24 | thrvalue <<0;
				// uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(uint32_t)(channel/16+1),0x10001}; //evtl. sind auch mehrere Kanäle auf einmal setzbar.
				std::array<uint32_t,18> c = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(uint32_t)(channel/16+1),0x10001}; //evtl. sind auch mehrere Kanäle auf einmal setzbar.

				TRBAccessMutex.Lock();
				ret=trb_register_write_mem(gBoardAddress,0xd400,0,c.data(),18);
				TRBAccessMutex.UnLock();
				if(ret!=-1) break;
				usleep(1000);
			}
			if(ret==-1) return ret;
			uint32_t ret_c[18];
			TRBAccessMutex.Lock();
			ret=trb_register_read(gBoardAddress,0xd412,ret_c,18);
			TRBAccessMutex.UnLock();
			if(check){
				uint16_t set_threshold;
				ret=ReadSingleThreshold(channel, set_threshold);
				if(ret==-1) return ret;
				if(gdirich_reporting_level>2) std::cout << "wanted " << thrvalue << " set " << set_threshold << std::endl;
				if(set_threshold!=thrvalue) return -1;
			}
			// std::cout << ret << std::endl;
			return ret;
		}
	}
}

int dirich::WriteThresholds(std::array<uint16_t,NRCHANNELS> thrarray, bool check)
{
	int ret;
	if(fsimulate){
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			fsim_current_threshold.at(ichannel) = thrarray[ichannel];
		}
	}
	else{
		std::array<uint32_t,18> cmd1;
		std::array<uint32_t,18> cmd2;
		if(gdirichver==1){
			uint16_t reg=0xa000+32-NRCHANNELS;
			uint32_t buffer[NRCHANNELS];
			for (int i=0;i<NRCHANNELS;i++) {
				buffer[i]=thrarray[NRCHANNELS-i-1];
			}
			TRBAccessMutex.Lock();
			ret=trb_register_write_mem(gBoardAddress,reg,0,buffer,NRCHANNELS);
			TRBAccessMutex.UnLock();
			if(ret==-1) return ret;
		}
		else{
			uint counter=0;
			for(int channel=0;channel<16;++channel){
				if(thrarray.at(channel)){
					uint8_t real_channel = channel%16+16*abs(gdirichver-3);
					cmd1.at(counter) = 0x8 << 20 | real_channel << 24 | thrarray.at(channel) <<0;
					counter++;
				}
			}
			cmd1.at(16)=1; 
			cmd1.at(17)=0x10000 | counter; 

			counter=0;
			for(int channel=16;channel<32;++channel){
				if(thrarray.at(channel)){
					uint8_t real_channel = channel%16+16*abs(gdirichver-3);
					cmd2.at(counter) = 0x8 << 20 | real_channel << 24 | thrarray.at(channel) <<0;
					counter++;
				}
			}
			cmd2.at(16)=2;
			cmd2.at(17)=0x10000 | counter;

			for(int failed=0;failed<100;++failed){
				TRBAccessMutex.Lock();
				ret=trb_register_write_mem(gBoardAddress,0xd400,0,cmd1.data(),18);
				TRBAccessMutex.UnLock();

				if(ret!=-1) break;
				usleep(1000);
			}
			if(ret==-1) return ret;

			std::array<uint32_t,18> ret_c1;
			TRBAccessMutex.Lock();
			ret=trb_register_read(gBoardAddress,0xd412,ret_c1.data(),18);
			TRBAccessMutex.UnLock();
			if(ret==-1) return ret;

			for(int failed=0;failed<100;++failed){
				TRBAccessMutex.Lock();
				ret=trb_register_write_mem(gBoardAddress,0xd400,0,cmd2.data(),18);
				TRBAccessMutex.UnLock();
				if(ret!=-1) break;
				usleep(1000);
			}
			if(ret==-1) return ret;

			std::array<uint32_t,18> ret_c2;
			TRBAccessMutex.Lock();
			ret=trb_register_read(gBoardAddress,0xd412,ret_c2.data(),18);
			TRBAccessMutex.UnLock();
			if(ret==-1) return ret;
		}

		if(gdirich_reporting_level>2){
			std::cout << "wanted" << std::endl;
			for(auto& c_iterator : cmd1){
				std::cout<< std::hex << (c_iterator & 0xffff) << " ";
			}
			for(auto& c_iterator : cmd2){
				std::cout<< std::hex << (c_iterator & 0xffff) << " ";
			}
			std::cout << std::endl;
		}

		if(check){
			// usleep(THRESHDELAY);
			std::array<uint16_t, NRCHANNELS> set_thresholds;
			ret=ReadThresholds(set_thresholds);
			if(ret==-1) return ret;
			if(gdirich_reporting_level>2){
			std::cout << "set_blabla" << std::endl;
				for(auto& set_threshold : set_thresholds){
					std::cout << std::hex << set_threshold << " ";
				}
				std::cout << std::endl;
			}					
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				if(thrarray.at(ichannel)!=0 && abs(thrarray.at(ichannel)-set_thresholds.at(ichannel))>0){
				// if(thrarray.at(ichannel)!=0 && abs(thrarray.at(ichannel)-set_thresholds.at(ichannel))>2){
					if(gdirich_reporting_level>2){
						std::cerr << "not the same Threshold " << std::dec << ichannel << " " << std::hex << thrarray.at(ichannel) << " " << set_thresholds.at(ichannel) << std::endl;
					}
					return -1;
				} 
			}
		}
	}
	return ret;
}

int dirich::WriteThresholds(uint16_t thrvalue, bool check)
{
	int ret;
	std::array<uint16_t,NRCHANNELS> thrarray;
	thrarray.fill(thrvalue);
	ret=WriteThresholds(thrarray, check);
	return ret;
}

void dirich::SetSingleThresholdmV(uint8_t channel ,double thrinmV=30.)
{
	if(channel>=NRCHANNELS){
		std::cerr << "Channel: " << std::dec << channel << " not specified" << std::endl;
		return;
	}
	int baseline=fbaseline[channel];
	if(baseline==0){
		std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s channel: " << std::dec << unsigned(channel) << " has no baseline! (baseline==0)" << std::endl;
		return;
	} 
	int thrinD=Thr_mVtoD(thrinmV);
	int newthreshold= thrinD==0 ? 0 : baseline+thrinD;

	int ret=WriteSingleThreshold(channel, newthreshold, true);
	if(ret<0){
			std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s setting Thresholds failed" << std::endl;
	} 
	fthresholdmV.at(channel)=thrinmV;

	// printf("Setting threshold for channel %i \n",channel);
	// printf("Baseline: %i \n",baseline);
	// printf("Threshold %3.0f in digits: %i \n",thrinmV, thrinD);
	// printf("New threshold value %i set \n\n",newthreshold);
	// DrawGraphs(dirichptr,cc.at(dirichptr.first),Form("Threshold: %3.0f mV",thrinmV));
}

void dirich::SetThresholdsmV(std::array<double,NRCHANNELS> thrarrayinmV)
{
	std::array<uint16_t,NRCHANNELS> thrarrayD;
	for(int ichannel=0; ichannel<NRCHANNELS;++ichannel){
		if(fbaseline.at(ichannel)==0){
			std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s channel: " << std::dec << unsigned(ichannel) << " has no baseline! (baseline==0)" << std::endl;
	 		thrarrayD.at(ichannel)=0;
	 	}
	 	else{
		 	if(thrarrayinmV.at(ichannel)!=0)
				thrarrayD.at(ichannel)=Thr_mVtoD(thrarrayinmV.at(ichannel))+fbaseline.at(ichannel);
			else 
				thrarrayD.at(ichannel)=0;
		}
	}
	int ret=0;
	for(int tries=0;tries<100;++tries){
		ret=WriteThresholds(thrarrayD, true);
		if(ret!=-1) break;
	}
	if(ret<0){
	 std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s setting Thresholds failed" << std::endl;
	 return;
	}	
	fthresholdmV = thrarrayinmV;
}


void dirich::SetThresholdsmV(double thrinmV=30.)
{
	std::array<double,NRCHANNELS> thrarray;
	thrarray.fill(thrinmV);
	SetThresholdsmV(thrarray);
}


int dirich::ReadSingleScaler(uint8_t channel, uint32_t& scalervalue, std::chrono::high_resolution_clock::time_point& access_time)
{
	int ret;
	if (channel>NRCHANNELS)
		return -1;
	uint16_t reg=0xc000+channel;
	if(fjans_readout) reg=0xdfc0+channel; //readout of Jans implementation!!!!
	// else reg=0xc000+channel;
	uint32_t buffer[2];
	TRBAccessMutex.Lock();
	access_time = std::chrono::high_resolution_clock::now();
	ret=trb_register_read(gBoardAddress,reg, buffer, 2);
	TRBAccessMutex.UnLock();
	if ( (gBoardAddress != buffer[0]) || (ret != 2) ) 
		return -1;

	scalervalue=buffer[1] & 0x7fffffff;
	return 0;
}

int dirich::ReadScalers(uint32_t* scalervalues, std::chrono::high_resolution_clock::time_point& access_time)
{
	int ret;
	uint16_t reg=0xc000+1;
	if(fjans_readout) reg=0xdfc0; //readout of Jans implementation!!!!
	// else reg=0xc000+1;
	uint32_t buffer[NRCHANNELS+1];
	for (int i=0;i<NRCHANNELS+1; i++) buffer[i]=0;
	TRBAccessMutex.Lock();
	access_time = std::chrono::high_resolution_clock::now();
	// ret=trb_register_read(gBoardAddress,reg,buffer,NRCHANNELS);
	ret=trb_register_read_mem(gBoardAddress,reg,0,NRCHANNELS,buffer,NRCHANNELS+1);
	// ret=trb_register_read_mem(gBoardAddress,reg,0,NRCHANNELS,buffer,NRCHANNELS+1);
	TRBAccessMutex.UnLock();
	// std::cout << std::hex << buffer[0] << " " << gBoardAddress << " " << fjans_readout << std::dec << " " << ret <<	std::endl;
	// if (ret != NRCHANNELS+1) 
	// if ( gBoardAddress != buffer[0]) 
	if ( (ret == NRCHANNELS+1) && (fjans_readout || gBoardAddress == (buffer[0] & 0xffff)) ){
		for (int i=0;i<NRCHANNELS;i++){
			scalervalues[i]=buffer[i+1] & 0xfffffff;
		}
		return 0;
	}
	else{
		return -1;
	}
}

double dirich::GetSingleRate(double delay, uint8_t channel)
{
	if(fsimulate){
		long long int counter = 0;
		for(unsigned long long int itime=0;itime<delay*1000000;++itime){
		// for(unsigned long long int itime=0;itime<delay*1000000000000;++itime){
			// if(itime%1000000==0) std::cout << std::dec << itime << " ";
			GetRateMutex.Lock();
			std::uniform_real_distribution<double> is_spp_dist(0,1);
			std::normal_distribution<double> noise_dist(fsim_baseline.at(channel),fsim_noise_sigma.at(channel));
			std::normal_distribution<double> spp_dist(fsim_singlephotonpeakposition.at(channel),2000);
			if(fsim_current_threshold.at(channel)>fsim_baseline.at(channel)){
				double spp_pulse_prob = is_spp_dist(rnd);
				// std::cout << "spp_pulse_prob" << spp_pulse_prob << std::endl;
				if(is_spp_dist(rnd)<=4E-9){
					if(noise_dist(rnd)+spp_dist(rnd)>fsim_current_threshold.at(channel)){
						counter++;
						itime+=20000-1;
					}
				}
			}
			else{
				double val_1 = noise_dist(rnd);
				if((val_1>fsim_current_threshold.at(channel) && fsim_current_threshold.at(channel)>fsim_baseline.at(channel)) || (val_1<fsim_current_threshold.at(channel) && fsim_current_threshold.at(channel)<fsim_baseline.at(channel))){
					counter++;
					itime+=20000-1;
					// std::cout << "yeah1" << std::endl;
				}
			}
		}
		GetRateMutex.UnLock();
		double counter_d = 0;
		// std::cout << "fsim_current_threshold.at(ichannel) " << std::dec << fsim_current_threshold.at(ichannel) << std::endl;
		// std::cout << "counter.at(ichannel) " << counter.at(ichannel) << std::endl;
		counter_d=counter*1000000/delay;
		// counter_d.at(ichannel)=1.*counter.at(ichannel)/delay;
		// std::cout << "counter_d.at(ichannel) " << counter_d.at(ichannel) << std::endl;
		// std::cout << std::endl;
		return counter_d;
	}
	else{
		int ret=-1;
		uint32_t scaler1;
		uint32_t scaler2;
		GetRateMutex.Lock();
		std::chrono::high_resolution_clock::time_point start1;
		for(int iterator=0;iterator<100;++iterator){
			if(ret>=0) break;
			ret=ReadSingleScaler(channel,scaler1,start1);
			// std::cout << iterator << "\t";
			if(iterator==99) printf("Error reading start_scalers !\n");
		}
		GetRateMutex.UnLock();
		// std::cout << scaler1 << "\n";

		usleep(1e6*delay);
		ret=-1;
		GetRateMutex.Lock();
		std::chrono::high_resolution_clock::time_point stop1;
		for(int iterator=0;iterator<100;++iterator){
			if(ret>=0) break;
			ret=ReadSingleScaler(channel,scaler2,stop1);
			// std::cout << iterator << "\t";
			if(iterator==99) printf("Error reading end_scalers !\n");
		}
		GetRateMutex.UnLock();
		// std::cout << scaler2 << "\n";

		double exactdelay1=std::chrono::duration_cast<std::chrono::microseconds>(stop1-start1).count() * 1e-6;
		uint32_t scaler_diff=scaler2<scaler1?
			(2<<27)+scaler2-scaler1 : scaler2-scaler1;
		double rate=1.*rate/exactdelay1;
		return rate;
	}
}

double* dirich::GetRates(double delay)
{
	if(fsimulate){
		std::array <long long int,NRCHANNELS> counter;
		counter.fill(0);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(unsigned long long int itime=0;itime<delay*1000000;++itime){
			// for(unsigned long long int itime=0;itime<delay*1000000000000;++itime){
				// if(itime%1000000==0) std::cout << std::dec << itime << " ";
				// GetRateMutex.Lock();
				std::uniform_real_distribution<double> is_spp_dist(0,1);
				std::normal_distribution<double> noise_dist(fsim_baseline.at(ichannel),fsim_noise_sigma.at(ichannel));
				std::normal_distribution<double> spp_dist(fsim_singlephotonpeakposition.at(ichannel),2000);
				if(fsim_current_threshold.at(ichannel)>fsim_baseline.at(ichannel)){
					double spp_pulse_prob = is_spp_dist(rnd);
					// std::cout << "spp_pulse_prob" << spp_pulse_prob << std::endl;
					if(is_spp_dist(rnd)<=4E-9){
						if(noise_dist(rnd)+spp_dist(rnd)>fsim_current_threshold.at(ichannel)){
							counter.at(ichannel)++;
							itime+=20000-1;
						}
					}
				}
				else{
					double val_1 = noise_dist(rnd);
					if((val_1>fsim_current_threshold.at(ichannel) && fsim_current_threshold.at(ichannel)>fsim_baseline.at(ichannel)) || (val_1<fsim_current_threshold.at(ichannel) && fsim_current_threshold.at(ichannel)<fsim_baseline.at(ichannel))){
						counter.at(ichannel)++;
						itime+=20000-1;
						// std::cout << "yeah1" << std::endl;
					}
				}
			}
			// GetRateMutex.UnLock();
		}
		std::array <double,NRCHANNELS> counter_d;
		counter_d.fill(0);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			// std::cout << "fsim_current_threshold.at(ichannel) " << std::dec << fsim_current_threshold.at(ichannel) << std::endl;
			// std::cout << "counter.at(ichannel) " << counter.at(ichannel) << std::endl;
			counter_d.at(ichannel)=counter.at(ichannel)*1000000/delay;
			// counter_d.at(ichannel)=1.*counter.at(ichannel)/delay;
			// std::cout << "counter_d.at(ichannel) " << counter_d.at(ichannel) << std::endl;
		}
		// std::cout << std::endl;
		return counter_d.data();
	}
	else{
		int ret=-1;
		uint32_t scaler1[NRCHANNELS];
		uint32_t scaler2[NRCHANNELS];
		double* ratevalues = (double*) calloc(NRCHANNELS, sizeof(double*));
		GetRateMutex.Lock();
		std::chrono::high_resolution_clock::time_point start1;
		for(int iterator=0;iterator<100;++iterator){
			if(ret>=0) break;
			ret=ReadScalers(scaler1,start1);
			// std::cout << iterator << "\t";
			if(iterator==99) printf("Error reading start_scalers !\n");
		}
		GetRateMutex.UnLock();

		usleep(1e6*delay);
		ret=-1;
		GetRateMutex.Lock();
		std::chrono::high_resolution_clock::time_point stop1;
		for(int iterator=0;iterator<100;++iterator){
			if(ret>=0) break;
			ret=ReadScalers(scaler2,stop1);
			// std::cout << iterator << "\t";
			if(iterator==99) printf("Error reading end_scalers !\n");
		}
		GetRateMutex.UnLock();

		double exactdelay1=std::chrono::duration_cast<std::chrono::microseconds>(stop1-start1).count() * 1e-6;
		for (int i=0; i<NRCHANNELS; i++) {
			uint32_t scaler_diff=scaler2[i]<scaler1[i]?
				(2<<27)+scaler2[i]-scaler1[i] : scaler2[i]-scaler1[i];
			double rate=1.*rate/exactdelay1;
			ratevalues[i]=rate;
		}
		return ratevalues;
	}
}

void dirich::DoThreshScan(){
	std::array<uint16_t,NRCHANNELS> gLowerEdge_array;
	for(auto& one_gLowerEdge_array : gLowerEdge_array)
		one_gLowerEdge_array=gLowerEdge;
	std::array<uint16_t,NRCHANNELS> gUpperEdge_array;
	for(auto& one_gUpperEdge_array : gUpperEdge_array)
		one_gUpperEdge_array=gUpperEdge;
	DoThreshScan(0, NRCHANNELS, gLowerEdge_array, gUpperEdge_array, gMeasureTime, gStepsize, gNrPasses, 1);
}
void dirich::DoFineThreshScan(){
	std::array<uint16_t,NRCHANNELS> gLowerEdge_array;
	std::array<uint16_t,NRCHANNELS> gUpperEdge_array;
	for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
		gLowerEdge_array.at(ichannel)=fbaseline.at(ichannel)-0.75*fnoisewidth.at(ichannel);
		gUpperEdge_array.at(ichannel)=fbaseline.at(ichannel)+0.75*fnoisewidth.at(ichannel);
	}
	DoThreshScan(0, NRCHANNELS, gLowerEdge_array, gUpperEdge_array, gMeasureTime, 10, gNrPasses, 0);
	// for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
	//	 gUpperEdge_array.at(ichannel)=fbaseline.at(ichannel)+0.75*fnoisewidth.at(ichannel);
	//	 gLowerEdge_array.at(ichannel)=fbaseline.at(ichannel)+0.25*fnoisewidth.at(ichannel);
	// }
	// DoThreshScan(0, NRCHANNELS, gLowerEdge_array, gUpperEdge_array, gMeasureTime, 10, gNrPasses, 0); 
}
void dirich::DoThreshScan(uint8_t FirstChannel, uint8_t LastChannel, std::array<uint16_t,NRCHANNELS> FromThr, std::array<uint16_t,NRCHANNELS> ToThr, double MeasureTime, uint16_t StepSize, int NrPasses, int clear_graph)
{
	if(gdirich_reporting_level>=1){
		std::cout << std::dec << (int)FirstChannel << " " << (int)LastChannel << " " << FromThr.at(0) << " " << ToThr.at(0) << " " << MeasureTime << " " << StepSize << " " << NrPasses << " " << clear_graph << std::endl;
	}
	int ret;
	if(fjans_readout){
		uint32_t scaler_switch[] = {0xffffffff};
		TRBAccessMutex.Lock();
		ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching on scaler for this dirich .... only needed when using jan's readout
		TRBAccessMutex.UnLock();
		if(ret<0) std::cerr << "Error switching on scalers" << std::endl;
	}
	// std::cout << "Set0" << std::endl;

	if(clear_graph==1){
		for (int ichannel=0; ichannel<NRCHANNELS; ichannel++){
			gRateGraphs[ichannel]->Set(0);
			gRateGraphsOverBase[ichannel]->Set(0);
		}
	}

	for(int ipass=0;ipass<NrPasses;++ipass){
		for(int tries=0;tries<100;++tries){
			ret=WriteThresholds(OFFTHRESH, true);
			if(ret!=-1) break;
			usleep(THRESHDELAY);
		}
		if(ret<0){
		 std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s setting Thresholds failed" << std::endl;
		 return;
		}
		usleep(THRESHDELAY);
		int max_diff=0;
		for(int i=0;i<NRCHANNELS;++i){
			int temp_diff = ToThr.at(i)-FromThr.at(i);
			max_diff = temp_diff>max_diff ? temp_diff : max_diff;
		}
		for (int addthresh=0; addthresh<=max_diff; addthresh+=StepSize){
			if(gdirich_reporting_level==1){
				// if((addthresh/StepSize*100)%(max_diff/StepSize)==0)
				std::cout << "\r" << std::setw(10) << std::setprecision(2) << std::fixed << 1.*(addthresh/StepSize*100)/(max_diff/StepSize) << "%" << std::flush;
			}
			std::array<uint16_t,NRCHANNELS> threshold_value;
			for (int ichannel=0; ichannel<LastChannel; ichannel++){
				threshold_value.at(ichannel) = (ichannel+ipass)%NrPasses==0 ? FromThr.at(ichannel)+addthresh : 0;
			}
			ret=WriteThresholds(threshold_value, false);
			usleep(THRESHDELAY);
			double* rates;
			rates = GetRates(MeasureTime);
			for (int ichannel=0; ichannel<LastChannel; ichannel++){
				if(threshold_value.at(ichannel)==0) continue;
				gRateGraphs[ichannel]->SetPoint(gRateGraphs[ichannel]->GetN(),1.*threshold_value.at(ichannel),1.*rates[ichannel]);
				if(gdirich_reporting_level>2){
					std::cout << 1.*threshold_value.at(ichannel) << " " << 1.*rates[ichannel] << std::endl;
				}
			}
		}
	}
	for (int ichannel=FirstChannel; ichannel<LastChannel; ++ichannel){
		gRateGraphs[ichannel]->Sort();
	}
	if(fjans_readout){
		usleep(THRESHDELAY);
		uint32_t scaler_switch[] = {0x0};
		TRBAccessMutex.Lock();
		ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching off scaler for this dirich .... only needed when using jan's readout
		TRBAccessMutex.UnLock();
		if(ret<0) std::cerr << "Error switching off scalers" << std::endl;
	}
	if(gdirich_reporting_level==1){
		std::cout << std::endl;
	}
}

void dirich::DoThreshScanOverBase(){
	if(std::all_of(fbaseline.cbegin(), fbaseline.cend(), [](int i){ return i==0; })){
		std::cerr << "dirich 0x" << std::hex << gBoardAddress << std::dec << " has no baseline yet. Please load or scan one (load_base, system_thr_scan)" << std::endl;
		return;
	}
	std::array<double, NRCHANNELS> UpperEdge_over_array; 
	for(auto& UpperEdge_over_array_element : UpperEdge_over_array){
		UpperEdge_over_array_element=gUpperEdge_over;
	}
	DoThreshScanOverBase(0, NRCHANNELS, UpperEdge_over_array, gMeasureTime_over, gStepsize_over, gNrPasses_over);
	// DoThreshScanOverBase(0, NRCHANNELS, Thr_DtomV(gUpperEdge-*min_element(fbaseline.begin(),fbaseline.end(),find_min_wo_zero)), gMeasureTime, Thr_DtomV(gStepsize), gNrPasses);
	// std::cout << "minimal baseline = " << *min_element(fbaseline.begin(),fbaseline.end(),find_min_wo_zero) << std::endl;
}
void dirich::DoThreshScanOverBase(uint8_t FirstChannel, uint8_t LastChannel, std::array<double,NRCHANNELS> ToThrmV, double MeasureTime, double StepSize, int NrPasses){
	if(gdirich_reporting_level>=1){
		std::cout << std::dec << (int)FirstChannel << " " << (int)LastChannel << " " << ToThrmV.at(0) << " " << MeasureTime << " " << StepSize << " " << NrPasses << std::endl;
	} 
	double gMeasureTime_over_temp=3.;
	int gMeasures = MeasureTime/gMeasureTime_over_temp;
	int ret;
	if(fjans_readout){
		uint32_t scaler_switch[] = {0xffffffff};
		TRBAccessMutex.Lock();
		ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching on scaler for this dirich .... only needed when using jan's readout
		TRBAccessMutex.UnLock();
		if(ret<0) std::cerr << "Error switching on scalers" << std::endl; 
	}

	for(int ipass=0;ipass<NrPasses;++ipass){
		for(int tries=0;tries<100;++tries){
			ret=WriteThresholds(OFFTHRESH, true);
			if(ret!=-1) break;
			usleep(THRESHDELAY);
		}
		if(ret<0){
		 std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s setting Thresholds failed" << std::endl;
		 return;
		}
		usleep(THRESHDELAY);
		int max=0;
		for(int i=0;i<NRCHANNELS;++i){
			int temp_diff = ToThrmV.at(i);
			max = ToThrmV.at(i)>max ? ToThrmV.at(i) : max;
		}

		for (int addthresh=0; addthresh<=max; addthresh+=StepSize){
			if(gdirich_reporting_level==1){
				std::cout << "\r" << std::setw(10) << std::setprecision(2) << std::fixed << 1.*(addthresh/StepSize*100)/(max/StepSize) << "%" << std::flush;
			}
			std::array<uint16_t,NRCHANNELS> threshold_value;
			for (int ichannel=0; ichannel<LastChannel; ichannel++){
				threshold_value.at(ichannel) = (ichannel+ipass)%NrPasses==0 ? fbaseline.at(ichannel)+fnoisewidth.at(ichannel)/4+Thr_mVtoD(addthresh) : 0;
			}
			for(int tries=0;tries<100;++tries){
				ret=WriteThresholds(threshold_value, true);
				if(ret!=-1) break;
				usleep(THRESHDELAY);
			}
			if(ret<0){
			 std::cerr << "dirich 0x" << std::hex << gBoardAddress << "'s setting Thresholds failed" << std::endl;
			 return;
			}			
			usleep(THRESHDELAY);
			for(int measures=0;measures<gMeasures;++measures){
				double* rates;
				rates = GetRates(gMeasureTime_over_temp);
				for (int ichannel=0; ichannel<LastChannel; ichannel++){
					if(threshold_value.at(ichannel)==0) continue;
					gRateGraphsOverBase[ichannel]->SetPoint(gRateGraphsOverBase[ichannel]->GetN(),1.*Thr_DtomV(threshold_value.at(ichannel)-fbaseline.at(ichannel)),1.*rates[ichannel]);
					if(gdirich_reporting_level>2){
						std::cout << 1.*threshold_value.at(ichannel) << " " << 1.*rates[ichannel] << std::endl;
					}
				}
			}

			int finish_counter=0;
			for (int ichannel=FirstChannel+ipass; ichannel<LastChannel; ichannel+=NrPasses){
				int number_of_points = gRateGraphsOverBase[ichannel]->GetN();
				if(number_of_points>3){
					if(gRateGraphsOverBase[ichannel]->GetY()[number_of_points-1]< 3. && gRateGraphsOverBase[ichannel]->GetY()[number_of_points-2] < 3. && gRateGraphsOverBase[ichannel]->GetY()[number_of_points-3] < 3.) finish_counter++;
				}
			}
			if(finish_counter==(LastChannel-(FirstChannel+ipass))/NrPasses){
				if(gdirich_reporting_level>=1){
					std::cout << "Stopped Scan above Threshold at threshold of " << addthresh << " as no larger Rate than 3 Hz was observed in any channel for the last three thresholds" << std::endl;
				} 
				break;
			} 
		}
	}

	for (int ichannel=FirstChannel; ichannel<=LastChannel; ++ichannel){
		gRateGraphsOverBase[ichannel]->Sort();
	}
	if(fjans_readout){
		usleep(THRESHDELAY);
		uint32_t scaler_switch[] = {0x0};
		TRBAccessMutex.Lock();
		ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching off scaler for this dirich .... only needed when using jan's readout
		TRBAccessMutex.UnLock();
		if(ret<0) std::cerr << "Error switching off scalers" << std::endl;
	} 
}

void dirich::MakeDiffGraphsOverBase(){
	// MakeDiffGraphsOverBase(gLowerEdge, gUpperEdge);
// }
// void dirich::MakeDiffGraphsOverBase(uint16_t FromThr, uint16_t ToThr){
	for(int ichannel;ichannel<NRCHANNELS;++ichannel){
		gDiffRateGraphsOverBase[ichannel]->Set(0);
		if(fbaseline[ichannel]==0) continue;

		std::map<double,std::vector<double>> points_per_x;
		double* x_val = gRateGraphsOverBase[ichannel]->GetX();
		double* y_val = gRateGraphsOverBase[ichannel]->GetY();

		std::vector<std::pair<double,double> > val_med;
		for (int ipoint=0; ipoint<gRateGraphsOverBase[ichannel]->GetN(); ++ipoint) {
			if(points_per_x.count(x_val[ipoint])==0){
				std::vector<double> temp;
				points_per_x.insert(std::pair<double,std::vector<double>>(x_val[ipoint],temp));
			}
			points_per_x.at(x_val[ipoint]).push_back(y_val[ipoint]);
		}
		for(auto& x_vals : points_per_x){
			std::sort(x_vals.second.begin(),x_vals.second.end());
			double med=0;
			if(x_vals.second.size()%2==0){
				med=(x_vals.second.at(x_vals.second.size()/2)+x_vals.second.at((x_vals.second.size()+1)/2))/2;
			}
			else{
				med=x_vals.second.at(x_vals.second.size()/2);
			}
			val_med.push_back(std::pair<double,double>(x_vals.first,med));
		}
		for (int ipoint=2;ipoint<val_med.size()-2;++ipoint){
			gDiffRateGraphsOverBase[ichannel]->SetPoint(gDiffRateGraphsOverBase[ichannel]->GetN(),
																					(val_med.at(ipoint-1).first+val_med.at(ipoint+1).first) * 0.5,
																					(-val_med.at(ipoint-2).second+8*val_med.at(ipoint-1).second-8*val_med.at(ipoint+1).second+val_med.at(ipoint+2).second) /
																					(12*fabs(val_med.at(ipoint-1).first-val_med.at(ipoint+1).first)));
		}
		// //smoothing the graph
		// double* x_values = gDiffRateGraphsOverBase[ichannel]->GetX();
		// double* y_values = gDiffRateGraphsOverBase[ichannel]->GetY();
		// int number = gDiffRateGraphsOverBase[ichannel]->GetN();
		// gDiffRateGraphsOverBase[ichannel]->Set(0);
		// for (int ipoint=1; ipoint<number-1; ipoint+=3) {
		// 	gDiffRateGraphsOverBase[ichannel]->SetPoint(ipoint-1,x_values[ipoint],1.*(y_values[ipoint-1]+y_values[ipoint]+y_values[ipoint+1])/3);
		// }
	}
}

void dirich::MakeGraphsOverBase(){
	MakeGraphsOverBase(gUpperEdge);
}
void dirich::MakeGraphsOverBase(uint16_t ToThr){
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		if(fbaseline[ichannel]==0) continue;
		for(int ipoint=0; ipoint<gRateGraphs[ichannel]->GetN(); ++ipoint){
			if(gRateGraphs[ichannel]->GetX()[ipoint]<fbaseline.at(ichannel)) continue;
			gRateGraphsOverBase[ichannel]->SetPoint(gRateGraphsOverBase[ichannel]->GetN(),Thr_DtomV((gRateGraphs[ichannel]->GetX()[ipoint])-fbaseline[ichannel]),gRateGraphs[ichannel]->GetY()[ipoint]);
		}
	}
}

void dirich::AnalyzeBaseline(){
	AnalyzeBaseline(50000);
}
void dirich::AnalyzeBaseline(uint32_t NoiseThreshold)
{
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		int NrBins=gRateGraphs[ichannel]->GetN();
		int noiseedgeleft=-1;
		int noiseedgeright=-1;
		double max=0;
		double max_x=0;

		for (int ibin=0; ibin<NrBins-1; ibin++) {
			if( 
				(noiseedgeleft==-1) &&
				(gRateGraphs[ichannel]->GetY()[ibin] < NoiseThreshold) &&
				(gRateGraphs[ichannel]->GetY()[ibin+1] >= NoiseThreshold)
			){
				// linear interpolation
				double y1=gRateGraphs[ichannel]->GetY()[ibin];
				double y2=gRateGraphs[ichannel]->GetY()[ibin+1];
				double x1=gRateGraphs[ichannel]->GetX()[ibin];
				double x2=gRateGraphs[ichannel]->GetX()[ibin+1];
				noiseedgeleft=x1+(NoiseThreshold-y1)/(y2-y1)*(x2-x1);
				// std::cout << ichannel << " LeftEdge " << ibin << " " << x1 << " " << x2 << " " << y1 << " " << y2 << " " << noiseedgeleft << std::endl;
			}

			if( 
				(noiseedgeright==-1) &&
				(gRateGraphs[ichannel]->GetY()[NrBins-ibin-1] < NoiseThreshold) &&
				(gRateGraphs[ichannel]->GetY()[NrBins-ibin-2] >= NoiseThreshold)
			){
				// linear interpolation
				double y1=gRateGraphs[ichannel]->GetY()[NrBins-ibin-2];
				double y2=gRateGraphs[ichannel]->GetY()[NrBins-ibin-1];
				double x1=gRateGraphs[ichannel]->GetX()[NrBins-ibin-2];
				double x2=gRateGraphs[ichannel]->GetX()[NrBins-ibin-1];
				noiseedgeright=x2-(NoiseThreshold-y2)/(y1-y2)*(x2-x1);
				// std::cout << ichannel << " RightEdge " << ibin << " " << x1 << " " << x2 << " " << y1 << " " << y2 << " " << noiseedgeright << std::endl;
			}
			if(gRateGraphs[ichannel]->GetY()[ibin]>max){
				max=gRateGraphs[ichannel]->GetY()[ibin];
				max_x=gRateGraphs[ichannel]->GetX()[ibin];
			}
			if(noiseedgeleft!=-1 && noiseedgeright!=-1) break;
		}
		if(noiseedgeleft==-1 || noiseedgeright==-1){
			if(max==0 || NrBins<2){
				fbaseline_old[ichannel] = fbaseline[ichannel];
				fbaseline[ichannel] = 0;
				fnoisewidth_old[ichannel] = fnoisewidth[ichannel];
				fnoisewidth[ichannel] = 0;
				std::cout << "No baseline found for channel " << ichannel << "on dirich 0x" << std::hex << gBoardAddress << std::dec << std::endl;
			}
			else{
				fbaseline_old[ichannel] = fbaseline[ichannel];
				fbaseline[ichannel] = max_x;
				fnoisewidth_old[ichannel] = fnoisewidth[ichannel];
				fnoisewidth[ichannel] = 2*(gRateGraphs[ichannel]->GetX()[0]-gRateGraphs[ichannel]->GetX()[1]);
			}
		}
		else {
			fbaseline_old[ichannel] = fbaseline[ichannel];
			fbaseline[ichannel]=(noiseedgeleft+noiseedgeright) / 2.;
			fnoisewidth_old[ichannel] = fnoisewidth[ichannel];
			fnoisewidth[ichannel]=(noiseedgeright-noiseedgeleft);
			TLine* baseline_line = new TLine(fbaseline[ichannel],gRateGraphs[ichannel]->GetHistogram()->GetMinimum(),fbaseline[ichannel],gRateGraphs[ichannel]->GetHistogram()->GetMaximum());
			baseline_line->SetLineWidth(2);
			baseline_line->SetLineColor(kRed);
			gRateGraphs[ichannel]->GetListOfFunctions()->Add(baseline_line);
		} 
	}
}
