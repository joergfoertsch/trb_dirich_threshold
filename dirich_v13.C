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
// 0-31: TDC input ichannels (same index as for threshold setting)

const size_t BUFFER_SIZE4mb = 4194304;	/* 4MByte */
static uint32_t buffer4mb[4194304];


static std::mt19937_64 rnd;

#ifndef NCH
	const int NRCHANNELS = 32;			 //Nr of TDC ichannels in dirich
	const int CHPCHAIN = 16;			 //Nr of TDC ichannels pre dirich-chain
	#define NCH
#endif

#ifndef THC
	const int OFFTHRESH = 1000; //Value to switch off ichannel
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
	std::array<int,NRCHANNELS> forientation;

	int ReadSingleThreshold(uint8_t ichannel, uint16_t &thrvalue);
	int ReadThresholds(std::array<uint16_t,NRCHANNELS>& thrarray);
	int WriteSingleThreshold(uint8_t ichannel, uint16_t thrvalue, bool check);
	int WriteThresholds(std::array<uint16_t,NRCHANNELS> thrarray, bool check);
	int WriteThresholds(uint16_t thrvalue, bool check);
	int ReadSingleScaler(
		uint8_t ichannel, 
		uint32_t &scalervalue, 
		std::chrono::system_clock::time_point& access_time
	);
	int ReadScalers(uint32_t* scalervalues, std::chrono::system_clock::time_point& access_time);
	// int GetRates(uint32_t* ratevalues, double delay=1);

	int gdirichver = 3;

public:
	// constructor and destructor
	dirich();
	dirich(uint16_t gBoardAddress);
	virtual ~dirich();


	// converter functions
	static double Thr_DtomV(uint32_t value) {return (double)value *2500. / 65536; }
	static uint32_t Thr_mVtoD(double value) {return value /2500. *65536+.5; }


	// setter functions
	//threshold
	void SetSingleThresholdmV(uint8_t ichannel, double thrinmV);
	void SetThresholdsmV(std::array<double,NRCHANNELS> thrarrayinmV);
	void SetThresholdsmV(double thrinmV);
	//baseline
	void SetSingleBaseline(uint8_t ichannel, uint16_t baseline) {
		if(ichannel<NRCHANNELS) 
			fbaseline[ichannel] = baseline; 
		else 
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl;
	}
	void SetSingleBaseline_old(uint8_t ichannel, uint16_t baseline) {
		if(ichannel<NRCHANNELS) 
			fbaseline_old[ichannel] = baseline; 
		else 
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl;
	}
	//noisewidth
	void SetSingleNoisewidth(uint8_t ichannel, uint16_t noisewidth) {
		if(ichannel<NRCHANNELS) 
			fnoisewidth[ichannel] = noisewidth; 
		else 
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl;
	}
	void SetSingleNoisewidth_old(uint8_t ichannel, uint16_t noisewidth) {
		if(ichannel<NRCHANNELS) 
			fnoisewidth_old[ichannel] = noisewidth; 
		else 
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl;
	}

	// getter functions

	uint16_t GetBoardAddress() {return gBoardAddress;} //board address
	uint64_t GetBoardUID() {return gBoardUID;} //board address

	double	GetSingleRate(double delay=1, uint8_t ichannel=0); //rates from scaler
	double* GetRates(double delay=1); //rates from scaler


	double GetSingleThresholdmV(uint8_t ichannel) {
		if(ichannel<NRCHANNELS) 
			return fthresholdmV[ichannel]; 
		else{
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl; 
			return 0.;
		}
	}//threshold
	std::array<double,NRCHANNELS> GetThresholdsmV() {return fthresholdmV;}

	uint16_t GetSingleBaseline(uint8_t ichannel) {
		if(ichannel<NRCHANNELS) 
			return fbaseline[ichannel]; 
		else{
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl; 
			return 0;
		}
	}//baseline
	std::array<uint16_t,NRCHANNELS> GetBaselines() {return fbaseline;}

	uint16_t GetSingleBaseline_old(uint8_t ichannel) {
		if(ichannel<NRCHANNELS) 
			return fbaseline_old[ichannel]; 
		else{
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl; 
			return 0;
		}
	}//baseline
	std::array<uint16_t,NRCHANNELS> GetBaselines_old() {return fbaseline_old;}

	uint16_t GetSingleNoisewidth(uint8_t ichannel) {
		if(ichannel<NRCHANNELS) 
			return fnoisewidth[ichannel]; 
		else{
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl; 
			return 0;
		}
	}//noisewidth
	std::array<uint16_t,NRCHANNELS> GetNoisewidths() {return fnoisewidth;}

	uint16_t GetSingleNoisewidth_old(uint8_t ichannel) {
		if(ichannel<NRCHANNELS) return fnoisewidth_old[ichannel]; 
		else{
			std::cerr << "Channel: " << ichannel << " not specified" << std::endl; 
			return 0;
		}
	}//noisewidth
	std::array<uint16_t,NRCHANNELS> GetNoisewidths_old() {return fnoisewidth_old;}

	// threshold functions
	void DoBaselineScan	( );
	void DoBaselineScan	( 
		uint32_t SearchedNoise, 
		uint16_t MaxStepSize, 
		double MeasureTime, 
		int NrPasses
	);
	void AnalyzeBaseline ( );
	void AnalyzeBaseline ( uint32_t NoiseThreshold);
	void DoThreshScan	( );
	void DoThreshScan	( 
		uint8_t FirstChannel, 
		uint8_t LastChannel, 
		std::array<uint16_t,NRCHANNELS> FromThr, 
		std::array<uint16_t,NRCHANNELS> ToThr, 
		double MeasureTime, 
		uint16_t StepSize, 
		int NrPasses, 
		int clear_graph
	);
	void DoFineThreshScan ( );
	void DoThreshSearch	( );
	void DoThreshSearch	( 
		double Perc, 
		bool SPP_SPV /*1==SPP*/, 
		double MeasureTime, 
		int16_t StepSize, 
		int NrPasses
	);
	void DoThreshScanOverBase ( );
	void DoThreshScanOverBase ( 
		uint8_t FirstChannel, 
		uint8_t LastChannel, 
		std::array<double,NRCHANNELS> ToThrmV, 
		double MeasureTime, 
		double StepSize, 
		int NrPasses
	);
	void MakeGraphsOverBase ( );
	void MakeDiffGraphsOverBase ( );
	void MakeDiffGraphsOverBase (int case_type);
	// void MakeDiffGraphsOverBase ( uint16_t FromThr, uint16_t ToThr);
	void FindMinThreshScanOverBase(double gThreshold_finding_method);

	int WhichDirichVersion (){return gdirichver;}

	// threshold visualization items
	std::array<TGraph*,NRCHANNELS> gRateGraphs;
	std::array<TGraph*,NRCHANNELS> gRateGraphsOverBase;
	std::array<TGraph*,NRCHANNELS> gDiffRateGraphsOverBase;
	// void ClearGraphs();


	// settings for threshold measurement
	double gMeasureTime; 	 											//Time [s] to determine rate
	std::array<uint16_t,NRCHANNELS> gLowerEdge;		 	//start of thresholdscan
	std::array<uint16_t,NRCHANNELS> gUpperEdge;		 	//end of thresholdscan
	int gStepsize; 	 														//stepsize fr thresholdscan
	int gNrPasses; 	 														//Nr of passes to scan all ichannels

	double	 gMeasureTime_over; //Time [s] to determine rate
	int gLowerEdge_over;	//end of thresholdscan over baseline
	int gUpperEdge_over;	//end of thresholdscan over baseline
	int gStepsize_over;	 	//stepsize for thresholdscan
	int gNrPasses_over;	 	//Nr of passes to scan all ichannels

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
		std::cerr 
			<< "No DiRICH found with Address:" 
			<< std::hex << BoardAddress 
			<< "\nNot adding DiRICH" 
			<< std::endl;
		return;
	}
	else if(ret!=4){
		std::cerr 
			<< "Too many DiRICH found with Address:" 
			<< std::hex << BoardAddress 
			<< "(Amount: "<< ret/4 << ")\nNot adding DiRICH" 
			<< std::endl;
		return;
	}
	else{
		gBoardAddress = buffer4mb[3];
		uint64_t temp_store = buffer4mb[0];
		uint64_t temp_store2 = buffer4mb[1];
		gBoardUID = temp_store << 32 | temp_store2;
	}

	uint32_t cmd = 0x0 | 0xff << 24;
	uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0x10001};

	std::array<uint32_t,CHPCHAIN+2> ret_c;
	int ret1=0, ret2=0;
	for(int failed=0;failed<100;++failed){
		TRBAccessMutex.Lock();
		ret1=trb_register_write_mem(gBoardAddress,0xd400,0,c,CHPCHAIN+2);
		TRBAccessMutex.UnLock();
		// if(ret==-1) continue;
		usleep(100);

		TRBAccessMutex.Lock();
		ret2=trb_register_read(gBoardAddress,0xd412,ret_c.data(),2);
		TRBAccessMutex.UnLock();	
		if(ret1!=-1 && ret2==2 && (ret_c.at(1) & 0xff00) == 0x100) break; 
		//check if correct version (for dirich 0x100) is on the side FPGA
	}
	if(ret1==-1 || ret2!=2 || (ret_c.at(1) & 0xff00) != 0x100){
			std::cerr 
			<< "No DiRICH2 Threshold FPGA of newest version detected (Version:0x" 
			<< std::hex << (ret_c.at(1) & 0xff00) 
			<< ")\nNot adding DiRICH" << std::dec << std::endl;
		return;
	}
	else{
		gdirichver = 3;
	}

	// gBoardAddress=BoardAddress;
	gdirich_reporting_level=0;
	gMeasureTime=.3;
	gLowerEdge.fill(28000);
	gUpperEdge.fill(32000);
	gStepsize=75;
	// gStepsize=25;
	gNrPasses=1;
	gMeasureTime_over=30.;
	gLowerEdge_over=0;
	gUpperEdge_over=600;
	gStepsize_over=10;
	gNrPasses_over=1;
	gThreshold_finding_method=0;
	for (int iichannel=0; iichannel<NRCHANNELS; iichannel++) {
		gRateGraphs.at(iichannel)=new TGraph();
		gRateGraphs.at(iichannel)->SetTitle(
			Form(
				"Rate graph of dirich 0x%x's ichannel %i;Threshold;Rate",
				gBoardAddress,
				iichannel
			)
		);
		gRateGraphs.at(iichannel)->SetName(
			Form(
				"Rate graph of dirich 0x%x's ichannel %i",
				gBoardAddress,
				iichannel
			)
		);

		gRateGraphsOverBase.at(iichannel)=new TGraph();
		gRateGraphsOverBase.at(iichannel)->SetTitle(
			Form(

			
				"Rate graph over baseline of dirich 0x%x's ichannel %i;Threshold in mV;Rate",
				gBoardAddress,
				iichannel
			)
		);
		gRateGraphsOverBase.at(iichannel)->SetName(
			Form(
				"Rate graph over baseline of dirich 0x%x's ichannel %i",
				gBoardAddress,
				iichannel
			)
		);

		gDiffRateGraphsOverBase.at(iichannel)=new TGraph();
		gDiffRateGraphsOverBase.at(iichannel)->SetTitle(
			Form(
				"Differentiated rate graph over baseline of dirich 0x%x's ichannel %i;"
				"Threshold in mV;Differentiated rate",
				gBoardAddress,
				iichannel
			)
		);
		gDiffRateGraphsOverBase.at(iichannel)->SetName(
			Form(
				"Differentiated rate graph over baseline of dirich 0x%x's ichannel %i",
				gBoardAddress,
				iichannel
			)
		);

		fbaseline.at(iichannel) = 0;
		fbaseline_old.at(iichannel) = 0;
		fnoisewidth.at(iichannel) = 0;
		fnoisewidth_old.at(iichannel) = 0;
		fthresholdmV.at(iichannel) = 0.;
		forientation.at(iichannel) = 
		iichannel%2==0 ?
		1 : 1; //to cope with different polarities in PADIWA
	}
	return;
}

dirich::~dirich()
{
	for (int iichannel=0; iichannel<NRCHANNELS; iichannel++) {
		if(gRateGraphs.at(iichannel)) delete gRateGraphs.at(iichannel);
		if(gDiffRateGraphsOverBase.at(iichannel)) delete gDiffRateGraphsOverBase.at(iichannel);
	}
}

int dirich::ReadSingleThreshold(uint8_t ichannel, uint16_t& thrvalue)
{
	if (ichannel>NRCHANNELS-1)
		return -1;
	int ret;
	int reg;
	uint8_t real_ichannel;
	if(gdirichver==0){
		reg=0xa000+31-ichannel; //old firwmare
	}
	if(gdirichver==1){
		reg=0xa000+ichannel; //new firwmare 
	}
	if(gdirichver==2){
		real_ichannel = ichannel%CHPCHAIN+CHPCHAIN;
	}
	if(gdirichver==3){
		real_ichannel = ichannel%CHPCHAIN;
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
		uint32_t cmd = 0x0 << 20 | real_ichannel << 24 | thrvalue << 0;
		//evtl. sind auch mehrere Kanäle auf einmal lesbar.
		uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(uint32_t)ichannel/CHPCHAIN+1,0x10001}; 
		for(int failed=0;failed<100;++failed){
			TRBAccessMutex.Lock();
			ret=trb_register_write_mem(gBoardAddress,0xd400,0,c,CHPCHAIN+2);
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
			std::cout 
				<< std::hex << gBoardAddress 
				<< " " << std::dec << (int)ichannel 
				<< " " << (int)real_ichannel << std::hex 
				<< " " << ret_c[0] 
				<< " " << thrvalue 
				<< std::endl;
		return 0;
	}
}

int dirich::ReadThresholds(std::array<uint16_t,NRCHANNELS>& thrarray){
	int ret=0;
	for(uint8_t iichannel=0;iichannel<NRCHANNELS;++iichannel){
		ret=ReadSingleThreshold(iichannel,thrarray.at(iichannel));
		// std::cout << std::hex << thrarray.at(iichannel) << std::endl;
		if(ret<0) return ret;
	}
	return ret;
}

int dirich::WriteSingleThreshold(uint8_t ichannel, uint16_t thrvalue, bool check) 
{
	if (ichannel>NRCHANNELS-1)
		return -1;
	int ret=0;

	if(gdirichver==1){
//	 int reg=0xa000+31-ichannel; old firwmare
		int reg=0xa000+ichannel; //new firwmare 
		TRBAccessMutex.Lock();
		ret=trb_register_write(gBoardAddress, reg, (uint32_t)thrvalue);
		TRBAccessMutex.UnLock();
		if(ret==-1) return ret;
	}
	else{
		for(int failed=0;failed<100;++failed){
			uint8_t real_ichannel = ichannel%CHPCHAIN+CHPCHAIN*abs(gdirichver-3);
			uint32_t cmd = 0x8 << 20 | real_ichannel << 24 | thrvalue <<0;
			std::array<uint32_t,CHPCHAIN+2> c = {
				cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				(uint32_t)(ichannel/CHPCHAIN+1),
				0x00001
			}; 

			TRBAccessMutex.Lock();
			ret=trb_register_write_mem(gBoardAddress,0xd400,0,c.data(),CHPCHAIN+2);
			TRBAccessMutex.UnLock();
			if(ret!=-1) break;
			usleep(1000);
		}
		if(ret==-1) return ret;
		if(check){
			uint16_t set_threshold;
			ret=ReadSingleThreshold(ichannel, set_threshold);
			if(ret==-1) return ret;
			if(gdirich_reporting_level>2) 
				std::cout 
					<< "wanted " << thrvalue 
					<< " set " << set_threshold 
				<< std::endl;
			if(set_threshold!=thrvalue) return -1;
		}
		return ret;
	}
}

int dirich::WriteThresholds(std::array<uint16_t,NRCHANNELS> thrarray, bool check)
{
	int ret;
	std::array<std::array<uint32_t,CHPCHAIN+2>,NRCHANNELS/CHPCHAIN> cmd;
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
		std::array<uint,NRCHANNELS/CHPCHAIN> counter;
		counter.fill(0);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			if(thrarray.at(ichannel)){
				uint8_t real_channel = ichannel%CHPCHAIN+CHPCHAIN*abs(gdirichver-3);
				cmd.at(ichannel/CHPCHAIN).at(counter.at(ichannel/CHPCHAIN)) = 
					0x8 << 20 | real_channel << 24 | thrarray.at(ichannel) <<0;
				counter.at(ichannel/CHPCHAIN)++;
			}
		}
		for(int ichain=0;ichain<NRCHANNELS/CHPCHAIN;++ichain){
			cmd.at(ichain).at(CHPCHAIN)=ichain+1; 
			cmd.at(ichain).at(CHPCHAIN+1)=0x00000 | counter.at(ichain); 

			for(int failed=0;failed<100;++failed){
				TRBAccessMutex.Lock();
				ret=trb_register_write_mem(gBoardAddress,0xd400,0,cmd.at(ichain).data(),CHPCHAIN+2);
				TRBAccessMutex.UnLock();

				if(ret!=-1) break;
				usleep(1000);
			}
			if(ret==-1) return ret;

		}
	}

	if(gdirich_reporting_level>2){
		std::cout << "wanted" << std::endl;
		for(auto& chain_it : cmd){
			for(auto& c_iterator : chain_it){
				std::cout << std::hex << (c_iterator & 0xffff) << " ";
			}
		}
		std::cout << std::endl;
	}

	if(check){
		std::array<uint16_t, NRCHANNELS> set_thresholds;
		ret=ReadThresholds(set_thresholds);
		if(ret==-1) return ret;
		if(gdirich_reporting_level>2){
		std::cout << "set_thr" << std::endl;
			for(auto& set_threshold : set_thresholds){
				std::cout << std::hex << set_threshold << " ";
			}
			std::cout << std::endl;
		}
		for(int iichannel=0;iichannel<NRCHANNELS;++iichannel){
			if(thrarray.at(iichannel)!=0 && abs(thrarray.at(iichannel)-set_thresholds.at(iichannel))>0){
			// if(thrarray.at(iichannel)!=0 && abs(thrarray.at(iichannel)-set_thresholds.at(iichannel))>2){
				if(gdirich_reporting_level>2){
					std::cerr 
						<< "not the same Threshold " 
						<< std::dec << iichannel 
						<< " " << std::hex << thrarray.at(iichannel) 
						<< " " << set_thresholds.at(iichannel) 
						<< std::endl;
				}
				return -1;
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

void dirich::SetSingleThresholdmV(uint8_t ichannel ,double thrinmV=30.)
{
	if(ichannel>=NRCHANNELS){
		std::cerr << "Channel: " << std::dec << ichannel << " not specified" << std::endl;
		return;
	}
	int baseline=fbaseline[ichannel];
	if(baseline==0){
		std::cerr 
			<< "dirich 0x" 
			<< std::hex << gBoardAddress 
			<< "'s ichannel: " << std::dec << unsigned(ichannel) 
			<< " has no baseline! (baseline==0)" 
			<< std::endl;
		return;
	} 
	int newthreshold = thrinmV==0 ? 0 : baseline+forientation.at(ichannel)+Thr_mVtoD(thrinmV);

	int ret=WriteSingleThreshold(ichannel, newthreshold, false);
	if(ret<0){
			std::cerr 
			<< "dirich 0x" << std::hex << gBoardAddress 
			<< "'s setting Thresholds failed" 
			<< std::endl;
	} 
	fthresholdmV.at(ichannel) = thrinmV;

}

void dirich::SetThresholdsmV(std::array<double,NRCHANNELS> thrarrayinmV)
{
	std::array<uint16_t,NRCHANNELS> thrarrayD;
	for(int iichannel=0; iichannel<NRCHANNELS;++iichannel){
		if(fbaseline.at(iichannel)==0){
			std::cerr 
			<< "dirich 0x" << std::hex << gBoardAddress 
			<< "'s ichannel: " << std::dec << unsigned(iichannel) 
			<< " has no baseline! (baseline==0)" 
			<< std::endl;
	 		thrarrayD.at(iichannel)=0;
	 	}
	 	else{
		 	if(thrarrayinmV.at(iichannel)!=0)
				thrarrayD.at(iichannel) = 
					forientation.at(iichannel)
						*Thr_mVtoD(thrarrayinmV.at(iichannel))
					+fbaseline.at(iichannel);
			else 
				thrarrayD.at(iichannel) = 0;
		}
	}
	int ret=0;
	for(int tries=0;tries<100;++tries){
		ret=WriteThresholds(thrarrayD, false);
		if(ret!=-1) break;
	}
	if(ret<0){
	 std::cerr 
	 	<< "dirich 0x" << std::hex << gBoardAddress 
	 	<< "'s setting Thresholds failed" 
	 	<< std::endl;
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


int dirich::ReadSingleScaler(
	uint8_t ichannel, 
	uint32_t& scalervalue, 
	std::chrono::system_clock::time_point& access_time
)
{
	int ret;
	if (ichannel>NRCHANNELS)
		return -1;
	uint16_t reg=0xc000+ichannel;
	// else reg=0xc000+ichannel;
	uint32_t buffer[2];
	TRBAccessMutex.Lock();
	// auto access_time_1 = std::chrono::system_clock::now();
	ret=trb_register_read(gBoardAddress,reg, buffer, 2);
	// auto access_time_2 = std::chrono::system_clock::now();
	// std::cout 
		// << std::chrono::duration_cast<std::chrono::microseconds>(access_time_1-access_time_2).count()
	// << std::endl;
	access_time = std::chrono::system_clock::now();
	TRBAccessMutex.UnLock();
	if ( (gBoardAddress != buffer[0]) || (ret != 2) ) 
		return -1;

	scalervalue=buffer[1] & 0x7fffffff;
	return 0;
}

int dirich::ReadScalers(uint32_t* scalervalues, std::chrono::system_clock::time_point& access_time)
{
	int ret;
	uint16_t reg=0xc000+1;
	uint32_t buffer[NRCHANNELS+1];
	for (int i=0;i<NRCHANNELS+1; i++) buffer[i]=0;
	TRBAccessMutex.Lock();
	// auto access_time_1 = std::chrono::system_clock::now();
	ret=trb_register_read_mem(gBoardAddress,reg,0,NRCHANNELS,buffer,NRCHANNELS+1);
	access_time = std::chrono::system_clock::now();
	// auto access_time_2 = std::chrono::system_clock::now();
	// std::cout 
		// << std::chrono::duration_cast<std::chrono::microseconds>(access_time_1-access_time_2).count()
	// << std::endl;
	TRBAccessMutex.UnLock();
	if ( (ret == NRCHANNELS+1) && gBoardAddress == (buffer[0] & 0xffff) ){
		if(gdirich_reporting_level>=5){
			std::cout << "scalers:" << std::endl;
		}
		for (int i=0;i<NRCHANNELS;i++){
			scalervalues[i]=buffer[i+1] & 0x7fffffff;
			if(gdirich_reporting_level>=5){
				std::cout << std::dec <<scalervalues[i] << " ";
			}
		}
		if(gdirich_reporting_level>=5){
			std::cout << std::endl;
		}
		return 0;
	}
	else{
		return -1;
	}
}

double dirich::GetSingleRate(double delay, uint8_t ichannel)
{
	int ret=-1;
	uint32_t scaler1;
	uint32_t scaler2;
	GetRateMutex.Lock();
	std::chrono::system_clock::time_point start1;
	for(int iterator=0;iterator<100;++iterator){
		if(ret>=0) break;
		ret=ReadSingleScaler(ichannel,scaler1,start1);
		// std::cout << iterator << "\t";
		if(iterator==99) printf("Error reading start_scalers !\n");
	}
	GetRateMutex.UnLock();
	// std::cout << scaler1 << "\n";

	usleep(1e6*delay);
	ret=-1;
	GetRateMutex.Lock();
	std::chrono::system_clock::time_point stop1;
	for(int iterator=0;iterator<100;++iterator){
		if(ret>=0) break;
		ret=ReadSingleScaler(ichannel,scaler2,stop1);
		// std::cout << iterator << "\t";
		if(iterator==99) printf("Error reading end_scalers !\n");
	}
	GetRateMutex.UnLock();
	// std::cout << scaler2 << "\n";

	double exactdelay1 =
	 	1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(stop1-start1).count();
	uint32_t scaler_diff=scaler2<scaler1?
		(2<<31)+scaler2-scaler1 : scaler2-scaler1;
	double rate=1.*scaler_diff/exactdelay1;
	if(gdirich_reporting_level>=5){
			std::cout 
				<< exactdelay1 << std::dec 
				<< " " << scaler2 
				<< " " << scaler1 
				<< " " << scaler_diff 
				<< " " << rate 
				<< std::endl;
	}	
	return rate;
}

double* dirich::GetRates(double delay)
{
	int ret=-1;
	uint32_t scaler1[NRCHANNELS];
	uint32_t scaler2[NRCHANNELS];
	double* ratevalues = (double*) calloc(NRCHANNELS, sizeof(double));
	GetRateMutex.Lock();
	std::chrono::system_clock::time_point start1;
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
	std::chrono::system_clock::time_point stop1;
	for(int iterator=0;iterator<100;++iterator){
		if(ret>=0) break;
		ret=ReadScalers(scaler2,stop1);
		// std::cout << iterator << "\t";
		if(iterator==99) printf("Error reading end_scalers !\n");
	}
	GetRateMutex.UnLock();

	double exactdelay1 =
	 	1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(stop1-start1).count();
	for (int i=0; i<NRCHANNELS; i++) {
		uint32_t scaler_diff=scaler2[i]<scaler1[i]?
			(1<<31)+scaler2[i]-scaler1[i] : scaler2[i]-scaler1[i];
		double rate=1.*scaler_diff/exactdelay1;
		ratevalues[i]=rate;
		if(gdirich_reporting_level>=4){
			std::cout 
				<< exactdelay1 << std::dec 
				<< " " << scaler2[i] 
				<< " " << scaler1[i] 
				<< " " << scaler_diff 
				<< " " << rate 
				<< std::endl;
		}
	}
	return ratevalues;
}

void dirich::DoThreshScan(){
	DoThreshScan(
		0, 
		NRCHANNELS, 
		gLowerEdge, 
		gUpperEdge, 
		gMeasureTime, 
		gStepsize, 
		gNrPasses, 
		1
	);
}
void dirich::DoFineThreshScan(){
	std::array<uint16_t,NRCHANNELS> gLowerEdge_array;
	std::array<uint16_t,NRCHANNELS> gUpperEdge_array;
	for(int iichannel=0;iichannel<NRCHANNELS;++iichannel){
		gLowerEdge_array.at(iichannel)=fbaseline.at(iichannel)-0.75*fnoisewidth.at(iichannel);
		gUpperEdge_array.at(iichannel)=fbaseline.at(iichannel)+0.75*fnoisewidth.at(iichannel);
	}
	DoThreshScan(
		0, 
		NRCHANNELS, 
		gLowerEdge_array, 
		gUpperEdge_array, 
		gMeasureTime, 
		10, 
		gNrPasses, 
		0
	);
}
void dirich::DoThreshScan(
	uint8_t FirstChannel, 
	uint8_t LastChannel, 
	std::array<uint16_t,NRCHANNELS> FromThr, 
	std::array<uint16_t,NRCHANNELS> ToThr, 
	double MeasureTime, 
	uint16_t StepSize, 
	int NrPasses, 
	int clear_graph
)
{
	if(gdirich_reporting_level>=1){
		std::cout 
			<< std::dec << (int)FirstChannel 
			<< " " << (int)LastChannel 
			<< " " << FromThr.at(0) 
			<< " " << ToThr.at(0) 
			<< " " << MeasureTime 
			<< " " << StepSize 
			<< " " << NrPasses 
			<< " " << clear_graph 
			<< std::endl;
	}
	int ret;

	if(clear_graph==1){
		for (int iichannel=0; iichannel<NRCHANNELS; iichannel++){
			gRateGraphs.at(iichannel)->Set(0);
			gRateGraphsOverBase.at(iichannel)->Set(0);
		}
	}

	for(int ipass=0;ipass<NrPasses;++ipass){
		for(int tries=0;tries<100;++tries){
			ret=WriteThresholds(OFFTHRESH, false);
			if(ret!=-1) break;
			usleep(THRESHDELAY);
		}
		if(ret<0){
		 std::cerr 
		 	<< "dirich 0x" 
		 	<< std::hex << gBoardAddress 
		 	<< "'s setting Thresholds failed" 
		 	<< std::endl;
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
				std::cout 
				<< "\r" 
				<< std::setw(10) << std::setprecision(2) << std::fixed 
				<< 1.*(addthresh/StepSize*100)/(max_diff/StepSize) 
				<< "%" << std::flush;
			}
			std::array<uint16_t,NRCHANNELS> threshold_value;
			for (int iichannel=0; iichannel<LastChannel; iichannel++){
				threshold_value.at(iichannel) = 
				(iichannel+ipass)%NrPasses==0 ? 
				FromThr.at(iichannel)+addthresh : 0;
			}
			ret=WriteThresholds(threshold_value, false);
			usleep(THRESHDELAY);
			double* rates;
			rates = GetRates(MeasureTime);
			for (int iichannel=0; iichannel<LastChannel; iichannel++){
				if(threshold_value.at(iichannel)==0) continue;
				gRateGraphs.at(iichannel)->SetPoint(
					gRateGraphs.at(iichannel)->GetN(),
					1.*threshold_value.at(iichannel),
					1.*rates[iichannel]
				);
				if(gdirich_reporting_level>2){
					std::cout 
						<< 1.*threshold_value.at(iichannel) 
						<< " " << 1.*rates[iichannel] 
						<< std::endl;
				}
			}
		}
	}
	for (int iichannel=FirstChannel; iichannel<LastChannel; ++iichannel){
		gRateGraphs.at(iichannel)->Sort();
	}
	if(gdirich_reporting_level==1){
		std::cout << std::endl;
	}
}

void dirich::DoThreshScanOverBase(){
	if(std::all_of(fbaseline.cbegin(), fbaseline.cend(), [](int i){ return i==0; })){
		std::cerr 
			<< "dirich 0x" << std::hex << gBoardAddress 
			<< std::dec << " has no baseline yet. Please load or scan one (load_base, system_thr_scan)" 
			<< std::endl;
		return;
	}
	std::array<double, NRCHANNELS> UpperEdge_over_array; 
	for(auto& UpperEdge_over_array_element : UpperEdge_over_array){
		UpperEdge_over_array_element=gUpperEdge_over;
	}
	DoThreshScanOverBase(
		0, 
		NRCHANNELS, 
		UpperEdge_over_array, 
		gMeasureTime_over, 
		gStepsize_over, 
		gNrPasses_over
	);
}
void dirich::DoThreshScanOverBase(
	uint8_t FirstChannel, 
	uint8_t LastChannel, 
	std::array<double,NRCHANNELS> ToThrmV, 
	double MeasureTime, 
	double StepSize, 
	int NrPasses
)
{
	if(gdirich_reporting_level>=1){
		std::cout 
			<< std::dec << (int)FirstChannel 
			<< " " << (int)LastChannel 
			<< " " << ToThrmV.at(0) 
			<< " " << MeasureTime 
			<< " " << StepSize 
			<< " " << NrPasses 
			<< std::endl;
	} 
	double gMeasureTime_over_temp=3.;
	int gMeasures = MeasureTime/gMeasureTime_over_temp;
	int ret;

	for(int ipass=0;ipass<NrPasses;++ipass){
		for(int tries=0;tries<100;++tries){
			ret=WriteThresholds(OFFTHRESH, false);
			if(ret!=-1) break;
			usleep(THRESHDELAY);
		}
		if(ret<0){
		 std::cerr 
		 	<< "dirich 0x" << std::hex << gBoardAddress 
		 	<< "'s setting Thresholds failed" 
		 	<< std::endl;
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
				std::cout 
					<< "\r" 
					<< std::setw(10) << std::setprecision(2) << std::fixed 
					<< 1.*(addthresh/StepSize*100)/(max/StepSize) 
					<< "%" << std::flush;
			}
			std::array<uint16_t,NRCHANNELS> threshold_value;
			for (int iichannel=0; iichannel<LastChannel; iichannel++){
				threshold_value.at(iichannel) = 
				(iichannel+ipass)%NrPasses==0 ? 
				(
					fbaseline.at(iichannel)
					+forientation.at(iichannel)
					*(
						fnoisewidth.at(iichannel)/4
						+Thr_mVtoD(addthresh)
					)
				) : 0;
			}
			for(int tries=0;tries<100;++tries){
				ret=WriteThresholds(threshold_value, false);
				if(ret!=-1) break;
				usleep(THRESHDELAY);
			}
			if(ret<0){
			 std::cerr 
			 	<< "dirich 0x" << std::hex << gBoardAddress 
			 	<< "'s setting Thresholds failed" 
			 	<< std::endl;
			 return;
			}			
			usleep(THRESHDELAY);
			for(int measures=0;measures<gMeasures;++measures){
				double* rates = GetRates(gMeasureTime_over_temp);
				for (int iichannel=0; iichannel<LastChannel; iichannel++){
					if(threshold_value.at(iichannel)==0) continue;
					gRateGraphsOverBase.at(iichannel)->SetPoint(
						gRateGraphsOverBase.at(iichannel)->GetN(),
						1.*Thr_DtomV(threshold_value.at(iichannel)-fbaseline.at(iichannel)),
						1.*rates[iichannel]
					);
					if(gdirich_reporting_level>2){
						std::cout 
							<< 1.*threshold_value.at(iichannel) 
							<< " " << 1.*rates[iichannel] 
							<< std::endl;
					}
				}
			}

			int finish_counter=0;
			for (int iichannel=FirstChannel+ipass; iichannel<LastChannel; iichannel+=NrPasses){
				int number_of_points = gRateGraphsOverBase.at(iichannel)->GetN();
				if(number_of_points>3){
					if(
						gRateGraphsOverBase.at(iichannel)->GetY()[number_of_points-1]< 3. 
						&& gRateGraphsOverBase.at(iichannel)->GetY()[number_of_points-2] < 3. 
						&& gRateGraphsOverBase.at(iichannel)->GetY()[number_of_points-3] < 3.
						) 
							finish_counter++;
				}
			}
			if(finish_counter==(LastChannel-(FirstChannel+ipass))/NrPasses){
				if(gdirich_reporting_level>=1){
					std::cout 
					<< "Stopped Scan above Threshold at threshold of " << addthresh 
					<< " as no larger Rate than 3 Hz was observed in any ichannel for the last three thresholds" 
					<< std::endl;
				} 
				break;
			} 
		}
	}

	for (int iichannel=FirstChannel; iichannel<=LastChannel; ++iichannel){
		gRateGraphsOverBase.at(iichannel)->Sort();
	}
}

void dirich::MakeDiffGraphsOverBase(){
	MakeDiffGraphsOverBase(1);
}
void dirich::MakeDiffGraphsOverBase(int case_type){
	auto get_med = [](std::multiset<double> input_set) {
		return 1.*(
			*std::next(input_set.begin(), (int)floor(1.*input_set.size()/2))
			+(*std::next(input_set.begin(), (int)ceil(1.*input_set.size()/2)-1))
			// input_set.at((int)floor(1.*input_set.size()/2))
			// +input_set.at((int)ceil(1.*input_set.size()/2)-1)
			)/2;
	};

	for(int iichannel;iichannel<NRCHANNELS;++iichannel){
		gDiffRateGraphsOverBase.at(iichannel)->Set(0);
		if(fbaseline.at(iichannel)==0) continue;
			
		if(gRateGraphsOverBase.at(iichannel)->GetN()<2) continue;
		std::map<double,std::multiset<double>> points_per_x;
		double* x_val = gRateGraphsOverBase.at(iichannel)->GetX();
		double* y_val = gRateGraphsOverBase.at(iichannel)->GetY();

		std::vector<std::pair<double,double> > val_med;
		for (int ipoint=0; ipoint<gRateGraphsOverBase.at(iichannel)->GetN(); ++ipoint) {
			// std::cout << std::dec << "\t\t" << ipoint << std::endl;
			if(points_per_x.count(x_val[ipoint])==0){
				std::multiset<double> temp;
				points_per_x.insert(std::pair<double,std::multiset<double>>(x_val[ipoint],temp));
				// std::cout << "\t\tnew" << std::endl;
			}
			points_per_x.at(x_val[ipoint]).insert(y_val[ipoint]);
			// std::cout << std::dec << "\t\t" << points_per_x.size() << std::endl;
		}

		switch(case_type){
			case 0:
			for (int ipoint=0;ipoint<points_per_x.size();++ipoint){
				auto iterator = std::next(points_per_x.begin(),ipoint);
				// std::cout << std::dec << "\t\t" << ipoint << std::endl;
				// std::cout << std::dec << "\t\t" << ipoint << " " << *((*std::prev(iterator,1)).second.begin()) << " " << get_med((*std::prev(iterator,1)).second) << std::endl;
				gDiffRateGraphsOverBase.at(iichannel)->SetPoint(
					gDiffRateGraphsOverBase.at(iichannel)->GetN(),
					(*iterator).first,
					(
						get_med((*iterator).second)
					)
				);
			}
			break;
			case 1:
			for (int ipoint=2;ipoint<points_per_x.size()-2;++ipoint){
				auto iterator = std::next(points_per_x.begin(),ipoint);
				// std::cout << std::dec << "\t\t" << ipoint << std::endl;
				// std::cout << std::dec << "\t\t" << ipoint << " " << *((*std::prev(iterator,1)).second.begin()) << " " << get_med((*std::prev(iterator,1)).second) << std::endl;
				gDiffRateGraphsOverBase.at(iichannel)->SetPoint(
					gDiffRateGraphsOverBase.at(iichannel)->GetN(),
					((*std::prev(iterator,1)).first+(*std::next(iterator,1)).first) * 0.5,
					(
						-1.*get_med((*std::prev(iterator,2)).second)
						+8.*get_med((*std::prev(iterator,1)).second)
						-8.*get_med((*std::next(iterator,1)).second)
						+1.*get_med((*std::next(iterator,2)).second)
					)/(
						6*fabs((*std::prev(iterator,1)).first-(*std::next(iterator,1)).first)
					)
				);
			}
			break;
			case 2:
			case 3:
			default:
			for (int ipoint=0;ipoint<points_per_x.size()-1;++ipoint){
				auto iterator = std::next(points_per_x.begin(),ipoint);
				// std::cout << std::dec << "\t\t" << ipoint << std::endl;
				// std::cout << std::dec << "\t\t" << ipoint << " " << *((*std::prev(iterator,1)).second.begin()) << " " << get_med((*std::prev(iterator,1)).second) << std::endl;
				gDiffRateGraphsOverBase.at(iichannel)->SetPoint(
					gDiffRateGraphsOverBase.at(iichannel)->GetN(),
					(*iterator).first,
					(
						-1.*get_med((*std::next(iterator,1)).second)
						+get_med((*iterator).second)
					)/(
						fabs((*iterator).first-(*std::next(iterator,1)).first)
					)
				);
			}
			if(case_type==3){
				//smoothing the graph
				double* x_values = gDiffRateGraphsOverBase.at(iichannel)->GetX();
				double* y_values = gDiffRateGraphsOverBase.at(iichannel)->GetY();
				int number = gDiffRateGraphsOverBase.at(iichannel)->GetN();
				gDiffRateGraphsOverBase.at(iichannel)->Set(0);
				for (int ipoint=1; ipoint<number-1; ipoint++) {
					gDiffRateGraphsOverBase.at(iichannel)->SetPoint(
						ipoint-1,x_values[ipoint],
						1.*(y_values[ipoint-1]+y_values[ipoint]+y_values[ipoint+1])/3
					);
				}
			}
			break;
		}
	}
}

void dirich::MakeGraphsOverBase(){
	for (int iichannel=0; iichannel<NRCHANNELS; iichannel++) {
		if(fbaseline.at(iichannel)==0) continue;
		for(int ipoint=0; ipoint<gRateGraphs.at(iichannel)->GetN(); ++ipoint){
			if(gRateGraphs.at(iichannel)->GetX()[ipoint]<fbaseline.at(iichannel)) continue;
			gRateGraphsOverBase.at(iichannel)->SetPoint(
				gRateGraphsOverBase.at(iichannel)->GetN(),
				Thr_DtomV((gRateGraphs.at(iichannel)->GetX()[ipoint])-fbaseline.at(iichannel)),
				gRateGraphs.at(iichannel)->GetY()[ipoint]
			);
		}
	}
}

void dirich::AnalyzeBaseline(){
	AnalyzeBaseline(50000);
}
void dirich::AnalyzeBaseline(uint32_t NoiseThreshold)
{
	for (int iichannel=0; iichannel<NRCHANNELS; iichannel++) {
		int NrBins=gRateGraphs.at(iichannel)->GetN();
		int noiseedgeleft=-1;
		int noiseedgeright=-1;
		double max=0;
		double max_x=0;

		for (int ibin=0; ibin<NrBins-1; ibin++) {
			if( 
				(noiseedgeleft==-1) &&
				(gRateGraphs.at(iichannel)->GetY()[ibin] < NoiseThreshold) &&
				(gRateGraphs.at(iichannel)->GetY()[ibin+1] >= NoiseThreshold)
			){
				// linear interpolation
				double y1=gRateGraphs.at(iichannel)->GetY()[ibin];
				double y2=gRateGraphs.at(iichannel)->GetY()[ibin+1];
				double x1=gRateGraphs.at(iichannel)->GetX()[ibin];
				double x2=gRateGraphs.at(iichannel)->GetX()[ibin+1];
				noiseedgeleft=x1+(NoiseThreshold-y1)/(y2-y1)*(x2-x1);
			}

			if( 
				(noiseedgeright==-1) &&
				(gRateGraphs.at(iichannel)->GetY()[NrBins-ibin-1] < NoiseThreshold) &&
				(gRateGraphs.at(iichannel)->GetY()[NrBins-ibin-2] >= NoiseThreshold)
			){
				// linear interpolation
				double y1=gRateGraphs.at(iichannel)->GetY()[NrBins-ibin-2];
				double y2=gRateGraphs.at(iichannel)->GetY()[NrBins-ibin-1];
				double x1=gRateGraphs.at(iichannel)->GetX()[NrBins-ibin-2];
				double x2=gRateGraphs.at(iichannel)->GetX()[NrBins-ibin-1];
				noiseedgeright=x2-(NoiseThreshold-y2)/(y1-y2)*(x2-x1);
			}
			if(gRateGraphs.at(iichannel)->GetY()[ibin]>max){
				max=gRateGraphs.at(iichannel)->GetY()[ibin];
				max_x=gRateGraphs.at(iichannel)->GetX()[ibin];
			}
			if(gRateGraphs.at(iichannel)->GetY()[NrBins-ibin-1]>max){
				max=gRateGraphs.at(iichannel)->GetY()[NrBins-ibin-1];
				max_x=gRateGraphs.at(iichannel)->GetX()[NrBins-ibin-1];
			}
			if(noiseedgeleft!=-1 && noiseedgeright!=-1) break;
		}
		if(noiseedgeleft==-1 || noiseedgeright==-1){
			if(max==0 || NrBins<2){
				fbaseline_old.at(iichannel) = fbaseline.at(iichannel);
				fbaseline.at(iichannel) = 0;
				fnoisewidth_old.at(iichannel) = fnoisewidth.at(iichannel);
				fnoisewidth.at(iichannel) = 0;
				std::cout 
					<< "No baseline found for ichannel " 
					<< iichannel << "on dirich 0x" << std::hex << gBoardAddress 
					<< std::dec << std::endl;
			}
			else{
				fbaseline_old.at(iichannel) = fbaseline.at(iichannel);
				fbaseline.at(iichannel) = max_x;
				fnoisewidth_old.at(iichannel) = fnoisewidth.at(iichannel);
				fnoisewidth.at(iichannel) = 
					2*(
						gRateGraphs.at(iichannel)->GetX()[0]-gRateGraphs.at(iichannel)->GetX()[1]
					);
			}
		}
		else {
			fbaseline_old.at(iichannel) = fbaseline.at(iichannel);
			fbaseline.at(iichannel)=(noiseedgeleft+noiseedgeright) / 2.;
			fnoisewidth_old.at(iichannel) = fnoisewidth.at(iichannel);
			fnoisewidth.at(iichannel)=(noiseedgeright-noiseedgeleft);
			TLine* baseline_line = new TLine(fbaseline.at(iichannel),
				gRateGraphs.at(iichannel)->GetHistogram()->GetMinimum(),
				fbaseline.at(iichannel),gRateGraphs.at(iichannel)->GetHistogram()->GetMaximum()
			);
			baseline_line->SetLineWidth(2);
			baseline_line->SetLineColor(kRed);
			gRateGraphs.at(iichannel)->GetListOfFunctions()->Add(baseline_line);
		} 
	}
}
