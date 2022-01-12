#include "trbnetcom.h"
#include "stdint.h"
#include "unistd.h"
#include <iostream>
#include <vector>
#include "TGraph.h"
#include "TGraphErrors.h"
#include <chrono>
#include <array>
#include <iomanip>
#include <random>
#include <math.h>
#include <thread>
#include <mutex>

//**************************************
//dirich handling routines
//**************************************
// Channel Numbers:
// 0-31: TDC input ichannels (same index as for threshold setting)

#ifndef MB4
	// const size_t BUFFER_SIZE4mb = 4194304;	/* 4MByte */
	const size_t BUFFER_SIZE4mb = 1048576;	// 1MByte holds space for 260k DiRICHes (UID-request)
	static uint32_t buffer4mb[BUFFER_SIZE4mb];
	#define MB4
#endif

static std::mt19937_64 rnd;

#ifndef NCH
	const int NRCHANNELS = 32;			 //Nr of TDC ichannels in dirich
	const int CHPCHAIN = 16;			 //Nr of TDC ichannels pre dirich-chain
	#define NCH
#endif

#ifndef THC
	const int OFFTHRESH_high = 65535; //Value to switch off ichannel
	const int OFFTHRESH_low = 1; //Value to switch off ichannel
	const int THRESHDELAY = 100000; //Delay [mus] for thresh change to succeed
	const int SPICOMDELAY = 30000; //Delay [mus] for std SPI request to be completed
	#define THC
#endif

#ifndef REFV
	const double REF_VOLT = 2500.;
	#define REFV
#endif

const bool self_check_threshold = false;
const uint16_t MAXTHROFFSET = 5; //Maximal difference between wanted and set threshold

std::mutex GetRateMutex;

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

	inline int ReadSingleThreshold(uint8_t ichannel, uint16_t &thrvalue);
	inline int ReadThresholds(std::array<uint16_t,NRCHANNELS>& thrarray);
	inline int WriteSingleThreshold(uint8_t ichannel, uint16_t thrvalue, bool check, int nof_checks);
	inline int WriteThresholds(std::array<uint16_t,NRCHANNELS> thrarray, bool check, int nof_checks);
	inline int WriteThresholds(uint16_t thrvalue, bool check, int nof_checks);
	inline int ReadSingleScaler(
		uint8_t ichannel, 
		uint32_t &scalervalue, 
		std::chrono::system_clock::time_point& access_time
	);
	inline int ReadScalers(uint32_t* scalervalues, std::chrono::system_clock::time_point& access_time);
	// int GetRates(uint32_t* ratevalues, double delay=1);

	int gdirichver;

public:
	// constructor and destructor
	dirich();
	dirich(uint16_t gBoardAddress);
	virtual ~dirich();


	// converter functions
	static double Thr_DtomV(uint32_t value) {return (double)value *REF_VOLT / 65536; }
	static uint32_t Thr_mVtoD(double value) {return value /REF_VOLT *65536+.5; }


	// setter functions
	//threshold
	void SetSingleThresholdmV(uint8_t ichannel, double thrinmV);
	void SetThresholdsmV(std::array<double,NRCHANNELS> thrarrayinmV);
	void SetThresholdsmV(double thrinmV);
	void SetPattern(long pattern);
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
	inline int SetTDCSetting();
	inline int SetTDCSetting(uint32_t setting);
	inline int SetTDCSetting(std::array<uint32_t,NRCHANNELS/(CHPCHAIN*2)> setting);

	void SetOrientation(std::array<int,NRCHANNELS> orientation) {
		forientation = orientation;
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
	
	int GetTDCSetting();

	std::array<int,NRCHANNELS> GetOrientation() { return forientation;}

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
	std::array<TGraphErrors*,NRCHANNELS> gRateGraphs;
	std::array<TGraphErrors*,NRCHANNELS> gRateGraphsOverBase;
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

	bool gTDC_read = false;
	std::array<uint32_t,NRCHANNELS/(CHPCHAIN*2)> gTDC_setting{{0x0}};

	std::array<uint16_t,NRCHANNELS> gCurrent_Threshold;
	std::mutex Current_Thr_Mutex;
	std::chrono::steady_clock::time_point gCurrent_Threshold_time;

};


dirich::dirich()
{
	dirich(0);
}

dirich::dirich(uint16_t BoardAddress)
	{

	int ret=0;
	for(int tries=0;tries<100;++tries){
		ret=Ttrb::read_uid(BoardAddress, buffer4mb, BUFFER_SIZE4mb);
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
	ret=0;
	for(int failed=0;failed<100;++failed){
		// std::cout << gBoardAddress << " " << failed << std::endl;
		ret=Ttrb::register_write_mem(gBoardAddress,0xd400,0,c,CHPCHAIN+2);
		std::this_thread::sleep_for(std::chrono::microseconds(SPICOMDELAY));
		ret=Ttrb::register_read(gBoardAddress,0xd412,ret_c.data(),2);
		if(ret<0) return;
		if(ret==2 && (ret_c.at(1) & 0xff00) == 0x100) break;
		//check if correct version (for dirich 0x100) is on the side FPGA
	}
	if(ret!=2 || (ret_c.at(1) & 0xff00) != 0x100){
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
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		gRateGraphs.at(ichannel)=new TGraphErrors();
		gRateGraphs.at(ichannel)->SetTitle(
			Form(
				"Rate graph of dirich 0x%x's ichannel %i;Threshold;Rate",
				gBoardAddress,
				ichannel
			)
		);
		gRateGraphs.at(ichannel)->SetName(
			Form(
				"Rate graph of dirich 0x%x's ichannel %i",
				gBoardAddress,
				ichannel
			)
		);

		gRateGraphsOverBase.at(ichannel)=new TGraphErrors();
		gRateGraphsOverBase.at(ichannel)->SetTitle(
			Form(

			
				"Rate graph over baseline of dirich 0x%x's ichannel %i;Threshold in mV;Rate",
				gBoardAddress,
				ichannel
			)
		);
		gRateGraphsOverBase.at(ichannel)->SetName(
			Form(
				"Rate graph over baseline of dirich 0x%x's ichannel %i",
				gBoardAddress,
				ichannel
			)
		);

		gDiffRateGraphsOverBase.at(ichannel)=new TGraph();
		gDiffRateGraphsOverBase.at(ichannel)->SetTitle(
			Form(
				"Differentiated rate graph over baseline of dirich 0x%x's ichannel %i;"
				"Threshold in mV;Differentiated rate",
				gBoardAddress,
				ichannel
			)
		);
		gDiffRateGraphsOverBase.at(ichannel)->SetName(
			Form(
				"Differentiated rate graph over baseline of dirich 0x%x's ichannel %i",
				gBoardAddress,
				ichannel
			)
		);

		fbaseline.at(ichannel) = 0;
		fbaseline_old.at(ichannel) = 0;
		fnoisewidth.at(ichannel) = 0;
		fnoisewidth_old.at(ichannel) = 0;
		fthresholdmV.at(ichannel) = 0.;
		forientation.at(ichannel) = 
		ichannel%2==0 ?
		1 : 1; //to cope with different polarities in PADIWA
	}
	return;
}

dirich::~dirich()
{
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		if(gRateGraphs.at(ichannel)) delete gRateGraphs.at(ichannel);
		if(gDiffRateGraphsOverBase.at(ichannel)) delete gDiffRateGraphsOverBase.at(ichannel);
	}
}

int dirich::ReadSingleThreshold(uint8_t ichannel, uint16_t& thrvalue)
{
	if(self_check_threshold){
		if (ichannel>NRCHANNELS-1)
			return -1;
		int ret;
		int reg=0;
		uint8_t real_ichannel=0;
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
			ret=Ttrb::register_read(gBoardAddress,reg, buffer, 2);

			if((gBoardAddress != buffer[0]) || (ret != 2)) return -1;

			thrvalue=(buffer[1] & 0xffff); 
			return 0;
		}	
		else{
			uint32_t cmd = 0x0 << 20 | real_ichannel << 24 | thrvalue << 0;
			//evtl. sind auch mehrere Kanäle auf einmal lesbar.
			uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(uint32_t)ichannel/CHPCHAIN+1,0x10001}; 
			ret=Ttrb::register_write_mem(gBoardAddress,0xd400,0,c,CHPCHAIN+2);
			std::this_thread::sleep_for(std::chrono::microseconds(SPICOMDELAY));
			uint32_t ret_c[2];
			ret=Ttrb::register_read(gBoardAddress,0xd412,ret_c,2);
			if((gBoardAddress != ret_c[0]) || (ret != 2)) return -1;
			thrvalue=(ret_c[1] & 0xffff);
			if(gdirich_reporting_level>3)
				std::cout 
					<< std::hex << gBoardAddress 
					<< " " << std::dec << (int)ichannel 
					<< " " << (int)real_ichannel << std::hex 
					<< " " << ret_c[0] 
					<< " " << ret_c[1] 
					<< " " << thrvalue 
					<< std::endl;
			return 0;
		}
	}
	else{
		auto temp_set_time = std::chrono::steady_clock::now();
		while(true){
			Current_Thr_Mutex.lock();
			if(gCurrent_Threshold_time>temp_set_time) 
				break;
			Current_Thr_Mutex.unlock();
			std::this_thread::sleep_for(std::chrono::microseconds(2*SPICOMDELAY));
		}		
		// std::this_thread::sleep_until(
		// 	gCurrent_Threshold_time
		// 	+std::chrono::microseconds((NRCHANNELS+1)*SPICOMDELAY)
		// );
		thrvalue = gCurrent_Threshold.at(ichannel);
		Current_Thr_Mutex.unlock();
		return 0;
	}
}

int dirich::ReadThresholds(std::array<uint16_t,NRCHANNELS>& thrarray){
	int ret=0;
	if(self_check_threshold){
		for(uint8_t ichannel=0;ichannel<NRCHANNELS;++ichannel){
			ret=ReadSingleThreshold(ichannel,thrarray.at(ichannel));
			// std::cout << std::hex << thrarray.at(ichannel) << std::endl;
			if(ret<0) return ret;
		}
		return ret;
	}
	else{
		auto temp_set_time = std::chrono::steady_clock::now();
		while(true){
			Current_Thr_Mutex.lock();
			if(gCurrent_Threshold_time>temp_set_time) 
				break;
			Current_Thr_Mutex.unlock();
			std::this_thread::sleep_for(std::chrono::microseconds(2*SPICOMDELAY));
			// std::cout << "NOPE" << std::endl;
		}		
		// std::this_thread::sleep_until(
		// 	gCurrent_Threshold_time
		// 	+std::chrono::microseconds((NRCHANNELS+1)*SPICOMDELAY)
		// );		
		thrarray = gCurrent_Threshold;
		Current_Thr_Mutex.unlock();
		return 0;
	}
}

int dirich::WriteSingleThreshold(uint8_t ichannel, uint16_t thrvalue, bool check, int nof_checks)
{
	if (ichannel>NRCHANNELS-1)
		return -1;
	int ret=0;

	if(gdirichver==1){
//	 int reg=0xa000+31-ichannel; old firwmare
		int reg=0xa000+ichannel; //new firwmare 
		ret=Ttrb::register_write(gBoardAddress, reg, (uint32_t)thrvalue);
		return ret;
	}
	else{
		uint8_t real_ichannel = ichannel%CHPCHAIN+CHPCHAIN*abs(gdirichver-3);
		uint32_t cmd = 0x8 << 20 | real_ichannel << 24 | thrvalue <<0;
		std::array<uint32_t,CHPCHAIN+2> c = {{
			cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			(uint32_t)( 1 << ((int)ichannel/CHPCHAIN)),
			0x00001 | ((unsigned int)check) << 16
		}};
		ret=Ttrb::register_write_mem(gBoardAddress,0xd400,0,c.data(),CHPCHAIN+2);
		if(ret==-1) return ret;
		std::this_thread::sleep_for(std::chrono::microseconds(THRESHDELAY));		
		if(check){
			uint32_t temp[18];
			ret=Ttrb::register_read(gBoardAddress,0xd412, temp,18);
			if(ret==-1) return ret;
		}
		int failed=0;
		for(failed=0;failed<nof_checks;++failed){
			uint16_t set_threshold;
			ret=ReadSingleThreshold(ichannel, set_threshold);
			if(ret==-1) break;
			if(gdirich_reporting_level>3) 
				std::cout 
					<< "wanted " << thrvalue 
					<< " set " << set_threshold 
				<< std::endl;
			if(abs(set_threshold-thrvalue) < MAXTHROFFSET) break;
			else{
				if(gdirich_reporting_level>2) 
				std::cout 
					<< "Thresholds not matching: wanted " << thrvalue 
					<< " set " << set_threshold 
				<< std::endl;
			}
			ret=Ttrb::register_write_mem(gBoardAddress,0xd400,0,c.data(),CHPCHAIN+2);
			if(ret==-1) break;
			std::this_thread::sleep_for(std::chrono::microseconds(THRESHDELAY));
			if(check){
				uint32_t temp[18];
				ret=Ttrb::register_read(gBoardAddress,0xd412, temp,18);
				if(ret==-1) break;
			}
		}
		if(failed==nof_checks-1) return -1;
		else return ret;
	}
}

int dirich::WriteThresholds(std::array<uint16_t,NRCHANNELS> thrarray, bool check, int nof_checks)
{
	int ret;
	std::array<std::array<uint32_t,CHPCHAIN+2>,NRCHANNELS/CHPCHAIN> cmd;
	if(gdirichver==1){
		uint16_t reg=0xa000+32-NRCHANNELS;
		uint32_t buffer[NRCHANNELS];
		for (int i=0;i<NRCHANNELS;i++) {
			buffer[i]=thrarray[NRCHANNELS-i-1];
		}
		ret=Ttrb::register_write_mem(gBoardAddress,reg,0,buffer,NRCHANNELS);
		return ret;
	}
	else{
		std::array<uint,NRCHANNELS/CHPCHAIN> counter;
		counter.fill(0);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			if(thrarray.at(ichannel)){
				uint8_t real_channel = ichannel%CHPCHAIN+CHPCHAIN*abs(gdirichver-3);
				cmd.at(ichannel/CHPCHAIN).at(counter.at(ichannel/CHPCHAIN)) = 
					0x8 << 20 | real_channel << 24 | thrarray.at(ichannel) << 0;
				counter.at(ichannel/CHPCHAIN)++;
			}
		}
		int failed=0;
		for(int ichain=0;ichain<NRCHANNELS/CHPCHAIN;++ichain){
			if(counter.at(ichain) == 0x0)
				continue;
			cmd.at(ichain).at(CHPCHAIN)= 1 << ichain;
			cmd.at(ichain).at(CHPCHAIN+1)=0x00000 | counter.at(ichain) | ((int)check) << 16; 
			ret=Ttrb::register_write_mem(gBoardAddress,0xd400,0,cmd.at(ichain).data(),CHPCHAIN+2);
			if(ret==-1) break;
			std::this_thread::sleep_for(std::chrono::microseconds(THRESHDELAY));		
			if(check){
				uint32_t temp[18];
				ret=Ttrb::register_read(gBoardAddress,0xd412, temp,18);
				if(ret==-1) break;
			}
			for(;failed<nof_checks;++failed){
				std::array<uint16_t, NRCHANNELS> set_thresholds;
				ret=ReadThresholds(set_thresholds);
				if(ret==-1) break;
				int equal_it=0;
				for(int ichannel=0;ichannel<CHPCHAIN;++ichannel){
					if(gdirich_reporting_level>3){
						std::cout 
							<< thrarray.at(ichannel+CHPCHAIN*ichain) << "/"
							<< set_thresholds.at(ichannel) << "\t";
					}
					if(
						thrarray.at(ichannel+CHPCHAIN*ichain)==0 
						|| abs(
							thrarray.at(ichannel+CHPCHAIN*ichain)
							-set_thresholds.at(ichannel+CHPCHAIN*ichain)
							)<MAXTHROFFSET
						){
						equal_it++;
					}
					else if(gdirich_reporting_level>2){
						std::cout 
							<< "Threshold not equal at channel " 
							<< ichannel << " "
							<< CHPCHAIN*ichain << " "
							<< thrarray.at(ichannel+CHPCHAIN*ichain) << "/"
							<< set_thresholds.at(ichannel+CHPCHAIN*ichain) << std::endl;
					}
				}
				if(gdirich_reporting_level>2) 
					std::cout << std::endl;
				if(equal_it==CHPCHAIN) break;
				ret=Ttrb::register_write_mem(gBoardAddress,0xd400,0,cmd.at(ichain).data(),CHPCHAIN+2);
				if(ret==-1) break;
				std::this_thread::sleep_for(std::chrono::microseconds(THRESHDELAY));
				if(check){
					uint32_t temp[18];
					ret=Ttrb::register_read(gBoardAddress,0xd412, temp,18);
					if(ret==-1) break;
				}
			}
		}
		if(failed==nof_checks-1) return -1;
		else return ret;
	}
}

int dirich::WriteThresholds(uint16_t thrvalue, bool check, int nof_checks)
{
	int ret;
	std::array<uint16_t,NRCHANNELS> thrarray;
	thrarray.fill(thrvalue);
	ret=WriteThresholds(thrarray, check, nof_checks);
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

	int ret=WriteSingleThreshold(ichannel, newthreshold, true, 100);
	if(ret<0){
			std::cerr 
			<< "dirich 0x" << std::hex << gBoardAddress 
			<< "'s setting Thresholds failed" 
			<< std::endl;
	} 
	if(thrinmV!=0){
		fthresholdmV.at(ichannel) = thrinmV;
	}

}

void dirich::SetThresholdsmV(std::array<double,NRCHANNELS> thrarrayinmV)
{
	std::array<uint16_t,NRCHANNELS> thrarrayD;
	for(int ichannel=0; ichannel<NRCHANNELS;++ichannel){
		if(fbaseline.at(ichannel)==0){
			std::cerr 
			<< "dirich 0x" << std::hex << gBoardAddress 
			<< "'s ichannel: " << std::dec << unsigned(ichannel) 
			<< " has no baseline! (baseline==0)" 
			<< std::endl;
	 		thrarrayD.at(ichannel)=0;
	 	}
	 	else{
		 	if(thrarrayinmV.at(ichannel)!=0)
				thrarrayD.at(ichannel) = 
					forientation.at(ichannel)
						*Thr_mVtoD(thrarrayinmV.at(ichannel))
					+fbaseline.at(ichannel);
			else 
				thrarrayD.at(ichannel) = 0;
		}
	}
	int ret=0;
	ret=WriteThresholds(thrarrayD, true, 100);
	if(ret<0){
	 std::cerr 
	 	<< "dirich 0x" << std::hex << gBoardAddress 
	 	<< "'s setting Thresholds failed" 
	 	<< std::endl;
	 return;
	}	
	for(uint i=0;i<thrarrayinmV.size();++i){
		if(thrarrayinmV.at(i)!=0){
			fthresholdmV.at(i) = thrarrayinmV.at(i);
		}
	}
}

void dirich::SetThresholdsmV(double thrinmV=30.)
{
	std::array<double,NRCHANNELS> thrarray;
	thrarray.fill(thrinmV);
	SetThresholdsmV(thrarray);
}

void dirich::SetPattern(long pattern)
{
	std::array<uint16_t,NRCHANNELS> thresholdvals;
	for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
		thresholdvals.at(ichannel) = 
		(pattern >> ichannel) % 2 == 1 ? 
		0 : OFFTHRESH_low;
	}
	int ret=0;
	ret=WriteThresholds(thresholdvals, true, 100);
	if(ret<0){
		std::cerr 
			<< "dirich 0x" << std::hex << gBoardAddress 
			<< "'s setting Thresholds failed" 
			<< std::endl;
		return;
	}
	for(uint i=0;i<thresholdvals.size();++i){
		if(thresholdvals.at(i)!=0){
			fthresholdmV.at(i) = 0;
		}
	}
}

int dirich::SetTDCSetting(){
	if(gTDC_read==false){
		std::cerr << "no setting for TDC saved" << std::endl;
		return -1;
	}
	else{
		return SetTDCSetting(gTDC_setting);
	}
}
int dirich::SetTDCSetting(uint32_t setting){
	std::array<uint32_t,NRCHANNELS/(2*CHPCHAIN)> temp_arr;
	temp_arr.fill(setting);
	return SetTDCSetting(temp_arr);
}
int dirich::SetTDCSetting(std::array<uint32_t,NRCHANNELS/(2*CHPCHAIN)> setting){
	int ret = 0;
	for(int i=0;i<NRCHANNELS/(2*CHPCHAIN);++i){
		ret=Ttrb::register_write(gBoardAddress, 0xc802+i, setting.at(i)); //switch off TDC
		if(ret==-1){
			std::cerr 
				<< "Setting TDCs status failed for dirich " 
				<< std::hex << gBoardAddress << std::dec 
				<< std::endl;
			return -1;
		}
	}
	return ret;
}

int dirich::ReadSingleScaler(
	uint8_t ichannel, 
	uint32_t& scalervalue, 
	std::chrono::system_clock::time_point& access_time
)
{
	int ret=-1;
	if (ichannel>NRCHANNELS)
		return -1;
	uint16_t reg=0xc000+ichannel;
	// else reg=0xc000+ichannel;
	uint32_t buffer[2];
	for(int tries=0;tries<Ttrb::NOFCOMTRIES;++tries){
		GetRateMutex.lock();
		ret=Ttrb::register_read(gBoardAddress,reg, buffer, 2);
		access_time = std::chrono::system_clock::now();
		GetRateMutex.unlock();
		if( ret==2 && gBoardAddress == buffer[0] ) 
			break;
		std::this_thread::sleep_for(std::chrono::milliseconds(Ttrb::faildelay(Ttrb::gen)));
	}
	if( ret==2 && gBoardAddress == buffer[0] ){
		scalervalue=buffer[1] & 0x7fffffff;
		return 0;
	}
	else
		return -1;
}

int dirich::ReadScalers(uint32_t* scalervalues, std::chrono::system_clock::time_point& access_time)
{
	// std::cout << "Getting Scalers" << std::endl;
	int ret;
	uint16_t reg=0xc000+1;
	uint32_t buffer[NRCHANNELS+1];
	// for (int i=0;i<NRCHANNELS+1; i++) buffer[i]=0;
	for(int tries=0;tries<Ttrb::NOFCOMTRIES;++tries){
		// std::cout << "Getting Scalers try " << tries << std::endl;
		GetRateMutex.lock();
		ret=Ttrb::register_read_mem(gBoardAddress,reg,0,NRCHANNELS,buffer,NRCHANNELS+1);
		// for(int i=0;i<ret;++i) std::cout << buffer[i] <<std::endl;
		access_time = std::chrono::system_clock::now();
		GetRateMutex.unlock();
		if( ret==NRCHANNELS+1 && gBoardAddress == (buffer[0] & 0xffff) ) 
			break;
		std::this_thread::sleep_for(std::chrono::milliseconds(Ttrb::faildelay(Ttrb::gen)));
	}	
	if( (ret == NRCHANNELS+1) && gBoardAddress == (buffer[0] & 0xffff) ){
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
	else
		return -1;
}

double dirich::GetSingleRate(double delay, uint8_t ichannel)
{
	int ret=-1;
	uint32_t scaler1=0;
	uint32_t scaler2=0;

	std::chrono::system_clock::time_point start1;
	ret=ReadSingleScaler(ichannel,scaler1,start1);
	if(ret<0){
		std::cerr << "Error reading start_scalers" << std::endl;
		return -1;
	}

	std::this_thread::sleep_for(std::chrono::microseconds((int)(1e6*delay)));

	ret=-1;
	std::chrono::system_clock::time_point stop1;
	ret=ReadSingleScaler(ichannel,scaler2,stop1);
	if(ret<0){
		std::cerr << "Error reading end_scalers" << std::endl;
		return -2;
	}	

	double exactdelay1 =
	 	1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(stop1-start1).count();
	uint32_t scaler_diff=scaler2<scaler1?
		(1<<31)+scaler2-scaler1 : scaler2-scaler1;
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
	// std::cout << "Getting rates" << std::endl;

	int ret=-1;
	uint32_t scaler1[NRCHANNELS];
	uint32_t scaler2[NRCHANNELS];
	double* ratevalues = (double*) calloc(NRCHANNELS, sizeof(double));

	std::chrono::system_clock::time_point start1;
	ret=ReadScalers(scaler1,start1);
	if(ret<0){
		std::cerr << "Error reading start_scalers" << std::endl;
		for(int i=0;i<NRCHANNELS;++i)
			ratevalues[i] = -1;
		return ratevalues;
	}
	// std::cout << "Going to sleep " << std::hex << gBoardAddress << std::endl;
	std::this_thread::sleep_for(std::chrono::microseconds((int)(1e6*delay)));

	ret=-1;
	std::chrono::system_clock::time_point stop1;
	ret=ReadScalers(scaler2,stop1);
	if(ret<0){
		std::cerr << "Error reading start_scalers" << std::endl;
		for(int i=0;i<NRCHANNELS;++i)
			ratevalues[i] = -1;
		return ratevalues;
	}

	double exactdelay1 =
	 	1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(stop1-start1).count();
	for (int i=0; i<NRCHANNELS; i++) {
		uint64_t scaler_diff=scaler2[i]<scaler1[i]?
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
	// std::cout << "DONE Getting rates" << std::endl;
	return ratevalues;
}

int dirich::GetTDCSetting(){
	int ret = 0;
	for(int i=0;i<NRCHANNELS/(2*CHPCHAIN);++i){
		uint32_t temp_tdc_setting[2];
		ret=Ttrb::register_read(gBoardAddress, 0xc802+i, temp_tdc_setting, 2); //switch off TDC
		// std::cout << std::hex << temp_tdc_setting[0] << "\t" << temp_tdc_setting[1] << std::endl;
		if(ret!=2 || temp_tdc_setting[0]!=gBoardAddress){
			std::cerr 
				<< "Reading TDCs status failed for dirich " 
				<< std::hex << gBoardAddress << std::dec 
				<< " -> TDC for that dirich will be left switched off" 
				<< std::endl;
			gTDC_read = false;
			return -1;
		}
		gTDC_setting.at(i) = temp_tdc_setting[1];
		gTDC_read = true;
	}
	return ret;	
	
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
	for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
		gLowerEdge_array.at(ichannel)=fbaseline.at(ichannel)-0.75*fnoisewidth.at(ichannel);
		gUpperEdge_array.at(ichannel)=fbaseline.at(ichannel)+0.75*fnoisewidth.at(ichannel);
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
		for (int ichannel=0; ichannel<NRCHANNELS; ichannel++){
			gRateGraphs.at(ichannel)->Set(0);
			gRateGraphsOverBase.at(ichannel)->Set(0);
		}
	}
	
	GetTDCSetting(); //Get TDC values
	SetTDCSetting(0x0); //SwitchTDC off

	for(int ipass=0;ipass<NrPasses;++ipass){
		ret=WriteThresholds(OFFTHRESH_low, true, 0);
		if(ret<0){
		 std::cerr 
		 	<< "dirich 0x" 
		 	<< std::hex << gBoardAddress 
		 	<< "'s setting Thresholds failed" 
		 	<< std::endl;
		 return;
		}
		int max_diff=0;
		for(int i=0+ipass;i<NRCHANNELS;i+=NrPasses){
			int temp_diff = ToThr.at(i)-FromThr.at(i);
			max_diff = temp_diff>max_diff ? temp_diff : max_diff;
		}
		for (int addthresh=0; addthresh<=max_diff; addthresh+=StepSize){
			if(gdirich_reporting_level==1){
				// if((addthresh/StepSize*100)%(max_diff/StepSize)==0)
				std::cout 
				<< "\r" 
				<< std::setw(10) << std::setprecision(2) << std::fixed 
				<< 1.*((addthresh/StepSize+ipass*max_diff/StepSize)*100)/(max_diff/StepSize*NrPasses)
				<< "%" << std::flush;
			}
			std::array<uint16_t,NRCHANNELS> threshold_value;
			for (int ichannel=0; ichannel<LastChannel; ichannel++){
				threshold_value.at(ichannel) = 
				(ichannel+ipass)%NrPasses==0 ? 
				FromThr.at(ichannel)+addthresh : 0;
			}
			ret=WriteThresholds(threshold_value, true, 0);
			if(ret<0){
			 std::cerr 
			 	<< "dirich 0x" 
			 	<< std::hex << gBoardAddress 
			 	<< "'s setting Thresholds failed" 
			 	<< std::endl;
			 return;
			}
			double* rates;
			rates = GetRates(MeasureTime);
			for (int ichannel=0; ichannel<LastChannel; ichannel++){
				if(threshold_value.at(ichannel)==0) continue;
				gRateGraphs.at(ichannel)->SetPoint(
					gRateGraphs.at(ichannel)->GetN(),
					1.*threshold_value.at(ichannel),
					1.*rates[ichannel]
				);
				gRateGraphs.at(ichannel)->SetPointError(
					gRateGraphs.at(ichannel)->GetN()-1,
					1,
					sqrt(1.*rates[ichannel])/sqrt(MeasureTime)
				);				
				if(gdirich_reporting_level>3){  
					std::cout 
						<< 1.*threshold_value.at(ichannel) 
						<< " " << 1.*rates[ichannel] 
						<< std::endl;
				}
			}
		}
	}

	SetTDCSetting(); //Set TDC back to original values

	for (int ichannel=FirstChannel; ichannel<LastChannel; ++ichannel){
		gRateGraphs.at(ichannel)->Sort();
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
	
	GetTDCSetting(); //Get TDC values
	SetTDCSetting(0x0); //SwitchTDC off

	for(int ipass=0;ipass<NrPasses;++ipass){
		ret=WriteThresholds(OFFTHRESH_low, true, 100);
		if(ret<0){
		 std::cerr 
		 	<< "dirich 0x" << std::hex << gBoardAddress 
		 	<< "'s setting Thresholds failed" 
		 	<< std::endl;
		 return;
		}
		int max=0;
		for(int i=0+ipass;i<NRCHANNELS;i+=NrPasses){
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
			for (int ichannel=0; ichannel<LastChannel; ichannel++){
				threshold_value.at(ichannel) = 
				(ichannel+ipass)%NrPasses==0 ? 
				(
					fbaseline.at(ichannel)
					+forientation.at(ichannel)
					*(
						fnoisewidth.at(ichannel)/4
						+Thr_mVtoD(addthresh)
					)
				) : 0;
			}
			ret=WriteThresholds(threshold_value, true, 100);
			if(ret<0){
			 std::cerr 
			 	<< "dirich 0x" << std::hex << gBoardAddress 
			 	<< "'s setting Thresholds failed" 
			 	<< std::endl;
			 return;
			}
			if(gMeasures)for(int measures=0;measures<gMeasures;++measures){
				double* rates = GetRates(gMeasureTime_over_temp);
				for (int ichannel=0; ichannel<LastChannel; ichannel++){
					if(threshold_value.at(ichannel)==0) continue;
					gRateGraphsOverBase.at(ichannel)->SetPoint(
						gRateGraphsOverBase.at(ichannel)->GetN(),
						1.*Thr_DtomV(threshold_value.at(ichannel)-fbaseline.at(ichannel)),
						1.*rates[ichannel]
					);
					if(gdirich_reporting_level>2){
						std::cout 
							<< 1.*threshold_value.at(ichannel) 
							<< " " << 1.*rates[ichannel] 
							<< std::endl;
					}
				}
			}
			else{
				double* rates = GetRates(MeasureTime);
				for (int ichannel=0; ichannel<LastChannel; ichannel++){
					if(threshold_value.at(ichannel)==0) continue;
					gRateGraphsOverBase.at(ichannel)->SetPoint(
						gRateGraphsOverBase.at(ichannel)->GetN(),
						1.*Thr_DtomV(threshold_value.at(ichannel)-fbaseline.at(ichannel)),
						1.*rates[ichannel]
					);
					if(gdirich_reporting_level>2){
						std::cout 
							<< 1.*threshold_value.at(ichannel) 
							<< " " << 1.*rates[ichannel] 
							<< std::endl;
					}
				}
			}
			int finish_counter=0;
			for (int ichannel=FirstChannel+ipass; ichannel<LastChannel; ichannel+=NrPasses){
				int number_of_points = gRateGraphsOverBase.at(ichannel)->GetN();
				if(number_of_points>3.*gMeasures){
					if(
						gRateGraphsOverBase.at(ichannel)->GetY()[number_of_points-1]< 3. 
						&& gRateGraphsOverBase.at(ichannel)->GetY()[number_of_points-2] < 3. 
						&& gRateGraphsOverBase.at(ichannel)->GetY()[number_of_points-3] < 3.
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

	SetTDCSetting(); //Set TDC back to original values

	for (int ichannel=FirstChannel; ichannel<LastChannel; ++ichannel){
		gRateGraphsOverBase.at(ichannel)->Sort();
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

	for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
		gDiffRateGraphsOverBase.at(ichannel)->Set(0);
		if(fbaseline.at(ichannel)==0) continue;
			
		if(gRateGraphsOverBase.at(ichannel)->GetN()<2) continue;
		std::map<double,std::multiset<double>> points_per_x;
		double* x_val = gRateGraphsOverBase.at(ichannel)->GetX();
		double* y_val = gRateGraphsOverBase.at(ichannel)->GetY();

		std::vector<std::pair<double,double> > val_med;
		for (int ipoint=0; ipoint<gRateGraphsOverBase.at(ichannel)->GetN(); ++ipoint) {
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
			for (int ipoint=0;ipoint<int(points_per_x.size());++ipoint){
				auto iterator = std::next(points_per_x.begin(),ipoint);
				// std::cout << std::dec << "\t\t" << ipoint << std::endl;
				// std::cout << std::dec << "\t\t" << ipoint << " " << *((*std::prev(iterator,1)).second.begin()) << " " << get_med((*std::prev(iterator,1)).second) << std::endl;
				gDiffRateGraphsOverBase.at(ichannel)->SetPoint(
					gDiffRateGraphsOverBase.at(ichannel)->GetN(),
					(*iterator).first,
					(
						get_med((*iterator).second)
					)
				);
			}
			break;
			case 1:
			for (int ipoint=2;ipoint<int(points_per_x.size())-2;++ipoint){
				auto iterator = std::next(points_per_x.begin(),ipoint);
				// std::cout << std::dec << "\t\t" << ipoint << std::endl;
				// std::cout << std::dec << "\t\t" << ipoint << " " << *((*std::prev(iterator,1)).second.begin()) << " " << get_med((*std::prev(iterator,1)).second) << std::endl;
				gDiffRateGraphsOverBase.at(ichannel)->SetPoint(
					gDiffRateGraphsOverBase.at(ichannel)->GetN(),
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
			for (int ipoint=0;ipoint<int(points_per_x.size())-1;++ipoint){
				auto iterator = std::next(points_per_x.begin(),ipoint);
				// std::cout << std::dec << "\t\t" << ipoint << std::endl;
				// std::cout << std::dec << "\t\t" << ipoint << " " << *((*std::prev(iterator,1)).second.begin()) << " " << get_med((*std::prev(iterator,1)).second) << std::endl;
				gDiffRateGraphsOverBase.at(ichannel)->SetPoint(
					gDiffRateGraphsOverBase.at(ichannel)->GetN(),
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
				double* x_values = gDiffRateGraphsOverBase.at(ichannel)->GetX();
				double* y_values = gDiffRateGraphsOverBase.at(ichannel)->GetY();
				int number = gDiffRateGraphsOverBase.at(ichannel)->GetN();
				gDiffRateGraphsOverBase.at(ichannel)->Set(0);
				for (int ipoint=1; ipoint<number-1; ipoint++) {
					gDiffRateGraphsOverBase.at(ichannel)->SetPoint(
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
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		if(fbaseline.at(ichannel)==0) continue;
		for(int ipoint=0; ipoint<gRateGraphs.at(ichannel)->GetN(); ++ipoint){
			if(gRateGraphs.at(ichannel)->GetX()[ipoint]<fbaseline.at(ichannel)) continue;
			gRateGraphsOverBase.at(ichannel)->SetPoint(
				gRateGraphsOverBase.at(ichannel)->GetN(),
				Thr_DtomV((gRateGraphs.at(ichannel)->GetX()[ipoint])-fbaseline.at(ichannel)),
				gRateGraphs.at(ichannel)->GetY()[ipoint]
			);
		}
	}
}

void dirich::AnalyzeBaseline(){
	AnalyzeBaseline(50000);
}
void dirich::AnalyzeBaseline(uint32_t NoiseThreshold)
{
	for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		int NrBins=gRateGraphs.at(ichannel)->GetN();
		int noiseedgeleft=-1;
		int noiseedgeright=-1;
		double max=0;
		double max_x=0;

		for (int ibin=0; ibin<NrBins-1; ibin++) {
			if( 
				(noiseedgeleft==-1) &&
				(gRateGraphs.at(ichannel)->GetY()[ibin] < NoiseThreshold) &&
				(gRateGraphs.at(ichannel)->GetY()[ibin+1] >= NoiseThreshold)
			){
				// linear interpolation
				double y1=gRateGraphs.at(ichannel)->GetY()[ibin];
				double y2=gRateGraphs.at(ichannel)->GetY()[ibin+1];
				double x1=gRateGraphs.at(ichannel)->GetX()[ibin];
				double x2=gRateGraphs.at(ichannel)->GetX()[ibin+1];
				noiseedgeleft=x1+(NoiseThreshold-y1)/(y2-y1)*(x2-x1);
			}

			if( 
				(noiseedgeright==-1) &&
				(gRateGraphs.at(ichannel)->GetY()[NrBins-ibin-1] < NoiseThreshold) &&
				(gRateGraphs.at(ichannel)->GetY()[NrBins-ibin-2] >= NoiseThreshold)
			){
				// linear interpolation
				double y1=gRateGraphs.at(ichannel)->GetY()[NrBins-ibin-2];
				double y2=gRateGraphs.at(ichannel)->GetY()[NrBins-ibin-1];
				double x1=gRateGraphs.at(ichannel)->GetX()[NrBins-ibin-2];
				double x2=gRateGraphs.at(ichannel)->GetX()[NrBins-ibin-1];
				noiseedgeright=x2-(NoiseThreshold-y2)/(y1-y2)*(x2-x1);
			}
			if(gRateGraphs.at(ichannel)->GetY()[ibin]>max){
				max=gRateGraphs.at(ichannel)->GetY()[ibin];
				max_x=gRateGraphs.at(ichannel)->GetX()[ibin];
			}
			// if(gRateGraphs.at(ichannel)->GetY()[NrBins-ibin-1]>max){
			// 	max=gRateGraphs.at(ichannel)->GetY()[NrBins-ibin-1];
			// 	max_x=gRateGraphs.at(ichannel)->GetX()[NrBins-ibin-1];
			// }
			if(noiseedgeleft!=-1 && noiseedgeright!=-1) break;
		}
		if(noiseedgeleft==-1 || noiseedgeright==-1){
			if(max==0 || NrBins<2){
				fbaseline_old.at(ichannel) = fbaseline.at(ichannel);
				fbaseline.at(ichannel) = 0;
				fnoisewidth_old.at(ichannel) = fnoisewidth.at(ichannel);
				fnoisewidth.at(ichannel) = 0;
				std::cout 
					<< "No baseline found for ichannel " 
					<< ichannel << "on dirich 0x" << std::hex << gBoardAddress 
					<< std::dec << std::endl;
			}
			else{
				fbaseline_old.at(ichannel) = fbaseline.at(ichannel);
				fbaseline.at(ichannel) = max_x;
				fnoisewidth_old.at(ichannel) = fnoisewidth.at(ichannel);
				fnoisewidth.at(ichannel) = 
					2*(
						abs(gRateGraphs.at(ichannel)->GetX()[0]-gRateGraphs.at(ichannel)->GetX()[1])
					);
			}
		}
		else {
			fbaseline_old.at(ichannel) = fbaseline.at(ichannel);
			fbaseline.at(ichannel)=(noiseedgeleft+noiseedgeright) / 2.;
			fnoisewidth_old.at(ichannel) = fnoisewidth.at(ichannel);
			fnoisewidth.at(ichannel)=(noiseedgeright-noiseedgeleft);
			TLine* baseline_line = new TLine(fbaseline.at(ichannel),
				gRateGraphs.at(ichannel)->GetHistogram()->GetMinimum(),
				fbaseline.at(ichannel),gRateGraphs.at(ichannel)->GetHistogram()->GetMaximum()
			);
			baseline_line->SetLineWidth(2);
			baseline_line->SetLineColor(kRed);
			gRateGraphs.at(ichannel)->GetListOfFunctions()->Add(baseline_line);
		} 
	}
}
