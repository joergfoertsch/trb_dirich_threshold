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

#include <random>
#include <math.h>

//**************************************
//dirich handling routines
//**************************************
// Channel Numbers:
// 0-31: TDC input channels (same index as for threshold setting)

const size_t BUFFER_SIZE4mb = 4194304;  /* 4MByte */
static uint32_t buffer4mb[4194304];


static std::mt19937_64 rnd;

#ifndef NCH
  const int NRCHANNELS = 32;       //Nr of TDC channels in dirich
  // const int NRCHANNELS = 4;       //Nr of TDC channels in dirich
  #define NCH
#endif
// const int OFFTHRESH =  28000;   //Value to switch off channel
const int OFFTHRESH =  1000;   //Value to switch off channel
// const int THRESHDELAY = 1000;  //Delay [mus] for thresh change to succeed
const int THRESHDELAY = 100000;  //Delay [mus] for thresh change to succeed
const int MAXWIDTH = 1000; //Maximal width of Noisepeak in DAC units
const int MINTHRESHOLD = 23000;  //Minimum Threshold to search baseline in
const int MAXTHRESHOLD = 34000;  //Maximum Threshold to search baseline in

TMutex TRBAccessMutex;
TMutex GetRateMutex;

// bool find_min_wo_zero(uint16_t i, uint16_t j) { return (i!=0 && i<j); }

class dirich 
{

private:
// public:
  uint16_t gBoardAddress;      //Board gBoardAddress
  uint64_t gBoardUID;          //UID of Board


  std::array<uint16_t,NRCHANNELS> fbaseline;
  std::array<uint16_t,NRCHANNELS> fbaseline_old;
  std::array<uint16_t,NRCHANNELS> fnoisewidth;
  std::array<uint16_t,NRCHANNELS> fnoisewidth_old;
  std::array<double,NRCHANNELS>   fthresholdmV;
      
  std::array<uint16_t,NRCHANNELS> fsim_baseline;
  std::array<uint16_t,NRCHANNELS> fsim_noise_sigma;
  std::array<uint16_t,NRCHANNELS> fsim_noise_width;
  std::array<uint16_t,NRCHANNELS> fsim_current_threshold;
  std::array<uint16_t,NRCHANNELS> fsim_singlephotonpeakposition;
  std::array<TH1*,NRCHANNELS> fsim_distro;

  int ReadSingleThreshold(uint8_t channel, uint16_t &thrvalue);
  int WriteSingleThreshold(uint8_t channel, uint16_t thrvalue);
  int ReadThresholds(uint16_t* thrarray);
  int WriteThresholds(uint16_t* thrarray);
  int WriteThresholds(uint16_t thrvalue);
  int ReadSingleScaler(uint8_t channel, uint32_t &scalervalue, std::chrono::high_resolution_clock::time_point& access_time);
  int ReadScalers(uint32_t* scalervalues, std::chrono::high_resolution_clock::time_point& access_time);
  // int GetRates(uint32_t* ratevalues, double delay=1);

  bool fsimulate = true;
  bool fjans_readout = false;
  int  gdirichver = 3;

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
  // void SetThresholdsmV(double* thrinmV);
  // void SetThresholdsmV(double thrinmV);
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
  
  double  GetSingleRate(double delay=1, uint8_t channel=0); //rates from scaler
  double* GetRates(double delay=1); //rates from scaler
  
  uint16_t ReadSingleThreshold(uint8_t channel){if(channel>=NRCHANNELS){std::cerr << "Channel: " << channel << " not specified" << std::endl; return 0;} uint16_t thrarray[32]={0}; if(ReadThresholds(thrarray)!=-1) return thrarray[channel]; else{std::cerr << "Threshold could not be read" << std::endl; return 0;}};//read threshold from !dirich!
  uint16_t* ReadThresholds() {uint16_t* thrarray = (uint16_t*) calloc(NRCHANNELS, sizeof(uint16_t*)); if(ReadThresholds(thrarray)!=-1) return thrarray; else{uint16_t* errthrarray = (uint16_t*) calloc(NRCHANNELS, sizeof(uint16_t*)); return errthrarray;}};
  
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
  void DoBaselineScan         ( );
  void DoBaselineScan         ( uint32_t SearchedNoise, uint16_t MaxStepSize, double MeasureTime, int NrPasses);
  void AnalyzeBaseline        ( );
  void AnalyzeBaseline        ( uint32_t NoiseThreshold);
  void DoThreshScan           ( );
  void DoThreshScan           ( uint8_t FirstChannel, uint8_t LastChannel, uint16_t FromThr, uint16_t ToThr, double MeasureTime, uint16_t StepSize, int NrPasses); 
  void DoThreshSearch         ( );
  void DoThreshSearch         ( double Perc, bool SPP_SPV /*1==SPP*/, double MeasureTime, int16_t StepSize, int NrPasses);
  void DoThreshScanOverBase   ( );
  void DoThreshScanOverBase   ( uint8_t FirstChannel, uint8_t LastChannel, double FromThrmV, double ToThrmV, double MeasureTime, double StepSize, int NrPasses);
  void MakeGraphsOverBase     ( );
  void MakeGraphsOverBase     ( uint16_t ToThr);
  void MakeDiffGraphsOverBase ( );
  // void MakeDiffGraphsOverBase ( uint16_t FromThr, uint16_t ToThr);
  void FindMinThreshScanOverBase(double gThreshold_finding_method);
  
  bool IsSim (){return fsimulate;}
  bool IsJansReadout (){return fjans_readout;}
  int  WhichDirichVersion (){return gdirichver;}

  // threshold visualization items
  std::array<TGraph*,NRCHANNELS> gRateGraphs;
  std::array<TGraph*,NRCHANNELS> gRateGraphsOverBase;
  std::array<TGraph*,NRCHANNELS> gDiffRateGraphsOverBase;
  // void ClearGraphs();


  // settings for threshold measurement
  double 	 gMeasureTime; 	 //Time [s] to determine rate
  int    	 gLowerEdge;   	 //start of thresholdscan
  int    	 gUpperEdge;   	 //end of thresholdscan
  int    	 gStepsize;    	 //end of thresholdscan
  int    	 gNrPasses;    	 //Nr of passes to scan all channels

  double   gMeasureTime_over;   //Time [s] to determine rate
  int      gLowerEdge_over;     //end of thresholdscan
  int      gUpperEdge_over;     //end of thresholdscan
  int      gStepsize_over;      //end of thresholdscan
  int      gNrPasses_over;      //Nr of passes to scan all channels

  double   gThreshold_finding_method; //Method to find perfect threshold: 
                                      //0: searches for the minimum in the differentiated spectrum or for the minimal gradient 
                                      //0<value<5: tries to find peak and sigma of the single photon distribution and sets the threshold to value*sigma
                                      //5<value<100: tries to find the single photon peak and sets the threshold to value% of the spp-position

};


dirich::dirich()
{
  dirich(0);
}

dirich::dirich(uint16_t BoardAddress)
  {
  
  int ret=0;
  for(int i=0;i<100;++i){
    TRBAccessMutex.Lock();
    ret=trb_read_uid(BoardAddress, buffer4mb, BUFFER_SIZE4mb);
    TRBAccessMutex.UnLock();
    if(ret>0) break;
  }
  if(ret<=0){
    std::cerr << "No DiRICH found with Address:" << BoardAddress << "\nNot adding DiRICH" << std::endl;
    return;
  }
  else if(ret!=4){
    std::cerr << "Too many DiRICH found with Address:" << BoardAddress << "(Amount: "<< ret << ")\nNot adding DiRICH" << std::endl;
    return;
  }
  else{
    gBoardAddress = buffer4mb[3];
    uint64_t temp_store = buffer4mb[0];
    uint64_t temp_store2 = buffer4mb[1];
    gBoardUID = temp_store << 32 | temp_store2;
  }
  // std::cout << ret << " " << std::hex << buffer4mb[0] << " " << buffer4mb[1] << " " << gBoardUID << " " << gBoardAddress << std::endl;

  uint32_t cmd = 0x0 << 20 | 0xff << 24;
  uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0x10001};

  for(int failed=0;failed<100;++failed){
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xd400,0,c,18);
    TRBAccessMutex.UnLock();
    if(ret!=-1) break;
    usleep(1000);
  }
  // std::cout << ret << std::endl;

  std::array<uint32_t,18> ret_c;
  TRBAccessMutex.Lock();
  ret=trb_register_read(gBoardAddress,0xd412,ret_c.data(),18);
  TRBAccessMutex.UnLock();
  if(ret!=2 || (ret_c.at(1) & 0xff00)!=0x100){
    std::cerr << "No DiRICH2 Threshold FPGA of newest version detected\nNot adding DiRICH" << std::endl;
    return;
  }
  else{
    gdirichver = 3;
  }

  // gBoardAddress=BoardAddress;
  gMeasureTime=.3;
  gLowerEdge=28000;
  gUpperEdge=32000;
  // gStepsize=10;
  gStepsize=50;
  gNrPasses=2;
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
    // if( BoardAddress == 0x1235 ||  BoardAddress == 0x1231 ||  BoardAddress == 0x1217 ||  BoardAddress == 0x1233 ||  BoardAddress == 0x1238 ||  BoardAddress == 0x1214 ||  BoardAddress == 0x1237){
    //   fjans_readout=true;
    // }
    // else{
    //   fjans_readout=false;
    // }
    // int ret=0;
    // //check if TDC is implemented!
    // uint32_t temp_tdc_setting[2];
    // TRBAccessMutex.Lock();
    // ret=trb_register_read(BoardAddress, 0xc802, temp_tdc_setting, 2); //get TDC-settings
    // TRBAccessMutex.UnLock();
    // std::cout << "jans readout check phase 1: ret=" << std::dec << ret << " temp_tdc_setting[0]=" << std::hex << temp_tdc_setting[0] << " temp_tdc_setting[2]=" << temp_tdc_setting[1] << std::endl; 
    // if(ret!=2 || temp_tdc_setting[0]!=BoardAddress){
    //   std::cout << "jans readout check phase 1: is true" << std::endl;
    //   fjans_readout=true;
    // }

    // if(fjans_readout!=true){ 
    //   TRBAccessMutex.Lock();
    //   ret=trb_register_write(BoardAddress, 0xc802, temp_tdc_setting[1]); //reset TDC-settings
    //   TRBAccessMutex.UnLock();
    //   std::cout << "jans readout check phase 2: ret=" << std::dec << ret << " temp_tdc_setting[0]=" << std::hex << temp_tdc_setting[0] << " temp_tdc_setting[2]=" << temp_tdc_setting[1] << std::endl; 
    //   if(ret==-1){
    //     std::cout << "jans readout check phase 2: is true" << std::endl;
    //     fjans_readout=true;
    //   }
    // }

  //   // if(fjans_readout!=true){ 
  //   //   ret=-1;
  //   //   uint32_t scaler1[NRCHANNELS];
  //   //   // std::cout << "GetSingleRate" << scaler1 << " " << scaler2 << " " << ratevalues << std::endl;
  //   //   GetRateMutex.Lock();
  //   //   // auto start = std::chrono::high_resolution_clock::now();
  //   //   std::chrono::high_resolution_clock::time_point start1;
  //   //   for(int iterator=0;iterator<100;++iterator){
  //   //     if(ret>=0) break;
  //   //     ret=ReadScalers(scaler1,start1);
  //   //     // std::cout << iterator << "\t";
  //   //     if(iterator==99)fjans_readout=true;
  //   //   }
  //   //   GetRateMutex.UnLock();
  //   // }
  //   if(fjans_readout==true){ 
  //     std::cout << "using jans scaler implementation for dirich "  << std::hex << BoardAddress << std::dec << std::endl;
  //   }

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

int dirich::ReadSingleThreshold(uint8_t channel, uint16_t &thrvalue)
{
  if (channel>NRCHANNELS-1)
    return -1;
  int ret;
  if(gBoardAddress<=0x1210){
//   int reg=0xa000+31-channel; old firwmare
    int reg=0xa000+channel; //new firwmare 
    uint32_t buffer[2];
    TRBAccessMutex.Lock();
    ret=trb_register_read(gBoardAddress,reg, buffer, 2);
    TRBAccessMutex.UnLock();

    if((gBoardAddress != buffer[0]) || (ret != 2)) return -1;
    
    thrvalue=buffer[1];  
    return 0;
  }
  else{
    // uint8_t real_channel = channel%16+16;
    uint8_t real_channel = channel%16;
    uint32_t cmd = 0x0 << 20 | real_channel << 24 | thrvalue << 0;

    uint32_t c[] = {cmd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(uint32_t)channel/16+1,0x10001}; //evtl. sind auch mehrere Kanäle auf einmal lesbar.
    for(int failed=0;failed<100;++failed){
      TRBAccessMutex.Lock();
      ret=trb_register_write_mem(gBoardAddress,0xd400,0,c,18);
      TRBAccessMutex.UnLock();
      if(ret!=-1) break;
      usleep(1000);
    }
    uint32_t ret_c[18];
    TRBAccessMutex.Lock();
    ret=trb_register_read(gBoardAddress,0xd412,ret_c,18);
    TRBAccessMutex.UnLock();
    // if((gBoardAddress != ret_c[0]) || (ret != 2)) return -1;
    
    thrvalue=ret_c[1];  
    return 0;
  }
}

int dirich::WriteSingleThreshold(uint8_t channel, uint16_t thrvalue) 
{
  if (channel>NRCHANNELS-1)
    return -1;
  if(fsimulate)
    fsim_current_threshold.at(channel) = thrvalue;
  else{
    int ret=0;

    if(gdirichver==1){
  //   int reg=0xa000+31-channel; old firwmare
      int reg=0xa000+channel; //new firwmare 
      TRBAccessMutex.Lock();
      ret=trb_register_write(gBoardAddress, reg, (uint32_t)thrvalue);
      TRBAccessMutex.UnLock();
      // std::cout << ret << std::endl;
      return ret;
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
      uint32_t ret_c[18];
      TRBAccessMutex.Lock();
      ret=trb_register_read(gBoardAddress,0xd412,ret_c,18);
      TRBAccessMutex.UnLock();

      return ret;
    }
  }
}

int dirich::ReadThresholds(uint16_t* thrarray)
{
  int ret;
  uint16_t reg=0xa000+32-NRCHANNELS;
  uint32_t buffer[NRCHANNELS+1];
  
  TRBAccessMutex.Lock();
  ret=trb_register_read_mem(gBoardAddress,reg,0,NRCHANNELS,buffer,NRCHANNELS+1);
  TRBAccessMutex.UnLock();
  if (ret != NRCHANNELS+1) return -1;

  for (int i=0;i<NRCHANNELS;i++) {
    thrarray[i]=buffer[NRCHANNELS-i-1];
  }
  return 0;
}

int dirich::WriteThresholds(uint16_t* thrarray)
{
  if(fsimulate){
    for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
      fsim_current_threshold.at(ichannel) = thrarray[ichannel];
    }
  }
  else{
    int ret;
    std::array<uint32_t,18> c1;
    std::array<uint32_t,18> c2;
    if(gdirichver==1){
  	  uint16_t reg=0xa000+32-NRCHANNELS;
  	  uint32_t buffer[NRCHANNELS];
  	  for (int i=0;i<NRCHANNELS;i++) {
  	    buffer[i]=thrarray[NRCHANNELS-i-1];
  	  }
  	  TRBAccessMutex.Lock();
  	  ret=trb_register_write_mem(gBoardAddress,reg,0,buffer,NRCHANNELS);
  	  TRBAccessMutex.UnLock();
  	  return ret;
    }
    else{
    std::array<uint32_t,32> cmd;
    	for(int channel=0;channel<NRCHANNELS;++channel){
        // uint8_t real_channel = channel%16+16;
      	uint8_t real_channel = channel%16+16*abs(gdirichver-3);
      	cmd.at(channel) = 0x8 << 20 | real_channel << 24 | thrarray[channel] <<0;
        // std::cout << std::hex << thrarray[channel] << " threhsold for channel " << channel << " gives: " << cmd.at(channel) << std::endl;
    	}

      c1 = {{cmd.at(0),cmd.at(1),cmd.at(2),cmd.at(3),cmd.at(4),cmd.at(5),cmd.at(6),cmd.at(7),cmd.at(8),cmd.at(9),cmd.at(10),cmd.at(11),cmd.at(12),cmd.at(13),cmd.at(14),cmd.at(15),(uint32_t)1,0x10010}};
      for(int failed=0;failed<100;++failed){
        TRBAccessMutex.Lock();
        ret=trb_register_write_mem(gBoardAddress,0xd400,0,c1.data(),18);
        TRBAccessMutex.UnLock();

        if(ret!=-1) break;
        usleep(1000);
      }

    	std::array<uint32_t,18> ret_c1;
      TRBAccessMutex.Lock();
      if(ret!=-1) ret=trb_register_read(gBoardAddress,0xd412,ret_c1.data(),18);
      TRBAccessMutex.UnLock();

      if(ret==-1) return ret;

      // std::cout << "done setting first fpga" << std::endl;
      // for(auto& c_iterator : c1){
      //   std::cout << c_iterator << " ";
      // }
      // std::cout << std::endl;
      // for(auto& c_iterator : ret_c1){
      //   std::cout << c_iterator << " ";
      // }
      // std::cout << std::endl;

      c2 = {{cmd.at(16+0),cmd.at(16+1),cmd.at(16+2),cmd.at(16+3),cmd.at(16+4),cmd.at(16+5),cmd.at(16+6),cmd.at(16+7),cmd.at(16+8),cmd.at(16+9),cmd.at(16+10),cmd.at(16+11),cmd.at(16+12),cmd.at(16+13),cmd.at(16+14),cmd.at(16+15),(uint32_t)2,0x10010}};
  		for(int failed=0;failed<100;++failed){
        TRBAccessMutex.Lock();
        ret=trb_register_write_mem(gBoardAddress,0xd400,0,c2.data(),18);
        TRBAccessMutex.UnLock();
        if(ret!=-1) break;
        usleep(1000);
      }

    	std::array<uint32_t,18> ret_c2;
      TRBAccessMutex.Lock();
      if(ret!=-1) ret=trb_register_read(gBoardAddress,0xd412,ret_c2.data(),18);
      TRBAccessMutex.UnLock();

      // if(ret==-1) return ret;
      // std::cout << "done setting second fpga" << std::endl;
      // for(auto& c_iterator : c2){
        // std::cout << c_iterator << " ";
      // }
      // std::cout << std::endl;    
      // if(ret_c1!=c1 || ret_c2!=c2){
      //   std::cout << "return values not equal!" << std::endl;
      //   ret=-1;
      // } 
      return ret;
    }
  }
}

int dirich::WriteThresholds(uint16_t thrvalue)
{
  int ret;
  uint16_t thrarray[NRCHANNELS];
  for (int i=0; i<NRCHANNELS; i++)
    thrarray[i]=thrvalue;
  ret=WriteThresholds(thrarray);
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
  int newthreshold=baseline+thrinD;

  WriteSingleThreshold(channel,newthreshold);
  fthresholdmV[channel]=thrinmV;

  // printf("Setting threshold for channel %i \n",channel);
  // printf("Baseline: %i \n",baseline);
  // printf("Threshold %3.0f in digits:  %i \n",thrinmV, thrinD);
  // printf("New threshold value %i set \n\n",newthreshold);
  // DrawGraphs(dirichptr,cc.at(dirichptr.first),Form("Threshold: %3.0f mV",thrinmV));
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
  // std::cout << std::hex << buffer[0] << " " << gBoardAddress << " " << fjans_readout << std::dec << " " << ret <<  std::endl;
  // if (ret != NRCHANNELS+1) 
  // if ( gBoardAddress != buffer[0]) 
  if ( (ret == NRCHANNELS+1) && (fjans_readout || gBoardAddress == (buffer[0] & 0xffff)) ){
    for (int i=0;i<NRCHANNELS;i++){
      scalervalues[i]=buffer[i+1] & 0x7fffffff;
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
      // GetRateMutex.Lock();
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
    // GetRateMutex.UnLock();
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
    double rate=scaler2<scaler1?
      (2<<23)+scaler2-scaler1 : scaler2-scaler1;
    rate=1.*rate/exactdelay1;
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
      double rate=scaler2[i]<scaler1[i]?
  			(2<<23)+scaler2[i]-scaler1[i] : scaler2[i]-scaler1[i];
      rate=1.*rate/exactdelay1;
      ratevalues[i]=rate;
    }
    return ratevalues;
  }
}

void dirich::DoBaselineScan(){
	DoBaselineScan(50000, 250, .3, 4);
}
void dirich::DoBaselineScan(uint32_t SearchedNoise, uint16_t StartStepSize, double MeasureTime, int NrPasses)
{
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
    gRateGraphs[ichannel]->Set(0);
  }  
  int ret=0;
  if(fjans_readout){
    uint32_t scaler_switch[] = {0xffffffff};
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching on scaler for this dirich .... only needed when using jan's readout
    TRBAccessMutex.UnLock();
    if(ret<0) std::cerr << "Error switching on scalers" << std::endl;
  }

  int number_of_checks=2;

  std::array<uint16_t, NRCHANNELS> low_edge;
  low_edge.fill(0);
  std::array<uint16_t, NRCHANNELS> high_edge;
  high_edge.fill(0);  

  for(int ipass=0;ipass<NrPasses;++ipass){
    // std::cout  << std::dec << "Baselinescan @ " << 50/NrPasses*ipass << "\%" << std::endl;
    int n_of_iterations=0;
    std::array<int16_t, NRCHANNELS> step_size;
    std::array<uint16_t, NRCHANNELS> threshold_value;
    std::array<double, NRCHANNELS> rate;
    std::array<double, NRCHANNELS> old_rate;
    std::array<int, NRCHANNELS> status;
    status.fill(0);
    rate.fill(0);
    old_rate.fill(0);    
    for(int ichannel = 0; ichannel<NRCHANNELS; ++ichannel ){
      step_size.at(ichannel)= ichannel%NrPasses==ipass ? StartStepSize : 0;
      threshold_value.at(ichannel)= ichannel%NrPasses==ipass ? MINTHRESHOLD : OFFTHRESH;
      status.at(ichannel)= ichannel%NrPasses==ipass ? 0 : 4*number_of_checks-2;
    }
    while(!std::all_of(status.cbegin(), status.cend(), [number_of_checks](int i){ return i==4*number_of_checks-2;}) 
      && !std::all_of(step_size.cbegin(), step_size.cend(), [](int i){ return i==0;}) 
      && n_of_iterations<500)
    {
    // while(!std::all_of(step_size.cbegin(), step_size.cend(), [](int i){ return i==0;}) && n_of_iterations<200){
      std::cout << std::dec << n_of_iterations << std::endl;
      for(auto& step_size_ch : step_size){
        std::cout << std::dec << step_size_ch << "\t";
      }
      std::cout << std::dec << std::endl;
      for(auto& threshold_value_ch : threshold_value){
        std::cout << std::dec << threshold_value_ch << "\t";
      }
      std::cout << std::dec << std::endl; 
      for(auto& rate_ch : rate){
        std::cout << std::dec << rate_ch << "\t";
      }
      std::cout << std::dec << std::endl;      
      for(auto& old_rate_ch : old_rate){
        std::cout << std::dec << old_rate_ch << "\t";
      }
      std::cout << std::dec << std::endl; 
      for(auto& status_ch : status){
        std::cout << std::dec << status_ch << "\t";
      }
      std::cout << std::dec << std::endl;             
      std::cout << std::dec << std::endl;             
      std::cout << std::dec << std::endl;             

      ++n_of_iterations;
      ret=WriteThresholds(threshold_value.data());
      // if(n_of_iterations==1)usleep(1000000);
      usleep(THRESHDELAY);
      double* temp_rate;
      temp_rate = GetRates(MeasureTime);
      memcpy(rate.data(),temp_rate,NRCHANNELS*sizeof(double));

      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
        gRateGraphs.at(ichannel)->SetPoint(gRateGraphs.at(ichannel)->GetN(),1.*threshold_value.at(ichannel),1.*rate.at(ichannel));
        if(status.at(ichannel)<number_of_checks){
          if(status.at(ichannel)%2==0 && old_rate.at(ichannel) < SearchedNoise && rate.at(ichannel) > SearchedNoise){
            status.at(ichannel)++;
            threshold_value.at(ichannel)-=2*step_size.at(ichannel);
          }
          else if(status.at(ichannel)%2==1 && old_rate.at(ichannel) > SearchedNoise && rate.at(ichannel) < SearchedNoise){
            status.at(ichannel)++;
            // threshold_value.at(ichannel)-=step_size.at(ichannel);
          }
          else{
            status.at(ichannel)=0;
          }
          if(status.at(ichannel)==number_of_checks*2-1){
            step_size.at(ichannel)/=3;
            if(step_size.at(ichannel)!=0){
              status.at(ichannel)=0;
            }
            else{
              step_size.at(ichannel)= ichannel%NrPasses==ipass ? -StartStepSize : 0;
              threshold_value.at(ichannel)= ichannel%NrPasses==ipass ? MAXTHRESHOLD : OFFTHRESH;
              rate.at(ichannel)==0;
              low_edge.at(ichannel)=threshold_value.at(ichannel);
            }
          }
        }
        if(status.at(ichannel)<4*number_of_checks-2){
          if(status.at(ichannel)%2==0 && old_rate.at(ichannel) < SearchedNoise && rate.at(ichannel) > SearchedNoise){
            status.at(ichannel)++;
            threshold_value.at(ichannel)-=2*step_size.at(ichannel);
          }
          else if(status.at(ichannel)%2==1 && old_rate.at(ichannel) > SearchedNoise && rate.at(ichannel) < SearchedNoise){
            status.at(ichannel)++;
            // threshold_value.at(ichannel)-=step_size.at(ichannel);
          }
          else{
            status.at(ichannel)=0;
          }
          if(status.at(ichannel)==number_of_checks*2-1){
            step_size.at(ichannel)/=3;
            if(step_size.at(ichannel)!=0){
              status.at(ichannel)=0;
            }
            else{
              high_edge.at(ichannel)=threshold_value.at(ichannel);
            }
          }
        }
        threshold_value.at(ichannel)+=step_size.at(ichannel);
      }
      old_rate=rate;
    }
  }

  if(fjans_readout){
    usleep(THRESHDELAY);
    uint32_t scaler_switch[] = {0x0};
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching off scaler for this dirich .... only needed when using jan's readout
    TRBAccessMutex.UnLock();
    if(ret<0) std::cerr << "Error switching off scalers" << std::endl;
  }
  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
    // std::cout << std::dec << "high_edge\t" << high_edge.at(ichannel) << "\tlow_edge\t" << low_edge.at(ichannel) << std::endl;
    fnoisewidth_old.at(ichannel) = fnoisewidth.at(ichannel);
    fnoisewidth.at(ichannel) = high_edge.at(ichannel)-low_edge.at(ichannel);
    fbaseline_old.at(ichannel) = fbaseline.at(ichannel);
    fbaseline.at(ichannel) = (high_edge.at(ichannel)+low_edge.at(ichannel))/2;
    // std::cout << std::dec << "fnoisewidth_old" << "\t" << fnoisewidth_old.at(ichannel) << "\t" << "fnoisewidth" << "\t" << fnoisewidth.at(ichannel) << "\t" << "fbaseline_old" << "\t" << fbaseline_old.at(ichannel) << "\t" << "fbaseline" << "\t" << fbaseline.at(ichannel) << std::endl;
  }
  for (int ichannel=0; ichannel<NRCHANNELS; ++ichannel){
    gRateGraphs.at(ichannel)->Sort();
  }  
}
  
void dirich::DoThreshScan(){
	DoThreshScan(0, NRCHANNELS, gLowerEdge, gUpperEdge, gMeasureTime, gStepsize, gNrPasses);
}
void dirich::DoThreshScan(uint8_t FirstChannel, uint8_t LastChannel, uint16_t FromThr, uint16_t ToThr, double MeasureTime, uint16_t StepSize, int NrPasses)
{
  // std::cout << std::dec << (int)FirstChannel << " " << (int)LastChannel << " " << FromThr << " " << ToThr << " " << MeasureTime << " " << StepSize << " " << NrPasses << std::endl;
  int ret;
  if(fjans_readout){
    uint32_t scaler_switch[] = {0xffffffff};
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching on scaler for this dirich .... only needed when using jan's readout
    TRBAccessMutex.UnLock();
    if(ret<0) std::cerr << "Error switching on scalers" << std::endl;
  }
  // std::cout << "Set0" << std::endl;
  // TGraph* temp[NRCHANNELS];
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
    // if(gRateGraphs[ichannel]){
      // std::cout << "deleting gRateGraphs[" << ichannel << "]" << std::endl;
      // delete gRateGraphs[ichannel];
      // std::cout << "deleted" << std::endl;
      // gRateGraphs[ichannel]=new TGraph();
    // }
    // gRateGraphs[ichannel]->SetTitle(Form("Rate graph of dirich 0x%x's channel %i;Threshold;Rate",gBoardAddress,ichannel));
    // gRateGraphs[ichannel]->SetName(Form("Rate graph of dirich 0x%x's channel %i",gBoardAddress,ichannel));
    // temp[ichannel] = new TGraph();
    gRateGraphs[ichannel]->Set(0);
  }
  for(int ipass=0;ipass<NrPasses;++ipass){
    // std::cout << "Pass Nr. " << ipass+1 << " of " << NrPasses << std::endl;
    ret=WriteThresholds(OFFTHRESH);
    usleep(THRESHDELAY);
    ret=WriteThresholds(OFFTHRESH);
    usleep(THRESHDELAY);
    for (int thresh=FromThr; thresh<=ToThr; thresh+=StepSize){
      // std::cout << std::dec << ((thresh-FromThr)/StepSize+1)*(ipass+1) << "\t" << ((ToThr-FromThr)/StepSize+1)*(NrPasses+1) << "\t" << int(1.*((thresh-FromThr)/StepSize+1)*(ipass+1)/(((ToThr-FromThr)/StepSize+1)*(NrPasses+1))*100) << "\t" << int(1.*((thresh-FromThr)/StepSize+1)*(ipass+1)/(((ToThr-FromThr)/StepSize+1)*(NrPasses+1))*100)%20 << std::endl;
      // if(int(((thresh-FromThr)+(ToThr-FromThr)*(ipass))/((ToThr-FromThr)*(NrPasses))*100.)%20==0) std::cout << "Baselinescan of dirich 0x" << std::hex << gBoardAddress << std::dec <<" is @ " << int(((thresh-FromThr)+(ToThr-FromThr)*(ipass))/((ToThr-FromThr)*(NrPasses))*100.) << "%" << std::endl;
      for (int ichannel=FirstChannel+ipass; ichannel<=LastChannel; ichannel+=NrPasses){
        ret=WriteSingleThreshold(ichannel,thresh);
      }
      usleep(THRESHDELAY);
      double* rates;
      rates = GetRates(MeasureTime);
      for (int ichannel=FirstChannel+ipass; ichannel<LastChannel; ichannel+=NrPasses){
        std::cout << std::dec << "Channel " << ichannel << " PointIndex " << gRateGraphs[ichannel]->GetN() << " thr " << thresh << " rate " << rates[ichannel] << std::endl;
        gRateGraphs[ichannel]->SetPoint(gRateGraphs[ichannel]->GetN(),1.*thresh,1.*rates[ichannel]);
        // temp[ichannel]->SetPoint(temp[ichannel]->GetN(),1.*thresh,1.*rates[ichannel]);
        // std::cout << std::dec << "2Channel " << ichannel << " PointIndex " << gRateGraphs[ichannel]->GetN() << " thr " << thresh << " rate " << rates[ichannel] << std::endl;
        // gRateGraphs[ichannel]->Sort();
      }
      std::cout << std::endl;
    }
  }
  for (int ichannel=FirstChannel; ichannel<LastChannel; ++ichannel){
    // std::cout << ichannel << " 1" << std::endl;
    // for (int ipoint=0; ipoint<temp[ichannel]->GetN(); ++ipoint){
    //   gRateGraphs[ichannel]->SetPoint(ipoint,temp[ichannel]->GetX()[ipoint],temp[ichannel]->GetY()[ipoint]); 
    // }
    // std::cout << ichannel << " 2" << std::endl;
    gRateGraphs[ichannel]->Sort();
    // std::cout << ichannel << " 3" << std::endl;
    // delete temp[ichannel];
    // std::cout << ichannel << " 4" << std::endl;
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

void dirich::DoThreshSearch(){
  if(std::all_of(fbaseline.cbegin(), fbaseline.cend(), [](int i){ return i==0; })){
    std::cerr << "dirich 0x" << std::hex << gBoardAddress << std::dec << " has no baseline yet. Please load or scan one (load_base, system_thr_scan)" << std::endl;
    return;
  }
  if(gThreshold_finding_method==0){
    DoThreshSearch(1., 0, 30., 100, 1);
  }
  else if(gThreshold_finding_method>0 && gThreshold_finding_method<5){
    DoThreshScanOverBase( 0, NRCHANNELS, 0., 999., 30., 10, 1);
    MakeDiffGraphsOverBase();
    FindMinThreshScanOverBase(gThreshold_finding_method);
    std::cout << "currently this method is not implemented" << std::endl;
  } 
  else if(gThreshold_finding_method>5 && gThreshold_finding_method<100){
    DoThreshSearch(gThreshold_finding_method, 1, 30., 400, 1);
  }  
}
void dirich::DoThreshSearch(double Perc, bool SPP_SPV /*1==SPP*/, double MeasureTime, int16_t StepSize, int NrPasses){
  int ret;
  if(fjans_readout){
    uint32_t scaler_switch[] = {0xffffffff};
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching on scaler for this dirich .... only needed when using jan's readout
    TRBAccessMutex.UnLock();
    if(ret<0) std::cerr << "Error switching on scalers" << std::endl;  
  }
    // std::cout << "Set0" << std::endl;
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
    gRateGraphsOverBase[ichannel]->Set(0);
    gDiffRateGraphsOverBase[ichannel]->Set(0);
  }

  std::array<double,NRCHANNELS> SPP_SPV_value;
  for(int ipass=0;ipass<NrPasses;++ipass){
    int n_of_iterations=0;
    std::array<uint16_t, NRCHANNELS> step_size;
    std::array<uint16_t, NRCHANNELS> threshold_value;
    std::array<std::vector<double>, NRCHANNELS> rate;
    std::array<std::vector<double>, NRCHANNELS> diff_rate;
    for(int ichannel = 0; ichannel<NRCHANNELS; ++ichannel ){
      step_size.at(ichannel)= ichannel%NrPasses==ipass ? StepSize : 0;
      // step_size.at(ichannel)= ichannel==0 ? StepSize : 0;
      // threshold_value.at(ichannel)= ichannel%NrPasses==ipass ? fbaseline.at(ichannel)+fnoisewidth.at(ichannel)/5+.5 : OFFTHRESH;
      threshold_value.at(ichannel)= ichannel%NrPasses==ipass ? fbaseline.at(ichannel)+fnoisewidth.at(ichannel)/5+.5 : OFFTHRESH;
      SPP_SPV_value.at(ichannel)=0;
    }
    while(
        n_of_iterations==0 
        || (
          !std::all_of(step_size.cbegin(), step_size.cend(), [](int i){ return i==0;}) 
          && n_of_iterations<200 
          && !(
            std::all_of(rate.cbegin(), rate.cend(), [](std::vector<double> i){ return i.back()< 3.;}) 
            && std::all_of(step_size.cbegin(), step_size.cend(), [StepSize](int i){ return i==StepSize || i==0;})
          )
        )
      ){
    
      ++n_of_iterations;
      std::cout << "n_of_iterations " << n_of_iterations << std::endl;
      ret=WriteThresholds(threshold_value.data());
      if(n_of_iterations==1)usleep(1000000);
      usleep(2*THRESHDELAY);
      double* temp_rate;
      temp_rate = GetRates(MeasureTime);
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
        if(step_size.at(ichannel)==0) continue;
        std::cout << "ichannel " << ichannel << std::endl;
        std::cout << "threshold_value.at(ichannel) " << threshold_value.at(ichannel) << "\t";
        std::cout << "temp_rate[ichannel] " << temp_rate[ichannel] << "\t";
        rate.at(ichannel).push_back(temp_rate[ichannel]);
        std::cout << "rate.at(ichannel) " << rate.at(ichannel).back() << std::endl;
        gRateGraphsOverBase.at(ichannel)->SetPoint(gRateGraphsOverBase.at(ichannel)->GetN(),Thr_DtomV(threshold_value.at(ichannel)-fbaseline.at(ichannel)),rate.at(ichannel).back());
        std::cout << "rate.at(ichannel).size() " << rate.at(ichannel).size() << std::endl;
        if(step_size.at(ichannel)==StepSize && rate.at(ichannel).size()>5){
          diff_rate.at(ichannel).push_back(-1.*rate.at(ichannel).at(rate.at(ichannel).size()-5)+16.*rate.at(ichannel).at(rate.at(ichannel).size()-4)+32.*rate.at(ichannel).at(rate.at(ichannel).size()-3)-16.*rate.at(ichannel).at(rate.at(ichannel).size()-2)+1.*rate.at(ichannel).at(rate.at(ichannel).size()-1));
          diff_rate.at(ichannel).back()/=(12.*step_size.at(ichannel)*step_size.at(ichannel));
        }
        if(rate.at(ichannel).size()>=15){
          double left_side_derivative=-1.*rate.at(ichannel).at(rate.at(ichannel).size()-15)+8.*rate.at(ichannel).at(rate.at(ichannel).size()-14)-8.*rate.at(ichannel).at(rate.at(ichannel).size()-12)+1.*rate.at(ichannel).at(rate.at(ichannel).size()-11);
          double middle_derivative=-1.*rate.at(ichannel).at(rate.at(ichannel).size()-10)+8.*rate.at(ichannel).at(rate.at(ichannel).size()-9)-8.*rate.at(ichannel).at(rate.at(ichannel).size()-7)+1.*rate.at(ichannel).at(rate.at(ichannel).size()-6);
          double right_side_derivative=-1.*rate.at(ichannel).at(rate.at(ichannel).size()-5)+8.*rate.at(ichannel).at(rate.at(ichannel).size()-4)-8.*rate.at(ichannel).at(rate.at(ichannel).size()-2)+1.*rate.at(ichannel).at(rate.at(ichannel).size()-1);
          std::cout << "threshold_value.at(ichannel) " << threshold_value.at(ichannel) << " threshold_value.at(ichannel)-13*step_size.at(ichannel)-fbaseline.at(ichannel) " << threshold_value.at(ichannel)-13*step_size.at(ichannel)-fbaseline.at(ichannel) << std::endl;
          gDiffRateGraphsOverBase.at(ichannel)->SetPoint(gDiffRateGraphsOverBase.at(ichannel)->GetN(),Thr_DtomV(threshold_value.at(ichannel)-13*step_size.at(ichannel)-fbaseline.at(ichannel)),left_side_derivative/(12.*step_size.at(ichannel)));
          std::cout << "left_side_derivative " << left_side_derivative << " middle_derivative " << middle_derivative << " right_side_derivative " << right_side_derivative << std::endl;
          if((SPP_SPV==1 && left_side_derivative<middle_derivative && right_side_derivative<middle_derivative) || (SPP_SPV==0 && left_side_derivative>middle_derivative && right_side_derivative>middle_derivative)){
            rate.at(ichannel).clear();
            if(step_size.at(ichannel)/10==0){
              gDiffRateGraphsOverBase.at(ichannel)->SetPoint(gDiffRateGraphsOverBase.at(ichannel)->GetN(),Thr_DtomV(threshold_value.at(ichannel)-8*step_size.at(ichannel)-fbaseline.at(ichannel)),middle_derivative/(12.*step_size.at(ichannel)));
              gDiffRateGraphsOverBase.at(ichannel)->SetPoint(gDiffRateGraphsOverBase.at(ichannel)->GetN(),Thr_DtomV(threshold_value.at(ichannel)-3*step_size.at(ichannel)-fbaseline.at(ichannel)),right_side_derivative/(12.*step_size.at(ichannel)));
              // SPP_SPV_value.at(ichannel)=Perc*Thr_DtomV(threshold_value.at(ichannel)-step_size.at(ichannel)*5.5-fbaseline.at(ichannel));
              SPP_SPV_value.at(ichannel)=Perc*Thr_DtomV(threshold_value.at(ichannel)-step_size.at(ichannel)*5.5-fbaseline.at(ichannel));
              threshold_value.at(ichannel)=OFFTHRESH;
            }
            else{
              threshold_value.at(ichannel)=threshold_value.at(ichannel)-step_size.at(ichannel)*6;
            }
            step_size.at(ichannel)/=10;
          }
        }
        threshold_value.at(ichannel)=threshold_value.at(ichannel)+step_size.at(ichannel);
      }
    }
    if(SPP_SPV==0){
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
        if(step_size.at(ichannel)==StepSize){
          double minimum=1E6;
          for(int diff_rate_it=0;diff_rate_it<diff_rate.at(ichannel).size();++diff_rate_it){
            if(diff_rate.at(ichannel).at(diff_rate_it)<minimum){
              minimum=diff_rate.at(ichannel).at(diff_rate_it);
              diff_rate_it*step_size.at(ichannel)+fbaseline.at(ichannel)+fnoisewidth.at(ichannel)/5+.5;
            } 
          }
        }
      }
    }    
  }

  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
    std::cout << Thr_mVtoD(SPP_SPV_value.at(ichannel))+fbaseline.at(ichannel) << "\t";
    gDiffRateGraphsOverBase.at(ichannel)->Sort();
    gRateGraphsOverBase.at(ichannel)->Sort();
    SetSingleThresholdmV(ichannel,-1*SPP_SPV_value.at(ichannel));
  }
  std::cout << std::endl;
  if(fjans_readout){
    usleep(THRESHDELAY);
    uint32_t scaler_switch[] = {0x0};
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching off scaler for this dirich .... only needed when using jan's readout
    TRBAccessMutex.UnLock();  
    if(ret<0) std::cerr << "Error switching off scalers" << std::endl;
  }
}
void dirich::DoThreshScanOverBase(){
  if(std::all_of(fbaseline.cbegin(), fbaseline.cend(), [](int i){ return i==0; })){
    std::cerr << "dirich 0x" << std::hex << gBoardAddress << std::dec << " has no baseline yet. Please load or scan one (load_base, system_thr_scan)" << std::endl;
    return;
  }
  DoThreshScanOverBase(0, NRCHANNELS, gLowerEdge_over, gUpperEdge_over, gMeasureTime_over, gStepsize_over, gNrPasses_over);
  // DoThreshScanOverBase(0, NRCHANNELS, Thr_DtomV(gUpperEdge-*min_element(fbaseline.begin(),fbaseline.end(),find_min_wo_zero)), gMeasureTime, Thr_DtomV(gStepsize), gNrPasses);
  // std::cout << "minimal baseline = " << *min_element(fbaseline.begin(),fbaseline.end(),find_min_wo_zero) << std::endl;
}
void dirich::DoThreshScanOverBase(uint8_t FirstChannel, uint8_t LastChannel, double FromThrmV, double ToThrmV, double MeasureTime, double StepSize, int NrPasses){
  int ret;
  if(fjans_readout){
    uint32_t scaler_switch[] = {0xffffffff};
    TRBAccessMutex.Lock();
    ret=trb_register_write_mem(gBoardAddress,0xdf80,0,scaler_switch,1); //switching on scaler for this dirich .... only needed when using jan's readout
    TRBAccessMutex.UnLock();
    if(ret<0) std::cerr << "Error switching on scalers" << std::endl;  
  }

    // std::cout << "Set0" << std::endl;
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
    gRateGraphsOverBase[ichannel]->Set(0);
  }
  for(int ipass=0;ipass<NrPasses;++ipass){
    ret=WriteThresholds(OFFTHRESH);
    usleep(THRESHDELAY);
    for (double thresh=FromThrmV; fabs(thresh)<=ToThrmV; thresh+=StepSize){
      // if(int(1.*((thresh-0)+(ToThrmV-0)*(ipass))/((ToThrmV-0)*(NrPasses))*100)%20==0) std::cout << "Thresholdscan over Noiseband of dirich 0x" << std::hex << gBoardAddress << std::dec <<" is @ " << int(1.*((thresh-0)+(ToThrmV-0)*(ipass))/((ToThrmV-0)*(NrPasses))*100) << "%" << std::endl;
      for (int ichannel=FirstChannel+ipass; ichannel<=LastChannel; ichannel+=NrPasses){
        if(fbaseline[ichannel]==0) continue;
        // std::cout << thresh << " " << Thr_DtomV(fnoisewidth[ichannel]) << std::endl;
        if(fabs(thresh)<Thr_DtomV(fnoisewidth[ichannel])/4){
        // if(thresh<-1*fthresholdmV[ichannel]){
          // std::cout << "skipping" << std::endl;
          continue;
        }
        ret=WriteSingleThreshold(ichannel,Thr_mVtoD(thresh)+fbaseline[ichannel]);
      }
      usleep(THRESHDELAY);
      double* rates;
      rates = GetRates(MeasureTime);
      for (int ichannel=FirstChannel+ipass; ichannel<LastChannel; ichannel+=NrPasses){
        if(fbaseline[ichannel]==0) continue;
        if(fabs(thresh)<Thr_DtomV(fnoisewidth[ichannel])/4) continue;
        // if(thresh<-1*fthresholdmV[ichannel]) continue;
        // std::cout << std::dec << "Channel " << ichannel << " thr " << thresh << " rate " << rates[ichannel] << std::endl;
        gRateGraphsOverBase[ichannel]->SetPoint(gRateGraphsOverBase[ichannel]->GetN(),thresh,1.*rates[ichannel]);
      }

      int finish_counter=0;
      for (int ichannel=FirstChannel+ipass; ichannel<LastChannel; ichannel+=NrPasses){
        int number_of_points = gRateGraphsOverBase[ichannel]->GetN();
        if(number_of_points>3){
          if(gRateGraphsOverBase[ichannel]->GetY()[number_of_points-1]< 3. && gRateGraphsOverBase[ichannel]->GetY()[number_of_points-2] < 3. && gRateGraphsOverBase[ichannel]->GetY()[number_of_points-3] < 3.) finish_counter++;
          // std::cout << "gRateGraphsOverBase[ichannel]->GetY()[number_of_points-1] " << gRateGraphsOverBase[ichannel]->GetY()[number_of_points-1] << " gRateGraphsOverBase[ichannel]->GetY()[number_of_points-2] " << gRateGraphsOverBase[ichannel]->GetY()[number_of_points-2] << " gRateGraphsOverBase[ichannel]->GetY()[number_of_points-3] " << gRateGraphsOverBase[ichannel]->GetY()[number_of_points-3] << std::endl;
        }
      }
      std::cout << "thresh " << thresh << " finish_counter " << finish_counter << " (FirstChannel-(FirstChannel+ipass))/NrPasses " << std::dec << (LastChannel-(FirstChannel+ipass))/NrPasses << std::endl;
      if(finish_counter==(LastChannel-(FirstChannel+ipass))/NrPasses) break;

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

void dirich::FindMinThreshScanOverBase(double gThreshold_finding_method){
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++){
    TGraph* temp = new TGraph("temp","temp");
    for (int ipoint=2; ipoint<gDiffRateGraphsOverBase[ichannel]->GetN()-2; ++ipoint){
          temp->SetPoint(temp->GetN(),
                        (gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint-1]+gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint+1]) * 0.5,
                        (-gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint-2]+8*gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint-1]-8*gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint+1]+gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint+2]) /
                        (12*fabs(gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint-1]-gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint+1])));
    }
    double minimum_value=temp->GetY()[0];
    double minimum=temp->GetX()[0];
    int counter=0;
    for (int ipoint=1; ipoint<temp->GetN(); ++ipoint){
      if(temp->GetY()[ipoint-1]<minimum_value){
        minimum_value = temp->GetY()[ipoint];
        minimum = temp->GetX()[ipoint];
      }
      if(temp->GetY()[ipoint-1]>0 && temp->GetY()[ipoint]<0){
        minimum_value = -9999999;
        minimum = temp->GetX()[ipoint-1];
        counter++;
      }
    }
    std::cout << "perfect threshold for dirich 0x" << std::hex << gBoardAddress << "'s channel " << ichannel << " is at " << minimum << "lying in the" << counter << "rd valley." << std::endl;
    SetSingleThresholdmV(ichannel,minimum);
  }

}

void dirich::MakeDiffGraphsOverBase(){
	// MakeDiffGraphsOverBase(gLowerEdge, gUpperEdge);
// }
// void dirich::MakeDiffGraphsOverBase(uint16_t FromThr, uint16_t ToThr){
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
		// if(FromThr==0){
		// 	FromThr=fbaseline[ichannel]+fnoisewidth[ichannel]/2;
		// }  	
    // std::cout << "Set0" << std::endl;
    gDiffRateGraphsOverBase[ichannel]->Set(0);
    // std::cout << "Set" << gDiffRateGraphsOverBase[ichannel]->GetN() << std::endl;
    // doing the derivative (5 point stencil)
    for (int ipoint=2; ipoint<gRateGraphsOverBase[ichannel]->GetN()-2; ++ipoint) {
      // gDiffRateGraphsOverBase[ichannel]->SetPoint(gDiffRateGraphsOverBase[ichannel]->GetN(),
      //                                     (gRateGraphsOverBase[ichannel]->GetX()[ipoint-1]+gRateGraphsOverBase[ichannel]->GetX()[ipoint+1]) * 0.5,
      //                                     (gRateGraphsOverBase[ichannel]->GetY()[ipoint-1]-gRateGraphsOverBase[ichannel]->GetY()[ipoint+1]) /
      //                                     fabs(gRateGraphsOverBase[ichannel]->GetX()[ipoint-1]-gRateGraphsOverBase[ichannel]->GetX()[ipoint+1]));
      gDiffRateGraphsOverBase[ichannel]->SetPoint(gDiffRateGraphsOverBase[ichannel]->GetN(),
                                          (gRateGraphsOverBase[ichannel]->GetX()[ipoint-1]+gRateGraphsOverBase[ichannel]->GetX()[ipoint+1]) * 0.5,
                                          (-gRateGraphsOverBase[ichannel]->GetY()[ipoint-2]+8*gRateGraphsOverBase[ichannel]->GetY()[ipoint-1]-8*gRateGraphsOverBase[ichannel]->GetY()[ipoint+1]+gRateGraphsOverBase[ichannel]->GetY()[ipoint+2]) /
                                          (12*fabs(gRateGraphsOverBase[ichannel]->GetX()[ipoint-1]-gRateGraphsOverBase[ichannel]->GetX()[ipoint+1])));
    }
    //smoothing the graph
    double* x_values = gDiffRateGraphsOverBase[ichannel]->GetX();
    double* y_values = gDiffRateGraphsOverBase[ichannel]->GetY();
    for (int ipoint=1; ipoint<gDiffRateGraphsOverBase[ichannel]->GetN()-1; ipoint+=3) {
      gDiffRateGraphsOverBase[ichannel]->SetPoint(ipoint-1,x_values[ipoint],1.*(y_values[ipoint-1]+y_values[ipoint]+y_values[ipoint+1])/3);
    }
  }
}

void dirich::MakeGraphsOverBase(){
  MakeGraphsOverBase(gUpperEdge);
}
void dirich::MakeGraphsOverBase(uint16_t ToThr){
  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
    // std::cout << "Set0" << std::endl;
    gRateGraphsOverBase[ichannel]->Set(0);
    // std::cout << "Set" << gRateGraphsOverBase[ichannel]->GetN() << std::endl;
    if(fbaseline[ichannel]==0) continue;
    for(int ipoint=0; ipoint<gRateGraphs[ichannel]->GetN(); ++ipoint){
      if(gRateGraphs[ichannel]->GetX()[ipoint]<fbaseline[ichannel]+1.*fnoisewidth[ichannel]/2) continue;
      // std::cout << ipoint << " " << gRateGraphsOverBase[ichannel]->GetN() << " " << Thr_DtomV((gRateGraphs[ichannel]->GetX()[ipoint])-fbaseline[ichannel]) << " " << gRateGraphs[ichannel]->GetY()[ipoint] << std::endl;
      gRateGraphsOverBase[ichannel]->SetPoint(gRateGraphsOverBase[ichannel]->GetN(),Thr_DtomV((gRateGraphs[ichannel]->GetX()[ipoint])-fbaseline[ichannel]),gRateGraphs[ichannel]->GetY()[ipoint]);
    }
    // std::cout << gBoardAddress << " " << ichannel << " " << gRateGraphsOverBase[ichannel]->GetN() << " " << gRateGraphs[ichannel]->GetN() << " " << fbaseline[ichannel] << " " << gRateGraphs[ichannel]->GetX()[gRateGraphs[ichannel]->GetN()-1] << std::endl;
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

      if(noiseedgeleft!=-1 && noiseedgeright!=-1) break;
    }
    if(noiseedgeleft==-1 || noiseedgeright==-1){
      fbaseline_old[ichannel] = fbaseline[ichannel];
      fbaseline[ichannel] = 0;
      fnoisewidth_old[ichannel] = fnoisewidth[ichannel];
      fnoisewidth[ichannel] = 0;
      std::cout << "No baseline found for channel " << ichannel << "on dirich 0x" << std::hex << gBoardAddress << std::dec << std::endl;
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
