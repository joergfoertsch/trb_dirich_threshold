#include "trbnet.h"
// #include "dirich_sim.C"

#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TFile.h"
#include "TText.h"
#include "TMath.h"
#include "TThread.h"
#include <map>
#include <array>
#include <sstream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/range.hpp>
// #include <iomanip>

#include "dirich_v11.C"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

// #define 0 0
// #define LASTCHANNEL 31

#ifndef NCH
  const int NRCHANNELS = 32;       //Nr of TDC channels in dirich
  // const int NRCHANNELS = 4;       //Nr of TDC channels in dirich
  #define NCH
#endif

std::map<uint16_t,std::shared_ptr<dirich>> dirichlist ={};

std::map<uint16_t,TCanvas*> canvaslist;

std::vector<TCanvas*> canvasvector;

TH2* get_2D_rate_histo(std::shared_ptr<dirich>  dirichptr){
	TH2D* histo;
	gStyle->SetOptStat(0);
  if(dirichptr==NULL){
    histo = new TH2D("2D Rate vs. Threshold of all diriches","2D Rate vs. Threshold of all diriches",dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5,(dirichlist.begin()->second->gUpperEdge-dirichlist.begin()->second->gLowerEdge)/dirichlist.begin()->second->gStepsize,dirichlist.begin()->second->gLowerEdge,dirichlist.begin()->second->gUpperEdge);
		int idirich=0;
		// std::map<uint16_t,dirich*>::iterator dirichlistiterator = dirichlist.begin();
		TLine* dirich_line_left = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()));
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
    for (auto& dirichitem : dirichlist){
			TLine* dirich_line_right = new TLine(histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()));
			dirich_line_right->SetLineWidth(2);
      dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
    	TText* dirich_name = new TText(histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1)-(0.1*(histo->GetYaxis()->GetBinLowEdge(1)-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()))),Form("0x%x",dirichitem.first));
    	// TText* dirich_name = new TText(histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1)+500,Form("0x%x",dirichitem.first));
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				// std::cout << (int)dirichlist.size() << std::endl;
      	// std::cout << "dirich: " << std::hex << dirichlistiterator.first << std::dec << "GetN " << dirichlistiterator->second->gRateGraphs[ichannel]->GetN() << std::endl;
      	for(int ipoint=0;ipoint<dirichitem.second->gRateGraphs[ichannel]->GetN();++ipoint){
      		// std::cout << idirich << " " << std::hex << dirichlistiterator.first << std::dec << " " <<  ichannel << " " << ipoint << " " << dirichlistiterator->second->gRateGraphs[ichannel]->GetX()[ipoint] << " " << dirichlistiterator->second->gRateGraphs[ichannel]->GetY()[ipoint] << std::endl;
        	histo->Fill(idirich*NRCHANNELS+ichannel,dirichitem.second->gRateGraphs[ichannel]->GetX()[ipoint],dirichitem.second->gRateGraphs[ichannel]->GetY()[ipoint]);
        	// std::cout << ichannel << " " << int(dirichlist.size()*NRCHANNELS/20+1) << std::endl;
					if(ichannel%8==0) histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel));
					// if(ichannel%int(dirichlist.size()*NRCHANNELS/20+1)==0) histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel));
					else histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,"");
      	}
				TLine* baseline_line = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),dirichitem.second->GetSingleBaseline(ichannel),histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),dirichitem.second->GetSingleBaseline(ichannel));
				baseline_line->SetLineColor(kRed);
				baseline_line->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(baseline_line);
        TLine* baseline_line_old = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),dirichitem.second->GetSingleBaseline_old(ichannel),histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),dirichitem.second->GetSingleBaseline_old(ichannel));
        baseline_line_old->SetLineColor(kBlack);
        baseline_line_old->SetLineWidth(2);
        histo->GetListOfFunctions()->Add(baseline_line_old);        
      }
    	// ++dirichlistiterator;
    	++idirich;
    }
  }
  else{
	  histo = new TH2D(Form("2D Rate vs. Threshold of %x",dirichptr->GetBoardAddress()),Form("2D Rate vs. Threshold of %x",dirichptr->GetBoardAddress()),NRCHANNELS,-.5,NRCHANNELS-.5,(dirichptr->gUpperEdge-dirichptr->gLowerEdge)/dirichptr->gStepsize,dirichptr->gLowerEdge,dirichptr->gUpperEdge);
	  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(int ipoint=0;ipoint<dirichptr->gRateGraphs[ichannel]->GetN();++ipoint){
				histo->Fill(ichannel,dirichptr->gRateGraphs[ichannel]->GetX()[ipoint],dirichptr->gRateGraphs[ichannel]->GetY()[ipoint]);
			}
			TLine* baseline_line = new TLine(histo->GetXaxis()->GetBinLowEdge(ichannel+1),dirichptr->GetSingleBaseline(ichannel),histo->GetXaxis()->GetBinUpEdge(ichannel+1),dirichptr->GetSingleBaseline(ichannel));
			baseline_line->SetLineColor(kRed);
			baseline_line->SetLineWidth(2);
			histo->GetListOfFunctions()->Add(baseline_line);
      TLine* baseline_line_old = new TLine(histo->GetXaxis()->GetBinLowEdge(ichannel+1),dirichptr->GetSingleBaseline_old(ichannel),histo->GetXaxis()->GetBinUpEdge(ichannel+1),dirichptr->GetSingleBaseline_old(ichannel));
      baseline_line_old->SetLineColor(kBlack);
      baseline_line_old->SetLineWidth(2);
      histo->GetListOfFunctions()->Add(baseline_line_old);      
	  }
	}
	// std::cout << "finished histo" << std::endl;
  
	histo->SetMinimum(0);
	histo->GetXaxis()->SetTitle("Channel Nr");
	// histo->GetXaxis()->SetTitleOffset();
	histo->GetYaxis()->SetTitle("Threshold");
	histo->GetZaxis()->SetTitle("Rate");
  return histo;
}

TMultiGraph* get_2D_mgr_diff_over_thr_histo(std::shared_ptr<dirich>  dirichptr){
	TMultiGraph* multig = new TMultiGraph();
	if(dirichptr==NULL){
    // multig->SetTitle("Differentiated rate graph over baseline of all dirich;Channel Nr;Threshold;Differentiated rate");
		multig->SetTitle("Differentiated rate graph over baseline of all dirich;Threshold;Differentiated rate");
		multig->SetName("Differentiated rate graph over baseline of all dirich (Mutligraph)");
	  for (auto& dirichitem : dirichlist){
	  	for(auto& gDiffRateGraphsOverBaseIT : dirichitem.second->gDiffRateGraphsOverBase){
				multig->Add(gDiffRateGraphsOverBaseIT,"PL");
			}
		}
	}
	else{
    // multig->SetTitle(Form("Differentiated rate graph over baseline of dirich 0x%x;Channel Nr;Threshold;Differentiated rate",dirichptr->GetBoardAddress()));
		multig->SetTitle(Form("Differentiated rate graph over baseline of dirich 0x%x;Threshold;Differentiated rate",dirichptr->GetBoardAddress()));
		multig->SetName(Form("Differentiated rate graph over baseline of dirich 0x%x (Multigraph)",dirichptr->GetBoardAddress()));
  	for(auto& gDiffRateGraphsOverBaseIT : dirichptr->gDiffRateGraphsOverBase){
			multig->Add(gDiffRateGraphsOverBaseIT,"PL");
		}
	}
  // multig->SetMinimum(0);
  // multig->GetHistogram()->GetYaxis()->SetRangeUser(0,100);
	return multig;
}

TGraph2D* get_2D_gr_diff_over_thr_histo(std::shared_ptr<dirich>  dirichptr){
	TGraph2D* g2d = new TGraph2D();
	if(dirichptr==NULL){
		g2d->SetTitle("Differentiated rate graph over baseline of all dirich;Channel Nr;Threshold;Differentiated rate");
		g2d->SetName("Differentiated rate graph over baseline of all dirich (2D_Graph)");
		int idirich=0;
	  for (auto& dirichitem : dirichlist){
	  	int ichannel=0;
	  	for(auto& gDiffRateGraphsOverBaseIT : dirichitem.second->gDiffRateGraphsOverBase){
	  		// std::cout << " " << idirich << " " << dirichitem.first << " " << ichannel << " " << gDiffRateGraphsOverBaseIT->GetN() << std::endl;
	    	for(int ipoint=0;ipoint<gDiffRateGraphsOverBaseIT->GetN();++ipoint){
	    		// std::cout << " " << idirich << " " << dirichitem.first << " " << ichannel << " " << ipoint << " " << gDiffRateGraphsOverBaseIT->GetX()[ipoint] << " " << gDiffRateGraphsOverBaseIT->GetY()[ipoint] << std::endl;
					g2d->SetPoint(g2d->GetN(),idirich*NRCHANNELS+ichannel,gDiffRateGraphsOverBaseIT->GetX()[ipoint],gDiffRateGraphsOverBaseIT->GetY()[ipoint]);
				}
				ichannel++;
			}
			idirich++;
		}
	}
	else{
		g2d->SetTitle(Form("Differentiated rate graph over baseline of dirich 0x%x;Channel Nr;Threshold;Differentiated rate",dirichptr->GetBoardAddress()));
		g2d->SetName(Form("Differentiated rate graph over baseline of dirich 0x%x (2D_Graph)",dirichptr->GetBoardAddress()));
  	int ichannel=0;
  	for(auto& gDiffRateGraphsOverBaseIT : dirichptr->gDiffRateGraphsOverBase){
  		// std::cout << " " << idirich << " " << dirichitem.first << " " << ichannel << " " << gDiffRateGraphsOverBaseIT->GetN() << std::endl;
    	for(int ipoint=0;ipoint<gDiffRateGraphsOverBaseIT->GetN();++ipoint){
    		// std::cout << " " << idirich << " " << dirichitem.first << " " << ichannel << " " << ipoint << " " << gDiffRateGraphsOverBaseIT->GetX()[ipoint] << " " << gDiffRateGraphsOverBaseIT->GetY()[ipoint] << std::endl;
				g2d->SetPoint(g2d->GetN(),ichannel,gDiffRateGraphsOverBaseIT->GetX()[ipoint],gDiffRateGraphsOverBaseIT->GetY()[ipoint]);
			}
			ichannel++;
		}
	}
  g2d->SetMinimum(0);
  g2d->GetZaxis()->SetRangeUser(0,100);
	return g2d;
}

TH2* get_2D_diff_over_thr_histo(std::shared_ptr<dirich>  dirichptr){
	TH2D* histo;
	TH2D* divided_histo;
  // divided_histo->SetDirectory(0);
	gStyle->SetOptStat(0);
  if(dirichptr==NULL){
    double max_value=-9999;
  	double min_value=9999;
  	double min_width=1000;
    for (auto& dirichitem : dirichlist){
	  	for(auto& gDiffRateGraphsOverBaseIT : dirichitem.second->gDiffRateGraphsOverBase){
        if(gDiffRateGraphsOverBaseIT->GetN()!=0 && max_value<gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1]) max_value = gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1];
	  		if(gDiffRateGraphsOverBaseIT->GetN()!=0 && max_value<gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1]) max_value = gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1];
	  		if(gDiffRateGraphsOverBaseIT->GetN()>=2 && min_width>abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1])) min_width = abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1]);
	  	}
	  }
    histo = new TH2D("2D Differentiated Rate vs. Threshold over baseline of all diriches","2D Differentiated Rate vs. Threshold over baseline of all diriches",dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5,max_value/min_width/2,0,max_value);
    divided_histo = new TH2D("temp_diff","temp_diff",dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5,max_value/min_width/2,0,max_value);
		int idirich=0;
		TLine* dirich_line_left = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()));
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
    for (auto& dirichitem : dirichlist){
			TLine* dirich_line_right = new TLine(histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()));
			dirich_line_right->SetLineWidth(2);    	
      dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
    	TText* dirich_name = new TText(histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1)-(0.1*(histo->GetYaxis()->GetBinLowEdge(1)-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()))),Form("0x%x",dirichitem.first));
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				// std::cout << (int)dirichlist.size() << std::endl;
      	// std::cout << "dirich: " << std::hex << dirichitem.first << std::dec << "GetN " << dirichitem.second->gRateGraphs[ichannel]->GetN() << std::endl;
      	for(int ipoint=0;ipoint<dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetN();++ipoint){
      		// std::cout << idirich << " dirich " << std::hex << dirichitem.first << std::dec << " C " <<  ichannel << " P " << ipoint << " X " << dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint] << " Y " << dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint] << std::endl;
        	histo->Fill(idirich*NRCHANNELS+ichannel,dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint],dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint]);
        	divided_histo->Fill(idirich*NRCHANNELS+ichannel,dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint]);
        	// std::cout << ichannel << " " << int(dirichlist.size()*NRCHANNELS/20+1) << std::endl;
					if(ichannel%8==0) histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel));
					else histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,"");
      	}
				TLine* thr_line = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),-1*dirichitem.second->GetSingleThresholdmV(ichannel),histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),-1*dirichitem.second->GetSingleThresholdmV(ichannel));
				thr_line->SetLineColor(kRed);
				thr_line->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(thr_line);
      }
    	++idirich;
    }
    idirich=0;
    // for (auto& dirichitem : dirichlist){
    //   for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
    //   	histo->SetBinContent(idirich*NRCHANNELS+ichannel+1,int(-1*dirichitem.second->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,-1);
    //   	divided_histo->SetBinContent(idirich*NRCHANNELS+ichannel+1,int(-1*dirichitem.second->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,1);
    //   }
    //   idirich++;
    // }
  }
  else{
  	double max_value=0;
  	double min_width=10000;
		for(auto& gDiffRateGraphsOverBaseIT : dirichptr->gDiffRateGraphsOverBase){
			if(gDiffRateGraphsOverBaseIT->GetN()!=0 && max_value<gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1]) max_value = gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1];
			if(gDiffRateGraphsOverBaseIT->GetN()>=2 && min_width>abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1])) min_width = abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1]);
	  }
    histo = new TH2D(Form("2D Differentiated Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),Form("2D Differentiated Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),NRCHANNELS,-.5,NRCHANNELS-.5,max_value/min_width/2,0,max_value);
    divided_histo = new TH2D("temp_diff","temp_diff",NRCHANNELS,-.5,NRCHANNELS-.5,max_value/min_width/2,0,max_value);
	  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(int ipoint=0;ipoint<dirichptr->gDiffRateGraphsOverBase[ichannel]->GetN();++ipoint){
      	histo->Fill(ichannel,dirichptr->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint],dirichptr->gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint]);
      	divided_histo->Fill(ichannel,dirichptr->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint]);
			}
      TLine* thr_line = new TLine(histo->GetXaxis()->GetBinLowEdge(ichannel+1),-1*dirichptr->GetSingleThresholdmV(ichannel),histo->GetXaxis()->GetBinUpEdge(ichannel+1),-1*dirichptr->GetSingleThresholdmV(ichannel));
      thr_line->SetLineColor(kRed);
      thr_line->SetLineWidth(2);
      histo->GetListOfFunctions()->Add(thr_line);
	  }
    // for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
    // 	histo->SetBinContent(ichannel+1,int(-1*dirichptr->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,-1);
    // 	divided_histo->SetBinContent(ichannel+1,int(-1*dirichptr->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,1);
    // }
	}
	// std::cout << "finished histo" << std::endl;

  histo->Divide(divided_histo);
  
  histo->Divide(divided_histo);
  for(int ibin=1;ibin<(histo->GetNbinsX()+2)*(histo->GetNbinsY()+2);++ibin){
    // if(histo->GetBinContent(ibin)<-1.)histo->SetBinContent(ibin,0);
    histo->SetBinError(ibin,0);
  }
  histo->SetMinimum(0.);
	// histo->GetZaxis()->SetRangeUser(0.,30.);
	histo->GetXaxis()->SetTitle("Channel Nr");
	// histo->GetXaxis()->SetTitleOffset();
	histo->GetYaxis()->SetTitle("Threshold");
	histo->GetZaxis()->SetTitle("Differentiated rate");
  return histo;
}

TH2* get_2D_rate_over_thr_histo(std::shared_ptr<dirich>  dirichptr){
	TH2D* histo;
	TH2D* divided_histo;
  // divided_histo->SetDirectory(0);
	gStyle->SetOptStat(0);
  if(dirichptr==NULL){
    double max_value=-9999;
  	double min_value=9999;
  	double min_width=1000;
    for (auto& dirichitem : dirichlist){
	  	for(auto& gRateGraphsOverBaseIT : dirichitem.second->gRateGraphsOverBase){
        if(gRateGraphsOverBaseIT->GetN()!=0 && max_value<gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1]) max_value = gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1];
	  		if(gRateGraphsOverBaseIT->GetN()!=0 && max_value<gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1]) min_value = gRateGraphsOverBaseIT->GetX()[0];
	  		if(gRateGraphsOverBaseIT->GetN()>=2 && min_width>abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1])) min_width = abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1]);
	  	}
	  }
    histo = new TH2D("2D Rate vs. Threshold over baseline of all diriches","2D Rate vs. Threshold over baseline of all diriches",dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5,max_value/min_width/2,0,max_value);
    divided_histo = new TH2D("temp_diff","temp_diff",dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5,max_value/min_width/2,0,max_value);
		int idirich=0;
		TLine* dirich_line_left = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()));
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
    for (auto& dirichitem : dirichlist){
			TLine* dirich_line_right = new TLine(histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()));
			dirich_line_right->SetLineWidth(2);    	
      dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
    	TText* dirich_name = new TText(histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1)-(0.1*(histo->GetYaxis()->GetBinLowEdge(1)-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY()))),Form("0x%x",dirichitem.first));
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				// std::cout << (int)dirichlist.size() << std::endl;
      	// std::cout << "dirich: " << std::hex << dirichitem.first << std::dec << "GetN " << dirichitem.second->gRateGraphs[ichannel]->GetN() << std::endl;
      	for(int ipoint=0;ipoint<dirichitem.second->gRateGraphsOverBase[ichannel]->GetN();++ipoint){
      		// std::cout << idirich << " dirich " << std::hex << dirichitem.first << std::dec << " C " <<  ichannel << " P " << ipoint << " X " << dirichitem.second->gRateGraphsOverBase[ichannel]->GetX()[ipoint] << " Y " << dirichitem.second->gRateGraphsOverBase[ichannel]->GetY()[ipoint] << std::endl;
        	histo->Fill(idirich*NRCHANNELS+ichannel,dirichitem.second->gRateGraphsOverBase[ichannel]->GetX()[ipoint],dirichitem.second->gRateGraphsOverBase[ichannel]->GetY()[ipoint]);
        	divided_histo->Fill(idirich*NRCHANNELS+ichannel,dirichitem.second->gRateGraphsOverBase[ichannel]->GetX()[ipoint]);
        	// std::cout << ichannel << " " << int(dirichlist.size()*NRCHANNELS/20+1) << std::endl;
					if(ichannel%8==0) histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel));
					else histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,"");
      	}
				TLine* thr_line = new TLine(histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),-1*dirichitem.second->GetSingleThresholdmV(ichannel),histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),-1*dirichitem.second->GetSingleThresholdmV(ichannel));
				thr_line->SetLineColor(kRed);
				thr_line->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(thr_line);
      }
    	++idirich;
    }
    idirich=0;
    // for (auto& dirichitem : dirichlist){
    //   for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
    //   	histo->SetBinContent(idirich*NRCHANNELS+ichannel+1,int(-1*dirichitem.second->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,-1);
    //   	divided_histo->SetBinContent(idirich*NRCHANNELS+ichannel+1,int(-1*dirichitem.second->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,1);
    //   }
    //   idirich++;
    // }
  }
  else{
  	double max_value=0;
  	double min_width=1000;
		for(auto& gRateGraphsOverBaseIT : dirichptr->gRateGraphsOverBase){
			if(gRateGraphsOverBaseIT->GetN()!=0 && max_value<gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1]) max_value = gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1];
			if(gRateGraphsOverBaseIT->GetN()>=2 && min_width>abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1])) min_width = abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1]);
	  }
    histo = new TH2D(Form("2D Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),Form("2D Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),NRCHANNELS,-.5,NRCHANNELS-.5,max_value/min_width/2,0,max_value);
    divided_histo = new TH2D("temp_diff","temp_diff",NRCHANNELS,-.5,NRCHANNELS-.5,max_value/min_width/2,0,max_value);
	  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(int ipoint=0;ipoint<dirichptr->gRateGraphsOverBase[ichannel]->GetN();++ipoint){
      	histo->Fill(ichannel,dirichptr->gRateGraphsOverBase[ichannel]->GetX()[ipoint],dirichptr->gRateGraphsOverBase[ichannel]->GetY()[ipoint]);
      	divided_histo->Fill(ichannel,dirichptr->gRateGraphsOverBase[ichannel]->GetX()[ipoint]);
			}
      TLine* thr_line = new TLine(histo->GetXaxis()->GetBinLowEdge(ichannel+1),-1*dirichptr->GetSingleThresholdmV(ichannel),histo->GetXaxis()->GetBinUpEdge(ichannel+1),-1*dirichptr->GetSingleThresholdmV(ichannel));
      thr_line->SetLineColor(kRed);
      thr_line->SetLineWidth(2);
      histo->GetListOfFunctions()->Add(thr_line);
	  }
    // for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
    // 	histo->SetBinContent(ichannel+1,int(-1*dirichptr->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,-1);
    // 	divided_histo->SetBinContent(ichannel+1,int(-1*dirichptr->GetSingleThresholdmV(ichannel)/histo->GetYaxis()->GetBinWidth(1))+1,1);
    // }
	}
	// std::cout << "finished histo" << std::endl;

  histo->Divide(divided_histo);

	histo->SetMinimum(0);
	histo->GetXaxis()->SetTitle("Channel Nr");
	// histo->GetXaxis()->SetTitleOffset();
	histo->GetYaxis()->SetTitle("Threshold");
	histo->GetZaxis()->SetTitle("Rate");
  return histo;
}

TH1* get_noisewidth_histo(std::shared_ptr<dirich>  dirichptr){
  TH1* histo;
  if(dirichptr==NULL){
    histo = new TH1D("Noisewidthhistogram of all diriches","Noisewidthhistogram of all diriches",dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5);
		int idirich=0;
    for (auto& dirichitem : dirichlist){
    	TText* dirich_name = new TText(histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1)-(0.1*histo->GetYaxis()->GetBinLowEdge(1)-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())),Form("0x%x",dirichitem.first));
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
        histo->SetBinContent(idirich*NRCHANNELS+ichannel+1,dirich::Thr_DtomV(dirichitem.second->GetSingleNoisewidth(ichannel)));
				if(ichannel%8==0) histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel));
				else histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,"");
      }
  	++idirich;
    }
		TLine* dirich_line_left = new TLine(histo->GetXaxis()->GetBinLowEdge(1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge(1),histo->GetMaximum()*1.05);
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
    for(int i=0;i<idirich;++i){
			TLine* dirich_line_right = new TLine(histo->GetXaxis()->GetBinLowEdge((i+1)*NRCHANNELS+1),histo->GetYaxis()->GetBinLowEdge(1),histo->GetXaxis()->GetBinLowEdge((i+1)*NRCHANNELS+1),histo->GetMaximum()*1.05);
			dirich_line_right->SetLineWidth(2);    	
      dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);    	
    }
  }
  else{
	  histo = new TH1D(Form("Noisewidthhistogram of %x",dirichptr->GetBoardAddress()),Form("Noisewidthhistogram of %x",dirichptr->GetBoardAddress()),NRCHANNELS,-.5,NRCHANNELS-.5);
	  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
	    histo->SetBinContent(ichannel+1,(dirich::Thr_DtomV(dirichptr->GetSingleNoisewidth(ichannel))));
	  }
	}

  histo->GetYaxis()->SetTitle("NoisewidthinmV");
  histo->SetMinimum(0);
  histo->GetXaxis()->SetTitle("Channel Nr");
  return histo;
}

TH1* get_diff_histo(std::shared_ptr<dirich>  dirichptr, bool baseline1_noisewidth0){
	TH1* histo;
  if(dirichptr==NULL){
  	if(baseline1_noisewidth0==1) histo = new TH1D("Difference in baseline of all diriches","Difference in baseline of all diriches",dirichlist.size()*200,-300,+300);
  	else histo = new TH1D("Difference in noisewidth of all diriches","Difference in noisewidth of all diriches",dirichlist.size()*200,-300,+300);
		int idirich=0;
    for (auto& dirichitem : dirichlist){
      for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
        if(baseline1_noisewidth0==1) histo->Fill(dirichitem.second->GetSingleBaseline(ichannel)-dirichitem.second->GetSingleBaseline_old(ichannel));
        else histo->Fill(dirichitem.second->GetSingleNoisewidth(ichannel)-dirichitem.second->GetSingleNoisewidth_old(ichannel));
      }
    }
  	++idirich;
  }
  else{
  	if(baseline1_noisewidth0==1) histo = new TH1D(Form("Difference in baseline of %x",dirichptr->GetBoardAddress()),Form("Difference in baseline of %x",dirichptr->GetBoardAddress()),200,-300,+300);
  	else histo = new TH1D(Form("Difference in noisewidth of %x",dirichptr->GetBoardAddress()),Form("Difference in noisewidth of %x",dirichptr->GetBoardAddress()),200,-300,+300);
	  for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
      if(baseline1_noisewidth0==1) histo->Fill(dirichptr->GetSingleBaseline(ichannel)-dirichptr->GetSingleBaseline_old(ichannel));
      else histo->Fill(dirichptr->GetSingleNoisewidth(ichannel)-dirichptr->GetSingleNoisewidth_old(ichannel));
	  }
	}

  histo->GetYaxis()->SetTitle("Number of");
  if(baseline1_noisewidth0==1) histo->GetXaxis()->SetTitle("Difference between old and new baseline");
  else histo->GetXaxis()->SetTitle("Difference between old and new noisewidth");
  return histo;
}

void clear_canvas_vector(){
	while(canvasvector.size()!=0){
		if(canvasvector.back()==NULL){
			std::cout <<  "1" << std::endl;
			canvasvector.pop_back();
		}
		else{
			std::cout <<  "2" << std::endl;
			canvasvector.back()->Clear();
			std::cout <<  "3" << std::endl;
			canvasvector.back()->Close();
			std::cout <<  "4" << std::endl;
			canvasvector.back()->Closed();
			std::cout <<  "5" << std::endl;
			delete canvasvector.back();
			std::cout <<  "6" << std::endl;
			canvasvector.back()=NULL;
			std::cout <<  "7" << std::endl;
			canvasvector.pop_back();
			std::cout <<  "8" << std::endl;
		}
	}
}

void draw_multigraph2D(TMultiGraph* multigraph,TCanvas* canvas){
	if(canvas==0){
		canvasvector.emplace_back(new TCanvas(Form("Canvas%i",(int)canvasvector.size()),Form("Canvas%i",(int)canvasvector.size()),1920,1080));
		canvasvector.back()->cd(0);
	}
	else{
	  canvas->cd(0);
	}
  multigraph->Draw("a fb l3d");
  gPad->SetTheta(0);
  gPad->SetPhi(-90);
  gPad->Update();
  if(canvas==0){
    canvasvector.back()->Modified();
    canvasvector.back()->Update();
  }
  else{
    canvas->Modified();
    canvas->Update();
  }
}

void draw_multigraph(TMultiGraph* multigraph,TCanvas* canvas){
  if(canvas==0){
    canvasvector.emplace_back(new TCanvas(Form("Canvas%i",(int)canvasvector.size()),Form("Canvas%i",(int)canvasvector.size()),1920,1080));
    canvasvector.back()->cd(0);
  }
  else{
    canvas->cd(0);
  }
  multigraph->Draw("alp");
  double max_value = 0;
  for(auto&& graph : (*multigraph->GetListOfGraphs())){
    for(int i=10 ; i < ((TGraph*)graph)->GetN() ; ++i){
      max_value = ((TGraph*)graph)->GetY()[i] > max_value ? ((TGraph*)graph)->GetY()[i] : max_value;
    }
  }
  // std::cout << "max_value" << max_value << std::endl;
  multigraph->GetHistogram()->GetYaxis()->SetRangeUser(0,max_value==0 ? 100 : max_value*1.05);
  if(canvas==0){
    canvasvector.back()->Modified();
    canvasvector.back()->Update();
  }
  else{
    canvas->Modified();
    canvas->Update();
  }
}

void draw_graph2D(TGraph2D* graph2d,TCanvas* canvas){
	if(canvas==0){
		canvasvector.emplace_back(new TCanvas(Form("Canvas%i",(int)canvasvector.size()),Form("Canvas%i",(int)canvasvector.size()),1920,1080));
		canvasvector.back()->cd(0);
	}
	else{
	  canvas->cd(0);
	}
  graph2d->Draw("surf1");
  gPad->SetTheta(0);
  gPad->SetPhi(-90);
  gPad->Update();  
  if(canvas==0){
    canvasvector.back()->Modified();
    canvasvector.back()->Update();
  }
  else{
    canvas->Modified();
    canvas->Update();
  }
}

void draw_histo(TH1* histo,TCanvas* canvas){
	if(canvas==0){
		canvasvector.emplace_back(new TCanvas(Form("Canvas%i",(int)canvasvector.size()),Form("Canvas%i",(int)canvasvector.size()),1920,1080));
		canvasvector.back()->cd(0);
	}
	else{
	  canvas->cd(0);
	}
  histo->Draw();
  if(canvas==0){
    canvasvector.back()->Modified();
    canvasvector.back()->Update();
  }
  else{
    canvas->Modified();
    canvas->Update();
  }
}

void draw_histo(TH2* histo,TCanvas* canvas){
	if(canvas==0){
		canvasvector.emplace_back(new TCanvas(Form("Canvas%i",(int)canvasvector.size()),Form("Canvas%i",(int)canvasvector.size()),1920,1080));
		canvasvector.back()->cd(0);
	}
	else{
	  canvas->cd(0);
	}
  histo->Draw("COLZ");
  if(canvas==0){
    canvasvector.back()->Modified();
    canvasvector.back()->Update();
  }
  else{
    canvas->Modified();
    canvas->Update();
  }
}

void set_thresholds(std::shared_ptr<dirich> dirichptr, double thrinmV=30.)
{
  if(dirichptr==0){
  	// std::cout << "setting threshold for all diriches: " << std::endl;
    for (auto& dirichlistitem: dirichlist){
      set_thresholds(dirichlistitem.second, thrinmV);
    }
  }
  else if(dirichlist.find(dirichptr->GetBoardAddress())!=dirichlist.end()){
    // if(thrinmV>0.){
      // std::cerr << "positive thresholds are not \"allowed\"!\ninverting value" << std::endl;
      // thrinmV = -1*thrinmV;
    // }
  	// std::cout << "setting threshold for dirich: " << std::hex << dirichptr->GetBoardAddress() << std::endl;
    for (int ichannel=0; ichannel<NRCHANNELS; ++ichannel) {
      dirichptr->SetSingleThresholdmV(ichannel ,thrinmV); 
    }
  // draw_graphs(dirichptr, canvaslist.at(dirichptr->GetBoardAddress()), makechannelvector(0,31), 1, 1, 1, 1);
  }
}

void set_thresholds(std::shared_ptr<dirich> dirichptr, double* thrinmV)
{
  if(dirichptr==0){
  	std::cout << "setting threshold for all diriches: " << std::endl;
    for (auto& dirichlistitem: dirichlist){
      set_thresholds(dirichlistitem.second, thrinmV);
    }
  }
  else if(dirichlist.find(dirichptr->GetBoardAddress())!=dirichlist.end()){
    if (!std::all_of(thrinmV, thrinmV+32, [](double i){ return i > 0.; })) {
      std::cerr << "positive thresholds are not \"allowed\"!\ncanceling set_thresholds" << std::endl;
      return;
    }
  	std::cout << "setting threshold for dirich: " << std::hex << dirichptr->GetBoardAddress() << std::endl;
    for (int ichannel=0; ichannel<NRCHANNELS; ++ichannel) {
  		dirichptr->SetSingleThresholdmV(ichannel ,thrinmV[ichannel]); 
  	}
    // draw_graphs(dirichptr, canvaslist.at(dirichptr->GetBoardAddress()), makechannelvector(0,31), 1, 1, 1, 1);
  }
}

void save_base(std::shared_ptr<dirich>  dirichptr, std::string filename, bool append){
	std::ofstream file;
  if(append) file.open(filename+".thr", std::ios_base::app);
	else file.open(filename+".thr");

	if(!file) std::cerr << "File for saving (" << filename+".thr" << ") could not be opened!" << std::endl;

	if(dirichptr==NULL){
  	for (auto& dirichlistitem: dirichlist){
			file << "# Scan-Settings for 0x" << std::hex << dirichlistitem.first << std::dec << "\n# gMeasureTime\tgLowerEdge\tgUpperEdge\tgStepsize\tgNrPasses\tgMeasureTime_over\tgUpperEdge_over\tgStepsize_over\tgNrPasses_over" << std::endl;
			file << "# " << dirichlistitem.second->gMeasureTime << "\t" << dirichlistitem.second->gLowerEdge << "\t" << dirichlistitem.second->gUpperEdge << "\t" << dirichlistitem.second->gStepsize << "\t" << dirichlistitem.second->gNrPasses << "\t" << dirichlistitem.second->gMeasureTime_over << "\t" << dirichlistitem.second->gUpperEdge_over << "\t" << dirichlistitem.second->gStepsize_over << "\t" << dirichlistitem.second->gNrPasses_over << std::endl;
			file << "# Scan-Data\n# dirich\tchannel\tbaseline\twidth in mV\tthreshold in mV over baseline" << std::endl;
		  for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
  			file  
  			<< std::hex << dirichlistitem.first << std::dec << "\t" 
  			<< ichannel << "\t" 
  			<< dirichlistitem.second->GetSingleBaseline(ichannel) << "\t" 
  			<< dirich::Thr_DtomV(dirichlistitem.second->GetSingleNoisewidth(ichannel)) << "\t" 
  			<< dirichlistitem.second->GetSingleThresholdmV(ichannel) 
  			<< std::endl;
  		}
  	}
	}
	else{
		file << "# Scan-Settings for 0x" << std::hex << dirichptr->GetBoardAddress() << std::dec << "\n# gMeasureTime\tgLowerEdge\tgUpperEdge\tgStepsize\tgNrPasses\tgMeasureTime_over\tgUpperEdge_over\tgStepsize_over\tgNrPasses_over" << std::endl;
		file << "# " << dirichptr->gMeasureTime << "\t" << dirichptr->gLowerEdge << "\t" << dirichptr->gUpperEdge << "\t" << dirichptr->gStepsize << "\t" << dirichptr->gNrPasses << "\t" << dirichptr->gMeasureTime_over << "\t" << dirichptr->gUpperEdge_over << "\t" << dirichptr->gStepsize_over << "\t" << dirichptr->gNrPasses_over << std::endl;
		file << "# Scan-Data\n# dirich\tchannel\tbaseline\twidth in mV\tthreshold in mV over baseline" << std::endl;
		for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
			file  
			<< std::hex << dirichptr->GetBoardAddress() << std::dec << "\t" 
			<< ichannel << "\t" 
			<< dirichptr->GetSingleBaseline(ichannel) << "\t" 
			<< dirich::Thr_DtomV(dirichptr->GetSingleNoisewidth(ichannel)) << "\t" 
			<< dirichptr->GetSingleThresholdmV(ichannel) 
			<< std::endl;
		}
	}
}

void load_base(std::shared_ptr<dirich>  dirichptr, std::string filename, bool uselast, bool set_base, bool set_thr){
	std::string dirichaddress_string="";
	uint16_t dirichaddress=0;
  int channel=0;
  int baseline=0;
  double width=0;
  double thresholdinmV=0;
	for (auto& dirichlistitem: dirichlist){
		if(dirichlistitem.second==NULL){
			std::cerr << "dirich 0x" << std::hex << dirichlistitem.first << std::dec << " was found uninitialized\nRun initialize_diriches(1/0) first!" << std::endl;
			return;
		}
	}
	std::ifstream file;
	file.open(filename);
	if(!file) std::cerr << "File for loading (" << filename << ") could not be opened!" << std::endl;
	if(dirichptr==NULL){
		while(!file.eof()){
			std::string line;
			std::getline(file, line);
		  std::istringstream iss(line);
		  iss >> dirichaddress_string;
		  if(dirichaddress_string=="#"){
		  	std::string dummy;
		  	std::getline(iss,dummy);
		  	continue;
		  }
		  iss >> channel >> baseline >> width >> thresholdinmV;
		  if(iss.tellg()!=-1) std::cerr << "Error reading line:\n" << line << "\nRead in:" << "\ndirichaddress:0x" << dirichaddress_string << "\nchannel:" << channel << "\nbaseline:" << baseline << "\nwidth:" << width << "\nthresholdinmV:" << thresholdinmV << std::endl;
			else{
      	dirichaddress = (uint16_t)stoi(dirichaddress_string,0,16);
        // std::cout << "0x" << dirichaddress_string << " " << std::hex << dirichaddress << std::dec <<  " " << channel << " " << baseline << std::endl;
				if(dirichlist.count(dirichaddress)!=0){
          if(set_base!=0){
            dirichlist.at(dirichaddress)->SetSingleBaseline_old(channel, dirichlist.at(dirichaddress)->GetSingleBaseline(channel));
            dirichlist.at(dirichaddress)->SetSingleBaseline(channel, baseline);

            dirichlist.at(dirichaddress)->SetSingleNoisewidth_old(channel, dirichlist.at(dirichaddress)->GetSingleNoisewidth(channel));
            dirichlist.at(dirichaddress)->SetSingleNoisewidth(channel, dirich::Thr_mVtoD(width));
          }
          if(set_thr!=0) dirichlist.at(dirichaddress)->SetSingleThresholdmV(channel, thresholdinmV);
        }
				else{ 
					std::cerr << "dirich 0x" << std::hex << dirichaddress << std::dec << " was not found in list of initialized diriches" << std::endl;
					continue;
				}
			}
		}
  	for (auto& dirichlistitem: dirichlist){
  		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
  			if(dirichlistitem.second->GetSingleBaseline(ichannel)==0) std::cerr << "No Baseline for dirich 0x" << std::hex << dirichlistitem.first << std::dec << "'s channel " << ichannel << " found in loading-file" << std::endl;
  		}
  	}
	}
	else if(uselast){
    if(set_base!=0){
      for(int ichannel=0;ichannel<32;++ichannel){
        dirichptr->SetSingleBaseline_old(ichannel, dirichptr->GetSingleBaseline(ichannel));
        dirichptr->SetSingleNoisewidth_old(ichannel, dirichptr->GetSingleNoisewidth(ichannel));
      }
    }
		while(!file.eof()){
			std::string line;
			std::getline(file, line);
		  std::istringstream iss(line);
		  iss >> dirichaddress_string;
		  if(dirichaddress_string=="#"){
		  	std::string dummy;
		  	std::getline(iss,dummy);
		  	continue;
		  }
		  iss >> channel >> baseline >> width >> thresholdinmV;
		  if(iss.tellg()!=-1) std::cerr << "Error reading line:\n" << line << "\nRead in:" << "\ndirichaddress:0x" << dirichaddress_string << "\nchannel:" << channel << "\nbaseline:" << baseline << "\nwidth:" << width << "\nthresholdinmV:" << thresholdinmV << std::endl;
			else{
        if(set_base!=0){
				  dirichptr->SetSingleBaseline(channel, baseline);
				  dirichptr->SetSingleNoisewidth(channel, dirich::Thr_mVtoD(width));
        }
        if(set_thr!=0) dirichptr->SetSingleThresholdmV(channel, thresholdinmV);
			}
		}
	}
  else{
    if(set_base!=0){
      for(int ichannel=0;ichannel<32;++ichannel){
        dirichptr->SetSingleBaseline_old(ichannel, dirichptr->GetSingleBaseline(ichannel));
        dirichptr->SetSingleNoisewidth_old(ichannel, dirichptr->GetSingleNoisewidth(ichannel));
      } 
    }   
    while(!file.eof()){
      std::string line;
      std::getline(file, line);
      std::istringstream iss(line);
		  iss >> dirichaddress_string;
		  if(dirichaddress_string=="#"){
		  	std::string dummy;
		  	std::getline(iss,dummy);
		  	continue;
		  }
		  iss >> channel >> baseline >> width >> thresholdinmV;
    	dirichaddress = (uint16_t)stoi(dirichaddress_string,0,16);
      if(iss.tellg()!=-1) std::cerr << "Error reading line:\n" << line << "\nRead in:" << "\ndirichaddress:0x" << dirichaddress_string << "\nchannel:" << channel << "\nbaseline:" << baseline << "\nwidth:" << width << "\nthresholdinmV:" << thresholdinmV << std::endl;
      else if(dirichaddress==dirichptr->GetBoardAddress()){
        if(set_base!=0){
          dirichptr->SetSingleBaseline(channel, baseline);
				  dirichptr->SetSingleNoisewidth(channel, dirich::Thr_mVtoD(width));
        }
        if(set_thr!=0) dirichptr->SetSingleThresholdmV(channel, thresholdinmV);
      }
    }
  }
}

void save_graphs(std::shared_ptr<dirich>  dirichptr, std::string filename){
	TFile* file=new TFile(Form("%s.root", filename.c_str()),"RECREATE");
  if(dirichptr==NULL){
  	file->cd();
		get_noisewidth_histo(0)->Write();
		get_2D_rate_histo(0)->Write();
		get_2D_rate_over_thr_histo(0)->Write();
		get_2D_diff_over_thr_histo(0)->Write();
		get_2D_gr_diff_over_thr_histo(0)->Write();
		get_2D_mgr_diff_over_thr_histo(0)->Write();
    for (auto& dirichlistitem: dirichlist) {
    	std::cout << "saving graphs of 0x" << std::hex << dirichlistitem.first << std::dec << std::endl;
			TDirectory *dirich_dir = file->mkdir(Form("dirich_0x%x",dirichlistitem.first));
			dirich_dir->cd();
			get_noisewidth_histo(dirichlistitem.second)->Write();
			get_2D_rate_histo(dirichlistitem.second)->Write();
			get_2D_rate_over_thr_histo(dirichlistitem.second)->Write();			
			get_2D_diff_over_thr_histo(dirichlistitem.second)->Write();
			get_2D_gr_diff_over_thr_histo(dirichlistitem.second)->Write();
			get_2D_mgr_diff_over_thr_histo(dirichlistitem.second)->Write();			
			TDirectory *channels = dirich_dir->mkdir("channels");
			channels->cd();
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				TDirectory *ch = channels->mkdir(Form("ch:%i",ichannel));
				ch->cd();
				dirichlistitem.second->gRateGraphs[ichannel]->Write();
				dirichlistitem.second->gRateGraphsOverBase[ichannel]->Write();
				dirichlistitem.second->gDiffRateGraphsOverBase[ichannel]->Write();
			}
		}
	}
  else{
  	file->cd();
  	TDirectory *dirich_dir = file->mkdir(Form("dirich_0x%x",dirichptr->GetBoardAddress()));
		dirich_dir->cd();
		get_noisewidth_histo(dirichptr)->Write();
		get_2D_rate_histo(dirichptr)->Write();
		get_2D_rate_over_thr_histo(dirichptr)->Write();			
		get_2D_diff_over_thr_histo(dirichptr)->Write();
		get_2D_gr_diff_over_thr_histo(dirichptr)->Write();
		get_2D_mgr_diff_over_thr_histo(dirichptr)->Write();			
		TDirectory *channels = dirich_dir->mkdir("channels");
		channels->cd();
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			TDirectory *ch = channels->mkdir(Form("ch:%i",ichannel));
			ch->cd();
			dirichptr->gRateGraphs[ichannel]->Write();
			dirichptr->gRateGraphsOverBase[ichannel]->Write();
			dirichptr->gDiffRateGraphsOverBase[ichannel]->Write();
		}
  }
  file->Close();
}

void save(){
  std::array<char, 64> buffer;
  buffer.fill(0);
  time_t rawtime;
  time(&rawtime);
  const auto timeinfo = localtime(&rawtime);
  strftime(buffer.data(), sizeof(buffer), "%Y%m%d_%H%M%S", timeinfo);
  std::string str = std::string(buffer.data()) + "_std_save";
  // std::string str = "./save/"+ string(buffer.data()) + "_std_save";

  // auto t = std::time(nullptr);
  // auto tm = *std::localtime(&t);

  // std::ostringstream oss;
  // oss << std::put_time(&tm, "%Y%m%d_%H%M%S") << "_std_save.dico";
  // auto str = oss.str();
  save_base(NULL,str,0);
  save_graphs(NULL,str);
}

void* scanthread_baseline(void* dirichptr) //Argument is pointer to DiRICH class instance
{
  TThread::Printf("Starting baseline for Dirich at address 0x%x",((dirich*)dirichptr)->GetBoardAddress());
	((dirich*)dirichptr)->DoBaselineScan();
  TThread::Printf("Threshold scan for Dirich at address 0x%x done ! ",((dirich*)dirichptr)->GetBoardAddress()); 
  return 0;
}

void* scanthread_nrml(void* dirichptr) //Argument is pointer to DiRICH class instance
{
  TThread::Printf("Starting threshscan for Dirich at address 0x%x",((dirich*)dirichptr)->GetBoardAddress());
	std::cout << "DoThreshScan" << std::endl;
	((dirich*)dirichptr)->DoThreshScan();
	std::cout << "AnalyzeBaseline" << std::endl;
	((dirich*)dirichptr)->AnalyzeBaseline();
	std::cout << "MakeGraphsOverBase" << std::endl;
	((dirich*)dirichptr)->MakeGraphsOverBase();
  std::cout << "MakeDiffGraphsOverBase" << std::endl;
  ((dirich*)dirichptr)->MakeDiffGraphsOverBase();
  TThread::Printf("Threshold scan for Dirich at address 0x%x done ! ",((dirich*)dirichptr)->GetBoardAddress()); 
  return 0;
}

void* scanthread_over(void* dirichptr) //Argument is pointer to DiRICH class instance
{
  TThread::Printf("Starting threshscan_over for Dirich at address 0x%x",((dirich*)dirichptr)->GetBoardAddress());
	((dirich*)dirichptr)->DoThreshScanOverBase();
  ((dirich*)dirichptr)->MakeDiffGraphsOverBase();
  TThread::Printf("Threshold scan for Dirich at address 0x%x done ! ",((dirich*)dirichptr)->GetBoardAddress()); 
  return 0;
}

void* scanthread_thr(void* dirichptr) //Argument is pointer to DiRICH class instance
{
  TThread::Printf("Starting threshsearch for Dirich at address 0x%x",((dirich*)dirichptr)->GetBoardAddress());
  ((dirich*)dirichptr)->DoThreshSearch();
  TThread::Printf("Threshold search for Dirich at address 0x%x done ! ",((dirich*)dirichptr)->GetBoardAddress()); 
  return 0;
}

// void* scanthread_std(void* dirichptr) //Argument is pointer to DiRICH class instance
// {
//   TThread::Printf("Starting threshscan for Dirich at address 0x%x",((dirich*)dirichptr)->GetBoardAddress());
// 	((dirich*)dirichptr)->DoThreshScan();
//   // ((dirich*)dirichptr)->MakeDiffGraphsOverBase();
//   TThread::Printf("Threshold scan for Dirich at address 0x%x done ! ",((dirich*)dirichptr)->GetBoardAddress()); 
//   return 0;
// }

// void single_thr_scan(int type=0, std::shared_ptr<dirich>  dirich_to_scan=NULL){
// 	int ret=0;
// 	if(dirich_to_scan==NULL){
//     std::cerr << "DiRICH not initialized!" << std::endl;
//     return;
//   } 

//   uint32_t TDC_setting[2];
//   ret=trb_register_read(dirich_to_scan->GetBoardAddress(), 0xc802, TDC_setting, 2); //switch off TDC
//   if(ret==-1){
//     std::cerr << "Reading TDCs status failed for dirich " << std::hex << dirich_to_scan->GetBoardAddress() << std::dec << " -> TDC for that dirich will be left switched off" << std::endl;
//   }

//   ret=trb_register_write(dirich_to_scan->GetBoardAddress(), 0xc802, 0x00000000); //switch off TDC
//   if(ret==-1){
//     std::cerr << "Switching off TDCs failed for dirich " << std::hex << dirich_to_scan->GetBoardAddress() << std::dec << " -> Interupting baselinescan" << std::endl;
//     // return;
//   }

// 	TThread* thread;
//   switch(type){
//   case 1:
//   	thread = new TThread(Form("Thread_%i",(int)dirich_to_scan->GetBoardAddress()), scanthread_baseline, (void*) dirich_to_scan);
//   	break;
//   case 2:
//   	thread = new TThread(Form("Thread_%i",(int)dirich_to_scan->GetBoardAddress()), scanthread_thr, (void*) dirich_to_scan);
//   	break;
//   case 3:
//   	thread = new TThread(Form("Thread_%i",(int)dirich_to_scan->GetBoardAddress()), scanthread_over, (void*) dirich_to_scan);
//   	break;
//   case 0:
//   default:
//   	thread = new TThread(Form("Thread_%i",(int)dirich_to_scan->GetBoardAddress()), scanthread_nrml, (void*) dirich_to_scan);
//   	break;
//   }
//   usleep(1000);
//   thread->Run(); 
//   printf("Waiting: \n");
// 	usleep(1000);
  
//   thread->Join();
//   thread->Delete();

//   printf("System scan done ! \n");
//   switch(type){
//   case 1:
//     save();
//   	break;
//   case 2:
//     save();
//   	break;
//   case 3:
//     save();
//   	break;    	
//   case 0:
//   default:
//   	save();
//   	break;
//   }

//   ret=trb_register_write(dirich_to_scan->GetBoardAddress(), 0xc802, TDC_setting[1]);
//   if(ret==-1){
//     std::cerr << "Switching on TDCs failed for dirich: " << std::hex << dirich_to_scan->GetBoardAddress() << std::dec << std::endl;
//   }

// }

void system_thr_scan(int type=0)
{
	int ret=0;

  std::map<uint16_t,uint32_t> TDC_setting;
  std::map<uint16_t,int> TDC_set;
  for (auto& dirichlistitem: dirichlist) {
    if(dirichlistitem.second==NULL){
      std::cerr << "DiRICH " << std::hex << dirichlistitem.first << std::dec << " not initialized! Not switching off TDC" << std::endl;
      TDC_set.insert(std::pair<uint16_t,int>(dirichlistitem.first,1));
      continue;
    }
    if(dirichlistitem.second->IsJansReadout()){
      std::cerr << "DiRICH " << std::hex << dirichlistitem.first << std::dec << " has no TDC" << std::endl;
      TDC_set.insert(std::pair<uint16_t,int>(dirichlistitem.first,2));
      continue;
    } 
    if(dirichlistitem.second->IsSim()){
      std::cerr << "DiRICH " << std::hex << dirichlistitem.first << std::dec << " is only simulated" << std::endl;
      TDC_set.insert(std::pair<uint16_t,int>(dirichlistitem.first,2));
      continue;
    }   
    uint32_t temp_tdc_setting[2];
  	ret=trb_register_read(dirichlistitem.first, 0xc802, temp_tdc_setting, 2); //switch off TDC
    // std::cout << std::hex << temp_tdc_setting[0] << "\t" << temp_tdc_setting[1] << std::endl;
    if(ret!=2 || temp_tdc_setting[0]!=dirichlistitem.first){
      std::cerr << "Reading TDCs status failed for dirich " << std::hex << dirichlistitem.first << std::dec << " -> TDC for that dirich will be left switched off" << std::endl;
      temp_tdc_setting[1] = 0x0;
    }

  	ret=trb_register_write(dirichlistitem.first, 0xc802, 0x00000000); //switch off TDC
    if(ret==-1){
      TDC_set.insert(std::pair<uint16_t,int>(dirichlistitem.first,3));
      continue;
    }
    TDC_set.insert(std::pair<uint16_t,int>(dirichlistitem.first,4));
    TDC_setting.insert(std::pair<uint16_t,uint32_t>(dirichlistitem.first,temp_tdc_setting[1]));
	}

  std::vector <TThread*> threadlist;  
  // Initialize instances of dirich class for each module
  for (auto& dirichlistitem: dirichlist){
    if(dirichlistitem.second==NULL){
      std::cerr << "DiRICH " << std::hex << dirichlistitem.first << std::dec << " not initialized!" << std::endl;
      continue;
    } 
    if(TDC_set.at(dirichlistitem.first)==3){
      std::cerr << "Switching off TDCs failed for dirich " << std::hex << dirichlistitem.first << std::dec << " -> Skipping" << std::endl;
      continue;
    }    
    switch(type){
    case 1:
    	threadlist.push_back(new TThread(Form("Thread_%i",(int)dirichlistitem.first), scanthread_baseline, (void*) dirichlistitem.second.get()));
    	break;
    case 2:
    	threadlist.push_back(new TThread(Form("Thread_%i",(int)dirichlistitem.first), scanthread_thr, (void*) dirichlistitem.second.get()));
    	break;
    case 3:
    	threadlist.push_back(new TThread(Form("Thread_%i",(int)dirichlistitem.first), scanthread_over, (void*) dirichlistitem.second.get()));
    	break;
    case 0:
    default:
    	threadlist.push_back(new TThread(Form("Thread_%i",(int)dirichlistitem.first), scanthread_nrml, (void*) dirichlistitem.second.get()));
    	break;
    }
    usleep(1000);
    threadlist.back()->Run(); 
  }
  // cout << threadlist.size() << std::endl;
  printf("Waiting: \n");
  usleep(1000);

  for(auto& thread : threadlist){
      // cout << thread->GetState() << std::endl;
      // thread.second->Join();
      thread->Join();
      thread->Delete();
  } 
  // threadlist.clear();
  printf("System scan done ! \n");
  switch(type){
  case 1:
    save();
  	break;
  case 2:
  	break;
  case 3:
    save();
  	break;    	
  case 0:
  default:
  	save();
  	break;
  }

	for (auto& TDC_setting_item: TDC_setting) {
    if(TDC_set.at(TDC_setting_item.first)!=4) continue;
		ret=trb_register_write(TDC_setting_item.first, 0xc802, TDC_setting_item.second); //switch off TDC
		if(ret==-1){
			std::cerr << "Switching on TDCs failed for dirich: " << std::hex << TDC_setting_item.first << std::dec << std::endl;
		}
	}
}

void initialize_diriches(bool search_dirich){
// void initialize_diriches(bool search_dirich, std::vector<int> ranges, int NrPasses, double meas_time){
  TH1::AddDirectory(0);
  int ret=0;
  ret=init_ports();

  if(ret==-1){
    std::cerr << "failed to initialize trb-net ports" << std::endl;
  }
	dirichlist.clear();

  if(search_dirich){
    int dirich_counter=0;
    // const size_t size4mb = 4194304;
    const size_t size4mb = 8000; //sufficient for 2000 DiRICHes
    uint32_t buffer[size4mb];
    for(int i=0;i<100;++i){
      TRBAccessMutex.Lock();
      ret=trb_read_uid(0xfe51, buffer, size4mb);
      TRBAccessMutex.UnLock();
      if(ret>0) break;
    }
    if(ret<0){
      std::cerr << "No TRB3 Modules found!!!" << std::endl;
      return;
    }
    for(int i=0;i<ret;i+=4){
      // if(buffer[i+3]>0x1200 && buffer[i+3]<0x1200)
      dirichlist.insert(std::make_pair(uint16_t(buffer[i+3]),std::shared_ptr<dirich>(new dirich(uint16_t(buffer[i+3])))));
      ++dirich_counter;
      std::cout << "Created DiRICH-Object for DiRICH with address: 0x" << std::hex << buffer[i+3] << std::endl;
      // std::cout << dirichlist.at(uint16_t(buffer[i+3]))->GetBoardUID() << std::endl;
      if(dirichlist.at(uint16_t(buffer[i+3]))->gMeasureTime!=.3){ //pls change it according to your initialization... Sure one should rather throw during init... but well I am lazy
      	std::cout << "DiRICH 0x" << std::hex << uint16_t(buffer[i+3]) << " not correclty initialized. Deleting!" << std::endl;
        // delete dirichlist.at(uint16_t(buffer[i+3]));
        dirichlist.erase(uint16_t(buffer[i+3]));
      }
      // std::cout << dirichlist.at(uint16_t(buffer[i+3]))->gLowerEdge << std::endl;
      // std::cout << dirichlist.at(uint16_t(buffer[i+3]))->gUpperEdge << std::endl;
    }
    std::cout << "Found " << std::dec << dirich_counter << " different diriches\nInitialized " << dirichlist.size() << " out of those" << std::endl;
  }
}

void setup_scan_parameters(std::shared_ptr<dirich>  dirichptr, double gMeasureTime, int gLowerEdge, int gUpperEdge, int gStepsize, int gNrPasses){
	if(dirichptr==NULL){
		for (auto& dirichlistitem: dirichlist) {
			if(dirichlistitem.second==NULL){
				std::cerr << "dirich 0x" << std::hex << dirichlistitem.first << std::dec << " not initialized" << std::endl;
				continue;
			}
			dirichlistitem.second->gMeasureTime = gMeasureTime;
			dirichlistitem.second->gLowerEdge = gLowerEdge;
			dirichlistitem.second->gUpperEdge = gUpperEdge;
			dirichlistitem.second->gStepsize = gStepsize;
			dirichlistitem.second->gNrPasses = gNrPasses;
		}
	}
	else{
		dirichptr->gMeasureTime = gMeasureTime;
		dirichptr->gLowerEdge = gLowerEdge;
		dirichptr->gUpperEdge = gUpperEdge;
		dirichptr->gStepsize = gStepsize;
		dirichptr->gNrPasses = gNrPasses;
	}
}

void setup_scan_parameters_over_thr_mV(std::shared_ptr<dirich>  dirichptr, double gMeasureTime, double gUpperEdgemV, double gStepsizemV, int gNrPasses){
	if(dirichptr==NULL){
		for (auto& dirichlistitem: dirichlist) {
			if(dirichlistitem.second==NULL){
				std::cerr << "dirich 0x" << std::hex << dirichlistitem.first << std::dec << " not initialized" << std::endl;
				continue;
			}
			dirichlistitem.second->gMeasureTime_over = gMeasureTime;
			dirichlistitem.second->gUpperEdge_over = gUpperEdgemV;
			dirichlistitem.second->gStepsize_over = gStepsizemV;
			dirichlistitem.second->gNrPasses_over = gNrPasses;
		}
	}
	else{
		dirichptr->gMeasureTime_over = gMeasureTime;
		dirichptr->gUpperEdge_over = gUpperEdgemV;
		dirichptr->gStepsize_over = gStepsizemV;
		dirichptr->gNrPasses_over = gNrPasses;
	}
}


int main(int argc, char* argv[]){
  std::string loading_file = "";
  std::string save_file = "";
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("info,i", "only initialize diriches")
      // ("use-dirich,u", po::value<std::vector<std::string>>()->multitoken(), "Add diriches to the dirichlist. If dirich does not exist it will be simulated")
      // ("dont-search-dirich", "Dont't search for active diriches during startup.")
      // ("scan-baseline-new,n", "Do new baselinescan (Michaels method). No parameters need to be given!")
      ("scan-baseline,b", po::value<std::vector<std::string>>()->multitoken(), "Do standard baselinescan. Six parameters need to be given:dirich (if 0, all diriches), measure-time (s), threshold-start-value, threshold-end-value, threshold-step-width, number of cycles (two refers to every second channel measured at a time. Values will be set for all diriches!")
      ("find-threshold,r", "find the perfect threshold for the given dirich/maptm-channel-combination. No parameters need to be given.")
      // ("find-threshold,i", po::value<double>(),"find the perfect threshold for the given dirich/maptm-channel-combination. The parameter specifies the method to find the perfect threshold:\n0: searches for the minimum in the differentiated spectrum or for the minimal gradient\n0<value<5: tries to find peak and sigma of the single photon distribution and sets the threshold to value*sigma (!!!!currently not implemented!!!)\n5<value<100: tries to find the single photon peak and sets the threshold to value% of the spp-position")
      ("scan-above-noise,a", po::value<std::vector<std::string>>()->multitoken(), "Do scan for threshold-values greater than the diriches noiseband. Five parameters need to be given:dirich (if 0, all diriches), measure-time (s), threshold-end-value (mV), threshold-step-width (mV), number of cycles (two refers to every second channel measured at a time")
      ("load-baseline,l", po::value<std::vector<std::string>>()->multitoken(), "This option loads the baseline from the file specified in --loading-file. If no file was specified, the latest produced file is choosen. One can specify a certain dirich by using this options parameter. Be aware that this option overwrites the baseline retreived from the baselinescan")
      ("load-threshold", po::value<std::vector<std::string>>()->multitoken(), "This option loads the threshold from the file specified in --loading-file. If no file was specified, the latest produced file is choosen. One can specify a certain dirich by using this options parameter. Be aware that the thresholds are overwriten by --set-threshold")
      ("loading-file,f", po::value<std::string>(&loading_file)->default_value(""), "File to load thresholds and/or baseline from")
      ("save,s", po::value<std::vector<std::string>>()->multitoken(), "Save histograms and data of specified dirich after everything else is executed! Autosaves will be still produced and saved via \"DATE_std_save{.thr,.root}\". Savefile can be set via --save-file")
      ("save-file", po::value<std::string>(&save_file)->default_value(""), "Save histograms and data. If no file specified, a std. filename will be produced")
      ("draw-scan-baseline,d", po::value<std::vector<std::string>>()->multitoken(), "Draw the results of the baselinescan. Dirich can be specified using this options parameter. Obviously this function fails if no scan was done!")
      ("draw-scan-above-noise", po::value<std::vector<std::string>>()->multitoken(), "Draw the results of the thresholdscan above the diriches noiseband. Dirich can be specified using this options parameter. Obviously this function fails if no scan was done!")
      ("draw-scan-above-noise-diff-gr", po::value<std::vector<std::string>>()->multitoken(), "Draw the results of the baselinescan above the diriches noiseband as differential plot. Dirich can be specified using this options parameter. Obviously this function fails if no scan was done!")
      ("draw-noisewidth,w", po::value<std::vector<std::string>>()->multitoken(), "Draw the noisewidth. Dirich can be specified using this options parameter. Obviously this function fails if neither a scan was done nor a threshold-setting was loaded!")
      ("set-threshold,t", po::value<std::vector<std::string>>()->multitoken(), "Set threshold for specified diriches in mV. First Parameter specifies the dirich (0 equals all dirichs), the second the threshold. Only positive threshold values are accepted, as the minus-sign induces errors.")
  ;
// implicit_value(std::vector<std::string>{"0"},"0")
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }
  if(loading_file==""){
    fs::path latest;
    std::time_t latest_tm {};
    for (auto&& entry : boost::make_iterator_range(fs::directory_iterator("."), {})) {
      fs::path p = entry.path();
      if (is_regular_file(p) && p.extension() == ".thr") 
      {
        std::time_t timestamp = fs::last_write_time(p);
        if (timestamp > latest_tm) {
          latest = p;
          latest_tm = timestamp;
        }
      }
    }
    loading_file = latest.filename().string();
  }
  // if(vm.count("use-diriches")){
  //   for(auto& use_diriches_options : vm["use-diriches"].as<std::vector<std::string>>()){
  //     dirichlist.emplace(std::stoi(use_diriches_options.substr(use_diriches_options.find("0x")!=std::string::npos ? use_diriches_options.find("0x")+2 : 0),NULL,16),(dirich*)NULL);
  //   }
  // }

  // if(!vm.count("dont-search-dirich")) initialize_diriches(1);
  // else initialize_diriches(0);
  
  initialize_diriches(1);
  if(vm.count("info")){
    return 0;
  }
  // if(vm.count("scan-baseline-new")){
  //   system_thr_scan(1);

  //   if(vm.count("draw-scan-baseline")){
  //     // std::cout << "draw-scan-baseline" << std::endl;
  //     if(vm["draw-scan-baseline"].empty()){
  //       draw_histo(get_2D_rate_histo(NULL),NULL);
  //     }
  //     for(auto& draw_scan_baseline_options : vm["draw-scan-baseline"].as<std::vector<std::string>>()){
  //       // std::cout << "input for draw-scan-baseline: " << draw_scan_baseline_options << std::endl;
  //       if(draw_scan_baseline_options=="0") draw_histo(get_2D_rate_histo(NULL),NULL);
  //       else draw_histo(get_2D_rate_histo(dirichlist.at(std::stoi(draw_scan_baseline_options.substr(draw_scan_baseline_options.find("0x")!=std::string::npos ? draw_scan_baseline_options.find("0x")+2 : 0),NULL,16))),NULL);
  //     }
  //   }
  //   if(vm.count("draw-noisewidth")){
  //     if(vm["draw-noisewidth"].empty()){
  //       draw_histo(get_noisewidth_histo(NULL),NULL);
  //     }
  //     for(auto& draw_noisewidth_options : vm["draw-noisewidth"].as<std::vector<std::string>>()){
  //       if(draw_noisewidth_options=="0") draw_histo(get_noisewidth_histo(NULL),NULL);
  //       else draw_histo(get_noisewidth_histo(dirichlist.at(std::stoi(draw_noisewidth_options.substr(draw_noisewidth_options.find("0x")!=std::string::npos ? draw_noisewidth_options.find("0x")+2 : 0),NULL,16))),NULL);
  //     }
  //   }
  // }

  if(vm.count("scan-baseline")){
    if(vm["scan-baseline"].empty() || (vm["scan-baseline"].as<std::vector<std::string>>()).size() < 6){
      std::cout << "no or less than six arguments were provided for option --scan-baseline:\nrunning scan with std. parameters" << std::endl;    
    }
    else{
      std::vector<std::vector<std::string>> each_scan_base_opt;
      std::vector<std::string> temp_vec;
      for(auto& scan_base_opt : vm["scan-baseline"].as<std::vector<std::string>>()){
        std::cout << scan_base_opt << std::endl;
        if(scan_base_opt.find("0x")!=std::string::npos || scan_base_opt=="0"){
          // std::cout << "0 or 0x" <<std::endl;
          if(temp_vec.size()==6){
            each_scan_base_opt.push_back(temp_vec);
          }
          temp_vec.clear();
        }
        temp_vec.push_back(scan_base_opt);
      }
      if(temp_vec.size()==6){
        each_scan_base_opt.push_back(temp_vec);
      }      
      for(auto& one_scan_base_opt : each_scan_base_opt){
        std::cout << std::stoi(one_scan_base_opt.at(0).substr(one_scan_base_opt.at(0).find("0x")!=std::string::npos ? one_scan_base_opt.at(0).find("0x")+2 : 0),NULL,16) << "\t"
                  << std::stod(one_scan_base_opt.at(1)) << "\t"
                  << std::stoi(one_scan_base_opt.at(2)) << "\t"
                  << std::stoi(one_scan_base_opt.at(3)) << "\t"
                  << std::stoi(one_scan_base_opt.at(4)) << "\t"
                  << std::stoi(one_scan_base_opt.at(5)) 
                  << std::endl;
        setup_scan_parameters(one_scan_base_opt.at(0) == "0" ? 0 : dirichlist.at(std::stoi(one_scan_base_opt.at(0).substr(one_scan_base_opt.at(0).find("0x")!=std::string::npos ? one_scan_base_opt.at(0).find("0x")+2 : 0),NULL,16))
          , std::stod(one_scan_base_opt.at(1))
          , std::stoi(one_scan_base_opt.at(2))
          , std::stoi(one_scan_base_opt.at(3))
          , std::stoi(one_scan_base_opt.at(4))
          , std::stoi(one_scan_base_opt.at(5))
        );
      }
    }
    system_thr_scan(0);

    if(vm.count("draw-scan-baseline")){
      // std::cout << "draw-scan-baseline" << std::endl;
      if(vm["draw-scan-baseline"].empty()){
        draw_histo(get_2D_rate_histo(NULL),NULL);
      }
      for(auto& draw_scan_baseline_options : vm["draw-scan-baseline"].as<std::vector<std::string>>()){
        // std::cout << "input for draw-scan-baseline: " << draw_scan_baseline_options << std::endl;
        if(draw_scan_baseline_options=="0") draw_histo(get_2D_rate_histo(NULL),NULL);
        else draw_histo(get_2D_rate_histo(dirichlist.at(std::stoi(draw_scan_baseline_options.substr(draw_scan_baseline_options.find("0x")!=std::string::npos ? draw_scan_baseline_options.find("0x")+2 : 0),NULL,16))),NULL);
      }
    }
    if(vm.count("draw-noisewidth")){
      if(vm["draw-noisewidth"].empty()){
        draw_histo(get_noisewidth_histo(NULL),NULL);
      }
      for(auto& draw_noisewidth_options : vm["draw-noisewidth"].as<std::vector<std::string>>()){
        if(draw_noisewidth_options=="0") draw_histo(get_noisewidth_histo(NULL),NULL);
        else draw_histo(get_noisewidth_histo(dirichlist.at(std::stoi(draw_noisewidth_options.substr(draw_noisewidth_options.find("0x")!=std::string::npos ? draw_noisewidth_options.find("0x")+2 : 0),NULL,16))),NULL);
      }
    }
  }

  if(vm.count("load-baseline")){
    // std::cout << "inside load baseline" << std::endl;
    if(loading_file==""){
      std::cout << "no loading-file found!\n! aborting !" << std::endl;
    }
    else{
    	std::cout << "loading_file: " << loading_file << std::endl;
      for(auto& load_baseline_opt : vm["load-baseline"].as<std::vector<std::string>>()){
        if(load_baseline_opt=="0") load_base(NULL, loading_file, 0, 1, 0);
        else load_base(dirichlist.at(std::stoi(load_baseline_opt.substr(load_baseline_opt.find("0x")!=std::string::npos ? load_baseline_opt.find("0x")+2 : 0),NULL,16)), loading_file, 0, 1, 0);
      }
      if(vm.count("draw-noisewidth")){
        if(vm["draw-noisewidth"].empty()){
          draw_histo(get_noisewidth_histo(NULL),NULL);
        }        
        for(auto& draw_noisewidth_options : vm["draw-noisewidth"].as<std::vector<std::string>>()){
          if(draw_noisewidth_options=="0") draw_histo(get_noisewidth_histo(NULL),NULL);
          else draw_histo(get_noisewidth_histo(dirichlist.at(std::stoi(draw_noisewidth_options.substr(draw_noisewidth_options.find("0x")!=std::string::npos ? draw_noisewidth_options.find("0x")+2 : 0),NULL,16))),NULL);
        }
      }
    }
  }

  if(vm.count("scan-above-noise")){
    if(vm["scan-above-noise"].empty() || (vm["scan-above-noise"].as<std::vector<std::string>>()).size() < 5){
      std::cout << "no or less than five arguments were provided for option --scan-above-noise:\nrunning scan with std. parameters" << std::endl;    
    }
    else{
      std::vector<std::vector<std::string>> each_scan_above_noise_opt;
      std::vector<std::string> temp_vec;
      for(auto& scan_above_noise_opt : vm["scan-above-noise"].as<std::vector<std::string>>()){
        if(scan_above_noise_opt.find("0x")!=std::string::npos || scan_above_noise_opt=="0"){
          if(temp_vec.size()==5){
            each_scan_above_noise_opt.push_back(temp_vec);
          }
          temp_vec.clear();
        }
        temp_vec.push_back(scan_above_noise_opt);
      }
      if(temp_vec.size()==5){
        each_scan_above_noise_opt.push_back(temp_vec);
      }   
      for(auto& one_scan_above_noise_opt : each_scan_above_noise_opt){
        std::cout << std::stoi(one_scan_above_noise_opt.at(0).substr(one_scan_above_noise_opt.at(0).find("0x")!=std::string::npos ? one_scan_above_noise_opt.at(0).find("0x")+2 : 0),NULL,16) << "\t"
                  << std::stod(one_scan_above_noise_opt.at(1)) << "\t"
                  << std::stod(one_scan_above_noise_opt.at(2)) << "\t"
                  << std::stod(one_scan_above_noise_opt.at(3)) << "\t"
                  << std::stoi(one_scan_above_noise_opt.at(4)) 
                  << std::endl;
        setup_scan_parameters_over_thr_mV(one_scan_above_noise_opt.at(0) == "0" ? 0 : dirichlist.at(std::stoi(one_scan_above_noise_opt.at(0).substr(one_scan_above_noise_opt.at(0).find("0x")!=std::string::npos ? one_scan_above_noise_opt.at(0).find("0x")+2 : 0),NULL,16))
          , std::stod(one_scan_above_noise_opt.at(1))
          , std::stod(one_scan_above_noise_opt.at(2))
          , std::stod(one_scan_above_noise_opt.at(3))
          , std::stoi(one_scan_above_noise_opt.at(4))
        );
      }
    }
    system_thr_scan(3);
  }

  if(vm.count("find-threshold")){
    for(auto& dirich : dirichlist){
      dirich.second->gThreshold_finding_method = 0;
      // dirich.second->gThreshold_finding_method = vm["find-threshold"].as<double>();
      std::cout << dirich.second->gThreshold_finding_method << std::endl;
    }
    system_thr_scan(2);
  }

  if(vm.count("load-threshold")){
    for(auto& load_threshold_options : vm["load-threshold"].as<std::vector<std::string>>()){
      if(load_threshold_options=="0") load_base(NULL, loading_file,0, 0, 1);
      else load_base(dirichlist.at(std::stoi(load_threshold_options.substr(load_threshold_options.find("0x")!=std::string::npos ? load_threshold_options.find("0x")+2 : 0),NULL,16)), loading_file, 0, 0, 1);
    }
  }

  if(vm.count("set-threshold")){
    if(vm["set-threshold"].empty() || (vm["set-threshold"].as<std::vector<std::string>>()).size() < 2){
      std::cout << "no or less than two arguments were provided for option --set-threshold:\nno thresholds will be set" << std::endl;    
    }
    else{
      std::vector<std::vector<std::string>> each_set_threshold_opt;
      std::vector<std::string> temp_vec;
      for(auto& set_threshold_opt : vm["set-threshold"].as<std::vector<std::string>>()){
        if(set_threshold_opt.find("0x")!=std::string::npos || set_threshold_opt=="0"){
          if(temp_vec.size()==2){
            each_set_threshold_opt.push_back(temp_vec);
          }
          temp_vec.clear();
        }
        temp_vec.push_back(set_threshold_opt);
      }
      if(temp_vec.size()==2){
        each_set_threshold_opt.push_back(temp_vec);
      }   
      for(auto& one_set_threshold_opt : each_set_threshold_opt){
        std::cout << std::stoi(one_set_threshold_opt.at(0).substr(one_set_threshold_opt.at(0).find("0x")!=std::string::npos ? one_set_threshold_opt.at(0).find("0x")+2 : 0),NULL,16) << "\t"
                  << std::stod(one_set_threshold_opt.at(1))
                  << std::endl;
        set_thresholds(one_set_threshold_opt.at(0) == "0" ? 0 : dirichlist.at(std::stoi(one_set_threshold_opt.at(0).substr(one_set_threshold_opt.at(0).find("0x")!=std::string::npos ? one_set_threshold_opt.at(0).find("0x")+2 : 0),NULL,16)), std::stod(one_set_threshold_opt.at(1)));
      }
    }
  }

  if(vm.count("scan-above-noise")){
    if(vm.count("draw-scan-above-noise")){
      if(vm["draw-scan-above-noise"].empty()){
        draw_histo(get_2D_rate_over_thr_histo(NULL),NULL);
      }      
      for(auto& draw_scan_above_noise_options : vm["draw-scan-above-noise"].as<std::vector<std::string>>()){
        if(draw_scan_above_noise_options=="0") draw_histo(get_2D_rate_over_thr_histo(NULL),NULL);
        else draw_histo(get_2D_rate_over_thr_histo(dirichlist.at(std::stoi(draw_scan_above_noise_options.substr(draw_scan_above_noise_options.find("0x")!=std::string::npos ? draw_scan_above_noise_options.find("0x")+2 : 0),NULL,16))),NULL);
      }
    }
    if(vm.count("draw-scan-above-noise-diff-gr")){
      if(vm["draw-scan-above-noise-diff-gr"].empty()){
        draw_multigraph(get_2D_mgr_diff_over_thr_histo(NULL),NULL);
      }      
      for(auto& draw_scan_above_noise_diff_options : vm["draw-scan-above-noise-diff-gr"].as<std::vector<std::string>>()){
        if(draw_scan_above_noise_diff_options=="0") draw_multigraph(get_2D_mgr_diff_over_thr_histo(NULL),NULL);
        else draw_multigraph(get_2D_mgr_diff_over_thr_histo(dirichlist.at(std::stoi(draw_scan_above_noise_diff_options.substr(draw_scan_above_noise_diff_options.find("0x")!=std::string::npos ? draw_scan_above_noise_diff_options.find("0x")+2 : 0),NULL,16))),NULL);
      }
    }
  }

  if(save_file==""){
    std::array<char, 64> buffer;
    buffer.fill(0);
    time_t rawtime;
    time(&rawtime);
    const auto timeinfo = localtime(&rawtime);
    strftime(buffer.data(), sizeof(buffer), "%Y%m%d_%H%M%S", timeinfo);
    save_file = std::string(buffer.data()) + "_std_save";      
  }
  else{
    save_file=vm["save-file"].as<std::string>();
  }

  if(vm.count("save")){
    for(auto& save_options : vm["save"].as<std::vector<std::string>>()){
      if(save_options=="0"){
        save_base(NULL,save_file,1);
        save_graphs(NULL,save_file);
      }
      else{
        save_base(dirichlist.at(std::stoi(save_options.substr(save_options.find("0x")!=std::string::npos ? save_options.find("0x")+2 : 0),NULL,16)),save_file,1);
        save_graphs(dirichlist.at(std::stoi(save_options.substr(save_options.find("0x")!=std::string::npos ? save_options.find("0x")+2 : 0),NULL,16)),save_file);
      }
    }
  }

  std::string str = save_file + "_all_canvases.pdf";
  uint counter=0;
  for(auto& canvases : canvasvector){
    if(counter==0 && canvasvector.size()>1) canvases->Print(Form("%s(",str.c_str()));
    else if(counter==canvasvector.size()-1) canvases->Print(Form("%s)",str.c_str()));
    else canvases->Print(str.c_str());
    counter++;
  }

  return 0;
}
