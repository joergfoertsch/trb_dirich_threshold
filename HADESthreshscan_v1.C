// #include "trbnetcom.h"
// #include "dirich_sim.C"

#include "TROOT.h"
#include "TError.h"
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

#include <iostream>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <map>
#include <array>
#include <sstream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <future>
// #include <functional>

#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/range.hpp>
// #include <iomanip>

#include "dirich_v13.C"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int gcheck_thresholds = 1;
std::mutex gcheck_thresholds_mutex;
// #define 0 0
// #define LASTCHANNEL 31

// #ifndef NCH
// 	const int NRCHANNELS = 32;			 //Nr of TDC ichannels in dirich
// 	const int CHPCHAIN = 16;			 //Nr of TDC ichannels pre dirich-chain
// 	#define NCH
// #endif

// #ifndef THC
// 	const int OFFTHRESH =	1;	 //Value to switch off channel
// 	const int THRESHDELAY = 100000;	//Delay [mus] for thresh change to succeed
// 	#define THC
// #endif

const uint16_t BROADCAST = 0xfe51;

std::map<uint16_t,std::shared_ptr<dirich>> dirichlist ={};

std::map<uint16_t,TCanvas*> canvaslist;

std::vector<TCanvas*> canvasvector;

TH2* get_2D_rate_histo(std::shared_ptr<dirich>	dirichptr)
{
	TH2D* histo;
	gStyle->SetOptStat(0);
	if(dirichptr==NULL){
		histo = new TH2D(
			"2D Rate vs. Threshold of all diriches","2D Rate vs. Threshold of all diriches",
			dirichlist.size()*NRCHANNELS,-.5,dirichlist.size()*NRCHANNELS-.5,
			(
				dirichlist.begin()->second->gUpperEdge.at(0)
				-dirichlist.begin()->second->gLowerEdge.at(0)
			)/dirichlist.begin()->second->gStepsize,
			dirichlist.begin()->second->gLowerEdge.at(0),
			dirichlist.begin()->second->gUpperEdge.at(0)
		);
		int idirich=0;
		// std::map<uint16_t,dirich*>::iterator dirichlistiterator = dirichlist.begin();
		TLine* dirich_line_left = new TLine(
			histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),
			histo->GetYaxis()->GetBinLowEdge(1),
			histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),
			histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
		);
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
		for (auto& dirichitem : dirichlist){
			TLine* dirich_line_right = new TLine(
				histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1),
				histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
			);
			dirich_line_right->SetLineWidth(2);
			dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
			TText* dirich_name = new TText(
				histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1)
				-(0.1*(
					histo->GetYaxis()->GetBinLowEdge(1)
					-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
				)),
				Form("0x%x",dirichitem.first)
			);
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				for(int ipoint=0;ipoint<dirichitem.second->gRateGraphs[ichannel]->GetN();++ipoint){
					histo->Fill(
						idirich*NRCHANNELS+ichannel,
						dirichitem.second->gRateGraphs[ichannel]->GetX()[ipoint],
						dirichitem.second->gRateGraphs[ichannel]->GetY()[ipoint]
					);
					if(ichannel%8==0) 
						histo->GetXaxis()->SetBinLabel(
							idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel)
						);
					else histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,"");
				}
				TLine* baseline_line = new TLine(
					histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),
					dirichitem.second->GetSingleBaseline(ichannel),
					histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),
					dirichitem.second->GetSingleBaseline(ichannel)
				);
				baseline_line->SetLineColor(kRed);
				baseline_line->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(baseline_line);
				TLine* baseline_line_old = new TLine(
					histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),
					dirichitem.second->GetSingleBaseline_old(ichannel),
					histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),
					dirichitem.second->GetSingleBaseline_old(ichannel)
				);
				baseline_line_old->SetLineColor(kBlack);
				baseline_line_old->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(baseline_line_old);
			}
			// ++dirichlistiterator;
			++idirich;
		}
	}
	else{
		histo = new TH2D(
			Form("2D Rate vs. Threshold of %x",dirichptr->GetBoardAddress()),
			Form("2D Rate vs. Threshold of %x",dirichptr->GetBoardAddress()),
			NRCHANNELS,
			-.5,
			NRCHANNELS-.5,
			(dirichptr->gUpperEdge.at(0)-dirichptr->gLowerEdge.at(0))/dirichptr->gStepsize,
			dirichptr->gLowerEdge.at(0),
			dirichptr->gUpperEdge.at(0)
		);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(int ipoint=0;ipoint<dirichptr->gRateGraphs[ichannel]->GetN();++ipoint){
				histo->Fill(
					ichannel,
					dirichptr->gRateGraphs[ichannel]->GetX()[ipoint],
					dirichptr->gRateGraphs[ichannel]->GetY()[ipoint]
				);
			}
			TLine* baseline_line = new TLine(
				histo->GetXaxis()->GetBinLowEdge(ichannel+1),
				dirichptr->GetSingleBaseline(ichannel),
				histo->GetXaxis()->GetBinUpEdge(ichannel+1),
				dirichptr->GetSingleBaseline(ichannel)
			);
			baseline_line->SetLineColor(kRed);
			baseline_line->SetLineWidth(2);
			histo->GetListOfFunctions()->Add(baseline_line);
			TLine* baseline_line_old = new TLine(
				histo->GetXaxis()->GetBinLowEdge(ichannel+1),
				dirichptr->GetSingleBaseline_old(ichannel),
				histo->GetXaxis()->GetBinUpEdge(ichannel+1),
				dirichptr->GetSingleBaseline_old(ichannel)
			);
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

TMultiGraph* get_2D_mgr_diff_over_thr_histo(std::shared_ptr<dirich>	dirichptr)
{
	TMultiGraph* multig = new TMultiGraph();
	if(dirichptr==NULL){
		multig->SetTitle(
			"Differentiated rate graph over baseline of all dirich;Threshold;Differentiated rate"
		);
		multig->SetName(
			"Differentiated rate graph over baseline of all dirich (Mutligraph)"
		);
		for (auto& dirichitem : dirichlist){
			for(auto& gDiffRateGraphsOverBaseIT : dirichitem.second->gDiffRateGraphsOverBase){
				multig->Add(gDiffRateGraphsOverBaseIT,"PL");
			}
		}
	}
	else{
		multig->SetTitle(Form(
			"Differentiated rate graph over baseline of dirich 0x%x;Threshold;Differentiated rate"
			,dirichptr->GetBoardAddress()
		));
		multig->SetName(Form(
			"Differentiated rate graph over baseline of dirich 0x%x (Multigraph)"
			,dirichptr->GetBoardAddress()
		));
		for(auto& gDiffRateGraphsOverBaseIT : dirichptr->gDiffRateGraphsOverBase){
			multig->Add(gDiffRateGraphsOverBaseIT,"PL");
		}
	}
	// multig->SetMinimum(0);
	// multig->GetHistogram()->GetYaxis()->SetRangeUser(0,100);
	return multig;
}

TGraph2D* get_2D_gr_diff_over_thr_histo(std::shared_ptr<dirich> dirichptr)
{
	TGraph2D* g2d = new TGraph2D();
	if(dirichptr==NULL){
		g2d->SetTitle(
			"Differentiated rate graph over baseline of all dirich;"
			"Channel Nr;"
			"Threshold;"
			"Differentiated rate"
		);
		g2d->SetName("Differentiated rate graph over baseline of all dirich (2D_Graph)");
		int idirich=0;
		for (auto& dirichitem : dirichlist){
			int ichannel=0;
			for(auto& gDiffRateGraphsOverBaseIT : dirichitem.second->gDiffRateGraphsOverBase){
				for(int ipoint=0;ipoint<gDiffRateGraphsOverBaseIT->GetN();++ipoint){
					g2d->SetPoint(
						g2d->GetN(),
						idirich*NRCHANNELS+ichannel,
						gDiffRateGraphsOverBaseIT->GetX()[ipoint],
						gDiffRateGraphsOverBaseIT->GetY()[ipoint]
					);
				}
				ichannel++;
			}
			idirich++;
		}
	}
	else{
		g2d->SetTitle(Form(
			"Differentiated rate graph over baseline of dirich 0x%x;"
			"Channel Nr;"
			"Threshold;"
			"Differentiated rate"
			,dirichptr->GetBoardAddress())
		);
		g2d->SetName(Form(
			"Differentiated rate graph over baseline of dirich 0x%x (2D_Graph)"
			,dirichptr->GetBoardAddress())
		);
		int ichannel=0;
		for(auto& gDiffRateGraphsOverBaseIT : dirichptr->gDiffRateGraphsOverBase){
			for(int ipoint=0;ipoint<gDiffRateGraphsOverBaseIT->GetN();++ipoint){
				g2d->SetPoint(
					g2d->GetN(),
					ichannel,
					gDiffRateGraphsOverBaseIT->GetX()[ipoint],
					gDiffRateGraphsOverBaseIT->GetY()[ipoint]
				);
			}
			ichannel++;
		}
	}
	g2d->SetMinimum(0);
	g2d->GetZaxis()->SetRangeUser(0,100);
	return g2d;
}

TH2* get_2D_diff_over_thr_histo(std::shared_ptr<dirich>	dirichptr)
{
	TH2D* histo;
	TH2D* divided_histo;
	// divided_histo->SetDirectory(0);
	gStyle->SetOptStat(0);
	if(dirichptr==NULL){
		double max_value=-9999;
		// double min_value=9999;
		double min_width=1000;
		for (auto& dirichitem : dirichlist){
			for(auto& gDiffRateGraphsOverBaseIT : dirichitem.second->gDiffRateGraphsOverBase){
				if(
					gDiffRateGraphsOverBaseIT->GetN()!=0 
					&& max_value<gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1]
				) 
					max_value = gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1];
				if(
					gDiffRateGraphsOverBaseIT->GetN()!=0 
					&& max_value<gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1]
				) 
					max_value = gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1];
				if(
					gDiffRateGraphsOverBaseIT->GetN()>=2 
					&& min_width>abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1])
				) 
					min_width = abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1]);
			}
		}
		histo = new TH2D(
			"2D Differentiated Rate vs. Threshold over baseline of all diriches",
			"2D Differentiated Rate vs. Threshold over baseline of all diriches",
			dirichlist.size()*NRCHANNELS,
			-.5,
			dirichlist.size()*NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		divided_histo = new TH2D(
			"temp_diff",
			"temp_diff",
			dirichlist.size()*NRCHANNELS,
			-.5,dirichlist.size()*NRCHANNELS-.5,
			max_value/min_width/2,
			0,
			max_value
		);
		int idirich=0;
		TLine* dirich_line_left = new TLine(
			histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),
			histo->GetYaxis()->GetBinLowEdge(1),
			histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),
			histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
		);
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
		for (auto& dirichitem : dirichlist){
			TLine* dirich_line_right = new TLine(
				histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1),
				histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
			);
			dirich_line_right->SetLineWidth(2);
			dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
			TText* dirich_name = new TText(
				histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),
					histo->GetYaxis()->GetBinLowEdge(1)
					-(
						0.1*(histo->GetYaxis()->GetBinLowEdge(1)
						-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
					)),
				Form("0x%x",dirichitem.first)
			);
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				for(int ipoint=0;ipoint<dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetN();++ipoint){
					histo->Fill(
						idirich*NRCHANNELS+ichannel,
						dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint],
						dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint]
					);
					divided_histo->Fill(
						idirich*NRCHANNELS+ichannel,
						dirichitem.second->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint]
					);
					// std::cout << ichannel << " " << int(dirichlist.size()*NRCHANNELS/20+1) << std::endl;
					if(ichannel%8==0) histo->GetXaxis()->SetBinLabel(
						idirich*NRCHANNELS+ichannel+1,Form("%i",ichannel)
					);
					else histo->GetXaxis()->SetBinLabel(idirich*NRCHANNELS+ichannel+1,"");
				}
				TLine* thr_line = new TLine(
					histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),
					-1*dirichitem.second->GetSingleThresholdmV(ichannel),
					histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),
					-1*dirichitem.second->GetSingleThresholdmV(ichannel)
				);
				thr_line->SetLineColor(kRed);
				thr_line->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(thr_line);
			}
			++idirich;
		}
		idirich=0;
	}
	else{
		double max_value=0;
		double min_width=10000;
		for(auto& gDiffRateGraphsOverBaseIT : dirichptr->gDiffRateGraphsOverBase){
			if(
				gDiffRateGraphsOverBaseIT->GetN()!=0 
				&& max_value<gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1]
			) 
				max_value = gDiffRateGraphsOverBaseIT->GetX()[gDiffRateGraphsOverBaseIT->GetN()-1];
			if(
				gDiffRateGraphsOverBaseIT->GetN()>=2 
				&& min_width>abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1])
			) 
				min_width = abs(gDiffRateGraphsOverBaseIT->GetX()[0]-gDiffRateGraphsOverBaseIT->GetX()[1]);
		}
		histo = new TH2D(
			Form("2D Differentiated Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),
			Form("2D Differentiated Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),
			NRCHANNELS,
			-.5,
			NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		divided_histo = new TH2D(
			"temp_diff","temp_diff",
			NRCHANNELS,
			-.5,
			NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(int ipoint=0;ipoint<dirichptr->gDiffRateGraphsOverBase[ichannel]->GetN();++ipoint){
				histo->Fill(
					ichannel,
					dirichptr->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint],
					dirichptr->gDiffRateGraphsOverBase[ichannel]->GetY()[ipoint]
				);
				divided_histo->Fill(
					ichannel,
					dirichptr->gDiffRateGraphsOverBase[ichannel]->GetX()[ipoint]
				);
			}
			TLine* thr_line = new TLine(
				histo->GetXaxis()->GetBinLowEdge(ichannel+1),
				-1*dirichptr->GetSingleThresholdmV(ichannel),
				histo->GetXaxis()->GetBinUpEdge(ichannel+1),
				-1*dirichptr->GetSingleThresholdmV(ichannel)
			);
			thr_line->SetLineColor(kRed);
			thr_line->SetLineWidth(2);
			histo->GetListOfFunctions()->Add(thr_line);
		}
	}

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

TH2* get_2D_rate_over_thr_histo(std::shared_ptr<dirich>	dirichptr)
{
	TH2D* histo;
	TH2D* divided_histo;
	// divided_histo->SetDirectory(0);
	gStyle->SetOptStat(0);
	if(dirichptr==NULL){
		double max_value=-9999;
		// double min_value=9999;
		double min_width=1000;
		for (auto& dirichitem : dirichlist){
			for(auto& gRateGraphsOverBaseIT : dirichitem.second->gRateGraphsOverBase){
				if(
					gRateGraphsOverBaseIT->GetN()!=0 
					&& max_value<gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1]
				){
					max_value = gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1];
					// min_value = gRateGraphsOverBaseIT->GetX()[0];
				}
				if(
					gRateGraphsOverBaseIT->GetN()>=2 
					&& min_width>abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1])
				) 
					min_width = abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1]);
			}
		}
		histo = new TH2D(
			"2D Rate vs. Threshold over baseline of all diriches",
			"2D Rate vs. Threshold over baseline of all diriches",
			dirichlist.size()*NRCHANNELS,
			-.5,
			dirichlist.size()*NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		divided_histo = new TH2D(
			"temp_diff",
			"temp_diff",
			dirichlist.size()*NRCHANNELS,
			-.5,
			dirichlist.size()*NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		int idirich=0;
		TLine* dirich_line_left = new TLine(
			histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),
			histo->GetYaxis()->GetBinLowEdge(1),
			histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+1),
			histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
		);
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
		for (auto& dirichitem : dirichlist){
			TLine* dirich_line_right = new TLine(
				histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1),
				histo->GetXaxis()->GetBinLowEdge((idirich+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
			);
			dirich_line_right->SetLineWidth(2);
			dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
			TText* dirich_name = new TText(
				histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1)
				-(0.1*(
						histo->GetYaxis()->GetBinLowEdge(1)
						-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
				)),
				Form("0x%x",dirichitem.first)
			);
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				for(int ipoint=0;ipoint<dirichitem.second->gRateGraphsOverBase[ichannel]->GetN();++ipoint){
					histo->Fill(
						idirich*NRCHANNELS+ichannel,
						dirichitem.second->gRateGraphsOverBase[ichannel]->GetX()[ipoint],
						dirichitem.second->gRateGraphsOverBase[ichannel]->GetY()[ipoint]
					);
					divided_histo->Fill(
						idirich*NRCHANNELS+ichannel,
						dirichitem.second->gRateGraphsOverBase[ichannel]->GetX()[ipoint]
					);
					if(ichannel%8==0) 
						histo->GetXaxis()->SetBinLabel(
							idirich*NRCHANNELS+ichannel+1,
							Form("%i",ichannel)
						);
					else 
						histo->GetXaxis()->SetBinLabel(
							idirich*NRCHANNELS+ichannel+1,""
						);
				}
				TLine* thr_line = new TLine(
					histo->GetXaxis()->GetBinLowEdge(idirich*NRCHANNELS+ichannel+1),
					-1*dirichitem.second->GetSingleThresholdmV(ichannel),
					histo->GetXaxis()->GetBinUpEdge(idirich*NRCHANNELS+ichannel+1),
					-1*dirichitem.second->GetSingleThresholdmV(ichannel)
				);
				thr_line->SetLineColor(kRed);
				thr_line->SetLineWidth(2);
				histo->GetListOfFunctions()->Add(thr_line);
			}
			++idirich;
		}
		idirich=0;
	}
	else{
		double max_value=0;
		double min_width=1000;
		for(auto& gRateGraphsOverBaseIT : dirichptr->gRateGraphsOverBase){
			if(
				gRateGraphsOverBaseIT->GetN()!=0 
				&& max_value<gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1]
			) 
				max_value = gRateGraphsOverBaseIT->GetX()[gRateGraphsOverBaseIT->GetN()-1];
			if(
				gRateGraphsOverBaseIT->GetN()>=2 
				&& min_width>abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1])
			) 
				min_width = abs(gRateGraphsOverBaseIT->GetX()[0]-gRateGraphsOverBaseIT->GetX()[1]);
		}
		histo = new TH2D(
			Form("2D Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),
			Form("2D Rate vs. Threshold over baseline of %x",dirichptr->GetBoardAddress()),
			NRCHANNELS,
			-.5,
			NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		divided_histo = new TH2D(
			"temp_diff",
			"temp_diff",
			NRCHANNELS,
			-.5,
			NRCHANNELS-.5,
			max_value/min_width/2,0,
			max_value
		);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			for(int ipoint=0;ipoint<dirichptr->gRateGraphsOverBase[ichannel]->GetN();++ipoint){
				histo->Fill(
					ichannel,
					dirichptr->gRateGraphsOverBase[ichannel]->GetX()[ipoint],
					dirichptr->gRateGraphsOverBase[ichannel]->GetY()[ipoint]
				);
				divided_histo->Fill(ichannel,dirichptr->gRateGraphsOverBase[ichannel]->GetX()[ipoint]);
			}
			TLine* thr_line = new TLine(
				histo->GetXaxis()->GetBinLowEdge(ichannel+1),
				-1*dirichptr->GetSingleThresholdmV(ichannel),
				histo->GetXaxis()->GetBinUpEdge(ichannel+1),
				-1*dirichptr->GetSingleThresholdmV(ichannel)
			);
			thr_line->SetLineColor(kRed);
			thr_line->SetLineWidth(2);
			histo->GetListOfFunctions()->Add(thr_line);
		}
	}

	histo->Divide(divided_histo);

	histo->SetMinimum(0);
	histo->GetXaxis()->SetTitle("Channel Nr");
	// histo->GetXaxis()->SetTitleOffset();
	histo->GetYaxis()->SetTitle("Threshold");
	histo->GetZaxis()->SetTitle("Rate");
	return histo;
}

TH1* get_noisewidth_histo(std::shared_ptr<dirich>	dirichptr)
{
	TH1* histo;
	if(dirichptr==NULL){
		histo = new TH1D(
			"Noisewidthhistogram of all diriches",
			"Noisewidthhistogram of all diriches",
			dirichlist.size()*NRCHANNELS,
			-.5,
			dirichlist.size()*NRCHANNELS-.5
		);
		int idirich=0;
		for (auto& dirichitem : dirichlist){
			TText* dirich_name = new TText(
				histo->GetXaxis()->GetBinCenter((idirich+1./2)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1)
				-(
					0.1*histo->GetYaxis()->GetBinLowEdge(1)
					-histo->GetYaxis()->GetBinUpEdge(histo->GetNbinsY())
				),
				Form("0x%x",dirichitem.first)
			);
			dirich_name->SetTextAlign(22);
			dirich_name->SetTextColor(kRed+2);
			dirich_name->SetTextFont(43);
			dirich_name->SetTextSize(20);
			histo->GetListOfFunctions()->Add(dirich_name);
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				histo->SetBinContent(
					idirich*NRCHANNELS+ichannel+1,
					dirich::Thr_DtomV(dirichitem.second->GetSingleNoisewidth(ichannel))
				);
				if(ichannel%8==0) 
					histo->GetXaxis()->SetBinLabel(
						idirich*NRCHANNELS+ichannel+1,
						Form("%i",ichannel)
					);
				else 
					histo->GetXaxis()->SetBinLabel(
						idirich*NRCHANNELS+ichannel+1,
						""
					);
			}
		++idirich;
		}
		TLine* dirich_line_left = new TLine(
			histo->GetXaxis()->GetBinLowEdge(1),
			histo->GetYaxis()->GetBinLowEdge(1),
			histo->GetXaxis()->GetBinLowEdge(1),
			histo->GetMaximum()*1.05
		);
		dirich_line_left->SetLineWidth(2);
		dirich_line_left->SetLineColor(kRed);
		histo->GetListOfFunctions()->Add(dirich_line_left);
		for(int i=0;i<idirich;++i){
			TLine* dirich_line_right = new TLine(
				histo->GetXaxis()->GetBinLowEdge((i+1)*NRCHANNELS+1),
				histo->GetYaxis()->GetBinLowEdge(1),
				histo->GetXaxis()->GetBinLowEdge((i+1)*NRCHANNELS+1),
				histo->GetMaximum()*1.05
			);
			dirich_line_right->SetLineWidth(2);
			dirich_line_right->SetLineColor(kRed);
			histo->GetListOfFunctions()->Add(dirich_line_right);
		}
	}
	else{
		histo = new TH1D(
			Form("Noisewidthhistogram of %x",dirichptr->GetBoardAddress()),
			Form("Noisewidthhistogram of %x",dirichptr->GetBoardAddress()),
			NRCHANNELS,
			-.5,
			NRCHANNELS-.5
		);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			histo->SetBinContent(
				ichannel+1,
				(dirich::Thr_DtomV(dirichptr->GetSingleNoisewidth(ichannel)))
			);
		}
	}

	histo->GetYaxis()->SetTitle("NoisewidthinmV");
	histo->SetMinimum(0);
	histo->GetXaxis()->SetTitle("Channel Nr");
	return histo;
}

TH1* get_diff_histo(std::shared_ptr<dirich>	dirichptr, bool baseline1_noisewidth0)
{
	TH1* histo;
	if(dirichptr==NULL){
		if(baseline1_noisewidth0==1) histo = new TH1D(
			"Difference in baseline of all diriches",
			"Difference in baseline of all diriches",
			dirichlist.size()*200,
			-300,
			+300
		);
		else 
			histo = new TH1D(
				"Difference in noisewidth of all diriches",
				"Difference in noisewidth of all diriches",
				dirichlist.size()*200,
				-300,
				+300
			);
		int idirich=0;
		for (auto& dirichitem : dirichlist){
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				if(baseline1_noisewidth0==1) 
					histo->Fill(
						dirichitem.second->GetSingleBaseline(ichannel)
						-dirichitem.second->GetSingleBaseline_old(ichannel)
					);
				else 
					histo->Fill(
						dirichitem.second->GetSingleNoisewidth(ichannel)
						-dirichitem.second->GetSingleNoisewidth_old(ichannel)
					);
			}
		}
		++idirich;
	}
	else{
		if(baseline1_noisewidth0==1) 
			histo = new TH1D(
				Form("Difference in baseline of %x",dirichptr->GetBoardAddress()),
				Form("Difference in baseline of %x",dirichptr->GetBoardAddress()),
				200,
				-300,
				+300
			);
		else 
			histo = new TH1D(
				Form("Difference in noisewidth of %x",dirichptr->GetBoardAddress()),
				Form("Difference in noisewidth of %x",dirichptr->GetBoardAddress()),
				200,
				-300,
				+300
			);
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			if(baseline1_noisewidth0==1) 
				histo->Fill(
					dirichptr->GetSingleBaseline(ichannel)
					-dirichptr->GetSingleBaseline_old(ichannel)
				);
			else 
				histo->Fill(
					dirichptr->GetSingleNoisewidth(ichannel)
					-dirichptr->GetSingleNoisewidth_old(ichannel)
				);
		}
	}

	histo->GetYaxis()->SetTitle("Number of");
	if(baseline1_noisewidth0==1) 
		histo->GetXaxis()->SetTitle("Difference between old and new baseline");
	else 
		histo->GetXaxis()->SetTitle("Difference between old and new noisewidth");
	return histo;
}

void clear_canvas_vector()
{
	while(canvasvector.size()!=0){
		if(canvasvector.back()==NULL){
			std::cout <<	"1" << std::endl;
			canvasvector.pop_back();
		}
		else{
			std::cout <<	"2" << std::endl;
			canvasvector.back()->Clear();
			std::cout <<	"3" << std::endl;
			canvasvector.back()->Close();
			std::cout <<	"4" << std::endl;
			canvasvector.back()->Closed();
			std::cout <<	"5" << std::endl;
			delete canvasvector.back();
			std::cout <<	"6" << std::endl;
			canvasvector.back()=NULL;
			std::cout <<	"7" << std::endl;
			canvasvector.pop_back();
			std::cout <<	"8" << std::endl;
		}
	}
}

void draw_multigraph2D(TMultiGraph* multigraph,TCanvas* canvas)
{
	if(canvas==0){
		canvasvector.emplace_back(
			new TCanvas(
				Form("Canvas%i",(int)canvasvector.size()),
				Form("Canvas%i",(int)canvasvector.size()),
				1920,
				1080
			)
		);
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

void draw_multigraph(TMultiGraph* multigraph,TCanvas* canvas)
{
	if(canvas==0){
		canvasvector.emplace_back(
			new TCanvas(
				Form("Canvas%i",(int)canvasvector.size()),
				Form("Canvas%i",(int)canvasvector.size()),
				1920,
				1080
			)
		);
		canvasvector.back()->cd(0);
	}
	else{
		canvas->cd(0);
	}
	multigraph->Draw("alp");
	double max_value = 0;
	for(auto&& graph : (*multigraph->GetListOfGraphs())){
		for(int i=10 ; i < ((TGraph*)graph)->GetN() ; ++i){
			max_value = 
			((TGraph*)graph)->GetY()[i] > max_value ? 
			((TGraph*)graph)->GetY()[i] : max_value;
		}
	}
	// std::cout << "max_value" << max_value << std::endl;
	multigraph->GetHistogram()->GetYaxis()->SetRangeUser(
		0,
			max_value==0 ? 
			100 : max_value*1.05
	);
	if(canvas==0){
		canvasvector.back()->Modified();
		canvasvector.back()->Update();
	}
	else{
		canvas->Modified();
		canvas->Update();
	}
}

void draw_graph2D(TGraph2D* graph2d,TCanvas* canvas)
{
	if(canvas==0){
		canvasvector.emplace_back(
			new TCanvas(
				Form("Canvas%i",(int)canvasvector.size()),
				Form("Canvas%i",(int)canvasvector.size()),
				1920,
				1080
			)
		);
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

void draw_histo(TH1* histo,TCanvas* canvas)
{
	if(canvas==0){
		canvasvector.emplace_back(
			new TCanvas(
				Form("Canvas%i",(int)canvasvector.size()),
				Form("Canvas%i",(int)canvasvector.size()),
				1920,
				1080
			)
		);
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

void draw_histo(TH2* histo,TCanvas* canvas)
{
	if(canvas==0){
		canvasvector.emplace_back(
			new TCanvas(
				Form("Canvas%i",(int)canvasvector.size()),
				Form("Canvas%i",(int)canvasvector.size()),
				1920,
				1080
			)
		);
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
	if(thrinmV<0.){
		std::cerr 
			<< "negative thresholds are not \"allowed\"!\ninverting value" 
			<< std::endl;
		thrinmV = -1*thrinmV;
	}	
	if(dirichptr==0){
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 2;
		gcheck_thresholds_mutex.unlock();		

		std::vector<std::thread> threads;
		for(auto& dirichlistitem : dirichlist)
			threads.push_back(
				std::thread(
					[&dirichlistitem, thrinmV](){
						dirichlistitem.second->SetThresholdsmV(thrinmV);
					}
				)
			);

		for(auto& thread : threads)
			thread.join();

		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 1;
		gcheck_thresholds_mutex.unlock();
	}
	else if(dirichlist.find(dirichptr->GetBoardAddress())!=dirichlist.end()){
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 2;
		gcheck_thresholds_mutex.unlock();		

		dirichptr->SetThresholdsmV(thrinmV);

		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 1;
		gcheck_thresholds_mutex.unlock();		
	}
	else{
		std::cerr 
			<< "No DiRICH 0x" << std::hex << dirichptr->GetBoardAddress() 
			<< " found" 
			<< std::endl;
	}
}

void set_thresholds_to_noise(std::shared_ptr<dirich> dirichptr, double part_of_noisewidth=1.5)
{
	if(dirichptr==nullptr){
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 2;
		gcheck_thresholds_mutex.unlock();		

		std::vector<std::thread> threads;
		for(auto& dirichlistitem : dirichlist){
			std::array<double,NRCHANNELS> thresholdvals;
			std::array<uint16_t,NRCHANNELS> noisevalues = dirichlistitem.second->GetNoisewidths();
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				thresholdvals.at(ichannel) = 
					part_of_noisewidth*dirich::Thr_DtomV(.5*noisevalues.at(ichannel));
			}

			threads.push_back(
				std::thread(
					[&dirichlistitem, thresholdvals](){
						dirichlistitem.second->SetThresholdsmV(thresholdvals);
					}
				)
			);
		}

		for(auto& thread : threads)
			thread.join();
		
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 1;
		gcheck_thresholds_mutex.unlock();		
	}
	else if(dirichlist.find(dirichptr->GetBoardAddress())!=dirichlist.end()){
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 2;
		gcheck_thresholds_mutex.unlock();

		std::array<double,NRCHANNELS> thresholdvalues;
		std::array<uint16_t,NRCHANNELS> noisevalues = dirichptr->GetNoisewidths();
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			thresholdvalues.at(ichannel) = 
				part_of_noisewidth*dirich::Thr_DtomV(.5*noisevalues.at(ichannel));
		}

		dirichptr->SetThresholdsmV(thresholdvalues);
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 1;
		gcheck_thresholds_mutex.unlock();		
	}
	else{
		std::cerr 
			<< "No DiRICH 0x" << std::hex << dirichptr->GetBoardAddress() 
			<< " found" 
			<< std::endl;
	}
}

void set_pattern(std::shared_ptr<dirich> dirichptr, uint32_t pattern=4294967295)
{
	gcheck_thresholds_mutex.lock();
	gcheck_thresholds = 2;
	gcheck_thresholds_mutex.unlock();		

	if(dirichptr==0){
		std::vector<std::thread> threads;
		for(auto& dirichlistitem : dirichlist){
			threads.push_back(
				std::thread(
					[&dirichlistitem,&pattern](){
						dirichlistitem.second->SetPattern(pattern);
					}
				)
			);
		}

		for(auto& thread : threads)
			thread.join();
	}
	else if(dirichlist.find(dirichptr->GetBoardAddress())!=dirichlist.end()){
		dirichptr->SetPattern(pattern);
	}
	else{
		std::cerr << "No DiRICH 0x" << std::hex << dirichptr->GetBoardAddress() 
			<< " found" 
			<< std::endl;
	}

	gcheck_thresholds_mutex.lock();
	gcheck_thresholds = 1;
	gcheck_thresholds_mutex.unlock();		
}

void measure_rate(std::shared_ptr<dirich>	dirichptr, std::string filename, double measure_time)
{
	if(dirichptr==NULL){
		std::ofstream file;
		file.open(filename, std::ios_base::app);
		if(!file) std::cerr << "File for saving (" << filename << ") could not be opened!" << std::endl;
		
		std::unordered_map<uint16_t,std::future<double*>> rates;
		for(auto& dirich : dirichlist){
			rates.insert(std::pair<uint16_t,std::future<double*>>(
				dirich.first, 
				std::async(std::launch::async,
					&dirich::GetRates, dirich.second.get(), measure_time
				)
			));
		}
		for(auto& one_rates : rates){
			file 
				<< "# Scan-Data\n# dirich\tchannel\trate\terror\t" 
				<< std::endl;
				one_rates.second.wait();
			double* temp_rate_arr = one_rates.second.get();
			for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
				file
					<< std::hex << one_rates.first << std::dec << "\t" 
					<< ichannel << "\t" 
					<< temp_rate_arr[ichannel] << "\t" 
					<< sqrt(temp_rate_arr[ichannel])/sqrt(measure_time) << "\t" 
					<< std::endl;
			}
		}
		if(file)
			file.close();
	}
	else{
		std::ofstream file;
		file.open(filename, std::ios_base::app);
		if(!file) std::cerr << "File for saving (" << filename << ") could not be opened!" << std::endl;
		
		double* rates = dirichptr->GetRates(measure_time);

		file 
			<< "# Scan-Data\n# dirich\tchannel\trate\terror\t" 
			<< std::endl;
		for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
			file
				<< std::hex << dirichptr->GetBoardAddress() << std::dec << "\t" 
				<< ichannel << "\t" 
				<< rates[ichannel] << "\t" 
				<< sqrt(rates[ichannel])/sqrt(measure_time) << "\t" 
				<< std::endl;
		}
		if(file)
			file.close();			
	}
}

void save_base(std::shared_ptr<dirich>	dirichptr, std::string filename, bool append)
{
	std::ofstream file;
	if(append) file.open(filename+".thr", std::ios_base::app);
	else file.open(filename+".thr");

	if(!file) std::cerr << "File for saving (" << filename+".thr" << ") could not be opened!" << std::endl;

	if(dirichptr==NULL){
		std::cout << "saving minimal_data" << std::endl;
		int counter=0;
		for (auto& dirichlistitem: dirichlist){
			std::cout 
				<< "\r" 
				<< std::setw(10) << std::setprecision(2) << std::fixed 
				<< 1.*((counter+1)*100)/(dirichlist.size()) 
				<< "%" << std::flush;
			file 
				<< "# Scan-Settings for 0x" 
				<< std::hex << dirichlistitem.first 
				<< std::dec 
				<< "\n# gMeasureTime\tgLowerEdge(0)\tgUpperEdge(0)\tgStepsize\tgNrPasses\tgMeasureTime_over"
					"\tgUpperEdge_over\tgStepsize_over\tgNrPasses_over" 
				<< std::endl;
			file 
				<< "# " 
				<< dirichlistitem.second->gMeasureTime 
				<< "\t" << dirichlistitem.second->gLowerEdge.at(0) 
				<< "\t" << dirichlistitem.second->gUpperEdge.at(0) 
				<< "\t" << dirichlistitem.second->gStepsize 
				<< "\t" << dirichlistitem.second->gNrPasses 
				<< "\t" << dirichlistitem.second->gMeasureTime_over 
				<< "\t" << dirichlistitem.second->gUpperEdge_over 
				<< "\t" << dirichlistitem.second->gStepsize_over 
				<< "\t" << dirichlistitem.second->gNrPasses_over 
				<< std::endl;
			file 
				<< "# Scan-Data\n# dirich\tchannel\tbaseline\twidth in mV\tthreshold in mV over baseline" 
				<< std::endl;
			for (int ichannel=0; ichannel<NRCHANNELS; ichannel++) {
				file
					<< std::hex << dirichlistitem.first << std::dec << "\t" 
					<< ichannel << "\t" 
					<< dirichlistitem.second->GetSingleBaseline(ichannel) << "\t" 
					<< dirich::Thr_DtomV(dirichlistitem.second->GetSingleNoisewidth(ichannel)) << "\t" 
					<< dirichlistitem.second->GetSingleThresholdmV(ichannel) 
					<< std::endl;
			}
			counter++;
		}
		std::cout << std::endl;
	}
	else{
		file 
			<< "# Scan-Settings for 0x" << std::hex << dirichptr->GetBoardAddress() 
			<< std::dec 
			<< "\n# gMeasureTime\tgLowerEdge(0)\tgUpperEdge(0)\tgStepsize\tgNrPasses\tgMeasureTime_over"
				"\tgUpperEdge_over\tgStepsize_over\tgNrPasses_over" 
			<< std::endl;
		file 
			<< "# " 
			<< dirichptr->gMeasureTime 
			<< "\t" << dirichptr->gLowerEdge.at(0) 
			<< "\t" << dirichptr->gUpperEdge.at(0)
			<< "\t" << dirichptr->gStepsize 
			<< "\t" << dirichptr->gNrPasses 
			<< "\t" << dirichptr->gMeasureTime_over 
			<< "\t" << dirichptr->gUpperEdge_over 
			<< "\t" << dirichptr->gStepsize_over 
			<< "\t" << dirichptr->gNrPasses_over 
			<< std::endl;
		file 
			<< "# Scan-Data\n# dirich\tchannel\tbaseline\twidth in mV\tthreshold in mV over baseline" 
			<< std::endl;
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

void load_base(std::shared_ptr<dirich>	dirichptr, 
	std::string filename, 
	bool uselast, 
	bool set_base, 
	bool set_thr
	)
{
	std::string dirichaddress_string="";
	uint16_t dirichaddress=0;
	int channel=0;
	int baseline=0;
	double width=0;
	double thresholdinmV=0;
	for (auto& dirichlistitem: dirichlist){
		if(dirichlistitem.second==NULL){
			std::cerr 
				<< "dirich 0x" << std::hex << dirichlistitem.first 
				<< std::dec << " was found uninitialized\nRun initialize_diriches(1/0) first!" 
				<< std::endl;
			return;
		}
	}
	std::ifstream file;
	file.open(filename);
	if(!file) std::cerr << "File for loading (" << filename << ") could not be opened!" << std::endl;
	if(dirichptr==nullptr){
		std::unordered_map<uint16_t,std::array<double,32>> thresholds;
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
			// if(dirichaddress_string=="") continue;
			iss >> channel >> baseline >> width >> thresholdinmV;
			if(iss.tellg()!=-1) 
				std::cerr 
					<< "Error reading line:\n" << line 
					<< "\nRead in:" 
					<< "\ndirichaddress:0x" << dirichaddress_string 
					<< "\nchannel:" << channel 
					<< "\nbaseline:" << baseline 
					<< "\nwidth:" << width 
					<< "\nthresholdinmV:" << thresholdinmV 
					<< std::endl;
			else{
				dirichaddress = (uint16_t)stoi(dirichaddress_string,0,16);
				if(dirichlist.count(dirichaddress)!=0){
					if(set_base!=0){
						dirichlist.at(dirichaddress)->SetSingleBaseline_old(
							channel, 
							dirichlist.at(dirichaddress)->GetSingleBaseline(channel)
						);
						dirichlist.at(dirichaddress)->SetSingleBaseline(
							channel, 
							baseline
						);

						dirichlist.at(dirichaddress)->SetSingleNoisewidth_old(
							channel, 
							dirichlist.at(dirichaddress)->GetSingleNoisewidth(channel)
						);
						dirichlist.at(dirichaddress)->SetSingleNoisewidth(
							channel, 
							dirich::Thr_mVtoD(width)
						);
					}
					if(set_thr!=0){
						if(thresholds.find(dirichaddress)==thresholds.end()){
							std::array <double,32> temp_array;
							temp_array.fill(0);
							thresholds.insert(std::make_pair(dirichaddress,temp_array));
						}
						thresholds.at(dirichaddress).at(channel) = thresholdinmV;
					}
				}
				else{ 
					std::cerr 
						<< "dirich 0x" 
						<< std::hex << dirichaddress 
						<< std::dec << " was not found in list of initialized diriches" 
						<< std::endl;
					continue;
				}
			}
		}
		std::vector<std::thread> threads;
		if(set_thr!=0){
			gcheck_thresholds_mutex.lock();
			gcheck_thresholds = 2;
			gcheck_thresholds_mutex.unlock();

			for(auto& one_threshold : thresholds)
				threads.push_back(std::thread(
					[&one_threshold](){
						dirichlist.at(one_threshold.first)->SetThresholdsmV(
							one_threshold.second
						);
					}
				));
			for(auto& one_thread : threads)
				one_thread.join();

			gcheck_thresholds_mutex.lock();
			gcheck_thresholds = 1;
			gcheck_thresholds_mutex.unlock();
		}
		for (auto& dirichlistitem: dirichlist){
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				if(dirichlistitem.second->GetSingleBaseline(ichannel)==0) 
					std::cerr 
						<< "No Baseline for dirich 0x" << std::hex << dirichlistitem.first 
						<< std::dec << "'s channel " << ichannel 
						<< " found in loading-file" 
						<< std::endl;
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
		std::unordered_map<uint16_t,std::array<double,32>> thresholds;
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
			if(iss.tellg()!=-1) 
				std::cerr 
					<< "Error reading line:\n" << line 
					<< "\nRead in:" 
					<< "\ndirichaddress:0x" << dirichaddress_string 
					<< "\nchannel:" << channel 
					<< "\nbaseline:" << baseline 
					<< "\nwidth:" << width 
					<< "\nthresholdinmV:" << thresholdinmV 
					<< std::endl;
			else{
				if(set_base!=0){
					dirichptr->SetSingleBaseline(channel, baseline);
					dirichptr->SetSingleNoisewidth(channel, dirich::Thr_mVtoD(width));
				}
				if(set_thr!=0){
					if(thresholds.find(dirichaddress)==thresholds.end()){
						std::array <double,32> temp_array;
						temp_array.fill(0);
						thresholds.insert(std::make_pair(dirichaddress,temp_array));
					}
					thresholds.at(dirichaddress).at(channel) = thresholdinmV;
				}
			}
		}
		std::vector<std::thread> threads;
		if(set_thr!=0){
			gcheck_thresholds_mutex.lock();
			gcheck_thresholds = 2;
			gcheck_thresholds_mutex.unlock();

			for(auto& one_threshold : thresholds)
				threads.push_back(std::thread(
					[&one_threshold](){
						dirichlist.at(one_threshold.first)->SetThresholdsmV(
							one_threshold.second
						);
					}
				));
			for(auto& one_thread : threads)
				one_thread.join();
		
			gcheck_thresholds_mutex.lock();
			gcheck_thresholds = 1;
			gcheck_thresholds_mutex.unlock();
		}
	}
	else{
		if(set_base!=0){
			for(int ichannel=0;ichannel<32;++ichannel){
				dirichptr->SetSingleBaseline_old(
					ichannel, 
					dirichptr->GetSingleBaseline(ichannel)
				);
				dirichptr->SetSingleNoisewidth_old(
					ichannel, 
					dirichptr->GetSingleNoisewidth(ichannel)
				);
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
			if(iss.tellg()!=-1) 
				std::cerr 
					<< "Error reading line:\n" << line 
					<< "\nRead in:" 
					<< "\ndirichaddress:0x" << dirichaddress_string 
					<< "\nchannel:" << channel 
					<< "\nbaseline:" << baseline 
					<< "\nwidth:" << width 
					<< "\nthresholdinmV:" << thresholdinmV 
					<< std::endl;
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

void save_graphs(std::shared_ptr<dirich>	dirichptr, std::string filename){
	TFile* file=new TFile(Form("%s.root", filename.c_str()),"RECREATE");
	if(dirichptr==NULL){
		file->cd();
		// get_noisewidth_histo(0)->Write();
		// get_2D_rate_histo(0)->Write();
		// get_2D_rate_over_thr_histo(0)->Write();
		// get_2D_diff_over_thr_histo(0)->Write();
		// get_2D_gr_diff_over_thr_histo(0)->Write();
		// get_2D_mgr_diff_over_thr_histo(0)->Write();
		std::cout << "saving graphs" << std::endl;
		int counter=0;
		for (auto& dirichlistitem: dirichlist) {
			std::cout 
				<< "\r" 
				<< std::setw(10) << std::setprecision(2) << std::fixed 
				<< 1.*((counter+1)*100)/(dirichlist.size()) 
				<< "%" << std::flush;
			TDirectory *dirich_dir = file->mkdir(Form("dirich_0x%x",dirichlistitem.first));
			dirich_dir->cd();
			get_noisewidth_histo(dirichlistitem.second)->Write();
			get_2D_rate_histo(dirichlistitem.second)->Write();
			get_2D_rate_over_thr_histo(dirichlistitem.second)->Write();
			get_2D_diff_over_thr_histo(dirichlistitem.second)->Write();
			// get_2D_gr_diff_over_thr_histo(dirichlistitem.second)->Write();
			// get_2D_mgr_diff_over_thr_histo(dirichlistitem.second)->Write();
			// TDirectory *channels = dirich_dir->mkdir("channels");
			// channels->cd();
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				// TDirectory *ch = channels->mkdir(Form("ch:%i",ichannel));
				// ch->cd();
				dirichlistitem.second->gRateGraphs[ichannel]->Write();
				dirichlistitem.second->gRateGraphsOverBase[ichannel]->Write();
				dirichlistitem.second->gDiffRateGraphsOverBase[ichannel]->Write();
			}
			counter++;
		}
	std::cout << std::endl;
	}
	else{
		file->cd();
		TDirectory *dirich_dir = file->mkdir(Form("dirich_0x%x",dirichptr->GetBoardAddress()));
		dirich_dir->cd();
		get_noisewidth_histo(dirichptr)->Write();
		get_2D_rate_histo(dirichptr)->Write();
		get_2D_rate_over_thr_histo(dirichptr)->Write();
		get_2D_diff_over_thr_histo(dirichptr)->Write();
		// get_2D_gr_diff_over_thr_histo(dirichptr)->Write();
		// get_2D_mgr_diff_over_thr_histo(dirichptr)->Write();
		// TDirectory *channels = dirich_dir->mkdir("channels");
		// channels->cd();
		for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
			// TDirectory *ch = channels->mkdir(Form("ch:%i",ichannel));
			// ch->cd();
			dirichptr->gRateGraphs[ichannel]->Write();
			dirichptr->gRateGraphsOverBase[ichannel]->Write();
			dirichptr->gDiffRateGraphsOverBase[ichannel]->Write();
		}
	}
	file->Close();
	// gROOT->GetListOfFiles()->Remove(file); // to get a faster closing time 
}

void save()
{
	std::array<char, 64> buffer;
	buffer.fill(0);
	time_t rawtime;
	time(&rawtime);
	const auto timeinfo = localtime(&rawtime);
	strftime(buffer.data(), sizeof(buffer), "%Y%m%d_%H%M%S", timeinfo);
	std::string str = std::string(buffer.data()) + "_std_save";

	save_base(NULL,str,0);
	save_graphs(NULL,str);
}

void* scanthread_nrml(void* dirichptr) //Argument is pointer to DiRICH class instance
{
	if(((dirich*)dirichptr)->gdirich_reporting_level>=1) 
		std::cout 
			<< "Starting threshscan for Dirich at address 0x"
			<< std::hex << ((dirich*)dirichptr)->GetBoardAddress()
		<< std::endl;
	// std::cout << "DoThreshScan" << std::endl;
	((dirich*)dirichptr)->DoThreshScan();
	// std::cout << "AnalyzeBaseline" << std::endl;
	((dirich*)dirichptr)->AnalyzeBaseline();
	// std::cout << "FineScan" << std::endl;
	((dirich*)dirichptr)->DoFineThreshScan();
	// std::cout << "AnalyzeBaseline" << std::endl;
	((dirich*)dirichptr)->AnalyzeBaseline();
	// std::cout << "MakeGraphsOverBase" << std::endl;
	((dirich*)dirichptr)->MakeGraphsOverBase();
	// std::cout << "MakeDiffGraphsOverBase" << std::endl;
	((dirich*)dirichptr)->MakeDiffGraphsOverBase();
	if(((dirich*)dirichptr)->gdirich_reporting_level>=1) 
		std::cout 
			<< "Threshscan for Dirich at address 0x"
			<< std::hex << ((dirich*)dirichptr)->GetBoardAddress()
			<< " done"
		<< std::endl;
	return 0;
}

void* scanthread_over(void* dirichptr) //Argument is pointer to DiRICH class instance
{
	if(((dirich*)dirichptr)->gdirich_reporting_level>=1) 
	std::cout 
		<< "Starting threshscan_over for Dirich at address 0x"
		<< std::hex << ((dirich*)dirichptr)->GetBoardAddress()
	<< std::endl;
	((dirich*)dirichptr)->DoThreshScanOverBase();
	((dirich*)dirichptr)->MakeDiffGraphsOverBase();
	if(((dirich*)dirichptr)->gdirich_reporting_level>=1) 
		std::cout 
			<< "Threshscan_over for Dirich at address 0x"
			<< std::hex << ((dirich*)dirichptr)->GetBoardAddress()
			<< " done"
		<< std::endl;
	return 0;
}

void system_thr_scan(int type=0)
{
	if(type==0){
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 1;
		gcheck_thresholds_mutex.unlock();
	}
	else{
		gcheck_thresholds_mutex.lock();
		gcheck_thresholds = 2;
		gcheck_thresholds_mutex.unlock();		
	}
	std::vector <std::thread*> threadlist;
	// Initialize instances of dirich class for each module
	for (auto& dirichlistitem: dirichlist){
		if(dirichlistitem.second==NULL){
			std::cerr 
				<< "DiRICH " << std::hex << dirichlistitem.first 
				<< std::dec << " not initialized!" 
				<< std::endl;
			continue;
		}
		switch(type){
		case 1:
			threadlist.push_back(
				new std::thread(
					scanthread_over, 
					(void*) dirichlistitem.second.get()
				)
			);
			break;
		case 0:
		default:
			threadlist.push_back(
				new std::thread(
					scanthread_nrml, 
					(void*) dirichlistitem.second.get()
					)
				);
			break;
		}
	}
	usleep(1000);

	for(auto& thread : threadlist){
			thread->join();
			delete thread;
	}
	gcheck_thresholds_mutex.lock();
	gcheck_thresholds = 1;
	gcheck_thresholds_mutex.unlock();
	// threadlist.clear();
	// printf("System scan done ! \n");
	// switch(type){
	// case 1:
	// 	save();
	// 	break;
	// case 0:
	// default:
	// 	save();
	// 	break;
	// }

}

void* check_thresholds(){
	// return 0;
	if(self_check_threshold==true) return 0;
	int ret = 0;
	int break_counter=0;
	std::map<uint16_t,std::array<uint16_t,NRCHANNELS>> thresholds;
	std::array<uint16_t,NRCHANNELS> temp_array;
	temp_array.fill(0);
	for(auto& dirichlistitem : dirichlist){
		thresholds.insert(std::make_pair(dirichlistitem.first,temp_array));
	}
	uint32_t temp_buffer4mb[BUFFER_SIZE4mb];
	while(true){
		gcheck_thresholds_mutex.lock();
		int temp_gcheck_thresholds = gcheck_thresholds;
		gcheck_thresholds_mutex.unlock();
		switch(temp_gcheck_thresholds){
			case 0:
			default:
			return 0;
			case 1:
			std::this_thread::sleep_for(std::chrono::microseconds(10*THRESHDELAY));
			break;
			case 2:
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){
				uint32_t real_ichannel = ichannel%CHPCHAIN;
				uint32_t c[] = {
					(0x0 << 20 | real_ichannel << 24),
					0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					(uint32_t)ichannel/CHPCHAIN+1,
					0x10001
				}; 
				ret=Ttrb_register_write_mem(BROADCAST,0xd400,0,c,CHPCHAIN+2);
				if(ret<0){
					std::cerr << "Can't retreive Thresholds (1)!!!" << std::endl;
					break_counter++;
					if(break_counter>100) return 0;
					else break;
				}
				std::this_thread::sleep_for(std::chrono::microseconds(SPICOMDELAY));
				ret=Ttrb_register_read(BROADCAST,0xd412,temp_buffer4mb,BUFFER_SIZE4mb);
				if(ret<0){
					std::cerr << "Can't retreive Thresholds (2)!!!" << std::endl;
					break_counter++;
					if(break_counter>100) return 0;
					else break;
				}
				for(int i=0;i<ret;i+=2){
					try{
						thresholds.at(temp_buffer4mb[i]).at(ichannel) = uint16_t(temp_buffer4mb[i+1] & 0xffff); }
					catch(...){}
				}
			}
			auto curr_time = std::chrono::steady_clock::now();
			for(auto& threshold : thresholds){
				auto& dirich = dirichlist.at(threshold.first);

				dirich->Current_Thr_Mutex.lock();
				dirich->gCurrent_Threshold = threshold.second;
				dirich->gCurrent_Threshold_time = curr_time;
				dirich->Current_Thr_Mutex.unlock();
			}
			std::this_thread::sleep_for(std::chrono::microseconds(SPICOMDELAY));
			break;
		}
	}
	return 0;
}

void initialize_diriches(std::vector<uint16_t> diriches = {})
{
// void initialize_diriches(bool search_dirich, std::vector<int> ranges, int NrPasses, double meas_time){
	TH1::AddDirectory(0);
	gErrorIgnoreLevel = kError;
	
	int ret=0;
	ret=init_ports();
	if(ret==-1){
		std::cerr << "failed to initialize trb-net ports" << std::endl;
	}

	dirichlist.clear();

	std::unordered_map<uint16_t,std::future<dirich*>> inited_diriches;
	int dirich_counter=0;

	if(diriches.empty()){
		ret=Ttrb_read_uid(BROADCAST, buffer4mb, BUFFER_SIZE4mb);
		if(ret<4){
			std::cerr << "No TRB3 Modules found!!!" << std::endl;
			return;
		}
		for(int i=0;i<ret;i+=4){
			// if(buffer[i+3]>0x1200 && buffer[i+3]<0x1200)
			inited_diriches.insert(
				std::make_pair(
					uint16_t(buffer4mb[i+3]),
					std::async(
						std::launch::async,
						[](uint16_t uid) ->dirich* {
							return new dirich(uid);
						},
						uint16_t(buffer4mb[i+3])
					)
				)
			);
		++dirich_counter;
		}
	}
	else{
		for(auto& dirich_uid : diriches){
			// if(buffer[i+3]>0x1200 && buffer[i+3]<0x1200)
			inited_diriches.insert(
				std::make_pair(
					dirich_uid,
					std::async(
						std::launch::async,
						[](uint16_t uid) ->dirich* {
							return new dirich(uid);
						},
						dirich_uid
					)
				)
			);
		++dirich_counter;
		}
	}
	for(auto& one_dirich : inited_diriches){
		one_dirich.second.wait();
		dirich* temp_dirich_prt = one_dirich.second.get();
		if(temp_dirich_prt->WhichDirichVersion()!=3){ 
			//pls change it according to your initialization... Sure one should rather throw during init... but well I am lazy
			std::cerr 
				<< "DiRICH 0x" << std::hex << one_dirich.first 
				<< " not correclty initialized. Deleting!" 
				<< std::endl;
			delete temp_dirich_prt;
		}
		else
			dirichlist.insert(
				std::make_pair(one_dirich.first,std::shared_ptr<dirich>(temp_dirich_prt))
			);
	}
	if(diriches.empty()){
		std::cout 
			<< "Found " << std::dec << inited_diriches.size() << " different diriches\n"
			<< "Initialized " << dirichlist.size() << " out of those" << std::endl;
	}
	else{
		std::cout 
			<< "Looked for " << std::dec << diriches.size() << " different diriches\n"
			<< "Initialized " << dirichlist.size() << " out of those" << std::endl;
	}
	if(dirichlist.size()==0) exit(EXIT_FAILURE);
	// for(auto& dirichlistitem : dirichlist){
	// 	dirichlistitem.second->gdirich_reporting_level=3;	
	// }
	dirichlist.begin()->second->gdirich_reporting_level=1;
}

void setup_scan_parameters(
	std::shared_ptr<dirich>	dirichptr, 
	double gMeasureTime, 
	int gLowerEdge, 
	int gUpperEdge, 
	int gStepsize, 
	int gNrPasses
)
{
	if(dirichptr==NULL){
		for (auto& dirichlistitem: dirichlist) {
			if(dirichlistitem.second==NULL){
				std::cerr 
					<< "dirich 0x" << std::hex << dirichlistitem.first 
					<< std::dec << " not initialized" 
					<< std::endl;
				continue;
			}
			dirichlistitem.second->gMeasureTime = gMeasureTime;
			dirichlistitem.second->gLowerEdge.fill(gLowerEdge);
			dirichlistitem.second->gUpperEdge.fill(gUpperEdge);
			dirichlistitem.second->gStepsize = gStepsize;
			dirichlistitem.second->gNrPasses = gNrPasses;
		}
	}
	else{
		dirichptr->gMeasureTime = gMeasureTime;
		dirichptr->gLowerEdge.fill(gLowerEdge);
		dirichptr->gUpperEdge.fill(gUpperEdge);
		dirichptr->gStepsize = gStepsize;
		dirichptr->gNrPasses = gNrPasses;
	}
}

void setup_scan_parameters_over_thr_mV(
	std::shared_ptr<dirich>	dirichptr, 
	double gMeasureTime, 
	double gUpperEdgemV, 
	double gStepsizemV, 
	int gNrPasses
)
{
	if(dirichptr==NULL){
		for (auto& dirichlistitem: dirichlist) {
			if(dirichlistitem.second==NULL){
				std::cerr 
					<< "dirich 0x" << std::hex << dirichlistitem.first 
					<< std::dec << " not initialized" 
					<< std::endl;
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
	std::string loading_file_threshold = "";
	std::string save_file = "";

  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

	// Declare the supported options.	
	po::options_description desc("Allowed options", (unsigned)w.ws_col, 
      (unsigned)w.ws_col/3);
	desc.add_options()
		(
			"help,h",
			"produce this help message"
		)
		(
			"verbosity,v", 
			po::value<int>(), 
			"Set verbosity level"
		)
		(
			"use-dirich,u", 
			po::value<std::vector<std::string>>()->multitoken(), 
			"If this option is specified, "
			"only diriches named in its parametes will be used during "
			"this runtime"
		)
		(
			"scan-baseline,b", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Do standard baselinescan for all initialized diriches! "
				"Additionally groups of six parameters can to be given, "
				"to alter the scan parameters for one/all diriches : "
				"dirich (if 0, all diriches), "
				"measure-time (s), "
				"threshold-start-value, "
				"threshold-end-value, "
				"threshold-step-width, "
				"number of cycles (two refers to every second channel measured at a time)."
		)
		(
			"draw-scan-baseline,d", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Draw the results of the baselinescan. "
				"Dirich can be specified using this options parameter. "
				"Obviously this function fails if no scan was done!"
		)
		(
			"draw-scan-above-noise", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Draw the results of the thresholdscan above the diriches noiseband. "
				"Dirich can be specified using this options parameter. "
				"Obviously this function fails if no scan was done!"
			)
		(
			"draw-scan-above-noise-diff-gr", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Draw the results of the baselinescan above the diriches noiseband as differential plot. "
				"Dirich can be specified using this options parameter. "
				"Obviously this function fails if no scan was done!"
		)
		(
			"draw-noisewidth,w", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Draw the noisewidth. "
				"Dirich can be specified using this options parameter. "
				"Obviously this function fails if neither a scan was done nor a threshold-setting was loaded!"
		)
		(
			"measure-rate,r",
			po::value<double>()->implicit_value(10.),
				"measure the rate for all initialized diriches. "
				"Parameter is the measure time in seconds. If non is given, the standard value is 10s. "
				"Results are saved in an \"_rate.dat\"-file with corresponding date/time."
		)
		// (
		// 	"find-threshold,i", 
		// 	po::value<double>(),
		// 	"
		// 		Find the perfect threshold for the given dirich/maptm-channel-combination. 
		// 		The parameter specifies the method to find the perfect threshold:
		// 		\n0: searches for the minimum in the differentiated spectrum or for the minimal gradient
		// 		\n0<value<5: tries to find peak and sigma of the single photon distribution and 
		// 			sets the threshold to value*sigma (!!!!currently not implemented!!!)
		// 		\n5<value<100: tries to find the single photon peak and 
		// 			sets the threshold to value% of the spp-position
		// 	"
		// )
		(
			"invert,i", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
			"Inverts all following threshold-settings. "
			"Meaning the inversion of the following commands: "
			"\"--scan-above-noise\",\"--load-threshold\",\"--set-to-noise\",\"--set-threshold\""
		)
		(
			"scan-above-noise,a", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Do scan for threshold-values greater than the diriches noiseband " 
				"for all initialized diriches! "
				"Additionally groups of five parameters can to be given, "
				"to alter the scan parameters for one/all diriches : "
				"dirich (if 0, all diriches), "
				"measure-time (s), "
				"threshold-end-value (mV), "
				"threshold-step-width (mV), "
				"number of cycles (two refers to every second channel measured at a time)."
		)
		(
			"load-baseline,l", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"This option loads the baseline from the file specified in --loading-file. "
				"If no file was specified, the latest produced file is choosen. "
				"One can specify a certain dirich by using this options parameter. "
				"Be aware that this option overwrites the baseline retreived from the baselinescan"
		)
		(
			"loading-file,f", 
			po::value<std::string>(&loading_file)->default_value(""), 
			"File to load thresholds and/or baseline from"
		)
		(
			"load-threshold", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"This option loads the threshold from the file specified in --loading-file-threshold. "
				"If no file was specified, the latest produced file is choosen. "
				"One can specify a certain dirich by using this options parameter. "
				"Be aware that the thresholds are overwriten by --set-threshold"
		)
		(
			"loading-file-threshold", 
			po::value<std::string>(&loading_file_threshold)->default_value(""), 
			"File to load thresholds from. If no file is specified, the normal loading-file is used"
		)
		(
			"set-to-noise,n", 
			po::value<std::vector<std::string>>()->multitoken(), 
				"Set threshold to a certain distance in terms of noisewidth for specified diriches. "
				"First parameter specifies the dirich (0 equals all DiRICHes), "
				"the second the part of the half-noisebandwidth. "
				"If only one parameter is given. All initialized diriches are set to the part of the half-noisebandwidth"
				" specified as parameter."
				"Only positive threshold values are accepted, as the minus-sign induces errors."
				"If you want to set negative thresholds use this command with a positive threshold"
				" and invert it using \"--invert\". However a threshold of 0 is impossible!"
		)
		(
			"set-threshold,t", 
			po::value<std::vector<std::string>>()->multitoken(), 
				"Set threshold for specified diriches in mV after pre amplification. "
				"First parameter specifies the dirich (0 equals all DiRICHes), "
				"the second the threshold. "
				"If only one parameter is given. The thresholds of all initialized diriches are set to the"
				" specified as parameter."
				"Only positive threshold values are accepted, as the minus-sign induces errors."
				"If you want to set negative thresholds use this command with a positive threshold"
				" and invert it using \"--invert\". However a threshold of 0 is impossible!"
		)
		(
			"set-pattern,p", 
			po::value<std::vector<std::string>>()->multitoken(), 
				"Set pattern for specified diriches. "
				"First parameter specifies the dirich (0 equals all DiRICHes), the second the pattern. "
				"The pattern is derived by interpreting the second parameter "
				"If only one parameter is given. All initialized diriches are set to the pattern"
				" specified as parameter."
				"	as bitpattern and disabling each channel where the corresponding bit equals 0. "
				"To disable one or many diriches completely you need to put the pattern 00!. "
				"And... What you are searching for is 1431655765/2863311530"
		)
		(
			"save-file", 
			po::value<std::string>(&save_file)->default_value(""), 
			"Names the name and path of the file to be saved to. "
			"If not specified, a std. filename will be produced according to: "
			"\"(current_date_time)_std_save{.thr,.root}\""
		)
		(
			"no-autosave", 
			"Disables the autosave feature. "
			"Use this and the --save-base/-s option to only save baseline-data."
		)
		(
			"save-base,s", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Save data of all or (via parameter) specified dirich after everything else is executed! "
				"Autosaves will be still produced and saved via \"DATE_std_save{.thr,.root}\". "
				"Savefile can be set via --save-file"
		)
		(
			"save-graphs,g", 
			po::value<std::vector<std::string>>()->multitoken()->zero_tokens(), 
				"Save histograms of all or (via parameter) specified dirich after everything else is executed! "
				"Autosaves will be still produced and saved via \"DATE_std_save{.thr,.root}\". "
				"Savefile can be set via --save-file"
		)	
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << "This program can be used to perform the the following tasks for a DiRICH-FEB of the TRB family:\n"
		<< "-- measure the baseline position and noise width\n"
		<< "-- set thresholds according to basline and noisewidth\n"
		<< "-- save and load baslines and thresholds. If not disabled (--no-autosave), a save is made after each scan/change of thresholds\n"
		<< "-- measure rates w.r.t threshold above baseline\n"
		<< "-- measure rates at a certain threshold\n"
		<< "-- visualize baselines and signal-shapes above the noise band\n"
		<< "\n\nA standard command-chain used would be:\n"
		<< "\t\t./HADESthreshscan_v1 -b -t 50\n Makes a baselinescan (-b) and sets the threshold to 50mV above baseline for all initialized diriches (-t)\n\n"
		<< "\t\t./HADESthreshscan_v1 -l -n 5\n Loads the baselines from the latest produced threshold-file (*.thr) (-l) and sets the threshold to 5*noiswidth/2 for all initialized diriches (-n)\n\n"
		<< "\t\t./HADESthreshscan_v1 -b --loading-file-threshold path/to/threshold.thr --load-threshold\n Makes a baseline scan and sets thresholds (--load-threshold) according to the specified file (--loading-file-threshold)\n"
		<< std::endl;
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

	if(vm.count("use-dirich")){
		std::vector<uint16_t> temp_diriches;
		for(auto& use_diriches_options : vm["use-dirich"].as<std::vector<std::string>>()){
			temp_diriches.push_back(
				std::stoi(
					use_diriches_options.substr(use_diriches_options.find("0x")!=std::string::npos ? 
					use_diriches_options.find("0x")+2 : 0),NULL,16
				)
			);
		}
		initialize_diriches(temp_diriches);
	}
	else{
		initialize_diriches();
	}

	std::cout << "All set and done" << std::endl;
	std::thread* threshold_checker= new std::thread(check_thresholds);
	gcheck_thresholds_mutex.lock();
	gcheck_thresholds = 1;
	gcheck_thresholds_mutex.unlock();

	if(vm.count("verbosity")){
		std::cout << "Setting Verbosity to level " << vm["verbosity"].as<int>() << std::endl;
		for(auto& dirich : dirichlist){
			dirich.second->gdirich_reporting_level=vm["verbosity"].as<int>();
		}
	}

	if(vm.count("invert")){
		if(vm["invert"].as<std::vector<std::string>>().size() < 1 
			|| vm["invert"].as<std::vector<std::string>>().at(0)=="0") {
			std::cout << "Inverting all diriches"<< std::endl;
			for(auto& dirich : dirichlist){
				auto curr_orientation = dirich.second->GetOrientation();
				for(auto& one_curr_orientation : curr_orientation) {
					one_curr_orientation*=-1;
				}
				dirich.second->SetOrientation(curr_orientation);
			}
		}
		else{
			for(auto& invert_opt : vm["invert"].as<std::vector<std::string>>()){
				uint16_t dirichnr = std::stoi(
					invert_opt.substr(invert_opt.find("0x")!=std::string::npos ?
					invert_opt.find("0x")+2 : 0),NULL,16
				);
				std::cout << "Inverting 0x" << std::hex << dirichnr << std::dec << std::endl;
				try{
					auto curr_orientation = dirichlist.at(dirichnr)->GetOrientation();
					for(auto& one_curr_orientation : curr_orientation) {
						one_curr_orientation*=-1;
					}
					dirichlist.at(dirichnr)->SetOrientation(curr_orientation);
				}
				catch(...) {
					std::cout << "Dirich 0x" << std::hex << dirichnr << std::dec 
						<< "does not exist" << std::endl;
				}
			}
		}
	}

	if(vm.count("scan-baseline")){
		if(vm["scan-baseline"].as<std::vector<std::string>>().size() < 6){
			std::cout 
				<< "no or less than six arguments were provided for option --scan-baseline:"
				"\nrunning scan with std. parameters" 
				<< std::endl;
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
				std::cout 
					<< std::stoi(
						one_scan_base_opt.at(0).substr(one_scan_base_opt.at(0).find("0x")!=std::string::npos ? 
						one_scan_base_opt.at(0).find("0x")+2 : 0),NULL,16
					) 
					<< "\t" << std::stod(one_scan_base_opt.at(1)) 
					<< "\t" << std::stoi(one_scan_base_opt.at(2)) 
					<< "\t" << std::stoi(one_scan_base_opt.at(3)) 
					<< "\t" << std::stoi(one_scan_base_opt.at(4)) 
					<< "\t" << std::stoi(one_scan_base_opt.at(5)) 
					<< std::endl;
				setup_scan_parameters(
					one_scan_base_opt.at(0) == "0" ? 
					0 : dirichlist.at(
							std::stoi(
								one_scan_base_opt.at(0).substr(one_scan_base_opt.at(0).find("0x")!=std::string::npos ?
						 		one_scan_base_opt.at(0).find("0x")+2 : 0),NULL,16
					 		)
						)
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
			if(vm["draw-scan-baseline"].as<std::vector<std::string>>().empty()){
				draw_histo(get_2D_rate_histo(NULL),NULL);
			}
			for(auto& draw_scan_baseline_options : vm["draw-scan-baseline"].as<std::vector<std::string>>()){
				// std::cout << "input for draw-scan-baseline: " << draw_scan_baseline_options << std::endl;
				if(draw_scan_baseline_options=="0") 
					draw_histo(get_2D_rate_histo(NULL),NULL);
				else 
					draw_histo(
						get_2D_rate_histo(
							dirichlist.at(
								std::stoi(
									draw_scan_baseline_options.substr(draw_scan_baseline_options.find("0x")!=std::string::npos ? 
										draw_scan_baseline_options.find("0x")+2 : 0
									),
								NULL,
								16
								)
							)
						),
						NULL
					);
			}
		}
		if(vm.count("draw-noisewidth")){
			if(vm["draw-noisewidth"].as<std::vector<std::string>>().empty()){
				draw_histo(get_noisewidth_histo(NULL),NULL);
			}
			for(auto& draw_noisewidth_options : vm["draw-noisewidth"].as<std::vector<std::string>>()){
				if(draw_noisewidth_options=="0") 
					draw_histo(get_noisewidth_histo(NULL),NULL);
				else 
					draw_histo(
						get_noisewidth_histo(
							dirichlist.at(
								std::stoi(
									draw_noisewidth_options.substr(draw_noisewidth_options.find("0x")!=std::string::npos ? 
										draw_noisewidth_options.find("0x")+2 : 0
									),
									NULL,
									16
								)
							)
						),
					NULL
				);
			}
		}
	}

	if(vm.count("load-baseline")){
		if(loading_file==""){
			std::cout << "no loading-file found!\n! Loading no baselines !" << std::endl;
		}
		else{
			std::cout << "loading_file: " << loading_file << std::endl;
			if(vm["load-baseline"].as<std::vector<std::string>>().empty()){ 
				load_base(nullptr, loading_file, 0, 1, 0);
			}
			else for(auto& load_baseline_opt : vm["load-baseline"].as<std::vector<std::string>>()){
				if(load_baseline_opt=="0"){
					load_base(nullptr, loading_file, 0, 1, 0);
				}
				else{
					load_base(
						dirichlist.at(
							std::stoi(
								load_baseline_opt.substr(load_baseline_opt.find("0x")!=std::string::npos ? 
									load_baseline_opt.find("0x")+2 : 0),
								NULL,
								16)
							), 
						loading_file, 
						0, 1, 0);
				}
			}
			if(vm.count("draw-noisewidth")){
				if(vm["draw-noisewidth"].as<std::vector<std::string>>().empty()){
					draw_histo(get_noisewidth_histo(NULL),NULL);
				}
				for(auto& draw_noisewidth_options : vm["draw-noisewidth"].as<std::vector<std::string>>()){
					if(draw_noisewidth_options=="0") 
						draw_histo(get_noisewidth_histo(NULL),NULL);
					else 
						draw_histo(
							get_noisewidth_histo(
								dirichlist.at(
									std::stoi(draw_noisewidth_options.substr(draw_noisewidth_options.find("0x")!=std::string::npos ? 
										draw_noisewidth_options.find("0x")+2 : 0),
									NULL,
									16
									)
								)
							),
							NULL
						);
				}
			}
		}
	}

	if(vm.count("scan-above-noise")){
		if(vm["scan-above-noise"].as<std::vector<std::string>>().size() < 5){
			std::cout 
			<< "no or less than five arguments were provided for option --scan-above-noise:\nrunning scan with std. parameters" 
			<< std::endl;
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
				std::cout 
				<< std::stoi(one_scan_above_noise_opt.at(0).substr(
					one_scan_above_noise_opt.at(0).find("0x")!=std::string::npos ? 
						one_scan_above_noise_opt.at(0).find("0x")+2 : 0),
					NULL,
					16
				) << "\t"
				<< std::stod(one_scan_above_noise_opt.at(1)) << "\t"
				<< std::stod(one_scan_above_noise_opt.at(2)) << "\t"
				<< std::stod(one_scan_above_noise_opt.at(3)) << "\t"
				<< std::stoi(one_scan_above_noise_opt.at(4)) 
					<< std::endl;
				setup_scan_parameters_over_thr_mV(
					one_scan_above_noise_opt.at(0) == "0" ? 
						0 : dirichlist.at(
							std::stoi(one_scan_above_noise_opt.at(0).substr(one_scan_above_noise_opt.at(0).find("0x")!=std::string::npos ? 
								one_scan_above_noise_opt.at(0).find("0x")+2 : 0),
							NULL,
							16
						)
					)
					, std::stod(one_scan_above_noise_opt.at(1))
					, std::stod(one_scan_above_noise_opt.at(2))
					, std::stod(one_scan_above_noise_opt.at(3))
					, std::stoi(one_scan_above_noise_opt.at(4))
				);
			}
		}
		system_thr_scan(1);
	}

	if(vm.count("load-threshold")){
		if(loading_file_threshold==""){
			loading_file_threshold=loading_file;
		}
		if(loading_file_threshold==""){
			std::cout << "no loading-file found!\n! Loading no thresholds !" << std::endl;
		}
		else{
			if(vm["load-threshold"].as<std::vector<std::string>>().empty()){
				std::cout << "loading_file_threshold: " << loading_file_threshold << std::endl;
				load_base(NULL, loading_file_threshold,0, 0, 1);
			}
			for(auto& load_threshold_options : vm["load-threshold"].as<std::vector<std::string>>()){
				std::cout << "loading_file_threshold: " << loading_file_threshold << std::endl;
				if(load_threshold_options=="0"){
					load_base(NULL, loading_file_threshold,0, 0, 1);
				}
				else{
					load_base(
						dirichlist.at(
								std::stoi(
								load_threshold_options.substr(
									load_threshold_options.find("0x")!=std::string::npos ? load_threshold_options.find("0x")+2 : 0
									),
								NULL,
								16
							)
						), 
						loading_file_threshold, 
						0, 0, 1);
				}
			}
		}
	}

	if(vm.count("scan-above-noise")){
		if(vm.count("draw-scan-above-noise")){
			if(vm["draw-scan-above-noise"].as<std::vector<std::string>>().empty()){
				draw_histo(get_2D_rate_over_thr_histo(NULL),NULL);
			}
			for(auto& draw_scan_above_noise_options : vm["draw-scan-above-noise"].as<std::vector<std::string>>()){
				if(draw_scan_above_noise_options=="0") 
					draw_histo(get_2D_rate_over_thr_histo(NULL),NULL);
				else 
					draw_histo(
						get_2D_rate_over_thr_histo(
							dirichlist.at(
								std::stoi(
									draw_scan_above_noise_options.substr(draw_scan_above_noise_options.find("0x")!=std::string::npos ? 
										draw_scan_above_noise_options.find("0x")+2 : 0),
									NULL,
									16
								)
							)
						),
						NULL
					);
			}
		}
		if(vm.count("draw-scan-above-noise-diff-gr")){
			if(vm["draw-scan-above-noise-diff-gr"].as<std::vector<std::string>>().empty()){
				draw_multigraph(get_2D_mgr_diff_over_thr_histo(NULL),NULL);
			}
			for(auto& draw_scan_above_noise_diff_options : vm["draw-scan-above-noise-diff-gr"].as<std::vector<std::string>>()){
				if(draw_scan_above_noise_diff_options=="0") 
					draw_multigraph(get_2D_mgr_diff_over_thr_histo(NULL),NULL);
				else 
					draw_multigraph(
						get_2D_mgr_diff_over_thr_histo(
							dirichlist.at(
								std::stoi(
									draw_scan_above_noise_diff_options.substr(draw_scan_above_noise_diff_options.find("0x")!=std::string::npos ? 
										draw_scan_above_noise_diff_options.find("0x")+2 : 0),
									NULL,
									16
								)
							)
						),
						NULL
					);
			}
		}
	}

	if(vm.count("set-to-noise")){
		if(vm["set-to-noise"].as<std::vector<std::string>>().size() < 2){
			std::cout << "Setting Threshold of all diriches to: " 
				<< std::stod(vm["set-to-noise"].as<std::vector<std::string>>().at(0))
				<< " times the half-noisebandwidth" << std::endl;
			set_thresholds_to_noise(nullptr, std::stod(vm["set-to-noise"].as<std::vector<std::string>>().at(0)));
		}
		else{
			std::vector<std::vector<std::string>> each_set_threshold_noise_opt;
			std::vector<std::string> temp_vec;
			for(auto& set_threshold_noise_opt : vm["set-to-noise"].as<std::vector<std::string>>()){
				if(set_threshold_noise_opt.find("0x")!=std::string::npos || set_threshold_noise_opt=="0"){
					if(temp_vec.size()==2){
						each_set_threshold_noise_opt.push_back(temp_vec);
					}
					temp_vec.clear();
				}
				temp_vec.push_back(set_threshold_noise_opt);
			}
			if(temp_vec.size()==2){
				each_set_threshold_noise_opt.push_back(temp_vec);
			}	 
			for(auto& one_set_threshold_noise_opt : each_set_threshold_noise_opt){
				std::cout << "Setting Threshold of: "
									<< std::stoi(
										one_set_threshold_noise_opt.at(0).substr(one_set_threshold_noise_opt.at(0).find("0x")!=std::string::npos ? 
											one_set_threshold_noise_opt.at(0).find("0x")+2 : 0),
										NULL,
										16
									) << "\t"
									<< " to: "
									<< std::stod(one_set_threshold_noise_opt.at(1))
									<< " times the noisebandwidth" << std::endl;
				set_thresholds_to_noise(
					one_set_threshold_noise_opt.at(0) == "0" ? 
					0 : dirichlist.at(
						std::stoi(
							one_set_threshold_noise_opt.at(0).substr(one_set_threshold_noise_opt.at(0).find("0x")!=std::string::npos ? 
								one_set_threshold_noise_opt.at(0).find("0x")+2 : 0),
							NULL,
							16
						)
					), 
					std::stod(one_set_threshold_noise_opt.at(1))
				);
			}
		}
	}

	if(vm.count("set-threshold")){
		if(vm["set-threshold"].as<std::vector<std::string>>().size() < 2){
			std::cout << "Setting Threshold of all diriches to: "
				<< std::stod(vm["set-threshold"].as<std::vector<std::string>>().at(0))
				<< " mV" << std::endl;
			set_thresholds(0, std::stoi(vm["set-threshold"].as<std::vector<std::string>>().at(0)));
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
				std::cout << "Setting Threshold of: "
									<< std::stoi(
											one_set_threshold_opt.at(0).substr(one_set_threshold_opt.at(0).find("0x")!=std::string::npos ? 
												one_set_threshold_opt.at(0).find("0x")+2 : 0),
											NULL,
											16
										) 
									<< "\t" << " to: "
									<< std::stod(one_set_threshold_opt.at(1))
									<< std::endl;
				set_thresholds(
					one_set_threshold_opt.at(0) == "0" ? 
						0 : dirichlist.at(
							std::stoi(
								one_set_threshold_opt.at(0).substr(one_set_threshold_opt.at(0).find("0x")!=std::string::npos ? 
									one_set_threshold_opt.at(0).find("0x")+2 : 0),
								NULL,
								16
							)
						), 
					std::stoi(one_set_threshold_opt.at(1))
				);
			}
		}
	}

	if(vm.count("set-pattern")){
		if(vm["set-pattern"].as<std::vector<std::string>>().size() < 2){
			std::cout << "Setting all diriches to a pattern of: ";
			for(int ichannel=0;ichannel<NRCHANNELS;++ichannel){ 
				std::cout << (std::stol(vm["set-pattern"].as<std::vector<std::string>>().at(0)) >> ichannel) % 2;
			}
			std::cout << std::endl;
			set_pattern(0, std::stol(vm["set-pattern"].as<std::vector<std::string>>().at(0)));
		}
		else{
			std::vector<std::vector<std::string>> each_set_pattern_opt;
			std::vector<std::string> temp_vec;
			for(auto& set_pattern_opt : vm["set-pattern"].as<std::vector<std::string>>()){
				if(set_pattern_opt.find("0x")!=std::string::npos || set_pattern_opt=="0"){
					if(temp_vec.size()==2){
						each_set_pattern_opt.push_back(temp_vec);
					}
					temp_vec.clear();
				}
				temp_vec.push_back(set_pattern_opt);
			}
			if(temp_vec.size()==2){
				each_set_pattern_opt.push_back(temp_vec);
			}	 
			for(auto& one_set_pattern_opt : each_set_pattern_opt){
				std::cout << "Setting "
									<< std::stoi(
										one_set_pattern_opt.at(0).substr(one_set_pattern_opt.at(0).find("0x")!=std::string::npos ? 
										one_set_pattern_opt.at(0).find("0x")+2 : 0),NULL,16) 
									<< "\t"	<< " to a pattern of: ";
									for(int ichannel=0;ichannel<NRCHANNELS;++ichannel) 
										std::cout << (std::stol(one_set_pattern_opt.at(1)) >> ichannel) % 2;
									std::cout << std::endl;
				set_pattern(
					one_set_pattern_opt.at(0) == "0" ? 
					0 : dirichlist.at(
						std::stoi(
							one_set_pattern_opt.at(0).substr(one_set_pattern_opt.at(0).find("0x")!=std::string::npos ? 
								one_set_pattern_opt.at(0).find("0x")+2 : 0),
							NULL,
							16
						)
					), 
					std::stol(one_set_pattern_opt.at(1))
				);
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

	if(!vm.count("no-autosave") 
			&& (vm.count("scan-baseline") || vm.count("scan-above-noise") 
					|| vm.count("set-to-noise") || vm.count("set-threshold")
					|| vm.count("set-pattern") )
	){
		save_base(NULL,save_file,1);
		save_graphs(NULL,save_file);
	}

	if(vm.count("save-base")){
		if(vm["save-base"].as<std::vector<std::string>>().empty()){
			save_base(NULL,save_file,1);
		}
		for(auto& save_options : vm["save-base"].as<std::vector<std::string>>()){
			if(save_options=="0"){
				save_base(NULL,save_file,1);
			}
			else{
				save_base(
					dirichlist.at(
						std::stoi(
							save_options.substr(save_options.find("0x")!=std::string::npos ? 
								save_options.find("0x")+2 : 0),
							NULL,
							16
						)
					),
					save_file,
					1
				);
			}
		}
	}
	if(vm.count("save-graphs")){
		if(vm["save-graphs"].as<std::vector<std::string>>().empty()){
			save_graphs(NULL,save_file);
		}
		for(auto& save_options : vm["save-graphs"].as<std::vector<std::string>>()){
			if(save_options=="0"){
				save_graphs(NULL,save_file);
			}
			else{
				save_graphs(
					dirichlist.at(
						std::stoi(
							save_options.substr(save_options.find("0x")!=std::string::npos ? 
								save_options.find("0x")+2 : 0),
							NULL,
							16
						)
					),
					save_file
				);
			}
		}
	}

	std::string rate_file = save_file + "_rate.dat";
	if(vm.count("measure-rate")){
		double rate_measure_time = vm["measure-rate"].as<double>();
		std::cout << "Measuring rate of all diriches over " 
				<< rate_measure_time << " seconds" << std::endl;
		measure_rate(
			0,
			rate_file,
			rate_measure_time
		);
	}

	std::string str = save_file + "_all_canvases.pdf";
	uint counter=0;
	for(auto& canvases : canvasvector){
		if(counter==0 && canvasvector.size()>1) canvases->Print(Form("%s(",str.c_str()));
		else if(counter==canvasvector.size()-1) canvases->Print(Form("%s)",str.c_str()));
		else canvases->Print(str.c_str());
		counter++;
	}

	usleep(100000);
	gcheck_thresholds_mutex.lock();
	gcheck_thresholds = 0;
	gcheck_thresholds_mutex.unlock();	
	threshold_checker->join();

	return 0;
}
