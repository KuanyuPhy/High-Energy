/**
 * @file	exercise_04.cc
 * @date	Aug 15, 2011
 * @author	mprim
 * @brief	brief description
 *
 * long description
 */

// ROOT general
#include "TROOT.h"
#include "TMath.h"
// ROOT data
#include "TFile.h"
#include "TChain.h"
#include "TLeaf.h"
// ROOT plot
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "THStack.h"
// RooFit
#include "RooAbsPdf.h"
#include "RooAddModel.h"
#include "RooAddPdf.h"
#include "RooArgusBG.h"
#include "RooBifurGauss.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooGlobalFunc.h"
#include "RooHist.h"
#include "RooMCStudy.h"
#include "RooMsgService.h"
#include "RooNumIntConfig.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
// program options
#include "options.h"
// Utilities
#include "utility.h"
// std library
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace RooFit;


/**
 * @brief Sets some reasonable global ROOT Style settings
 */
inline void SetRootStyle()
{
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.04, "xy");
	gStyle->SetLabelOffset(0.006, "y");
	gStyle->SetTitleSize(0.04, "xy");
	gStyle->SetTitleOffset(1.0, "x");
	gStyle->SetTitleOffset(1.3, "y");
	gStyle->SetNdivisions(505, "x");

	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.10);
	gStyle->SetPadTopMargin(0.05);

	gStyle->SetFillColor(0);
	gStyle->SetMarkerSize(0.8);
	gStyle->SetLineColor(kBlack);
	gStyle->SetLineWidth(1);

	gStyle->SetLegendBorderSize(0);
}
void ProduceToyMCTuples(const RooAbsPdf &pdf_signal, const RooAbsPdf &pdf_background, const RooArgSet &argset, Options *opt,
		unsigned int nsig_events = 5000, unsigned int nbkg_events =5000 ) {

	RooDataSet *toy_mcdata_sig = pdf_signal.generate(argset,nsig_events);
	RooDataSet *toy_mcdata_bkg = pdf_background.generate(argset,nbkg_events);

	TFile *output_file = new TFile((opt->GetOutputDir()+"toymc_basic.root").c_str(),"RECREATE");
	TTree *output_tree = new TTree("h1","h1");
	float Mb_c,de;
	int mcinfo;
	output_tree->Branch("Mb_c",&Mb_c,"Mb_c/F");
	output_tree->Branch("de",&de,"de/F");
	output_tree->Branch("mcinfo",&mcinfo,"mcinfo/I");

	for(unsigned int i = 0; i < nsig_events; ++i) {
		const RooArgSet* values = toy_mcdata_sig->get(i);
		Mb_c = values->getRealValue("Mb_c");
		de = values->getRealValue("de");
		mcinfo = 1;
		output_tree->Fill();
	}

	for(unsigned int i = 0; i < nbkg_events; ++i) {
		const RooArgSet* values = toy_mcdata_bkg->get(i);
		Mb_c = values->getRealValue("Mb_c");
		de = values->getRealValue("de");
		mcinfo = 0;
		output_tree->Fill();
	}

	output_tree->Write();
	delete output_tree;
	output_file->Close();
	delete output_file;
}

/**
 * @brief Main method
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[])
{
	/**
	 * Parse program options and setup some global settings
	 */
	// parse program options
	Options *opt = new Options();
	if (!(opt->ParseOptions(argc, argv)))
	{
		std::exit(EXIT_FAILURE);
	}
	else
	{
		opt->PrintOptions();
	}
	// create output directory
	util::execute_command((util::to_string("mkdir -p ") + opt->GetOutputDir()).c_str());
	// set the global ROOT style
	SetRootStyle();
	// create a chain with all input files
	TChain *chain = new TChain(opt->GetTreename().c_str());
	for (unsigned int i = 0; i < opt->GetFilenames().size(); ++i)
	{
		chain->Add(opt->GetFilenames()[i].c_str());
	}
	// create an empty canvas to draw all the plots on
	TCanvas *canvas = new TCanvas("plot", "plot", 1280, 1024);

	//----------資料是如何讀近來的呢？？？
	/**
	 * Part 1: Read input data and define variables
	 */
	// define the different variables
	RooRealVar Mb_c("Mb_c", "M_{bc}", 5.24, 5.29, "GeV");
	RooRealVar de("de", "#Delta E", -0.2, 0.2, "GeV");
	RooRealVar contsupp("contsupp", "Continuum Suppression", -15, 15);
	RooCategory mcinfo("mcinfo", "MC Info"); //no use
	mcinfo.defineType("background", 0);
	mcinfo.defineType("signal", 1);

	// fill the RooDataSet
	RooArgSet *varset = new RooArgSet(Mb_c, de); //這裡的最要和root裡的依樣 , Convert n-tuple  to RooDataSet
	RooDataSet data("data", "data", chain, *varset);
	// Print some information of the RooDataSet (e.g. number of events)
	data.Print();

	// define full range of input data
	Mb_c.setRange("full_range", 5.24, 5.29);
	de.setRange
	("full_range", -0.2, 0.2);
	contsupp.setRange("full_range", -15, 15);

	/**i
	 * Part 2: Define 1-dimensional PDFs
	 */
	// ----- 1-dim Mb_c PDF -----
	// Gaussian 紅色
	//RooRealVar mbc_signal_mean("mbc_signal_mean", "Mb_c signal mean", 5.289, 5.27, 5.29);			//5.28中間值 下限值 上限值    minuit演算法
	//RooRealVar mbc_signal_width("mbc_signal_width", "Mb_c signal width",0.5,0.0,1.0);  //---------這些寬度之後要自己設嗎？為什麼title要寫兩次
	//RooGaussian mbc_signal("mbc_signal", "Mb_c signal", Mb_c, mbc_signal_mean, mbc_signal_width); //後半部再把前面宣告的變數抓近來去fit gaussian
	//RooGaussian mbc("mbc_signal", "Mb_c signal", Mb_c, mbc_signal_mean, mbc_signal_width);//tttttttttttttttttttttttt
	
        //one Gasussian
       //RooRealVar mga1_mean("mga1_mean", "mga1_mean", 5.28504, 5.26999,5.29);
       //RooRealVar mga1_width("mga1_width", "mga1_width",0.0025999,0,0.003599);
        //RooGaussian mga1("mga1", "ga1",Mb_c, mga1_mean, mga1_width);
        //RooProdPdf mbc_pdf("mmbc_signal_ga","mmbc signal_ga",RooArgList(mga1,mbc));

     	//RooRealVar mbc_signal_fraction("mbc_signal_fraction", "Mb_c signal fraction", 0.5, 0.0, 1.0);
	//RooAddPdf mbc_pdf("mbc_pdf","Mb_c pdf", mbc,mbc_signal_fraction);//testtttt
	
        //Christal Ball
        RooRealVar mbc_signal_mean("mbc_signal_mean", "mbc signal mean",5.28,5.265,5.29);
        RooRealVar mbc_sigma_delta("mbc_sigma_delta","mbc sigma_delta",0.1,0.0,0.2);
        RooRealVar mbc_alpha("mbc_alpha","alpha",1,0.0,2.0);
        RooRealVar mbc_n("n","mbc_n",5.0,0.0,1000.0);
        RooCBShape mbc_pdf("mbc_signal", "mbc signal",Mb_c, mbc_signal_mean,mbc_sigma_delta,mbc_alpha,mbc_n);

	// Argus 藍色
	/*
	RooRealVar mbc_background_m0("mbc_background_m0", "Mb_c background m0", 5.289);
	RooRealVar mbc_background_c("mbc_background_c", "Mb_c background c", -60.0, -200.0, -1.0);
	RooRealVar mbc_background_p("mbc_background_p", "Mb_c background p", 0.5);
	RooArgusBG mbc_background("mbc_background", "Mb_c background", Mb_c, mbc_background_m0, mbc_background_c, mbc_background_p);
	*/
	// Gaussian +++ Argus 兩個加起來
	//RooRealVar mbc_signal_fraction("mbc_signal_fraction", "Mb_c signal fraction", 0.5, 0.0, 1.0); //比例  這個gaussain 和 argus的比例
	//RooAddPdf mbc_pdf("mbc_pdf", "Mb_c pdf", mbc_signal, mbc_background, mbc_signal_fraction);	//把兩個pdf加起來
	//RooAddPdf mbc_pdf("mbc_pdf","Mb_c pdf", mbc_signal);//testtttt

	//////////////////////////////////
	// ----- 1-dim deltaE PDF -----
	//////////////////////////////////
   	//Christal Ball
   	RooRealVar deltae_signal_mean("deltae_signal_mean", "de signal mean",0.18, -0.2, 0.25);
	RooRealVar sigma_delta("sigma_delta","sigma_delta",0.04,0,3);//0.04,0,0.3
	RooRealVar alpha("alpha","alpha",0.6,0.0,5);
	RooRealVar n("n","n",4.0,0.0,500);
    	RooCBShape deltae_pdf("deltae_signal", "de signal",de, deltae_signal_mean,sigma_delta,alpha,n);	
	//RooProdPdf deltae_signal("deltae_signal","de signal", RooArgList(gaussx,deltae_signal1)) ;
	
	//Three Gasussian
	//RooRealVar ga1_mean("ga1_mean", "ga1_mean", 0.1,-0.2,0.1);			
	//RooRealVar ga1_width("ga1_width", "ga1_width",0.0001,0.0,0.06); 

	//RooRealVar ga2_mean("ga2_mean", "ga2_mean", 0.0, -0.2, 0.3);			
	//RooRealVar ga2_width("ga2_width", "ga2_width",0.05, 0.0, 1.0);  

	//RooRealVar ga3_mean("ga3_mean", "ga3_mean", 0.1, -0.2, -.2);			
	//RooRealVar ga3_width("ga3_width", "ga3_width", 0.03, 0.0, 0.06);   

	//RooGaussian ga1("ga1", "ga1", de, ga1_mean, ga1_width); 
    	//RooGaussian ga2("ga2", "ga2", de, ga2_mean, ga2_width); 
    	//RooGaussian ga3("ga3", "ga3", de, ga3_mean, ga3_width); 
	//RooProdPdf deltae_pdf("de_signal_pdf", "de signal_pdf",RooArgList(ga2,deltae)); 
	//RooGaussian deltae_background("deltae_background", "deltae_background",de,ga1_mean,ga1_width); 

	//Gaussian +++ Christal Ball
	//RooRealVar deltae_signal_fraction("deltae_signal_fraction", "de signal fraction", 0.9,0.0,1.0);
	//RooAddPdf deltae_pdf("deltae_pdf", "de pdf", deltae_signal, deltae_background, deltae_signal_fraction);
	
	//////
	//Bifur-Gaussian
	//////
        //RooRealVar dde_signal_mean("dde_signal_mean","dde signal mean",0.1,-0.2,0.3);
        //RooRealVar dde_signal_width_left("dde_signal_width_left","dde signal width left",0.01,0.0,10.0);
        //RooRealVar dde_signal_width_right("dde_signal_width_right","dde signal width right",0.01,0.0,3.0);
        //RooBifurGauss dde_signal("dde_signal","dde signal",de,dde_signal_mean,dde_signal_width_left,dde_signal_width_right);
	//RooProdPdf deltae_pdf("de_signal_pdf", "de signal_pdf",RooArgList(dde_signal,deltae));
	/**
	 * Part 3: Fit 1-dimensional PDFs to data
	 */
	if (opt->GetSkipFit() == false)
	{ //-------這裡在幹嘛呢？
		// fit to Mb_c  執行fit
		RooFitResult *mbc_result = mbc_pdf.fitTo(data, NumCPU(opt->GetNumcpu()), Timer(true), Save(true));
		// fit to de
		RooFitResult *deltae_result = deltae_pdf.fitTo(data, NumCPU(opt->GetNumcpu()), Timer(true), Save(true));
		
	

		// in this example we dont't use the RooFitResults later, so we clean them up directly
		delete mbc_result;
		delete deltae_result;
		
	}

	/**
	 * Part 4: Visualize the 1-dimensional fit results
	 */
	// visualize the Mb_c fit result
	canvas->Clear();          //---這些BINS是在幹嘛呢？？
	RooPlot *plot_mbc = Mb_c.frame(Bins(opt->GetBins()), Name("plot_mbc"), Title("plot_mbc"));
	data.plotOn(plot_mbc, Name("data_hist"), CutRange("full_range"), MarkerColor(kBlack));
	mbc_pdf.plotOn(plot_mbc, Name("mbc_fit_result"), ProjectionRange("full_range"), Components("mbc_pdf"), LineColor(kBlue), LineStyle(kSolid));
	mbc_pdf.plotOn(plot_mbc, Name("mbc_sig_fit_result"), ProjectionRange("full_range"), Components("mbc_signal"), LineColor(kRed), LineStyle(kDashDotted));
	mbc_pdf.plotOn(plot_mbc, Name("mbc_bkg_fit_result"), ProjectionRange("full_range"), Components("mbc_background"), LineColor(kBlue), LineStyle(kDashed));
	plot_mbc->Draw("");
	canvas->SaveAs((opt->GetOutputDir() + "1D_plot_mbc.eps").c_str());
	delete plot_mbc;

	// visualize the de fit result
	canvas->Clear();
	RooPlot *plot_deltae = de.frame(Bins(opt->GetBins()), Name("plot_deltae"), Title("plot_deltae"));
	data.plotOn(plot_deltae, Name("data_hist"), CutRange("full_range"), MarkerColor(kBlack));
	deltae_pdf.plotOn(plot_deltae, Name("deltae_fit_result"), ProjectionRange("full_range"), Components("deltae_pdf"), LineColor(kBlue), LineStyle(kSolid));
	//------
	deltae_pdf.plotOn(plot_deltae, Name("deltae_sig_fit_result"), ProjectionRange("full_range"), Components("deltae_signal"), LineColor(kRed), LineStyle(kDashDotted));
	deltae_pdf.plotOn(plot_deltae, Name("deltae_bkg_fit_result"), ProjectionRange("full_range"), Components("deltae_background"), LineColor(kBlue), LineStyle(kDashed));
	
	//put information on
	deltae_pdf.paramOn(plot_deltae,Layout(0.6));
	//deltae_pdf.paramOn(plot_deltae, Layout(0.6), Format("NEU", AutoPrecision(2)));
	//data.statOn(plot_deltae);
	// TODO: Add signal and background component to the deltaE plot
	plot_deltae->Draw("");
	canvas->SaveAs((opt->GetOutputDir() + "1D_plot_deltae.eps").c_str());
	//RooAbsPdf::paramOn(frame)

	//float xxx=deltae_signal_mean.getVal();
	//std::cout<<"Mean"<<xxx;

	delete plot_deltae;

	
	/**
	 * Part 5: Define multi-dimensionale PDFs
	 */
	// construct final signal and background PDF as RooProdPdf
	//RooProdPdf final_signal("final_signal","final signal pdf",RooArgList(mbc_signal,deltae_signal));
	RooProdPdf final_signal("final_signal","final signal pdf",RooArgList(mbc_pdf,deltae_pdf));//////////////

	//RooProdPdf final_background("final_background","final background pdf",RooArgList(mbc_background,deltae_background));tttttttttttttttt

	// variables for extended likelihood fit
	unsigned long nentries = (long)chain->GetEntries();
	// set reasonable starting values for signal and background events (for extened maximum likelihood fit)
	RooRealVar nsig("nsig","number signal events",nentries*0.6,0,nentries*1.2);
	//RooRealVar nbkg("nbkg","number background events",nentries*0.4,0,nentries*1.2);ttttttttttttttttttttttttttttttt

	// final PDF as RooAddPdf of signal+background
	//RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(final_signal,final_background),RooArgList(nsig,nbkg));ttttttttttttttttttttttttt
	RooAddPdf final_pdf("final_pdf","final pdf",final_signal,nsig);//tttttttttttttttttttttt
	/**
	 * Part 6: Fit multi-dimensional PDFs
	 */
	if(opt->GetSkipFit() == false) {
		// fit the multi-dimensional pdf to the data distribution
		RooFitResult* final_result = final_pdf.fitTo(data, NumCPU(opt->GetNumcpu()), Timer(true), Extended(true), Save(true));

		// in this example we dont't use the RooFitResult later, so we clean it up directly
		delete final_result;
	}
	/**
	 * Part 7: Visualize the multi-dimensional fit results
	 */
	// define signal box ranges for each variable, with cuts on other variables applied
	Mb_c.setRange("signal_box_mbc",5.2,5.29);
	de.setRange("signal_box_mbc",-0.1,0.1);
	//contsupp.setRange("signal_box_mbc",0,15);

	Mb_c.setRange("signal_box_deltae",5.27,5.29);
	de.setRange("signal_box_deltae",-0.2,0.2);
	//contsupp.setRange("signal_box_deltae",0,15);

	//Mb_c.setRange("signal_box_contsupp",5.27,5.29);
	//de.setRange("signal_box_contsupp",-0.1,0.1);
	//contsupp.setRange("signal_box_contsupp",-15,15);

	// visualize the mbc fit result by cutting on other dimensions via the ProjectionRange argument
	canvas->Clear();
	RooPlot* plot_mbc_projection = Mb_c.frame(Bins(opt->GetBins()), Name("plot_mbc_projection"), Title("plot_mbc_projection"));
	data.plotOn(plot_mbc_projection, Name("data_hist"), CutRange("signal_box_mbc"), MarkerColor(kBlack));
	final_pdf.plotOn(plot_mbc_projection, Name("mbc_fit_result"), ProjectionRange("signal_box_mbc"), Components("final_pdf"), LineColor(kBlue), LineStyle(kSolid));
	final_pdf.plotOn(plot_mbc_projection, Name("mbc_sig_fit_result"), ProjectionRange("signal_box_mbc"), Components("final_signal"), LineColor(kRed), LineStyle(kDashDotted));
	final_pdf.plotOn(plot_mbc_projection, Name("mbc_bkg_fit_result"), ProjectionRange("signal_box_mbc"), Components("final_background"), LineColor(kBlue), LineStyle(kDashed));
	plot_mbc_projection->Draw("");
	canvas->SaveAs((opt->GetOutputDir()+"plot_mbc_projection.eps").c_str());
	delete plot_mbc_projection;

	// visualize the deltae fit result by cutting on other dimensions via the ProjectionRange argument
	canvas->Clear();
	RooPlot* plot_deltae_projection = de.frame(Bins(opt->GetBins()), Name("plot_deltae_projection"), Title("plot_deltae_projection"));
	data.plotOn(plot_deltae_projection, Name("data_hist"), CutRange("signal_box_deltae"), MarkerColor(kBlack));
	final_pdf.plotOn(plot_deltae_projection, Name("deltae_fit_result"), ProjectionRange("signal_box_deltae"), Components("final_pdf"), LineColor(kBlue), LineStyle(kSolid));
	final_pdf.plotOn(plot_deltae_projection, Name("deltae_sig_fit_result"), ProjectionRange("signal_box_deltae"), Components("final_signal"), LineColor(kRed), LineStyle(kDashDotted));
	final_pdf.plotOn(plot_deltae_projection, Name("deltae_bkg_fit_result"), ProjectionRange("signal_box_deltae"), Components("final_background"), LineColor(kBlue), LineStyle(kDashed));
	final_pdf.paramOn(plot_deltae_projection,Layout(0.6));
	plot_deltae_projection->Draw("");
	canvas->SaveAs((opt->GetOutputDir()+"plot_deltae_projection.eps").c_str());
	delete plot_deltae_projection;

	// visualize the contsupp fit result by cutting on other dimensions via the ProjectionRange argument

	 

	/**
	 * Clean up and exit program
	 */
	delete varset;
	delete canvas;
	delete chain;
	delete opt;

	return EXIT_SUCCESS;
}
