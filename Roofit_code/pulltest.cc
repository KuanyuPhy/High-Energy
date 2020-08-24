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
#include "RooHistPdf.h"
#include "TNtuple.h"
// program options
#include "options.h"
// Utilities
#include "utility.h"
// std library
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "RooCategory.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
using namespace std;
#include <fstream>
#include <iomanip>
#include <TRandom3.h>
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
		unsigned int nsig_events = 2000, unsigned int nbkg_events =3000) {

	RooDataSet *toy_mcdata_sig = pdf_signal.generate(argset,nsig_events);
	RooDataSet *toy_mcdata_bkg = pdf_background.generate(argset,nbkg_events);

	TFile *output_file = new TFile((opt->GetOutputDir()+"Ega_toymc_basic.root").c_str(),"RECREATE");
	TTree *output_tree = new TTree("h1","h1");
	double de;
	int mcinfo;
	//output_tree->Branch("Mb_c",&Mb_c,"Mb_c/F");
	output_tree->Branch("de",&de,"de/D");
	output_tree->Branch("mcinfo",&mcinfo,"mcinfo/I");

	for(unsigned int i = 0; i < nsig_events; ++i) {
		const RooArgSet* values = toy_mcdata_sig->get(i);
		//Mb_c = values->getRealValue("Mb_c");
		de = values->getRealValue("de");
		mcinfo = 1;
		output_tree->Fill();
	}

	for(unsigned int i = 0; i < nbkg_events; ++i) {
		const RooArgSet* values = toy_mcdata_bkg->get(i);
		//Mb_c = values->getRealValue("Mb_c");
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
	TCanvas *canvas = new TCanvas("plot", "plot", 1680, 1024);

	//----------資料是如何讀近來的呢？？？
	/**
	 * Part 1: Read input data and define variables
	 */
	// define the different variables
	//RooRealVar Mb_c("Mb_c", "M_{bc}", 5.2, 5.29, "GeV");
	RooRealVar de("de", "#Delta E", -0.2, 0.2, "GeV");
	//RooRealVar contsupp("contsupp", "Continuum Suppression", -15, 15);
	//RooCategory mcinfo("mcinfo", "MC Info"); //no use
	//mcinfo.defineType("background", 0);
	//mcinfo.defineType("signal", 1);

	// fill the RooDataSet
	RooArgSet *varset = new RooArgSet(de); //這裡的最要和root裡的依樣 , Convert n-tuple  to RooDataSet
	RooDataSet data("data", "data", chain, *varset);
	// Print some information of the RooDataSet (e.g. number of events)
	data.Print();

	// define full range of input data
	//Mb_c.setRange("full_range", 5.2, 5.29);
	de.setRange("full_range", -0.2, 0.2);
	//contsupp.setRange("full_range", -15, 15);

	// ----- 1-dim deltaE PDF -----
        RooRealVar mean("deltae_signal_mean", "de signal mean", 0.0);
        RooRealVar sigma_delta1("sigma_delta","sigma_delta",0.02718);//0.04,0,0.3
        RooRealVar alpha1("alpha","alpha",0.5946);
        RooRealVar n1("n","n",130);
        RooCBShape sig1("sig1", "Signal component 1",de,mean,sigma_delta1,alpha1,n1);

	RooRealVar mga1_mean("mga1_mean", "mga1_mean", 0.0);
        RooRealVar mga1_width("mga1_width", "mga1_width",0.01);
        RooGaussian sig2("sig2", "sig2",de, mga1_mean, mga1_width);
	RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal",0.9961);

        RooAddPdf deltae_signal("deltae_signal", "de signal_pdf", RooArgList(sig1,sig2), sig1frac);	
    
    // -------- bb de background ------------
    
        RooRealVar deltae_a0("deltae_a0", "deltae_a0", -0.5314);
        RooRealVar deltae_a1("deltae_a1", "deltae_a1", 0.023);
        //RooRealVar deltae_a2("deltae_a1", "deltae_a1", 0.1, -1.0, 1.0);
        RooChebychev deltae_bbackground("deltae_bbackground", "de b background",de, RooArgList(deltae_a0,deltae_a1));
    // --------- qq de background ------------
	 	RooRealVar deltae_a2("deltae_a2", "deltae_a2", -0.4406);
        RooRealVar deltae_a3("deltae_a3", "deltae_a3", 0.132);
        //RooRealVar deltae_a2("deltae_a1", "deltae_a1", 0.1, -1.0, 1.0);
        RooChebychev deltae_qbackground("deltae_qbackground", "de q background",de, RooArgList(deltae_a0,deltae_a1));

	//unsigned long nentries = (long)chain->GetEntries();
	// set reasonable starting values for signal and background events (for extened maximum likelihood fit)
	RooRealVar nsig("nsig","number signal events",3.0,-100.0,100);
	RooRealVar nbbkg("nbbkg","number b background events",123,-100.0,6000.0);//,0.0,10000.0);	
	RooRealVar nqbkg("nqbkg","number q background events",123,-100.0,6000.0);//,0.0,10000.0);
	RooAddPdf deltae_pdf("deltae_pdf", "deltae_pdf",RooArgList(deltae_signal,deltae_bbackground,deltae_qbackground),RooArgList(nsig,nbbkg,nqbkg));

	/**
	 * Part 3: Fit 1-dimensional PDFs to data
	 */
	//if (opt->GetSkipFit() == false)
	 //-------這裡在幹嘛呢？
		
		// fit to de
		RooFitResult *deltae_result = deltae_pdf.fitTo( data, NumCPU(opt->GetNumcpu()),RooFit::Minos(kTRUE), Timer(true), Extended(true), Save(true));
		
	
std::cout<<"1"<<std::endl;
		// in this example we dont't use the RooFitResults later, so we clean them up directly
		
		delete deltae_result;
		
	

	/**
	 * Part 4: Visualize the 1-dimensional fit results
	 */
	// visualize the Mb_c fit result
	
	// visualize the de fit result
	std::cout<<"2"<<std::endl;

	canvas->Clear();
	RooPlot *plot_deltae = de.frame(Bins(opt->GetBins()), Name("plot_deltae"), Title("plot_deltae"));
	data.plotOn(plot_deltae, Name("data_hist"), CutRange("full_range"), MarkerColor(kBlack));
	deltae_pdf.plotOn(plot_deltae, Name("deltae_fit_result"), ProjectionRange("full_range"), Components("deltae_pdf"), LineColor(kBlue), LineStyle(kSolid));
	//------
	deltae_pdf.plotOn(plot_deltae, Name("deltae_sig_fit_result"), ProjectionRange("full_range"), Components("deltae_signal"), LineColor(kRed), LineStyle(kDashDotted));
	deltae_pdf.plotOn(plot_deltae, Name("deltae_bbkg_fit_result"), ProjectionRange("full_range"), Components("deltae_bbackground"), LineColor(kBlue), LineStyle(kDashed));
	deltae_pdf.plotOn(plot_deltae, Name("deltae_qbkg_fit_result"), ProjectionRange("full_range"), Components("deltae_qbackground"), LineColor(kGreen), LineStyle(kDashed));
	
	//put information on
	deltae_pdf.paramOn(plot_deltae,Layout(0.6));
	//deltae_pdf.paramOn(plot_deltae, Layout(0.6), Format("NEU", AutoPrecision(2)));
	//data.statOn(plot_deltae);
	// TODO: Add signal and background component to the deltaE plot
	plot_deltae->Draw("");
	canvas->SaveAs((opt->GetOutputDir() + "1D_plot_deltae.pdf").c_str());
	//RooAbsPdf::paramOn(frame)

	//float xxx=deltae_signal_mean.getVal();
	//std::std::cout<<"Mean"<<xxx;

	delete plot_deltae;

	std::cout<<"3"<<std::endl;

	
	/**
	 * Part 8: Produce a ToyMC from a PDF
	 */
	// this is a method that gets two PDFs, a list of variables, options objcet and numbers of signal and background
	// events to generate a ToyMC from it.
	//Christal Ball
   

	// Gaussian +++ Christal Ball
	

	//ProduceToyMCTuples(deltae_signal,deltae_background,de,opt,2,983);

	/**
	 * Part 9: Do an ensemble test
	 */
/*
	if(opt->GetSkipEnsembleTest() == false) 
	{
		// define the ensemble test object by passing PDF that is used for generation and fitting (not necessarily the same,
		// check documentation/examples for the other case). Also pass variables that are generated and whether it is a
		// Extended fit or not (depends on your PDF). Silence(true) makes it less noisy.
		// FitOptions can be everything you would pass to a normal fit. e.g. that it is a extended fit.
		RooMCStudy *ensemble_test = new RooMCStudy(deltae_pdf,RooArgSet(de),Extended(true),Silence(true),
				FitOptions(NumCPU(opt->GetNumcpu()), Timer(true), Extended(true), Save(true)));
		// now run the ensemble test 100 times
		ensemble_test->generateAndFit(1000);

		// create an empty canvas to draw all the ensemble tests plots on
		TCanvas *canvas_ensemble = new TCanvas("ensemble_plot","ensemble_plot",1280*2,1024*2);
		canvas_ensemble->Divide(2,2);
		canvas_ensemble->cd(1);
		// this plot contains the 100 fitted results of deltaE signal mean
		RooPlot* frame1 = ensemble_test->plotParam(nsig,Bins(opt->GetBins()));
		frame1->Draw();
		canvas_ensemble->cd(2);
		// this plot contains the 100 fitted errors of deltaE signal mean
		RooPlot* frame2 = ensemble_test->plotError(nsig,Bins(opt->GetBins()));
		frame2->Draw();
		canvas_ensemble->cd(3);
		// this plot contains the 100 fitted pulls of deltaE signal mean
		RooPlot* frame3 = ensemble_test->plotPull(nsig,Bins(opt->GetBins()),FitGauss(true));
		frame3->Draw();
		canvas_ensemble->cd(4);
		// this plot contains the negative log likelihood for each of the 100 fits
		RooPlot* frame4 = ensemble_test->plotNLL(Bins(opt->GetBins()));
		frame4->Draw();
		canvas_ensemble->SaveAs((opt->GetOutputDir()+"ensemble_test.eps").c_str());
		delete canvas_ensemble;

TCanvas *canvas_ensemble1 = new TCanvas("ensemble_plot","ensemble_plot",1680,1024);
RooPlot* frame = ensemble_test->plotPull(nsig,Bins(opt->GetBins()),FitGauss(true));
		frame->Draw();

canvas_ensemble1->SaveAs((opt->GetOutputDir() + "pull.pdf").c_str());
delete canvas_ensemble;
	}
*/
///////////////////////////////////////////////////////
//--------------------

//toy start
 
 RooRealVar pulls( "pulls", "pulls", -5,5 );
 RooDataSet pulldata("pulldata","pulldata",RooArgSet(pulls)) ;
 //double fittedsig,fitted;
int runnum=2000;
TRandom3 a;
TRandom3 b;
TRandom3 c;
for(int i=0;i<runnum;i++)
{
RooRealVar N_sig("N_sig", "N_sig",200,-1000.0,5000.0);
RooRealVar N_bb("N_bb", "N_bb",983,-1000.0, 20000.0);
RooRealVar N_qq("N_qq", "N_qq",983,-1000.0, 20000.0);
      
        
 
 float toynumsig=a.PoissonD(1.74),toynumbb=b.PoissonD(918),toynumqq=c.PoissonD(1916);// my expect value out put with poisson distribution
 
 RooDataSet* toysig = deltae_signal.generate(RooArgList(de),toynumsig);  //generate toy MC by my model 
 RooDataSet* toybb = deltae_bbackground.generate(RooArgList(de),toynumbb);
 RooDataSet* toyqq = deltae_qbackground.generate(RooArgList(de),toynumqq);

 RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(deltae_signal,deltae_bbackground,deltae_qbackground),RooArgList(N_sig,N_bb,N_qq));

 toysig->append(*toybb);
 toysig->append(*toyqq);
 
RooDataSet data2=*toysig;
 RooFitResult* ToyfitRes = final_pdf.fitTo(*toysig,Extended(),RooFit::Minos(kTRUE),RooFit::Save(kTRUE));
 float fitsig=N_sig.getVal(); //get N signal value
 //float fitbb=N_bb.getVal();

delete ToyfitRes;

canvas->Clear();
RooPlot *plot_deltae = de.frame(Bins(opt->GetBins()), Name("plot_deltae"), Title("plot_deltae"));
	data2.plotOn(plot_deltae, Name("data_hist"), CutRange("full_range"), MarkerColor(kBlack));
	final_pdf.plotOn(plot_deltae, Name("deltae_fit_result"), ProjectionRange("full_range"), Components("deltae_pdf"), LineColor(kBlue), LineStyle(kSolid));
	//------
	final_pdf.plotOn(plot_deltae, Name("deltae_sig_fit_result"), ProjectionRange("full_range"), Components("deltae_signal"), LineColor(kRed), LineStyle(kDashDotted));
	final_pdf.plotOn(plot_deltae, Name("deltae_bbkg_fit_result"), ProjectionRange("full_range"), Components("deltae_bbackground"), LineColor(kBlue), LineStyle(kDashed));
	final_pdf.plotOn(plot_deltae, Name("deltae_qbkg_fit_result"), ProjectionRange("full_range"), Components("deltae_qbackground"), LineColor(kGreen), LineStyle(kDashed));

	final_pdf.paramOn(plot_deltae,Layout(0.6));
	plot_deltae->Draw("");
	char xxx[32];
	sprintf (xxx,"%s%d%s","pull",i,".pdf"); // print each picture
	canvas->SaveAs((opt->GetOutputDir() + xxx).c_str());


	delete plot_deltae;

 //改成正和負的error，fit出來和當初丟進去的差,如過是正的error就存正的error，如果是負的就存負的error

 //float sigerror=(N_sig.getAsymErrorHi()-N_sig.getAsymErrorLo())/2;
 float sigerror1=abs((N_sig.getAsymErrorLo()));
 float sigerror2=(N_sig.getAsymErrorHi());
//float bberror=(N_bb.getAsymErrorHi()-N_bb.getAsymErrorLo())/2;
if (sigerror1==0)
{continue;}

else {
if (fitsig>toynumsig)
    {
        pulls=(fitsig-toynumsig)/sigerror1;
    }
    else if (fitsig<toynumsig)
    {
        pulls=(fitsig-toynumsig)/sigerror2;
    }  

std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxFitSig "<<i<<" =  "<<fitsig<<std::endl;
    std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxtoynumsig "<<i<<" =  "<<toynumsig<<std::endl;
  std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxsigerror "<<i<<" =  "<<sigerror1<<std::endl;
  std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxN_sig.getAsymErrorHi "<<i<<" =  "<<N_sig.getAsymErrorHi()<<std::endl;
  std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxN_sig.getAsymErrorLo()"<<i<<" =  "<<N_sig.getAsymErrorLo()<<std::endl;

 std::cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxpull =   "<<(fitsig-toynumsig)/sigerror2<<std::endl;
 pulldata.add(RooArgSet(pulls));
}
}
  RooRealVar gmean("gmean","mean of gaussian",-0.05,-0.15,0.15) ;
  RooRealVar sigma("sigma","width of gaussian",0.91,0,1.5) ;
  RooGaussian gauss("gauss","gaussian PDF",pulls,gmean,sigma) ;
  gauss.fitTo(pulldata);
 canvas->Clear();
      RooPlot *pullplot_result = pulls.frame(Bins(40));
    pulldata.plotOn( pullplot_result, Name("Pull") , LineColor(kBlack), LineStyle(kSolid));
   gauss.plotOn( pullplot_result, Name("gauss") ,Components("gauss"), LineColor(kBlue), LineStyle(kSolid));
gauss.paramOn(pullplot_result,Layout(0.6));
    pullplot_result->Draw("");
canvas->SaveAs((opt->GetOutputDir() + "Ega_1D_plot_pull.pdf").c_str());

//toy end----------------------------------------



//--------------------

	/**
	 * Clean up and exit program
	 */
	delete varset;
	delete canvas;
	delete chain;
	delete opt;

	return EXIT_SUCCESS;


}


