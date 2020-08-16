#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
void roofit_qq()
{
    TChain *chain = new TChain("h1");
    chain->Add("bbmd.root");
    using namespace RooFit;
    RooRealVar mbc("mbc", "M_{bc}", 5.24, 5.29, "GeV");
    RooRealVar de("de", "#Delta E", -0.2, 0.2, "GeV");
    //RooRealVar masspsi("masspsi","masspsi",3.0,3.5,"GeV");
    RooCategory mcinfo("mcinfo", "MCIN");
    mcinfo.defineType("background", 0);
    mcinfo.defineType("signal", 1);
    RooArgSet *varset = new RooArgSet(de, mbc);
    RooDataSet data("data", "data", chain, *varset);
    data.Print();
    mbc.setRange("full_range",5.24,5.29);
    de.setRange("full_range", -0.2, 0.2);

    ///////////////////////
    //de signal
    ///////////////////////
	RooRealVar de_a0("de_a0", "de_a0",0, -100, 100);
        RooRealVar de_a1("de_a1", "de_a1",0, -5.0, 10.0);
        RooChebychev deltae_pdf("deltae_pdf", "de signal_pdf",de, RooArgList(de_a0,de_a1));
	//RooPolynomial deltae_pdf("deltae_pdf","de signal_pdf", de,RooArgList(de_a0,de_a1));
	//RooChebychev deltae_pdf("deltae_pdf", "de signal_pdf",de, de_a0);
	///////////////////////
        //mbc signal
        ///////////////////////
	RooRealVar argpar("argpar","argus shape parameter",-20.0,-100.,-1.);
	RooArgusBG mbc_signal("mbc_signal","Argus PDF",mbc,RooConst(5.29),argpar);

	RooFitResult *de_result = new RooFitResult("de_result", "de_result");
    de_result = deltae_pdf.fitTo(data, Timer(true), Save(true));
    TCanvas *canvas = new TCanvas("canvas", "c1");
    canvas->Clear();
    RooPlot *plot_deltae = de.frame(Bins(100), Name("plot_deltae"), Range("full_range"), Title("plot_deltae"));
    data.plotOn(plot_deltae, Name("data_hist"), CutRange("full_range"), MarkerColor(kBlack));
    deltae_pdf.plotOn(plot_deltae, Name("deltae_fit_result"), ProjectionRange("full_range"), Components("deltae_pdf"), LineColor(kBlue), LineStyle(kSolid));
    //put information on
    deltae_pdf.paramOn(plot_deltae, Layout(0.6));   
    plot_deltae->Draw("");
    canvas->SaveAs("bb1D_plot_de.pdf");
    delete plot_deltae;

    canvas->Clear();
    RooFitResult *mbc_result = new RooFitResult("mbc_result", "mbc_result");
    mbc_result = mbc_signal.fitTo(data, Timer(true), Save(true));
    TCanvas *canvas1 = new TCanvas("canvas1", "c2");
    
    RooPlot *plot_mbc = mbc.frame(Bins(100), Name("plot_mbc"), Range("full_range"), Title("plot_mbc"));
    data.plotOn(plot_mbc, Name("data_hist"), CutRange("full_range"), MarkerColor(kBlack));
    mbc_signal.plotOn(plot_mbc, Name("mbc_fit_result"), ProjectionRange("full_range"), Components("mbc_signal"), LineColor(kRed), LineStyle(kSolid));
    //------
    //deltae_pdf.plotOn(plot_deltae, Name("deltae_sig_fit_result"), ProjectionRange("full_range"), Components("deltae_signal"), LineColor(kRed), LineStyle(kDashDotted));
    //deltae_pdf.plotOn(plot_deltae, Name("deltae_bkg_fit_result"), ProjectionRange("full_range"), Components("deltae_background"), LineColor(kBlue), LineStyle(kDashed));
 
    mbc_signal.paramOn(plot_mbc, Layout(0.6));
    //deltae_pdf.paramOn(plot_deltae, Layout(0.6), Format("NEU", AutoPrecision(2)));
    //data.statOn(plot_deltae);
    // TODO: Add signal and background component to the deltaE plot
 plot_mbc->Draw("");
 canvas1->SaveAs("bb1D_plot_mbc.pdf");
    
}
