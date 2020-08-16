#include "TROOT.h"
#include "TMath.h"
// ROOT data
 #include "TFile.h"
 #include "TChain.h"
 #include "TLeaf.h"
// // ROOT plot
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
 #include "TH2.h"
 #include "TH3.h"
 #include "RooCategory.h"
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 #include <fstream>
 #include <iomanip>
#include <fstream>
#include <iomanip>


        using namespace RooFit;
int main() {
        TChain*te = new TChain("h1");
        TChain*cs = new TChain("h1");
        TChain*bg = new TChain("h1");
        TChain*tot = new TChain("h1");
        te->Add("true.root");
        cs->Add("cross.root");
        bg->Add("bg.root");
        tot->Add("tot.root");
            
        RooRealVar mbc("mbc","mbc",5.24,5.29,"GeV");
        RooRealVar de("de","de",-0.5,0.5,"GeV");
        RooRealVar massrho0("massrho0","massrho0",0.6,1.4,"GeV");
        RooRealVar nbtr("nbtr","nbtr",-10,10,"GeV");
        RooArgSet *varset = new RooArgSet(mbc,de,massrho0,nbtr);
        RooDataSet data("te","te",te,*varset);
        RooDataSet data2("cs","cs",cs,*varset);
        RooDataSet data3("bg","bg",bg,*varset);
        RooDataSet data4("tot","tot",tot,*varset);
        TH2F *m1v = new TH2F("m1v","m1v",100,5.24,5.29,100,-0.5,0.5);        
        TH2F *d12 = new TH2F("d12", "d12" ,100,5.24,5.29,100,-0.5,0.5);
        TH3F *d23 = new TH3F("d22", "d23" ,100,5.24,5.29,100,-0.5,0.5,100,0.6,1.4);
        TH1F *d31 = new TH1F("d31", "d31" ,100,0.6,1.4);
       // TH3F*3d3 = new TH3F(100,5.24,5.29,100,-0.5,0.5,100,0.6,1.4);
        te->Draw("mbc:de>>d12");
        cs->Draw("mbc:de:massrho0>>d23");
        bg->Draw("massrho0>>d31");
      //  bg->Draw("mbc:de:massrho0>>3d3");
        RooDataHist teh("teh","teh",RooArgList(de,mbc),2d1);
        RooDataHist csh("csh","csh",RooArgList(massrho0,de,mbc),3d2);
        RooDataHist bgh("bgh","bgh",RooArgList(massrho0),1d3);
       // RooDataHist bgh("bgh","bgh",RooArgList(massrho0,de,mbc),3d3);
        RooHistPdf tep("tep","tep",RooArgList(de,mbc),teh,0);
        RooHistPdf csp("csp","csp",RooArgList(massrho0,de,mbc),csh,0);
        RooHistPdf bgp("bgp","bgp",RooArgList(massrho0),bgh,0);
       // RooHistPdf bgp("bgp","bgp",RooArgList(massrho0,de,mbc),bgh,0);

      //  data.Print();
        mbc.setRange("full_range",5.24,5.29);
        RooRealVar mbc_bg_m0 ("mbc_bg_m0","mbc_bg_m0",5.289);
        RooRealVar mbc_bg_c ("mbc_bg_c","mbc_bg_c",-0.106);
        RooRealVar mbc_bg_p ("mbc_bg_p","mbc_bg_p",0.5);
        RooArgusBG mbc_bg("mbc_bg","mbc_bg",mbc,mbc_bg_m0,mbc_bg_c,mbc_bg_p); 

     
        de.setRange("full_range2",-0.5,0.5);
       // RooRealVar de_signal_mean ("de_signal_mean","de_signal_mean",0.0);
       // RooRealVar de_signal_width ("de_signal_width","de_signal_width",0.0092);
       // RooGaussian de_signal("de_signal","de_signal",de,de_signal_mean,de_signal_width);
       // RooRealVar de_background_width ("de_background_width","de_background_width",0.13);
        RooRealVar de_a ("de_a","de_a",0.24);
        RooRealVar de_b ("de_b","de_b",-0.677);
        RooRealVar de_c ("de_c","de_c",-0.061);
        RooChebyShev de_bg("de_bg","de bg",de,RooArgList(de_b,de_a,de_c));            
      // RooGaussian de_signal("de_signal","de_signal",de,de_signal_mean,de_signal_width);
      // RooGaussian de_background("de_background","de background",de,de_signal_mean,de_background_width);
      // RooRealVar de_signal_fraction("de_signal_fraction","de_signal_fraction",0.94);
      // RooAddPdf de_true("de_true","de_true",de_signal,de_background,de_signal_fraction);


        massrho0.setRange("full_range3",0.6,1.4);
        RooRealVar rho0_mean("rho0_mean","rho0_mean",0.765);
        RooRealVar rho0_width("rho0_width","rho0_width",0.132);
        RooBreitWigner rho0("rho0","rho0",massrho0,rho0_mean,rho0_width);




        nbtr.setRange("full_range4",-10,10);
        RooRealVar nbtr_true_mean("nbtr_true_mean","nbtr_true_mean",1.817);
        RooRealVar nbtr_true_width("nbtr_true_width","nbtr_true_width",1.8617);
        RooGaussian nbtr_true("nbtr_true","nbtr_true",nbtr,nbtr_true_mean,nbtr_true_width);
        RooRealVar nbtr_cross_signal_mean ("nbtr_cross_signal_mean","nbtr_cross_signal_mean",1.21);
        RooRealVar nbtr_cross_background_mean ("nbtr_cross_background_mean","nbtr_cross_background_mean",0.551);
        RooRealVar nbtr_cross_signal_width ("nbtr_cross_signal_width","nbtr_cross_signal_width",1.711);
        RooRealVar nbtr_cross_background_width ("nbtr_cross_background_width","nbtr_cross_background_width",2.54);
        RooGaussian nbtr_cross_signal("nbtr_cross_signal","nbtr_cross_signal",nbtr,nbtr_cross_signal_mean,nbtr_cross_signal_width);
        RooGaussian nbtr_cross_background("nbtr_cross_background","nbtr_cross_background",nbtr,nbtr_cross_background_mean,nbtr_cross_background_width);
        RooRealVar nbtr_cross_fraction("nbtr_cross_fraction","nbtr_cross_fraction",0.8);
        RooAddPdf nbtr_cross("nbtr_cross","nbtr_cross",RooArgList(nbtr_cross_signal,nbtr_cross_background),RooArgList(nbtr_cross_fraction));
        RooRealVar nbtr_signal_mean ("nbtr_signal_mean","nbtr_signal_mean",-0.19);
        RooRealVar nbtr_background_mean ("nbtr_background_mean","nbtr_background_mean",-0.26);
        RooRealVar nbtr_signal_width ("nbtr_signal_width","nbtr_signal_width",2.21);
        RooRealVar nbtr_background_width ("nbtr_background_width","nbtr_background_width",1.39);
        RooBreitWigner nbtr_signal("nbtr_signal","nbtr_signal",nbtr,nbtr_signal_mean,nbtr_signal_width);
        RooRealVar nbtr_sec_mean("nbtr_sec_mean","nbtr_sec_mean",-1.19);
        RooRealVar nbtr_sec_width("nbtr_sec_width","nbtr_sec_width",2.21);
        RooGaussian nbtr_background("nbtr_background","nbtr background",nbtr,nbtr_background_mean,nbtr_background_width);
        RooGaussian nbtr_sec("nbtr_sec","nbtr_sec",nbtr,nbtr_sec_mean,nbtr_sec_width);
        RooRealVar nbtr_signal_fraction("nbtr_signal_fraction","nbtr_signal_fraction",0.165);
        RooRealVar nbtr_sec_fraction("nbtr_sec_fraction","nbtr_sec_fraction",0.78);
        RooAddPdf nbtr_sec_pdf("nbtr_sec_pdf","nbtr_sec_pdf",RooArgList(nbtr_signal,nbtr_background),RooArgList(nbtr_signal_fraction));
        RooAddPdf nbtr_bg("nbtr_bg","nbtr_bg",RooArgList(nbtr_sec_pdf,nbtr_sec),RooArgList(nbtr_sec_fraction));
       










        RooProdPdf final_true ("final_true","final_true_pdf",RooArgList(tep,rho0,nbtr_true));
        RooProdPdf final_cross ("final_cross","final_cross_pdf",RooArgList(csp,nbtr_cross));
        RooProdPdf final_bg ("final_bg","final_bg_pdf",RooArgList(bgp,mbc_bg,de_bg,nbtr_bg));
        unsigned long nentries=(long)tot->GetEntries();
        RooRealVar nsig ("nsig","signal events",nentries*0.4,0,nentries*1.2);
        RooRealVar nbg ("nbg","bg events ",nentries*0.6,0,nentries*1.2);
        RooAddPdf final_pdf("final_pdf","final_pdf",RooArgList(final_true,final_cross,final_bg),RooArgList(0.915*nsig,0.085*nsig,nbg));
        RooFitResult *final_result = new RooFitResult("mbc_final_result","mbc_final_result");
        final_result = final_pdf.fitTo(data4, Timer(true), Save(true));
        mbc.setRange("signal_box_mbc",5.2,5.29);
        de.setRange("signal_box_de",-0.5,0.5);
        massrho0.setRange("signal_box_massrho0",0.6,1.4);
        nbtr.setRange("full_range4",-10,10);
       
        TCanvas*canvas= new TCanvas("canvas","c1");
        canvas->Clear();
        RooPlot*plot_mbc_projection=mbc.frame(Bins(100),Name("plot_mbc_projection"),Title("plot_mbc_projection"));
        data4.plotOn(plot_mbc_projection,Name("data_hist"),CutRange("signal_box_mbc"),MarkerColor(kBlack));
        final_pdf.plotOn(plot_mbc_projection,Name("mbc_fit_result"),ProjectionRange("signal_box_mbc"),Components("final_pdf"),LineColor(kBlue),LineStyle(kSolid));
        final_pdf.plotOn(plot_mbc_projection,Name("mbc_fit_result"),ProjectionRange("signal_box_mbc"),Components("final_true"),LineColor(kRed),LineStyle(kSolid));
        final_pdf.plotOn(plot_mbc_projection,Name("mbc_fit_result"),ProjectionRange("signal_box_mbc"),Components("final_cross"),LineColor(kGreen),LineStyle(kSolid));
        final_pdf.plotOn(plot_mbc_projection,Name("mbc_fit_result"),ProjectionRange("signal_box_mbc"),Components("final_bg"),LineColor(kViolet),LineStyle(kSolid));
        plot_mbc_projection->Draw("");

        TCanvas*canvas2= new TCanvas("canvas2","c2");
        canvas2->Clear();
        RooPlot*plot_de_projection=de.frame(Bins(100),Name("plot_de_projection"),Title("plot_de_projection"));
        data4.plotOn(plot_de_projection,Name("data_hist"),CutRange("signal_box_de"),MarkerColor(kBlack));
        final_pdf.plotOn(plot_de_projection,Name("de_fit_result"),ProjectionRange("signal_box_de"),Components("final_pdf"),LineColor(kBlue),LineStyle(kSolid));
        final_pdf.plotOn(plot_de_projection,Name("de_fit_result"),ProjectionRange("signal_box_de"),Components("final_true"),LineColor(kRed),LineStyle(kSolid));
        final_pdf.plotOn(plot_de_projection,Name("de_fit_result"),ProjectionRange("signal_box_de"),Components("final_cross"),LineColor(kGreen),LineStyle(kSolid));
        final_pdf.plotOn(plot_de_projection,Name("de_fit_result"),ProjectionRange("signal_box_de"),Components("final_bg"),LineColor(kViolet),LineStyle(kSolid));
        plot_de_projection->Draw("");
       
        TCanvas*canvas3= new TCanvas("canvas3","c3");
        canvas3->Clear();
        RooPlot*plot_massrho0_projection=massrho0.frame(Bins(100),Name("plot_massrho0_projection"),Title("plot_massrho0_projection"));
        data4.plotOn(plot_massrho0_projection,Name("data_hist"),CutRange("signal_box_massrho0"),MarkerColor(kBlack));
        final_pdf.plotOn(plot_massrho0_projection,Name("massrho0_fit_result"),ProjectionRange("signal_box_massrho0"),Components("final_pdf"),LineColor(kBlue),LineStyle(kSolid));
        final_pdf.plotOn(plot_massrho0_projection,Name("massrho0_fit_result"),ProjectionRange("signal_box_massrho0"),Components("final_true"),LineColor(kRed),LineStyle(kSolid));
        final_pdf.plotOn(plot_massrho0_projection,Name("massrho0_fit_result"),ProjectionRange("signal_box_nassrho0"),Components("final_cross"),LineColor(kGreen),LineStyle(kSolid));
        final_pdf.plotOn(plot_massrho0_projection,Name("massrho0_fit_result"),ProjectionRange("signal_box_massrho0"),Components("final_bg"),LineColor(kViolet),LineStyle(kSolid));
        plot_massrho0_projection->Draw("");

        TCanvas*canvas4= new TCanvas("canvas4","c4");
        canvas4->Clear();
        RooPlot*plot_nbtr_projection=nbtr.frame(Bins(100),Name("plot_nbtr_projection"),Title("plot_nbtr_projection"));
        data4.plotOn(plot_nbtr_projection,Name("data_hist"),CutRange("signal_box_nbtr"),MarkerColor(kBlack));
        final_pdf.plotOn(plot_nbtr_projection,Name("nbtr_fit_result"),ProjectionRange("signal_box_nbtr"),Components("final_pdf"),LineColor(kBlue),LineStyle(kSolid));
        final_pdf.plotOn(plot_nbtr_projection,Name("nbtr_fit_result"),ProjectionRange("signal_box_nbtr"),Components("final_true"),LineColor(kRed),LineStyle(kSolid));
        final_pdf.plotOn(plot_nbtr_projection,Name("nbtr_fit_result"),ProjectionRange("signal_box_nbtr"),Components("final_cross"),LineColor(kGreen),LineStyle(kSolid));
        final_pdf.plotOn(plot_nbtr_projection,Name("nbtr_fit_result"),ProjectionRange("signal_box_nbtr"),Components("final_bg"),LineColor(kViolet),LineStyle(kSolid));
        plot_nbtr_projection->Draw("");


        cout<<"nentries = "<<nentries<<endl;
}
