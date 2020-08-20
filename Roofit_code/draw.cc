#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void eta();
void chic0();
void B0Mbc();
void B0dE();
etaadd()

void draw()
{
    //eta();
	etaadd()
    //chic0();
    //B0Mbc();
    //B0dE();
}

void eta()
{	
    TString title = "#eta(#gamma#gamma)";
	TString xaxistitle = "Mass(GeV/c^{2})";
	TString yaxistitle = "events/0.0009GeV/c^{2}";
	TString variables = "etam";
    Int_t bins = 180;
	Float_t Lower = 0.45;
	Float_t Upper = 0.62;
    TCanvas *c = new TCanvas;
	//TFile *file = new TFile("gSKIM_B");
	TFile *file = new TFile("cod010701.root");
	//	TFile *file = new TFile("gSKIM_B0_Jpsi_ee_gamma.root");
	TTree *tree = (TTree *)file->Get("h1");
	Float_t M;
	Int_t entries = (Int_t)tree->GetEntries();
	tree->SetBranchAddress(variables, &M);
    
    TH1F *hM = new TH1F(variables, title, bins, Lower, Upper);
	for (Int_t i = 0; i < entries; i++)
	{
		tree->GetEntry(i);
		hM->Fill(M);
	}
	hM->GetXaxis()->SetTitle(xaxistitle);
	hM->GetYaxis()->SetTitle(yaxistitle);
	hM->Draw();
	gSystem->ProcessEvents();
	TImage *img = TImage::Create();
	img->FromPad(c);
	img->WriteImage("M_Jpsi_B02eegg.png");
}
void etaadd()
{	
    TString title = "#eta(#gamma#gamma)";
	TString xaxistitle = "Mass(GeV/c^{2})";
	TString yaxistitle = "events/0.0009GeV/c^{2}";
	TString variables = "etam";
    Int_t bins = 180;
	Float_t Lower = 0.45;
	Float_t Upper = 0.62;
    TCanvas *c = new TCanvas;
	//TFile *file = new TFile("gSKIM_B");
	TFile *file = new TFile("cod010701.root");
	TChain *qq=new TChain("h1");
	TChain *bb=new TChain("h1")
	//	TFile *file = new TFile("gSKIM_B0_Jpsi_ee_gamma.root");
	TTree *tree = (TTree *)file->Get("h1");
	Float_t M;
	Int_t entries = (Int_t)tree->GetEntries();
	tree->SetBranchAddress(variables, &M);
    
    TH1F *hM = new TH1F(variables, title, bins, Lower, Upper);
	for (Int_t i = 0; i < entries; i++)
	{
		tree->GetEntry(i);
		hM->Fill(M);
	}
	hM->GetXaxis()->SetTitle(xaxistitle);
	hM->GetYaxis()->SetTitle(yaxistitle);
	hM->Draw();
	gSystem->ProcessEvents();
	TImage *img = TImage::Create();
	img->FromPad(c);
	img->WriteImage("M_Jpsi_B02eegg.png");
}
