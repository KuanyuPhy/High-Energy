#include <TFile.h>

void Jpsi();
void kstar0();
void B0Mbc();
void B0dE();

void eegg()
{

	Jpsi();
	kstar0();
	B0Mbc();
	B0dE();
}

void Jpsi()
{
	//TString title = "MC Signal B0 #rightarrow J/#psi(e+e-) #k(#gamma#gamma) J/psi Mass distribution";
	TString title = "J/#psi(e+e-)";
	TString xaxistitle = "Mass(GeV/c^2)";
	TString yaxistitle = "Number of events(counts)";
	TString variables = "jpsimass"; ////////////
	Int_t bins = 300;
	Float_t Lower = 2.95;
	Float_t Upper = 3.15;
	TCanvas *c = new TCanvas;
	//TFile *file = new TFile("gSKIM_B");
	TFile *file = new TFile("B0_Jpsi_Kstar0_test5.root");
	//	TFile *file = new TFile("gSKIM_B0_Jpsi_ee_gamma.root");
	TTree *tree = (TTree *)file->Get("h1");
	Float_t M;
	Int_t entries = (Int_t)tree->GetEntries();
	tree->SetBranchAddress(variables, &M);
	//eta -> gamma m 0.54~0.5
	//TH1F *hM = new TH1F("M","M distribution with signal",100,0.545,0.55);

	//Jpsi -> e e
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

void kstar0()
{
	//TString title = "MC Signal B0 #rightarrow J/#psi(e+e-) #eta(#gamma#gamma) eta Mass distribution";
	TString title = "kstar0";
	TString xaxistitle = "Mass(GeV/c^2)";
	TString yaxistitle = "Number of events(counts)";
	TString variables = "k_star0m";
	Int_t bins = 100;
	Float_t Lower = 0.6;
	Float_t Upper = 1.6;
	TCanvas *c = new TCanvas;
	//TFile *file = new TFile("gSKIM_B0_eta_gg.root");
	TFile *file = new TFile("B0_Jpsi_Kstar0_test5.root");
	//	TFile *file = new TFile("gSKIM_B0_Jpsi_ee_gamma.root");
	TTree *tree = (TTree *)file->Get("h1");
	Float_t M;
	Int_t entries = (Int_t)tree->GetEntries();
	tree->SetBranchAddress(variables, &M);
	//eta -> gamma m 0.54~0.5
	//TH1F *hM = new TH1F("M","M distribution with signal",100,0.545,0.55);

	//Jpsi -> e e
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
	img->WriteImage("k_star0m.png");
}

void B0Mbc()
{
	//TString title = "MC Signal B0 #rightarrow J/#psi(e+e-) #eta(#gamma#gamma) Mbc distribution";
	TString title = "mbc";
	TString xaxistitle = "Mass(GeV/c^2)";
	TString yaxistitle = "Number of events(counts)";
	TString variables = "Mb_c";
	Int_t bins = 300;
	Float_t Lower = 5.2;
	Float_t Upper = 5.3;
	TCanvas *c = new TCanvas;
	//TFile *file = new TFile("gSKIM_B0_eta_gg.root");
	TFile *file = new TFile("B0_Jpsi_Kstar0_test5.root");
	//	TFile *file = new TFile("gSKIM_B0_Jpsi_ee_gamma.root");
	TTree *tree = (TTree *)file->Get("h1");
	Float_t M;
	Int_t entries = (Int_t)tree->GetEntries();
	tree->SetBranchAddress(variables, &M);
	//eta -> gamma m 0.54~0.5
	//TH1F *hM = new TH1F("M","M distribution with signal",100,0.545,0.55);

	//Jpsi -> e e
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
	img->WriteImage("Mbc.png");
}

void B0dE()
{
	//TString title = "MC Signal B0 #rightarrow J/#psi(e+e-) #eta(#gamma#gamma) #DeltaE distribution";
	TString title = "#DeltaE";
	TString xaxistitle = "#DeltaE(GeV)";
	TString yaxistitle = "Number of events(counts)";
	TString variables = "de";
	Int_t bins = 100;
	Float_t Lower = -0.2;
	Float_t Upper = 0.2;
	TCanvas *c = new TCanvas;
	TFile *file = new TFile("B0_Jpsi_Kstar0_test5.root");
	TTree *tree = (TTree *)file->Get("h1");
	Float_t M;
	Int_t entries = (Int_t)tree->GetEntries();
	tree->SetBranchAddress(variables, &M);
	//eta -> gamma m 0.54~0.5
	//TH1F *hM = new TH1F("M","M distribution with signal",100,0.545,0.55);

	//Jpsi -> e e
	
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
	//img->WriteImage("dE_B0_B02eegg.png");
}
