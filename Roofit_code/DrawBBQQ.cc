void DrawBBQQ(){ 
	TChain *BB = new TChain("h1");
	//BB->Add("~/Desktop/rootfiles/x3872gammamumuBB_st0_cut1.root");
	//BB->Add("~/Desktop/rootfiles/x3872gammamumuBB_st1-9_cut1.root");
	//BB->Add("~/Desktop/rootfiles/step1_x3872gammamumuBB_st1-9_cut.root");
	//BB->Add("~/Desktop/rootfiles/step1_x3872gammamumuBB_st1-9_cutdm.root");
	//BB->Add("~/Desktop/rootfiles/step1_x3872gammaeeBB_cut.root");
	//BB->Add("~/Desktop/rootfiles/step1_x3872gammaeeBB_cutdm.root");
	BB->Add("~/fs02/bgMC/k0pipigammaBB.root");
	//BB->Add("~/BAS_2012May/neurobayes_kppg_qq/step1_k0pipigammaBB_0*.root");
	//BB->Add("~/Desktop/rootfiles/step1_k0pipigammaBB_cutnb.root");
	//BB->Add("~/Desktop/rootfiles/step1_k0pipigammaBB_cutdm.root");

	TChain *QQ = new TChain("h1");
	//QQ->Add("~/Desktop/rootfiles/x3872gammamumuQQ_cut1.root");
	//QQ->Add("~/Desktop/rootfiles/step1_x3872gammamumuQQ_cut.root");
	//QQ->Add("~/Desktop/rootfiles/step1_x3872gammamumuQQ_cutdm.root");
	//QQ->Add("~/Desktop/rootfiles/step1_x3872gammaeeQQ_cut.root");
	//QQ->Add("~/Desktop/rootfiles/step1_x3872gammaeeQQ_cutdm.root");
	QQ->Add("~/fs02/bgMC/k0pipigammaQQ.root");
	//QQ->Add("~/BAS_2012May/neurobayes_kppg_qq/step1_k0pipigammaQQ_02.root");
	//QQ->Add("~/Desktop/rootfiles/step1_k0pipigammaQQ_cut02nb.root");
	//QQ->Add("~/Desktop/rootfiles/step1_k0pipigammaQQ_cut02dm.root");

    TChain *sig = new TChain("h1");
	//sig->Add("~/Desktop/rootfiles/mcx3872gamma_cut1.root");
	//sig->Add("~/Desktop/rootfiles/step1_mcx3872gammaee_cut_02.root");
	//sig->Add("~/Desktop/rootfiles/step1_mcx3872gamma_cut02dm.root");
	//sig->Add("~/Desktop/rootfiles/step1_mcx3872gammaee_cut02dm.root");
	sig->Add("~/Desktop/rootfiles/mck0pipigamma_cut.root");
	//sig->Add("~/BAS_2012May/neurobayes_kppg_qq/step1_mck0pipigamma_cut_02.root");
	//sig->Add("~/Desktop/rootfiles/step1_mck0pipigamma_cut02nb.root");
	//sig->Add("~/Desktop/rootfiles/step1_mck0pipigamma_cut02dm.root");

	TChain *jpsi = new TChain("h1");
	jpsi->Add("~/Desktop/rootfiles/k0pipigammaRarec.root");
	jpsi->Add("~/Desktop/rootfiles/k0pipigammaRarem.root");
	//jpsi->Add("~/Desktop/rootfiles/step1_jpsixee_02.root");
	//jpsi->Add("~/Desktop/rootfiles/step1_jpsixee_cut02dm.root");
	//jpsi->Add("~/Desktop/rootfiles/step1_k0pipigammaRare_cutnb.root");

	TChain *data = new TChain("h1");
	data->Add("~/Desktop/rootfiles/k0pipigammaData_cut.root");
	//data->Add("~/Desktop/rootfiles/step1_k0pipigammaData_cutnb.root");

	double bmin=5.2;
	double bmax=5.29;
	//char figtitle [101]= "M_{K#pi#pi} after NB selection (Signal Region)";
	char figtitle [101]= "Modified M_{bc} before NB selection";
	//char figtitle [101]= "#DeltaE after NB selection";
	//char figtitle [101]= "Energy ratio of #pi^{+} & #pi^{-} from B^{0} (CM frame, After NB selection)";
	//char figtitle [101]= "Cosine Angle Between (K^{0}_{S} #pi^{+} #pi^{-}) and K^{0}_{S} (CM Frame)";

	double beamnum=100;
	TH1F* bb1 = new TH1F("bb1",figtitle,beamnum,bmin,bmax);
	TH1F* qq1 = new TH1F("qq1",figtitle,beamnum,bmin,bmax);
	TH1F* sig1 = new TH1F("sig1",figtitle,beamnum,bmin,bmax);
	TH1F* scf1 = new TH1F("scf1",figtitle,beamnum,bmin,bmax);
	TH1F* jpsi1 = new TH1F("jpsi1",figtitle,beamnum,bmin,bmax);
	TH1F* data1 = new TH1F("data1",figtitle,beamnum,bmin,bmax);

	TCanvas *canvas = new TCanvas(figtitle,figtitle,900,600);
  	canvas->cd(1);

  	int sigr=0;
  	//TCut selection = "cosxkcm>0&&kpp0mass<4&&0<=chisqexk&&chisqexk<100&&e9oe25>0.78&&e3>1.5&&sqrt(px1*px1+py1*py1+pz1*pz1)>0.1&&sqrt(px2*px2+py2*py2+pz2*pz2)>0.1&&abs(dz1)<5&&abs(dr1)<2&&abs(dz2)<5&&abs(dr2)<2&&cosllcm>0&&massks0<0.516&&massks0>0.478&&vchisq>=0&&vchisq<40";
  	TCut yeskppg1_1 = "(bpd1==ggmo01&&bpd2==22&&np==1&&nbpd==2)||(bnd1==ggmo01&&bnd2==22&&np==-1&&nbnd==2)";
	TCut yeskppg1_21 = "((abs(bkd1)==311&&abs(bkd2)==211&&bkd3==-bkd2)||(abs(bkd3)==311&&abs(bkd2)==211&&bkd1==-bkd2))&&nkd==3";
	TCut yeskppg1_22 = "((abs(bkd1)==311&&bkd2==113)||(abs(bkd2)==311&&bkd1==113))&&nkd==2&&abs(hi1)==211&&hi2==-hi1";
	TCut yeskppg1 = yeskppg1_1&&(yeskppg1_21||yeskppg1_22);
	
	TCut yeskppg2_1 = "(bpd2==22&&np==1&&nbpd==2)||(bnd2==22&&np==-1&&nbnd==2)";
	TCut yeskppg2_21 = "((abs(bkd1)==323&&abs(bkd2)==211&&ggmo01==bkd1)||(abs(bkd2)==323&&abs(bkd1)==211&&ggmo01==bkd2))&&nkd==2";
	TCut yeskppg2_22 = "((abs(bkd1)==10321&&abs(bkd2)==211&&ggmo01==bkd1)||(abs(bkd2)==10321&&abs(bkd1)==211&&ggmo01==bkd2))&&nkd==2";
	TCut yeskppg2 = yeskppg2_1&&(yeskppg2_21||yeskppg2_22);

  	TCut yeskppg=yeskppg1||yeskppg2;

  	if(sigr==1){
  		BB -> Draw("mbc>>bb1", "mbc>5.27&&mbc<5.29&&de>-0.15&&de<0.1&&costhr>-2");//bg1
		QQ -> Draw("mbc>>qq1", "mbc>5.27&&mbc<5.29&&de>-0.15&&de<0.1&&costhr>-2");//bg2
		sig-> Draw("mbc>>sig1","mbc>5.27&&mbc<5.29&&de>-0.15&&de<0.1&&costhr>-2&&hindex!=0&&dp1lab<0.05&&dp2lab<0.05");//sig &&hindex!=0 //&&hindex!=0&&deltap01<0.05&&deltap02<0.05
		sig-> Draw("mbc>>scf1","mbc>5.27&&mbc<5.29&&de>-0.15&&de<0.1&&costhr>-2&&(hindex==0||dp1lab>0.05||dp2lab>0.05)");
		jpsi->Draw("mbc>>jpsi1",(!yeskppg)&&"mbc>5.27&&mbc<5.29&&de>-0.15&&de<0.1&&costhr>-2");//jpsi
		data->Draw("mbc>>data1", "mbc>5.27&&mbc<5.29&&de>-0.15&&de<0.1&&costhr>-2");//data
  	}else{
		BB -> Draw("mbc>>bb1", "de<0.2&&costhr>-2");//bg1
		QQ -> Draw("mbc>>qq1", "de<0.2&&costhr>-2");//bg2
		sig-> Draw("mbc>>sig1","de<0.2&&costhr>-2&&hindex!=0&&dp1lab<0.05&&dp2lab<0.05");//sig &&hindex!=0 //&&hindex!=0&&deltap01<0.05&&deltap02<0.05
		sig-> Draw("mbc>>scf1","de<0.2&&costhr>-2&&(hindex==0||dp1lab>0.05||dp2lab>0.05)");
		//BB->Draw("mbc>>bb1", "costhr>-2&&((np==1&&bpd1!=443&&bpd2!=443)||(np==-1&&bnd1!=443&&bnd2!=443)||abs(np)!=1)&&bxd1!=443&&bxd2!=443&&mo000!=443&&mo001!=443");//bg1
		jpsi->Draw("mbc>>jpsi1",(!yeskppg)&&"de<0.2&&costhr>-2");//jpsi
		data->Draw("mbc>>data1", "de<0.2&&costhr>-2");//data
  	}
	
	sig1->SetLineColor(kAzure+2);
	bb1->SetLineColor(kPink-1);
	qq1->SetLineColor(kGreen+1);
	scf1->SetLineColor(kOrange-6);
	jpsi1->SetLineColor(kViolet-3);
	data1->SetLineColor(kBlack);

    sig1->SetLineWidth(3);
	bb1->SetLineWidth(3);
	qq1->SetLineWidth(3);
	scf1->SetLineWidth(3);
	jpsi1->SetLineWidth(3);
	data1->SetLineWidth(3);

	//sig1->Scale(771581000*1.99e-5*1e-6);
    sig1->Scale(771581000*1.99e-5*1e-6*10./9);
    bb1->Scale(1./10);
    qq1->Scale(930./929./6);
    //qq1->Scale(1./6);
    scf1->Scale(771581000*1.99e-5*1e-6*10./9);
    jpsi1->Scale(1./50);
    //scf1->Scale(771581000*1.99e-5*1e-6);

    data1->Draw();
    bb1->Draw("same");
    qq1->Draw("same");
    scf1->Draw("same");
    sig1->Draw("same");
    jpsi1->Draw("same");

	//TH1* norm = sig1->DrawNormalized();
	//bb1-> DrawNormalized("same");
	//qq1-> DrawNormalized("same");
	//scf1-> DrawNormalized("same");

	auto legend = new TLegend(0.7,0.75,0.98,0.95); //right legend
	//auto legend = new TLegend(0.12,0.72,0.34,0.9); //left legend
	
	//legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry("sig1","True signal","l");
	legend->AddEntry("bb1","BB background","l");// (without J/#psi)
	legend->AddEntry("qq1","QQ background","l");
	legend->AddEntry("scf1","Self cross feed","l");
	legend->AddEntry("jpsi1","Rare Decay","l");
	legend->AddEntry("data1","Data","l");
	//legend->AddEntry("jpsi1","J/#psi + inclusive","l");
	//legend->AddEntry("sig1","B^{0} #rightarrow X(3872) K^{0}","l");
	//legend->AddEntry("bb1","B^{0} #rightarrow X(3872) K*(892)^{0}","l");
	//legend->AddEntry("qq1","B^{0} #rightarrow X(3872) K^{0} #pi^{0}","l");
	legend->Draw();

	//data1->GetYaxis()->SetTitle("Counts");
	//data1->GetXaxis()->SetTitle("M_{K#pi#pi} (GeV/c^{2})");
	data1->GetXaxis()->SetTitle("M_{bc} (GeV/c^{2})");
	//data1->GetXaxis()->SetTitle("#DeltaE (GeV)");
	data1->SetStats(0);

	char CanvasTitle[200];
  	sprintf(CanvasTitle,"~/Desktop/20180426/BBQQ_kppg_mbc_beforenb.png");
  	canvas->SaveAs(CanvasTitle); 
  	/*
	TCut cut = "costhr>-2&&mbc>5.2";
	TCut cutsig = "hindex!=0&&deltap01<0.05&&deltap02<0.05";
	TH2F* hist = new TH2F("hist", "histogram", 100, 5.2, 5.3, 100, -0.5, 0.5);
	sig -> Project("hist", "de:mbc", cut&&cutsig);
	hist->GetEntries()
	hist->GetCorrelationFactor()
	hist->Draw()
  	*/
} 
