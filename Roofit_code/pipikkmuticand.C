#include "TFile.h"
#include "TTree.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
     
void pipikkmuticand()
{ 

TChain *data = new TChain("h1");
//TChain *qq = new TChain("h1");
//TChain *bg = new TChain("h1");
//TChain *bb = new TChain("h1");

data->Add("qqnvcut.root");

//qq->Add("step2_kfsw_1003_qq_sb_fit_t5st_v3.root");
//qq->Add("step2_kfsw_0903_qq_sb_fit_st2.root");
//qq->Add("step2_kfsw_0903_qq_sb_fit_st3.root");

//bb->Add("step2_kfsw_0903_bb_sb_fit_st1.root");
//bb->Add("step2_kfsw_0903_bb_sb_fit_st2.root");
//bb->Add("genemc_k_1209_pre1.root");

//bg->Add("step2_kfsw_1001_bb_sb_fit_t5st.root");
//bg->Add("step2_kfsw_1001_qq_sb_fit_t5st.root");


float nb_vlike,hi6,hi5,imode,flavor,qr,chix,chiy,chiz,chie,chim,mbc,de,ebeam,etapx,etapy,etapz,etae,etam,kpx,kpy,kpz,ke,km,drk1,dzk1,kpx1,kpy1,kpz1,ke1,km1,drk2,dzk2,ppx,ppy,ppz,pe,pm,drpi,dzpi,ppx1,ppy1,ppz1,pe1,pm1,drpi1,dzpi1,chisqexk,evtcount,bvertexx,bvertexy,bvertexz,overtexx,overtexy,overtexz,exfit,chisqExK,pmiss,emiss,vtflag,deltaz,R2,costhr,costhp,spher,cosb,R2s,R1so,R2so,R3so,R4so,R1gso,R2gso,R3gso,R4gso,R1oo,R2oo,R3oo,R4oo,Expid,Runid,Evtid,eventid,farmid,imm2,mmm2,et,H0oo,H1oo,H2oo,H3oo,H4oo,H0son,H1son,H2son,H3son,H4son,H0soc,H1soc,H2soc,H3soc,H4soc,H0som,H1som,H2som,H3som,H4som,hi1,hi2,hi3,hi4,mo1,gmo1,ggmo1,mo2,gmo2,ggmo2,mo3,gmo3,mo4,gmo4,mo5,gmo5,mo6,gmo6,bnd1,bnd2,bnd3,bnd4,bnd5,bnd6,nbnd ,bpd1,bpd2,bpd3,bpd4,bpd5,bpd6,nbpd,hindex,pxb,pyb,pzb,eb,ptb,gax,gay,gaz,gae,gam,cos1,cos1cm,cosh1,gax1,gay1,gaz1,gae1,gam1,pi0v,pi0p,pi0v1,pi0p1,e9oe25,e9oe251,cos0,cos0cm,cosh0,pidk1,chrgk1,vchisq,ipx,ipy,ipz,charged,mgridp,mgridn,gas;



  int a;

  //int nev = 0;
  int gev = 0;
 

data->SetBranchAddress("gas",&gas);
data->SetBranchAddress("de",&de);
data->SetBranchAddress("mbc",&mbc);
data->SetBranchAddress("ebeam",&ebeam);
data->SetBranchAddress("pxb",&pxb);
data->SetBranchAddress("pyb",&pyb);
data->SetBranchAddress("pzb",&pzb);
data->SetBranchAddress("eb",&eb);
data->SetBranchAddress("ptb",&ptb);
data->SetBranchAddress("chix",&chix);
data->SetBranchAddress("chiy",&chiy);
data->SetBranchAddress("chiz",&chiz);
data->SetBranchAddress("chie",&chie);
data->SetBranchAddress("chim",&chim);
data->SetBranchAddress("gax",&gax);
data->SetBranchAddress("gay",&gay);
data->SetBranchAddress("gaz",&gaz);
data->SetBranchAddress("gae",&gae);
data->SetBranchAddress("gam",&gam);
data->SetBranchAddress("cos1",&cos1);
data->SetBranchAddress("cos1cm",&cos1cm);
data->SetBranchAddress("cosh1",&cosh1);
data->SetBranchAddress("gax1",&gax1);
data->SetBranchAddress("gay1",&gay1);
data->SetBranchAddress("gaz1",&gaz1);
data->SetBranchAddress("gae1",&gae1);
data->SetBranchAddress("gam1",&gam1);
data->SetBranchAddress("pi0v",&pi0v);
data->SetBranchAddress("pi0p",&pi0p);
data->SetBranchAddress("pi0v1",&pi0v1);
data->SetBranchAddress("pi0p1",&pi0p1);
data->SetBranchAddress("e9oe25",&e9oe25);
data->SetBranchAddress("e9oe251",&e9oe251);
data->SetBranchAddress("etapx",&etapx);
data->SetBranchAddress("etapy",&etapy);
data->SetBranchAddress("etapz",&etapz);
data->SetBranchAddress("etae",&etae);
data->SetBranchAddress("etam",&etam);
data->SetBranchAddress("cos0",&cos0);
data->SetBranchAddress("cos0cm",&cos0cm);
data->SetBranchAddress("cosh0",&cosh0);
data->SetBranchAddress("kpx",&kpx);
data->SetBranchAddress("kpy",&kpy);
data->SetBranchAddress("kpz",&kpz);
data->SetBranchAddress("ke",&ke);
//data->SetBranchAddress("kp",&kp);
//data->SetBranchAddress("km",&km);
data->SetBranchAddress("kpx1",&kpx1);
data->SetBranchAddress("kpy1",&kpy1);
data->SetBranchAddress("kpz1",&kpz1);
data->SetBranchAddress("ke1",&ke1);
//data->SetBranchAddress("kp1",&kp1);
//data->SetBranchAddress("km1",&km1);
data->SetBranchAddress("pidk1",&pidk1);
data->SetBranchAddress("chrgk1",&chrgk1);
data->SetBranchAddress("drk1",&drk1);
data->SetBranchAddress("dzk1",&dzk1);
data->SetBranchAddress("drk2",&drk2);
data->SetBranchAddress("dzk2",&dzk2);
data->SetBranchAddress("ppx",&ppx);
data->SetBranchAddress("ppy",&ppy);
data->SetBranchAddress("ppz",&ppz);
data->SetBranchAddress("pe",&pe);
data->SetBranchAddress("pm",&pm);
data->SetBranchAddress("drpi",&drpi);
data->SetBranchAddress("dzpi",&dzpi);
data->SetBranchAddress("ppx1",&ppx1);
data->SetBranchAddress("ppy1",&ppy1);
data->SetBranchAddress("ppz1",&ppz1);
data->SetBranchAddress("pe1",&pe1);
data->SetBranchAddress("pm1",&pm1);
data->SetBranchAddress("drpi1",&drpi1);
data->SetBranchAddress("dzpi1",&dzpi1);
data->SetBranchAddress("vchisq",&vchisq);
data->SetBranchAddress("ipx",&ipx);
data->SetBranchAddress("ipy",&ipy);
data->SetBranchAddress("ipz",&ipz);
data->SetBranchAddress("chisqexk",&chisqexk);
data->SetBranchAddress("evtcount",&evtcount);
data->SetBranchAddress("bvertexx",&bvertexx);
data->SetBranchAddress("bvertexy",&bvertexy);
data->SetBranchAddress("bvertexz",&bvertexz);
data->SetBranchAddress("overtexx",&overtexx);
data->SetBranchAddress("overtexy",&overtexy);
data->SetBranchAddress("overtexz",&overtexz);
data->SetBranchAddress("exfit",&exfit);
data->SetBranchAddress("pmiss",&pmiss);
data->SetBranchAddress("emiss",&emiss);
data->SetBranchAddress("vtflag",&vtflag);
data->SetBranchAddress("deltaz",&deltaz);
data->SetBranchAddress("R2",&R2);
data->SetBranchAddress("costhr",&costhr);
data->SetBranchAddress("costhp",&costhp);
data->SetBranchAddress("spher",&spher);
data->SetBranchAddress("cosb",&cosb);
data->SetBranchAddress("R2s",&R2s);
data->SetBranchAddress("R1so",&R1so);
data->SetBranchAddress("R2so",&R2so);
data->SetBranchAddress("R3so",&R3so);
data->SetBranchAddress("R4so",&R4so);
data->SetBranchAddress("R1gso",&R1gso);
data->SetBranchAddress("R2gso",&R2gso);
data->SetBranchAddress("R3gso",&R3gso);
data->SetBranchAddress("R4gso",&R4gso);
data->SetBranchAddress("R1oo",&R1oo);
data->SetBranchAddress("R2oo",&R2oo);
data->SetBranchAddress("R3oo",&R3oo);
data->SetBranchAddress("R4oo",&R4oo);
data->SetBranchAddress("Evtid",&Evtid);
data->SetBranchAddress("Runid",&Runid);
data->SetBranchAddress("Expid",&Expid);
data->SetBranchAddress("eventid",&eventid);
data->SetBranchAddress("farmid",&farmid);
data->SetBranchAddress("imm2",&imm2);
data->SetBranchAddress("mmm2",&mmm2);
data->SetBranchAddress("et",&et);
data->SetBranchAddress("H0oo",&H0oo);
data->SetBranchAddress("H1oo",&H1oo);
data->SetBranchAddress("H2oo",&H2oo);
data->SetBranchAddress("H3oo",&H3oo);
data->SetBranchAddress("H4oo",&H4oo);
data->SetBranchAddress("H0son",&H0son);
data->SetBranchAddress("H1son",&H1son);
data->SetBranchAddress("H2son",&H2son);
data->SetBranchAddress("H3son",&H3son);
data->SetBranchAddress("H4son",&H4son);
data->SetBranchAddress("H0soc",&H0soc);
data->SetBranchAddress("H1soc",&H1soc);
data->SetBranchAddress("H2soc",&H2soc);
data->SetBranchAddress("H3soc",&H3soc);
data->SetBranchAddress("H4soc",&H4soc);
data->SetBranchAddress("H0som",&H0som);
data->SetBranchAddress("H1som",&H1som);
data->SetBranchAddress("H2som",&H2som);
data->SetBranchAddress("H3som",&H3som);
data->SetBranchAddress("H4som",&H4som);
data->SetBranchAddress("qr",&qr);
data->SetBranchAddress("flavor",&flavor);
data->SetBranchAddress("imode",&imode);
data->SetBranchAddress("hi1",&hi1);
data->SetBranchAddress("hi2",&hi2);
data->SetBranchAddress("hi3",&hi3);
data->SetBranchAddress("hi4",&hi4);
data->SetBranchAddress("hi5",&hi5);
data->SetBranchAddress("hi6",&hi6);
data->SetBranchAddress("mo1",&mo1);
data->SetBranchAddress("gmo1",&gmo1);
data->SetBranchAddress("ggmo1",&ggmo1);
data->SetBranchAddress("mo2",&mo2);
data->SetBranchAddress("gmo2",&gmo2);
data->SetBranchAddress("ggmo2",&ggmo2);
data->SetBranchAddress("mo3",&mo3);
data->SetBranchAddress("gmo3",&gmo3);
data->SetBranchAddress("mo4",&mo4);
data->SetBranchAddress("gmo4",&gmo4);
data->SetBranchAddress("mo5",&mo5);
data->SetBranchAddress("gmo5",&gmo5);
data->SetBranchAddress("mo6",&mo6);
data->SetBranchAddress("gmo6",&gmo6);
data->SetBranchAddress("charged",&charged);
data->SetBranchAddress("bnd1",&bnd1);
data->SetBranchAddress("bnd2",&bnd2);
data->SetBranchAddress("bnd3",&bnd3);
data->SetBranchAddress("bnd4",&bnd4);
data->SetBranchAddress("bnd5",&bnd5);
data->SetBranchAddress("bnd6",&bnd6);
//data->SetBranchAddress("bnd7",&bnd7);
//data->SetBranchAddress("bnd8",&bnd8);
//data->SetBranchAddress("bnd9",&bnd9);
data->SetBranchAddress("nbnd",&nbnd);
data->SetBranchAddress("bpd1",&bpd1);
data->SetBranchAddress("bpd2",&bpd2);
data->SetBranchAddress("bpd3",&bpd3);
data->SetBranchAddress("bpd4",&bpd4);
data->SetBranchAddress("bpd5",&bpd5);
data->SetBranchAddress("bpd6",&bpd6);
//data->SetBranchAddress("bpd7",&bpd7);
//data->SetBranchAddress("bpd8",&bpd8);
//data->SetBranchAddress("bpd9",&bpd9);
data->SetBranchAddress("nbpd",&nbpd);
data->SetBranchAddress("mgridp",&mgridp);
data->SetBranchAddress("mgridn",&mgridn);
data->SetBranchAddress("hindex",&hindex);
data->SetBranchAddress("nb_vlike",&nb_vlike);

  TFile *f1 = new TFile("qqmd.root","recreate");
  TTree *h1 = new TTree("h1", "h1");
h1->Branch("gas",&gas,"gas");
h1->Branch("de",&de,"de");
h1->Branch("mbc",&mbc,"mbc");
h1->Branch("ebeam",&ebeam,"ebeam");
h1->Branch("pxb",&pxb,"pxb");
h1->Branch("pyb",&pyb,"pyb");
h1->Branch("pzb",&pzb,"pzb");
h1->Branch("eb",&eb,"eb");
h1->Branch("ptb",&ptb,"ptb");
h1->Branch("chix",&chix,"chix");
h1->Branch("chiy",&chiy,"chiy");
h1->Branch("chiz",&chiz,"chiz");
h1->Branch("chie",&chie,"chie");
h1->Branch("chim",&chim,"chim");
h1->Branch("gax",&gax,"gax");
h1->Branch("gay",&gay,"gay");
h1->Branch("gaz",&gaz,"gaz");
h1->Branch("gae",&gae,"gae");
h1->Branch("gam",&gam,"gam");
h1->Branch("cos1",&cos1,"cos1");
h1->Branch("cos1cm",&cos1cm,"cos1cm");
h1->Branch("cosh1",&cosh1,"cosh1");
h1->Branch("gax1",&gax1,"gax1");
h1->Branch("gay1",&gay1,"gay1");
h1->Branch("gaz1",&gaz1,"gaz1");
h1->Branch("gae1",&gae1,"gae1");
h1->Branch("gam1",&gam1,"gam1");
h1->Branch("pi0v",&pi0v,"pi0v");
h1->Branch("pi0p",&pi0p,"pi0p");
h1->Branch("pi0v1",&pi0v1,"pi0v1");
h1->Branch("pi0p1",&pi0p1,"pi0p1");
h1->Branch("e9oe25",&e9oe25,"e9oe25");
h1->Branch("e9oe251",&e9oe251,"e9oe251");
h1->Branch("etapx",&etapx,"etapx");
h1->Branch("etapy",&etapy,"etapy");
h1->Branch("etapz",&etapz,"etapz");
h1->Branch("etae",&etae,"etae");
h1->Branch("etam",&etam,"etam");
h1->Branch("cos0",&cos0,"cos0");
h1->Branch("cos0cm",&cos0cm,"cos0cm");
h1->Branch("cosh0",&cosh0,"cosh0");
h1->Branch("kpx",&kpx,"kpx");
h1->Branch("kpy",&kpy,"kpy");
h1->Branch("kpz",&kpz,"kpz");
h1->Branch("ke",&ke,"ke");
//h1->Branch("kp",&kp,"kp");
h1->Branch("km",&km,"km");
h1->Branch("kpx1",&kpx1,"kpx1");
h1->Branch("kpy1",&kpy1,"kpy1");
h1->Branch("kpz1",&kpz1,"kpz1");
h1->Branch("ke1",&ke1,"ke1");
//h1->Branch("kp1",&kp1,"kp1");
h1->Branch("km1",&km1,"km1");
h1->Branch("pidk1",&pidk1,"pidk1");
h1->Branch("chrgk1",&chrgk1,"chrgk1");
h1->Branch("drk1",&drk1,"drk1");
h1->Branch("dzk1",&dzk1,"dzk1");
h1->Branch("drk2",&drk2,"drk2");
h1->Branch("dzk2",&dzk2,"dzk2");
h1->Branch("ppx",&ppx,"ppx");
h1->Branch("ppy",&ppy,"ppy");
h1->Branch("ppz",&ppz,"ppz");
h1->Branch("pe",&pe,"pe");
h1->Branch("pm",&pm,"pm");
h1->Branch("drpi",&drpi,"drpi");
h1->Branch("dzpi",&dzpi,"dzpi");
h1->Branch("ppx1",&ppx1,"ppx1");
h1->Branch("ppy1",&ppy1,"ppy1");
h1->Branch("ppz1",&ppz1,"ppz1");
h1->Branch("pe1",&pe1,"pe1");
h1->Branch("pm1",&pm1,"pm1");
h1->Branch("drpi1",&drpi1,"drpi1");
h1->Branch("dzpi1",&dzpi1,"dzpi1");
h1->Branch("vchisq",&vchisq,"vchisq");
h1->Branch("ipx",&ipx,"ipx");
h1->Branch("ipy",&ipy,"ipy");
h1->Branch("ipz",&ipz,"ipz");
h1->Branch("chisqexk",&chisqexk,"chisqexk");
h1->Branch("evtcount",&evtcount,"evtcount");
h1->Branch("bvertexx",&bvertexx,"bvertexx");
h1->Branch("bvertexy",&bvertexy,"bvertexy");
h1->Branch("bvertexz",&bvertexz,"bvertexz");
h1->Branch("overtexx",&overtexx,"overtexx");
h1->Branch("overtexy",&overtexy,"overtexy");
h1->Branch("overtexz",&overtexz,"overtexz");
h1->Branch("exfit",&exfit,"exfit");
h1->Branch("pmiss",&pmiss,"pmiss");
h1->Branch("emiss",&emiss,"emiss");
h1->Branch("vtflag",&vtflag,"vtflag");
h1->Branch("deltaz",&deltaz,"deltaz");
h1->Branch("R2",&R2,"R2");
h1->Branch("costhr",&costhr,"costhr");
h1->Branch("costhp",&costhp,"costhp");
h1->Branch("spher",&spher,"spher");
h1->Branch("cosb",&cosb,"cosb");
h1->Branch("R2s",&R2s,"R2s");
h1->Branch("R1so",&R1so,"R1so");
h1->Branch("R2so",&R2so,"R2so");
h1->Branch("R3so",&R3so,"R3so");
h1->Branch("R4so",&R4so,"R4so");
h1->Branch("R1gso",&R1gso,"R1gso");
h1->Branch("R2gso",&R2gso,"R2gso");
h1->Branch("R3gso",&R3gso,"R3gso");
h1->Branch("R4gso",&R4gso,"R4gso");
h1->Branch("R1oo",&R1oo,"R1oo");
h1->Branch("R2oo",&R2oo,"R2oo");
h1->Branch("R3oo",&R3oo,"R3oo");
h1->Branch("R4oo",&R4oo,"R4oo");
h1->Branch("Evtid",&Evtid,"Evtid");
h1->Branch("Runid",&Runid,"Runid");
h1->Branch("Expid",&Expid,"Expid");
h1->Branch("eventid",&eventid,"eventid");
h1->Branch("farmid",&farmid,"farmid");
h1->Branch("imm2",&imm2,"imm2");
h1->Branch("mmm2",&mmm2,"mmm2");
h1->Branch("et",&et,"et");
h1->Branch("H0oo",&H0oo,"H0oo");
h1->Branch("H1oo",&H1oo,"H1oo");
h1->Branch("H2oo",&H2oo,"H2oo");
h1->Branch("H3oo",&H3oo,"H3oo");
h1->Branch("H4oo",&H4oo,"H4oo");
h1->Branch("H0son",&H0son,"H0son");
h1->Branch("H1son",&H1son,"H1son");
h1->Branch("H2son",&H2son,"H2son");
h1->Branch("H3son",&H3son,"H3son");
h1->Branch("H4son",&H4son,"H4son");
h1->Branch("H0soc",&H0soc,"H0soc");
h1->Branch("H1soc",&H1soc,"H1soc");
h1->Branch("H2soc",&H2soc,"H2soc");
h1->Branch("H3soc",&H3soc,"H3soc");
h1->Branch("H4soc",&H4soc,"H4soc");
h1->Branch("H0som",&H0som,"H0som");
h1->Branch("H1som",&H1som,"H1som");
h1->Branch("H2som",&H2som,"H2som");
h1->Branch("H3som",&H3som,"H3som");
h1->Branch("H4som",&H4som,"H4som");
h1->Branch("qr",&qr,"qr");
h1->Branch("flavor",&flavor,"flavor");
h1->Branch("imode",&imode,"imode");
h1->Branch("hi1",&hi1,"hi1");
h1->Branch("hi2",&hi2,"hi2");
h1->Branch("hi3",&hi3,"hi3");
h1->Branch("hi4",&hi4,"hi4");
h1->Branch("hi5",&hi5,"hi5");
h1->Branch("hi6",&hi6,"hi6");
h1->Branch("mo1",&mo1,"mo1");
h1->Branch("gmo1",&gmo1,"gmo1");
h1->Branch("ggmo1",&ggmo1,"ggmo1");
h1->Branch("mo2",&mo2,"mo2");
h1->Branch("gmo2",&gmo2,"gmo2");
h1->Branch("ggmo2",&ggmo2,"ggmo2");
h1->Branch("mo3",&mo3,"mo3");
h1->Branch("gmo3",&gmo3,"gmo3");
h1->Branch("mo4",&mo4,"mo4");
h1->Branch("gmo4",&gmo4,"gmo4");
h1->Branch("mo5",&mo5,"mo5");
h1->Branch("gmo5",&gmo5,"gmo5");
h1->Branch("mo6",&mo6,"mo6");
h1->Branch("gmo6",&gmo6,"gmo6");
h1->Branch("charged",&charged,"charged");
h1->Branch("bnd1",&bnd1,"bnd1");
h1->Branch("bnd2",&bnd2,"bnd2");
h1->Branch("bnd3",&bnd3,"bnd3");
h1->Branch("bnd4",&bnd4,"bnd4");
h1->Branch("bnd5",&bnd5,"bnd5");
h1->Branch("bnd6",&bnd6,"bnd6");
//h1->Branch("bnd7",&bnd7,"bnd7");
//h1->Branch("bnd8",&bnd8,"bnd8");
//h1->Branch("bnd9",&bnd9,"bnd9");
h1->Branch("nbnd",&nbnd,"nbnd");
h1->Branch("bpd1",&bpd1,"bpd1");
h1->Branch("bpd2",&bpd2,"bpd2");
h1->Branch("bpd3",&bpd3,"bpd3");
h1->Branch("bpd4",&bpd4,"bpd4");
h1->Branch("bpd5",&bpd5,"bpd5");
h1->Branch("bpd6",&bpd6,"bpd6");
//h1->Branch("bpd7",&bpd7,"bpd7");
//h1->Branch("bpd8",&bpd8,"bpd8");
//h1->Branch("bpd9",&bpd9,"bpd9");
h1->Branch("nbpd",&nbpd,"nbpd");
h1->Branch("mgridp",&mgridp,"mgridp");
h1->Branch("mgridn",&mgridn,"mgridn");
h1->Branch("hindex",&hindex,"hindex");
h1->Branch("nb_vlike",&nb_vlike,"nb_vlike");


int nevs = 0;
        int nev = 0;
        int evtid = 0;
        float chisq = 0.0;
      	data->GetEntry(0, 0);
        h1->Fill();
        while (data->GetEntry(nev, 0))
        {
                if (Evtid == evtid)
                {
                        if (chisqexk < chisq)
                        {
                                nevs = nev;
                                chisq = chisqexk;
                        }
                }
                else
                {
                        if (nevs == 0)
                        {
                                nevs = 1;
                        }
                        else
                        {
                                data->GetEntry(nevs, 0);
                                h1->Fill();
                                data->GetEntry(nev, 0);
                        }
                        nevs = nev;
                        evtid = Evtid;
                        chisq = chisqexk;
                }
                nev++;
        }
        h1->Write();

        cout << "# of writen events: " << h1->GetEntries() << endl;


return;


}
