
// FILE : chi_c0_eta_pipipipigg.cc
// Author : Tzuying
// Date : 2009 Nov. 29 (charged particle polar angle cos_theta)
//        2009 Dec. 21 (deltaZ & sigMC gamma energy)
//        2009 DEC. 25 (mod. gamma energy_default value & nbnd nbpd for pipi KK & hindex for pipi KK)
//        2010 Jan. 13 (add chisq for ExKfitter)
//        2010 Jan. 26 (correct hindex for kk pipi mode & modify ExKfitter vertexing quality chisqExK & B0 B0Bar ptype change)
//	  2010 Feb. 12 (correct Beam energy)
//        2010 Mar. 17 (pipi kk correct loop)
//	  2010 Jun. 14 (add kid eff.) & correct the idhep of KK mode
//	  2010 Oct. 31 kid tight cut K>0.6 pi<0.4 && add kid calibration for KK mode
//	  2011 Jan. 16 TOF is used.

// Description: for 4s B0 to h anti-h mode,from jelov
//
//

#include "belle.h"
#include "particle/Particle.h"
#include "particle/PID.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "benergy/BeamEnergy.h"
#include "shape.h"
#include "findLambda.h"
#include "UserUtility.h"
#include "toolbox/FoxWolfr.h"
#include "toolbox/Thrust.h"
#include "brutus/brutus_f.h"
#include "hamlet/Hamlet.h"	//tagging
#include "hamlet/AnaBrecon.h" // to use AnaBrecon class ***
#include "tagv/TagV.h"
#include "ExKFitter.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "kid_eff_06.h"  // sysmetic error
#include "kid/atc_pid.h" // sysmetic error
#include "eid/eid.h"
#include "mdst/Muid_mdst.h"
#include "mdst/mdst.h"
#include HEPEVT_H
#include EVTCLS_H
#include MDST_H
#include BELLETDF_H
#include "mdst/findKs.h"
#include "kfitter/kvertexfitter.h" // to use vertexfitter.
#include "kfitter/kmassfitter.h"   // to use massfitter
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/khelix2xyz.h"
#include "kfitter/kfitterparticle.h"
#include "helix/Helix.h"
#include "k_sfw.h"
#include "CLHEP/Vector/LorentzVector.h" // /afs/afs11/belle/belle/b20090127_0910/src/util/belleCLHEP/Vector
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <iostream>
#include "userinfo.h"
#include "panther/panther.h"
#include "ip/IpProfile.h"
//#include "./mdst3.icc"

#include MDST_H
#include BELLETDF_H
#include HEPEVT_H

using namespace std;
//the unit is GeV
//#define E_HER 7.998213 // high energy ring i.e electron ring
//#define E_LER 3.499218    // low energy ring i.e. positron ring

#define BTUPLE

#if defined(BELLE_NAMESPACE)
namespace Belle
{
#endif
//utilities
int checkMultiUse(Particle &i, Particle &j);
///////////問學長這裡需要改class的名稱嗎？
// Module class
class chi_c0_eta_pipipipigg : public Module
{ //start class chi_c0_eta_pipipipigg
  public:
	chi_c0_eta_pipipipigg(void){};
	~chi_c0_eta_pipipipigg(void){};
	void init(int *);
	void term(void);
	void disp_stat(const char *){};
	void hist_def(void);
	void event(BelleEvent *, int *);
	void begin_run(BelleEvent *, int *);
	void end_run(BelleEvent *, int *){};
	void other(int *, BelleEvent *, int *){};
	//void shape(Particle &, float &, float &, float &, float &, float par_sfw[17], Vector4 &);
	void shape(Particle &, float &, float &, float &, float &, float par_sfw[17], Vector4 &, HepPoint3D &, int &, HepPoint3D &);
	void GetImpactParameters(const Mdst_charged *, double *, double *, int); //double*傳址呼叫,  double&傳參考呼叫
	double deltaZ(Particle, double &, double &, double &, double &);
	///////////這裡有需要改嗎？

  public:
  private:
	BelleHistogram *m_pi0, *m_b, *m_b_beam;
	BelleHistogram *e_gamm;
	BelleHistogram *dec_lenx, *dec_alenx;
	BelleHistogram *dec_leny, *dec_aleny;
	BelleHistogram *m_lamv, *m_lamppi, *m_lamkf;
	BelleHistogram *p_elec, *p_muon;
	BelleTuple *b_tpl;
	int evtcount;
	brutus_f Fisher_fw;
	brutus_f Fisher_ksfw[7];
	KID_eff_06 eff_s1_kpi_kp;
	KID_eff_06 eff_s1_kpi_km;
	KID_eff_06 eff_s1_kpi_kp_fake;
	KID_eff_06 eff_s1_kpi_km_fake;
	KID_eff_06 eff_s1_kpi_pip;
	KID_eff_06 eff_s1_kpi_pim;
	KID_eff_06 eff_s1_kpi_pip_fake;
	KID_eff_06 eff_s1_kpi_pim_fake;
	KID_eff_06 eff_s1_pipi_pip;
	KID_eff_06 eff_s1_pipi_pim;
	KID_eff_06 eff_s1_pipi_pip_fake;
	KID_eff_06 eff_s1_pipi_pim_fake;
	KID_eff_06 eff_s1_kk_kp;
	KID_eff_06 eff_s1_kk_km;
	KID_eff_06 eff_s1_kk_kp_fake;
	KID_eff_06 eff_s1_kk_km_fake;

	KID_eff_06 eff_s2_kpi_kp;
	KID_eff_06 eff_s2_kpi_km;
	KID_eff_06 eff_s2_kpi_kp_fake;
	KID_eff_06 eff_s2_kpi_km_fake;
	KID_eff_06 eff_s2_kpi_pip;
	KID_eff_06 eff_s2_kpi_pim;
	KID_eff_06 eff_s2_kpi_pip_fake;
	KID_eff_06 eff_s2_kpi_pim_fake;
	KID_eff_06 eff_s2_pipi_pip;
	KID_eff_06 eff_s2_pipi_pim;
	KID_eff_06 eff_s2_pipi_pip_fake;
	KID_eff_06 eff_s2_pipi_pim_fake;
	KID_eff_06 eff_s2_kk_kp;
	KID_eff_06 eff_s2_kk_km;
	KID_eff_06 eff_s2_kk_kp_fake;
	KID_eff_06 eff_s2_kk_km_fake;

}; //end class chi_c0_eta_pipipipigg
extern "C" Module_descr *mdcl_chi_c0_eta_pipipipigg()
{
	chi_c0_eta_pipipipigg *module = new chi_c0_eta_pipipipigg;
	Module_descr *dscr = new Module_descr("chi_c0_eta_pipipipigg", module);
	IpProfile::define_global(dscr);
	BeamEnergy::define_global(dscr);
	return dscr;
} //end extern "C" Module_descr
void chi_c0_eta_pipipipigg::hist_def(void)
{ //start void chi_c0_eta_pipipipigg::hist_def(void)
	extern BelleTupleManager *BASF_Histogram;
	BelleTupleManager &tm = *BASF_Histogram;
#ifdef BTUPLE

	//ntuple name length--> 8 characters(upper limit)//////宣告變數//////有需要改的嗎？****
	b_tpl = BASF_Histogram->ntuple("B0 --> hh", "btag ggmo1 ggmo2 Etar2He Etar1He exk_eta vchisq vchisq_2 Mb_c dE qr flavor imode bnd1 bnd2 bnd3 bnd4 bnd5 bnd6 bnd7 bnd8 bnd9 nbnd charged enerbpd1 enerbpd2 enerbpd3 enerbpd4 enerbpd5 enerbpd6 enerbpd7 enerbpd8 enerbpd9 bvertexx bvertexy bvertexz overtexx overtexy overtexz exfit chisqExK vtflag deltaz px2 py2 pz2 e2 mass2 pid2 chrg2 dr2 dz2 cos_2polar pkid2 ppiid2 pid3 chrg3 dr3 dz3 cos_3polar pkid3 ppiid3 pid4 chrg4 dr4 dz4 icos_4polar pkid4 ppiid4  pid5 chrg5 dr5 dz5 cos_5polar pkid5 ppiid5  pid6 chrg6 dr6 dz6 cos_6polar pkid6 ppiid6 ipx ipy ipz R2 costhr costhp spher cosb R2s R1so R2so R3so R4so R1gso R2gso R3gso R4gso R1oo R2oo R3oo R4oo evtflag EVTID RUNID EXPID eventid farmid ebeam mode hi1 hi2 hi3 hi4 hi5 hi6 mo1 gmo1 mo2 gmo2 mo3 gmo3 mo4 gmo4 mo5 gmo5 mo6 gmo6 hindex hrepeat pmiss emiss et imm2 mmm2 H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son H0soc H1soc H2soc H3soc H4soc H0som H1som H2som H3som H4som CHI0 chi_c0px chi_c0py chi_c0pz chi_c0e chi_c0mass ETA etapx etapy etapz etae etamass gapx gapy gapz gae gamass gapx1 gapy1 gapz1 gae1 gamass1 gapx1 gapy1 gapz1 gae1 gamass1 etapx etapy etapz etae etamass Kpx Kpy Kpz Ke Kp Kmass Kpx1 Kpy1 Kpz1 Ke1 Kp1	Kmass1 ppx ppy ppz pe pp pmass p1px1 p1py1 p1pz1 p1e1 p1p1 p1mass1 ");
#endif

	e_gamm = BASF_Histogram->histogram("Energy of Photon", 100, 0.0, 3.5);
	m_pi0 = BASF_Histogram->histogram("Mass(Pi0 - GeV)", 100, 0.115, 0.145);
	m_b = BASF_Histogram->histogram("Mass(B - GeV)", 50, 5.2, 5.3);
	m_b_beam = BASF_Histogram->histogram("Massfit(B - GeV)", 50, 5.2, 5.3);
	dec_lenx = BASF_Histogram->histogram("LAM decay length x cm", 40, 0., 20.);
	dec_leny = BASF_Histogram->histogram("LAM decay length y cm", 40, 0., 20.);
	dec_alenx = BASF_Histogram->histogram("ALAM decay length x cm", 40, 0., 20.);
	dec_aleny = BASF_Histogram->histogram("ALAM decay length y cm", 40, 0., 20.);
	m_lamv = BASF_Histogram->histogram("Mass(LamdaVee2- GeV)", 100, 1.05, 1.15);
	m_lamppi = BASF_Histogram->histogram("Mass(Lamda ppi- GeV)", 100, 1.05, 1.15);
	m_lamkf = BASF_Histogram->histogram("Mass(Lamda fit- GeV)", 100, 1.05, 1.15);
	Fisher_fw.histogram(10, "Super Fox-Wolfram");
	k_sfw::initialize(Fisher_ksfw);
	p_elec = BASF_Histogram->histogram("P_electron", 50, 0.0, 5.0);
	p_muon = BASF_Histogram->histogram("P_muon", 50, 0.0, 5.0);
} //void chi_c0_eta_pipipipigg::hist_def(void)

////////////////要不要改class
void chi_c0_eta_pipipipigg::init(int *)
{ //start void chi_c0_eta_pipipipigg::init(int*)
	evtcount = 0;

	Hamlet::init();
	////////////////////////////////////這裡在算系統誤差，一開始做分析不用理他，也可以另外做
	// 1st argument : probability cut value (0.1,0.2,...,0.9)
	// 2nd argument : 1 = prob(K/pi)>0.X for kaon (kaon eff)
	//                2 = prob(K/pi)>0.X for pion (pion fake)
	//                3 = prob(pi/K)>0.X for pion (pion eff)
	//                4 = prob(pi/K)>0.X for kaon (kaon fake)
	//                where 0.X is the number given as 1st argument
	// 3rd argument : name (anything is OK. but it should be different  with each other)

	// 1st argument : probability cut value (0.1,0.2,...,0.9)
	// 2nd argument : 1 = prob(K/pi)>0.X for kaon (kaon eff)
	//                2 = prob(K/pi)>0.X for kaon candidate (pion fake to Kaon)
	//                3 = prob(pi/K)>0.X for pion (pion eff)
	//                4 = prob(pi/K)>0.X for pion candidate (kaon fake to pion)
	//                where 0.X is the number given as 1st argument
	// 3rd argument : name (anything is OK. but it should be different  with each other)

	/////////////////////////////這裡也是在做系統誤差嗎？可以刪掉嗎？////////////////////
	eff_s1_kpi_kp.init(.6, 1, "track_s1kpi_kp", "kideff-2006-svd1-pos.dat");		   //K plus eff.
	eff_s1_kpi_km.init(.6, 1, "track_s1kpi_km", "kideff-2006-svd1-neg.dat");		   //K minus eff.
	eff_s1_kpi_kp_fake.init(.6, 2, "track_s1kpi_kp_fake", "kideff-2006-svd1-pos.dat"); //K fake plus eff.
	eff_s1_kpi_km_fake.init(.6, 2, "track_s1kpi_km_fake", "kideff-2006-svd1-neg.dat"); //K fake minus eff.

	eff_s1_kpi_pip.init(.6, 3, "track_s1kpi_pip", "kideff-2006-svd1-pos.dat");			 //K plus eff.
	eff_s1_kpi_pim.init(.6, 3, "track_s1kpi_pim", "kideff-2006-svd1-neg.dat");			 //K minus eff.
	eff_s1_kpi_pip_fake.init(.6, 4, "track_s1kpi_pip_fake", "kideff-2006-svd1-pos.dat"); //K fake plus eff.
	eff_s1_kpi_pim_fake.init(.6, 4, "track_s1kpi_pim_fake", "kideff-2006-svd1-neg.dat"); //K fake minus eff.

	// just test for comparing the previous study
	//    eff_s1_kpi_kp.init( .6, 1, "track_s1kpi_kp", "kideff-2005-svd1-pos.dat" );//K plus eff.
	//    eff_s1_kpi_km.init( .6, 1, "track_s1kpi_km", "kideff-2005-svd1-neg.dat" );//K minus eff.
	//    eff_s1_kpi_kp_fake.init( .6, 2, "track_s1kpi_kp_fake", "kideff-2005-svd1-pos.dat" );//K fake plus eff.
	//    eff_s1_kpi_km_fake.init( .6, 2, "track_s1kpi_km_fake", "kideff-2005-svd1-neg.dat" );//K fake minus eff.

	//    eff_s1_kpi_pip.init( .6, 3, "track_s1kpi_pip", "kideff-2006-svd1-pos.dat" );//K plus eff.
	//    eff_s1_kpi_pim.init( .6, 3, "track_s1kpi_pim", "kideff-2006-svd1-neg.dat" );//K minus eff.
	//    eff_s1_kpi_pip_fake.init( .6, 4, "track_s1kpi_pip_fake", "kideff-2006-svd1-pos.dat" );//K fake plus eff.
	//    eff_s1_kpi_pim_fake.init( .6, 4, "track_s1kpi_pim_fake", "kideff-2006-svd1-neg.dat" );//K fake minus eff.

	eff_s1_pipi_pip.init(.6, 3, "track_s1pipi_pip", "kideff-2006-svd1-pos.dat");		   //pi plus eff.
	eff_s1_pipi_pip_fake.init(.6, 4, "track_s1pipi_pip_fake", "kideff-2006-svd1-pos.dat"); //pi fake plus eff.
	eff_s1_pipi_pim.init(.6, 3, "track_s1pipi_pim", "kideff-2006-svd1-neg.dat");		   //pi minus eff.
	eff_s1_pipi_pim_fake.init(.6, 4, "track_s1pipi_pim_fake", "kideff-2006-svd1-neg.dat"); //pi fake minus eff.

	// just test for comparing the previous study
	//    eff_s1_pipi_pip.init( .6, 3, "track_s1pipi_pip", "kideff-2005-svd1-pos.dat" );//pi plus eff.
	//    eff_s1_pipi_pip_fake.init( .6, 4, "track_s1pipi_pip_fake", "kideff-2005-svd1-pos.dat" );//pi fake plus eff.
	//    eff_s1_pipi_pim.init( .6, 3, "track_s1pipi_pim", "kideff-2005-svd1-neg.dat" );//pi minus eff.
	//    eff_s1_pipi_pim_fake.init( .6, 4, "track_s1pipi_pim_fake", "kideff-2005-svd1-neg.dat" );//pi fake minus eff.

	eff_s1_kk_kp.init(.6, 1, "track_s1kk_kp", "kideff-2006-svd1-pos.dat"); //k plus eff.
	eff_s1_kk_km.init(.6, 1, "track_s1kk_km", "kideff-2006-svd1-neg.dat"); //k minus eff.

	eff_s1_kk_kp_fake.init(.6, 2, "track_s1kk_kp_fake", "kideff-2006-svd1-pos.dat"); //k fake plus eff.
	eff_s1_kk_km_fake.init(.6, 2, "track_s1kk_km_fake", "kideff-2006-svd1-neg.dat"); //k fake minus eff.

	// just test for comparing the previous study
	//    eff_s1_kk_kp.init( .6, 1, "track_s1kk_kp", "kideff-2005-svd1-pos.dat" );//k plus eff.
	//    eff_s1_kk_km.init( .6, 1, "track_s1kk_km", "kideff-2005-svd1-neg.dat" );//k minus eff.

	//    eff_s1_kk_kp_fake.init( .6, 2, "track_s1kk_kp_fake", "kideff-2005-svd1-pos.dat" );//k fake plus eff.
	//    eff_s1_kk_km_fake.init( .6, 2, "track_s1kk_km_fake", "kideff-2005-svd1-neg.dat" );//k fake minus eff.

	eff_s2_kpi_kp.init(.6, 1, "track_s2kpi_kp", "kideff-2010-svd2-pos.dat");		   //K plus eff.
	eff_s2_kpi_km.init(.6, 1, "track_s2kpi_km", "kideff-2010-svd2-neg.dat");		   //K minus eff.
	eff_s2_kpi_kp_fake.init(.6, 2, "track_s2kpi_kp_fake", "kideff-2010-svd2-pos.dat"); //K fake plus eff.
	eff_s2_kpi_km_fake.init(.6, 2, "track_s2kpi_km_fake", "kideff-2010-svd2-neg.dat"); //K fake minus eff.

	eff_s2_kpi_pip.init(.6, 3, "track_s2kpi_pip", "kideff-2010-svd2-pos.dat");			 //K plus eff.
	eff_s2_kpi_pim.init(.6, 3, "track_s2kpi_pim", "kideff-2010-svd2-neg.dat");			 //K minus eff.
	eff_s2_kpi_pip_fake.init(.6, 4, "track_s2kpi_pip_fake", "kideff-2010-svd2-pos.dat"); //K fake plus eff.
	eff_s2_kpi_pim_fake.init(.6, 4, "track_s2kpi_pim_fake", "kideff-2010-svd2-neg.dat"); //K fake minus eff.

	eff_s2_pipi_pip.init(.6, 3, "track_s2pipi_pip", "kideff-2010-svd2-pos.dat");		   //pi plus eff.
	eff_s2_pipi_pip_fake.init(.6, 4, "track_s2pipi_pip_fake", "kideff-2010-svd2-pos.dat"); //pi fake plus eff.
	eff_s2_pipi_pim.init(.6, 3, "track_s2pipi_pim", "kideff-2010-svd2-neg.dat");		   //pi minus eff.
	eff_s2_pipi_pim_fake.init(.6, 4, "track_s2pipi_pim_fake", "kideff-2010-svd2-neg.dat"); //pi fake minus eff.

	eff_s2_kk_kp.init(.6, 1, "track_s2kk_kp", "kideff-2010-svd2-pos.dat"); //k plus eff.
	eff_s2_kk_km.init(.6, 1, "track_s2kk_km", "kideff-2010-svd2-neg.dat"); //k minus eff.

	eff_s2_kk_kp_fake.init(.6, 2, "track_s2kk_kp_fake", "kideff-2010-svd2-pos.dat"); //k fake plus eff.
	eff_s2_kk_km_fake.init(.6, 2, "track_s2kk_km_fake", "kideff-2010-svd2-neg.dat"); //k fake minus eff.

} //end void chi_c0_eta_pipipipigg::init(int*)

////////////////////////////////////設這個物件要幹嘛呢？把他//會不會怎麼樣？////////////////////////////////////
void chi_c0_eta_pipipipigg::term(void)
{   //start void chi_c0_eta_pipipipigg::term(void)
	/*
		   eff_s2_kpi_kp.calculate();
		   eff_s2_kpi_km.calculate();
		   eff_s2_kpi_kp_fake.calculate();
		   eff_s2_kpi_km_fake.calculate();
		   eff_s2_kpi_pip.calculate();
		   eff_s2_kpi_pim.calculate();
		   eff_s2_kpi_pip_fake.calculate();
		   eff_s2_kpi_pim_fake.calculate();
		   eff_s2_pipi_pip.calculate();
		   eff_s2_pipi_pim.calculate();
		   eff_s2_pipi_pip_fake.calculate();
		   eff_s2_pipi_pim_fake.calculate();
		   eff_s2_kk_kp.calculate();
		   eff_s2_kk_km.calculate();
		   eff_s2_kk_kp_fake.calculate();
		   eff_s2_kk_km_fake.calculate();


		   eff_s1_kpi_kp.calculate();
		   eff_s1_kpi_km.calculate();
		   eff_s1_kpi_kp_fake.calculate();
		   eff_s1_kpi_km_fake.calculate();
		   eff_s1_kpi_pip.calculate();
		   eff_s1_kpi_pim.calculate();
		   eff_s1_kpi_pip_fake.calculate();
		   eff_s1_kpi_pim_fake.calculate();
		   eff_s1_pipi_pip.calculate();
		   eff_s1_pipi_pim.calculate();
		   eff_s1_pipi_pip_fake.calculate();
		   eff_s1_pipi_pim_fake.calculate();
		   eff_s1_kk_kp.calculate();
		   eff_s1_kk_km.calculate();
		   eff_s1_kk_kp_fake.calculate();
		   eff_s1_kk_km_fake.calculate();


		   eff_s2_kpi_kp.dump();
		   eff_s2_kpi_km.dump();
		   eff_s2_kpi_kp_fake.dump();
		   eff_s2_kpi_km_fake.dump();
		   eff_s2_kpi_pip.dump();
		   eff_s2_kpi_pim.dump();
		   eff_s2_kpi_pip_fake.dump();
		   eff_s2_kpi_pim_fake.dump();
		   eff_s2_pipi_pip.dump();
		   eff_s2_pipi_pim.dump();
		   eff_s2_pipi_pip_fake.dump();
		   eff_s2_pipi_pim_fake.dump();
		   eff_s2_kk_kp.dump();
		   eff_s2_kk_km.dump();
		   eff_s2_kk_kp_fake.dump();
		   eff_s2_kk_km_fake.dump();


		   eff_s1_kpi_kp.dump();
		   eff_s1_kpi_km.dump();
		   eff_s1_kpi_kp_fake.dump();
		   eff_s1_kpi_km_fake.dump();
		   eff_s1_kpi_pip.dump();
		   eff_s1_kpi_pim.dump();
		   eff_s1_kpi_pip_fake.dump();
		   eff_s1_kpi_pim_fake.dump();
		   eff_s1_pipi_pip.dump();
		   eff_s1_pipi_pim.dump();
		   eff_s1_pipi_pip_fake.dump();
		   eff_s1_pipi_pim_fake.dump();
		   eff_s1_kk_kp.dump();
		   eff_s1_kk_km.dump();
		   eff_s1_kk_kp_fake.dump();
		   eff_s1_kk_km_fake.dump();

		*/
} //end void chi_c0_eta_pipipipigg::term(void)

// Set IP location   ////////////// HepPoint3D是class嗎？為什麼不用::////////////////////////
HepPoint3D IP(0, 0, 0);
HepSymMatrix IPerr(3, 0);

void chi_c0_eta_pipipipigg::begin_run(BelleEvent *evptr, int *status)
{					  //start void chi_c0_eta_pipipipigg::begin_run
	eid::init_data(); // available in the new lib (after b199907*)

	//beam energy
	BeamEnergy::begin_run();

	// Get IP profile data from $BELLE_POSTGRES_SERVER
	IpProfile::begin_run();

	// Dump IP profile data to STDOUT (optional)
	IpProfile::dump();

	// Set IP and error
	IP = IpProfile::position();
	IPerr = IpProfile::position_err();

	//tagging
	// To initialize LH tables by EvtGen MC
	Hamlet::begin_run(Hamlet::MULT_DIM_LH);
} //end void chi_c0_eta_pipipipigg::begin_run

///////////////////////////重點開始！////////
///////////////////////////////////一開始的class需要換嗎？//////////////////////////
void chi_c0_eta_pipipipigg::event(BelleEvent *evptr, int *status) //～～～～～～～～為什麼逗號前面的不用宣告呢～～～～～～～～～～
{																		   //start void chi_c0_eta_pipipipigg::event(BelleEvent *evptr, int *status)
	const HepPoint3D &ip = IpProfile::position();						   ////////什麼意思呢？
	Belle_event_Manager &bevt_mgr = Belle_event_Manager::get_manager();
	Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
	Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
	Mdst_pi0_Manager &pi0_mag = Mdst_pi0_Manager::get_manager(); //~~~~~~~這裏需要加我們的particle進去嗎？~~~~~~~~~
													Mdst_vee2_Manager &Vee2Mgr = Mdst_vee2_Manager::get_manager(); //~~~~~~~get_manager是什麼意思呢？~~~~~~~~~~~
	Mdst_event_add_Manager &mevtmgr = Mdst_event_add_Manager::get_manager();
	Mdst_vee_daughters_Manager &veedmgr = Mdst_vee_daughters_Manager::get_manager();
	Mdst_klm_mu_ex_Manager &klmmgr = Mdst_klm_mu_ex_Manager::get_manager();
	Mdst_gamma_Manager &gamma_Mgr = Mdst_gamma_Manager::get_manager();
	Mdst_ecl_Manager &ecl_mag = Mdst_ecl_Manager::get_manager();
	Mdst_ecl_trk_Manager &ecltrk_mag = Mdst_ecl_trk_Manager::get_manager();

	Mdst_ecl_aux_Manager &eclaux_mag = Mdst_ecl_aux_Manager::get_manager();
	Evtcls_hadronic_flag_Manager &evtcls = Evtcls_hadronic_flag_Manager::get_manager();
	Gen_hepevt_Manager &gen_mgr = Gen_hepevt_Manager::get_manager();
	//////////////////////////////////////////////這些綠掉的不用理他嗎？
	//for reprocessed data
	//  remove_duplicates();
	//scale_momenta(1.00246, 1.0);
	//for reprocessed exp7 data
	//scale_momenta(1.00221, 1.0);
	//for reprocessed exp9 data
	//scale_momenta(1.00149, 1.0);

	//cout << "void all mdst" << endl;

	evtcount++;
	//cout << "event = " << evtcount << endl;

	double E_HER = BeamEnergy::E_HER();
	double E_LER = BeamEnergy::E_LER();
	double cross_angle = BeamEnergy::Cross_angle(); //radian
	static Vector4 P_BEAM(-E_HER * sin(cross_angle), 0.0, E_LER - E_HER * cos(cross_angle), E_HER + E_LER);
	//static 靜態局部變數，讓變數執行玩函式後會保存下載，下次繼續使用
	int Run = 0;
	int Evt = 0;
	int Eventid;
	int Farmid;
	int Exp;

	Run = bevt_mgr[0].RunNo();
	Evt = bevt_mgr[0].EvtNo();
	Exp = bevt_mgr[0].ExpNo();
	Eventid = (bevt_mgr[0].EvtNo() & 0x0FFFFFFF);
	Farmid = bevt_mgr[0].EvtNo() >> 28;

	//cout << "static all var" << endl;

	// set IP and error
	int IPUsable = 0;
	if (IpProfile::usable())
	{ //482
		IP = IpProfile::position(1);
		IPerr = IpProfile::position_err_b_life_smeared(1);
		IPUsable = 1;
	}
	else
	{
		IP = HepPoint3D(0, 0, 0);
		IPerr = HepSymMatrix(3, 0);

		//cout << "set IP and error" << endl;
	}

	//=============================================
	// Various cut values declared as const double.
	//=============================================
	// for J/Psi, pi0, Xc1, Xc2, J/Psi(2s) mass
	const double MPSI_PDG = 3.0969160;			 // J/Psi mass in 2009 PDG average.
	const double META_PDG = 0.5478530;			 // eta mass in 2009 PDG average.
	const double METAP_PDG = 0.95778;			 // eta' mass in 2009 PDG average.
	const double MPI0_PDG = 0.1349766;			 // pi0 mass in 2009 PDG average.
	const double MRHO_PDG = 0.77549;			 // rho mass in 2009 PDG average.
	const double M_ee_ub = 0.036;				 // "Exclusive" HP recommendateion.    3.132916
	const double M_ee_lb = 0.15;				 // After including radiated photons.  2.946916
	const double M_mumu_ub = 0.036;				 // "Exclusive" HP recommendateion.  3.132916
	const double M_mumu_lb = 0.06;				 // "Exclusive" HP recommendateion.  3.036916
	const double M_pi0_min = 0.134 - 0.0054 * 3; //pi0 mass, 3 sigma      0.1178
	const double M_pi0_max = 0.134 + 0.0054 * 3; //pi0 mass, 3 sigma      0.1502
	const double MKs_PDG = 0.497614;			 //Ks mass in 2010 PDG average.
   	const double MCHIC0_PDG = 3.4150;												 //==============
												 //宣告你的變數，自己定義，可自行增加
												 //===============================
	//list of Particle
	std::vector<Particle> gamma;
	std::vector<Particle> pi_0;
	std::vector<Particle> pi_plus;
	std::vector<Particle> pi_minus;
	std::vector<Particle> k_plus;
	std::vector<Particle> k_minus;
	std::vector<Particle> k_short;
	std::vector<Particle> p_plus;
	std::vector<Particle> p_minus;
	//std::vector<Particle> Lamda;
	//std::vector<Particle> Lamdabar;
	//std::vector<Particle> Sigma_plus;
	std::vector<Particle> mu_plus;
	std::vector<Particle> mu_minus;
	std::vector<Particle> e_plus;
	std::vector<Particle> e_minus;
	std::vector<Particle> tk_plus;
	std::vector<Particle> tk_minus;
	std::vector<Particle> ttk_plus;
	std::vector<Particle> ttk_minus;
	std::vector<Particle> tpi_plus;
	std::vector<Particle> tpi_minus;
	std::vector<Particle> B_cand1;
	//std::vector<Particle> B_cand2;
	//std::vector<Particle> B_cand3;
	//std::vector<Particle> B_cand4;
	//std::vector<Particle> B_cand5;
	std::vector<Particle> eta; //自己增加加定義
	std::vector<Particle> chi_c0;
	//std::vector<Particle> chi_c1;
	//std::vector<Particle> chi_c2;
	//std::vector<Particle> Jpsi;
	std::vector<Particle> gamma1;
	std::vector<Particle> gammass;

	//cout << "void your own var" << endl;

	//儲存特定particle資訊，ＩＤ標準質量
	//Particle Type(/sw/belle/belle/b20040727_1143/share/data-files/qq98/decay.dec)  寫法去這裡找
	Ptype ptype_gamma("GAMM");
	Ptype ptype_pi_0("PI0");
	Ptype ptype_pi_plus("PI+");
	Ptype ptype_pi_minus("PI-");
	Ptype ptype_k_plus("K+");
	Ptype ptype_k_minus("K-");
	//	Ptype ptype_k_short("K0S");
	Ptype ptype_p_plus("P+");
	Ptype ptype_p_minus("AP+");
	//Ptype ptype_Lamda("LAM");
	Ptype ptype_Lamdabar("ALAM");
	//Ptype ptype_Sigma_plus("SIG+");
	Ptype ptype_mu_plus("MU+");
	Ptype ptype_mu_minus("MU-");
	Ptype ptype_e_plus("E+");
	Ptype ptype_e_minus("E-");
	//Ptype ptype_Bu_cand("B+");
	Ptype ptype_Bd_cand("B0");
	//Ptype ptype_B0B_cand("B0B");
	//Ptype ptype_Bs_cand("BS0");
	Ptype ptype_eta_cand("ETA");
	Ptype ptype_chi_c0_cand("CHI0");
	//Ptype ptype_chi_c1_cand("CHI1");
	//Ptype ptype_chi_c2_cand("CHI2");
	//Ptype ptype_jpsi("PSI");

	//cout << "store special particle info" << endl;

	gamma.clear();
	pi_0.clear();
	pi_plus.clear();
	pi_minus.clear();
	k_plus.clear();
	k_minus.clear();
	//k_short.clear();
	//Lamda.clear();
	//Lamdabar.clear();
	e_plus.clear();
	e_minus.clear();
	mu_plus.clear();
	mu_minus.clear();
	tpi_plus.clear();
	tpi_minus.clear();
	tk_plus.clear();
	tk_minus.clear();
	ttk_plus.clear();
	ttk_minus.clear();
	B_cand1.clear();
	eta.clear(); //自己增加加定義
	chi_c0.clear();
	gamma1.clear();
	gammass.clear();
	/*atc_pid (particle identification)
		  1. include "kid/atc_pid.h"
		  2.atc_pid  selkpi(accq,tofq,cdcq,ids,idb)
		  3.particle=(0,1,2,3,4)=(e,mu,pi,k,proton)accq:
		  3 -->Probability is calculated based on measured Npe and PDF.
		  This option gives higher eff. than accq=0, especially in the high momentum region.
tofq: 1 -->Used with checking Z hit position w.r.t the track extrapolation.
0 -->Used without checking Z hit position w.r.t the track extrapolation.
cdcq: 5 -->Correct dE/dx run-dependent shift in 2001 summer reprocesing (exp11,13).
Also correct shift in the high momentum pion (all exp.No.).
*/

	//cout << "act info" << endl;

	//  atc_pid selKpi(0,1,0,3,2);//for the reprocessed exp7 data
	atc_pid selkpi(3, 1, 5, 3, 2);
	atc_pid selkp(3, 1, 5, 3, 4);
	atc_pid selppi(3, 1, 5, 4, 2);
	atc_pid selpk(3, 1, 5, 4, 3);
	//  atc_pid selmuk(3,1,5,1,3);
	//  atc_pid selmupi(3,1,5,1,2);

	// distinguish MC from data,MC is MCstatus == 1
	int MCstatus = 0;

	//cout << "MCstatus=0" << endl;

	for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin(); //i 本身就是物件
		 i != gen_mgr.end(); i++)
	{ 
		MCstatus = 1;
		//	cout << "MCstatus=1" << endl;
	}
	//fill pi or k lists from MDST_Charged Data Base
	for (std::vector<Mdst_charged>::iterator i = charged_mag.begin(); i != charged_mag.end(); i++)

	//cout<<"55555" <<endl;

	{
		Muid_mdst muon(*i);
		Mdst_charged &ch = *i;
		int mu_sta = muon.Status();
		int outcome = muon.Outcome();
		int mu_level = muon.Level();
		int reject = muon.Prerejection();
		double mu_like = muon.Muon_likelihood();
		double chi2 = muon.Chi_2();
		//for eid
		eid sel_e(ch);
		//float eid_prob = sel_e.prob(0,-1,0);
		float eid_prob = sel_e.prob(3, -1, 5);
		//Select k/pi candidates

		//cout << "select k/pi cand" << endl;

		if ((reject || mu_like < 0.9) && eid_prob < 0.9)
		{							   //603
			if (selkpi.prob(*i) > 0.6) //select k pi
									   //if(selkpi.prob(*i)>0.)
			{						   //606
				//cout << "sel k/pi ing" << endl;
				if ((*i).charge() > 0)
				{								 //608
					if (ch.trk().quality() == 0) // Obtain "Good track"好的定義為0
					{							 //610
						//cout << "push back k+" << endl;
						Particle tmp(*i, ptype_k_plus);
						k_plus.push_back(tmp);

					} //610
				}	 //608
				else
				{								 //616
					if (ch.trk().quality() == 0) // Obtain "Good track"
					{							 //618
						Particle tmp(*i, ptype_k_minus);
						k_minus.push_back(tmp);
						//cout << "push back k-" << endl;
					}						//618
				}							//616
			}								//606
			else if (selkpi.prob(*i) < 0.4) //prob 機率 跟likelyhood 類似
											//if (selkpi.prob(*i)<1.0)
			{
				if ((*i).charge() > 0)
				{								 //628
					if (ch.trk().quality() == 0) // Obtain "Good track"
					{							 //630
						Particle tmp(*i, ptype_pi_plus);
						pi_plus.push_back(tmp);
						//cout << "push back pi+" << endl;
					} //630
				}	 //628
				else
				{								 //636
					if (ch.trk().quality() == 0) // Obtain "Good track"
					{							 //638
						Particle tmp(*i, ptype_pi_minus);
						pi_minus.push_back(tmp);
						//cout << "push back pi-" << endl;
					} //638
				}	 //636
			}
		}

		//Select mu+/mu- candidates

		//cout<<" select mu" <<endl;

		else if (!reject && mu_like > 0.8)
		{
			if ((*i).charge() > 0)
			{
				//cout << "selet mu ing" << endl;
				Particle tmp(*i, ptype_mu_plus);
				mu_plus.push_back(tmp);
				//cout << "push back mu+" << endl;
			}
			else
			{
				Particle tmp(*i, ptype_mu_minus);
				mu_minus.push_back(tmp);
				//cout << "push back mu-" << endl;
			}
		}

		//std::cout<<"end select mu"<<endl;

		//cout<<"select electron"<<endl;
		//Select electron
		else if (eid_prob > 0.6) //每部都改成else if 不然會重複選取
		{
			if ((*i).charge() > 0)
			{
				Particle tmp(*i, ptype_e_plus);
				e_plus.push_back(tmp);
			}
			else
			{
				Particle tmp(*i, ptype_e_minus);
				e_minus.push_back(tmp);
			}
		}

		//Select porton
		else
		{
			if ((*i).charge() > 0)
			{
				//cout << "sel porton ing" << endl;
				Particle tmp(*i, ptype_p_plus);
				p_plus.push_back(tmp);
				//cout << "push back p+" << endl;
			}
			else
			{
				Particle tmp(*i, ptype_p_minus);
				p_minus.push_back(tmp);
				//cout << "push back p-" << endl;
			}

		} //681
		  //=======================
		  //以上選完charge particle
		  //=======================

		//end of filling the charged  particles
	}
	double alldr = 1;
	double alldz = 5;
	for (std::vector<Particle>::iterator i = pi_plus.begin(); i != pi_plus.end(); i++)
	{
		double dr, dz;
		const Mdst_charged *ch1 = &(*i).mdstCharged();
		GetImpactParameters(ch1, &dr, &dz, 2);
		if (abs(dr) < alldr && abs(dz) < alldz)
		{
			Particle tmp(*i);
			tpi_plus.push_back(tmp);
		}
	}
	for (std::vector<Particle>::iterator i = pi_minus.begin(); i != pi_minus.end(); i++)
	{
		double dr, dz;
		const Mdst_charged *ch1 = &(*i).mdstCharged();
		GetImpactParameters(ch1, &dr, &dz, 2);
		if (abs(dr) < alldr && abs(dz) < alldz)
		{
			Particle tmp(*i);
			tpi_minus.push_back(tmp);
		}
	}
	for (std::vector<Particle>::iterator i = k_plus.begin(); i != k_plus.end(); i++)
	{
		double dr, dz;
		const Mdst_charged *ch1 = &(*i).mdstCharged();
		GetImpactParameters(ch1, &dr, &dz, 3);
		if (abs(dr) < alldr && abs(dz) < alldz)
		{
			Particle tmp(*i);
			tk_plus.push_back(tmp);
		}
	}
	for (std::vector<Particle>::iterator i = k_minus.begin(); i != k_minus.end(); i++)
	{
		double dr, dz;
		const Mdst_charged *ch1 = &(*i).mdstCharged();
		GetImpactParameters(ch1, &dr, &dz, 3);
		if (abs(dr) < alldr && abs(dz) < alldz)
		{
			Particle tmp(*i);
			tk_minus.push_back(tmp);
		}
	}
	// fill Ks list from Mdst_vee2 bank
	for (std::vector<Mdst_vee2>::iterator i = Vee2Mgr.begin(); i != Vee2Mgr.end(); i++)
	{
		int goodVee;
		FindKs findks;
		if ((*i).kind() == 1) // 1 for Ks, 2 for Lambda, 3 for anti-Lambda , 4 for gamma -> ee
		{					 
			Particle tmp(*i);
			k_short.push_back(tmp);
		}
	}	 
	
	//==================
	//select gamma/pi_0
	//==================

	//cout << "start select gamma/pi0" << endl;

	//gamma.clear();
	for (std::vector<Mdst_gamma>::iterator it = gamma_Mgr.begin(); it != gamma_Mgr.end(); it++)
	{
		//cout << "start gamma" << endl;
		int flags1 = -1;
		int flags2 = -1;
		int flags3 = -1;
		int gammamatch = (*it).ecl().match(); //宣告gammamatch  去ecl的match函數找

		Mdst_ecl_aux_Manager &eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
		Mdst_ecl_aux &aux = eclaux_mgr(Panther_ID((*it).ecl().get_ID()));

		double e9oe25 = aux.e9oe25(); //energy 以中心3*3(個晶格所量到的能量) /5*5 (個晶格所量到的能量),ecl是以晶格組成的量能器，
		//e9oe25m越接近1，帶表能量集中

		//	cout << "gamma selection debug" << endl;
		//cout << "(*it)= " << (*it) << endl;

		for (std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s != eclaux_mag.end(); s++)
		{

			Mdst_ecl_aux &ch = *s;
			for (std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p != ecl_mag.end(); p++)
			{

				if ((*s).get_ID() == (*p).get_ID())
				{
					//forward e>0.1
					if ((*p).theta() < 0.547 && (*p).energy() > 0.1)
						flags1 = 1;
					//forward e>0.05
					if ((*p).theta() > 0.562 && (*p).theta() < 2.246 && (*p).energy() > 0.05)
						flags2 = 1;
					//foreard e>0.15
					if ((*p).theta() > 2.246 && (*p).energy() > 0.15)
						flags3 = 1;
				}
			}

		} //end ecl

		//cout << "end ecl" << endl;

		if ((flags1 || flags2 || flags3) && gammamatch == 0) //多加一個matchgamma＝0代表沒有track
		{
			Particle tmp(*it);
			gamma1.push_back(tmp);
		}
	}
	//cout << "end gammal selection" << endl;

	for (std::vector<Mdst_pi0>::iterator i = pi0_mag.begin(); i != pi0_mag.end(); i++)
	{
		Particle pi0(*i);
		//std::cout << "start pi0" << std::endl;
		int pi_n;
		double e1, e2, a, p1, p2, p, cos;
		e1 = (pi0).child(0).e();
		e2 = (pi0).child(1).e();
		a = abs(e1 - e2) / (e1 + e2);
		/////////////////////////////////

		Mdst_pi0 pi0yo = (pi0).mdstPi0();

		double pi0mass = pi0yo.mass();
		if (e1 > 0.05 && e2 > 0.05 && pi0mass > 0.1178 && pi0mass < 0.1502 && a < 0.9) //pi0 select mass
		{																			
			int flags1 = -1;
			int flags2 = -1;
			int flags3 = -1;

			for (std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s != eclaux_mag.end(); s++)
			{ 
				//cout << "ecl begin" << endl;
				Mdst_ecl_aux &ch = *s;
				for (std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p != ecl_mag.end(); p++)
				{
					if ((*s).get_ID() == (*p).get_ID())
					{ 
						//forward e>0.1
						if ((*p).theta() < 0.547 && (*p).energy() > 0.1)
							flags1 = 1;
						//barrel e>0.05
						if ((*p).theta() > 0.562 && (*p).theta() < 2.246 && (*p).energy() > 0.05)
							flags2 = 1;
						//backward e>0.15
						if ((*p).theta() > 2.246 && (*p).energy() > 0.15)
							flags3 = 1;
					} 

				} 
				//cout << "end ecl begin" << endl;
			} //end ecl
			if (flags1)
			{ 
				Particle tmp(*i);
				pi_0.push_back(tmp);
			} 
			if (flags2)
			{ 
				Particle tmp(*i);
				pi_0.push_back(tmp);
			} 
			if (flags3)
			{ 
				Particle tmp(*i);
				pi_0.push_back(tmp);
			} 
		}
		//end Gamma selection
		//b_tpl->dumpData();
		//*status = 1;
	} //end pi_0 
	//cout << "end pi0 selection " << endl;
	for (std::vector<Particle>::iterator i2 = gamma1.begin(); i2 != gamma1.end(); i2++)
	{

		int flag = -1;
		for (std::vector<Particle>::iterator i = pi_0.begin(); i != pi_0.end(); i++)
		{

			for (int k = 0; k < 2; k++)
			{
				Particle child = (*i).child(k);
				if ((*i2).mdstGamma().ecl().get_ID() == child.mdstGamma().ecl().get_ID())
				{
					flag = 1;
				}
			}
		}
		//cout << "end pi0 begin" << endl;
		if (flag == 1)
		{

			continue;
		}
		else
		{
			Particle tmp(*i2);
			gammass.push_back(tmp);
			//cout << "end gamma begin" << endl;
		}
	}
	//========================
	//reconstruction chic0mass
	//========================
	HepLorentzVector chic_pipi_lv;
	for (std::vector<Particle>::iterator i = tpi_plus.begin(); i != tpi_plus.end(); i++)
	{
		Particle &pi1_plus_chic = *i; //

		for (std::vector<Particle>::iterator i2 = i + 1; i2 != tpi_plus.end(); i2++)
		{
			Particle &pi2_plus_chic = *i2;
			for (std::vector<Particle>::iterator i3 = tpi_minus.begin(); i3 != tpi_minus.end(); i3++)
			{
				Particle &pi1_minus_chic = *i3;
				for (std::vector<Particle>::iterator i4 = i3 + 1; i4 != tpi_minus.end(); i4++)
				{
					Particle &pi2_minus_chic = *i4;
					chic_pipi_lv = pi1_plus_chic.p() + pi2_plus_chic.p() + pi1_minus_chic.p() + pi2_minus_chic.p();

					Particle particle_chic_pipi(chic_pipi_lv, ptype_chi_c0_cand);
					particle_chic_pipi.relation().append(pi1_plus_chic);
					particle_chic_pipi.relation().append(pi2_plus_chic);
					particle_chic_pipi.relation().append(pi1_minus_chic);
					particle_chic_pipi.relation().append(pi2_minus_chic);
					chi_c0.push_back(particle_chic_pipi);
				}
			}
		}
	}
	//cout << "start combination" << endl;

	combination(eta, ptype_eta_cand, gammass, gammass,0.0973,0.0527);
	combination(B_cand1, ptype_Bd_cand, chi_c0, eta);

	//cout << "bcandsize" << B_cand1.size() << endl;
	//cout << "etasize" << eta.size() << endl;
	//cout << "chic0size" << chi_c0.size() << endl;

	//(*i)chi_c0(0):chi_c0
	//(*i)chi_c0(0).child(0):pi_plus
	int btag =0;
	int run = 0;
	float vchisq2;
	float vchisq;
	float gae;
	float gae1;
	double lastchic0mass=1;
	double lastetamass = 0.1;
	double lastchoose=500;
	for (std::vector<Particle>::iterator i = B_cand1.begin(); i != B_cand1.end(); i++)
	{
		HepLorentzVector b_p((*i).p());
		HepLorentzVector boost_vector(-E_HER * sin(cross_angle), 0.0, E_LER - E_HER * cos(cross_angle), E_HER + E_LER);
		b_p.boost(boost_vector.boostVector());
		double ebeam = BeamEnergy::E_beam_corr();
		double mass_sqr = ebeam * ebeam - b_p.vect().mag2();			   //Ｂmass 質量平方
		double mass = (mass_sqr > 0.) ? sqrt(mass_sqr) : -sqrt(-mass_sqr); // true  --> sqrt(mass_sqr)
		float de = b_p.e() - ebeam;										   //delta E 理論上是0
		double etamass = (*i).child(1).p().mag();
		double chi_c0mass = (*i).child(0).p().mag();
		double gammaenergy1=(*i).child(1).child(0).e();
        double gammaenergy2=(*i).child(1).child(1).e();
		double choose = abs(chi_c0mass-MCHIC0_PDG)+abs(etamass - META_PDG);
		run += 1;
		kvertexfitter kvf;
		addTrack2fit(kvf, (*i).child(0).child(2)); // add "pi+" to kfitter
		addTrack2fit(kvf, (*i).child(0).child(3)); // add "pi-" to kfitter
		addTrack2fit(kvf, (*i).child(0).child(0));
		addTrack2fit(kvf, (*i).child(0).child(1));
		unsigned err = kvf.fit(); // do "fitting"
		//float vchisq;
		if (err == 0)
		{
			vchisq = kvf.chisq();
		}
		else
		{
			vchisq = -1.;
		}
		//cout<<"vchisq= "<<vchisq<<endl;
		if (((*i).child(1).child(0).px() == (*i).child(1).child(1).px()) && ((*i).child(1).child(0).py() == (*i).child(1).child(1).py()) &&
			((*i).child(1).child(0).pz() == (*i).child(1).child(1).pz()))
		{
			continue;
		}
		else
		{	
			//if(gammaenergy1 > 0 && gammaenergy1 > 0.05 && gammaenergy2 > 0 && gammaenergy2 > 0.05)
			//{
    				if(vchisq > 0 && vchisq <= 100)
  				{
       					if(mass >= 5.2 && mass <=5.29)
        				{
            					if(abs(de)<0.2)
            					{
                					if(etamass > 0.49 && etamass < 0.58)
                					{
                    						if(chi_c0mass > 3.2 && chi_c0mass < 3.7)
                    						{
                        						if(choose < lastchoose)
                        						{
					                            		lastchoose = choose;
                            							lastchic0mass = chi_c0mass;
         					                   			lastetamass = etamass;
  					                          			btag = run;
                            							vchisq = vchisq2;
                					            		gae = gammaenergy1;
                				           				gae1 = gammaenergy2;
                        						}
                        						else
                        						{
                            							continue;
                        						}
                    						}
                					}
            					}
        				}	
    				}
			//}	
		}
	}
	//	cout<<"end first B loop"<<endl;
	int runi = 0;
	//cout << "start scend B combination" << endl;
	//B_cand1 -> k pi
	for (std::vector<Particle>::iterator i = B_cand1.begin(); i != B_cand1.end(); i++)
	{ 
		runi += 1;
		//cout<<"input bug0 info"<<endl;
		HepLorentzVector b_p((*i).p());
		//cout<<"input bug1 info"<<endl;
		//boost_vector 很重要，對撞的四動量(Px,Py,Px,|P|),固定Ｙ方向是0              右邊beam的能量和左邊beam的能量 400多行有定義
		HepLorentzVector boost_vector(-E_HER * sin(cross_angle), 0.0, E_LER - E_HER * cos(cross_angle), E_HER + E_LER);
		//cout<<"input bug2 info"<<endl;
		b_p.boost(boost_vector.boostVector());
		//cout<<"input bug3 info"<<endl;
		//直接boost到 upsilon 4s 的mass fram(CM fram) boost就是單純座標轉換
		//轉換到beam上，動量空間係，upsilon4s是靜止的，從這裡bb飛出去
		//平常抓到得值都是實驗室座標系得直，轉換到beam的坐標系，兩邊的beam合起來，站在beam上看，這樣看某西些特性看得比較清楚
		//double ebeam = Benergy();
		double ebeam = BeamEnergy::E_beam_corr();
		//cout<<"input bug4 info"<<endl;
		double mass_sqr = ebeam * ebeam - b_p.vect().mag2();			   //Ｂmass 質量平方
		//cout<<"input bug5 info"<<endl;
		double mass = (mass_sqr > 0.) ? sqrt(mass_sqr) : -sqrt(-mass_sqr); // true  --> sqrt(mass_sqr)
		//這就是Mbc B meson constrain mass // false --> -sqrt(-mass_sqr)
		//HepLorentzVector etainfo1((*i).child(1).child(0));
		//etainfo1.boost(boost_vector.boostVector());
		//HepLorentzVector etainfo2((*i).child(1).child(1));
		//etainfo2.boost(boost_vector.boostVector());
		//cout<<"input bug6 info"<<endl;
		float de = b_p.e() - ebeam; //delta E 理論上是0
		//cout<<"input bug7 info"<<endl;
		if (btag == runi)
		{
			//#ifdef BTUPLE //這行可有可無
			b_tpl->column("vchisq", vchisq2);//上面宣告過的變數，和算出來的值  
                        //cout<<"input vchisq info"<<endl;
			b_tpl->column("Mb_c", mass);
			//cout<<"input B mass info"<<endl;
			b_tpl->column("dE", de);
			//cout<<"input de info"<<endl;
			b_tpl->column("ebeam", ebeam);
			//cout<<"input ebeam info"<<endl;
			b_tpl->column("btag", btag);
			//cout<<"input btag info"<<endl;
			//cout <<"bcm:"<<b_mass<<" de:"<<de<<endl;
			////fill ntuple of chi_c0
			b_tpl->column("chi_c0px", (*i).child(0).px());
			b_tpl->column("chi_c0py", (*i).child(0).py());
			b_tpl->column("chi_c0pz", (*i).child(0).pz());
			b_tpl->column("chi_c0e", (*i).child(0).e());
			b_tpl->column("chi_c0mass", (*i).child(0).p().mag());
			//cout<<"input chic0 info"<<endl;
			////fill ntuple of gamma
			b_tpl->column("gapx", (*i).child(1).child(0).px());
			b_tpl->column("gapy", (*i).child(1).child(0).py());
			b_tpl->column("gapz", (*i).child(1).child(0).pz());
			b_tpl->column("gae",(*i).child(1).child(0).e());
			b_tpl->column("gamass", (*i).child(1).child(0).p().mag());
			//分辨gamma&&gamma1
			b_tpl->column("gapx1", (*i).child(1).child(1).px());
			b_tpl->column("gapy1", (*i).child(1).child(1).py());
			b_tpl->column("gapz1", (*i).child(1).child(1).pz());
			//b_tpl->column("gae1",gae1);
			b_tpl->column("gae1",(*i).child(1).child(1).e());
			b_tpl->column("gamass1", (*i).child(1).child(1).p().mag());
			////fill ntuple of eta
			b_tpl->column("etapx", (*i).child(1).px());
			b_tpl->column("etapy", (*i).child(1).py());
			b_tpl->column("etapz", (*i).child(1).pz());
			b_tpl->column("etae", (*i).child(1).e());
			b_tpl->column("etamass", (*i).child(1).p().mag());
			//cout<<"input eta info"<<endl;
			//fill ntuple of K+
		//	b_tpl->column("Kpx", (*i).child(0).child(2).px());
		//	b_tpl->column("Kpy", (*i).child(0).child(2).py());
		//	b_tpl->column("Kpz", (*i).child(0).child(2).pz());
		//	b_tpl->column("Ke", (*i).child(0).child(2).e());
			//b_tpl->column("Kp", (*i).child(1).child(2).mag());//four vector transfers to 3D vector , and mag is magtitude//總動量
		//	b_tpl->column("Kmass1", (*i).child(0).child(2).p().mag());
			//cout<<"input k+ info"<<endl;
			//fill ntuple of K-
		//	b_tpl->column("Kpx1", (*i).child(0).child(3).px());
		//	b_tpl->column("Kpy1", (*i).child(0).child(3).py());
		//	b_tpl->column("Kpz1", (*i).child(0).child(3).pz());
		//	b_tpl->column("Ke1", (*i).child(0).child(3).e());
			//b_tpl->column("Kp1", (*i).child(1).child(3).mag());// four vector transfers to 3D vector , and mag is magtitude
		//	b_tpl->column("Kmass1", (*i).child(0).child(3).p().mag());
			//cout<<"input k- info"<<endl;
			//fill ntuple of pi+
		//b_tpl->column("ppx", (*i).child(0).child(0).px());
		//	b_tpl->column("ppy", (*i).child(0).child(0).py());
		//	b_tpl->column("ppz", (*i).child(0).child(0).pz());
		//	b_tpl->column("pe", (*i).child(0).child(0).e());
			//b_tpl->column("pp", (*i).child(1).child(0).mag());// four vector transfers to 3D vector , and mag is magtitude
		//	b_tpl->column("pmass", (*i).child(0).child(0).p().mag());
			//cout<<"input pi+ info"<<endl;
			//----------------------------------- information of pi- from B (*i).child(0).child(1) -------------------------
			//fill ntuple of pi-
			//b_tpl->column("p1px1", (*i).child(0).child(1).px());
			//b_tpl->column("p1py1", (*i).child(0).child(1).py());
			//b_tpl->column("p1pz1", (*i).child(0).child(1).pz());
			//b_tpl->column("p1e1", (*i).child(0).child(1).e());
			//b_tpl->column("p1p1", (*i).child(1).child(1).mag());//four vector transfers to 3D vector , and mag is magtitude
			//b_tpl->column("p1mass1", (*i).child(0).child(1).p().mag());
			//cout<<"input pi- info"<<endl;
			//cout<<"save all particle info"<<endl;
/*
			const Mdst_charged *chk3 = &(*i).child(0).child(0).mdstCharged(); //pi+info
			float k3pi_id = selkpi.prob(chk3);
			//float k3pi_id_pk = selpk.prob(chk3);
			//float k3pi_id_ppi = selppi.prob(chk3);
			b_tpl->column("pid3", k3pi_id);
			//b_tpl->column("pkid1", k3pi_id_pk);
			//b_tpl->column("ppiid", k1pi_id_ppi);
			float chrg3 = chk3->charge();
			b_tpl->column("chrg3", chrg3);
			double dr3, dz3;						  //0,1,2,3,4 e mu pi k p
			GetImpactParameters(chk3, &dr3, &dz3, 2); //算dr dz的函數
			b_tpl->column("dr3", float(dr3));
			b_tpl->column("dz3", float(dz3));
			//polar angle for charge particle
			double cos_polar_3;
			Vector3 ch3_momentum((*i).child(0).child(0).px(), (*i).child(0).child(0).py(), (*i).child(0).child(0).pz());
			Vector3 P_BEAM_polar(E_HER * sin(cross_angle), 0.0, -E_LER + E_HER * cos(cross_angle));
			cos_polar_3 = ch3_momentum.dot(P_BEAM_polar);
			cos_polar_3 = cos_polar_3 / (ch3_momentum.mag() * P_BEAM_polar.mag());
			b_tpl->column("cos_3polar", float(cos_polar_3));
			//end of pi+(-)
			//b_tpl->column("p2",(*i).child(1).vect().mag());
			//b_tpl->column("mass2", (*i).child(1).momentum().p().mag());
			const Mdst_charged *chk4 = &(*i).child(0).child(1).mdstCharged(); //pi- info
			float k4pi_id = selkpi.prob(chk4);
			float k4pi_id_pk = selpk.prob(chk4);
			//float k4pi_id_ppi = selppi.prob(chk4);
			b_tpl->column("pkid4", k4pi_id_pk);
			//b_tpl->column("ppiid2", k2pi_id_ppi);
			float chrg4 = chk4->charge();
			b_tpl->column("chrg4", chrg4);
			//cout<<"start GetImpactParameters"<<endl;
			double dr4, dz4;
			GetImpactParameters(chk4, &dr4, &dz4, 2); //
			b_tpl->column("dr4", float(dr4));
			b_tpl->column("dz4", float(dz4));
			//polar angle for charge particle
			double cos_polar_4;
			Vector3 ch4_momentum((*i).child(0).child(1).px(), (*i).child(0).child(1).py(), (*i).child(0).child(1).pz());
			cos_polar_4 = ch4_momentum.dot(P_BEAM_polar);
			cos_polar_4 = cos_polar_4 / (ch4_momentum.mag() * P_BEAM_polar.mag());
			b_tpl->column("cos_4polar", cos_polar_4);
			//end of pi-(+)

			const Mdst_charged *chk5 = &(*i).child(0).child(2).mdstCharged(); //k+ info
			float k5pi_id = selkpi.prob(chk5);
			//float k5pi_id_pk = selpk.prob(chk5);
			//float k5pi_id_ppi = selppi.prob(chk5);
			b_tpl->column("pkid5", k5pi_id);
			//b_tpl->column("ppiid5", k5pi_id_ppi);
			float chrg5 = chk5->charge();
			b_tpl->column("chrg5", chrg5);
			//cout<<"start GetImpactParameters"<<endl;
			double dr5, dz5;
			GetImpactParameters(chk5, &dr5, &dz5, 3); //
			b_tpl->column("dr5", float(dr5));
			b_tpl->column("dz5", float(dz5));
			//polar angle for charge particle
			double cos_polar_5;
			Vector3 ch5_momentum((*i).child(0).child(2).px(), (*i).child(0).child(2).py(), (*i).child(0).child(2).pz());
			cos_polar_5 = ch5_momentum.dot(P_BEAM_polar);
			cos_polar_5 = cos_polar_5 / (ch5_momentum.mag() * P_BEAM_polar.mag());
			b_tpl->column("cos_5polar", cos_polar_5);
			//end of k+(-)

			const Mdst_charged *chk6 = &(*i).child(0).child(3).mdstCharged(); //k- info
			float k6pi_id = selkpi.prob(chk6);
			//float k6pi_id_pk = selpk.prob(chk6);
			//float k6pi_id_ppi = selppi.prob(chk6);
			b_tpl->column("pkid6", k6pi_id);
			//b_tpl->column("ppiid6", k6pi_id_ppi);
			float chrg6 = chk6->charge();
			b_tpl->column("chrg6", chrg6);
			//cout<<"start GetImpactParameters"<<endl;
			double dr6, dz6;
			GetImpactParameters(chk6, &dr6, &dz6, 3); //
			b_tpl->column("dr6", float(dr6));
			b_tpl->column("dz6", float(dz6));
			//polar angle for charge particle
			double cos_polar_6;
			Vector3 ch6_momentum((*i).child(0).child(3).px(), (*i).child(0).child(3).py(), (*i).child(0).child(3).pz());
			cos_polar_6 = ch6_momentum.dot(P_BEAM_polar);
			cos_polar_6 = cos_polar_6 / (ch6_momentum.mag() * P_BEAM_polar.mag());
			b_tpl->column("cos_6polar", cos_polar_6);
			//end of k-(+)
*/
			//cout << "input ip" << endl;
			b_tpl->column("ipx", IP.x());
			b_tpl->column("ipy", IP.y());
			b_tpl->column("ipz", IP.z());
			//cout << "end input ip" << endl;
			//#endif
			//for shape variables
			Vector4 OtherBVector;
			float R2, spher, cos_thr, cos_thp, par_sfw[13];
			//Particle &ch3 = (*i).child(0).child(0);
			//Particle &ch4 = (*i).child(0).child(1);
			//Particle &ch5 = (*i).child(0).child(2);
			//Particle &ch6 = (*i).child(0).child(3);
			Particle &ch7 = (*i).child(1).child(0);
			Particle &ch8 = (*i).child(1).child(1);
			// for kid calibration!   k pi
			//double plabchild3 = ch3.ptot();
			//double costhechild3 = ch3.momentum().p().pz() / ch3.ptot();

			//double plabchild4 = ch4.ptot();
			//double costhechild4 = ch4.momentum().p().pz() / ch4.ptot();

			//double plabchild5 = ch5.ptot();
			//double costhechild5 = ch5.momentum().p().pz() / ch5.ptot();

			//double plabchild6 = ch6.ptot();
			//double costhechild6 = ch6.momentum().p().pz() / ch6.ptot();
			//=========================================
			// Vertex Constrain
			//=========================================
			//Mdst_charged Mdst_ch3 = ch3.mdstCharged();
			//Mdst_charged Mdst_ch4 = ch4.mdstCharged();
			//Mdst_charged Mdst_ch5 = ch5.mdstCharged();
			//Mdst_charged Mdst_ch6 = ch6.mdstCharged();
			Mdst_gamma Mdst_ch7 = ch7.mdstGamma();
			Mdst_gamma Mdst_ch8 = ch8.mdstGamma();
			//ExKFitterParticle KFch3(Mdst_ch3, 2); //KFitter 回推碰撞點
			//ExKFitterParticle KFch4(Mdst_ch4, 2);
			//ExKFitterParticle KFch5(Mdst_ch5, 3); //KFitter 回推碰撞點
			//ExKFitterParticle KFch6(Mdst_ch6, 3);
			ExKFitterParticle KF_ch7(Mdst_ch7, IP, IPerr);
			ExKFitterParticle KF_ch8(Mdst_ch8, IP, IPerr);

			ExKFitterVertex eta_Vertex(IP, IPerr);
			ExKFitterMass eta_Mass(META_PDG);

			/*
				ExKFitterParticle B;
				B.LinkParticle(&KFch3);
				B.LinkParticle(&KFch4);
				B.LinkParticle(&KFch5);
				B.LinkParticle(&KFch6);
				B.LinkVertex(&Bpm_Vertex);
			*/
			ExKFitterConstrain con;
			con.SetMassConstrain();
			con.LinkParticle(&KF_ch8);
			con.LinkParticle(&KF_ch7);
			con.LinkMass(&eta_Mass);
			con.LinkVertex(&eta_Vertex);

			ExKFitter Core;
			Core.LinkConstrain(&con);

			int ret = Core.Minimize();

			//double Xi();
			ExKFitterParticle exk_eta;
			exk_eta.LinkParticle(&KF_ch7);
			exk_eta.LinkParticle(&KF_ch8);
			exk_eta.LinkVertex(&eta_Vertex);

			float bvtxx = eta_Vertex.Vertex().x();
			float bvtxy = eta_Vertex.Vertex().y();
			float bvtxz = eta_Vertex.Vertex().z();
			//float chisqExK = Core.Xi();
			HepPoint3D bvertex = eta_Vertex.Vertex();
			if (ret == 0)
			{
				// only B.Update is really needed, the others are done automatically after minimization.
				// the composited vitrual particle track position is at its correspondent vertex, indeed.
				//	B.Update();
				//	KFch3.Update();
				//	KFch4.Update();
				//	KFch5.Update();
				//	KFch6.Update();
				exk_eta.Update();
			}
			b_tpl->column("exk_eta", exk_eta.Mass());
			//cout<<"end save exk  info"<<endl;
			int vertexflag;
			HepPoint3D overtex;
			//for shape variables
			shape((*i), R2, spher, cos_thr, cos_thp, par_sfw, OtherBVector, bvertex, vertexflag, overtex);
			//test****
			b_tpl->column("bvertexx", bvertex.x());
			b_tpl->column("bvertexy", bvertex.y());
			b_tpl->column("bvertexz", bvertex.z());
			b_tpl->column("overtexx", overtex.x());
			b_tpl->column("overtexy", overtex.y());
			b_tpl->column("overtexz", overtex.z());
			b_tpl->column("exfit", ret);
			//b_tpl->column("chisqExK", chisqExK);

			Vector4 ExpectedOtherBVector(0., 0., 0., ebeam);
			Vector4 P_miss = ExpectedOtherBVector - OtherBVector;
			b_tpl->column("pmiss", P_miss.rho());
			b_tpl->column("emiss", P_miss.e());

			// cosine of B and beam dir. the angle between B beam and P_beam but in rest frame or lab frame?
			float cosb = Vector3(b_p.vect()).dot(P_BEAM.vect());
			cosb = cosb / (Vector3(b_p.vect()).mag() * P_BEAM.vect().mag());

			//#ifdef BTUPLE
			//fill ntuple of shape variables
			b_tpl->column("vtflag", vertexflag);
			if (!vertexflag)
			{
				double deltaZ = bvertex.z() - overtex.z();
				b_tpl->column("deltaz", deltaZ);
			}
			b_tpl->column("R2", R2);
			b_tpl->column("costhr", float(cos_thr));
			b_tpl->column("costhp", float(cos_thp));
			b_tpl->column("spher", float(spher));
			b_tpl->column("cosb", cosb);

			b_tpl->column("R2s", par_sfw[0]);
			b_tpl->column("R1so", par_sfw[1]);
			b_tpl->column("R2so", par_sfw[2]);
			b_tpl->column("R3so", par_sfw[3]);
			b_tpl->column("R4so", par_sfw[4]);
			b_tpl->column("R1gso", par_sfw[9]);
			b_tpl->column("R2gso", par_sfw[10]);
			b_tpl->column("R3gso", par_sfw[11]);
			b_tpl->column("R4gso", par_sfw[12]);
			b_tpl->column("R1oo", par_sfw[5]);
			b_tpl->column("R2oo", par_sfw[6]);
			b_tpl->column("R3oo", par_sfw[7]);
			b_tpl->column("R4oo", par_sfw[8]);

			//b_tpl->column("evtflag", evt_flag);

			b_tpl->column("EVTID", Evt);
			b_tpl->column("RUNID", Run);
			b_tpl->column("EXPID", Exp);
			b_tpl->column("eventid", Eventid);
			b_tpl->column("farmid", Farmid);
			//    k_sfw variables  //可以不用動，把訊號跟背景分到兩邊，把兩邊差別變數差異極大化
			k_sfw ksfw_obj((*i));
			const double ksfw(ksfw_obj.fd());
			const int iksfw(ksfw_obj.i_mm2());
			b_tpl->column("imm2", iksfw);
			double miss(ksfw_obj.mm2());
			b_tpl->column("mmm2", miss);
			double et(ksfw_obj.e_t());
			b_tpl->column("et", et);
			double H0oo(ksfw_obj.Hoo(0));
			double H1oo(ksfw_obj.Hoo(1));
			double H2oo(ksfw_obj.Hoo(2));
			double H3oo(ksfw_obj.Hoo(3));
			double H4oo(ksfw_obj.Hoo(4));
			b_tpl->column("H0oo", H0oo);
			b_tpl->column("H1oo", H1oo);
			b_tpl->column("H2oo", H2oo);
			b_tpl->column("H3oo", H3oo);
			b_tpl->column("H4oo", H4oo);
			double H0son(ksfw_obj.Hso_n(0));
			double H1son(ksfw_obj.Hso_n(1));
			double H2son(ksfw_obj.Hso_n(2));
			double H3son(ksfw_obj.Hso_n(3));
			double H4son(ksfw_obj.Hso_n(4));
			b_tpl->column("H0son", H0son);
			b_tpl->column("H1son", H1son);
			b_tpl->column("H2son", H2son);
			b_tpl->column("H3son", H3son);
			b_tpl->column("H4son", H4son);
			double H0soc(ksfw_obj.Hso_c(0));
			double H1soc(ksfw_obj.Hso_c(1));
			double H2soc(ksfw_obj.Hso_c(2));
			double H3soc(ksfw_obj.Hso_c(3));
			double H4soc(ksfw_obj.Hso_c(4));
			b_tpl->column("H0soc", H0soc);
			b_tpl->column("H1soc", H1soc);
			b_tpl->column("H2soc", H2soc);
			b_tpl->column("H3soc", H3soc);
			b_tpl->column("H4soc", H4soc);
			double H0som(ksfw_obj.Hso_m(0));
			double H1som(ksfw_obj.Hso_m(1));
			double H2som(ksfw_obj.Hso_m(2));
			double H3som(ksfw_obj.Hso_m(3));
			double H4som(ksfw_obj.Hso_m(4));
			b_tpl->column("H0som", H0som);
			b_tpl->column("H1som", H1som);
			b_tpl->column("H2som", H2som);
			b_tpl->column("H3som", H3som);
			b_tpl->column("H4som", H4som);

			//cout<<"save shape info"<<endl;

			// B tagging
			// Ntuple entries:
			// r
			//Hamlet 再組另外一個B的方法
			//Hamlet ham;
			//ham.setBcp(B_cand2);
			//b_tpl->column("r", ham.fbtg_mult_dim_likelihood().fq());

			Hamlet hamlet;
			hamlet.setBcp(*i);
			hamlet.setTagMethod(Hamlet::MULT_DIM_LH);
			const double qr_evtgen = hamlet.q();
			//const Fbtag_MultiDimLikelihood0 & mdlh_evtgen = hamlet,fbtag_mult_dim_likelihood();
			b_tpl->column("qr", qr_evtgen);

			int flavor = hamlet.flavor(); // +1 : for tag-side being B0 or B+
			// -1 : for tag-side being B0-bar or B-
			//  0 : does not specify a b-flavor
			// -2 : does not specify a b-flavor (no leptons and kaons are found)
			b_tpl->column("flavor", flavor);
			int imode = hamlet.imode(); // 1 : tag a b-flavor with high-electron method
			// 2 : tag a b-flavor with high-muon method
			// 4 : tag a b-flavor with kaon method
			// others : does not specify a b-flavor
			b_tpl->column("imode", imode);
			//****
			//cout << "start to get Mdstinfo " << endl;

			//------------------------- MC truth matching --------------------------------------------------------------------------------------------------

			// the record of HEP_ID (only for Sig_MC)_from ywdoung
			if (MCstatus == 1)
			{

				const Mdst_gamma ch1 = (*i).child(1).child(0).mdstGamma();
				const Mdst_gamma ch2 = (*i).child(1).child(1).mdstGamma();
				const Mdst_charged ch3 = (*i).child(0).child(0).mdstCharged(); //pi+
				const Mdst_charged ch4 = (*i).child(0).child(1).mdstCharged(); //pi+
				const Mdst_charged ch5 = (*i).child(0).child(2).mdstCharged(); //chic02pi-
				const Mdst_charged ch6 = (*i).child(0).child(3).mdstCharged(); //chic02pi-

				//cout << "1" << endl;

				Gen_hepevt Evtch1 = get_hepevt(ch1); //連到table，去 get_hepevt他是一個class
				Gen_hepevt Evtch2 = get_hepevt(ch2);
				Gen_hepevt Evtch3 = get_hepevt(ch3);
				Gen_hepevt Evtch4 = get_hepevt(ch4);
				Gen_hepevt Evtch5 = get_hepevt(ch5);
				Gen_hepevt Evtch6 = get_hepevt(ch6);

				//cout << "2" << endl;

				b_tpl->column("hi1", Evtch1.idhep()); //call partical的編號，並把它存起來
				b_tpl->column("hi2", Evtch2.idhep());
				b_tpl->column("hi3", Evtch3.idhep());
				b_tpl->column("hi4", Evtch4.idhep());
				b_tpl->column("hi5", Evtch5.idhep());
				b_tpl->column("hi6", Evtch6.idhep());

				//cout << "3" << endl;

				Gen_hepevt EvtP1;
				Gen_hepevt EvtP2;
				Gen_hepevt EvtP3;
				Gen_hepevt EvtP4;
				Gen_hepevt EvtP5;
				Gen_hepevt EvtP6;
				Gen_hepevt EvtGP1;
				Gen_hepevt EvtGP2;
				Gen_hepevt EvtGP3;
				Gen_hepevt EvtGP4;
				Gen_hepevt EvtGP5;
				Gen_hepevt EvtGP6;
                Gen_hepevt EvtGGP1;
                Gen_hepevt EvtGGP2;	
				//cout << "4" << endl;
				if (Evtch1.mo(0))					   //mo=mother particle/如果有媽媽才會有東西，才會執行
				{									   //1279
					EvtP1 = gen_mgr[Evtch1.mo(0) - 1]; //gen_mgr 就是去抓這個值所對應的particle的資訊存進去
					b_tpl->column("mo1", EvtP1.idhep());
					//cout << "mo1 ok" << endl;
					if (EvtP1.mo(0))
					{
						EvtGP1 = gen_mgr[EvtP1.mo(0) - 1];
						b_tpl->column("gmo1", EvtGP1.idhep());
						//cout << "gmo1 ok" << endl;
                        if(EvtGP1.mo(0))
                        {
                            EvtGGP1 = gen_mgr[EvtGP1.mo(0) - 1];
                            b_tpl->column("ggmo1", EvtGGP1.idhep());
                        }
					}
				}
				//cout << "end all mo1*" << endl;

				if (Evtch2.mo(0))
				{
					EvtP2 = gen_mgr[Evtch2.mo(0) - 1];
					b_tpl->column("mo2", EvtP2.idhep());
					//cout << "mo2 ok " << endl;
					if (EvtP2.mo(0))
					{
						EvtGP2 = gen_mgr[EvtP2.mo(0) - 1];
						b_tpl->column("gmo2", EvtGP2.idhep());
						//cout << "gmo2 ok" << endl;
						if(EvtGP2.mo(0))
                                                {
						EvtGGP2 = gen_mgr[EvtGP2.mo(0) - 1];
                                                b_tpl->column("ggmo2", EvtGGP2.idhep());
						}

					}
				}

				//cout << "end all mo2*" << endl;

				if (Evtch3.mo(0))
				{
					EvtP3 = gen_mgr[Evtch3.mo(0) - 1];
					b_tpl->column("mo3", EvtP3.idhep());
					//cout << "mo3 ok" << endl;
					if (EvtP3.mo(0))
					{
						EvtGP3 = gen_mgr[EvtP3.mo(0) - 1];
						b_tpl->column("gmo3", EvtGP3.idhep());
						//cout << "gmo3 ok" << endl;
					}
				}

				//cout << "end all mo3*" << endl;

				if (Evtch4.mo(0))
				{
					EvtP4 = gen_mgr[Evtch4.mo(0) - 1];
					b_tpl->column("mo4", EvtP4.idhep());
					//cout << "mo4 ok" << endl;
					if (EvtP4.mo(0))
					{
						EvtGP4 = gen_mgr[EvtP4.mo(0) - 1];
						b_tpl->column("gmo4", EvtGP4.idhep());
						//cout << "gmo4 ok" << endl;
					}
				}

				//cout << "end all mo4*" << endl;

				if (Evtch5.mo(0))
				{
					EvtP5 = gen_mgr[Evtch5.mo(0) - 1];
					b_tpl->column("mo5", EvtP5.idhep());
					//cout << "mo5 ok" << endl;
					if (EvtP5.mo(0))
					{
						EvtGP5 = gen_mgr[EvtP5.mo(0) - 1];
						b_tpl->column("gmo5", EvtGP5.idhep());
						//cout << "gmo5 ok" << endl;
					}
				}
				//cout << "end all mo5*" << endl;
				if (Evtch6.mo(0))
				{
					EvtP6 = gen_mgr[Evtch6.mo(0) - 1];
					b_tpl->column("mo6", EvtP6.idhep());
					//cout << "mo6 ok" << endl;
					if (EvtP6.mo(0))
					{
						EvtGP6 = gen_mgr[EvtP6.mo(0) - 1];
						b_tpl->column("gmo6", EvtGP6.idhep());
						//cout << "gmo6 ok" << endl;
					}
				}
				//cout << "end all mo6*" << endl;
			//cout << "end mo particle identification" << endl;
				//------------------------------- trace from top to down (from B meson)_from poyuan--------------------	
					   int nbpd = 0, nbnd = 0;
					   for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin(); i != gen_mgr.end(); i++)
					   {
					   double energy = -1;//default value

					   if ((*i).idhep() == 521 || (*i).idhep() == 511)
					   {
					   int da1 = (*i).da(0), da2 = (*i).da(1);
					   if ((*i).idhep() == 521) b_tpl->column("charged", 1);
					   else b_tpl->column("charged", 0);
					   nbpd = (da2 - da1 + 1);
					//da2,da1 is defined in mdst file. nbpd --> the totle child
					//                    cout<<da2<<"-"<<da1<<(da2-da1+1)<<endl;
					b_tpl->column("nbpd", nbpd);
					for (int start = 0; start < (da2 - da1 + 1); start++)
					{
					Gen_hepevt Evda = gen_mgr[(*i).da(0) - 1 + start];
					char bpdnumber[32];
					char bpdnumber_energy[72];
					energy = Evda.P(3);
					sprintf(bpdnumber, "%s%d", "bpd", start + 1);
					sprintf(bpdnumber_energy, "%s%d", "enerbpd", start + 1);
					b_tpl->column(bpdnumber, Evda.idhep());
					b_tpl->column(bpdnumber_energy, energy);
					}
					}
					
					else if ((*i).idhep() == -521 || (*i).idhep() == -511)
					{
					int da1 = (*i).da(0), da2 = (*i).da(1);

					//cout << da2 << "-" << da1 << (da2 - da1 + 1) << endl;
					if ((*i).idhep() == -521) b_tpl->column("charged", 1);
					else b_tpl->column("charged", 0);

					nbnd = (da2 - da1 + 1);
					b_tpl->column("nbnd", nbnd);

					for (int start = 0; start<(da2 - da1 + 1); start++)
					{
					Gen_hepevt Evda = gen_mgr[(*i).da(0) - 1 + start];
					char bndnumber[32];
					char bndnumber_energy[72];
					energy = Evda.P(3);
					sprintf(bndnumber, "%s%d", "bnd", start + 1);//bnd1,bnd2,bnd3....
					sprintf(bndnumber_energy, "%s%d", "enerbnd", start + 1);
					b_tpl->column(bndnumber, Evda.idhep());
					b_tpl->column(bndnumber_energy, energy);
					}
					}
					}//trace from top to down (from B meson) loop


				//					cout << "start to get idhep number" << endl;
				//======================
				//	evt.pdl info
				//======================
				//	eta=221
				//
				//	gamma=22
				//
				//	pi+=211
				//
				//	pi-=-211
				//
				//	K+=321
				//
				//	K-=-321
				//
				//	chic0=10441
				//
				//	B0=511
				//=====================
				int hindex = 0;
				double hrepeat = 0;
				//查表，每一個particle對應一個數字，並確保給一個particle 和他的mother particle 還有granma particle 都依樣

				if (Evtch1.idhep() == 22 && Evtch2.idhep() == 22 && Evtch3.idhep() == 211 && Evtch4.idhep() == 211 && Evtch5.idhep() == -211 && Evtch6.idhep() == -211) //k+ pi-  用代號
				{

					if (EvtP1.idhep() == 221 && EvtP2.idhep() == 221 && EvtP3.idhep() == 10441 && EvtP4.idhep() == 10441 && EvtP5.idhep() == 10441 && EvtP6.idhep() == 10441 &&
						abs(EvtGP1.idhep()) == 511 && abs(EvtGP2.idhep()) == 511 && abs(EvtGP3.idhep()) == 511 && abs(EvtGP4.idhep()) == 511 && abs(EvtGP5.idhep()) == 511 && abs(EvtGP6.idhep()) == 511)

					{
						if (Evtch1.mo(0) == Evtch2.mo(0))
						{
							hindex = 1; // hrepeat = hindex*0.1 + 1 ;
										//cout << "Evtch1.mo==Evtch2.mo" << endl;
						}
						if (Evtch3.mo(0) == Evtch4.mo(0) && Evtch3.mo(0) == Evtch5.mo(0) && Evtch3.mo(0) == Evtch6.mo(0))
						{
							hindex = 1;
							//cout << "Evtch3.mo == Evtch4.mo == Evtch5.mo == Evtch6.mo " << endl;
						}
						else
						{
							hindex = 2;
						}
						if (EvtP3.mo(0) == EvtP4.mo(0) && EvtP3.mo(0) == EvtP5.mo(0) && EvtP3.mo(0) == EvtP6.mo(0) && EvtP3.mo(0) == EvtP1.mo(0) && EvtP3.mo(0) == EvtP2.mo(0))
						{
							hindex = 1;
							//cout << "EvtP3.mo == EvtP4.mo== EvtP5.mo== EvtP6.mo == EvtP1.mo== EvtP2.mo" << endl;
							// hrepeat = hindex*0.1 + 1 ;
						}
						/*
							   		else if (EvtGP3.mo(0) == EvtGP1.mo(0))
							   		{
							  		 hindex = 1;
							   		}	*/
						else
						{
							hindex = 2;
							//hrepeat = hindex*0.1 + 1 ;
							//cout << "hindex=2" << endl;
						}
					}
					else
					{
						hindex = 2;
					}
				}
			else if(abs(Evtch1.idhep()) == 11 && abs(Evtch2.idhep()) == 11 && EvtP1.idhep() == 22 && EvtP2.idhep() == 22 && Evtch3.idhep() == 211 && Evtch4.idhep() == 211 && Evtch5.idhep() == -211 && Evtch6.idhep() == -211)//case1
				{
					if(EvtGP1.idhep() == 221 && EvtGP2.idhep() == 221&& abs(EvtGGP1.idhep()) == 511 && abs(EvtGGP2.idhep()) == 511&& EvtP3.idhep() == 10441 && EvtP4.idhep() == 10441 && EvtP5.idhep() == 10441 && EvtP6.idhep() == 10441&& abs(EvtGP3.idhep()) == 511 && abs(EvtGP4.idhep()) == 511 && abs(EvtGP5.idhep()) == 511 && abs(EvtGP6.idhep()) == 511)
					{
						if( EvtP1.mo(0) == EvtP2.mo(0))
						{
							hindex = 1;
						}
						if(Evtch3.mo(0) == Evtch4.mo(0) && Evtch3.mo(0) == Evtch5.mo(0) && Evtch3.mo(0) == Evtch6.mo(0))
						{
							hindex =1;
						}
						if(EvtP3.mo(0) == EvtP4.mo(0) && EvtP3.mo(0) == EvtP5.mo(0) && EvtP3.mo(0) == EvtP6.mo(0)&& EvtP3.mo(0) == EvtGP1.mo(0)&& EvtP3.mo(0) == EvtGP2.mo(0))
						{
							hindex = 1;
						}
						else
						{
							hindex=2 ;
						}	
					}
					else
					{
						hindex=2 ;
					}
				}
				else if(abs(Evtch1.idhep()) == 11 && Evtch2.idhep() == 22 && EvtP1.idhep() == 22 && EvtP2.idhep() == 221 && Evtch3.idhep() == 211 && Evtch4.idhep() == 211 && 
        Evtch5.idhep() == -211 && Evtch6.idhep() == -211)//case2
				{
					if(EvtGP1.idhep() == 221 && abs(EvtGP2.idhep()) == 511 && abs(EvtGGP1.idhep()) == 511 && abs(EvtGP2.idhep()) == 511&& EvtP3.idhep() == 10441 && EvtP4.idhep() == 10441 && EvtP5.idhep() == 10441 && EvtP6.idhep() == 10441&& abs(EvtGP3.idhep()) == 511 && abs(EvtGP4.idhep()) == 511 && abs(EvtGP5.idhep()) == 511 && abs(EvtGP6.idhep()) == 511)
					{
						if( EvtP1.mo(0) == Evtch2.mo(0))
						{
							hindex = 1;
						}
						if(Evtch3.mo(0) == Evtch4.mo(0) && Evtch3.mo(0) == Evtch5.mo(0) && Evtch3.mo(0) == Evtch6.mo(0))
						{
							hindex =1;
						}
						if(EvtP3.mo(0) == EvtP4.mo(0) && EvtP3.mo(0) == EvtP5.mo(0) && EvtP3.mo(0) == EvtP6.mo(0)&& EvtP3.mo(0) == EvtGP1.mo(0)&& EvtP3.mo(0) == EvtP2.mo(0))
						{
							hindex = 1;
						}
						else
						{
							hindex=2 ;
						}	
					}
					else
					{
						hindex=2 ;
					}
				}
                else if(abs(Evtch2.idhep()) == 11 && Evtch1.idhep() == 22 && EvtP2.idhep() == 22 && EvtP1.idhep() == 221 && Evtch3.idhep() == 211 && Evtch4.idhep() == 211 && 
        Evtch5.idhep() == -211 && Evtch6.idhep() == -211)//case3
				{
					if(EvtGP2.idhep() == 221 && abs(EvtGP1.idhep()) == 511 && abs(EvtGGP2.idhep()) == 511 && abs(EvtGP1.idhep()) == 511&& EvtP3.idhep() == 10441 && EvtP4.idhep() == 10441 && EvtP5.idhep() == 10441 && EvtP6.idhep() == 10441&& abs(EvtGP3.idhep()) == 511 && abs(EvtGP4.idhep()) == 511 && abs(EvtGP5.idhep()) == 511 && abs(EvtGP6.idhep()) == 511)
					{
						if( EvtP2.mo(0) == Evtch1.mo(0))
						{
							hindex = 1;
						}
						if(Evtch3.mo(0) == Evtch4.mo(0) && Evtch3.mo(0) == Evtch5.mo(0) && Evtch3.mo(0) == Evtch6.mo(0))
						{
							hindex =1;
						}
						if(EvtP3.mo(0) == EvtP4.mo(0) && EvtP3.mo(0) == EvtP5.mo(0) && EvtP3.mo(0) == EvtP6.mo(0)&& EvtP3.mo(0) == EvtGP2.mo(0)&& EvtP3.mo(0) == EvtP1.mo(0))
						{
							hindex = 1;
						}
						else
						{
							hindex=2 ;
						}	
					}
					else
					{
						hindex=2 ;
					}
				}
				else
				{
					hindex=2;	
				}
				b_tpl->column("hindex", hindex);
				//b_tpl->column("hrepeat", hrepeat);now use
				//cout << "end to get idhep number" << endl;
				//#endif
			} //end of idhep(MCstatus)
			b_tpl->dumpData();
			*status = 1; //結束後要加，之後才會從新的一行開始
		}
		//cout << "end debug" << endl;

		cout << "(*i)= " << *i << endl;
	}
} //end class
//end of fill pi or k lists from MDST_Charged Data Base
//======================================================
/////////////////以下都是宣告function////////////////////
//======================================================

//void chi_c0_eta_pipipipigg::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P)
void chi_c0_eta_pipipipigg::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P, HepPoint3D &bvertex, int &vertexflag, HepPoint3D &overtex)
{
	TagVK ver; // Vertex Reconstructor with kvertexfitter
	ver.setdefault(b, bvertex);
	//ver.useKsVeto();
	ver.dontUseIPforFit();
	ver.useTubeforFit();
          //cout<<"input bug6 info"<<endl;
	double E_HER = BeamEnergy::E_HER();
	double E_LER = BeamEnergy::E_LER();
	double cross_angle = BeamEnergy::Cross_angle();

	//Particle Type
	Ptype ptype_gamma("GAMM");
	static Vector4 P_BEAM(-E_HER * sin(cross_angle), 0.0, E_LER - E_HER * cos(cross_angle),
						  E_HER + E_LER);

	std::vector<Particle> char_plus;
	std::vector<Particle> char_minus;
	std::vector<Particle> g; //gamma particle

	Ptype ptype_g("GAMM");
	Ptype ptype_char_plus("PI+");
	Ptype ptype_char_minus("PI-");

	std::vector<Particle> candiBfinalParticle;
	std::vector<Vector4> allFinal;
	std::vector<Vector4> candiBfinal;
	std::vector<Vector4> candiBfinalgamm;
	std::vector<Vector4> otherBfinal;
	std::vector<Vector4> Pmiss;

	//fill all charged particle to charged bank
	//into particle list
	Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
	for (std::vector<Mdst_charged>::iterator i = charged_mag.begin();
		 i != charged_mag.end(); i++)
	{
		if ((*i).charge() > 0.)
		{
			Particle tmp(*i, ptype_char_plus);
			char_plus.push_back(tmp);
		}
		else
		{
			Particle tmp(*i, ptype_char_minus);
			char_minus.push_back(tmp);
		}
	}

	//fill gamma list from Mdst_gamma
	//into particle list
	Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
	for (std::vector<Mdst_gamma>::iterator i = gamma_mag.begin();
		 i != gamma_mag.end(); i++)
	{
		Particle tmp(*i);
		g.push_back(tmp);
	}
	//for super fox-wolfram
	//cout << "no of B final particle = " << b.relation().nFinalStateParticles()<<endl;

	for (int i = 0; i < b.relation().nFinalStateParticles();
		 i++)
	{
		candiBfinalParticle.push_back(b.relation().finalStateParticle(i));
		//cout << " i = " << i << " mass = " <<  b.relation().finalStateParticle(i).mass() <<  endl;
		Vector4 temp_P4 = b.relation().finalStateParticle(i).p();
		temp_P4.boost(P_BEAM.boostVector());
		candiBfinal.push_back(temp_P4);
	}

	//by one daughter from b (now omega)
	//cout << "omega final state no ="
	//     << b.relation().child(0).relation().nFinalStateParticles() << endl;
	for (int i = 0; i < b.relation().child(0).relation().nFinalStateParticles(); i++)
	{
		Vector4 temp_P4 =
			b.relation().child(0).relation().finalStateParticle(i).p();
		temp_P4.boost(P_BEAM.boostVector());
		candiBfinalgamm.push_back(temp_P4);
	}

	int numver = 0;
	for (std::vector<Particle>::iterator it1 = char_plus.begin();
		 it1 != char_plus.end(); it1++)
	{
		Vector4 temp_P4 = (*it1).p();
		temp_P4.boost(P_BEAM.boostVector());
		allFinal.push_back(temp_P4);
		int notBDauFlag = 1;
		for (std::vector<Particle>::iterator j = candiBfinalParticle.begin();
			 j != candiBfinalParticle.end(); j++)
			if (checkMultiUse(*it1, *j))
			{
				notBDauFlag = 0;
				break;
			}
		if (notBDauFlag)
		{
			otherBfinal.push_back(temp_P4);
			otherB_P = otherB_P + temp_P4;
			//cout <<" "<<(*it1).relation().mdstCharged().get_ID();
			Particle &tmp = *it1;
			ver.push_back(&tmp);
			numver++;
		}
	}
	for (std::vector<Particle>::iterator it1 = char_minus.begin();
		 it1 != char_minus.end(); it1++)
	{
		Vector4 temp_P4 = (*it1).p();
		temp_P4.boost(P_BEAM.boostVector());
		allFinal.push_back(temp_P4);
		int notBDauFlag = 1;
		for (std::vector<Particle>::iterator j = candiBfinalParticle.begin();
			 j != candiBfinalParticle.end(); j++)
			if (checkMultiUse(*it1, *j))
			{
				notBDauFlag = 0;
				break;
			}
		if (notBDauFlag)
		{
			otherBfinal.push_back(temp_P4);
			otherB_P = otherB_P + temp_P4;
			//cout <<" "<<(*it1).relation().mdstCharged().get_ID();
			Particle &tmp = *it1;
			ver.push_back(&tmp);
			numver++;
		}
	}

	vertexflag = ver.fit();
	if (!vertexflag)
	{
		overtex = ver.vtx();
		//cout<<numver<<":"<<ver.vtx()<<" "<<ver.used_particles().size()<<endl;
	}

	//cout << endl;
	//cout << "non B candi ga_id = ";
	for (std::vector<Particle>::iterator it1 = g.begin();
		 it1 != g.end(); it1++)
	{
		Vector4 temp_P4 = (*it1).p();
		temp_P4.boost(P_BEAM.boostVector());
		allFinal.push_back(temp_P4);
		int notBDauFlag = 1;
		for (std::vector<Particle>::iterator j = candiBfinalParticle.begin();
			 j != candiBfinalParticle.end(); j++)
			if (checkMultiUse(*it1, *j))
			{
				notBDauFlag = 0;
				break;
			}
		if (notBDauFlag)
		{
			otherBfinal.push_back(temp_P4);
			otherB_P = otherB_P + temp_P4;
			//cout <<" "<<(*it1).relation().mdstGamma().get_ID();
		}
	}

	//float ebeam = Benergy();
	double ebeam = BeamEnergy::E_beam_corr();
	HepLorentzVector boost_vector(E_HER * sin(cross_angle), 0.0, -E_LER + E_HER * cos(cross_angle), E_HER + E_LER);
	Vector4 ExpectedOtherBVector(0., 0., 0., ebeam);
	Vector4 OtherBVector;
	OtherBVector.boost(boost_vector.boostVector());
	Vector4 temp_Pmiss = ExpectedOtherBVector - OtherBVector;
	Pmiss.push_back(temp_Pmiss);

	//cout << endl;
	//cout << " no of ch+ = "<< char_plus.size();
	//cout << " no of ch- = "<< char_minus.size();
	//cout << " no of gam = "<< g.size() << endl;

	if ((candiBfinal.size() + otherBfinal.size()) != allFinal.size())
		//cout << candiBfinal.size() << "+" << otherBfinal.size() << " !=" << allFinal.size() << endl;

		//Start to fill shape variables
		Vector3 candiBthr = thrust(candiBfinal);
	Vector3 otherBthr = thrust(otherBfinal);
	//Fox Wolfram Moment
	SuperFoxWolfram foxWolfram;
	SuperFoxWolfram sfwSO;
	SuperFoxWolfram sfwS1O;
	SuperFoxWolfram sfwOO;
	SuperFoxWolfram sfwS3O;
	//SuperFoxWolfram sfwSS;
	foxWolfram.fill(allFinal, allFinal);
	sfwSO.fill(candiBfinal, otherBfinal);
	//sfwSS.fill(candiBfinal,candiBfinal);

	par_sfw[0] = (float)foxWolfram.R(2);
	par_sfw[1] = (float)sfwSO.R(1);
	par_sfw[2] = (float)sfwSO.R(2);
	par_sfw[3] = (float)sfwSO.R(3);
	par_sfw[4] = (float)sfwSO.R(4);

	sfwOO.fill(otherBfinal, otherBfinal);
	par_sfw[5] = (float)sfwOO.R(1);
	par_sfw[6] = (float)sfwOO.R(2);
	par_sfw[7] = (float)sfwOO.R(3);
	par_sfw[8] = (float)sfwOO.R(4);

	sfwS1O.fill(candiBfinalgamm, otherBfinal);
	par_sfw[9] = (float)sfwS1O.R(1);
	par_sfw[10] = (float)sfwS1O.R(2);
	par_sfw[11] = (float)sfwS1O.R(3);
	par_sfw[12] = (float)sfwS1O.R(4);

	//for shape variables
	float H0, R1, R3, R4;
	int ntrk = char_plus.size() + char_minus.size() + g.size(), trkc = 0;
	//int ntrk=g_mag.size()+charged_mag.size(),trkcount=0;
	float ptrk[ntrk * 3];
	int itrk[ntrk];
	//float ebeam = Benergy();

	for (std::vector<Particle>::iterator i = g.begin();
		 i != g.end(); i++)
	{
		Vector4 g_P = (*i).p();
		g_P.boost(P_BEAM.boostVector());
		//Vector3 g_boost = g_P;

		ptrk[trkc * 3] = g_P.px();
		ptrk[trkc * 3 + 1] = g_P.py();
		ptrk[trkc * 3 + 2] = g_P.pz();
		trkc++;
	}

	for (std::vector<Particle>::iterator i = char_plus.begin();
		 i != char_plus.end(); i++)
	{
		Vector4 charged_P = (*i).p();
		charged_P.boost(P_BEAM.boostVector());
		//Vector3 charged_boost = charged_P;

		ptrk[trkc * 3] = charged_P.px();
		ptrk[trkc * 3 + 1] = charged_P.py();
		ptrk[trkc * 3 + 2] = charged_P.pz();
		trkc++;
	} //end for Particle_List

	for (std::vector<Particle>::iterator i = char_minus.begin();
		 i != char_minus.end(); i++)
	{
		Vector4 charged_P = (*i).p();
		charged_P.boost(P_BEAM.boostVector());
		//Vector3 charged_boost = charged_P;

		ptrk[trkc * 3] = charged_P.px();
		ptrk[trkc * 3 + 1] = charged_P.py();
		ptrk[trkc * 3 + 2] = charged_P.pz();
		trkc++;
	} //end for Particle_List
	//end of shape variables

	fwjet2(ntrk, ptrk, ebeam, &H0, &R1, &R2, &R3, &R4);

	// variables for thrust

	std::vector<Vector4> ptl;
	std::vector<Vector4> ptlb;
	Vector3 Bthr;

	for (std::vector<Particle>::iterator j = g.begin();
		 j != g.end(); j++)
	{
		int mask = 0;
		for (int i = 0; i < b.relation().nFinalStateParticles(); i++)
		{
			if (b.relation().finalStateParticle(i).pType() ==
				ptype_gamma)
			{
				if ((*j).relation().mdstGamma().get_ID() ==
					b.relation().finalStateParticle(i).mdstGamma().get_ID())
					mask = 1;
			}
		}
		Vector4 g_P = (*j).p();
		g_P.boost(P_BEAM.boostVector());
		Vector4 tmpptl(g_P.px(), g_P.py(), g_P.pz(),
					   g_P.e());

		if (mask == 1)
		{
			ptlb.push_back(tmpptl);
			//cout << "gamma " << j <<" is been masked" << endl;
		}
		else
		{
			ptl.push_back(tmpptl);
		}
	}

	for (std::vector<Particle>::iterator j = char_plus.begin();
		 j != char_plus.end(); j++)
	{
		int mask = 0;
		for (int i = 0; i < b.relation().nFinalStateParticles(); i++)
		{
			if (b.relation().finalStateParticle(i).pType() !=
				ptype_gamma)
			{
				if ((*j).relation().mdstCharged().get_ID() ==
					b.relation().finalStateParticle(i).mdstCharged().get_ID())
					mask = 1;
			}
		}

		Vector4 charged_P = (*j).p();
		charged_P.boost(P_BEAM.boostVector());
		Vector4 tmpptl(charged_P.px(), charged_P.py(), charged_P.pz(),
					   charged_P.e());

		if (mask == 1)
		{
			ptlb.push_back(tmpptl);
			//cout << "track " << j <<" is been masked" << endl;
		}
		else
		{
			ptl.push_back(tmpptl);
		}
	}
	for (std::vector<Particle>::iterator j = char_minus.begin();
		 j != char_minus.end(); j++)
	{
		int mask = 0;
		for (int i = 0; i < b.relation().nFinalStateParticles(); i++)
		{
			if (b.relation().finalStateParticle(i).pType() !=
				ptype_gamma)
			{
				if ((*j).relation().mdstCharged().get_ID() ==
					b.relation().finalStateParticle(i).mdstCharged().get_ID())
					mask = 1;
			}
		}

		Vector4 charged_P = (*j).p();
		charged_P.boost(P_BEAM.boostVector());
		Vector4 tmpptl(charged_P.px(), charged_P.py(), charged_P.pz(),
					   charged_P.e());

		if (mask == 1)
		{
			ptlb.push_back(tmpptl);
			//cout << "track " << j <<" is been masked" << endl;
		}
		else
		{

			ptl.push_back(tmpptl);
		}
	}
	//end for Particle_List

	//one of the B daughters, used to calculate the thrust angle

	Vector4 ome_P = b.relation().child(0).p();
	ome_P.boost(P_BEAM.boostVector());

	Vector4 b_P = b.p();
	b_P.boost(P_BEAM.boostVector());

	Vector3 thr_axis = thrust(ptl);
	Vector3 thr_axis_b = thrust(ptlb);

	cos_thr = (thr_axis.dot(thr_axis_b));
	cos_thr = cos_thr / (thr_axis.mag() * thr_axis_b.mag());

	//Vector3 tmpjet(b_P.px(),b_P.py(),b_P.pz());
	Vector3 tmpjet(ome_P.px(), ome_P.py(), ome_P.pz());
	spher = Sper(ptl, tmpjet);
}
//for int t--> 0:e 1:mu 2:pi 3:k 4:p
void chi_c0_eta_pipipipigg::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t)
{

	*dr = 1000.0;
	*dz = 1000.0;

	if (charged->trk())
	{
		Mdst_trk_fit &trkFit = charged->trk().mhyp(t);
		if (trkFit)
		{
			HepVector a(5, 0);
			a[0] = trkFit.helix(0);
			a[1] = trkFit.helix(1);
			a[2] = trkFit.helix(2);
			a[3] = trkFit.helix(3);
			a[4] = trkFit.helix(4);

			HepPoint3D pivot(trkFit.pivot(0), trkFit.pivot(1), trkFit.pivot(2));
			Helix ltrk(pivot, a);
			ltrk.pivot(IP);

			*dr = ltrk.dr();
			*dz = ltrk.dz();
		}
	}
}
//#if defined(BELLE_NAMESPACE)
} // namespace Belle
//#endif
