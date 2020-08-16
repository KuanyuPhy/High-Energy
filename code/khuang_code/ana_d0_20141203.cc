

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
#include "hamlet/Hamlet.h" //tagging
#include "tagv/TagV.h"

#include "ExKFitter.h"

#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"

#include "kid_eff_06.h"
#include "kid/atc_pid.h"  
#include "eid/eid.h"       
#include "mdst/Muid_mdst.h"
#include "mdst/mdst.h"
#include HEPEVT_H
#include EVTCLS_H
#include MDST_H
#include BELLETDF_H
#include "mdst/findKs.h"
#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/khelix2xyz.h"
#include "kfitter/kfitterparticle.h"
#include "helix/Helix.h"
#include "k_sfw.h"

#include "CLHEP/Vector/LorentzVector.h"  // /afs/afs11/belle/belle/b20090127_0910/src/util/belleCLHEP/Vector
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <iostream.h>
#include "userinfo.h"
#include "panther/panther.h"
#include "ip/IpProfile.h"
//#include "./mdst3.icc"

#include MDST_H
#include BELLETDF_H
#include HEPEVT_H 



//the unit is GeV
//#define E_HER 7.998213 // high energy ring i.e electron ring
//#define E_LER 3.499218    // low energy ring i.e. positron ring



#define BTUPLE

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//utilities

int checkMultiUse(Particle& i, Particle& j);


// Module class
class ana_d0_20141203 : public Module 
{
public:
  ana_d0_20141203(void){};
  ~ana_d0_20141203(void){};
  void init(int*);
  void term(void);
  void disp_stat(const char*){};
  void hist_def(void);
  void event(BelleEvent*, int*); 
  void begin_run(BelleEvent*, int*);
  void end_run(BelleEvent*, int*){};
  void other(int*, BelleEvent*, int*){};
//  void shape(Particle &, float &, float &, float &, float &, float par_sfw[17], Vector4 &);
  void shape(Particle &, float &, float &, float &, float &, float par_sfw[17],Vector4 &,HepPoint3D &,int &, HepPoint3D &);

  void GetImpactParameters(const Mdst_charged*, double* ,double* ,int);

  double deltaZ(Particle, double&, double&, double&, double&);


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

};


extern "C" Module_descr *mdcl_ana_d0_20141203()
{
  ana_d0_20141203 *module = new ana_d0_20141203;
  Module_descr *dscr = new Module_descr ( "ana_d0_20141203", module );
  IpProfile::define_global(dscr);
  BeamEnergy::define_global(dscr);
  return dscr;
}



void 
ana_d0_20141203::hist_def(void)
{ 
  extern BelleTupleManager *BASF_Histogram;  
  BelleTupleManager& tm = *BASF_Histogram;
#ifdef BTUPLE

//ntuple name length--> 8 characters(upper limit)
  b_tpl = BASF_Histogram->ntuple("B0 --> D0 pi","massd0 pxd0 pyd0 pzd0 ed0 mass0 px0 py0 pz0 e0 pid0 chrg0 dr0 dz0 cos0 mass1 px1 py1 pz1 e1 pid1 chrg1 dr1 dz1 cos1 vchisq");

#endif

  e_gamm   = BASF_Histogram->histogram("Energy of Photon", 100, 0.0, 3.5);
  m_pi0    = BASF_Histogram->histogram("Mass(Pi0 - GeV)",100, 0.115, 0.145);
  m_b      = BASF_Histogram->histogram("Mass(B - GeV)",50, 5.2, 5.3);
  m_b_beam = BASF_Histogram->histogram("Massfit(B - GeV)",50, 5.2, 5.3);
  dec_lenx = BASF_Histogram->histogram("LAM decay length x cm",40, 0., 20.);
  dec_leny = BASF_Histogram->histogram("LAM decay length y cm",40, 0., 20.);
  dec_alenx= BASF_Histogram->histogram("ALAM decay length x cm",40, 0., 20.);
  dec_aleny= BASF_Histogram->histogram("ALAM decay length y cm",40, 0., 20.);
  m_lamv   = BASF_Histogram->histogram("Mass(LamdaVee2- GeV)",100, 1.05, 1.15);
  m_lamppi = BASF_Histogram->histogram("Mass(Lamda ppi- GeV)",100, 1.05, 1.15);
  m_lamkf  = BASF_Histogram->histogram("Mass(Lamda fit- GeV)",100, 1.05, 1.15);
  Fisher_fw.histogram(10,"Super Fox-Wolfram");
  k_sfw::initialize(Fisher_ksfw);

  p_elec = BASF_Histogram->histogram("P_electron", 50, 0.0, 5.0);
  p_muon = BASF_Histogram->histogram("P_muon", 50, 0.0, 5.0);
} 


void ana_d0_20141203::init(int*)
{
  evtcount=0;

  Hamlet::init();

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


    eff_s1_kpi_kp.init( .6, 1, "track_s1kpi_kp", "kideff-2006-svd1-pos.dat" );//K plus eff.
    eff_s1_kpi_km.init( .6, 1, "track_s1kpi_km", "kideff-2006-svd1-neg.dat" );//K minus eff.
    eff_s1_kpi_kp_fake.init( .6, 2, "track_s1kpi_kp_fake", "kideff-2006-svd1-pos.dat" );//K fake plus eff.
    eff_s1_kpi_km_fake.init( .6, 2, "track_s1kpi_km_fake", "kideff-2006-svd1-neg.dat" );//K fake minus eff.

    eff_s1_kpi_pip.init( .6, 3, "track_s1kpi_pip", "kideff-2006-svd1-pos.dat" );//K plus eff.
    eff_s1_kpi_pim.init( .6, 3, "track_s1kpi_pim", "kideff-2006-svd1-neg.dat" );//K minus eff.
    eff_s1_kpi_pip_fake.init( .6, 4, "track_s1kpi_pip_fake", "kideff-2006-svd1-pos.dat" );//K fake plus eff.
    eff_s1_kpi_pim_fake.init( .6, 4, "track_s1kpi_pim_fake", "kideff-2006-svd1-neg.dat" );//K fake minus eff.

// just test for comparing the previous study
//    eff_s1_kpi_kp.init( .6, 1, "track_s1kpi_kp", "kideff-2005-svd1-pos.dat" );//K plus eff.
//    eff_s1_kpi_km.init( .6, 1, "track_s1kpi_km", "kideff-2005-svd1-neg.dat" );//K minus eff.
//    eff_s1_kpi_kp_fake.init( .6, 2, "track_s1kpi_kp_fake", "kideff-2005-svd1-pos.dat" );//K fake plus eff.
//    eff_s1_kpi_km_fake.init( .6, 2, "track_s1kpi_km_fake", "kideff-2005-svd1-neg.dat" );//K fake minus eff.
  
//    eff_s1_kpi_pip.init( .6, 3, "track_s1kpi_pip", "kideff-2006-svd1-pos.dat" );//K plus eff.
//    eff_s1_kpi_pim.init( .6, 3, "track_s1kpi_pim", "kideff-2006-svd1-neg.dat" );//K minus eff.
//    eff_s1_kpi_pip_fake.init( .6, 4, "track_s1kpi_pip_fake", "kideff-2006-svd1-pos.dat" );//K fake plus eff.
//    eff_s1_kpi_pim_fake.init( .6, 4, "track_s1kpi_pim_fake", "kideff-2006-svd1-neg.dat" );//K fake minus eff.


    eff_s1_pipi_pip.init( .6, 3, "track_s1pipi_pip", "kideff-2006-svd1-pos.dat" );//pi plus eff.
    eff_s1_pipi_pip_fake.init( .6, 4, "track_s1pipi_pip_fake", "kideff-2006-svd1-pos.dat" );//pi fake plus eff.
    eff_s1_pipi_pim.init( .6, 3, "track_s1pipi_pim", "kideff-2006-svd1-neg.dat" );//pi minus eff.
    eff_s1_pipi_pim_fake.init( .6, 4, "track_s1pipi_pim_fake", "kideff-2006-svd1-neg.dat" );//pi fake minus eff.


// just test for comparing the previous study
//    eff_s1_pipi_pip.init( .6, 3, "track_s1pipi_pip", "kideff-2005-svd1-pos.dat" );//pi plus eff.
//    eff_s1_pipi_pip_fake.init( .6, 4, "track_s1pipi_pip_fake", "kideff-2005-svd1-pos.dat" );//pi fake plus eff.
//    eff_s1_pipi_pim.init( .6, 3, "track_s1pipi_pim", "kideff-2005-svd1-neg.dat" );//pi minus eff.
//    eff_s1_pipi_pim_fake.init( .6, 4, "track_s1pipi_pim_fake", "kideff-2005-svd1-neg.dat" );//pi fake minus eff.


    eff_s1_kk_kp.init( .6, 1, "track_s1kk_kp", "kideff-2006-svd1-pos.dat" );//k plus eff.
    eff_s1_kk_km.init( .6, 1, "track_s1kk_km", "kideff-2006-svd1-neg.dat" );//k minus eff.

    eff_s1_kk_kp_fake.init( .6, 2, "track_s1kk_kp_fake", "kideff-2006-svd1-pos.dat" );//k fake plus eff.
    eff_s1_kk_km_fake.init( .6, 2, "track_s1kk_km_fake", "kideff-2006-svd1-neg.dat" );//k fake minus eff.

// just test for comparing the previous study
//    eff_s1_kk_kp.init( .6, 1, "track_s1kk_kp", "kideff-2005-svd1-pos.dat" );//k plus eff.
//    eff_s1_kk_km.init( .6, 1, "track_s1kk_km", "kideff-2005-svd1-neg.dat" );//k minus eff.

//    eff_s1_kk_kp_fake.init( .6, 2, "track_s1kk_kp_fake", "kideff-2005-svd1-pos.dat" );//k fake plus eff.
//    eff_s1_kk_km_fake.init( .6, 2, "track_s1kk_km_fake", "kideff-2005-svd1-neg.dat" );//k fake minus eff.



    eff_s2_kpi_kp.init( .6, 1, "track_s2kpi_kp", "kideff-2010-svd2-pos.dat" );//K plus eff.
    eff_s2_kpi_km.init( .6, 1, "track_s2kpi_km", "kideff-2010-svd2-neg.dat" );//K minus eff.
    eff_s2_kpi_kp_fake.init( .6, 2, "track_s2kpi_kp_fake", "kideff-2010-svd2-pos.dat" );//K fake plus eff.
    eff_s2_kpi_km_fake.init( .6, 2, "track_s2kpi_km_fake", "kideff-2010-svd2-neg.dat" );//K fake minus eff.

    eff_s2_kpi_pip.init( .6, 3, "track_s2kpi_pip", "kideff-2010-svd2-pos.dat" );//K plus eff.
    eff_s2_kpi_pim.init( .6, 3, "track_s2kpi_pim", "kideff-2010-svd2-neg.dat" );//K minus eff.
    eff_s2_kpi_pip_fake.init( .6, 4, "track_s2kpi_pip_fake", "kideff-2010-svd2-pos.dat" );//K fake plus eff.
    eff_s2_kpi_pim_fake.init( .6, 4, "track_s2kpi_pim_fake", "kideff-2010-svd2-neg.dat" );//K fake minus eff.
  
    eff_s2_pipi_pip.init( .6, 3, "track_s2pipi_pip", "kideff-2010-svd2-pos.dat" );//pi plus eff.
    eff_s2_pipi_pip_fake.init( .6, 4, "track_s2pipi_pip_fake", "kideff-2010-svd2-pos.dat" );//pi fake plus eff.
    eff_s2_pipi_pim.init( .6, 3, "track_s2pipi_pim", "kideff-2010-svd2-neg.dat" );//pi minus eff.
    eff_s2_pipi_pim_fake.init( .6, 4, "track_s2pipi_pim_fake", "kideff-2010-svd2-neg.dat" );//pi fake minus eff.
  
    eff_s2_kk_kp.init( .6, 1, "track_s2kk_kp", "kideff-2010-svd2-pos.dat" );//k plus eff.
    eff_s2_kk_km.init( .6, 1, "track_s2kk_km", "kideff-2010-svd2-neg.dat" );//k minus eff.

    eff_s2_kk_kp_fake.init( .6, 2, "track_s2kk_kp_fake", "kideff-2010-svd2-pos.dat" );//k fake plus eff.
    eff_s2_kk_km_fake.init( .6, 2, "track_s2kk_km_fake", "kideff-2010-svd2-neg.dat" );//k fake minus eff.


}


void ana_d0_20141203::term (void)
{
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

}

// Set IP location
HepPoint3D IP(0,0,0);
HepSymMatrix IPerr(3,0);


void ana_d0_20141203::begin_run(BelleEvent *evptr, int *status)
{
  eid::init_data();  // available in the new lib (after b199907*)  

  //beam energy
  BeamEnergy::begin_run();

  // Get IP profile data from $BELLE_POSTGRES_SERVER
  IpProfile::begin_run();
  
  // Dump IP profile data to STDOUT (optional)
  IpProfile::dump();

  // Set IP and error
  IP    = IpProfile::position();
  IPerr = IpProfile::position_err();

  //tagging 
  // To initialize LH tables by EvtGen MC 
  Hamlet::begin_run(Hamlet::MULT_DIM_LH);
}
  

void 
ana_d0_20141203::event(BelleEvent *evptr, int *status)
{ 
       const HepPoint3D             &ip     = IpProfile::position();

       Belle_event_Manager& bevt_mgr        = Belle_event_Manager::get_manager();
       Mdst_charged_Manager &charged_mag    = Mdst_charged_Manager::get_manager();
       Mdst_gamma_Manager &gamma_mag        = Mdst_gamma_Manager::get_manager();
       Mdst_pi0_Manager &pi0_mag            = Mdst_pi0_Manager::get_manager();
       Mdst_vee2_Manager& Vee2Mgr           = Mdst_vee2_Manager::get_manager();
       Mdst_event_add_Manager& mevtmgr      = Mdst_event_add_Manager::get_manager();
       Mdst_vee_daughters_Manager& veedmgr  = Mdst_vee_daughters_Manager::get_manager();
       Mdst_klm_mu_ex_Manager& klmmgr       = Mdst_klm_mu_ex_Manager::get_manager();
       Mdst_ecl_aux_Manager &eclaux_mag     = Mdst_ecl_aux_Manager::get_manager();
       Evtcls_hadronic_flag_Manager &evtcls = Evtcls_hadronic_flag_Manager::get_manager();
       Gen_hepevt_Manager &gen_mgr          = Gen_hepevt_Manager::get_manager();

  //for reprocessed data
  //  remove_duplicates();
  //scale_momenta(1.00246, 1.0);
  //for reprocessed exp7 data
  //scale_momenta(1.00221, 1.0);
  //for reprocessed exp9 data
  //scale_momenta(1.00149, 1.0);


  evtcount++;
  //cout << "event = " << evtcount << endl;

  double E_HER=BeamEnergy::E_HER();
  double E_LER=BeamEnergy::E_LER();
  double cross_angle=BeamEnergy::Cross_angle();//radian
  static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle), E_HER+E_LER);


  int Run=0;
  int Evt=0;
  int Eventid;
  int Farmid;
  int Exp;

  Run = bevt_mgr[0].RunNo();
  Evt = bevt_mgr[0].EvtNo();
  Exp = bevt_mgr[0].ExpNo();
  Eventid = (bevt_mgr[0].EvtNo() & 0x0FFFFFFF );
  Farmid = bevt_mgr[0].EvtNo() >> 28;


// set IP and error
	int IPUsable = 0;
	if(IpProfile::usable())
	{
         IP    = IpProfile::position(1);
         IPerr = IpProfile::position_err_b_life_smeared(1);
         IPUsable = 1;
    	}
	else 
	{
         IP    = HepPoint3D(0,0,0);
         IPerr = HepSymMatrix(3,0);
    	}



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
  std::vector<Particle> Lamda;
  std::vector<Particle> Lamdabar;
  std::vector<Particle> Sigma_plus;
  std::vector<Particle> mu_plus;
  std::vector<Particle> mu_minus;
  std::vector<Particle> e_plus;
  std::vector<Particle> e_minus;
  std::vector<Particle> B_cand1;
  std::vector<Particle> B_cand2;
  std::vector<Particle> B_cand3;
  std::vector<Particle> B_cand4;
  std::vector<Particle> B_cand5;
std::vector<Particle> D0;
 std::vector<Particle> D0B;


  //Particle Type(/belle/belle/b20040727_1143/share/data-files/qq98/decay.dec)
  Ptype ptype_gamma("GAMM");
Ptype ptype_D0("D0");
Ptype ptype_D0B("D0B");
  Ptype ptype_pi_0("PI0");
  Ptype ptype_pi_plus("PI+");
  Ptype ptype_pi_minus("PI-");
  Ptype ptype_k_plus("K+"); 
  Ptype ptype_k_minus("K-");
  Ptype ptype_k_short("K0S");
  Ptype ptype_p_plus("P+");
  Ptype ptype_p_minus("AP+");
  Ptype ptype_Lamda("LAM");
  Ptype ptype_Lamdabar("ALAM");
  Ptype ptype_Sigma_plus("SIG+");
  Ptype ptype_mu_plus("MU+");
  Ptype ptype_mu_minus("MU-");
  Ptype ptype_e_plus("E+");
  Ptype ptype_e_minus("E-");
  Ptype ptype_Bu_cand("B+");
  Ptype ptype_Bd_cand("B0");
  Ptype ptype_B0B_cand("B0B");
  Ptype ptype_Bs_cand("BS0");


 

//  atc_pid selKpi(0,1,0,3,2);//for the reprocessed exp7 data
    atc_pid selkpi(3,1,5,3,2); 
    atc_pid selkp(3,1,5,3,4);
//  atc_pid selppi(3,1,5,4,2); 
//  atc_pid selpk(3,1,5,4,3);
//  atc_pid selmuk(3,1,5,1,3);
//  atc_pid selmupi(3,1,5,1,2); 


  // distinguish MC from data,MC is MCstatus == 1
  int MCstatus=0;
  for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();
       i != gen_mgr.end(); i++)
  {
    MCstatus=1;
  }


  //fill all partilce from Gen_hepevt Data Base



  //fill pi or k lists from MDST_Charged Data Base
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin(); i != charged_mag.end(); i++)
  {

        Muid_mdst muon( *i );

        Mdst_charged& ch = *i;

        int mu_sta = muon.Status();
        int outcome = muon.Outcome();
        int mu_level = muon.Level();
        int reject = muon.Prerejection();
        double mu_like = muon.Muon_likelihood();
        double chi2 = muon.Chi_2();


      //for eid
      eid sel_e(ch);
      //float eid_prob = sel_e.prob(0,-1,0);
      float eid_prob = sel_e.prob(3,1,5);


   

    if((reject || mu_like < 0.95) && eid_prob < 0.95 )
	{
          if(selkpi.prob(*i)>0.6)
	//if(selkpi.prob(*i)>0.)
	  {
              	if ((*i).charge()>0)
		{
                  if (ch.trk().quality() == 0)  // Obtain "Good track"
                  {
                  Particle tmp(*i, ptype_k_plus);
                  k_plus.push_back(tmp);
              	  }
                }
		else
		{
                 if (ch.trk().quality() == 0)  // Obtain "Good track"
                  {
                  Particle tmp(*i, ptype_k_minus);
                  k_minus.push_back(tmp);
                  }
              	}
          }
	  else if (selkpi.prob(*i)<0.4)
	 //if (selkpi.prob(*i)<1.0)
	  {
              if ((*i).charge()>0)
		{
		 if (ch.trk().quality() == 0)  // Obtain "Good track"
                  {
                  Particle tmp(*i, ptype_pi_plus);
                  pi_plus.push_back(tmp);   
		  }
	        }
	      else
		{
                 if (ch.trk().quality() == 0)  // Obtain "Good track"
                  {
                  Particle tmp(*i, ptype_pi_minus);
                  pi_minus.push_back(tmp);
		  }
              	}	
          }
	}//end loop of e, mu rejection


}
//end of filling the charged  particles



    combination(D0,ptype_D0B,k_plus,pi_minus);
    combination(D0,ptype_D0,k_minus,pi_plus);

//for B candicate 2:

  for(std::vector<Particle>::iterator i =D0.begin(); i != D0.end(); i++)
  { 

#ifdef BTUPLE

//------------------- information of D0 (*i) -------------------------

	b_tpl->column("massd0", (*i).p().mag());
	b_tpl->column("pxd0", (*i).px());
	b_tpl->column("pyd0", (*i).py());
	b_tpl->column("pzd0", (*i).pz());
	b_tpl->column("ed0", (*i).e());



//------------------- information of K (*i).child(0) -------------------------

        b_tpl->column("mass0", (*i).child(0).p().mag());
        b_tpl->column("px0", (*i).child(0).px());
        b_tpl->column("py0", (*i).child(0).py());
        b_tpl->column("pz0", (*i).child(0).pz());
        b_tpl->column("e0", (*i).child(0).e());

        const Mdst_charged* chk0 = &(*i).child(0).mdstCharged();
        float k0pi_id = selkpi.prob(chk0);
        b_tpl->column("pid0", k0pi_id);
          
        float chrg0 = chk0->charge();
        b_tpl->column("chrg0", chrg0);

        double dr0,dz0;
        GetImpactParameters(chk0,&dr0,&dz0,3);
        b_tpl->column("dr0", float(dr0));
        b_tpl->column("dz0", float(dz0));

        double cos_polar_0;
        Vector3 ch0_momentum((*i).child(0).px(), (*i).child(0).py(), (*i).child(0).pz());
        Vector3 P_BEAM_polar ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));

        cos_polar_0=ch0_momentum.dot(P_BEAM_polar);
        cos_polar_0=cos_polar_0/(ch0_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos0", cos_polar_0);



//------------------- information of pi (*i).child(1) -------------------------

        b_tpl->column("mass1", (*i).child(1).p().mag());
        b_tpl->column("px1", (*i).child(1).px());
        b_tpl->column("py1", (*i).child(1).py());
        b_tpl->column("pz1", (*i).child(1).pz());
        b_tpl->column("e1", (*i).child(1).e());

        const Mdst_charged* chk1 = &(*i).child(1).mdstCharged();
        float k1pi_id = selkpi.prob(chk1);
        b_tpl->column("pid1", k1pi_id);

        float chrg1 = chk1->charge();
        b_tpl->column("chrg1", chrg1);

        double dr1,dz1;
        GetImpactParameters(chk1,&dr1,&dz1,2);
        b_tpl->column("dr1", float(dr1));
        b_tpl->column("dz1", float(dz1));

        double cos_polar_1;
        Vector3 ch1_momentum((*i).child(0).px(), (*i).child(0).py(), (*i).child(0).pz());

        cos_polar_1=ch1_momentum.dot(P_BEAM_polar);
        cos_polar_1=cos_polar_1/(ch1_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos1", cos_polar_1);

//------------------- D0 vertex fit with K fitter  -------------------------

        kvertexfitter kvf;
        addTrack2fit(kvf, (*i).child(0)); // add "K+(-)" to kfitter
        addTrack2fit(kvf, (*i).child(1)); // add "pi-(+)" to kfitter
        unsigned err = kvf.fit(); // do "fitting"
        float vchisq;
        if(err == 0) // success.
              vchisq = kvf.chisq();
        else vchisq = -1.;

	b_tpl->column("vchisq", vchisq);


	b_tpl->dumpData(); 
        *status = 1;
#endif

}//end of D0




//B_cand1 for k pi



}//end of fill pi or k lists from MDST_Charged Data Base

//void ana_d0_20141203::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P)
void ana_d0_20141203::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex)
{
   TagVK ver; // Vertex Reconstructor with kvertexfitter
   ver.setdefault(b,bvertex);
   //ver.useKsVeto();
   ver.dontUseIPforFit();
   ver.useTubeforFit();

  double E_HER=BeamEnergy::E_HER();
  double E_LER=BeamEnergy::E_LER();
  double cross_angle=BeamEnergy::Cross_angle();   
 
  //Particle Type
  Ptype ptype_gamma("GAMM");  
  static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle),
                          E_HER+E_LER );
  
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
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); i++)
  {
    if ( (*i).charge() > 0.)
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
  
  for( int i=0;i< b.relation().nFinalStateParticles();
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
  for(int i=0; i<b.relation().child(0).relation().nFinalStateParticles();i++) 
  {
    Vector4 temp_P4 =
      b.relation().child(0).relation().finalStateParticle(i).p();
    temp_P4.boost(P_BEAM.boostVector());
    candiBfinalgamm.push_back(temp_P4);
  }
  
  int numver=0;
  for( std::vector<Particle>::iterator it1=char_plus.begin();
       it1!=char_plus.end(); it1++ ) 
  {
    Vector4 temp_P4 = (*it1).p();
    temp_P4.boost(P_BEAM.boostVector());
    allFinal.push_back(temp_P4);
    int notBDauFlag=1;
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag) 
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //cout <<" "<<(*it1).relation().mdstCharged().get_ID();
      Particle& tmp=*it1;
      ver.push_back(&tmp);
      numver++;
      }
  }
  for ( std::vector<Particle>::iterator it1=char_minus.begin();
        it1!=char_minus.end(); it1++ ) 
  {
    Vector4 temp_P4 = (*it1).p();
    temp_P4.boost(P_BEAM.boostVector());
    allFinal.push_back(temp_P4);
    int notBDauFlag=1;
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag) 
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //cout <<" "<<(*it1).relation().mdstCharged().get_ID();
      Particle& tmp=*it1;
      ver.push_back(&tmp);
      numver++;
      }  
  }

  vertexflag=ver.fit();
  if(!vertexflag) 
  {
    overtex = ver.vtx();
    //cout<<numver<<":"<<ver.vtx()<<" "<<ver.used_particles().size()<<endl;
  }

  //cout << endl;
  //cout << "non B candi ga_id = ";
  for ( std::vector<Particle>::iterator it1=g.begin();
        it1!= g.end(); it1++ ) 
  {
    Vector4 temp_P4 = (*it1).p();
    temp_P4.boost(P_BEAM.boostVector());
    allFinal.push_back(temp_P4);
    int notBDauFlag=1;
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag) 
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //cout <<" "<<(*it1).relation().mdstGamma().get_ID();
      }
  }

//  float ebeam = Benergy();
      double ebeam = BeamEnergy::E_beam_corr();
  HepLorentzVector boost_vector(E_HER*sin(cross_angle), 0.0,-E_LER+E_HER*cos(cross_angle),E_HER+E_LER);
  Vector4 ExpectedOtherBVector(0.,0.,0.,ebeam);
  Vector4 OtherBVector;
  OtherBVector.boost(boost_vector.boostVector());
  Vector4 temp_Pmiss = ExpectedOtherBVector - OtherBVector;
  Pmiss.push_back(temp_Pmiss);

  //cout << endl;
  //cout << " no of ch+ = "<< char_plus.size();
  //cout << " no of ch- = "<< char_minus.size();
  //cout << " no of gam = "<< g.size() << endl;

  if((candiBfinal.size()+otherBfinal.size()) != allFinal.size())
    cout << candiBfinal.size()<<"+"<<otherBfinal.size()<<" !="<< allFinal.size()<<endl;


  //Start to fill shape variables
  Vector3 candiBthr=thrust(candiBfinal);
  Vector3 otherBthr=thrust(otherBfinal);
  //Fox Wolfram Moment
  SuperFoxWolfram foxWolfram;
  SuperFoxWolfram sfwSO;
  SuperFoxWolfram sfwS1O;
  SuperFoxWolfram sfwOO;
  SuperFoxWolfram sfwS3O;
  //SuperFoxWolfram sfwSS;
  foxWolfram.fill(allFinal,allFinal);
  sfwSO.fill(candiBfinal,otherBfinal);
  //sfwSS.fill(candiBfinal,candiBfinal);

  par_sfw[0] = (float)foxWolfram.R(2);
  par_sfw[1] = (float)sfwSO.R(1);
  par_sfw[2] = (float)sfwSO.R(2);
  par_sfw[3] = (float)sfwSO.R(3);
  par_sfw[4] = (float)sfwSO.R(4);

  sfwOO.fill(otherBfinal,otherBfinal);
  par_sfw[5] = (float)sfwOO.R(1);
  par_sfw[6] = (float)sfwOO.R(2);
  par_sfw[7] = (float)sfwOO.R(3);
  par_sfw[8] = (float)sfwOO.R(4);

  sfwS1O.fill(candiBfinalgamm,otherBfinal);
  par_sfw[9] = (float)sfwS1O.R(1);
  par_sfw[10] = (float)sfwS1O.R(2);
  par_sfw[11] = (float)sfwS1O.R(3);
  par_sfw[12] = (float)sfwS1O.R(4);

  //for shape variables
  float H0,R1,R3,R4;
  int ntrk=char_plus.size()+char_minus.size()+g.size(),trkc=0;
  //int ntrk=g_mag.size()+charged_mag.size(),trkcount=0;
  float ptrk[ntrk*3];
  int itrk[ntrk];
  //float ebeam = Benergy();

  for (std::vector<Particle>::iterator i = g.begin();
       i != g.end(); i++)
  {
    Vector4 g_P= (*i).p();
    g_P.boost(P_BEAM.boostVector());
    //Vector3 g_boost = g_P;

    ptrk[trkc*3]=g_P.px();
    ptrk[trkc*3+1]=g_P.py();
    ptrk[trkc*3+2]=g_P.pz();
    trkc++;
  }

  for (std::vector<Particle>::iterator i = char_plus.begin();
       i != char_plus.end(); i++) 
  {
    Vector4 charged_P = (*i).p();
    charged_P.boost(P_BEAM.boostVector());
    //Vector3 charged_boost = charged_P;

    ptrk[trkc*3]=charged_P.px();
    ptrk[trkc*3+1]=charged_P.py();
    ptrk[trkc*3+2]=charged_P.pz();
    trkc++;
  } //end for Particle_List

  for (std::vector<Particle>::iterator i = char_minus.begin();
       i != char_minus.end(); i++) 
  {
    Vector4 charged_P = (*i).p();
    charged_P.boost(P_BEAM.boostVector());
    //Vector3 charged_boost = charged_P;

    ptrk[trkc*3]=charged_P.px();
    ptrk[trkc*3+1]=charged_P.py();
    ptrk[trkc*3+2]=charged_P.pz();
    trkc++;
  } //end for Particle_List
  //end of shape variables

  fwjet2(ntrk,ptrk,ebeam,&H0,&R1,&R2,&R3,&R4);


  // variables for thrust

  std::vector<Vector4> ptl;
  std::vector<Vector4> ptlb;
  Vector3 Bthr;


  for (std::vector<Particle>::iterator j = g.begin();
       j != g.end(); j++)
  {
    int mask=0;
    for (int i=0; i<b.relation().nFinalStateParticles(); i++)
    {
      if (b.relation().finalStateParticle(i).pType() ==
          ptype_gamma) 
      {
        if ( (*j).relation().mdstGamma().get_ID() ==
             b.relation().finalStateParticle(i).mdstGamma().get_ID() )
          mask=1;

      }
    }
      Vector4 g_P= (*j).p();
      g_P.boost(P_BEAM.boostVector());
      Vector4 tmpptl(g_P.px(),g_P.py(),g_P.pz(),
                     g_P.e());

    if (mask==1) 
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
    int mask=0;
    for (int i=0; i<b.relation().nFinalStateParticles(); i++)
    {
      if (b.relation().finalStateParticle(i).pType() !=
          ptype_gamma) 
      {
        if ( (*j).relation().mdstCharged().get_ID() ==
             b.relation().finalStateParticle(i).mdstCharged().get_ID() )
          mask=1;
      }
    }

    Vector4 charged_P = (*j).p();
    charged_P.boost(P_BEAM.boostVector());
    Vector4 tmpptl(charged_P.px(),charged_P.py(),charged_P.pz(),
                   charged_P.e());

    if (mask==1) 
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
    int mask=0;
    for (int i=0; i<b.relation().nFinalStateParticles(); i++)
    {
      if (b.relation().finalStateParticle(i).pType() !=
          ptype_gamma) 
      {
        if ( (*j).relation().mdstCharged().get_ID() ==
             b.relation().finalStateParticle(i).mdstCharged().get_ID() )
          mask=1;
      }
    }

    Vector4 charged_P = (*j).p();
    charged_P.boost(P_BEAM.boostVector());
    Vector4 tmpptl(charged_P.px(),charged_P.py(),charged_P.pz(),
                   charged_P.e());

    if (mask==1) 
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

  Vector4 b_P= b.p();
  b_P.boost(P_BEAM.boostVector());

  Vector3 thr_axis = thrust(ptl);
  Vector3 thr_axis_b = thrust(ptlb);

  cos_thr = (thr_axis.dot(thr_axis_b));
  cos_thr = cos_thr / (thr_axis.mag()*thr_axis_b.mag());

//Vector3 tmpjet(b_P.px(),b_P.py(),b_P.pz());
  Vector3 tmpjet(ome_P.px(),ome_P.py(),ome_P.pz());
  spher = Sper(ptl, tmpjet);

}


                      


//for int t--> 0:e 1:mu 2:pi 3:k 4:p 
void ana_d0_20141203::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t) 
     {


        *dr = 1000.0;
        *dz = 1000.0;

        if(charged->trk())
	{
            Mdst_trk_fit& trkFit = charged->trk().mhyp(t);
            if(trkFit)
		{
                HepVector a(5,0);
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

        
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


              


