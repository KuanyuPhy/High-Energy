

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
#include "kfitter/kmassfitter.h"#include "kfitter/kmassvertexfitter.h"
#include "kfitter/khelix2xyz.h"
#include "kfitter/kfitterparticle.h"
#include "helix/Helix.h"
#include "k_sfw.h"

#include "CLHEP/Vector/LorentzVector.h"  // /afs/afs11/belle/belle/b20090127_0910/src/util/belleCLHEP/Vector
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <iostream>
#include "userinfo.h"
#include "panther/panther.h"
#include "ip/IpProfile.h"
#include "koppenburg/pi0eta_prob.h"
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


using namespace std;
//utilities

int checkMultiUse(Particle& i, Particle& j);


// Module class
class ana_psigmag_0809 : public Module 
{
public:
  ana_psigmag_0809(void){};
  ~ana_psigmag_0809(void){};
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


extern "C" Module_descr *mdcl_ana_psigmag_0809()
{
  ana_psigmag_0809 *module = new ana_psigmag_0809;
  Module_descr *dscr = new Module_descr ( "ana_psigmag_0809", module );
  IpProfile::define_global(dscr);
  BeamEnergy::define_global(dscr);
  return dscr;
}



void 
ana_psigmag_0809::hist_def(void)
{ 
  extern BelleTupleManager *BASF_Histogram;  
  BelleTupleManager& tm = *BASF_Histogram;
  
#ifdef BTUPLE

  //uplimit of branch number is 256
//branch name length--> 8 characters(upper limit)
  b_tpl = BASF_Histogram->ntuple("Btopsigmag","mbc modmbc de ebeam EVTID RUNID EXPID eventid farmid \
  pxb pyb pzb eb \
  massp1 pxp1 pyp1 pzp1 ep1 pidp1kp pidp1ppi pidp1pk chrgp1 drp1 dzp1 cosp1 \
  massSigm pxSigm  pySigm  pzSigm  epSigm \
  massLam pxLam  pyLam  pzLam  epLam Lam_dr Lam_dz Lam_zdst Lam_vx Lam_vy Lam_vz Lam_vr chisqL goodLam \
  mass100 px100 py100 pz100 e100 id100kp id100ppi id100kpi id100pk chrg100 dr100 dz100 cos100 \
  mass101 px101 py101 pz101 e101 id101kp id101ppi id101kpi id101pk chrg101 dr101 dz101 cos101 \
  pxGam2  pyGam2  pzGam2  epGam2 pi0veto2 cospig2 \
  pxGam1  pyGam1  pzGam1  epGam1 pi0veto1 pi0prob1 \
  m_ps chisqexk \
  hindex hi100 hi101 hi11 hi2 hi0 mo100 mo101 mo11 mo2 mo0 gmo100 gmo101 gmo11 gmo0 gmo2 ggmo100 ggmo101 ggmo11 gggmo100 gggmo101 gggmo11 \
  charged mgridp mgridn nbnd nbpd bpd1 bpd2 bpd3 bpd4 bnd1 bnd2 bnd3 bnd4 \
  bpnumber bvertexx bvertexy bvertexz overtexx overtexy overtexz \
  pmiss emiss R2 costhr costhp spher cosb vtflag deltaz \
  R2s R1so R2so R3so R4so R1gso R2gso R3gso R4gso R1oo R2oo R3oo R4oo \
  imm2 mmm2 et H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son \
  H0soc H1soc H2soc H3soc \
  H0som H1som H2som H3som H4som H4soc");


#endif

  //
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


void ana_psigmag_0809::init(int*)
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


void ana_psigmag_0809::term (void)
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


void ana_psigmag_0809::begin_run(BelleEvent *evptr, int *status)
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
ana_psigmag_0809::event(BelleEvent *evptr, int *status)
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

  std::vector<Particle> Lam_list;
  std::vector<Particle> Lambar_list;

  std::vector<Particle> Sigma_list;
  std::vector<Particle> Sigma_bar_list;

  std::vector<Particle> Sigma_plus;
  std::vector<Particle> mu_plus;
  std::vector<Particle> mu_minus;
  std::vector<Particle> e_plus;
  std::vector<Particle> e_minus;
  std::vector<Particle> B_cand1;
  std::vector<Particle> D0;
  std::vector<Particle> D0B;


  //define Particle Type
  //reference in(/sw/belle/belle/b20040727_1143/share/data-files/qq98/decay.dec)
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
  
  Ptype ptype_Lambda("LAM");
  Ptype ptype_Lambdabar("ALAM");

  Ptype ptype_Sigma0("SIG0");
  Ptype ptype_Sigma0_bar("ASIG0");


  //  Ptype ptype_Sigma_plus("SIG+");
  //Ptype ptype_Sigma_minus("ASIG+");
  Ptype ptype_mu_plus("MU+");
  Ptype ptype_mu_minus("MU-");
  Ptype ptype_e_plus("E+");
  Ptype ptype_e_minus("E-");
  Ptype ptype_Bplus("B+");
  Ptype ptype_Bminus("B-");
  
  //Declare the likelihood function
  //  atc_pid selKpi(0,1,0,3,2);//for the reprocessed exp7 data
  atc_pid selkpi(3,1,5,3,2); 
  atc_pid selkp(3,1,5,3,4);
  atc_pid selppi(3,1,5,4,2); 
  atc_pid selpk(3,1,5,4,3);
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
      
      //std::cout<<"Mdst charged...."<<std::endl;
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
      
      
      //lepton veto
      if((reject || mu_like < 0.95) && eid_prob < 0.95 )
	{
	  //Selection of K
          if(selkpi.prob(*i)>0.6 && selkp.prob(*i)>0.6)
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
	  //Selection of pi
	  else if (selkpi.prob(*i)<0.4 && selppi.prob(*i)<0.4)
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
	  //Selection of proton
	  else if (selpk.prob(*i)>0.6 && selppi.prob(*i)>0.6)
	    {
	      if ((*i).charge()>0)
		{
		  if (ch.trk().quality() == 0)  // Obtain "Good track"
		    {
		      Particle tmp(*i, ptype_p_plus);
		      p_plus.push_back(tmp);   
		    }
	        }
	      else
		{
		  if (ch.trk().quality() == 0)  // Obtain "Good track"
		    {
		      Particle tmp(*i, ptype_p_minus);
		      p_minus.push_back(tmp);
		    }
              	}	
	    }
	}//end loop of e, mu rejection
      
      
    }
  //end of filling the charged  particles
  
  //filling Lambda list
  //Vee2 selection: kind=1 for k_short, kind=2 for Lambda, kind=3 for Lambdabar
  for(std::vector<Mdst_vee2>::iterator i =Vee2Mgr.begin(); i != Vee2Mgr.end(); i++){
    //std::cout<<"Lambda particle"<<std::endl;
    if((*i).kind()==2){
      float masslam = (*i).p(3)*(*i).p(3)-(*i).p(0)*(*i).p(0)-(*i).p(1)*(*i).p(1)-(*i).p(2)*(*i).p(2);
      if(masslam > 0) masslam = sqrt(masslam);
      int lamflag =0;
      FindLambda find_ilam;
      find_ilam.candidates(*i,IP);
      lamflag = find_ilam.goodLambda();
      if(masslam > 1.05 && masslam < 1.18 && lamflag > 0){
	Particle tmp(*i);
	Lam_list.push_back(tmp);
      }  
    }
    if((*i).kind()==3){
      float masslam = (*i).p(3)*(*i).p(3)-(*i).p(0)*(*i).p(0)-(*i).p(1)*(*i).p(1)-(*i).p(2)*(*i).p(2);
      if(masslam > 0) masslam = sqrt(masslam);
      int lamflag = 0;
      FindLambda find_ilam;
      find_ilam.candidates(*i,IP);
      lamflag = find_ilam.goodLambda();
      if(masslam > 1.05 && masslam < 1.18 && lamflag > 0){
	Particle tmp(*i);
	Lambar_list.push_back(tmp);
      }
    }
  }//end of filling Lambda list
  
  //filling gamma list
  for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin();i != gamma_mag.end();i++){
    Particle tmp(*i);
    //HepLorentzVector Pcm_gamma(tmp.p());
    //HepLorentzVector boost(E_HER*sin(0.022), 0.0, -E_LER+E_HER*cos(0.022), E_HER+E_LER);
    //Pcm_gamma.boost(-boost.boostVector());
    if(tmp.e()>0.01){
      gamma.push_back(tmp);
      //ecm_mdst_gamma->accumulate(Pcm_gamma.e(),1);
    }
  }//end of filling gamma list
  
  
  combination(Sigma_list,ptype_Sigma0,Lam_list,gamma);
  combination(Sigma_bar_list,ptype_Sigma0_bar,Lambar_list,gamma);
  combination(B_cand1,ptype_Bminus,p_minus,Sigma_list,gamma);
  combination(B_cand1,ptype_Bplus,p_plus,Sigma_bar_list,gamma);
  
  //for B candicate 21
  //For B->p Sigma gamma decay
  //B-->p1 Sigm Gam1
  //       Sigm --> Lam Gam2
  //                Lam --> p2 pi
  
  //for checking if we found the Sigma
  /*
  for(std::vector<Particle>::iterator i =Sigma_list.begin(); i != Sigma_list.end(); i++)
    { 
      std::cout<<"Sigma0 start"<<std::endl;
    }
  */

  
  for(std::vector<Particle>::iterator i =B_cand1.begin(); i != B_cand1.end(); i++)
    { 
      
#ifdef BTUPLE
      //------------------- dE and Mbc of B meson -------------------------
      
      //std::cout<<"B ntuple start...."<<std::endl;  
      
      HepLorentzVector b_p((*i).p() );
      
      HepLorentzVector boost_vector(-E_HER*sin(cross_angle), 0.0,
                                    E_LER-E_HER*cos(cross_angle),
                                    E_HER+E_LER);
      b_p.boost( boost_vector.boostVector() );

      double ebeam = BeamEnergy::E_beam_corr();
      double mass_sqr = ebeam*ebeam - b_p.vect().mag2();
      double mass  = (mass_sqr > 0.) ? sqrt(mass_sqr) :  -sqrt(-mass_sqr);
      float de = b_p.e()-ebeam;


      if (mass < 5.2 || de < -0.5 || de > 0.5) continue;

      
      b_tpl->column("mbc", mass);
      b_tpl->column("de", de);
      b_tpl->column("ebeam", BeamEnergy::E_beam_corr());
      
      b_tpl->column("EVTID", Evt);
      b_tpl->column("RUNID", Run);
      b_tpl->column("EXPID", Exp);
      b_tpl->column("eventid", Eventid);
      b_tpl->column("farmid", Farmid);


//------------------- information of B (*i) -------------------------

      b_tpl->column("pxb", (*i).px());
      b_tpl->column("pyb", (*i).py());
      b_tpl->column("pzb", (*i).pz());
      b_tpl->column("eb", (*i).e());
      
//------------------- information of p+ (p1) (*i).child(0) -----------------------
      b_tpl->column("massp1", (*i).child(0).p().mag());
      b_tpl->column("pxp1", (*i).child(0).px());
      b_tpl->column("pyp1", (*i).child(0).py());
      b_tpl->column("pzp1", (*i).child(0).pz());
      b_tpl->column("ep1", (*i).child(0).e());

      //boosted to CMS frame
      HepLorentzVector p1_p((*i).child(0).p() );
      p1_p.boost( boost_vector.boostVector() );

      
      const Mdst_charged* chp1 = &(*i).child(0).mdstCharged();
      float kpi_idp1_kp = selkp.prob(chp1);
      b_tpl->column("pidp1kp", kpi_idp1_kp);
      
      float kpi_idp1_ppi = selppi.prob(chp1);
      b_tpl->column("pidp1ppi", kpi_idp1_ppi);
      
      float kpi_idp1_pk = selpk.prob(chp1);
      b_tpl->column("pidp1pk", kpi_idp1_pk);

      
      float chrgp1 = chp1->charge();
      b_tpl->column("chrgp1", chrgp1);
      
      double drp1,dzp1;
      GetImpactParameters(chp1,&drp1,&dzp1,3);
      b_tpl->column("drp1", float(drp1));
      b_tpl->column("dzp1", float(dzp1));
      
      double cos_polar_p1;
      Vector3 chp1_momentum((*i).child(0).px(), (*i).child(0).py(), (*i).child(0).pz());
      Vector3 P_BEAM_polarp1 ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));
      
      cos_polar_p1=chp1_momentum.dot(P_BEAM_polarp1);
      cos_polar_p1=cos_polar_p1/(chp1_momentum.mag()*P_BEAM_polarp1.mag());
      b_tpl->column("cosp1", cos_polar_p1);
      
      //------------------- information of Sigma0 (*i).child(1) -------------------------
      b_tpl->column("massSigm", (*i).child(1).p().mag());
      b_tpl->column("pxSigm", (*i).child(1).px());
      b_tpl->column("pySigm", (*i).child(1).py());
      b_tpl->column("pzSigm", (*i).child(1).pz());
      b_tpl->column("epSigm", (*i).child(1).e());
      
      
      //------------------- information of Lambda0 from Sigma (*i).child(1).child(0) ---------------
      b_tpl->column("massLam", (*i).child(1).child(0).p().mag());
      b_tpl->column("pxLam", (*i).child(1).child(0).px());
      b_tpl->column("pyLam", (*i).child(1).child(0).py());
      b_tpl->column("pzLam", (*i).child(1).child(0).pz());
      b_tpl->column("epLam", (*i).child(1).child(0).e());

      //boosted to CMS frame
      HepLorentzVector Lam_p((*i).child(1).child(0).p() );
      Lam_p.boost( boost_vector.boostVector() );

      
      float fd = sqrt(((*i).child(1).child(0).mdstVee2().v(0)-IP.x())*((*i).child(1).child(0).mdstVee2().v(0)-IP.x())+
		      ((*i).child(1).child(0).mdstVee2().v(1)-IP.y())*((*i).child(1).child(0).mdstVee2().v(1)-IP.y()));                                  
      b_tpl->column("Lam_dr",fd);           
      fd = abs((*i).child(1).child(0).mdstVee2().v(2)-IP.z());
      b_tpl->column("Lam_dz",fd);
      b_tpl->column("Lam_zdst", (*i).child(1).child(0).mdstVee2().z_dist());
      b_tpl->column("Lam_vx", (*i).child(1).child(0).mdstVee2().v(0));
      b_tpl->column("Lam_vy", (*i).child(1).child(0).mdstVee2().v(1));
      b_tpl->column("Lam_vz", (*i).child(1).child(0).mdstVee2().v(2));  
      b_tpl->column("Lam_vr", sqrt((*i).child(1).child(0).mdstVee2().v(0)*(*i).child(1).child(0).mdstVee2().v(0) +(*i).child(1).child(0).mdstVee2().v(1)*(*i).child(1).child(0).mdstVee2().v(1)) ); 
      b_tpl->column("chisqL", (*i).child(1).child(0).mdstVee2().chisq());                                                                                
      FindLambda find_ilam;
      find_ilam.candidates((*i).child(1).child(0).mdstVee2(),IP);
      b_tpl->column("goodLam",find_ilam.goodLambda());
      
      
      //------------------- information of pi+ from Lambda-bar or p+ from Lambda (*i).child(1).child(0).child(0) -------
      b_tpl->column("mass100", (*i).child(1).child(0).child(0).p().mag());
      b_tpl->column("px100", (*i).child(1).child(0).child(0).px());
      b_tpl->column("py100", (*i).child(1).child(0).child(0).py());
      b_tpl->column("pz100", (*i).child(1).child(0).child(0).pz());
      b_tpl->column("e100", (*i).child(1).child(0).child(0).e());
      
      const Mdst_charged* chp2 = &(*i).child(1).child(0).child(0).mdstCharged();
      float kpi_idp2_kp = selkp.prob(chp2);
      b_tpl->column("id100kp", kpi_idp2_kp);
           
      float kpi_idp2_ppi = selppi.prob(chp2);
      b_tpl->column("id100ppi", kpi_idp2_ppi);
      
      float kpi_idp2_kpi = selkpi.prob(chp2);
      b_tpl->column("id100kpi", kpi_idp2_kpi);

      float kpi_idp2_pk = selpk.prob(chp2);
      b_tpl->column("id100pk", kpi_idp2_pk);
      

      
      float chrgp2 = chp2->charge();
      b_tpl->column("chrg100", chrgp2);
      
      double drp2,dzp2;
      GetImpactParameters(chp2,&drp2,&dzp2,3);
      b_tpl->column("dr100", float(drp2));
      b_tpl->column("dz100", float(dzp2));
      
      double cos_polar_p2;
      Vector3 chp2_momentum((*i).child(1).child(0).child(0).px(), (*i).child(1).child(0).child(0).py(), (*i).child(1).child(0).child(0).pz());
      Vector3 P_BEAM_polarp2 ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));
      
      cos_polar_p2=chp2_momentum.dot(P_BEAM_polarp2);
      cos_polar_p2=cos_polar_p2/(chp2_momentum.mag()*P_BEAM_polarp2.mag());
      b_tpl->column("cos100", cos_polar_p2);
      
      //-------------------- information of p- from Lambda-bar or pi- from Lambda (*i).child(1).child(0).child(1) -----------
      b_tpl->column("mass101", (*i).child(1).child(0).child(1).p().mag());
      b_tpl->column("px101", (*i).child(1).child(0).child(1).px());
      b_tpl->column("py101", (*i).child(1).child(0).child(1).py());
      b_tpl->column("pz101", (*i).child(1).child(0).child(1).pz());
      b_tpl->column("e101", (*i).child(1).child(0).child(1).e());
      
      const Mdst_charged* chpi = &(*i).child(1).child(0).child(1).mdstCharged();
      float kpi_idpi = selkpi.prob(chpi);
      b_tpl->column("id101kpi", kpi_idpi);
      
      float ppi_idpi = selppi.prob(chpi);
      b_tpl->column("id101ppi", ppi_idpi);

      float pk_idpi = selpk.prob(chpi);
      b_tpl->column("id101pk", pk_idpi);

      float kp_idpi = selkp.prob(chpi);
      b_tpl->column("id101kp", kp_idpi);

      
      
      float chrgpi = chpi->charge();
      b_tpl->column("chrg101", chrgpi);
      
      double drpi,dzpi;
      GetImpactParameters(chpi,&drpi,&dzpi,2);
      b_tpl->column("dr101", float(drpi));
      b_tpl->column("dz101", float(dzpi));
      
      double cos_polar_pi;
      Vector3 chpi_momentum((*i).child(1).child(0).child(1).px(), (*i).child(1).child(0).child(1).py(), (*i).child(1).child(0).child(1).pz());
      
      cos_polar_pi=chpi_momentum.dot(P_BEAM_polarp1);
      cos_polar_pi=cos_polar_pi/(chpi_momentum.mag()*P_BEAM_polarp1.mag());
      b_tpl->column("cos101", cos_polar_pi);
      
      
      //------------------- information of Gamma from Sigma (*i).child(1).child(1) -----------------
      //b_tpl->column("massGam2", (*i).child(1).child(1).p().mag());
      b_tpl->column("pxGam2", (*i).child(1).child(1).px());
      b_tpl->column("pyGam2", (*i).child(1).child(1).py());
      b_tpl->column("pzGam2", (*i).child(1).child(1).pz());
      b_tpl->column("epGam2", (*i).child(1).child(1).e());

      double pigam2px = (*i).child(1).child(1).px();
      double pigam2py = (*i).child(1).child(1).py();
      double pigam2pz = (*i).child(1).child(1).pz();
      double psqrtg = pigam2px*pigam2px + pigam2py*pigam2py + pigam2pz*pigam2pz;
      double cospig2 = pigam2pz/sqrt(psqrtg);

      b_tpl->column("cospig2", cospig2);

      //boosted to CMS frame
      HepLorentzVector Gam2_p((*i).child(1).child(1).p() );
      float piveto2;
      piveto2=0;
      //pizero veto loop for all possible gamma     
      for(std::vector<Particle>::iterator k = gamma.begin();k != gamma.end();k++){
        Particle tmp(*k);
        //HepLorentzVector boost(E_HER*sin(0.022), 0.0, -E_LER+E_HER*cos(0.022), E_HER+E_LER);
        //Pcm_gamma.boost(-boost.boostVector());
        //When the gamma energy is larger than 30 MeV
	if(tmp.e()>0.03){
          HepLorentzVector P_gamma(tmp.p());
          HepLorentzVector Tpi0_p(Gam2_p+P_gamma);
          double pimass_sq = Tpi0_p.vect().mag2();
          double piinvmass = (pimass_sq > 0.) ? sqrt(pimass_sq) :  -sqrt(-pimass_sq);
          //pi0 mass 134.97 MeV +-18MeV (116.9, 153 MeV)                                                                    
          if (piinvmass<0.1530 && piinvmass>0.1169){
            piveto2=1.0;}
          //ecm_mdst_gamma->accumulate(Pcm_gamma.e(),1);  
        }
      }//end of filling gamma list                                                                                          
      b_tpl->column("pi0veto2",piveto2);
      Gam2_p.boost( boost_vector.boostVector() );

      
      //------------------- information of psedoparticle that combinated proton and sigma0 --------
      HepLorentzVector p_0((*i).child(0).p());  
      HepLorentzVector p_1((*i).child(1).p());        
      HepLorentzVector p_dib(p_0 + p_1);	
      b_tpl->column("m_ps",p_dib.mag()); // Save the magnitude of the 4 monentum -->mass 
      
      //------------------- information of Gamma (*i).child(2) -----------------
      //b_tpl->column("massGam1", (*i).child(2).p().mag());
      b_tpl->column("pxGam1", (*i).child(2).px());
      b_tpl->column("pyGam1", (*i).child(2).py());
      b_tpl->column("pzGam1", (*i).child(2).pz());
      b_tpl->column("epGam1", (*i).child(2).e());
      

      //boosted to CMS frame
      HepLorentzVector Gam1_p((*i).child(2).p() );
      //Gam1_p.boost( boost_vector.boostVector() );

      float piveto1;
      double cos_polar_g;
      piveto1=0;
      double maxprob = 0.0;
      double pi0prob;
      //pizero veto loop for all possible gamma     
      for(std::vector<Particle>::iterator k = gamma.begin();k != gamma.end();k++){
        Particle tmp(*k);
        //HepLorentzVector boost(E_HER*sin(0.022), 0.0, -E_LER+E_HER*cos(0.022), E_HER+E_LER);
        //Pcm_gamma.boost(-boost.boostVector());
        //When the gamma energy is larger than 30 MeV
	if(tmp.e()>0.03){
          HepLorentzVector P_gamma(tmp.p());
          HepLorentzVector Tpi0_p(Gam1_p+P_gamma);
          double pimass_sq = Tpi0_p.vect().mag2();
          double piinvmass = (pimass_sq > 0.) ? sqrt(pimass_sq) :  -sqrt(-pimass_sq);
          //pi0 mass 134.97 MeV +-18MeV (116.9, 153 MeV)                                                                    
          if (piinvmass<0.1530 && piinvmass>0.1169){
            piveto1=1.0;}
          //ecm_mdst_gamma->accumulate(Pcm_gamma.e(),1);  
	  Vector3 chg_momentum((*k).px(), (*k).py(), (*k).pz());
	  
	  cos_polar_g = chg_momentum.dot(P_BEAM_polarp1);
	  cos_polar_g  = cos_polar_g/(chg_momentum.mag()*P_BEAM_polarp1.mag());
	  double polar_ang = acos(cos_polar_g);
	  pi0prob = Pi0_Prob(piinvmass, tmp.e(), polar_ang);
	  if (maxprob < pi0prob) maxprob = pi0prob;
	  
        }
      }//end of filling gamma list                                                                                          
      b_tpl->column("pi0veto1",piveto1);
      b_tpl->column("pi0prob1", maxprob);
      Gam1_p.boost( boost_vector.boostVector() );
	

      //------------------- Modyfied Mbc calculation --------------------------
      double E_diff, one_epGam12, one_epGam2, scpx, scpy, scpz, sce;
      //std::cout<<"Modyfied Mbc"<<std::endl;
      E_diff = ebeam - p1_p.e() - Lam_p.e();
      one_epGam12 = E_diff/(Gam1_p.e()+Gam2_p.e());
      //one_epGam2 = E_diff/Gam2_p.e();
      
      scpx= Gam1_p.px()*one_epGam12;
      scpy= Gam1_p.py()*one_epGam12;
      scpz= Gam1_p.pz()*one_epGam12;
      sce = Gam1_p.e()*one_epGam12;
      Gam1_p.setPx(scpx);
      Gam1_p.setPy(scpy);
      Gam1_p.setPz(scpz);
      Gam1_p.setE(sce);
      scpx= Gam2_p.px()*one_epGam12;
      scpy= Gam2_p.py()*one_epGam12;
      scpz= Gam2_p.pz()*one_epGam12;
      sce = Gam2_p.e()*one_epGam12;
      Gam2_p.setPx(scpx);
      Gam2_p.setPy(scpy);
      Gam2_p.setPz(scpz);
      Gam2_p.setE(sce);

      
      HepLorentzVector TGam_p(Gam1_p+Gam2_p);
      HepLorentzVector b_pMody(p1_p +Lam_p+ TGam_p);
      
      double modmass_sqr = ebeam*ebeam - b_pMody.vect().mag2();
      double modmass  = (modmass_sqr > 0.) ? sqrt(modmass_sqr) :  -sqrt(-modmass_sqr);

      b_tpl->column("modmbc", modmass);

      //------------------- B vertex fit with ExKfitter  -------------------------
      //std::cout<<"vertex fitter"<<std::endl;
      Mdst_charged Mdst_100 =(*i).child(1).child(0).child(0).mdstCharged(); //p+ from Lambda-bar or pi+ from Lambda
      Mdst_charged Mdst_101 =(*i).child(1).child(0).child(1).mdstCharged(); //pi- from Lambda-bar or p- from Lambda
      Mdst_charged Mdst_0   =(*i).child(0).mdstCharged(); //proton
      

      int pid_100=0;
      int pid_101=0;
      if ((*i).child(1).child(0).child(0).p().mag()>0.90 && (*i).child(1).child(0).child(1).p().mag()>0.10 && (*i).child(1).child(0).child(1).p().mag() <0.20 ){
	pid_100=4;
	pid_101=2;
      }
      else if ((*i).child(1).child(0).child(0).p().mag()>0.10 && (*i).child(1).child(0).child(0).p().mag()<0.20 && (*i).child(1).child(0).child(1).p().mag() <0.90 ){
	pid_100=2;
	pid_101=4; 
      }

      ExKFitterParticle KF_100(Mdst_100, pid_100); 
      ExKFitterParticle KF_101(Mdst_101, pid_101);
      ExKFitterParticle KF_0(Mdst_0, 4); //proton
      
      HepPoint3D Sigma0_init;
      

      //Don't use gamma to reconstruct the vertex, as the uncertainty is large
      
      Sigma0_init.setX(IP.x() + (*i).child(1).p().px()/(*i).child(1).p().rho());
      Sigma0_init.setY(IP.y() + (*i).child(1).p().py()/(*i).child(1).p().rho());
      Sigma0_init.setZ(IP.z() + (*i).child(1).p().pz()/(*i).child(1).p().rho());
      ExKFitterVertex Sigma0_Vertex(Sigma0_init);
      
      HepPoint3D Lambda0_init;
      Lambda0_init.setX((*i).child(1).child(0).mdstVee2().v(0));
      Lambda0_init.setY((*i).child(1).child(0).mdstVee2().v(1));
      Lambda0_init.setZ((*i).child(1).child(0).mdstVee2().v(2));
      ExKFitterVertex Lambda0_Vertex(Lambda0_init);
      
      
      ExKFitterVertex B_Vertex(IP,IPerr);
      
      ExKFitterParticle Lambda0;
      Lambda0.LinkParticle(&KF_100);
      Lambda0.LinkParticle(&KF_101);
      Lambda0.LinkVertex(&Lambda0_Vertex);
      ExKFitterConstrain con1;
      con1.SetVertexConstrain();
      con1.LinkParticle(&KF_100);
      con1.LinkParticle(&KF_101);
      con1.LinkVertex(&Lambda0_Vertex);
      
      ExKFitterParticle B;
      //B.LinkParticle(&Sigma0);
      B.LinkParticle(&Lambda0);
      B.LinkParticle(&KF_0);
      //B.LinkParticle(&KF_2);
      B.LinkVertex(&B_Vertex);
      ExKFitterConstrain con3;
      con3.SetVertexConstrain();
      //con3.LinkParticle(&Sigma0);
      con3.LinkParticle(&Lambda0);
      con3.LinkParticle(&KF_0);
      //con3.LinkParticle(&KF_2);
      con3.LinkVertex(&B_Vertex);
      
      ExKFitter Core;
      Core.LinkConstrain(&con1);
      //Core.LinkConstrain(&con2);
      Core.LinkConstrain(&con3);
      
      int ret = Core.Minimize();
      float chisqExK = Core.Chisq();
      float dof_exk = Core.N_DegreeOfFreedom();
      
      HepPoint3D bvertex = B_Vertex.Vertex();
      if(ret==0)
        {
	  B.Update();
        }
      
      b_tpl->column("chisqexk",chisqExK/dof_exk);
      
      //------------------------- MC truth matching ------------------------
      
      int hindex = 0;
      //std::cout<<"MC truth matching..."<<std::endl;
      if (MCstatus == 1){
	const Mdst_charged ch_100 =(*i).child(1).child(0).child(0).mdstCharged(); //pi+ from Lambda-bar or p+ from Lambda
	const Mdst_charged ch_101 =(*i).child(1).child(0).child(1).mdstCharged(); //p- from Lambda-bar or pi- from Lambda
	const Mdst_gamma ch_11  =(*i).child(1).child(1).mdstGamma(); //gamma from sigma
	const Mdst_gamma ch_2   =(*i).child(2).mdstGamma(); //gamma
	const Mdst_charged ch_0   =(*i).child(0).mdstCharged(); //proton
	
	
        Gen_hepevt Evtch_100=get_hepevt(ch_100);
        Gen_hepevt Evtch_101=get_hepevt(ch_101);
	Gen_hepevt Evtch_11=get_hepevt(ch_11);
	Gen_hepevt Evtch_2=get_hepevt(ch_2);
	Gen_hepevt Evtch_0=get_hepevt(ch_0);

        Gen_hepevt EvtP_100;
        Gen_hepevt EvtP_101;
	Gen_hepevt EvtP_11;
	Gen_hepevt EvtP_2;
        Gen_hepevt EvtP_0;
	
	Gen_hepevt EvtGP_100;
	Gen_hepevt EvtGP_101;
	Gen_hepevt EvtGP_11;
	Gen_hepevt EvtGP_0;
	Gen_hepevt EvtGP_2;

	Gen_hepevt EvtGGP_100;
	Gen_hepevt EvtGGP_101;
	Gen_hepevt EvtGGP_11;
	Gen_hepevt EvtGGP_2;

	Gen_hepevt EvtGGGP_100;
	Gen_hepevt EvtGGGP_101;
	Gen_hepevt EvtGGGP_11;

	Gen_hepevt EvtGGGGP_11;

	b_tpl->column("hi100", Evtch_100.idhep());
	b_tpl->column("hi101", Evtch_101.idhep());
	b_tpl->column("hi11", Evtch_11.idhep());
	b_tpl->column("hi2", Evtch_2.idhep());
	b_tpl->column("hi0", Evtch_0.idhep());

	int mo11, gmo11, ggmo11, gggmo11, ggggmo11;
	int mo2, gmo2, ggmo2;
	
        if (Evtch_100.mo(0))
        {
	  //parent of pi+, should be Lambda-bar or p+ from Lambda
	  EvtP_100=gen_mgr[Evtch_100.mo(0)-1]; //-1 because of index of fotran
	  b_tpl->column("mo100", EvtP_100.idhep());
	  if (EvtP_100.mo(0))
	    {
	      //grandparent of pi+, should be sigma0
	      EvtGP_100=gen_mgr[EvtP_100.mo(0)-1];
	      b_tpl->column("gmo100", EvtGP_100.idhep());
	    
	      if (EvtGP_100.mo(0))
		{
		  //three level up of pi+, shold be the psudoparticle
		  EvtGGP_100=gen_mgr[EvtGP_100.mo(0)-1];
		  b_tpl->column("ggmo100", EvtGGP_100.idhep());
		  if (EvtGGP_100.mo(0))
		    {
		      //four level up of pi+, shold be B+-
		      EvtGGGP_100=gen_mgr[EvtGGP_100.mo(0)-1];
		      b_tpl->column("gggmo100", EvtGGGP_100.idhep());
		    }
		}
	    }
        }
	
        if (Evtch_101.mo(0))
        {
	  //parent of p-, should be Lambda, or parent of pi- from Lambda-bar
	  EvtP_101=gen_mgr[Evtch_101.mo(0)-1]; //-1 because of index of fotran
	  b_tpl->column("mo101", EvtP_101.idhep());
	  if (EvtP_101.mo(0))
	    {
	      //grandparent of p-, should be sigma0
	      EvtGP_101=gen_mgr[EvtP_101.mo(0)-1];
	      b_tpl->column("gmo101", EvtGP_101.idhep());
	    
	      if (EvtGP_101.mo(0))
		{
		  //three level up of p-, shold be the psudoparticle
		  EvtGGP_101=gen_mgr[EvtGP_101.mo(0)-1];
		  b_tpl->column("ggmo101", EvtGGP_101.idhep());
		  if (EvtGGP_101.mo(0))
		    {
		      //four level up of p-, shold be B+-
		      EvtGGGP_101=gen_mgr[EvtGGP_101.mo(0)-1];
		      b_tpl->column("gggmo101", EvtGGGP_101.idhep());
		    }
		}
	    }
        }
	
        if (Evtch_11.mo(0))
	  {
	    //parent of gamma2, should be sigma0
            EvtP_11=gen_mgr[Evtch_11.mo(0)-1];
	    b_tpl->column("mo11", EvtP_11.idhep());
	    mo11 = EvtP_11.idhep();
	    if (EvtP_11.mo(0))
	      {
		//grandparent of gamma2, should be psudoparticle
                EvtGP_11=gen_mgr[EvtP_11.mo(0)-1];
		b_tpl->column("gmo11", EvtGP_11.idhep());
		gmo11 = EvtGP_11.idhep();
		if (EvtGP_11.mo(0))
		  {
		    //three level up than gamma2, shold be B+-
		    EvtGGP_11=gen_mgr[EvtGP_11.mo(0)-1];
		    b_tpl->column("ggmo11", EvtGGP_11.idhep());
		    ggmo11 =  EvtGGP_11.idhep();
		    if (EvtGGP_11.mo(0))
		      {
			//three level up than gamma2, shold be B+-
			EvtGGGP_11=gen_mgr[EvtGGP_11.mo(0)-1];
			b_tpl->column("gggmo11", EvtGGGP_11.idhep());
			gggmo11 =  EvtGGGP_11.idhep();
			if (EvtGGGP_11.mo(0))
			  {
			    //three level up than gamma2, shold be B+-
			    EvtGGGGP_11=gen_mgr[EvtGGGP_11.mo(0)-1];
			    //b_tpl->column("gggmo11", EvtGGGP_11.idhep());
			    ggggmo11 =  EvtGGGGP_11.idhep();
			
			  }
		      }
		  }
	      }
	  }
	
	if (Evtch_2.mo(0))
	  {
	    //parent of gamma 1, should be B+-
            EvtP_2=gen_mgr[Evtch_2.mo(0)-1];
            b_tpl->column("mo2", EvtP_2.idhep());
	    mo2 = EvtP_2.idhep();
	    if (EvtP_2.mo(0))
	      {
		//three level up than gamma2, shold be B+-
		EvtGP_2=gen_mgr[EvtP_2.mo(0)-1];
		b_tpl->column("gmo2", EvtGP_2.idhep());
		gmo2 = EvtGP_2.idhep();
		if (EvtGP_2.mo(0))
		  {
		    //four level up than gamma2, shold be B+-
		    EvtGGP_2=gen_mgr[EvtGP_2.mo(0)-1];
		    //b_tpl->column("gmo2", EvtGP_2.idhep());
		    ggmo2 = EvtGP_2.idhep();
		  }
	      }
	  }
	
	if (Evtch_0.mo(0))
	  {
	    //parent of proton 1, should be psudoparticle
            EvtP_0=gen_mgr[Evtch_0.mo(0)-1];
            b_tpl->column("mo0", EvtP_0.idhep());
	    if (EvtP_0.mo(0))
	      {
		//grandparent of proton 1, should be B+-
                EvtGP_0=gen_mgr[EvtP_0.mo(0)-1];
                b_tpl->column("gmo0", EvtGP_0.idhep());
	      }
	  }
	

	//B^{+} decay, Lambda-bar decay to p+ (100) and pi- (101)
        if(Evtch_100.idhep()== 211 && EvtP_100.idhep()== -3122 && EvtGP_100.idhep()== -3212 && EvtGGP_100.idhep()== 426 && EvtGGGP_100.idhep()== 521
	   && Evtch_101.idhep()== -2212 && EvtP_101.idhep()== -3122 && EvtGP_101.idhep()== -3212 && EvtGGP_101.idhep()== 426 && EvtGGGP_101.idhep()== 521
	   && Evtch_0.idhep()== 2212 && EvtP_0.idhep()== 426 && EvtGP_0.idhep()== 521 ){
	  //list all possible situations for gamma ray id
	  if ( ((ggmo11==521 || gggmo11==521|| ggggmo11==521)&& (mo11==-3212|| gmo11==-3212|| ggmo11==-3212) && 
		(abs(mo11)!=3122 && abs(gmo11)!=3122 && abs(ggmo11)!=3122 && abs(mo11)!=2212 && abs(gmo11)!=2212 && abs(ggmo11)!=2212))
	       &&( mo2==521 || gmo2==521 || ggmo2==521)
	       ){ hindex = 1;}
	}
	//B^{-} decay, Lambda decya to p-(101) and pi+ (100)
        else if(Evtch_100.idhep()== 2212 && EvtP_100.idhep()== 3122 && EvtGP_100.idhep()== 3212 && EvtGGP_100.idhep()== -426 && EvtGGGP_100.idhep()== -521
		&& Evtch_101.idhep()== -211 && EvtP_101.idhep()== 3122 && EvtGP_101.idhep()== 3212 && EvtGGP_101.idhep()== -426 && EvtGGGP_101.idhep()== -521
		&& Evtch_0.idhep()== -2212 && EvtP_0.idhep()== -426 && EvtGP_0.idhep()== -521 ){
	  if( ((ggmo11== -521 || gggmo11== -521|| ggggmo11== -521)&&(mo11== 3212||gmo11== 3212||ggmo11== 3212) && 
	       (abs(mo11)!=3122 && abs(gmo11)!=3122 && abs(ggmo11)!=3122 && abs(mo11)!=2212 && abs(gmo11)!=2212 && abs(ggmo11)!=2212))
	       &&( mo2== -521 || gmo2== -521 || ggmo2== -521)
	     ){ hindex = 1;}
	}
        else
	  { hindex = 0;}

	b_tpl->column("hindex", hindex);	
	//--------------------- trace from top to down (from B meson)------------------------------------
        int nbpd=0,nbnd=0;
	//std::cout<<"Trace from top to down...."<<std::endl;
	for (std::vector<Gen_hepevt>::iterator k = gen_mgr.begin();k != gen_mgr.end(); k++)
	  {
	    int mgr_id;
	    const Gen_hepevt &h_Brec = (*k);
	    //const Gen_hepevt &h_Brec = (*i).genHepevt();
	    
	    mgr_id = h_Brec ? h_Brec.idhep(): 0;
	    
	    //b_tpl->column("mgrid", mgr_id);
	    if (mgr_id == 521||mgr_id == 511)
	      {
		int da1=h_Brec.da(0),da2=h_Brec.da(1);

		if (mgr_id == 521) b_tpl->column("charged",1);
		else b_tpl->column("charged",0);
		nbpd=(da2-da1+1);
		b_tpl->column("nbpd",nbpd); 
		b_tpl->column("mgridp",mgr_id);
		for(int start=0;start<(da2-da1+1);start++)
		  {                                                                                     
		    Gen_hepevt Evda=gen_mgr[da1-1+start];
		    char bpdnumber[32];                                                                 
		    sprintf (bpdnumber,"%s%d","bpd",start+1);    
		    b_tpl->column(bpdnumber,Evda.idhep());
		  } 
	      }
	    else if (mgr_id == -521||mgr_id == -511)
	      {
		int da1=h_Brec.da(0),da2=h_Brec.da(1);
	    
		if (mgr_id == -521) b_tpl->column("charged",-1);
		else b_tpl->column("charged",0);
		nbnd=(da2-da1+1);
		b_tpl->column("nbnd",nbnd);
		b_tpl->column("mgridn",mgr_id);
		for(int start=0;start<(da2-da1+1);start++)                                               
		  {                                   
		    Gen_hepevt Evda=gen_mgr[da1-1+start]; 
		    char bndnumber[32];           
		    sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....  
		    b_tpl->column(bndnumber,Evda.idhep());
		  }
	      }	
	  }	
      }
      //-----------      continuum suppresion  ------------------                         
      //std::cout<<"continuum suppresion..."<<std::endl;
      int vertexflag;                                                                                    
      HepPoint3D overtex;                                                                                
      
      Vector4 OtherBVector;                                                                                                  
      float R2, spher, cos_thr, cos_thp, par_sfw[13];                                                                        
      
      shape((*i), R2,  spher, cos_thr, cos_thp, par_sfw, OtherBVector, bvertex, vertexflag, overtex);                        
    
      b_tpl->column("bvertexx",bvertex.x());                                                                                 
      b_tpl->column("bvertexy",bvertex.y());                                                                                 
      b_tpl->column("bvertexz",bvertex.z());                                                                                 
      b_tpl->column("overtexx",overtex.x());     
      b_tpl->column("overtexy",overtex.y());                                                                                 
      b_tpl->column("overtexz",overtex.z());                                                                                 
      
      Vector4 ExpectedOtherBVector(0.,0.,0.,ebeam);                                                                          
      Vector4 P_miss = ExpectedOtherBVector - OtherBVector;                                                                  
      b_tpl->column("pmiss",P_miss.rho());                                                                                   
      b_tpl->column("emiss",P_miss.e());                                                                                     
      
      // cosine of B and beam dir.                                                                                           
      float cosb = Vector3(b_p.vect()).dot(P_BEAM.vect());                                                                   
      cosb = cosb/(Vector3(b_p.vect()).mag()*P_BEAM.vect().mag());                                                           
      
      b_tpl->column("R2", R2);                                                                                               
      b_tpl->column("costhr", float(cos_thr));                                                                               
      b_tpl->column("costhp", float(cos_thp));                                                                               
      b_tpl->column("spher",  float(spher));   
      b_tpl->column("cosb", cosb);                                                                                           
      
      b_tpl->column("vtflag",vertexflag);                                                                                    
      if( !vertexflag)                                                                                                       
        {                                                                                                                    
	  double deltaZ=bvertex.z()-overtex.z();
	  b_tpl->column("deltaz",deltaZ);
        }
      
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
      
      //k_sfw variables                               
      k_sfw ksfw_obj((*i));         
      const double ksfw(ksfw_obj.fd());        
      const int iksfw(ksfw_obj.i_mm2());       
      b_tpl->column("imm2",iksfw);   
      double miss(ksfw_obj.mm2());   
      b_tpl->column("mmm2",miss);       
      double et(ksfw_obj.e_t());     
      b_tpl->column("et",et);        
      double H0oo(ksfw_obj.Hoo(0));  
      double H1oo(ksfw_obj.Hoo(1));      
      double H2oo(ksfw_obj.Hoo(2));  
      double H3oo(ksfw_obj.Hoo(3));
      double H4oo(ksfw_obj.Hoo(4));
      b_tpl->column("H0oo",H0oo); 
      b_tpl->column("H1oo",H1oo);
      b_tpl->column("H2oo",H2oo);
      b_tpl->column("H3oo",H3oo);
      b_tpl->column("H4oo",H4oo);
      double H0son(ksfw_obj.Hso_n(0));
      double H1son(ksfw_obj.Hso_n(1));
      double H2son(ksfw_obj.Hso_n(2));
      double H3son(ksfw_obj.Hso_n(3));
      double H4son(ksfw_obj.Hso_n(4));
      b_tpl->column("H0son",H0son);
      b_tpl->column("H1son",H1son);
      b_tpl->column("H2son",H2son);
      b_tpl->column("H3son",H3son);
      b_tpl->column("H4son",H4son);
      double H0soc(ksfw_obj.Hso_c(0));
      double H1soc(ksfw_obj.Hso_c(1));
      double H2soc(ksfw_obj.Hso_c(2));
      double H3soc(ksfw_obj.Hso_c(3));
      double H4soc(ksfw_obj.Hso_c(4));
      b_tpl->column("H0soc",H0soc);
      b_tpl->column("H1soc",H1soc);
      b_tpl->column("H2soc",H2soc);
      b_tpl->column("H3soc",H3soc);
      b_tpl->column("H4soc",H4soc);

      double H0som(ksfw_obj.Hso_m(0));
      double H1som(ksfw_obj.Hso_m(1)); 
      double H2som(ksfw_obj.Hso_m(2));
      double H3som(ksfw_obj.Hso_m(3));
      double H4som(ksfw_obj.Hso_m(4));
      b_tpl->column("H0som",H0som);
      b_tpl->column("H1som",H1som);
      b_tpl->column("H2som",H2som);
      b_tpl->column("H3som",H3som);
      b_tpl->column("H4som",H4som);

      std::cout<<"Dump data...."<<std::endl;
      b_tpl->dumpData(); 
      *status = 1;
#endif
    	
    }//end of B candidates


}//end of events

//void ana_psigmag_0809::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P)
void ana_psigmag_0809::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex)
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
//std::cout << "no of B final particle = " << b.relation().nFinalStateParticles()<<endl;
  
  for( int i=0;i< b.relation().nFinalStateParticles();
       i++)
  {
    candiBfinalParticle.push_back(b.relation().finalStateParticle(i));
    //std::cout << " i = " << i << " mass = " <<  b.relation().finalStateParticle(i).mass() <<  endl;
    Vector4 temp_P4 = b.relation().finalStateParticle(i).p();
    temp_P4.boost(P_BEAM.boostVector());
    candiBfinal.push_back(temp_P4);
  }
  
  //by one daughter from b (now omega)
  //std::cout << "omega final state no ="
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
    //Check if this track is used by the rec. signal B
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag) 
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //std::cout <<" "<<(*it1).relation().mdstCharged().get_ID();
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
    //Check if this track is used by the rec. signal B
    for(std::vector<Particle>::iterator j= candiBfinalParticle.begin();
        j!= candiBfinalParticle.end();j++)
      if(checkMultiUse(*it1,*j)){ notBDauFlag=0; break;}
      if(notBDauFlag) 
      {
      otherBfinal.push_back(temp_P4);
      otherB_P=otherB_P+temp_P4;
      //std::cout <<" "<<(*it1).relation().mdstCharged().get_ID();
      Particle& tmp=*it1;
      ver.push_back(&tmp);
      numver++;
      }  
  }

  vertexflag=ver.fit();
  if(!vertexflag) 
  {
    overtex = ver.vtx();
    //std::cout<<numver<<":"<<ver.vtx()<<" "<<ver.used_particles().size()<<endl;
  }

  //std::cout << endl;
  //std::cout << "non B candi ga_id = ";
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
      //std::cout <<" "<<(*it1).relation().mdstGamma().get_ID();
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

  //std::cout << endl;
  //std::cout << " no of ch+ = "<< char_plus.size();
  //std::cout << " no of ch- = "<< char_minus.size();
  //std::cout << " no of gam = "<< g.size() << endl;

  if((candiBfinal.size()+otherBfinal.size()) != allFinal.size())
    std::cout << candiBfinal.size()<<"+"<<otherBfinal.size()<<" !="<< allFinal.size()<<endl;


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
      //std::cout << "gamma " << j <<" is been masked" << endl;
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
      //std::cout << "track " << j <<" is been masked" << endl;
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
      //std::cout << "track " << j <<" is been masked" << endl;
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
void ana_psigmag_0809::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t) 
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
