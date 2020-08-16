

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
//#define STUPLE
//#define PTUPLE
//#define PPTUPLE

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace std;
//utilities

int checkMultiUse(Particle& i, Particle& j);


// Module class
class ana_psigplusg_1211 : public Module 
{
public:
  ana_psigplusg_1211(void){};
  ~ana_psigplusg_1211(void){};
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
  BelleTuple *b_tpl, *s_tpl, *p_tpl, *pp_tpl;
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


extern "C" Module_descr *mdcl_ana_psigplusg_1211()
{
  ana_psigplusg_1211 *module = new ana_psigplusg_1211;
  Module_descr *dscr = new Module_descr ( "ana_psigplusg_1211", module );
  IpProfile::define_global(dscr);
  BeamEnergy::define_global(dscr);
  return dscr;
}



void 
ana_psigplusg_1211::hist_def(void)
{ 
  extern BelleTupleManager *BASF_Histogram;  
  BelleTupleManager& tm = *BASF_Histogram;
  
#ifdef BTUPLE

  //uplimit of branch number is 256
  //branch name length--> 8 characters(upper limit)
  b_tpl = BASF_Histogram->ntuple("Btopsigplusg","mbc modmbc de ebeam EVTID RUNID EXPID eventid farmid \
  pxb pyb pzb eb \
  massp1 pxp1 pyp1 pzp1 ep1 pidp1kp pidp1ppi pidp1pk chrgp1 drp1 dzp1 cosp1 cmp1 \
  massSigm pxSigm  pySigm  pzSigm  epSigm \
  mass10 px10  py10 pz10 e10 id10kp id10ppi id10kpi id10pk chrg10 dr10 dz10 cos10 lab10 cmp10 \
  mass11 px11 py11 pz11 e11 cos11 cmpi11 cospig0 cospig1 pig0e pig1e \
  pxGam1  pyGam1  pzGam1  epGam1 piveto pi0prob \
  m_ps chisqexk chisqem \
  hindex hi10 hi110 hi111 hi2 hi0 mo10 mo110 mo111 mo2 mo0 gmo10 gmo110 gmo111 gmo2 gmo0 ggmo10 ggmo110 ggmo111 gggmo110 gggmo111 \
  charged mgridp mgridn nbnd nbpd bpd1 bpd2 bpd3 bpd4 bnd1 bnd2 bnd3 bnd4 \
  bpnumber bvertexx bvertexy bvertexz overtexx overtexy overtexz	\
  pmiss emiss R2 costhr costhp spher cosb vtflag deltaz \
  nbpd nbnd bnd1 bnd2 bnd3 bnd4 bnd5 \
  R2s R1so R2so R3so R4so R1gso R2gso R3gso R4gso R1oo R2oo R3oo R4oo \
  imm2 mmm2 et H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son \
  H0soc H1soc H2soc H3soc \
  H0som H1som H2som H3som H4som H4soc");
 
#endif

#ifdef STUPLE

  //uplimit of branch number is 256
  //branch name length--> 8 characters(upper limit)
  s_tpl = BASF_Histogram->ntuple("Sigmatoppi0","EVTID RUNID EXPID mSigm pxSigm pySigm \
  pzSigm epSigm abspi0 egdif shindex");
#endif

#ifdef PTUPLE
  //uplimit of branch number is 256
  //branch name length--> 8 characters(upper limit)
  p_tpl = BASF_Histogram->ntuple("pi0togg","EVTID RUNID EXPID mpi0 pxpi0 pypi0 \
  pzpi0 epi0 cos11 pig0e pig1e abspi0 egdif phindex");
#endif

#ifdef PPTUPLE
  //uplimit of branch number is 256
  //branch name length--> 8 characters(upper limit)
  pp_tpl = BASF_Histogram->ntuple("p_plus","EVTID RUNID EXPID mp pxp pyp \
  pzp ep pmon lpk lppi phindex");
#endif

  //hindex hi10 hi11 hi2 hi0 mo10 mo11 mo2 mo0 gmo10 gmo11 gmo2 gmo0 ggmo10 ggmo11 \
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


void ana_psigplusg_1211::init(int*)
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


void ana_psigplusg_1211::term (void)
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


void ana_psigplusg_1211::begin_run(BelleEvent *evptr, int *status)
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
ana_psigplusg_1211::event(BelleEvent *evptr, int *status)
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
  std::cout<<"event:"<<evtcount<<std::endl;

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
  
  //Ptype ptype_Lambda("LAM");
  //Ptype ptype_Lambdabar("ALAM");

  //Ptype ptype_Sigma0("SIG0");
  //Ptype ptype_Sigma0_bar("ASIG0");


  Ptype ptype_Sigma_plus("SIG+");
  //Ptype ptype_Sigma_minus("ASIG+");
  Ptype ptype_mu_plus("MU+");
  Ptype ptype_mu_minus("MU-");
  Ptype ptype_e_plus("E+");
  Ptype ptype_e_minus("E-");
  Ptype ptype_Bplus("B+");
  Ptype ptype_Bminus("B-");
  Ptype ptype_B0("B0");
  Ptype ptype_antiB0("B0B");
  
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
	  //Selection of proton else //else if (selpk.prob(*i)>0.3 && selppi.prob(*i)>0.3)
	  else  if (selpk.prob(*i)>0.6 && selppi.prob(*i)>0.6)
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

  
  //filling gamma list
  for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin();i != gamma_mag.end();i++){
    Particle tmp(*i);
    //HepLorentzVector Pcm_gamma(tmp.p());
    //HepLorentzVector boost(E_HER*sin(0.022), 0.0, -E_LER+E_HER*cos(0.022), E_HER+E_LER);
    //Pcm_gamma.boost(-boost.boostVector());
    if(tmp.e()>0.04){
      gamma.push_back(tmp);
      //ecm_mdst_gamma->accumulate(Pcm_gamma.e(),1);
    }
  }//end of filling gamma list


  for(std::vector<Mdst_pi0>::iterator i = pi0_mag.begin();i!=pi0_mag.end();i++){

    Particle tmp(*i);
    Particle g1((*i).gamma(0));
    Particle g2((*i).gamma(1));

    float pi0mass =  (*i).mass(); //(*i).p(3)*(*i).p(3)-(*i).p(0)*(*i).p(0)-(*i).p(1)*(*i).p(1)-(*i).p(2)*(*i).p(2);
    float abspi = (*i).p(0)*(*i).p(0)+(*i).p(1)*(*i).p(1)+(*i).p(2)*(*i).p(2);
    float egdif = abs(g1.e()-g2.e())/(g1.e()+g2.e());
    
    if(abspi >0) abspi = sqrt(abspi);
    //abs momentum larger than 100MeV/c
    //if(abspi>0.01 ){
    if(abspi>0.1 ){
      //energy of two daughter gammas should have enegy larger than 50MeV and the diff of two gamma is less than 0.9 of the summation of two energies
      //if (g1.e()>0.02 && g2.e()>0.02){
      if (g1.e()>0.05 && g2.e()>0.05 && egdif <0.9){
	//if (g1.e()>0.05 && g2.e()>0.05 ){
	if (pi0mass>0.118 && pi0mass<0.150){
	  //std::cout<<"pi0 particle"<<std::endl;
	pi_0.push_back(tmp);
	}	
      }
    }
  }
  
  combination(Sigma_plus,ptype_Sigma_plus,p_plus,pi_0);
  //std::cout<<"combine Sigma+"<<std::endl;
  combination(B_cand1,ptype_antiB0,p_minus,Sigma_plus,gamma);
  //std::cout<<"combine B0-bar"<<std::endl;
  
  //for B candicate 21
  //For B->p Sigma gamma decay
  //B-->p- Sigm Gam1
  //       Sigm --> p+ pi0
  //                   pi0 --> gamma gamma
  
  //for checking if we found the Sigma

  //p+ loop
  for(std::vector<Particle>::iterator i =p_plus.begin(); i != p_plus.end(); i++)
    {
      //std::cout<<"proton start"<<std::endl;
      
#ifdef PPTUPLE
      pp_tpl->column("EVTID", Evt);
      pp_tpl->column("RUNID", Run);
      pp_tpl->column("EXPID", Exp);
      //---------p_plus mass
      pp_tpl->column("mp", (*i).p().mag());
      
      pp_tpl->column("pxp", (*i).px());
      pp_tpl->column("pyp", (*i).py());
      pp_tpl->column("pzp", (*i).pz());
      pp_tpl->column("ep", (*i).e());
      
      float pmon, psqre;
      psqre = (*i).px()*(*i).px()+(*i).py()*(*i).py()+(*i).pz()*(*i).pz();
      pmon = sqrt(psqre);

      pp_tpl->column("pmon", pmon);

      const Mdst_charged* chp1 = &(*i).mdstCharged();
      
      float kpi_idp_ppi = selppi.prob(chp1);
      pp_tpl->column("lppi", kpi_idp_ppi);
      
      float kpi_idp_pk = selpk.prob(chp1);
      pp_tpl->column("lpk", kpi_idp_pk);
	
      int phindex=0;
      //p_plus MC truth
      if (MCstatus == 1){
	const Mdst_charged chs_10 =(*i).mdstCharged(); //p+ from Sigma+

	Gen_hepevt Evtchs_10=get_hepevt(chs_10);

	Gen_hepevt EvtPs_10;


	if (Evtchs_10.mo(0))
	  {
	    //parent of p-, should be Sigma+
	    EvtPs_10=gen_mgr[Evtchs_10.mo(0)-1]; //-1 because of index of fotran
	  }
	if(Evtchs_10.idhep()== 2212 && EvtPs_10.idhep()== 3222){
	  phindex=1;}
	else{
	  phindex=0;}
	pp_tpl->column("phindex", phindex);	
      }
      std::cout<<"Dump proton+ data...."<<std::endl;
      pp_tpl->dumpData(); 
#endif
    }

  

  
  //pi0 loop
  for(std::vector<Particle>::iterator i =pi_0.begin(); i != pi_0.end(); i++)
    {
      //std::cout<<"pi zero start"<<std::endl;
      
#ifdef PTUPLE
      p_tpl->column("EVTID", Evt);
      p_tpl->column("RUNID", Run);
      p_tpl->column("EXPID", Exp);
      //---------pi0 mass
      //p_tpl->column("mpi0", (*i).p().mass();

      p_tpl->column("pxpi0", (*i).px());
      p_tpl->column("pypi0", (*i).py());
      p_tpl->column("pzpi0", (*i).pz());
      p_tpl->column("epi0", (*i).e());
      
      //If pi0  
      Mdst_pi0 pi0y = (*i).mdstPi0();
      float mass11 = pi0y.mass();
      p_tpl->column("mpi0", mass11);

      Vector3 chpi0_momentum((*i).px(), (*i).py(), (*i).pz());
      Vector3 P_BEAM_polarpi0 ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));
      double cos_polar_pi0=chpi0_momentum.dot(P_BEAM_polarpi0);
      cos_polar_pi0=cos_polar_pi0/(chpi0_momentum.mag()*P_BEAM_polarpi0.mag());
      //Looking into its two daughters gamma rays
      Mdst_gamma& gam0 = pi0y.gamma(0);
      Mdst_gamma& gam1 = pi0y.gamma(1);
      
      double pigam0px, pigam0py, pigam0pz, cospig0;
      double pigam1px, pigam1py, pigam1pz, cospig1;
      double pig0e, pig1e;

      pigam0px = gam0.px();
      pigam0py = gam0.py();
      pigam0pz = gam0.pz();
      double psqrtg0 = pigam0px*pigam0px + pigam0py*pigam0py + pigam0pz*pigam0pz;
      cospig0 = pigam0pz/sqrt(psqrtg0);
      pig0e = sqrt(psqrtg0);

      pigam1px = gam1.px();
      pigam1py = gam1.py();
      pigam1pz = gam1.pz();
      double psqrtg1 = pigam1px*pigam1px + pigam1py*pigam1py + pigam1pz*pigam1pz;
      cospig1 = pigam1pz/sqrt(psqrtg1);
      pig1e = sqrt(psqrtg1); //gam1.e();

      
      p_tpl->column("cos11", cos_polar_pi0);
      //p_tpl->column("cospig0", cospig0);
      //p_tpl->column("cospig1", cospig1);
      p_tpl->column("pig0e", pig0e);
      p_tpl->column("pig1e", pig1e);

      float gp1E = (*i).child(0).e();
      float gp2E = (*i).child(1).e();
      
      float abspip = (*i).px()*(*i).px()+(*i).py()*(*i).py()+(*i).pz()*(*i).pz();
      //float egdifp = abs(gp1E-gp2E)/(gp1E+gp2E);
      float egdifp = abs(pig0e-pig1e)/(pig1e+pig1e);

      p_tpl->column("abspi0", abspip);
      p_tpl->column("egdif", egdifp);

      int phindex=0;
      std::cout<<"MC truth matching for pi0..."<<std::endl;
      //--------MC truth for sigma and its daughters-----
      if (MCstatus == 1){
	const Mdst_gamma chp_110  =(*i).child(0).mdstGamma(); //gamma from pi0
        const Mdst_gamma chp_111  =(*i).child(1).mdstGamma(); //gamma from pi0
	
	Gen_hepevt Evtchp_110=get_hepevt(chp_110);
	Gen_hepevt Evtchp_111=get_hepevt(chp_111);
	
        Gen_hepevt EvtPp_110;
        Gen_hepevt EvtPp_111;
	
	Gen_hepevt EvtGPp_110;
        Gen_hepevt EvtGPp_111;
	
	Gen_hepevt EvtGGPp_110;
        Gen_hepevt EvtGGPp_111;
	
	Gen_hepevt EvtGGGPp_110;
        Gen_hepevt EvtGGGPp_111;

	Gen_hepevt EvtGGGGPp_110;
        Gen_hepevt EvtGGGGPp_111;
	

	int mop110=0;
	int mop111=0;
	int gmop110=0;
	int gmop111=0;
	int ggmop110=0; 
	int ggmop111=0;
	int gggmop110=0; 
	int gggmop111=0;
	int ggggmop110=0; 
	int ggggmop111=0;

	if (Evtchp_110.mo(0))
	  {
	    //parent of gamma, should be pi0
	    EvtPp_110=gen_mgr[Evtchp_110.mo(0)-1]; //-1 because of index of fotran
	    mop110 = EvtPp_110.idhep();
	    if (EvtPp_110.mo(0))
	      {
		//grandparent of gamma, should be Sigma+
		EvtGPp_110=gen_mgr[EvtPp_110.mo(0)-1];
		gmop110 = EvtGPp_110.idhep();
		if (EvtGPp_110.mo(0))
		  {
		    EvtGGPp_110=gen_mgr[EvtGPp_110.mo(0)-1];
		    ggmop110 = EvtGGPp_110.idhep();
		    if (EvtGGPp_110.mo(0))
		      {
			EvtGGGPp_110=gen_mgr[EvtGGPp_110.mo(0)-1];
			gggmop110 = EvtGGGPp_110.idhep();	
			if (EvtGGGPp_110.mo(0))
			  {
			    EvtGGGGPp_110=gen_mgr[EvtGGGPp_110.mo(0)-1];
			    ggggmop110 = EvtGGGGPp_110.idhep();	
			  }
		      }
		  }
	      }
	  }
	
	if (Evtchp_111.mo(0))
	  {
	    //parent of gamma, should be pi0
	    EvtPp_111=gen_mgr[Evtchp_111.mo(0)-1]; //-1 because of index of fotran
	    mop111 = EvtPp_111.idhep();
	    if (EvtPp_111.mo(0))
	      {
		//grandparent of gamma, should be Sigma+
		EvtGPp_111=gen_mgr[EvtPp_111.mo(0)-1];
		gmop111 = EvtGPp_111.idhep();
		if (EvtGPp_111.mo(0))
		  {
		    //three level up of gamma, shold be the psudoparticle
		    EvtGGPp_111=gen_mgr[EvtGPp_111.mo(0)-1];
		    ggmop111 = EvtGGPp_111.idhep();
		    if (EvtGGPp_111.mo(0))
		      {
			EvtGGGPp_111=gen_mgr[EvtGGPp_111.mo(0)-1];
			gggmop111 = EvtGGGPp_111.idhep();
			if (EvtGGGPp_111.mo(0))
			  {
			    EvtGGGGPp_111=gen_mgr[EvtGGGPp_111.mo(0)-1];
			    ggggmop111 = EvtGGGGPp_111.idhep();	
			  }
		      }
		  }
	      }
	  }

	
	//gamma, gamma from pi0, if any of its parents, graparents, or gra-grapaerents is Sigma
	if(mop110==3222||gmop110==3222||ggmop110==3222||gggmop110==3222){
	  //if(gmop110==-511||ggmop110==-511||gggmop110==-511||ggggmop110==-511){
	  phindex = 1;
	}else {
	  phindex = 0;}
	if(mop111==3222||gmop111==3222||ggmop111==3222||gggmop111==3222){
	  //if(gmop111==-511||ggmop111==-511||gggmop111==-511||ggggmop111==-511){
	  phindex *= 1;
	}else {
	  phindex *= 0;}
	
	p_tpl->column("phindex", phindex);
      }
      std::cout<<"Dump Sigma data...."<<std::endl;
      p_tpl->dumpData(); 
#endif
    }
  

  //Sigma plus loop
  for(std::vector<Particle>::iterator i =Sigma_plus.begin(); i != Sigma_plus.end(); i++)
    { 
      //std::cout<<"Sigma+ start"<<std::endl;
      
#ifdef STUPLE
      s_tpl->column("EVTID", Evt);
      s_tpl->column("RUNID", Run);
      s_tpl->column("EXPID", Exp);
      //---------sigma mass 
      s_tpl->column("mSigm", (*i).p().mag());
      s_tpl->column("pxSigm", (*i).px());
      s_tpl->column("pySigm", (*i).py());
      s_tpl->column("pzSigm", (*i).pz());
      s_tpl->column("epSigm", (*i).e());

      //Particle g1((*i).gamma(0));
      //Particle g2((*i).gamma(1));
      float g1E = (*i).child(1).child(0).e();
      float g2E = (*i).child(1).child(1).e();
      

      //float pi0mass =  (*i).p(3)*(*i).p(3)-(*i).p(0)*(*i).p(0)-(*i).p(1)*(*i).p(1)-(*i).p(2)*(*i).p(2);
      float abspiS = (*i).child(1).px()*(*i).child(1).px()+(*i).child(1).py()*(*i).child(1).py()+(*i).child(1).pz()*(*i).child(1).pz();
      float egdifS = abs(g1E-g2E)/(g1E+g2E);
    
      s_tpl->column("abspi0", abspiS);
      s_tpl->column("egdif", egdifS);

      int shindex = 0;
      std::cout<<"MC truth matching for Sigma+..."<<std::endl;
      //--------MC truth for sigma and its daughters-----
      if (MCstatus == 1){
	const Mdst_charged chs_10 =(*i).child(0).mdstCharged(); //p+ from Sigma+
	const Mdst_gamma chs_110  =(*i).child(1).child(0).mdstGamma(); //gamma from pi0
        const Mdst_gamma chs_111  =(*i).child(1).child(1).mdstGamma(); //gamma from pi0

	Gen_hepevt Evtchs_10=get_hepevt(chs_10);
	Gen_hepevt Evtchs_110=get_hepevt(chs_110);
	Gen_hepevt Evtchs_111=get_hepevt(chs_111);

	Gen_hepevt EvtPs_10;
        Gen_hepevt EvtPs_110;
        Gen_hepevt EvtPs_111;

	Gen_hepevt EvtGPs_110;
        Gen_hepevt EvtGPs_111;

	Gen_hepevt EvtGGPs_110;
        Gen_hepevt EvtGGPs_111;

	Gen_hepevt EvtGGGPs_110;
        Gen_hepevt EvtGGGPs_111;

	int mos110=0;
	int mos111=0;
	int gmos110=0;
	int gmos111=0;
	int ggmos110=0; 
	int ggmos111=0;
	int gggmos110=0; 
	int gggmos111=0;

	// int ismo110, ismo111, isgmo110, isgmo111, isggmo110, isggmo111;
	if (Evtchs_10.mo(0))
	  {
	    //parent of p-, should be Sigma+
	    EvtPs_10=gen_mgr[Evtchs_10.mo(0)-1]; //-1 because of index of fotran
	  }

	if (Evtchs_110.mo(0))
	  {
	    //parent of gamma, should be pi0
	    EvtPs_110=gen_mgr[Evtchs_110.mo(0)-1]; //-1 because of index of fotran
	    mos110 = EvtPs_110.idhep();
	    if (EvtPs_110.mo(0))
	      {
		//grandparent of gamma, should be Sigma+
		EvtGPs_110=gen_mgr[EvtPs_110.mo(0)-1];
		gmos110 = EvtGPs_110.idhep();
		if (EvtGPs_110.mo(0))
		  {
		    EvtGGPs_110=gen_mgr[EvtGPs_110.mo(0)-1];
		    ggmos110 = EvtGGPs_110.idhep();
		    if (EvtGGPs_110.mo(0))
		      {
			EvtGGGPs_110=gen_mgr[EvtGGPs_110.mo(0)-1];
			gggmos110 = EvtGGGPs_110.idhep();
			
		      }
		  }
	      }
	  }
	
	if (Evtchs_111.mo(0))
	  {
	    //parent of gamma, should be pi0
	    EvtPs_111=gen_mgr[Evtchs_111.mo(0)-1]; //-1 because of index of fotran
	    mos111 = EvtPs_111.idhep();
	    if (EvtPs_111.mo(0))
	      {
		//grandparent of gamma, should be Sigma+
		EvtGPs_111=gen_mgr[EvtPs_111.mo(0)-1];
		gmos111 = EvtGPs_111.idhep();
		if (EvtGPs_111.mo(0))
		  {
		    //three level up of gamma, shold be the psudoparticle
		    EvtGGPs_111=gen_mgr[EvtGPs_111.mo(0)-1];
		    ggmos111 = EvtGGPs_111.idhep();
		    if (EvtGGPs_111.mo(0))
		      {
			EvtGGGPs_111=gen_mgr[EvtGGPs_111.mo(0)-1];
			gggmos111 = EvtGGGPs_111.idhep();
		      }
		  }
	      }
	  }

	//Simga+ decay to p+ (10) and pi0->gamma, gamma (110,111)
	// p+
	
        if(Evtchs_10.idhep()== 2212 && EvtPs_10.idhep()== 3222){
	  shindex=1;}
	else{
	  shindex=0;}
	
	//gamma, gamma from pi0, if any of its parents, graparents, or gra-grapaerents is Sigma
	//if(EvtPs_110.idhep()==3222||EvtGPs_110.idhep()==3222||EvtGGPs_110.idhep()==3222){
	if(mos110==3222||gmos110==3222||ggmos110==3222||gggmos110==3222){
	  shindex *= 1;
	}else {
	  shindex *= 0;}
	
	if(mos111==3222||gmos111==3222||ggmos111==3222||gggmos111==3222){
	  shindex *= 1;
	}else {
	  shindex *= 0;}
	
	s_tpl->column("shindex", shindex);
	
      }

      
      std::cout<<"Dump Sigma data...."<<std::endl;
      s_tpl->dumpData(); 
      //*status = 1;

#endif
    }
  

  
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
      
//------------------- information of p- (p1) (*i).child(0) -----------------------
      b_tpl->column("massp1", (*i).child(0).p().mag());
      b_tpl->column("pxp1", (*i).child(0).px());
      b_tpl->column("pyp1", (*i).child(0).py());
      b_tpl->column("pzp1", (*i).child(0).pz());
      b_tpl->column("ep1", (*i).child(0).e());

      
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

      //boosted to CMS frame
      HepLorentzVector p1_p((*i).child(0).p() );
      p1_p.boost( boost_vector.boostVector() );

      double cmp1_sq = p1_p.px()*p1_p.px() + p1_p.py()*p1_p.py() + p1_p.pz()*p1_p.pz();
      
      double cmp1 = (cmp1_sq > 0.) ? sqrt(cmp1_sq) : -sqrt(-cmp1_sq);
      b_tpl->column("cmp1", cmp1); //pi momentum at cms frame
      
      //------------------- information of Sigma+ (*i).child(1) -------------------------
      b_tpl->column("massSigm", (*i).child(1).p().mag());
      b_tpl->column("pxSigm", (*i).child(1).px());
      b_tpl->column("pySigm", (*i).child(1).py());
      b_tpl->column("pzSigm", (*i).child(1).pz());
      b_tpl->column("epSigm", (*i).child(1).e());
      
      
      //------------------- information of p+ from Sigma+  -------
      b_tpl->column("mass10", (*i).child(1).child(0).p().mag());
      b_tpl->column("px10", (*i).child(1).child(0).px());
      b_tpl->column("py10", (*i).child(1).child(0).py());
      b_tpl->column("pz10", (*i).child(1).child(0).pz());
      b_tpl->column("e10", (*i).child(1).child(0).e());
      
      const Mdst_charged* chp2 = &(*i).child(1).child(0).mdstCharged();
      float kpi_idp2_kp = selkp.prob(chp2);
      b_tpl->column("id10kp", kpi_idp2_kp);
           
      float kpi_idp2_ppi = selppi.prob(chp2);
      b_tpl->column("id10ppi", kpi_idp2_ppi);
      
      float kpi_idp2_kpi = selkpi.prob(chp2);
      b_tpl->column("id10kpi", kpi_idp2_kpi);

      float kpi_idp2_pk = selpk.prob(chp2);
      b_tpl->column("id10pk", kpi_idp2_pk);
      
      
      float chrgp2 = chp2->charge();
      b_tpl->column("chrg10", chrgp2);
      
      double drp2,dzp2;
      GetImpactParameters(chp2,&drp2,&dzp2,3);
      b_tpl->column("dr10", float(drp2));
      b_tpl->column("dz10", float(dzp2));
      
      double cos_polar_p2;
      Vector3 chp2_momentum((*i).child(1).child(0).px(), (*i).child(1).child(0).py(), (*i).child(1).child(0).pz());
      Vector3 P_BEAM_polarp2 ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));

      cos_polar_p2=chp2_momentum.dot(P_BEAM_polarp2);
      cos_polar_p2=cos_polar_p2/(chp2_momentum.mag()*P_BEAM_polarp2.mag());
      b_tpl->column("cos10", cos_polar_p2);
      
      //boosted to CMS frame 
      HepLorentzVector p10_p((*i).child(1).child(0).p());
      double lab10_sq = p10_p.px()*p10_p.px() + p10_p.py()*p10_p.py() + p10_p.pz()*p10_p.pz();
      double lab10 = (lab10_sq > 0.) ? sqrt(lab10_sq) : -sqrt(-lab10_sq);
      b_tpl->column("lab10", lab10);

      p10_p.boost( boost_vector.boostVector() );
      double cmp10_sq = p10_p.px()*p10_p.px() + p10_p.py()*p10_p.py() + p10_p.pz()*p10_p.pz();
      
      double cmp10 = (cmp10_sq > 0.) ? sqrt(cmp10_sq) : -sqrt(-cmp10_sq);
      b_tpl->column("cmp10", cmp10); //pi momentum at cms frame

      //------------------- information of pi0 from Sigma+ (*i).child(1).child(1) -----------
      // b_tpl->column("mass11", (*i).child(0).child(1).p().mag());      
      b_tpl->column("px11", (*i).child(1).child(1).px());
      b_tpl->column("py11", (*i).child(1).child(1).py());
      b_tpl->column("pz11", (*i).child(1).child(1).pz());
      b_tpl->column("e11", (*i).child(1).child(1).e());
      double mass11; // 
      
      //If pi0  
      Mdst_pi0 pi0y = (*i).child(1).child(1).mdstPi0();
      mass11 = pi0y.mass();
      Vector3 chpi0_momentum((*i).child(1).child(1).px(), (*i).child(1).child(1).py(), (*i).child(1).child(1).pz());
      Vector3 P_BEAM_polarpi0 ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));
      double cos_polar_pi0=chpi0_momentum.dot(P_BEAM_polarpi0);
      cos_polar_pi0=cos_polar_pi0/(chpi0_momentum.mag()*P_BEAM_polarpi0.mag());
      //Looking into its two daughters gamma rays
      Mdst_gamma& gam0 = pi0y.gamma(0);
      Mdst_gamma& gam1 = pi0y.gamma(1);
      
      double pigam0px, pigam0py, pigam0pz, cospig0;
      double pigam1px, pigam1py, pigam1pz, cospig1;
      double pig0e, pig1e;
      
      pigam0px = gam0.px();
      pigam0py = gam0.py();
      pigam0pz = gam0.pz();
      double psqrtg0 = pigam0px*pigam0px + pigam0py*pigam0py + pigam0pz*pigam0pz;
      cospig0 = pigam0pz/sqrt(psqrtg0);
      pig0e = sqrt(psqrtg0);
      
      pigam1px = gam1.px();
      pigam1py = gam1.py();
      pigam1pz = gam1.pz();
      double psqrtg1 = pigam1px*pigam1px + pigam1py*pigam1py + pigam1pz*pigam1pz;
      cospig1 = pigam1pz/sqrt(psqrtg1);
      pig1e = sqrt(psqrtg1); //gam1.e();
      
      //Filling the variables, here
      b_tpl->column("mass11", mass11);
      b_tpl->column("cos11", cos_polar_pi0);
      b_tpl->column("cospig0", cospig0);
      b_tpl->column("cospig1", cospig1);
      b_tpl->column("pig0e", pig0e);
      b_tpl->column("pig1e", pig1e);
      
      //boosted to CMS frame
      HepLorentzVector pi11_p((*i).child(1).child(1).p());
      pi11_p.boost( boost_vector.boostVector() );
      
      double cmpi11_sq = pi11_p.px()*pi11_p.px() + pi11_p.py()*pi11_p.py() + pi11_p.pz()*pi11_p.pz();
      
      double cmpi11 = (cmpi11_sq > 0.) ? sqrt(cmpi11_sq) : -sqrt(-cmpi11_sq);
      b_tpl->column("cmpi11", cmpi11); //pi momentum at cms frame
      
      HepLorentzVector pig1_p((*i).child(1).child(1).child(0).p());
      HepLorentzVector pig2_p((*i).child(1).child(1).child(1).p());
      
      pig1_p.boost( boost_vector.boostVector() );
      pig2_p.boost( boost_vector.boostVector() );
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

      float piveto;
      piveto=0;
      float maxprob = -10.0;
      float pi0prob, cos_polar_g;
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
          //pi0 mass 134.7 MeV +-18MeV (116.9, 153 MeV)
          if (piinvmass<0.1530 && piinvmass>0.1169){
            piveto=1.0;}
          //ecm_mdst_gamma->accumulate(Pcm_gamma.e(),1);                      


          Vector3 chg_momentum((*k).px(), (*k).py(), (*k).pz());

          cos_polar_g = chg_momentum.dot(P_BEAM_polarp1);
          cos_polar_g  = cos_polar_g/(chg_momentum.mag()*P_BEAM_polarp1.mag());
          double polar_ang = acos(cos_polar_g);
          pi0prob = Pi0_Prob(piinvmass, tmp.e(), polar_ang);
          if (maxprob < pi0prob) maxprob = pi0prob;

        }
      }//end of filling gamma list
      b_tpl->column("piveto",piveto);
      b_tpl->column("pi0prob",maxprob);
      Gam1_p.boost( boost_vector.boostVector() );
      
      //------------------- Modyfied Mbc calculation --------------------------
      double E_diff, one_epGam12, one_epGam2, scpx, scpy, scpz, sce;
      //std::cout<<"Modyfied Mbc"<<std::endl;
      E_diff = ebeam - p1_p.e() - p10_p.e() ;
      one_epGam12 = E_diff/(Gam1_p.e()+pig1_p.e()+pig2_p.e());
      //one_epGam2 = E_diff/Gam2_p.e();
      
      scpx= Gam1_p.px()*one_epGam12;
      scpy= Gam1_p.py()*one_epGam12;
      scpz= Gam1_p.pz()*one_epGam12;
      sce = Gam1_p.e()*one_epGam12;
      Gam1_p.setPx(scpx);
      Gam1_p.setPy(scpy);
      Gam1_p.setPz(scpz);
      Gam1_p.setE(sce);
      scpx= pi11_p.px()*one_epGam12;
      scpy= pi11_p.py()*one_epGam12;
      scpz= pi11_p.pz()*one_epGam12;
      sce = pi11_p.e()*one_epGam12;
      pi11_p.setPx(scpx);
      pi11_p.setPy(scpy);
      pi11_p.setPz(scpz);
      pi11_p.setE(sce);
      
      HepLorentzVector TGam_p(Gam1_p+pi11_p);
      HepLorentzVector b_pMody(p1_p + p10_p + TGam_p);
      
      double modmass_sqr = ebeam*ebeam - b_pMody.vect().mag2();
      double modmass  = (modmass_sqr > 0.) ? sqrt(modmass_sqr) :  -sqrt(-modmass_sqr);

      b_tpl->column("modmbc", modmass);

      
      //------------------- B vertex fit with ExKfitter  -------------------------
      //std::cout<<"vertex fitter"<<std::endl;
      Mdst_charged Mdst_10 = (*i).child(1).child(0).mdstCharged(); //p+ from Sigma+
      //Mdst_charged Mdst_11 =(*i).child(1).child(0).child(1).mdstCharged(); //pi0 from Sigma+
      Mdst_charged Mdst_0 = (*i).child(0).mdstCharged(); //proton
      
      int pid_10=4;
      
      ExKFitterParticle KF_10(Mdst_10, pid_10); //proton
      //ExKFitterParticle KF_101(Mdst_101, pid_101);
      ExKFitterParticle KF_0(Mdst_0, 4); //proton
      
      HepPoint3D Sigma_init;
      
      //Don't use gamma to reconstruct the vertex, as the uncertainty is large
      
      Sigma_init.setX(IP.x() + (*i).child(1).p().px()/(*i).child(1).p().rho());
      Sigma_init.setY(IP.y() + (*i).child(1).p().py()/(*i).child(1).p().rho());
      Sigma_init.setZ(IP.z() + (*i).child(1).p().pz()/(*i).child(1).p().rho());
      ExKFitterVertex Sigma_Vertex(Sigma_init);
      
      ExKFitterVertex B_Vertex(IP,IPerr);
      
      ExKFitterParticle Sigma_plus;
      Sigma_plus.LinkParticle(&KF_10);
      Sigma_plus.LinkVertex(&Sigma_Vertex);
      ExKFitterConstrain con1;
      con1.SetVertexConstrain();
      con1.LinkParticle(&KF_10);
      con1.LinkVertex(&Sigma_Vertex);
      
      ExKFitterParticle B;
      B.LinkParticle(&Sigma_plus);
      B.LinkParticle(&KF_0);
      B.LinkVertex(&B_Vertex);
      ExKFitterConstrain con3;
      con3.SetVertexConstrain();
      con3.LinkParticle(&Sigma_plus);
      con3.LinkParticle(&KF_0);
      con3.LinkVertex(&B_Vertex);
      
      ExKFitter Core;
      Core.LinkConstrain(&con1);
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
      
      //------------------  Sigma+ mass fit --------------------------------------
      //std::cout<<"Sigma mass fitter"<<std::endl;
      //Mdst bank
      Mdst_gamma Mdst_110 =(*i).child(1).child(1).child(0).mdstGamma(); //gamma from pi0
      Mdst_gamma Mdst_111 =(*i).child(1).child(1).child(1).mdstGamma(); //gamma from pi0
      Mdst_charged Mdst_10p = (*i).child(1).child(0).mdstCharged(); //p+ from Sigma+
      
      ExKFitterParticle ExPi0_g1(Mdst_110);
      ExKFitterParticle ExPi0_g2(Mdst_111);

      ExKFitterParticle ExP_10(Mdst_10p, pid_10); //proton
      
      Sigma_init.setX(IP.x() + (*i).child(1).p().px()/(*i).child(1).p().rho());
      Sigma_init.setY(IP.y() + (*i).child(1).p().py()/(*i).child(1).p().rho());
      Sigma_init.setZ(IP.z() + (*i).child(1).p().pz()/(*i).child(1).p().rho());
      ExKFitterVertex Sigma_MVertex(Sigma_init);
      
      ExKFitterMass Pi_Mass(0.134976); //pi0 mass
      ExKFitterMass Sigma_Mass(1.18937); //sigma mass
      
      ExKFitterParticle Sigma_fitter;
      Sigma_fitter.LinkParticle(&ExPi0_g1);
      Sigma_fitter.LinkParticle(&ExPi0_g2);
      Sigma_fitter.LinkParticle(&ExP_10);
      Sigma_fitter.LinkVertex(&Sigma_MVertex);

      ExKFitterConstrain conm1;
      conm1.SetMassConstrain();
      conm1.LinkParticle(&ExPi0_g1);
      conm1.LinkParticle(&ExPi0_g2);
      conm1.LinkVertex(&Sigma_MVertex);
      conm1.LinkMass(&Pi_Mass);
      
      ExKFitterConstrain conm2;
      conm2.SetMassConstrain();
      conm2.LinkParticle(&ExPi0_g1);
      conm2.LinkParticle(&ExPi0_g2);
      conm2.LinkParticle(&ExP_10);
      conm2.LinkVertex(&Sigma_MVertex);
      conm2.LinkMass(&Sigma_Mass);

      ExKFitterConstrain conm3;
      conm3.SetVertexConstrain();
      conm3.LinkParticle(&Sigma_fitter);
      conm3.LinkVertex(&B_Vertex);
      
      ExKFitter CoreM;
      CoreM.LinkConstrain(&conm1);
      CoreM.LinkConstrain(&conm2);
      CoreM.LinkConstrain(&conm3);
      int retm = CoreM.Minimize();
      float chisqExKm = CoreM.Chisq();
      float dof_exkm = CoreM.N_DegreeOfFreedom();
      
      b_tpl->column("chisqem",chisqExKm/dof_exkm);

      //------------------------- MC truth matching ------------------------
      
      int hindex = 0;
      //std::cout<<"MC truth matching..."<<std::endl;
      if (MCstatus == 1){
	const Mdst_charged ch_10 =(*i).child(1).child(0).mdstCharged(); //p+ from Sigma+
	const Mdst_gamma ch_110  =(*i).child(1).child(1).child(0).mdstGamma(); //gamma from pi0
	const Mdst_gamma ch_111  =(*i).child(1).child(1).child(1).mdstGamma(); //gamma from pi0
	const Mdst_gamma ch_2   =(*i).child(2).mdstGamma(); //gamma
	const Mdst_charged ch_0   =(*i).child(0).mdstCharged(); //proton
	

        Gen_hepevt Evtch_10=get_hepevt(ch_10);

        Gen_hepevt Evtch_110=get_hepevt(ch_110);
        Gen_hepevt Evtch_111=get_hepevt(ch_111);
	Gen_hepevt Evtch_2=get_hepevt(ch_2);
	Gen_hepevt Evtch_0=get_hepevt(ch_0);

        Gen_hepevt EvtP_10;
        Gen_hepevt EvtP_110;
	Gen_hepevt EvtP_111;
	Gen_hepevt EvtP_2;
        Gen_hepevt EvtP_0;
	
	Gen_hepevt EvtGP_10;
	Gen_hepevt EvtGP_110;
	Gen_hepevt EvtGP_111;
	Gen_hepevt EvtGP_0;
	Gen_hepevt EvtGP_2;

	Gen_hepevt EvtGGP_2;
	Gen_hepevt EvtGGP_10;
	Gen_hepevt EvtGGP_110;
	Gen_hepevt EvtGGP_111;
	
	Gen_hepevt EvtGGGP_110;
	Gen_hepevt EvtGGGP_111;
	
	Gen_hepevt EvtGGGGP_110;
	Gen_hepevt EvtGGGGP_111;
	
	int hi110, hi111;
	int is110, is111;
	int hi2;

	int mo2, gmo2, ggmo2;
	int mo110,mo111,gmo110,gmo111, ggmo110, ggmo111;
	int gggmo110, gggmo111, ggggmo110, ggggmo111;
	int ismo110, ismo111, isgmo110, isgmo111, isggmo110, isggmo111;

	hi110 = Evtch_110.idhep();
	hi111 = Evtch_111.idhep();
	is110 = Evtch_110.isthep();
	is111 = Evtch_111.isthep();
	hi2 = Evtch_2.idhep();

	b_tpl->column("hi10", Evtch_10.idhep());
	b_tpl->column("hi110", hi110);
	b_tpl->column("hi111", hi111);
	b_tpl->column("hi2", Evtch_2.idhep());
	b_tpl->column("hi0", Evtch_0.idhep());
	
        if (Evtch_10.mo(0))
        {
	  //parent of p-, should be Sigma+
	  EvtP_10=gen_mgr[Evtch_10.mo(0)-1]; //-1 because of index of fotran
	  b_tpl->column("mo10", EvtP_10.idhep());
	  if (EvtP_10.mo(0))
	    {
	      //grandparent of p-, should be the psudoparticle
	      EvtGP_10=gen_mgr[EvtP_10.mo(0)-1];
	      b_tpl->column("gmo10", EvtGP_10.idhep());
	    
	      if (EvtGP_10.mo(0))
		{
		  //three level up of p-, shold be B0
		  EvtGGP_10=gen_mgr[EvtGP_10.mo(0)-1];
		  b_tpl->column("ggmo10", EvtGGP_10.idhep());
		}
	    }
        }
	
        if (Evtch_110.mo(0))
        {
	  //parent of gamma, should be pi0
	  EvtP_110=gen_mgr[Evtch_110.mo(0)-1]; //-1 because of index of fotran
	  mo110 = EvtP_110.idhep();
	  ismo110 = EvtP_110.isthep();
	  b_tpl->column("mo110", mo110);
	  if (EvtP_110.mo(0))
	    {
	      //grandparent of gamma, should be Sigma+
	      EvtGP_110=gen_mgr[EvtP_110.mo(0)-1];
	      gmo110 = EvtGP_110.idhep();
	      isgmo110 = EvtGP_110.isthep();
	      b_tpl->column("gmo110", gmo110);
	      if (EvtGP_110.mo(0))
		{
		  //three level up of gamma, shold be the psudoparticle
		  EvtGGP_110=gen_mgr[EvtGP_110.mo(0)-1];
		  ggmo110 = EvtGGP_110.idhep();
		  isggmo110 = EvtGGP_110.isthep();
		  b_tpl->column("ggmo110", ggmo110);
		  if (EvtGGP_110.mo(0))
		    {
		      //four level up of gamma, shold be B0
		      EvtGGGP_110=gen_mgr[EvtGGP_110.mo(0)-1];
		      gggmo110 = EvtGGGP_110.idhep();
		      b_tpl->column("gggmo110", EvtGGGP_110.idhep());
		      if (EvtGGGP_110.mo(0))
			{
			  //fifth level up of gamma, one more parent generation to check 
			  EvtGGGGP_110=gen_mgr[EvtGGGP_110.mo(0)-1];
			  ggggmo110 = EvtGGGGP_110.idhep();
			}
		    }
		}
	    }
        }

        if (Evtch_111.mo(0))
        {
	  //parent of gamma, should be pi0
	  EvtP_111=gen_mgr[Evtch_111.mo(0)-1]; //-1 because of index of fotran
	  mo111 = EvtP_111.idhep();
	  ismo111 = EvtP_111.isthep();
	  b_tpl->column("mo111", mo111);
	  if (EvtP_111.mo(0))
	    {
	      //grandparent of gamma, should be Sigma+
	      EvtGP_111=gen_mgr[EvtP_111.mo(0)-1];
	      gmo111 = EvtGP_111.idhep();
	      isgmo111 = EvtGP_111.isthep();
	      b_tpl->column("gmo111", gmo111);
	    
	      if (EvtGP_111.mo(0))
		{
		  //three level up of gamma, shold be the psudoparticle
		  EvtGGP_111=gen_mgr[EvtGP_111.mo(0)-1];
		  ggmo111 = EvtGGP_111.idhep();
		  isggmo111 = EvtGGP_111.isthep();
		  b_tpl->column("ggmo111", ggmo111);
		  if (EvtGGP_111.mo(0))
		    {
		      //four level up of gamma, shold be B0
		      EvtGGGP_111=gen_mgr[EvtGGP_111.mo(0)-1];
		      gggmo111 = EvtGGGP_111.idhep();
		      b_tpl->column("gggmo111", EvtGGGP_111.idhep());
		      if (EvtGGGP_111.mo(0))
			{
			  //fifth level up of gamma, one more parent generation to check 
			  EvtGGGGP_111=gen_mgr[EvtGGGP_111.mo(0)-1];
			  ggggmo111 = EvtGGGGP_111.idhep();
			}
		    }
		}
	    }
        }
	
	
	if (Evtch_2.mo(0))
	  {
	    //parent of gamma 1, should be B0
            EvtP_2=gen_mgr[Evtch_2.mo(0)-1];
            b_tpl->column("mo2", EvtP_2.idhep());
	    mo2 = EvtP_2.idhep();
	    if (EvtP_2.mo(0))
	      {
		//two level up than gamma1, shold be B0
		EvtGP_2=gen_mgr[EvtP_2.mo(0)-1];
		b_tpl->column("gmo2", EvtGP_2.idhep());
		gmo2 = EvtGP_2.idhep();
		if (EvtGP_2.mo(0))
		  {
		    //three level up than gamma1, shold be B0
		    EvtGGP_2=gen_mgr[EvtGP_2.mo(0)-1];
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
		//grandparent of proton 1, should be B0
                EvtGP_0=gen_mgr[EvtP_0.mo(0)-1];
                b_tpl->column("gmo0", EvtGP_0.idhep());
	      }
	  }
	
	//B0 decay, Simga+ decay to p+ (10) and pi0->gamma, gamma (110,111)
        if(Evtch_10.idhep()== 2212 && EvtP_10.idhep()== 3222 && EvtGP_10.idhep()== 9487 && EvtGGP_10.idhep()== -511 
	   && Evtch_0.idhep()== -2212 && EvtP_0.idhep()== 9487 && EvtGP_0.idhep()== -511 ){
	  if ((mo2==-511||gmo2==-511||ggmo2==-511)&&
	      (hi2!=9487||mo2!=9487||gmo2!=9487||ggmo2!=9487||abs(mo2)!=2212||abs(gmo2)!=2212||abs(ggmo2)!=2212) )
	    { hindex = 1;}	  
	}

	//From here, to look up for two gamma rays from pi0 decay
	if ((ggmo110==-511||gggmo110==-511||ggggmo110==-511)&&(hi110==111 ||mo110==111 || gmo110==111|| ggmo110==111)){
	  hindex*=1;
	}
	else{
	  hindex*=0;}

	if ((ggmo111==-511||gggmo111==-511||ggggmo111==-511)&&(hi110==111 ||mo110==111 || gmo110==111|| ggmo110==111)){
	  hindex*=1;
	}
	else{
	  hindex*=0;}
		
	b_tpl->column("hindex", hindex);
	
	//--------------------- trace from top to down (from B meson)------------------------------------
        int nbpd=0,nbnd=0;
	//std::cout<<"Trace from top to down...."<<std::endl;
	for (std::vector<Gen_hepevt>::iterator k = gen_mgr.begin();k != gen_mgr.end(); k++)
          {
            int mgr_id;
            const Gen_hepevt &h_Brec = (*k);
            mgr_id = h_Brec ? h_Brec.idhep(): 0;

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

//---------------------------------------------

      //std::cout<<"Dump data...."<<std::endl;
      b_tpl->dumpData(); 
      *status = 1;
#endif
    } //end of B candidates

}//end of events

//void ana_psigplusg_1211::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P)
void ana_psigplusg_1211::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex)
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


  if (otherBfinal.size()<50)
    {
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
}


                      


//for int t--> 0:e 1:mu 2:pi 3:k 4:p 
void ana_psigplusg_1211::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t) 
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
