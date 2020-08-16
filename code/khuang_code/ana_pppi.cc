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
class ana_pppi : public Module 
{
public:
  ana_pppi(void){};
  ~ana_pppi(void){};
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
  double helic4(HepLorentzVector& childp, HepLorentzVector& motherp);


public:
private:
  BelleTuple *b_tpl;
  int evtcount;
  brutus_f Fisher_fw;
  brutus_f Fisher_ksfw[7];

};


extern "C" Module_descr *mdcl_ana_pppi()
{
  ana_pppi *module = new ana_pppi;
  Module_descr *dscr = new Module_descr ( "ana_pppi", module );
  IpProfile::define_global(dscr);
  BeamEnergy::define_global(dscr);
  return dscr;
}



void 
ana_pppi::hist_def(void)
{ 
  extern BelleTupleManager *BASF_Histogram;  
  BelleTupleManager& tm = *BASF_Histogram;
#ifdef BTUPLE
//ntuple name length--> 8 characters(upper limit)
  b_tpl = BASF_Histogram->ntuple("B -->pp_minus piPM ","de mbc massb pxb pyb pzb eb\
  massrho pxrho pyrho pzrho erho runid evtid putid expid eventid farmid \
  massdp massdm chisqexk\
  mass0 mass0_re px0 py0 pz0 e0 pid0 chrg0 dr0 dz0 cos0 cos0_bb\
  mass1 mass1_re px1 py1 pz1 e1 pid1 chrg1 dr1 dz1 cos1 cos1_bb\
  mass2 mass2_re px2 py2 pz2 e2 pid2 chrg2 dr2 dz2 cos2 cos2_bb\
  mpp\
  cos02_0 cos12_0 cos01_0 cos01_1 cos02_1 cos12_1 cos02_01 cos12_01 cos01_01 cos02_bb cos12_bb cos01_bb\
  hindex hi0 mo0 gmo0 hi1 mo1 gmo2 hi2 mo2 \
  bpd1 bpd2 bpd3 bpd4 bpd5 bpd6 bpd7 bpd8 bpd9 nbpd \
  bnd1 bnd2 bnd3 bnd4 bnd5 bnd6 bnd7 bnd8 bnd9 nbnd charged \
  bvertexx bvertexy bvertexz overtexx overtexy overtexz \
  pmiss emiss R2 costhr costhp spher cosb vtflag deltaz R2s \
  R1so R2so R3so R4so R1gso R2gso R3gso R4gso R1oo R2oo R3oo R4oo \
  imm2 mmm2 et H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son \
  H4son H0som H1som H2som H3som H4som H0soc H1soc H2soc H3soc H4soc \
  qr flavor imode");

#endif
  Fisher_fw.histogram(10,"Super Fox-Wolfram");
  k_sfw::initialize(Fisher_ksfw);

} 


void ana_pppi::init(int*)
{
  evtcount=0;

  Hamlet::init();

}


void ana_pppi::term (void)
{

}

// Set IP location
HepPoint3D IP(0,0,0);
HepSymMatrix IPerr(3,0);


void ana_pppi::begin_run(BelleEvent *evptr, int *status)
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
ana_pppi::event(BelleEvent *evptr, int *status)
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
       Mdst_ecl_Manager &ecl_mag            = Mdst_ecl_Manager::get_manager();
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

  std::vector<Particle> Lam_list;
  std::vector<Particle> Lambar_list;

  std::vector<Particle> Sigma_list;
  std::vector<Particle> Sigmabar_list;

  std::vector<Particle> Sigma_plus;
  std::vector<Particle> mu_plus;
  std::vector<Particle> mu_minus;
  std::vector<Particle> e_plus;
  std::vector<Particle> e_minus;
  std::vector<Particle> B_cand1;
  std::vector<Particle> D0;
  std::vector<Particle> D0B;
  std::vector<Particle> rho0;
  std::vector<Particle> rho_plus;
  std::vector<Particle> rho_minus;

  //Particle Type(/sw/belle/belle/b20040727_1143/share/data-files/qq98/decay.dec)
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


  Ptype ptype_Sigma_plus("SIG+");
  Ptype ptype_mu_plus("MU+");
  Ptype ptype_mu_minus("MU-");
  Ptype ptype_e_plus("E+");
  Ptype ptype_e_minus("E-");
  Ptype ptype_Bplus("B+");
  Ptype ptype_Bminus("B-");
  Ptype ptype_rho0("RHO0");
  Ptype ptype_rho_plus("RHO+");
  Ptype ptype_rho_minus("RHO-");
  Ptype ptype_B0("B0");
  Ptype ptype_B0b("B0B");
//  atc_pid selKpi(0,1,0,3,2);//for the reprocessed exp7 data
    atc_pid selkpi(3,1,5,3,2); 
    atc_pid selkp(3,1,5,3,4);
  atc_pid selppi(3,1,5,4,2); 
  atc_pid selpk(3,1,5,4,3);
//  atc_pid selmuk(3,1,5,1,3);
//  atc_pid selmupi(3,1,5,1,2); 

  // distinguish MC from data,MC is MCstatus == 1
  int MCstatus=0; 
  for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin(); i != gen_mgr.end(); i++){
	MCstatus=1;
  }



  //fill pi or k lists from MDST_Charged Data Base
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin(); i != charged_mag.end(); i++){
  	
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
          if(selppi.prob(*i)>0.6 && selpk.prob(*i)>0.6)
	//if(selkpi.prob(*i)>0.)
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
	  else if (selkpi.prob(*i)<0.4)
	 //if (selppi.prob(*i)<1.0)
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

	combination(B_cand1, ptype_Bplus, p_plus, p_minus, pi_plus);
	combination(B_cand1, ptype_Bminus, p_plus, p_minus, pi_minus);	
//for B candicate 21

  for(std::vector<Particle>::iterator i =B_cand1.begin(); i != B_cand1.end(); i++)
{ 

#ifdef BTUPLE
//------------------- dE and Mbc of B meson -------------------------
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

        b_tpl->column("evtid", Evt);
        b_tpl->column("runid", Run);
        b_tpl->column("expid", Exp);
        b_tpl->column("eventid", Eventid);
        b_tpl->column("farmid", Farmid);
        
//------------------- information of B (*i) -------------------------
        b_tpl->column("pxb", (*i).px());
        b_tpl->column("pyb", (*i).py());
        b_tpl->column("pzb", (*i).pz());
        b_tpl->column("eb", (*i).e());
        
//------------------- information of p from pp_minus (*i).child(0) -------------------------

		float px0 = (*i).child(0).px();
		float py0 = (*i).child(0).py();
		float pz0 = (*i).child(0).pz();
		float e0 = (*i).child(0).e();
		float mass0_re = sqrt(pow(e0,2)-(px0*px0+py0*py0+pz0*pz0));
	
        b_tpl->column("mass0", (*i).child(0).p().mag());
        b_tpl->column("mass0_re", mass0_re);        
        b_tpl->column("px0", px0);
        b_tpl->column("py0", py0);
        b_tpl->column("pz0", pz0);
        b_tpl->column("e0", e0);
        

        const Mdst_charged* chk0 = &(*i).child(0).mdstCharged();
        float ppi_id0 = selppi.prob(chk0);
        b_tpl->column("pid0", ppi_id0);
          
        float chrg0 = chk0->charge();
        b_tpl->column("chrg0", chrg0);

        double dr0,dz0;
        GetImpactParameters(chk0,&dr0,&dz0,4);
        b_tpl->column("dr0", float(dr0));
        b_tpl->column("dz0", float(dz0));

//------------------- information of cos0 -------------------------	
	
        double cos_polar_0;
        Vector3 ch0_momentum((*i).child(0).px(), (*i).child(0).py(), (*i).child(0).pz());
        Vector3 P_BEAM_polar ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));

        cos_polar_0=ch0_momentum.dot(P_BEAM_polar);
        cos_polar_0=cos_polar_0/(ch0_momentum.mag()*P_BEAM_polar.mag());
        
        b_tpl->column("cos0", cos_polar_0);

   
//------------------- information of p from PP- (*i).child(1) -------------------------
		float px1 = (*i).child(1).px();
		float py1 = (*i).child(1).py();
		float pz1 = (*i).child(1).pz();
		float e1 = (*i).child(1).e();
		float mass1_re = sqrt(pow(e1,2)-(px1*px1+py1*py1+pz1*pz1));
			
        b_tpl->column("mass1", (*i).child(1).p().mag());
        b_tpl->column("mass1_re", mass1_re);        
        b_tpl->column("px1", px1);
        b_tpl->column("py1", py1);
        b_tpl->column("pz1", pz1);
        b_tpl->column("e1", e1);
        
        const Mdst_charged* chk1 = &(*i).child(1).mdstCharged();
        float ppi_id1 = selppi.prob(chk1);
        b_tpl->column("pid1", ppi_id1);

        float chrg1 = chk1->charge();
        b_tpl->column("chrg1", chrg1);

        double dr1,dz1;
        GetImpactParameters(chk1,&dr1,&dz1,4);
        b_tpl->column("dr1", float(dr1));
        b_tpl->column("dz1", float(dz1));

        double cos_polar_1;
        Vector3 ch1_momentum((*i).child(1).px(), (*i).child(1).py(), (*i).child(1).pz());

        cos_polar_1=ch1_momentum.dot(P_BEAM_polar);
        cos_polar_1=cos_polar_1/(ch1_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos1", cos_polar_1);
        

//------------------- information of pp_minus (*i).child(0)/(1) -------------------------
		float px2 = (*i).child(2).px();
		float py2 = (*i).child(2).py();
		float pz2 = (*i).child(2).pz();
		float e2 = (*i).child(2).e();
		float mass2_re = sqrt(pow(e2,2)-(px2*px2+py2*py2+pz2*pz2));
		
        b_tpl->column("mass2", (*i).child(2).p().mag());
        b_tpl->column("mass2_re", mass2_re);        
        b_tpl->column("px2", px2);
        b_tpl->column("py2", py2);
        b_tpl->column("pz2", pz2);
        b_tpl->column("e2", e2);

        const Mdst_charged* chk2 = &(*i).child(2).mdstCharged();
        float ppi_id2 = selppi.prob(chk2);
        b_tpl->column("pid2", ppi_id2);
        
        float chrg2 = chk2->charge();
        b_tpl->column("chrg2", chrg2); 
        
        double dr2,dz2;
        GetImpactParameters(chk2,&dr2,&dz2,2);
        b_tpl->column("dr2", float(dr2));
        b_tpl->column("dz2", float(dz2));

        double cos_polar_2;
        Vector3 ch2_momentum((*i).child(2).px(), (*i).child(2).py(), (*i).child(2).pz());

        cos_polar_2=ch2_momentum.dot(P_BEAM_polar);
        cos_polar_2=cos_polar_2/(ch2_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos2", cos_polar_2);

//------------------- mpp -------------------------
        b_tpl->column("mpp", sqrt(pow(mass0_re,2)+pow(mass1_re,2)+2*(e0*e1-px0*px1-py0*py1-pz0*pz1)));
        
//begin to calculate the angle between ppbar and pi+-
//declare boost vector
        HepLorentzVector boostVectorCh0((*i).child(0).p());
        HepLorentzVector boostVectorCh1((*i).child(1).p());
        HepLorentzVector boostVectorCh01 = boostVectorCh0+boostVectorCh1;       
        HepLorentzVector boost_vector_beam(-E_HER*sin(cross_angle), 0.0,
                                    E_LER-E_HER*cos(cross_angle),
                                    E_HER+E_LER);
									        
//declare LorentzVector
        HepLorentzVector ltzVectorCh0((*i).child(0).p());
        HepLorentzVector ltzVectorCh1((*i).child(1).p());
        HepLorentzVector ltzVectorCh2((*i).child(2).p());
        HepLorentzVector beamDir(E_HER*sin(cross_angle), 0.0,
                                    -E_LER+E_HER*cos(cross_angle),
                                    E_HER+E_LER);
                                    
//boost to child 0
       	ltzVectorCh0.boost(-boostVectorCh0.boostVector());
       	ltzVectorCh1.boost(-boostVectorCh0.boostVector());
       	ltzVectorCh2.boost(-boostVectorCh0.boostVector());       
       	
        double cos02_0 = ((ltzVectorCh0.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
        double cos12_0 = ((ltzVectorCh1.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh1.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
        double cos01_0 = ((ltzVectorCh0.vect().dot(ltzVectorCh1.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh1.vect().mag2())));
         
  		        
        b_tpl->column("cos02_0", cos02_0);				   
        b_tpl->column("cos12_0", cos12_0);
       b_tpl->column("cos01_0", cos01_0);
       
//boost to child 1

        ltzVectorCh0 = ((*i).child(0).p());
        ltzVectorCh1 = ((*i).child(1).p());  
        ltzVectorCh2 = ((*i).child(2).p());   		      

       	ltzVectorCh0.boost(-boostVectorCh1.boostVector());
       	ltzVectorCh1.boost(-boostVectorCh1.boostVector());
       	ltzVectorCh2.boost(-boostVectorCh1.boostVector());       
		   
	     double cos02_1 = ((ltzVectorCh0.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
        double cos12_1 = ((ltzVectorCh1.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh1.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
        double cos01_1 = ((ltzVectorCh0.vect().dot(ltzVectorCh1.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh1.vect().mag2())));
                
        b_tpl->column("cos02_1", cos02_1);				   
        b_tpl->column("cos12_1", cos12_1);			
        b_tpl->column("cos01_1", cos01_1);
		        
//boost to ppbar
        ltzVectorCh0 = ((*i).child(0).p());
        ltzVectorCh1 = ((*i).child(1).p());  
        ltzVectorCh2 = ((*i).child(2).p());   	
       	ltzVectorCh0.boost(-boostVectorCh01.boostVector());
       	ltzVectorCh1.boost(-boostVectorCh01.boostVector());
       	ltzVectorCh2.boost(-boostVectorCh01.boostVector());       

        double cos02_01 = ((ltzVectorCh0.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
        double cos12_01 = ((ltzVectorCh1.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh1.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
        double cos01_01 = ((ltzVectorCh0.vect().dot(ltzVectorCh1.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh1.vect().mag2())));
                 		        
        b_tpl->column("cos02_01", cos02_01);				   
        b_tpl->column("cos12_01", cos12_01);        
        b_tpl->column("cos01_01", cos01_01);        

//boost to BBbar
        ltzVectorCh0 = ((*i).child(0).p());
        ltzVectorCh1 = ((*i).child(1).p());  
        ltzVectorCh2 = ((*i).child(2).p());   	
     	ltzVectorCh0.boost( boost_vector_beam.boostVector() );
     	ltzVectorCh1.boost( boost_vector_beam.boostVector() );
     	ltzVectorCh2.boost( boost_vector_beam.boostVector() );
		beamDir.boost( boost_vector_beam.boostVector() );
		

        double cos0_bb = ((ltzVectorCh0.vect().dot(beamDir.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(beamDir.vect().mag2())));
        double cos1_bb =  ((ltzVectorCh1.vect().dot(beamDir.vect()))/(sqrt(ltzVectorCh1.vect().mag2())*sqrt(beamDir.vect().mag2())));
        double cos2_bb =  ((ltzVectorCh2.vect().dot(beamDir.vect()))/(sqrt(ltzVectorCh2.vect().mag2())*sqrt(beamDir.vect().mag2())));    

		double cos02_bb = (ltzVectorCh0.vect().dot(ltzVectorCh2.vect())/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
     	double cos12_bb = ((ltzVectorCh1.vect().dot(ltzVectorCh2.vect()))/(sqrt(ltzVectorCh1.vect().mag2())*sqrt(ltzVectorCh2.vect().mag2())));
      	double cos01_bb = ((ltzVectorCh0.vect().dot(ltzVectorCh1.vect()))/(sqrt(ltzVectorCh0.vect().mag2())*sqrt(ltzVectorCh1.vect().mag2())));
		        	
        b_tpl->column("cos02_bb", cos02_bb);				   
        b_tpl->column("cos12_bb", cos12_bb);        
        b_tpl->column("cos01_bb", cos01_bb);       
        b_tpl->column("cos0_bb", cos0_bb);      
        b_tpl->column("cos1_bb", cos1_bb);      
        b_tpl->column("cos2_bb", cos2_bb);              
//end calculate the angle between ppbar and pi+-        
			


//------------------- B vertex fit with ExKfitter  -------------------------
       // ExKFitterVertex B_Vertex(IP,IPerr);
       // HepPoint3D bvertex = B_Vertex.Vertex();

        Mdst_charged Mdst_0 =(*i).child(0).mdstCharged();
        Mdst_charged Mdst_1 =(*i).child(1).mdstCharged();
        Mdst_charged Mdst_2 =(*i).child(2).mdstCharged();
        

        ExKFitterParticle KF_0(Mdst_0, 4);
        ExKFitterParticle KF_1(Mdst_1, 4);
		ExKFitterParticle KF_2(Mdst_2, 2);

		ExKFitterVertex B_Vertex(IP,IPerr);

		ExKFitterParticle B;
        B.LinkParticle(&KF_0);
		B.LinkParticle(&KF_1);
		B.LinkParticle(&KF_2);
		B.LinkVertex(&B_Vertex);
	
		ExKFitterConstrain con1;
		con1.SetVertexConstrain();
	    con1.LinkParticle(&KF_0);
		con1.LinkParticle(&KF_1);
		con1.LinkParticle(&KF_2);
		con1.LinkVertex(&B_Vertex);

        ExKFitter Core;
        Core.LinkConstrain(&con1);
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
if (MCstatus == 1)
{
    const Mdst_charged ch_0 = (*i).child(0).mdstCharged();
    const Mdst_charged ch_1 = (*i).child(1).mdstCharged();
    const Mdst_charged ch_2 = (*i).child(2).mdstCharged();
    
    Gen_hepevt Evtch_0=get_hepevt(ch_0);
    Gen_hepevt Evtch_1=get_hepevt(ch_1);
     Gen_hepevt Evtch_2=get_hepevt(ch_2);   

    Gen_hepevt EvtP_0;
    Gen_hepevt EvtP_1;
    Gen_hepevt EvtP_2;
    
	Gen_hepevt EvtGP_0;
	Gen_hepevt EvtGP_1;

	b_tpl->column("hi0", Evtch_0.idhep());
	b_tpl->column("hi1", Evtch_1.idhep());
	b_tpl->column("hi2", Evtch_2.idhep());	

	//find the ID of the parent 
        if (Evtch_0.mo(0))
        {
            EvtP_0=gen_mgr[Evtch_0.mo(0)-1];
	    b_tpl->column("mo0", EvtP_0.idhep());
		if (EvtP_0.mo(0))
		{
		EvtGP_0=gen_mgr[EvtP_0.mo(0)-1];
		 b_tpl->column("gmo0", EvtGP_0.idhep());
		}
        }

        if (Evtch_1.mo(0))
        {
            EvtP_1=gen_mgr[Evtch_1.mo(0)-1];
	    b_tpl->column("mo1", EvtP_1.idhep());
		if (EvtP_1.mo(0))
                {
                EvtGP_1=gen_mgr[EvtP_1.mo(0)-1];
		b_tpl->column("gmo1", EvtGP_1.idhep());
                }
        }

        if (Evtch_2.mo(0))
        {
            EvtP_2=gen_mgr[Evtch_2.mo(0)-1];
	    b_tpl->column("mo2", EvtP_2.idhep());
        }

        if(Evtch_0.idhep() == 2212 && Evtch_1.idhep() == -2212\ 
		&& EvtP_0.idhep() == 4431 && EvtP_1.idhep() == 4431\
    	&& EvtGP_0.idhep() == 521 && EvtGP_1.idhep() == 521\
		&& Evtch_2.idhep() == 211 && EvtP_2.idhep() == 521)
		{ hindex = 1;}
        else if(Evtch_0.idhep() == 2212 && Evtch_1.idhep() == -2212\ 
		&& EvtP_0.idhep() == 4431 && EvtP_1.idhep() == 4431\
    	&& EvtGP_0.idhep() == -521 && EvtGP_1.idhep() == -521\
		&& Evtch_2.idhep() == -211 && EvtP_2.idhep() == -521)
        { hindex = 1;}
        else
        { hindex = 0;}
        
	b_tpl->column("hindex", hindex);
	//--------------------- trace from top to down (from B meson)---------move this part outside the B candidate loop will lead to error
        int nbpd=0,nbnd=0;
//        int cout1=0;
//        int cout2=0;  
        for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();i != gen_mgr.end(); i++)
        {               
                  if ((*i).idhep()==521||(*i).idhep()== 511)
                  {       
//				  printf("cout1 is %d \n", cout1);
//                  cout1 ++;
                  
                    int da1=(*i).da(0),da2=(*i).da(1);
                           
                    if ((*i).idhep() == 521) b_tpl->column("charged",1);
                    else b_tpl->column("charged",0);
                    nbpd=(da2-da1+1);
                    b_tpl->column("nbpd",nbpd);
                      for(int start=0;start<(da2-da1+1);start++)
                      {
                        Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
                        char bpdnumber[32];
                        sprintf (bpdnumber,"%s%d","bpd",start+1);
                        b_tpl->column(bpdnumber,Evda.idhep());
                      }
                   }
                   else if((*i).idhep()== -521 || (*i).idhep()== -511)
                   {
//                printf("cout2 is %d \n", cout2);
 //                 cout2 ++;
                  
                        int da1=(*i).da(0),da2=(*i).da(1);

                          if ((*i).idhep() == -521) b_tpl->column("charged",1);
                          else b_tpl->column("charged",0);
                        nbnd=(da2-da1+1);
                        b_tpl->column("nbnd",nbnd);
                         for(int start=0;start<(da2-da1+1);start++)
                         {
                           Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
                           char bndnumber[32];
                           sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....
                           b_tpl->column(bndnumber,Evda.idhep());
                         }
                   }
        }
}//------------------------- MC truth matching ------------------------

//  	continuum suppresion
        int vertexflag;
        HepPoint3D overtex;

        Vector4 OtherBVector;
        float R2, spher, cos_thr, cos_thp, par_sfw[13];
        
        shape((*i), R2,  spher, cos_thr, cos_thp, par_sfw, OtherBVector,bvertex, vertexflag ,overtex);


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

        Hamlet hamlet;
        hamlet.setBcp(*i);
        hamlet.setTagMethod(Hamlet::MULT_DIM_LH);
        const double qr_evtgen = hamlet.q();
       // const Fbtag_MultiDimLikelihood0 & mdlh_evtgen = hamlet,fbtag_mult_dim_likelihood();
        b_tpl->column("qr",qr_evtgen);
        
        int flavor = hamlet.flavor();
        // +1 : for tag-side being B0 or B+
        // -1 : for tag-side being B0-bar or B-
        //  0 : does not specify a b-flavor
        // -2 : does not specify a b-flavor (no leptons and kaons are found)
        b_tpl->column("flavor",flavor);
       int imode = hamlet.imode();// 1 : tag a b-flavor with high-electron method
                                  // 2 : tag a b-flavor with high-muon method
                                  // 4 : tag a b-flavor with kaon method
                                  // others : does not specify a b-flavor
        b_tpl->column("imode",imode);
                                                                                                      
	b_tpl->dumpData(); 
        *status = 1;
#endif
} //end of PP-


//B_cand1 for k pi


}//end of fill pi or k lists from MDST_Charged Data Base

//void ana_pppi::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P)
void ana_pppi::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex)
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
    std::cout << candiBfinal.size()<<"+"<<otherBfinal.size()<<" !="<< allFinal.size()<<std::endl;


  //Start to fill shape variables
  Vector3 candiBthr=thrust(candiBfinal);
if (otherBfinal.size()<50)
{

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

}
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
if (otherBfinal.size()<50)
{

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
void ana_pppi::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t) 
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
double ana_pppi::helic4(HepLorentzVector& childp, HepLorentzVector& motherp){
  double hel;
  HepLorentzVector nchildp = childp;
  HepLorentzVector nmotherp = motherp;
  nchildp.boost(-(nmotherp.boostVector()));
  hel = nchildp.px()*nmotherp.px()+nchildp.py()*nmotherp.py()+nchildp.pz()*nmotherp.pz();
  hel = hel/sqrt(nchildp.px()*nchildp.px()+nchildp.py()*nchildp.py()+nchildp.pz()*nchildp.pz());
  hel = hel/sqrt(nmotherp.px()*nmotherp.px()+nmotherp.py()*nmotherp.py()+nmotherp.pz()*nmotherp.pz());
  return hel;
}        
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif