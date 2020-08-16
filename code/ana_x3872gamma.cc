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


#include "koppenburg/pi0eta_prob.h"
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

int checkMultiUse(Particle& i, Particle& j);


// Module class
class ana_x3872gamma : public Module
{
public:
    ana_x3872gamma(void) {};
    ~ana_x3872gamma(void) {};
    void init(int*);
    void term(void);
    void disp_stat(const char*) {};
    void hist_def(void);
    void event(BelleEvent*, int*);
    void begin_run(BelleEvent*, int*);
    void end_run(BelleEvent*, int*) {};
    void other(int*, BelleEvent*, int*) {};
//  void shape(Particle &, float &, float &, float &, float &, float par_sfw[17], Vector4 &);
    void shape(Particle &, float &, float &, float &, float &, float par_sfw[17],Vector4 &,HepPoint3D &,int &, HepPoint3D &);

    void GetImpactParameters(const Mdst_charged*, double* ,double* ,int);

    double deltaZ(Particle, double&, double&, double&, double&);


public:
private:
    BelleTuple *b_tpl;
    int evtcount;
    brutus_f Fisher_ksfw[7];
};


extern "C" Module_descr *mdcl_ana_x3872gamma()
{
    ana_x3872gamma *module = new ana_x3872gamma;
    Module_descr *dscr = new Module_descr ( "ana_x3872gamma", module );
    IpProfile::define_global(dscr);
    BeamEnergy::define_global(dscr);
    return dscr;
}



void
ana_x3872gamma::hist_def(void)
{
    extern BelleTupleManager *BASF_Histogram;
    BelleTupleManager& tm = *BASF_Histogram;
#ifdef BTUPLE

//ntuple name length--> 8 characters(upper limit)
    b_tpl = BASF_Histogram->ntuple("B0 -> x3872ga","de mbc0 mbc massb pxb pyb pzb eb ebeam\
      massx3872 pxx3872 pyx3872 pzx3872 ex3872 sumecl e9oe25\
      mass00 px00 py00 pz00 e00 pxi01 pyi01 pzi01 pxi02 pyi02 pzi02\
      mass01 px01 py01 pz01 e01 pid01 chrg01 dr01 dz01 cos01 cos01cm\
      mass02 px02 py02 pz02 e02 pid02 chrg02 dr02 dz02 cos02 cos02cm\
      mass000 px000 py000 pz000 e000 pid000 chrg000 dr000 dz000 cos000 cos000cm\
      mass001 px001 py001 pz001 e001 pid001 chrg001 dr001 dz001 cos001 cos001cm\
      pt00 pt01 pt02 pt000 pt001 ptb ptx3872 pt1 coshij coshir coshipxx coshipxyx coshipzx\
      mass1 px1 py1 pz1 e1 pid1 chrg1 dr1 dz1 cos1 cos1cm cos0 cos0cm cos00 cos00cm evtcount\
      vchisq vchisq1 chisqexk hindex hi00 mo00 gmo00 hi01 mo01 gmo01 hi02 mo02 gmo02 hi1 mo1 gmo1 \
      ggmo1 gggmo1 rho0mass rho0e rho0cos rhocoscm ptrho0 pxrho0 pyrho0 pzrho0\
      hi000 mo000 gmo000 ggmo000 hi001 mo001 gmo001 ggmo001\
      nmo000 nmo001 ngmo000 ngmo001 nggmo000 nggmo001 nmo01 nmo02 ngmo01 ngmo02 nmo1 ngmo1 nggmo1 ngggmo1\
      imm2 mmm2 et H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son\
      H0soc H1soc H2soc H3soc H4soc H0som H1som H2som H3som H4som\
      bvertexx bvertexy bvertexz overtexx overtexy overtexz pmiss emiss\
      R2 costhr spher cosb vtflag deltaz ggmo01 nggmo01 ggmo02 nggmo02 coshrho0 drrho0 dzrho0 chrgrho0 pidrho0\
      R2s R1so R2so R3so R4so R1gso R2gso R3gso R4gso R1oo R2oo R3oo R4oo\
      bpd1 bpd2 bpd3 bpd4 bpd5 nbpd bpdid1 bpdid2 bpdid3 bpdid4 bpdid5\
      npx3872 npjpsi npxd npjd nnx3872 nnjpsi nnxd nnjd deltap02 deltap01 jindex\
      bnd1 bnd2 bnd3 bnd4 bnd5 nbnd charged bndid1 bndid2 bndid3 bndid4 bndid5\
      bpxdid1 bpxdid2 bpxdid3 bpxdid4 bpxdid5 bpxd1 bpxd2 bpxd3 bpxd4 bpxd5\
      bpjdid1 bpjdid2 bpjdid3 bpjdid4 bpjdid5 bpjd1 bpjd2 bpjd3 bpjd4 bpjd5\
      bnxdid1 bnxdid2 bnxdid3 bnxdid4 bnxdid5 bnxd1 bnxd2 bnxd3 bnxd4 bnxd5\
      bnjdid1 bnjdid2 bnjdid3 bnjdid4 bnjdid5 bnjd1 bnjd2 bnjd3 bnjd4 bnjd5\
      cosll cosjl cos2pi cosrpi cosrj cosxj cosxg cospgg cosegg\
      cosllcm cosjlcm cos2picm cosrpicm cosrjcm cosxjcm cosxgcm cospggcm coseggcm\
      probpi0 probeta0 mpi0 meta0 pie2 etae2 np qr flavor imode\
      cosh0 cosh00 cosh000 cosh001 cosh01 cosh02 cosh1 expnum runnum\
      coshipxj coshipyj coshipzj coshiej coshipxr coshipyr coshipzr coshier");//ebeam added,20171210
/*
    b_tpl = BASF_Histogram->ntuple("B0 -> x3872ga","de mbc0 mbc massb pxb pyb pzb eb ebeam\
      massx3872 pxx3872 pyx3872 pzx3872 ex3872 sumecl e9oe25\
      mass00 px00 py00 pz00 e00 pxi01 pyi01 pzi01 pxi02 pyi02 pzi02\
      mass01 px01 py01 pz01 e01 pid01 chrg01 dr01 dz01 cos01 \
      mass02 px02 py02 pz02 e02 pid02 chrg02 dr02 dz02 cos02 \
      mass000 px000 py000 pz000 e000 pid000 chrg000 dr000 dz000 cos000 \
      mass001 px001 py001 pz001 e001 pid001 chrg001 dr001 dz001 cos001 \
      pt00 pt01 pt02 pt000 pt001 ptb ptx3872 pt1\
      mass1 px1 py1 pz1 e1 pid1 chrg1 dr1 dz1 cos1 cos1cm cos0 cos0cm cos00 cos00cm evtcount\
      vchisq vchisq1 chisqex0 chisqe00 chisqexk hindex hi00 mo00 gmo00 hi01 mo01 gmo01 hi02 mo02 gmo02 hi1 mo1 gmo1 \
      ggmo1 gggmo1 ggggmo1 gggggmo1 rho0mass rho0e rho0cos rhocoscm\
      hi000 mo000 gmo000 ggmo000 hi001 mo001 gmo001 ggmo001\
      nmo000 nmo001 ngmo000 ngmo001 nggmo000 nggmo001 nmo01 nmo02 ngmo01 ngmo02 nmo1 ngmo1 nggmo1 ngggmo1 nggggmo1 ngggggm1\
      imm2 mmm2 et H0oo H1oo H2oo H3oo H4oo H0son H1son H2son H3son H4son\
      H0soc H1soc H2soc H3soc H4soc H0som H1som H2som H3som H4som\
      bvertexx bvertexy bvertexz overtexx overtexy overtexz pmiss emiss\
      R2 costhr spher cosb vtflag deltaz\
      R2s R1so R2so R3so R4so R1gso R2gso R3gso R4gso R1oo R2oo R3oo R4oo\
      bpd1 bpd2 bpd3 bpd4 bpd5 bpd6 bpd7 bpd8 bpd9 nbpd \
      bpdid1 bpdid2 bpdid3 bpdid4 bpdid5 bpdid6 bpdid7 bpdid8 bpdid9\
      nx3872 njpsi ngamma nxd njd deltap02 deltap01\
      bndid1 bndid2 bndid3 bndid4 bndid5 bndid6 bndid7 bndid8 bndid9\
      bxdid1 bxdid2 bxdid3 bxdid4 bxdid5 bxdid6 bxdid7 bxdid8 bxdid9\
      bjdid1 bjdid2 bjdid3 bjdid4 bjdid5 bjdid6 bjdid7 bjdid8 bjdid9\
      bxd1 bxd2 bxd3 bxd4 bxd5 bxd6 bxd7 bxd8 bxd9\
      bjd1 bjd2 bjd3 bjd4 bjd5 bjd6 bjd7 bjd8 bjd9\
      bnd1 bnd2 bnd3 bnd4 bnd5 bnd6 bnd7 bnd8 bnd9 nbnd charged");*/
 k_sfw::initialize(Fisher_ksfw);

#endif
}


void ana_x3872gamma::init(int*)
{
    evtcount=0;

    Hamlet::init();

}


void ana_x3872gamma::term (void)
{
}

// Set IP location
HepPoint3D IP(0,0,0);
HepSymMatrix IPerr(3,0);


void ana_x3872gamma::begin_run(BelleEvent *evptr, int *status)
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
ana_x3872gamma::event(BelleEvent *evptr, int *status)
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
    Mdst_ecl_Manager &ecl_mag            = Mdst_ecl_Manager::get_manager();
    Evtcls_hadronic_flag_Manager &evtcls = Evtcls_hadronic_flag_Manager::get_manager();
    Gen_hepevt_Manager &gen_mgr          = Gen_hepevt_Manager::get_manager();

    //Mdst_klong_Manager &klong_mag        = Mdst_klong_Manager::get_manager();

    //for reprocessed data
    //  remove_duplicates();
    //scale_momenta(1.00246, 1.0);
    //for reprocessed exp7 data
    //scale_momenta(1.00221, 1.0);
    //for reprocessed exp9 data
    //scale_momenta(1.00149, 1.0);


    evtcount++;
    //cout << "event = " << evtcount << endl;
//debugging
    double E_HER=BeamEnergy::E_HER();
    double E_LER=BeamEnergy::E_LER();
    double cross_angle=BeamEnergy::Cross_angle();//radian
    static Vector4 P_BEAM ( -E_HER*sin(cross_angle), 0.0, E_LER-E_HER*cos(cross_angle), E_HER+E_LER);

    //int var1[5]={0},var2[5]={0},var3[5]={0};//STAR

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
    std::vector<Particle> gamma1;
    //std::vector<Particle> pi_0;
    std::vector<Particle> pi_plus;
    std::vector<Particle> pi_minus;
    
    //std::vector<Particle> k_plus;
    //std::vector<Particle> k_minus;
    //std::vector<Particle> k_short;
    //std::vector<Particle> p_plus;
    //std::vector<Particle> p_minus;
    //std::vector<Particle> Lamda;
    //std::vector<Particle> Lamdabar;
    //std::vector<Particle> Sigma_plus;
    std::vector<Particle> mu_plus;
    std::vector<Particle> mu_minus;
    std::vector<Particle> e_plus;
    std::vector<Particle> e_minus;
    std::vector<Particle> tte_plus;
    std::vector<Particle> tte_minus;
    std::vector<Particle> B_cand1;
    //std::vector<Particle> D0;
    //std::vector<Particle> D0B;
    std::vector<Particle> B0;
    std::vector<Particle> B0B;
    std::vector<Particle> jpsi;
    std::vector<Particle> X3872;
    std::vector<Particle> rho0;

    //Particle Type(/belle/belle/b20040727_1143/share/data-files/qq98/decay.dec)
    Ptype ptype_gamma("GAMM");
    Ptype ptype_gamma1("GAMM");
    //Ptype ptype_D0("D0");
    //Ptype ptype_D0B("D0B");
    //Ptype ptype_pi_0("PI0");
    Ptype ptype_pi_plus("PI+");
    Ptype ptype_pi_minus("PI-");
    //Ptype ptype_pi1("PI+");
    //Ptype ptype_pi2("PI-");
    //Ptype ptype_k_plus("K+");
    //Ptype ptype_k_minus("K-");
    //Ptype ptype_k_short("K0S");
    //Ptype ptype_p_plus("P+");
    //Ptype ptype_p_minus("AP+");
    //Ptype ptype_Lamda("LAM");
    //Ptype ptype_Lamdabar("ALAM");
    //Ptype ptype_Sigma_plus("SIG+");
    Ptype ptype_mu_plus("MU+");
    Ptype ptype_mu_minus("MU-");
    Ptype ptype_e_plus("E+");
    Ptype ptype_e_minus("E-");
    Ptype ptype_tte_plus("E+");
    Ptype ptype_tte_minus("E-");
    //Ptype ptype_Bplus("B+");
    //Ptype ptype_Bminus("B-");

    Ptype ptype_jpsi("PSI");
    Ptype ptype_X3872("PSI");//QUESTIONCHI1
    Ptype ptype_B0("B0");
    Ptype ptype_B0B("B0B");
    Ptype ptype_rho0("RHO0");
//QUESTION



//  atc_pid selKpi(0,1,0,3,2);//for the reprocessed exp7 data
    atc_pid selkpi(3,1,5,3,2);
    //atc_pid selkp(3,1,5,3,4);
//  atc_pid selppi(3,1,5,4,2);
//  atc_pid selpk(3,1,5,4,3);
    //atc_pid selmuk(3,1,5,1,3);
    //atc_pid selmupi(3,1,5,1,2);


    // distinguish MC from data,MC is MCstatus == 1
    int MCstatus=0;
    for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin();
            i != gen_mgr.end(); i++)
    {
        MCstatus=1;
    }
    
    // I add gamma here //STAR
    for(std::vector<Mdst_gamma>::iterator it = gamma_mag.begin();it !=gamma_mag.end(); it++)
    {
      int flag=0;
      /*
      int nfinal = Btag.relation().nFinalStateParticles();
      for(int k=0; k < nfinal; k++)
      {
        Particle child = Btag.relation().finalStateParticle(k);
        if(child.mdstGamma())
        {
          if(child.mdstGamma().get_ID()==(*it).get_ID())
            flag=1;
        }
        if(flag==1) continue;

      }
      */
      //cout<<"line908";

      int match= (*it).ecl().match();
      if(match!=0) continue;
      
      Particle tmp(*it);
      gamma.push_back(tmp);

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

//QUESTION
        // muon
        if( !reject && mu_like>0.8 )
        {
            if ((*i).charge()>0)
            {
                Particle tmp(*i, ptype_mu_plus);
                mu_plus.push_back(tmp);//mu+ ++
            }
            else
            {
                Particle tmp(*i, ptype_mu_minus);
                mu_minus.push_back(tmp);
            }
        }
        // electron
        if( eid_prob>0.9 )//STAR
        {
            if ((*i).charge()>0)////like e+
            {
                if (ch.trk().quality() == 0)// Obtain "Good track"
                {
                    Particle tmp(*i, ptype_e_plus);
                    e_plus.push_back(tmp);
                }
            }
            else ////like e-
            {
                if (ch.trk().quality() == 0)// Obtain "Good track"
                {
                    Particle tmp(*i, ptype_e_minus);
                    e_minus.push_back(tmp);
                }
            }
        }
//===========================================
        if((reject || mu_like < 0.95) && eid_prob < 0.95 )
        {
            if(selkpi.prob(*i)>0.6)// like k
                //if(selkpi.prob(*i)>0.)
            {
                if ((*i).charge()>0)
                {
                    if (ch.trk().quality() == 0)  // Obtain "Good track"
                    {
                        //Particle tmp(*i, ptype_k_plus);
                        //k_plus.push_back(tmp);
                    }
                }
                else
                {
                    if (ch.trk().quality() == 0)  // Obtain "Good track"
                    {
                        //Particle tmp(*i, ptype_k_minus);
                        //k_minus.push_back(tmp);
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

    /*
    X(3872) 771*10^6    1.4e-6~2.11e-7 for J/psi gamma
      gamma psi(2S) (>3%)
         e+e-(7.89e-3)  -> 0.0002367     3rd
         mu+mu-(7.9e-3)  -> 0.0002370    4th
      pi+ pi- J/psi(1S) (>2.6%)
         e+e-(5.971%) -> 0.001552    1st
         mu+ mu-(5.961%) -> 0.001550     2nd
    */

//cout << "combining..." << endl;
    //combination(jpsi,ptype_jpsi,tte_plus,tte_minus);

    combination(jpsi,ptype_jpsi,mu_plus,mu_minus);
    combination(rho0,ptype_rho0,pi_plus,pi_minus);
    //cout << "combining1..." << endl;
    //combination(X3872,ptype_X3872,jpsi,pi_plus,pi_minus);//3 body decay?
    combination(X3872,ptype_X3872,jpsi,rho0);
    //cout << "combining2..." << endl;
    combination(B_cand1,ptype_B0,X3872,gamma);
    //cout << "combined" << endl;
    //combination(B_cand1,ptype_B0,X3872,gamma1);

//for B candicate 21
//cout << "bcand1begin=" << B_cand1<<endl;
//cout << "bcand1begin=" << B_cand1.begin() << "bcand1end" << B_cand1.end() <<endl;
    for(std::vector<Particle>::iterator i =B_cand1.begin(); i != B_cand1.end(); i++)
    {//cout<<"hi"<<endl;
#ifdef BTUPLE
//------------------- dE and Mbc of B meson -------------------------
    //cout<<"dE and Mbc of B meson"<<endl;

		b_tpl->column("expnum",Exp);
        b_tpl->column("runnum",Run);
        
        HepLorentzVector b_p((*i).p() );

        HepLorentzVector boost_vector(-E_HER*sin(cross_angle), 0.0,
                E_LER-E_HER*cos(cross_angle),
                E_HER+E_LER);
        b_p.boost( boost_vector.boostVector() );

        double ebeam = BeamEnergy::E_beam_corr();
        double mass_sqr = ebeam*ebeam - b_p.vect().mag2();
        double mass  = (mass_sqr > 0.) ? sqrt(mass_sqr) :  -sqrt(-mass_sqr);
        float de = b_p.e()-ebeam;


        if (de < -0.5 || de > 0.2) continue;

        b_tpl->column("mbc0", mass);
        b_tpl->column("de", de);
        b_tpl->column("ebeam", BeamEnergy::E_beam_corr());

//------------------- information of B (*i) -------------------------
     //cout<<"information of B"<<endl;
        HepLorentzVector b_pb((*i).p());
        b_pb.boost(boost_vector.boostVector());

        b_tpl->column("pxb", b_pb.px());
        b_tpl->column("pyb", b_pb.py());
        b_tpl->column("pzb", b_pb.pz());
        b_tpl->column("eb" , b_pb.e());
        b_tpl->column("ptb", sqrt((*i).px()*(*i).px()+(*i).py()*(*i).py()));

        Vector3 chb_momentum_cm(b_pb.px(), b_pb.py(), b_pb.pz());


//------------------- information of X(3872) (*i).child(0) -------------------------
     //cout<<"information of x3872"<<endl;
        HepLorentzVector b_p0((*i).child(0).p() );
        b_p0.boost( boost_vector.boostVector() );
        double pxx3872=b_p0.px();
        double pyx3872=b_p0.py();
        double pzx3872=b_p0.pz();
        double massx3872 = (*i).child(0).p().mag();
        //if (massx3872 > 3.95 || massx3872 < 3.7) continue;
        b_tpl->column("massx3872", (*i).child(0).p().mag());
        b_tpl->column("pxx3872", (*i).child(0).px());
        b_tpl->column("pyx3872", (*i).child(0).py());
        b_tpl->column("pzx3872", (*i).child(0).pz());
        double ex3872 = b_p0.e();
        b_tpl->column("ex3872", (*i).child(0).e());

        double cos_polar_0;
        Vector3 ch0_momentum((*i).child(0).px(), (*i).child(0).py(), (*i).child(0).pz());
        Vector3 P_BEAM_polar ( E_HER*sin(cross_angle), 0.0, -E_LER+E_HER*cos(cross_angle));

        cos_polar_0=ch0_momentum.dot(P_BEAM_polar);
        cos_polar_0=cos_polar_0/(ch0_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos0", cos_polar_0);

        double cos_polar_0_cm;
        Vector3 ch0_momentum_cm(b_p0.px(), b_p0.py(), b_p0.pz());

        cos_polar_0_cm=ch0_momentum_cm.dot(P_BEAM_polar);
        cos_polar_0_cm=cos_polar_0_cm/(ch0_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos0cm", cos_polar_0_cm);

        b_tpl->column("ptx3872", sqrt((*i).child(0).px()*(*i).child(0).px()+(*i).child(0).py()*(*i).child(0).py()));

        HepLorentzVector b_ph0((*i).child(0).p());
        b_ph0.boost(-(*i).p().boostVector());
        Vector3 ch0_momentum_h(b_ph0.px(), b_ph0.py(), b_ph0.pz());
        double cos_hel_0;
        cos_hel_0 = ch0_momentum_h.dot(chb_momentum_cm);
        cos_hel_0 = cos_hel_0/(ch0_momentum_h.mag()*chb_momentum_cm.mag());
        b_tpl->column("cosh0", cos_hel_0);

//------------------- information of J/psi from X(3872) (*i).child(0).child(0) -------------------------
    //cout<<"information of jpsi"<<endl;
        HepLorentzVector b_p00((*i).child(0).child(0).p() );
        b_p00.boost( boost_vector.boostVector() );
        double mass00 = (*i).child(0).child(0).p().mag();
        //if (/*massx3872-mass00 > 0.795 || massx3872-mass00 < 0.755 ||*/mass00<3.03||mass00>3.13) continue;
        b_tpl->column("mass00", (*i).child(0).child(0).p().mag());
        b_tpl->column("px00", (*i).child(0).child(0).px());
        b_tpl->column("py00", (*i).child(0).child(0).py());
        b_tpl->column("pz00", (*i).child(0).child(0).pz());
        b_tpl->column("e00", (*i).child(0).child(0).e());


        double cos_polar_00;
        Vector3 ch00_momentum((*i).child(0).child(0).px(), (*i).child(0).child(0).py(), (*i).child(0).child(0).pz());

        cos_polar_00=ch00_momentum.dot(P_BEAM_polar);
        cos_polar_00=cos_polar_00/(ch00_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos00", cos_polar_00);

        double cos_polar_00_cm;
        Vector3 ch00_momentum_cm(b_p00.px(), b_p00.py(), b_p00.pz());

        cos_polar_00_cm=ch00_momentum_cm.dot(P_BEAM_polar);
        cos_polar_00_cm=cos_polar_00_cm/(ch00_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos00cm", cos_polar_00_cm);

        b_tpl->column("pt00", sqrt((*i).child(0).child(0).px()*(*i).child(0).child(0).px()+(*i).child(0).child(0).py()*(*i).child(0).child(0).py()));

        HepLorentzVector b_ph00((*i).child(0).child(0).p());
        b_ph00.boost(-(*i).child(0).p().boostVector());
        Vector3 ch00_momentum_h(b_ph00.px(), b_ph00.py(), b_ph00.pz());
        double cos_hel_00;
        cos_hel_00 = ch00_momentum_h.dot(ch0_momentum_cm);
        cos_hel_00 = cos_hel_00/(ch00_momentum_h.mag()*ch0_momentum_cm.mag());
        b_tpl->column("cosh00", cos_hel_00);

//------------------- information of rho0 from X(3872) (*i).child(0).child(1) -------------------------
    //cout<<"information of rho0"<<endl;
        double rho0mass = (*i).child(0).child(1).p().mag();
        b_tpl->column("rho0mass", rho0mass);
        //if (rho0mass-massx3872+mass00<-0.15) continue;
        double pxrho0 = (*i).child(0).child(1).px();
        double pyrho0 = (*i).child(0).child(1).py();
        double pzrho0 = (*i).child(0).child(1).pz();
        b_tpl->column("pxrho0", pxrho0);
        b_tpl->column("pyrho0", pyrho0);
        b_tpl->column("pzrho0", pzrho0);
        b_tpl->column("rho0e", (*i).child(0).child(1).e());
        b_tpl->column("ptrho0", sqrt((*i).child(0).child(1).px()*(*i).child(0).child(1).px()+(*i).child(0).child(1).py()*(*i).child(0).child(1).py()));

        double cos_polar_rho0;
        Vector3 rho0_momentum1((*i).child(0).child(1).px(), (*i).child(0).child(1).py(), (*i).child(0).child(1).pz());
  
        cos_polar_rho0=rho0_momentum1.dot(P_BEAM_polar);
        cos_polar_rho0=cos_polar_rho0/(rho0_momentum1.mag()*P_BEAM_polar.mag());
        b_tpl->column("rho0cos", cos_polar_rho0);

        HepLorentzVector b_prho0((*i).child(0).child(1).p() );
        b_prho0.boost( boost_vector.boostVector() );
        double cos_polar_rho0_cm;
        Vector3 rho0_momentum_cm1(b_prho0.px(), b_prho0.py(), b_prho0.pz());

        cos_polar_rho0_cm=rho0_momentum_cm1.dot(P_BEAM_polar);
        cos_polar_rho0_cm=cos_polar_rho0_cm/(rho0_momentum_cm1.mag()*P_BEAM_polar.mag());
        b_tpl->column("rhocoscm", cos_polar_rho0_cm);

        HepLorentzVector b_phrho0((*i).child(0).child(1).p());
        b_phrho0.boost(-(*i).child(0).p().boostVector());
        Vector3 chrho0_momentum_h(b_phrho0.px(), b_phrho0.py(), b_phrho0.pz());
        double cos_hel_rho0;
        cos_hel_rho0 = chrho0_momentum_h.dot(ch0_momentum_cm);
        cos_hel_rho0 = cos_hel_rho0/(chrho0_momentum_h.mag()*ch0_momentum_cm.mag());
        b_tpl->column("coshrho0", cos_hel_rho0);

//------------------- information of pi+ from rho0 (*i).child(0).child(1).child(0) -------------------------
    //cout<<"information of pi+"<<endl;
        b_tpl->column("mass01", (*i).child(0).child(1).child(0).p().mag());
        double px01 = (*i).child(0).child(1).child(0).px();
        double py01 = (*i).child(0).child(1).child(0).py();
        double pz01 = (*i).child(0).child(1).child(0).pz();
        b_tpl->column("px01", px01);
        b_tpl->column("py01", py01);
        b_tpl->column("pz01", pz01);

        b_tpl->column("e01", (*i).child(0).child(1).child(0).e());

        b_tpl->column("pt01", sqrt((*i).child(0).child(1).child(0).px()*(*i).child(0).child(1).child(0).px()+(*i).child(0).child(1).child(0).py()*(*i).child(0).child(1).child(0).py()));
        const Mdst_charged* chk01 = &(*i).child(0).child(1).child(0).mdstCharged();
        float kpi_id01 = selkpi.prob(chk01);
        b_tpl->column("pid01", kpi_id01);
        //float kpi_idp01_muk = selmuk.prob(chk01);
        //b_tpl->column("id01muk", kpi_idp01_muk);
        //     
        //float kpi_idp01_mupi = selmupi.prob(chk01);
        //b_tpl->column("id01mupi", kpi_idp01_mupi);

        float chrg01 = chk01->charge();
        b_tpl->column("chrg01", chrg01);

        double dr01,dz01;
        GetImpactParameters(chk01,&dr01,&dz01,2);//STAR
        b_tpl->column("dr01", float(dr01));
        b_tpl->column("dz01", float(dz01));

        double cos_polar_01;
        Vector3 ch01_momentum((*i).child(0).child(1).child(0).px(), (*i).child(0).child(1).child(0).py(), (*i).child(0).child(1).child(0).pz());
    
        cos_polar_01=ch01_momentum.dot(P_BEAM_polar);
        cos_polar_01=cos_polar_01/(ch01_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos01", cos_polar_01);

        HepLorentzVector b_p01((*i).child(0).child(1).child(0).p() );
        b_p01.boost( boost_vector.boostVector() );
        double cos_polar_01_cm;
        Vector3 ch01_momentum_cm(b_p01.px(), b_p01.py(), b_p01.pz());

        cos_polar_01_cm=ch01_momentum_cm.dot(P_BEAM_polar);
        cos_polar_01_cm=cos_polar_01_cm/(ch01_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos01cm", cos_polar_01_cm);

        HepLorentzVector b_ph01((*i).child(0).child(1).child(0).p());
        b_ph01.boost(-(*i).child(0).child(1).p().boostVector());
        Vector3 ch01_momentum_h(b_ph01.px(), b_ph01.py(), b_ph01.pz());
        double cos_hel_01;
        cos_hel_01 = ch01_momentum_h.dot(rho0_momentum_cm1);
        cos_hel_01 = cos_hel_01/(ch01_momentum_h.mag()*rho0_momentum_cm1.mag());
        b_tpl->column("cosh01", cos_hel_01);

//------------------- information of pi- from X(3872) (*i).child(0).child(1).child(1) -------------------------
    //cout<<"information of pi-"<<endl;
        b_tpl->column("mass02", (*i).child(0).child(1).child(1).p().mag());
        double px02 = (*i).child(0).child(1).child(1).px();
        double py02 = (*i).child(0).child(1).child(1).py();
        double pz02 = (*i).child(0).child(1).child(1).pz();
        b_tpl->column("px02", px02);
        b_tpl->column("py02", py02);
        b_tpl->column("pz02", pz02);
        b_tpl->column("e02", (*i).child(0).child(1).child(1).e());

        b_tpl->column("pt02", sqrt((*i).child(0).child(1).child(1).px()*(*i).child(0).child(1).child(1).px()+(*i).child(0).child(1).child(1).py()*(*i).child(0).child(1).child(1).py()));
        const Mdst_charged* chk02 = &(*i).child(0).child(1).child(1).mdstCharged();
        float kpi_id02 = selkpi.prob(chk02);
        b_tpl->column("pid02", kpi_id02);
        //float kpi_idp02_muk = selmuk.prob(chk02);
        //b_tpl->column("id02muk", kpi_idp02_muk);
        //     
        //float kpi_idp02_mupi = selmupi.prob(chk02);
        //b_tpl->column("id02mupi", kpi_idp02_mupi);
        float chrg02 = chk02->charge();
        b_tpl->column("chrg02", chrg02);

        double dr02,dz02;
        GetImpactParameters(chk02,&dr02,&dz02,2);//STAR
        b_tpl->column("dr02", float(dr02));
        b_tpl->column("dz02", float(dz02));

        double cos_polar_02;
        Vector3 ch02_momentum((*i).child(0).child(1).child(1).px(), (*i).child(0).child(1).child(1).py(), (*i).child(0).child(1).child(1).pz());

        cos_polar_02=ch02_momentum.dot(P_BEAM_polar);
        cos_polar_02=cos_polar_02/(ch02_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos02", cos_polar_02);
        
        HepLorentzVector b_p02((*i).child(0).child(1).child(1).p() );
        b_p02.boost( boost_vector.boostVector() );
        double cos_polar_02_cm;
        Vector3 ch02_momentum_cm(b_p02.px(), b_p02.py(), b_p02.pz());

        cos_polar_02_cm=ch02_momentum_cm.dot(P_BEAM_polar);
        cos_polar_02_cm=cos_polar_02_cm/(ch02_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos02cm", cos_polar_02_cm);

        HepLorentzVector b_ph02((*i).child(0).child(1).child(1).p());
        b_ph02.boost(-(*i).child(0).child(1).p().boostVector());
        Vector3 ch02_momentum_h(b_ph02.px(), b_ph02.py(), b_ph02.pz());
        double cos_hel_02;
        cos_hel_02 = ch02_momentum_h.dot(rho0_momentum_cm1);
        cos_hel_02 = cos_hel_02/(ch02_momentum_h.mag()*rho0_momentum_cm1.mag());
        b_tpl->column("cosh02", cos_hel_02);


        //if (abs(dr01)>2||abs(dr02)>2||abs(dz01)>5||abs(dz02)>5||px01*px01+py01*py01+pz01*pz01<0.01||px02*px02+py02*py02+pz02*pz02<0.01) continue;
//------------------- information of e+ (mu+) from Jpsi (*i).child(0).child(0).child(0) -------------------------
    //cout<<"information of mu+"<<endl;
        b_tpl->column("mass000", (*i).child(0).child(0).child(0).p().mag());
        b_tpl->column("px000", (*i).child(0).child(0).child(0).px());
        b_tpl->column("py000", (*i).child(0).child(0).child(0).py());
        b_tpl->column("pz000", (*i).child(0).child(0).child(0).pz());
        b_tpl->column("e000", (*i).child(0).child(0).child(0).e());

        b_tpl->column("pt000", sqrt((*i).child(0).child(0).child(0).px()*(*i).child(0).child(0).child(0).px()+(*i).child(0).child(0).child(0).py()*(*i).child(0).child(0).child(0).py()));
        const Mdst_charged* chk000 = &(*i).child(0).child(0).child(0).mdstCharged();
        float kpi_id000 = selkpi.prob(chk000);
        b_tpl->column("pid000", kpi_id000);
        //float kpi_idp000_muk = selmuk.prob(chk000);
        //b_tpl->column("id000muk", kpi_idp000_muk);
        //     
        //float kpi_idp000_mupi = selmupi.prob(chk000);
        //b_tpl->column("id000mpi", kpi_idp000_mupi);
        float chrg000 = chk000->charge();
        b_tpl->column("chrg000", chrg000);

        double dr000,dz000;
        GetImpactParameters(chk000,&dr000,&dz000,1);//STAR
        b_tpl->column("dr000", float(dr000));
        b_tpl->column("dz000", float(dz000));

        double cos_polar_000;
        Vector3 ch000_momentum((*i).child(0).child(0).child(0).px(), (*i).child(0).child(0).child(0).py(), (*i).child(0).child(0).child(0).pz());

        cos_polar_000=ch000_momentum.dot(P_BEAM_polar);
        cos_polar_000=cos_polar_000/(ch000_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos000", cos_polar_000);

        HepLorentzVector b_p000((*i).child(0).child(0).child(0).p() );
        b_p000.boost( boost_vector.boostVector() );
        double cos_polar_000_cm;
        Vector3 ch000_momentum_cm(b_p000.px(), b_p000.py(), b_p000.pz());

        cos_polar_000_cm=ch000_momentum_cm.dot(P_BEAM_polar);
        cos_polar_000_cm=cos_polar_000_cm/(ch000_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos000cm", cos_polar_000_cm);

        HepLorentzVector b_ph000((*i).child(0).child(0).child(0).p());
        b_ph000.boost(-(*i).child(0).child(0).p().boostVector());
        Vector3 ch000_momentum_h(b_ph000.px(), b_ph000.py(), b_ph000.pz());
        double cos_hel_000;
        cos_hel_000 = ch000_momentum_h.dot(ch00_momentum_cm);
        cos_hel_000 = cos_hel_000/(ch000_momentum_h.mag()*ch00_momentum_cm.mag());
        b_tpl->column("cosh000", cos_hel_000);

//------------------- information of e- (mu-) from Jpsi (*i).child(0).child(0).child(1) -------------------------
    //cout<<"information of mu-"<<endl;
        b_tpl->column("mass001", (*i).child(0).child(0).child(1).p().mag());
        b_tpl->column("px001", (*i).child(0).child(0).child(1).px());
        b_tpl->column("py001", (*i).child(0).child(0).child(1).py());
        b_tpl->column("pz001", (*i).child(0).child(0).child(1).pz());
        b_tpl->column("e001", (*i).child(0).child(0).child(1).e());

        b_tpl->column("pt001", sqrt((*i).child(0).child(0).child(1).px()*(*i).child(0).child(0).child(1).px()+(*i).child(0).child(0).child(1).py()*(*i).child(0).child(0).child(1).py()));
        const Mdst_charged* chk001 = &(*i).child(0).child(0).child(1).mdstCharged();
        float kpi_id001 = selkpi.prob(chk001);
        b_tpl->column("pid001", kpi_id001);
        //float kpi_idp001_muk = selmuk.prob(chk001);
        //b_tpl->column("id001muk", kpi_idp001_muk);
        //     
        //float kpi_idp001_mupi = selmupi.prob(chk001);
        //b_tpl->column("id001mpi", kpi_idp001_mupi);
        float chrg001 = chk001->charge();
        b_tpl->column("chrg001", chrg001);

        double dr001,dz001;
        GetImpactParameters(chk001,&dr001,&dz001,1);//STAR
        b_tpl->column("dr001", float(dr001));
        b_tpl->column("dz001", float(dz001));
        //if (abs(dr000)>0.2||abs(dr001)>0.2||abs(dz000)>2||abs(dz001)>2) continue;
        double cos_polar_001;
        Vector3 ch001_momentum((*i).child(0).child(0).child(1).px(), (*i).child(0).child(0).child(1).py(), (*i).child(0).child(0).child(1).pz());

        cos_polar_001=ch001_momentum.dot(P_BEAM_polar);
        cos_polar_001=cos_polar_001/(ch001_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos001", cos_polar_001);

        HepLorentzVector b_p001((*i).child(0).child(0).child(1).p() );
        b_p001.boost( boost_vector.boostVector() );
        double cos_polar_001_cm;
        Vector3 ch001_momentum_cm(b_p001.px(), b_p001.py(), b_p001.pz());

        cos_polar_001_cm=ch001_momentum_cm.dot(P_BEAM_polar);
        cos_polar_001_cm=cos_polar_001_cm/(ch001_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos001cm", cos_polar_001_cm);

        HepLorentzVector b_ph001((*i).child(0).child(0).child(1).p());
        b_ph001.boost(-(*i).child(0).child(0).p().boostVector());
        Vector3 ch001_momentum_h(b_ph001.px(), b_ph001.py(), b_ph001.pz());
        double cos_hel_001;
        cos_hel_001 = ch001_momentum_h.dot(ch00_momentum_cm);
        cos_hel_001 = cos_hel_001/(ch001_momentum_h.mag()*ch00_momentum_cm.mag());
        b_tpl->column("cosh001", cos_hel_001);
//------------------- information of gamma from B (*i).child(1) -------------------------
    //cout<<"information of gamma"<<endl;
        HepLorentzVector b_p1((*i).child(1).p() );
        b_p1.boost( boost_vector.boostVector() );
        //if (b_p1.e()<0.6) continue;
        b_tpl->column("mass1", (*i).child(1).p().mag());
        double px1=b_p1.px();
        double py1=b_p1.py();
        double pz1=b_p1.pz();
        b_tpl->column("px1", (*i).child(1).px());
        b_tpl->column("py1", (*i).child(1).py());
        b_tpl->column("pz1", (*i).child(1).pz());
        b_tpl->column("e1", (*i).child(1).e());
        b_tpl->column("pt1", sqrt((*i).child(1).px()*(*i).child(1).px()+(*i).child(1).py()*(*i).child(1).py()));

        //Gamma cannot get charge and impact parameter (dr, dz).

        double cos_polar_1;
        Vector3 ch1_momentum((*i).child(1).px(), (*i).child(1).py(), (*i).child(1).pz());

        cos_polar_1=ch1_momentum.dot(P_BEAM_polar);
        cos_polar_1=cos_polar_1/(ch1_momentum.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos1", cos_polar_1);

        double cos_polar_1_cm;
        Vector3 ch1_momentum_cm(b_p1.px(), b_p1.py(), b_p1.pz());

        cos_polar_1_cm=ch1_momentum_cm.dot(P_BEAM_polar);
        cos_polar_1_cm=cos_polar_1_cm/(ch1_momentum_cm.mag()*P_BEAM_polar.mag());
        b_tpl->column("cos1cm", cos_polar_1_cm);

        HepLorentzVector b_ph1((*i).child(1).p());
        b_ph1.boost(-(*i).p().boostVector());
        Vector3 ch1_momentum_h(b_ph1.px(), b_ph1.py(), b_ph1.pz());
        double cos_hel_1;
        cos_hel_1 = ch1_momentum_h.dot(chb_momentum_cm);
        cos_hel_1 = cos_hel_1/(ch1_momentum_h.mag()*chb_momentum_cm.mag());
        b_tpl->column("cosh1", cos_hel_1);

//------------------- Modify Mbc -------------------------
        double temp1 = (ebeam-ex3872)/sqrt(px1*px1+py1*py1+pz1*pz1);
        double temp2 = (pxx3872+px1*temp1)*(pxx3872+px1*temp1)+(pyx3872+py1*temp1)*(pyx3872+py1*temp1)+(pzx3872+pz1*temp1)*(pzx3872+pz1*temp1);
        double mass_correct = sqrt(ebeam*ebeam-temp2);
        b_tpl->column("mbc", mass_correct);
        if (mass_correct<5.2) continue;

//------------------- E9/E25 -------------------------
        Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
        Mdst_ecl_aux &aux = eclaux_mgr(Panther_ID((*i).child(1).mdstGamma().ecl().get_ID()));

        double e9oe25 = aux.e9oe25();
        //if (e9oe25<0.87) continue;
        b_tpl->column("e9oe25",e9oe25); 

//------------------- sumecl -------------------------
        double sum_ecl = 0;

        for(std::vector<Mdst_ecl_aux>::iterator s = eclaux_mag.begin(); s !=eclaux_mag.end();s++ )
        {
            Mdst_ecl_aux& ch = *s;
            for(std::vector<Mdst_ecl>::iterator p = ecl_mag.begin(); p!=ecl_mag.end(); p++)
            {

                if ((*s).get_ID()==(*p).get_ID())
                {
                  //forward e>0.1
                  if ((*p).theta() < 0.562 && (*p).energy() > 0.1 )
                      sum_ecl = sum_ecl + (*p).energy();
                  //barrel e>0.05
                  if ((*p).theta() > 0.562 && (*p).theta() < 2.246 && (*p).energy()>0.05 )
                      sum_ecl = sum_ecl + (*p).energy();
                  //backward e>0.15
                  if ((*p).theta() > 2.246 && (*p).energy() > 0.15 )
                      sum_ecl = sum_ecl + (*p).energy();
                }
            }
        }//end ecl
        b_tpl->column("sumecl",sum_ecl); 

//------------------- J/psi vertex fit with K fitter  -------------------------
   //cout<<"J/psi vertex fit with K fitter"<<endl;
        kvertexfitter kvf;
        addTrack2fit(kvf, (*i).child(0).child(0).child(0)); // add "mu+(-)" to kfitter
        addTrack2fit(kvf, (*i).child(0).child(0).child(1)); // add "mu-(+)" to kfitter
        unsigned err = kvf.fit(); // do "fitting"
        float vchisq;
        if(err == 0) // success.
              vchisq = kvf.chisq();
        else vchisq = -1.;
        //if (vchisq>20||vchisq<0) continue;
        b_tpl->column("vchisq", vchisq);

//------------------- PiPi vertex fit with K fitter  -------------------------
   //cout<<"PiPi vertex fit with K fitter"<<endl;
        kvertexfitter kvf1;
        addTrack2fit(kvf1, (*i).child(0).child(1).child(0)); // add "pi+(-)" to kfitter
        addTrack2fit(kvf1, (*i).child(0).child(1).child(1)); // add "pi-(+)" to kfitter
        unsigned err1 = kvf1.fit(); // do "fitting"
        float vchisq1;
        if(err1 == 0) // success.
              vchisq1 = kvf1.chisq();
        else vchisq1 = -1.;
        //if (vchisq1>80||vchisq1<0) continue;
        b_tpl->column("vchisq1", vchisq1);

//------------------- X(3872) vertex fit with ExKfitter  -------------------------
        
   //cout<<"X(3872) vertex fit with ExKfitter"<<endl;
        Mdst_charged Mdst_01 =(*i).child(0).child(1).child(0).mdstCharged();
        Mdst_charged Mdst_02 =(*i).child(0).child(1).child(1).mdstCharged();
        Mdst_charged Mdst_000 =(*i).child(0).child(0).child(0).mdstCharged();
        Mdst_charged Mdst_001 =(*i).child(0).child(0).child(1).mdstCharged();
        //Mdst_gamma Mdst_1 =(*i).child(1).mdstGamma();

        ExKFitterParticle KF_000(Mdst_000, 1);
        ExKFitterParticle KF_001(Mdst_001, 1);//STAR e mu pi k b
        ExKFitterParticle KF_02(Mdst_02, 2);
        ExKFitterParticle KF_01(Mdst_01, 2);
        //ExKFitterParticle KF_1(Mdst_1);
        
        //ExKFitterMass B_Mass(5.27958); <- Mass Constraint Fit

        HepPoint3D X3872_init;
        HepPoint3D JPSI_init;
        HepPoint3D rho0_init;

        X3872_init.setX(IP.x() + (*i).child(0).p().px()/(*i).child(0).p().rho());
        X3872_init.setY(IP.y() + (*i).child(0).p().py()/(*i).child(0).p().rho());
        X3872_init.setZ(IP.z() + (*i).child(0).p().pz()/(*i).child(0).p().rho());
        ExKFitterVertex X3872_Vertex(X3872_init);
        JPSI_init.setX(IP.x() + (*i).child(0).child(0).p().px()/(*i).child(0).child(0).p().rho());
        JPSI_init.setY(IP.y() + (*i).child(0).child(0).p().py()/(*i).child(0).child(0).p().rho());
        JPSI_init.setZ(IP.z() + (*i).child(0).child(0).p().pz()/(*i).child(0).child(0).p().rho());
        ExKFitterVertex JPSI_Vertex(JPSI_init);
        rho0_init.setX(IP.x() + (*i).child(0).child(1).p().px()/(*i).child(0).child(1).p().rho());
        rho0_init.setY(IP.y() + (*i).child(0).child(1).p().py()/(*i).child(0).child(1).p().rho());
        rho0_init.setZ(IP.z() + (*i).child(0).child(1).p().pz()/(*i).child(0).child(1).p().rho());
        ExKFitterVertex RHO0_Vertex(rho0_init);

        ExKFitterVertex B_Vertex(IP,IPerr);
        
        //cout<<"con1"<<endl;
        ExKFitterParticle JPSI;
        JPSI.LinkParticle(&KF_000);
        JPSI.LinkParticle(&KF_001);
        JPSI.LinkVertex(&JPSI_Vertex);

        ExKFitterConstrain con1;
        con1.SetVertexConstrain();
        con1.LinkParticle(&KF_000);
        con1.LinkParticle(&KF_001);
        con1.LinkVertex(&JPSI_Vertex);

        ExKFitterParticle RHO0;
        RHO0.LinkParticle(&KF_01);
        RHO0.LinkParticle(&KF_02);
        RHO0.LinkVertex(&RHO0_Vertex);

        ExKFitterConstrain con11;
        con11.SetVertexConstrain();
        con11.LinkParticle(&KF_01);
        con11.LinkParticle(&KF_02);
        con11.LinkVertex(&RHO0_Vertex);

        //cout<<"con2"<<endl;
        ExKFitterParticle X3872;
        X3872.LinkParticle(&JPSI);
        X3872.LinkParticle(&RHO0);
        //X3872.LinkParticle(&KF_02);
        X3872.LinkVertex(&X3872_Vertex);

        ExKFitterConstrain con2;
        con2.SetVertexConstrain();
        con2.LinkParticle(&JPSI);
        con2.LinkParticle(&RHO0);
        //con2.LinkParticle(&KF_02);
        con2.LinkVertex(&X3872_Vertex);

        //cout<<"con3"<<endl;
        //ExKFitterParticle B;
        //B.LinkParticle(&X3872);
        //B.LinkParticle(&KF_1);
        //B.LinkVertex(&B_Vertex);
        

        //ExKFitterConstrain con3;
//
        //con3.SetMassConstrain();
        //con3.LinkParticle(&X3872);
        //con3.LinkParticle(&KF_1);
        //con3.LinkVertex(&B_Vertex);
        //con3.LinkMass(&B_Mass);

        //ExKFitterConstrain con11=con1;
        //ExKFitterConstrain con12=con1;
        //ExKFitterConstrain con21=con2;
        //cout<<"Core"<<endl;
        ExKFitter Core;
        Core.LinkConstrain(&con1);
        Core.LinkConstrain(&con11);
        Core.LinkConstrain(&con2);
        //Core.LinkConstrain(&con3);
        int ret = Core.Minimize();
        float chisqExK = Core.Chisq();
        float dof_exk = Core.N_DegreeOfFreedom();

        //ExKFitter Core0;
        //Core0.LinkConstrain(&con11);
        //
        //int ret0 = Core0.Minimize();
        //float chisqExK0 = Core0.Chisq();
        //float dof_exk0 = Core0.N_DegreeOfFreedom();

        //ExKFitter Core00;
        //Core00.LinkConstrain(&con12);
        //Core00.LinkConstrain(&con21);

        //int ret00 = Core00.Minimize();
        //float chisqExK00 = Core00.Chisq();
        //float dof_exk00 = Core00.N_DegreeOfFreedom();

        //HepPoint3D jvertex = JPSI_Vertex.Vertex();
        HepPoint3D bvertex = B_Vertex.Vertex();
        if(ret==0){X3872.Update();}
        //if(ret0==0){JPSI.Update();}
        //if(ret00==0){JPSI.Update();}
        //if (chisqExK/dof_exk>100||chisqExK/dof_exk<0) continue;
        b_tpl->column("chisqexk",chisqExK/dof_exk);
        //b_tpl->column("chisqex0",chisqExK0/dof_exk0);
        //b_tpl->column("chisqe00",chisqExK00/dof_exk00);
        b_tpl->column("evtcount",evtcount);


        //ExKFitterVertex B_Vertex(IP,IPerr);
       // HepPoint3D bvertex = B_Vertex.Vertex();
/*
        //------------------- pipi_combination (write by myself)  -------------------------
        //cout<<"pipi_combination"<<endl;
        std::vector<Particle> pi1;
        std::vector<Particle> pi2;

        Particle tmp1((*i).child(0).child(1));
        Particle tmp2((*i).child(0).child(2));
        pi1.push_back(tmp1);
        pi2.push_back(tmp2);

        combination(rho0,ptype_rho0,pi1,pi2);

        double massrho0,erho0,cos_polar_rho0,cos_polar_rho0_cm;
        Vector3 rho0_momentum1,rho0_momentum_cm1;

        for(std::vector<Particle>::iterator i1 =rho0.begin(); i1 != rho0.end(); i1++)
        {
            HepLorentzVector rho0_p((*i1).p() );
            rho0_p.boost( boost_vector.boostVector() );
            massrho0 = rho0_p.mag();
            erho0 = rho0_p.e();

            Vector3 rho0_momentum((*i1).px(), (*i1).py(), (*i1).pz());
            rho0_momentum1=rho0_momentum;

            cos_polar_rho0=rho0_momentum.dot(P_BEAM_polar);
            cos_polar_rho0=cos_polar_rho0/(rho0_momentum.mag()*P_BEAM_polar.mag());
    
            Vector3 rho0_momentum_cm(rho0_p.px(), rho0_p.py(), rho0_p.pz());
            rho0_momentum_cm1=rho0_momentum_cm;

            cos_polar_rho0_cm=rho0_momentum_cm.dot(P_BEAM_polar);
            cos_polar_rho0_cm=cos_polar_rho0_cm/(rho0_momentum_cm.mag()*P_BEAM_polar.mag());

        }
        b_tpl->column("rho0mass", massrho0);
        b_tpl->column("rho0e", erho0);
        b_tpl->column("rho0cos", cos_polar_rho0);
        b_tpl->column("rhocoscm", cos_polar_rho0_cm);
*/
//------------------- pi0 & eta0 veto  -------------------------
   //cout<<"pi0 & eta0 veto"<<endl;
        //double Closest_Pi0_Probability(Mdst_gamma& IN, Mdst_gamma& OUT, double& m);
        //double Closest_Eta_Probability(Mdst_gamma& IN,  Mdst_gamma& OUT, double& m); 
        //out<<"ping1"<<endl;
        Mdst_gamma gamma2,gamma3;
        double mpi0,meta0;
        double probpi0 = Closest_Pi0_Probability((*i).child(1).mdstGamma(), gamma2, mpi0);
        //if (probpi0>0.3) continue;
        double probeta0 = Closest_Eta_Probability((*i).child(1).mdstGamma(),  gamma3, meta0); 
        b_tpl->column("probpi0", probpi0);
        b_tpl->column("probeta0", probeta0);
        b_tpl->column("mpi0", mpi0);
        b_tpl->column("meta0", meta0);

        //cout<<"ping2"<<endl;
        //cout<<"probpi0 = "<<probpi0<<" , probeta0 = "<<probeta0<<endl;
        //cout<<"mpi0 = "<<mpi0<<" , meta0 = "<<meta0<<endl;
        //cout<<"gamma2 = "<<gamma2<<" , gamma3 = "<<gamma3<<endl;
        Vector3 pi0_gamma2_momentum,pi0_gamma2_momentum_cm,eta0_gamma2_momentum,eta0_gamma2_momentum_cm;
        if (gamma2){
            double pie2 = gamma2.ecl().energy();
            b_tpl->column("pie2", pie2);
            HepLorentzVector Ppi0_gamma2( gamma2.px(), gamma2.py(), gamma2.pz());
            Ppi0_gamma2.boost( boost_vector.boostVector() );
            Vector3 pi0_gamma2_momentum1( gamma2.px(), gamma2.py(), gamma2.pz());
            pi0_gamma2_momentum=pi0_gamma2_momentum1;
            Vector3 pi0_gamma2_momentum_cm1(Ppi0_gamma2.px(), Ppi0_gamma2.py(), Ppi0_gamma2.pz());
            pi0_gamma2_momentum_cm=pi0_gamma2_momentum_cm1;
        }else{
            b_tpl->column("pie2", -1);
        }
        
        if (gamma3){
            double etae2 = gamma3.ecl().energy();
            b_tpl->column("etae2", etae2);
            HepLorentzVector Peta0_gamma2( gamma3.px(), gamma3.py(), gamma3.pz());
            Peta0_gamma2.boost( boost_vector.boostVector() );
            Vector3 eta0_gamma2_momentum1( gamma3.px(), gamma3.py(), gamma3.pz());
            eta0_gamma2_momentum=eta0_gamma2_momentum1;
            Vector3 eta0_gamma2_momentum_cm1(Peta0_gamma2.px(), Peta0_gamma2.py(), Peta0_gamma2.pz());
            eta0_gamma2_momentum_cm=eta0_gamma2_momentum_cm1;
        }else{
            b_tpl->column("etae2", -1);
        }
        //cout<<"ping3"<<endl;


//------------------- Cosine Angles  -------------------------
   //cout<<"Cosine Angles"<<endl;
        double cosll, cosjl, cos2pi, cosrpi, cosrj, cosxj, cosxg, cospgg, cosegg;
        double cosllcm, cosjlcm, cos2picm, cosrpicm, cosrjcm, cosxjcm, cosxgcm, cospggcm, coseggcm;

        //Angle of mu+ mu- vertex
        //cout<<"pong1"<<endl;
        cosll=ch000_momentum.dot(ch001_momentum)/(ch000_momentum.mag()*ch001_momentum.mag());
        b_tpl->column("cosll",cosll);
        cosllcm=ch000_momentum_cm.dot(ch001_momentum_cm)/(ch000_momentum_cm.mag()*ch001_momentum_cm.mag());
        b_tpl->column("cosllcm",cosllcm);

        //cout<<"pong2"<<endl;
        cosjl=ch000_momentum.dot(ch00_momentum)/(ch000_momentum.mag()*ch00_momentum.mag());
        b_tpl->column("cosjl",cosjl);
        cosjlcm=ch000_momentum_cm.dot(ch00_momentum_cm)/(ch000_momentum_cm.mag()*ch00_momentum_cm.mag());
        b_tpl->column("cosjlcm",cosjlcm);

        //Angle of pi+ pi- vertex
        //cout<<"pong3"<<endl;
        cos2pi=ch01_momentum.dot(ch02_momentum)/(ch01_momentum.mag()*ch02_momentum.mag());
        b_tpl->column("cos2pi",cos2pi);
        cos2picm=ch01_momentum_cm.dot(ch02_momentum_cm)/(ch01_momentum_cm.mag()*ch02_momentum_cm.mag());
        b_tpl->column("cos2picm",cos2picm);

        //cout<<"pong4"<<endl;
        cosrpi=ch01_momentum.dot(rho0_momentum1)/(ch01_momentum.mag()*rho0_momentum1.mag());
        b_tpl->column("cosrpi",cosrpi);
        cosrpicm=ch01_momentum_cm.dot(rho0_momentum_cm1)/(ch01_momentum_cm.mag()*rho0_momentum_cm1.mag());
        b_tpl->column("cosrpicm",cosrpicm);

        //Angle of J/psi rho0 vertex
        //cout<<"pong5"<<endl;
        cosrj=rho0_momentum1.dot(ch00_momentum)/(rho0_momentum1.mag()*ch00_momentum.mag());
        b_tpl->column("cosrj",cosrj);
        cosrjcm=rho0_momentum_cm1.dot(ch00_momentum_cm)/(rho0_momentum_cm1.mag()*ch00_momentum_cm.mag());
        b_tpl->column("cosrjcm",cosrjcm);

        //cout<<"pong6"<<endl;
        cosxj=ch0_momentum.dot(ch00_momentum)/(ch0_momentum.mag()*ch00_momentum.mag());
        b_tpl->column("cosxj",cosxj);
        cosxjcm=ch0_momentum_cm.dot(ch00_momentum_cm)/(ch0_momentum_cm.mag()*ch00_momentum_cm.mag());
        b_tpl->column("cosxjcm",cosxjcm);
        //if (cosxjcm<0.65||cosllcm>-0.5) continue;
        //Angle of X(3872) gamma vertex
        //cout<<"pong7"<<endl;
        cosxg=ch0_momentum.dot(ch1_momentum)/(ch0_momentum.mag()*ch1_momentum.mag());
        b_tpl->column("cosxg",cosxg);
        cosxgcm=ch0_momentum_cm.dot(ch1_momentum_cm)/(ch0_momentum_cm.mag()*ch1_momentum_cm.mag());
        b_tpl->column("cosxgcm",cosxgcm);

        //// Angle between X(3872) and B0 ignored.

        //Angle of pi0 vertex
        //cout<<"pong8"<<endl;
        if (gamma2){
            cospgg=ch1_momentum.dot(pi0_gamma2_momentum)/(ch1_momentum.mag()*pi0_gamma2_momentum.mag());
            b_tpl->column("cospgg",cospgg);
            cospggcm=ch1_momentum_cm.dot(pi0_gamma2_momentum_cm)/(ch1_momentum_cm.mag()*pi0_gamma2_momentum_cm.mag());
            b_tpl->column("cospggcm",cospggcm);
        }else{
            b_tpl->column("cospgg",-2);
            b_tpl->column("cospggcm",-2);
        }

        //Angle of eta0 vertex
        //cout<<"pong9"<<endl;
        if (gamma3){
            cosegg=ch1_momentum.dot(eta0_gamma2_momentum)/(ch1_momentum.mag()*eta0_gamma2_momentum.mag());
            b_tpl->column("cosegg",cosegg);
            coseggcm=ch1_momentum_cm.dot(eta0_gamma2_momentum_cm)/(ch1_momentum_cm.mag()*eta0_gamma2_momentum_cm.mag());
            b_tpl->column("coseggcm",coseggcm);
        }else{
            b_tpl->column("cosegg",-2);
            b_tpl->column("coseggcm",-2);
        }


//------------------- k_sfw variables  -------------------------
 //cout<<"k_sfw variables"<<endl;
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

//------------------------- MC truth matching ------------------------
        //cout<<"MC truth matching"<<endl;
        int hindex = 0;
        if (MCstatus == 1)
        {
         //cout<<"MCstatus=1"<<endl;
            const Mdst_charged ch_01 = (*i).child(0).child(1).child(0).mdstCharged();//pi+
            const Mdst_charged ch_02 = (*i).child(0).child(1).child(1).mdstCharged();//pi-
            const Mdst_charged ch_000 = (*i).child(0).child(0).child(0).mdstCharged();//e+
            const Mdst_charged ch_001 = (*i).child(0).child(0).child(1).mdstCharged();//e-
            const Mdst_gamma ch_1 = (*i).child(1).mdstGamma();//gamma
            Gen_hepevt Evtch_02=get_hepevt(ch_02);
            Gen_hepevt Evtch_01=get_hepevt(ch_01);
            Gen_hepevt Evtch_000=get_hepevt(ch_000);
            Gen_hepevt Evtch_001=get_hepevt(ch_001);
            Gen_hepevt Evtch_1=get_hepevt(ch_1);
            Gen_hepevt EvtP_02 ;
            Gen_hepevt EvtP_01 ;
            Gen_hepevt EvtP_000 ;
            Gen_hepevt EvtP_001 ;
            Gen_hepevt EvtP_1 ;
            Gen_hepevt EvtGP_02 ;
            Gen_hepevt EvtGP_01 ;
            Gen_hepevt EvtGP_000 ;
            Gen_hepevt EvtGP_001 ;
            Gen_hepevt EvtGP_1 ;
            Gen_hepevt EvtGGP_1 ;
            Gen_hepevt EvtGGGP_1 ;
            Gen_hepevt EvtGGGGP_1 ;
            Gen_hepevt EvtGGGGGP_1 ;
            Gen_hepevt EvtGGP_000 ;
            Gen_hepevt EvtGGP_001 ;
            Gen_hepevt EvtGGP_02 ;
            Gen_hepevt EvtGGP_01 ;
            
            b_tpl->column("hi02", Evtch_02.idhep());
            b_tpl->column("hi01", Evtch_01.idhep());
            b_tpl->column("hi000", Evtch_000.idhep());
            b_tpl->column("hi001", Evtch_001.idhep());
            b_tpl->column("hi1", Evtch_1.idhep());

            Gen_hepevt EvtY4S = gen_mgr[0];

            double pxi02 = Evtch_02.PX(), pyi02 = Evtch_02.PY(), pzi02 = Evtch_02.PZ();
            double pxi01 = Evtch_01.PX(), pyi01 = Evtch_01.PY(), pzi01 = Evtch_01.PZ();
            b_tpl->column("pxi02", pxi02);
            b_tpl->column("pyi02", pyi02);
            b_tpl->column("pzi02", pzi02);
            b_tpl->column("pxi01", pxi01);
            b_tpl->column("pyi01", pyi01);
            b_tpl->column("pzi01", pzi01);
            double pxi000 = Evtch_000.PX(), pyi000 = Evtch_000.PY(), pzi000 = Evtch_000.PZ();
            double pxi001 = Evtch_001.PX(), pyi001 = Evtch_001.PY(), pzi001 = Evtch_001.PZ();
            
            double deltap02 = sqrt((px02-pxi02)*(px02-pxi02)+(py02-pyi02)*(py02-pyi02)+(pz02-pzi02)*(pz02-pzi02));
            double deltap01 = sqrt((px01-pxi01)*(px01-pxi01)+(py01-pyi01)*(py01-pyi01)+(pz01-pzi01)*(pz01-pzi01));
            
            b_tpl->column("deltap02", deltap02);
            b_tpl->column("deltap01", deltap01);


            //b_tpl->column("vx02", Evtch_02.VX()-EvtY4S.VX());
            //b_tpl->column("vy02", Evtch_02.VY()-EvtY4S.VY());
            //b_tpl->column("vz02", Evtch_02.VZ()-EvtY4S.VZ());
            //b_tpl->column("vx01", Evtch_01.VX()-EvtY4S.VX());
            //b_tpl->column("vy01", Evtch_01.VY()-EvtY4S.VY());
            //b_tpl->column("vz01", Evtch_01.VZ()-EvtY4S.VZ());
            //b_tpl->column("vx000", Evtch_000.VX()-EvtY4S.VX());
            //b_tpl->column("vy000", Evtch_000.VY()-EvtY4S.VY());
            //b_tpl->column("vz000", Evtch_000.VZ()-EvtY4S.VZ());
            //b_tpl->column("vx001", Evtch_001.VX()-EvtY4S.VX());
            //b_tpl->column("vy001", Evtch_001.VY()-EvtY4S.VY());
            //b_tpl->column("vz001", Evtch_001.VZ()-EvtY4S.VZ());
            //b_tpl->column("vx1", Evtch_1.VX()-EvtY4S.VX());
            //b_tpl->column("vy1", Evtch_1.VY()-EvtY4S.VY());
            //b_tpl->column("vz1", Evtch_1.VZ()-EvtY4S.VZ());
            

            //b_tpl->column("vxY", EvtY4S.VX());
            //b_tpl->column("vyY", EvtY4S.VY());
            //b_tpl->column("vzY", EvtY4S.VZ());

            if (Evtch_02.mo(0))
            {
                EvtP_02=gen_mgr[Evtch_02.mo(0)-1];
                b_tpl->column("nmo02", Evtch_02.mo(0));
                b_tpl->column("mo02", EvtP_02.idhep());
                if (EvtP_02.mo(0))
                {
                    EvtGP_02=gen_mgr[EvtP_02.mo(0)-1];
                    b_tpl->column("ngmo02", EvtP_02.mo(0));
                    b_tpl->column("gmo02", EvtGP_02.idhep());
                    //cout<<"gmo02="<<EvtGP_02.idhep()<<endl;
                    if (EvtGP_02.mo(0))
                	{
                	    EvtGGP_02=gen_mgr[EvtGP_02.mo(0)-1];
                	    b_tpl->column("ggmo02", EvtGGP_02.idhep());
                	    b_tpl->column("nggmo02", EvtGP_02.mo(0));
                	    //cout<<"ggmo02="<<EvtGP_02.idhep()<<endl;
                	}
                }
            }

            if (Evtch_01.mo(0))
            {
                EvtP_01=gen_mgr[Evtch_01.mo(0)-1];
                b_tpl->column("mo01", EvtP_01.idhep());
                b_tpl->column("nmo01", Evtch_01.mo(0));
                if (EvtP_01.mo(0))
                {
                    EvtGP_01=gen_mgr[EvtP_01.mo(0)-1];
                    b_tpl->column("gmo01", EvtGP_01.idhep());
                    b_tpl->column("ngmo01", EvtP_01.mo(0));
                    //cout<<"gmo01="<<EvtGP_01.idhep()<<endl;
                    if (EvtGP_01.mo(0))
                	{
                	    EvtGGP_01=gen_mgr[EvtGP_01.mo(0)-1];
                	    b_tpl->column("ggmo01", EvtGGP_01.idhep());
                	    b_tpl->column("nggmo01", EvtGP_01.mo(0));
                	    //cout<<"ggmo01="<<EvtGP_01.idhep()<<endl;
                	}
                }
            }

            if (Evtch_000.mo(0))
            {
                EvtP_000=gen_mgr[Evtch_000.mo(0)-1];
                b_tpl->column("mo000", EvtP_000.idhep());
                b_tpl->column("nmo000", Evtch_000.mo(0));
                if (EvtP_000.mo(0))
                {
                    EvtGP_000=gen_mgr[EvtP_000.mo(0)-1];
                    b_tpl->column("gmo000", EvtGP_000.idhep());
                    b_tpl->column("ngmo000", EvtP_000.mo(0));
                    if (EvtGP_000.mo(0)){
                        EvtGGP_000=gen_mgr[EvtGP_000.mo(0)-1];
                        b_tpl->column("ggmo000", EvtGGP_000.idhep());
                        b_tpl->column("nggmo000", EvtGP_000.mo(0));
                    }
                    //cout<<"gmo000="<<EvtGP_000.idhep()<<endl;
                }
            }

            if (Evtch_001.mo(0))
            {
                EvtP_001=gen_mgr[Evtch_001.mo(0)-1];
                b_tpl->column("mo001", EvtP_001.idhep());
                b_tpl->column("nmo001", Evtch_001.mo(0));
                if (EvtP_001.mo(0))
                {
                    EvtGP_001=gen_mgr[EvtP_001.mo(0)-1];
                    b_tpl->column("gmo001", EvtGP_001.idhep());
                    b_tpl->column("ngmo001", EvtP_001.mo(0));
                    if (EvtGP_001.mo(0)){
                        EvtGGP_001=gen_mgr[EvtGP_001.mo(0)-1];
                        b_tpl->column("ggmo001", EvtGGP_001.idhep());
                        b_tpl->column("nggmo001", EvtGP_001.mo(0));
                    }
                    //cout<<"gmo001="<<EvtGP_001.idhep()<<endl;
                }
            }

           //cout<<"Ready to test cosh of generator"<<endl;
            if(Evtch_000.mo(0)&&Evtch_01.mo(0)&&EvtP_000.mo(0)){
             //cout<<"Evtch_000.mo(0)&&Evtch_01.mo(0)&&EvtP_000.mo(0)"<<endl;
            	HepLorentzVector boost_x3872(EvtGP_000.PX(),EvtGP_000.PY(),EvtGP_000.PZ(),EvtGP_000.E());
            	
            	HepLorentzVector pijpsihep(EvtP_000.PX(),EvtP_000.PY(),EvtP_000.PZ(),EvtP_000.E());
            	HepLorentzVector pirho0hep(EvtP_01.PX(),EvtP_01.PY(),EvtP_01.PZ(),EvtP_01.E());
            	pijpsihep.boost( -boost_x3872.boostVector() );
            	pirho0hep.boost( -boost_x3872.boostVector() );
                boost_x3872.boost( boost_vector.boostVector() );
            	Vector3 pijpsi(pijpsihep.px(), pijpsihep.py(), pijpsihep.pz());
            	Vector3 pirho0(pirho0hep.px(), pirho0hep.py(), pirho0hep.pz());
            	Vector3 pix3872(boost_x3872.px(),boost_x3872.py(),boost_x3872.pz());
	
           //cout<<"cos_hel_jpsi"<<endl;
              double cos_hel_jpsi;
        	  cos_hel_jpsi = pijpsi.dot(pix3872);
        	  cos_hel_jpsi = cos_hel_jpsi/(pijpsi.mag()*pix3872.mag());
        	  b_tpl->column("coshij", cos_hel_jpsi);
        	  b_tpl->column("coshipxj", pijpsihep.px());
              b_tpl->column("coshipyj", pijpsihep.py());
              b_tpl->column("coshipzj", pijpsihep.pz());
        	  b_tpl->column("coshiej", pijpsihep.e());

           //cout<<"cos_hel_rho0i"<<endl;
        	  double cos_hel_rho0i;
        	  cos_hel_rho0i = pirho0.dot(pix3872);
        	  cos_hel_rho0i = cos_hel_rho0i/(pirho0.mag()*pix3872.mag());
        	  b_tpl->column("coshir", cos_hel_rho0i);
        	  b_tpl->column("coshipxr", pirho0hep.px());
              b_tpl->column("coshipyr", pirho0hep.py());
              b_tpl->column("coshipzr", pirho0hep.pz());
        	  b_tpl->column("coshier", pirho0hep.e());

           //cout<<"coshipxx"<<endl;
              b_tpl->column("coshipxx", boost_x3872.px());
              b_tpl->column("coshipyx", boost_x3872.py());
              b_tpl->column("coshipzx", boost_x3872.pz());
        	  }else{
             //cout<<"else-100"<<endl;
              b_tpl->column("coshij", -100);
              b_tpl->column("coshipxj", -100);
              b_tpl->column("coshipyj", -100);
              b_tpl->column("coshipzj", -100);
              b_tpl->column("coshiej", -100);
    
              b_tpl->column("coshir", -100);
              b_tpl->column("coshipxr", -100);
              b_tpl->column("coshipyr", -100);
              b_tpl->column("coshipzr", -100);
              b_tpl->column("coshier", -100);

              b_tpl->column("coshipxx",  -100);
              b_tpl->column("coshipyx",  -100);
              b_tpl->column("coshipzx",  -100);
            }
        //cout<<"Done test cosh of generator"<<endl;
            if (Evtch_1.mo(0))
            {
                EvtP_1=gen_mgr[Evtch_1.mo(0)-1];
                b_tpl->column("mo1", EvtP_1.idhep());
                b_tpl->column("nmo1", Evtch_1.mo(0));
                if (EvtP_1.mo(0))
                {
                    EvtGP_1=gen_mgr[EvtP_1.mo(0)-1];
                    b_tpl->column("gmo1", EvtGP_1.idhep());
                    b_tpl->column("ngmo1", EvtP_1.mo(0));
                    //cout<<"gmo1="<<EvtGP_1.idhep()<<endl;

                    if (EvtGP_1.mo(0))
                    {
                        //cout<<"ping20"<<endl;
                        EvtGGP_1=gen_mgr[EvtGP_1.mo(0)-1];
                        b_tpl->column("ggmo1", EvtGGP_1.idhep());
                        b_tpl->column("nggmo1", EvtGP_1.mo(0));
                        //cout<<"ping2"<<endl;
                        if (EvtGGP_1.mo(0))
                        {
                            //cout<<"ping30"<<endl;
                            EvtGGGP_1=gen_mgr[EvtGGP_1.mo(0)-1];
                            b_tpl->column("gggmo1", EvtGGGP_1.idhep());
                            b_tpl->column("ngggmo1", EvtGGP_1.mo(0));
                            //cout<<"ping3"<<endl;
                            if (EvtGGGP_1.mo(0))
                            {
                                //cout<<"ping40"<<endl;
                                EvtGGGGP_1=gen_mgr[EvtGGGP_1.mo(0)-1];
                                //b_tpl->column("ggggmo1", EvtGGGGP_1.idhep());
                                //b_tpl->column("nggggmo1", EvtGGGP_1.mo(0));
                                //cout<<"ping4"<<endl;
                                if (EvtGGGGP_1.mo(0))
                                {
                                    //cout<<"ping50"<<endl;
                                    EvtGGGGGP_1=gen_mgr[EvtGGGGP_1.mo(0)-1];
                                    //b_tpl->column("gggggmo1", EvtGGGGGP_1.idhep());
                                    //b_tpl->column("ngggggm1", EvtGGGGP_1.mo(0));
                                    //cout<<"ping5"<<endl;
                                }
                            }
                        }
                    }
                }
            }
    //cout<<"truth metric"<<endl;

         //cout <<"ch_000 = "<<Evtch_000.idhep()<<", ch_001 = "<<Evtch_001.idhep()<<endl;
         //cout <<"ch_1 = "<< Evtch_1.idhep()<<endl;
         //cout <<", P_000 = "<<EvtP_000.idhep()<<", P_001 = "<<EvtP_001.idhep()<<endl;
         //cout <<"GP_000 = "<<EvtGP_000.idhep()<<", GP_001 = "<<EvtGP_001.idhep();
         //cout<<", GGP_000 = "<<EvtGGP_000.idhep()<<", GGP_001 = "<<EvtGGP_001.idhep()<<endl;
         //cout <<"ch_01 = "<<Evtch_01.idhep()<<", ch_02 = "<<Evtch_02.idhep()<<", P_01 = "<<EvtP_01.idhep()<<", P_02 = "<<EvtP_02.idhep()<<endl;
         //cout <<"GP_01 = "<<EvtGP_01.idhep()<<", GP_02 = "<<EvtGP_02.idhep();
         //cout <<", GGP_01 = "<<EvtGGP_01.idhep()<<", GGP_02 = "<<EvtGGP_02.idhep()<<endl;
         //",  P_1 = "<<EvtP_1.idhep()<<endl;
         if(Evtch_000.mo(0)&&Evtch_001.mo(0)&&Evtch_01.mo(0)&&Evtch_02.mo(0)&&Evtch_1.mo(0)){
         	//cout<<"Evtch_000.mo(0)&&Evtch_001.mo(0)&&Evtch_01.mo(0)&&Evtch_02.mo(0)&&Evtch_1.mo(0)"<<endl;
           if (EvtP_000.mo(0)&&EvtP_001.mo(0)&&EvtP_01.mo(0)&&EvtP_02.mo(0)){
           	//cout<<"EvtP_000.mo(0)&&EvtP_001.mo(0)&&EvtP_01.mo(0)&&EvtP_02.mo(0)"<<endl;
           	if(EvtGP_01.mo(0)&&EvtGP_02.mo(0)&&EvtGP_000.mo(0)&&EvtGP_001.mo(0)){
           		//cout<<"EvtGP_01.mo(0)&&EvtGP_02.mo(0)&&EvtGP_000.mo(0)&&EvtGP_001.mo(0)"<<endl;

        //cout <<"ch_000 = "<<Evtch_000.idhep()<<", ch_001 = "<<Evtch_001.idhep()<<", P_000 = "<<EvtP_000.idhep()<<", P_001 = "<<EvtP_001.idhep()<<endl;
        //cout <<"GP_000 = "<<EvtGP_000.idhep()<<", GP_001 = "<<EvtGP_001.idhep();
        //cout<<", GGP_000 = "<<EvtGGP_000.idhep()<<", GGP_001 = "<<EvtGGP_001.idhep()<<endl;
        //cout <<"ch_01 = "<<Evtch_01.idhep()<<", ch_02 = "<<Evtch_02.idhep()<<", P_01 = "<<EvtP_01.idhep()<<", P_02 = "<<EvtP_02.idhep()<<endl;
        //cout <<"GP_01 = "<<EvtGP_01.idhep()<<", GP_02 = "<<EvtGP_02.idhep();
        //cout <<", GGP_01 = "<<EvtGGP_01.idhep()<<", GGP_02 = "<<EvtGGP_02.idhep()<<endl;
        //cout <<"ch_1 = "<< Evtch_1.idhep();
        //cout <<",  P_1 = "<<EvtP_1.idhep()<<endl;

            if(abs(Evtch_000.idhep()) == 13 && Evtch_001.idhep() == -Evtch_000.idhep() && abs(Evtch_01.idhep()) == 211 && Evtch_02.idhep() == -Evtch_01.idhep() 
                && EvtP_000.idhep() == 443 && EvtP_001.idhep() == 443 && EvtGP_01.idhep() == 120443 && EvtGP_02.idhep() == 120443 && Evtch_02.mo(0) == Evtch_01.mo(0) 
                && EvtGP_000.idhep() == 120443 && EvtGP_001.idhep() == 120443 && abs(EvtGGP_01.idhep()) == 511 && EvtGGP_02.idhep() == EvtGGP_01.idhep()
                && Evtch_1.idhep() == 22 && EvtP_1.idhep() == EvtGGP_01.idhep() && EvtGGP_000.idhep() == EvtGGP_01.idhep() && EvtGGP_001.idhep() == EvtGGP_01.idhep()
                && EvtP_01.idhep() == 113 && EvtP_02.idhep() == 113 )
            {
                hindex = 1;
            //cout<<"CASE1!!!"<<endl;
            }
            else
            {
                if (abs(Evtch_000.idhep()) == 13 && Evtch_001.idhep() == -Evtch_000.idhep() && abs(Evtch_01.idhep()) == 211 && Evtch_02.idhep() == -Evtch_01.idhep() 
                && EvtP_000.idhep() == 443 && EvtP_001.idhep() == 443 && EvtGP_01.idhep() == 120443 && EvtGP_02.idhep() == 120443 && Evtch_02.mo(0) == Evtch_01.mo(0) 
                && EvtGP_000.idhep() == 120443 && EvtGP_001.idhep() == 120443 && abs(EvtGGP_01.idhep()) == 511 && EvtGGP_02.idhep()==EvtGGP_01.idhep()
                && EvtGGP_001.idhep()==EvtGGP_000.idhep() && EvtGGP_001.idhep()==EvtGGP_01.idhep() && EvtP_01.idhep() == 113 && EvtP_02.idhep() == 113 )
                {
                //cout<<"if gamma false"<<endl;
                    if ((abs(Evtch_1.idhep())==11||Evtch_1.idhep()==22) && (abs(EvtP_1.idhep())==11||EvtP_1.idhep()==22))
                    {
                    //cout<<"if 22->11->22"<<endl;
                        if (EvtP_1.idhep() == 22 && EvtGP_1.idhep() == EvtGGP_01.idhep())
                        {
                        //cout<<"if"<<endl; hindex = 2;
                        }
                        else if (abs(EvtGP_1.idhep())==11||EvtGP_1.idhep()==22)
                        {
                        //cout<<"else if"<<endl;
                            if (EvtGP_1.idhep()==22 && EvtGGP_1.idhep() == EvtGGP_01.idhep()){
                                hindex = 3;
                            }else if (EvtGGP_1.idhep()==22 && EvtGGGP_1.idhep() == EvtGGP_01.idhep()){
                                hindex = 4;
                            }else if ((abs(EvtGGP_1.idhep())==11||EvtGGP_1.idhep()==22) && EvtGGGP_1.idhep()==22 && EvtGGGGP_1.idhep() == EvtGGP_01.idhep()){
                                hindex = 5;
                            }else if ((abs(EvtGGP_1.idhep())==11||EvtGGP_1.idhep()==22) && (abs(EvtGGGP_1.idhep())==11||EvtGGGP_1.idhep()==22) && EvtGGGGP_1.idhep()==22 && EvtGGGGGP_1.idhep() == EvtGGP_01.idhep()){
                                hindex = 6;
                            }else{
                                hindex = 0;
                            }
                        }
                        else
                        {
                        //cout<<"else"<<endl; 
                          hindex = 0;
                        }
                    //cout<<"ping02"<<endl;
                    }
                    else
                    {
                        hindex = 0;
                    //cout<<"ping00"<<endl;
                    }
                }else{
                    hindex = 0; 
                  //cout<<"ping000"<<endl;
                }
                
            //cout<<"NOCASE!!!"<<endl;
            }//End truth metric
          }else{hindex = 0;}//End EvtGP
          }else{hindex = 0;}//End EvtP
          }else{hindex = 0;}//End Evtch
            b_tpl->column("hindex", hindex);

            //if(hindex==0) continue;
//421: D0, 321: K+, 211:Pi+, 521:B+, 120443:X(3872)
            // Jpsi:443, e-:11, gamma: 22
            //(Evtch_00.idhep() == 321 (K+) && Evtch_01.idhep() == -211 (pi-) && EvtP_00.idhep() == -421 (D0Bar) && EvtP_01.idhep() == -421 (D0Bar) &&
            //Evtch_1.idhep() == 211 (PI+) && EvtP_1.idhep() == 521 (B+) && EvtGP_00.idhep() == 521 (B+) && EvtGP_01.idhep() == 521 (B+))
            //B->D0pi+ ->K+pi-
            //B->x3872gamma ->e+e-
//--------------------- trace from top to down (from B meson)------------------------------------
 //cout<<"Topdown"<<endl;
            int nbpd=0,nbnd=0,nxd=0,njd=0;
            int jindex = 0;
            for (std::vector<Gen_hepevt>::iterator i = gen_mgr.begin(); i != gen_mgr.end(); i++)
            {
                if((*i).idhep()==443) jindex=1;
                //b_tpl->column("bid",i);
                if ((*i).idhep()==521||(*i).idhep()== 511)
                {
                    b_tpl->column("np",1);
                    int da1=(*i).da(0),da2=(*i).da(1);

                    if ((*i).idhep() == 521)
                        {
                            b_tpl->column("charged",1);
                            //cout<<"charged521"<<endl;
                        }
                    else b_tpl->column("charged",0);
                    nbpd=(da2-da1+1);
                    b_tpl->column("nbpd",nbpd);

                    //cout<<"nbpd!!!"<<endl;
                    for(int start=0; start<(da2-da1+1); start++)
                    {
                        Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];
                        char bpdnumber[32];char bpdid[32];
                        sprintf (bpdnumber,"%s%d","bpd",start+1);
                        sprintf (bpdid,"%s%d","bpdid",start+1);
                        b_tpl->column(bpdnumber,Evda.idhep());
                        b_tpl->column(bpdid,(*i).da(0)+start);
                        //cout<<"bpd!!!"<<endl;
                    }

                    Gen_hepevt Evdax3872=gen_mgr[(*i).da(0)-1];
                    b_tpl->column("npx3872", (*i).da(0));
                    int dax1=Evdax3872.da(0),dax2=Evdax3872.da(1);
                    nxd=(dax2-dax1+1);
                    b_tpl->column("npxd",nxd);
                    //cout<<"nbpd!!!"<<endl;
                    if (Evdax3872.da(0)){
                        for(int start1=0; start1<(dax2-dax1+1); start1++)
                        {
                            Gen_hepevt Evdax=gen_mgr[Evdax3872.da(0)-1+start1];
                            char bxdnumber[32];char bxdid[32];
                            sprintf (bxdnumber,"%s%d","bpxd",start1+1);
                            sprintf (bxdid,"%s%d","bpxdid",start1+1);
                            b_tpl->column(bxdnumber,Evdax.idhep());
                            b_tpl->column(bxdid,Evdax3872.da(0)+start1);
                            //cout<<"bpd!!!"<<endl;
                        }
                    
                        Gen_hepevt Evdajpsi=gen_mgr[Evdax3872.da(0)-1];
                        b_tpl->column("npjpsi", Evdax3872.da(0));
                        int daj1=Evdajpsi.da(0),daj2=Evdajpsi.da(1);
                        njd=(daj2-daj1+1);
                        b_tpl->column("npjd",njd);
                        //cout<<"nbpd!!!"<<endl;
                        if (Evdajpsi.da(0)){
                            for(int start2=0; start2<(daj2-daj1+1); start2++)
                            {
                                Gen_hepevt Evdaj=gen_mgr[Evdajpsi.da(0)-1+start2];
                                char bjdnumber[32];char bjdid[32];
                                sprintf (bjdnumber,"%s%d","bpjd",start2+1);
                                sprintf (bjdid,"%s%d","bpjdid",start2+1);
                                b_tpl->column(bjdnumber,Evdaj.idhep());
                                b_tpl->column(bjdid,Evdajpsi.da(0)+start2);
                                //cout<<"bpd!!!"<<endl;
                            }
                        }
                    }
                    //Gen_hepevt Evdagamma=gen_mgr[(*i).da(0)];
                    //b_tpl->column("ngamma", (*i).da(0)+1);
                }
                else if((*i).idhep()== -521 || (*i).idhep()== -511)
                {
                    b_tpl->column("np",-1);
                    int da1=(*i).da(0),da2=(*i).da(1);

                    if ((*i).idhep() == -521)
                        {
                            b_tpl->column("charged",1);
                            //cout<<"charged-521"<<endl;
                        }
                    else b_tpl->column("charged",0);
                    nbnd=(da2-da1+1);
                    b_tpl->column("nbnd",nbnd);
                    //cout<<"nbnd!!!"<<endl;
                    for(int start=0; start<(da2-da1+1); start++)
                    {
                        Gen_hepevt Evda=gen_mgr[(*i).da(0)-1+start];

                        char bndnumber[32];char bndid[32];
                        sprintf (bndnumber,"%s%d","bnd",start+1);//bnd1,bnd2,bnd3....
                        sprintf (bndid,"%s%d","bndid",start+1);
                        b_tpl->column(bndnumber,Evda.idhep());
                        b_tpl->column(bndid,(*i).da(0)+start);
                        //cout<<"bnd!!!"<<endl;
                    }
                    Gen_hepevt Evdax3872=gen_mgr[(*i).da(0)-1];
                    b_tpl->column("nnx3872", (*i).da(0));
                    int dax1=Evdax3872.da(0),dax2=Evdax3872.da(1);
                    nxd=(dax2-dax1+1);
                    b_tpl->column("nnxd",nxd);
                    //cout<<"nbpd!!!"<<endl;
                    if (Evdax3872.da(0)){
                        for(int start1=0; start1<(dax2-dax1+1); start1++)
                        {
                            Gen_hepevt Evdax=gen_mgr[Evdax3872.da(0)-1+start1];
                            char bxdnumber[32];char bxdid[32];
                            sprintf (bxdnumber,"%s%d","bnxd",start1+1);
                            sprintf (bxdid,"%s%d","bnxdid",start1+1);
                            b_tpl->column(bxdnumber,Evdax.idhep());
                            b_tpl->column(bxdid,Evdax3872.da(0)+start1);
                            //cout<<"bpd!!!"<<endl;
                        }
                    
                        Gen_hepevt Evdajpsi=gen_mgr[Evdax3872.da(0)-1];
                        b_tpl->column("nnjpsi", Evdax3872.da(0));
                        int daj1=Evdajpsi.da(0),daj2=Evdajpsi.da(1);
                        njd=(daj2-daj1+1);
                        b_tpl->column("nnjd",njd);
                        //cout<<"nbpd!!!"<<endl;
                        if (Evdajpsi.da(0)){
                            for(int start2=0; start2<(daj2-daj1+1); start2++)
                            {
                                Gen_hepevt Evdaj=gen_mgr[Evdajpsi.da(0)-1+start2];
                                char bjdnumber[32];char bjdid[32];
                                sprintf (bjdnumber,"%s%d","bnjd",start2+1);
                                sprintf (bjdid,"%s%d","bnjdid",start2+1);
                                b_tpl->column(bjdnumber,Evdaj.idhep());
                                b_tpl->column(bjdid,Evdajpsi.da(0)+start2);
                                //cout<<"bpd!!!"<<endl;
                            }
                        }
                    }
                    //Gen_hepevt Evdagamma=gen_mgr[(*i).da(0)];
                    //b_tpl->column("ngamma", (*i).da(0)+1);
                }
            }
            b_tpl->column("jindex",jindex);
        }


//-----------      continuum suppresion  ------------------                                                                                                
 //cout<<"continuum suppresion..."<<std::endl;
      int vertexflag;                                                                                    
      HepPoint3D overtex;                                                                                
      
      Vector4 OtherBVector;                                                                                                  
      float R2, spher, cos_thr, cos_thp, par_sfw[13];                                                                        
      
      
      shape((*i), R2,  spher, cos_thr, cos_thp, par_sfw, OtherBVector, bvertex, vertexflag, overtex);                        
      //if (abs(cos_thr)<=1){}else{continue;}
      
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
      //b_tpl->column("costhp", float(cos_thp));         
      b_tpl->column("spher",  float(spher));   
      b_tpl->column("cosb", cosb);         
      b_tpl->column("vtflag",vertexflag);         
      if( !vertexflag)                                
      {                                                                 
        double deltaZ1=bvertex.z()-overtex.z();
        b_tpl->column("deltaz",deltaZ1);
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

    }//end of D0

//B_cand1 for k pi
/*
cout<<"============================="<<endl;
cout<<"Anglecase1:"<<endl;
cout<<var1[0]<<","<<var1[1]<<","<<var1[2]<<","<<var1[3]<<","<<var1[4]<<","<<endl;
cout<<"Anglecase2:"<<endl;
cout<<var2[0]<<","<<var2[1]<<","<<var2[2]<<","<<var2[3]<<","<<var2[4]<<","<<endl;
cout<<"Anglecase3:"<<endl;
cout<<var3[0]<<","<<var3[1]<<","<<var3[2]<<","<<var3[3]<<","<<var3[4]<<","<<endl;
cout<<"============================="<<endl;*/
}//end of fill pi or k lists from MDST_Charged Data Base

//void ana_d0pi::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13], Vector4 &otherB_P)
void ana_x3872gamma::shape(Particle &b, float &R2, float &spher, float &cos_thr, float &cos_thp, float par_sfw[13],Vector4 &otherB_P, HepPoint3D &bvertex, int & vertexflag, HepPoint3D & overtex)
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

    for( int i=0; i< b.relation().nFinalStateParticles();
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
    for(int i=0; i<b.relation().child(0).relation().nFinalStateParticles(); i++)
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
                j!= candiBfinalParticle.end(); j++)
            if(checkMultiUse(*it1,*j))
            {
                notBDauFlag=0;
                break;
            }
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
                j!= candiBfinalParticle.end(); j++)
            if(checkMultiUse(*it1,*j))
            {
                notBDauFlag=0;
                break;
            }
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
                j!= candiBfinalParticle.end(); j++)
            if(checkMultiUse(*it1,*j))
            {
                notBDauFlag=0;
                break;
            }
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
    /*
      if((candiBfinal.size()+otherBfinal.size()) != allFinal.size())
        //cout << candiBfinal.size()<<"+"<<otherBfinal.size()<<" !="<< allFinal.size()<<endl;
    */

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
void ana_x3872gamma::GetImpactParameters(const Mdst_charged *charged, double *dr, double *dz, int t)
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