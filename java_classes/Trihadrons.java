/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analyzers;

/**
 *
 * @author tbhayward
 */

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;


public class Trihadrons {
    
    protected byte helicity;
    protected int runnum;
    
    protected double test;
    
    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    
    // labels are unnumbered if they refer to the trihadron (perhaps a meson) and numbered for individual
    // hadrons. Convention is ordered by mass, then charge. For example in pi+pi- pi+ is hadron 1
    // in proton+pi+ the proton is p1, in k+pi- the kaon is p1.
    protected double Q2, W, gamma, nu, x, y, z, z1, z2, z3, z12, z13, z23;
    protected double Mx, Mx1, Mx2, Mx3, Mx12, Mx13, Mx23; // Mx is the Mx(ep1p2p3), Mx1 is Mx(e[p1]p2p3), etc.
    protected double Mh, Mh12, Mh13, Mh23;
    protected double pT, pT1, pT2, pT3, pT12, pT13, pT23;
    protected double xF, xF1, xF2, xF3, xF12, xF13, xF23;
    protected double zeta, zeta1, zeta2, zeta3, zeta12, zeta13, zeta23;
    protected double eta, eta1, eta2, eta3, eta12, eta13, eta23;
    protected double eta_gN, eta1_gN, eta2_gN, eta3_gN, eta12_gN, eta13_gN, eta23_gN;
    // eta is the rapidity, preferred by theorists in the Breit frame (e.g. eta1 is in Breit) 
    // eta_gN is the rapidity in the gamma*-nucleon COM frame
    // the difference between two rapidities is Lorentz invariant, i.e.
    // eta1-eta2 = eta1_COM - eta2_COM

    // theta is defined as the angle between the hadron p1 in the pair center-of-mass frame 
    // and the direction of the pair, Ph, in the photon-target rest frame.
    // phih and phiR are defined from the vectors Ph = p1+p2 and 2R = p1-p2, see 
    // hep-ex:2101.04842
    protected double theta, phih, phiR, phi1, phi2, phi3, phi12, phi13, phi23, Delta_phi12, Delta_phi13, Delta_phi23;
    
    // depolarization vectors defining the polarization lost during the transfer from beam to 
    // the virtual photon. 
    // in ALU BSAs the twist 2 terms are modified by C/A and the twist 3 terms by W/A
    // B and V come in AUL
    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;
    
    protected double e_px, e_py, e_pz, e_p, e_e, e_theta, e_phi; // electron kinematics
    protected double p_px, p_py, p_pz, p_p, p_e; // dihadron kinematics, mass is Mh
    protected double p1_px, p1_py, p1_pz, p1_p, p1_e, p1_theta, p1_phi; // p1 kinematics
    protected double p1_px_unc, p1_py_unc, p1_pz_unc, p1_p_unc, mass2, mass2_unc;
    protected double p2_px, p2_py, p2_pz, p2_p, p2_e, p2_theta, p2_phi; // p2 kinematics
    protected double p3_px, p3_py, p3_pz, p3_p, p3_e, p3_theta, p3_phi; // p3 kinematics
    protected double vx_e, vx_p1, vx_p2, vx_p3, vy_e, vy_p1, vy_p2, vy_p3, vz_e, vz_p1, vz_p2, vz_p3;
    
    protected double p_Breit_pz, p1_Breit_pz, p2_Breit_pz, p3_Breit_pz, p12_Breit_pz, p13_Breit_pz, p23_Breit_pz;
    protected double p_gN_pz, p1_gN_pz, p2_gN_pz, p3_gN_pz, p12_gN_pz, p13_gN_pz, p23_gN_pz;
    
    
    public static boolean channel_test(Trihadrons variables) {
//        if (variables.helicity==0){ return false; }
        if (variables.Q2()<1) { return false; } 
        else if (variables.W()<2) { return false; } 
        else if (variables.y()>0.75) { return false; } 
//        else if (variables.Mx()<1.5) { return false; } 
//        else if (variables.xF1()<0.0) { return false; } 
//        else if (variables.xF2()<0.0) { return false; } 
//        else if (variables.Mx1()<0.0) { return false; } 
//        else if (variables.Mx2()<0.0) { return false; } 
	return true;
    }
    
    public Trihadrons(DataEvent event, PhysicsEvent recEvent, int p1PID, int p1Index, int p2PID, int p2Index, 
            int p3PID, int p3Index) {
        // provide the PDG PID of the three hadrons
        
        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle"); // load particle bank
        
        helicity = eventBank.getByte("helicity", 0);
//        helicity = eventBank.getByte("helicityRaw", 0);
        runnum = configBank.getInt("run",0); // used for beam energy and polarization
        
    
        // default beam energy set to rga fall 2018
        double Eb = 10.6041; // rga fall 2018
        if (runnum >= 5032 && runnum <= 5666) { Eb = 10.6041; } // RGA Fall 2018
        else if (runnum >= 6616 && runnum <= 6783) { Eb = 10.1998; } // RGA Spring 2019
        else if (runnum >= 6120 && runnum <= 6399) { Eb = 10.5986; } // RGB Spring 2019
        else if (runnum >= 6409 && runnum <= 6604) { Eb = 10.1998; } // RGB Spring 2019
        else if (runnum >= 11093 && runnum <= 11283) { Eb = 10.4096; } // RGB Fall 2019
        else if (runnum >= 11284 && runnum <= 11300) { Eb = 4.17179; } // RGB Fall 2019
        else if (runnum >= 11323 && runnum <= 11571) { Eb = 10.3894; } // RGB Spring 2020
//        Eb = 22.00;
        
        num_elec = recEvent.countByPid(11); // returns number of electrons
	num_piplus = recEvent.countByPid(211); 
	num_piminus = recEvent.countByPid(-211);
	num_kplus = recEvent.countByPid(321);
	num_kminus = recEvent.countByPid(-321);
        num_protons = recEvent.countByPid(2212);
        num_particles = num_elec+num_piplus+num_piminus+num_kplus+num_kminus+num_protons;
        
        // Set up Lorentz vectors
        // beam electron
        LorentzVector lv_beam = new LorentzVector();
	lv_beam.setPxPyPzM(0, 0, Math.pow(Eb*Eb-particle_mass(11)*particle_mass(11),0.5), 
                particle_mass(11));
        LorentzVector lv_target = new LorentzVector();
        // target, proton for RGA... what mass to use for RGB (deuterium target)?
	lv_target.setPxPyPzM(0, 0, 0, particle_mass(2212));
        // pull from rec banks for outgoing particles
        // electron
        String electron_index = "[11,0]"; // highest p, kinematic fitter should require FD etc
	Particle scattered_electron = recEvent.getParticle(electron_index); //
        LorentzVector lv_e = new LorentzVector();
	lv_e.setPxPyPzM(scattered_electron.px(), scattered_electron.py(), 
            scattered_electron.pz(), particle_mass(11));
        // hadrons set up below (to allow for iteration over more than two hadrons in an event)
        
        // kinematics of electron
        e_px = lv_e.px(); e_py = lv_e.py(); e_pz = lv_e.pz(); e_p = lv_e.p(); e_e = lv_e.e(); 
        e_theta = scattered_electron.theta();
        e_phi = scattered_electron.phi();
        if (e_phi < 0) { e_phi = 2*Math.PI + e_phi; }
                
        // DIS variables
        LorentzVector lv_q = new LorentzVector(lv_beam); lv_q.sub(lv_e);
	Q2 = -lv_q.mass2();
	nu = lv_beam.e()-lv_e.e();
	x  = Q2 / (2 * particle_mass(2212) * nu);
	W  = Math.pow(Math.pow(particle_mass(2212),2)+2*particle_mass(2212)*nu - Q2, 0.5);
	y = nu/lv_beam.e();
        gamma = 2*particle_mass(2212)*x/Math.pow(Q2, 0.5);
        
        // Depolarization variables
        Depolarization_A = 1/(1+gamma*gamma)*(1-y+y*y/2+y*y*gamma*gamma/4);
        Depolarization_B = 1/(1+gamma*gamma)*(1-y-y*y*gamma*gamma/4);
        Depolarization_C = (y/Math.pow(1+gamma*gamma, 0.5))*(1-y/2);
        Depolarization_V = (2-y)/(1+gamma*gamma)*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
        Depolarization_W = y/(Math.pow(1+gamma*gamma, 0.5))*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
    
        
        // set up boost to gamma*-nucleon center of mass frame
        LorentzVector gN = new LorentzVector(lv_q);
	gN.add(lv_target);
	Vector3 gNBoost = gN.boostVector();
	gNBoost.negative();
        
        // set up boost to Breit frame, this needs to be cross checked
        LorentzVector Breit = new LorentzVector(lv_q);
        LorentzVector Breit_target = new LorentzVector();
        Breit_target.setPxPyPzM(0, 0, 0, 2*x*particle_mass(2212));
        Breit.add(Breit_target);
        Vector3 BreitBoost = Breit.boostVector();
        BreitBoost.negative();
        
        // set up hadrons and dihadron
        String p1Index_string = "["+p1PID+","+p1Index+"]";
        Particle p1 = recEvent.getParticle(p1Index_string);
        String p2Index_string = "["+p2PID+","+p2Index+"]";
        Particle p2 = recEvent.getParticle(p2Index_string);
        String p3Index_string = "["+p3PID+","+p3Index+"]";
        Particle p3 = recEvent.getParticle(p3Index_string);
        
        String combined_index_string = "["+p1PID+","+p1Index+"]+["+p2PID+","+p2Index+"]+["+p3PID+","+p3Index+"]";
        Particle trihadron = recEvent.getParticle(combined_index_string);
        Mh = trihadron.mass();
        combined_index_string = "["+p1PID+","+p1Index+"]+["+p2PID+","+p2Index+"]";
        Particle dihadron12 = recEvent.getParticle(combined_index_string);
        Mh12 = dihadron12.mass();
        combined_index_string = "["+p1PID+","+p1Index+"]+["+p3PID+","+p3Index+"]";
        Particle dihadron13 = recEvent.getParticle(combined_index_string);
        Mh13 = dihadron13.mass();
        combined_index_string = "["+p2PID+","+p2Index+"]+["+p3PID+","+p3Index+"]";
        Particle dihadron23 = recEvent.getParticle(combined_index_string);
        Mh23 = dihadron23.mass();
        
        vx_e = scattered_electron.vx();
        vx_p1 = p1.vx();
        vx_p2 = p2.vx();
        vx_p3 = p3.vx();
        vy_e = scattered_electron.vy();
        vy_p1 = p1.vy();
        vy_p2 = p2.vy();
        vy_p3 = p3.vy();
        vz_e = scattered_electron.vz();
        vz_p1 = p1.vz();
        vz_p2 = p2.vz();
        vz_p3 = p3.vz();
        
        LorentzVector lv_p = new LorentzVector();
        lv_p.setPxPyPzM(trihadron.px(), trihadron.py(), trihadron.pz(), trihadron.mass());
        LorentzVector lv_p1 = new LorentzVector();
        lv_p1.setPxPyPzM(p1.px(), p1.py(), p1.pz(), p1.mass());
        LorentzVector lv_p2 = new LorentzVector();
        lv_p2.setPxPyPzM(p2.px(), p2.py(), p2.pz(), p2.mass());
        LorentzVector lv_p3 = new LorentzVector();
        lv_p3.setPxPyPzM(p3.px(), p3.py(), p3.pz(), p3.mass());
        LorentzVector lv_p12 = new LorentzVector();
        lv_p12.setPxPyPzM(dihadron12.px(), dihadron12.py(), dihadron12.pz(), dihadron12.mass());
        LorentzVector lv_p13 = new LorentzVector();
        lv_p13.setPxPyPzM(dihadron13.px(), dihadron13.py(), dihadron13.pz(), dihadron13.mass());
        LorentzVector lv_p23 = new LorentzVector();
        lv_p23.setPxPyPzM(dihadron23.px(), dihadron23.py(), dihadron23.pz(), dihadron23.mass());
        
        
        // kinematics of hadrons
        p1_px = lv_p1.px(); p1_py = lv_p1.py(); p1_pz = lv_p1.pz(); p1_p = lv_p1.p(); p1_e = p1.e();
        p1_theta = p1.theta();
        p1_phi = p1.phi();
        if (p1_phi < 0) { p1_phi = 2*Math.PI + p1_phi; }
        
        p2_px = lv_p2.px(); p2_py = lv_p2.py(); p2_pz = lv_p2.pz(); p2_p = lv_p2.p(); p2_e = p2.e();
        p2_theta = p2.theta();
        p2_phi = p2.phi();
        if (p2_phi < 0) { p2_phi = 2*Math.PI + p2_phi; }
        
        p3_px = lv_p3.px(); p3_py = lv_p3.py(); p3_pz = lv_p3.pz(); p3_p = lv_p3.p(); p3_e = p3.e();
        p3_theta = p3.theta();
        p3_phi = p3.phi();
        if (p3_phi < 0) { p3_phi = 2*Math.PI + p3_phi; }
        
        p_px = lv_p.px(); p_py = lv_p.py(); p_pz = lv_p.pz(); p_p = lv_p.p(); p_e = lv_p.e();

        
        z = lv_p.e()/lv_q.e();
        z1 = lv_p1.e()/lv_q.e();
	z2 = lv_p2.e()/lv_q.e();
        z3 = lv_p3.e()/lv_q.e();
        z12 = lv_p12.e()/lv_q.e();
        z13 = lv_p13.e()/lv_q.e();
        z23 = lv_p23.e()/lv_q.e();
        
        // missing mass calculations
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target); lv_Mx.sub(lv_p1); lv_Mx.sub(lv_p2);
            lv_Mx.sub(lv_p3); Mx = lv_Mx.mass();  // missing mass with all observed
        LorentzVector lv_Mx1 = new LorentzVector(lv_q); lv_Mx1.add(lv_target); lv_Mx1.sub(lv_p1); 
        Mx1 = lv_Mx1.mass(); // missing mass with p1 observed
        LorentzVector lv_Mx2 = new LorentzVector(lv_q); lv_Mx2.add(lv_target); lv_Mx2.sub(lv_p2);
        Mx2 = lv_Mx2.mass(); // missing mass with p2 observed
        LorentzVector lv_Mx3 = new LorentzVector(lv_q); lv_Mx3.add(lv_target); lv_Mx3.sub(lv_p3);
        Mx3 = lv_Mx3.mass(); // missing mass with p3 observed
        LorentzVector lv_Mx12 = new LorentzVector(lv_q); lv_Mx12.add(lv_target); lv_Mx12.sub(lv_p1); lv_Mx12.sub(lv_p2); 
        Mx12 = lv_Mx12.mass(); // missing mass with p1 and p2 observed
        LorentzVector lv_Mx13 = new LorentzVector(lv_q); lv_Mx13.add(lv_target); lv_Mx13.sub(lv_p1); lv_Mx13.sub(lv_p3); 
        Mx13 = lv_Mx13.mass(); // missing mass with p1 and p3 observed
        LorentzVector lv_Mx23 = new LorentzVector(lv_q); lv_Mx23.add(lv_target); lv_Mx13.sub(lv_p2); lv_Mx13.sub(lv_p3); 
        Mx23 = lv_Mx23.mass(); // missing mass with p2 and p3 observed
        
        // boost to gamma*-nucleon center of mass frame
        LorentzVector lv_p_gN = new LorentzVector(lv_p); lv_p_gN.boost(gNBoost);
        LorentzVector lv_p1_gN = new LorentzVector(lv_p1); lv_p1_gN.boost(gNBoost);
        LorentzVector lv_p2_gN = new LorentzVector(lv_p2); lv_p2_gN.boost(gNBoost);
        LorentzVector lv_p3_gN = new LorentzVector(lv_p3); lv_p3_gN.boost(gNBoost);
        LorentzVector lv_p12_gN = new LorentzVector(lv_p12); lv_p12_gN.boost(gNBoost);
        LorentzVector lv_p13_gN = new LorentzVector(lv_p13); lv_p13_gN.boost(gNBoost);
        LorentzVector lv_p23_gN = new LorentzVector(lv_p23); lv_p23_gN.boost(gNBoost);
        LorentzVector lv_e_gN = new LorentzVector(lv_e); lv_e_gN.boost(gNBoost);
        LorentzVector lv_target_gN = new LorentzVector(lv_target); lv_target_gN.boost(gNBoost);
        LorentzVector lv_Mx_gN = new LorentzVector(lv_Mx); lv_Mx_gN.boost(gNBoost);
        Vector3 lv_e_gN_unit = new Vector3();
        lv_e_gN_unit.setMagThetaPhi(1, lv_e_gN.theta(), lv_e_gN.phi());
        LorentzVector lv_q_gN = new LorentzVector(lv_q); lv_q_gN.boost(gNBoost);
        Vector3 lv_q_gN_unit = new Vector3();
        lv_q_gN_unit.setMagThetaPhi(1, lv_q_gN.theta(), lv_q_gN.phi());
        Vector3 lv_y_gN = new Vector3(); // I guess this isn't really a lorenzvector so the lv name is wrong 
        lv_y_gN = lv_q_gN.vect().cross(lv_e_gN.vect()); 
        
        // boost to Breit infinite momentum frame
        LorentzVector lv_p_Breit = new LorentzVector(lv_p); lv_p_Breit.boost(BreitBoost);
        LorentzVector lv_p1_Breit = new LorentzVector(lv_p1); lv_p1_Breit.boost(BreitBoost);
        LorentzVector lv_p2_Breit = new LorentzVector(lv_p2); lv_p2_Breit.boost(BreitBoost);
        LorentzVector lv_p3_Breit = new LorentzVector(lv_p3); lv_p3_Breit.boost(BreitBoost);
        LorentzVector lv_p12_Breit = new LorentzVector(lv_p12); lv_p12_Breit.boost(BreitBoost);
        LorentzVector lv_p13_Breit = new LorentzVector(lv_p13); lv_p13_Breit.boost(BreitBoost);
        LorentzVector lv_p23_Breit = new LorentzVector(lv_p23); lv_p23_Breit.boost(BreitBoost);
        LorentzVector lv_e_Breit = new LorentzVector(lv_e); lv_e_Breit.boost(BreitBoost);
        Vector3 lv_e_Breit_unit = new Vector3();
        lv_e_Breit_unit.setMagThetaPhi(1, lv_e_Breit.theta(), lv_e_Breit.phi());
        LorentzVector lv_q_Breit = new LorentzVector(lv_q); lv_q_Breit.boost(BreitBoost);
        Vector3 lv_q_Breit_unit = new Vector3();
        lv_q_Breit_unit.setMagThetaPhi(1, lv_q_Breit.theta(), lv_q_Breit.phi());
        // note that in the Breit frame +z is antialigned to the direction of q
        
  
        // set up boost to dihadron rest frame (p2 and p3 center of mass frame)
        // just used to calculate theta
        Vector3 pCOMBoost = lv_p.boostVector(); pCOMBoost.negative();
        LorentzVector lv_p23_COM = new LorentzVector(lv_p); lv_p23_COM.boost(pCOMBoost);
        LorentzVector lv_p2_COM = new LorentzVector(lv_p2); lv_p2_COM.boost(pCOMBoost);
        LorentzVector lv_p3_COM = new LorentzVector(lv_p3); lv_p3_COM.boost(pCOMBoost);
//        theta = Math.acos(lv_p1_COM.vect().dot(lv_p_gN.vect()) / 
//                (lv_p1_COM.vect().mag()*lv_p_gN.vect().mag()));
        theta = Math.acos(lv_p2_COM.vect().dot(lv_p.vect()) / 
                (lv_p2_COM.vect().mag()*lv_p.vect().mag()));
  
        
        pT = lv_q_gN_unit.cross(lv_p_gN.vect()).mag();
        pT1 = lv_q_gN_unit.cross(lv_p1_gN.vect()).mag();
        pT2 = lv_q_gN_unit.cross(lv_p2_gN.vect()).mag();
        pT3 = lv_q_gN_unit.cross(lv_p3_gN.vect()).mag();
        pT12 = lv_q_gN_unit.cross(lv_p12_gN.vect()).mag();
        pT13 = lv_q_gN_unit.cross(lv_p13_gN.vect()).mag();
        pT23 = lv_q_gN_unit.cross(lv_p23_gN.vect()).mag();
    
        xF =  2*(lv_p_gN.vect() .dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF1 = 2*(lv_p1_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF2 = 2*(lv_p2_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF3 = 2*(lv_p3_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF12 = 2*(lv_p12_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF13 = 2*(lv_p13_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF23 = 2*(lv_p23_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        
        zeta = lv_p_gN.e()/lv_target_gN.e(); 
        zeta1 = lv_p1_gN.e()/lv_target_gN.e(); 
        zeta2 = lv_p2_gN.e()/lv_target_gN.e(); 
        zeta3 = lv_p3_gN.e()/lv_target_gN.e(); 
        zeta12 = lv_p12_gN.e()/lv_target_gN.e(); 
        zeta13 = lv_p13_gN.e()/lv_target_gN.e(); 
        zeta23 = lv_p23_gN.e()/lv_target_gN.e(); 

        
        p_gN_pz = lv_p_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p1_gN_pz = lv_p1_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p2_gN_pz = lv_p2_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p3_gN_pz = lv_p3_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p12_gN_pz = lv_p12_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p13_gN_pz = lv_p13_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p23_gN_pz = lv_p23_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p_Breit_pz = lv_p_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p1_Breit_pz = lv_p1_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p2_Breit_pz = lv_p2_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p3_Breit_pz = lv_p3_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p12_Breit_pz = lv_p12_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p13_Breit_pz = lv_p13_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p23_Breit_pz = lv_p23_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        
        // Breit frame rapidity
        eta = -0.5*Math.log((lv_p_Breit.e()+p_Breit_pz) /  (lv_p_Breit.e()-p_Breit_pz));
        eta1 = -0.5*Math.log((lv_p1_Breit.e()+p1_Breit_pz) /  (lv_p1_Breit.e()-p1_Breit_pz));
        eta2 = -0.5*Math.log((lv_p2_Breit.e()+p2_Breit_pz) /   (lv_p2_Breit.e()-p2_Breit_pz));
        eta3 = -0.5*Math.log((lv_p3_Breit.e()+p3_Breit_pz) /  (lv_p3_Breit.e()-p3_Breit_pz));
        eta12 = -0.5*Math.log((lv_p12_Breit.e()+p12_Breit_pz) /  (lv_p12_Breit.e()-p12_Breit_pz));
        eta13 = -0.5*Math.log((lv_p13_Breit.e()+p13_Breit_pz) /  (lv_p13_Breit.e()-p13_Breit_pz));
        eta23 = -0.5*Math.log((lv_p23_Breit.e()+p23_Breit_pz) /  (lv_p23_Breit.e()-p23_Breit_pz));
        
        // gamma*-nucleon frame rapidity
        eta_gN = 0.5*Math.log((lv_p_gN.e()+p_gN_pz) / (lv_p_gN.e()-p_gN_pz));
        eta1_gN = 0.5*Math.log((lv_p1_gN.e()+p1_gN_pz) /  (lv_p1_gN.e()-p1_gN_pz));
        eta2_gN = 0.5*Math.log((lv_p2_gN.e()+p2_gN_pz) /  (lv_p2_gN.e()-p2_gN_pz));
        eta3_gN = 0.5*Math.log((lv_p3_gN.e()+p3_gN_pz) /  (lv_p3_gN.e()-p3_gN_pz));
        eta12_gN = 0.5*Math.log((lv_p12_gN.e()+p12_gN_pz) /  (lv_p12_gN.e()-p12_gN_pz));
        eta13_gN = 0.5*Math.log((lv_p13_gN.e()+p13_gN_pz) /  (lv_p13_gN.e()-p13_gN_pz));
        eta23_gN = 0.5*Math.log((lv_p23_gN.e()+p23_gN_pz) /  (lv_p23_gN.e()-p23_gN_pz));
        
        Vector3 vecH = new Vector3();
        vecH.setMagThetaPhi(lv_p3_gN.vect().mag()/z3, lv_p3_gN.vect().theta(), lv_p3_gN.vect().phi());
        Vector3 vecR = new Vector3(lv_p2_gN.vect()); // not R yet
        vecR.setMagThetaPhi(lv_p2_gN.vect().mag()/z2, lv_p2_gN.vect().theta(), lv_p2_gN.vect().phi());
        vecR.sub(vecH); // this is really R now that the subtraction is done
        
	Vector3 vectRt = new Vector3();
	Vector3 R_Q = new Vector3();
							
	R_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());
	vectRt=vecR;
	vectRt.sub(R_Q);

	Vector3 vectPh = new Vector3(lv_p_gN.vect());
        Vector3 vectPh1 = new Vector3(lv_p1_gN.vect());
        Vector3 vectPh2 = new Vector3(lv_p2_gN.vect());
        Vector3 vectPh3 = new Vector3(lv_p3_gN.vect());
        Vector3 vectPh12 = new Vector3(lv_p12_gN.vect());
        Vector3 vectPh13 = new Vector3(lv_p13_gN.vect());
        Vector3 vectPh23 = new Vector3(lv_p23_gN.vect());
	Vector3 Pt_Q = new Vector3();
	Pt_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());
	Vector3 vectPhT = new Vector3(vectPh);
	vectPhT.sub(Pt_Q);
        Vector3 vectPhT1 = new Vector3(vectPh1);
	vectPhT1.sub(Pt_Q);
        Vector3 vectPhT2 = new Vector3(vectPh2);
	vectPhT2.sub(Pt_Q);
        Vector3 vectPhT3 = new Vector3(vectPh3);
	vectPhT3.sub(Pt_Q);
        Vector3 vectPhT12 = new Vector3(vectPh12);
	vectPhT12.sub(Pt_Q);
        Vector3 vectPhT13 = new Vector3(vectPh13);
	vectPhT13.sub(Pt_Q);
        Vector3 vectPhT23 = new Vector3(vectPh23);
	vectPhT23.sub(Pt_Q);
        
        Vector3 vT = new Vector3(lv_q_gN_unit.cross(lv_e_gN_unit)); vT.unit();
	Vector3 vTR = new Vector3(lv_q_gN_unit.cross(vectRt)); vTR.unit();
	Vector3 vTH = new Vector3(lv_q_gN_unit.cross(vectPhT)); vTH.unit();
        Vector3 vTH1 = new Vector3(lv_q_gN_unit.cross(vectPhT1)); vTH1.unit();
        Vector3 vTH2 = new Vector3(lv_q_gN_unit.cross(vectPhT2)); vTH2.unit();
        Vector3 vTH3 = new Vector3(lv_q_gN_unit.cross(vectPhT3)); vTH3.unit();
        Vector3 vTH12 = new Vector3(lv_q_gN_unit.cross(vectPhT12)); vTH12.unit();
        Vector3 vTH13 = new Vector3(lv_q_gN_unit.cross(vectPhT13)); vTH13.unit();
        Vector3 vTH23 = new Vector3(lv_q_gN_unit.cross(vectPhT23)); vTH23.unit();
        
        double cosPhiR = vT.dot(vTR);
	double sinPhiR = lv_e_gN.vect().cross(vectRt).dot(lv_q_gN_unit);
	double cosPhiH = vT.dot(vTH);
	double sinPhiH = lv_e_gN.vect().cross(vectPhT).dot(lv_q_gN_unit);
        double cosPhiH1 = vT.dot(vTH1);
	double sinPhiH1 = lv_e_gN.vect().cross(vectPhT1).dot(lv_q_gN_unit);
        double cosPhiH2 = vT.dot(vTH2);
	double sinPhiH2 = lv_e_gN.vect().cross(vectPhT2).dot(lv_q_gN_unit);
        double cosPhiH3 = vT.dot(vTH3);
	double sinPhiH3 = lv_e_gN.vect().cross(vectPhT3).dot(lv_q_gN_unit);
        double cosPhiH12 = vT.dot(vTH12);
	double sinPhiH12 = lv_e_gN.vect().cross(vectPhT12).dot(lv_q_gN_unit);
        double cosPhiH13 = vT.dot(vTH13);
	double sinPhiH13 = lv_e_gN.vect().cross(vectPhT13).dot(lv_q_gN_unit);
        double cosPhiH23 = vT.dot(vTH23);
	double sinPhiH23 = lv_e_gN.vect().cross(vectPhT23).dot(lv_q_gN_unit);

        // scaling
	double rScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectRt).mag();
	sinPhiR = sinPhiR/rScale;
	double hScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh).mag();
	sinPhiH = sinPhiH/hScale;
        double h1Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh1).mag();
	sinPhiH1 = sinPhiH1/h1Scale;
        double h2Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh2).mag();
	sinPhiH2 = sinPhiH2/h2Scale;
        double h3Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh3).mag();
	sinPhiH3 = sinPhiH3/h3Scale;
        double h12Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh12).mag();
	sinPhiH12 = sinPhiH12/h12Scale;
        double h13Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh13).mag();
	sinPhiH13 = sinPhiH13/h13Scale;
        double h23Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh23).mag();
	sinPhiH23 = sinPhiH23/h23Scale;
        
        phih = Math.acos(cosPhiH); 
        phi1 = Math.acos(cosPhiH1); 
        phi2 = Math.acos(cosPhiH2);
        phi3 = Math.acos(cosPhiH3);
        phi12 = Math.acos(cosPhiH12);
        phi13 = Math.acos(cosPhiH13);
        phi23 = Math.acos(cosPhiH23);
	phiR = Math.acos(cosPhiR); 
        
        if (sinPhiR < 0.0) {  phiR = 2*Math.PI - phiR; }
	if (sinPhiH < 0.0) { phih = 2*Math.PI - phih; }
        if (sinPhiH1 < 0.0) { phi1 = 2*Math.PI - phi1; }
        if (sinPhiH2 < 0.0) { phi2 = 2*Math.PI - phi2; }
        if (sinPhiH3 < 0.0) { phi3 = 2*Math.PI - phi3; }
        if (sinPhiH12 < 0.0) { phi12 = 2*Math.PI - phi12; }
        if (sinPhiH13 < 0.0) { phi13 = 2*Math.PI - phi13; }
        if (sinPhiH23 < 0.0) { phi23 = 2*Math.PI - phi23; }
        
        Delta_phi12 = phi2 - phi1;
        if (Delta_phi12 < 0) { Delta_phi12+=2*Math.PI; }
        Delta_phi13 = phi3 - phi1;
        if (Delta_phi13 < 0) { Delta_phi13+=2*Math.PI; }
        Delta_phi23 = phi3 - phi2;
        if (Delta_phi23 < 0) { Delta_phi23+=2*Math.PI; }
         
    }
   
    
    public int get_helicity() { // -1, 0, or 1. 0 equals unassigned by EventBuilder
        if (runnum <= 5666) {
            return -1*helicity;
        } else if ( runnum >= 6616 && runnum <= 6783) {
            return -1*helicity;
        } else if ( runnum >= 6120 && runnum <= 6604) { 
            return -1*helicity;
        } else if ( runnum >= 11093 && runnum <= 11283) {
            return helicity;
        } else if ( runnum >= 11284 && runnum < 11300) {
            return -1*helicity;
        } else if ( runnum >= 11323 && runnum <= 11571) {
            return helicity;
        }
        System.out.println("runnum not found, assigning helicity flip");
        return -1*helicity;
    }
    
    public int get_runnum() { return runnum; }; // returns run number for polarizations and energy
    
    public int num_elec() { return num_elec; } // returns number of electrons
    
    public int num_piplus() { return num_piplus; } // returns number of piplus
    
    public int num_piminus() { return num_piminus; } // returns number of piminus
    
    public int num_kplus() { return num_kplus; }// returns number of kplus
    
    public int num_kminus() { return num_kminus; } // returns number of kminus
    
    public int num_protons() { return num_protons; } // returns number of protons
    
    public double test() { return test; } // returns test var

    public double Q2() { return Double.valueOf(Math.round(Q2*100000))/100000; } // returns Q2
    
    public double W() { return Double.valueOf(Math.round(W*100000))/100000; }// returns W
    
    public double gamma() { return Double.valueOf(Math.round(gamma*100000))/100000; } // returns gamma
    
    public double nu() { return Double.valueOf(Math.round(nu*100000))/100000; }// returns nu
    
    public double x() { return Double.valueOf(Math.round(x*100000))/100000; }// returns x
    
    public double y() { return Double.valueOf(Math.round(y*100000))/100000; }// returns y
    
    public double z() { return Double.valueOf(Math.round(z*100000))/100000; }// returns z
    
    public double z1() { return Double.valueOf(Math.round(z1*100000))/100000; }// returns z1
    
    public double z2() { return Double.valueOf(Math.round(z2*100000))/100000; }// returns z2
    
    public double z3() { return Double.valueOf(Math.round(z3*100000))/100000; }// returns z3
    
    public double z12() { return Double.valueOf(Math.round(z12*100000))/100000; }// returns z12
    
    public double z13() { return Double.valueOf(Math.round(z13*100000))/100000; }// returns z13
    
    public double z23() { return Double.valueOf(Math.round(z23*100000))/100000; }// returns z23

    public double Mx() { return Double.valueOf(Math.round(Mx*100000))/100000; }
    
    public double Mx1() { return Double.valueOf(Math.round(Mx1*100000))/100000; }
    
    public double Mx2() { return Double.valueOf(Math.round(Mx2*100000))/100000; }
    
    public double Mx3() { return Double.valueOf(Math.round(Mx3*100000))/100000; }
    
    public double Mx12() { return Double.valueOf(Math.round(Mx12*100000))/100000; }
    
    public double Mx13() { return Double.valueOf(Math.round(Mx13*100000))/100000; }
    
    public double Mx23() { return Double.valueOf(Math.round(Mx23*100000))/100000; }
    
    public double Mh() { return Double.valueOf(Math.round(Mh*100000))/100000; }// returns Mh
    
    public double Mh12() { return Double.valueOf(Math.round(Mh12*100000))/100000; }// returns Mh12
    
    public double Mh13() { return Double.valueOf(Math.round(Mh13*100000))/100000; }// returns Mh13
    
    public double Mh23() { return Double.valueOf(Math.round(Mh23*100000))/100000; }// returns Mh23
    
    public double pT() { return Double.valueOf(Math.round(pT*100000))/100000; }// returns pT
    
    public double pT1() { return Double.valueOf(Math.round(pT1*100000))/100000; }// returns pT1
    
    public double pT2() { return Double.valueOf(Math.round(pT2*100000))/100000; }// returns pT2
    
    public double pT3() { return Double.valueOf(Math.round(pT3*100000))/100000; }// returns pT3
    
    public double pT12() { return Double.valueOf(Math.round(pT12*100000))/100000; }// returns pT12
    
    public double pT13() { return Double.valueOf(Math.round(pT13*100000))/100000; }// returns pT13
    
    public double pT23() { return Double.valueOf(Math.round(pT23*100000))/100000; }// returns pT23
    
    public double pTpT() { return Double.valueOf(Math.round(pT1*pT2*100000))/100000; } // returns ptpt
    
    public double xF() { return Double.valueOf(Math.round(xF*100000))/100000; }// returns xF
    
    public double xF1() { return Double.valueOf(Math.round(xF1*100000))/100000; }// returns xF1
    
    public double xF2() { return Double.valueOf(Math.round(xF2*100000))/100000; }// returns xF2
    
    public double xF3() { return Double.valueOf(Math.round(xF3*100000))/100000; }// returns xF3
    
    public double xF12() { return Double.valueOf(Math.round(xF12*100000))/100000; }// returns xF12
    
    public double xF13() { return Double.valueOf(Math.round(xF13*100000))/100000; }// returns xF13
    
    public double xF23() { return Double.valueOf(Math.round(xF23*100000))/100000; }// returns xF23
    
    public double zeta() { return Double.valueOf(Math.round(zeta*100000))/100000; }// returns zeta1
    
    public double zeta1() { return Double.valueOf(Math.round(zeta1*100000))/100000; }// returns zeta1
    
    public double zeta2() { return Double.valueOf(Math.round(zeta2*100000))/100000; }// returns zeta2
    
    public double zeta3() { return Double.valueOf(Math.round(zeta3*100000))/100000; }// returns zeta3
    
    public double zeta12() { return Double.valueOf(Math.round(zeta12*100000))/100000; }// returns zeta12
    
    public double zeta13() { return Double.valueOf(Math.round(zeta13*100000))/100000; }// returns zeta13
    
    public double zeta23() { return Double.valueOf(Math.round(zeta23*100000))/100000; }// returns zeta23
    
    public double p1_Breit_pz() { return Double.valueOf(Math.round(p1_Breit_pz*100000))/100000; } 
    // returns p1 pz in Breit
    
    public double p2_Breit_pz() { return Double.valueOf(Math.round(p2_Breit_pz*100000))/100000; } 
    // returns p2 pz in Breit
    
    public double p3_Breit_pz() { return Double.valueOf(Math.round(p3_Breit_pz*100000))/100000; } 
    // returns p3 pz in Breit
    
    public double p12_Breit_pz() { return Double.valueOf(Math.round(p12_Breit_pz*100000))/100000; } 
    // returns p12 pz in Breit
    
    public double p13_Breit_pz() { return Double.valueOf(Math.round(p13_Breit_pz*100000))/100000; } 
    // returns p13 pz in Breit
    
    public double p23_Breit_pz() { return Double.valueOf(Math.round(p23_Breit_pz*100000))/100000; } 
    // returns p23 pz in Breit
    
    public double p1_gN_pz() { return Double.valueOf(Math.round(p1_gN_pz*100000))/100000; } // returns p1 pz in gN
    
    public double p2_gN_pz() { return Double.valueOf(Math.round(p2_gN_pz*100000))/100000; } // returns p2 pz in gN
    
    public double p3_gN_pz() { return Double.valueOf(Math.round(p3_gN_pz*100000))/100000; } // returns p3 pz in gN
    
    public double p12_gN_pz() { return Double.valueOf(Math.round(p12_gN_pz*100000))/100000; } // returns p12 pz in gN
    
    public double p13_gN_pz() { return Double.valueOf(Math.round(p13_gN_pz*100000))/100000; } // returns p13 pz in gN
    
    public double p23_gN_pz() { return Double.valueOf(Math.round(p23_gN_pz*100000))/100000; } // returns p23 pz in gN
    
    public double eta() { return Double.valueOf(Math.round(eta*100000))/100000; }// returns eta in the Breit frame
    
    public double eta1() { return Double.valueOf(Math.round(eta1*100000))/100000; }// returns eta1 in the Breit frame
    
    public double eta2() { return Double.valueOf(Math.round(eta2*100000))/100000; }// returns eta2 in the Breit frame
    
    public double eta3() { return Double.valueOf(Math.round(eta3*100000))/100000; }// returns eta3 in the Breit frame
    
    public double eta12() { return Double.valueOf(Math.round(eta12*100000))/100000; }// returns eta12 in the Breit frame
    
    public double eta13() { return Double.valueOf(Math.round(eta13*100000))/100000; }// returns eta13 in the Breit frame
    
    public double eta23() { return Double.valueOf(Math.round(eta23*100000))/100000; }// returns eta23 in the Breit frame
    
    public double eta_gN() { return Double.valueOf(Math.round(eta_gN*100000))/100000; }// returns eta_gN
    
    public double eta1_gN() { return Double.valueOf(Math.round(eta1_gN*100000))/100000; }// returns eta1_gN
    
    public double eta2_gN() { return Double.valueOf(Math.round(eta2_gN*100000))/100000; }// returns eta2_gN
    
    public double eta3_gN() { return Double.valueOf(Math.round(eta3_gN*100000))/100000; }// returns eta3_gN
    
    public double Delta_eta12() { return Double.valueOf(Math.round((eta2-eta1)*100000))/100000; }// returns Delta_eta12,
    
    public double Delta_eta13() { return Double.valueOf(Math.round((eta3-eta1)*100000))/100000; }// returns Delta_eta13,
    
    public double Delta_eta23() { return Double.valueOf(Math.round((eta3-eta2)*100000))/100000; }// returns Delta_eta32,
    // should be Lorentz invariant
    // maybe set a check for difference of Delta_eta in Breit and gN frames in future...
    // make sure the sign is what you really want
    
    public double theta() { return Double.valueOf(Math.round(theta*100000))/100000; }
    // returns theta, the dihadron "decay angle"
    
    public double phih() { return Double.valueOf(Math.round(phih*100000))/100000; }// returns phih
    
    public double phiR() { return Double.valueOf(Math.round(phiR*100000))/100000; }// returns phiR
    
    public double phi1() { return Double.valueOf(Math.round(phi1*100000))/100000; }
    // returns phi1 (gamma*-nucleon frame)
    
    public double phi2() { return Double.valueOf(Math.round(phi2*100000))/100000; }
    // returns phi2 (gamma*-nucleon frame)
    
    public double phi3() { return Double.valueOf(Math.round(phi3*100000))/100000; }
    // returns phi3 (gamma*-nucleon frame)
    
    public double phi12() { return Double.valueOf(Math.round(phi12*100000))/100000; }
    // returns phi12 (gamma*-nucleon frame)
    
    public double phi13() { return Double.valueOf(Math.round(phi13*100000))/100000; }
    // returns phi13 (gamma*-nucleon frame)
    
    public double phi23() { return Double.valueOf(Math.round(phi23*100000))/100000; }
    // returns phi23 (gamma*-nucleon frame)
    
    public double Delta_phi12() {return Double.valueOf(Math.round(Delta_phi12*100000))/100000; }//returns Delta_phi (p2-p1)
    
    public double Delta_phi13() {return Double.valueOf(Math.round(Delta_phi13*100000))/100000; }//returns Delta_phi (p3-p1)
    
    public double Delta_phi23() {return Double.valueOf(Math.round(Delta_phi23*100000))/100000; }//returns Delta_phi (p3-p2)
     
    public double Depolarization_A() { return Double.valueOf(Math.round(Depolarization_A*100000))/100000; }
    // returns Depolarization_A
    
    public double Depolarization_B() { return Double.valueOf(Math.round(Depolarization_B*100000))/100000; }
    // returns Depolarization_B
    
    public double Depolarization_C() { return Double.valueOf(Math.round(Depolarization_C*100000))/100000; }
    // returns Depolarization_C
    
    public double Depolarization_V() { return Double.valueOf(Math.round(Depolarization_V*100000))/100000; }
    // returns Depolarization_V
    
    public double Depolarization_W() { return Double.valueOf(Math.round(Depolarization_W*100000))/100000; }
    // returns Depolarization_W
    
    public double e_px() { return Double.valueOf(Math.round(e_px*100000))/100000; }// returns electron lab frame px
    
    public double e_py() { return Double.valueOf(Math.round(e_py*100000))/100000; }// returns electron lab frame py
    
    public double e_pz() { return Double.valueOf(Math.round(e_pz*100000))/100000; }// returns electron lab frame pz
    
    public double e_p() { return Double.valueOf(Math.round(e_p*100000))/100000; }// returns electron lab frame p
    
    public double e_e() { return Double.valueOf(Math.round(e_e*100000))/100000; }// returns electron lab frame energy
    
    public double e_theta() { return Double.valueOf(Math.round(e_theta*100000))/100000; } // returns electron lab 
    // frame polar angle
    
    public double e_phi() { return Double.valueOf(Math.round(e_phi*100000))/100000; } // returns electron lab 
    // frame polar angle
    
    public double p_px() { return Double.valueOf(Math.round(p_px*100000))/100000; }// returns dihadron lab frame px
    
    public double p_py() { return Double.valueOf(Math.round(p_py*100000))/100000; }// returns dihadron lab frame py
    
    public double p_pz() { return Double.valueOf(Math.round(p_pz*100000))/100000; }// returns dihadron lab frame pz
    
    public double p_p() { return Double.valueOf(Math.round(p_p*100000))/100000; }// returns dihadron lab frame p
    
    public double p_e() { return Double.valueOf(Math.round(p_e*100000))/100000; }// returns dihadron lab frame energy
    
    public double p1_px() { return Double.valueOf(Math.round(p1_px*100000))/100000; }// returns hadron 1 lab frame px
    
    public double p1_py() { return Double.valueOf(Math.round(p1_py*100000))/100000; }// returns hadron 1 lab frame py
    
    public double p1_pz() { return Double.valueOf(Math.round(p1_pz*100000))/100000; }// returns hadron 1 lab frame pz
    
    public double p1_p() { return Double.valueOf(Math.round(p1_p*100000))/100000; }// returns hadron 1 lab frame p
    
    public double p1_e() { return Double.valueOf(Math.round(p1_e*100000))/100000; }// returns hadron 1 lab frame energy
    
    public double p1_theta() { return Double.valueOf(Math.round(p1_theta*100000))/100000; } // returns p1 lab 
    // frame polar angle
    
    public double p1_phi() { return Double.valueOf(Math.round(p1_phi*100000))/100000; } // returns p1 lab 
    // frame polar angle
    
    public double p2_px() { return Double.valueOf(Math.round(p2_px*100000))/100000; }// returns hadron 2 lab frame px
    
    public double p2_py() { return Double.valueOf(Math.round(p2_py*100000))/100000; }// returns hadron 2 lab frame py
    
    public double p2_pz() { return Double.valueOf(Math.round(p2_pz*100000))/100000; }// returns hadron 2 lab frame pz
    
    public double p2_p() { return Double.valueOf(Math.round(p2_p*100000))/100000; }// returns hadron 2 lab frame p
    
    public double p2_e() { return Double.valueOf(Math.round(p2_e*100000))/100000; }// returns hadron 2 lab frame energy
    
    public double p2_theta() { return Double.valueOf(Math.round(p2_theta*100000))/100000; } // returns p2 lab frame polar angle
    
    public double p2_phi() { return Double.valueOf(Math.round(p2_phi*100000))/100000; } // returns p2 lab frame polar angle
    
    public double p3_px() { return Double.valueOf(Math.round(p3_px*100000))/100000; }// returns hadron 3 lab frame px
    
    public double p3_py() { return Double.valueOf(Math.round(p3_py*100000))/100000; }// returns hadron 3 lab frame py
    
    public double p3_pz() { return Double.valueOf(Math.round(p3_pz*100000))/100000; }// returns hadron 3 lab frame pz
    
    public double p3_p() { return Double.valueOf(Math.round(p3_p*100000))/100000; }// returns hadron 3 lab frame p
    
    public double p3_e() { return Double.valueOf(Math.round(p3_e*100000))/100000; }// returns hadron 3 lab frame energy
    
    public double p3_theta() { return Double.valueOf(Math.round(p3_theta*100000))/100000; } // returns p3 lab frame polar angle
    
    public double p3_phi() { return Double.valueOf(Math.round(p3_phi*100000))/100000; } // returns p3 lab frame polar angle
    
    public double vx_e() { return Double.valueOf(Math.round(vx_e*100000))/100000; }// returns electron x vertex
    
    public double vx_p1() { return Double.valueOf(Math.round(vx_p1*100000))/100000; }// returns p1 x vertex
    
    public double vx_p2() { return Double.valueOf(Math.round(vx_p2*100000))/100000; }// returns p2 x vertex
    
    public double vx_p3() { return Double.valueOf(Math.round(vx_p3*100000))/100000; }// returns p3 x vertex
    
    public double vy_e() { return Double.valueOf(Math.round(vy_e*100000))/100000; }// returns electron y vertex
    
    public double vy_p1() { return Double.valueOf(Math.round(vy_p1*100000))/100000; }// returns p1 y vertex
    
    public double vy_p2() { return Double.valueOf(Math.round(vy_p2*100000))/100000; }// returns p2 y vertex
    
    public double vy_p3() { return Double.valueOf(Math.round(vy_p3*100000))/100000; }// returns p3 y vertex
    
    public double vz_e() { return Double.valueOf(Math.round(vz_e*100000))/100000; }// returns electron z vertex
    
    public double vz_p1() { return Double.valueOf(Math.round(vz_p1*100000))/100000; }// returns p1 z vertex
    
    public double vz_p2() { return Double.valueOf(Math.round(vz_p2*100000))/100000; }// returns p2 z vertex
    
    public double vz_p3() { return Double.valueOf(Math.round(vz_p3*100000))/100000; }// returns p3 z vertex
    
    private static double particle_mass (int pid) {
	if (pid==11||pid==-11) { // electron is pid=11, positron is pid=-11
            return 0.0005109989461;
        } else if (pid == 111) {
            return 0.1349768;
	} else if (pid==211||pid==-211) { // pions
            return 0.139570;
	} else if (pid==321||pid==-321) { // kaons
            return 0.493677;
	} else if (pid==2212||pid==-2212) { // protons
            return 0.938272;
//            return 1.875;
	} else if (pid==113) { // rho0
            return 0.7754;
        }
            return -1;
    }
}
