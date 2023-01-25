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


public class Inclusive {
    
    protected byte helicity;
    protected int runnum;
    
    protected double test;
    
    protected int num_elec, num_particles;
    
    protected double Q2, W, gamma, nu, x, y, Mx;
    
    // depolarization vectors defining the polarization lost during the transfer from beam to 
    // the virtual photon. 
    // in ALU BSAs the twist 2 terms are modified by C/A and the twist 3 terms by W/A
    // B and V come in AUL
    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;
    
    protected double e_px, e_py, e_pz, e_p, e_e, e_theta, e_phi, vz_e; // electron kinematics
    
    
    public static boolean channel_test(Inclusive variables) {
//        if (variables.helicity==0){ 
//            System.out.println("You're returning false because helicity = 0. Is this data or MC?");
//            return false; }
        if (variables.Q2()<1) { return false; } 
        if (variables.W()<2) { return false; } 
//        if (variables.xF()<0.0) { return false; } 
//        else if (variables.y()>0.80) { return false; } 
//        else if (variables.Mx()<1.5) { return false; } 
	return true;
    }
    
    public Inclusive(DataEvent event, PhysicsEvent recEvent) {
        // provide the PDG PID of the two hadrons
        
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
//        Eb = 22;
        
        num_elec = recEvent.countByPid(11); // returns number of electrons
        num_particles = num_elec;
        
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
    
        
        vz_e = scattered_electron.vz();
        
        
        // missing mass calculations
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target);  
        Mx = lv_Mx.mass();     
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
        } else if ( runnum >= 11323 && runnum < 11571) {
            return helicity;
        }
        return -1*helicity;
    }
    
    public int get_runnum() { return runnum; }; // returns run number for polarizations and energy
    
    public int num_elec() { return num_elec; } // returns number of electrons
    
    public double test() { return test; } // returns test var

    public double Q2() { return Double.valueOf(Math.round(Q2*100000))/100000; } // returns Q2
    
    public double W() { return Double.valueOf(Math.round(W*100000))/100000; }// returns W
    
    public double gamma() { return Double.valueOf(Math.round(gamma*100000))/100000; } // returns gamma
    
    public double nu() { return Double.valueOf(Math.round(nu*100000))/100000; }// returns nu
    
    public double x() { return Double.valueOf(Math.round(x*100000))/100000; }// returns x
    
    public double y() { return Double.valueOf(Math.round(y*100000))/100000; }// returns y
    
    public double Mx() { return Double.valueOf(Math.round(Mx*100000))/100000; }// returns Mx(ep1p2)
    
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
  
    public double vz_e() { return Double.valueOf(Math.round(vz_e*100000))/100000; }// returns electron z vertex
    
    
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
