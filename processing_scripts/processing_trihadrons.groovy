/*
 * author Timothy B. Hayward
 * 
 * SIDIS dihadron 
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import analyzers.*;

// dilks CLAS QA analysis
import clasqa.QADB

// filetype for gathering files in directory
import groovy.io.FileType;


public class processing_trihadrons {

	public static void main(String[] args) {
		File[] hipo_list;
		if (args.length == 0) {
			// exits program if input directory not specified 
    	   	println("ERROR: Please enter a hipo file directory as the first argument");
    	  	System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list = directory.listFiles();
    	}

    	def list = []

		def dir = new File(args[0])
		dir.eachFileRecurse (FileType.FILES) { file ->
		  list << file
		  // println(file.toString()); println(); 
		}

		println(); println(); println();

		String p1_Str;
		if (args.length < 2) {
			// assigns pi+ to p1
			println("WARNING: Specify a PDG PID for p1! Set to proton. \n");
			p1_Str = "2212";
		} else {
			p1_Str = args[1];
			println("Set p1 PID = "+p1_Str+"\n");
		}

		String p2_Str;
		if (args.length < 3) {
			// assigns pi- to p2
			println("WARNING: Specify a PDG PID for p2! Set to pi+. \n");
			p2_Str = "211";
		} else {
			p2_Str = args[2];
			println("Set p2 PID = "+p2_Str+"\n");
		}

		String p3_Str;
		if (args.length < 3) {
			// assigns pi- to p2
			println("WARNING: Specify a PDG PID for p3! Set to pi-. \n");
			p3_Str = "-211";
		} else {
			p3_Str = args[3];
			println("Set p3 PID = "+p3_Str+"\n");
		}

		String output_file;
		if (args.length < 5) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"trihadron_dummy_out.txt\".\n");
			output_file = "trihadron_dummy_out.txt"
		} else {
			output_file = args[4];
		}

		int n_files;
		if ((args.length < 6)||(Integer.parseInt(args[5])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[5]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		println("hello world");

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041); // load my kinematic fitter/PID
		// GenericKinematicFitter research_fitter = new proton_energy_loss_corrections_fitter(10.6041);
		// GenericKinematicFitter research_fitter = new b2b_PRD_fitter(10.6041); // load my kinematic fitter/PID
		// GenericKinematicFitter research_fitter = new event_builder_fitter(10.6041); // load my kinematic fitter/PID
		EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":"+p3_Str+":X+:X-:Xn"); 
		// set filter for final states

		// setup QA database
		QADB qa = new QADB();

		int num_events = 0;
		int current_file = 0;
		// for (int current_file; current_file<n_files; current_file++) {
		while (current_file < n_files) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event;

			while(reader.hasEvent()==true){
				num_events++; 
				if (num_events%5000 == 0) { // not necessary
					print("Processed: "+num_events+" events. On file "+Integer.toString(current_file)
					+" out of "+n_files+". ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0); // collect info for QA
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

			    boolean process_event = false;
			    if (runnum == 11) { // if run number = 11 then it is MC and we don't use QA
			    	process_event = filter.isValid(research_Event);
			    } else {
			    	process_event = (filter.isValid(research_Event) && qa.OkForAsymmetry(runnum,evnum));
			    }

				if (process_event) {
					int num_p1 = research_Event.countByPid(p1_Str.toInteger());  // get # of particles w/ pid1
					int num_p2 = research_Event.countByPid(p2_Str.toInteger()); // get # of particles w/ pid2
					int num_p3 = research_Event.countByPid(p3_Str.toInteger()); // get # of particles w/ pid3

					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { // cycle over all combinations
						for (int current_p2 = 0; current_p2 < num_p2; current_p2++) {
							for (int current_p3 = 0; current_p3 < num_p3; current_p3++) {

								if (current_p1 == current_p2 && p1_Str.toInteger() == p2_Str.toInteger()) {continue; }
								if (current_p1 == current_p3 && p1_Str.toInteger() == p3_Str.toInteger()) {continue; }
								if (current_p2 == current_p3 && p2_Str.toInteger() == p3_Str.toInteger()) {continue; }

								Trihadrons variables = new Trihadrons(event, research_Event, 
									p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2, 
									p3_Str.toInteger(), current_p3);
								// this is my class for defining all relevant kinematic variables

								if (variables.channel_test(variables)) {
									int helicity = variables.get_helicity(); // helicity of event, might be 0

									// lab kinematics
									double e_p = variables.e_p();
									double e_theta = variables.e_theta();
									double e_phi = variables.e_phi();
									double p1_p = variables.p1_p();
									double p1_theta = variables.p1_theta();
									double p1_phi = variables.p1_phi();
									double p2_p = variables.p2_p();
									double p2_theta = variables.p2_theta();
									double p2_phi = variables.p2_phi();
									double p3_p = variables.p3_p();
									double p3_theta = variables.p3_theta();
									double p3_phi = variables.p3_phi();

									// DIS variables
									double Q2 = variables.Q2();
									double W = variables.W();
									double y = variables.y();
									double Mx = variables.Mx();
									double Mx1 = variables.Mx1();
									double Mx2 = variables.Mx2();
									double Mx3 = variables.Mx3();
									double Mx12 = variables.Mx12();
									double Mx23 = variables.Mx23();
									double Mx13 = variables.Mx13();

									// SIDIS variables
									double x = variables.x();
									double z = variables.z();
									double zeta = variables.zeta();
									double xF = variables.xF();
									double pT = variables.pT();
									double eta = variables.eta();
									double eta_gN = variables.eta_gN();

									// SIDIS trihadron variables
									double z1 = variables.z1();
									double z2 = variables.z2();
									double z3 = variables.z3();
									double z12 = variables.z12();
									double z13 = variables.z13();
									double z23 = variables.z23();

									double xF1 = variables.xF1();
									double xF2 = variables.xF2();
									double xF3 = variables.xF3();
									double xF12 = variables.xF12();
									double xF13 = variables.xF13();
									double xF23 = variables.xF23();

									double zeta1 = variables.zeta1();
									double zeta2 = variables.zeta2(); 
									double zeta3 = variables.zeta3();
									double zeta12 = variables.zeta12();
									double zeta13 = variables.zeta13();
									double zeta23 = variables.zeta23();

									double Mh = variables.Mh();
									double Mh12 = variables.Mh12();
									double Mh13 = variables.Mh13();
									double Mh23 = variables.Mh23();

									double pT1 = variables.pT1();
									double pT2 = variables.pT2();
									double pT3 = variables.pT3();
									double pT12 = variables.pT12();
									double pT13 = variables.pT13();
									double pT23 = variables.pT23();
									double pTpT = variables.pTpT();

									double eta1 = variables.eta1();
									double eta2 = variables.eta2();
									double eta3 = variables.eta3();
									double eta12 = variables.eta12();
									double eta13 = variables.eta13();
									double eta23 = variables.eta23();

									// angles 
									double phi1 = variables.phi1();
									double phi2 = variables.phi2();
									double phi3 = variables.phi3();
									double phi12 = variables.phi12();
									double phi13 = variables.phi13();
									double phi23 = variables.phi23();
									double Delta_phi12 = variables.Delta_phi12();
									double Delta_phi13 = variables.Delta_phi13();
									double Delta_phi23 = variables.Delta_phi23();
									double phih = variables.phih();
									double phiR = variables.phiR();
									double theta = variables.theta();

									// vertices 
									double vz_e = variables.vz_e();
									double vz_p1 = variables.vz_p1();
									double vz_p2 = variables.vz_p2();
									double vz_p3 = variables.vz_p3();

									// depolarization factors
									double Depolarization_A = variables.Depolarization_A();
									double Depolarization_B = variables.Depolarization_B();
									double Depolarization_C = variables.Depolarization_C();
									double Depolarization_V = variables.Depolarization_V();
									double Depolarization_W = variables.Depolarization_W();

									// append event to next line of the text file
									file.append(runnum+" "+evnum+" "+helicity+" ");
									file.append(e_p+" "+e_theta+" "+e_phi+" "+vz_e+" ");
									file.append(p1_p+" "+p1_theta+" "+p1_phi+" "+vz_p1+" ");
									file.append(p2_p+" "+p2_theta+" "+p2_phi+" "+vz_p2+" ");
									file.append(p3_p+" "+p3_theta+" "+p3_phi+" "+vz_p3+" ");
									file.append(Q2+" "+W+" ");
									file.append(Mx+" "+Mx1+" "+Mx2+" "+Mx3+" "+Mx12+" "+Mx13+" "+Mx23+" ");
									file.append(x+" "+y+" ");
									file.append(z+" "+z1+" "+z2+" "+z3+" "+z12+" "+z13+" "+z23+" ");
									file.append(zeta+" "+zeta1+" "+zeta2+" "+zeta3+" "+zeta12+" "+zeta13+" "+zeta23+" ");
									file.append(pT+" "+pT1+" "+pT2+" "+pT3+" "+pT12+" "+pT13+" "+pT23+" ");
									file.append(Mh+" "+Mh12+" "+Mh13+" "+Mh23+" ");
									file.append(xF+" "+xF1+" "+xF2+" "+xF3+" "+xF12+" "+xF13+" "+xF23+" ");
									file.append(eta+" "+eta1+" "+eta2+" "+eta3+" "+eta12+" "+eta13+" "+eta23+" ");
									file.append(phi1+" "+phi2+" "+phi3+" "+phi12+" "+phi13+" "+phi23+" "+phih+" "+phiR+" "+theta+" ");
									file.append(Delta_phi12+" "+Delta_phi13+" "+Delta_phi23+" ");
									file.append(Depolarization_A+" "+Depolarization_B+" "+Depolarization_C+" ");
									file.append(Depolarization_V+" "+Depolarization_W+"\n");
								}
							}
						}
					}
				}
			}
			println(); println();
			print("1:runnum, 2:evnum, 3:helicity, ");
			print("4: e_p, 5: e_theta, 6: e_phi, 7: vz_e, ");
			print("8: p1_p, 9: p1_theta, 10: p1_phi, 11: vz_p1, ");
			print("12: p2_p, 13: p2_theta, 14: p2_phi, 15: vz_p2, ");
			print("16: p3_p, 17: p3_theta, 18: p3_phi, 19: vz_p3 ");
			print("20: Q2, 21: W");
			print("22: Mx, 23: Mx1, 24: Mx2, 25: Mx3, 26: Mx12, 27: Mx13, 28: Mx23, ");
			print("29: x, 30: y, ");
			print("31: z, 32: z1, 33: z2, 34: z3, 35: z12, 36: z13, 37: z23, ");
			print("38: zeta, 39: zeta1, 40: zeta2, 41: zeta3, 42: zeta12, 43: zeta13, 44: zeta23, ");
			print("45: pT, 46: pT1, 47: pT2, 48: pT3, 49: pT12, 50: pT13, 51: pT23, ");
			print("52: Mh, 53: Mh12, 54: Mh13, 55: Mh23, ");
			print("56: xF, 57: xF1, 58: xF2, 59: xF3, 60: xF12, 61: xF13, 62: xF23, ");
			print("63: eta, 64: eta1, 65: eta2, 66: eta3, 67: eta12, 68: eta13, 69: eta23, ");
			print("70: phi1, 71: phi2, 72: phi3, 73: phi12, 74: phi13, 75: phi23, 76: phih, 77: phiR, 78: theta, ");
			print("79: Delta_phi12, 80: Delta_phi13, 81: Delta_phi23, ")
			print("82: DepA, 83: DepB, 84: DepC, 85: DepV, 86: DepW.");

			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("Set p2 PID = "+p2_Str+"\n");
			println("output file is: "+file);
		}

	}
}