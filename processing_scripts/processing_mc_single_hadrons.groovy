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

// filetype for gathering files in directory
import groovy.io.FileType;

public class processing_mc {

	public static double phi_calculation (double x, double y) {
		// tracks are given with Cartesian values and so must be converted to cylindrical
		double phi = Math.toDegrees(Math.atan2(x,y));
		phi = phi - 90;
		if (phi < 0) {
			phi = 360 + phi;
		}
		phi = 360 - phi;
		return phi;	
	}

	public static double theta_calculation (double x, double y, double z) {
		// convert cartesian coordinates to polar angle
		double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
		return (double) (180/Math.PI)*Math.acos(z/r);
	}



	public static void main(String[] args) {

		double scale = 2.0;


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

		// list.each {
		// 	println it.path
		// }

		println(); println(); println();

		String p1_Str;
		if (args.length < 2) {
			// assigns pi+ to p1
			println("WARNING: Specify a PDG PID for p1! Set to pi+. \n");
			p1_Str = "211";
		} else {
			p1_Str = args[1];
			println("Set p1 PID = "+p1_Str+"\n");
		}


		String output_file;
		if (args.length < 3) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"dihadron_dummy_out.txt\".\n");
			output_file = "dihadron_dummy_out.txt"
		} else {
			output_file = args[2];
		}

		int n_files;
		if ((args.length < 4)||(Integer.parseInt(args[3])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[3]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041);
		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
		EventFilter filter = new EventFilter("11:"+p1_Str+":"+":X+:X-:Xn");


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
			HipoDataEvent event = reader.getNextEvent(); 

			while(reader.hasEvent()==true){
				num_events++; 
				if (num_events%100000 == 0) { 
					print("processed: "+num_events+" events.     ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0);
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);
			    PhysicsEvent mc_Event = mc_fitter.getPhysicsEvent(event);

				if (filter.isValid(research_Event)) {

					HipoDataBank recBank = (HipoDataBank) event.getBank("REC::Event");
					HipoDataBank lundBank = (HipoDataBank) event.getBank("MC::Lund");
					HipoDataBank mcBank = (HipoDataBank) event.getBank("MC::Particle");

					int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 

					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {

						Particle exp_e = research_Event.getParticleByPid(11,0);
						Particle exp_p1 = research_Event.getParticleByPid(p1_Str.toInteger(),current_p1);

						Hadron variables = new Hadron(event, research_Event, 
							p1_Str.toInteger(), current_p1);
						Hadron mc_variables = new Hadron(event, mc_Event, 
							p1_Str.toInteger(), current_p1);

						if (variables.channel_test(variables)) {

							// lab kinematics data
							double e_p = variables.e_p();
							double e_theta = variables.e_theta();
							double p1_p = variables.p_p();
							double p1_theta = variables.p_theta();

							// lab kinematics MC
							double mc_e_p = mc_variables.e_p();
							double mc_e_theta = mc_variables.e_theta();
							double mc_p1_p = mc_variables.p_p();
							double mc_p1_theta = mc_variables.p_theta();

							// DIS variables data
							double Q2 = variables.Q2();
							double W = variables.W();
							double y = variables.y();
							double Mx = variables.Mx();

							// DIS variables MC
							double mc_Q2 = mc_variables.Q2();
							double mc_W = mc_variables.W(); 
							double mc_y = mc_variables.y();
							double mc_Mx = mc_variables.Mx();

							// SIDIS variables data
							double x = variables.x();
							double z = variables.z();
							double xF = variables.xF();
							double pT = variables.pT();
							double eta = variables.eta();
							double eta_gN = variables.eta_gN();
							double zeta = variables.zeta();

							// SIDIS variables MC
							double mc_x = mc_variables.x();
							double mc_z = mc_variables.z();
							double mc_xF = mc_variables.xF();
							double mc_pT = mc_variables.pT();
							double mc_eta = mc_variables.eta();
							double mc_eta_gN = mc_variables.eta_gN();
							double mc_zeta = mc_variables.zeta();

							// angles data
							double trento_phi = variables.phi();

							// angles MC
							double mc_trento_phi = mc_variables.phi();

							boolean matching_e = false;
							boolean matching_p1 = false;

							int matching_p1_pid = 0;
							int mc_p1_parent_index = 0;
							for (int current_part = 0; current_part < lundBank.rows(); current_part++) {
								int pid = lundBank.getInt("pid", current_part);
								if (matching_p1) { continue; }
								double mc_px = lundBank.getFloat("px", current_part);
								double mc_py = lundBank.getFloat("py", current_part);
								double mc_pz = lundBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
								double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
									exp_p1.pz());

								matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_p1) {
									matching_p1_pid = pid;
									mc_p1_parent_index = lundBank.getInt("parent", current_part)-1;
								}
								// if (Mx > 1.5 && xF > 0) {
								// 	if (lundBank.getInt("pid", mc_p1_parent_index)==2212 && matching_p1_pid == 211) {
								// 		println("PROTON DETECTED "+num_p1+" "+(current_part+1)+" "+mc_p1_parent_index+" "+matching_p1_pid+" ")
								// 		lundBank.show();
								// 	}
								// }
							}
							int mc_p1_parent = lundBank.getInt("pid", mc_p1_parent_index);


							matching_e = false;
							matching_p1 = false;

							int matching_e_pid = 0;
							int mc_e_parent_index = 0;
							for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
								int pid = mcBank.getInt("pid", current_part);
								if (matching_e) { continue; }
								double mc_px = mcBank.getFloat("px", current_part);
								double mc_py = mcBank.getFloat("py", current_part);
								double mc_pz = mcBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_e.px(), exp_e.py());
								double exp_theta = theta_calculation(exp_e.px(), exp_e.py(), exp_e.pz());

								matching_e = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_e) {
									matching_e_pid = pid;
								}
							}

							matching_p1_pid = 0;
							mc_p1_parent_index = 0;
							for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
								int pid = mcBank.getInt("pid", current_part);
								if (matching_p1) { continue; }
								double mc_px = mcBank.getFloat("px", current_part);
								double mc_py = mcBank.getFloat("py", current_part);
								double mc_pz = mcBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
								double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
									exp_p1.pz());

								matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_p1) {
									matching_p1_pid = pid;
								}
							}

							// if (Mx > 1.55 && Mx < 1.65) {
							// 	lundBank.show();
							// }

							file.append(e_p+" "+mc_e_p+" "+e_theta+" "+mc_e_theta+" ");
							file.append(p1_p+" "+mc_p1_p+" "+p1_theta+" "+mc_p1_theta+" ");
							file.append(Q2+" "+mc_Q2+" "+W+" "+mc_W+" ");
							file.append(Mx+" "+mc_Mx+" ");
							file.append(x+" "+mc_x+" "+y+" "+mc_y+" "+z+" "+mc_z+" ");
							file.append(xF+" "+mc_xF+" ");
							file.append(pT+" "+mc_pT+" ");
							file.append(zeta+" "+mc_zeta+" ");
							file.append(eta+" "+mc_eta+" ");
							file.append(trento_phi+" "+mc_trento_phi+" ");
							file.append(matching_e_pid+" "+matching_p1_pid+" ");
							file.append(mc_p1_parent+"\n");

						}
					}
				}
			}

			println(); println();
			print("1: e_p, 3: e_theta, 5: p1_p, 7: p1_theta, ");
			print("9: Q2, 11: W, 13: Mx, 15: x, 17: y, 19: z, 21: xF, 23: pT, 25: zeta, ");
			print("27: eta, 29: phi, ");
			print("30: matching e pid, 31: matching p1 pid, 32: p1 parent id\n");

			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("output file is: "+file);
			println();
		}
	}
}