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


public class processing_single_hadrons_RICH_studies {

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
			println("WARNING: Specify a PDG PID for p1! Set to pi+ (211). \n");
			p1_Str = "211";
		} else {
			p1_Str = args[1];
			println("Set p1 PID = "+p1_Str+"\n");
		}

		String output_file;
		if (args.length < 3) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".\n");
			output_file = "hadron_dummy_out.txt"
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
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041); // load my kinematic fitter/PID
		// GenericKinematicFitter research_fitter = new monte_carlo_fitter(10.6041);
		// GenericKinematicFitter research_fitter = new event_builder_fitter(10.6041); // load my kinematic fitter/PID
		// GenericKinematicFitter RICH_fitter = new RICH_fitter(10.6041); // load fitter using RICH
		// GenericKinematicFitter research_fitter = new proton_energy_loss_corrections_fitter(10.6041); // energy loss
		EventFilter filter = new EventFilter("11:"+p1_Str+":X+:X-:Xn"); // set filter for final states
		// setup QA database
		QADB qa = new QADB();

		int num_events = 0;
		int num_hadrons = 0;
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
				if (num_events%10000 == 0) { // not necessary, just updates output
					print("processed: "+num_events+" events. ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0); // collect info for QA
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

			    boolean process_event = false;
			    if (runnum == 11 || runnum >= 11571) { // if run number = 11 then it is MC and we don't use QA
			    	process_event = filter.isValid(research_Event);
			    } else {
			    	process_event = (filter.isValid(research_Event) && qa.OkForAsymmetry(runnum,evnum));
			    }
			    // println("Hello");
				if (process_event) {

					int num_p1 = research_Event.countByPid(p1_Str.toInteger());  // get # of particles w/ pid1

					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { // cycle over all combinations

						Hadron_RICH_studies variables = new Hadron_RICH_studies(event, research_Event, 
							p1_Str.toInteger(), current_p1);
						// this is my class for defining all relevant kinematic variables

						if (variables.channel_test(variables)) {
							int helicity = variables.get_helicity(); // helicity of event, might be 0
							
							// lab kinematics
							double e_p = variables.e_p(); // lab frame momentum
							double e_theta = variables.e_theta(); // lab polar angle
							double e_phi = variables.e_phi();  // lab azimuthal angle
							double p_phi = variables.p_phi(); // lab azimuthal angle
							double p_p = variables.p_p(); // lab momentum
							double p_theta = variables.p_theta(); // lab polar angle

							// DIS variables
							double Q2 = variables.Q2(); // my jar cuts on Q2 > 1
							double W = variables.W(); // my jar cuts on W > 1
							double y = variables.y();
							double Mx = variables.Mx();
							double Mx2 = variables.Mx2(); // missing mass square

							// SIDIS variables
							double x = variables.x(); // Bjorken-x
							double z = variables.z();
							double xF = variables.xF(); // Feynman-x
							double pT = variables.pT();
							double eta = variables.eta(); // rapidity 
							double zeta = variables.zeta();

							// angles 
							double phi = variables.phi(); // trento phi of the hadron

							// vertices 
							double vz_e = variables.vz_e();
							double vz_p = variables.vz_p();

							// rich tests
							int best_PID = variables.best_PID();
							int pindex = variables.pindex();
							int hindex = variables.hindex();
							int emilay = variables.emilay();
							int emico = variables.emico();
							int enico = variables.enico();
							int emqua = variables.emqua();

							double RQ_prob = variables.RQ_prob();
							double ReQ_prob = variables.ReQ_prob();
							double el_prob = variables.el_prob();
							double pi_prob = variables.pi_prob();
							double k_prob = variables.k_prob();
							double pr_prob = variables.pr_prob();
							double mchi2 = variables.mchi2();
							double best_etaC = variables.best_etaC();
							double best_c2 = variables.best_c2();
							double best_RL = variables.best_RL();
							double best_ntot = variables.best_ntot();
							double best_mass = variables.best_mass();

							double chi2pid = variables.chi2pid(); 
							double beta = variables.beta();  

							// append event to next line of the text file
							file.append(runnum+" "+evnum+" "+helicity+" ");
							file.append(e_p+" "+e_theta+" "+e_phi+" "+vz_e+" ");
							file.append(p_p+" "+p_theta+" "+p_phi+" "+vz_p+" ");
							file.append(Q2+" "+W+" "+Mx+" "+Mx2+" "+x+" "+y+" "+z+" "+xF+" ");
							file.append(pT+" "+zeta+" ");
							file.append(eta+" ");
							file.append(phi+" ");
							file.append(hindex+" "+pindex+" "+emilay+" "+emico+" "+enico+" ");
							file.append(emqua+" "+mchi2+" "+best_PID+" "+RQ_prob+" "+ReQ_prob+" ");
							file.append(el_prob+" "+pi_prob+" "+k_prob+" "+pr_prob+" ");
							file.append(best_etaC+" "+best_c2+" "+best_RL+" "+best_ntot+" "+best_mass+" ");
							file.append(chi2pid+" "+beta);
							file.append("\n");
						}
					}
				}
			}
			println(); println();
			print("1:runnum, 2:evnum, 3:helicity, ");
			print("4:e_p, 5:e_theta, 6:e_phi, 7:vz_e, ")
			print("8:p_p, 9:p_theta, 10:p_phi, 11:vz_p, ");
			print("12:Q2, 13:W, 14:Mx, 15:Mx2, 16:x, 17:y, 18:z, ");
			print("19:xF, ");
			print("20:pT, 21:zeta ");
			print("22:eta, ");
			print("23:phi, "); // this is the trento phi of the Hadron
			print("24:hindex, 25:pindex, 26:emilay, 27: emico, 28: enico, ");
			print("29:emqua, 30:mchi2, 31:best_PID, 32:RQ_prob, 33:ReQ_prob, ");
			print("34:el_prob, 35:pi_prob, 36:k_prob, 37:pr_prob, ");
			print("38:best_etaC, 39:best_c2, 40:best_RL, 41:best_ntot, 42:best_mass, ");
			print("43:chi2pid, 44:beta.");
			print("\n");

			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("output file is: "+file);
		}

	}
}