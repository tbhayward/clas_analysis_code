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


public class processing_inclusive {

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

		String output_file;
		if (args.length < 2) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"inclusive_dummy_out.txt\".\n");
			output_file = "inclusive_dummy_out.txt"
		} else {
			output_file = args[1];
		}

		int n_files;
		if ((args.length < 3)||(Integer.parseInt(args[2])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[2]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041); // load my kinematic fitter/PID
		// GenericKinematicFitter research_fitter = new monte_carlo_fitter(10.6041);
		// GenericKinematicFitter research_fitter = new event_builder_fitter(10.6041); // load my kinematic fitter/PID
		// GenericKinematicFitter RICH_fitter = new RICH_fitter(10.6041); // load fitter using RICH
		// GenericKinematicFitter research_fitter = new proton_energy_loss_corrections_fitter(10.6041); // energy loss
		EventFilter filter = new EventFilter("11:X+:X-:Xn"); // set filter for final states
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
				if (process_event) {

					Inclusive variables = new Inclusive(event, research_Event);
					// this is my class for defining all relevant kinematic variables

					if (variables.channel_test(variables)) {
						int helicity = variables.get_helicity(); // helicity of event, might be 0
						
						// lab kinematics
						double e_p = variables.e_p(); // lab frame momentum
						double e_theta = variables.e_theta(); // lab polar angle
						double e_phi = variables.e_phi();  // lab azimuthal angle
						double vz_e = variables.vz_e();

						// DIS variables
						double Q2 = variables.Q2(); // my jar cuts on Q2 > 1
						double W = variables.W(); // my jar cuts on W > 1
						double x = variables.x(); // Bjorken-x
						double y = variables.y();
						double Mx = variables.Mx();

						// append event to next line of the text file
						file.append(runnum+" "+evnum+" "+helicity+" ");
						file.append(e_p+" "+e_theta+" "+e_phi+" "+vz_e+" ");
						file.append(Q2+" "+W+" "+Mx+" "+x+" "+y);
						file.append("\n");
					}
				}
			}
			println(); println();
			print("1:runnum, 2:evnum, 3:helicity, ");
			print("4:e_p, 5:e_theta, 6:e_phi, 7:vz_e, ")
			print("8:Q2, 9:W, 10:Mx, 11:x, 12:y");
			print("\n");

			println(); println();
			println("output file is: "+file);
		}

	}
}