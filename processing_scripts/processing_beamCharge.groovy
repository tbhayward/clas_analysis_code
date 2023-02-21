/*
 * author Timothy B. Hayward
 * 
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


public class processing_beamCharge {

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
		}

		println(); println(); 

		int n_files;
		if ((args.length < 2)||(Integer.parseInt(args[1])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[1]);
		}

		int num_events = 0;
		int current_file = 0;
		String beamChargeList = '';
		float beamChargeMax = 0;
		float posHelbeamChargeTotal = 0;
		float negHelbeamChargeTotal = 0;
		int runnum;
		while (current_file < n_files) {
			beamChargeMax = 0;
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event = reader.getNextEvent(); 

			double charge_for_positive_helicities, charge_for_negative_helicities;
			double charge_for_run;
			while(reader.hasEvent()==true){
				num_events++; 
				if (num_events%1000000 == 0) { // not necessary, just updates output
					print("processed: "+num_events+" events, max beamCharge of current ");
					print("run = "+beamChargeMax+".\n");
				}

				if (event.hasBank("RUN::scaler")) {
					float beamCharge = event.getBank("RUN::scaler").getFloat("fcupgated",0);
					if (beamCharge > beamChargeMax) { beamChargeMax = beamCharge; }
    			}

    			if (event.hasBank("HEL::scaler")) {
    				if (event.getBank("HEL::scaler").getInt("helicity",0) == 1) {
    					float beamCharge = event.getBank("HEL::scaler").getFloat("fcupgated",0);
    					posHelbeamChargeTotal+=beamCharge;
					} else if (event.getBank("HEL::scaler").getInt("helicity",0) == -1) {
						float beamCharge = event.getBank("HEL::scaler").getFloat("fcupgated",0);
    					negHelbeamChargeTotal+=beamCharge;
					}
    			}

				// get run and event numbers
				event = reader.getNextEvent();
			    runnum = event.getBank("RUN::config").getInt('run',0); 

			}
			beamChargeList+="{"+runnum.toString()+", "+beamChargeMax.toString()+", ";
			beamChargeList+=posHelbeamChargeTotal.toString()+", "+negHelbeamChargeTotal.toString()+"}, "
			println(); println(); println();
			print("{"+beamChargeList+"};");

			beamChargeMax = 0;
			posHelbeamChargeTotal = 0;
			negHelbeamChargeTotal = 0;
		}

	}
}