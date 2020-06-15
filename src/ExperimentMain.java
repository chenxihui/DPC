import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultEdge;


public class ExperimentMain {
	public static void main(String[] args) throws IOException {
		
		/*
		 * Re-identification attack
		 */
		//reIdentificationAttack attack = new reIdentificationAttack();
		//attack.attackSingleFile("epinions8000", "epinions8000_w_triangle_final.txt");
		//attack.attackSingleFile("facebook", "facebook_w_triangle_test_final.txt");
		//attack.attackSingleFile("petster", "petster_w_triangle.txt");

		//attack.attackInBatchForExperiment(args[0], args[1], args[2]);

		//attack.attackInBatchForExperiment("facebook_combined_experiment", "facebook", "2.0");

		ExperimentMain expMain = new ExperimentMain();
		//expMain.sampleSubgraph("epinions_combined_mutual.txt", "epinions", 10000);
		//expMain.testDiffSenarios(args);
		
		//expMain.averageThetaFCompare("ipython/epinions8000_DCSBM", "epinions8000", "DCSBM");
		//expMain.averageThetaFCompare("ipython/epinions8000_combined_experiment", "epinions8000", "");
		//expMain.averageThetaFCompare("ipython/epinions8000_TriCycle", "epinions8000", "TriCycle");
		//expMain.averageThetaFCompare("ipython/epinions_TriCycle_combined_experiment", "epinions", "TriCycle");
		//expMain.averageThetaFCompare("ipython/epinions_combined_experiment", "epinions", "");
		//expMain.averageThetaFCompare("ipython/epinions_DCSBM_combined_experiment", "epinions", "DCSBM");
		//expMain.averageThetaFCompare("ipython/facebook_TriCycle_combined_experiment", "facebook", "TriCycle");
		//expMain.averageThetaFCompare("ipython/facebook_DCSBM_combined_experiment", "facebook", "DCSBM");
		//expMain.averageThetaFCompare("ipython/facebook_combined_experiment", "facebook", "");
		//expMain.averageThetaFCompare("epinions_DCSBM_combined_experiment", "epinions", "DCSBM");

		/*
		String dataname = args[0];
		int start_no= Integer.parseInt(args[1]);
		int end_no = Integer.parseInt(args[2]);
		int method = Integer.parseInt(args[3]);
		*/
		
		String dataname = "facebook";
		int start_no= 0;
		int end_no = 0;
		int method = 1;
		
		switch(method) {
		case 1: 
			expMain.generateDPCAGMGraphInBatch(dataname, start_no, end_no);
			break;
		case 2:
			expMain.generateDPDCSBMGraphInBatch(dataname, start_no, end_no);
			break;
		case 3:
			expMain.generateDPTriCycleInBatch(dataname, start_no, end_no);
			break;
			
		}
		
//		expMain.generateDPDCSBMGraphInBatch("facebook", 0, 40);
		
	}
	

	public void generateDPCAGMGraphInBatch(String dataname, int start_run_no, int end_run_no) throws IOException {
		double [] elist = {2.0, 3.0, 4.0, 5.0};

		//double [] elist = {0.0};

		String filedir =  dataname + "_combined_experiment" ;
	    File directory = new File(filedir);
	    if (! directory.exists()){
	        directory.mkdir();
	    }
		for (int num_run = start_run_no; num_run <= end_run_no; num_run ++ ) {
			for (double e : elist) {
				double part_epsilon = Math.round(1.0 * e* 100/2) * 1.0 / 100;

				String pFile = filedir+"/" + dataname + "_partition_mod_" + e + "_"+ num_run + ".txt";
				String inputFileName = filedir +"/" + dataname + "InputFile_" + e + "_" + num_run + ".txt";	

				if (e - 0.0 > 0.0001) {
					File f = new File(pFile);
					if (!f.exists()) {
						Partition rp;
						do {
							//rp = this.run_moddivisive(pFile, dataname, part_epsilon);
							rp = this.run_moddivisive_abt(pFile, dataname, part_epsilon);
						}while (rp.partition.size()<=1); 
					}
				}else
					pFile = filedir+"/" + dataname + "_partition_mod_" + e + "_0"+ ".txt";

				String graphFile = dataname + "_combined.txt";
				String abtFile = "attribute_" + dataname + "_combined.txt";

				dpCAGM dpc ;
				if(!new File(inputFileName).exists()) {
					if (e - 0.0 > 0.0001)
						dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e , true, true, false, false);
					else
						dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e , false, true, false, false);
					dpc.generateCAGMInputFile(inputFileName);
				}
			}
		}
		for (int i = start_run_no; i<= end_run_no; i++) {
			for (double e:elist) {
				String filename = String.format("%s/%sInputFile_%d.0_%d.txt", filedir, dataname, (int)e, i);
				if (!new File(filename).exists() ) {
					System.out.println(filename );
					continue;
				}
				
				String file_wt_file = filedir + "/" + dataname + "_wt_triangles_" +(int)e + ".0_final_" + i +".txt";
				String file_w_file = filedir + "/" + dataname + "_w_triangles_" +(int)e + ".0_final_" + i +".txt";
				if(!new File(file_w_file).exists()) {
					CAGMGenerator dgen = new CAGMGenerator(dataname, filename);
					dgen.Generate_graph_triangles_iterative(file_wt_file, file_w_file, 0.99);
				}
			}
		}

		
		String graphFile = dataname + "_combined.txt"; 
		String abtFile = "attribute_"+dataname + "_combined.txt"; 
		
		for (int i = start_run_no; i<= end_run_no; i++) 
			for (double e : elist) {
				String randomGraphFile = filedir + "/" + dataname + "_w_triangles_" +(int)e + ".0_final_" + i +".txt";
				
				String pFile = filedir + "/" + dataname + "_partition_mod_" + e + "_"
						+i+".txt";
				if (e-0.0 < 0.0000001)
					pFile = filedir + "/" + dataname + "_partition_mod_" + e + "_0" +".txt";

				
				File f = new File(randomGraphFile);
				if (!f.exists())
					continue;
				f = new File(pFile);
				if(!f.exists())
					continue;
				dpCAGM dpc;
				if (e-0.0 < 0.000001)
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, false , true, false, false);
				else
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, true, true , false, false);

				dpc.genDpCAGM(randomGraphFile, i);
			}
		
		for (int i = start_run_no; i<= end_run_no; i++) 
			for (double e : elist) {
				String abt = filedir+"/" + dataname+"_attribute_" + e + "_" +i+ ".txt";
				String Ffile =filedir+"/"+  dataname + "_F_" + e + "_" + i + ".txt";
				String Fg1file =filedir+"/"+ dataname + "_F_g1_" + e + "_" + i + ".txt";
				String supFile =filedir+"/"+ dataname+"_sup_" + e + "_" + i + ".txt";
				
				if(!new File(abt).exists() || !new File(Ffile).exists() || 
						!new File(Fg1file).exists() || !new File(supFile).exists())
					continue;

				String inputFile = String.format("%s/%sInputFile_%d.0_%d.txt", filedir, dataname, (int)e, i);

				String file_wt= filedir + "/" + dataname + "_wt_triangles_" +(int)e + ".0_final_dp_" + i +".txt";
				String file_w= filedir + "/" + dataname + "_w_triangles_" +(int)e + ".0_final_dp_" + i +".txt";

				if(!new File(file_w).exists()) {
					DPCAGMGenerator dp = new DPCAGMGenerator(dataname, inputFile, abt, Ffile, Fg1file, supFile);
					dp.Generate_graph_triangles_iterative(file_wt, file_w, 0.99);
				}
		}
		
	}

	public void combined_process_to_inputFileGen(String dataname, double epsilon, 
			int start_num_run, int end_num_run) throws IOException {
		/*
		 * calculate the partition file
		 */

		for (int num_run = start_num_run; num_run < end_num_run; num_run++) {
			double part_epsilon = Math.round(1.0 * epsilon * 100/2) * 1.0 / 100;

			String pFile = dataname + "_partition_mod_" + epsilon + "_"+ num_run + ".txt";
			this.run_moddivisive(pFile, dataname, part_epsilon);

			/*
			 * generate the CAGM input file named by facebookInputFile.txt
			 */

			String graphFile = dataname + "_combined.txt";
			String abtFile = "attribute_" + dataname + "_combined.txt";

			dpCAGM dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, epsilon, true, true, false, false);
			String inputFileName = dataname + "InputFile_" + epsilon + "_" + num_run + ".txt";
			dpc.generateCAGMInputFile(inputFileName);
		}
		
	}

	public void generateDPDCSBMGraphInBatch(String dataname, int start_run_no, int end_run_no) throws IOException {
		double [] elist = {2.0, 3.0, 4.0, 5.0};

		String filedir =  dataname + "_DCSBM_combined_experiment" ;
	    File directory = new File(filedir);
	    if (! directory.exists()){
	        directory.mkdir();
	    }

		for (int num_run = start_run_no; num_run <= end_run_no; num_run ++ ) {
			for (double e : elist) {
				double part_epsilon = Math.round(1.0 * e* 100/2) * 1.0 / 100;

				String pFile = filedir+"/" + dataname + "_partition_mod_" + e + "_"+ num_run + ".txt";
				File f = new File(pFile);
				if (!f.exists()) {
					Partition rp;
					do {
						rp = this.run_moddivisive(pFile, dataname, part_epsilon);
					}while (rp.partition.size()==1); 
				}

				String graphFile = dataname + "_combined.txt";
				String abtFile = "attribute_" + dataname + "_combined.txt";
				
				
				String inputFileName = filedir +"/" + dataname + "_DCSBMInputFile_" + e + "_" + num_run + ".txt";
				String inputEdgeFile = filedir + "/" + dataname + "_DCSBMInputEdgeFile"+ e + "_" + num_run + ".txt";

				if (new File(inputFileName).exists() && new File(inputEdgeFile).exists())
					continue;
				dpCAGM dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e , true, false,false, true);
				dpc.GenerateInputFileDCSBM(inputFileName, inputEdgeFile);
			}
		}
		for (int i = start_run_no; i<= end_run_no; i++) {
			for (double e:elist) {
				String filename = String.format("%s/%s_DCSBMInputFile_%d.0_%d.txt", filedir, dataname, (int)e, i);
				String fileEdgename = String.format("%s/%s_DCSBMInputEdgeFile%d.0_%d.txt", filedir, dataname, (int)e, i);
				if (!new File(filename).exists() || !new File(fileEdgename).exists()) {
					System.out.println(filename + " " + fileEdgename);
					continue;
				}
				
				String fileOut = filedir + "/" + dataname + "_DCSBM_graph_" +(int)e + ".0_final_" + i +".txt";
				if (new File(fileOut).exists())
					continue;
				DCSBM_Generator dgen = new DCSBM_Generator();
				dgen.generateGraph(filename, fileEdgename, fileOut);
			}
		}

		String graphFile = dataname + "_combined.txt"; 
		String abtFile = "attribute_"+dataname + "_combined.txt"; 
		
		for (int i = start_run_no; i<= end_run_no; i++) 
			for (double e : elist) {
				String randomGraphFile = filedir +"/" + dataname + "_DCSBM_graph_" + e+
						"_final_"+i+".txt";
				
				String pFile = filedir + "/" + dataname + "_partition_mod_" + e + "_"
						+i+".txt";
				
				File f = new File(randomGraphFile);
				if (!f.exists())
					continue;
				f = new File(pFile);
				if(!f.exists())
					continue;
				dpCAGM dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, true, false, false, true);

				dpc.genDpCAGM(randomGraphFile, i);
			}
		
		for (int i = start_run_no; i<= end_run_no; i++) 
			for (double e : elist) {
				String abt =filedir+"/"+ dataname+"_DCSBM_attribute_" + e + "_" +i+ ".txt";
				String Ffile =filedir+"/"+ dataname + "_DCSBM_F_" + e + "_" + i + ".txt";
				String Fg1file =filedir+"/"+ dataname + "_DCSBM_F_g1_" + e + "_" + i + ".txt";
				String supFile =filedir+"/"+ dataname+"_DCSBM_sup_" + e + "_" + i + ".txt";
				if(!new File(abt).exists() || !new File(Ffile).exists() || 
						!new File(Fg1file).exists() || !new File(supFile).exists())
					continue;
	
				DPDCSBM_Generator dp = new DPDCSBM_Generator(dataname, abt, Ffile, Fg1file, supFile);
				String inputFile = String.format("%s/%s_DCSBMInputFile_%d.0_%d.txt", filedir, dataname, (int)e, i);
				String numEdgeFile = String.format("%s/%s_DCSBMInputEdgeFile%d.0_%d.txt", filedir, dataname, (int)e, i);
				String fileOut = filedir + "/" + dataname + "_DCSBM_graph_" +(int)e + ".0_final_dp_" + i +".txt";
				if (new File(fileOut).exists())
					continue;
				dp.generateGraph(inputFile, numEdgeFile, fileOut);
		}
		
	}
	public void generateRandomGraphInBatch(String dataname, String filedir, int start_run_no, int end_run_no) throws IOException {
		double[] elist = {0.0, 2.0, 3.0, 4.0, 5.0};
		
		for (int i = start_run_no; i<= end_run_no; i++) {
			for (double e:elist) {
				String filename = String.format("%s/%sInputFile_%d.0_%d.txt", filedir, dataname, (int)e, i);
				if (!new File(filename).exists())
					continue;
				
				String file_wt_file = dataname + "_wt_triangles_" +(int)e + ".0_final_" + i +".txt";
				String file_w_file = dataname + "_w_triangles_" +(int)e + ".0_final_" + i +".txt";
				
				CAGMGenerator gen = new CAGMGenerator(dataname, filename);
				
				if (gen.comList.keySet().size() ==1 )
					continue;
				gen.Generate_graph_triangles_iterative(file_wt_file, file_w_file, 0.98);
				//(file_wt_file, file_w_file);
				
			}
		}
	}
	
	public void run_dpcCAGM_batch(String dataname, int start_from, int end_at) throws IOException{
		String graphFile = dataname + "_combined.txt"; 
		String abtFile = "attribute_"+dataname + "_combined.txt"; 
		
		double[] elist = {2.0, 3.0, 4.0, 5.0};

		for (double e : elist)
			for (int i = start_from; i<= end_at; i++) {
				String randomGraphFile = dataname + "_combined_experiment/" + dataname + "_w_triangles_" + e+
						"_final_"+i+".txt";
				
				String pFile = dataname + "_combined_experiment/" + dataname + "_partition_mod_" + e + "_"
						+i+".txt";
				
				File f = new File(randomGraphFile);
				if (!f.exists())
					continue;
				f = new File(pFile);
				if(!f.exists())
					continue;
				dpCAGM dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, true, true, false, false);

				dpc.genDpCAGM(randomGraphFile, i);
			}
	}
	public void run_dpcCAGM_batch_otherModels(String dataname, int start_from, int end_at) throws IOException{
		String graphFile = dataname + "_combined.txt"; 
		String abtFile = "attribute_"+dataname + "_combined.txt"; 
		
		double[] elist = {0.0, 2.0, 3.0, 4.0, 5.0};
		

		for (double e : elist)
			for (int i = start_from; i<= end_at; i++) {
				String randomGraphFile = dataname + "_TriCycle_/" + dataname + "_TriCycle_w_triangles_" + e+
						"_final_"+i+".txt";
				
				String pFile = dataname + "_combined_experiment/" + dataname + "_partition_mod_" + e + "_"
						+i+".txt";
				
				File f = new File(randomGraphFile);
				if (!f.exists())
					continue;
				dpCAGM dpc;
				if (Math.abs(e-0.0)<0.0001)
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, false, false, true, false);
				else
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, true, false, true, false);
					
				dpc.genDpCAGM(randomGraphFile, i);
			}
	}
	public  void run_dpCAGM(String dataname, double epsilon) throws IOException {
		String graphFile = dataname + "_combined.txt"; 
		String abtFile = "attribute_"+dataname + "_combined.txt" ; 
		String pFile = dataname + "_partition_mod_0.37_6.txt"; 
		String randomGraphFile = dataname + "_graph_with_triangles.txt";
		dpCAGM dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, epsilon,true, true, false, false);

		String inputFileName = dataname + "InputFile_0.37_6.txt";
		//dpc.generateCAGMInputFile(inputFileName);
		dpc.genDpCAGM(randomGraphFile, 0);
	}
		
	public Partition run_moddivisive(String pFile, String dataname, double epsilon) throws IOException {
		graph g = new graph(); 
		double epsilon_m = 0.01;
		double global_sensitivity = 0.0003;

		int maxL = 5;
		int k = 2;
		int K = 500;
        /*
         * main experiments starts here
         */
		System.out.println("Start reading graph...");
		g.readGraphOnly(dataname + "_combined.txt");
		System.out.println("Graph read...");

		//double epsilon = Math.round(0.5 * e * 100) * 1.0 / 100;
//				0.3 * Math.log(g.getGraph().vertexSet().size());
		
		// create a ModDivisive object with epsilon, espilon_m, global sensitivity, maxL and k,
		ModDivisiv mod = new ModDivisiv(epsilon, epsilon_m, global_sensitivity, maxL, k, K);

		Partition rp = mod.ModDivisiveMain(g, dataname, 1);
		System.out.println("The final partition is " + rp.toString());	
		FileWriter fout = new FileWriter(pFile);
		fout.write(rp.toString());
		fout.close();
		return rp;
	}

    public Partition run_moddivisive_abt(String pFile, String dataname, double epsilon) throws IOException{
    	graph g = new graph(); 
		attributeGraph ag = new attributeGraph();
		double epsilon_m = 0.01;

		int maxL = 5;
		int k = 2;
		int K = 500;	
        /*
         * main experiments starts here
         */

		String dataName = "facebook";

		System.out.println("Start reading graph...");
	    g.readGraphOnly(dataName + "_combined.txt");

		System.out.println("Graph read...");
		System.out.println("Start reading attribute graph...");
		ag.readGraphOnly(dataName + "_combined_attribute.txt", g);
	    System.out.println("Attribute graph read...");
		double global_sensitivity = 0.98* 0.0003 + 0.02 *  (60.0/g.getGraph().vertexSet().size());

					
					// create a ModDivisive object with epsilon, espilon_m, global sensitivity, maxL
					// and k,
		ModDivisiv mod = new ModDivisiv(epsilon, epsilon_m, global_sensitivity, maxL, k, K, 0.98);

		Partition rp = mod.ModDivisiveMain_attribute(ag, g, dataName, 0.98, 1);
		System.out.println("The final partition is " + rp.toString());	
		FileWriter fout = new FileWriter(pFile);
		fout.write(rp.toString());
		fout.close();

		return rp;
				

    }

    public void generateDPTriCycleInBatch(String dataname, int start_run_no, int end_run_no) throws IOException {
		double [] elist = {0.0, 2.0, 3.0, 4.0, 5.0};
		
		/*
		 * calculate the partition file
		 */
		String filedir =  dataname + "_TriCycle_combined_experiment" ;

		for (int num_run = start_run_no; num_run <= end_run_no; num_run ++ ) {
			for (double e : elist) {
				String pFile = filedir+"/" + dataname + "_TriCycle_partition_mod_" + e + "_"+ num_run + ".txt";

				String graphFile = dataname + "_combined.txt";
				String abtFile = "attribute_" + dataname + "_combined.txt";
				dpCAGM dpc;
				if (Math.abs(e - 0.0)<0.0001)
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, false, false, true, false);
				else
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, true, false, true, false);

				String inputFileName = filedir + "/" + dataname + "_TriCycleInputFile" + e + "_" + num_run + ".txt";

				if (!new File(inputFileName).exists())
					dpc.generateCAGMInputFile(inputFileName);
			}
		}
		
		for (int i = start_run_no; i<= end_run_no; i++) {
			for (double e:elist) {
				String filename = String.format("%s/%s_TriCycleInputFile%d.0_%d.txt", filedir, dataname, (int)e, i);
				if (!new File(filename).exists() ) {
					System.out.println(filename );
					continue;
				}
				
				String file_wt_file = filedir + "/" + dataname + "_TriCycle_wt_triangles_" +(int)e + ".0_final_" + i +".txt";
				String file_w_file = filedir + "/" + dataname + "_TriCycle_w_triangles_" +(int)e + ".0_final_" + i +".txt";
				if(!new File(file_w_file).exists()) {
					TriCycle dgen = new TriCycle(dataname, filename);
					dgen.generate_graph_triangles(file_wt_file, file_w_file);
				}
			}
		}

		String graphFile = dataname + "_combined.txt"; 
		String abtFile = "attribute_"+dataname + "_combined.txt"; 
		
		for (int i = start_run_no; i<= end_run_no; i++) 
			for (double e : elist) {
				String randomGraphFile = filedir + "/" + dataname + "_TriCycle_w_triangles_" +(int)e + ".0_final_" + i +".txt";
				
				String pFile = filedir + "/" + dataname + "_partition_mod_" + e + "_" +i+".txt";
				if (e-0.0 < 0.0000001)
					pFile = filedir + "/" + dataname + "_partition_mod_" + e + "_0" + ".txt";

				if (!new File(randomGraphFile).exists())
					continue;
				dpCAGM dpc;
				if (Math.abs(e-0.0)<0.0001)
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, false, false, true, false);
				else
					dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, e, true, false, true, false);
						
				dpc.genDpCAGM(randomGraphFile, i);	
			}
		
		for (int i = start_run_no; i<= end_run_no; i++) 
			for (double e : elist) {
				String abt =filedir+"/"+ dataname+"_TriCycle_attribute_" + e + "_" +i+ ".txt";
				String Ffile =filedir+"/"+ dataname + "_TriCycle_F_" + e + "_" + i + ".txt";
				String Fg1file =filedir+"/"+ dataname + "_TriCycle_F_g1_" + e + "_" + i + ".txt";
				String supFile =filedir+"/"+ dataname+"_TriCycle_sup_" + e + "_" + i + ".txt";
				
				String inputFile = String.format("%s/%s_TriCycleInputFile%d.0_%d.txt", filedir, dataname, (int)e, i);
				DPTriCycle dp = new DPTriCycle(dataname, inputFile, abt, Ffile, Fg1file, supFile);

				String file_wt= filedir + "/" + dataname + "_TriCycle_wt_triangles_" +(int)e + ".0_final_dp_" + i +".txt";
				String file_w= filedir + "/" + dataname + "_TriCycle_w_triangles_" +(int)e + ".0_final_dp_" + i +".txt";
				
				if(!new File(file_w).exists())
					dp.generate_graph_triangles(file_wt, file_w);
		}
		
	}
	public void run_moddivisive(String [] args) throws IOException {
		graph g = new graph(); 
	
		String dataName = "epinions";
		double epsilon =  1.0;
		double epsilon_m = 0.1;
		double global_sensitivity = 0.00001;
		int maxL = 5;
		int k = 2;
		int K = 20;
        /*
         * main experiments starts here
         */
		System.out.println("Start reading graph...");
		g.readGraphOnly(dataName + "_combined.txt");
		System.out.println("Graph read...");

		// create a ModDivisive object with epsilon, espilon_m, global sensitivity, maxL and k,
		ModDivisiv mod = new ModDivisiv(epsilon, epsilon_m, global_sensitivity, maxL, k, K);

		
		Partition rp = mod.ModDivisiveMain(g, dataName, 1);
		System.out.println("The final partition is " + rp.toString());	
	}
	public static void test() {
		
double [] seq = {4.0, 5.0, 1.0, 9.0, 4.0, 3.0, 4.0, 6.0, 10.0, 11.0, 15};
		
		double [] s = new double[seq.length];
		double  max = 0;
		double [] sum = new double[seq.length];
		sum[0] = seq[0];
		for (int i = 1; i<seq.length; i++)
			sum[i] = sum[i-1] + seq[i]; 
		
		for (int k = 0; k< seq.length; k++) {
			double min = Double.MAX_VALUE;
			for(int j = k; j<seq.length; j++) {
				if (k>0)
					sum[j] -= seq[k-1];
				if(sum[j]/(j-k+1) < min)
					min = sum[j]/(j-k+1);
			}
			if(max < min)
				max = min;
			s[k] = max;
		}
			
		
		for (double i : s)
			System.out.println(i);
		
	}
	
	public void averageThetaFCompare(String dir, String dataname, String Model) throws IOException {
		double[] elist = {2.0, 3.0, 4.0, 5.0};
		/*if (dataname.contentEquals("epinions")) {
			elist[0] = 6.0;
			elist[1] = 7.0;
			elist[2] = 8.0;
			elist[3] = 9.0;
		}*/
			
		HashMap<Double, Double> result = new HashMap<Double, Double>();
		
		for (double e: elist) {
			int count = 0;
			double sum = 0;
			for(int i = 0; i<=99; i++) {
				String dp_final_file = dir + File.separator + dataname + "_w_triangles_"
						+ e + "_final_dp_" + i + ".txt";

				if(Model.contentEquals("TriCycle"))
					dp_final_file = dir + File.separator + dataname + "_TriCycle_w_triangles_"
						+ e + "_final_dp_" + i + ".txt";
				else if(Model.contentEquals("DCSBM"))
					dp_final_file = dir + File.separator + dataname + "_DCSBM_graph_"
						+ e + "_final_dp_" + i + ".txt";

				if (!new File(dp_final_file).exists())
					continue;
				double r = this.thetaFCompare(dir, dataname, Model,  e, i);
				sum += r;
				count ++ ;
			}
			result.put(e, sum/count);
		}
		
		for (double e: result.keySet())
			System.out.println(e + ": " + result.get(e) + "\n");
		
	}
	public double thetaFCompare(String resultdir, String dataname, String Model, double epsilon, int serialno) throws IOException {
		ResultProcessor rp = new ResultProcessor();
		
		String graph_file = dataname + "_combined.txt";
		String dp_final_file = resultdir + File.separator + dataname + "_w_triangles_"
				+ epsilon + "_final_dp_" + serialno + ".txt";
		String attribute_file = resultdir + File.separator + dataname + "_attribute_" + epsilon + "_" + serialno+ ".txt";
		String original_abt_file = "attribute_" + dataname + "_combined.txt"; 
	
		int num_abt = 50;
		if (Model.contentEquals("TriCycle")){
			dp_final_file = resultdir + File.separator + dataname + "_TriCycle_w_triangles_"
				+ epsilon + "_final_dp_" + serialno + ".txt";
			attribute_file = resultdir + File.separator + dataname + "_TriCycle_attribute_" + epsilon + "_" + serialno+ ".txt";
			
		
		} else if(Model.contentEquals("DCSBM")) {
			dp_final_file = resultdir + File.separator + dataname + "_DCSBM_graph_"
				+ epsilon + "_final_dp_" + serialno + ".txt";
			attribute_file = resultdir + File.separator + dataname + "_DCSBM_attribute_" + epsilon + "_" + serialno+ ".txt";
		}
				
		if (dataname.contentEquals("petster"))
			num_abt = 13;	

//		String partition_file_ori =  dataname + "_partition_no_privacy.txt";
		//String partition_file_gen =  partGenDir + File.separator + dataname + "_partition_cesna_" + epsilon + "_"
		//		+ serialno+ "cmtyvv.txt";
		//if(Model.contentEquals("TriCycle"))
		//	partition_file_ = "";
		String partition_file_ori =  dataname + "_partition_cesna.txt";

/*		double res = rp.compareThetaF_ratio("facebook_Graph_with_triangles_final_dp.txt", 
				"facebook_combined.txt", "facebook_attribute.txt", "attribute_facebook_combined.txt", 
				"facebook_partition_mod.txt", 50);
*/		
		
		System.out.println(attribute_file);
		double res = rp.compareThetaF_ratio(dp_final_file, graph_file, attribute_file, original_abt_file, 
				partition_file_ori, num_abt);
		System.out.println(res);
		return res;
		
	}

	public void combined_process_to_inputFileGen_otherModels(String dataname, double epsilon, 
			int num_run_total, int start_num_run) throws IOException {
		/*
		 * calculate the partition file
		 */
		for (int num_run = start_num_run; num_run < num_run_total; num_run++) {
			//double part_epsilon = Math.round(1.0 * epsilon * 100/2) * 1.0 / 100;

			String pFile = dataname + "_TriCycle_partition_mod_" + epsilon + "_"+ num_run + ".txt";
			//this.run_moddivisive(pFile, dataname, part_epsilon);

			/*
			 * generate the CAGM input file named by facebookInputFile.txt
			 */

			String graphFile = dataname + "_combined.txt";
			String abtFile = "attribute_" + dataname + "_combined.txt";
			dpCAGM dpc;
			if (Math.abs(epsilon-0.0)<0.0001)
			    dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, epsilon, false, false, true, false);
			else
				dpc = new dpCAGM(dataname, graphFile, abtFile, pFile, epsilon, true, false, true, false);

			String inputFileName = dataname + "InputFile_TriCycle_" + epsilon + "_" + num_run + ".txt";
			dpc.generateCAGMInputFile(inputFileName);
		}
	}
	
	public void sampleSubgraph(String graphFile, String dataname, int num_nodes) throws IOException {
		graph g = new graph();
		g.readGraphOnly(graphFile);
	
		int rand_node = new Random().nextInt(g.getGraph().vertexSet().size());
		String v = new ArrayList<String>(g.getGraph().vertexSet()).get(rand_node);
		ArrayList<String> list = new ArrayList<String>();
		
		list.add(v);
		int curPos = 0;
		while(list.size() < num_nodes && curPos <= list.size()) {
			String currV = list.get(curPos);
			for(DefaultEdge e: g.getGraph().edgesOf(currV)) {
				String neighbour = g.getGraph().getEdgeSource(e).contentEquals(currV)?g.getGraph().getEdgeTarget(e):
					g.getGraph().getEdgeSource(e);
				if(!list.contains(neighbour)) {
					System.out.println(list.size() + " nodes added.");
					
					list.add(neighbour);
					if (list.size() >= num_nodes)
						break;
				}
			}
			curPos ++ ;
		}
		
		AsSubgraph<String, DefaultEdge> subG = new AsSubgraph<String, DefaultEdge>(g.getGraph(), new HashSet<String>(list));
	    
		System.out.println(subG.vertexSet().size());
		System.out.println(subG.edgeSet().size());
		
		System.out.println(g.getGraph().vertexSet().size());
		System.out.println(g.getGraph().edgeSet().size());
		
		
		FileWriter fout = new FileWriter(dataname + num_nodes+"_combined.txt");
		for (DefaultEdge e: subG.edgeSet()) {
			String src = subG.getEdgeSource(e);
			String target = subG.getEdgeTarget(e);
			fout.write(src + " " + target + "\n");
		}
		fout.close();
		
	}
}
	
