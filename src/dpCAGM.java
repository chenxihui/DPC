import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import org.jgrapht.GraphMetrics;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;


public class dpCAGM {
	private double[] abt_dist;
	private Partition p;
	//private ArrayList<HashMap<Integer, Double>> Theta_x;
	private double[][] Theta_x ;
	private graph g;
	private double[][][] Theta_F;
	private double epsilon;
	private HashMap<String, ArrayList<Integer>> attribute;
	private int num_abt;
	private String dataname;
	private double[][] theta_F_btw;
	//ratio version
	private double[][] theta_F_ratio;
	private double[] theta_btw_ratio; 
	private double[][] theta_f_g1_ratio;
	private double[] theta_f_btw_g1_ratio;
	private boolean privacy_considered; 
	private boolean DPCModel; 
	private boolean TriModel;
	private boolean DCSBMModel;
	
	public dpCAGM() {
		
	}
	public dpCAGM(String dataname, String graphFile, String abtFile, String pFile, double epsilon,
			boolean privacy_considered, boolean DPCModel, boolean TriModel, boolean DCSBMModel) throws IOException {
		
		System.out.println("Start reading the original graph ...");
		this.dataname = dataname;
		this.g = new graph();
		g.readGraphOnly(graphFile);
		
		System.out.println("Done. Strat reading the attributes ...");
		readAttribute(abtFile);

		System.out.println("Done. Start reading the partition ...");
		if (DPCModel||DCSBMModel ) {
	
			BufferedReader in = new BufferedReader(new FileReader(pFile));
			String pString = in.readLine();
			this.p = new Partition(pString);
			in.close();
		}else {
			String pString = "[[";
			for (String node:this.attribute.keySet())
				pString = pString + node + ",";
			pString = pString.substring(0, pString.length()-1) + "]]";
			this.p = new Partition(pString);
		}
		p.calNode2Com();
		
		System.out.println("Initial setting reading finished.");
		System.out.println("");

		this.epsilon = epsilon;
		this.privacy_considered = privacy_considered;
		this.DPCModel = DPCModel;
		this.DCSBMModel = DCSBMModel;
		this.TriModel = TriModel;
	}
	
	public void getAbtDist(double epsilon) {
		
		this.abt_dist = new double[this.p.partition.size()];
		Arrays.fill(abt_dist, 0.0);
		
		for (String v : this.p.node2com.keySet()) {
			int com = this.p.node2com.get(v);
			abt_dist[com] += 1;
		}
		
		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, 2.0/ epsilon);
		double total = 0;
		for(int com = 0; com < abt_dist.length; com ++) {
			abt_dist[com] += lng.nextLaplacian();
			total += abt_dist[com];
		}
		
		for(int com = 0; com < abt_dist.length; com ++) 
			abt_dist[com] /= total;
	}

	public void calThetaCX(double epsilon) {
		
		this.Theta_x =  new double[this.p.partition.size()][num_abt];
		
		Iterator<ArrayList<String>> iter = p.partition.iterator();

		int com_index = 0;


		while (iter.hasNext()) {
			
			ArrayList<String> com = iter.next();
			
			for(String node : com) {
				if (!this.attribute.containsKey(node)) {
					System.out.println("node not existing");
					Random rng = new Random();
					ArrayList<Integer> temp = new ArrayList<Integer>();
					for (int i =0; i< num_abt; i++)
						temp.add(rng.nextInt(2));
					this.attribute.put(node, temp);
				}
					
				ArrayList<Integer> abt_node = this.attribute.get(node);
				//int hashcode = abt_node.hashCode();

				
				for (int i = 0; i < this.num_abt; i++) {
					if (abt_node.get(i) == 0)
						this.Theta_x[com_index][i] ++;
				}
			}
						
			for ( int i = 0; i< num_abt; i++) {
				if (this.privacy_considered) {
					LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, 1.0/epsilon);
					this.Theta_x[com_index][i] += lng.nextLaplacian();
				}
				if (this.Theta_x[com_index][i]<0) 
					this.Theta_x[com_index][i] = 1;
				if (this.Theta_x[com_index][i] > com.size())
					this.Theta_x[com_index][i] = com.size();
				this.Theta_x[com_index][i] /= com.size(); 
				System.out.println(com_index +  " " + i + "th attribute:" + this.Theta_x[com_index][i]);
			}

			com_index ++;
		}
	}
	
	public double calProbThetaX(double prob, double lambda) {
		return prob * lambda + (1.0/this.g.getGraph().vertexSet().size()) * (1-lambda);
	}
	
	public void calThetaF(double epsilon, int k) {
		//this.Theta_F = new ArrayList<HashMap<Integer, Double>>();

		double [][][] theta_F = new double[this.p.partition.size()][num_abt][3]; 
		
		//this.abtPairMap = new HashMap<Integer, String>();

		LaplacianNoiseGenerator lgn = new LaplacianNoiseGenerator(0, 2*k/epsilon);

		//HashMap<Integer, Double> temp_theta_F_btw = new HashMap<Integer, Double>();
		this.theta_F_btw = new double[num_abt][3];
		
		Iterator<ArrayList<String>> iter = this.p.partition.iterator();
		
		int num_node_btw_com = 0;
		
		int com_index = 0;
		while(iter.hasNext()) {
			//HashMap<Integer, Double> temp_Theta_F = new HashMap<Integer, Double>();

			ArrayList<String> com = iter.next();
			for(String node : com) {
				if(!this.g.getGraph().containsVertex(node))
					continue;
				for (DefaultEdge e:this.g.getGraph().edgesOf(node)) {
					String vs = g.getGraph().getEdgeSource(e);
					if(vs.contentEquals(node))
						vs = g.getGraph().getEdgeTarget(e);

					ArrayList<Integer> abt_node = this.attribute.get(node);
					ArrayList<Integer> abt_vs = this.attribute.get(vs);

					if (!com.contains(vs)) {
						num_node_btw_com++;
						for (int i = 0; i < num_abt; i++) {
							if (abt_node.get(i) + abt_vs.get(i) == 0)
								theta_F_btw[i][0] += 1;
							else if (abt_node.get(i) + abt_vs.get(i) == 1)
								theta_F_btw[i][1]++;
							else
								theta_F_btw[i][2]++;
						}
					} else {
						for (int i = 0; i < num_abt; i++) {
							if (abt_node.get(i) + abt_vs.get(i) == 0)
								theta_F[com_index][i][0]++;
							else if (abt_node.get(i) + abt_vs.get(i) == 1)
								theta_F[com_index][i][1]++;
							else
								theta_F[com_index][i][2]++;
						}
					}
				}
			}

			for (int i=0; i< num_abt; i++) {

				double sum = 0.0 ;
				for (int j = 0; j < 3; j++) {
					double noisy = theta_F[com_index][i][j]/2 + lgn.nextLaplacian();
					if (noisy <= 0)
						noisy = 1;
					sum += noisy;
					theta_F[com_index][i][j] = noisy;
				}
				for (int j = 0; j<3; j++)
					theta_F[com_index][i][j] /= sum;
			}

			this.Theta_F = theta_F;
			com_index ++;
		}
		
		for(int i = 0; i < num_abt; i++) { 
			double sum = 0.0;
			for (int j = 0; j< 3; j++) {
				double noisy = theta_F_btw[i][j]/2 + lgn.nextLaplacian();
				if (noisy <= 0)
					noisy = 1;
				sum += noisy;
				theta_F_btw[i][j] = noisy;
			}
			
			for (int j = 0; j< 3; j++) {
				theta_F_btw[i][j] /= sum;
			}
		}
	}
	/*
	 * calcualte the differentially private ordered degree sequence
	 */
	public double calM_ij(double[] seq, int i, int j) {
		double m = 0.0;

		if(i>j) {
			return Double.MIN_VALUE;
		}
		/*if (i>=seq.length)
			i = seq.length -1;
		if (j>=seq.length)
			j = seq.length -1 ;*/
		for (int index = i; index <=j; index ++) {
			m += seq[index];
		}
		return m/(j-i+1);
	}
	
	public int [] calOrderedDegreeSequence(double[] seq) {
		int[] s = new int[seq.length];
		
		LinkedList<Integer> J = new LinkedList<Integer>();
		J.push(seq.length-1);
		for (int k = seq.length-1; k>=0; k-- ) {
			int j_a = k; 
			int j = J.getLast();
			
			while (!J.isEmpty() && calM_ij(seq, j_a +1, j) <= calM_ij(seq, k, j_a)) {
				j_a = j;
				J.pop();
				if (!J.isEmpty())
					j = J.getLast();
			}
			J.push(j_a);
		}
		
		int b = 0;
		while (!J.isEmpty()) {
			int j_a = J.pop();
			for(int k = b; k<=j_a; k++) {
				s[k] = (int)Math.round(calM_ij(seq, b, j_a));	
				//if (s[k] <= 0) 
					//s[k] = 1;
			}
			b = j_a +1;
		}
		
		return s;
	}
	
	
	public int [] calOrderedDegreeSequence_naive_version(double[] seq) {
		int [] s = new int [seq.length];
		if (seq.length==0)
			return s;
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
			s[k] =(int) Math.round(max);
		}
			
		return s;
	}
	
	public int[] graphicalProcess(int[] seq ) {
		
		/*
		 * seq is in a non-decreasing order
		 * s is also non-decreasing
		 * a is non-increasing
		 */
		int [] s = seq.clone();
		
		Arrays.sort(s);

		ArrayList<Integer> a = new ArrayList<Integer>() ;
		

		for (int i=seq.length-1; i>=0; i--)
			a.add(seq[i]);

		SimpleGraph<Integer, DefaultEdge> g = new SimpleGraph<Integer, DefaultEdge> (DefaultEdge.class);
		
		for (int i=0; i < seq.length; i++) {
			g.addVertex(i);
		}
		
		for (int i = 0; i < seq.length; i++) {
			int pos = 0;
			for (int j = i+1; j <= seq.length-1; j++) {
				if(a.get(j) != 0)
					pos ++ ;
			}
			int zi = s[s.length-i-1];
			int h = Math.min(zi, pos);
		 
			for (int j = i+1; j<=i+h; j++) {
				g.addEdge(i, j );
				a.set(j, a.get(j)-1);
			}
		}

		int[] s_final = new int[seq.length];
		for (int node:g.vertexSet()) {
			int d = g.degreeOf(node);
			s_final[node] = d;
		}
		return s_final;
	}
	
	public HashMap<Integer, ArrayList<NodeDegree>> calDegreeSequenceInBtw(double epsilon) {
	
		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, 2/epsilon);
		HashMap<Integer, ArrayList<NodeDegree>> degreeList_com = new HashMap<Integer, ArrayList<NodeDegree>> ();
		HashMap<String, Integer> degree_btw_global = new HashMap<String, Integer>();

		Iterator<ArrayList<String>> iter = p.partition.iterator();
		int com_index = 0;
	
		
		while(iter.hasNext()) {

			ArrayList<String> com = iter.next();
			if (com.size()==1 && com.get(0).isEmpty())
				 continue;
			if (com.size()==0)
					continue;
			ArrayList<NodeDegree> degreeList = new ArrayList<NodeDegree>(com.size());
		
			for (String v : com) {
				int degree_in = 0;
				int degree_btw = 0;

				if (!this.g.getGraph().containsVertex(v))
					continue;
				for (DefaultEdge e: this.g.getGraph().edgesOf(v)) {
					String vt = g.getGraph().getEdgeSource(e);
					if (vt.contentEquals(v))
						vt = g.getGraph().getEdgeTarget(e);
					if (com.contains(vt))
						degree_in ++;
					else
						degree_btw ++;
				}
				

				degreeList.add(new NodeDegree(v, degree_in, degree_btw, com_index));
			}
			
			
			
			//sort the list according to the degree_in values

			Collections.sort(degreeList);
			
			if (this.privacy_considered) {

				//add noise to the list
				double [] seq_in = new double[degreeList.size()];

				for (int i = 0; i< degreeList.size(); i++) {
					NodeDegree nd = degreeList.get(i);
					double noise = lng.nextLaplacian();
					nd.degree_in += noise;
					if (nd.degree_in < 0) 
						nd.degree_in =0;

					noise = lng.nextLaplacian();
					nd.degree_btw += noise;
					if(nd.degree_btw < 0)
						nd.degree_btw = 0;

					seq_in[i] = nd.degree_in;
					degree_btw_global.put(nd.node, (int)Math.round(nd.degree_btw));
				}

				double sum = 0 ;
				for (double d : seq_in)
					sum += d;

				int [] s = this.calOrderedDegreeSequence_naive_version(seq_in);

				sum = 0 ;
				for (double d : s)
					sum += d;


				for(int i = 0; i < degreeList.size(); i++) {
					NodeDegree nd = degreeList.get(i);
					nd.degree_in = s[i];
				}
			}
			
			degreeList_com.put(com_index, degreeList);
			com_index ++;
		}

		if (this.privacy_considered) {
			//post-processing to ensure the graphicalality of degree_btw
			String [] nodeList  = new String[degree_btw_global.keySet().size()];
			int [] degree_btw_list = new int[degree_btw_global.keySet().size()];

			int i = 0;
			for (String node : degree_btw_global.keySet()) {
				nodeList[i] = node;
				degree_btw_list[i] = degree_btw_global.get(node);
				i++;
			}

			//		int[] s = this.graphicalProcess(degree_btw_list);

			double sum = 0 ;
			for (double d : degree_btw_list)
				sum += d;


			for (int com_id : degreeList_com.keySet()) {
				ArrayList<NodeDegree> comND = degreeList_com.get(com_id);
				for (NodeDegree nd : comND) {
					nd.degree_btw = degree_btw_global.get(nd.node);
				}
			}
		}
		return degreeList_com;
	}
	
	public void generateCAGMInputFile(String outputFileName) throws IOException {
		

		System.out.println("Generating CAGM input file ...");

		HashMap<String, Double> epsilon_dist= SplitEpsilon();

		//calculate degree sequence 
		System.out.println("Calculating the degree sequence both intra- and inter-community ...");
		HashMap<Integer, ArrayList<NodeDegree>> dseq = this.calDegreeSequenceInBtw(epsilon_dist.get("theta_degree"));
		System.out.println("Degree sequence calculated.");
		System.out.println();

		//calculate triangle_in, triangle_btw

		System.out.println("Calculating the number of intra-community triangles ...");
		int tri_in = 0;
		int tri_btw = 0;
		if (this.privacy_considered && this.TriModel == false) {
			tri_in = LadderFunction.countTriangle_normalDP(this.g, epsilon_dist.get("theta_tria_in"), this.p);
			if(this.dataname.contentEquals("petster"))
				tri_in = LadderFunction.countTriangleDP(this.g, epsilon_dist.get("theta_tria_in"), this.p);


			System.out.println("    " + tri_in + " intra-community triangles calcualted.");

			System.out.println("  Calculating the number of inter-community triangles ...");
			tri_btw = LadderFunction.countTriangle_normalDP(this.g, epsilon_dist.get("theta_tria")) - tri_in;
			if(this.dataname.contentEquals("petster"))
				tri_btw =LadderFunction.countTriangleDP(this.g, epsilon_dist.get("theta_tria")) - tri_in;
			System.out.println("    " + tri_btw + " inter-community triangles calcualted.");

		}else if(!this.privacy_considered && this.DPCModel) {
			tri_in = LadderFunction.calNumIntraTriangles_v2(this.g, this.p);

			System.out.println("    " + tri_in + " intra-community triangles calcualted.");

			System.out.println("  Calculating the number of inter-community triangles ...");
			
			tri_btw = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph())- tri_in;

			System.out.println("    " + tri_btw + " inter-community triangles calcualted.");


		} else if (this.TriModel) {
			tri_in = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
			tri_btw = 0;
			System.out.println("    " + tri_in + " intra-community triangles calcualted.");
		}
	
		System.out.println("  Generating the output file.");
		FileWriter fout = new FileWriter(outputFileName);
		fout.write(tri_in + "\t" + tri_btw +"\n");
		for (int com_index : dseq.keySet()) {
			ArrayList<NodeDegree> ndlist = dseq.get(com_index);
			for (NodeDegree nd:ndlist) {
				String line = nd.node + "\t" + nd.com_index + "\t" + nd.degree_in + "\t" + nd.degree_btw +"\n";
				fout.append(line);
			}
		}
		fout.close();
		System.out.println("CAGM Input file generated successfully.");
	}
	
	public g_theta_F_value calculateThetaF(graph g1, HashMap<String, ArrayList<Integer>> abt_rnd, double [][][] theta_f_g1, double[][] theta_btw_g1){
		
		theta_f_g1 = new double[this.p.partition.size()][num_abt][3]; 
		
		theta_btw_g1 = new double[num_abt][3];
		
		Iterator<ArrayList<String>> iter = this.p.partition.iterator();
		
		int com_index = 0;

		while(iter.hasNext()) {

			ArrayList<String> com = iter.next();

			for(String node : com) {
				if(!g1.getGraph().containsVertex(node))
					continue;
				for (DefaultEdge e : g1.getGraph().edgesOf(node)) {
					String vs = g1.getGraph().getEdgeSource(e);

					if(vs.contentEquals(node))
						vs = g1.getGraph().getEdgeTarget(e);
					

					ArrayList<Integer> abt_node = abt_rnd.get(node);
					ArrayList<Integer> abt_vs = abt_rnd.get(vs);

					if (!com.contains(vs)) {
						for (int i = 0; i < num_abt; i++) {
							if (abt_node.get(i) + abt_vs.get(i) == 0)
								theta_btw_g1[i][0] += 1;
							else if (abt_node.get(i) + abt_vs.get(i) == 1)
								theta_btw_g1[i][1]++;
							else
								theta_btw_g1[i][2]++;
						}
					} else {
						for (int i = 0; i < num_abt; i++) {
							if (abt_node.get(i) + abt_vs.get(i) == 0)
								theta_f_g1[com_index][i][0]++;
							else if (abt_node.get(i) + abt_vs.get(i) == 1)
								theta_f_g1[com_index][i][1]++;
							else
								theta_f_g1[com_index][i][2]++;
						}
					}
				}
			}

			for (int i=0; i< num_abt; i++) {

				double sum = 0.0 ;
				for (int j = 0; j < 3; j++) 
					sum += theta_f_g1[com_index][i][j];
				
				for (int j = 0; j<3; j++) {
					theta_f_g1[com_index][i][j]= (theta_f_g1[com_index][i][j]*0.5+1)/(sum/2 + 1);
					theta_f_g1[com_index][i][j] = theta_f_g1[com_index][i][j] <0.001?0.001:theta_f_g1[com_index][i][j];
				}
				
				sum = 0.0 ;
				for (int j = 0; j < 3; j++) 
					sum += theta_f_g1[com_index][i][j];
				for (int j = 0; j<3; j++) 
					theta_f_g1[com_index][i][j]= theta_f_g1[com_index][i][j]/sum;

				    
			}
			com_index ++;
		}
		
		for(int i = 0; i < num_abt; i++) { 
			double sum = 0.0;
			for (int j = 0; j< 3; j++) 
				sum += theta_btw_g1[i][j];
			
			for (int j = 0; j< 3; j++) {
				theta_btw_g1[i][j] = (theta_btw_g1[i][j] * 0.5 +1)/ (sum/2 +1);
				theta_btw_g1[i][j] = theta_btw_g1[i][j]<0.001? 0.001:theta_btw_g1[i][j];
			}
		}
		return new g_theta_F_value(theta_f_g1, theta_btw_g1);
	}

	public void genDpCAGM(String randomGraphFile, int no) throws IOException {
		
		HashMap<String, Double> epsilon_dist= SplitEpsilon();
		
		//calculate Theta_C which is the community division
		Partition C = this.p;
	
		//calculate Theta_X which is this.Theta_x 
		System.out.println("Calculating Theta_xc ...");

		this.calThetaCX(epsilon_dist.get("theta_xc")); 	
		
		//Calculate Theta_F which is this.Theta_F
		System.out.println("Theta_xc calculated. Start calculating theta_f ...");
		this.calThetaF_v2(epsilon_dist.get("theta_f"), 100);
			
		System.out.println("Theta_F calculated. Start calculating theta_f_g1 ...");
		HashMap<String, ArrayList<Integer>> attribute_rnd = RandomAttributes();
	
		graph g1 = new graph();
		g1.readGraphOnly(randomGraphFile);
		
		/*
		double[][][] theta_f_g1 = null;
		double [][] theta_f_btw_g1 = null; 

		g_theta_F_value g_f_v = calculateThetaF(g1, attribute_rnd, theta_f_g1, theta_f_btw_g1);
		theta_f_g1= g_f_v.theta_f;
		theta_f_btw_g1 = g_f_v.theta_f_btw;
		*/

		calculateThetaF_ratio(g1, attribute_rnd);
	
		//ArrayList<HashMap<Integer, Double>> A = new ArrayList<HashMap<Integer, Double>>(); 

		double [] sup = new double[this.p.partition.size()+1];


		for (int com_index = 0 ; com_index<this.p.partition.size(); com_index ++) {
			double ms = 0.0;
			for (int i = 0; i < 11; i++) {
				System.out.println(this.theta_F_ratio[com_index][i] +  ": " +theta_F_ratio[com_index][i]);
				double r0 = this.theta_F_ratio[com_index][i]/theta_f_g1_ratio[com_index][i];

				if (ms < r0)
					ms = r0;
			}
			sup[com_index] = Math.min(ms, 25.0);
		}

		double ms = 0.0;
		for (int i = 0; i < 11; i++) {
			double r0 = this.theta_btw_ratio[i]/this.theta_f_btw_g1_ratio[i];

			if (ms < r0) ms = r0;
		}
		sup[this.p.partition.size()] = Math.min(ms, 25.0);

		//Store the attribute generated for the new graph
	
		FileWriter fabt;
		if (this.DPCModel) {
			String filedir = this.dataname + "_combined_experiment/";
			fabt = new FileWriter(filedir + this.dataname+"_attribute_" + this.epsilon + "_" +no+ ".txt");
		}else if (this.TriModel) {
			String filedir = this.dataname + "_TriCycle_combined_experiment/";
			fabt = new FileWriter(filedir + this.dataname+"_TriCycle_attribute_" + this.epsilon + "_" +no+ ".txt");
		}else { 
			String filedir = this.dataname + "_DCSBM_combined_experiment/";
			fabt = new FileWriter(filedir+this.dataname+"_DCSBM_attribute_" + this.epsilon + "_" +no+ ".txt");
		}			
		
		for (String node : attribute_rnd.keySet()) {
			String outStr= node;
			ArrayList<Integer> abt = attribute_rnd.get(node);
			for (int i: abt)
				outStr = outStr +" " + i;  
			fabt.write(outStr+"\n");
		}
		fabt.close();

		// Store theta_F to file	
		FileWriter fF ;
		if (this.DPCModel) {
			String filedir = this.dataname + "_combined_experiment/";
			fF = new FileWriter(filedir + this.dataname + "_F_" + this.epsilon + "_" + no + ".txt");
		}else if (this.TriModel) {
			String filedir = this.dataname + "_TriCycle_combined_experiment/";
			fF = new FileWriter(filedir+this.dataname + "_TriCycle_F_" + this.epsilon + "_" + no + ".txt");
		}else {
			String filedir = this.dataname + "_DCSBM_combined_experiment/";
			fF = new FileWriter(filedir + this.dataname + "_DCSBM_F_" + this.epsilon + "_" + no + ".txt");
		}	
		for (int com_index = 0; com_index< this.p.partition.size(); com_index ++ ) {
			String outStr = "" + com_index;
			for (int i = 0 ; i< 11; i++) {
				outStr = outStr +  " " + this.theta_F_ratio[com_index][i];
			}
			fF.write(outStr + "\n");
		}
		String outStr = "" + "-1";

		for (int i = 0 ; i< 11; i++) {
			outStr = outStr +  " " + this.theta_btw_ratio[i];
		}
		fF.write(outStr+"\n");

		fF.close();
		// output theta_f_g1 to file 

		if(this.DPCModel) {
			String filedir = this.dataname + "_combined_experiment/";
			fF = new FileWriter(filedir + this.dataname + "_F_g1_" + this.epsilon + "_" + no + ".txt");
		}else if (this.TriModel) {
			String filedir = this.dataname + "_TriCycle_combined_experiment/";
			fF = new FileWriter(filedir + this.dataname + "_TriCycle_F_g1_" + this.epsilon + "_" + no + ".txt");
		}else {
			String filedir = this.dataname + "_DCSBM_combined_experiment/";
			fF = new FileWriter(filedir + this.dataname + "_DCSBM_F_g1_" + this.epsilon + "_" + no + ".txt");
		}
		for (int com_index = 0; com_index< this.p.partition.size(); com_index ++ ) {
			outStr = "" + com_index;
			for (int i = 0 ; i< 11; i++) {
				outStr = outStr +  " " + theta_f_g1_ratio[com_index][i] ;
			}
			fF.write(outStr + "\n");
		}
		outStr = "" + "-1";
		for (int i = 0 ; i< 11; i++) {
			outStr = outStr +  " " + this.theta_f_btw_g1_ratio[i] ;
		}

		fF.write(outStr+"\n");

		fF.close();
		FileWriter fout;
		if (this.DPCModel) {
			String filedir = this.dataname + "_combined_experiment/";
			fout = new FileWriter(filedir + this.dataname+"_sup_" + this.epsilon + "_" + no + ".txt");
		}else if (this.TriModel){
			String filedir = this.dataname + "_TriCycle_combined_experiment/";
			fout = new FileWriter(filedir + this.dataname+"_TriCycle_sup_" + this.epsilon + "_" + no + ".txt");
		}else {
			String filedir = this.dataname + "_DCSBM_combined_experiment/";
			fout = new FileWriter(filedir + this.dataname+"_DCSBM_sup_" + this.epsilon + "_" + no + ".txt");
		}
			
		for (int com_index = 0; com_index <= this.p.partition.size() ; com_index ++ ) {
			if (com_index < this.p.partition.size())
				fout.write(com_index + " " + sup[com_index]+"\n");
			else
				fout.write(-1 + " " + sup[com_index]+"\n");
		}

		fout.close();
	}
	

	HashMap<String, ArrayList<Integer>> RandomAttributes(){
		
		HashMap<String, ArrayList<Integer>> sample_res = new HashMap<String, ArrayList<Integer>> ();
	
		Random rng = new Random();

		for(int com_index = 0; com_index < this.p.partition.size(); com_index ++) {
			ArrayList<String> com = this.p.partition.get(com_index);

			for (String node : com) {
				ArrayList<Integer> abt_node = new ArrayList<Integer>();

				for (int i = 0; i<num_abt; i++) {
					double r = rng.nextDouble();
					if(r < this.Theta_x[com_index][i])
						abt_node.add(0);
					else
						abt_node.add(1);
				}
				
				sample_res.put(node, abt_node);
			}
		}
		return sample_res;
	}
	
	HashMap<String, Double> SplitEpsilon(){
		HashMap<String, Double> dist = new HashMap<String, Double>();
		double e = this.epsilon;
		if(this.TriModel)
			e *= 2.0;
		dist.put("theta_c", 0.5 * e);
		dist.put("theta_f", 1.0 * e/6);
		dist.put("theta_xc", e/12);
		dist.put("theta_degree", e/12);
		dist.put("theta_tria_in", e/12);
		dist.put("theta_tria", e/12);
		
		return dist;
	}
	
	ArrayList<Integer> getArrayListFromString(String ptn){
		ArrayList<Integer> r = new ArrayList<Integer>();
		
		ptn = ptn.trim();
		ptn = ptn.substring(1, -1);
		String[] fds = ptn.split(", ");
		for (int i = 0; i<fds.length; i++)
			r.add(Integer.parseInt(fds[i]));
		return r;
	}
	public void readAttribute(String file) throws IOException{
		this.attribute = new HashMap<String, ArrayList<Integer>>();
		
		BufferedReader in = new BufferedReader(new FileReader(file));
		
		String line = "";

		while ((line = in.readLine()) != null) {
				
			String [] fds = line.split(" ");
			ArrayList<Integer> av = new ArrayList<Integer>(fds.length-1);
			
			for(int i =1; i<fds.length; i++) 
				av.add((int)Double.parseDouble(fds[i]));
		
			this.attribute.put(fds[0], av);
			this.num_abt = 50;
			switch(this.dataname) {
				case "facebook": this.num_abt = 50; break;
				case "petster": this.num_abt = 13; break;
				case "epinions": this.num_abt = 50;
			}
			//this.num_abt = 100;
			//this.num_abt = av.size();
		}
		in.close();
	}	
	
	public void calThetaF_v2(double epsilon, int k) {
		double[][] theta_f = new double[this.p.partition.size()][11];
		double [] theta_btw = new double[11];

		//HashMap<Integer, Double> temp_theta_F_btw = new HashMap<Integer, Double>();
		
		Iterator<ArrayList<String>> iter = this.p.partition.iterator();
		
		int num_node_btw_com = 0;
		
		int com_index = 0;

		while(iter.hasNext()) {
			//HashMap<Integer, Double> temp_Theta_F = new HashMap<Integer, Double>();

			ArrayList<String> com = iter.next();
			for(String node : com) {
				if(!this.g.getGraph().containsVertex(node))
					continue;
				for (DefaultEdge e:this.g.getGraph().edgesOf(node)) {
					String vs = g.getGraph().getEdgeSource(e);
					if(vs.contentEquals(node))
						vs = g.getGraph().getEdgeTarget(e);

					//System.out.println(node + " " + vs + " " + this.attribute.containsKey(node) + " " + this.attribute.containsKey(vs));
					
					if (!this.attribute.containsKey(node) || !this.attribute.containsKey(vs))
						continue;
					ArrayList<Integer> abt_node = this.attribute.get(node);
					ArrayList<Integer> abt_vs = this.attribute.get(vs);
					int ratio = (int) Math.floor(this.abt_distance(abt_node,abt_vs)/0.1);
					if (!com.contains(vs)) {
						num_node_btw_com++;
						theta_btw[ratio] += 1.0;
					}	
					else{
						theta_f[com_index][ratio] +=1;
					}
				}
			}

			double sum = 0.0 ;
			double noisy;
			for (int i=0; i< 11; i++) {
				if (this.privacy_considered) {
					LaplacianNoiseGenerator lgn = new LaplacianNoiseGenerator(0, 2*k/epsilon);

					noisy = theta_f[com_index][i]*1.0/2 + lgn.nextLaplacian();
				}else
					noisy = theta_f[com_index][i]*1.0/2;
				if (noisy <= 0)
					noisy = 1;
				sum += noisy;
				theta_f[com_index][i] = noisy;
			}
			for (int i = 0; i<11; i++)
					theta_f[com_index][i] /= sum;
			
			com_index ++;
		}
		
		this.theta_F_ratio = theta_f;
		this.theta_btw_ratio = theta_btw;
		double sum = 0.0;
		double noisy;
		for(int i = 0; i < 11; i++) { 
			if(this.privacy_considered) {
				LaplacianNoiseGenerator lgn = new LaplacianNoiseGenerator(0, 2*k/epsilon);
				noisy = theta_btw[i]/2 + lgn.nextLaplacian();
			}else
				noisy = theta_btw[i]/2;
			if (noisy <= 0)
				noisy = 1;
			sum += noisy;
			theta_btw[i] = noisy;
		}
			
		for (int i = 0; i< 11;  i++) {
				theta_btw[i] /= sum;
		}	
		
	}

	public void GenerateDCSBMInputFile(String inputFilename) throws IOException {

		Iterator<ArrayList<String>> iter = p.partition.iterator();
		FileWriter fout = new FileWriter(this.dataname+"InputFile_DCSBM.txt");
		
		int[] nDegreeCom = new int[p.partition.size()];
		
		int com_index = 0;
		while(iter.hasNext()) {
			int nDegree= 0;
			ArrayList<String> com = iter.next();
			for (String node:com) {
				if (g.getGraph().containsVertex(node)) {
					fout.append(node + " " + com_index + " " + g.getGraph().degreeOf(node) + "\n");
					nDegree += g.getGraph().degreeOf(node);
				}
			}

			nDegreeCom[com_index] = nDegree;
			com_index ++;
		}
		
		fout.close();
		
		int[][] nEdgeBtwCom = new int[p.partition.size()][p.partition.size()];
		
		for (String vs:g.getGraph().vertexSet()) {
			int com_vs = p.node2com.get(vs);
			for (DefaultEdge e : g.getGraph().edgesOf(vs)) {
				String vt = g.getGraph().getEdgeSource(e);
				if (vt.contentEquals(vs))
					vt = g.getGraph().getEdgeTarget(e);
				int com_vt = p.node2com.get(vt);
				nEdgeBtwCom[com_vs][com_vt] ++;
			}
		}
						
		fout = new FileWriter(this.dataname+"InputNumberEdges_DCSBM.txt");
		for (com_index = 0; com_index <p.partition.size(); com_index++ ) {
			String output = com_index + " " + nDegreeCom[com_index] ;
			for(int j =0 ; j <p.partition.size(); j++ ) {
				output = output + " " + nEdgeBtwCom[com_index][j]; 
			}
			output += "\n";
			fout.write(output);
		}
		fout.close();
	}

    double abt_distance(ArrayList<Integer> abt_vs, ArrayList<Integer> abt_vt) {
    	double distance = 0 ;

    
    	double prod = 0;
    	double ds = 0 ;
    	double dt = 0; 
    	
    	for (int i = 0; i<this.num_abt; i++ ) {
    	   prod += (abt_vs.get(i) * abt_vt.get(i));
    	   ds += Math.pow(abt_vs.get(i), 2);
    	   dt += Math.pow(abt_vt.get(i), 2);
    	}
    	distance = prod *1.0 /(Math.sqrt(ds) * Math.sqrt(dt));

    	return distance;
    }
    public void calculateThetaF_ratio(graph g1, HashMap<String, ArrayList<Integer>> abt_rnd ){
    	
		this.theta_f_g1_ratio = new double[this.p.partition.size()][11]; 
		
		this.theta_f_btw_g1_ratio = new double[11];
		
		Iterator<ArrayList<String>> iter = this.p.partition.iterator();
		
		int com_index = 0;

		while(iter.hasNext()) {

			ArrayList<String> com = iter.next();

			for(String node : com) {
				if(!g1.getGraph().containsVertex(node))
					continue;
				for (DefaultEdge e : g1.getGraph().edgesOf(node)) {
					String vs = g1.getGraph().getEdgeSource(e);

					if(vs.contentEquals(node))
						vs = g1.getGraph().getEdgeTarget(e);
					

					ArrayList<Integer> abt_node = abt_rnd.get(node);
					ArrayList<Integer> abt_vs = abt_rnd.get(vs);

					int ratio = (int)Math.floor(this.abt_distance(abt_node, abt_vs)/0.1);

					if (!com.contains(vs)) {
						theta_f_btw_g1_ratio[ratio] += 1;
							
					} else {
						theta_f_g1_ratio[com_index][ratio] += 1;
					}
					
				}
			}

			double sum = 0.0 ;
			for (int i=0; i< 11; i++) 
				sum += theta_f_g1_ratio[com_index][i];
				
			for (int i = 0; i<11; i++) {
					theta_f_g1_ratio[com_index][i] = (theta_f_g1_ratio[com_index][i] *0.5 +1)/(sum/2 + 1);
			}
				
			com_index ++;
		}
		
		double sum = 0.0;
		for(int i = 0; i < 11; i++)  
			sum += theta_f_btw_g1_ratio[i];
			
		for (int i = 0; i< 11; i++) {
				theta_f_btw_g1_ratio[i] = (theta_f_btw_g1_ratio[i] * 0.5 +1)/ (sum/2 +1);
		}
	}

    public void GenerateInputFileDCSBM(String inputFileName, String inputEdgeFile) throws IOException {
    	System.out.println("Generating DCSBM input file ...");

		HashMap<String, Double> epsilon_dist= SplitEpsilon();

		//calculate degree sequence 
		System.out.println("Calculating the degree sequence both intra- and inter-community ...");
		HashMap<Integer, ArrayList<NodeDegree>> dseq = this.calDegreeSequenceDCSBM(epsilon_dist.get("theta_degree"));

		System.out.println("Degree sequence calculated.");
		System.out.println();

    	
		FileWriter fout = new FileWriter(inputFileName);
		
		int[] nDegreeCom = new int[p.partition.size()];
		
		for (int com_index :dseq.keySet()){
			int nDegree= 0;
			ArrayList<NodeDegree> com = dseq.get(com_index);
			
			for (NodeDegree nd:com) {
				String node = nd.node;
				int degree = nd.degree_in;
				fout.append(node + " " + com_index + " " + degree  + "\n");
					nDegree += degree; 
				}

			nDegreeCom[com_index] = nDegree;
			com_index ++;
		}
		
		fout.close();
		
		int[][] nEdgeBtwCom = new int[p.partition.size()][p.partition.size()];
		
		for (String vs:g.getGraph().vertexSet()) {
			int com_vs = p.node2com.get(vs);
			for (DefaultEdge e : g.getGraph().edgesOf(vs)) {
				String vt = g.getGraph().getEdgeSource(e);
				if (vt.contentEquals(vs))
					vt = g.getGraph().getEdgeTarget(e);
				int com_vt = p.node2com.get(vt);
				nEdgeBtwCom[com_vs][com_vt] ++;
			}
		}
		
		for (int com_index = 0; com_index <p.partition.size(); com_index++ ) {
			for(int j =0 ; j <= com_index; j++ ) {
				LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, 2/(this.epsilon/6));
				nEdgeBtwCom[com_index][j] += (int)lng.nextLaplacian(); 
				
				if (com_index != j)
					nEdgeBtwCom[j][com_index] = nEdgeBtwCom[com_index][j];
			}
		}
		
						
		//fout = new FileWriter(this.dataname+"InputNumberEdges_DCSBM.txt");
		fout = new FileWriter(inputEdgeFile);
		for (int com_index = 0; com_index <p.partition.size(); com_index++ ) {
			String output = com_index + " " + nDegreeCom[com_index] ;
			for(int j =0 ; j <p.partition.size(); j++ ) {
				output = output + " " + nEdgeBtwCom[com_index][j]; 
			}
			output += "\n";
			fout.write(output);
		}
		fout.close();
	}
    public HashMap<Integer, ArrayList<NodeDegree>> calDegreeSequenceDCSBM(double epsilon) {
    	
		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, 2/epsilon);
		HashMap<Integer, ArrayList<NodeDegree>> degreeList_com = new HashMap<Integer, ArrayList<NodeDegree>> ();

		Iterator<ArrayList<String>> iter = p.partition.iterator();
		int com_index = 0;
	
		int total_degree = 0;
		
		while(iter.hasNext()) {

			ArrayList<String> com = iter.next();
			ArrayList<NodeDegree> degreeList = new ArrayList<NodeDegree>(com.size());
		
			for (String v : com) {

				if (!this.g.getGraph().containsVertex(v))
					continue;
				int degree = this.g.getGraph().degreeOf(v);

				total_degree += degree;
				degreeList.add(new NodeDegree(v, degree,0, com_index));
			}
			
			if (degreeList.size() <= 0)
				continue;
			//sort the list according to the degree_in values

			Collections.sort(degreeList);
			
			if (this.privacy_considered) {

				//add noise to the list
				double [] seq_in = new double[degreeList.size()];

				for (int i = 0; i< degreeList.size(); i++) {
					NodeDegree nd = degreeList.get(i);
					double noise = lng.nextLaplacian();
					//System.out.println(noise);
					nd.degree_in += noise;
					if (nd.degree_in < 0) 
						nd.degree_in =0;

					seq_in[i] = nd.degree_in;
				}

				double sum = 0 ;
				for (double d : seq_in)
					sum += d;

				int [] s = this.calOrderedDegreeSequence_naive_version(seq_in);

				sum = 0 ;
				for (double d : s)
					sum += d;


				for(int i = 0; i < degreeList.size(); i++) {
					NodeDegree nd = degreeList.get(i);
					nd.degree_in = s[i];
				}
			}
			
			degreeList_com.put(com_index, degreeList);
			com_index ++;
		}

		
		return degreeList_com;
	}
	
}

class g_theta_F_value {
	double[][][] theta_f;
	double [][] theta_f_btw;
	g_theta_F_value(double[][][] f, double[][] btw){
	  this.theta_f = f;
	  this.theta_f_btw = btw;
	}
}

class NodeDegree implements Comparable<NodeDegree>{
	String node;
	int degree_in;
	int degree_btw;
	int com_index;
	int degree;
	
	NodeDegree(){
		node = "";
		degree_in = 0 ;
		degree_btw = 0;
	}
	NodeDegree(String v, int degree, int com_index){
		this.degree = degree; 
		this.com_index = com_index;
		this.node = v;
	}
	
	NodeDegree(String node, int degree_in, int degree_btw, int com_index){
		this.node = node; 
		this.degree_btw = degree_btw;
		this.degree_in = degree_in;
		this.com_index = com_index;
	}
	public int compareTo(NodeDegree nd) {
		return (this.degree_in < nd.degree_in? -1:(this.degree_in == nd.degree_in? 0:1));
	}
}

