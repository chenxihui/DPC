import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import org.jgrapht.graph.DefaultEdge;

public class DPDCSBM_Generator {
	
	private String dataname;
	Partition p;
	graph g;
	HashMap<String, int[]> nodeAbtMap;
	HashMap<Integer, double[]> theta_F;
	HashMap<Integer, double[]> theta_F_g1;
	HashMap<Integer, Double> sup;
	HashMap<String, Integer> comMap;
	
	public DPDCSBM_Generator() {
		
	}
	
	public DPDCSBM_Generator(String dataname, String abtFile, String Ffile, 
			String Fg1File, String supFile) throws IOException {
		nodeAbtMap = new HashMap<String, int[]>();
		theta_F = new HashMap<Integer, double[]>();
		theta_F_g1 = new HashMap<Integer, double[]>();
		sup= new HashMap<Integer, Double>() ;
		this.readAttribute(abtFile, Ffile, Fg1File, supFile);
		this.dataname = dataname;
	}
	
	void readAttribute(String abtFile, String Ffile, String Fg1File, String supFile) throws IOException {
		/*
		 * read attributes from abtFile
		 */
		BufferedReader fin = new BufferedReader(new FileReader(abtFile));
		String line = "";
		while ((line = fin.readLine())!=null) {
			String[] fds = line.split(" ");
			int [] abt = new int[fds.length-1];

			for (int i = 1; i< fds.length; i++)
				abt[i-1] = Integer.parseInt(fds[i]);
			this.nodeAbtMap.put(fds[0], abt);
		}
		fin.close();
	
		/*
		 * read Theta_F of the originla graph
		 */
		
		fin = new BufferedReader(new FileReader(Ffile));
		while ((line = fin.readLine())!= null) {
			String[] fds = line.split(" ");
			int com_index = Integer.parseInt(fds[0]);
			double [] f = new double[fds.length-1];
			for (int i = 1; i<fds.length; i++) {
				f[i-1] = Double.parseDouble(fds[i]);
			}
			this.theta_F.put(com_index, f);
				
		}
		fin.close();
		
		/*
		 * read theta_f_g1
		 */
		
		fin = new BufferedReader(new FileReader(Fg1File));
		while ((line = fin.readLine())!= null) {
			String[] fds = line.split(" ");
			int com_index = Integer.parseInt(fds[0]);
			double [] f = new double[fds.length-1];
			for (int i = 1; i<fds.length; i++) {
				f[i-1] = Double.parseDouble(fds[i]);
			}
			this.theta_F_g1.put(com_index, f);
		}
		
		fin.close();
		
		/*
		 * read sup file
		 */
		
		fin = new BufferedReader(new FileReader(supFile));
		sup = new HashMap<Integer, Double>();

		while((line = fin.readLine())!=null) {
			String[] fds = line.split(" ");
			int com_index = Integer.parseInt(fds[0]);
			double sup_v = Double.parseDouble(fds[1]);
			
			sup.put(com_index, sup_v);
		}
	}

	

	public void GenerateInputFile() throws IOException {

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

	public void generateGraph(String inputFile, String numEdgeFile, String fileOut) throws NumberFormatException, IOException {
		HashMap<String, Integer> degreeMap = new HashMap<String, Integer> ();
		comMap = new HashMap<String, Integer>();

		BufferedReader fin = new BufferedReader(new FileReader(inputFile));
		String line;
		int total_degree = 0;
		
		System.out.println("Reading degree input file ...");
		
		while ((line = fin.readLine()) != null) {
			String[] fds = line.split(" ");
			int com_index = Integer.parseInt(fds[1]);
			String node = fds[0];
			int degree = Integer.parseInt(fds[2]);
			total_degree += degree;
			
			degreeMap.put(node, degree);
			comMap.put(node, com_index);
		}

		fin.close();
		
		System.out.println("Finished Reading degree input file. Starting reading edge file ...");
		fin = new BufferedReader(new FileReader(numEdgeFile));

		int [] degreeCom = null ;
	    int [][] edgeBtwCom = null;

	    while((line = fin.readLine())!=null) {
	    	String[] fds = line.split(" ");
	    	int com_index = Integer.parseInt(fds[0]);
	    	int num_com = fds.length - 2;
	    	
	    	if (degreeCom == null) {
	    		degreeCom = new int[num_com];
	    		edgeBtwCom = new int[num_com][num_com];
	    	}
	    	degreeCom[com_index] = Integer.parseInt(fds[1]);
	    	
	    	for (int i =2; i<fds.length; i++){
	    		edgeBtwCom[com_index][i-2] = Integer.parseInt(fds[i]);
	    	}
	    }
	    
	    System.out.println("Starting calculating degree probability distribution ...");
	    HashMap<String, Double> probDist = new HashMap<String, Double>();
	    for (String v : degreeMap.keySet()) {
	    	probDist.put(v, degreeMap.get(v)*1.0/total_degree);
	    }
	    
	    
	    System.out.println("Adding nodes to the graph ...");
	    graph g = new graph();
	    
	    for (String v : probDist.keySet())
	    	if(!g.getGraph().containsVertex(v))
	    		g.getGraph().addVertex(v);
	    
	    int m = 0;
	    int num_edge = total_degree/2;
	   
	    System.out.println("Start generating edges ...");
	    while (m < total_degree/2) {
	    	String v1 =  sampleVertex(probDist);
	    	String v2 = sampleVertex(probDist);
	    	if(!comMap.containsKey(v1) || !comMap.containsKey(v2))
	    		continue;
	        if (v1.contentEquals(v2))
	        	continue;

	    	double p = (2.0 * num_edge / degreeCom[comMap.get(v1)])
	    		* (edgeBtwCom[comMap.get(v1)][comMap.get(v2)]*1.0/degreeCom[comMap.get(v2)]);
	    	
	    	//System.out.println(v1 + "  " + v2 + " " + p);
	    	double r = new Random().nextDouble();
	    	if (r<=p) {
	    		if (!g.getGraph().containsEdge(v1, v2) && !v1.contentEquals(v2)) {
					double Af = this.calculateAf(v1, v2);
					double u = new Random().nextDouble();
					if (u<Af) {
						g.getGraph().addEdge(v1, v2);
						m ++;
						System.out.println(m + "/" + num_edge + ": edge ("+ 
			    					v1+","+v2 +") added with prob "+ p+".");
					}
	    		}
	    	}
	    }
	    
	   g.printGraph(fileOut);
		
		
	}
	
	public double calculateAf(String v1, String v2) {
		int com_index = -1;

		int num_abt = 50;
		if (this.dataname.contentEquals("epinions"))
			num_abt = 50;
		else if (this.dataname.contentEquals("facebook"))
			num_abt = 50;
		else if (this.dataname.contentEquals("petster"))
			num_abt = 13;
		
		if (this.comMap.get(v1) == this.comMap.get(v2))
			com_index = this.comMap.get(v1);
		
		double supc = this.sup.get(com_index);
		int[] abt_v1 = this.nodeAbtMap.get(v1);
		int[] abt_v2 = this.nodeAbtMap.get(v2);
	
		double prod_v1 = 0, prod_v2 = 0;
		double prod_cross = 0;
		for (int i = 0; i< num_abt; i++) {
			prod_v1 += (abt_v1[i] * abt_v1[i]);
			prod_v2 += (abt_v2[i] * abt_v2[i]);
			prod_cross += (abt_v1[i] * abt_v2[i]);
		}
		double division = Math.sqrt(prod_v1) * Math.sqrt(prod_v2);
	
		double distance = 0;
		if (division - 0.0 <0.00000001)
			distance = 0;
		else {
			distance = prod_cross*1.0/division;
		}
		
		int ratio =(int)Math.floor(distance/0.1);
		
		double[] F = this.theta_F.get(com_index);
		double f = F[ratio];
		F = this.theta_F_g1.get(com_index);
		double fg1 = F[ratio];
		double R = f/fg1;
		return R/supc;
	}
	public String sampleVertex(HashMap<String, Double> probDist) {
		Random rng = new Random();
		double r = rng.nextDouble();
		double sum = 0;
		String sampledV = "";
		for(String v : probDist.keySet()) {
			sum += probDist.get(v);
			
			if (sum>r) {
				sampledV = v;
				break;
			}
		}
		return sampledV;
	}
	public void generateGraph_not_used(String inputFile, String numEdgeFile) throws NumberFormatException, IOException {
		HashMap<Integer, HashMap<String ,Integer>> comNodeDegree = new HashMap<Integer, HashMap<String, Integer>>();
		
		BufferedReader fin = new BufferedReader(new FileReader(inputFile));
		String line;
		while ((line = fin.readLine()) != null) {
			String[] fds = line.split(" ");
			int com_index = Integer.parseInt(fds[1]);
			String node = fds[0];
			int degree = Integer.parseInt(fds[2]);
			
			if (!comNodeDegree.containsKey(com_index)) {
				comNodeDegree.put(com_index, new HashMap<String, Integer>());
			}
			comNodeDegree.get(com_index).put(node, degree);
		}
	}
	public void GenerateInputFile(double epsilon) throws IOException {

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
}
