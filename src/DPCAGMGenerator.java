import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.jgrapht.GraphMetrics;
import org.jgrapht.alg.connectivity.ConnectivityInspector;
import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultEdge;

public class DPCAGMGenerator {
	HashMap<String, Integer> indegreeMap ;
	HashMap<String, Integer> btwdegreeMap ;
	HashMap<String, Integer> comMap;
	HashMap<Integer, ArrayList<String>> comList;
	HashMap<Integer, Integer> degreeInComMap;
	HashMap<Integer, Integer> degreeBtwComMap;
	LinkedList<edge> intraEdgeList;
	LinkedList<edge> btwEdgeList;
	
	int in_triangle;
	int btw_triangle; 
	graph g; 
	
	String dataname;

	HashMap<Integer, Double> pi_btwDegC; 
	HashMap<Integer, HashMap<String, Double>> pi_btwDegInC; 
	HashMap<Integer, HashMap<String, Double>> pi_inDegInC; 
	
	HashMap<String, int[]> nodeAbtMap;
	HashMap<Integer, double[]> theta_F;
	HashMap<Integer, double[]> theta_F_g1;
	HashMap<Integer, Double> sup;
	
	public DPCAGMGenerator(String dataname, String inputFile, String abtFile, String Ffile, 
			String Fg1File, String supFile) throws IOException {
		readInputFile(inputFile);
		this.degreeInComMap = new HashMap<Integer, Integer>() ;
		this.degreeBtwComMap = new HashMap<Integer, Integer>() ;
		this.intraEdgeList = new LinkedList<edge>();
		this.btwEdgeList = new LinkedList<edge>();
		
		for (int com_index: comList.keySet()) {
			int in_degree = 0;
			int btw_degree = 0;
			
			for (String v : comList.get(com_index)) {
				in_degree += this.indegreeMap.get(v);
				btw_degree += this.btwdegreeMap.get(v);
			}
			this.degreeInComMap.put(com_index, in_degree);
			this.degreeBtwComMap.put(com_index, btw_degree);
		}
		this.dataname = dataname;
		
		nodeAbtMap = new HashMap<String, int[]>();
		theta_F = new HashMap<Integer, double[]>();
		theta_F_g1 = new HashMap<Integer, double[]>();
		sup= new HashMap<Integer, Double>() ;
		this.readAttribute(abtFile, Ffile, Fg1File, supFile);
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

	public double calculateAf(String v1, String v2) {
		int num_abt = 50;
		if (this.dataname.contains("epinions"))
			num_abt = 50;
		else if (this.dataname.contentEquals("facebook"))
			num_abt = 50;
		else if (this.dataname.contentEquals("petster"))
			num_abt = 13;

		int com_index = -1;
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
	public void readInputFile(String inputFile) throws IOException {
		
		this.g = new graph();
		this.indegreeMap = new HashMap<String, Integer> ();
		this.btwdegreeMap = new HashMap<String, Integer> ();
		this.comMap = new HashMap<String, Integer>();

		this.comList = new HashMap<Integer, ArrayList<String>>();
		BufferedReader fin = new BufferedReader(new FileReader(inputFile));
		String line = fin.readLine();
		String[] fds = line.split("\t");
	    
		this.in_triangle = (int)Double.parseDouble(fds[0]);
		this.btw_triangle = (int) Double.parseDouble(fds[1]);
		
		while ((line= fin.readLine())!= null) {
			fds = line.split("\t");
			
			String node = fds[0];
			
			if (!g.getGraph().containsVertex(node))
				g.getGraph().addVertex(node);

			int com_index = Integer.parseInt(fds[1]);
			int in_degree = Integer.parseInt(fds[2]);
			int btw_degree = Integer.parseInt(fds[3]);
			
			indegreeMap.put(node, in_degree);
			btwdegreeMap.put(node, btw_degree);
			comMap.put(node, com_index);
			
			if (!comList.keySet().contains(com_index)) {
				ArrayList<String> vlist = new ArrayList<String>();
				vlist.add(node);
				comList.put(com_index, vlist);
			}else
				comList.get(com_index).add(node);
		}
	}
	public void generate_graph_no_triangles(String file_wt_tri) throws IOException {
		
		
		this.pi_inDegInC = new HashMap<Integer, HashMap<String, Double>>();
		for(int com_index : this.comList.keySet()) {
			ArrayList<String> com = this.comList.get(com_index);
			int m_in_c = Math.round(this.degreeInComMap.get(com_index)/2);
			
			int m_gen = 0;
			HashMap<String, Double> pi_in_degree = new HashMap<String, Double>(); 
			for (String v : com) {
				pi_in_degree.put(v, this.indegreeMap.get(v)*1.0/this.degreeInComMap.get(com_index));
			}
			
			this.pi_inDegInC.put(com_index, pi_in_degree);
			//int com_size = com.size();
			int n_try = 0;
			while (m_gen < m_in_c && n_try < com.size()) {
				String v1 = sampleVertex(pi_in_degree);
				int num_tries = 0;
				while (num_tries <1000) {
					String v2 = sampleVertex(pi_in_degree);
					if(!v1.contentEquals(v2) && !this.g.getGraph().containsEdge(v1, v2)) {
						double Af = this.calculateAf(v1, v2);
						double u = new Random().nextDouble();
						if (u<Af) {
							this.g.getGraph().addEdge(v1,  v2);
							m_gen ++;
							num_tries = 0;
							this.intraEdgeList.add(new edge(v1, v2));
							System.out.println("Community " + com_index + ": intra-edge " + m_gen + " of " + m_in_c
								+ ": (" + v1 + "," + v2+") is added.1");
							break;
						}
					}
					num_tries ++;
				}
				if(num_tries == 1000)
					n_try ++;
				else
					n_try = 0;
				//System.out.println(n_try);
			}
		}
		
		int sum_btw_degree = 0 ;
		for (int com_index: this.degreeBtwComMap.keySet())
			sum_btw_degree += this.degreeBtwComMap.get(com_index);
	
		int m_btw =Math.round( sum_btw_degree/2);
		
		this.pi_btwDegC = new HashMap<Integer, Double> ();
		this.pi_btwDegInC = new HashMap<Integer, HashMap<String, Double>> ();
		
		for (int com_index: this.degreeBtwComMap.keySet()) {
			pi_btwDegC.put(com_index, this.degreeBtwComMap.get(com_index)*1.0/sum_btw_degree);
			HashMap<String, Double> pi_btw_nodes_inCom = new HashMap<String, Double>();
			
			for (String v : this.comList.get(com_index)) {
				pi_btw_nodes_inCom.put(v, this.btwdegreeMap.get(v)*1.0/this.degreeBtwComMap.get(com_index));
			}
			pi_btwDegInC.put(com_index, pi_btw_nodes_inCom);
		}
		
		int m_gen = 0;
		
		int n_tries = 0;

		while (m_gen < m_btw && n_tries < 2* m_btw) {
			//int c1 = this.sampleCom(pi_btwDegC);
			//int c2 = this.sampleCom(pi_btwDegC);
			twoCom temp = this.sampleTwoComs(pi_btwDegC);
			int c1 = temp.c1;
			int c2 = temp.c2;
			
			if(c1 == c2)
				continue;
			
			String v1 = sampleVertex (pi_btwDegInC.get(c1));
			String v2 = sampleVertex (pi_btwDegInC.get(c2));
			
			if (!g.getGraph().containsEdge(v1, v2) && !v1.contentEquals(v2)) {
				double Af = this.calculateAf(v1, v2);
				double u = new Random().nextDouble();
						
				if (u < Af) {
					this.g.getGraph().addEdge(v1,  v2);
					m_gen ++;
					this.btwEdgeList.add(new edge(v1, v2));
					System.out.println("inter-edge " + m_gen + " of " + m_btw + " : (" + v1 + "," + v2 + ") is added.");
					n_tries = 0;
				}
			}else
				n_tries ++;
			
		}
		if (file_wt_tri.isEmpty())
				this.g.printGraph(this.dataname + "_graph_without_triangles.txt");
			else
				g.printGraph(file_wt_tri);

	}

	public int sampleCom(HashMap<Integer, Double> probDist) {
		Random rng = new Random();
		double r = rng.nextDouble();
		double sum = 0;
		int sampledV =0 ;
		for(int com_index : probDist.keySet()) {
			sum += probDist.get(com_index);
			
			if (sum>r) {
				sampledV = com_index;
				break;
			}
		}
		return sampledV;
	}

	public twoCom sampleTwoComs(HashMap<Integer, Double> probDist) {
		Random rng = new Random();
		double r = rng.nextDouble();
		double sum = 0;
		int sampledC1 =0 ;
		int sampledC2 =0 ;
		for(int com_index : probDist.keySet()) {
			sum += probDist.get(com_index);
			
			if (sum>r) {
				sampledC1 = com_index;
				break;
			}
		}
		
		r = new Random().nextDouble() * (1.0 - probDist.get(sampledC1));
		for(int com_index : probDist.keySet()) {
			if (com_index != sampledC1)
				sum += probDist.get(com_index);
			
			if (sum>r) {
				sampledC2 = com_index;
				break;
			}
		}
		
		return new twoCom(sampledC1, sampledC2);
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
	
	String sampleNeighbourInCom(String v) {
		String sampledV = "";
		if (!this.comMap.containsKey(v)) {
			System.out.println("error!");
			return "";
		}
		int com_v = this.comMap.get(v);
		ArrayList<String> neiList = new ArrayList<String>();
		for(DefaultEdge e : this.g.getGraph().edgesOf(v)) {
			String nei = this.g.getGraph().getEdgeSource(e);
			if (nei.contentEquals(v))
				nei = this.g.getGraph().getEdgeTarget(e);
			if(this.comMap.get(nei) == com_v)
				neiList.add(nei);
		}
		if (neiList.isEmpty())
			return "";
		else {
			int i = new Random().nextInt(neiList.size());
			return neiList.get(i);
		}
	}
	String sampleNeighbourBtwCom(String v) {
		String sampledV = "";
		int com_v = this.comMap.get(v);
		ArrayList<String> neiList = new ArrayList<String>();
		for(DefaultEdge e : this.g.getGraph().edgesOf(v)) {
			String nei = this.g.getGraph().getEdgeSource(e);
			if (nei.contentEquals(v))
				nei = this.g.getGraph().getEdgeTarget(e);
			if(this.comMap.get(nei) != com_v)
				neiList.add(nei);
		}
		if (neiList.isEmpty())
			return "";
		else {
			int i = new Random().nextInt(neiList.size());
			return neiList.get(i);
		}
	}
	void generate_graph_triangle(String file_w) throws IOException {
		//this.generate_graph_no_triangles(file_wt);
		int num_in_triangles = this.calNumIntraTriangles_v2();
		System.out.println("The number of intra_triangles is " + num_in_triangles);
	
		while(num_in_triangles<this.in_triangle) {
			int c = new Random().nextInt(this.comList.keySet().size());
			String v1 = this.sampleVertex(this.pi_inDegInC.get(c));
			String v2 = this.sampleNeighbourInCom(v1);
			if (v2.isEmpty()) 
				continue;
			String v3 = this.sampleNeighbourInCom(v2);
			if (v3.isEmpty()) 
				continue;
			
			if (!v1.contentEquals(v3) && !this.g.getGraph().containsEdge(v1, v3)) {
				int n_tries = 0;
				edge e = null;
				String vs = "", vt = "";

				while (n_tries < 5 * this.intraEdgeList.size()) {
					e = this.intraEdgeList.poll();
					vs = e.v1;
					vt = e.v2;
					if (g.getGraph().edgesOf(vs).size() == 1 ||
							g.getGraph().edgesOf(vt).size() ==1 || 
							this.num_common_neighbours(vs, vt) ==0) {
						this.intraEdgeList.add(e);
						n_tries ++;
					}else
						break;
				}
				if (n_tries == 5 * this.intraEdgeList.size()) 
					continue;
				
				int num_tri_prev = this.num_common_neighbours_in(vs, vt);
				this.g.getGraph().removeEdge(vs, vt);
				
				int num_tri_after = this.num_common_neighbours_in(v1,  v3);
				
				if (num_tri_prev < num_tri_after) {
					double Af = this.calculateAf(v1, v3);
					double u = new Random().nextDouble();
					if (u<Af) {
						this.g.getGraph().addEdge(v1, v3);
						this.intraEdgeList.add(new edge(v1, v3));
						num_in_triangles += (num_tri_after - num_tri_prev);
						System.out.printf("(%s, %s) n_prev =%d, n_after=%d, intra_triangle=%d of %d \n", 
								v1, v3, num_tri_prev, num_tri_after, num_in_triangles, this.in_triangle);
					}else {
						this.g.getGraph().addEdge(vs,  vt);
						this.intraEdgeList.add(e);
					}
				}else {
					this.g.getGraph().addEdge(vs,  vt);
					this.intraEdgeList.add(e);
				}
			}
			
		}
		
		int num_btw_triangles = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph())
				 - num_in_triangles;
		
		System.out.println("Start increasing inter-community triangles ...");
		
		while (num_btw_triangles < this.btw_triangle) {
			int c = this.sampleCom(this.pi_btwDegC);
			String v1 = this.sampleVertex(this.pi_btwDegInC.get(c));
			String v2 = this.sampleNeighbourBtwCom(v1);
			
			if(v2.isEmpty()) continue;
			
			String v3 = this.sampleNeighbourInCom(v2);
			
			edge e = this.btwEdgeList.poll();
			String vs = e.v1, vt = e.v2;
			
			int num_tri_prev = this.num_common_neighbours(vs, vt);
			this.g.getGraph().removeEdge(vs, vt);
			
			int num_tri_after = this.num_common_neighbours(v1, v3);
			
			if (num_tri_prev < num_tri_after && !this.g.getGraph().containsEdge(v1, v3)) {
				double Af = this.calculateAf(v1, v3);
				double u = new Random().nextDouble();
				if (u<Af) {
					this.g.getGraph().addEdge(v1,  v3);
					this.btwEdgeList.add(new edge(v1,v3));
					num_btw_triangles += (num_tri_after - num_tri_prev);
					System.out.printf("(%s, %s) n_prev =%d, n_after=%d, inter_triangle=%d of %d \n", 
						v1, v3, num_tri_prev, num_tri_after, num_btw_triangles, this.btw_triangle);
				}else {
					this.g.getGraph().addEdge(vs, vt);
					this.btwEdgeList.add(e);
				}
			}else {
				this.g.getGraph().addEdge(vs, vt);
				this.btwEdgeList.add(e);
			}
			
		}
		if (file_w.isEmpty())
			this.g.printGraph(this.dataname + "_graph_with_triangles.txt");
		else {
			int beginIndex = file_w.lastIndexOf("_final");
			this.g.printGraph(file_w.substring(0, beginIndex)+".txt");
		}

		this.postProcessing();
	}
	

	public int num_common_neighbours(String vi, String vj) {
		Set<DefaultEdge> n_i = g.getGraph().edgesOf(vi);
		int a = 0;
		for(DefaultEdge e : n_i) {
			String neighbour = g.getGraph().getEdgeTarget(e);
			if (vi.contentEquals(neighbour))
				neighbour = g.getGraph().getEdgeSource(e);
			
			if(g.getGraph().getEdge(vj, neighbour)!=null)
				a++;
		}
		return a;

	}

	public int num_common_neighbours_in(String vi, String vj) {
		ArrayList<String> com = this.comList.get(this.comMap.get(vi));
		Set<DefaultEdge> n_i = g.getGraph().edgesOf(vi);

		int a = 0;
		for(DefaultEdge e : n_i) {
			String neighbour = g.getGraph().getEdgeTarget(e);
			if (vi.contentEquals(neighbour))
				neighbour = g.getGraph().getEdgeSource(e);
			if (!com.contains(neighbour))
				continue;
			
			if(g.getGraph().getEdge(vj, neighbour)!=null)
				a++;
		}
		return a;
	}
	int calNumIntraTriangles() {
		int num_tria = 0;
		
		for (int com_index : this.comList.keySet()) {
			ArrayList<String> com = this.comList.get(com_index);
			if(com.size()==1 && com.get(0).isEmpty()) {
				continue;
			}
			if(com.size()==0)
				continue;
			for (int i = 0; i<com.size(); i++)
				for (int j = i+1; j<com.size(); j++)
					for (int k = j+1; k<com.size(); k++) {
						String v1 = com.get(i);
						String v2 = com.get(j);
						String v3 = com.get(k);
						
						if(g.getGraph().containsEdge(v1, v2) && g.getGraph().containsEdge(v2,  v3) && 
								g.getGraph().containsEdge(v3, v1)) 
							num_tria ++ ;
					}
		}
		return num_tria;
	}
	
	int calNumIntraTriangles_v2() {
		int num_tria = 0;

		for (int com_index : this.comList.keySet()) {
			ArrayList<String> com = this.comList.get(com_index);
			
			Set<String> set = new HashSet<String>(com);
			AsSubgraph sg = new AsSubgraph(g.getGraph(), set);
			int n_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(sg);
			num_tria += n_tri;
			System.out.printf("com_index=%d, n_tri=%d\n", com_index, n_tri);
		}
		return num_tria;
	}
	
	public void Generate_graph_triangles_iterative(String file_wt, String file_w, double threshold) throws IOException {
		this.generate_graph_no_triangles(file_wt);
	
		int total_num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		if (this.btw_triangle ==0)
			return;
		
		int num_tries = 0;
		while (total_num_tri < (this.btw_triangle + this.in_triangle) * threshold && num_tries <5) {
			num_tries ++;
			generate_graph_triangle(file_w);
			total_num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		}

		this.g.printGraph(file_w);
	
	}
	void postProcessing() {
		ConnectivityInspector<String, DefaultEdge> cinspect = new ConnectivityInspector<String, DefaultEdge>(g.getGraph());
		System.out.println("Start post-processing ...");

		int num_tries = 0; 
		while(!cinspect.isConnected() && num_tries <=10) {
			num_tries ++ ;
			List<Set<String>> S =  cinspect.connectedSets();
			
			int maxSize = 0;
			int main_sid = -1;
			for (int sid = 0; sid <S.size(); sid++) {
				if (maxSize < S.get(sid).size()) {
					maxSize = S.get(sid).size();
					main_sid = sid;
				}
			}
			System.out.println("Main componet calculated.");

			HashMap<String, Double> pi_DegMainS= new HashMap<String, Double>();
			
			int sum_degree = 0;
			for (String v : S.get(main_sid)) {
				int degree = this.btwdegreeMap.get(v) + this.indegreeMap.get(v);
				pi_DegMainS.put(v, degree*1.0);
				sum_degree += degree;
			}
			for (String v : S.get(main_sid)) {
				double prob = pi_DegMainS.get(v)/sum_degree; 
				pi_DegMainS.put(v, prob);
			}
			
			for (int sid =0; sid <S.size(); sid ++) {
				System.out.println("sid = " + sid);

				if (sid == main_sid) continue;
				
				Set<String> s = S.get(sid);
				for (String v : s) {
					if (v.isEmpty())
						continue;
					String vt = this.sampleVertex(pi_DegMainS);
					if (!this.g.getGraph().containsEdge(v, vt)) {
						int randomEdge = new Random().nextInt(this.intraEdgeList.size() + this.btwEdgeList.size());
						if(randomEdge < this.intraEdgeList.size()) {
							edge e = this.intraEdgeList.remove(randomEdge);
							this.g.getGraph().removeEdge(e.v1, e.v2);
						}else {
							edge e = this.btwEdgeList.remove(randomEdge - this.intraEdgeList.size());
							this.g.getGraph().removeEdge(e.v1, e.v2);
						}
						this.g.getGraph().addEdge(v, vt);
						if(this.comMap.get(v) == this.comMap.get(vt))
							this.intraEdgeList.add(new edge(v, vt));
						else
							this.btwEdgeList.add(new edge(v,vt));
					}
				}
			}
			cinspect = new ConnectivityInspector<String, DefaultEdge>(g.getGraph());

		}
		
	}
	
}

