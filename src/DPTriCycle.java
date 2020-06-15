import java.util.List;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;
import java.util.Set;

import org.jgrapht.GraphMetrics;
import org.jgrapht.alg.connectivity.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;

public class DPTriCycle {
	HashMap<String, Integer> degreeMap;
	LinkedList<edge> edgeList;
	graph g;
	String dataname;
	int num_triangle;
	HashMap<String, Double> pi;
	
	HashMap<String, int[]> nodeAbtMap;
	HashMap<Integer, double[]> theta_F;
	HashMap<Integer, double[]> theta_F_g1;
	HashMap<Integer, Double> sup;
	
	DPTriCycle(String dataname, String inputFilename,  String abtFile, String Ffile, 
			String Fg1File, String supFile) throws IOException{
		this.g = new graph();
		degreeMap = new HashMap<String, Integer>();
		edgeList = new LinkedList<edge>();
		this.dataname = dataname;
		readInputFile(inputFilename);

		pi = new HashMap<String, Double>();

		int total_degree = 2*calNumEdge();

		for(String v : this.degreeMap.keySet())
			pi.put(v, this.degreeMap.get(v)*1.0/total_degree);
		
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
	
	void readInputFile(String inputFilename) throws IOException {
		
		BufferedReader fin = new BufferedReader(new FileReader(inputFilename));
		String line = fin.readLine();
		String [] fds = line.split("\t");
		
		//this.num_triangle = Integer.parseInt(fds[0]) + Integer.parseInt(fds[1]);
		this.num_triangle = (int)(Double.parseDouble(fds[0])+ Double.parseDouble(fds[1]));
		
		while ((line = fin.readLine())!= null) {
			fds = line.split("\t");
			String node = fds[0];
			degreeMap.put(node, Integer.parseInt(fds[2])+Integer.parseInt(fds[3]));
			if (!g.getGraph().containsVertex(node))
				g.getGraph().addVertex(node);
		}
	}
	
	int calNumEdge() {
		int num_edge = 0;
		for (String v : this.degreeMap.keySet())
			num_edge += this.degreeMap.get(v);
		return num_edge/2;
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
	void generate_graph_no_triangles(String file_wt) throws IOException{
		int m = calNumEdge();
		int m_gen = 0 ;
				
		while (m_gen < m) {
			String v1 = this.sampleVertex(pi);
			String v2 = this.sampleVertex(pi);
			
			if (!v1.contentEquals(v2) && !this.g.getGraph().containsEdge(v1, v2)) {
				double Af = this.calculateAf(v1, v2);
				double u = new Random().nextDouble();
				if (u < Af) {
					this.g.getGraph().addEdge(v1, v2);
					this.edgeList.add(new edge(v1, v2));
					m_gen ++;
					System.out.printf("edge %d/%d: (%s, %s) is added.\n", m_gen, m, v1, v2);
				}
			}
		}
		
		if (file_wt.isEmpty())
			g.printGraph(this.dataname + "_TriCycle_graph_without_triangles.txt");
		else
			g.printGraph(file_wt);
		
	}
	
	void generate_graph_triangles(String file_wt, String file_w) throws IOException {
		this.generate_graph_no_triangles(file_wt);
		
		int n_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		
		System.out.println("Start increasing triangles ...");
		while (n_tri < this.num_triangle) {
			String v1 = sampleVertex(pi);
			String v2 = this.sampleNeighbour(v1);
			
			if (v2.isEmpty())
				continue;
			String v3 = this.sampleNeighbour(v2);
			if(v3.isEmpty())
				continue;
			
			if (!v1.contentEquals(v3) && !this.g.getGraph().containsEdge(v1, v3)) {
				edge e = this.edgeList.pop();
				String vs = e.v1, vt = e.v2;
				
				int n_prev = this.num_common_neighbours(vs, vt);
				this.g.getGraph().removeEdge(vs,  vt);
				int n_after = this.num_common_neighbours(v1,  v3);
				if (n_prev < n_after) {
					double Af = this.calculateAf(v1, v3);
					double u = new Random().nextDouble();
					if (u < Af) {
						this.g.getGraph().addEdge(v1,  v3);
						this.edgeList.add(new edge(v1, v3));
					
						n_tri += (n_after - n_prev);
					
						System.out.printf("(%s, %s) is added. n_prev=%d, n_after=%d, n_triangles=%d/%d\n", v1, v3, n_prev,
							n_after, n_tri, this.num_triangle);
					}else {
						this.g.getGraph().addEdge(vs,  vt);
						this.edgeList.add(new edge(vs, vt));
					}
				}else {
					this.g.getGraph().addEdge(vs,  vt);
					this.edgeList.add(new edge(vs, vt));
				}
			}
		}
		this.postProcessing();
		if (file_w.isEmpty())
			g.printGraph(this.dataname+"_TriCycle_graph_with_griangles_final.txt");
		else
			g.printGraph(file_w);
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
	String sampleNeighbour(String v) {
		ArrayList<String> neiList = new ArrayList<String>();
		for(DefaultEdge e : this.g.getGraph().edgesOf(v)) {
			String nei = this.g.getGraph().getEdgeSource(e);
			if (nei.contentEquals(v))
				nei = this.g.getGraph().getEdgeTarget(e);
				neiList.add(nei);
		}
		if (neiList.isEmpty())
			return "";
		else {
			int i = new Random().nextInt(neiList.size());
			return neiList.get(i);
		}
	}
	ArrayList<String> getNeighbourList(String v){
		ArrayList<String> neiList = new ArrayList<String>();
		for(DefaultEdge e : this.g.getGraph().edgesOf(v)) {
			String nei = this.g.getGraph().getEdgeSource(e);
			if (nei.contentEquals(v))
				nei = this.g.getGraph().getEdgeTarget(e);
				neiList.add(nei);
		}
		return neiList;
	}
	void postProcessing() {
		int m_expected = this.calNumEdge();
		
		ConnectivityInspector<String, DefaultEdge> cinspect = new ConnectivityInspector<String, DefaultEdge>(g.getGraph());
		System.out.println("Start post-processing ...");

		int n_times = 0;
		while(!cinspect.isConnected() && n_times <=50) {
			n_times ++;
			int m = g.getGraph().edgeSet().size();
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
			
			for (int sid =0; sid <S.size(); sid ++) {
				System.out.println("sid = " + sid);

				if (sid == main_sid) continue;
				
				Set<String> s = S.get(sid);
				
				for (String v : s) {
					ArrayList<String> neighbourList = getNeighbourList(v);
					if (!neighbourList.isEmpty()) {
						for(String nv : neighbourList) {
							g.getGraph().removeEdge(v,nv);
						}
					}
					m -= neighbourList.size();
					
					for (int i= 0; i< this.degreeMap.get(v); i++) {
						int n_tries = 0;
						String target = "";
						
						while (n_tries < S.get(main_sid).size()) {
							target = sampleVertex(pi);
							if (S.get(main_sid).contains(target))
								break;
							else
								n_tries ++;
						}
						if (n_tries >= S.get(main_sid).size()) {
							ArrayList<String> list = new ArrayList<String>();
							list.addAll(S.get(main_sid));
							target =list.get(new Random().nextInt(S.get(main_sid).size()));
						}
						if(!target.isEmpty()) {
							g.getGraph().addEdge(v, target);
							m ++ ;
						}
						if(m >= m_expected) {
							Set<DefaultEdge> eSet = g.getGraph().edgeSet();
							ArrayList<DefaultEdge> eList = new ArrayList<DefaultEdge>();
							eList.addAll(eSet);
							int pos = new Random().nextInt(eSet.size());
							g.getGraph().removeEdge(eList.get(pos));
							m --;
						}
					}
				}
			}	
			cinspect = new ConnectivityInspector<String, DefaultEdge>(g.getGraph());
							
		}
	}
	public double calculateAf(String v1, String v2) {
		int com_index = 0;
		int num_abt = 50;
		if (this.dataname.contentEquals("epinions"))
			num_abt = 50;
		else if (this.dataname.contentEquals("facebook"))
			num_abt = 50;
		else if (this.dataname.contentEquals("petster"))
			num_abt = 13;
		
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

		
}
