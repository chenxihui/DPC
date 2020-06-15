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

public class TriCycle {
	HashMap<String, Integer> degreeMap;
	LinkedList<edge> edgeList;
	graph g;
	String dataname;
	int num_triangle;
	HashMap<String, Double> pi;
	
	TriCycle(String dataname, String inputFilename) throws IOException{
		this.g = new graph();
		degreeMap = new HashMap<String, Integer>();
		edgeList = new LinkedList<edge>();
		this.dataname = dataname;
		readInputFile(inputFilename);

		pi = new HashMap<String, Double>();

		int total_degree = 2*calNumEdge();

		for(String v : this.degreeMap.keySet())
			pi.put(v, this.degreeMap.get(v)*1.0/total_degree);

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
				this.g.getGraph().addEdge(v1, v2);
				this.edgeList.add(new edge(v1, v2));
				m_gen ++;
				System.out.printf("edge %d/%d: (%s, %s) is added.\n", m_gen, m, v1, v2);
				
			}
		}
		
		if (file_wt.isEmpty())
			g.printGraph(this.dataname + "_TriCycle_graph_without_triangles.txt");
		else
			g.printGraph(file_wt);
		
	}
	
	void generate_graph_triangles(String file_wt, String file_w) throws IOException {
		if (!file_wt.isEmpty())
			this.generate_graph_no_triangles(file_wt);
		
		int n_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		
		System.out.println("Start increasing triangles ...");
		while (n_tri < this.num_triangle) {
			String v1 = sampleVertex(pi);
			if (!g.getGraph().containsVertex(v1))
				continue;
			String v2 = this.sampleNeighbour(v1);
			
			if(!g.getGraph().containsVertex(v2))
				continue;
			
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
					this.g.getGraph().addEdge(v1,  v3);
					this.edgeList.add(new edge(v1, v3));
					
					n_tri += (n_after - n_prev);
					
					System.out.printf("(%s, %s) is added. n_prev=%d, n_after=%d, n_triangles=%d/%d\n", v1, v3, n_prev,
							n_after, n_tri, this.num_triangle);
					
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
	public void Generate_graph_triangles_iterative(String file_wt, String file_w, double threshold) throws IOException {
		this.generate_graph_no_triangles(file_wt);
	
		int total_num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());

		while (total_num_tri < this.num_triangle * threshold) {
			this.generate_graph_triangles("", file_w);
			total_num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		}

		this.g.printGraph(file_w);
	
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

		int ntimes = 0;
		while(!cinspect.isConnected() && ntimes < 16) {
			ntimes ++ ;
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
							for (edge e:this.edgeList) {
								String v1 = e.v1;
								String v2 = e.v2;
								if ((v1.contentEquals(v) && v2.contentEquals(nv)) || (v1.contentEquals(nv) && v2.contentEquals(v))) {
									this.edgeList.remove(e);
									break;
								}
							}
						}
					}
					//m -= neighbourList.size();
					m = this.g.getGraph().edgeSet().size();
					System.out.println(m);
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
							this.edgeList.add(new edge(v, target));
							m ++ ;
						}
						while(m >= m_expected) {
							Set<DefaultEdge> eSet = g.getGraph().edgeSet();
							ArrayList<DefaultEdge> eList = new ArrayList<DefaultEdge>();
							eList.addAll(eSet);
							int pos = new Random().nextInt(eSet.size());
							String mv = this.g.getGraph().getEdgeSource(eList.get(pos));
							String nv = this.g.getGraph().getEdgeTarget(eList.get(pos));
							g.getGraph().removeEdge(eList.get(pos));
							m --;
							for (edge e:this.edgeList) {
								String v1 = e.v1;
								String v2 = e.v2;
								if ((v1.contentEquals(mv) && v2.contentEquals(nv)) || (v1.contentEquals(nv) && v2.contentEquals(mv))) {
									this.edgeList.remove(e);
									break;
								}
							}

						}
					}
				}
			}	
			cinspect = new ConnectivityInspector<String, DefaultEdge>(g.getGraph());
							
		}
	}
		
}
