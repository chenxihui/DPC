import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.jgrapht.GraphMetrics;
import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultEdge;

public class LadderFunction {
	
	public static int num_common_neighbours(graph g, String vi, String vj) {
		Set<DefaultEdge> n_i = g.getGraph().edgesOf(vi);
		//Set<DefaultEdge> n_j = g.getGraph().edgesOf(vj);
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

	public static int num_common_neighbours_in(graph g, String vi, String vj, ArrayList<String> com) {
		Set<DefaultEdge> n_i = g.getGraph().edgesOf(vi);

		//Set<DefaultEdge> n_j = g.getGraph().edgesOf(vj);
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
	/**
	 * Compute the number of nodes that connect to either node vi or nj but not both
	 * @param g
	 * @param vi
	 * @param vj
	 * @return
	 */
	public static int get_degree_com(graph g, String v, ArrayList<String> com) {

		Set<DefaultEdge> e_v = g.getGraph().edgesOf(v);
		
		int a = 0 ;
		for (DefaultEdge e : e_v) {
			String neighbour = g.getGraph().getEdgeTarget(e);
			if (neighbour.contentEquals(v))
				neighbour = g.getGraph().getEdgeSource(e);
			if (com.contains(neighbour))
				a ++;
		}
		
		return a;
	}

	public static int num_node_connect_to_one_node_in(graph g, String vi, String vj, int a, ArrayList<String> com){
		int b=0;
		/**
		 * degree of vi and vj
		 */
		int di = get_degree_com(g, vi, com);
		
		int dj = get_degree_com(g, vj, com); 
		
		int x= g.getGraph().containsEdge(vi, vj)?1:0;
		
		b = di + dj - 2*a - 2*x;
		
		return b;
	}

	public static int num_node_connect_to_one_node(graph g, String vi, String vj, int a) {
		int b=0;
		/**
		 * degree of vi and vj
		 */
		int di = g.getGraph().edgesOf(vi).size();
		int dj = g.getGraph().edgesOf(vj).size();
		
		int x= g.getGraph().containsEdge(vi, vj)?1:0;
		
		b = di + dj - 2*a - 2*x;
		
		return b;
	}
	public static int LS_g_t(graph g, long t, int[][] a, int[][] b, String[] vlist) {
		
		String [] V = vlist;
		
	/**
	 * n is the number of nodes in graph g
	 */
	    int n = V.length;
	    /**
	     * Compute LS(g,t)
	     */
	    int LS = 0; 
		for (int i=0; i<V.length; i++) {
			for(int j=i +1 ; j<V.length; j++) {
				
				int LSij = (int) Math.min(n-2, a[i][j]+Math.floor((t+Math.min(t,b[i][j]))/2));
				
				if(LS<LSij) LS = LSij;

				if (LS == n-2) return LS;
			}
		}
		return LS;
	}
	
	public static int[] LS_g_t_in(graph g, Partition p, int M) {
		
		int n = g.getGraph().vertexSet().size();

				
		Iterator<ArrayList<String>> iter = p.partition.iterator();
		
		int[] res = new int[M+1];
	
		Arrays.fill(res, 0);

		
		int com_index = 0; 
		while(iter.hasNext()) {
			System.out.println("current com=" + com_index);
			
			ArrayList<String> com = iter.next();
			for (int i = 0; i < com.size(); i ++) {
				

				if(!g.getGraph().containsVertex(com.get(i)))
					continue;
				for(int j = i+1; j<com.size(); j++) {
					if(!g.getGraph().containsVertex(com.get(j)))
						continue;
					System.out.println("current com=" + com_index + " i=" + i +"/" + com.size() + " j=" + j);
					int a = num_common_neighbours_in(g, com.get(i), com.get(j), com);
					int b = num_node_connect_to_one_node_in(g, com.get(i), com.get(j), a, com);
					
					for(int k = 0; k<=M; k++) {
						if (res[k] == n-2 )
							continue;
						int LSij_k = (int) Math.min(n-2, a + Math.floor((k+Math.min(k,b))/2));
						if(res[k] < LSij_k)
							res[k] = LSij_k;
					}
				}
			}
			com_index ++ ;
		}
	    
		return res;
	}
	
	public static int countTriangle_normalDP(graph g, double epsilon) {
		int num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		
		System.out.println("the correct number of triangles: "+num_tri);
		
	    int n = g.getGraph().vertexSet().size();

		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, (n - 2.0)/ epsilon);
	
		return (int) (num_tri + Math.floor(lng.nextLaplacian()));
		
	}
	public static int countTriangleDP(graph g, double epsilon) {
		/**
		 * dp_num_tri is the number of triangles under DP
		 * num_tri is the true number of triangles in g
		 * n is the number of nodes in g
		 */
		int dp_num_tri = 0;
		int num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(g.getGraph());
		
		System.out.println("the correct number of triangles: "+num_tri);
		
	    int n = g.getGraph().vertexSet().size();

	    int M = 2 * n; 
	    int[] leftB = new int[M+1];
	    int[] leftE = new int[M+1];
	    int[] rightB = new int[M+1];
	    int[] rightE = new int[M+1];
	    double [] weight = new double[M+2];
	    
	    leftB[0]= num_tri;
	    leftE[0]= num_tri;
	    rightB[0]= num_tri;
	    rightB[0]= num_tri;
	    
	    //weight[0] = 1 * exp(e/2*deltaq *0)
	    weight[0] = 1;

	    /**
	     * dst is the distance of I_{t-1} from num_tri
	     */

	    int dst = 0;
	    
	    int [][] a = new int[n][n];
	    int [][] b = new int[n][n];
	    
	    String[] vlist =  g.getGraph().vertexSet().toArray(new String[0]);

	    for (int i = 0; i < vlist.length; i ++)
				for(int j = i+1; j<vlist.length; j++) {
					a[i][j] = num_common_neighbours(g, vlist[i], vlist[j]);
					b[i][j] = num_node_connect_to_one_node(g, vlist[i], vlist[j], a[i][j]);
				}
	    
	    for (int t=1; t<=M; t++) {
	    	System.out.println("\t LS(g, t) with t=" + t);
	    	/**
	    	 * I_t_mois_1 is the I_{t-1}
	    	 */
	    	int I_t_mois_1 = (int) LS_g_t(g,t-1,a,b, vlist);
	    	/**
	    	 * Calculate the bucket of level t
	    	 * (LeftB[t], leftE[t]) (rightB[t], rightE[t])
	    	 */
	    	leftB[t] = num_tri - dst - I_t_mois_1;
	    	leftE[t] = num_tri - dst;
	    	rightB[t] = num_tri + dst;
	    	rightE[t] = num_tri + dst + I_t_mois_1;
	    	
	    	weight[t] = 2*I_t_mois_1 * Math.exp(-1*(epsilon/2) * t);
	    	dst += I_t_mois_1;
	    }
	    weight[M+1] = 2*(n-2) * Math.exp(-1* (epsilon/2)*(M+1))/(1-Math.exp(-1*epsilon/2));
	    /**
	     * calcualte the probability distribution over t by normalising weight
	     */
	    double[] prob = new double[M+2];
	    double totalWeight = 0; 
	    for (int t=0; t<=M+1; t++) 
	    	totalWeight += weight[t];
	    
	    for(int t=0; t<=M+1; t++) {
	    	prob[t] = weight[t]/totalWeight;
	    }
	    
	    /**
	     * sample t based on prob
	     */
	    double rand = Math.random();
	    
	    totalWeight = prob[0] ; int t=0;
	    while (totalWeight<rand && t<=M+1) {
	    	t++;
	    	totalWeight+=prob[t];
	    }
	    //Sample the count of triangles from range[t]
	    if(t<=M) {
	    	int range = leftE[t] - leftB[t]+1;
	    	dp_num_tri=-1;
			while (dp_num_tri < 0) {
				rand = Math.random();
				
				if (rand < 0.5)
					dp_num_tri = leftB[t] + (int) Math.round(rand / (0.5 / range));
				else
					dp_num_tri = rightB[t] + (int) Math.round((rand-0.5) / (0.5 / range));
			}
	    }else {
	    	rand = Math.random();
	    	int i = (int) Math.round(Math.log(1-rand)/Math.log(Math.exp(-1*epsilon/2)));
		    	
	    	int range = n-2;
	    	dp_num_tri=-1;
			while (dp_num_tri < 0) {
				rand = Math.random();
				if (rand < 0.5)
					dp_num_tri = num_tri- (dst + i *(n-2)+(int) Math.round(rand / (0.5 / range)));
				else
					dp_num_tri = num_tri + dst + i*(n-2)+(int) Math.round((rand-0.5) / (0.5 / range));
			}
	    }
		return dp_num_tri;
	}

	public static int countTriangle_normalDP(graph g, double epsilon, Partition p) {
		int num_tri_in = calNumIntraTriangles_v2(g, p);
		System.out.println("the correct number of intra-community triangles: "+ num_tri_in);
		
	    int n = g.getGraph().vertexSet().size();

		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, (n - 2.0)/ epsilon);
	
		return (int) (num_tri_in + Math.floor(lng.nextLaplacian()));
	}
	
	static int calNumIntraTriangles_v2(graph g, Partition p) {
		int num_tria = 0;
		int com_index = 0;
		Iterator<ArrayList<String>> iter = p.partition.iterator();
		while(iter.hasNext()) {
			ArrayList<String> com = iter.next();
			
			Set<String> set = new HashSet<String>(com);
			AsSubgraph sg = new AsSubgraph(g.getGraph(), set);
			int n_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles(sg);
			num_tria += n_tri;
			System.out.printf("com_index=%d, n_tri=%d\n", com_index, n_tri);
			com_index ++;
		}
		return num_tria;
	}
	
	public static int countTriangleDP(graph g, double epsilon, Partition p) {

		int num_tri_in = calNumIntraTriangles_v2(g, p);
		/**
		 * dp_num_tri is the number of triangles under DP
		 * num_tri is the true number of triangles in g
		 * n is the number of nodes in g
		 */
		int dp_num_tri = 0;
		//int num_tri = (int) GraphMetrics.<String, DefaultEdge>getNumberOfTriangles_com(g.getGraph());
		int num_tri = num_tri_in;
		System.out.println("\t The correct number of triangles: "+num_tri);
		
	    int n = g.getGraph().vertexSet().size();
	    
	    int M = 2 * n; 
	    int[] leftB = new int[M+1];
	    int[] leftE = new int[M+1];
	    int[] rightB = new int[M+1];
	    int[] rightE = new int[M+1];
	    double [] weight = new double[M+2];
	    
	    leftB[0]= num_tri;
	    leftE[0]= num_tri;
	    rightB[0]= num_tri;
	    rightB[0]= num_tri;
	    
	    //weight[0] = 1 * exp(e/2*deltaq *0)
	    weight[0] = 1;
	    /**
	     * dst is the distance of I_{t-1} from num_tri
	     */
	    int dst = 0;
		
	    int [] LSgt = LS_g_t_in(g, p, M);
	    
	    for (int t=1; t<=M; t++) {
	    	/**
	    	 * I_t_mois_1 is the I_{t-1}
	    	 */
	    	int I_t_mois_1 = (int) LSgt[t-1];
	    	/**
	    	 * Calculate the bucket of level t
	    	 * (LeftB[t], leftE[t]) (rightB[t], rightE[t])
	    	 */
	    	leftB[t] = num_tri - dst - I_t_mois_1;
	    	leftE[t] = num_tri - dst;
	    	rightB[t] = num_tri + dst;
	    	rightE[t] = num_tri + dst + I_t_mois_1;
	    	
	    	weight[t] = 2*I_t_mois_1 * Math.exp(-1*(epsilon/2) * t);
	    	dst += I_t_mois_1;
	    }
	    weight[M+1] = 2*(n-2) * Math.exp(-1* (epsilon/2)*(M+1))/(1-Math.exp(-1*epsilon/2));
	    /**
	     * calcualte the probability distribution over t by normalising weight
	     */
	    double[] prob = new double[M+2];
	    double totalWeight = 0; 
	    for (int t=0; t<=M+1; t++) 
	    	totalWeight += weight[t];
	    
	    for(int t=0; t<=M+1; t++) {
	    	prob[t] = weight[t]/totalWeight;
	    }
	    
	    /**
	     * sample t based on prob
	     */
	    double rand = Math.random();
	    
	    totalWeight = prob[0] ; int t=0;
	    while (totalWeight<rand && t<=M+1) {
	    	t++;
	    	totalWeight+=prob[t];
	    }
	    //Sample the count of triangles from range[t]
	    if(t<=M) {
	    	int range = leftE[t] - leftB[t]+1;
	    	dp_num_tri=-1;
			while (dp_num_tri < 0) {
				rand = Math.random();
				
				if (rand < 0.5)
					dp_num_tri = leftB[t] + (int) Math.round(rand / (0.5 / range));
				else
					dp_num_tri = rightB[t] + (int) Math.round((rand-0.5) / (0.5 / range));
			}
	    }else {
	    	rand = Math.random();
	    	int i = (int) Math.round(Math.log(1-rand)/Math.log(Math.exp(-1*epsilon/2)));
		    	
	    	int range = n-2;
	    	dp_num_tri=-1;
			while (dp_num_tri < 0) {
				rand = Math.random();
				if (rand < 0.5)
					dp_num_tri = num_tri- (dst + i *(n-2)+(int) Math.round(rand / (0.5 / range)));
				else
					dp_num_tri = num_tri + dst + i*(n-2)+(int) Math.round((rand-0.5) / (0.5 / range));
			}
	    }
	    
		return dp_num_tri;
	}
	static int calNumIntraTriangles(graph g, Partition p) {
		int num_tria = 0;
		Iterator<ArrayList<String>> iter = p.partition.iterator();
		
		while(iter.hasNext()) {
			ArrayList<String> com = iter.next();
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
}
