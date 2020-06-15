import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;

public class attributeGraph {
	
	private SimpleWeightedGraph<String, DefaultWeightedEdge> g;
	private HashMap<String, ArrayList<Integer>> attribute;
	private HashMap<String, Double> adegree;
 	HashMap<String, HashMap<String, Double>> weight;
 	int num_abt;
  
 	public HashMap<String, ArrayList<Integer>> getAttributeMatrix() {
 		return this.attribute;
 	}
	public SimpleWeightedGraph<String, DefaultWeightedEdge> getGraph(){
		return this.g;
	}
	
	public attributeGraph() {
		this.g = new SimpleWeightedGraph<String, DefaultWeightedEdge> (DefaultWeightedEdge.class);
	}
	public attributeGraph(graph graph, String attributeFile, double threashold) throws IOException {
		this.g = new SimpleWeightedGraph<String, DefaultWeightedEdge> (DefaultWeightedEdge.class);
		this.attribute= new HashMap<String, ArrayList<Integer>>();
		this.adegree= new HashMap<String, Double>();

		this.readAttributeOnly(attributeFile);

		for (String v : graph.getGraph().vertexSet()) {
			if (!this.g.containsVertex(v))
				this.g.addVertex(v);
		}
		
		ArrayList<String> vList = new ArrayList<String>(graph.getGraph().vertexSet());
        		
        for (int i = 0; i < vList.size(); i++)		
        	for (int j = i+1; i<vList.size(); j++) {
        		String vs = vList.get(i);
        		String vt = vList.get(j);

        		double weight = 1 - getDistance(vs, vt);
        		if (weight > threashold) {
        			DefaultWeightedEdge edge = new DefaultWeightedEdge();
        		
        			this.g.addEdge(vs, vt, edge);
        			this.g.setEdgeWeight(edge, weight);
        			this.adegree.put(vs, this.adegree.get(vs)+weight);
        			this.adegree.put(vt, this.adegree.get(vt)+weight);
        		}
        	}
	}

	public double getDistance(String vs, String vt) {
		ArrayList<Integer> Avs = this.attribute.get(vs);
		ArrayList<Integer> Avt = this.attribute.get(vt);
	    double total = 0; 
		for (int i = 0; i< Avs.size(); i++) {
		   total += Math.abs(Avs.get(i) - Avt.get(i));
		}
		return total*1.0/Avs.size();
	}
	
	public int[] getNumberDegreeInComs(String v, Partition p) {
		int [] degInCom = new int[p.partition.size()];
		for (int i = 0 ; i<p.partition.size(); i++)
			degInCom[i] = 0;
		for(DefaultWeightedEdge e : g.edgesOf(v)) {
			String source = g.getEdgeSource(e);
			if (source.equals(v)) 
				source = g.getEdgeTarget(e);
			if(p.node2com.containsKey(source))
				degInCom[p.node2com.get(source)] ++;
		}
		return degInCom;
	}
	//get the number of edges between a given list of nodes.
	
	public int getNumberOfEdges(ArrayList<String> vlist) {
		
		if(vlist.isEmpty())
			return 0;
		
		int nEdges = 0;
		/*
		for (String v : vlist) {
			for(DefaultEdge e :graph.edgesOf(v)) {
				String source = graph.getEdgeSource(e);
				String target = graph.getEdgeTarget(e);
				
				if(vlist.contains(source) && vlist.contains(target))
					nEdges ++;
			}
		}
		*/
		for (DefaultWeightedEdge e : this.g.edgeSet()) {
			String source = this.g.getEdgeSource(e);
			String target = this.g.getEdgeTarget(e);

			if (vlist.contains(source) && vlist.contains(target))
				nEdges ++;
		}
		return nEdges;

	}

	//return the total number of edges of a given list of nodes
	
	public int getSumOfDegrees(ArrayList<String> vlist) {
		
		if(vlist.size() ==0)
			return 0;
		
		Graph<String, DefaultWeightedEdge>  graph =this.g;
		int nDegrees = 0;

		Iterator<String> iter = vlist.iterator();
		while (iter.hasNext()) {
			String v = iter.next();
			//System.out.println(vlist.toString()+ "the size of vlist is " + vlist.size() + ". node is x"+v + "x");
			if(!v.isEmpty())
				nDegrees += graph.degreeOf(v);
		}
		return nDegrees;
	}

	public double getAttributeModularity(Partition p) {
		int nEdges = this.g.edgeSet().size();

		double Q = 0;

		Iterator<ArrayList<String>> iter = p.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> c = iter.next();
			
			if(c.isEmpty())
				continue;
			
			double cEdges = getNumberOfEdges(c);
			int cDegrees = getSumOfDegrees(c);
			Q += ((double) cEdges / nEdges - Math.pow((double) cDegrees / (2 * nEdges), 2));
			// System.out.println(cEdges + " " +nEdges+ " "+cDegrees+ " "+Q);
		}
		return Q;
	}
	
	public void readAttributeOnly(String file) throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(file));
		
		String line = "";
		while ((line = in.readLine()) != null) {
			String [] fds = line.split(" ");
			ArrayList<Integer> av = new ArrayList<Integer>(fds.length-1);
			
			for(int i =1; i<fds.length; i++) {
				av.add(Integer.parseInt(fds[i]));
			}
			
			this.attribute.put(fds[0], av);
			this.num_abt = av.size();
		}
		in.close();
	}	
 	
	public double getAttributePartialModularity(Partition p) {
		int nEdges = p.nEdges_abt;
		
		double Q = 0;

		Iterator<ArrayList<String>> iter = p.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> c = iter.next();
			
			if(c.isEmpty())
				continue;
			
			double cEdges = getNumberOfEdges(c);
			//int cDegrees = getSumOfDegrees(c);
			double cDegrees = 0;
			for (String v :c) {
				cDegrees += p.degreeP_abt.get(v);
			}
			Q += ((double) cEdges / nEdges - Math.pow((double) cDegrees*1.0 / (2 * nEdges), 2));
		}
		return Q;
	}
	
	public void readGraphOnly(String file, graph g) throws IOException {
		
		for (String node : g.getGraph().vertexSet())
			this.g.addVertex(node);
		
		
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line = "";
		while((line = in.readLine()) != null) {
			String[] edge_vertice = line.trim().split(" ");
			String source = edge_vertice[0];
			String target = edge_vertice[1];

			if(!source.equals(target) && this.g.containsVertex(source) && this.g.containsVertex(target)) {

				this.g.addEdge(source,  target);
			}
		}
		in.close();
	}
	public double getAttributeEntropy(Partition p) {
		double ent = 0 ;
		int N_nodes = this.g.vertexSet().size();
	
		
		Iterator<ArrayList<String>> iter = p.partition.iterator();

	    while (iter.hasNext()) {
	    	ArrayList<String> com = iter.next();
	    	if (com.size() == 0)
	    		continue;
	    	int[] abt_com = new int[this.num_abt];
	    	int n_node_com = com.size();

	    	for (String node:com){
	    		for(int i=0; i<this.attribute.get(node).size(); i++) {
	    			if (this.attribute.get(node).get(i) == 1) {
	    				abt_com[i] ++;
	    			}
	    		}
	    		for(int i=0; i<abt_com.length; i++) {
	    			double ent_com = - abt_com[i]*1.0/N_nodes * (Math.log(abt_com[i]*1.0/N_nodes)/Math.log(2.0))
	    					- (n_node_com - abt_com[i])*1.0/N_nodes * (Math.log((n_node_com - abt_com[i])*1.0/N_nodes)/Math.log(2.0));
	    			ent += ent_com ;
	    		}
	    	}	
	    }
	    ent /= this.num_abt;
	    
	    return ent;

		
	}
	public void readGraphOnly(String file) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line = "";
		while((line = in.readLine()) != null) {
			String[] edge_vertice = line.trim().split(" ");
			String source = edge_vertice[0];
			String target = edge_vertice[1];
			
			this.g.addVertex(source);
			this.g.addVertex(target);

			if(!source.equals(target) && this.g.containsVertex(source) && this.g.containsVertex(target)) {

				this.g.addEdge(source,  target);
			}
		}
		in.close();
	}
}
