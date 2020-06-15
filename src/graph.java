import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.function.Supplier;

import org.jgrapht.Graph;
import org.jgrapht.generate.GnpRandomGraphGenerator;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.util.SupplierUtil;

public class graph {
	private Graph<String,DefaultEdge> g;
	private Graph<String,DefaultWeightedEdge> g1;
	//global attribute entropy
	
	public graph() {
		
		this.g = new SimpleGraph<String, DefaultEdge> (DefaultEdge.class);

	}
	public graph(int i) {
		
		this.g1 = new SimpleGraph<String, DefaultWeightedEdge> (DefaultWeightedEdge.class);
	}

	public Graph<String, DefaultWeightedEdge> getWeightedGraph(){
		return this.g1;
	}
	public void readWeightedGraphOnly(String file) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line = "";
		while((line = in.readLine()) != null) {
			String[] edge_vertice = line.trim().split(" ");
			String source = edge_vertice[0];
			String target = edge_vertice[1];
			

			if(! this.g1.containsVertex(source))
				this.g1.addVertex(source);
			if(! this.g1.containsVertex(target))
				this.g1.addVertex(target);

			if(!source.equals(target)) {
				this.g1.addEdge(source,  target);
			}
		}
		in.close();
	}
	public void readGraphOnly(String file) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line = "";
		while((line = in.readLine()) != null) {
			String[] edge_vertice = line.trim().split(" ");
			String source = edge_vertice[0];
			String target = edge_vertice[1];
			

			if(! this.g.containsVertex(source))
				this.g.addVertex(source);
			if(! this.g.containsVertex(target))
				this.g.addVertex(target);

			if(!source.equals(target)) {
				this.g.addEdge(source,  target);
			}
		}
		in.close();
	}
	//create a graph with g
	public graph(Graph<String,DefaultEdge> g) {
		this.g=g;
	}
	
	// Read an existing graph from the file located at graphFile
	public graph(String graphFile) {
		g = new SimpleGraph<String, DefaultEdge>(DefaultEdge.class);
		
		try {
			String attribute_file = graphFile.substring(0, graphFile.indexOf("."));
			
			BufferedReader file = new BufferedReader(new FileReader(attribute_file));
			String line = file.readLine();
			
			//read the number of nodes, edges and attributes
			
			
			file.close();
			
			
			//read edges
			BufferedReader gFile = new BufferedReader(new FileReader(graphFile));
			while((line=gFile.readLine())!=null){
				String [] edge = line.split(" ");
				g.addEdge(edge[0], edge[1]);
			}
			gFile.close();
		}catch(IOException e){}
		
	}
	
	//Generate a random graph with node as AVertex and assign the generated graph to g
	public Graph<AVertex, DefaultEdge> generateRandomGraph(int numberVertex, double probability, String filepath)
			throws IOException {

		Supplier<AVertex> vSupplier = new Supplier<AVertex>() {
			private int id = 0;

			@Override
			public AVertex get() {

				String vid = "v" + id++;

				AVertex v = new AVertex(vid);
				v.generate_random_attribute();
				return v;
			}
		};

		Graph<AVertex, DefaultEdge> gnpGraph = new SimpleGraph<AVertex, DefaultEdge>(vSupplier,
				SupplierUtil.createDefaultEdgeSupplier(), false);

		// Create the CompleteGraphGenerator object
		GnpRandomGraphGenerator<AVertex, DefaultEdge> gnpGenerator = new GnpRandomGraphGenerator<AVertex, DefaultEdge>(
				numberVertex, probability);

		gnpGenerator.generateGraph(gnpGraph);

		System.out.println(gnpGraph.vertexSet().size());
		File file = new File(filepath);
		if (!file.exists())
			file.createNewFile();

		FileWriter writer = new FileWriter(file);

		// Print out the graph to be sure it's really complete

		Iterator<AVertex> iter = gnpGraph.vertexSet().iterator();
		while (iter.hasNext()) {
			AVertex vertex = iter.next();
			String str = vertex.getID() + ": " + vertex.getAttributeInString();
			System.out.println(str);
			writer.write(str + "\n");
		}

		Iterator<DefaultEdge> edgeiter = gnpGraph.edgeSet().iterator();
		while (edgeiter.hasNext()) {
			DefaultEdge edge = edgeiter.next();
			System.out.println(edge.toString());
			writer.write(edge.toString() + "\n");

		}
		writer.close();
		return gnpGraph;
	}

	//Generate a random graph with string vertices and assign the result to g. 
	//The generated graph will stored as a dat file at filepath.
	public Graph<String, DefaultEdge> generateRandomAttributeGraph(int numberVertex, double probability, int attribute_number,
			String filepath) throws IOException {

		Supplier<String> vSupplier = new Supplier<String>() {
			private int id = 0;

			@Override
			public String get() {

				String vid = Integer.toString(id);
				id++;

				return vid;
			}
		};

		Graph<String, DefaultEdge> gnpGraph = new SimpleGraph<String, DefaultEdge>(vSupplier,
				SupplierUtil.createDefaultEdgeSupplier(), false);

		// Create the CompleteGraphGenerator object
		GnpRandomGraphGenerator<String, DefaultEdge> gnpGenerator = new GnpRandomGraphGenerator<String, DefaultEdge>(
				numberVertex, probability);

		gnpGenerator.generateGraph(gnpGraph);

		String attribute_file = filepath.substring(0, filepath.indexOf(".")) +"_attribute.txt";
	    File fwrite = new File(attribute_file);
	    if(fwrite.exists())
	    	fwrite.delete();
		
		FileWriter writer = new FileWriter(attribute_file, true);

		// Print out the graph to be sure it's really complete
		writer.write(gnpGraph.vertexSet().size() + " " + gnpGraph.edgeSet().size() + " "+attribute_number+ "\n");

		Iterator<String> iter = gnpGraph.vertexSet().iterator();
		while (iter.hasNext()) {
			String vertex = iter.next();
			String str = vertex + ": [";

			ArrayList<Integer> attribute = new ArrayList<Integer>(attribute_number);

			for (int i = 0; i < attribute_number; i++) {
				Random rng = new Random();
				int value_to_add = 0;
				if (rng.nextDouble() <= 0.5)
					value_to_add = 1;

				attribute.add(value_to_add);
				if (i != attribute_number - 1)
					str = str + value_to_add + ", ";
				else
					str = str + value_to_add + "]";
			}

		//	System.out.println(str);
			writer.write(str + "\n");
		}
		writer.close();
		fwrite = new File(filepath);
	    if(fwrite.exists())
	    	fwrite.delete();

		writer = new FileWriter(filepath, true);
		
		
		Iterator<DefaultEdge> edgeiter = gnpGraph.edgeSet().iterator();
		while (edgeiter.hasNext()) {
			DefaultEdge edge = edgeiter.next();
			
			String source = gnpGraph.getEdgeSource(edge);
			String target = gnpGraph.getEdgeTarget(edge);

			// System.out.println(edge.toString());
			writer.write(source + " " + target + "\n");
		}
		writer.close();
		
		this.g = gnpGraph;
		System.out.println("Graph is generated.");
		return gnpGraph;
	}

	public int[] getNumberDegreeInComs(String v, Partition p) {
		int [] degInCom = new int[p.partition.size()];
		for (int i = 0 ; i<p.partition.size(); i++)
			degInCom[i] = 0;
		for(DefaultEdge e : g.edgesOf(v)) {
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
		
		Graph<String, DefaultEdge> graph = this.g;
		int nEdges = 0;
		for (String v : vlist) {
			for(DefaultEdge e :graph.edgesOf(v)) {
				String vs = graph.getEdgeSource(e);
				if(vs.equalsIgnoreCase(v))
					vs = graph.getEdgeTarget(e);
				
				if(vlist.contains(vs))
					nEdges ++;
			}
		}
		/*
		for (DefaultEdge e : graph.edgeSet()) {
			String source = graph.getEdgeSource(e);
			String target = graph.getEdgeTarget(e);

			if (vlist.contains(source) && vlist.contains(target))
				nEdges++;
		}

		return nEdges;
		*/
		return nEdges/2;

	}

	//return the total number of edges of a given list of nodes
	
	public int getSumOfDegrees(ArrayList<String> vlist) {
		
		if(vlist.size() ==0)
			return 0;
		
		Graph<String, DefaultEdge>  graph =this.g;
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
	
	//Calculate the modularity of partition in g
	public double getModularityAttribute(attributeGraph ag, Partition p, double w_s, double w_a) {
		//Graph<String, DefaultEdge> graph = this.g;
		
		int nEdges = g.edgeSet().size();

		double Q_s = 0;
		double Q_a = 0;

		Iterator<ArrayList<String>> iter = p.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> c = iter.next();
			
			if(c.isEmpty())
				continue;
			
			int cEdges = getNumberOfEdges(c);
			int cDegrees = getSumOfDegrees(c);
			Q_s += ((double) cEdges / nEdges - Math.pow((double) cDegrees / (2 * nEdges), 2));
			// System.out.println(cEdges + " " +nEdges+ " "+cDegrees+ " "+Q);
			
		}
		Q_a = ag.getAttributeModularity(p);
	
		return w_s*Q_s + w_a*Q_a;
	}
	
   	public double getPartialModularity(Partition p) {
		//Graph<String, DefaultEdge> graph = this.g;
		
		//int nEdges = g.edgeSet().size();
   		int nEdges = p.nEdges;
		double Q = 0;

		Iterator<ArrayList<String>> iter = p.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> c = iter.next();
			
			if(c.isEmpty())
				continue;
			
			int cEdges = getNumberOfEdges(c);
			//int cDegrees = getSumOfDegrees(c);
			int cDegrees = 0;
			for (String v :c) {
				cDegrees += p.degreeP.get(v);
			}
			Q += ((double) cEdges / nEdges - Math.pow((double) cDegrees*1.0 / (2 * nEdges), 2));
			// System.out.println(cEdges + " " +nEdges+ " "+cDegrees+ " "+Q);
		}
		return Q;
	}
	public double getPartialModularity_abt(attributeGraph ag, Partition p, double w_s, double w_a) {
		//Graph<String, DefaultEdge> graph = this.g;
		
		//int nEdges = g.edgeSet().size();
   		int nEdges = p.nEdges;
		double Q_s = 0;

		Iterator<ArrayList<String>> iter = p.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> c = iter.next();
			
			if(c.isEmpty())
				continue;
			
			int cEdges = getNumberOfEdges(c);
			//int cDegrees = getSumOfDegrees(c);
			int cDegrees = 0;
			for (String v :c) {
				cDegrees += p.degreeP.get(v);
			}
			Q_s += ((double) cEdges / nEdges - Math.pow((double) cDegrees*1.0 / (2 * nEdges), 2));
			// System.out.println(cEdges + " " +nEdges+ " "+cDegrees+ " "+Q);
		}
		double Q_a = ag.getAttributePartialModularity(p);
		return w_s*Q_s + w_a * Q_a;
	}
	public double getModularity(Partition p) {
		//Graph<String, DefaultEdge> graph = this.g;
		
		int nEdges = g.edgeSet().size();
//		int nEdges = this.getNumberOfEdges(p.getVertexList());

		double Q = 0;

		Iterator<ArrayList<String>> iter = p.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> c = iter.next();
			
			if(c.isEmpty())
				continue;
			
			int cEdges = getNumberOfEdges(c);
			int cDegrees = getSumOfDegrees(c);
			Q += ((double) cEdges / nEdges - Math.pow((double) cDegrees / (2 * nEdges), 2));
		}
		return Q;
	}
	
	
	public Graph<String, DefaultEdge> getGraph() {
		return this.g;
	}
	
	public void printGraph(String outputFile) throws IOException {
		FileWriter fout = new FileWriter(outputFile);
		for (DefaultEdge e : this.g.edgeSet()) {
			String vs = this.g.getEdgeSource(e);
			String vt = this.g.getEdgeTarget(e);
			fout.write(vs + " " + vt + "\n");
		}
		fout.close();
	}
}
