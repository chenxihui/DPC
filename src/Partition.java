
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;
import java.util.Iterator;
import java.util.Set;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultWeightedEdge;

public class Partition {
	ArrayList<ArrayList<String>> partition;
	HashMap<String, Integer> node2com;
	int nEdges; 
	int nEdges_abt;
	HashMap<String, Integer> degreeP;
	HashMap<String, Integer> degreeP_abt;
	
	HashMap<Integer, ArrayList<String>> partition_map;

	public Partition() {
		this.partition = new ArrayList<ArrayList<String>>();
	    this.node2com= new HashMap<String, Integer> ();
	}
	
	public void readPartition_Map (String s) {
		if(s.isEmpty()) {
			System.out.println("the sampled partition in McMc is empty!");
			System.exit(0);
		}
			
		this.partition = new ArrayList<ArrayList<String>>();
		String[] blocks = s.substring(1, s.length()-1).split(", ");
		for(String block :blocks) {
			String[] vertices = block.substring(1, block.length()-1).split(",");
			if(vertices.length == 0 )
				continue;
			ArrayList<String> av = new ArrayList<String>(vertices.length);
			for (int x = 0; x<vertices.length; x++) {
				av.add(vertices[x]);
			}
			this.partition.add(av);
		}
	    this.node2com= new HashMap<String, Integer> ();

	}
	// Create a partition from a string like [[v0,v3,v7,v9], [v1,v4,v5,v6], [v2,v8]]
	public Partition(String s) {
		if(s.isEmpty()) {
			System.out.println("the sampled partition in McMc is empty!");
			System.exit(0);
		}
			
		this.partition = new ArrayList<ArrayList<String>>();
		String[] blocks = s.substring(1, s.length()-1).split(", ");
		for(String block :blocks) {
			String[] vertices = block.substring(1, block.length()-1).split(",");
			if(vertices.length == 0 )
				continue;
			ArrayList<String> av = new ArrayList<String>(vertices.length);
			for (int x = 0; x<vertices.length; x++) {
				av.add(vertices[x]);
			}
			this.partition.add(av);
		}
	    this.node2com= new HashMap<String, Integer> ();
	}

	public void calculatePartitionParas_abt(attributeGraph g) {
		ArrayList<String> vlist = this.getVertexList();

		this.degreeP_abt = new HashMap<String, Integer>();

		for (String v : vlist) {
			int degreeV = 0;
			for (DefaultWeightedEdge e: g.getGraph().edgesOf(v)) {
				String vs = g.getGraph().getEdgeTarget(e);
				if (vs.equalsIgnoreCase(v))
					vs = g.getGraph().getEdgeSource(e);
				if (vs.equalsIgnoreCase(v))
					continue;
				if (vlist.contains(vs)) {
					degreeV ++ ;
					this.nEdges_abt ++;
				}	
			}
			this.degreeP_abt.put(v, degreeV);
		}			
		this.nEdges_abt = (int) this.nEdges_abt/2 ;
	}

	public void calculatePartitionParas(graph g) {
		ArrayList<String> vlist = this.getVertexList();
    	this.degreeP = new HashMap<String, Integer>();
    	
		for (String v : vlist) {
			int degreeV = 0;
			for (DefaultEdge e: g.getGraph().edgesOf(v)) {
				String vs = g.getGraph().getEdgeTarget(e);
				if (vs.equalsIgnoreCase(v))
					vs = g.getGraph().getEdgeSource(e);
				if (vs.equalsIgnoreCase(v))
					continue;
				if (vlist.contains(vs)) {
					degreeV ++ ;
					this.nEdges ++;
				}	
			}
			this.degreeP.put(v, degreeV);
				
		}
		this.nEdges = (int) this.nEdges/2;
	
	}

	int getSumOfDegrees(ArrayList<String> com) {
		int total = 0;
		for (String v : com) {
			total += this.degreeP.get(v);
		}
		return total;
	}

	int getSumOfDegrees_abt(ArrayList<String> com) {
		int total = 0;
		for (String v : com) {
			total += this.degreeP_abt.get(v);
		}
		return total;
	}

	Partition(ArrayList<ArrayList<String>> partition) {
		this.partition = partition;
	}

	public Partition clone() {
		Partition p = new Partition();
		Iterator<ArrayList<String>> iter = this.partition.iterator();

		while (iter.hasNext()) {
			ArrayList<String> block = iter.next();
			ArrayList<String> cloned_block = new ArrayList<String>(block.size());
			Iterator<String> iter_element = block.iterator();

			while (iter_element.hasNext()) {
				cloned_block.add(iter_element.next());
			}
			p.partition.add(cloned_block);
		}

		return p;
	}

	public boolean addBlock(ArrayList<String> block) {
		this.partition.add(block);
		return true;
	}
	public boolean contains(String v) {
		Iterator<ArrayList<String>> iter = this.partition.iterator();
		while (iter.hasNext()) {
			ArrayList<String> block = iter.next();
			if (block.contains(v))
				return true;
		}
		return false;
	}

	public ArrayList<String> getVertexList() {
		ArrayList<String> combined_list = new ArrayList<String>();
		Iterator<ArrayList<String>> iter = this.partition.iterator();
		while (iter.hasNext()) {
			
			Iterator<String> iter_node = iter.next().iterator();
			while(iter_node.hasNext()) {
				combined_list.add(iter_node.next());
			}
		}

		return combined_list;
	}

	public Partition nextMCMCPartition() {
		Partition p = this.clone();
		Random rng = new Random();

		int from_block = rng.nextInt(p.partition.size());
		while (p.partition.get(from_block).size() == 0)
			from_block = rng.nextInt(p.partition.size());

		if (p.partition.get(from_block).size() > 0) {
			int index_in_block = rng.nextInt(p.partition.get(from_block).size());

			int move_to_block = rng.nextInt(p.partition.size());
			while (move_to_block == from_block)
				move_to_block = rng.nextInt(p.partition.size());

			String vertex_to_move = (String) p.partition.get(from_block).remove(index_in_block);
			p.partition.get(move_to_block).add(vertex_to_move);
		}
		return p;

	}

	public String toString() {
		String str = "[";
		Iterator<ArrayList<String>> iter = partition.iterator();
		while (iter.hasNext()) {
			ArrayList<String> block = iter.next();
			str = str + "[";
			Iterator<String> iterElement = block.iterator();
			while (iterElement.hasNext()) {
				String s = iterElement.next();
				if(s.trim().isEmpty())
					continue;
				else
					str = str + s.trim();

				if (iterElement.hasNext())
					str = str + ",";
			}
			if (iter.hasNext())
				str = str + "], ";
			else
				str = str + "]]";
		}
		return str;
	}

	public void randomPartition(Graph<String, DefaultEdge> g, int K) {
		Set<String> vertices = g.vertexSet();
		ArrayList<String> vertexSet = new ArrayList<String>(vertices.size());

		Iterator<String> iter = vertices.iterator();

		while (iter.hasNext()) {
			vertexSet.add(iter.next());
		}
		randomPartition(vertexSet, K);
		calNode2Com();

	}

	public void calNode2Com() {
		this.node2com.clear();
        for(int com_index = 0; com_index<this.partition.size(); com_index++) {
        	ArrayList<String> com = this.partition.get(com_index); 
        	for (String node : com) {
        		this.node2com.put(node, com_index);
        	}
		}
	}
	public void randomPartition(ArrayList<String> vertexSet, int K) {
		Random rng = new Random();

		double total = 0.0;
		double[] r = new double[K];

		for (int i = 0; i < K; i++) {
			r[i] = rng.nextDouble();
			total += r[i];
		}

		int currentPos = 0;
		for (int step = 0; step < K; step++) {

			int blocksize = (int) Math.floor(r[step] / total * vertexSet.size());
			if (blocksize <=0)
				continue;
			if (step == K - 1)
				blocksize = vertexSet.size() - currentPos;

			// System.out.println(blocksize+ " "+ total);

			ArrayList<String> block = new ArrayList<String>(blocksize);

			for (int i = currentPos; i < currentPos + blocksize; i++) {
				block.add(vertexSet.get(i));
			}

			this.partition.add(block);

			currentPos += blocksize;

		}
		if(this.partition.size()<=1 && vertexSet.size()>=K) {
			int blocksize = Math.floorDiv(vertexSet.size(),K);
			currentPos = 0;
			for (int step = 0; step < K; step++) {
				if (step == K - 1)
					blocksize = vertexSet.size() - currentPos;

			// System.out.println(blocksize+ " "+ total);

				ArrayList<String> block = new ArrayList<String>(blocksize);

				for (int i = currentPos; i < currentPos + blocksize; i++) {
					block.add(vertexSet.get(i));
				}

				this.partition.add(block);

				currentPos += blocksize;

			}

		}
		
		this.calNode2Com();
	}

	public void sort() {
		Iterator<ArrayList<String>> iter = this.partition.iterator();
		while (iter.hasNext()) {
			iter.next().sort(new SortVertexIdentity());
		}

		this.partition.sort(new sortPartition());

	}
	
	public void toFile(String filename) throws IOException {
		FileWriter out = new FileWriter(filename);
		Iterator<ArrayList<String>> iter = this.partition.iterator();
		int com_id = 0;
	    while(iter.hasNext()) {
	    	String outString = "circle" + com_id +"\t";
	    	Iterator<String> iter_vertex = iter.next().iterator();
	    	
	    	while(iter_vertex.hasNext()) {
	    		outString.concat(iter_vertex.next()+"\t");
	    	}
	    	outString.concat("\n");
	    	out.write(outString);
	    	
	    	com_id ++;
	    }
	    out.close();
	}
}

class sortPartition implements Comparator<ArrayList<String>> {
	public int compare(ArrayList<String> list1, ArrayList<String> list2) {
		if (list1.size() == 0 || list2.size() == 0)
			return list1.size() - list2.size();

		//list1.sort(new SortVertexIdentity());
		//list2.sort(new SortVertexIdentity());

		String v1 = list1.get(0);
		String v2 = list2.get(0);
		
		return v1.compareTo(v2);

/*		if (v1.isEmpty() || v2.isEmpty())
			return v1.compareTo(v2);

		int id1 = Integer.parseInt(v1.substring(1));
		int id2 = Integer.parseInt(v2.substring(1));
		// System.out.println(id1+" "+id2);
		return id1 - id2;*/

	}
}

class SortVertexIdentity implements Comparator<String> {
	public int compare(String v1, String v2) {
		if (v1.isEmpty() || v2.isEmpty())
			return v1.compareTo(v2);

		//int id1 = Integer.parseInt(v1.substring(1));
		//int id2 = Integer.parseInt(v2.substring(1));
		// System.out.println(id1+" "+id2);

		return v1.compareTo(v2);
	}
}