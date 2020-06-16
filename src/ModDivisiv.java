import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.Set;

//@example:class:begin
public class ModDivisiv {
	// private static final int SIZE = 10;

	double global_sensitivity;
	double epsilon;
	//budget for bestcut
	double epsilon_m;
	double lambda;
	// maximum Level considered.
	int maxL;
	// number of branches, i.e., k
	int k;

	private final double DEFAULT_LAMBDA = 2.0;
	//number of iterations in a partition calculation 
	private int K; 
	
	private double w_s;

	public ModDivisiv(double epsilon, double epsilon_m, double global_sensitivity, int maxL, int k, int K) {
		this.epsilon = epsilon;
		this.global_sensitivity = global_sensitivity;
		this.maxL = maxL;
		this.k = k;
		this.lambda = DEFAULT_LAMBDA;
		this.epsilon_m =epsilon_m;
		this.K = K ;
	}

	public ModDivisiv(double epsilon, double epsilon_m, double global_sensitivity, int maxL, int k, int K, double w_s) {
		this.epsilon = epsilon;
		this.global_sensitivity = global_sensitivity;
		this.maxL = maxL;
		this.k = k;
		this.lambda = DEFAULT_LAMBDA;
		this.epsilon_m =epsilon_m;
		this.K = K ;
		this.w_s = w_s;
	}

	public HashMap<String, Double> getProbDistribution(HashMap<String, Integer> hmap) {
		HashMap<String, Double> result = new HashMap<String, Double>();
		Iterator<String> iter = hmap.keySet().iterator();
		int sum = 0;
		while (iter.hasNext()) {
			sum += (hmap.get(iter.next()));
		}

		iter = hmap.keySet().iterator();
		while (iter.hasNext()) {
			String key = iter.next();
			result.put(key, hmap.get(key) / (sum * 1.0));
		}

		return result;
	}
	
	public HashMap<Integer, Double> getProbDistribution_hash(HashMap<Integer, Integer> hmap) {
		HashMap<Integer, Double> result = new HashMap<Integer, Double>();
		Iterator<Integer> iter = hmap.keySet().iterator();
		int sum = 0;
		while (iter.hasNext()) {
			sum += (hmap.get(iter.next()));
		}

		iter = hmap.keySet().iterator();
		while (iter.hasNext()) {
			Integer key = iter.next();
			result.put(key, hmap.get(key) / (sum * 1.0));
		}

		return result;
	}


	public double getDPmass(double epsilon, double global_sensitivity, double score) {
		return Math.exp(epsilon * score / (0.2 * global_sensitivity));
	}

	public double getJumpProb(double epsilon, double global_sensitivity,  double score_incr) {
		return Math.min(1.0, Math.exp(epsilon * score_incr/(2*global_sensitivity)));
	}



/**
 * ModMCMC with attributes
 * @param p
 * @param g
 * @param windowsize
 * @param levelepsilon
 * @return
 */
	
	public Partition ModMCMC_attribute(Partition p, attributeGraph ag, graph g, double levelepsilon, double w_s, long rS) {
		
		double w_a = 1 - w_s;
		

		if (p.partition.size() < this.k)
			return p;

		double next_modularity = 0;

		System.out.println("Start computing partition paras...");               
		p.calculatePartitionParas(g);                                           
		System.out.println("Finish computing partition paras...");  	
		
		System.out.println("Start computing attribute partition paras...");               
		p.calculatePartitionParas_abt(ag);                                           
		System.out.println("Finish computing attribute partition paras...");  	
		
		if (p.nEdges == 0) 
			return p;
		
		if (p.nEdges_abt - 0<0.001) {
			return this.ModMCMC_v4(p, g, levelepsilon, rS);
		}
		
		int nEdges = p.nEdges;
		int nEdgesAbt = p.nEdges_abt;
		
		System.out.println("Start computing CurrModularity...");
		double curr_mod = g.getPartialModularity(p);
		double curr_modularity = g.getPartialModularity_abt(ag, p, w_s, w_a);
		
		System.out.println("o_mod =" + curr_mod + " curr_modularity="+curr_modularity);
		
		System.out.println("Finish computing CurrentModularity.");
	
		System.out.println("Start computing ComDegree ...");                    
		HashMap<Integer, Integer> ComDegree = new HashMap<Integer, Integer>();  
		HashMap<Integer, Integer> ComDegree_abt= new HashMap<Integer, Integer>();  
		for (int i = 0; i < p.partition.size(); i++) {                          
		     ArrayList<String> com = p.partition.get(i);                         
		     ComDegree.put(i, p.getSumOfDegrees(com));                           
		     ComDegree_abt.put(i, p.getSumOfDegrees_abt(com));                           

		 }                                                                       
		 System.out.println("Finish computing ComDegree ...");

		long count = 0 ;
		while (count < this.K * rS) {
			Random rng = new Random();
			int from_block = rng.nextInt(p.partition.size());
			while (p.partition.get(from_block).size() == 0)
			    from_block = rng.nextInt(p.partition.size());
			
			if (!p.partition.get(from_block).isEmpty()) {
				int index_in_block = rng.nextInt(p.partition.get(from_block).size());
				ArrayList<String> com_from =(ArrayList<String>) p.partition.get(from_block); 

				String vertex_to_remove = p.partition.get(from_block).get(index_in_block);
			    int[] degInCom = g.getNumberDegreeInComs(vertex_to_remove, p);
			    int[] degInComAbt = ag.getNumberDegreeInComs(vertex_to_remove, p);
			    			    
				double remove_cost = - (degInCom[from_block])
						+ (ComDegree.get(from_block) - p.degreeP.get(vertex_to_remove))*
						p.degreeP.get(vertex_to_remove)*1.0/(2 * nEdges);

				double remove_costAbt = - (degInComAbt[from_block])
						+ (ComDegree_abt.get(from_block) - p.degreeP_abt.get(vertex_to_remove))*
						p.degreeP_abt.get(vertex_to_remove)/(2 * nEdgesAbt);

				int move_to_block = rng.nextInt(p.partition.size());
				while (move_to_block == from_block)
					move_to_block = rng.nextInt(p.partition.size());

				ArrayList<String> com_to =(ArrayList<String>) p.partition.get(move_to_block); 

				double incr = remove_cost + (degInCom[move_to_block]) - 
						ComDegree.get(move_to_block)*p.degreeP.get(vertex_to_remove)*1.0/(2 * nEdges);
				
				double incrAbt = remove_costAbt + (degInComAbt[move_to_block]) - 
						ComDegree_abt.get(move_to_block)*p.degreeP_abt.get(vertex_to_remove)*1.0/(2 * nEdgesAbt);

				double increased_mod = w_s* incr/nEdges + w_a * incrAbt/nEdgesAbt;
				next_modularity = curr_modularity + increased_mod; 
				
				double next_mod = curr_modularity + incr/nEdges;
				
				double move_probability = this.getJumpProb(levelepsilon, this.global_sensitivity,  increased_mod);

				if (Math.random() < move_probability) {
				    com_from.remove(vertex_to_remove);	
				    com_to.add(vertex_to_remove);
					curr_modularity = next_modularity;
					curr_mod = next_mod;
					p.node2com.put(vertex_to_remove, move_to_block);

					ComDegree.put(from_block, ComDegree.get(from_block)-p.degreeP.get(vertex_to_remove));
					ComDegree.put(move_to_block, ComDegree.get(move_to_block)+p.degreeP.get(vertex_to_remove));

					ComDegree_abt.put(from_block, ComDegree_abt.get(from_block)-p.degreeP_abt.get(vertex_to_remove));
					ComDegree_abt.put(move_to_block, ComDegree_abt.get(move_to_block)+p.degreeP_abt.get(vertex_to_remove));
				}

				
				/* Sort partition so that the hash value of the partition is unique */
				if (count % 10000 == 0)
					System.out.println("count="+count+ ", current_modularity=" +
						 curr_modularity + " curr_mod="+curr_mod+" ,"+ "increased_mod=" + increased_mod + ", " + move_probability + " epsilon= " + this.epsilon +  " ws =" + this.w_s); 
			}
			count ++;
		}	
		return p;

	}
	/**
	 * divide the privacy budget given \lambda and \epsilon
	 */
	double[] dividePrivacyBudget(double epsilon) {
		double[] eA = new double[maxL];
		eA[maxL - 1] = epsilon * (1 - this.lambda) / (1 - Math.pow(this.lambda, this.maxL));
		for (int i = this.maxL - 2; i >= 0; i--) {
			eA[i] = this.lambda * eA[i + 1];
		}
		return eA;
	}

	/**
	 * The main function of ModDivisive
	 * @throws IOException 
	 */
	public Partition ModDivisiveMain(graph g, String dataName, int num_iter) throws IOException {
	    FileWriter out = new FileWriter(dataName+ "_" + this.epsilon + "_" + this.epsilon_m+"_" + this.maxL +"_"+ 
							this.k+"_" + this.K +"_partition.txt",true)	;
		double[] eA = dividePrivacyBudget(this.epsilon - this.maxL*this.epsilon_m);

		
		/* initiate root with all vertices in g*/
		MD_node root = new MD_node(0);
		Set<String> vertices = g.getGraph().vertexSet();
		Iterator<String> iter = vertices.iterator();
		while (iter.hasNext()) {
			root.getVertexList().add(iter.next());
		}
		
		/* initiate the queue with root */
		LinkedList<MD_node> queue = new LinkedList<MD_node>();
		queue.add(root);
		out.append(root.toString());
		out.close();
		System.out.println("First node added. Start the main loop ...");
		
		/* the main loop start here */
		while (!queue.isEmpty()) {
			
			MD_node r = queue.pollFirst();
			if(r.getVertexList().size()<30)
				continue;
			
			if (r.getLevel() < maxL) {
				Partition p = new Partition();
				
				System.out.println("  Start generating the first random partition ...");

				p.randomPartition(r.getVertexList(), this.k);

				System.out.println("  Finish generating the first random partition.");
				
				//System.out.println(r.getVertexList().toString()+ r.getVertexList().size());
				//System.out.println("The randomised partition: " + p.toString());
				
				/* Sample a new partition from a random one following differential privacy */

				//Partition new_partition = this.ModMCMC_v2(p, g, WINDOW_SIZE, eA[r.getLevel()],r.getVertexList().size());
				Partition new_partition = this.ModMCMC_v4(p, g,  eA[r.getLevel()],r.getVertexList().size());
				

				//System.out.println("level: "+r.getLevel() + " new_partition: " + new_partition.toString());

				/* Add the blocks in the new partition to the queue */
				Iterator<ArrayList<String>> iter_p = new_partition.partition.iterator();
				
				int count = 0;
				while (iter_p.hasNext()) {
					MD_node node = new MD_node(r.getLevel() + 1);
					node.setVertexList(iter_p.next());
					node.setParent(r);

					// set the parent of the node
					r.getChildren().add(node);
					// add the new node to the queue
					queue.add(node);

					/*out = new FileWriter(dataName+ "_" + this.epsilon + "_" + this.epsilon_m+"_" + this.maxL +"_"+ 
							this.k+"_" + this.K +"_partition.txt",true);
					out.write("count="+count+" " + node.toString());
					
					out.close();*/
				}
			}
		}
		// calculate the final partition with the last parameter as the privacy budget.
		Partition rp = bestCut(g, root);
		/*out = new FileWriter(dataName +"_" + this.epsilon + "_" + this.epsilon_m+"_" + this.maxL +"_"+ 
							this.k+"_" + this.K + "_"+ num_iter + "_best_cut.txt", true);
		out.write(rp.toString());
		out.close();*/
		return rp;
	}

	
	public Partition ModDivisiveMain_attribute(attributeGraph ag, graph g, String dataname, double w_s, int num_iter) throws IOException {
		if (Math.abs(w_s - 1.0) <0.00001)
			return this.ModDivisiveMain(g, dataname, num_iter);
		
		 FileWriter out = new FileWriter(dataname+ "_" + this.epsilon + "_" + this.epsilon_m+"_" + this.maxL +"_"+ 
					this.k+"_" + this.K +"_partition.txt",true);
				 
		double w_a = 1 - w_s;
		/* Calculate the privacy budget on each level */
		double[] eA = dividePrivacyBudget(this.epsilon - this.maxL*this.epsilon_m);

		/* initiate the queue by adding the root node with all vertices in g*/
		MD_node root = new MD_node(0);
		Set<String> vertices = g.getGraph().vertexSet();
		Iterator<String> iter = vertices.iterator();
		while (iter.hasNext()) {
			root.getVertexList().add(iter.next());
		}
		LinkedList<MD_node> queue = new LinkedList<MD_node>();
		queue.add(root);

		while (!queue.isEmpty()) {
			MD_node r = queue.pollFirst();
			if(r.getVertexList().size()<30)
				continue;
			if (r.getLevel() < maxL) {
				Partition p = new Partition();
				
				System.out.println("  Start generating the first random partition ...");
				p.randomPartition(r.getVertexList(), this.k);
				System.out.println("  Finish generating the first random partition.");

				Partition new_partition = this.ModMCMC_attribute(p,ag, g, eA[r.getLevel()], 
						w_s, r.getVertexList().size());

				//System.out.println("level: "+r.getLevel() + " new_partition: " + new_partition.toString());

				Iterator<ArrayList<String>> iter_p = new_partition.partition.iterator();
				while (iter_p.hasNext()) {
					MD_node node = new MD_node(r.getLevel() + 1);
					node.setVertexList(iter_p.next());
					node.setParent(r);

					// set the parent of the node
					r.getChildren().add(node);
					// add the new node to the queue
					queue.add(node);
				}
			}
		}
		// calculate the final partition with the last parameter as the privacy budget.
		System.out.println("Start calculating best_cut ... ");
		Partition rp = bestCut(ag, g, root, w_s, w_a);
		out.write(rp.toString());
		out.close();
		return rp;

		
	}

	/**
	 * calculate the best cut given the k-nary tree and the graph g queue: <-
	 * xxxxxx<- stack: xxxxxx<->
	 **/
	public Partition bestCut(graph g, MD_node r) {
		LinkedList<MD_node> stack = new LinkedList<MD_node>();
		LinkedList<MD_node> queue = new LinkedList<MD_node>();

		// get the stack with leaf nodes closer to the top
		queue.addLast(r);
		while (!queue.isEmpty()) {
			MD_node node = queue.removeFirst();
			stack.addLast(node);
			Iterator<MD_node> iter = node.getChildren().iterator();

			while (iter.hasNext()) {
				queue.addLast(iter.next());
			}
		}

		// Create the laplacian noise generator for modularity noise
		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, this.global_sensitivity / epsilon_m);

		while (!stack.isEmpty()) {
			MD_node node = stack.removeLast();
			// calcualte the modularity of the block defined in node
			ArrayList<ArrayList<String>> temp = new ArrayList<ArrayList<String>>();
			temp.add(node.getVertexList());
			Partition ptn = new Partition(temp);
			double mod = g.getModularity(ptn) + lng.nextLaplacian();

			/* update the modularity in the node */
			node.setModularity(mod);

			/**
			 * compare the node's modularity with the sum of all its children's modularity
			 * determine whether the block is selected as a block.
			 */

			if (node.getChildren().isEmpty()) {
				node.setSelf(true);
			} else {
				/**
				 * calculate the sum of modularity of all the children of the node
				 */
				double sum_m = 0;
				Iterator<MD_node> iter = node.getChildren().iterator();
				while (iter.hasNext())
					sum_m += iter.next().getModularity();

				if (node.getModularity() < sum_m) {
					node.setSelf(false);
				} else
					node.setSelf(true);
			}
		}

		// transform the result to a partition
		Partition p = new Partition();

		queue.clear();

		queue.add(r);
		while (!queue.isEmpty()) {
			MD_node node = queue.removeFirst();
			if (node.getSelf() == true) {
				p.addBlock(node.getVertexList());
			} else {
				Iterator<MD_node> iter = node.getChildren().iterator();
				while (iter.hasNext())
					queue.addLast(iter.next());
			}
		}
		return p;
	}

	
	/**
	 * bestcut with attributes
	 */
	public Partition bestCut(attributeGraph ag, graph g, MD_node r, double w_s, double w_a) {
		LinkedList<MD_node> stack = new LinkedList<MD_node>();
		LinkedList<MD_node> queue = new LinkedList<MD_node>();

		// get the stack with leaf nodes closer to the top
		queue.addLast(r);
		while (!queue.isEmpty()) {
			MD_node node = queue.removeFirst();
			stack.addLast(node);
			Iterator<MD_node> iter = node.getChildren().iterator();

			while (iter.hasNext()) {
				queue.addLast(iter.next());
			}
		}

		// Create the laplacian noise generator for modularity noise
		LaplacianNoiseGenerator lng = new LaplacianNoiseGenerator(0, this.global_sensitivity / epsilon_m);

		while (!stack.isEmpty()) {
			MD_node node = stack.removeLast();
			// calcualte the modularity of the block defined in node
			ArrayList<ArrayList<String>> temp = new ArrayList<ArrayList<String>>();
			temp.add(node.getVertexList());
			Partition ptn = new Partition(temp);
			double mod = g.getModularityAttribute(ag, ptn, w_s, w_a) + lng.nextLaplacian();
			
			//double mod1 = g.getModularity(ptn);
			
			//System.out.println(mod+ " " + mod1);

			/* update the modularity in the node */
			node.setModularity(mod);

			/**
			 * compare the node's modularity with the sum of all its children's modularity
			 * determine whether the block is selected as a block.
			 */

			if (node.getChildren().isEmpty()) {
				node.setSelf(true);
			} else {
				/**
				 * calculate the sum of modularity of all the children of the node
				 */
				double sum_m = 0;
				Iterator<MD_node> iter = node.getChildren().iterator();
				while (iter.hasNext())
					sum_m += iter.next().getModularity();

				if (node.getModularity() < sum_m) {
					node.setSelf(false);
				} else
					node.setSelf(true);
			}
		}

		// transform the result to a partition
		Partition p = new Partition();

		queue.clear();

		queue.add(r);
		while (!queue.isEmpty()) {
			MD_node node = queue.removeFirst();
			if (node.getSelf() == true) {
				p.addBlock(node.getVertexList());
			} else {
				Iterator<MD_node> iter = node.getChildren().iterator();
				while (iter.hasNext())
					queue.addLast(iter.next());
			}
		}
		return p;
	}

	public String updatePDigest(String pDigest, String vertex_to_remove, String com_to_item0) {
		int vtr_pos = pDigest.indexOf("," + vertex_to_remove + ",");
		if (vtr_pos >= 0) {
			pDigest = pDigest.substring(0, vtr_pos + 1) + pDigest.substring(vtr_pos + vertex_to_remove.length() + 2);
			int pos_bracket = pDigest.indexOf("]", pDigest.indexOf(com_to_item0));
			pDigest = pDigest.substring(0, pos_bracket) + "," + vertex_to_remove + pDigest.substring(pos_bracket);
		} else {
			vtr_pos = pDigest.indexOf("," + vertex_to_remove + "]");
			if (vtr_pos >= 0) {
				pDigest = pDigest.substring(0, vtr_pos) + pDigest.substring(vtr_pos + vertex_to_remove.length() + 1);
				int pos_bracket = pDigest.indexOf("]", pDigest.indexOf(com_to_item0));
				pDigest = pDigest.substring(0, pos_bracket) + "," + vertex_to_remove + pDigest.substring(pos_bracket);
			} else {
				vtr_pos = pDigest.indexOf("[" + vertex_to_remove + ",");
				if (vtr_pos >= 0) {
					pDigest = pDigest.substring(0, vtr_pos + 1)
							+ pDigest.substring(vtr_pos + 1 + vertex_to_remove.length() + 1);
					int pos_bracket = pDigest.indexOf("]", pDigest.indexOf(com_to_item0));
					pDigest = pDigest.substring(0, pos_bracket) + "," + vertex_to_remove
							+ pDigest.substring(pos_bracket);
				} else {
					vtr_pos = pDigest.indexOf("[" + vertex_to_remove + "]");
					if (vtr_pos >= 0) {
						pDigest = pDigest.substring(0, vtr_pos + 1)
								+ pDigest.substring(vtr_pos + 1 + vertex_to_remove.length());
						int pos_bracket = pDigest.indexOf("]", pDigest.indexOf(com_to_item0));
						pDigest = pDigest.substring(0, pos_bracket) + "," + vertex_to_remove
								+ pDigest.substring(pos_bracket);
					} else {
						System.out.println("error");
					}
				}
			}
		}
		return pDigest; 
	}

	

	public Partition ModMCMC_v4(Partition p, graph g,  double levelepsilon, long rS) {
		
		if(p.partition.size()<this.k)
			return p;
		/* hist stores the histogram of samples generated */


		/* lastMeanModularity stores the mean modularity in the last window
		 * and meanDifference stores the difference of mean modularities between current and the last window
		 */
		double next_modularity = 0;
		

		System.out.println("Start computing partition paras...");
		p.calculatePartitionParas(g);
		System.out.println("Finish computing partition paras...");

		//int nEdges = g.getNumberOfEdges(p.getVertexList());
		if (p.nEdges == 0)
			return p;
		int nEdges = p.nEdges;
		System.out.println("nEdges = " + nEdges);

		//int nEdges = g.getGraph().edgeSet().size();

		System.out.println("Start computing CurrModularity...");

		double curr_modularity = g.getPartialModularity(p);
//		double curr_modularity_1 = g.getModularity(p);
		
		if (curr_modularity >1 || curr_modularity <-1)
			System.out.println("error");

		long count = 0 ;
		
		System.out.println("Start computing ComDegree ...");
	    HashMap<Integer, Integer> ComDegree = new HashMap<Integer, Integer>();
	    for (int i = 0; i < p.partition.size(); i++) {
	    	ArrayList<String> com = p.partition.get(i);
	    	ComDegree.put(i, p.getSumOfDegrees(com));
	    }
		System.out.println("Finish computing ComDegree ...");
	    
	    
		while (count < this.K * rS) {
			
				/*
				 * get the vertex to move and the destination block
				 */
			Random rng = new Random();

			int from_block = rng.nextInt(p.partition.size());
			while (p.partition.get(from_block).size() == 0)
			    from_block = rng.nextInt(p.partition.size());
			
			if (!p.partition.get(from_block).isEmpty()) {
				int index_in_block = rng.nextInt(p.partition.get(from_block).size());
				ArrayList<String> com_from =(ArrayList<String>) p.partition.get(from_block); 

				String vertex_to_remove = p.partition.get(from_block).get(index_in_block);
			    int[] degInCom = g.getNumberDegreeInComs(vertex_to_remove, p);
				
			    /*
				double remove_cost = - (degInCom[from_block])
						+ (ComDegree.get(from_block) - g.getGraph().degreeOf(vertex_to_remove))*
						g.getGraph().degreeOf(vertex_to_remove)/(2 * nEdges);
					*/
				double remove_cost = - (degInCom[from_block])
						+ (ComDegree.get(from_block) - p.degreeP.get(vertex_to_remove))*
						p.degreeP.get(vertex_to_remove)*1.0/(2 * nEdges);


				int move_to_block = rng.nextInt(p.partition.size());
				while (move_to_block == from_block)
					move_to_block = rng.nextInt(p.partition.size());

				ArrayList<String> com_to =(ArrayList<String>) p.partition.get(move_to_block); 

				/*
				double incr = remove_cost + (degInCom[move_to_block]) - 
						ComDegree.get(move_to_block)*g.getGraph().degreeOf(vertex_to_remove)/(2 * nEdges);
			*/	
				
				double incr = remove_cost + (degInCom[move_to_block]) - 
						ComDegree.get(move_to_block)*p.degreeP.get(vertex_to_remove)*1.0/(2 * nEdges);
				
				next_modularity = curr_modularity + incr/nEdges; 
				
				double move_probability = this.getJumpProb(levelepsilon, this.global_sensitivity, incr/nEdges);
				if (Math.random() < move_probability) {
				    com_from.remove(vertex_to_remove);	
				    com_to.add(vertex_to_remove);
					curr_modularity = next_modularity;
					p.node2com.put(vertex_to_remove, move_to_block);
					ComDegree.put(from_block, ComDegree.get(from_block)-p.degreeP.get(vertex_to_remove));
					ComDegree.put(move_to_block, ComDegree.get(move_to_block)+p.degreeP.get(vertex_to_remove));
					
				}

				
				/* Sort partition so that the hash value of the partition is unique */
				
				if (count % 10000 == 0)
				 System.out.println("count="+count+ ", current_modularity=" +
						 curr_modularity + " , incr="+ incr/nEdges + ", " + move_probability + " vertex_to_remove=" + vertex_to_remove); 
				 //+ "\n 		 Current partition is " + current_p.toString());
			}
			count ++;
			
		}
		
		/* sample a new partition from the distribution calculated from hist */
	// convert sample_partition which is a string to a partition
		return p;

	}
}

