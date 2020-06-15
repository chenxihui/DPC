import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import org.jgrapht.alg.scoring.ClusteringCoefficient;
import org.jgrapht.graph.DefaultEdge;

public class reIdentificationAttack {
	public graph gOri;
	public graph gGen;
	double[][] weightMgOri;
	double[][] weightMgGen;
	
	public void attackInBatchForExperiment(String folderName, String datasetName, String e) throws IOException {

		String originalFileName;
		if(datasetName.contentEquals("facebook"))
			originalFileName = datasetName + "_combined_no_ego.txt";
		else
			originalFileName = datasetName + "_combined.txt";
				
		String originalPartName = "partition_" + datasetName + "_combined.txt";
	
		File dir = new File(folderName);
		  File[] directoryListing = dir.listFiles();
		  if (directoryListing != null) {
		    for (File partfile: directoryListing) {
		    	
		    	if(!partfile.getName().contains("parition_" + datasetName)||!partfile.getName().contains(e))
		    		continue;
		    	String partGenFile = folderName + File.separator + partfile.getName();
		    	String genGraphFileName = folderName + File.separator + 
		    			partfile.getName().substring(partfile.getName().indexOf("_")+1);
		    	
		    	System.out.println("Start attacking file " + genGraphFileName);
		    	reIdentificationAttack attack = new reIdentificationAttack(originalFileName, genGraphFileName);
		    	String outMapFile = folderName + File.separator + "vertexMapResult/Map_" + 
		    			partfile.getName().substring(partfile.getName().indexOf("_")+1);
		    	if (new File(outMapFile).exists())
		    		continue;
		    	FileWriter outMap = new FileWriter(folderName + File.separator + "vertexMapResult/Map_" + 
		    			partfile.getName().substring(partfile.getName().indexOf("_")+1));

		    	HashMap<String, String> vertexMap = attack.theAttack(originalPartName, partGenFile);
		    	for (String v : vertexMap.keySet()) {
		    		outMap.write(v + " " + vertexMap.get(v)+ "\n");
		    	}
		    	outMap.close();
		    }
		  } 
	}
	public void attackSingleFile(String datasetName, String genGraphFileName ) throws IOException {

		String originalFileName;
		if(datasetName.contentEquals("facebook"))
			originalFileName = datasetName + "_combined_no_ego.txt";
		else
			originalFileName = datasetName + "_combined.txt";
				
		String originalPartName =  datasetName + "PartGen.txt";
	
		String partGenFile = "parition_"+ genGraphFileName;
		reIdentificationAttack attack = new reIdentificationAttack(originalFileName, genGraphFileName);
		String outMapFile = "vertexMapResult/Map_" + genGraphFileName;
		FileWriter outMap = new FileWriter(outMapFile);
		HashMap<String, String> vertexMap = attack.theAttack(originalPartName, partGenFile);
		for (String v : vertexMap.keySet()) {
		    outMap.write(v + " " + vertexMap.get(v)+ "\n");
		}
		outMap.close();
		
		double successrate = 0; 
	    for (String v : vertexMap.keySet()) {
	    	if (v.contentEquals(vertexMap.get(v)))
	    		successrate +=1.0d;
	    }
	    successrate = successrate/attack.gOri.getGraph().vertexSet().size();
	    System.out.println("The success rate is " + successrate);
	}

	reIdentificationAttack(){
		
	}
	reIdentificationAttack(String graphFileOri, String graphFileGen) throws IOException{
		gOri = new graph();
		gOri.readGraphOnly(graphFileOri);

		gGen = new graph();
		gGen.readGraphOnly(graphFileGen);
	}

	HashMap<String, String> theAttack(String oriPartFile, String genPartFile) 
			throws IOException {
		BufferedReader fin = new BufferedReader(new FileReader(oriPartFile));
		String line = fin.readLine();
		Partition partOri = new Partition(line);
		fin.close();
		
		fin = new BufferedReader(new FileReader(genPartFile));
		line = fin.readLine();
		Partition partGen = new Partition(line);
		fin.close();
		
		partOri.calNode2Com();
		partGen.calNode2Com();
		
		ArrayList<String> seedlist = this.getSeedFrom4Cliques(16);
		HashMap<Integer, Integer> comMap = this.communityMapping(partOri, partGen, seedlist);
		
		weightMgOri = new double[partOri.partition.size()][partOri.partition.size()];
		weightMgGen = new double[partGen.partition.size()][partGen.partition.size()];
		graph mgOri= getMacroGraph(partOri, "Original");
		graph mgGen= getMacroGraph(partGen, "Generated");
		
		
		int num_iteration = 0;
		while(num_iteration< 5 * Math.max(mgOri.getGraph().vertexSet().size(), mgGen.getGraph().vertexSet().size())) {
			num_iteration ++ ;
			Random rng = new Random();
			int comPos = rng.nextInt(comMap.size());
			Iterator<Integer> iter = comMap.keySet().iterator();
			int i = 0 ;
			int comOri=0; 
			while(iter.hasNext() && i<comPos) {
				comOri= iter.next();
				i++;
			}

			if (!comMap.containsKey(comOri))
				continue;
			int comGen = comMap.get(comOri);
			
			if(mgOri.getGraph().degreeOf(comOri+"")>0) {
				Iterator<DefaultEdge> iterE = mgOri.getGraph().edgesOf(comOri+"").iterator();
				int rngNei = -1;
				while(iterE.hasNext() ) {
					DefaultEdge e = iterE.next();
					
					rngNei = Integer.parseInt(Integer.toString(comOri).contentEquals(mgOri.getGraph().getEdgeSource(e))?mgOri.getGraph().getEdgeTarget(e):
						mgOri.getGraph().getEdgeSource(e));
					if(!comMap.containsKey(rngNei)) 
						break;
				}
				if(rngNei == -1 || comMap.containsKey(rngNei))
					continue;
				double maxScore = 0;
				int comGenMapped = -1;
				for(DefaultEdge e : mgGen.getGraph().edgesOf(comGen+"")) {
					int neiGen = Integer.parseInt(Integer.parseInt(mgGen.getGraph().getEdgeSource(e))== comGen? 
							mgGen.getGraph().getEdgeTarget(e):mgGen.getGraph().getEdgeSource(e));
					if(comMap.containsValue(neiGen))
						continue;
					double s = getSimScore(rngNei,neiGen , mgOri, mgGen, comMap);
					if(maxScore < s) {
						comGenMapped = neiGen;
						maxScore = s;
					}
				}
				if(comGenMapped!= -1 && Double.compare(maxScore, 0) > 0 ) {
					comMap.put(rngNei, comGenMapped);
				}
			}
		}
		
		
		HashMap<String, String> vertexMap = new HashMap<String, String>();
		for (String v : seedlist)
			vertexMap.put(v, v);
		System.out.println("Start mapping in individual communities ... ");
		for (int comOri : comMap.keySet()) {

			int comGen = comMap.get(comOri);
			ArrayList<String> vComOri = partOri.partition.get(comOri);
			ArrayList<String> vComGen = partGen.partition.get(comGen);
			

			System.out.println("Identifying seeds in community " + comOri);
			identifySeedInCom(vComOri, vComGen, vertexMap);
			System.out.println("Propagating in community " + comOri);
			propagateInCom(vComOri, vComGen, vertexMap);
		}
		
		System.out.println("Start global propagation ...");
		propagationGlobal(vertexMap);
		System.out.println(vertexMap.toString());
		return vertexMap;
	}
	

	
	public HashMap<Integer, Integer> communityMapping(Partition partOri, Partition partGen, ArrayList<String> seedlist) throws IOException {
		
		HashMap<Integer, Integer> comMap = new HashMap<Integer, Integer>();

		for(int com_index = 0; com_index< partOri.partition.size(); com_index ++) {
			ArrayList<String> seedinCom = new ArrayList<String>();
			HashMap<Integer, Integer> comCount = new HashMap<Integer, Integer> ();
			for (String v: seedlist) {
				if(partOri.node2com.get(v) == com_index){
					seedinCom.add(v);
					int mapCom = partGen.node2com.get(v);
					if (comMap.keySet().contains(mapCom))
						continue;
					if(comCount.containsKey(mapCom)) {
						int m = comCount.get(mapCom) + 1 ;
						comCount.put(mapCom, m);
					}else
						comCount.put(mapCom, 1);
				}
			}
			
			int selectCom = -1;
			if(comCount.size()>0) {
				int maxMapping = 0;
				for (int i : comCount.keySet()) {
					if (maxMapping < comCount.get(i)) {
						maxMapping = comCount.get(i);
						selectCom = i;
					}
				}

			}
			if (selectCom >-1)
				comMap.put(com_index, selectCom);
		}
		
		return comMap;
	}
	
		private void propagationGlobal(HashMap<String, String> vertexMap) {
		// TODO Auto-generated method stub

		int num_iteration = 0;

		while(num_iteration< 2 * gOri.getGraph().vertexSet().size()) {
			num_iteration ++ ;
			Random rng = new Random();
			int pos = rng.nextInt(vertexMap.size());
			Iterator<String> iter = vertexMap.keySet().iterator();
			int i = 0 ;
			String vertexOri = ""; 
			while(iter.hasNext() && i<pos) {
				vertexOri= iter.next();
				i++;
			}

			String vertexGen = vertexMap.get(vertexOri);
			
			if(!gOri.getGraph().vertexSet().contains(vertexOri)) {
				System.out.println("Vertex Ori not existed during global propagation: " + vertexOri);
				continue;
			}
			
			if(gOri.getGraph().degreeOf(vertexOri)>0) {
				Iterator<DefaultEdge> iterE = gOri.getGraph().edgesOf(vertexOri).iterator();
				String rngNei = "";
				while(iterE.hasNext() && (vertexMap.keySet().contains(rngNei))) {
					DefaultEdge e = iterE.next();
					rngNei = vertexOri.contentEquals(gOri.getGraph().getEdgeSource(e))?
							gOri.getGraph().getEdgeTarget(e): 
								gOri.getGraph().getEdgeSource(e);
				}
				if(rngNei.isEmpty()  || vertexMap.keySet().contains(rngNei))
					continue;
				double maxScore = 0;
				String vertexGenMapped = "";
				for(DefaultEdge e : gGen.getGraph().edgesOf(vertexGen)) {
					String neiGen = gGen.getGraph().getEdgeSource(e).contentEquals(vertexGen)? 
							gGen.getGraph().getEdgeTarget(e):gGen.getGraph().getEdgeSource(e);
					if(vertexMap.containsValue(neiGen))
						continue;
					double s = getSimScore(rngNei,neiGen,  vertexMap);
					 
					if(maxScore < s) {
						vertexGenMapped = neiGen;
						maxScore = s;
					}
				}
				if(vertexGenMapped.isEmpty() && Double.compare(maxScore, 0) > 0 ) {
					vertexMap.put(rngNei, vertexGenMapped);
				}
			}
		}
	}


	private double getSimScore(String rngNei, String neiGen, HashMap<String, String> vertexMap) {
		// TODO Auto-generated method stub

		double s = 0;
	
		ArrayList<String> neighbourOfNeiGen = new ArrayList<String> ();
		
		for(DefaultEdge e : gGen.getGraph().edgesOf(neiGen)) {
			String neiStr = gGen.getGraph().getEdgeSource(e).contentEquals(neiGen)?
					gGen.getGraph().getEdgeTarget(e): gGen.getGraph().getEdgeSource(e);
			
			neighbourOfNeiGen.add(neiStr);
		}
		
		for(DefaultEdge e : gOri.getGraph().edgesOf(rngNei)) {
			String neiStr = gOri.getGraph().getEdgeSource(e).contentEquals(rngNei)?
					gOri.getGraph().getEdgeTarget(e): gOri.getGraph().getEdgeSource(e);

			if(vertexMap.containsKey(neiStr)) {
				String mappedVertex= vertexMap.get(neiStr) ;
				if(neighbourOfNeiGen.contains(mappedVertex)) {
					s = s + 1 ;
				}
			}
		}
		
		s = s/Math.sqrt(gOri.getGraph().degreeOf(rngNei) * gGen.getGraph().degreeOf(neiGen));
		
		return s;
	}


	private void propagateInCom(ArrayList<String> vComOri, ArrayList<String> vComGen,
			HashMap<String, String> vertexMap) {

		int num_iteration = 0;
		ArrayList<String> mappedVertices = new ArrayList<String> ();
		for(String v: vertexMap.keySet()) {
			if(vComOri.contains(v))
				mappedVertices.add(v);
		}

		while(num_iteration< 2 * Math.max(vComOri.size(), vComGen.size())) {
			num_iteration ++ ;
			Random rng = new Random();
			int pos = rng.nextInt(mappedVertices.size());
			//Iterator<String> iter = mappedVertices.iterator();
			//int i = 0 ;
			String vertexOri = mappedVertices.get(pos); 
			//while(iter.hasNext() && i<pos) {
			//	vertexOri= iter.next();
			//	i++;
			//}

			String vertexGen = vertexMap.get(vertexOri);
			
			if(!gOri.getGraph().vertexSet().contains(vertexOri)) {
				System.out.println("Vertex Ori not existed: " + vertexOri);
				continue;
			}
				
			if(gOri.getGraph().degreeOf(vertexOri)>0) {
				Iterator<DefaultEdge> iterE = gOri.getGraph().edgesOf(vertexOri).iterator();
				String rngNei = "";
				while(iterE.hasNext() ) {
					DefaultEdge e = iterE.next();
					rngNei = vertexOri.contentEquals(gOri.getGraph().getEdgeSource(e))?
							gOri.getGraph().getEdgeTarget(e): 
								gOri.getGraph().getEdgeSource(e);
					if((vertexMap.containsKey(rngNei)|| !vComOri.contains(rngNei)))
						continue;
				}
				if(rngNei.isEmpty()  || vertexMap.containsKey(rngNei))
					continue;
				double maxScore = 0;
				String vertexGenMapped = "";
				for(DefaultEdge e : gGen.getGraph().edgesOf(vertexGen)) {
					String neiGen = gGen.getGraph().getEdgeSource(e).contentEquals(vertexGen)? 
							gGen.getGraph().getEdgeTarget(e):gGen.getGraph().getEdgeSource(e);
					if(vertexMap.containsValue(neiGen)||!vComGen.contains(neiGen))
						continue;
					double s = getSimScoreInCom(rngNei,neiGen , vComOri, vComGen, vertexMap);
					 
					if(maxScore < s) {
						vertexGenMapped = neiGen;
						maxScore = s;
					}
				}
				if(vertexGenMapped.isEmpty() && Double.compare(maxScore, 0) > 0 ) {
					vertexMap.put(rngNei, vertexGenMapped);
					mappedVertices.add(rngNei);
				}
			}
		}
	
	}


	private double getSimScoreInCom(String rngNei, String neiGen,ArrayList<String> vComOri,
			ArrayList<String> vComGen,
			HashMap<String, String> vertexMap) {

		double s = 0;
	
		ArrayList<String> neighbourOfNeiGen = new ArrayList<String> ();
		
		for(DefaultEdge e : gGen.getGraph().edgesOf(neiGen)) {
			String neiStr = gGen.getGraph().getEdgeSource(e).contentEquals(neiGen)?
					gGen.getGraph().getEdgeTarget(e): gGen.getGraph().getEdgeSource(e);
			
			neighbourOfNeiGen.add(neiStr);
		}
		
		for(DefaultEdge e : gOri.getGraph().edgesOf(rngNei)) {
			String neiStr = gOri.getGraph().getEdgeSource(e).contentEquals(rngNei+"")?
					gOri.getGraph().getEdgeTarget(e): gOri.getGraph().getEdgeSource(e);

			if(vertexMap.containsKey(neiStr)) {
				String mappedVertex= vertexMap.get(neiStr) ;
				if(neighbourOfNeiGen.contains(mappedVertex)) {
					s = s + 1 ;
				}
			}
		}
		
		s = s/Math.sqrt(gOri.getGraph().degreeOf(rngNei) * gGen.getGraph().degreeOf(neiGen));
		
		return s;
	}


	void identifySeedInCom(ArrayList<String> vComOri, ArrayList<String> vComGen, HashMap<String, String> vertexMap) {
		
		ClusteringCoefficient<String, DefaultEdge> gCCOri = new ClusteringCoefficient(gOri.getGraph());
		ClusteringCoefficient<String, DefaultEdge> gCCGen = new ClusteringCoefficient(gGen.getGraph());
		
		for(String v1: vComOri) {
			if(vertexMap.keySet().contains(v1))
				continue;
			ArrayList<Double> degL= new ArrayList<Double>();
			ArrayList<Double> ccL= new ArrayList<Double>();
			scoreEntity maxDegScore = new scoreEntity("", 0); 
			scoreEntity sndDegScore = new scoreEntity("", 0);
			scoreEntity maxCCScore = new scoreEntity("", 0); 
			scoreEntity sndCCScore = new scoreEntity("", 0);

			for(String v2: vComGen) {
				if(vertexMap.values().contains(v2))
					continue;
				double degScore = 1-getDegreeScore(v1, v2);
				double ccScore = 1-getCCScore(gCCOri.getVertexScore(v1), gCCGen.getVertexScore(v2));
				
				degL.add(degScore);
				ccL.add(ccScore);

				if(maxDegScore.getScore()< degScore) {
					sndDegScore.setVertex(maxDegScore.getVertex());
					sndDegScore.setScore(maxDegScore.getScore());
					maxDegScore.setVertex(v2);
					maxDegScore.setScore(degScore);
				}
				if(maxCCScore.getScore()< ccScore) {
					sndCCScore.setVertex(maxCCScore.getVertex());
					sndCCScore.setScore(maxCCScore.getScore());
					maxCCScore.setVertex(v2);
					maxCCScore.setScore(ccScore);
				}
			}
			
			double degVar = variance(degL);
			double ccVar = variance(ccL);
			
			double Dd = Math.abs(maxDegScore.getScore() - sndDegScore.getScore())/degVar;
			double Dcc = Math.abs(maxCCScore.getScore() - sndCCScore.getScore())/ccVar;
			
			if(Dd >= 0.1)
				vertexMap.put(v1,maxDegScore.getVertex());
			else if (Dcc >= 0.1)
				vertexMap.put(v1, maxCCScore.getVertex());
		}
	}
	
	private double getCCScore(double cc1, double cc2) {
		
		return Math.abs(cc1-cc2)/Math.max(cc1, cc2);
	}


	private double getDegreeScore(String v1, String v2) {
		double deg1 = this.gOri.getGraph().degreeOf(v1);
		double deg2 = this.gGen.getGraph().degreeOf(v2);
		
		return Math.abs(deg1-deg2)/Math.max(deg1, deg2);
	}


	double getSimScore(int rngNei, int neiGen, graph mgOri, graph mgGen, HashMap<Integer, Integer> comMap) {
		double s = 0;
	
		ArrayList<String> neighbourOfNeiGen = new ArrayList<String> ();
		
		for(DefaultEdge e : mgGen.getGraph().edgesOf(neiGen+"")) {
			String neiStr = mgGen.getGraph().getEdgeSource(e).contentEquals(neiGen+"")?mgGen.getGraph().getEdgeTarget(e):
				mgGen.getGraph().getEdgeSource(e);
			neighbourOfNeiGen.add(neiStr);
		}
		
		for(DefaultEdge e : mgOri.getGraph().edgesOf(rngNei+"")) {
			String neiStr = mgOri.getGraph().getEdgeSource(e).contentEquals(rngNei+"")?mgOri.getGraph().getEdgeTarget(e):
				mgOri.getGraph().getEdgeSource(e);
			if(comMap.containsKey(Integer.parseInt(neiStr))) {
				String mappedCom = Integer.toString(comMap.get(Integer.parseInt(neiStr)));
				if(neighbourOfNeiGen.contains(mappedCom)) {
					double w1 = weightMgOri[Integer.parseInt(neiStr)][rngNei];
					//		mgOri.getGraph().getEdgeWeight(mgOri.getGraph().getEdge(neiStr, rngNei+""));
					double w2 = weightMgGen[Integer.parseInt(mappedCom)][neiGen];
					//		mgGen.getGraph().getEdgeWeight(mgGen.getGraph().getEdge(mappedCom, neiGen+""));
					s = s + (1 - Math.sqrt(w1) - Math.sqrt(w2));
				}
			}
		}
		
		s = s/Math.sqrt(mgOri.getGraph().degreeOf(rngNei+"") * mgGen.getGraph().degreeOf(neiGen + ""));
		
		return s;
	}
	
	graph getMacroGraph(Partition part, String graphType) {
		
		graph g;
		double [][] weight;
		if(graphType.contentEquals("Original")) {
			g = gOri;
			weight = this.weightMgOri;
		} else {
			g = gGen;
			weight = this.weightMgGen;
		}

		graph mg = new graph();
		
		for (int i=0; i<part.partition.size(); i++) {
			mg.getGraph().addVertex(""+i);
		}
		
		for (DefaultEdge e: g.getGraph().edgeSet()) {
			String v1 = g.getGraph().getEdgeSource(e);
			String v2 = g.getGraph().getEdgeTarget(e);
			
			int com1 = part.node2com.get(v1);
			int com2 = part.node2com.get(v2);
			
			if(com1 != com2) {
				if (!mg.getGraph().containsEdge(com1+"", com2+"")) {
					mg.getGraph().addEdge(com1 + "", com2 + "");
					weight[com1][com2] = 1;
					weight[com2][com1] = 1;
					//mg.getGraph().setEdgeWeight(com1+"", com2+"", 1);
				}else {
					double preWeight = mg.getGraph().getEdgeWeight(mg.getGraph().getEdge(com1+"", com2+""));
					weight[com1][com2] ++;
					weight[com2][com1] ++;
					//mg.getGraph().setEdgeWeight(com1+"", com2+"", preWeight + 1);
				}
			}
		}
		for(int i = 0 ; i< weight.length; i++) {
			int sum = 0 ;
			for (int j=0; j<weight.length; j++)
				sum += weight[i][j];
			
			for (int j=0; j<weight.length; j++)
				 weight[i][j]/= sum;
		}
		
		return mg;
	}
	
	ArrayList<String> getSeedFrom4Cliques(int num_clique){
		graph g = gOri;

		ArrayList<String> nodeList = new ArrayList<String>(g.getGraph().vertexSet());
		ArrayList<String> seedlist = new ArrayList<String>();
		Random rng = new Random();
		
		int numLoop=0;
	
		int totalClique = 0;
		while(totalClique < num_clique && numLoop <=100) {
			int index = rng.nextInt(nodeList.size());
			String v = nodeList.get(index);
			while (seedlist.contains(v)) {
				index = rng.nextInt(nodeList.size());
				v = nodeList.get(index);
			}
			
			numLoop ++;
			if (g.getGraph().edgesOf(v).size() < 3)
				continue;
			ArrayList<String> neighbours = new ArrayList<String> ();
			for (DefaultEdge e : g.getGraph().edgesOf(v)) {
				String n = g.getGraph().getEdgeSource(e).contentEquals(v)? g.getGraph().getEdgeTarget(e):
					g.getGraph().getEdgeSource(e);
				if (g.getGraph().degreeOf(n) >=3 && !seedlist.contains(n))
					neighbours.add(n);
			}
			if(neighbours.size()>=3) {
				outLoop:
				for (int i=0; i< neighbours.size()-2; i++) {
					String v1= neighbours.get(i);
					for(int j = i+1; j<neighbours.size()-1; j++) {
						String v2= neighbours.get(j);
						if(g.getGraph().containsEdge(v1, v2)) {
							for (int k= j+1; k<neighbours.size(); k++) {
								String v3= neighbours.get(k);
								if(g.getGraph().containsEdge(v1, v3) && g.getGraph().containsEdge(v2, v3)) {
									seedlist.add(v1);
									seedlist.add(v2);
									seedlist.add(v3);
									seedlist.add(v);
									totalClique ++;
									numLoop = 0;
									break outLoop;
								}
							}
						}
					}
				}
			}
			
		}
		
		return seedlist;
	}
	public double sum(ArrayList<Double> list) {
        double sum = 0;
        
        for (double line:list) {
            sum = sum + line;
     
        }
        return sum;
    }
    
    public double average(ArrayList<Double> list) {
      double average = (double)sum(list) / (double) list.size(); 
        return average;
    }

    public double variance(ArrayList<Double> list) {
        double z = 0;
        int y = 0;
        double x = 0;
        for (double word : list) {
        	x = (double) list.get(y)* list.get(y);
        	z = z + x;
            
        	y++;
        }
        
        
        double var = (z - (sum(list) * sum(list)) / list.size()) / (list.size()-1);
//        System.out.println(var);
        return var;
    }
}
class scoreEntity{
	private String vertex;
	private double score;
	
	scoreEntity(String vertex, double score){
		this.score = score;
		this.vertex = vertex;
	}
	void setVertex(String vertex) {
		this.vertex = vertex;
	}
	void setScore(double score) {
		this.score = score;
	}
	double getScore() {
		return this.score;
	}
	String getVertex() {
		return this.vertex;
	}
}