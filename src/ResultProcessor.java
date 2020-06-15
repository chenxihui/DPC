import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.jgrapht.graph.DefaultEdge;

public class ResultProcessor {
	
	String file ; 
	
	ResultProcessor(){
		
	}
	
	double compareThetaF_ratio(String generatedGraph, String originalGraph, 
			String abtGenerated, String abtOriginal, String partitionFileOri, int num_abt) throws IOException {
		double distance = 0.0;
		this.file = generatedGraph;
		graph g_gen = new graph();
		g_gen.readGraphOnly(generatedGraph);
		
		graph g_ori = new graph();
		g_ori.readGraphOnly(originalGraph);
		
		HashMap<String, ArrayList<Integer>> attribute_gen = readAttribute(abtGenerated);
		HashMap<String, ArrayList<Integer>> attribute_ori = readAttribute(abtOriginal);
	
		System.out.println("Done. Start reading the partition ... " + generatedGraph + " " + originalGraph);
		Partition pOri;
		BufferedReader in = new BufferedReader(new FileReader(partitionFileOri));

		String pString = "[";
		String line;

		while((line = in.readLine())!=null) {
			pString = pString + "[";
			String [] fds = line.split("\t");
			for(int i = 0; i<fds.length; i++)
				pString = pString + fds[i] +",";
			pString = pString +"], ";
		}
		pString = pString.substring(0, pString.length()-2)+"]";
		
//		String pString = in.readLine();
		pOri = new Partition(pString);
		in.close();

		Theta_F_Ratio_Store ratio_gen = this.calculateThetaF_ratio(g_gen, attribute_gen,attribute_ori, pOri, num_abt);
//		System.out.println("original");
		Theta_F_Ratio_Store ratio_ori = this.calculateThetaF_ratio(g_ori, attribute_ori,attribute_gen, pOri, num_abt);
	
	    //distance = compareRatioKL(ratio_gen, ratio_ori);	
	    distance = compareRatioHellinger(ratio_gen, ratio_ori);
		
		return distance;
	}

	double compareRatio(Theta_F_Ratio_Store ratio_gen, Theta_F_Ratio_Store ratio_ori) {
		
		double [][] theta_f_gen = ratio_gen.theta_f_ratio;
		double[][] theta_f_ori = ratio_ori.theta_f_ratio;
		
		double[] theta_btw_gen = ratio_gen.theta_f_btw_ratio;
		double[] theta_btw_ori = ratio_ori.theta_f_btw_ratio;
		
		double avg = 0;
		for (int com_index = 0; com_index < theta_f_gen.length; com_index ++) {
			double total_sum = 0;
			for(int i = 0; i< theta_f_gen[com_index].length; i++) 
				total_sum += (Math.abs(theta_f_gen[com_index][i] - theta_f_ori[com_index][i]));
			avg += (total_sum/theta_f_gen[0].length);
		}

		double total_sum = 0;
		for(int i = 0; i< theta_btw_gen.length; i++) 
			total_sum += (Math.abs(theta_btw_gen[i] - theta_btw_ori[i]));
		avg += (total_sum/theta_btw_gen.length);
		
		avg /= (theta_f_gen.length +1);
			
		return avg;
	}
	double compareRatioKL(Theta_F_Ratio_Store ratio_gen, Theta_F_Ratio_Store ratio_ori) {
		
		double [][] theta_f_gen = ratio_gen.theta_f_ratio;
		double[][] theta_f_ori = ratio_ori.theta_f_ratio;
		
		double[] theta_btw_gen = ratio_gen.theta_f_btw_ratio;
		double[] theta_btw_ori = ratio_ori.theta_f_btw_ratio;
		double avg = 0;
		for (int com_index = 0; com_index < theta_f_gen.length; com_index ++) {
			avg += (KL(theta_f_ori[com_index], theta_f_ori[com_index]));
		}
		avg += (KL(theta_btw_gen, theta_btw_ori));
		
		avg /= (theta_f_gen.length +1);
			
		return avg;

	
	}
	
	double compareRatioHellinger(Theta_F_Ratio_Store ratio_gen, Theta_F_Ratio_Store ratio_ori) {
		
		double [][] theta_f_gen = ratio_gen.theta_f_ratio;
		double[][] theta_f_ori = ratio_ori.theta_f_ratio;
		
		double[] theta_btw_gen = ratio_gen.theta_f_btw_ratio;
		double[] theta_btw_ori = ratio_ori.theta_f_btw_ratio;
		double avg = 0;

		for (int com_index = 0; com_index < theta_f_gen.length; com_index ++) {
			double d = Hellinger(theta_f_gen[com_index], theta_f_ori[com_index]);
//			System.out.println(d);
			avg += d;

		}
		double d = (Hellinger(theta_btw_gen, theta_btw_ori));
		if (d>1) d=1;
//			System.out.println(d);
		avg+=d;
		avg /= (theta_f_ori.length +1);
			
		return avg;
	}
	
	double Hellinger(double[] p, double [] q) {
		double sum = 0;
		for (int i =0; i<p.length; i++) {
			//System.out.println(p[i] + " "+ q[i]);
			sum += Math.pow(Math.sqrt(p[i])-Math.sqrt(q[i]), 2);
		}
		
		sum = Math.sqrt(sum);
		return sum/Math.sqrt(2);
		
	}
	
	double KL(double[] p, double [] q) {
		double sum = 0;
		for (int i =0; i<p.length; i++)
			sum += (p[i] * Math.log(p[i]/q[i])/Math.log(2.0));
		for (int i =0; i<p.length; i++)
			sum += (q[i] * Math.log(q[i]/p[i])/Math.log(2.0));
		return sum;
		
	}
	public HashMap<String, ArrayList<Integer>> readAttribute(String file) throws IOException{

		HashMap<String, ArrayList<Integer>> attribute = new HashMap<String, ArrayList<Integer>>();
		
		BufferedReader in = new BufferedReader(new FileReader(file));
		
		String line = "";

		while ((line = in.readLine()) != null) {
				
			String [] fds = line.split(" ");
			ArrayList<Integer> av = new ArrayList<Integer>(fds.length-1);
			
			for(int i =1; i<fds.length; i++) 
				av.add((int)Double.parseDouble(fds[i]));
		
			attribute.put(fds[0], av);

			//this.num_abt = av.size();
		}
		in.close();

		return attribute;
	}	
	
	
	public Theta_F_Ratio_Store calculateThetaF_ratio(graph g1, HashMap<String, ArrayList<Integer>> abt_rnd, 
			HashMap<String, ArrayList<Integer>> abt_rnd_2,
			Partition p, int num_abt ){
    	
		double[][] theta_f_g1_ratio = new double[p.partition.size()][11]; 
		
		double[] theta_f_btw_g1_ratio = new double[11];
		
		Iterator<ArrayList<String>> iter = p.partition.iterator();
		
		int com_index = 0;

		while(iter.hasNext()) {

			ArrayList<String> com = iter.next();

			for(String node : com) {
				if(!g1.getGraph().containsVertex(node)) {
//					System.out.println("xxx");
					continue;
				}
				for (DefaultEdge e : g1.getGraph().edgesOf(node)) {
					String vs = g1.getGraph().getEdgeSource(e);

					if(vs.contentEquals(node))
						vs = g1.getGraph().getEdgeTarget(e);
					

					ArrayList<Integer> abt_node = abt_rnd.get(node);
					ArrayList<Integer> abt_vs = abt_rnd.get(vs);
					if (abt_rnd.get(node) == null) {
						//System.out.println("no attribute found");
						abt_node = abt_rnd_2.get(node);
					//	continue;
					}
					if (abt_vs == null) {
						abt_vs= abt_rnd_2.get(vs);
					//	continue;
					}
					if (abt_rnd.get(node) == null || abt_vs == null) {
						continue;
					}
					int ratio = (int)Math.floor(abt_distance(abt_node, abt_vs, num_abt)/0.1);

					if (!com.contains(vs)) {
						theta_f_btw_g1_ratio[ratio] += 1;
							
					} else {
						theta_f_g1_ratio[com_index][ratio] += 1;
					}
					
				}
			}
			
			double sum = 0.0 ;
			for (int i=0; i< 11; i++) 
				sum += theta_f_g1_ratio[com_index][i];
				
			for (int i = 0; i<11; i++) {
					theta_f_g1_ratio[com_index][i] = (theta_f_g1_ratio[com_index][i]*1.0)/(sum + 1);
			}
				
			com_index ++;
		}
		
		double sum = 0.0;
		for(int i = 0; i < 11; i++)  
			sum += theta_f_btw_g1_ratio[i];
			
		for (int i = 0; i< 11; i++) {
				theta_f_btw_g1_ratio[i] = (theta_f_btw_g1_ratio[i]*1.0)/ (sum +1);
		}
		
		Theta_F_Ratio_Store res = new Theta_F_Ratio_Store(theta_f_g1_ratio, theta_f_btw_g1_ratio);
		return res; 
	}
	
    double abt_distance(ArrayList<Integer> abt_vs, ArrayList<Integer> abt_vt, int num_abt) {
    	double distance = 0 ;

    
    	double prod = 0;
    	double ds = 0 ;
    	double dt = 0; 
    	
    	for (int i = 0; i<num_abt; i++ ) {
    	   prod += (abt_vs.get(i) * abt_vt.get(i));
    	   ds += Math.pow(abt_vs.get(i), 2);
    	   dt += Math.pow(abt_vt.get(i), 2);
    	}
    	distance = prod *1.0 /(Math.sqrt(ds) * Math.sqrt(dt));

    	return distance;
    }
}

class Theta_F_Ratio_Store{
	double [][] theta_f_ratio;
	double [] theta_f_btw_ratio;
	Theta_F_Ratio_Store(double [][] theta_f_ratio, double[] theta_f_btw_ratio){
		this.theta_f_btw_ratio = theta_f_btw_ratio;
		this.theta_f_ratio = theta_f_ratio;
	}
}