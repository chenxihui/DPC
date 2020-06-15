import java.util.ArrayList;

public class MD_node {
	private ArrayList<String> vertexList;
	private int level;
	private ArrayList<MD_node> children;
	private MD_node parent;
	
	/*variable self is used in the bestcut to indicate whether the node 
	is selected as a block in the final partition. set to false by default
	modularity is also used in best cut to store the modularity of the block 
	set to 0 by default*/
	
	private boolean self;
	private double modularity;
	
	public MD_node(int level) {
		this.vertexList = new ArrayList<String>();
		this.children = new ArrayList<MD_node>();
		this.setLevel(level);
		this.parent = null;
		this.setSelf(false);
		this.modularity = 0;
	}
	public MD_node(ArrayList<String> vertexList, int level) {
		this.vertexList = vertexList;
		this.children = new ArrayList<MD_node>();
		this.setLevel(level);
		this.parent = null;
		this.setSelf(false);
		this.modularity = 0;
	}
	public String toString() {
		String str = "";
		str = str.concat("level:"+this.level+";").concat("vertextlist:"
		    +this.vertexList.toString()+";").concat("modularity:"
		    +this.modularity+";").concat("self:" 
		    + this.self+";").concat("num_of_children:" 
		    + this.children.size()+"\n");
		return str;
	}
	public ArrayList<String> getVertexList() {
		return vertexList;
	}

	public void setVertexList(ArrayList<String> vertexList) {
		this.vertexList = vertexList;
	}

	public ArrayList<MD_node> getChildren() {
		return children;
	}

	public void setChildren(ArrayList<MD_node> children) {
		this.children = children;
	}
	public int getLevel() {
		return level;
	}
	public void setLevel(int level) {
		this.level = level;
	}
	public MD_node getParent() {
		return parent;
	}
	public void setParent(MD_node parent) {
		this.parent = parent;
	}
	public boolean getSelf() {
		return self;
	}
	public void setSelf(boolean self) {
		this.self = self;
	}
	public double getModularity() {
		return modularity;
	}
	public void setModularity(double modularity) {
		this.modularity = modularity;
	}
	
}
