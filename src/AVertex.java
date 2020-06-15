


import java.util.Random;
import java.util.Vector;

public class AVertex {
	final String vid;
	final int num_attribute = 20;
	final Vector<Integer> X;

	AVertex(String vid, Vector<Integer> X) {
		this.vid = vid;
		this.X = X;
	}

	AVertex(String vid) {
		this.vid = vid;
		X = new Vector<Integer>(num_attribute);
		for (int i = 0; i < num_attribute; i++)
			X.add(0);
	}

	public void generate_random_attribute() {
		Random rng = new Random();

		for (int i = 0; i < X.size(); i++) {
			double r = rng.nextDouble();
			if (r > 0.5)
				X.set(i, 1);
			else
				X.set(i, 0);

		}

	}

	public String getID() {
		return this.vid;
	}

	public String getAttributeInString() {
		String str = "[";
		for (int i = 0; i < this.X.size(); i++) {
			str = str + X.get(i);
			if (i != X.size() - 1)
				str = str + ", ";
		}
		str = str + "]";
		return str;
	}

	public Vector<Integer> getAttributes() {
		return X;
	}

	public String toString() {
		return vid;
	}

	public int hashCode() {
		return this.toString().hashCode();
	}

	public boolean equals(AVertex v) {
		return (this.vid.equals(v.getID()));
	}
}
