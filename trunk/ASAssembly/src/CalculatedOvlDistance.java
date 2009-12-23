// store all the calculated overlap distance
import java.util.TreeSet;

public class CalculatedOvlDistance {
	TreeSet<Dist> distances;
	
	public CalculatedOvlDistance() {
		distances = new TreeSet<Dist> ();
	}
	
	public void addDistance(int i1, int i2, int ovlDis, int ovlLen) {
		distances.add(new Dist(i1, i2, ovlDis, ovlLen));
	}
	
	/*
	 * search in distances to find if there is information between i1 and i2.
	 * If does, return int[2] which is same as the return from "getOVLDistance" in OvlDistance.java.
	 * int[0] - overlap length; int[2] - overlap distance.
	 * If not, return null.
	 */
	public int[] searchDistance(int i1, int i2) {
		int[] ret = new int[2];
		
		if (distances.contains(new Dist(i1, i2, 0 ,0))) {
			Dist tmp = (Dist) distances.tailSet(new Dist(i1, i2, 0, 0), true).first();
			ret[0] = tmp.ovlLen;
			ret[1] = tmp.ovlDis;
		} else if (distances.contains(new Dist(i2, i1, 0, 0))) {
			Dist tmp = (Dist) distances.tailSet(new Dist(i2, i1, 0, 0), true).first();
			ret[0] = (-1) * tmp.ovlLen;
			ret[1] = (-1) * tmp.ovlDis;
		} else {
			ret = null;
		}
		return ret;
	}
	
	public static void main(String args[]) {
		CalculatedOvlDistance in = new CalculatedOvlDistance();
		in.addDistance(1048575, 1048574, 11, 20);
		in.addDistance(1048575, 1048578, -12, -20);
		in.addDistance(1048575, 1048572, 13, 20);
		in.addDistance(1048575, 5, -14, -20);
		in.addDistance(1, 6, 15, 20);
		
		int[] ret = in.searchDistance(1048575, 1048578);
		if (ret != null) {
			System.out.println(ret[0] + " " + ret[1]);
		} else {
			System.out.println("Not found!");
		}
		
	}



	class Dist implements Comparable<Dist>{
		long indexOfNodes;
		int ovlDis;
		int ovlLen;
		
		Dist(int idx1, int idx2, int dis, int len) {
			indexOfNodes = (idx1 << 20) + idx2;
			ovlDis = dis;
			ovlLen = len;
		}
		
		public int compareTo(Dist other) {
			//Returns 0 if the argument is equal to this; 			
			//a value less than 0 if the argument is greater than this; 
			//and a value greater than 0 if the argument is less than this. 
			if (this.indexOfNodes == other.indexOfNodes) {
				return 0;
			} else if (this.indexOfNodes > other.indexOfNodes) {
				return 1;
			} else {
				return -1;
			}
		}
	}

}

