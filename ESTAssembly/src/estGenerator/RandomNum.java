package estGenerator;

import java.util.Random;

public class RandomNum {
	private Random ran;

	// Initialize with random seed
    public RandomNum() {
    	ran = new Random();
    }
 
    // Initialize with user specified seed
    public RandomNum(long s) {
    	ran = new Random(s);
    }

    // Returns a uniformly distributed random value from the interval [lower,upper)
	public int unifRan(int lower, int upper) {
		int reVal = ran.nextInt(upper-lower) + lower;
		return reVal;
	}
	
	//Returns a random value from an Exponential distribution with the given mean
	public int expoRan(int mean, int lower, int upper) {
		int reVal;
		do {
			reVal= (int)(-mean * Math.log(unif01()));
		} while ((reVal<lower) || (reVal>upper));
		return reVal;
	}

	// Generate uniform 0-1 random number
    private double unif01() {
        double r;
        do {
            r = ran.nextDouble();
        } while (r == 0);
        return r;
    }

    // generate error in the input string. 0.4base/100bases error.
	public String errEst(String s) {
		int len = s.length();
		String retStr = s;
		String tmp = null;
		double prob = len / 100 * 0.4;
		prob = 1;
		for (int i=0; i<5; i++){
		if (ran.nextDouble() < prob) {	//generate errors in the string
			int pos = unifRan(0, len);
			
			if (ran.nextDouble() < 1/3) { //delete
				if (pos < len-1) { //not the last character
					retStr = s.substring(0, pos) + s.substring(pos+1);
				} else {
					retStr = s.substring(0, pos);
				}
			} else if (ran.nextDouble() < 2/3) { //insert
				tmp = genBase();
				retStr = s.substring(0, pos) + tmp + s.substring(pos);
			} else { //change
				do {
					tmp = genBase();
				} while (s.charAt(pos) == tmp.charAt(0));
				if (pos < len-1) { //not the last character
					retStr = s.substring(0, pos) + tmp + s.substring(pos+1);
				} else {
					retStr = s.substring(0, pos)+ tmp;
				}
			}
			
		}
		}
		return retStr;
	}

	public String genBase() {
		if (ran.nextDouble() < 0.25) {
			return "A";
		} else if (ran.nextDouble() < 0.5) {
			return "G";
		} else if (ran.nextDouble() < 0.75) {
			return "C";
		} else {
			return "T";
		} 
	}

	/*
	 * generate a different base than the input one with a uniform probability.
	 */
	public char genDiffBase(char c) {
		char tmp;
		do {
			tmp = genBase().charAt(0);
		} while (c == tmp);
		return tmp;
	}
	/*
	 * decide if one event will happen according to the probability
	 * @param p the probability of the event happening
	 * @return true or false
	 */
	public boolean happenOrNot(double p) {
		if (ran.nextDouble() < p) {
			return true;
		} else {
			return false;
		}
	}
	
	/*
	 * Given a list of events, decide which one will happen
	 * @param p an double array which store the probability of each event
	 * @return int the index of the event in the input array
	 */
	public int whichHappen(double[] prob) {
		double[] p = new double[prob.length];
		for (int i=0; i<prob.length; i++) {
			p[i] = prob[i];
		}
		
		double tmp = ran.nextDouble();
		for (int i=1; i<p.length; i++) {
			p[i] = p[i-1] + p[i]; 
		}
		int index = 0;
		for (; index<p.length-1; index++) {
			if (tmp < p[index]) {
				return index;
			}
		}
		return p.length-1;
	}
}
