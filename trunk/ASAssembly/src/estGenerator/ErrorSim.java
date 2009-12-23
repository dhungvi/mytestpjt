package estGenerator;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

/*
 * This class simulate two kinds of errors: single base error and stutter.
 */
public class ErrorSim {
	RandomNum rand;
	double paraP1;
	double paraP2;
	double paraP3;
	double para1P4;
	double para2P4;
	double[] singleErrorProb; //substitution, deletion or insertion, each with probability 1/3
														  //here [0]:substitution, [1]:deletion, [2]:insertion.
	String oriEst;	//original est
	int[][] errorFlag;	//record error type. The first dimension is the index of the base; the second dimension is the flag of error.
						//errorFlag[][0]: 0-no single base error; 1-substitution; 2-deletion; 3-insertion.
						//errorFlag[][1]: 0-no stuttering error; n-number of added bases.
	int slashEta = 20;  //a parameter to calculate the probability of stuttering error.
	
	/*
	 * constructor
	 */
	public ErrorSim(RandomNum r, Properties props) {
		rand = r;
		paraP1 = Double.parseDouble(props.getProperty("paraP1"));
		paraP2 = Double.parseDouble(props.getProperty("paraP2"));
		paraP3 = Double.parseDouble(props.getProperty("paraP3"));
		para1P4 = Double.parseDouble(props.getProperty("para1P4"));
		para2P4 = Double.parseDouble(props.getProperty("para2P4"));
		double d1 = Double.parseDouble(props.getProperty("singleErrorProbSub"));
		double d2 = Double.parseDouble(props.getProperty("singleErrorProbIns"));
		double d3 = Double.parseDouble(props.getProperty("singleErrorProbDel"));
		singleErrorProb = new double[]{d1,d2,d3};
	}
	
	/*
	 * generating errors to an EST and return the modified est
	 * @param s the original est
	 * @return String an est with simulated errors bases on the original one
	 */
	public String genErrors(String s) {
		oriEst = s;
		errorFlag = new int[s.length()][2];
		
		for (int i=0; i<s.length(); i++) {
			errorFlag[i][0] = singleBaseError(i+1);
			errorFlag[i][1] = stutterError(i);
		}
		
		/*
		 * generate errors
		 */
		StringBuffer errorEst = new StringBuffer();
		for (int i=0; i<oriEst.length(); i++) {
			//0-no error; 1-substitution; 2-deletion; 3-insertion
			switch(errorFlag[i][0]) {
			case 0:
				errorEst.append(oriEst.charAt(i));
				break;
			case 1:
				errorEst.append(rand.genDiffBase(oriEst.charAt(i)));
				break;
			case 2:
				break;
			case 3:
				errorEst.append(rand.genBase());
				errorEst.append(oriEst.charAt(i));
				break;
			}
			
			//0-no error; n-number of repetition. We have put the current character into errorEst, so we just need add repeated
			//characters here.
			if (errorFlag[i][1] != 0) {
				for (int j=0; j<errorFlag[i][1]; j++) {
					errorEst.append(oriEst.charAt(i));
				}
			}
		}

		return errorEst.toString();
	}
	
	/*
	 *This method simulates single-base error.
	 *	1 Allow it to change at a probability of 0.01.
	 *	2 If it did not change in (1), then give a second change to change based on its distance from the 3’ end.  This probability here 
	 *		is 0.02x/l, where x is the index number of the base and l is the length of the est.  (So the first base has a probability of 
	 *		0.02/l of changing for this reason, and the last base has a probability of 0.02.)
	 *  3 At the end of the est things get even worse. So if there was no change in step (1) or (2), let it change with a probability of 
	 *  	exp(-0.5sqrt(l-x)). (So the first base has a probability of exp(-0.5sqrt(l)) -- very small, while the second-to-last base has
	 *  	a probability 0.606 and the last base has a probability of 1.
	 *  4 There is also extra problems at the beginning. So an error can happen with probability 0.5(1 – tanh(x – 50)). 
	 *  If a base is subjected to an error: it could be a substitution, deletion or insertion, each with probability 1/3.  (Here we treat 
	 *  “insertion” as inserting a new base right before it.)  For substitution and insertion, we pick from the other letters with a uniform
	 *   probability. 
	 *   
	 *  In this method, we decide if there is an error for the base with the index by calculating "1 - (1 – P1)(1 – P2)(1 – P3) (1-p4)".
	 *
	 *  @param index index of the base, starting from 1.
	 *  @return error type. 0-no error; 1-substitution; 2-deletion; 3-insertion.
	 */
	private int singleBaseError(int index) {
		int len = oriEst.length(); //length of the est
		double p1 = paraP1;
		double p2 = paraP2*index/len;
		double p3 = Math.exp(paraP3*Math.sqrt(1-index));
		double p4 = para1P4*(1-Math.tanh(index-para2P4));
		double p = 1 - (1-p1)*(1-p2)*(1-p3)*(1-p4); //the probability of an error happening
		
		if (rand.happenOrNot(p)) { //an error happens
			int i = rand.whichHappen(singleErrorProb); //0:substitution, 1:deletion, 2:insertion.
			return i+1;
		} else {
			return 0;
		}
	}
	
	/*
	 * This method simulates stuttering error.
	 * Stuttering: This is when a subsequence gets copied.  This is apparently more likely to happen at a T or a G.
	 * This error is only allowed when there is a sequence of multiple Ts or multiple Gs. If we are at position x, let r_G(x) 
	 * be the number of contiguous Gs that occur up to and including x. (So for the sequence AGGGA, r_G(0)=0, r_G(1) = 1, 
	 * r_G(2) = 2, r_G(3) = 3, and r_G(4) = 0.)  Let r_T(x) be the same for Ts. Then the probability of a stutter (in terms of
	 * the parameter \eta) is:   (0.5(1-cos(r_G(x)/ \eta)))^2 + (0.5(1-cos(r_T(x)/ \eta)))^2
	 * If stuttering does occur at index x, then the length of the insert sequence is chosen from a uniform distribution over the 
	 * range (0...2*r_G(x)) or (0...2*r_T(x)).
	 * 
	 *  In this method, we decide if there is an stuttering error for the base with the index, if no, we will return 0. If there 
	 *  is an error, we will return the number of repetition.
	 *
	 *  @param index index of the base, starting from 0; len length of the est.
	 *  @return error type. 0-no error; n-number of repetition.
	 */
	private int stutterError(int index) {
		int numOfKey = 0; //the number of contiguous 'T' or 'G'
		if (oriEst.charAt(index) == 'T') {
			for (int i=index; i>=0; i--) {
				if (oriEst.charAt(i) == 'T') {
					numOfKey++;
				} else {
					break;
				}
			}
		} else if (oriEst.charAt(index) == 'G') {
			for (int i=index; i>=0; i--) {
				if (oriEst.charAt(i) == 'G') {
					numOfKey++;
				} else {
					break;
				}
			}
		}
		
		int retVal = 0;
		if (numOfKey != 0) {
			double p = Math.pow((0.5*(1-Math.cos(numOfKey/slashEta))), 2);
			if (rand.happenOrNot(p)) {
				retVal = rand.unifRan(0, 2*numOfKey+1);
			}
		}
		
		return retVal;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
	    	return;
		}

		RandomNum ran = new RandomNum();
		ErrorSim es = new ErrorSim(ran, props);
		String input = "ATCGATCTTTTTGGGACTG";
		System.out.println(es.genErrors(input));
	}
	
	//only used for test by main
	private static Properties getProperties(String fName) throws IOException {
		Properties props = new Properties();
		File f = new File(fName);
        
        if (!f.exists()) {
        	return props;
        }
        
        props.load(new FileInputStream(f)); 
        return props;
    }

}
