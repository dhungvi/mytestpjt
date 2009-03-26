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

}
