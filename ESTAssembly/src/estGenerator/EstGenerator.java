package estGenerator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class EstGenerator {
	private final String SourceFile = "sq.fa";	//only has one line, no char '>' at start
	private final String OutFile = "est.fa";
	private int uniLower = 0;	//lower bound of uniform distribution	
	private int uniUpper = 280;  //upper bound of uniform distribution,uniUpper=lenOfGene-expoMean
	private int expoMean = 20; 	//mean of exponential distribution
	private int expoLower = 17; //lower bound of exponential distribution
	private int expoUpper = 23; //upper bound of exponential distribution
	private int numEsts = 60;	//number of ests
	private String oriStr;	//store the original string from SourceFile
	private RandomNum ran;
	
	public EstGenerator() {
		ran = new RandomNum();
	}
	
	public EstGenerator(long seed) {
		ran = new RandomNum(seed);
	}

	private void readSourceFile() {
		boolean bExists = false;
		
		try{ 
			File f = (new File(SourceFile));
			bExists = f.exists();
			if (!bExists) {
				System.out.println("SourceFile does not exist!");
				return;
			}

			BufferedReader in = new BufferedReader(new FileReader(f));
			oriStr = in.readLine();
			in.close();	
		}catch(IOException e){ 
			System.out.println(e.toString());
			return;
		}
	}

	private void writeToFile() {
		try {
			File f = (new File(OutFile));
			if (f.exists()) {
				f.delete();
			}
			
			BufferedWriter out = new BufferedWriter(new FileWriter(f));
	        
			for (int i=0; i<numEsts; i++) {
				/*
				 * format is: 
				 * >startPos.lenEst
				 * est
				 * >startPos.lenEst
				 * est
				 */
				out.write(">");
				int startPos = ran.unifRan(uniLower, uniUpper);
				int lenEst = ran.expoRan(expoMean, expoLower, expoUpper);
				out.write(Integer.toString(startPos));
				out.write(".");
		        out.write(Integer.toString(lenEst));
		        out.write("\n");
		        int endIndex = startPos + lenEst;
		        String inStr = "";
		        if (endIndex > oriStr.length()) {
		        	inStr = oriStr.substring(startPos);
		        } else {
		        	inStr = oriStr.substring(startPos, endIndex);
		        }
		        //out.write(ran.errEst(inStr));
		        out.write(inStr);
		        out.write("\n");
			}
	        out.flush();	
	        out.close();
	    } catch (IOException e) {
	    	System.out.println(e.toString());
	    }
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		EstGenerator gen = new EstGenerator();
		gen.readSourceFile();
		gen.writeToFile();
	}

}
