package estGenerator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;

public class EstGenerator {
	private String SourceFile;	//only has one line, no char '>' at start
	private String OutFile;
	private int uniLower;	//lower bound of uniform distribution	
	private int uniUpper;  //upper bound of uniform distribution,uniUpper=lenOfGene-expoMean
	private int expoMean; 	//mean of exponential distribution
	private int expoLower; //lower bound of exponential distribution
	private int expoUpper; //upper bound of exponential distribution
	private int numEsts;	//number of ests
	
	private String oriStr;	//store the original string from SourceFile
	private RandomNum ran;
	private ErrorSim errSim;
	
	public EstGenerator(Properties props) {
		ran = new RandomNum();
		errSim = new ErrorSim(ran, props);
		SourceFile = props.getProperty("SourceFile");
		OutFile = props.getProperty("OutFile");
		uniLower = Integer.parseInt(props.getProperty("uniLower"));
		uniUpper = Integer.parseInt(props.getProperty("uniUpper"));
		expoMean = Integer.parseInt(props.getProperty("expoMean"));
		expoLower = Integer.parseInt(props.getProperty("expoLower"));
		expoUpper = Integer.parseInt(props.getProperty("expoUpper"));
		numEsts = Integer.parseInt(props.getProperty("numEsts"));
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
				//int lenEst = ran.expoRan(expoMean, expoLower, expoUpper);
				int lenEst = ran.unifRan(expoLower, expoUpper);	//use uniform distribution temporarily for real gene
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
		        out.write(errSim.genErrors(inStr));
		        //out.write(ran.errEst(inStr));
		        //out.write(inStr);
		        out.write("\n");
			}
	        out.flush();	
	        out.close();
	    } catch (IOException e) {
	    	System.out.println(e.toString());
	    }
	}

	protected static Properties getProperties(String fName) throws IOException {
		Properties props = new Properties();
		File f = new File(fName);
        
        if (!f.exists()) {
        	return props;
        }
        
        props.load(new FileInputStream(f)); 
        return props;
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

		EstGenerator gen = new EstGenerator(props);
		gen.readSourceFile();
		gen.writeToFile();
	}

}
