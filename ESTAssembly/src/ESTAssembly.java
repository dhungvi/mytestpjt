import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import com.mhhe.clrs2e.*;

public class ESTAssembly {
	protected final int INT_MAX = Integer.MAX_VALUE;
	protected final int INT_MIN = Integer.MIN_VALUE;
	Properties props;
	Graph g;	//graph to store all the ests. It is generated in 'readEstFile' function.
	SixTuplesGeneration gen;
	Reconstruction rec;
	InclusionNodes incNodes;

	
	public ESTAssembly(Properties p) {
		props = p;
		incNodes = new InclusionNodes();
		g = new Graph(props, incNodes);
		gen = null;
		rec = null;
	}
	
	public void assemble() {
		readEstFile(props.getProperty("EstFile"));
		readMST(props.getProperty("MSTFile"));
		gen = new SixTuplesGeneration(props, g, incNodes);
		rec = new Reconstruction(props, g, gen.getAlignArray(), gen.getLeftEnds(), incNodes);
		rec.getConsensus();
	}
	
	/*
	 * Read ests from the input file;
	 * Generate a Graph object, all the ests are considered to be one node in the graph;
	 * No edge in the graph. Edges will be added in "createAlignArray" function.
	 * 
	 * FASTA format:
	 * A sequence in FASTA format begins with a single-line description, followed by lines of 
	 * sequence data. The description line is distinguished from the sequence data by a greater-than
	 * (">") symbol in the first column. The word following the ">" symbol is the identifier of the 
	 * sequence, and the rest of the line is the description (both are optional). There should be no 
	 * space between the ">" and the first letter of the identifier. It is recommended that all lines 
	 * of text be shorter than 80 characters. The sequence ends if another line starting with a ">" 
	 * appears; this indicates the start of another sequence.
	 * 
	 * EST file comment format:
	 * >g001_001169_001679: first one is number of gene, second is index of tarting position, third is 
	 * 						index of ending position. Index starts from 0.
	 */
	private void readEstFile(String inFileName) {
		boolean bExists = false;
		ArrayList<String> ests;	//store all the ests. 
		ests = new ArrayList<String> ();
		
		try{ 
			File f = (new File(inFileName));
			bExists = f.exists();
			if (!bExists) {
				System.out.println("File does not exist!");
				return;
			}

			BufferedReader in = new BufferedReader(new FileReader(f));
			String str = in.readLine();
			while (str != null) {
				str = str.trim();
				// first line is comment line which begins from '>'
				if (str.charAt(0) == '>') {	//comment line begins from '>'
					String[] paras = str.split("_");
					ests.add(paras[1]);
					//ests.add("0"); //starting position
					ests.add(str); //comment
					
					//get est in the next lines
					str = in.readLine();
					StringBuffer estStr = new StringBuffer();
					while (str != null) {
						str = str.trim();
						if (str.compareTo("") != 0)	{
							if (str.charAt(0) != '>') {
								estStr.append(str.trim());
							} else  {
								ests.add(estStr.toString());
								break;
							}
						}
						str = in.readLine();
					}
					if (str == null) {
						ests.add(estStr.toString());
					}
				} 
			}
			in.close();			
		}catch(IOException e){ 
			System.out.println(e.toString());
			return;
		}
		
		//generate a graph from the input file
		int i=0;
		while (i<ests.size()) {
			//the sequence is upper-case
			g.addNode(new Node(ests.get(i), ests.get(i+1), ests.get(i+2).toUpperCase()));
			i = i+3;
		}
	}
	
	/*
	 * read a minimum spanning tree from the input MST file
	 */
	private void readMST(String inFileName) {
		int nOfNodes = g.graphNodes.size();
		int[][] nodes = new int[nOfNodes-1][3];	//store edges in MST, there are n-1 edges, n is number of nodes.

		//read mst from the input file
		boolean bExists = false;
		try{ 
			File f = (new File(inFileName));
			bExists = f.exists();
			if (!bExists) {
				System.out.println("File does not exist!");
				return;
			}

			BufferedReader in = new BufferedReader(new FileReader(f));
			String str = in.readLine();
			int curIndex = 0;
			while (str != null) {
				str = str.trim();
				if (str.charAt(0) != '#') {	//comment line begins from '#'
					String[] paras = str.split(",");
					if (Integer.parseInt(paras[0]) != -1) {	//-1 means root of MST
						int i0 = Integer.parseInt(paras[0]);
						int i1 = Integer.parseInt(paras[1]);
						int i2 = Integer.parseInt(paras[2]);
						nodes[curIndex][0] = i0;
						nodes[curIndex][1] = i1;
						nodes[curIndex][2] = i2;
						curIndex++;
					}
				} 
				str = in.readLine();
			}
			in.close();			
		}catch(IOException e){ 
			System.out.println(e.toString());
			return;
		}
		
		// Make a undirected MST.
		WeightedAdjacencyListGraph mst = 
			new WeightedAdjacencyListGraph(nOfNodes, false);
		for (int i=0; i<nOfNodes; i++) {	//i is the index of the node in graph
			mst.addVertex(i, Integer.toString(i));
		}
		for (int j=0; j<nodes.length; j++) {
			mst.addEdge(nodes[j][0], nodes[j][1], nodes[j][2]);
		}
		
		g.setMst(mst);
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

	public static void main(String[] args) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
	    	return;
		}

		ESTAssembly assemble = new ESTAssembly(props);

		assemble.assemble();
	}
}
