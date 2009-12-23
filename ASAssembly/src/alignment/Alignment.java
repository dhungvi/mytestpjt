package alignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

// Test all seven alignment algorithms

public class Alignment {
  public static void main(String[] args) {
		String s1 = "GCAacGCTgTGCAGGAATCTAgTTTTTACGCTTTGAGGTGCAAACCCCCCCACAGCCCTTCAtCGCCGGAGTTTTTTTTTTTATGCCGCTATgGGGAGCGGAGACGCCCAGCANCCTGGGGCCGTCATATCTTGTGgNCCCATGGCCCGCTATGGAGCGGACTAaCTTNCAGATTACcAGCAGCAaGCCTGGTGCtGCgGCGGACGGCAGCaTCGGCGGCTgGGGcgATGGtGGCGaGCAAGGTCGCTgACAGTCGCTGCCAGAACCGCCAGTCATTAAGCAcCAGAAAcTATCCGCTCGTGGGaAACGcAGTGCCCCAGCAGCAGTTATcNGTGcAGCTCtCGcCACGCAGGGAGCGCCGGAtaCCAGCCCCaGCAggTACAGCTATAAATTTTCAGCAACAGAtTCTCCcTTTtGGGCTTTCTtTNaTTTGTTGCTTACCGtCTAGTCGCCGcCTGGAAGcCTGCGCTAATCGGCCTGACACCtTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGAGGcN";
		String s2 = "CAcCCCTaCAgCGCCGGcGtNGCCGgTATGGGNNNNNNNGAGCGGAGcCGCCGGGGGGGGCAGCAGCCTTTGGGGCCGacNNATATCTTGTGCCCCCATGgtCCCGCTATGGAGCgGGACACCTTACAGATTACAGCAGCATGtTCTGGTGCTGCCGCGtACGGCAGCCTCGGCGGCTGGGGTCNTGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAtATTATCCGCTCGTGcGAACGCAGTGCCCCAcCaAGCaGTTATCCGTGGAGCTCGCGCACGNAGGGtAGCGCCGGAGGCCAGCCCCAGCaCCTACAGCTATAAATTTTCAGCAAaAaATTCTgCATTTTGGGCTTCTTTGGTTTGtTGCTTACCgcCTAGtCGCCgTCTGGAAGGCTGCGCTAATCGGtCTGACACCCTTGGgGCGTTGCCTCTCAtTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCcCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCaGAAGTCTCaCTTCTGCCaTg";
	s1 = s1.toUpperCase();
	s2 = s2.toUpperCase();

    Substitution sub = new ForGene();
    SWSmart smart = (new SWSmart (sub, 2, s1, s2));
    String[] strs = smart.getMatch();
    System.out.println(smart.getScore());
 	try{ 
		File outFile = new File("analysis.out");
		boolean bExists = outFile.exists();
		if (bExists) {
			outFile.delete();
		}
		BufferedWriter out1 = new BufferedWriter(new FileWriter(outFile, true));
		for (int i=0; i<strs.length; i++) {
			out1.write(strs[i]);
			out1.write("\n");
		}
		out1.flush();
		out1.close();
	}catch(IOException e){ 
		System.out.println(e.toString());
	} 
	
	String consensus = "CCCGATAGTCTTTGACCGTTCTATTTTTAGGGACCCGTNCCTAGTGGACGTTAAACTCTNGAGGCAACCTTTGACCGTTCTATTGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGATCGTCTTGCGCCACCCTCTCGCTCGCTACGCATGCGCTCCCGGCCGCCTACGCCGTGGCCAACACGGCCAGTTGCAGAGGCCGACATAACAGATCTGCAAACATCGGAATGAGAGACCCGGGAGACGTTCTGCCATCGCCGACACCACCATGATCAGTCGCTCTTTGCTCTGGCCGATCTGCCTAGACCATGGAGCTATGTTCGGACCTATCTGACTGTGATCGATATTATCCTGAGCAGATCGCGACTACCTAGCGAGCGGAGTAGACGGCCGCTGGCGCTCAAATCATTCTAACGGAAGAATGCCCACAGAGCGTTGGGGAGTTTAAGCCAGCTTGAGGACGGAGATCGCCATGATGCGGACGTTGCCGAGTGCCGGAGCGAAACGCCGCGGTGGAGAGCATCTGCAAAATGCGAGAACTCGCGCAGTCGAAGGTGTTGACGAGCTCTGTAGGTTCATACGTGGCAATCGTTGAAGTTGCATTTGTCGTGCTCATGTCTTTCAGTGTTTGAGACGCACATGTTGCATCGATGTTGGAGGCGCTGTATGCCCCAGCTCGCCCAGGAGCTCATAGGAGTGGTGCCCTGAATGCCCTGAGGCTCATTTGGTTGTGTTATGTCTGAGGTTCGCCATTAGCCGCCGTGAACGGAGCGACAGGCTCAACGCGCATGGCTATCACGCCCGCGAGCCCTCTCCATTGCGGTTCACACCAAGGCAGTGTAGGCATGGTTTGTGAGACTTTCAGACGCATGACTGCACTGAGACAGACGTACACGCACGCAGACTGGCAATATGCATACGGCACACCCTGGCACAGGCCCAGACGTTGCAGACGAGCGCCGCATTGACATTGCGCCGCACGCTGGGGGGGGAGTCGCCGTCAGGCACACGAAACGTTCATTGTCCTATCATAGGGCGCGACAGCTGTCGTGACCTGGCGTGTGACCGTGCATATCGTGTGACTGAGTGATTTGGGTGTGGTTTTGGCGTGTGGGAAGAGGGGAACGTGCCCGAGGTTGGTTGGCTCGACACGGTGTATTGCATGGGAGGCCGGAGTACGGGTGATCTCGAGAGAGGATATAGGTCGCGTCTGAGAGCGGCTAGAGCCGTCTCTCATTCCGATCCTAGGCGAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGAAAAACAAAT";
	String gene = "GAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGATCGTCTTGCGCCACCCTCTCGCTCGCTACGCATGCGCACCCGGCCGCCTACACCGTGGCCAACACGGCCAGTTGCAGAGGCCGACATAACAGATCTGCAAACATCGGAATGAGAGACCCGGGAGACGTTCCGCCATCGCCGACCACCATGATCAGTCGCTCTTTGCTCTGGCCGATCTGCCTAGACCATGGAGCTATGTTCGGACCTATCTGACTGTGATCGATATTATCCTGAGCAGATCGCGACTACTAGCGAGCGGAGTAGACGGCCGCTGGCGCTCTTATCATTCTAACGGAAGAATGCCCACAGAGCGTTGGGGAGTTTAAGCCAGCTTGAGGACGGAGATCGCCATGATGCGGACGTTGCCGAGTGCCGGAGCGAAACGCCGCGGTGGAGAGCATCTGCAAAATGCGAGAACTCGCACAGTCGAAGGTGTTGACGAGTTCTGTAGGTTCATACATGGCAGTCGTTGAAGTTGCGTTTGTCGTGCTCATGTCTTTCAGTGTTTGAGACGCACATGTTGCATCGATGTTGGAGGCGCTGTATGCCCCAGCTCGCCCAGGAGCTCATAGGAGTGGTGCCCTGAATGCCCTGAGGCTCATTTGGTTGTGTTATGTCTGAGGTTCGCCATTAGCCGCCGTGAACGGAGCGACAGGCTCAACGCGCATGGCTATCACGCCCGCGAGCCCTCTCCATTGCGGTTCACACCAAGGCAGTGTAGGCATGGTTTGTGAGACTTTCAGACGCATGACTGCACTGAGACAGACGTACACGCACGCAGACTGGCAATATGCATACGGCACACCCTGCACAGGCCCAGACGTTGCAGACGAGCGCCGCATTGACATTGCGCCGCACGCTGGGGGGGGAGTCGCCGTCAGGCACACGAAACGTTCATTGTCCTATCATAGGGCGCGACAGCTGTCATGACCTGGCGTGTGACCGTGCATATCGTGTGACTGAGTGATTTGGATGTGGTTTTGGCGTGTGGGAAGAGGGGAACGTGCCCGAGGTTGGTTGGCTCGACACGGTGTATTGCATGGGAGGCCGGAGTACGGGTGATCCGAGAGAGGATATAGGTCGCGTCTGAGAGCGGCTAGAGCCGTCTCTCATTCCGATCCTAGGCGAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTG";
    Substitution sub2 = new ForAScore();
    AScore ascore = (new AScore (sub2, 13, consensus, gene));
    System.out.println(ascore.getScore());
    

 }
}