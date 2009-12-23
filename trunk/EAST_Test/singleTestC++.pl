#!/usr/bin/perl
use strict;

my $base_dir = "/home/zhangy9";

my $oriGeneFile = $ARGV[0];
my $numOfEst = $ARGV[1];
my $errorRate = $ARGV[2];
my $id = $ARGV[3];
my $useQualityFile = $ARGV[4]; #1-use quality file for Cap3, 0-not use.
my $output; #store all the analysis results

my $geneFile = "$oriGeneFile\_$id";
my $estFile = "estFile_$id.fa";
my $mstFile = "mstFile_$id.fa";
my $eastConsensusFile = "$estFile.east.contigs";
my $eastSingletonFile = "$estFile.east.singlets";
my $eastNum = "$estFile.east.num"; #number of used ests
my $eastDebugFile = "$estFile.east.debug";

my $cap3ConsensusFile = "$estFile.cap.contigs";
my $cap3SingletonFile = "$estFile.cap.singlets";
my $cap3DebugFile = "$estFile.cap.debug";

	my @oriGene;
	$oriGene[0] = ">genefile";
	$oriGene[1] = getGene($oriGeneFile);
	chdir "$base_dir/peace";

	writeContigToFile($geneFile, @oriGene);

	my @contigs;
	my $inContigFile;

	chdir "$base_dir/peace";
	system("python wrapper.py $geneFile $numOfEst $errorRate $mstFile $estFile > peace.debug");
	chdir "$base_dir/ESTAssemblyC++";
	system("mv $base_dir/peace/$geneFile .");
	system("cp -f $base_dir/peace/$estFile .");
	system("mv $base_dir/peace/$mstFile .");
	my $time1 = time();
	#must ouput to a file, or else it will be printed to the screen and be written into the final output file.
	system("./Main $estFile $mstFile $eastConsensusFile $eastSingletonFile $eastNum > $eastDebugFile");
	my $time2 = time();
	my @time;
	$time[0] = $time2 - $time1; #EAST time
	
	chdir "$base_dir/cap3";
	system("mv $base_dir/peace/$estFile .");
	if ($useQualityFile == 1) {
		system("perl writeQualityScore.pl $estFile $errorRate");
	}
	$time1 = time();
	system("cap3 $estFile > $cap3DebugFile"); #must ouput to a file, or else it will be printed to the screen and be written into the final output file.
	$time2 = time();
	$time[1] = $time2 - $time1; #CAP3 time
	
	#put contigs from CAP3 into a file, one contig in one line.
	chdir "$base_dir/cap3";
	$inContigFile = $cap3ConsensusFile;
	my $outContigCap = "cap3.$id.out";
	@contigs = getConsensus($inContigFile);
	writeContigToFile($outContigCap, @contigs); #write contigs
	system("mv $outContigCap $base_dir/ESTAssemblyC++/");
	$inContigFile = $cap3SingletonFile;
	my $outSingCap = "cap3Singleton.$id.out";
	@contigs = getConsensus($inContigFile);
	writeContigToFile($outSingCap, @contigs); #write singletons
	system("mv $outSingCap $base_dir/ESTAssemblyC++/");
	
	#put contigs from EAST into a file, one contig in one line.
	chdir "$base_dir/ESTAssemblyC++";
	$inContigFile = $eastConsensusFile;
	my $outContigEast = "east.$id.out";
	@contigs = getConsensus($inContigFile);
	writeContigToFile($outContigEast, @contigs); #write contigs

	$inContigFile = $eastSingletonFile;
	my $outSingEast = "eastSingleton.$id.out";
	@contigs = getConsensus($inContigFile);
	writeContigToFile($outSingEast, @contigs); #write singletons

	$inContigFile = $eastNum;
	@contigs = getConsensus($inContigFile);
	my $numOfUsedEsts = $contigs[0];
	
	#compare cap3 to EAST, put the analysis into a file.
	chdir "$base_dir/ESTAssemblyC++";
	my @gene = getConsensus($geneFile);
	writeContigToFile($geneFile, @gene);
	my $cap3numOfUsedEsts = 0;
	my $configFile = "analysis.properties";
	my $analysisOut = "analysis.$id.out";
	system("java eSTAssembly.ResultAnalysis $time[1] $time[0] $cap3numOfUsedEsts $numOfUsedEsts eSTAssembly/$configFile $geneFile $outContigEast $outSingEast $outContigCap $outSingCap $analysisOut >> $eastDebugFile");

	
	#copy the intermediate files for debugging later
	#system("cp -f $geneFile $base_dir/analysis/.");
	#system("cp -f $eastDebugFile $base_dir/analysis/.");
	#system("cp -f $mstFile $base_dir/analysis/.");
	#system("cp -f $estFile $base_dir/analysis/.");
	#system("cp -f $eastConsensusFile $base_dir/analysis/.");
	#system("cp -f $base_dir/cap3/$cap3ConsensusFile $base_dir/analysis/.");
	
	my $analysis = (getConsensus($analysisOut))[0];
	print $analysis;

	#remove the intermediate files. Notice: CANNOT remove the original gene file $oriGeneFile, it will be used by other jobs.
	system("rm $geneFile");
	system("rm $estFile");
	system("rm $mstFile");
	system("rm $eastConsensusFile");
	system("rm $eastSingletonFile");
	system("rm $eastNum");
	system("rm $eastDebugFile");
	system("rm $outContigEast");
	system("rm $outSingEast");
	system("rm $analysisOut");
	system("rm $outContigCap");
	system("rm $outSingCap");
	system("rm $base_dir/cap3/$estFile*");

# read gene from the input gene file, put the gene in one line.
sub getGene {
	my $fileName = $_[0];
	open INPUTFILE, "<$fileName" or die ("Input file $fileName can not be opened!");
	my $gene;
	my $line;
	
	while (defined($line = <INPUTFILE>)) {
		chomp $line;
		if (!($line =~ />(.)*/i)) {
			#print ("$line\n");
			$gene = $gene . trim($line);
		}
	};

	close INPUTFILE;
	return $gene;
}

# read consensus from the input file, one consensus in one line.
# one parameter: the name of the input file. The file is in fasta format.
# output: an array which includes all the consensus read from the input file in fasta format.
sub getConsensus {
	my $contigFile = $_[0];
	open INPUTFILE, "<$contigFile" or die ("Input file $contigFile can not be opened!");
	my @genes;
	my $contig = "";
	my $line;
	my $index = 0;
	
	while (defined($line = <INPUTFILE>)) {
		chomp $line;
		if ($line =~ />(.)*/i) { #a line of comment
			if (!$contig eq "") {
				$genes[$index++] = $contig;
				$contig = "";
			}
		} else {
			$contig = $contig . trim($line); #put one contig in one line
		}
	};
	if (!$contig eq "") {
		$genes[$index] = $contig;
	}
	
	close INPUTFILE;
	return @genes;
}

# write all the contigs into an output file, one contig in  one line.
# two parameter: output file name, an array which includes all the contigs.
sub writeContigToFile {
	my $fileName = shift @_;
	my @contigs = @_;
	my $size = @contigs;
	open OUTPUTFILE, ">$fileName" or die ("Fail to open the output file $fileName!"); 
	for (my $i=0; $i<$size; $i++) {
		print OUTPUTFILE "$contigs[$i]\n";
	}
	
	close OUTPUTFILE;
}

#delete preceding and trailing white spaces
sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
