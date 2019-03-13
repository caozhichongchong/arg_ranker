#!/usr/bin/perl -w
use strict;

##Authur Xiao-Tao Jiang
##Email biofuture.jiang@gmail.com.
##Date: 2016/3/6
##Modified 2017-07-19 
##Description
##This pipeline is designed to process multisamples ARG identification, this is the part two pipeline

die " perl $0 <Blast6out> <Meta_data_info> <lenth> <e-value> <identity> <out_prefix>  <SARG_STRUCTURE> <SARG_FA> Authur: Xiao-Tao Jiang\n Email: biofuture.jiang\@gmail.com\n" unless (@ARGV == 8);
#-------------------------------------------------------------------------------------------------
#blastx aginst ARG database for accurately identification of reads for antibiotic resistence gene
my $efa = $ARGV[5];
my $blast6out = "$ARGV[0]";
#my $ARDB_STRUCTURE = "/workhome/JXT/GitProject/updation_20161224/DB/structure_20161119.list";
my $ARDB_STRUCTURE = $ARGV[6];
#my $ARDBFA = "/workhome/JXT/GitProject/updation_20161224/DB/SARG-20161119.fasta";
my $ARDBFA = $ARGV[7];
my $lenmatch = $ARGV[2];
my $evaluematch = $ARGV[3];
my $identitymatch = "$ARGV[4]";

##gene abundance ppm, against 16S and Cell number three tables  
my $geneppm = "$efa.normalize_ppm.gene.tab.txt";
my $subtypeppm = "$efa.normalize_ppm.subtype.tab.txt";
my $typeppm = "$efa.normalize_ppm.type.tab.txt";

#for 16s normalization
my $gene16s = "$efa.normalize_16s.gene.tab.txt";
my $subtype16s = "$efa.normalize_16s.subtype.tab.txt";
my $type16s = "$efa.normalize_16s.type.tab.txt";

#for cell number normalization
my $genecellnumber = "$efa.normalize_cellnumber.gene.tab.txt";
my $subtypecellnumber = "$efa.normalize_cellnumber.subtype.tab.txt";
my $typecellnumber = "$efa.normalize_cellnumber.type.tab.txt";

my $rlen = 100;  ## This reads length for pair-end sequencing 



##process blastx results and the structure information of arg database 

##process meta data-------------------------------------------------------------------------------
die "$!\n" unless open(Meta,"$ARGV[1]");
my %sample2reads;
my %sample216s;
my %sample2cellnumber;

my $headmeta = <Meta>;
my @hmeta = split(/\t/,$headmeta);

while(<Meta>){
	chomp;
	my @tt = split(/\t/,$_);
	$sample2reads{$tt[1]} = $tt[-3];
	#print "$tt[-3]\n";
	$sample216s{$tt[1]} = $tt[-2];
	$sample2cellnumber{$tt[1]} = $tt[-1];
}
close Meta;

#process ARDB to get the length information 
my %len;
die "$!\n" unless open(LEN, "$ARDBFA");
while(my $name = <LEN>){
	chomp($name);
	$name =~ s/^>//;
	my $seq = <LEN>; chomp($seq);
	my $idsarg = (split(/\s+/,$name))[0];
	my $le = length($seq);
	$len{$idsarg} = $le;
}
close LEN;

##process ARDB structure files------------------------------------------------------------------- 
die "$!\n" unless open(STRU, "$ARDB_STRUCTURE");
my %type;
my %subtype;
my %typelist;
my %subtypelist;
my %genelist;


<STRU>;
while(<STRU>){
	chomp;
	my @tem = split /\t/;
	my @stem = split("__", $tem[0]);
	$tem[1] =~ s/^\[//;
	$tem[1] =~ s/\]$//;
	my @ids = split(", ", $tem[1]);
	##for each ids identify their type and subtype
	for(my $i = 0; $i <=$#ids; $i++){
		$ids[$i] =~ s/^\'//;
		$ids[$i] =~ s/\'$//;
		#print "$ids[$i]\n";
		$subtype{$ids[$i]} = $tem[0];
		$type{$ids[$i]} = $stem[0]; 
		$genelist{$ids[$i]} = 1;
	}

	#including all type and subtype
	$typelist{$stem[0]} = 1;
	$subtypelist{$tem[0]} = 1;
}
close STRU;

##parse blast6out results-----------------------------------------------------------------------
die "$!\n" unless open(BLAST6, "$blast6out"); 
my %samplehit; #Hash->Hash  sample->ARGs type/subtype->number of this ARG
my %argssamplereads;
my $upper = "";
while(<BLAST6>){
	chomp;
	my @tem = split /\t/;
	if($tem[3] >= $lenmatch && $tem[2] >= $identitymatch && $evaluematch <= 1e-5 ){
		my $seqn = $tem[0];

		if($type{$tem[1]} && $subtype{$tem[1]}){
		}else{
		print "$tem[0] $tem[1]\t $!\n"; next;
		}
		$tem[0] =~ s/\_(\d+)$//g;
		#print "$tem[0]\n" unless(exists $len{$tem[1]});
		my $ratio = 1 * $rlen / ($len{$tem[1]} * 3);   ##alignment length ratio to the total length of that ARG
		next if ($seqn eq $upper);		
	
		if(exists $samplehit{$tem[0]}){

			$samplehit{$tem[0]}{$tem[1]} += $ratio;
			$samplehit{$tem[0]}{$type{$tem[1]}} += $ratio;
			$samplehit{$tem[0]}{$subtype{$tem[1]}} += $ratio;
		}else{
			##gene
			$samplehit{$tem[0]}{$tem[1]} = $ratio;
			##type
			$samplehit{$tem[0]}{$type{$tem[1]}} = $ratio;
			##subtype
			$samplehit{$tem[0]}{$subtype{$tem[1]}} = $ratio;

		}

		if(exists $argssamplereads{$tem[0]}){

			$argssamplereads{$tem[0]}{$tem[1]} += 1;
			$argssamplereads{$tem[0]}{$type{$tem[1]}} += 1;
			$argssamplereads{$tem[0]}{$subtype{$tem[1]}} += 1;
		}else{
			##gene
			$argssamplereads{$tem[0]}{$tem[1]} = 1;
			##type
			$argssamplereads{$tem[0]}{$type{$tem[1]}} = 1;
			##subtype
			$argssamplereads{$tem[0]}{$subtype{$tem[1]}} = 1;

		}
		$upper = $seqn;
	}

}
close BLAST6;

#--------------------------------------------Gene normalization-----------16S Normalization-------Cell Number Normalization--------
##For each ARG type subtype generate mothor tables------------------------------------------------
##Hash -> Hash


die "$!\n" unless open(GENEP, ">$geneppm");
die "$!\n" unless open(SUBP, ">$subtypeppm");
die "$!\n" unless open(TYPEP, ">$typeppm");

die "$!\n" unless open(GENEM, ">$gene16s");
die "$!\n" unless open(SUBM, ">$subtype16s");
die "$!\n" unless open(TYPEM, ">$type16s");
#------

die "$!\n" unless open(GENEC, ">$genecellnumber");
die "$!\n" unless open(SUBC, ">$subtypecellnumber");
die "$!\n" unless open(TYPEC, ">$typecellnumber");
#------

print GENEP "Gene abundance to million reads gene ppm\n\tSUBTYPE\tTYPE";
print GENEM "Gene abundance to 16s  ppm\n\tSUBTYPE\tTYPE";
print GENEC "Gene abundance to cellnumber ppm\n\tSUBTYPE\tTYPE";

print SUBP "subtype abundance to million reads\n";
print SUBM "subtype abundance to 16s reads\n";
print SUBC "subtype abundance to cellnumber\n";

print TYPEP "Type abundance to million reads\n";
print TYPEM "Type abundance to 16S\n";
print TYPEC "Type abundance to cellnumber\n";

for my $id (sort keys %sample2reads){
	
	
	print GENEP "\t$id";
	print GENEM "\t$id";
	print GENEC "\t$id";
	


	print SUBP "\t$id";
	print SUBM "\t$id";
	print SUBC "\t$id";

	print TYPEP "\t$id";
	print TYPEM "\t$id";
	print TYPEC "\t$id";

}
print GENEP "\n";
print GENEM "\n";
print GENEC "\n";



print SUBP "\n";
print SUBM "\n";
print SUBC "\n";

print TYPEP "\n";
print TYPEM "\n";
print TYPEC "\n";

##output gene mother table 

for my $gene (sort keys %genelist){
	print GENEP "$gene\t$subtype{$gene}\t$type{$gene}";
	print GENEM "$gene\t$subtype{$gene}\t$type{$gene}";
	print GENEC "$gene\t$subtype{$gene}\t$type{$gene}";
	for my $sam(sort keys %sample2reads){
		if(exists $samplehit{$sam}{$gene}){

			my $value = $samplehit{$sam}{$gene} /  $sample216s{$sam};

			my $valuecls = $samplehit{$sam}{$gene} / $sample2cellnumber{$sam};

			my $valuep = $argssamplereads{$sam}{$gene} * 1000000 / $sample2reads{$sam};

			print GENEM "\t$value";
			print GENEC "\t$valuecls";
			print GENEP "\t$valuep";
		}else{
			print GENEM "\t0";
			print GENEC "\t0";
			print GENEP "\t0";
		}
	}
	print GENEM "\n";
	print GENEC "\n";
	print GENEP "\n";
}
close GENEM;
close GENEC;
close GENEP;



##output subtype mothor table
for my $sub (sort keys %subtypelist){
	print SUBM "$sub";
	print SUBC "$sub";
	print SUBP "$sub";
	for my $sam(sort keys %sample2reads){
		if(exists $samplehit{$sam}{$sub}){

			my $value = $samplehit{$sam}{$sub} /  $sample216s{$sam};

			my $valuecls = $samplehit{$sam}{$sub} / $sample2cellnumber{$sam};

			my $valuep = $argssamplereads{$sam}{$sub} * 1000000 / $sample2reads{$sam};

			print SUBM "\t$value";
			print SUBC "\t$valuecls";
			print SUBP "\t$valuep";
		}else{
			print SUBM "\t0";
			print SUBC "\t0";
			print SUBP "\t0";
		}
	}
	print SUBM "\n";
	print SUBC "\n";
	print SUBP "\n";
}
close SUBM;
close SUBC;
close SUBP;
##output type mother table 
for my $ty (sort keys %typelist){

	print TYPEM "$ty";
	print TYPEC "$ty";
	print TYPEP "$ty";

	for my $sam(sort keys %sample2reads){
		if(exists $samplehit{$sam}{$ty}){
			my $value = $samplehit{$sam}{$ty}  / $sample216s{$sam};
			print TYPEM "\t$value";

			my $valuecls = $samplehit{$sam}{$ty} / $sample2cellnumber{$sam};
			print TYPEC "\t$valuecls";

			my $valuep = $argssamplereads{$sam}{$ty} * 1000000 / $sample2reads{$sam};
			print TYPEP "\t$valuep";

		}else{
			print TYPEM "\t0";
			print TYPEC "\t0";
			print TYPEP "\t0";
		}
	}
	print TYPEM "\n";
	print TYPEC "\n";
	print TYPEP "\n";
}
close TYPEM;
close TYPEC;
close TYPEP;

1;
