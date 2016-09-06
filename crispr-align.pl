#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $fastq_file_path;
my $library_file_path;
my $project_name;
my $cutadapt_seq;
my $resources_file;
my $email;
my $args = scalar @ARGV;
my $help = undef;
GetOptions (
  "fastq_file_path=s" => \$fastq_file_path,
  "library_file_path=s" => \$library_file_path,
  "project_name=s" => \$project_name,
  "cutadapt_seq=s" => \$cutadapt_seq,
  "resources_file=s" => \$resources_file,
  "email=s" => \$email,
  "help" => \$help,
  );

if(!$args || $help){
	&usage;
	exit(0);
}

my %resources;
open RES, "< $resources_file" or die "Unable to read resources from $resources_file: $!\n";
while(<RES>){
	next unless /^[^\t]+\t[^\t]+$/;
	my $line = $_;
	chomp($line);
	my ($key, $value) = split(/\t/, $line);
	$resources{$key} = $value;
}
close RES;

my @fastq_files_on_data;
# find the fastq files
if($fastq_file_path =~ /,/){
	# we have multiple fastq files joined by commas
	@fastq_files_on_data = split(/,/,$fastq_file_path);
}
elsif($fastq_file_path =~ /\.gz$/){
	# we have a single fastq file
	push(@fastq_files_on_data,$fastq_file_path);
}
else{
	# attempt to find fastq files in subdirectories below the path given
	@fastq_files_on_data = glob "$fastq_file_path/*/*.fastq.gz";
}
my $fastq_file_count = scalar(@fastq_files_on_data);
print "found $fastq_file_count fastq files\n";

# check if any fastq files are in a common folder.
# If so join them and fix the lane info.

my %fastq_file_sets_on_data; # use path to dir as key and comma-sep file path(s) as values
foreach my $fastq_file_on_data (@fastq_files_on_data){
	my $fastq_dir_on_data = $fastq_file_on_data;
	$fastq_dir_on_data =~ s/[^\/]+$//;
	if(exists $fastq_file_sets_on_data{$fastq_dir_on_data}){
		$fastq_file_sets_on_data{$fastq_dir_on_data} = $fastq_file_sets_on_data{$fastq_dir_on_data} . ',' . $fastq_file_on_data;
	}
	else{
		$fastq_file_sets_on_data{$fastq_dir_on_data} = $fastq_file_on_data;
	}
}

foreach my $fastq_dir_on_data (keys %fastq_file_sets_on_data){
	
	# check if analysis folder exists and create if not
	my $sample_name = $fastq_dir_on_data;
	$sample_name =~ s/\/?$//;
	$sample_name =~ s/^.*\/+([^\/]+)/$1/;
	
	print "folder name: $fastq_dir_on_data\nsample name: $sample_name\n\n";
	
	my $analysis_dir = $resources{'crisp_screen_analysis_dir'};
	my $project_dir = "$resources{'crisp_screen_analysis_dir'}/$project_name/";
	my $sample_dir = "$resources{'crisp_screen_analysis_dir'}/$project_name/$sample_name/";
	
	mkdir($analysis_dir, 0775) unless -e $analysis_dir;
	mkdir($project_dir, 0775) unless -e $project_dir;
	mkdir($sample_dir, 0775) unless -e $sample_dir;
	
	# store the first part of the job sumbission script in a string
	my $script_header = <<"END";
#!/bin/bash
#BSUB -u $email
#BSUB -o $sample_dir/run_crispr_shalign.out
#BSUB -e $sample_dir/run_crispr_shalign.err
#BSUB -P  $resources{'budget_code'}
#BSUB -J crispr_shalign
#BSUB -n 8

module load CPAN/perl5.10

END
	
	# if we have multiple fastq files in a folder, copy them to
	# the new analysis folder, merge and gunzip. Otherwise just
	# copy and gunzip
	my @fastq_files_on_data = split(
		",",
		$fastq_file_sets_on_data{$fastq_dir_on_data}
		);
	
	my $script_body;
	
	if($#fastq_files_on_data > 0){ # we have more than one fastq file - merge
		my @path_to_copied_fastq;
		foreach my $fastq_file_on_data (@fastq_files_on_data){
			$script_body .= "cp $fastq_file_on_data $sample_dir\n";
			my $fastq_file_name = $fastq_file_on_data;
			$fastq_file_name =~ s/.+\///;
			$script_body .= "gunzip $sample_dir/$fastq_file_name\n";
			$fastq_file_name =~ s/\.gz$//; #remove .gz from end to match the name after gunzip
			push(@path_to_copied_fastq, "$sample_dir/$fastq_file_name");
		}
		my $fastq_files_to_merge_string = join(',',@path_to_copied_fastq);
		$script_body .= <<"END"

perl /scratch/breakthr/jamesc/scripts/merge-fastqs.pl \\
--fastqs_to_merge $fastq_files_to_merge_string

perl /scratch/breakthr/jamesc/scripts/shALIGN/shalign.pl \\
--project 01 \\
--fastqFiles \\
$sample_dir/merged.fastq \\
--libraries \\
$library_file_path \\
--offset 0 \\
--length 19

END
	}
	else{
		$script_body .= "cp $fastq_files_on_data[0] $sample_dir\n";
		my $fastq_file_name = $fastq_files_on_data[0];
		$fastq_file_name =~ s/^.+\///;
		$script_body .= "gunzip $sample_dir/$fastq_file_name\n";
		$fastq_file_name =~ s/\.gz$//; # remove .gz to match file name after gunzip
		$script_body .= <<"END"

perl /scratch/breakthr/jamesc/scripts/shALIGN/shalign.pl \\
--project 01 \\
--fastqFiles \\
$sample_dir/$fastq_file_name \\
--libraries \\
$library_file_path \\
--offset 0 \\
--length 19

END
		
	}
	
	open BSUB, "> $sample_dir/run_crispr_shalign.sh" or die "Unable to write job submission script to $sample_dir: $!\n";
	print BSUB "$script_header\n$script_body\n";
	close BSUB;
	
	my $bsub_message = `bsub < $sample_dir/run_crispr_shalign.sh`;
	print "$bsub_message\n";
		
}


sub usage{

	my $usage = <<"END";
Usage:
	perl crispr-align.pl [options]
	
Options:
	--fastq_file_path		path to a fastq file or a directory containing subdirectories with fastq files
	--library_file_path		path to the CRISPR guide library file
	--project_name			a string to use in the output directory name
	--cutadapt_seq			a sequence to pass to cutadapt if needed
	--resources_file		path to a file listing information such as budget codes 

END

	print $usage;

}


