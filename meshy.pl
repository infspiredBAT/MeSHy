#!/usr/bin/perl -w

$| = 1;

use Getopt::Long;
use Statistics::Descriptive;
use XML::LibXML::Reader;
use LWP::Simple;
use XML::Validate;

$version              = localtime((stat($0))[9]);

$usage = "# MeSHy - unanticipated knowledge discovery through statistical ranking of MeSH term pairs
# Theodosiou T. and Darzentas N.
# $version
# http://bat.infspire.org | bat\@infspire.org

welcome

before you run MeSHy, be sure to:
- download the latest PubMed's MeSH term tree and frequencies (see code for details)
- update the total number of PubMed articles in the code
- have installed the following Perl modules:
\tGetopt::Long
\tLWP::Simple
\tXML::Validate
\tXML::LibXML::Reader
\tStatistics::Descriptive

meshy.pl

    either
--query       query you want to perform in PubMed, enclose it in double quotes
    or
--xml         results from PubMed in XML format

    other (optional) arguments
--threshold   if log_odds is between threshold, remove MeSH term [default=5.0]
--chem        use chemical terms [default=on]
--max_doc     define the maximum number of documents you want to handle [default=10,000]

";

my $pubmed_articles = 26000000; # total number of PubMed articles

# https://www.nlm.nih.gov/mesh/filelist.html
# ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/meshtrees/
# ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/meshtrees/mtrees2016.bin
my $mtree_file      = "mtrees.bin";

# https://mbr.nlm.nih.gov/Downloads.shtml
# wget https://mbr.nlm.nih.gov/Download/2016/FreqCounts/Chemical_freq_alpha.gz
# wget https://mbr.nlm.nih.gov/Download/2016/FreqCounts/MH_freq_alpha.gz
# wget https://mbr.nlm.nih.gov/Download/2016/FreqCounts/MH_SH_freq_alpha.gz
my $mh_freq_file    = "MH_freq_alpha";
my $mh_sh_freq_file = "MH_SH_freq_alpha";
my $ch_freq_file    = "Chemical_freq_alpha";

# https://www.nlm.nih.gov/cgi/mesh/2016/MB_cgi
my %mesh_cat = (
		"A"=>"Anatomy",
		"B"=>"Organisms",
		"C"=>"Diseases",
		"D"=>"Chemical and Drugs",
		"E"=>"Analytical, Diagnostic and Therapeutic Techniques and Equipment",
		"F"=>"Psychiatry and Physchology",
		"G"=>"Phenomena and Processes",
		"H"=>"Disciplines and Occupations",
		"I"=>"Anthropology, Education, Sociology and Social Phenomena",
		"J"=>"Technology, Industry, Agriculture",
		"K"=>"Humanities",
		"L"=>"Information Science",
		"M"=>"Named Groups",
		"N"=>"Health Care",
		"V"=>"Publication Characteristics",
		"Z"=>"Geographicals");

my %ids_mesh        = ();
my %score           = ();
my $cluster_pmid    = ();
my %expected        = ();
my $verbose         = 1;    # print progress info of script
my $cluster_file    = '';   # file name of the cluster file
my $query;
my $xml_filename;
my $threshold       = 5;    # log odds threshold for filtering MeSH terms
my $quali_switch    = 0;    # read and use MeSH qualifiers https://www.nlm.nih.gov/mesh/2016/QualiferTree2016.html
my $chem_switch     = 1;    # flag for reading or not the chemical terms
my $small           = 0;    # if set then compute only the first thousand articles
my $max_documents   = 10000;
my $pwd             = '.';
my $query_size      = 0;    # total number of articles in the query


my $options_res = GetOptions(
    "query:s"       =>\$query,
    "xml:s"         =>\$xml_filename,
    "max_doc:f"     =>\$max_documents,

    "chem"          =>\$chem_switch,
    "qualifiers!"   =>\$quali_switch,
    "threshold:f"   =>\$threshold,
    "small"         =>\$small,

    "mtree:s"       =>\$mtree_file,
    "mhfreqfile:s"  =>\$mh_freq_file,
    "mhshfreqfile:s"=>\$mh_sh_freq_file,
    "chfreqfile:s"  =>\$ch_freq_file,

    "clusterfile:s" =>\$cluster_file,

    "verbose!"      =>\$verbose,
    "pwd:s"         =>\$pwd
    );

if (!defined($query) && !defined($xml_filename) || (defined($query) && defined($xml_filename))) {
    die $usage;
}

$randomer = int(rand(1001)) . $$;
$randomer = sprintf "%010s", $randomer;

# check for a clustering file and get the pmids in each cluster
# fhe format should be cluster_name\tpmid
if($cluster_file ne ''){
    open(CLUSTER,"<$cluster_file") or die "Cannot read $cluster_file\n";
    while(<CLUSTER>){
        chomp;
        ($cluster, $pmid) = split(/\t/);
    #	$cluster_pmid{$cluster}{$pmid} = 1;
        $pmid_cluster{$pmid} = "$cluster";
    #	++$cluster_size{$cluster};
    }
}

# run query to PubMed using eutils
if (defined($query) =~ /\w/) {
    # Define library for the 'get' function used in the next section.
    # $utils contains route for the utilities.
    # $db, $query, and $report may be supplied by the user when prompted; 
    # if not answered, default values, will be assigned as shown below.
    
    print "(?) querying PubMed with: $query\n";

    my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
    
    my $db     = "Pubmed";
    my $report = "xml";
    my $esearch = "$utils/esearch.fcgi?" .
                  "db=$db&usehistory=y&term=";
    my $retstart;
    my $retmax=1000;
    my $results = $pwd."/meshy.query$randomer.xml";
    
    my $esearch_result = get($esearch . $query);
    
    $esearch_result =~ 
      m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;
    
    my $count    = $1;
    my $queryKey = $2;
    my $webEnv   = $3;
    
    if ($count > $max_documents){
        print "(!) we found $count documents but only the first $max_documents will be downloaded and validated\n";
        $count = $max_documents;
    } elsif ($count > 0) {
        print "(=) we found $count documents - downloading and validating\n";
    } else {
        print "(=) we found no documents - exiting\n";
        exit;
    }
    
    # this area defines a loop which will display $retmax citation results from 
    # Efetch each time the the Enter Key is pressed, after a prompt.
    for($retstart = 0; $retstart < $count; $retstart += $retmax) {
        my $left = $count - $retstart;
        print " $left\n" unless $retstart == 0;
        my $efetch = "$utils/efetch.fcgi?" .
        "retmode=$report&retstart=$retstart&retmax=$retmax&" .
        "db=$db&query_key=$queryKey&WebEnv=$webEnv";
        
        my $efetch_result = get($efetch);
        $efetch_results .= $efetch_result;
        sleep(2);
    }
    print " ok\n" if $retstart > $retmax;
    @split_efetch_results = (split /\n/, $efetch_results);
    foreach $line (@split_efetch_results) {
        unless ($line =~ /(^\<\?xml version|^\<\!DOCTYPE|^\<PubmedArticleSet|^\<\/PubmedArticleSet)/) {
            $cleaned_efetch_results .= "$line\n";
        }
    }
    $tmp2validate = <<XMLCODE;
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2015//EN" "http://www.ncbi.nlm.nih.gov/corehtml/query/DTD/pubmed_150101.dtd">
<PubmedArticleSet>
$cleaned_efetch_results
</PubmedArticleSet>
XMLCODE
    open(RES,">$results") or die "Cannot create file query_res.xml: $!\n";    
    binmode(RES, ":utf8");
    print RES $tmp2validate;
    close RES;
    $xml_filename = $results;
    $tmp2validate = do {local (@ARGV,$/) = "$xml_filename"; <>}; # file slurp!
    $validator = new XML::Validate(Type => 'LibXML');  
    if ($validator->validate($tmp2validate)) {
        print "(?) download complete and XML is valid\n";
    } else {
        print "(?) download complete but XML is invalid, this can happens over the net - please try again later, exiting\n";
        my $message = $validator->last_error()->{message};
        my $line = $validator->last_error()->{line};
        my $column = $validator->last_error()->{column};
        print "error: $message at line $line, column $column\n";
        exit;
    }    
    undef $tmp2validate;
    
} elsif (defined($xml_filename)) {

    print "(?) loading and validating PubMed XML\n";
    $tmp2validate = do {local (@ARGV,$/) = "$xml_filename"; <>}; # file slurp!
    $validator = new XML::Validate(Type => 'LibXML');  
    if ($validator->validate($tmp2validate)) {
        print "(?) XML is valid\n";
    } else {
        print "(?) XML is invalid - please try again with another file, exiting\n";
        my $message = $validator->last_error()->{message};
        my $line = $validator->last_error()->{line};
        my $column = $validator->last_error()->{column};
        print "error: $message at line $line, column $column\n";
        exit;
    }    
    undef $tmp2validate;
}

print "(?) reading PubMed XML file\n" unless $verbose == 0;

$pmid_pos = 1;
%ori_rank = ();
#$year_flag = 0;
$pub_years = ();

if ($chem_switch) { print "(?) including chemical terms\n"; }

# Read the PubMed XML file and all the relevant info
my $reader = new XML::LibXML::Reader(location => $xml_filename) or die "Cannot read $xml_filename\n";

while ($reader->read) {
    if ($reader->name eq "PMID" and $reader->depth == 3 and $reader->nodeType == 1) {
	    $reader->read;
        $pmid = $reader->value;
        $ori_rank{$pmid} = $pmid_pos;
        ++$pmid_pos;
    } elsif($reader->name eq "DateCreated"){
        $reader->read;
        $reader->read;
        if($reader->name eq "Year"){
            $reader->read;
            $pub_years{$pmid}= $reader->value;
        }
    } elsif ($reader->name eq "PublicationType" and $reader->nodeType == 1) {
    	$reader->read;
        unless ($reader->value =~ /Corrected|Duplicate|English Abstract|Journal Article|Research Support/i) {
            $pubtypes{$pmid} .= $reader->value . ";";
        }
    } elsif ($reader->name eq "DescriptorName" and $reader->nodeType == 1) {
    	$reader->read;
        $descriptor = $reader->value;
        $descriptors{$descriptor} = 1;
    } elsif ($reader->name eq "QualifierName" and $reader->nodeType == 1 and $quali_switch){
        $reader->read;
        $qualifiers{$descriptor}{$reader->value} = 1;
    } elsif ($reader->name eq "ChemicalList" and $chem_switch == 1 and $reader->nodeType == 1) {
        $chemflag = 1;
    }
    if ($reader->name eq "NameOfSubstance" && $chemflag and $reader->nodeType == 1) {
        $reader->read;        
        $descriptor = $reader->value;
        $descriptors{$descriptor} = 1;
        $mesh_tree{$descriptor}{'D'} = $mesh_cat{'D'};
    }
    if ($reader->name eq "PubmedArticle" and $reader->nodeType == 15) {
        if (defined($pubtypes{$pmid})) {
            chop $pubtypes{$pmid};
        }
        foreach $descriptor (keys %descriptors) {
            if (defined($qualifiers{$descriptor})) {
                foreach $qualifier (keys %{$qualifiers{$descriptor}}) {
                    $mem{"$descriptor/$qualifier"} = 1;
                }
            } else {
                $mem{$descriptor} = 1;
            }
        }
    
        $chemflag = 0;
        ++$query_size{all};
	    ++$query_size{$pmid_cluster{$pmid}} unless !defined($pmid_cluster{$pmid});
        # print "Processing article number: $query_size{all}\r" unless $verbose == 0; 
        foreach $descriptor_a (keys %mem) {
            ++$term_stats{all}{$descriptor_a}{in_pmids};
            # hold the number of occurences of each MeSH per cluster
	        ++$term_stats{$pmid_cluster{$pmid}}{$descriptor_a}{in_pmids} unless !defined($pmid_cluster{$pmid});
            foreach $descriptor_b (keys %mem) {
                if ($descriptor_a ne $descriptor_b && !defined($done{all}{$descriptor_b}{$descriptor_a})) {
                    if (defined($printmem{all}{$descriptor_a}{$descriptor_b})) {
                        ++$printmem{all}{$descriptor_a}{$descriptor_b};
                        $ids_mesh{all}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                    } elsif (defined($printmem{all}{$descriptor_b}{$descriptor_a})) {
                        ++$printmem{all}{$descriptor_b}{$descriptor_a};
                        $ids_mesh{all}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                    } else {
                        ++$printmem{all}{$descriptor_a}{$descriptor_b};
                        $ids_mesh{all}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                    }
                    $done{all}{$descriptor_a}{$descriptor_b} = 1;
                    $ids_mesh{all}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                }
                if (defined($pmid_cluster{$pmid}) && $descriptor_a ne $descriptor_b && !defined($done{$pmid_cluster{$pmid}}{$descriptor_b}{$descriptor_a})) {
		            if (defined($printmem{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b})) {
                        ++$printmem{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b};
			            $ids_mesh{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                    } elsif (defined($printmem{$pmid_cluster{$pmid}}{$descriptor_b}{$descriptor_a})) {
                        ++$printmem{$pmid_cluster{$pmid}}{$descriptor_b}{$descriptor_a};
			            $ids_mesh{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                    } else {
                        ++$printmem{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b};
			            $ids_mesh{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                    }
                    $done{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b} = 1;
                    $ids_mesh{$pmid_cluster{$pmid}}{$descriptor_a}{$descriptor_b}{$pmid} = 1;
                }
            }
        }
    	%mem = ();
    	%qualifiers = ();
    	%descriptors = ();
	    %done = ();
        # print only the first thousand if the small flag is set
        if ($query_size{all} > 1000 && $small == 1) {
            last;
        }
    }
}
$reader->finish;
#close XML;
#close MOD;
#print  "\nFinished reading PubMed file\n" unless $verbose == 0;

# read MeSH tree structure
open (MESHTREE,"<$mtree_file") or die "Cannot read $mtree_file\n";
$mesh_tree = ();
print "(?) reading MeSH tree structure\n";
while(<MESHTREE>){
    chomp;
    my ($mh, $pos) = split(/;/);
    my $mcategory =$mesh_cat{substr($pos,0,1)};
	$mesh_tree{$mh}{$pos} = $mcategory;
}
close MESHTREE;

# compute probability of each MeSH term in the query
foreach $group (keys %printmem){
    foreach $term (keys %{ $printmem{$group} }) {
	    $rel_freq{query}{$group}{$term} = $term_stats{$group}{$term}{in_pmids}/$query_size{$group};
	    if      ($rel_freq{query}{$group}{$term} >= 0.01)  { $pop{query}{$group}{$term} = 'pop';
	    } elsif ($rel_freq{query}{$group}{$term} <= 0.001) { $pop{query}{$group}{$term} = 'rare';
	    } else {                                    	     $pop{query}{$group}{$term} = 'norm'; }
	#    print LOG "prob query\t$group\t$rel_freq{query}{$group}{$term}\t$term\n";
    }
}

# compute the probability for each MeSH term for PubMed
open (MHFREQ, "$mh_freq_file") or die "Cannot read $mh_freq_file\n";
while (<MHFREQ>) {
    chomp;
    @split_line = (split /\|/);
    next unless defined($term_stats{all}{$split_line[4]}{in_pmids});
    $term_stats{all}{$split_line[4]}{pubmed} = $split_line[0];
    # the probability of each MeSH term to be found to PubMed and not in our query
    if ($term_stats{all}{$split_line[4]}{in_pmids} > $split_line[0]){
        #print "\tMeSH term \"$split_line[4]\" has larger freq ($term_stats{all}{$split_line[4]}{in_pmids})" .
        #" in the query than in whole PubMed ($split_line[0]) -- equalising with in-query frequency\n";
        #$remove{$split_line[4]} = 1; #remove MeSH term if we do not have a count for the whole PubMed
        $split_line[0] = $term_stats{all}{$split_line[4]}{in_pmids}; # instead of removing the MeSH term, equalise it...
    }
    $rel_freq{pubmed}{all}{minus}{$split_line[4]} = ( abs($split_line[0] - $term_stats{all}{$split_line[4]}{in_pmids}) ) / ($pubmed_articles - $query_size{all});
}
close MHFREQ;

# read qualifiers only if quali_switch is on
if($quali_switch){
open (MHSHFREQ, "$mh_sh_freq_file") or die "Cannot read $mh_sh_freq_file\n";
while (<MHSHFREQ>) {
    chomp;
    @split_line = (split /\|/);
    next unless defined($term_stats{all}{$split_line[2]}{in_pmids});
    $term_stats{all}{$split_line[2]}{pubmed} = $split_line[0];
    # the probability of each MeSH term to be found to PubMed and not in our query
    if ($term_stats{all}{$split_line[2]}{in_pmids} > $split_line[0]){
        #print "\tMeSH term \"$split_line[2]\" has larger freq ($term_stats{all}{$split_line[2]}{in_pmids})" .
        #"in the query than in whole PubMed ($split_line[0]) -- equalising with in-query frequency\n";
        #$remove{$split_line[2]} = 1; #remove MeSH term if we do not have a count for the whole PubMed
        $split_line[0] = $term_stats{all}{$split_line[2]}{in_pmids}; # instead of removing the MeSH term, equalise it...
    }
    $rel_freq{pubmed}{all}{minus}{$split_line[2]} = ( abs($split_line[0] - $term_stats{all}{$split_line[2]}{in_pmids}) ) / ($pubmed_articles - $query_size{all});
}
close MHSHFREQ;
}

open (CHFREQ, "$ch_freq_file") or die "Cannot read $ch_freq_file\n";
while (<CHFREQ>) {
    chomp;
    @split_line = (split /\|/);
    next unless defined($term_stats{all}{$split_line[2]}{in_pmids});
    $term_stats{all}{$split_line[2]}{pubmed} = $split_line[0];
    # the probability of each MeSH term to be found to PubMed and not in our query
    if ($term_stats{all}{$split_line[2]}{in_pmids} > $split_line[0]){
        #print "\tMeSH term \"$split_line[2]\" has larger freq ($term_stats{all}{$split_line[2]}{in_pmids})" .
        #"in the query than in whole PubMed ($split_line[0]) -- equalising with in-query frequency\n";
        #$remove{$split_line[2]} = 1; #remove MeSH term if we do not have a count for the whole PubMed
        $split_line[0] = $term_stats{all}{$split_line[2]}{in_pmids}; # instead of removing the MeSH term, equalise it...
    }
    $rel_freq{pubmed}{all}{minus}{$split_line[2]} = ( abs($split_line[0] - $term_stats{all}{$split_line[2]}{in_pmids}) ) / ($pubmed_articles - $query_size{all});
}
close CHFREQ;

# compute logOdds for each MeSH term and put in remove hash the ones with around 0 logOdds
print "(?) computing log odds for each MeSH term and filtering\n" unless $verbose == 0;
foreach $group (keys %printmem){
    foreach $term (keys %{ $printmem{$group} }){
        unless (defined($rel_freq{pubmed}{all}{minus}{$term})) {
            # equalising again, probably overestimating...
            $rel_freq{pubmed}{all}{minus}{$term} = ( abs($term_stats{all}{$term}{in_pmids} - $term_stats{all}{$term}{in_pmids}) ) / ($pubmed_articles - $query_size{all});
            #print "\tMeSH term \"$term\" has no freq in whole PubMed -- equalising with in-query frequency ($term_stats{all}{$term}{in_pmids})\n";
        }
        if (!defined($rel_freq{pubmed}{all}{minus}{$term})) {
            $log_odds = 0;
        } elsif ($rel_freq{pubmed}{all}{minus}{$term} == 0) { 
            $log_odds = 0;
        } else {
	        $log_odds = log($rel_freq{query}{$group}{$term} / $rel_freq{pubmed}{all}{minus}{$term});
	    }
        # if log_odds is between the threshold then remove the MeSH term
    	$remove{$term} = 1 if( (-$threshold < $log_odds and 0 > $log_odds) or ($threshold > $log_odds and 0 < $log_odds) );
    }
}


print "(?) scoring each MeSH pair\n" unless $verbose == 0;
foreach $group (keys %printmem) {
    foreach $a (keys %{ $printmem{$group} }) {
        next unless !defined($remove{$a}); # skip if logOdds near zero
        foreach $b (keys %{ $printmem{$group}{$a} }) {
            next unless !defined($remove{$b}); # skip if logOdds near zero

            $expected_score = $rel_freq{query}{$group}{$a} * $rel_freq{query}{$group}{$b}; # probability to find the MeSH pair in the query
            #$observed_score = $printmem{$group}{$a}{$b}; # the actual number of times we find the MeSH pair
            $observed_score = $printmem{$group}{$a}{$b} / $query_size{all}; # the relative number of times we find the MeSH pair
            #$expected{$group}{$a}{$b} = $expected_score;
            $new_score = log($observed_score / $expected_score); # the score for each MeSH pair
            if($new_score > 1){
                $annotation = "Observed less than expected";
            }
            elsif($new_score < -1) {
                $annotation = "Observed more than expected";
            }
            else {
                $annotation = "Observed as much as expected";
            }
            $mesh1only = (split '/', $a)[0];
            $mesh2only = (split '/', $b)[0];
            $semantic_score = sim($mesh1only,$mesh2only,%mesh_tree);

            $combined_score = $new_score * (1 - $semantic_score);
            $score{$group}{$combined_score}{$a}{$b} = $annotation;
            #foreach $pmid (keys %{ $ids_mesh{$group}{$a}{$b} }) {
            #    $ranked{$group}{$pmid}{$combined_score} = 1; 
            #}
            ++$term_stats{$group}{$a}{in_pairs};
            ++$term_stats{$group}{$b}{in_pairs};
            $pos1 = '';
            foreach $pos (sort keys %{$mesh_tree{$mesh1only}}) {
                $pos1 .= "$pos;";
            }
            chop $pos1;
            $pos2 = '';
            foreach $pos (sort keys %{$mesh_tree{$mesh2only}}) {
                $pos2 .= "$pos;";
            }
            chop $pos2;
        }
    }

}

# print MeSH pairs in categories. Use html format.
my $html_filename = $xml_filename . ".html";
open (HTML,">$html_filename") or die "Can not create file $html_filename";
print HTML <<HTMLCODE;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
HTMLCODE
    my $keyword;
my $type_of_input = "query";
if(defined($query)){
    $keyword = $query;
}
else {
    if($xml_filename =~ /\/meshy\/(.*)\.xml/){
	$keyword = $1;
    }
    else {
	$keyword = $xml_filename;
	$type_of_input = "PubMed XML file";
    }
}
print HTML "<title>MeSHy results for $type_of_input $keyword</title>";
# http://tablefilter.free.fr/dwn.php
print HTML <<HTMLCODE;
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<style type="text/css" media="screen, projection">
\@import "TableFilter/TF_Themes/Default/TF_Default.css";
</style>
<script src="TableFilter/tablefilter_all.js" language="javascript" type="text/javascript"></script>
<script src="TableFilter/sortabletable.js" language="javascript" type="text/javascript"></script>
<script src="TableFilter/tfAdapter.sortabletable.js" language="javascript" type="text/javascript"></script>
<style type="text/css" media="screen">
\@import "TableFilter/filtergrid.css";
div#navmenu li a#lnk01{ 
color:#333; font-weight:bold;
border-top:2px solid #ff9900;
background:#fff;
}
</style>
</head>
<body>
HTMLCODE
    my $number_documents = $pmid_pos - 1;
print HTML "<h1>MeSHy results for $type_of_input $keyword ($number_documents documents)</h1>";
print HTML <<HTMLCODE;
<table id="table1" class="mytable" cellspacing="0" cellpadding="0">
<thead>
<tr>
<th>Pair</th>
<th>MeSH 1</th><th>MeSH category</th>
<th>MeSH 2</th><th>MeSH category</th>
<th>MeSH 1 - MeSH 2</th>
<th>Score</th><th>PMIDs</th>
</tr>
</thead>
<tbody>
HTMLCODE

$count = 1;
foreach $group (sort keys %score) {
    foreach $scre (sort { $b <=> $a } keys % { $score{$group} }) {
    	foreach $mesh1 (sort keys %{ $score{$group}{$scre} }) {
    	    $mesh1only = (split '/', $mesh1)[0];
	        foreach $mesh2 (sort keys %{ $score{$group}{$scre}{$mesh1} }) {
        	    $mesh2only = (split '/', $mesh2)[0];
		        foreach $pmid (keys %{ $ids_mesh{$group}{$mesh1}{$mesh2} }){
                    push(@tmp,$pmid);
                    if (defined($pubtypes{$pmid})) {
                        $pubtypes2print = ":$pubtypes{$pmid}";
                    } else {
                        $pubtypes2print = "";
                    }
                    $pmid_html = "<a href=\"http://www.ncbi.nlm.nih.gov/pubmed?term=" . $pmid 
                    . "[uid]\" target=\"_blank\">$pmid($ori_rank{$pmid}):$pub_years{$pmid}$pubtypes2print</a>";
                    push(@tmp_html,$pmid_html);
        		}
                foreach $pos1 (keys %{ $mesh_tree{$mesh1only} }){
                    if ($pos1 =~ /\./) {
                        $depth1 = ($pos1 =~ tr/\./\./) + 1;
                    } else {
                        $depth1 = 1;
                    }
                    my $pos_str1 = "$mesh_tree{$mesh1only}{$pos1}";
                    $tmp_pos1{$pos_str1}{$depth1} = 1;
                }
                foreach $pos2 (keys %{ $mesh_tree{$mesh2only} }){
                    if ($pos2 =~ /\./) {
                        $depth2 = ($pos2 =~ tr/\./\./) + 1;
                    } else {
                        $depth2 = 1;
                    }
                    my $pos_str2 = "$mesh_tree{$mesh2only}{$pos2}";
                    $tmp_pos2{$pos_str2}{$depth2} = 1;
                }
                $pmids_html = join(" ",@tmp_html);
                foreach $pos1 (keys %tmp_pos1){
                    my $depth1 = join("-",sort {$a<=>$b} keys %{ $tmp_pos1{$pos1} });
                    $tree1 = "$pos1 ($depth1)";
                }
                foreach $pos2 (keys %tmp_pos2){
                    my $depth2 = join("-",sort {$a<=>$b} keys %{ $tmp_pos2{$pos2} });
                    $tree2 = "$pos2 ($depth2)";
                }
                $mesh1_html = "<a href=\"http://www.ncbi.nlm.nih.gov/mesh?term=" . $mesh1only . "\" target=\"_blank\">$mesh1</a>   [$pop{query}{$group}{$mesh1}]";
                $mesh2_html = "<a href=\"http://www.ncbi.nlm.nih.gov/mesh?term=" . $mesh2only . "\" target=\"_blank\">$mesh2</a>   [$pop{query}{$group}{$mesh2}]";
        
                print HTML "<tr>";
                printf HTML ("<td>$count</td>" . 
                    "<td>$mesh1_html</td><td>$tree1</td>" . 
                    "<td>$mesh2_html</td><td>$tree2</td>" . 
                    "<td>$mesh1_html - $mesh2_html</td>" .
                         "<td>%.3f</td><td>$pmids_html</td>\n",$scre);
                print HTML "</tr>\n";
                @tmp = ();
                @tmp_html=();
                %tmp_pos1 = ();
                %tmp_pos2 = ();
                $tree1 = "";
                $tree2 = "";
                ++$count;
            }
    	}
    }
}
$count = 0;
print HTML "</tbody>\n";
print HTML "</table>\n";
print HTML <<HTMLCODE;
<script language="javascript" type="text/javascript">
//<![CDATA[
var props = {  
        sort: true,  
        sort_config: {
                      sort_types:['Number','String','String','String','String','String','Number','String']
                     },  
        remember_grid_values: true,  
        status: true,
        highlight_keywords: true,
        alternate_rows: true,  
        rows_counter: true,  
        rows_counter_text: "Displayed pairs: ",  
        btn_reset: true,  
        btn_reset_text: "Clear",  
        status_bar: true,  
        col_1: "select",
        col_2: "multiple",  
        col_3: "select",
        col_4: "multiple",
        display_all_text: "< Show all >",  
        custom_slc_options: {
              cols:[2,4],
              texts: [['A:Anatomy','B:Organisms','C:Diseases','D:Chemical and Drugs','E:Analytical, Diagnostic and Therapeutic Techniques and Equipment','F:Psychiatry and Physchology','G:Phenomena and Processes','H:Disciplines and Occupations','I:Anthropology, Education, Sociology and Social Phenomena','J:Technology, Industry, Agriculture','K:Humanities','L:Information Science','M:Named Groups','N:Health Care','V:Publication Characteristics','Z:Geographicals'],['A:Anatomy','B:Organisms','C:Diseases','D:Chemical and Drugs','E:Analytical, Diagnostic and Therapeutic Techniques and Equipment','F:Psychiatry and Physchology','G:Phenomena and Processes','H:Disciplines and Occupations','I:Anthropology, Education, Sociology and Social Phenomena','J:Technology, Industry, Agriculture','K:Humanities','L:Information Science','M:Named Groups','N:Health Care','V:Publication Characteristics','Z:Geographicals']],
              values: [['Anatomy*','Organisms*','Diseases*','Chemical and Drugs*','Analytical, Diagnostic and Therapeutic Techniques and Equipment*','Psychiatry and Physchology*','Phenomena and Processes*','Disciplines and Occupations*','Anthropology, Education, Sociology and Social Phenomena*','Technology, Industry, Agriculture*','Humanities*','Information Science*','Named Groups*','Health Care*','Publication Characteristics*','Geographicals*'],['Anatomy*','Organisms*','Diseases*','Chemical and Drugs*','Analytical, Diagnostic and Therapeutic Techniques and Equipment*','Psychiatry and Physchology*','Phenomena and Processes*','Disciplines and Occupations*','Anthropology, Education, Sociology and Social Phenomena*','Technology, Industry, Agriculture*','Humanities*','Information Science*','Named Groups*','Health Care*','Publication Characteristics*','Geographicals*']],
              sorts: [false,false]
},
        col_width: ["15px","150px","200px","150px","200px","150px","70px","100px"],

extensions: { 

name:['ColsVisibility'], 
src:['TableFilter/TFExt_ColsVisibility/TFExt_ColsVisibility.js'], 
description:['Columns visibility manager'], 
initialize:[function(o){o.SetColsVisibility();}] 
},

}

var t01 = new TF('table1',props);
t01.AddGrid();
//]]>
</script>
HTMLCODE
print HTML "</body>\n";
print HTML "</html>\n";
close HTML;

print "(=) done - $xml_filename\n";

exit(0);


sub sim {
    my ($mesh1, $mesh2, %mtree) = @_;
    my $score = 0;
    my $simab = 0;
    
    my $stats = Statistics::Descriptive::Full->new();
    
    foreach $mid1 (keys %{ $mtree{$mesh1} }) {
        my $mid_stats = Statistics::Descriptive::Full->new();
        my $cca_max = 0;
        my $root1 = substr($mid1,0,1);
    
MID2:   foreach $mid2 (keys %{ $mtree{$mesh2} }) {
    
            my $root2 = substr($mid2,0,1);
            if($root1 ne $root2){
                $simab = 0;
                $stats->add_data($simab);
                next; 
            }
            else {
                my $count = 0;
                $mid1 =~ s/^\D//;
                $mid2 =~ s/^\D//;
                my @numbers1 = split(/\./,$mid1);
                my @numbers2 = split(/\./,$mid2);
                my $depth1 = scalar(@numbers1);
                my $depth2 = scalar(@numbers2);
                foreach $part1 (@numbers1){
                    if(defined($numbers2[$count]) and $part1 == $numbers2[$count]){
                        ++$count;
                    } else {
                        last;
                    }
                }
            SCORE:
                $maxdepth = max($depth1,$depth2);
                $simab = ($count+1) / ($maxdepth+1);
                $stats->add_data($simab);
                next MID2;
            }
        }
        if(scalar(keys %{ $mtree{$mesh2} }) == 0){ # if mesh2 has no mesh tree number, like Female
            $stats->add_data(0);
        }
    }
    
    my $total_ids1 = scalar(keys %{ $mtree{$mesh1} });
    my $total_ids2 = scalar(keys %{ $mtree{$mesh2} });
    if($total_ids1 == 0 or $total_ids2 == 0){
        $score = 0;
    } else {
        $score = $stats->sum() / ($total_ids1 * $total_ids2);
    }
    return $score;
}


sub max {
    my ($a, $b) = @_;
    my $maxab = 0;
    if ($a >= $b) { $maxab = $a;
    } else {        $maxab = $b; }
    return $maxab;
}
