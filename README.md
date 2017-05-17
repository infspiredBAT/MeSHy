# MeSHy
unanticipated knowledge discovery through statistical ranking of MeSH term pairs

</br>
server: http://bat.infspire.org/tools/meshy/</br>
citation: http://www.ncbi.nlm.nih.gov/pubmed/21684350</br>
code: we provide the main Perl script, and the server's CGI script and HTML page - you'll need a few more things...</br>
usage:
<pre>
before you run MeSHy, be sure to:
- download the latest PubMed's MeSH term tree and frequencies (see code for details)
- update (in the code) the total (approximate) number of PubMed articles
- have installed the following Perl modules:
	Getopt::Long
	LWP::Simple
	XML::LibXML::Reader
	Statistics::Descriptive
+ we use TableFilter (now: https://koalyptus.github.io/TableFilter/) for HTML output

meshy.pl

   either
--query       query you want to perform in PubMed, enclose it in double quotes
   or
--xml         results from PubMed in XML format

   other (optional) arguments
--threshold   if log_odds is between threshold, remove MeSH term [default=5.0]
--chem        use chemical terms [default=on]
--max_doc     define the maximum number of documents you want to handle [default=10,000]
</pre>
