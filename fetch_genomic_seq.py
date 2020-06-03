from Bio.Seq import Seq
from xml.dom.minidom import parseString
import urllib.request
import os, re
from Bio.Alphabet import generic_dna
__author__ = "Shiran Abadi"


def download_file(url, filename):
     #print >> sys.stderr, "filename:" + filename
	if not os.path.exists(filename):
		 f = urllib.request.urlretrieve(url,filename)
	return filename

	

def get_xml_filename(genome, chrom, startpos, endpos, cache_dir):

	current_filename = "{0}/{1}_{2}_{3}_{4}.xml".format(cache_dir, genome, chrom, str(startpos), str(endpos))
	return current_filename


def get_dna_coordinates_xmlfile(genome, chromosome, startpos, endpos, cache_dir):
	url = "http://genome.ucsc.edu/cgi-bin/das/{0}/dna?segment={1}:{2},{3}"
	print (url)

	#chr15:65637530,65637553"
	current_filename = get_xml_filename(genome, chromosome, startpos, endpos, cache_dir)
	#current_filename='user_files/results.csv'
	current_url = url.format(genome, chromosome, startpos, endpos)
	print(current_url)
	download_file(current_url, current_filename)
	return current_filename


def fetch_dna_coordinates(genome, chrom, startpos, endpos, cache_dir):
	current_filename = get_dna_coordinates_xmlfile(genome, chrom, startpos, endpos, cache_dir)
	with open(current_filename) as fp:
		xmldata = parseString(fp.read())
	seq = re.sub("\s", "", xmldata.childNodes[1].childNodes[1].childNodes[1].childNodes[0].data.upper())
	return seq


def get_seq_by_orientation(seq, strand):
	"""
		returns the original sequence if it's the plus strand, o/w returns the reverse-complement
	"""
	if strand == "-":
		seqobject = Seq(seq, generic_dna)
		seq = str(Seq.reverse_complement(seqobject))
	return seq

