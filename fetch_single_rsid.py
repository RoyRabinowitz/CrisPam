import copy
from Bio import Entrez, SeqIO
from Scripts.fetch_genomic_seq import *

__author__ = "Shiran Abadi"


def get_snp(snp_id, temp_dir):
	"""
	This function takes as input a snp identifiers and a temporary directory to save temp files to
	and returns a parsed dictionary of their data from Entrez.
	"""
	Entrez.email = 'A.N.Other@example.com'

	try:
		handle = Entrez.efetch(db='SNP', id=snp_id)
		snp_info = handle.read()
		handle.close()
	except:
		print("Failed to fetch SNP.")
		return None

	r = {} # Return dictionary variable
	field_names = ["SNP_ID", "CLINICAL_SIGNIFICANCE", "GENE_ID", "CHRPOS", "FXN_CLASS", "DOCSUM"]
	for field_name in field_names:
		try:
			r[field_name] = re.search("<" + field_name + ">(.*?)</" + field_name + ">", snp_info).group(1)
		except:
			r[field_name] = None

	res = {}
	res["formal_id"] = r["SNP_ID"]
	res["ClinicalSignificance"] = r["CLINICAL_SIGNIFICANCE"]
	res["GENE_ID"] = r["GENE_ID"]
	try:
		loc_attributes = re.search("(.*):(.*)", r["CHRPOS"])
	except:
		return []
	res["chromosome"] = loc_attributes.group(1)
	res["coordinate"] = int(loc_attributes.group(2))
	res["fxnClass"] = r["FXN_CLASS"]

	res["genomic_flanking"] = \
			fetch_dna_coordinates("hg38", res["chromosome"],
			                              res["coordinate"] - 25, res["coordinate"] + 25, temp_dir)

	## DOCSUM may have multiple mutations due to alternative splicing, and might look like:
	## <DOCSUM>HGVS=NC_000011.10:g.13492506G&gt;A,NC_000011.10:g.13492506G&gt;T,NC_000011.9:g.13514053G&gt;
	# A,NC_000011.9:g.13514053G&gt;T,NG_008962.1:g.8515C&gt;T,NG_008962.1:g.8515C&gt;A,NM_000315.4:c.247C&gt;
	# T,NM_000315.4:c.247C&gt;A,NM_000315.3:c.247C&gt;T,NM_000315.3:c.247C&gt;A,NM_000315.2:c.247C&gt;
	# T,NM_000315.2:c.247C&gt;A,NM_001316352.1:c.343C&gt;T,NM_001316352.1:c.343C&gt;
	# A,NP_000306.1:p.Arg83Ter,NP_001303281.1:p.Arg115Ter|SEQ=[G/A/T]|GENE=PTH:5741</DOCSUM>
	## So I take every record that begins with a nucleotide and then NM (mRNA record) and compute the reading frame according to the position that follows.
	## For example, A,NM_001316352.1:c.343C&gt means that there was an SNV to A at position 343 of the mRNA

	m = re.finditer("([ACGT]),(NM_.*?):c\.([0-9\-]+)", r["DOCSUM"])
	orig_nuc = re.search("SEQ\=\[([ACGT])", r["DOCSUM"]).group(1)
	mutations = {}
	for match in m:
		try:
			mut = (match.group(1), int(match.group(3)))
			mutations[mut] = mutations.get(mut, []) + [match.group(2)]
		except:
			pass # if we reached here, the coordinate is in the UTRs (noncoding part that might be in the NM query)

	# retrieve a valid codon from one of the mrna sequences
	all_snps = []
	for k, v in mutations.items():
		try:
			handle = Entrez.efetch(db='nucleotide', id=",".join(v), rettype="gb", retmode="text")
			response_item = SeqIO.parse(handle, "gb").__next__() #only need one
			for feature in response_item.features:
				if feature.type == "CDS":
					cds_start_index = feature.location.start
					cds_end_index = feature.location.end
					break

			mrna = response_item.seq[cds_start_index:cds_end_index]
			snp_idx = k[1]-1

			all_snps.append(copy.deepcopy(res))
			all_snps[-1]["5UTR"] = k[1] < 0
			all_snps[-1]["3UTR"] = k[1] > (cds_end_index - cds_start_index)

			if not (all_snps[-1]["5UTR"] or all_snps[-1]["3UTR"]):
				orig_codon = mrna[(snp_idx // 3) * 3:(snp_idx // 3 + 1) * 3]
				all_snps[-1]["original_codon"] = str(orig_codon)
				all_snps[-1]["reading_frame"] = snp_idx % 3 + 1
				all_snps[-1]["orientation"] = "+" if mrna[snp_idx]==orig_nuc else "-"
				all_snps[-1]["SNP"] = mrna[snp_idx] + ">" + k[0]
				all_snps[-1]["CDS_flanking"] = mrna[max(0,snp_idx-25):min(snp_idx+27, len(mrna))]
			else:
				all_snps[-1]["SNP"] = response_item.seq[cds_start_index+snp_idx+1] + ">" + k[0]

			handle.close()
		except:
			pass
	return all_snps


if __name__ == '__main__':
	snp_id = "rs397515468"
	temp_dir = ""
	snp_info = get_snp(snp_id, temp_dir)
	print (snp_info)
