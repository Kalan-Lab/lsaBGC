import os
import sys

expected_header = set(['accession', 'ambiguous_bases', 'checkm_completeness', 'checkm_contamination', 'checkm_marker_count',
				   'checkm_marker_lineage', 'checkm_marker_set_count', 'checkm_strain_heterogeneity', 'coding_bases',
				   'coding_density', 'contig_count', 'gc_count', 'gc_percentage', 'genome_size',
				   'gtdb_genome_representative', 'gtdb_representative', 'gtdb_taxonomy',
				   'gtdb_type_designation_ncbi_taxa', 'gtdb_type_designation_ncbi_taxa_sources',
				   'gtdb_type_species_of_genus', 'l50_contigs', 'l50_scaffolds', 'longest_contig', 'longest_scaffold',
				   'lsu_23s_contig_len', 'lsu_23s_count', 'lsu_23s_length', 'lsu_23s_query_id', 'lsu_5s_contig_len',
				   'lsu_5s_count', 'lsu_5s_length', 'lsu_5s_query_id', 'lsu_silva_23s_blast_align_len',
				   'lsu_silva_23s_blast_bitscore', 'lsu_silva_23s_blast_evalue', 'lsu_silva_23s_blast_perc_identity',
				   'lsu_silva_23s_blast_subject_id', 'lsu_silva_23s_taxonomy', 'mean_contig_length',
				   'mean_scaffold_length', 'mimag_high_quality', 'mimag_low_quality', 'mimag_medium_quality',
				   'n50_contigs', 'n50_scaffolds', 'ncbi_assembly_level', 'ncbi_assembly_name', 'ncbi_assembly_type',
				   'ncbi_bioproject', 'ncbi_biosample', 'ncbi_contig_count', 'ncbi_contig_n50', 'ncbi_country',
				   'ncbi_date', 'ncbi_genbank_assembly_accession', 'ncbi_genome_category', 'ncbi_genome_representation',
				   'ncbi_isolate', 'ncbi_isolation_source', 'ncbi_lat_lon', 'ncbi_molecule_count', 'ncbi_ncrna_count',
				   'ncbi_organism_name', 'ncbi_protein_count', 'ncbi_refseq_category', 'ncbi_rrna_count',
				   'ncbi_scaffold_count', 'ncbi_scaffold_l50', 'ncbi_scaffold_n50', 'ncbi_scaffold_n75',
				   'ncbi_scaffold_n90', 'ncbi_seq_rel_date', 'ncbi_spanned_gaps', 'ncbi_species_taxid',
				   'ncbi_ssu_count', 'ncbi_strain_identifiers', 'ncbi_submitter', 'ncbi_taxid', 'ncbi_taxonomy',
				   'ncbi_taxonomy_unfiltered', 'ncbi_total_gap_length', 'ncbi_total_length', 'ncbi_translation_table',
				   'ncbi_trna_count', 'ncbi_type_material_designation', 'ncbi_ungapped_length', 'ncbi_unspanned_gaps',
				   'ncbi_wgs_master', 'protein_count', 'scaffold_count', 'ssu_contig_len', 'ssu_count',
				   'ssu_gg_blast_align_len', 'ssu_gg_blast_bitscore', 'ssu_gg_blast_evalue',
				   'ssu_gg_blast_perc_identity', 'ssu_gg_blast_subject_id', 'ssu_gg_taxonomy', 'ssu_length',
				   'ssu_query_id', 'ssu_silva_blast_align_len', 'ssu_silva_blast_bitscore', 'ssu_silva_blast_evalue',
				   'ssu_silva_blast_perc_identity', 'ssu_silva_blast_subject_id', 'ssu_silva_taxonomy',
				   'total_gap_length', 'trna_aa_count', 'trna_count', 'trna_selenocysteine_count'])

file1 = sys.argv[1] # provide bac120_metadata_rXXX.tsv file

print('GenBank_ID\tGTDB_Genus\tGTDB_Species')
with open(file1) as of:
	for i, line in enumerate(of):
		line = line.strip('\n')
		ls = line.split('\t')
		if i == 0:
			curr_header = set(ls)
			if not (len(curr_header.difference(expected_header)) == 0 and len(expected_header.difference(curr_header)) == 0):
				sys.stderr.write('Headers not matching expectations for input file!')
				sys.exit(1)
		else:
			genbank_acc = ls[54]
			genus = ls[16].split(';g__')[1].split(';s__')[0].strip()
			species = ls[16].split(';s__')[1].strip()
			assert(genus != '' and species != '')
			print(genbank_acc + '\t' + genus + '\t' + species)