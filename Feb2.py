#!/usr/bin/env python3
import regex as re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import subprocess
from multiprocessing import Pool


def filelines_to_list(file):
    """Makes a list of each new line in a file. """
    with open(file, 'rU') as file_handle:
        file_list = [ind.rstrip() for ind in file_handle.readlines()]
    return file_list


def out_dir_maker(project_path):
    """Creates output directories for the stages of the pipeline. """
    dirs_to_make = [f"{project_path}/Blast_databases", f"{project_path}/Blast_results", f"{project_path}/Txtms_TrimtoTargets_alignments",
                    f"{project_path}/Blast_by_gene", f"{project_path}/Exonerate_out", f"{project_path}/Exonerate_clean"]
    for directory in dirs_to_make:
        try:
            os.mkdir(directory)
        except:
            pass
    return dirs_to_make


def blast_run(evalue, threads, in_path, query_file, out_path):
    """Run tblastn on each file in in_path with set evalue cutoff on x threads.
    Only saves top blast hit. Infile must end with '.nt'.
    Consecutive runs will overwrite the blast_output file.
    """
    filenames = os.listdir(in_path)
    for file in filenames:
        if file.endswith(".nt"):
            file_prefix = file.replace('.Trinity.annotated.nt', '')
            blast_output = subprocess.getoutput(
                f"tblastn -query {query_file} -db {in_path}/{file} -evalue {evalue} -outfmt '6 sseqid qseqid' -max_target_seqs 1 -max_hsps 1 -num_threads {threads}")
            with open(f"{out_path}/{file_prefix}.blast.txt", "w") as out_handle:
                out_handle.write(blast_output)
    return


def seq_fetcher(list_dir, list_file_endings, seq_dir, seq_file_ending, org_file):
    """Loads the blast results of each organism and writes the results by gene.
    Appends to existing files if present."""
    orglist = filelines_to_list(org_file)
    for org in orglist:
        hit_records = []
        record_dict = SeqIO.to_dict(SeqIO.parse(
            f"{seq_dir}/{org}{seq_file_ending}", "fasta"))
        tofetch = filelines_to_list(f"{list_dir}/{org}{list_file_endings}")
        for id in tofetch:
            record_id, gene_id = id.split()
            record = record_dict[record_id]
            record_i = SeqRecord(
                Seq(f"{record.seq}", SingleLetterAlphabet()), id=gene_id)
            hit_records.append(record_i)
        with open(f'{list_dir}/{org}.fa', 'w') as outhandle:
            SeqIO.write(hit_records, outhandle, "fasta")
    return


def cat_by_gene(org_file, loci_list, in_path, file_ending, out_path, cutoff):
    """This concatenates sequence files for all organisms by loci.
    Sequence files must have fasta heading of '>org-loci'
    The name of the alignments with fewer sequences than 'cutoff' are returned in check_list."""
    check_list = []
    org_list = filelines_to_list(org_file)
    org_dict = {}
    for org in org_list:
        exonerate_dict = SeqIO.to_dict(SeqIO.parse(
            f"{in_path}/{org}{file_ending}", "fasta"))
        org_dict[org] = exonerate_dict
    spurgenes = loci_list
    # print(f"Number of genes: {len(spurgenes)}")
    # print(f"Number of transcriptomes: {len(org_list)}")
    for spurgene in spurgenes:
        out_file_path=f"{out_path}/{spurgene}.FNA"
        with open(out_file_path, 'w') as out_handle:
            for org in org_list:
                record_id = f"{org}-{spurgene}"
                try:
                    record_i = SeqRecord(Seq(f"{org_dict[org][spurgene].seq}", SingleLetterAlphabet(
                    )), id=record_id, description="")
                    SeqIO.write(record_i, out_handle, "fasta")
                except:
                    pass
        with open(out_file_path, 'r') as out_handle:
            records=list(SeqIO.parse(out_handle,"fasta"))
            if len(records) < cutoff:
                check_list.append(spurgene)
    return check_list


def exoneratetor_bygene(query_file, in_path, in_file_ending, out_path, loci_file, org_file, num_threads):
    """Makes a list of commands to send to exonerator() and parallizes the run. """
    orglist = filelines_to_list(org_file)
    loclist = filelines_to_list(loci_file)
    query_dict = SeqIO.to_dict(SeqIO.parse(query_file, "fasta"))
    command_args = []
    for loc in loclist:
        in_file = f"{in_path}/{loc}{in_file_ending}"
        command_args.append(
            [query_file, in_file, in_file_ending, out_path, loc])
    p = Pool(num_threads)
    p.map(exonerator, command_args)
    return


def exonerator(command_args):
    """Run exonerate. This is a function to allow parallelization with multiprocessing.Pool()"""
    query_file, in_file, in_file_ending, out_path, loc = command_args
    exonerate_output = subprocess.getoutput(
        f"exonerate --model protein2genome -q {query_file} -t {in_file} -Q protein -T dna --showvulgar F --showalignment F --verbose 0 --fsmmemory 20G --ryo '>%ti %qi\n%tas\n'")
    with open(f"{out_path}/{loc}{in_file_ending}", 'w') as out_handle:
        out_handle.write(exonerate_output)
    return


def exonerateout_cleaner(exonerateout_path, file_ending, loci_file, exonerateclean_path):
    """This parses the exonerate output from exonerator so that the query hit and the target hit are the same loci."""
    loclist = filelines_to_list(loci_file)
    fasta_paths_to_dedupe = []
    for loc in loclist:
        to_write = []
        in_file = f"{exonerateout_path}/{loc}{file_ending}"
        records = SeqIO.parse(in_file, "fasta")
        for record in records:
            qloci = record.description.split()[1]
            tloci = record.id.split('-')[1]
            if qloci == tloci:
                to_write.append(record)
        out_path = f"{exonerateclean_path}/{loc}{file_ending}"
        with open(out_path, 'w') as out_handle:
            SeqIO.write(to_write, out_handle, 'fasta')
            fasta_paths_to_dedupe.append(os.path.abspath(out_path))
    return fasta_paths_to_dedupe


def fasta_deduper(fasta_file):
    """Rewrites over fasta file with only one seq per name. Keeps the longest sequence. Removes seq description"""
    u_records = []
    u_record_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        record.description=''
        if record.id in u_record_dict:
            if u_record_dict[record.id][0] < len(record.seq):
                u_record_dict[record.id] = [len(record.seq), record]
            else:
                pass
        else:
            u_record_dict[record.id] = [len(record.seq), record]
    u_records = [u_record_dict[record][1] for record in u_record_dict.keys()]
    with open(fasta_file, 'w') as out_handle:
        SeqIO.write(u_records, out_handle, "fasta")
    return

def fasta_subseter(subset_list, in_fasta_path, out_fasta_path):
    """Given a list of fasta headers, it will extract the seqs in that list from in_fasta_path
    and write them to out_fasta_path
    """
    input_dict = SeqIO.to_dict(SeqIO.parse(in_fasta_path, "fasta"))
    out_records = []
    for name in subset_list:
        out_records.append(input_dict[name])
    with open(out_fasta_path, "w") as out_handle:
        SeqIO.write(out_records, out_handle, 'fasta')
    return

def maffter(fasta_file_path,out_file_path,num_threads):
    """ Runs and auto mafft alignemnt on a fasta file"""
    mafft_output = subprocess.getoutput(
        f"mafft --thread {num_threads} --quiet --auto {fasta_file_path} ")
    with open(out_file_path, "w") as out_handle:
        out_handle.write(mafft_output)
    return


def txtm_fishing_pipe(project_path, num_threads, org_list_path, loci_list_path, blast_query_path, exonerate_query_path, e_vals, min_num_seqs):
    """Run pipeline in Python! Shell scripts are #OldSckool.
    This will iterate the blast searches over the e_vals.
    If a gene has no hits for a given evalue, it is lowered and searched again.
    The results from all of the blast searches are concatenated and fed into exonerate.
    Exonerate picks out the exons that we targeted.
    """
    dirs_to_make = out_dir_maker(project_path)
    org_list = filelines_to_list(org_list_path)
    blast_databases, blast_results, mafft_alignments, blast_by_gene, exonerate_out, exonerate_clean = dirs_to_make
    to_blast_list = filelines_to_list(loci_list_path)
#### Start the pipeline:
    for e_val in e_vals:
        blast_run(e_val, num_threads, blast_databases, blast_query_path, blast_results)
        seq_fetcher(blast_results, '.blast.txt', blast_databases, '.Trinity.annotated.nt', org_list_path)
        check_list = cat_by_gene(org_list_path, to_blast_list, blast_results, '.fa', blast_by_gene, min_num_seqs)
        if len(check_list) < 20:
            print(f"The folowing loci had no hits for {e_val}:{' '.join(check_list)}")
        elif len(check_list) == 0:
            print("All loci had hits at {e_val}")
        else:
            print(f"There were {len(check_list)} genes with no hits for {e_val}")
        to_blast_list = check_list
        new_query_path=f'{project_path}/SpurGenes_aa_{e_val}.fasta'
        fasta_subseter(to_blast_list,blast_query_path,new_query_path)
        blast_query_path=new_query_path
    exoneratetor_bygene(exonerate_query_path, blast_by_gene, '.FNA', exonerate_out, loci_list_path, org_list_path, num_threads)
    fasta_paths_to_dedupe = exonerateout_cleaner(exonerate_out, '.FNA', loci_list_path, exonerate_clean)
    for fasta_file in fasta_paths_to_dedupe:
        fasta_deduper(fasta_file)
        maffter(fasta_file,os.path.join(mafft_alignments,os.path.basename(fasta_file)),num_threads)
    return


project_path = "/Users/josec/Desktop/Crinoid_capture/Feb_hybTxCrinoid2"
num_threads = 7
min_num_seqs = 6
e_vals = ['1e-100','1e-90','1e-80','1e-70','1e-60','1e-50','1e-40','1e-30']
blast_query_path = f'{project_path}/SpurGenes_aa.fasta'
org_list_path = f'{project_path}/TranscriptomeList.txt'
loci_list_path = f'{project_path}/Trim_Loc_List.txt'
exonerate_query_path = f'{project_path}/SpurGenes_TrimtoTargets_AA.fa'

# Run the pipeline
txtm_fishing_pipe(project_path, num_threads, org_list_path, loci_list_path, blast_query_path, exonerate_query_path, e_vals, min_num_seqs)
