#! /usr/bin/env python3
import argparse
import os
import time
import subprocess
import multiprocessing


def parse_args():
    parser = argparse.ArgumentParser(prog='GSEA', description='gene set enrichment analysis.',
        epilog="""
                Requirement:
                python3
                HISAT2
                stringtie2
                cufflinks
                samtools
                R
                DESeq2(R packages)
               """)
    parser.add_argument('data_dir', help='data set.')
    parser.add_argument('-ref_fasta', help='reference genome sequence.', required=True)
    parser.add_argument('-ref_gtf', help='reference genome gtf/gff3 file.', required=True)
    parser.add_argument('-hisat2_index', help='hisat2 genome index file, if this\
        is not provided, prog will make index file auto.', default=None)
    parser.add_argument('-p', help='number of threads used.', default=int(os.cpu_count() * 0.8), type=int)
    args = parser.parse_args()
    return args.data_dir, args.ref_fasta, args.ref_gtf, args.hisat2_index, args.p


def struc_sample_info(data_dir):
    dt_out = {}
    for ele in os.listdir(data_dir):
        dt_out[ele] = {}
        for ele2 in os.listdir(os.path.join(data_dir, ele)):
            dt_out[ele][ele2] = []
            for ele3 in os.listdir(os.path.join(data_dir, ele, ele2)):
                dt_out[ele][ele2].append(os.path.join(data_dir, ele, ele2, ele3))
    return dt_out


def hisat2_indexing(reference_fasta, reference_gtf, threads, working_dir):
    cmdname_1 = 'hisat2_extract_splice_sites.py'
    ss_file = os.path.join(working_dir, 'ss_file')
    subprocess.run([cmdname_1, reference_gtf], stdout=open(ss_file, 'w'), stderr=open(os.devnull, 'w'))
    cmdname_2 = 'hisat2_extract_exons.py'
    exon_file = os.path.join(working_dir, 'exon_file')
    subprocess.run([cmdname_2, reference_gtf], stdout=open(exon_file, 'w'), stderr=open(os.devnull, 'w'))
    cmdname_3 = 'hisat2-build'
    index_file = os.path.join(working_dir, 'hisat2_index_file')
    subprocess.run([cmdname_3, '--ss', ss_file, '--exon', exon_file, '-p', str(threads), reference_fasta, index_file],
        stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w'))
    return index_file


def assign_thread(thread, works):
    if (thread / works) < 1:
        return 0
    else:
        dt_out = [0] * works
        i = 0
        for ele in range(thread):
            dt_out[i] += 1
            i += 1
            i = i % works
        return dt_out


def map_reads_hisat2(hisat2_index_file, sample_info, thread, working_dir):
    """
        Using HISAT to map reads to reference genome.
        1) index reference genome, if it is not
        provided by HISAT website or you just want to index by yourself, by the way, this step
        is quired RAM consuming.
            first, extract splice-site and exon information form gtf file. you need convert gff
            into gtf file if gtf file not available.
                extract_splice_site.py gtf_file >ss_file
                extract_exons.py gtf_file >exon_file
            second, build HISAT2 index.
                hisat2-build -p threads_num --ss ss_file --exon exon_file fasta_file_indexing output_file_name

        2)map reads to reference genome.
            hisat2 -p treads_num --dta -x reference_genome_index_prefix -1 fastq1 -2 fastq2 -S out.sam
    """

    def process_worker(data_point, hisat2_index_file, thread):
        reads_1 = data_point[0]
        reads_2 = data_point[1]
        sam_file = data_point[2]
        cmdname = 'hisat2'
        subprocess.run([cmdname, '-p', str(thread), '--dta', '-x', hisat2_index_file, '-1', reads_1, '-2', reads_2, '-S', sam_file],
            capture_output=True)
        return 0

    all_args = []
    works = 0
    for ele in sample_info:
        for ele2 in sample_info[ele]:
            sample_info[ele][ele2].append(os.path.join(working_dir, ele2 + '.sam'))
            works += 1
    # DEBUG: print(sample_info)
    thread_per_process = assign_thread(thread, works)
    if thread_per_process == 0:
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                all_args.append([sample_info[ele][ele2], hisat2_index_file, 1])
        i = 0
        while True:
            process = []
            for ele in range(thread):
                process.append(multiprocessing.Process(target=process_worker, args=all_args[i]))
                i += 1
                if i > works:
                    break
            for p in process:
                p.start()
            for p in process:
                p.join()
            if i > works:
                break
    else:
        i = 0
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                all_args.append([sample_info[ele][ele2], hisat2_index_file, thread_per_process[i]])
                i += 1
        process = []
        for ele in all_args:
            process.append(multiprocessing.Process(target=process_worker, args=ele))
        for p in process:
            p.start()
        for p in process:
            p.join()
    return 0


def sort_sam_samtools(sample_info, thread, working_dir):
    """
    samtools sort -@ thread -o bam_sorted_res sam_file
    """
    def process_worker(sam_file, bam_result, thread):
        cmdname = 'samtools'
        subprocess.run([cmdname, 'sort', '-@', str(thread), '-o', bam_result, sam_file],
            capture_output=True)
        return 0

    works = 0
    process_args = []
    for ele in sample_info:
        for ele2 in sample_info[ele]:
            sample_info[ele][ele2].append(os.path.join(working_dir, ele2 + '_sorted.bam'))
            works += 1
    thread_per_process = assign_thread(thread, works)
    if thread_per_process == 0:
        i = 0
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                process_args.append([sample_info[ele][ele2][-2], sample_info[ele][ele2][-1], 1])
        while True:
            process = []
            for ele in range(thread):
                process.append(multiprocessing.Process(target=process_worker, args=process_args[i]))
                i += 1
                if i > works:
                    break
            for p in process:
                p.start()
            for p in process:
                p.join()
            if i > works:
                break
    else:
        i = 0
        process = []
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                process_args.append([sample_info[ele][ele2][-2], sample_info[ele][ele2][-1], thread_per_process[i]])
        for ele in process_args:
            process.append(multiprocessing.Process(target=process_worker, args=ele))
        for p in process:
            p.start()
        for ele in process:
            p.join()
    return 0


def assemble_transcripts(sample_info, thread, reference_gtf, working_dir):
    """
        using stringtie to assemble transcripts.
        1)assembly
            stringtie -p threads_num -G reference_gtf/gff_file -o out_gtf_file -l name_prefix_for_output_transcripts bam_file
        2)merger transcripts
            stringtie --merge -p threads_num  -G reference_gtf/gff_file -o result_gtf merger_file_list
    """
    def process_worker(reference_gtf, bam_file, name_prefix, out_gtf_file, thread):
        cmdname = 'stringtie'
        subprocess.run([cmdname, '-p', str(thread), '-G', reference_gtf, '-o', out_gtf_file, '-l', name_prefix, bam_file],
            capture_output=True)
        return 0

    def merge_transcripts(sample_info, reference_gtf, thread, working_dir):
        merge_list = os.path.join(working_dir, 'merge_list')
        merged_transcripts = os.path.join(working_dir, 'merged_transcripts.gtf')
        tmp = []
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                tmp.append(sample_info[ele][ele2][-1])
        f_merge_list = open(merge_list, 'w')
        for ele in tmp:
            print(ele, file=f_merge_list)
        f_merge_list.close()
        cmdname = 'stringtie'
        subprocess.run([cmdname, '--merge', '-p', str(thread), '-G', reference_gtf, '-o', merged_transcripts, merge_list],
            capture_output=True)
        return merged_transcripts

    works = 0
    for ele in sample_info:
        for ele2 in sample_info[ele]:
            sample_info[ele][ele2].append(os.path.join(working_dir, ele2 + '.gtf'))
            works += 1
    thread_per_process = assign_thread(thread, works)
    process_args = []
    if thread_per_process == 0:
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                process_args.append([reference_gtf, sample_info[ele][ele2][-2], ele2, sample_info[ele][ele2][-1], 1])
        i = 0
        while True:
            process = []
            for ele in range(thread):
                process.append(multiprocessing.Process(target=process_worker, args=process_args[i]))
                i += 1
                if i > works:
                    break
            for p in process:
                p.start()
            for p in process:
                p.join()
            if i > works:
                break
    else:
        process = []
        i = 0
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                process_args.append([reference_gtf, sample_info[ele][ele2][-2], ele2, sample_info[ele][ele2][-1], thread_per_process[i]])
                i += 1
        for ele in process_args:
            process.append(multiprocessing.Process(target=process_worker, args=ele))
        for p in process:
            p.start()
        for p in process:
            p.join()
    merged_gtf = merge_transcripts(sample_info, reference_gtf, thread, working_dir)
    return merged_gtf


def estimate_abundance(sample_info, merged_gtf, thread, working_dir):
    """
    1)Estimate transcript abundances
        stringtie  -p -G  thread merged_gtf -o abundence_gtf sorted_bam
    """
    def process_worker(merged_gtf, abundance_gtf, bam_file, thread):
        cmdname = 'stringtie'
        subprocess.run([cmdname, '-p', str(thread), '-e', '-G', merged_gtf, '-o', abundance_gtf, bam_file],
            capture_output=True)
        return 0

    works = 0
    for ele in sample_info:
        for ele2 in sample_info[ele]:
            sample_info[ele][ele2].append(os.path.join(working_dir, ele2 + '_abundance.gtf'))
            works += 1
    print(sample_info)
    thread_per_process = assign_thread(thread, works)
    process_args = []
    if thread_per_process == 0:
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                process_args.append([merged_gtf, sample_info[ele][ele2][-1], sample_info[ele][ele2][3], 1])
        i = 0
        while True:
            process = []
            for ele in range(thread):
                process.append(multiprocessing.Process(target=process_worker, args=process_args[i]))
                i += 1
                if i > works:
                    break
            for p in process:
                p.start()
            for p in process:
                p.join()
            if i > works:
                break
    else:
        i = 0
        process = []
        for ele in sample_info:
            for ele2 in sample_info[ele]:
                process_args.append([merged_gtf, sample_info[ele][ele2][-1], sample_info[ele][ele2][3], thread_per_process[i]])
                i += 1
        for ele in process_args:
            process.append(multiprocessing.Process(target=process_worker, args=ele))
        for p in process:
            p.start()
        for p in process:
            p.join()
    return 0


def count_reads(sample_info, working_dir, reads_len=76):
    """
    using prepDE.py count reads from gtf file. read = (coverage * length)/reads_length.
    prepDE.py -i gtf.list -l reads_len
    """
    reads_abundance_list = os.path.join(working_dir, 'gtf.list')
    abundances_list = open(reads_abundance_list, 'w')
    gene_reads_count = os.path.join(working_dir, 'gene_count_matrix.csv')
    transcript_reads_count = os.path.join(working_dir, 'transcript_count_matrix.csv')
    for ele in sample_info:
        for ele2 in sample_info[ele]:
            print(ele2 + '\t' + sample_info[ele][ele2][5], file=abundances_list)
    abundances_list.close()
    subprocess.run(['prepDE.py', '-i', reads_abundance_list, '-l', str(reads_len), '-g', gene_reads_count, '-t', transcript_reads_count],
        capture_output=True)
    return gene_reads_count, transcript_reads_count


def reform_reads_count_file(file_name, sample_order):
    a_line = [line.rstrip().split(',') for line in open(file_name)]
    new_line = []
    head = a_line[0]
    index_order = []
    for ele in sample_order:
        index_order.append(head.index(ele))
    for line in a_line:
        tmp = []
        tmp.append(line[0])
        for index in index_order:
            tmp.append(line[index])
        new_line.append(tmp)
    f_out = open(file_name, 'w')
    for line in new_line:
        print(','.join(line), file=f_out)
    f_out.close()
    return 0


def call_different_transcripts(sample_info, working_dir, gene_reads_count, transcript_reads_count):
    """
    using DESeq2 R packages to call differently express genes/transcripts.
    """
    meta_file = os.path.join(working_dir, 'meta.csv')
    meta_f = open(meta_file, 'w')
    print('sample,group', file=meta_f)
    sample_order = []
    for ele in sample_info:
        for ele2 in sample_info[ele]:
            print(ele2 + ',' + ele, file=meta_f)
            sample_order.append(ele2)
    meta_f.close()
    reform_reads_count_file(gene_reads_count, sample_order)
    reform_reads_count_file(transcript_reads_count, sample_order)
    gene_compare_res = os.path.join(working_dir, 'gene_diff.csv')
    transcript_compare_res = os.path.join(working_dir, 'transcript_diff.csv')

    subprocess.run(['./call_differ_transcripts.R', gene_reads_count, transcript_reads_count, meta_file,
        gene_compare_res, transcript_compare_res], capture_output=False)
    return gene_compare_res, transcript_compare_res


def function_interpration():
    return 0


def main():
    data_dir, reference_fasta, reference_gtf, hisat2_index_file, thread = parse_args()
    working_dir = 'working_dir_' + time.strftime('%Y%m%d%H%M%S')
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    if reference_gtf.split('.')[-1] == 'gff':
        gft_file = os.path.join(working_dir, 'reference.gtf')
        cmdname = 'gffread'
        subprocess.run([cmdname, reference_gtf, '-T', '-o', gft_file], capture_output=True)
        reference_gtf = gft_file
    sample_info = struc_sample_info(data_dir)
    # DEBUG: print(sample_info)
    if not hisat2_index_file:
        hisat2_index_file = hisat2_indexing(reference_fasta, reference_gtf, thread, working_dir)
    map_reads_hisat2(hisat2_index_file, sample_info, thread, working_dir)
    sort_sam_samtools(sample_info, thread, working_dir)
    merged_gtf = assemble_transcripts(sample_info, thread, reference_gtf, working_dir)
    print("transcripts assembled.")
    estimate_abundance(sample_info, merged_gtf, thread, working_dir)
    print('transcripts abundance reestimated.')
    gene_reads_count, transcript_reads_count = count_reads(sample_info, working_dir, 76)
    gene_compare_res, transcript_compare_res = call_different_transcripts(sample_info, working_dir,
        gene_reads_count, transcript_reads_count)
    return 0


if __name__ == '__main__':
    main()
