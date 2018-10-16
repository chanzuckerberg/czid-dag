import json
import os
import shelve
import threading
import time
import traceback
from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count
import idseq_dag.util.s3 as s3
import idseq_dag.util.m8 as m8

from idseq_dag.steps.run_assembly import PipelineStepRunAssembly


MIN_READS_PER_GENUS = 100
MIN_CONTIG_SIZE = 50

class PipelineStepReclassifyReads(PipelineStep):
    '''
        Reclassify reads after alignment
    '''
    def run(self):
        '''
            1. organize the reads into different genus and get the relevant accessions from first pass
                1.1 decide based on # reads in a genus if we want to reclassify. min(100, 1% nonhost)  reads
                1.2 decide based on # species under a genus. > 1
            2a. download the actual reference sequences based on the accession data.
                build blastindex per genus
            2b. assemble  per genus at the same time from 1
            3. balst the contigs from 2b to 2a. get the one with the best bitscore
            4. build the bowtie index from 2b
            5. align the short reads from 1 to 4.
            6. for reads that align to the contigs, we assign the m8 entry
               based on where the best m8 entry the contig blast to and the species accordingly
            7. for reads that do not align to the contigs, we blast it and try to re-assign it

        '''
        (align_m8, deduped_m8, hit_summary, _orig_counts) = self.input_files_local[0]
        input_fasta = self.input_files_local[1][-1]
        (blast_m8, refined_m8, refined_hit_summary, refined_counts, genus_to_accession_list) = self.output_files_local()
        loc_db = s3.fetch_from_s3(
            self.additional_files["loc_db"],
            self.ref_dir_local,
            allow_s3mi=True)
        db_s3_path = self.additional_attributes["db"]
        db_type = self.additional_attributes["db_type"]
        lineage_db = s3.fetch_from_s3(
            self.additional_files["lineage_db"],
            self.ref_dir_local,
            allow_s3mi=True)

        (read_dict, accession_dict, selected_genera) = self.summarize_hits(hit_summary)
        output_basedir = os.path.join(self.output_dir_local, f"assembled_{db_type}")

        # Start a non-blocking function for download ref sequences
        genus_references = {} # genus -> accession list dir
        download_thread = threading.Thread(
            target = self.download_ref_sequences,
            args = [selected_genera, output_basedir, loc_db, db_s3_path, genus_references])
        download_thread.start()


        genus_fasta_files = self.group_reads_by_genus(input_fasta, read_dict,
                                                      selected_genera, output_basedir)
        genus_assembled = self.assemble_all(genus_fasta_files)
        download_thread.join()

        genus_blast_m8 = self.blast_all(genus_assembled, genus_references, db_type)
        (consolidated_dict, read2blastm8) = self.consolidate_all_results(genus_assembled,
                                                                          genus_blast_m8,
                                                                          read_dict,
                                                                          accession_dict)

        # output the results
        self.output_genus_accession_file(selected_genera, genus_to_accession_list)
        self.output_blast_m8_file(genus_blast_m8, blast_m8)
        self.generate_m8_and_hit_summary(consolidated_dict, read2blastm8,
                                         hit_summary, deduped_m8,
                                         refined_hit_summary, refined_m8)
        deuterostome_db = None
        evalue_type = 'raw'
        if self.additional_files.get("deuterostome_db"):
            deuterostome_db = s3.fetch_from_s3(self.additional_files["deuterostome_db"],
                                               self.ref_dir_local, allow_s3mi=True)
        m8.generate_taxon_count_json_from_m8(refined_m8, refined_hit_summary,
                                             evalue_type, db_type.upper(),
                                             lineage_db, deuterostome_db, refined_counts)
        # additional files to upload
        # assembled contigs. bowtie2.sam, blast.m8, top_blast.m8 for read to contig mapping
        for genus_taxid, assembled in genus_assembled.items():
            (contig_file, scaffold_file, read2contig, bowtie_sam) = assembled
            if contig_file:
                self.additional_files_to_upload.append(contig_file)
                self.additional_files_to_upload.append(scaffold_file)
                self.additional_files_to_upload.append(bowtie_sam)
        for genus_taxid, m8s in genus_blast_m8.items():
            (m8_raw, m8_top, refseq_file) = m8s
            self.additional_files_to_upload.append(m8_raw)
            self.additional_files_to_upload.append(m8_top)
            self.additional_files_to_upload.append(refseq_file)



    @staticmethod
    def output_genus_accession_file(selected_genera, genus_to_accession_list):
        with open(genus_to_accession_list, 'w') as gaf:
            for genus_taxid, accession_list in selected_genera.items():
                outstr = "\t".join([genus_taxid, ",".join(accession_list)]) + "\n"
                gaf.write(outstr)

    @staticmethod
    def output_blast_m8_file(genus_blast_m8, blast_m8):
        ''' output one file that include all the m8 content. append genus_taxid at front '''
        with open(blast_m8, 'w') as gbmf:
            for genus_taxid, m8s in genus_blast_m8.items():
                raw_m8 = m8s[0]
                for _contig_id, _accession_id, _percent_id, _alignment_length, e_value, _bitscore, line in m8.iterate_m8(raw_m8):
                    # add the genus taxid to each line
                    gbmf.write(f"{genus_taxid}:{line}")

    @staticmethod
    def generate_m8_and_hit_summary(consolidated_dict, read2blastm8,
                                    hit_summary, deduped_m8,
                                    refined_hit_summary, refined_m8):
        ''' generate new m8 and hit_summary based on consolidated_dict and read2blastm8 '''
        # Generate new hit summary
        with open(refined_hit_summary, 'w') as rhsf:
            with open(hit_summary, 'r', encoding='utf-8') as hsf:
                for line in hsf:
                    read_id = line.rstrip().split("\t")[0]
                    read = consolidated_dict[read_id]
                    output_str = "\t".join(read)
                    rhsf.write(output_str + "\n")
        # Generate new M8
        with open(refined_m8, 'w') as rmf:
            with open(deduped_m8, 'r', encoding='utf-8') as mf:
                for line in mf:
                    read_id = line.rstrip().split("\t")[0]
                    m8_line = read2blastm8.get(read_id)
                    if m8_line:
                        m8_fields = m8_line.split("\t")
                        m8_fields[0] = read_id
                        rmf.write("\t".join(m8_fields))
                    else:
                        rmf.write(line)


    @staticmethod
    def consolidate_all_results(genus_assembled_contig, genus_blast_m8, read_dict, accession_dict):
        ''' based on the blast results, refine read_dict and note  in read2blastm8'''
        consolidated_dict = read_dict # replace read_dict with the assemebled results if suitable
        read2blastm8 = {}
        for genus_taxid, assembled in genus_assembled_contig.items():
            read2contig = assembled[2]
            blast_top_m8 = genus_blast_m8.get(genus_taxid, (None, None, None))[1]
            if read2contig and blast_top_m8:
                contig2accession = {}
                for contig_id, accession_id, _percent_id, _alignment_length, e_value, _bitscore, line in m8.iterate_m8(blast_top_m8):
                    contig2accession[contig_id] = (accession_id, line)
                for read_id, contig_id in read2contig.items():
                    (accession, m8_line) = contig2accession.get(contig_id, (None, None))
                    if accession:
                        (species_taxid, genus_taxid) = accession_dict[accession]
                        consolidated_dict[read_id] += [f"{genus_taxid}:{contig_id}", accession, species_taxid, genus_taxid]
                        consolidated_dict[read_id][2] = species_taxid
                    if m8_line:
                        read2blastm8[read_id] = m8_line
        return (consolidated_dict, read2blastm8)

    @staticmethod
    def reference_fasta(genus_dir):
        return os.path.join(genus_dir, "refseq.fasta")

    def blast_all(self, genus_assembled, genus_references, db_type):
        ''' blast the assembled contigs to downloaded accession list '''
        genus_blast_m8 = {} # genus => (raw_m8, top_entry_m8, reference_fasta)
        for genus_taxid, output in genus_assembled.items():
            if output[0]:
                (contig, scaffold, read2contig, bowtie_sam) = output
                genus_dir = os.path.dirname(contig)
                # Make blast index
                reference_fasta = self.reference_fasta(genus_dir)
                blast_index_path = os.path.join(genus_dir, "blastindex")
                if s3.check_s3_presence(self.s3_path(reference_fasta)):
                    reference_fasta = s3.fetch_from_s3(self.s3_path(reference_fasta), genus_dir)
                else:
                    accession_dir = genus_references[genus_taxid]
                    command.execute(f"find {accession_dir}/ -type f | xargs -n 32 -P 1 cat >> {reference_fasta}")
                blast_type = 'nucl'
                blast_command = 'blastn'
                if db_type == 'nr':
                    blast_type = 'prot'
                    blast_command = 'blastx'
                command.execute(f"makeblastdb -in {reference_fasta} -dbtype {blast_type} -out {blast_index_path}")
                # blast the contig to the blast index
                output_m8 = os.path.join(genus_dir, 'blast.m8')
                top_entry_m8 = os.path.join(genus_dir, 'blast_top.m8')
                command.execute(f"{blast_command} -query {contig} -db {blast_index_path} -out {output_m8} -outfmt 6 -num_alignments 5 -num_threads 32")
                # further processing of getting the top m8 entry for each contig.
                PipelineStepReclassifyReads.get_top_m8(output_m8, top_entry_m8)
                genus_blast_m8[genus_taxid] = (output_m8, top_entry_m8, reference_fasta)
        return genus_blast_m8

    @staticmethod
    def get_top_m8(orig_m8, top_entry_m8):
        ''' Get top m8 file entry for each read from orig_m8 and output to top_entry_m8 '''
        with open(top_entry_m8, 'w') as top_m8f:
            top_line = None
            top_bitscore = 0
            current_read_id = None
            for read_id, _accession_id, _percent_id, _alignment_length, e_value, bitscore, line in m8.iterate_m8(orig_m8):
                # Get the top entry of each read_id based on the bitscore
                if read_id != current_read_id:
                    # Different batch start
                    if current_read_id: # Not the first line
                        top_m8f.write(top_line)
                    current_read_id = read_id
                    top_line = line
                    top_bitscore = bitscore
                elif bitscore > top_bitscore:
                    top_bitscore = bitscore
                    top_line = line
            top_m8f.write(top_line)


    def download_ref_sequences(self, selected_genera, output_basedir,
                               loc_db, db_s3_path,
                               genus_references):
        ''' Download accessions specified in the selected_genera '''
        threads = []
        error_flags = {}
        semaphore = threading.Semaphore(64)
        mutex = threading.RLock()

        bucket, key = db_s3_path[5:].split("/", 1)
        loc_dict = shelve.open(loc_db.replace('.db', ''), 'r')
        for genus_taxid, accession_list in selected_genera.items():
            genus_dir = os.path.join(output_basedir, genus_taxid)
            reference_fasta = self.reference_fasta(genus_dir)
            if s3.check_s3_presence(self.s3_path(reference_fasta)):
                # ref seq data already exists. skip downloading
                continue
            accession_dir = os.path.join(genus_dir, 'accessions')
            command.execute(f"mkdir -p {accession_dir}")
            genus_references[genus_taxid] = accession_dir
            for accession in accession_list:
                accession_out_file = os.path.join(accession_dir, accession)
                semaphore.acquire()
                t = threading.Thread(
                    target=PipelineStepReclassifyReads.fetch_sequence_for_thread,
                    args=[
                        error_flags, accession, accession_out_file, loc_dict,
                        bucket, key, semaphore, mutex
                    ])
                t.start()
                threads.append(t)
        for t in threads:
            t.join()
        if error_flags:
            raise RuntimeError("Error in getting sequences by accession list.")

    @staticmethod
    def fetch_sequence_for_thread(error_flags, accession, accession_out_file,
                                  loc_dict, bucket, key,
                                  semaphore, mutex):
        ''' fetch sequence from S3 for the specific accession'''
        try:
            entry = loc_dict.get(accession)
            if entry:
                range_start, name_length, seq_len = entry
                range_end = range_start + name_length + seq_len - 1
                num_retries = 3
                for attempt in range(num_retries):
                    try:
                        s3.fetch_byterange(range_start, range_end, bucket, key, accession_out_file)
                        break
                    except Exception as e:
                        if attempt + 1 < num_retries:  # Exponential backoff
                            time.sleep(1.0 * (4**attempt))
                        else:
                            msg = f"All retries failed for getting sequence by accession ID {accession}: {e}"
                        raise RuntimeError(msg)

        except:
            with mutex:
                if not error_flags:
                    traceback.print_exc()
                error_flags["error"] = 1
        finally:
            semaphore.release()


    def assemble_all(self, genus_fasta_files):
        ''' assemble the individual fasta files by genus '''
        genus_assembled = {} # output genus => (contigs.fasta scaffolds.fasta, readid -> contig, bowtie_sam)
        for genus_taxid, fasta_file in genus_fasta_files.items():
            genus_dir = os.path.dirname(fasta_file)
            output = [None, None, None, None]
            assembled_contig = os.path.join(genus_dir, 'contigs.fasta')
            assembled_scaffold = os.path.join(genus_dir, 'scaffolds.fasta')
            read2config = {}
            contig_stats_json = os.path.join(genus_dir, 'contig_stats.json')
            bowtie_sam = os.path.join(genus_dir, 'read-contig.sam')

            if s3.check_s3_presence(self.s3_path(assembled_contig)) and \
                s3.check_s3_presence(self.s3_path(assembled_scaffold)) and \
                s3.check_s3_presence(self.s3_path(bowtie_sam)):
                # check if file already assembled before, if so, reuse.
                assembled_contig = s3.fetch_from_s3(self.s3_path(assembled_contig), genus_dir)
                assembled_scaffold = s3.fetch_from_s3(self.s3_path(assembled_scaffold), genus_dir)
                bowtie_sam = s3.fetch_from_s3(self.s3_path(bowtie_sam), genus_dir)
                _contig_stats = defaultdict(0)
                PipelineStepRunAssembly.generate_info_from_sam(bowtie_sam,
                                                               read2contig, _contig_stats)
            else:
                PipelineStepRunAssembly.assemble(fasta_file, assembled_contig, assembled_scaffold,
                                                 bowtie_sam, contig_stats_json, read2config)


            if len(read2config) > 0: # assemble success
                output = [assembled_contig, assembled_scaffold, read2contig, bowtie_sam]

            genus_assembled[genus_taxid] = output

        return genus_assembled

    @staticmethod
    def group_reads_by_genus(input_fasta, read_dict, selected_genera, output_basedir):
        ''' returns genus_id => fasta file which includes sequences from the genus from alignment '''
        genus_fasta_data = defaultdict(str)
        genera_list = set(selected_genera.keys())
        with open(input_fasta, 'r', encoding='utf-8') as input_fasta_f:
            sequence_name = input_fasta_f.readline()
            sequence_data = input_fasta_f.readline()
            while sequence_name and sequence_data:
                read_id = sequence_name.rstrip().lstrip('>')
                entry = read_dict.get(read_id)
                if entry and len(entry) > 5:
                    genus_taxid = entry[5]
                    if genus_taxid in genera_list:
                        genus_fasta_data[genus_taxid] += sequence_name + sequence_data
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()
        # output the fasta files and create work dir
        genus_fasta_files = {}
        for genus_taxid, fasta_content in genus_fasta_data.items():
            output_dir = os.path.join(output_basedir, genus_taxid)
            command.execute(f"mkdir -p {output_dir}")
            output_file = os.path.join(output_dir, 'input.fa')
            with open(output_file, 'w') as outf:
                outf.write(fasta_content)
            genus_fasta_files[genus_taxid] = output_file

        return genus_fasta_files

    @staticmethod
    def summarize_hits(hit_summary_file):
        ''' Parse the hit summary file and get the relevant into'''
        read_dict = {} # read_id => line
        accession_dict = {} # accession => (species, genus)
        genus_read_counts = defaultdict(int) # genus => read_counts
        genus_species = defaultdict(set) # genus => list of species
        genus_accessions = defaultdict(set) # genus => list of accessions
        total_reads = 0
        with open(hit_summary_file, 'r') as hsf:
            for line in hsf:
                read = line.rstrip().split("\t")
                accession_id, species_taxid, genus_taxid = read[3:6]
                read_dict[read[0]] = read
                total_reads += 1
                if accession_id == 'None' or accession_id == "" or int(genus_taxid) < 0:
                    continue
                accession_dict[accession_id] = (species_taxid, genus_taxid)
                genus_read_counts[genus_taxid] += 1
                genus_species[genus_taxid].add(species_taxid)
                genus_accessions[genus_taxid].add(accession_id)
        selected_genera = {} # genus => accession_list
        for genus_taxid, reads in genus_read_counts.items():
            if reads >= MIN_READS_PER_GENUS and len(genus_species[genus_taxid]) > 1:
                selected_genera[genus_taxid] = list(genus_accessions[genus_taxid])

        return (read_dict, accession_dict, selected_genera)


    def count_reads(self):
        pass

