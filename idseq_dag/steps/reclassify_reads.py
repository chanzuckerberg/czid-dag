import json
from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count


MIN_READS_PER_GENUS = 100

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
        (blast_m8, refined_m8, refined_hit_summary, refined_counts, genus_to_accession_list) =
        self.output_files_local()
        (read_dict, accession_dict, selected_genera) = self.summarize_hits(hit_summary)
        output_basedir = os.path.join(self.output_dir_local, "reclassify")
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

        # Start a non-blocking function for download ref sequences
        genus_references = {}
        download_thread = threading.Thread(
                target = self.download_ref_sequences,
                args = [selected_genera, output_basedir, loc_db, db_s3_path, genus_references])
        download_thread.start()
        genus_fasta_files = self.group_reads_by_genus(input_fasta, read_dict, selected_genera, output_basedir)
        genus_assembled = self.assemble_all(genus_fasta_files)
        download_thread.join()

        genus_blast_m8 = self.blast_all(genus_assembled_contig, genus_references, db_type)
        (consolidated_dict, read2blastm8)  = self.consolidate_all_results(genus_assembled_contig,
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
        evalue_type = 'log10' if db_type == 'nr' else 'raw'
        m8.generate_taxon_count_json_from_m8(refined_m8, refined_hit_summary,
                                             evalue_type, db_type.upper(),
                                             lineage_db, deuterostome_db, refined_counts)
        # additional files to upload
        # assembled contigs. bowtie2.sam, blast.m8, top_blast.m8 for read to contig mapping

    @staticmethod
    def consolidate_all_results(genus_assembled_contig, genus_blast_m8, read_dict, accession_dict):
        consolidated_dict = read_dict
        read2blastm8 = {}
        for genus_taxid, assembled in genus_assembled_contig.items():
            read2contig = assembled[2]
            blast_top_m8 = genus_blast_m8.get(genus_taxid, (None, None))[1]
            if read2contig and blast_top_m8:
                contig2accession = {}
                for contig_id, accession_id, _percent_id, _alignment_length, e_value, _bitscore, line in iterate_m8(blast_top_m8):
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
    def blast_all(genus_assembled_contig, genus_references, db_type):
        genus_blast_m8 = {}
        for genus_taxid, output in genus_assembled_contig.items():
            if output[0]:
                (contig, scaffold, read2contig) = output
                genus_dir = os.path.dirname(contig)
                # Make blast index
                accession_dir = genus_references[genus_taxid]
                reference_fasta = os.path.join(genus_dir, "refseq.fasta")
                blast_index_path = os.path.join(genus_dir, "blastindex")
                command.execute(f"cat {accession_dir}/* > {reference_fasta}")
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
                genus_blast_m8[genus_taxid] = (output_m8, top_entry_m8)
        return genus_blast_m8

    @staticmethod get_top_m8(orig_m8, top_entry_m8):
        with open(top_entry_m8, 'w') as top_m8f:
            top_line = None
            top_bitscore = 0
            current_read_id = None
            for read_id, _accession_id, _percent_id, _alignment_length, e_value, bitscore, line in iterate_m8(orig_m8):
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


    @staticmethod
    def download_ref_sequences(selected_genera, output_basedir, loc_db, db_s3_path, genus_references):
        threads = []
        error_flags = {}
        semaphore = threading.Semaphore(64)
        mutex = threading.RLock()

        bucket, key = db_s3_path[5:].split("/", 1)
        loc_dict = shelve.open(loc_db.replace('.db', ''), 'r')
        for genus_taxid, accession_list in selected_genera.items():
            accession_dir = os.path.join(output_basedir, genus_taxid, 'accessions')
            command.execute(f"mkdir -p {accession_dir}")
            genus_references[genus_taxid] = accession_dir
            for accession in accession_list:
                accession_out_file=os.path.join(accession_dir, accession)
                semaphore.acquire()
                t = threading.Thread(
                    target=self.fetch_sequence_for_thread,
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


    @staticmethod
    def assemble_all(genus_fasta_files):
        genus_assembled = {}
        for genus_taxid, fasta_file in genus_fasta_files:
            genus_dir = os.path.dirname(fasta_file)
            assembled_dir = os.path.join(genus_dir, 'spades')
            output = (None, None, None) # (contigs.fasta path, scaffolds.fasta path, readid -> contig mapping)
            command.execute(f"mkdir -p {assembled_dir}")

            try:
                command.execute(f"spades.py -s {fasta_file} -o {assembled_dir} -m 60 -t 32 --only-assembler")
                assembled_contig = os.path.join(assembled_dir, 'contigs.fasta')
                assembled_scaffold = os.path.join(assembled_dir, 'scaffolds.fasta')

                for i, assembled in enumerate([assembled_contig, assembled_scaffold]):
                    if os.path.exists(assembled) and os.path.getsize(assembled) > MIN_CONTIG_SIZE:
                        # move the file up to the main dir
                        final_path = assembled.replace('spades/', '')
                        command.execute(f"mv {assembled} {final_path}")
                        output[i] = final_path
                # build the bowtie index based on the contigs
                if output[0]:
                    output[2] = self.generate_read_to_contig_mapping(output[0], fasta_file)
            except:
                pass
            genus_assembled[genus_taxid] = output
            command.execute(f"rm -rf {assembled_dir}")

        return genus_assembled

    @staticmethod
    def generate_read_to_contig_mapping(assembled_contig, fasta_file):
        genus_dir = os.path.dirname(fasta_file)
        # build bowtie index based on assembled_contig
        bowtie_index_path = os.path.join(genus_dir, 'bowtie-contig')
        bowtie_sam = os.path.join(genus_dir, "read-output.sam")
        command.execute(f"bowtie2-build {assembled_contig} {bowtie_index_path}")
        command.execute(f"bowtie2 -x {bowtie_index_path} -f -U {fasta_file} --very-sensitive -p 32 > {bowtie_sam}")
        read2contig = {}
        with open(bowtie_sam, "r", encoding='utf-8') as samf:
            for line in samf:
                if line[0] == '@':
                    continue
                fields = line.split("\t")
                read = fields[0]
                contig = fields[2]
                if contig != '*':
                    read2contig[read] = contig
        return read2contig


    @staticmethod
    def group_reads_by_genus(input_fasta, read_dict, selected_genera, output_basedir):
        ''' returns genus_id => fasta file which includes sequences from the genus from alignment '''
        genus_fasta_data = defaultdict("")
        genera_list = set(selected_genera.keys())
        with open(input_fasta, 'r', encoding='utf-8') as input_fasta_f:
                sequence_name = input_fasta_f.readline()
                sequence_data = input_fasta_f.readline()
                while sequence_name and sequence_data:
                    read_id = sequence_name.rstrip().lstrip('>')
                    genus_taxid = read_dict[read_id][5]
                    if genus_taxid in genera_list:
                        genus_fasta_data[genus_taxid] += sequence_name + sequence_data
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
        read_dict = {}
        accession_dict = {}
        genus_read_counts = defaultdict(0)
        genus_species = defaultdict(set)
        genus_accessions = defaultdict(set)
        total_reads = 0
        with open(hit_summary_file, 'r') as hsf:
            for line in hsf:
                read = line.rstrip().split("\t")
                accession_id = read[3]
                species_taxid = read[4]
                genus_taxid = read[5]
                read_dict[read[0]] = read
                total_reads += 1
                if accession_id == 'None' or int(genus_taxid) < 0:
                    continue
                accession_dict[accession_id] = (species_taxid, genus_taxid)
                genus_read_counts[genus_taxid] += 1
                genus_species[genus_taxid].add(species_taxid)
                genus_accessions[genus_taxid].add(accession_id)
        selected_genera = {}
        for genus_taxid, reads in genus_read_counts.items():
            if reads >= MIN_READS_PER_GENUS and len(genus_species[genus_taxid]) > 1:
                selected_genera[genus_taxid] = list(genus_accessions[genus_taxid])

        return (read_dict, accession_dict, selected_genera)


    def count_reads(self):
        pass

