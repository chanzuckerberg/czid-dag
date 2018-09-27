import json
import threading
import re
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log

from Bio import Entrez

class PipelineStepFetchTaxInfo(PipelineStep):
    '''
        fetch tax info based on a list
    '''
    def run(self):
        '''
            1. fetch the taxid -> wikipedia link mappin
            2. fetch wikipedia content
            3. store everything
        '''
        taxid_list = self.input_files_local[0][0]
        (taxid2wiki, taxid2desc) = self.output_files_local()
        Entrez.email = self.additional_attributes.get("entrez_email", "yunfang@chanzuckerberg.com")
        num_threads = self.additional_attributes.get("threads", 16)
        batch_size = self.additional_attributes.get("batch_size", 100)
        taxid2wikidict = {}
        threads = []
        semaphore = threading.Semaphore(num_threads)
        mutex = threading.RLock()
        batch = []
        with open(taxid_list, 'r') as taxf:
            for line in taxf:
                taxid = line.rstrip()
                if taxid == 'taxid':
                    continue # header
                batch.append(taxid)
                if len(batch) >= batch_size:
                    semaphore.acquire()
                    t = threading.Thread(
                        target=PipelineStepFetchTaxInfo.
                        get_taxid_mapping_for_batch,
                        args=[batch, taxid2wikidict, mutex, semaphore]
                        )
                    t.start()
                    threads.append(t)
                    batch = []
        if len(batch) > 0:
            semaphore.acquire()
            t = threading.Thread(
                target=PipelineStepFetchTaxInfo.
                get_taxid_mapping_for_batch,
                args=[batch, taxid2wikidict, mutex, semaphore]
                )
            t.start()
            threads.append(t)
        for t in threads:
            t.join()
        # output the data
        with open(taxid2wiki, 'w') as taxidoutf:
            for taxid, wikiurl in taxid2wikidict.items():
                taxidoutf.write(f"{taxid}\t{wikiurl}\n")
        # output dummay for actual wiki content for now
        command.execute(f"echo 1234 > {taxid2desc}")


    @staticmethod
    def get_taxid_mapping_for_batch(taxids, taxid2wikidict, mutex, semaphore):
        taxid_str = ",".join(taxids)
        log.write(f"fetching batch {taxid_str}")
        handle = Entrez.elink(dbfrom="taxonomy", id=taxid_str, cmd="llinks")
        record = Entrez.read(handle)
        handle.close()

        parsed = {}
        results = record[0]['IdUrlList']['IdUrlSet']
        for result in results:
            taxid = result['Id']
            wikiurl = ""
            for link in result['ObjUrl']:
                url = str(link['Url'])
                if re.search('wikipedia.org', url):
                    wikiurl = url
                    break
            parsed[taxid] = wikiurl
        semaphore.release()
        with mutex:
            taxid2wikidict.update(parsed)


    def count_reads(self):
        pass

