from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import json


class PipelineStepGenerateTaxidLocator(PipelineStep):
    def run(self):
        input_fasta = ""
        taxid_field = ""
        hit_type = ""
        output_fasta = ""
        output_json = ""
        TEMP_DIR = ""

        taxid_field_num = PipelineStepGenerateTaxidLocator.get_taxid_field_num(
            taxid_field, input_fasta)
        # Put every 2-line fasta record on a single line with delimiter
        # ":lineseparator:":
        cmd = "awk 'NR % 2 == 1 { o=$0 ; next } { print o \":lineseparator:\" $0 }' " + input_fasta
        # Sort the records based on the field containing the taxids
        cmd += " | sort -T %s --key %s --field-separator ':' --numeric-sort" % (
            TEMP_DIR, taxid_field_num)
        # Split every record back over 2 lines
        cmd += " | sed 's/:lineseparator:/\\n/g' > %s" % output_fasta
        command.execute(cmd)

        # Make JSON file giving the byte range of the file corresponding to each
        # taxid
        taxon_sequence_locations = []
        f = open(output_fasta, 'rb')
        sequence_name = f.readline()
        sequence_data = f.readline()

        taxid = PipelineStepGenerateTaxidLocator.get_taxid(
            sequence_name, taxid_field)
        first_byte = 0
        end_byte = first_byte + len(sequence_name) + len(sequence_data)
        while len(sequence_name) > 0 and len(sequence_data) > 0:
            sequence_name = f.readline()
            sequence_data = f.readline()
            new_taxid = PipelineStepGenerateTaxidLocator.get_taxid(
                sequence_name, taxid_field)
            if new_taxid != taxid:
                # Note on boundary condition: when end of file is reached, then
                # sequence_name == '' => new_taxid == 'none' => new_taxid != taxid
                # so last record will be written to output correctly.
                taxon_sequence_locations.append({
                    'taxid': int(taxid),
                    'first_byte': first_byte,
                    'last_byte': end_byte - 1,
                    'hit_type': hit_type
                })
                taxid = new_taxid
                first_byte = end_byte
                end_byte = first_byte + len(sequence_name) + len(sequence_data)
            else:
                end_byte += len(sequence_name) + len(sequence_data)
        f.close()

        with open(output_json, 'wb') as f:
            json.dump(taxon_sequence_locations, f)

    @staticmethod
    def get_taxid_field_num(taxid_field, input_fasta):
        with open(input_fasta) as f:
            sequence_name = f.readline()
        return sequence_name.replace('>',
                                     ':').split(":").index(taxid_field) + 1

    @staticmethod
    def get_taxid(sequence_name, taxid_field):
        parts = sequence_name.replace('>', ':').split(":%s:" % taxid_field)
        if len(parts) <= 1:
            # Sequence_name empty or taxid_field not found
            return 'none'
        taxid = parts[1].split(":")[0]
        # Example sequence_name: ">nr:-100:nt:684552:NR::NT:LT629734.1:HWI-ST640
        # :828:H917FADXX:2:1101:1424:15119/1"
        return taxid
