import os
import subprocess

from pprint import pprint,pformat

from ReadsUtils.ReadsUtilsClient import ReadsUtils
#from KBaseReport.KBaseReportClient import KBaseReport



def log(message):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(message)



class CutadaptRunner:

    CUTADAPT = 'cutadapt'

    def __init__(self, scratch):
        self.scratch = scratch
        self.clear_options()

    def clear_options(self):
        self.interleaved = False
        self.input_filename = None
        self.output_filename = None
        self.five_prime = None
        self.three_prime = None
        self.err_tolerance = None
        self.overlap = None
        self.min_read_length = None
        self.discard_untrimmed = None


    def set_input_file(self, filename):
        self.input_filename = filename

    def set_output_file(self, filename):
        self.output_filename = filename


    def set_three_prime_option(self, sequence, anchored):
        if anchored == 1:
            sequence = sequence + '$'
        self.three_prime = sequence

    def set_five_prime_option(self, sequence, anchored):
        if anchored == 1:
            sequence = '^' + sequence
        self.five_prime = sequence


    def set_error_tolerance(self, tolerance):
        self.err_tolerance = float(tolerance)

    def set_min_overlap(self, overlap):
        self.overlap = int(overlap)

    def set_min_read_length(self, min_read_length):
        self.min_read_length = int(min_read_length)

    def set_discard_untrimmed(self, discard_untrimmed):
        self.discard_untrimmed = int(discard_untrimmed)

    def set_interleaved(self, interleaved):
        self.interleaved = interleaved

    def _build_adapter_removal_options(self, cmd):
        if self.interleaved:
            cmd.append('--interleaved')

        if self.three_prime:
            cmd.append('-a')
            cmd.append(self.three_prime)
            if self.interleaved:
                cmd.append('-A')
                cmd.append(self.three_prime)

        if self.five_prime:
            cmd.append('-g')
            cmd.append(self.five_prime)
            if self.interleaved:
                cmd.append('-G')
                cmd.append(self.five_prime)

        if self.err_tolerance:
            cmd.append('--error-rate=' + str(self.err_tolerance))

        if self.overlap:
            cmd.append('--overlap=' + str(self.overlap))

        if self.min_read_length:
            cmd.append('--minimum-length=' + str(self.min_read_length))

        if int(self.discard_untrimmed) == 1:
            cmd.append('--discard-untrimmed')


    def run(self):
        cmd = [self.CUTADAPT]

        self._build_adapter_removal_options(cmd)

        if self.output_filename:
            cmd.append('-o')
            cmd.append(self.output_filename)

        if not self.input_filename:
            raise ValueError('Input filename must be set to run cutadapt')
        cmd.append(self.input_filename)

        log('running cutadapt:')
        log('    ' + ' '.join(cmd))

        p = subprocess.Popen(cmd,
                             cwd=self.scratch,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)

        report = ''
        while True:
            line = p.stdout.readline()
            if not line:
                break
            report += line
            log(line.replace('\n', ''))

        p.stdout.close()
        p.wait()
        report += "\n\n"
        log('process return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running cutadapt, return code: ' +
                             str(p.returncode) + '\n')
        return report


class CutadaptUtil:

    def __init__(self, config):
        pprint(config)
        self.scratch = config['scratch']
        self.callbackURL = config['SDK_CALLBACK_URL']


    def remove_adapters(self, params):
        print ("\nPARAMS:\n"+pformat(params)+"\n")  # DEBUG

        self.validate_remove_adapters_parameters(params)

        ca = CutadaptRunner(self.scratch)
        input_file_info = self._stage_input_file(ca, params['input_reads'], params['reads_type'])
        output_file = os.path.join(self.scratch, params['output_object_name'] + '.fq')
        ca.set_output_file(output_file)
        self._build_run(ca, params)
        report = ca.run()

        return self._package_result(output_file,
                                    params['output_object_name'],
                                    params['output_workspace'],
                                    input_file_info,
                                    report)

    def validate_remove_adapters_parameters(self, params):
        # check for required parameters
        for p in ['input_reads', 'output_workspace', 'output_object_name']:
            if p not in params:
                raise ValueError('"' + p + '" parameter is required, but missing')

        adapter_found = False
        if 'five_prime' in params and params['five_prime'] != None:
            adapter_found = True
            if 'adapter_sequence_5P' not in params['five_prime']:
                raise ValueError('"five_prime.adapter_sequence_5P" was not defined')
            if 'anchored_5P' in params['five_prime']:
                if params['five_prime']['anchored_5P'] not in [0, 1]:
                    raise ValueError('"five_prime.anchored_5P" must be either 0 or 1')

        if 'three_prime' in params and params['three_prime'] != None:
            adapter_found = True
            if 'adapter_sequence_3P' not in params['three_prime']:
                raise ValueError('"three_prime.adapter_sequence_3P" was not defined')
            if 'anchored_3P' in params['three_prime']:
                if params['three_prime']['anchored_3P'] not in [0, 1]:
                    raise ValueError('"three_prime.anchored_3P" must be either 0 or 1')

        if not adapter_found:
            raise ValueError ("Must configure at least one of 5' or 3' adapter")

        # TODO: validate values of error_tolerance and min_overlap_length


    def _stage_input_file(self, cutadapt_runner, ref, reads_type):

        ru = ReadsUtils(self.callbackURL)
        if reads_type == 'KBaseFile.PairedEndLibrary' or 'KBaseAssembly.PairedEndLibrary':
            input_file_info = ru.download_reads({
                    'read_libraries': [ref],
                    'interleaved': 'true'
                    })['files'][ref]
        elif reads_type == 'KBaseFile.SingleEndLibrary' or 'KBaseAssembly.SingleEndLibrary':
            input_file_info = ru.download_reads({
                    'read_libraries': [ref]
                    })['files'][ref]
        else:
            raise ValueError ("Can't download_reads() for object type: '"+str(reads_type)+"'")
        input_file_info['input_ref'] = ref
        file_location = input_file_info['files']['fwd']

        # DEBUG
        #with open (file_location, 'r', 0)  as fasta_file:
        #    for line in fasta_file.readlines():
        #        print ("LINE: '"+line+"'\n")

        interleaved = False
        if input_file_info['files']['type'] == 'interleaved':
            interleaved = True
        cutadapt_runner.set_interleaved(interleaved)
        cutadapt_runner.set_input_file(file_location)
        return input_file_info


    def _build_run(self, cutadapt_runner, params):
        if 'five_prime' in params:
            seq = params['five_prime']['adapter_sequence_5P']
            if seq:
                anchored = 1
                if 'anchored_5P' in params['five_prime'] and params['five_prime'] != None:
                    anchored = params['five_prime']['anchored_5P']
                cutadapt_runner.set_five_prime_option(seq, anchored)

        if 'three_prime' in params and params['three_prime'] != None:
            seq = params['three_prime']['adapter_sequence_3P']
            if seq:
                anchored = 0
                if 'anchored_3P' in params['three_prime']:
                    anchored = params['three_prime']['anchored_3P']
                cutadapt_runner.set_three_prime_option(seq, anchored)

        if 'error_tolerance' in params:
            cutadapt_runner.set_error_tolerance(params['error_tolerance'])

        if 'min_overlap_length' in params:
            cutadapt_runner.set_min_overlap(params['min_overlap_length'])

        if 'min_read_length' in params:
            cutadapt_runner.set_min_read_length(params['min_read_length'])

        if 'discard_untrimmed' in params:
            cutadapt_runner.set_discard_untrimmed(params['discard_untrimmed'])


    def _package_result(self, output_file, output_name, ws_name_or_id, data_info, report):
        upload_params = {
            'fwd_file': output_file,
            'name': output_name
        }

        if str(ws_name_or_id).isdigit():
            upload_params['wsid'] = int(ws_name_or_id)
        else:
            upload_params['wsname'] = str(ws_name_or_id)

        fields = [
            'sequencing_tech',
            'strain',
            'source',
            'read_orientation_outward',
            'insert_size_mean',
            'insert_size_std_dev'
        ]

        if 'input_ref' in data_info and data_info['input_ref'] != None and data_info['sequencing_tech']:
            upload_params['source_reads_ref'] = data_info['input_ref']
        else:
            for f in fields:
                if f in data_info:
                    upload_params[f] = data_info[f]
            if 'single_genome' in data_info:
                if data_info['single_genome'] == 'true':
                    upload_params['single_genome'] = 1
                elif data_info['single_genome'] == 'false':
                    upload_params['single_genome'] = 0
            if 'sequencing_tech' not in upload_params:
                upload_params['sequencing_tech'] = 'unknown'
            if not upload_params['sequencing_tech']:
                upload_params['sequencing_tech'] = 'unknown'


        if data_info['files']['type'] == 'interleaved':
            upload_params['interleaved'] = 1

        ru = ReadsUtils(self.callbackURL)
        result = ru.upload_reads(upload_params)

        # THE REPORT MUST BE CREATED OUTSIDE SO THAT LIBS AND SETS ARE HANDLED
        """
        # create report
        kbreport = KBaseReport(self.callbackURL)
        rep = kbreport.create({
                              'report': {
                                  'text_message': report,
                                  'objects_created': [{
                                      "ref": str(ws_name_or_id) + '/' + upload_params['name'],
                                      "description": ''
                                  }]
                              },
                              "workspace_name": str(ws_name_or_id)
                              })

        return {
            'report_ref': rep['ref'],
            'report_name': rep['name'],
            'output_reads_ref': result['obj_ref']
        }
        """
        return {
            'report': report,
            'output_reads_ref': result['obj_ref']
        }
