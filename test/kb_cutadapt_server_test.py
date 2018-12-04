# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests
import shutil
requests.packages.urllib3.disable_warnings()

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from requests_toolbelt import MultipartEncoder

from kb_cutadapt.kb_cutadaptImpl import kb_cutadapt
from kb_cutadapt.kb_cutadaptServer import MethodContext

from ReadsUtils.ReadsUtilsClient import ReadsUtils
from kb_cutadapt.authclient import KBaseAuth as _KBaseAuth


class kb_cutadaptTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_cutadapt'):
            cls.cfg[nameval[0]] = nameval[1]
        authServiceUrl = cls.cfg.get('auth-service-url',
                "https://kbase.us/services/authorization/Sessions/Login")
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_cutadapt',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.shockURL = cls.cfg['shock-url']
        cls.handleURL = cls.cfg['handle-service-url']
        cls.serviceWizardURL = cls.cfg['service-wizard-url']

        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = kb_cutadapt(cls.cfg)

        # setup data at the class level for now (so that the code is run
        # once for all tests, not before each test case.  Not sure how to
        # do that outside this function..)
        suffix = int(time.time() * 1000)
        wsName = "test_SetAPI_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': wsName})
        cls.wsName = wsName

        # handle files in test
        """
        # copy test file to scratch area
        fq_filename = "interleaved.fastq"
        fq_path = os.path.join(cls.cfg['scratch'], fq_filename)
        shutil.copy(os.path.join("data", fq_filename), fq_path)

        ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'])
        cls.read1ref = ru.upload_reads({
                'fwd_file': fq_path,
                'sequencing_tech': 'tech1',
                'wsname': wsName,
                'name': 'reads1',
                'interleaved':1
            })['obj_ref']
        """

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'shock_ids'):
            for shock_id in cls.shock_ids:
                print('Deleting SHOCK node: '+str(shock_id))
                cls.delete_shock_node(shock_id)

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_cutadapt_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    @classmethod
    def upload_file_to_shock(cls, file_path):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        header = dict()
        header["Authorization"] = "Oauth {0}".format(cls.ctx['token'])

        if file_path is None:
            raise Exception("No file given for upload to SHOCK!")

        with open(os.path.abspath(file_path), 'rb') as dataFile:
            files = {'upload': dataFile}
            response = requests.post(
                cls.shockURL + '/node', headers=header, files=files,
                stream=True, allow_redirects=True, timeout=30)

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            shock_id = result['data']['id']
            if not hasattr(cls, 'shock_ids'):
                cls.shock_ids = []
            cls.shock_ids.append(shock_id)

            return result["data"]

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.ctx['token'])}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    def getSingleEndLibInfo(self, read_lib_basename, lib_i=0):
        if hasattr(self.__class__, 'singleEndLibInfo_list'):
            try:
                info = self.__class__.singleEndLibInfo_list[lib_i]
                name = self.__class__.singleEndLibName_list[lib_i]
                if info != None:
                    if name != read_lib_basename:
                        self.__class__.singleEndLib_SetInfo[lib_i] = None
                        self.__class__.singleEndLib_SetName[lib_i] = None
                    else:
                        return info
            except:
                pass

        # 1) upload files to shock
        token = self.ctx['token']
        forward_shock_file = self.upload_file_to_shock('data/'+read_lib_basename+'.fwd.fq')
        #pprint(forward_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        # 3) save to WS
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':'artificial reads'
        }

        new_obj_info = self.wsClient.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':'test-'+str(lib_i)+'.se.reads',
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'kb_trimmomatic',
                                        'method':'test_runTrimmomatic'
                                    }
                                ]
                            }]
                        })[0]

        # store it
        if not hasattr(self.__class__, 'singleEndLibInfo_list'):
            self.__class__.singleEndLibInfo_list = []
            self.__class__.singleEndLibName_list = []
        for i in range(lib_i+1):
            try:
                assigned = self.__class__.singleEndLibInfo_list[i]
            except:
                self.__class__.singleEndLibInfo_list.append(None)
                self.__class__.singleEndLibName_list.append(None)

        self.__class__.singleEndLibInfo_list[lib_i] = new_obj_info
        self.__class__.singleEndLibName_list[lib_i] = read_lib_basename
        return new_obj_info

    def getPairedEndLibInfo(self, read_lib_basename, lib_i=0):
        if hasattr(self.__class__, 'pairedEndLibInfo_list'):
            try:
                info = self.__class__.pairedEndLibInfo_list[lib_i]
                name = self.__class__.pairedEndLibName_list[lib_i]
                if info != None:
                    if name != read_lib_basename:
                        self.__class__.singleEndLibInfo_list[lib_i] = None
                        self.__class__.singleEndLibName_list[lib_i] = None
                    else:
                        return info
            except:
                pass

        # 1) upload files to shock
        token = self.ctx['token']
        forward_shock_file = self.upload_file_to_shock('data/'+read_lib_basename+'.fwd.fq')
        reverse_shock_file = self.upload_file_to_shock('data/'+read_lib_basename+'.rev.fq')
        #pprint(forward_shock_file)
        #pprint(reverse_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        reverse_handle = hs.persist_handle({
                                        'id' : reverse_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': reverse_shock_file['file']['name'],
                                        'remote_md5': reverse_shock_file['file']['checksum']['md5']})

        # 3) save to WS
        paired_end_library = {
            'lib1': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':forward_shock_file['file']['size']
            },
            'lib2': {
                'file': {
                    'hid':reverse_handle,
                    'file_name': reverse_shock_file['file']['name'],
                    'id': reverse_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':reverse_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':reverse_shock_file['file']['size']

            },
            'interleaved':0,
            'sequencing_tech':'artificial reads'
        }

        new_obj_info = self.wsClient.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseFile.PairedEndLibrary',
                                'data':paired_end_library,
                                'name':'test-'+str(lib_i)+'.pe.reads',
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'kb_trimmomatic',
                                        'method':'test_runTrimmomatic'
                                    }
                                ]
                            }]
                        })[0]

        # store it
        if not hasattr(self.__class__, 'pairedEndLibInfo_list'):
            self.__class__.pairedEndLibInfo_list = []
            self.__class__.pairedEndLibName_list = []
        for i in range(lib_i+1):
            try:
                assigned = self.__class__.pairedEndLibInfo_list[i]
            except:
                self.__class__.pairedEndLibInfo_list.append(None)
                self.__class__.pairedEndLibName_list.append(None)

        self.__class__.pairedEndLibInfo_list[lib_i] = new_obj_info
        self.__class__.pairedEndLibName_list[lib_i] = read_lib_basename
        return new_obj_info


    # call this method to get the WS object info of a Single End Library Set (will
    # upload the example data if this is the first time the method is called during tests)
    def getSingleEndLib_SetInfo(self, read_libs_basename_list, refresh=False):
        if hasattr(self.__class__, 'singleEndLib_SetInfo'):
            try:
                info = self.__class__.singleEndLib_SetInfo
                if info != None:
                    if refresh:
                        self.__class__.singleEndLib_SetInfo = None
                    else:
                        return info
            except:
                pass

        # build items and save each PairedEndLib
        items = []
        for lib_i,read_lib_basename in enumerate (read_libs_basename_list):
            label    = read_lib_basename
            lib_info = self.getSingleEndLibInfo (read_lib_basename, lib_i)
            lib_ref  = str(lib_info[6])+'/'+str(lib_info[0])
            print ("LIB_REF["+str(lib_i)+"]: "+lib_ref+" "+read_lib_basename)  # DEBUG

            items.append({'ref': lib_ref,
                          'label': label
                          #'data_attachment': ,
                          #'info':
                         })

        # save readsset
        desc = 'test ReadsSet'
        readsSet_obj = { 'description': desc,
                         'items': items
                       }
        name = 'TEST_READSET'

        new_obj_info = self.wsClient.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseSets.ReadsSet',
                                'data':readsSet_obj,
                                'name':name,
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'kb_trimmomatic',
                                        'method':'test_runTrimmomatic'
                                    }
                                ]
                            }]
                        })[0]

        # store it
        self.__class__.singleEndLib_SetInfo = new_obj_info
        return new_obj_info


    # call this method to get the WS object info of a Paired End Library Set (will
    # upload the example data if this is the first time the method is called during tests)
    def getPairedEndLib_SetInfo(self, read_libs_basename_list, refresh=False):
        if hasattr(self.__class__, 'pairedEndLib_SetInfo'):
            try:
                info = self.__class__.pairedEndLib_SetInfo
                if info != None:
                    if refresh:
                        self.__class__.pairedEndLib_SetInfo[lib_i] = None
                    else:
                        return info
            except:
                pass

        # build items and save each PairedEndLib
        items = []
        for lib_i,read_lib_basename in enumerate (read_libs_basename_list):
            label    = read_lib_basename
            lib_info = self.getPairedEndLibInfo (read_lib_basename, lib_i)
            lib_ref  = str(lib_info[6])+'/'+str(lib_info[0])
            print ("LIB_REF["+str(lib_i)+"]: "+lib_ref+" "+read_lib_basename)  # DEBUG

            items.append({'ref': lib_ref,
                          'label': label
                          #'data_attachment': ,
                          #'info':
                         })

        # save readsset
        desc = 'test ReadsSet'
        readsSet_obj = { 'description': desc,
                         'items': items
                       }
        name = 'TEST_READSET'

        new_obj_info = self.wsClient.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseSets.ReadsSet',
                                'data':readsSet_obj,
                                'name':name,
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'kb_trimmomatic',
                                        'method':'test_runTrimmomatic'
                                    }
                                ]
                            }]
                        })[0]

        # store it
        self.__class__.pairedEndLib_SetInfo = new_obj_info
        return new_obj_info


    ##############
    # UNIT TESTS #
    ##############

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa


    ### TEST 1: run Cutadapt against just one single end library
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_SE_Lib")
    def test_basic_options_SE_Lib(self):

        print ("\n\nRUNNING: test_basic_options_SE_Lib()")
        print ("====================================\n\n")

        input_libs = ['cutadapt_1']
        output_name = 'trim5p.SELib'

        se_lib_info = self.getSingleEndLibInfo(input_libs[0])
        se_lib_ref = str(se_lib_info[6])+'/'+str(se_lib_info[0])

        p1 = {
            'input_reads': se_lib_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': {
                'adapter_sequence_5P': 'TGCCCTGCAAAAACGTCTGGAAA',
                'anchored_5P': 1
            },
            'three_prime': None
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p1)
        pprint(ret)

        # check the output
        single_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':se_lib_info[7] + '/' + single_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],single_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseFile.SingleEndLibrary')


    ### TEST 2: run Cutadapt against just one paired end library
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_PE_Lib")
    def test_basic_options_PE_Lib(self):

        print ("\n\nRUNNING: test_basic_options_PE_Lib()")
        print ("====================================\n\n")

        input_libs = ['cutadapt_1']
        output_name = 'trim5p.PELib'

        pe_lib_info = self.getPairedEndLibInfo(input_libs[0])
        pe_lib_ref = str(pe_lib_info[6])+'/'+str(pe_lib_info[0])

        p2 = {
            'input_reads': pe_lib_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': {
                'adapter_sequence_5P': 'TGCCCTGCAAAAACGTCTGGAAA',
                'anchored_5P': 1
            },
            'three_prime': None
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p2)
        pprint(ret)

        # check the output
        paired_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + paired_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],paired_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseFile.PairedEndLibrary')


    ### TEST 3: run Cutadapt against single end reads set
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_SE_ReadsSet")
    def test_basic_options_SE_ReadsSet(self):

        print ("\n\nRUNNING: test_basic_options_SE_ReadsSet()")
        print ("=========================================\n\n")

        input_libs = ['cutadapt_1','cutadapt_2']
        output_name = 'trim5p.SERS'

        se_lib_set_info = self.getSingleEndLib_SetInfo(input_libs)
        se_lib_set_ref = str(se_lib_set_info[6])+'/'+str(se_lib_set_info[0])

        p3 = {
            'input_reads': se_lib_set_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': {
                'adapter_sequence_5P': 'TGCCCTGCAAAAACGTCTGGAAA',
                'anchored_5P': 1
            },
            'three_prime': None
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p3)
        pprint(ret)

        # check the output
        single_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':se_lib_set_info[7] + '/' + single_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],single_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseSets.ReadsSet')


    ### TEST 4: run Cutadapt against paired end reads set
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_PE_ReadsSet")
    def test_basic_options_PE_ReadsSet(self):

        print ("\n\nRUNNING: test_basic_options_PE_ReadsSet()")
        print ("=========================================\n\n")

        input_libs = ['cutadapt_1','cutadapt_2']
        output_name = 'trim5p.PERS'

        pe_lib_set_info = self.getPairedEndLib_SetInfo(input_libs)
        pe_lib_set_ref = str(pe_lib_set_info[6])+'/'+str(pe_lib_set_info[0])

        p4 = {
            'input_reads': pe_lib_set_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': {
                'adapter_sequence_5P': 'TGCCCTGCAAAAACGTCTGGAAA',
                'anchored_5P': 1
            },
            'three_prime': None
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p4)
        pprint(ret)

        # check the output
        paired_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_set_info[7] + '/' + paired_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],paired_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseSets.ReadsSet')


    ### TEST 5: run Cutadapt against just one paired end library with 3 prime adapter (anchored)
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_PE_Lib_threeprime_anchored")
    def test_basic_options_PE_Lib_threeprime_anchored(self):

        print ("\n\nRUNNING: test_basic_options_PE_Lib_threeprime_anchored()")
        print ("========================================================\n\n")

        input_libs = ['cutadapt_1']
        output_name = 'trim3p_anchored.PELib'

        pe_lib_info = self.getPairedEndLibInfo(input_libs[0])
        pe_lib_ref = str(pe_lib_info[6])+'/'+str(pe_lib_info[0])

        p2 = {
            'input_reads': pe_lib_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': None,
            'three_prime': {
                'adapter_sequence_3P': 'ACGTACGTACGTAAA',
                'anchored_3P': 1
            },
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p2)
        pprint(ret)

        # check the output
        paired_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + paired_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],paired_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseFile.PairedEndLibrary')


    ### TEST 6: run Cutadapt against just one paired end library with 3 prime adapter (unanchored)
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_PE_Lib_threeprime_UNanchored")
    def test_basic_options_PE_Lib_threeprime_UNanchored(self):

        print ("\n\nRUNNING: test_basic_options_PE_Lib_threeprime_UNanchored()")
        print ("==========================================================\n\n")

        input_libs = ['cutadapt_1']
        output_name = 'trim3p_unanchored.PELib'

        pe_lib_info = self.getPairedEndLibInfo(input_libs[0])
        pe_lib_ref = str(pe_lib_info[6])+'/'+str(pe_lib_info[0])

        p2 = {
            'input_reads': pe_lib_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': None,
            'three_prime': {
                'adapter_sequence_3P': 'ACGTACGTACGT',
                'anchored_3P': 0
            },
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p2)
        pprint(ret)

        # check the output
        paired_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + paired_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],paired_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseFile.PairedEndLibrary')


    ### TEST 7: run Cutadapt against just one paired end library with 5 prime unanchored
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_basic_options_PE_Lib_fiveprime_UNanchored")
    def test_basic_options_PE_Lib_fiveprime_UNanchored(self):

        print ("\n\nRUNNING: test_basic_options_PE_Lib_fiveprime_UNanchored()")
        print ("=========================================================\n\n")

        input_libs = ['cutadapt_1']
        output_name = 'trim5p.PELib'

        pe_lib_info = self.getPairedEndLibInfo(input_libs[0])
        pe_lib_ref = str(pe_lib_info[6])+'/'+str(pe_lib_info[0])

        p2 = {
            'input_reads': pe_lib_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 0,
            'five_prime': {
                'adapter_sequence_5P': 'TGCAAAAACGTCTGGAAA',
                'anchored_5P': 0
            },
            'three_prime': None
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p2)
        pprint(ret)

        # check the output
        paired_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + paired_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],paired_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseFile.PairedEndLibrary')



    ### TEST 8: run Cutadapt against just one paired end library with discard_untrimmed
    #
    # Uncomment to skip this test
    # @unittest.skip("skipped test_discard_untrimmed_option_PE_Lib")
    def test_discard_untrimmed_option_PE_Lib(self):

        print ("\n\nRUNNING: test_discard_untrimmed_option_PE_Lib()")
        print ("===============================================\n\n")

        input_libs = ['cutadapt_1']
        output_name = 'trim5p.PELib'

        pe_lib_info = self.getPairedEndLibInfo(input_libs[0])
        pe_lib_ref = str(pe_lib_info[6])+'/'+str(pe_lib_info[0])

        p2 = {
            'input_reads': pe_lib_ref,
            'output_workspace': self.getWsName(),
            'output_object_name': output_name,
            'min_read_length': 50,
            'discard_untrimmed': 1,
            'five_prime': {
                'adapter_sequence_5P': 'TGCCCTGCAAAAACGTCTGGAAA',
                'anchored_5P': 1
            },
            'three_prime': None
        }

        ret = self.getImpl().remove_adapters(self.getContext(), p2)
        pprint(ret)

        # check the output
        paired_output_name = output_name
        info_list = self.wsClient.get_object_info([{'ref':pe_lib_info[7] + '/' + paired_output_name}], 1)
        self.assertEqual(len(info_list),1)
        output_reads_info = info_list[0]
        self.assertEqual(output_reads_info[1],paired_output_name)
        self.assertEqual(output_reads_info[2].split('-')[0],'KBaseFile.PairedEndLibrary')

