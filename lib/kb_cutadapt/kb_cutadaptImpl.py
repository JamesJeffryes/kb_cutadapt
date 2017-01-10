# -*- coding: utf-8 -*-
#BEGIN_HEADER

import os

from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService

from kb_cutadapt.CutadaptUtil import CutadaptUtil
from SetAPI.SetAPIServiceClient import SetAPI

#END_HEADER


class kb_cutadapt:
    '''
    Module Name:
    kb_cutadapt

    Module Description:
    A KBase module: cutadapt
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.1"
    GIT_URL = "https://github.com/dcchivian/kb_cutadapt"
    GIT_COMMIT_HASH = "a53fd77dab5f6c39c25c47560375dedfb467bf7d"

    #BEGIN_CLASS_HEADER

    def log(self, target, message):
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = os.path.abspath(config['scratch'])
        self.data = os.path.abspath(config['data'])
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        if self.config['SDK_CALLBACK_URL'] == None:
            raise ValueError ("SDK_CALLBACK_URL not set in environment")
        #END_CONSTRUCTOR
        pass


    def remove_adapters(self, ctx, params):
        """
        :param params: instance of type "RemoveAdaptersParams" -> structure:
           parameter "output_workspace" of String, parameter
           "output_object_name" of String, parameter "input_reads" of type
           "ws_ref" (@ref ws), parameter "five_prime" of type
           "FivePrimeOptions" (unfortunately, we have to name the fields
           uniquely between 3' and 5' options due to the current
           implementation of grouped parameters) -> structure: parameter
           "adapter_sequence_3P" of String, parameter "anchored_3P" of type
           "boolean" (@range (0, 1)), parameter "three_prime" of type
           "ThreePrimeOptions" -> structure: parameter "adapter_sequence_5P"
           of String, parameter "anchored_5P" of type "boolean" (@range (0,
           1)), parameter "error_tolerance" of Double, parameter
           "min_overlap_length" of Long
        :returns: instance of type "RemoveAdaptersResult" -> structure:
           parameter "report_ref" of String, parameter "output_reads_ref" of
           String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN remove_adapters
        console = []
        self.log(console, 'Running remove_adapters() with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        wsClient = workspaceService(self.config['workspace-url'], token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        #SERVICE_VER = 'dev'  # DEBUG
        SERVICE_VER = 'release'

        # param checks
        required_params = ['output_workspace',
                           'input_reads', 
                           'output_object_name'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['input_reads'])]


        # RUN
        exec_remove_adapters_retVal = self.exec_remove_adapters (ctx, params)[0]


        # build report
        #
        reportObj = {'objects_created':[], 
                     'text_message':''}

        # text report
        try:
            reportObj['text_message'] = exec_remove_adapters_retVal['report']
        except:
            raise ValueError ("no report generated by exec_remove_adapters()")

        # output object
        if exec_remove_adapters_retVal['output_reads_ref'] != None:
            reportObj['objects_created'].append({'ref':exec_remove_adapters_retVal['output_reads_ref'],
                                                 'description':'Post Cutadapt Reads'})
        else:
            raise ValueError ("no output generated by exec_remove_adapters()")

        # save report object
        report = KBaseReport(self.config['SDK_CALLBACK_URL'], token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create({'report':reportObj, 'workspace_name':params['output_workspace']})

        result = { 'output_reads_ref': exec_remove_adapters_retVal['output_reads_ref'], 'report_ref': report_info['ref'] }
        #END remove_adapters

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method remove_adapters return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def exec_remove_adapters(self, ctx, params):
        """
        :param params: instance of type "RemoveAdaptersParams" -> structure:
           parameter "output_workspace" of String, parameter
           "output_object_name" of String, parameter "input_reads" of type
           "ws_ref" (@ref ws), parameter "five_prime" of type
           "FivePrimeOptions" (unfortunately, we have to name the fields
           uniquely between 3' and 5' options due to the current
           implementation of grouped parameters) -> structure: parameter
           "adapter_sequence_3P" of String, parameter "anchored_3P" of type
           "boolean" (@range (0, 1)), parameter "three_prime" of type
           "ThreePrimeOptions" -> structure: parameter "adapter_sequence_5P"
           of String, parameter "anchored_5P" of type "boolean" (@range (0,
           1)), parameter "error_tolerance" of Double, parameter
           "min_overlap_length" of Long
        :returns: instance of type "exec_RemoveAdaptersResult" -> structure:
           parameter "report" of String, parameter "output_reads_ref" of
           String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN exec_remove_adapters
        console = []
        self.log(console, 'Running exec_remove_adapters() with parameters: ')
        self.log(console, "\n"+pformat(params))
        report = ''
        returnVal = dict()
        returnVal['output_reads_ref'] = None

        token = ctx['token']
        wsClient = workspaceService(self.config['workspace-url'], token=token)
        headers = {'Authorization': 'OAuth '+token}
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # param checks
        required_params = ['output_workspace',
                           'input_reads', 
                           'output_object_name'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"'")

        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[str(params['input_reads'])]

        # Determine whether read library or read set is input object
        #
        try:
            # object_info tuple
            [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)

            input_reads_obj_info = wsClient.get_object_info_new ({'objects':[{'ref':params['input_reads']}]})[0]
            input_reads_obj_type = input_reads_obj_info[TYPE_I]
            input_reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", input_reads_obj_type)  # remove trailing version
            #input_reads_obj_version = input_reads_obj_info[VERSION_I]  # this is object version, not type version
        except Exception as e:
            raise ValueError('Unable to get read library object from workspace: (' + str(params['input_reads_ref']) +')' + str(e))

        acceptable_types = ["KBaseSets.ReadsSet", "KBaseFile.PairedEndLibrary", "KBaseFile.SingleEndLibrary"]
        if input_reads_obj_type not in acceptable_types:
            raise ValueError ("Input reads of type: '"+input_reads_obj_type+"'.  Must be one of "+", ".join(acceptable_types))


        # get set
        #
        readsSet_ref_list = []
        readsSet_names_list = []
        if input_reads_obj_type != "KBaseSets.ReadsSet":
            readsSet_ref_list = [params['input_reads']]
            readsSet_names_list = [params['output_object_name']]
        else:
            try:
                #setAPI_Client = SetAPI (url=self.config['SDK_CALLBACK_URL'], token=ctx['token'])  # for SDK local.  doesn't work for SetAPI
                setAPI_Client = SetAPI (url=self.config['service-wizard-url'], token=ctx['token'])  # for dynamic service
                input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':params['input_reads'],'include_item_info':1})

            except Exception as e:
                raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(params['input_reads'])+")\n" + str(e))
            for readsLibrary_obj in input_readsSet_obj['data']['items']:
                readsSet_ref_list.append(readsLibrary_obj['ref'])
                NAME_I = 1
                readsSet_names_list.append(readsLibrary_obj['info'][NAME_I])


        # Iterate through readsLibrary memebers of set
        #
        report = ''
        cutadapt_readsSet_ref  = None
        cutadapt_readsLib_refs = []

        for reads_item_i,input_reads_library_ref in enumerate(readsSet_ref_list):
            exec_remove_adapters_OneLibrary_params = { 'output_workspace': params['output_workspace'],
                                                        'input_reads': input_reads_library_ref
                                     }
            if input_reads_obj_type != "KBaseSets.ReadsSet":
                exec_remove_adapters_OneLibrary_params['output_object_name'] = params['output_object_name']
            else:
                exec_remove_adapters_OneLibrary_params['output_object_name'] = readsSet_names_list[reads_item_i]

            msg = "\n\nRUNNING exec_remove_adapters_OneLibrary() ON LIBRARY: "+str(input_reads_library_ref)+" "+str(readsSet_names_list[reads_item_i])+"\n"
            msg += "----------------------------------------------------------------------------\n"
            report += msg
            self.log (console, msg)

            # RUN
            exec_remove_adapters_OneLibrary_retVal = self.exec_remove_adapters_OneLibrary (ctx, exec_remove_adapters_OneLibrary_params)[0]

            report += exec_remove_adapters_OneLibrary_retVal['report']+"\n\n"
            cutadapt_readsLib_refs.append (exec_remove_adapters_OneLibrary_retVal['output_reads_ref'])


        # Just one Library
        if input_reads_obj_type != "KBaseSets.ReadsSet":

            # create return output object
            result = { 'report': report,
                       'output_reads_ref': cutadapt_readsLib_refs[0],
                       }
        # ReadsSet
        else:

            # save cutadapt readsSet
            some_cutadapt_output_created = False
            items = []
            for i,lib_ref in enumerate(cutadapt_readsLib_refs):

                if lib_ref == None:
                    #items.append(None)  # can't have 'None' items in ReadsSet
                    continue
                else:
                    some_cutadapt_output_created = True
                    try:
                        label = input_readsSet_obj['data']['items'][i]['label']
                    except:
                        NAME_I = 1
                        label = wsClient.get_object_info_new ({'objects':[{'ref':lib_ref}]})[0][NAME_I]
                    label = label + "_cutadapt"

                    items.append({'ref': lib_ref,
                                  'label': label
                                  #'data_attachment': ,
                                  #'info':
                                      })
            if some_cutadapt_output_created:
                reads_desc_ext = " + Cutadapt"
                reads_name_ext = "_cutadapt"
                output_readsSet_obj = { 'description': input_readsSet_obj['data']['description']+reads_desc_ext,
                                        'items': items
                                        }
                output_readsSet_name = str(params['output_object_name'])+reads_name_ext
                cutadapt_readsSet_ref = setAPI_Client.save_reads_set_v1 ({'workspace_name': params['output_workspace'],
                                                                          'output_object_name': output_readsSet_name,
                                                                          'data': output_readsSet_obj
                                                                          })['set_ref']
            else:
                raise ValueError ("No cutadapt output created")


            # create return output object
            result = { 'report': report,
                       'output_reads_ref': cutadapt_readsSet_ref
                       }
        #END exec_remove_adapters

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method exec_remove_adapters return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

    def exec_remove_adapters_OneLibrary(self, ctx, params):
        """
        :param params: instance of type "RemoveAdaptersParams" -> structure:
           parameter "output_workspace" of String, parameter
           "output_object_name" of String, parameter "input_reads" of type
           "ws_ref" (@ref ws), parameter "five_prime" of type
           "FivePrimeOptions" (unfortunately, we have to name the fields
           uniquely between 3' and 5' options due to the current
           implementation of grouped parameters) -> structure: parameter
           "adapter_sequence_3P" of String, parameter "anchored_3P" of type
           "boolean" (@range (0, 1)), parameter "three_prime" of type
           "ThreePrimeOptions" -> structure: parameter "adapter_sequence_5P"
           of String, parameter "anchored_5P" of type "boolean" (@range (0,
           1)), parameter "error_tolerance" of Double, parameter
           "min_overlap_length" of Long
        :returns: instance of type "exec_RemoveAdaptersResult" -> structure:
           parameter "report" of String, parameter "output_reads_ref" of
           String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN exec_remove_adapters_OneLibrary
        print('Running cutadapt.remove_adapters')
        pprint(params)

        cutadapt = CutadaptUtil(self.config)
        result = cutadapt.remove_adapters(params)
        #END exec_remove_adapters_OneLibrary

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method exec_remove_adapters_OneLibrary return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
