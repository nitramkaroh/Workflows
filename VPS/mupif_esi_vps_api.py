"""
#-------------------------------------------------------------------------------
# Python PACKAGE/MODULE
#-------------------------------------------------------------------------------
#   NAME        : MUPIF_VPS_API
#   MAIN AUTHOR : SLL
#   REVIEWER    :
#   LAST CHANGE : 2018-02-05
#   FUNCTIONS   : write_py_inc, write_ori, open_erf, read_THNOD
#   DESCRIPTION : I/O funtions for VPS solver files
#
#
#   History Changes:
#    . 2018-02-05 by SLL: Setup of Python Package
#
#-------------------------------------------------------------------------------
"""
#
from builtins import str
import sys
from mupif import *
import os
import logging
log = logging.getLogger()
import mupif.Physics.PhysicalQuantities as PQ
import re
import os
import numpy as np
import subprocess
import shutil
import h5py
import itertools
import time
import Pyro4
import vpsConfig
from scipy.interpolate import interp1d
from datetime import datetime, timedelta

re_newkey = re.compile('^[\w|\s]{6}\/')

@Pyro4.expose
class VPS_API(Application.Application):


    def __init__ (self,file='',workdir='',modelID=1):

        super(VPS_API, self).__init__()

        modelID=3

        model  = vpsConfig.VPS_MODEL(modelID=modelID)

        self.solverProp = vpsConfig.VPS_SOLVER()

        try:
            print ("Model name            : %s" % (model.modelName))
            self.modelName = model.modelName
        except:
            raise APIError.APIError('Model ID = %d is not valid.' % (modelID))
        try:
            print ("Working directory     : %s" % (workdir))
            self.workDir = workdir
        except:
            raise APIError.APIError('Working directory not defined in input data.')
        try:
            print ("Model main input file : %s" % (model.inputFile))

            self.inputFile  = os.path.normpath(workdir + '/' + os.path.basename(model.inputFile))
            self.outputFile = os.path.normpath(workdir + '/' + os.path.splitext(os.path.basename(model.inputFile))[0] + '.out')
            self.resultFile = os.path.normpath(workdir + '/' + os.path.splitext(os.path.basename(model.inputFile))[0] + '_RESULT.erfh5')
        except:
            raise APIError.APIError('Model main input file is not defined in input data.')

        print ("\nStart reading input files: \n")

        Properties,Orientation = read_pc(model.inputFile,workdir)

        print ("\nDone!\n")


        print('\nFound MUPIF Input Properties:')
        print('-----------------------------\n')
        print('         #      Name     Value\n')
        for PID in Properties:
            if Properties[PID]['IO']=='input':
                print('%10d%10s%10g'% (PID,Properties[PID]['name'],Properties[PID]['value']))

        print('\nFound MUPIF Output Properties:')
        print('-----------------------------\n')
        print('         #\n')
        for PID in Properties:
            if Properties[PID]['IO']=='output':
                print('%10d'% (PID))

        self.Properties = Properties

        print('\nFound Orientation files:')
        print('--------------------------\n')
        if Orientation:
            print('ORI file name : %s' %(Orientation))
        else:
            print('NONE')

        self.orientationFile = Orientation

        # Initialize solver status flag
        # 0: solver has not been exectuted
        # 1: solver is running
        # 2: solver terminated normally
        # 3: solver terminated with an error
        self.solverStatus = 0

    def getField(self, fieldID, time):
        """
        Returns the requested field at given time. Field is identified by fieldID.

        :param FieldID fieldID: Identifier of the field
        :param float time: Target time

        :return: Returns requested field.
        :rtype: Field
        """
        if fieldID==10001:
            return readERFfield(self.resultFile,'DISPLACEMENT',time)
        else:
            raise APIError.APIError ('Field ID is not supported.')

    def setField(self, field):
        """
        Registers the given (remote) field in application.

        :param Field field: Remote field to be registered by the application
        """


    def getProperty(self, propID, time, objectID=0):
        """
        Returns property identified by its ID evaluated at given time.

        :param PropertyID propID: property ID
        :param float time: Time when property should to be evaluated
        :param int objectID: Identifies object/submesh on which property is evaluated (optional, default 0)

        :return: Returns representation of requested property
        :rtype: Property
        """

        if not propID in self.Properties:
            raise APIError.APIError ('PropertyID is not valid for current model')
        else:
            if self.Properties[propID]['IO']=='input':
                return Property.ConstantProperty(self.Properties[propID]['value'], self.Properties[propID]['MUPIF_PID'], self.Properties[propID]['type'], self.Properties[propID]['unit'], None)
            elif self.Properties[propID]['IO']=='output':

                if propID==91007:
                    BCK_ERF = os.path.normpath(self.workDir + '/' + 'BUCKLING_EV_BKM.erfh5')
                    OUTPUTtype = 'BUCKL'
                    result,result_time = readERF(BCK_ERF,OUTPUTtype,None,None,None)
                    return Property.ConstantProperty(result, self.Properties[propID]['MUPIF_PID'], ValueType.Vector, self.Properties[propID]['unit'], None)

            else:
                raise APIError.APIError ('Can not return property, IO type is not defined.')



    def setProperty(self, Property, objectID=0):
        """
        Register given property in the application

        :param Property property: Setting property
        :param int objectID: Identifies object/submesh on which property is evaluated (optional, default 0)
        """
        PID = Property.getPropertyID().value

        if not PID in self.Properties:
            raise APIError.APIError ('PropertyID is not in the model property list.')
        else:
            if self.Properties[PID]['IO']=='input':
                incfile_path = self.Properties[PID]['incfile_path']
                incfile = open(incfile_path, 'r')
                inc_lines = incfile.readlines()
                incfile.close()

                line_number = self.Properties[PID]['incfile_line_mumber']
                PYVARname   = self.Properties[PID]['name']
                PYVARvalue  = Property.getValue()

                linesplit = inc_lines[line_number-1].split('#')
                if PYVARname in linesplit[0]:
                    linesplit[0] = ('%s = %g' %(PYVARname,PYVARvalue))
                    if len(linesplit[0])<23:
                        linesplit[0] = linesplit[0].ljust(23)
                    self.Properties[PID]['value'] = Property.getValue()
                    inc_lines[line_number-1] = '#'.join(linesplit)

                    # Write include file
                    f = open(incfile_path, 'w')
                    f.writelines(inc_lines)
                    f.close()

                    # Reset solver status
                    self.solverStatus = 0

            else:
                raise APIError.APIError ('PropertyID does not correspond to a input property.')

    def getFunction(self, funcID, objectID=0):
        """
        Returns function identified by its ID

        :param FunctionID funcID: function ID
        :param int objectID: Identifies optional object/submesh on which property is evaluated (optional, default 0)

        :return: Returns requested function
        :rtype: Function
        """
        raise APIError.APIError ('Unknown funcID')

    def setFunction(self, func, objectID=0):
        """
        Register given function in the application.

        :param Function func: Function to register
        :param int objectID: Identifies optional object/submesh on which property is evaluated (optional, default 0)
        """
        raise APIError.APIError ('Unknown funcID')



    def solveStep(self, runInBackground=False):
        """
        Solves the problem for given time step.

        Proceeds the solution from actual state to given time.
        The actual state should not be updated at the end, as this method could be
        called multiple times for the same solution step until the global convergence
        is reached. When global convergence is reached, finishStep is called and then
        the actual state has to be updated.
        Solution can be split into individual stages identified by optional stageID parameter.
        In between the stages the additional data exchange can be performed.
        See also wait and isSolved services.

        :param TimeStep tstep: Solution step
        :param int stageID: optional argument identifying solution stage (default 0)
        :param bool runInBackground: optional argument, defualt False. If True, the solution will run in background (in separate thread or remotely).

        """

        start = "\"{:}\" -fp {:d} -nt {:d} -wd \"{:}\" \"{:}\" > \"{:}\" 2>&1".format(\
                self.solverProp.solverPath,\
                self.solverProp.fp,\
                self.solverProp.nt,\
                os.path.normpath(self.workDir),\
                os.path.normpath(self.inputFile),\
                os.path.normpath(self.outputFile))

        if os.path.exists(os.path.normpath(self.workDir) + '\\signal'):
            os.remove(os.path.normpath(self.workDir) + '\\signal')


        print('\nStart VPS solver run:')
        print('---------------------\n')

        print('Solver path       : %s' % self.solverProp.solverPath)
        print('Precision         : %d' % self.solverProp.fp)
        print('Number of threads : %d\n' % self.solverProp.nt)

        solverRun = subprocess.Popen(start, shell = True)

        self.solverRun = solverRun

        poll = solverRun.poll()
        if poll == None:
            self.solverStatus = 1

        if runInBackground:
            return solverRun
        else:
            # Wait until solver has finished
            t_sleep = 0

            spinner = itertools.cycle(['-', '/', '|', '\\'])

            while solverRun.poll() == None:
                text = "\r{:13}: {:} ... {:}".format("Running job", self.modelName, next(spinner))
                sys.stdout.write(text)
                #sys.stdout.write(next(spinner))   # write the next character
                sys.stdout.flush()                # flush stdout buffer (actual character display)
                time.sleep(1)
                t_sleep += 1
                sys.stdout.write('\r')            # erase the last written char


#                if t_sleep in itertools.chain(range(0,61), [300, 600, 1800, 3600]):
#                    text = "\r{:13}: {:} ... {:>5d}s".format("Running job", self.modelName, t_sleep)
#                    sys.stdout.write(text)
#                    sys.stdout.flush()

            sec = timedelta(seconds=t_sleep)

            text = "\n{:14} {:} {:} {:} {:}".format(" ", self.modelName,"finished in ", str(sec),"[HH:MM:SS]")
            sys.stdout.write(text)
            sys.stdout.flush()

            self.isSolved()

            return True

#    def solverStatus(self):
#        """
#        Get the solver status
#        """

    def solverSignal(self,signal):
        """
        Get the solver status
        """
        if signal == 'PLOT': # Write the current model state to a mesh plot file
            f = open((self.workDir + 'signal'), 'w')
            f.write('PLOT')
            f.close()
        elif signal == 'QUIT': # QUIT solver run
            f = open((self.workDir + 'signal'), 'w')
            f.write('QUIT')
            f.close()
        elif signal == 'STOP': # Write RESTART file and QUIT solver run
            f = open((self.workDir + 'signal'), 'w')
            f.write('STOP')
            f.close()
        elif signal == 'EXIT': # Combines PLOT AND STOP
            f = open((self.workDir + 'signal'), 'w')
            f.write('QUIT')
            f.close()
        else:
            raise APIError.APIError ('Unsopported solver signal.')

    def isSolved(self):
        """
        Check whether solve has completed.

        :return: Returns true or false depending whether solve has completed when executed in background.
        :rtype: bool
        """

        if self.solverStatus==0:
            return False
        elif self.solverStatus==1:
            lastoutline = tail(os.path.normpath(self.outputFile))
            if 'NORMAL TERMINATION' in str(lastoutline):
                self.solverStatus = 2
                return True
            elif 'ERROR' in str(lastoutline):
                self.solverStatus = 3
                return False
            else:
                print('WARNING : Could not read proper termination message. Check the out file.\n')
                return False
        elif self.solverStatus==2:
            return True
        else:
            return False

    def getAPIVersion(self):
        """
        :return: Returns the supported API version
        :rtype: str, int
        """
        return 1
    def getApplicationSignature(self):
        """
        Get application signature.

        :return: Returns the application identification
        :rtype: str
        """
        return "ESI VPS COMPOSELECTOR Workstation in Stuttgart (Germany)"

    def terminate(self):
        """
        Terminates the application. Shutdowns daemons if created internally.
        """

        f = open((self.workDir + 'signal'), 'w')
        f.write('QUIT')
        f.close()

        if hasattr(self, 'pyroDaemon'):
            if self.pyroDaemon:
                self.pyroDaemon.unregister(self)
                if not self.externalDaemon:
                    self.pyroDaemon.shutdown()


    def getURI(self):
        """
        :return: Returns the application URI or None if application not registered in Pyro
        :rtype: str
        """
        return self.pyroURI

#%% VPS IO functions

def read_pc(pc_fpath,workDir):
    """
    Read PYVAR vlaues and put .
    :param pc_fpath: file path with pamcrash input to read
    :return: dictionary with supported keywords and their definition parameters, list with ignored keywords
    """

    if not os.path.exists(pc_fpath):
        error_msg = ('ERROR: Input file : %s does not exist!\n' % (pc_fpath))
        raise  APIError.APIError(error_msg)

    MUPIFprop = {}

    includes = [pc_fpath]

    act_inc_dir = os.path.dirname(pc_fpath)

    i = 0


    while i <= len(includes) - 1:
        line_no = 0
        card_line_no = 0
        act_keyword = ''
        act_include = includes[i]

        incfile = open(act_include, 'r')
        inclines = incfile.readlines()
        incfile.close()

        act_inc_dir      = os.path.dirname(act_include)
        act_inc_basename = os.path.basename(act_include)

        shutil.copy(act_include, workDir, follow_symlinks=True)

        act_include = os.path.normpath(workDir + '/' + act_inc_basename)

        print('Reading input file : %s' % act_include)

        for line in inclines:
            line_no += 1

            if re_newkey.match(line):
                card_line_no = 1
                cardlines = {}
                act_keyword = line[:7]
            else:
                card_line_no += 1

            if act_keyword == 'INCLU /':
                if line.startswith('$') or line.startswith('#'):
                    continue
                if card_line_no == 1:
                    incname = line[8:].strip()
                    incpath = os.path.join(act_inc_dir, incname)
                    includes.append(incpath)

            elif act_keyword == 'PYVAR /':
                if 'MUPIF_INPUT_PID' in line:
                    if line.startswith('$') or line.startswith('#'):
                        continue
                    try:
                        linesplit = line.partition('#')

                        PYVARname  = re.sub(r"\s+$", "", linesplit[0].split('=')[0], flags=re.UNICODE)
                        PYVARvalue = float(linesplit[0].split('=')[-1])
                        PIDnumber  = int(linesplit[-1].split('MUPIF_INPUT_PID=')[1])

                        MUPIFattributes = getMUPIFpropID(PIDnumber)
                        if MUPIFattributes:
                            MUPIFprop[PIDnumber] = dict(IO = 'input', type = ValueType.Scalar, value = PYVARvalue , name = PYVARname,incfile_path = act_include, incfile_line_mumber = line_no, **MUPIFattributes)
                    except:
                        error_msg = ('ERROR: Cant read MUPIF_PID from input line # %d :  %s  in file : %s  \n' % ( line_no,line,act_include))
                        raise  APIError.APIError(error_msg)
                elif 'MUPIF_OUTPUT' in line:
                    try:
                        OUTPUTpid  = False
                        OUTPUTtype = False
                        OUTPUTid   = False
                        OUTPUTcomp = False
                        linesplit = line.split('/')
                        for s in linesplit:
                            if 'PID' in s:
                                OUTPUTpid  = int(s.split('=')[1])
                            elif 'TYPE' in s:
                                OUTPUTtype = re.sub(r"\s+", "", s.split('=')[1], flags=re.UNICODE)
                            elif 'VPSID' in s:
                                vpsid = re.sub(r"\s+", "", s.split('=')[1], flags=re.UNICODE)
                                if vpsid.isdigit():
                                    OUTPUTid   = int(vpsid)
                                else:
                                    OUTPUTid   = vpsid
                            elif 'COMPONENT' in s:
                                if 'ALL' in s:
                                    OUTPUTcomp = re.sub(r"\s+", "", s.split('=')[1], flags=re.UNICODE)
                                else:
                                    OUTPUTcomp = int(s.split('=')[1])

                        if OUTPUTpid: # and OUTPUTtype and OUTPUTid and OUTPUTcomp:
                            MUPIFattributes = getMUPIFpropID(OUTPUTpid)
                            if MUPIFattributes:
                                #MUPIFprop[OUTPUTpid] = dict(IO = 'output', type = ValueType.Vector , name = OUTPUTtype, OUTPUTtype = OUTPUTtype, OUTPUTid = OUTPUTid, OUTPUTcomp = OUTPUTcomp, **MUPIFattributes)
                                MUPIFprop[OUTPUTpid] = dict(IO = 'output', OUTPUTid = OUTPUTid, **MUPIFattributes)
                    except:
                        error_msg = ('ERROR: Cant read MUPIF_OUTPUT from input line # %d :  %s  in file : %s  \n' % ( line_no,line,act_include))
                        raise  APIError.APIError(error_msg)

            elif act_keyword == 'ORTHF /':
                if line.startswith('$') or line.startswith('#'):
                    card_line_no -= 1
                    continue
                if card_line_no == 2:
                    ori = line[4:].strip()
                    shutil.copy((act_inc_dir + '//' + ori), workDir, follow_symlinks=True)

        # include counter
        i += 1

    if not 'ori' in locals():
        ori = False

    return MUPIFprop,ori

def meshgen(nodes,cells, element_type):
    """
    Generates a mesh from connectivity and node array
    """

    vertexlist=[];
    celllist=[];
    num = 0
    mesh = Mesh.UnstructuredMesh();

    # generate vertices
    for node in nodes:
#        if (debug):
#            print ("Adding vertex %d: %f %f %f " % (node_ids[ii], nodes_coord[ii,0], nodes_coord[ii,1], nodes_coord[ii,2]))
        vertexlist.append(Vertex.Vertex(nodes[node]['MUPIF_NUM'], nodes[node]['MUPIF_NUM'], nodes[node]['coord']))
        num +=1

    # generate cells
    for cell in cells:
        if (element_type == 'QUAD4'):
#            if (debug):
#                print ("Adding quad %d: %d %d %d %d" % (elem_ids[ii], elem_connect[ii,0], elem_connect[ii,1], elem_connect[ii,2], elem_connect[ii,3]))
            celllist.append(Cell.Quad_2d_lin(mesh, cells[cell]['MUPIF_NUM'], cells[cell]['MUPIF_NUM'], cells[cell]['vertices']))
        else:
            raise APIError.APIError('Unsupported element type for the mesh creation!')

    mesh.setup (vertexlist, celllist);
    return mesh

def read_ori(ori_fpath,model):

    orifile = open(ori_fpath, 'r')
    orilines = orifile.readlines()

    act_keyword = ''
    for line in orilines:
        if line.startswith('$') or line.startswith('#'):
            continue
        if re_newkey.match(line):
            card_line_no = 1
            act_keyword = line[:7]
        else:
            card_line_no += 1
        if act_keyword == 'SOLID /':
            linesplit = line[7:].split()
            elid = convert_or_set_zero(linesplit[0], 'int')
            ort1 = convert_or_set_zero(linesplit[1], 'float')
            ort2 = convert_or_set_zero(linesplit[2], 'float')
            ort3 = convert_or_set_zero(linesplit[3], 'float')
            ort4 = convert_or_set_zero(linesplit[4], 'float')
            ort5 = convert_or_set_zero(linesplit[5], 'float')
            ort6 = convert_or_set_zero(linesplit[6], 'float')
            if model['elements'].has_key(elid):
               model['elements'][elid]['ortho'] = [ort1, ort2, ort3, ort4, ort5, ort6]


    return True


def readERF(erf_fpath,OUTPUTtype,OUTPUTid,OUTPUTcomp,time):

    if not os.path.exists(erf_fpath):
        error_msg = ('ERROR: Result file : %s does not exist!\n' % (erf_fpath))
        raise  APIError.APIError(error_msg)

    erffile = h5py.File(erf_fpath, 'r')

    try:
        if OUTPUTtype == 'THNOD':

            title = erffile['/post/constant/attributes/NODE/erfblock/title']
            for ii in range(title.maxshape[0]):
                if OUTPUTid in str(title[ii]):
                    THNODid = erffile['/post/constant/attributes/NODE/erfblock/entid'][ii]

            entid = erffile['/post/multistate/TIMESERIES1/multientityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/entid']
            indexval =  np.array(erffile['/post/multistate/TIMESERIES1/multientityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/indexval'][:,0])
            for ii in range(entid.maxshape[0]):
                if entid[ii] == THNODid:
                    res = np.array(erffile['/post/multistate/TIMESERIES1/multientityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/res'][:,ii,:])

        elif OUTPUTtype == 'SECFO':
            entid = erffile['/post/multistate/TIMESERIES1/multientityresults/SECTION/Section_Force/ZONE1_set1/erfblock/entid']
            indexval =  np.array(erffile['/post/multistate/TIMESERIES1/multientityresults/SECTION/Section_Force/ZONE1_set1/erfblock/indexval'][:,0])
            for ii in range(entid.maxshape[0]):
                if entid[ii] == OUTPUTid:
                    res = np.array(erffile['/post/multistate/TIMESERIES1/multientityresults/SECTION/Section_Force/ZONE1_set1/erfblock/res'][:,ii,:])

        elif OUTPUTtype == 'BUCKL':
            value = []
            for mode in erffile['/post/MODAL_BASIS/'].keys():
                value.append(erffile['/post/MODAL_BASIS/' + mode + '/EIGEN_VALUE/erfblock/res'][0,0])


            result_time = None

        if OUTPUTtype == 'THNOD' or OUTPUTtype == 'SECFO':
            if not OUTPUTcomp == 'ALL':
                res = res[:,OUTPUTcomp-1]
            mintime = indexval[0]
            maxtime = indexval[-1]
            if time == -1:
                result_time = indexval
                value  = res
            else:
                if time<mintime or time>maxtime:
                    print('ERROR: Result time is not in the model time frame!')
                    raise
                result_time = time
                if time == mintime:
                    if res.ndim == 1:
                        value  = res[0]
                    else:
                        value  = res[0,:]
                elif time == maxtime:
                    if res.ndim == 1:
                        value  = res[-1]
                    else:
                        value  = res[-1,:]
                else:
                    if res.ndim == 1:
                        value = np.interp(time,indexval,res)
                    else:
                        value = np.zeros(3)
                        for ii in range(3):
                            value[ii] = np.interp(time,indexval,res[:,ii])
        erffile.close()
    except:
        erffile.close()
        raise APIError.APIError('ERROR: Reading result file.')

    return value,result_time

def readERFfield(erf_fpath,field_type,time):

    if not os.path.exists(erf_fpath):
        error_msg = ('ERROR: Result file : %s does not exist!\n' % (erf_fpath))
        raise  APIError.APIError(error_msg)

    erffile = h5py.File(erf_fpath, 'r')

#    try:
    if field_type == 'DISPLACEMENT':

        shell_connect = np.array(erffile['/post/constant/connectivities/SHELL/erfblock/ic'])
        shell_ids     = np.array(erffile['/post/constant/connectivities/SHELL/erfblock/idele'])

        node_ids   = np.array(erffile['/post/constant/entityresults/NODE/COORDINATE/ZONE1_set0/erfblock/entid'])
        node_coord = np.array(erffile['/post/constant/entityresults/NODE/COORDINATE/ZONE1_set0/erfblock/res'])

        nodes = {}
        num = 0
        for ii in range(node_ids.shape[0]):
            nodes[node_ids[ii]] = dict(MUPIF_NUM = num, coord = (node_coord[ii,0], node_coord[ii,1], node_coord[ii,2]))
            num += 1

        cells = {}
        num = 0
        for ii in range(shell_ids.shape[0]):
            vertices=(nodes[shell_connect[ii,0]]['MUPIF_NUM'],\
                      nodes[shell_connect[ii,1]]['MUPIF_NUM'],\
                      nodes[shell_connect[ii,2]]['MUPIF_NUM'],\
                      nodes[shell_connect[ii,3]]['MUPIF_NUM'])

            cells[shell_ids[ii]] = dict(MUPIF_NUM = num,vertices = vertices)
            num += 1

        mesh = meshgen(nodes,cells, 'QUAD4')

        states = erffile['/post/singlestate']
        for state in states.keys():

            state_time = erffile['/post/singlestate/'+state + '/entityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/indexval'][0]
            if state_time == time.inUnitsOf('ms').value:
                entid = np.array(erffile['/post/singlestate/'+state + '/entityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/entid'])
                res   = np.array(erffile['/post/singlestate/'+state + '/entityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/res'])

                resultField = Field.Field(mesh,FieldID.FID_ESI_VPS_Displacement,ValueType.Vector,'mm',time,None,1)

                for ii in range(entid.shape[0]):
                    resultField.setValue(nodes[entid[ii]]['MUPIF_NUM'], (res[ii,0],res[ii,1],res[ii,2]))
                break
            elif state_time > time.inUnitsOf('ms').value:
                entid    = np.array(erffile['/post/singlestate/'+ state + '/entityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/entid'])
                next_res = np.array(erffile['/post/singlestate/'+ state + '/entityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/res'])
                prev_res = np.array(erffile['/post/singlestate/'+ prev_state + '/entityresults/NODE/Translational_Displacement/ZONE1_set1/erfblock/res'])

                res = np.zeros_like(next_res)

                res = prev_res*(state_time-time)/(state_time-prev_state_time) + next_res*(time-prev_state_time)/(state_time-prev_state_time)

                break
            else:
                prev_state      = state
                prev_state_time = state_time

        if state_time < time.inUnitsOf('ms').value:
            print('ERROR: Result time is not in the model time frame!')
            raise

    erffile.close()
#    except:
#        erffile.close()
#        raise APIError.APIError('ERROR: Reading result file.')

    return resultField



def float_max_prec(floatnum, num_fields):
    """
    :param floatnum: floating point number
    :param num_fields: number of fields available for writing
    :return: string formatted as float with maximum possible precision
    """
    num_leading_digits = len(str(abs(floatnum)).split('.')[0])
    format_string = '{:%i.%if}' % (num_fields, num_fields - num_leading_digits - 3)  #
    float_formatted = format_string.format(floatnum)
    # REMOVE TRAILING ZEROS
    while float_formatted[-1] == '0':
        float_formatted = ' ' + float_formatted[:-1]
    return float_formatted


def convert_or_set_zero(entry, target_type):
    """
    Convert entry to target type, or set to zero, if entry is empty
    :param entry: entry string from card line
    :param target_type: target type to convert to, either int or float
    :return: converted entry or zero
    """
    entry = entry.strip()
    if entry:
        if target_type == 'int':
            entry = int(entry)
        elif target_type == 'float':
            entry = float(entry)
    else:
        entry = 0
    return entry


def getMUPIFpropID(ID):


    # ESI VPS properties
    ID_DB = {}
    ID_DB[90001] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_TEND     ,unit = 'ms'      ,time = None)
    ID_DB[90002] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_E0t1,unit = 'kN/mm**2',time = None)
    ID_DB[90003] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_E0t2,unit = 'kN/mm**2',time = None)
    ID_DB[90004] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_E0t3,unit = 'kN/mm**2',time = None)
    ID_DB[90005] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_G012,unit = 'kN/mm**2',time = None)
    ID_DB[90006] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_G023,unit = 'kN/mm**2',time = None)
    ID_DB[90007] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_G013,unit = 'kN/mm**2',time = None)
    ID_DB[90008] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_NU12,unit = None,      time = None)
    ID_DB[90009] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_NU23,unit = None,      time = None)
    ID_DB[90010] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_NU13,unit = None,      time = None)
    ID_DB[90011] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_E0c1,unit = 'kN/mm**2',time = None)
    ID_DB[90012] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_RHO ,unit = 'kg/mm**3',time = None)
    ID_DB[90013] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_hPLY     ,unit = 'mm'      ,time = None)
    ID_DB[90014] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_XT  ,unit = 'kN/mm**2',time = None)
    ID_DB[90015] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_XC  ,unit = 'kN/mm**2',time = None)
    ID_DB[90016] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_YT  ,unit = 'kN/mm**2',time = None)
    ID_DB[90017] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_YC  ,unit = 'kN/mm**2',time = None)
    ID_DB[90018] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_PLY1_S12 ,unit = 'kN/mm**2',time = None)
    ID_DB[90019] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_MOMENT   ,unit = 'kNmm'    ,time = None)
    ID_DB[90020] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_ROTATION ,unit = None      ,time = None)

    ID_DB[91001] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_FIRST_FAILURE_MOM  ,unit = 'kNmm',time = None)
    ID_DB[91002] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_FIRST_FAILURE_ROT  ,unit = None  ,time = None)
    ID_DB[91003] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_FIRST_FAILURE_LOC  ,unit = 'mm'  ,time = None)
    ID_DB[91004] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_FIRST_FAILURE_ELE  ,unit = None  ,time = None)
    ID_DB[91005] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_FIRST_FAILURE_PLY  ,unit = None  ,time = None)
    ID_DB[91006] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_FIRST_FAILURE_PART ,unit = None  ,time = None)
    ID_DB[91007] = dict(MUPIF_PID = PropertyID.PID_ESI_VPS_BUCKL_LOAD         ,unit = 'kN'  ,time = None)

    if ID in ID_DB:
        return ID_DB[ID]
    else:
        print('WARNING: PropertyID : %d is not defined in the database and will be neclected.' % ID )
        return False




def tail(filepath):
    with open(filepath, "rb") as f:
        first = f.readline()      # Read the first line.
        f.seek(-2, 2)             # Jump to the second last byte.
        while f.read(1) != b"\n": # Until EOL is found...
            try:
                f.seek(-2, 1)     # ...jump back the read byte plus one more.
            except IOError:
                f.seek(-1, 1)
                if f.tell() == 0:
                    break
        last = f.readline()       # Read last line.
    return last