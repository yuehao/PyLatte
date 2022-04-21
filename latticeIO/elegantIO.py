from .fileIO import LatticeFile

import ply.lex as lex
import ply.yacc as yacc


import numpy.random as rand
import numpy as np
import copy
import re
import subprocess as subp
import shlex
import io
import sdds
import os


class elegantLexer:
    literals = "=,-*:&()"
    # reserved={r'[Ll][Ii][Nn][Ee]':'KWLINE'}
    reserved = {'LINE': 'KWLINE', 'USE': 'KWUSE', 'RETURN': 'KWRETURN'}
    tokens = [ 'NUMBER', 'STRING', 'RET'] + list(reserved.values())



    def t_NUMBER(self, t):
        r'([0-9]*\.?[0-9]+|[0-9]+\.[0-9]*)(([eE][-+]?[0-9]+)?)'
        t.value = float(t.value)
        return t

    t_ignore = ' \t'

    def t_STRING(self, t):
        r'[-%]?[A-Za-z][A-Za-z0-9_\.]*'
        t.type = elegantLexer.reserved.get(t.value.upper(), 'STRING')
        return t

    def t_QUOTE(self, t):
        r'\"'
        pass

    def t_COMMENT(self, t):
        r'(\#|!).*'
        pass

    def t_contline(self, t):
        r'\&[ \t]*\n'
        # print('foundit')
        t.lexer.lineno += 1
        self.continuetor = 1

    def t_RET(self, t):
        r'\n+'
        t.lexer.lineno += len(t.value)
        t.value = 'RET'
        return t

    def t_EOF(self, t):
        r'<<EOF>>'
        t.value = 'EOF'
        return t

    def t_error(self, t):
        print("Illegal character {}".format(t.value[0]));
        t.lexer.skip(1)

    def __init__(self, **kwargs):
        self.lexer = lex.lex(module=self, optimize=1, reflags=re.IGNORECASE, **kwargs)
        self.continuetor = 0

    def lexit(self, data):
        self.lexer.input(data)
        while True:
            tok = self.lexer.token()
            if not tok: break
            # print tok


class elegantParser:
    # starting = "statement"
    tokens = elegantLexer.tokens

    def __init__(self, **kwargs):
        self.dicttemp = {}
        self.listtemp = []
        self.dicts = {}
        self.thisline = {}
        self.useline = ''
        self.parser = yacc.yacc(debug=0, module=self, **kwargs)

    def p_statements(self, p):
        '''statements : statement
                      | statements statement'''
        if self.thisline != {}:
            self.dicts[self.thisline['NAME']] = self.thisline
        self.thisline = {}

    def p_statement(self, p):
        '''statement : STRING ':' STRING properties RET
                     | STRING ':' KWLINE '=' '(' bline ')' RET
                     | KWUSE ',' STRING RET
                     | RET
                     | KWRETURN'''

        if (len(p) == 6):
            self.thisline['NAME'] = p[1]
            self.thisline['TYPE'] = p[3]
            self.thisline.update(self.dicttemp)
            dicttemp = {}
        elif (len(p) == 9):
            self.thisline['NAME'] = p[1]
            self.thisline['TYPE'] = 'BEAMLINE'
            self.thisline['SEQUENCE'] = self.listtemp
            listtemp = []
        elif (len(p) == 4):
            self.useline = p[3]

    def p_bline(self, p):
        '''bline : STRING
                | NUMBER '*' STRING
                | bline ',' STRING
                | bline ',' NUMBER '*' STRING
                | bline ',' '(' bline ')'
                '''
        if len(p) == 2:
            self.listtemp.append(p[1])
        elif len(p) == 4:
            self.listtemp.append(p[3])
        elif len(p) == 6 and p[3]=='(':
            pass
        elif len(p) == 6:
            for _ in range(int(p[3])):
                self.listtemp.append(p[5])

    def p_properties(self, p):
        '''properties :
                      | properties property'''
        pass

    def p_property(self, p):
        '''property : ',' STRING '=' values
                    | ',' STRING '=' '\"' values '\"'
        '''
        if len(p) == 5:
            self.dicttemp[p[2]] = p[4]
        elif len(p) == 7:
            self.dicttemp[p[2]] = p[6]

    def p_values(self, p):
        '''values : NUMBER
                  | STRING
                  | '-' NUMBER'''
        if (len(p) == 2):
            p[0] = p[1]
        elif len(p) == 3:
            p[0] = -1.0 * p[2]

    def p_error(self, p):
        raise Exception(p)


class elegantLatticeFile(LatticeFile):
    elegantTypes=[]
    elegantParameterNames={}
    elegantTypesRev={}
    elegantTypesAsymmetry=[]
    typeDict={}
    propertyDict={}
    elegant_dir=''
    sdds_dir=''

    def __init__(self, filename='', readfrom='', elegant_dir_input='', sdds_dir_input=''):
        LatticeFile.__init__(self)
        self.lattice_format='Elegant'
        self.filename=filename
        if len(elegantLatticeFile.elegantTypes)==0:
            self.resetIODict(elegant_dir_input, sdds_dir_input)
        self.useline= ''
        if readfrom!='':
            self.parseFrom(readfrom)
        self.codeName='ELEGENT'.upper()

    def resetIODict(self, elegant_dir_input='', sdds_dir_input=''):
        #status,_=subp.
        printdict = elegantCommandFile('tempdict.ele')
        printdict.addCommand('print_dictionary', note='', filename='tempdict.dat', SDDS_form=1)
        printdict.write()

        cmdstr = 'elegant tempdict.ele'
        subp.call(shlex.split(cmdstr), stdout=subp.PIPE)


        '''
        cmdstr = 'sddsprintout tempdict.dat -param=ElementType'
        p1 = subp.Popen(shlex.split(cmdstr), stdout=subp.PIPE)
        out = p1.communicate()[0]

        outfile = io.BytesIO(out)

        typeind = 1
        for line_b in outfile:
            line=line_b.decode()
            if 'ElementType' in line:
                elenametemp = line.split()[-1]
                elegantLatticeFile.elegantTypes[elenametemp] = None
                elegantLatticeFile.elegantParameterNames[elenametemp] = {}
                cmdstr = 'sddsprintout tempdict.dat -fromPage={ind} -toPage={ind} -col=ParameterName -col=Default'.format(
                    ind=typeind)
                p2 = subp.Popen(shlex.split(cmdstr), stdout=subp.PIPE)
                out2 = p2.communicate()[0]

                typeind += 1
                outfile2 = io.BytesIO(out2)
                flag = 0
                for subline_byte in outfile2:
                    subline = subline_byte.decode()
                    if flag == 1:
                        temp = subline.split()
                        if temp[0] in ['E1', 'H1', 'END1_FOCUS']:  # These parameter break the symmetry
                            elegantLatticeFile.elegantTypesAsymmetry.append(elenametemp)
                        if len(temp) == 1:
                            elegantLatticeFile.elegantParameterNames[elenametemp][temp[0]] = None
                        else:
                            elegantLatticeFile.elegantParameterNames[elenametemp][temp[0]] = temp[1]
                    elif '---' in subline:
                        flag = 1
        '''

        try:
            f = open("tempdict.dat")
        except FileNotFoundError :
            print("File not accessible")
            return
        finally:
            f.close()

        eledict=sdds.SDDS(0)
        eledict.load("tempdict.dat")
        elegantLatticeFile.elegantTypes=eledict.parameterData[0]
        
        for i in range(len(elegantLatticeFile.elegantTypes)):
            eletype=elegantLatticeFile.elegantTypes[i]
            paramdict={}
            for j in range(len(eledict.columnData[0][i])):
                param=eledict.columnData[0][i][j]
                dt=eledict.columnData[2][i][j].upper()
                defval=eledict.columnData[3][i][j]
                if dt == 'STRING':
                    pass
                elif dt == 'DOUBLE':
                    defval=float(defval)
                else:
                    defval=int(defval)
                paramdict[param]=defval    
            elegantLatticeFile.elegantParameterNames[eletype]=paramdict

            



        cmdstr = 'rm -f tempdict.ele tempdict.dat'
        subp.call(shlex.split(cmdstr))
        # for et in elements.element.elementTypes:
        #    print(et)
        #    en=getattr(elements, et.lower()).toPrograms['ELEGANT'].upper()
        #    if en in elegantLatticeFile.elegantTypes:
        #        elegantLatticeFile.elegantTypes[en]=et
        #        elegantLatticeFile.elegantTypesRev[et]=en
        #    else:
        #        print('Error, Elegant doesnot recognized the element type {} from defined class {}'.format(en, et))
        #        exit()


        # revd=dict([reversed(i) for i in elegantLatticeFile.pyerl_Types.items()])
        # elegantLatticeFile.pyerl_Types.update(revd)
        # print(elegantLatticeFile.elegantTypesAsymmetry)
        # print(elegantLatticeFile.elegantParameterNames['CSBEND'])
        # elementTypes=['DRIFT', 'DIPOLE', 'QUAD', 'SEXTRUPOLE ,OCTUPOLE', 'MULTIPOLE', 'MATRIX', 'TPSMAP', 'BEAMBEAM', 'RFCAVITY', 'MARKER', 'BPM', 'SOLENOID', 'KICKER', 'CENTER', 'MALIGN']
        #elegantElementType = ['EDRIFT', 'CSBEND', 'KQUAD', 'KSEXT', 'OCTU', 'MULT', 'EMATRIX', 'NULL', 'NULL', 'RFCW',
        #                      'MARK', 'MONI', 'SOLE', 'KICKER', 'CENTER', 'MALIGN']
        #elegantLatticeFile.typeDict = dict(zip(elements.Element.elementTypes, elegantElementType))
        # commonPropertyNames = ['NAME', 'L', 'TYPE', 'RMS_DISPLACEMENT', 'RMS_ANGLE_DISPLACEMENT', 'TILT', 'NSTEP', 'APERTURE', 'DISPLACEMENT', 'ANGLE_DISPLACEMENT','GROUP']
        #elegantcommonPropertyName = ['NAME', 'L', 'TYPE', 'NULL', 'NULL', 'TILT', 'N_KICKS', 'NULL', 'DX_DY_DZ_SPLIT',
        #                             'NULL', 'GROUP']
        #elegantLatticeFile.propertyDict = dict(zip(elements.Element.propertyNames, elegantcommonPropertyName))
        propertyNames = [
            'BN', 'AN', 'BERROR', 'AERROR',
            'DESIGN_ORDER', 'MAX_ORDER', 'ANGLE', 'E1', 'E2',
            'K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'K7', 'K8',
            'K1S', 'K2S', 'K3S', 'K4S', 'K5S', 'K6S', 'K7S', 'K8S',  # MPOLE
            'VOLT', 'FREQ', 'PHASE', 'PHASE_REF',  # RFCAVITY
            'XC', 'PXC', 'YC', 'PYC', 'DEC', 'CTC',  # CENTER
            'DX', 'DPX', 'DY', 'DPY', 'DEOE', 'DCT',  # MALIGN
            'WEIGHT',  # WEIGHT
        ]
        elegantPropertyName = [
            'NULL', 'NULL', 'NULL', 'NULL',
            'NULL', 'NULL', 'ANGLE', 'E1', 'E2',
            'K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'K7', 'K8',
            'K1S', 'K2S', 'K3S', 'K4S', 'K5S', 'K6S', 'K7S', 'K8S',
            'VOLT', 'FREQ', 'PHASE', 'PHASE_REF',
            'XC', 'PXC', 'YC', 'PYC', 'DEC', 'CTC',
            'DX', 'DPX', 'DY', 'DPY', 'DEOE', 'DCT',
            'WEIGHT',

        ]
        elegantLatticeFile.propertyDict.update(dict(zip(propertyNames, elegantPropertyName)))

    def applyErrorToAll(self, typeName, parameterName, rmsValue, namecontain='', mode='ABSOLUTE', theseed=0):
        '''
        :param typeName: Apply to all elements of this type
        :param parameterName: parameter to be changed
        :param rmsValue: the rms value
        :param mode: If Absolute the value will be addon, if not the percentage will be changed.
        :return: return the (list_of_name, list_of_value)
        '''
        ename_of_type = []
        assigned_error = []
        if theseed > 0:
            rand.seed(theseed)
        for ele in self.elementList:
            if self.getElementRootType(ele['NAME']) == typeName.upper() and namecontain in ele['NAME']:
                ename_of_type.append(ele['NAME'])
                rdn = rand.randn() * rmsValue
                assigned_error.append(rdn)
                if parameterName not in ele:
                    ele[parameterName] = rdn
                else:
                    if mode == 'ABSOLUTE':
                        ele[parameterName] += rdn
                    else:
                        ele[parameterName] *= (1 + rdn)

        return ename_of_type, assigned_error



    def checkType(self, typename, parameterName=None):
        if typename.upper() in elegantLatticeFile.elegantTypes:
            if parameterName is None:
                return True
            if (parameterName.upper() in elegantLatticeFile.elegantParameterNames[typename.upper()]):
                return True
        else:
            if typename.upper() in self.elementNameDict:
                self.checkType(self.getElementRootType(typename.upper()),parameterName)
        return False



    def isDrift(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'DRIF' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False

    def isDipole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'BEND' in parent_type:
            temp=self.getElementProperties(ele_name)
            if 'K1' not in temp:
                temp['K1']=0
            return temp
        else:
            return False

    def isQuadrupole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'QUAD' in parent_type:
            temp=self.getElementProperties(ele_name)
            if 'K1' not in temp:
                temp['K1']=0
            return temp
        else:
            return False
    def isSextrupole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'SEXT' in parent_type:
            temp=self.getElementProperties(ele_name)
            if 'K2' not in temp:
                temp['K2']=0
            return temp
        else:
            return False

    def isSolenoid(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'SOLE' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False
    def isCavity(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'RFC' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False

    '''def plotBeamline(self, plt_axis, beamline_name, colors=('DarkOrchid', 'Maroon', 'DeepSkyBlue', 'ForestGreen'),
                     heights=(0.7, 1, 0.8, 0.4), s_start=0):

        :param plt_axis: matplotlib axis variable
        :param beamline_name: name of beamline to be plotted.
        :param colors: The color for LINAC, dipole, Quad and multipoles
        :return:

        self.setUseLine()
        bl_pos, bl_list = self.elementPosInUseLine, self.useLineList
        print(beamline_name, len(bl_pos), len(bl_list))
        import matplotlib.pyplot as mpl
        from matplotlib.patches import Rectangle
        from matplotlib.patches import FancyBboxPatch

        plt_axis.set_xlim([s_start, s_start + bl_pos[-1]])

        plt_axis.set_ylim([-1, 1])
        plt_axis.set_yticks([])
        # plt_axis.xaxis.set_ticks_position('none')
        plt_axis.axhline(0, color='black')

        for i in range(len(bl_list)):
            start = bl_pos[i] + s_start
            ele = bl_list[i]
            indlist = self.getParentElements(ele)
            l_ele = bl_pos[i + 1] - bl_pos[i]
            lasttype = self.elementList[indlist[-1]]['TYPE']
            if l_ele > 0:
                if lasttype == 'RFCA' or lasttype == 'RFCW':
                    plt_axis.add_patch(
                        FancyBboxPatch((start, -heights[0] / 2.0), l_ele, heights[0], ec=colors[0], fc=colors[0]))
                elif 'BEND' in lasttype:
                    shift = 0
                    for ind in indlist:
                        if 'K1' in self.elementList[ind]:
                            shift = heights[1] * 0.2 * np.sign(self.elementList[ind]['K1'])
                            break
                    plt_axis.add_patch(
                        Rectangle((start, -heights[1] / 2.0 + shift), l_ele, heights[1], angle=0.0, ec=colors[1],
                                  fc='none'))

                elif 'QUAD' in lasttype:
                    shift = 0
                    for ind in indlist:
                        if 'K1' in self.elementList[ind]:
                            shift = heights[2] * 0.5 * np.sign(self.elementList[ind]['K1'])
                            break
                    plt_axis.add_patch(
                        Rectangle((start, -heights[2] / 2.0 + shift), l_ele, heights[2], angle=0.0, ec=colors[2],
                                  fc='none'))

                elif 'DRIF' not in lasttype:
                    plt_axis.add_patch(
                        Rectangle((start, -heights[3] / 2.0), l_ele, heights[3], angle=0.0, ec=colors[3], fc='none'))
'''
    def parseFrom(self, inputfile):

        # i=0
        bufferstr = ''
        for line in inputfile:
            # i+=1
            # print(i)
            # print(line)
            m = elegantLexer()
            m.lexit(line)
            bufferstr = ''.join([bufferstr.upper(), line])
            if m.continuetor == 1:
                continue
            else:
                p = elegantParser()
                p.parser.parse(bufferstr, lexer=m.lexer)
                bufferstr = ''
                for k, v in p.dicts.items():
                    if v['TYPE'] == 'BEAMLINE':
                        vline = v['SEQUENCE']

                        self.appendToBeamline(k, *vline)

                    else:
                        if self.checkType(v['TYPE']) or self.getElementIndex(v['TYPE']) is not None:

                            self.addElement(v['NAME'], v['TYPE'], **v)

                        else:
                            print ("The element {} with type {} is not recognized".format(v['NAME'], v['TYPE']))

        if p.useline != '':
            self.useline = p.useline
        else:
            self.useline = self.beamlineList[-1]['NAME']

    def downConvert(self, thelattice):
        '''
        Convert from general lattice file to elegant lattice file:
        '''

        self.__init__()

        for ele in thelattice.elementList:
            self.addElement(ele.properties['NAME'], elegantLatticeFile.typeDict[ele.properties['TYPE']])
            for k, v in ele.properties.items():
                if k in elegantLatticeFile.propertyDict:
                    kconv = elegantLatticeFile.propertyDict[k]
                    if kconv == 'NULL':
                        continue
                else:
                    continue

                if kconv[-6:] == '_SPLIT':
                    keys = kconv.split('_')
                    if len(keys) == len(v) + 1:
                        for i in range(len(v)):
                            thekey = keys[i]
                            self.addElement(ele.properties['NAME'], elegantLatticeFile.typeDict[ele.properties['TYPE']],
                                            thekey=v[i])
                else:
                    self.addElement(ele.properties['NAME'], elegantLatticeFile.typeDict[ele.properties['TYPE']],
                                    **{kconv: v})

        for lin in thelattice.beamlineName:
            self.appendToBeamline(lin, *thelattice.beamlineNameDict[lin])

    def toPyLatte(self):
        '''
        Convert to PyLatte lattice from elegant lattice file:
        '''
        lattice=[]
        the_line=self.beamlineList[self.beamlineNameDict[self.useline]]['LINE']
        for elename in the_line:
            test=self.isDrift(elename)
            if test is not False:
                from .. elements.drift import Drift
                lattice.append(Drift(elename, **test))

            test = self.isDipole(elename)
            if test is not False:
                from ..elements.dipole import Dipole
                lattice.append(Dipole(elename, **test))

            test = self.isQuadrupole(elename)
            if test is not False:
                from ..elements.quadrupole import Quadrupole
                lattice.append(Quadrupole(elename, **test))

            test = self.isSextrupole(elename)
            if test is not False:
                from ..elements.sextrupole import Sextrupole
                lattice.append(Sextrupole(elename, **test))

        return lattice



    def write(self, outputfilename='', mode='w'):
        if outputfilename != '':
            self.filename = outputfilename
        outfile = open(self.filename, mode)
        # write elements
        for i in range(len(self.elementList)):
            counter = outfile.tell()
            outfile.write('{}:{}'.format(self.elementList[i]['NAME'], self.elementList[i]['TYPE']))

            for k, v in self.elementList[i].items():

                if k != 'NAME' and k != 'TYPE' and '__USAGE_DICT' not in k and '__ID_IN_USELINE' not in k:
                    if outfile.tell() - counter > 50:
                        outfile.write('&\n\t')
                        counter = outfile.tell()
                    outfile.write(', {} = {}'.format(k, v))

            outfile.write('\n')

        # write lines
        for i in range(len(self.beamlineList)):
            counter = outfile.tell()

            outfile.write('{}:LINE=('.format(self.beamlineList[i]['NAME']))

            for ei in range(len(self.beamlineList[i]['LINE']) - 1):
                if outfile.tell() - counter > 50:
                    outfile.write('&\n\t')
                    counter = outfile.tell()
                outfile.write('{}, '.format(self.beamlineList[i]['LINE'][ei]))

            if len(self.beamlineList[i]['LINE']) > 0:
                outfile.write('{})\n'.format(self.beamlineList[i]['LINE'][-1]))
            else:
                outfile.write(')\n')


class elegantCommandFile:
    commandlist = (
    'alter_elements', 'amplification_factors', 'analyze_map', 'aperture_data', 'bunched_beam', 'change_particle',
    'chromaticity', 'closed_orbit', 'correct', 'correction_matrix_output', 'correct_tunes', 'coupled_twiss_output',
    'divide_elements', 'error_element', 'error_control', 'find_aperture', 'floor_coordinates', 'frequency_map',
    'global_settings', 'insert_elements ,insert_sceffects', 'linear_chromatic_tracking_setup', 'link_control',
    'link_elements', 'load_parameters', 'matrix_output', 'modulate_elements', 'moments_output', 'momentum_aperture',
    'optimize', 'optimization_constraint', 'optimization_covariable', 'optimization_setup',
    'parallel_optimization_setup',
    'optimization_term', 'optimization_variable', 'print_dictionary', 'ramp_elements', 'rf_setup', 'replace_elements',
    'rpn_expression', 'rpn_load', 'run_control', 'run_setup', 'sasefel', 'save_lattice', 'sdds_beam', 'semaphores',
    'slice_analysis',
    'subprocess', 'steering_element', 'touschek_scatter', 'transmute_elements', 'twiss_analysis', 'twiss_output',
    'track',
    'tune_shift_with_amplitude', 'vary_element')

    def __init__(self, filename):
        self.commandlist = []
        self.filename = filename

    def checkCommand(self, typename):
        for tn in elegantCommandFile.commandlist:
            if tn == typename.lower():
                return True
        return False

    def addCommand(self, command, note='', **params):
        if self.checkCommand(command) == False:
            print('The command {} is not recognized.'.format(command))
            return
        thiscom = {}
        thiscom['NAME'] = command.lower()
        thiscom['NOTE'] = note

        for k, v in params.items():
            thiscom[k] = v
        self.commandlist.append(thiscom)

    def modifyCommand(self, commandname, mode='last', **params):
        indlist = []
        i = 0
        for command in self.commandlist:
            if command['NAME'] == commandname.lower():
                indlist.append(i)
            i += 1
        if len(indlist) == 1:

            for k, v in params.items():
                self.commandlist[indlist[0]][k] = v
            return
        elif len(indlist) == 0:
            print("No such command {} can be found".format(commandname))
        else:
            if mode == 'last':
                for k, v in params.items():
                    self.commandlist[indlist[-1]][k] = v
                return
            elif mode == 'first':
                for k, v in params.items():
                    self.commandlist[indlist[0]][k] = v
                return

            elif isinstance(mode, int):
                for k, v in params.items():
                    self.commandlist[indlist[mode]][k] = v
                return

            else:
                print("The mode is invalid")

    def repeatCommand(self, commandname, mode='last'):
        indlist = []
        i = 0
        for command in self.commandlist:
            if command['NAME'] == commandname.lower():
                indlist.append(i)
            i += 1
        if len(indlist) == 1:
            self.commandlist.append(self.commandlist[indlist[0]])
            return
        elif len(indlist) == 0:
            print("No such command {} can be found".format(commandname))
        else:
            if mode == 'last':
                self.commandlist.append(self.commandlist[indlist[-1]])
            elif mode == 'first':
                self.commandlist.append(self.commandlist[indlist[0]])
            elif isinstance(mode, int):
                self.commandlist.append(self.commandlist[indlist[mode]])
            else:
                print("The mode is invalid")




    def optimize(self, lattice, ):
        pass

    def write(self, outputfilename='', mode='w'):
        if outputfilename != '':
            self.filename = outputfilename
        outfile = open(self.filename, mode)
        for command in range(len(self.commandlist)):
            outfile.write('&{}\n'.format(self.commandlist[command]['NAME']))
            if self.commandlist[command]['NOTE'] != '':
                outfile.write('! {}\n'.format(self.commandlist[command]['NOTE']))
            for k, v in self.commandlist[command].items():
                if k != 'NAME' and k != 'NOTE':
                    if len(k) <= 19:
                        outfile.write('\t{:20s}= {},\n'.format(k, v))
                    elif len(k) <= 29:
                        outfile.write('\t{:30s}= {},\n'.format(k, v))
                    else:
                        outfile.write('\t{:40s}= {},\n'.format(k, v))
            outfile.write('&end\n\n')





def get_SDDS_column(SDDSfile, column_name=[], convert_to_float=True):
    s = sdds.SDDS(0)
    import os.path
    if not os.path.isfile(SDDSfile):
        print('file {} does not exist'.format(SDDSfile))
        exit(-1)
    try:
        s.load(SDDSfile)
    except:
        print("the file {} can not be loaded".format(SDDSfile))
        exit()

    get_list = []
    for cn in column_name:

        if cn in s.columnName:
            ind = s.columnName.index(cn)
        elif cn.isdigit():
            ind = (int)(cn)
            if ind >= len(cn.columnName):
                print("Invalid column index {}".format(ind))
                exit(0)

        elif cn.startswith('-') and cn[1:].isdigit():
            ind = (int)(cn[1:])
            if ind >= len(cn.columnName):
                print("Invalid column index {}".format(-1 * ind))
                exit(0)
            ind *= (-1)
        else:
            print("Invalid column_name {}".format(cn))
            exit(0)

        get_list.append(ind)

    len_col = len(s.columnData[0][0])
    if convert_to_float == True:
        res = np.empty((len(get_list), len_col))
    else:
        res = np.empty((len(get_list), len_col), dtype='S16')
    nrow = 0
    for ind in get_list:
        res[nrow] = np.array(s.columnData[ind][0])
        nrow += 1
    cnlist = np.array(s.columnName)
    cnunit = np.array(s.columnDefinition)
    if convert_to_float:
        return res.astype(float), cnlist[get_list], cnunit[get_list, 1]
    return res, cnlist[get_list], cnunit[get_list, 1]

def read_elegant_SDDS_content(filename, column_names=None):
    try:
        f = open(filename)
    except FileNotFoundError :
        print("File not accessible")
        return
    finally:
        f.close()
    import pandas as pd
    indcn=[]
    sddsfile = sdds.SDDS(0)
    sddsfile.load(filename)
    if column_names is None:
        column_names=sddsfile.columnName
        indcn=[i for i in range(len(sddsfile.columnName))]
    else:
        column_names=list(set(column_names) & set(sddsfile.columnName))
        for name in column_names:
            ind = sddsfile.columnName.index(name)
            indcn.append(ind)
    column_list = (np.array(sddsfile.columnData)[indcn, 0, :].astype(float)).transpose()

    sdds_parameters = {}
    for para_name, para_value in zip(sddsfile.parameterName, sddsfile.parameterData):
        sdds_parameters[para_name] = para_value[0]
    return pd.DataFrame(data=column_list, columns=column_names), sdds_parameters
    
    


def elegant_findtwiss(lattice, beamline_to_use=None, rootname='temp', matched=1,
                      initial_optics=[1, 0, 0, 0, 1, 0, 0, 0],
                      alternate_element={}, gamma0=1.0e4 / 0.511,
                      divide_element=0,
                      closed_orbit = False,
                      twiss_columns=['s', 'betax', 'alphax', 'psix', 'etax', 'etaxp',
                                     'betay', 'alphay', 'psiy', 'etay', 'etayp']
                      ):
    '''
    Find twiss parameter using elegant
    :param lattice: The elegantLatticeFile object to be used
    :param beamline_to_use: The beamline to be used.  Default is None, then lattice.useline is used
    :param rootname: The rootname of the output file.  Default is 'temp'.
    :param matched: 1: seek for matched solution, 0: Use initial_optics as initial condition
    :param initial_optics: Initial optics to start with, default is [1, 0, 0, 0, 1, 0, 0, 0]
    :param alternate_element: No use for now
    :param closed_orbit:  No use for now
    :param divide_element: Number of division of elements
    :param gamma0: The reference lorentz factor of the particle
    :param twiss_columns: The output columns of the twiss parameters, default is ['s', 'betax', 'alphax', 'psix', 'etax', 'etaxp',
                                     'betay', 'alphay', 'psiy', 'etay', 'etayp']
    :return: return tulip of twiss numpy array and dictionary of parameters:(twiss array, twiss parameters)
    '''


    if beamline_to_use is None:
        beamline_to_use = lattice.useline
    lattice.write('{}.lte'.format(rootname))
    cmd_file = elegantCommandFile('{}.ele'.format(rootname))
    cmd_file.addCommand('run_setup',
                        lattice='{}.lte'.format(rootname),
                        use_beamline=beamline_to_use,
                        rootname=rootname,
                        p_central=np.sqrt(np.power(gamma0, 2.0) - 1),
                        centroid='%s.cen',
                        default_order=3,
                        concat_order=3,
                        element_divisions=divide_element,
                        )
    if closed_orbit:
        cmd_file.addCommand('closed_orbit',
                        output='%s.clo',
                        closed_orbit_iterations=1500,
                        closed_orbit_accuracy=1e-12,
                        iteration_fraction=0.3

                        )
    cmd_file.addCommand('twiss_output',
                        matched=matched,
                        output_at_each_step=1,
                        filename='%s.twi',
                        radiation_integrals=1,
                        beta_x=initial_optics[0],
                        alpha_x=initial_optics[1],
                        eta_x=initial_optics[2],
                        etap_x=initial_optics[3],
                        beta_y=initial_optics[4],
                        alpha_y=initial_optics[5],
                        eta_y=initial_optics[6],
                        etap_y=initial_optics[7]
                        )
    cmd_file.addCommand('run_control')
    cmd_file.addCommand('bunched_beam')
    cmd_file.addCommand('track')
    cmd_file.write()
    cmdstr = 'elegant {}.ele'.format(rootname)
    with open(os.devnull, "w") as f:
        subp.call(shlex.split(cmdstr), stdout=f)

    twifile = sdds.SDDS(0)
    twifile.load('{}.twi'.format(rootname))
    inds = []
    for name in twiss_columns:
        ind = twifile.columnName.index(name)
        inds.append(ind)
    twiss_list = (np.array(twifile.columnData)[inds, 0, :].astype(float)).transpose()
    import pandas as pd
    twiss_parameter = {}
    for para_name, para_value in zip(twifile.parameterName, twifile.parameterData):
        twiss_parameter[para_name] = para_value[0]

    if closed_orbit:
        clofile = sdds.SDDS(0)
        clofile.load('{}.clo'.format(rootname))
        inds = []
        for name in ['x', 'xp', 'y', 'yp']:
            ind = clofile.columnName.index(name)
            inds.append(ind)
        clo_list = (np.array(clofile.columnData)[inds, 0, :].astype(float)).transpose()
        twiss_list=np.hstack((twiss_list, clo_list))
        twiss_columns=twiss_columns+['x', 'xp', 'y', 'yp']


    return pd.DataFrame(data=twiss_list, columns=twiss_columns), twiss_parameter


def elegant_track(lattice, beamline_to_use=None, Npar=1, rootname='temp',
                  initial_optics=[1, 0, 0, 0, 1, 0, 0, 0], emit_x=0.0, emit_nx=0.0, emit_y=0.0, emit_ny=0.0,
                  centroid_columns=['s', 'Cx', 'Cxp', 'Cy', 'Cyp', 'Cs', 'Cdelta', 'pCentral'],
                  sigma_columns=['s', 'Sx', 'Sxp', 'Sy', 'Syp', 'Ss', 'Sdelta', 'ex', 'ecx', 'ey', 'ecy'],
                  gamma0=1.0e4 / 0.511):
    '''

    :param lattice: The elegantLatticeFile object to be used
    :param beamline_to_use: The beamline to be used.  Default is None, then lattice.useline is used
    :param Npar: Number of particles used in tracking, default is 1
    :param rootname: rootname: The rootname of the output file.  Default is 'temp'.
    :param initial_optics: Initial optics to start with, default is [1, 0, 0, 0, 1, 0, 0, 0]
    :param emit_x: Initial geometric emittance of horizontal direction
    :param emit_nx: Initial normalized emittance of horizontal direction, ignored if emit_x is given
    :param emit_y: Initial geometric emittance of vertical direction
    :param emit_ny: Initial normalized emittance of vertical direction, ignored if emit_y is given
    :param centroid_columns: The output columns of centroids
    :param sigma_columns: The output columns of sigmas
    :param gamma0: The reference lorentz factor of the particle
    :return: return tulip of two numpy array:(centroid_array, sigma_array)
    '''
    if beamline_to_use is None:
        beamline_to_use = lattice.useline
    lattice.write('temp.lte')
    cmd_file = elegantCommandFile('{}.ele'.format(rootname))
    cmd_file.addCommand('run_setup',
                        lattice='{}.lte'.format(rootname),
                        use_beamline=beamline_to_use,
                        rootname=rootname,
                        p_central=np.sqrt(np.power(gamma0, 2.0) - 1),
                        centroid='%s.cen',
                        sigma='%s.sig',
                        default_order=3,
                        concat_order=3,
                        )

    cmd_file.addCommand('run_control')
    cmd_file.addCommand('bunched_beam',
                        n_particles_per_bunch=Npar,
                        emit_x=emit_x,
                        emit_nx=emit_nx,
                        beta_x=initial_optics[0],
                        alpha_x=initial_optics[1],
                        eta_x=initial_optics[2],
                        etap_x=initial_optics[3],
                        emit_y=emit_y,
                        emit_ny=emit_ny,
                        beta_y=initial_optics[4],
                        alpha_y=initial_optics[5],
                        eta_y=initial_optics[6],
                        etap_y=initial_optics[7],
                        )
    cmd_file.addCommand('track')
    cmd_file.write()
    cmdstr = 'elegant {}.ele'.format(rootname)
    with open(os.devnull, "w") as f:
        subp.call(shlex.split(cmdstr), stdout=f)

    cenfile = sdds.SDDS(0)
    cenfile.load('{}.cen'.format(rootname))

    inds = []
    for name in centroid_columns:
        ind = cenfile.columnName.index(name)
        inds.append(ind)
    centroid = np.array(cenfile.columnData)[inds, 0, :].astype(float)

    sigfile = sdds.SDDS(0)
    sigfile.load('{}.sig'.format(rootname))
    inds = []
    for name in sigma_columns:
        ind = sigfile.columnName.index(name)
        inds.append(ind)
    sigma = np.array(sigfile.columnData)[inds, 0, :].astype(float)

    return centroid, sigma



if __name__ == '__main__':
    # m=elegantLexer()

    # m.lexit("a:LINE=(a,b,&    \nc,d)\n")

    # p=elegantParser()
    # dataline='''m1:test\n
    # m2:quad,l=2,k1=1\n
    # a:line=(m1,m2,m4,m3)\n'''
    #   p.parser.parse(dataline,lexer=m.lexer)
    #  print(p.dicts)
    test = elegantLatticeFile('test1.lte')
    test.parseFrom('test.lte')
    # test1=test.duplicate('test1.lte','arc1_')
    # test1.write()
    # print(test.elements)
    # print(test.beamlines)

