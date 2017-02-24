import copy
import numpy as np


class LatticeFile(object):
    '''
    It can be reused by all supported external files
    '''
    def __init__(self):
        self.elementList = []
        self.elementNameDict = {}
        self.beamlineList = []
        self.beamlineNameDict = {}
        self.beamlineElementSet = {}
        self.useLineList = []
        self.elementPosInUseLine = None

        self.useline= ''

    def checkType(self, typename, parameterName=None):
        return True


    def parseFrom(self, filename):
        pass

    def getElementIndex(self, elename):
        elename = elename.upper()
        if elename[0] == '-':
            return self.getElementIndex(elename[1:])
        elif elename in self.elementNameDict:
            return self.elementNameDict[elename]
        else:
            return None


    def addElement(self, _name, _eletype, **params):
        roottype = ''
        #if self.checkType(_eletype) == False:
        #    if (self.getElementIndex(_eletype) is None):
        #        print ("The element {} with type {} is not recognized when adding this element".format(_name, _eletype))
        #        raise KeyError
        #    else:
        #        roottype = self.getElementRootType(_eletype)
        if self.getElementIndex(_eletype) is None:
            if self.checkType(_eletype) == False:
                print ("The element {} with type {} is not recognized when adding this element".format(_name, _eletype))
                raise KeyError
        else:
            roottype = self.getElementRootType(_eletype)


        thisele = {}
        name = _name.upper()
        eletype = _eletype.upper()
        thisele['NAME'] = name
        thisele['TYPE'] = eletype
        thisele['__USAGE_DICT'] = {}
        thisele['__ID_IN_USELINE'] = []

        for k, v in params.items():
            if k.upper() == 'NAME':
                pass
            elif k.upper() == 'TYPE':
                pass
            elif self.checkType(eletype, k) or self.checkType(roottype, k):
                thisele[k.upper()] = v
            else:
                print('Unrecognized parameter name {} in element {} when adding'.format(k, name))
                raise KeyError

        ind = self.getElementIndex(name)
        if ind is None:
            cur_len = len(self.elementList)
            self.elementNameDict[name] = cur_len
            self.elementList.append(thisele)
        else:
            self.elementList[ind] = thisele
            print ("Warning, the element {} is redefined when adding this element, overwriting with newer one.".format(name))
        return

    def getElementRootType(self, elename):
        elename = elename.upper()
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when getting root type'.format(elename))
            raise KeyError
        cur_type = self.elementList[ind]['TYPE']
        if cur_type in self.elementNameDict:
            return self.getElementRootType(cur_type)
        else:
            return cur_type


    def getElementProperties(self, elename, keyname=None, partial_key=''):
        elename = elename.upper()
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when getting properties'.format(elename))
            raise KeyError
        if keyname is None:
            root_type = self.getElementRootType(elename)
            cur_type = self.elementList[ind]['TYPE']
            if cur_type == root_type:
                return {k: v for k, v in self.elementList[ind].items()
                        if partial_key.upper() in k and k not in ['NAME', '__USAGE_DICT', '__ID_IN_USELINE']}
            else:
                return self.getElementProperties(cur_type, None, partial_key).update(
                    {k: v for k, v in self.elementList[ind].items()
                     if partial_key.upper() in k and k not in ['TYPE', 'NAME', '__USAGE_DICT', '__ID_IN_USELINE']})
        else:
            keyname = keyname.upper()
            if keyname in self.elementList[ind]:
                return self.elementList[ind][keyname]
            elif self.elementList[ind]['TYPE'] in self.elementNameDict:
                return self.getElementProperties(self.elementList[ind]['TYPE'], keyname)
            else:
                return None



    def getParentElements(self, elename):
        pe_list = []
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when finding parents'.format(elename))
            raise KeyError
        t = self.elementList[ind]['TYPE']
        pe_list.append(t)
        if t in self.elementNameDict:
            pe_list += self.getParentElements(t)
        return pe_list

    def compareElements(self, other_lattice, elename, other_name=None):
        elename = elename.upper()
        if other_name is None:
            other_name = elename
        else:
            other_name = other_name.upper()
        ind_other = other_lattice.getElementIndex(other_name)
        ind = self.getElementIndex(elename)
        if ind is None or ind_other is None:
            return False
        return self.getElementProperties(elename) == other_lattice.getElementProperties(other_name)

    def modifyElement(self, elename, increment=False, **params):
        elename = elename.upper()
        ind = self.getElementIndex(elename)
        if ind is None:
            print('Unrecognized element name {} when modifying'.format(elename))
            raise KeyError
        for k, v in params.items():
            if self.checkType(self.elementList[ind]['TYPE'], k):
                if k in self.elementList[ind] and increment:
                    self.elementList[ind][k.upper()] += v
                else:
                    self.elementList[ind][k.upper()] = v
            else:
                print('Unrecognized parameter name {} in element {}'.format(k, elename))
                raise KeyError
        return

    def modifyAllElements(self, eletype, increment=False, name_contain='', **params):
        eletype = eletype.upper()

        for ele in self.elementList:
            if ele['TYPE']==eletype and name_contain.upper() in ele['NAME']:
                self.modifyElement(ele['NAME'], increment, **params)

                # parameterName = parameterName.upper()
                # if eletype in self.elementNameDict:
                #    roottype = self.getElementRootType(eletype)
                # else:
                #   roottype = eletype
                # if self.checkType(eletype):
                #    roottype = eletype
                # else:
                #    roottype = self.getElementRootType(eletype)
                # if not self.checkType(roottype, parameterName):
                #    print('The parameter name {} is not good for element type {}'.format(parameterName, eletype))
                #    raise KeyError
                #current_value = self.getElementProperties(ele['NAME'], parameterName)
                #future_value = parameterValue
                #if increment and current_value is not None:
                #    future_value += current_value
                #if ele['TYPE'] == eletype:
                #    ele[parameterName] = future_value
                #elif eletype in self.getParentElements(ele['NAME']) and current_value != future_value:
                #    ele[parameterName] = future_value

        #for ele in self.elementList:
        #    if eletype in self.getParentElements(ele['NAME']) and name_contain in ele['NAME'] and ele['TYPE']!=eletype \
        #            and self.getElementProperties(ele['NAME'],parameterName) != parameterValue:
        #        ele[parameterName] = parameterValue
            #elif self.elementList[i]['TYPE'] in self.elementNameDict:
            #    if self.elementList[self.elementNameDict[self.elementList[i]['TYPE']]] == eletype:
            #        self.elementList[i][parameterName]=parameterValue

    def isDrift(self, ele, parent_type=None):
        pass
    def isDipole(self, ele, parent_type=None):
        pass
    def isQuadrupole(self, ele, parent_type=None):
        pass
    def isSolenoid(self, ele, parent_type=None):
        pass
    def isCavity(self, ele, parent_type=None):
        pass


    def plotBeamline(self, plt_axis, colors=('DarkOrchid', 'Maroon', 'DeepSkyBlue', 'ForestGreen'),
                     heights=(1.2, 0.8, 0.5, 0.5), s_start=0, other_components = None):
        '''
        :param plt_axis: matplotlib axis variable
        :param beamline_name: name of beamline to be plotted.
        :param colors: The color for shapes that represent cavity, dipole, Quad and solenoid
        :param heights: The height for shapes that represent cavity, dipole, Quad and solenoid
        :param s_start: offset of the starting point, for display only
        :param other_components: Dictionary for other components to plot.
        Format:{components:{'real_length':True, 'type':'on_axis', 'height':0.2, 'color':'b'}}
        :return:
        '''
        if self.elementPosInUseLine is None:
            self.setUseLine()
        bl_pos, bl_list = self.elementPosInUseLine, self.useLineList
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        from matplotlib.patches import Ellipse

        plt_axis.set_xlim([s_start, s_start + bl_pos[-1]])

        plt_axis.set_ylim([-1, 1])
        plt_axis.set_yticks([])
        # plt_axis.xaxis.set_ticks_position('none')
        plt_axis.axhline(0, color='black')
        other_components = {k.upper(): v for k, v in other_components.items()}
        for i in range(len(bl_list)):
            start = bl_pos[i] + s_start
            ele_name = bl_list[i]
            ele=self.elementList[self.elementNameDict[ele_name]]
            typelist = self.getParentElements(ele_name)

            l_ele = bl_pos[i + 1] - bl_pos[i]
            #lasttype = self.elementList[indlist[-1]]['TYPE']
            lasttype=typelist[-1]
            if l_ele > 0:
                if self.isCavity(ele, lasttype):
                    tempdict=self.isCavity(ele, lasttype)

                    #plt_axis.add_patch(
                    #    Rectangle((start, -heights[0] / 2.0), tempdict['L'], heights[0], ec=colors[0], fc=colors[0]))
                    plt_axis.add_patch(
                        Ellipse((start+tempdict['L']/2.0, 0), tempdict['L'], heights[0] ,linewidth=1, color=colors[0]))

                elif self.isDipole(ele, lasttype):
                    tempdict=self.isDipole(ele, lasttype)
                    shift = heights[1] * 0.2 * np.sign(tempdict['K1'])
                    plt_axis.add_patch(
                        Rectangle((start, -heights[1] / 2.0 + shift), tempdict['L'], heights[1], angle=0.0, ec=colors[1],
                                  fc='none'))

                elif self.isQuadrupole(ele, lasttype):
                    tempdict = self.isQuadrupole(ele, lasttype)
                    shift = heights[2] * 0.5 * np.sign(tempdict['K1'])

                    plt_axis.add_patch(
                        Rectangle((start, -heights[2] / 2.0 + shift), tempdict['L'], heights[2], angle=0.0, ec=colors[2],
                                  fc='none'))

                elif self.isSolenoid( ele, lasttype):
                    tempdict = self.isSolenoid(ele, lasttype)
                    plt_axis.add_patch(
                        Rectangle((start, -heights[3] / 2.0), tempdict['L'], heights[3], angle=0.0, ec=colors[3], fc='none'))

            elif not self.isDrift(ele, lasttype) and lasttype in other_components:
                p = other_components[lasttype]

                if p.get('real_length', False):
                    pass
                else:
                    if p.get('type', 'on_axis') == 'on_axis':
                        h = p.get('height', 0.4)
                        c = p.get('color', 'y')
                        plt_axis.axvline(start + l_ele / 2.0, ymin=0.5 - h / 4.0, ymax=0.5 + h / 4.0, color=c)
                    else:
                        h = p.get('height', 0.3)
                        c = p.get('color', 'k')
                        plt_axis.axvline(start + l_ele / 2.0, ymin=1.0-h/2, ymax=1.0, color=c)
                        plt_axis.axvline(start + l_ele / 2.0, ymin=0 / 2, ymax=h/2, color=c)





    def getBeamlineIndex(self, linename):
        linename = linename.upper()
        if linename[0] == '-':
            return self.getBeamlineIndex(linename[1:])
        elif linename in self.beamlineNameDict:
            return self.beamlineNameDict[linename]
        else:
            return None

    def appendToBeamline(self, linename, *elenames):
        linename = linename.upper()

        ind = self.getBeamlineIndex(linename)
        if ind is None:
            addaline = {}
            addaline['NAME'] = linename
            addaline['LINE'] = []
            addaline['__USAGE_DICT']={}
            cur_len = len(self.beamlineList)
            self.beamlineList.append(addaline)
            self.beamlineNameDict[linename] = cur_len
            ind=cur_len

        for elename in elenames:
            ind_ele = self.getElementIndex(elename)  # element is really an element
            ind_line = self.getBeamlineIndex(elename)  # element is a line?
            if ind_ele is None and ind_line is None:
                print('No element or beamline named {} are defined yet in line {} when appending to a line.'.format(elename, linename))
                raise KeyError
            else:
                self.beamlineList[ind]['LINE'].append(elename.upper())
                cur_pos=len(self.beamlineList[ind]['LINE'])-1
                if ind_line is not None:
                    # self.beamlineNameDict[linename]=checkline
                    # self.beamlineNameDict[elenames]=ind
                    # self.beamlineList[ind],self.beamlineList[checkline]=self.beamlineList[checkline],self.beamlineList[ind]
                    if ind_line > ind:
                        temp = self.beamlineList.pop(ind)
                        self.beamlineList.append(temp)
                        for k, v in self.beamlineNameDict.items():
                            if v > ind:
                                self.beamlineNameDict[k] = v - 1
                        self.beamlineNameDict[linename] = len(self.beamlineList) - 1
                        ind_line = self.getBeamlineIndex(elename)
                        ind = len(self.beamlineList) - 1

                    if linename not in self.beamlineList[ind_line]['__USAGE_DICT']:
                        self.beamlineList[ind_line]['__USAGE_DICT'][linename]=[]
                    self.beamlineList[ind_line]['__USAGE_DICT'][linename].append(cur_pos)
                if ind_ele is not None:
                    if linename not in self.elementList[ind_ele]['__USAGE_DICT']:
                        self.elementList[ind_ele]['__USAGE_DICT'][linename]=[]
                    self.elementList[ind_ele]['__USAGE_DICT'][linename].append(cur_pos)



    def elementInLine(self, elename, linename):
        '''
        check if the element with name elename is in the line of linename
        '''
        ind_ele=self.getElementIndex(elename)
        ind_line=self.getBeamlineIndex(elename)
        if ind_ele is not None:
            ind = ind_ele
            dict_temp = self.elementList[ind]['__USAGE_DICT']
        elif ind_line is not None:
            ind = ind_line
            dict_temp = self.beamlineList[ind]['__USAGE_DICT']
        else:
            return False
        if linename.upper() in dict_temp:
            return True
        else:
            for k,v in dict_temp.items():
                if self.getBeamlineIndex(k) is not None:
                    return self.elementInLine(k, linename)

    def addReverseLine(self, newlinename, linename):
        '''
        Add reversed line with newlinename
        '''
        count = self.getBeamlineIndex(linename)
        if count is None:
            print('The line with name {} can not be found, no reverse line can be added'.format(linename))
        self.beamlineNameDict[newlinename] = len(self.beamlineList)
        self.beamlineList.append(copy.deepcopy(self.beamlineList[count]))
        self.beamlineList[-1]['NAME'] = newlinename
        self.beamlineList[-1]['LINE'].reverse()
        for i in range(len(self.beamlineList[-1]['LINE'])):
            elename = self.beamlineList[-1]['LINE'][i]
            if elename[0] == '-':
                self.beamlineList[-1]['LINE'][i] = elename[1:]
            #elif self.checkAsymmetryType(self.elementList[self.elementNameDict[elename]]):
            #    self.beamlineList[-1]['LINE'][i] = '-' + elename
            #else:
            #    pass

    def _expandLine(self, line_ind):

        expandedLine=[]
        briefline = self.beamlineList[line_ind]['LINE']
        for elename in briefline:
            if elename in self.beamlineNameDict:
                expandedLine += self._expandLine(self.beamlineNameDict[elename])
            else:
                expandedLine.append(elename)
        return expandedLine

    def setUseLine(self, linename=-1):
        if linename == -1:
            linename = self.beamlineList[-1]['NAME']
            line_ind = len(self.beamlinelist)-1
        else:
            linename=linename.upper()
            line_ind = self.getBeamlineIndex(linename)
        self.useline=linename

        self.elementPosInUseLine = [0.0, ]

        if line_ind is None:
            print('The beamline {} does no exist, can not prepare the line to be used'.format(linename))
            raise KeyError

        self.useLineList = self._expandLine(line_ind)

        ele_ind=0
        for ele in self.useLineList:
            ind_ele = self.getElementIndex(ele)
            self.elementList[ind_ele]['__ID_IN_USELINE'].append(ele_ind)
            l_ele = self.getElementProperties(ele, 'L')
            if l_ele==None:
                l_ele=0.0
            self.elementPosInUseLine.append(self.elementPosInUseLine[-1] + l_ele)
            ele_ind += 1


    def loadAnElement(self, fromlattice, elename, prefix=''):
        elename = elename.upper()
        prefix = prefix.upper()
        ind = fromlattice.getElementIndex(elename)
        if ind is not None:
            ele = copy.deepcopy(fromlattice.elementList[ind])
            ind_this = self.getElementIndex(prefix + elename)
            if ind_this is not None:
                compare=self.compareElements(fromlattice, prefix + ele['NAME'], other_name=ele['NAME'])
                if compare:
                    return
                else:
                    print('Warning, the element {} has different definition'.format(prefix+elename))
                    return

            if self.checkType(ele['TYPE']):
                ele['NAME'] = prefix + ele['NAME']
                ele['__USAGE_DICT'] = {}
                ele['__ID_IN_USELINE'] = []
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)
            elif self.getElementIndex(prefix + ele['TYPE']) is None:
                self.loadAnElement(fromlattice, ele['TYPE'], prefix)
                ele['NAME'] = prefix + ele['NAME']
                ele['TYPE'] = prefix + ele['TYPE']
                ele['__USAGE_DICT'] = {}
                ele['__ID_IN_USELINE'] = []
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)
            elif self.compareElements(fromlattice, prefix + ele['TYPE'], other_name=ele['TYPE']):
                ele['NAME'] = prefix + ele['NAME']
                ele['TYPE'] = prefix + ele['TYPE']
                ele['__USAGE_DICT'] = {}
                ele['__ID_IN_USELINE'] = []
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)
                pass
            else:
                print('Warning, the root element {} has different definition'.format(ele['TYPE']))
                exit(-1)
                '''self.loadAnElement(fromlattice, ele['TYPE'], prefix)
                ele['NAME'] = prefix + ele['NAME']
                ele['TYPE'] = prefix + ele['TYPE']
                self.elementNameDict[ele['NAME']] = len(self.elementList)
                self.elementList.append(ele)'''


        else:
            print('Can not load element {} from lattice'.format(elename))
            exit(-1)

    def loadALine(self, fromlattice, linename, reverse=False, prefix='', newname=''):
        linename = linename.upper()
        if newname == '':
            newname = linename
        else:
            newname = newname.upper()
        prefix = prefix.upper()
        ind = fromlattice.getBeamlineIndex(linename)
        ind_this = self.getBeamlineIndex(prefix + newname)
        templine = []
        if ind is not None and ind_this is None:
            theline = copy.deepcopy(fromlattice.beamlineList[ind]['LINE'])

            if reverse:

                theline.reverse()
                for i in range(len(theline)):
                    if theline[i][0] == '-':
                        theline[i] = theline[i][1:]
                    else:
                        pass
                        # print(theline)
                        # for elename in reversed(fromlattice.beamlineList[pos]['LINE']):
                        #    if fromlattice.checkElementName(elename) is not None:
                        #        self.loadAnElement(fromlattice, elename, prefix)
                        #        #self.addToBeamLine(prefix+newname, prefix+elename)
                        #        templine.append(prefix+elename)
                        #    elif fromlattice.checkBeamlineName(elename)[0]>0:
                        #        self.loadALine(fromlattice, elename, reverse=True, prefix=prefix)
                        # self.addToBeamLine(prefix+newname, prefix+elename)
                        #        templine.append(prefix+elename)

            for elename in theline:
                if fromlattice.getElementIndex(elename) is not None:
                    if elename[0] == '-':
                        self.loadAnElement(fromlattice, elename[1:], prefix)
                        templine.append(prefix + elename[1:])

                    else:

                        self.loadAnElement(fromlattice, elename, prefix)
                        templine.append(prefix + elename)


                elif fromlattice.getBeamlineIndex(elename) is not None:

                    if elename[0] == '-':
                        self.loadALine(fromlattice, elename[1:], reverse=False, prefix=prefix)
                        templine.append('-' + prefix + elename[1:])
                    else:
                        self.loadALine(fromlattice, elename, reverse=False, prefix=prefix)
                        templine.append(prefix + elename)
            self.appendToBeamline(prefix + newname, *templine)
        elif ind is None:
                print("The line {} doesnot exist in the source".format(linename))
                exit()
        elif ind_this is not None:
                print("The line {} is already in the target".format(newname))
        else:
                print("Something weird happend in loadALine when loading {}".format(linename))
                exit()
        return prefix + newname

