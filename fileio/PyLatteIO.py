from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ply.lex as lex
import ply.yacc as yacc

import numpy as np
import copy
from .fileIO import LatticeFile
from .. import element
import re

class latteLexer:
    literals = "=,-:&()*/"
    # reserved={r'[Ll][Ii][Nn][Ee]':'KWLINE'}
    reserved = {'LINE': 'KWLINE', 'USE': 'KWUSE'}
    tokens = ['NUMBER', 'STRING', 'RET'] + list(reserved.values())

    def t_NUMBER(self, t):
        r'([0-9]*\.?[0-9]+|[0-9]+\.[0-9]*)(([eE][-+]?[0-9]+)?)'
        t.value = float(t.value)
        return t

    t_ignore = ' \t'

    def t_STRING(self, t):
        r'-?[A-Za-z][A-Za-z0-9_\.]*'
        t.type = latteLexer.reserved.get(t.value, 'STRING')
        return t

    def t_QUOTE(self, t):
        r'\"'
        pass

    def t_COMMENT(self, t):
        r'\#.*'
        pass

    def t_contline(self, t):
        r'\&[ \t]*\n'
        # print('foundit')
        t.lexer.lineno += 1
        self.continuetor = 1
        pass

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


class latteParser:
    # starting = "statement"
    tokens = latteLexer.tokens

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
                | bline ',' STRING '''
        if len(p) == 2:
            self.listtemp.append(p[1])
        elif len(p) == 4:
            self.listtemp.append(p[3])

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


class PyLatteLatticeFile(LatticeFile):
    def __init__(self):
        LatticeFile.__init__(self)

    def checkType(self, typename, parameterName=None):
        if typename.upper() in element.Element.elementTypes:
            if parameterName is None:
                return True
            else:
                parameterList=element.__dict__[typename.title()].propertyNames
                return parameterName.upper() in parameterList
        else:
            if typename.upper() in self.elementNameDict:
                self.checkType(self.getElementRootType(typename.upper()), parameterName)
        return False

    def isDrift(self, ele, parent_type=None):
        if ele['TYPE']=='DRIFT' or parent_type=='DRIFT':
            return {'L':ele['L']}
        else:
            return False

    def isDipole(self, ele, parent_type=None):
        if ele['TYPE'] == 'DIPOLE' or parent_type=='DIPOLE':
            return {'L': ele['L'], 'ANGLE':ele['ANGLE'], 'K1':ele['K1']}
        else:
            return False

    def isQuadrupole(self, ele, parent_type=None):
        if ele['TYPE'] == 'QUADRUPOLE' or parent_type=='QUADRUPOLE':
            return {'L': ele['L'], 'K1':ele['K1']}
        else:
            return False
    def isSolenoid(self,ele, parent_type=None):
        if ele['TYPE'] == 'SOLENOID' or parent_type=='SOLENOID':
            return {'L': ele['L'], 'B':ele['B']}
        else:
            return False
    def isCavity(self,ele, parent_type=None):
        if ele['TYPE'] == 'CAVITY' or parent_type=='CAVITY':
            return {'L': ele['L']}
        else:
            return False
