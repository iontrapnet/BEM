from __future__ import print_function
import sys, traceback, re, tokenize, struct
from ctypes import *

X64 = struct.calcsize('P') == 8
PY3 = sys.version_info[0] == 3

if X64:
    ml = windll.ml64i4
else:
    ml = windll.ml32i4
ml.MLInitialize.restype = c_void_p
ml.MLOpenString.restype = c_void_p
ml.MLErrorMessage.restype = c_char_p

ml.ILLEGALPKT = 0
ml.CALLPKT = 7
ml.EVALUATEPKT = 13
ml.RETURNPKT = 3

ml.MLTKERR = 0
ml.MLTKSTR = 34
ml.MLTKSYM = 35
ml.MLTKREAL = 42
ml.MLTKINT = 43
ml.MLTKFUNC = 70

if PY3:
    integer = (int,)
    string = (str,)
    long = int
    import builtins as __builtin__
    #from past.builtins import execfile
    def execfile(filename, myglobals=None, mylocals=None):
        """
        Read and execute a Python script from a file in the given namespaces.
        The globals and locals are dictionaries, defaulting to the current
        globals and locals. If only globals is given, locals defaults to it.
        """
        if myglobals is None:
            # There seems to be no alternative to frame hacking here.
            caller_frame = inspect.stack()[1]
            myglobals = caller_frame[0].f_globals
            mylocals = caller_frame[0].f_locals
        elif mylocals is None:
            # Only if myglobals is given do we set mylocals to it.
            mylocals = myglobals
        if not isinstance(myglobals, Mapping):
            raise TypeError('globals must be a mapping')
        if not isinstance(mylocals, Mapping):
            raise TypeError('locals must be a mapping')
        with open(filename, "rbU") as fin:
             source = fin.read()
        code = compile(source, filename, "exec")
        exec_(code, myglobals, mylocals)
else:
    integer = (int, long)
    string = (str, unicode)
    import __builtin__

class Expr(tuple):
    def __new__ (cls, *args):
        return super(Expr, cls).__new__(cls, tuple(args))
        
    def __str__(self):
        head = self[0]
        listq = head == 'List'
        if listq:
            r = '{'
        else:
            r = str(head) + '['
        if len(self) > 1:
            r += str(self[1])
            for i in self[2:]:
                r += ', ' + str(i)
        return r + ('}' if listq else ']')

slots = {}
def slot(key, *args):
    if type(key) not in integer:
        return key
    if len(args) == 0:
        return slots.get(key, None)
    else:
        slots[key] = args[0]
        return args[0]

PyObject = slot
PyFunction = slot

def pyeval(expr):
    t = type(expr)
    if t == Expr:
        size = len(expr)
        t = [pyeval(i) for i in expr]
        if callable(t[0]):
            try:
                r = t[0](*t[1:])
            except Exception as e:
                r = Expr(u'Error', *traceback.format_exc().splitlines()[3:])
            finally:
                return r
        else:
            return Expr(*t)
    elif t in string:
        if expr[0] == '"':
            return expr[1:-1]
        else:
            return globals().get(expr, getattr(__builtin__, expr, expr))
    else:
        return expr

def pyget(*args):
    r = args[0]
    for v in args[1:]:
        if type(v) in string:
            if v[0] == '"':
                v = v[1:-1]
            r = getattr(r, v, None)
        else:
            r = None
    return r
Part = pyget
MessageName = pyget

def pyset(*args):
    r = args[0]
    size = len(args)
    for i in range(1,len-1):
        v = args[i]
        if type(v) in string:
            if v[0] == '"':
                v = v[1:-1]
            if i == len - 2:
                try:
                    setattr(r, v, args[len-1])
                except Exception as e:
                    traceback.print_exc()
            r = getattr(r, v, None)
    return r
    
class PythonLink(object):
    @classmethod
    def init(cls):
        cls.env = ml.MLInitialize(c_void_p())
    
    @classmethod
    def deinit(cls):
        ml.MLDeinitialize(c_void_p(cls.env))
            
    def __init__(self):
        pass
        
    def install(self):
        err = (c_int*1)(0)
        mldir = br'C:\Program Files\Wolfram Research\Mathematica\11.0'
        #mlarg = b"-linklaunch -linkname '" + mldir + br"\math.exe'"
        mlarg = (' '.join(sys.argv[1:])).encode('utf-8')
        self.mlp = ml.MLOpenString(c_void_p(self.env), mlarg, byref(err))
        ml.MLConnect(c_void_p(self.mlp))
        for i in range(len(self.defs)):
            v = self.defs[i]
            self.put_expr(Expr(u'DefineExternal',u'"' + v[1] + u'"', u'"' + v[2] + u'"', i))
        self.put_expr(Expr(u'ToExpression', u'''"
        ClearAll[PyObject, PyFunction];
        Format[PyObject[id_Integer]] := Py[repr[slot[id]]];
        Format[PyFunction[id_Integer]] := Py[repr[slot[id]]];
        PyObject[id_Integer][(f:_[___])[args___]] := PyObject[id][f][args];
        PyObject[id_Integer][(f:Except[List])[args___]] := PyFunction[id, f][args];
        PyObject[id_Integer][k:Except[_List]] := Py[pyget[slot[id], k]];
        PyObject[id_Integer][{k__}] := Py[pyget[slot[id], k]];
        PyObject[id_Integer][k:Except[_List], v_] := Py[pyset[slot[id], k, v]];
        PyObject[id_Integer][{k__}, v_] := Py[pyset[slot[id], k, v]];
        PyFunction[id_Integer, k___][(f:Except[PyObject|PyFunction|List|Association])[args___]] := PyFunction[id, k, f][args];
        PyFunction[id_Integer][args___] := Py[slot[id][args]];
        PyFunction[id_Integer, k___, $][i_] := Py[pyget[slot[id], k, "__getitem__"][i]];
        PyFunction[id_Integer, k___, $][i_, v_] := Py[pyget[slot[id], k, "__setitem__"][i,v]];
        PyFunction[id_Integer, k__][args___] := Py[pyget[slot[id], k][args]];
        PyObject[PyFunction[id_Integer]] := PyObject[id];
        PyFunction[PyObject[id_Integer]] := PyFunction[id];
        PyObject[obj:Except[_Integer|_Blank|_Pattern]] := PyObject@@Py[obj];
        PyFunction[obj:Except[_Integer|_Blank|_Pattern]] := PyFunction@@Py[obj];
        $[PyObject[id_Integer]] := Py@mexpr@slot@id;
        "'''))
        ml.MLPutSymbol(c_void_p(self.mlp), b'End')
        ml.MLFlush(c_void_p(self.mlp))
    
    def uninstall(self):
        ml.MLClose(c_void_p(self.mlp))
    
    def get_error(self):
        r = ml.MLErrorMessage(c_void_p(self.mlp))
        r = b'"' + r + b'"'
        r = r.decode('utf-8')
        ml.MLClearError(c_void_p(self.mlp))
        ml.MLNewPacket(c_void_p(self.mlp))
        return Expr('Error', r)
    
    def get_packet(self):
        waitResult = ml.MLWaitForLinkActivity(c_void_p(self.mlp))
        print('waitResult = ',waitResult)
        if not waitResult:
            return 0
        pkt = ml.MLNextPacket(c_void_p(self.mlp))
        print('pkt = ',pkt)
        return pkt
    
    def get_expr(self):
        t = ml.MLGetNext(c_void_p(self.mlp))
        if t == ml.MLTKINT:
            r = (c_longlong*1)(0)
            return r[0] if ml.MLGetInteger64(c_void_p(self.mlp), byref(r)) else self.get_error()
        elif t == ml.MLTKREAL:
            r = (c_double*1)(0)
            ml.MLGetReal64(c_void_p(self.mlp), byref(r))
            return r[0]
        elif t == ml.MLTKSYM:
            bufp = (POINTER(c_char)*1)()
            ml.MLGetSymbol(c_void_p(self.mlp), byref(bufp))
            r = cast(bufp[0], c_char_p).value
            ml.MLReleaseString(c_void_p(self.mlp), bufp[0])
            return r.decode('utf-8')
        elif t == ml.MLTKSTR:
            bufp = (POINTER(c_char)*1)()
            size = (c_int*1)(0)
            ml.MLGetByteString(c_void_p(self.mlp), byref(bufp), byref(size), 0)
            r = cast(bufp[0], POINTER(c_char * size[0])).contents[:]
            r = b'"' + r + b'"'
            ml.MLReleaseByteString(c_void_p(self.mlp), bufp[0], size[0])
            return r.decode('utf-8')
        elif t == ml.MLTKFUNC:
            n = (c_int*1)(0)
            ml.MLGetArgCount(c_void_p(self.mlp), byref(n))
            n = n[0]
            return Expr(*(self.get_expr() for i in range(n+1)))
    
    def put_expr(self, expr):
        t = type(expr)
        if t == Expr:
            head = expr[0]
            if type(head) in string:
                ml.MLPutFunction(c_void_p(self.mlp), head.encode('utf-8'), len(expr) - 1)
            else:
                ml.MLPutNext(c_void_p(self.mlp), ml.MLTKFUNC)
                ml.MLPutArgCount(c_void_p(self.mlp), len(expr) - 1)
                self.put_expr(head)
            for i in expr[1:]:
                self.put_expr(i)
        elif t in string:
            if expr[0] == '"':
                expr = expr[1:-1].encode('utf-8')
                ml.MLPutByteString(c_void_p(self.mlp), expr, len(expr))
            #elif expr.isidentifier():
            elif re.match(tokenize.Name + '$', expr):
                ml.MLPutSymbol(c_void_p(self.mlp), expr.encode('utf-8'))
            else:
                expr = expr.encode('utf-8')
                ml.MLPutByteString(c_void_p(self.mlp), expr, len(expr))
        elif t in integer:
            ml.MLPutInteger64(c_void_p(self.mlp), c_longlong(expr))
        elif t == float:
            ml.MLPutReal64(c_void_p(self.mlp), c_double(expr))
        elif t == type(None):
            ml.MLPutSymbol(c_void_p(self.mlp), b'Null')
        else:
            key = long(id(expr))
            slot(key, expr)
            if callable(expr):
                ml.MLPutFunction(c_void_p(self.mlp), b'PyFunction', 1)
            else:
                ml.MLPutFunction(c_void_p(self.mlp), b'PyObject', 1)
            ml.MLPutInteger64(c_void_p(self.mlp), c_longlong(key))
    
    def answer(self):
        pkt = self.get_packet()
        if pkt == ml.CALLPKT:
            t = self.get_expr()
            expr = self.get_expr()
            f = self.defs[t]
            if type(f) == tuple:
                f = f[0]
            r = f(self,expr)
            self.put_expr(r)
            ml.MLEndPacket(c_void_p(self.mlp))
        elif pkt == ml.EVALUATEPKT:
            expr = (self.get_expr(),)
            expr = self.py(expr)
            self.put_expr(expr)
            ml.MLEndPacket(c_void_p(self.mlp))
        elif pkt == ml.RETURNPKT:
            expr = self.get_expr()
            return expr
        else:
            print(self.get_expr())
        return pkt != 0
    
    def py(self,expr):
        if len(expr) == 2:
            r = expr[1]
            if type(r) in string and r[0] == '"':
                r = r[1:-1]
                if r[:4] == 'def ' or r[:5] == 'from ' or r[:7] == 'import ' or r.find('=') > 0:
                    # exec
                    try:
                        r = compile(r, '<string>', 'exec')
                    except Exception as e:
                        return Expr(u'Error', *traceback.format_exc().splitlines()[3:])
                try:
                    r = eval(r, globals())
                except Exception as e:
                    return Expr(u'Error', *traceback.format_exc().splitlines()[3:])
            return pyeval(r)
        elif len(expr) == 3:
            r = expr[1]
            if type(r) in string:
                if r[0] == '"':
                    r = r[1:-1]
                t = self.py(Expr(None,expr[2])) 
                g = globals()
                g[r] = t
                return t
            elif type(r) == Expr:
                return pyset(*(r[1:] + (self.py(Expr(None,expr[2])),)))
        
    defs = [
        (py, u'Py[expr__]', u'{expr}')
    ]
    
def meval(*args):
    expr = args[-1]
    if type(expr) in string:
        expr = Expr(u'ToExpression', '"' + expr + '"')
    if len(args) == 2:
        expr = Expr(u'Set', args[0], expr)
    pl.put_expr(Expr(u'EvaluatePacket', expr))
    ml.MLEndPacket(c_void_p(pl.mlp))
    return pl.answer()

def mfunc(expr):
    def f(*args):
        return meval(pl,Expr(expr, *args))
    return f

def mexpr(expr):
    t = type(expr)
    if t in (tuple, list):
        return Expr('List', *expr)
    elif t == dict:
        r = []
        for k,v in expr.items():
            r.append(Expr('Rule', k, mexpr(v)))
        return Expr('Association', *r)
    else:
        return expr
        
#globals()['$'] = mexpr

def List(*args):
    return list(args)

def Association(*args):
    r = {}
    for v in args:
        if type(v) == Expr:
            r[v[1]] = v[2]
    return r

if __name__ == '__main__':    
    PythonLink.init()
            
    pl = PythonLink()
    pl.install()

    r = pl.answer()
    while r:
        r = pl.answer()
        
    #pl.answer()
    #pl.put_expr(Expr('EvaluatePacket', Expr('ToExpression', '"100!"')))
    #print(pl.answer())
    #pl.put_expr(Expr('EvaluatePacket', Expr('ToExpression', '"10!"')))
    #print(pl.answer())

    pl.uninstall()
    
    PythonLink.deinit()