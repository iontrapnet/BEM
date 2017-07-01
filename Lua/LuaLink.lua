require 'mathlink_h'
ffi = require('ffi')
--ml = ffi.load('ml64i4.dll')
ml = ffi.load('ml32i4.dll')

int64 = ffi.typeof(0LL)
function is_integer(x)
    return type(x) == 'cdata' and ffi.istype(x, int64)
end

local function to_string(self)
    if is_expr(self) then
        local len = #self
        local r = to_string(self[0])..'['
        for i=1,len-1 do
            r = r..to_string(self[i])..','
        end
        if len > 0 then
            r = r..to_string(self[len])
        end
        r = r..']'
        return r
    elseif is_integer(self) then
        return tostring(self):sub(1,-3)
    else
        return tostring(self)
    end
end
Expr = {
    __index = Expr,
    __tostring = to_string,
    __call = function (self, ...)
        if #self > 0 then
            return Expr{[0] = self, ...}
        else
            return Expr{[0] = self[0], ...}
        end
    end
}
setmetatable(Expr, {__call = function (self, obj)
    if type(obj) ~= 'table' or is_expr(obj) then
        obj = {[0] = obj}
    end
    setmetatable(obj, Expr)
    return obj
end})
function is_expr(self)
    return getmetatable(self) == Expr
end

Error = function (str)
    return Expr('Error')('"'..tostring(str)..'"')
end
local function get_error()
    local r = ml.MLErrorMessage(mlp)
    r = '"'..ffi.string(r)..'"'
    ml.MLClearError(mlp)
    ml.MLNewPacket(mlp)
    return Error(r)
end

local function get_packet()
    local waitResult = ml.MLWaitForLinkActivity(mlp)
    print('waitResult = ',waitResult)
    if not waitResult then
        return 0
    end
    local pkt = ml.MLNextPacket(mlp)
    print('pkt = ',pkt)
    return pkt
end

local function get_expr()
    local t = ml.MLGetNext(mlp)
    if t == ml.MLTKSYM then
        local bufp = ffi.new('const char*[1]')
        ml.MLGetSymbol(mlp, bufp)
        local r = ffi.string(bufp[0])
        ml.MLReleaseString(mlp, bufp[0])
        return r
    elseif t == ml.MLTKSTR then
        local bufp = ffi.new('const char*[1]')
        local len = ffi.new('int[1]')
        ml.MLGetByteString(mlp, bufp, len, 0)
        r = '"'..ffi.string(bufp[0],len[0])..'"'
        ml.MLReleaseByteString(mlp, bufp[0], len[0])
        return r
    elseif t == ml.MLTKINT then
        local r = ffi.new('mlint64[1]')
        return (ml.MLGetInteger64(mlp, r) ~= 0 and r[0]) or get_error()
    elseif t == ml.MLTKREAL then
        local r = ffi.new('double[1]')
        ml.MLGetReal64(mlp, r)
        return r[0]
    elseif t == ml.MLTKFUNC then
        local n = ffi.new('int[1]')
        ml.MLGetArgCount(mlp, n)
        local r = {}
        for i=0,n[0] do
            r[i] = get_expr()
        end
        return Expr(r)
    end
    return t
end

local slots = {}
function slot(key, ...)
    if type(key) ~= 'string' then
        return key
    end
    local args = {...}
    if #args == 0 then
        return slots[key]
    else
        slots[key] = args[1]
        return args[1]
    end
end
LuaFunction = slot
LuaObject = slot

local function put_expr(expr)
    if is_expr(expr) then
        local head = expr[0]
        local len = #expr
        if type(head) == 'string' then
            ml.MLPutFunction(mlp, head, len)
        else
            ml.MLPutNext(mlp, ml.MLTKFUNC)
            ml.MLPutArgCount(mlp, len)
            put_expr(head)
        end
        for i=1,len do
            put_expr(expr[i])
        end
    else
        local t = type(expr)
        if t == 'string' then
            if expr:sub(1,1) == '"' then
                expr = expr:sub(2,-2)
                ml.MLPutByteString(mlp, expr, #expr)
            elseif expr:match('^[a-zA-Z$][a-zA-Z$0-9]*$') then
                ml.MLPutSymbol(mlp, expr)
            else
                ml.MLPutByteString(mlp, expr, #expr)
            end
        elseif t == 'number' then
            ml.MLPutReal64(mlp, expr)
        elseif t == 'nil' then
            ml.MLPutSymbol(mlp, 'Null')
        elseif t == 'cdata' and ffi.istype(expr, int64) then
            ml.MLPutInteger64(mlp, expr)
        else
            local key = tostring(expr)
            slot(key, expr)
            if t == 'function' then
                ml.MLPutFunction(mlp, 'LuaFunction', 1)
            else
                ml.MLPutFunction(mlp, 'LuaObject', 1)
            end
            ml.MLPutString(mlp, key)
        end
    end
end

List = Expr('List')
function list(expr)
    if type(expr) == 'table' then
        if is_expr(expr) then
            return expr
        else
            --return List(unpack(expr))
            local r = Expr('List')
            for k,v in ipairs(expr) do
                r[#r+1] = list(v)
            end
            return r
        end
    else
        return expr
    end
end

Rule = Expr('Rule')
function dict(expr)
    if type(expr) == 'table' then
        if is_expr(expr) then
            return expr
        else
            local r = Expr('Association')
            for k,v in pairs(expr) do
                r[#r+1] = Rule(k, dict(v))
            end
            return r
        end
    else
        return expr
    end
end

dump = string.dump

local function leval(expr)
    local t = type(expr)
    if t == 'table' then
        if is_expr(expr) then
            local len = #expr
            local t = Expr{}
            for i=0,len do
                t[i] = leval(expr[i])
            end
            if t[0] == '$' then
                if type(t[1]) == 'table' then
                    t[0] = list
                elseif type(t[1]) == 'function' then
                    t[0] = dump
                end
            elseif t[0] == '$$' then
                if type(t[1]) == 'table' then
                    t[0] = dict
                end
            end
            local r = t
            if not (type(t[0]) == 'string' or type(t[0]) == 'number') then
                r = {pcall(t[0], unpack(t))}
                print(expr[0],t[0],r[1],r[2])
                if r[1] then
                    r = r[2]
                else
                    r = Error(r[2])
                end
            end
            return r
        else
            return expr
        end
    elseif t == 'string' then
        if expr:sub(1,1) == '"' then
            return expr:sub(2,-2)
        else
            return _G[expr] or expr
        end
    else
        return expr
    end
end

function lget(...)
    local args = {...}
    local r = args[1]
    for i=2,#args do
        local t = type(r)
        if t == 'nil' or t == 'number' or t == 'string' then
            break
        else
            local v = args[i]
            if type(v) == 'string' then
                if v:sub(1,1) == '"' then
                    v = v:sub(2,-2)
                end
            else
                v = tonumber(v)
            end
            r = r[v]
        end
    end
    return r
end
Part = lget
MessageName = lget

function lset(...)
    local args = {...}
    local r = args[1]
    local len = #args
    for i=2,len-1 do
        local t = type(r)
        if t == 'nil' or t == 'number' or t == 'string' then
            break
        else
            local v = args[i]
            if type(v) == 'string' then
                if v:sub(1,1) == '"' then
                    v = v:sub(2,-2)
                end
            else
                v = tonumber(v)
            end
            if i == len-1 then
                r[v] = args[len]
            end
            r = r[v]
        end
    end
    return r
end

function lua(expr)
    if #expr == 1 then
        local r = expr[1]
        if type(r) == 'string' and r:sub(1,1) == '"' then
            if not (r:match('^"[^{]+=') or r:match('^"\x1BLJ')) then
                r = '"return '..r:sub(2,-1)
            end
            --print(r)
            if r:match('^"\x1BLJ') then
                r = Expr('loadstring')(r)
            else
                r = Expr('loadstring')(r)()
            end
        end
        return leval(r)
    elseif #expr == 2 then
        local r = expr[1]
        if type(r) == 'string' then
            if r:sub(1,1) == '"' then
               r = r:sub(2,-2)
            end
            local t = _G
            --if r:find(':') then t = slots end
            return lset(t, r, lua{expr[2]})
        elseif is_expr(r) then
            r = Expr(lset)(unpack(r))
            r[#r+1] = lua{expr[2]}
            return leval(r)
        end
    end
end

local defs = {
    {lua, 'Lua[expr__]', '{expr}'},
}
local function answer()
    local pkt = get_packet()
    if pkt == ml.CALLPKT then
        local t = tonumber(get_expr())
        local expr = get_expr()
        local f = defs[t+1]
        if type(f) == 'table' then
            f = f[1]
        end
        local r = {pcall(f, expr)}
        if r then
            r = r[2]
        else
            r = Error(r[2])
        end
        put_expr(r)
        ml.MLEndPacket(mlp)
    elseif pkt == ml.EVALUATEPKT then
        local expr = {get_expr()}
        expr = lua(expr)
        put_expr(expr)
        ml.MLEndPacket(mlp)
    elseif pkt == ml.RETURNPKT then
        local expr = get_expr()
        return expr
    else
        print(get_expr())
    end
    return (pkt ~= 0)
end

function meval(...)
    local args = {...}
    local len = #args
    local expr = args[len]
    if type(expr) == 'string' then
        expr = Expr('ToExpression')('"'..expr..'"')
    end
    if len == 2 then
        expr = Expr('Set')(args[1], expr)
    end
    put_expr(Expr('EvaluatePacket')(expr))
    ml.MLEndPacket(mlp)
    return answer()
end

function mfunc(expr)
    if not is_expr(expr) then
        expr = Expr(expr)
    end
    return function(...)
        return meval(expr(...))
    end
end

function Association(...)
    local args = {...}
    local r = {}
    for k,v in ipairs(args) do
        if is_expr(v) then
            lset(r, v[1], v[2])
        end
    end
    return r
end

local function init()
    mlenv = ml.MLInitialize(ffi.cast('void*',0))
    local err = ffi.new('int[1]',0)
    --mldir = [[C:\Program Files\Wolfram Research\Mathematica\11.0\]]
    --mlarg = "-linklaunch -linkname '"..mldir.."math.exe'"
    mlarg = table.concat(arg, ' ')
    mlp = ml.MLOpenString(mlenv, mlarg, err)
    
    ml.MLConnect(mlp)
    for i=1,#defs do
        local v = defs[i]
        put_expr(Expr('DefineExternal')('"'..v[2]..'"', '"'..v[3]..'"', i - 1))
    end
    put_expr(Expr('ToExpression')([=["
    ClearAll[LuaObject, LuaFunction];
    Format[LuaObject[id_String]] := id;
    Format[LuaFunction[id_String]] := id;
    LuaObject[id_String][(f:_[___])[args___]] := LuaObject[id][f][args];
    LuaObject[id_String][(f:Except[List])[args___]] := LuaFunction[id, f][args];
    LuaObject[id_String][k:Except[_List]] := Lua[lget[slot[id], k]];
    LuaObject[id_String][{k__}] := Lua[lget[slot[id], k]];
    LuaObject[id_String][k:Except[_List], v_] := Lua[lset[slot[id], k, v]];
    LuaObject[id_String][{k__}, v_] := Lua[lset[slot[id], k, v]];
    LuaFunction[id_String, k___][(f:Except[LuaObject|LuaFunction|List|Association])[args___]] := LuaFunction[id, k, f][args];
    LuaFunction[id_String][args___] := Lua[slot[id][args]];
    LuaFunction[id_String, k__][args___] := Lua[lget[slot[id], k][args]];
    LuaObject[LuaFunction[id_String]] := LuaObject[id];
    LuaFunction[LuaObject[id_String]] := LuaFunction[id];
    $[LuaObject[id_String]] := Lua@list@slot@id;
    $$[LuaObject[id_String]] := Lua@dict@slot@id;
    $[LuaFunction[id_String]] := Lua@dump@slot@id;
    "]=]))
    ml.MLPutSymbol(mlp, 'End')
    ml.MLFlush(mlp)
end

init()

--[[answer()
print(meval(meval(Expr('Factorial')(100))))
print(meval(Expr('foo')(Expr('Integrate')('x', 'x'), 3, '"y"')))
print(meval('foo[Integrate[x,x],3,"y"]'))]]

local r = answer()
while r do
    r = answer()
end
ml.MLClose(mlp)
ml.MLDeinitialize(mlenv)