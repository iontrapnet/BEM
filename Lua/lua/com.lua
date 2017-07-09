require 'CLRPackage'
import 'System'
import 'System.Reflection'

local comtype = luanet.ctype(__ComObject)
local get_flags = luanet.enum(BindingFlags,'GetProperty,IgnoreCase,Public')
local set_flags = luanet.enum(BindingFlags,'SetProperty,IgnoreCase,Public')
local call_flags = luanet.enum(BindingFlags,'InvokeMethod,IgnoreCase,Public')
local empty = luanet.make_array(Object,{})

local com = {}

function com.get(obj, member)
	local ok,res = pcall(comtype.InvokeMember,comtype,member,get_flags,nil,obj,empty)
	return res,ok
end

function com.set(obj, member, value)
	local ok,res = pcall(comtype.InvokeMember,comtype,member,set_flags,nil,obj,luanet.make_array(Object,{value}))
	return res,ok
end

function com.invoke(obj, member, args, mods)
	local ok,res
	if mods then
		ok,res = pcall(comtype.InvokeMember,comtype,member,call_flags,nil,obj,args,mods,nil,nil)
	else
		ok,res = pcall(comtype.InvokeMember,comtype,member,call_flags,nil,obj,luanet.make_array(Object,args))
	end
	return res,ok
end

function com.wrapper(obj)
	return setmetatable({com=obj},{
	__index = function (self,key)
		return com.get(self.com,key)
	end,
	__newindex = function (self,key,value)
		return com.set(self.com,key,value)
	end,
	__call = function (self,key,...)
		return com.invoke(self.com,key,{...})
	end
	})
end

function com.CreateObject(progid)
    local ft = Type.GetTypeFromProgID(progid)
    local f = Activator.CreateInstance(ft)
    return com.wrapper(f)
end

return com