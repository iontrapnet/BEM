ffi = require('ffi')
local file = io.open('BEM.lua.h')
bem = ffi.load('BEM')
ffi.cdef(file:read('*all'))
file:close()

local argv = {...}
local argc = #argv
if argc < 1 then return -1 end

local path = argv[1]
--local base, ext = string.match(path, ".-([^\\]-).([^\\%.]+)$")
local base, ext = string.match(path, "(.-).([^\\%.]+)$")
local wr = bem.new_world(base, 0.000001, 64, 6, 6, 0.001, 1000)
if ext == 'csv' then
    bem.load_csv(wr, path, 44)
elseif ext == 'dxf' then
    if argc < 2 then return -1 end
    bem.load_dxf(wr, path, tonumber(argv[2]))
end
bem.solve(wr)
bem.del_world(wr)