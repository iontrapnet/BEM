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
local wr = bem.new_world(base, 0.001, 64, 6, 6, 1, 1000)
if ext == 'csv' then
    bem.load_csv(wr, path, 44)
elseif ext == 'dxf' then
    if argc < 2 then return -1 end
    bem.load_dxf(wr, path, tonumber(argv[2]))
end
bem.solve(wr)
local X0 = tonumber(argv[2])
local X1 = tonumber(argv[3])
local nX = tonumber(argv[4])
local Y0 = tonumber(argv[5])
local Y1 = tonumber(argv[6])
local nY = tonumber(argv[7])
local Z0 = tonumber(argv[8])
local Z1 = tonumber(argv[9])
local nZ = tonumber(argv[10])
bem.region(wr, X0, X1, nX, Y0, Y1, nY, Z0, Z1, nZ)
bem.del_world(wr)