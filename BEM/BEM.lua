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
local X0 = bem.X0
local X1 = bem.X1
local Y0 = bem.Y0
local Y1 = bem.Y1
local Z0 = bem.Z0
local Z1 = bem.Z1
local Xn = bem.Xn
local Yn = bem.Yn
local Zn = bem.Zn
bem.region(wr, X0, X1, Xn, Y0, Y1, Yn, Z0, Z1, Zn)
bem.del_world(wr)