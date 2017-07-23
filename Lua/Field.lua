local ffi = require 'ffi'
local mat = require 'matio.ffi'
local md5 = require 'md5'

local function loadmat(path, ...)
    local f = mat.open(path, mat.ACC_RDONLY)
    local args = {...}
    local r = {}
    for _,name in ipairs(args) do
        r[#r+1] = mat.varRead(f, name)
    end
    mat.close(f)
    return unpack(r)
end

local function savemat(path, ...)
    local f = mat.createVer(path, nil, mat.FT_MAT5)
    local args = {...}
    for _,arg in ipairs(args) do
        local var = mat.varCreate(arg.name, mat.C_DOUBLE, mat.T_DOUBLE, arg.rank, arg.dims, arg.data, 0)
        mat.varWrite(f, var, mat.COMPRESSION_ZLIB)
        mat.varFree(var)
    end
    mat.close(f)
end

local function strcast(x)
    if type(x) == 'number' then
        return ffi.string(ffi.new('double[1]',x),8)
    elseif type(x) == 'table' then
        local r = ''
        for _,v in ipairs(x) do
            r = r..strcast(v)
        end
        return r
    else
        return tostring(x)
    end
end

local function DataHash(xr,yr,zr)
    local m = md5.new()
    m:update(strcast({'double',2,1,9}))
    m:update(strcast{xr[1],xr[2],xr[3],yr[1],yr[2],yr[3],zr[1],zr[2],zr[3]})
    return md5.tohex(m:finish())
end

local function FileExists(path)
  local file = io.open(path)
  if file then
    return file:close() and true
  end
end

local function curdir()
    return io.popen'cd':read'*l'
end

PATH_SPHERE = '..\\Model\\sphere\\922019645250603132'
PATH_4ROD = '..\\Model\\4rod\\167634622912717531'
PATH_SURFACE = '..\\Model\\surface\\6157822360140778246'

local FieldSession = {}

local function MEngine()
    local m = FieldSession['eng']
    if not m then
        local com = require 'com'
        --m = com.CreateObject('matlab.application')
        local MLApp = CLRPackage('MLApp', 'MLApp')
        m = MLApp.MLAppClass()
        m:Execute('cd '..curdir()..'\\..\\MATLAB')
        FieldSession['eng'] = m
    end
    return m
end

function FieldInit(path, xr, yr, zr, data)
    FieldSession['xr'] = xr
    FieldSession['yr'] = yr
    FieldSession['zr'] = zr
    local pb
    if not data then
        data = DataHash(xr,yr,zr)
    end
    local matfile = path..'-'..data..'.mat'
    if FileExists(matfile) then
        pb = loadmat(matfile, 'pb')
    else
        if not FileExists(path..'.mat') then
            local err, ret = MEngine():Feval('ChargeBasis',0,path)
        end
        local cb, triangles = loadmat(path..'.mat','cb','triangles')
        local dx = (xr[2]-xr[1])/xr[3]
        local dy = (yr[2]-yr[1])/yr[3]
        local dz = (zr[2]-zr[1])/zr[3]
        local nx = xr[3]+1
        local ny = yr[3]+1
        local nz = zr[3]+1
        local n = nx*ny*nz
        local ct = 'double['..n..']'
        local x = ffi.new(ct)
        local y = ffi.new(ct)
        local z = ffi.new(ct)
        for i=0,nx-1 do
            for j=0,ny-1 do
                for k=0,nz-1 do
                    local idx = i*ny*nz+j*nz+k
                    x[idx] = xr[1]+i*dx
                    y[idx] = yr[1]+j*dy
                    z[idx] = zr[1]+k*dz
                end
            end 
        end
        local cbdata = ffi.cast('double*', cb.data)
        local pbdata = ffi.new('double['..4*n*cb.dims[1]..']')
        local file = io.open('..\\BEM\\BEM.lua.h')
        local bem = ffi.load('BEM')
        ffi.cdef(file:read('*all'))
        file:close()
        for i=0,cb.dims[1]-1 do
            bem.triangles_potential_field(cb.dims[0],triangles.data,cbdata+cb.dims[0]*i,n,x,y,z,pbdata+4*n*i)
        end
        pb = {name = 'pb', rank = 3, dims = ffi.new('int[3]',{4,n,cb.dims[1]}), data = pbdata}
        savemat(matfile, pb)
    end
    FieldSession['pb'] = pb
end

function Field(voltages, points, ratio)
    if ratio == None then ratio = 1 end
    local xr = FieldSession['xr']
    local yr = FieldSession['yr']
    local zr = FieldSession['zr']
    local pb = FieldSession['pb']
    local data = ffi.cast('double*',pb.data)
    local b0 = pb.dims[0]
    local b1 = b0 * pb.dims[1]
    local s2 = 1
    if pb.rank > 2 then s2 = tonumber(pb.dims[2]) end
    local nx = xr[3] + 1
    local ny = yr[3] + 1
    local nz = zr[3] + 1
    local r = {}
    for _, point in ipairs(points) do
        local dx = (nx-1)*(point[1]-xr[1])/(xr[2]-xr[1])
        local dy = (ny-1)*(point[2]-yr[1])/(yr[2]-yr[1])
        local dz = (nz-1)*(point[3]-zr[1])/(zr[2]-zr[1])
        local ix = math.floor(dx)
        local iy = math.floor(dy)
        local iz = math.floor(dz)
        dx = dx - ix
        dy = dy - iy
        dz = dz - iz
        local coe = {(1-dx)*(1-dy)*(1-dz),(1-dx)*(1-dy)*dz,(1-dx)*dy*(1-dz),(1-dx)*dy*dz,dx*(1-dy)*(1-dz),dx*(1-dy)*dz,dx*dy*(1-dz),dx*dy*dz}
        local idx = {iz+iy*nz+ix*ny*nz,1+iz+iy*nz+ix*ny*nz,iz+(1+iy)*nz+ix*ny*nz,1+iz+(1+iy)*nz+ix*ny*nz,iz+iy*nz+(1+ix)*ny*nz,1+iz+iy*nz+(1+ix)*ny*nz,iz+(1+iy)*nz+(1+ix)*ny*nz,1+iz+(1+iy)*nz+(1+ix)*ny*nz}
        local fx = 0
        local fy = 0
        local fz = 0
        for i=1,s2 do
            for j=1,8 do
                local b = (i-1)*b1+idx[j]*b0+1
                local c = voltages[i]*coe[j]
                fx = fx + data[b]*c
                fy = fy + data[b+1]*c
                fz = fz + data[b+2]*c
            end
        end
        r[#r+1] = {ratio*fx,ratio*fy,ratio*fz}
    end
    return r
end

--[[local tic = os.clock()
FieldInit(PATH_SPHERE,{-2,2,10},{-2,2,10},{-2,2,10})
print(os.clock() - tic)
os.exit()]]

local tic = os.clock()
local xr = {-5e-5,5e-5,100}
local yr = {-5e-5,5e-5,100}
local zr = {-5e-5,5e-5,100}
FieldInit(PATH_4ROD,xr,yr,zr,'5d49bb9d6704e98ea598bb82a952b646')
print(os.clock() - tic)

local e = 1.602176565e-19
local m = 171 * 1.6605402e-27                       
local field_coeff = 100 * e / m
local Ke = 8.9875517873681e9                  --Coulomb constant
local coulomb_coeff = e^2 * Ke / m
local lam = 369.5e-9                          --cooling beam
local kB = 1.38064852e-23                     --Boltzmann constant

local u_dc = 100
local u_ac = 1000
local f_rf = 12e6
local w_rf = 2*math.pi*f_rf
local dc_z = 0.122321e4
local w_z = (2 * dc_z * e * u_dc / m)^.5
local gap_z = (e^2 * Ke / (w_z^2 * m))^(1/3)
local v0 = math.sqrt(2 * 500 * kB / 3 / m) * 1e-1

local T_total = 5e-4
local dt = (1 / f_rf) / 1000
local n_T = math.floor(T_total / dt)
local step = math.floor(n_T / 1e4)
if step < 1 then step = 1 end

local n_ions = 5

local disp = require ('display')
local win = {}
local plt = {}
local r = {{},{}}
local v = {{},{}}
for i=1,n_ions do
    local config = {title = 'Ion '..tostring(i), ylabel = 'Position (um)', labels = {"t (us)",'X','Y','Z'}}
    win[i] = disp.plot({}, config)
    plt[i] = {}
    r[1][i] = {}
    r[1][i][1] = (2 * math.random() - 1) * 2e-6
    r[1][i][2] = (2 * math.random() - 1) * 2e-6
    r[1][i][3] = (2 * math.random() - 1) * 2e-6 + (2*i - n_ions - 1) * gap_z / 2
    v[1][i] = {}
    v[1][i][1] = (2 * math.random() - 1) * v0
    v[1][i][2] = (2 * math.random() - 1) * v0
    v[1][i][3] = (2 * math.random() - 1) * v0
end

tic = os.clock()
for k = 1,n_T do
    local t = dt * (k - 1)
    local k1 = (k - 1) % 2 + 1
    if (k - 1) % step == 0 then
        local j = (k - 1)/step + 1
        for i=1,n_ions do
            plt[i][j] = {1e6*t,1e6*r[k1][i][1],1e6*r[k1][i][2],1e6*r[k1][i][3]}          
        end
    end
    local u1 = u_dc + 0*math.cos(w_rf*t)
    local u2 = u_ac*math.cos(w_rf*t)
    local voltages = {u1,u2,0,u2,0,u1}
    local a = Field(voltages,r[k1],field_coeff)
    local k2 = k % 2 + 1
    for i=1,n_ions do
        r[k2][i] = {}
        v[k2][i] = {}
        for j=i+1,n_ions do
            local dx = r[k1][i][1] - r[k1][j][1]
            local dy = r[k1][i][2] - r[k1][j][2]
            local dz = r[k1][i][3] - r[k1][j][3]
            local dist = (dx^2 + dy^2 + dz^2)^1.5
            dx = dx * coulomb_coeff / dist
            dy = dy * coulomb_coeff / dist
            dz = dz * coulomb_coeff / dist
            a[i][1] = a[i][1] + dx
            a[i][2] = a[i][2] + dy
            a[i][3] = a[i][3] + dz
            a[j][1] = a[j][1] - dx
            a[j][2] = a[j][2] - dy
            a[j][3] = a[j][3] - dz
        end
        for j=1,3 do
            r[k2][i][j] = r[k1][i][j] + v[k1][i][j] * dt
            v[k2][i][j] = v[k1][i][j] + a[i][j] * dt
        end
    end 
end
print(os.clock() - tic)

tic = os.clock()
for i=1,n_ions do
    disp.plot(plt[i], {win = win[i]})
end
print(os.clock() - tic)