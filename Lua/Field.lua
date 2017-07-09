local ffi = require 'ffi'
local mat = require 'matio.ffi'

local function loadmat(path, name)
    local f = mat.open(path, mat.ACC_RDONLY)
    local r = mat.varRead(f, name)
    mat.close(f)
    return r
end

PATH_SPHERE = '..\\Model\\sphere\\922019645250603132'
PATH_4ROD = '..\\Model\\4rod\\167634622912717531'

local FieldSession = {}

function FieldInit(path, xr, yr, zr, data)
    FieldSession['xr'] = xr
    FieldSession['yr'] = yr
    FieldSession['zr'] = zr
    local matfile = path..'-'..data..'.mat'
    FieldSession['pb'] = loadmat(matfile, 'pb')
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
    if pb.rank > 2 then s2 = pb.dims[2] end
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
local d0 = 500e-6
local d1 = 1000e-6
local r0 = ( d1 / math.sqrt(2) - d0 / 2)
local v0 = (math.sqrt(2 * 500 * kB / 3 / m))*1e-30

local u_dc = 100
local u_ac = 500
local w_drive = 2*math.pi*12e6
local dc_z = 0.122321e4
local w_z = (2 * dc_z * e * u_dc / m)^.5
local gap_z = (e^2 * Ke / (w_z^2 * m))^(1/3)

local n_ions = 5
local disp = require ('display')
local win = {}
local r = {{}}
local v = {{}}
for i=1,n_ions do
    local config = {title = 'Ion '..tostring(i), ylabel = 'Position (um)', labels = {"t (us)",'X','Y','Z'}}
    win[i] = disp.plot({}, config)
    r[1][i] = {}
    r[1][i][1] = (2 * math.random() - 1) * r0 / 200
    r[1][i][2] = (2 * math.random() - 1) * r0 / 200
    r[1][i][3] = (2 * math.random() - 1) * 2e-6 + (2*i - n_ions - 1) * gap_z / 2
    v[1][i] = {}
    v[1][i][1] = (2 * math.random() - 1) * v0
    v[1][i][2] = (2 * math.random() - 1) * v0
    v[1][i][3] = (2 * math.random() - 1) * v0
end

tic = os.clock()
local T_total = 1e-4                          --time of simulation
local dt = 8e-10
local n_T = math.floor(T_total / dt)
local t = {}
for k = 1,n_T do
    t[k] = dt * (k - 1)
    local u = u_ac*math.cos(w_drive*t[k])
    local voltages = {u_dc,u,0,u,0,u_dc}
    local a = Field(voltages,r[k],field_coeff)
    r[k+1] = {}
    v[k+1] = {}
    for i=1,n_ions do
        r[k+1][i] = {}
        v[k+1][i] = {}
        for j=i+1,n_ions do
            local dx = r[k][i][1] - r[k][j][1]
            local dy = r[k][i][2] - r[k][j][2]
            local dz = r[k][i][3] - r[k][j][3]
            local dist = (dx^2 + dy^2 + dz^2)^1.5
            --[[if dist == 0 then
                print(k,i,j,dx,dy,dz)
                return
            end]]
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
            r[k+1][i][j] = r[k][i][j] + v[k][i][j] * dt
            v[k+1][i][j] = v[k][i][j] + a[i][j] * dt
        end
    end
    --[[if math.mod(k,1e4) == 0 then
        for i=1,n_ions do
            local data = {}
            for j=1,k do
                data[j] = {1e6*t[j],1e6*r[j][i][1],1e6*r[j][i][2],1e6*r[j][i][3]} 
            end
            disp.plot(data, {win = win[i]})
        end
    end]]
end
print(os.clock() - tic)

tic = os.clock()
for i=1,n_ions do
    local data = {}
    for k=1,n_T do
        data[k] = {1e6*t[k],1e6*r[k][i][1],1e6*r[k][i][2],1e6*r[k][i][3]} 
    end
    disp.plot(data, {win = win[i]})
end
print(os.clock() - tic)

--[[local nplot = require 'nplot'
for i=1,n_ions do
    local x = {}
    local y = {}
    local z = {}
    for k=1,n_T do
        x[k] = r[k][1][1]
        y[k] = r[k][1][2]
        z[k] = r[k][1][3]
    end
    nplot(t,x,y,z)
end]]