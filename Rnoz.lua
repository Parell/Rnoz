local R = 8314.46261815324 --universal gas constant
local g0 = 9.80665 --standard gravity

local P0 = 7000000
local Pe = 102325
local Pamb = 101325
local At = 0.672
local Ae = 10.75
local T0 = 3558.24
local k = 1.22
local M = 22.186

local network = {
    {name = "grauna",  IP = "210.26.30.34"},
    {name = "arraial", IP = "210.26.30.23"},
    {name = "lua",     IP = "210.26.23.12"},
    {name = "derain",  IP = "210.26.23.20"},
  }

local engine = {P0, Pe, Pamb, At, Ae, T0, k, M}

Set()


local Rgas = R / M

local mdot = ((At * P0) / math.sqrt(T0)) * math.sqrt(k/Rgas) * ((k + 1)/2)^ -((k + 1)/(k - 1)/2)
local ve = math.sqrt(((T0 * R) / M) * ((2 * k) / (k - 1)) * (1 - (Pe / P0) ^ ((k - 1) / k)))

local thrust = mdot * ve + (Pe - Pamb) * Ae


local isp = thrust / (mdot * g0)
local twr = thrust / (2965000 + 140000)

--local deltaV = isp * g0 * math.log(((123000 + 49620 + 2290000) + 140000) / (175300 + 140000))

local a = math.sqrt(k * Rgas * T0)
local mache = ve / a


function Set()
    P0 = 2
end







print(Rgas .. "   " .. mdot .. "   " .. ve .. "   " .. thrust)

print("\n" .. isp .. "   " .. twr .. "   " .. mache)

tables = 2

print("\n" )


--note: after one engine stats are calculated the next on calculated until complete and final delta v stats are posted

--add units to all

--add license and logo in Ascii