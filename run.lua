function string:split(delimiter)
    local result = {}
    for part in self:gmatch("[^"..delimiter.."]+") do
        result[#result+1] = part
    end
    return result
end

local function read_csv(path)
  local f = io.open(path)
  local rows = {}
  for line in f:lines() do
    local row = {}
    for _, cell in ipairs(line:split(',')) do
      row[#row + 1] = tonumber(cell)
    end
    rows[#rows + 1] = row
  end
  f:close()
  return rows
end

local function preprocessing(nodes, elems, title)
  local alpha = 1
  local beta = -1
  local f = io.open('input.dat', 'w')
  f:write(title..'\n')
  f:write(string.format('%d   %d   %d\n', #elems, #nodes, 0))
  f:write(string.format('%10.6e  %10.6e\n', alpha, beta))
  f:write('#Nodes:\n')
  for i, node in ipairs(nodes) do
        f:write(string.format(' %d     %20.12e     %20.12e  %20.12e  \n', i, node[1], node[2], node[3]))
  end
  f:write('#Elements:\n')
  for i, elem in ipairs(elems) do
        f:write(string.format(' %d     %d    %d     %d  %d  %20.12e \n', i, elem[1], elem[2], elem[3], elem[4], elem[5]))
  end
  f:close()
end

local argv = {...}
local argc = #argv
title = 'input'
if argc > 0 then title = argv[1] end
local nodes = read_csv'nodes.csv'
local elements = read_csv'elements.csv'
preprocessing(nodes, elements, title)