--
-- A torch client for `display` graphics server
-- Based heavily on https://github.com/clementfarabet/gfx.js/blob/master/clients/torch/js.lua
--

local mime = require 'mime'
local http = require 'socket.http'
local ltn12 = require 'ltn12'
local json = require 'cjson'
local ffi = require 'ffi'
local torch = require 'torch'

M = {
  url = 'http://localhost:8000/events'
}

local function uid()
  return 'pane_' .. (os.time() .. math.random()):gsub('%.', '')
end

local function send(command)
  -- TODO: make this asynchronous, don't care about result, but don't want to block execution
  command = json.encode(command)
  http.request({
    url = M.url,
    method = 'POST',
    headers = { ['content-length'] = #command, ['content-type'] = 'application/text' },
    source = ltn12.source.string(command),
  })
end

local function pane(type, win, title, content)
  win = win or uid()
  send({ command='pane', type=type, id=win, title=title, content=content })
  return win
end

-- Set the URL of the listening server
function M.configure(config)
  local port = config.port or 8000
  local hostname = config.hostname or '127.0.0.1'
  M.url = 'http://' .. hostname .. ':' .. port ..'/events'
end

-- data is either a 2-d torch.Tensor, or a list of lists
-- opts.labels is a list of series names, e.g.
-- plot({ { 1, 23 }, { 2, 12 } }, { labels={'iteration', 'score'} })
-- first series is always the X-axis
-- See http://dygraphs.com/options.html for supported options
function M.plot(data, opts)
  opts = opts or {}

  local dataset = {}
  if torch.typename(data) then
    for i = 1, data:size(1) do
      local row = {}
      for j = 1, data:size(2) do
        table.insert(row, data[{i, j}])
      end
      table.insert(dataset, row)
    end
  else
    dataset = data
  end

  -- clone opts into options
  options = {}
  for k, v in pairs(opts) do
    options[k] = v
  end

  options.file = dataset
  if options.labels then
    options.xlabel = options.xlabel or options.labels[1]
  end

  -- Don't pass our options to dygraphs. 'title' is ok
  options.win = nil

  return pane('plot', opts.win, opts.title, options)
end

function M.text(text, opts)
  opts = opts or {}

  return pane('text', opts.win, opts.title, text)
end

return M
