require "CLRPackage"
require "CLRForm"
import "System"
import "System.Windows.Forms"
import "System.Drawing"
import 'NPlot.dll'

local NPW = CLRPackage('NPlot.dll','NPlot.Windows')

local function doubles(t)
    return luanet.make_array(Double,t)
end

local colors = {Color.Black, Color.Blue, Color.Brown, Color,Cyan, Color.Gold, Color.Gray, Color.Green, Color.Orange, Color.Pink, Color.Purple, Color.Red, Color.Silver, Color.Wheat, Color.Yellow}

local function lineplot(s,x,y,i)
    x = doubles(x)
    y = doubles(y)
    local c = colors[i]
    local lp = LinePlot()
    lp.AbscissaData = x
    lp.DataSource = y
    lp.Color = c
    lp.Label = i --nil
    s:Add(lp)
    --[[local pp = PointPlot()
    pp.AbscissaData = x
    pp.DataSource = y
    pp.Marker.Type = Marker.MarkerType.FilledCircle
    pp.Marker.Color = c
    pp.Label = i
    s:Add(pp)]]
end

local function plot(...)
    local f = Form()
    local s = NPW.PlotSurface2D()
    s.Dock = DockStyle.Fill
    f.Controls:Add(s)
    local args = {...}
    for i=1,#args/2 do
        lineplot(s,args[2*i-1],args[2*i],i)
    end
    local legend = Legend()
	legend.VerticalEdgePlacement = Legend.Placement.Inside
	legend.HorizontalEdgePlacement = Legend.Placement.Outside
	legend.BorderStyle = LegendBase.BorderType.Line
	s.Legend = legend
    f:ShowDialog()
end

return plot