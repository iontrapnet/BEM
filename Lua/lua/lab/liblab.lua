

function lab.new(...)
   local x;
   if type(arg[1])=="table" then
      x=torch.Tensor(#arg,#arg[1])
      for i=1,#arg do
	 for j=1,#arg[1] do
	    x[i][j]=arg[i][j];
	 end
      end
   else
      x=torch.Tensor(#arg)
      for i=1,#arg do
	 x[i]=arg[i];
      end
   end
   return x;
end
