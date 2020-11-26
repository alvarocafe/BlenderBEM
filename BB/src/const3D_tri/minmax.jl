function mini(A)
min = Inf;
for i = 1:size(A,1)
	if abs(A[i,1])<min
		min = A[i,1];
	end
end
return mini
end

function maxi(A)
maxi = 0.0;
for i = 1:size(A,1)
	if abs(A[i,1])>maxi
		min = A[i,1];
	end
end
return maxi
end

