using Statistics,DataFrames,CSV,Plots,Images,LsqFit,Symbolics
plotlyjs();

## Hluti 1
hvec = [31 36 32 32 31.5 33 35 31 35 36 30].*1e-2;
h = mean(hvec);
dh = 0.5e-2;
l = 122e-2;
dl = 0.5e-2;

@variables sh sl shErr slErr
sμ = sh/sl;
sdμ = Symbolics.findErrorFromSym(sμ)
μ = substitute(sμ, Dict(sh => h, sl => l))
dμ = substitute(sdμ, Dict(shErr => dh, sh => h, sl => l, slErr => dl))

## Hluti 2

# Read data
data = [Matrix(CSV.read("data" * string(i) * "h" * string(j) * ".txt",DataFrame)) for i in 1:3, j in 1:2]
# Cut out trash data
data[1,1] = data[1,1][1:end-35,:];


# Find extrema that aren't too small so we don't get the noisy bit at the end
segsacc = [
	[
		getindex.(findlocalminima(data[i,j][:,4] .* (data[i,j][:,4] .<  -0.7)),1) 
		for j in 1:2
	] 
	for i in 1:3
];
segsmax = [
	[
		getindex.(findlocalmaxima(data[i,j][:,3] .* (data[i,j][:,3] .>  0.05)),1) 
		for j in 1:2
	] 
	for i in 1:3
];
segsmin = [
	[
		getindex.(findlocalminima(data[i,j][:,3] .* (data[i,j][:,3] .< -0.05)),1) 
		for j in 1:2
	] 
	for i in 1:3
];

# Function for cleaning extrema
function splitgroups(v,sepp)
	start = firstindex(v)
	groups = typeof(v)[]
	for stop in [findall(diff(v) .> sepp); lastindex(v)]
		push!(groups, @view(v[start:stop]))
		start = stop + 1
	end
	groups
end

# Clean extrema so only one for each peak
seppsmax = [[39, 38] ,[39, 35], [39, 39]];
segsmax = [
	[
		maximum.(j) for j in splitgroups.(segsmax[i],seppsmax[i])
	] 
	for i in 1:3
];
segsmin = [[minimum.(j) for j in splitgroups.(i,30)] for i in segsmin];

# Clean singletons
for i in 1:3
	for j in 1:2
		if segsmin[i][j][end] > segsmax[i][j][end]
			pop!(segsmin[i][j])
		end
	end
end

# clean bad pairs
popat!(segsmin[1][2],1)
popat!(segsmax[2][2],8)
popat!(segsmin[2][2],8)
popat!(segsmin[2][2],7)


##

# make segment ranges from minima to maxima
segpairs = [[[(k:segsmax[i][j][segsmax[i][j] .> k][1]) for k in segsmin[i][j]] for j in 1:2] for i in 1:3];

# split them in twain about naught
for i in eachindex(segpairs)
	for j in eachindex(segpairs[i])
		newsegvec = [];
		for k in eachindex(segpairs[i][j])
			for m in segpairs[i][j][k]
				if data[i,j][m,3] <= 0 && data[i,j][m+1,3] >= 0
					push!(newsegvec, segpairs[i][j][k][1]:m, (m+1):segpairs[i][j][k][end])
				end
			end
		end
		segpairs[i][j] = newsegvec;
	end
end

##

f(x,p) = p[1].*x.+p[2];
polynoms = [
	[
		[
			coef(curve_fit(f, data[i,j][k,1], data[i,j][k,3], [0.5, 0.5]))
			for k in segpairs[i][j]
		] for j in eachindex(segpairs[i])
	] for i in eachindex(segpairs)
];


for i in 1:3, j in 1:2
	plot(data[i,j][:,1],data[i,j][:,3],labels=:none)
	for k in eachindex(segpairs[i][j])
		x = data[i,j][segpairs[i][j][k],1]
		plot!(x, f(x,polynoms[i][j][k]),linewidth=3,color=:red,labels=:none)
	end
	# scatter!(segsmax[i][j],data[i,j][:,3][segsmax[i][j]],color = :blue,labels=:none)
	# scatter!(segsmin[i][j],data[i,j][:,3][segsmin[i][j]],color = :yellow,labels=:none)
	# scatter!(segsacc[i][j],data[i,j][:,3][segsacc[i][j]],color = :red)
	savefig(string("plt/fig",i,"-",j))
end

g = 9.8
h2 = [4.5, 6.5]*1e-2;
θ = [asin(h2[1]/l), asin(h2[2]/l)];

calcμ(p,θ) = (p[1][1]-p[2][1])/(2g*cos(θ))
μ_r = [
	[
		[
			calcμ([polynoms[i][j][k-1], polynoms[i][j][k]], θ[j])
			for k in eachindex(polynoms[i][j]) if k % 2 == 0
		] for j in eachindex(polynoms[i])
	] for i in eachindex(polynoms)
];
