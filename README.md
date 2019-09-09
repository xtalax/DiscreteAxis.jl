# DiscreteAxis
A small utility module containing some useful types and methods for dealing with axies and collections of axies, for use with discretized functions stored in `AbstractArrays`. They can be iterated over and indexed in to like normal arrays, and have math operations `+, -, *,\` defined for `Number` types.

#Usage
This package implements a wrapper around a `StepRangeLen`, making it convenient and concise to retrieve important information about its step size and length. This was mostly implemented since I find that it makes code more readable.

It also implements a `Space2D/3D` which have the fields `x, y`\ `z`, each containing a DiscreteAxis. These make it convenient to pass around multiple axies when working on multidimensional code.

#Create a Linear Axis
```
x = LinearAxis(xmin, xmax, Δx) # compensating axis to the new units
y = LinearAxis(ymin, ymax, Δy)
> x.N == length(x)
true
> y.Δ == abs(y[2] - y[1])
true

space = Space(x,y,z) #Coming soon, space with any number of arguments and configurable names

>size(space) == (s.x.N, s.y.N, s.z.N)
true

midpoint = [div(x.N,2) for x in space]
```
#Coming soon:
- Iterate over all the points in a space
- Broadcast over a space, returning an Array of the correct shape
- A `DiscreteFunction` Type, containing a discretized function and its axies
