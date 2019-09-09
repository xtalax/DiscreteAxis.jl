module DiscreteAxis
export LinearAxis, LogAxis, FuncAxis, CoordinateSystem, interior_range
export DAxis, Space3D, Space2D, ProductAxis
export flatten, domain

    using StaticArrays, Rotations
    import Base.iterate
    import Base.convert
    import Base.length
    import Base.getindex
    import Base./
    import Base.*
    import Base.+
    import Base.-
    import Base.size
    import Base.first
    import Base.last
    abstract type DAxis{T} end
    abstract type ProductAxis{T,N} <: FieldVector{N, DAxis{T}} end


    struct LinearAxis{T<:Number} <: DAxis{T}
        pts::AbstractVector{T}
        Δ::T
        N::Int64
        i::typeof(1:10)
        function LinearAxis(min::T,max::T, N::Int64) where {T<:Number}
            points = range(min, stop=max, length=N)
            Δ = abs(points[1]-points[2])

            new{T}(points,Δ,N,1:N)
        end
        function LinearAxis(min::T, max::T, Δ::T) where {T<:Number}
            points = range(min, stop = max,  step = Δ)
            N = length(points)
            new{T}(points,Δ,N,1:N)
        end
        LinearAxis(pts, Δ, N) = new{typeof(Δ)}(pts, Δ, N, 1:N)
    end
    gettype(axis::DAxis{T}) where {T} = T
    Base.first(x::LinearAxis) = x.pts[1]
    Base.last(x::LinearAxis) = x.pts[end]

    struct Space3D{T} <: ProductAxis{T,3}
        x::LinearAxis{T}
        y::LinearAxis{T}
        z::LinearAxis{T}
        Space3D(x::LinearAxis{T}, y::LinearAxis{T}, z::LinearAxis{T})  where {T<:Number} = new{T}(x,y,z)
    end
    struct Space2D{T} <: ProductAxis{T,2}
        x::LinearAxis{T}
        y::LinearAxis{T}
        Space2D(x::LinearAxis{T}, y::LinearAxis{T})  where {T<:Number} = new{T}(x,y)
    end
    struct DiscreteSpace{T<:Number, N}
        dims::SVector{N, LinearAxis{T}}
        dimtype::SVector{N, String}
        DiscreteSpace(Axies...) = new{gettype(Axies[1]), length(axies)}(SVector(Axies...))
    end

    function interior_range(space::ProductAxis, nPML::Int)
        return Tuple([nPML+1:x.N-nPML for x in space])
    end


    gettype(s::ProductAxis{T}) where {T} = T
    Base.size(s::Space3D) = (s.x.N, s.y.N, s.z.N)
    Base.size(s::Space2D) = (s.x.N, s.y.N)

    function originindex(s::Space3D)
        y = [0, 0, 0]
        for (i,x) in space
            index, minval = findmin(abs.(x.pts))
            if minval < x.Δ/2
                y[i] = index
            else
                throw(ArgumentError("Origin not contained in space"))
            end
        end
        return y
    end

    #Base.getindex(s::Space3D, i::Int, j::Int, k::Int) = [s.x.pts[i], s.y.pts[j],s.z.pts[k]]
    #Base.getindex(s::Space3D, I::CartesianIndex) = [s.x.pts[I[1]], s.y.pts[I[2]], s.z.pts[I[3]]]
    Base.CartesianIndices(s::ProductAxis) = CartesianIndices((x.i for x in s))

    function transform_flat(s::Space3D; origin = zeros(gettype(s), 3), azimuth = 0.0, elevation = 0.0)
        T = gettype(s)
        R = Matrix(RotZXY(azimuth, elevation, 0.0)) #generate rotation matrix
        flat = fill(zeros(T,3), size(s)...)
        for (i,x) in enumerate(s.x.pts), (j,y) in enumerate(s.y.pts), (k,z) in enumerate(s.z.pts)
            flat[i,j,k] = R * ([x,y,z] .- origin)
        end
        return flat::Array{returntype, 3}
    end




    iterate(x::LinearAxis, n::Int) = iterate(x.pts,n)
    iterate(x::LinearAxis) = iterate(x.pts)
    length(x::LinearAxis) = x.N

    getindex(x::LinearAxis, n::Int) = getindex(x.pts, n)
    /(x::LinearAxis, a) = LinearAxis(x.pts ./ a, x.Δ/a, x.N)
    *(x::LinearAxis, a) = LinearAxis(x.pts .* a, x.Δ*a, x.N)
    *(a, x::LinearAxis) = *(x, a)
    +(x::LinearAxis, a) = LinearAxis(x.pts .+ a, x.Δ, x.N)
    +(a, x::LinearAxis) = +(x::LinearAxis, a)
    -(x::LinearAxis, a) = LinearAxis(x.pts .- a, x.Δ, x.N)
    -(a, x::LinearAxis) = -(x::LinearAxis, a)

    function domain(y::DAxis)
        maximum(y.pts) - minimum(y.pts)
    end

    function domain(ymin,ymax)
        ymax-ymin
    end


#=
    struct LogAxis <: DAxis
        pts::SVector{Number}
        Δ::SVector{Number}
        N::Int
        base::String
        function LogAxis(min, max, N::Int, base::String)
            pts = logspace(min,max,N)
            Δ = SVector[points[i+1]-points[i] for i in 1:(N-1)]
            new(points,Δ,N,base)
        end
    end

    struct FuncAxis <: Axis
        pts::SVector{Number}
        Δ::SVector{Number}
        N::Int
        func::function
        inverse::function

        function FuncAxis(min, max, N::Int, func, inverse)
            SVector(func.(range(inverse(min),stop=inverse(max),N)))
            Δ = SVector[points[i+1]-points[i] for i in 1:(N-1)]
            new(points,Δ,N,base)
        end
    end

    function convert(T, x::DAxis)
        convert(T, x.pts)
    end

    struct CoordinateSystem
        axes::SVector{DAxis}
        metric::SVector{Char}
        CoordinateSystem(axes, metric) = new(axes,metric)
    end
    =#
end
