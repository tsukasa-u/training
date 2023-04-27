function test()
    # a::Array{Array{Int, 1}, 1} = Array{Int, 1}[]
    # push!(a, Real[1])
    # push!(a, Real[1, 2])
    # push!(a, Real[1, 2, 3])
    # # println(typeof(a))
    # # println(typeof(a[1]))
    # println(a[2][2])

    # f(x, y) = x + y
    # g(x, y) = x - y
    # a = (f, g)
    # print(a[1](1, 2))

    a = (1.0, 1.0, NaN,)
    println(a)
    println(typeof(a))
end

test()