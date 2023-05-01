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

    # a = (1.0, 2.0, NaN,)
    # b = ((3.0), (5.0), (7.0))
    # println(b.+a)
    # println(typeof(a))
    # print(isnan.(a))

    # f1(x) = x
    # i = 1
    # a = :f$1
    # print()
    # b = [1, 2, 3]
    # a = [[1 for i in 1:5] for j in 1:2]
    # println(typeof(a))
    # print([ (begin
    #     a = 1 + 1
    #     b = a + 1
    #     return a
    # end) for i in 1:5])
    # for i in 1:1
    #     print(1)
    # end

    # a = [1 for i in 1:5, j in 1:8, k in 1:10]
    # println(a)
    # println(a[1, 1])
end

test()