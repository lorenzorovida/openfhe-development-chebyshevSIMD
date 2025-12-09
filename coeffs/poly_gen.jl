using Polynomials

n = 2^13

inpt = rand(-1:0.000001:1, n)

open("input.txt", "w") do f
    for x in inpt
        println(f, x)
    end
end

coeffs = []

for i in 1:n
   push!(coeffs, rand(-0.01:0.0000001:0.01, 495))
end

for i in 1:n
    open("coeffs/p$(i).txt", "w") do f
        for x in coeffs[i]
            println(f, x)
        end
    end
end

for i in 1:n
    p = ChebyshevT(coeffs[i])
    sum[i] = p(inpt[i])
end

open("expected.txt", "w") do f
    for x in sum
        println(f, x)
    end
end
