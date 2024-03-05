#Julia 1.10
#conjugate gradient calculating script
#C Mancuso 2019

using LinearAlgebra

function return_m(A, y, iterations)

    M, N = size(A)[2], size(A)[1]
    x = zeros(M , 1) #starting values

    s = y - (A * x)
    p = A' * s
	r = p
	q = A * p
	old = r' * r

	for i = 1:iterations
		α =((r' * r) / (q' * q))
		x = x +  (p * α)
		s = s - (q * α)
		r = A' * s
		new = r' * r
		β = new / old
		old = new
		p = r + (p * β)
		q = A * p
	end
    return x

end
