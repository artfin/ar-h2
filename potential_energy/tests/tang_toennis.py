#!/usr/bin/python3

from numpy import exp, arange, int32, array

# The pre-computed factorial numbers.
factorials = {
	2 : array([1, 1, 2]),
	3 : array([1, 1, 2, 6]),
	4 : array([1, 1, 2, 6, 24]),
	5 : array([1, 1, 2, 6, 24, 120]),
	6 : array([1, 1, 2, 6, 24, 120, 720]),
	7 : array([1, 1, 2, 6, 24, 120, 720, 5040]),
	8 : array([1, 1, 2, 6, 24, 120, 720, 5040, 40320])
}


def tang_toennies(n, r, E):
	""" The Tang-Toennies damping function. The equation is:
		f_{n}(r) = 1.0 - \exp{(-Er)}\sum_{i=0}^{k}{\frac{(Er)^k}{k!}}
	where n is the order, r is the damping target and E is the coefficient.
	Parameters
	----------
	n : int
		The order of the damping function.
	r : float
		The damping target.
	E : float
		The damping coefficient.
	Returns
	-------
	v : float
		The damping result.
	"""
	Er = E * r
	k = arange(0, n+1, dtype=int32)
	return 1.0 - exp(-Er) * ((Er**k / factorials[n]).sum())

if __name__ == '__main__':
    print("Value of tang_toennis damping function: ", tang_toennies(8, 8.0, 1.2))
