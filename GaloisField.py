class GaloisField:
    def __init__(self, p):
        """
        Initializes the galois field GF(p)

        :param p: prime number (characteristic of galois field)
        """
        self.p = p
        assert is_prime(p)

    def element(self, a):
        """
        Converts a number to a galois field element

        :param a: number to convert
        :return equivalent galois field element
        """
        result = a % self.p
        return result

    def add(self, a, b):
        """
        Applies the addition operation
        to two elements of the galois field

        :param a: 1st addend
        :param b: 2nd addend
        :return sum (galois field element)
        """
        result = self.element(a + b)
        return result

    def sub(self, a, b):
        """
        Applies the subtraction operation
        to two elements of the galois field

        :param a: minuend
        :param b: subtrahend
        :return: difference (galois field element)
        """
        result = self.element(a - b)
        return result

    def mul(self, a, b):
        """
        Applies the multiplication operation
        to two elements of the galois field

        :param a: 1st factor
        :param b: 2nd factor
        :return: product (galois field element)
        """
        result = self.element(a * b)
        return result

    def div(self, a, b):
        """
        Applies the division operation
        to two elements of the galois field

        :param a: dividend
        :param b: divisor
        :return: quotient (galois field element)
        """
        b_inverse = self.mul_inverse(b)
        result = self.mul(a, b_inverse)
        return result

    def pow(self, a, b):
        """
        Applies the power operation
        to two elements of the galois field

        :param a: base
        :param b: power
        :return: value (galois field element)
        """
        if b < 0:
            a = self.mul_inverse(a)
            b = -b

        b %= self.p - 1

        result = 1
        current_multiplier = a
        current_degree = 1
        while b != 0:
            b_potential = b - current_degree
            next_degree = current_degree * 2
            if b_potential % next_degree == 0:
                b = b_potential
                result = self.mul(result, current_multiplier)
            current_multiplier = self.mul(current_multiplier, current_multiplier)
            current_degree = next_degree
        return result

    def mul_inverse(self, a):
        """
        Takes a multiplicative inverse of number

        :param a: number
        :return: multiplicative inverse (galois field element)
        """
        assert a != 0
        a = self.element(a)
        _, result, _ = extended_euclidean_algorithm(a, self.p)
        result = self.element(result)
        return result

    def add_inverse(self, a):
        """
        Takes a additive inverse of number

        :param a: number
        :return: additive inverse (galois field element)
        """
        result = self.sub(0, a)
        return result


def is_prime(a):
    """
    Checks if number is prime

    :param a: number to check
    :return: is the number prime (boolean indicator)
    """
    if a < 2:
        return False

    for i in range(2, int(a ** 0.5) + 1):
        if a % i == 0:
            return False
    return True


def extended_euclidean_algorithm(a, b):
    """
    Calculates the greatest common divisor and coefficients of
    the Diophantine equation ax + by = gcd(a, b)

    :param a: 1st number for gcd
    :param b: 2nd number for gcd
    :return: gcd(a, b), x, y : ax + by = gcd(a, b)
    """
    u = (a, 1, 0)
    v = (b, 0, 1)
    while v[0] != 0:
        q = u[0] // v[0]
        u, v = v, (u[0] - q * v[0], u[1] - q * v[1], u[2] - q * v[2])
    return u