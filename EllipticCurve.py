from GaloisField import *


class Point:
    def __init__(self, x=0, y=0, is_neutral=False):
        """
        Initializes the point

        :param x: x coordinate of the point
        :param y: y coordinate of the point
        :param is_neutral: is point the neutral element (boolean indicator)
        """
        self.is_neutral = is_neutral
        self.x = x
        self.y = y

    def __eq__(self, other):
        """
        :param other: point to check equivalence
        :return: boolean indicator of points equivalence
        """
        x_equivalence = self.x == other.x
        y_equivalence = self.y == other.y
        equivalence = x_equivalence and y_equivalence
        return equivalence

    def __str__(self):
        """
        :return: simple string representation for points
        """
        if self.is_neutral:
            return '(neutral element)'
        else:
            return '(' + str(self.x) + ', ' + str(self.y) + ')'


class EllipticCurve:
    def __init__(self, prime_number, x1_coefficient, x0_coefficient):
        """
        Initializes the elliptic curve y ** 2 = x ** 3 + a * x + b (mod p)

        :param prime_number: characteristic of galois field
        :param x1_coefficient: a from elliptic curve equation
        :param x0_coefficient: b from elliptic curve equation
        """
        self.galois_field = GaloisField(prime_number)
        self.x1_coefficient = self.galois_field.element(x1_coefficient)
        self.x0_coefficient = self.galois_field.element(x0_coefficient)

        discriminant = self.__calc_discriminant()
        assert discriminant != 0

    def is_point_on_curve(self, point):
        """
        Checks if point is on the curve

        :param point: point to check
        :return: is point on the curve (boolean indicator)
        """
        if point.is_neutral:
            return True

        return self.are_coordinates_on_curve(point.x, point.y)

    def are_coordinates_on_curve(self, x, y):
        """
        Checks if coordinates are on the curve

        :param x: x coordinate to check
        :param y: y coordinate to check
        :return: are coordinates on the curve (boolean indicator)
        """
        x = self.galois_field.element(x)
        y = self.galois_field.element(y)
        right = self.__calc_elliptic_right(x)
        left = self.__calc_elliptic_left(y)
        return left == right

    def add(self, point1, point2):
        """
        Applies the addition operation
        to two points of the elliptic curve

        :param point1: 1st point
        :param point2: 2nd point
        :return: sum (elliptic curve point)
        """
        assert self.is_point_on_curve(point1)
        assert self.is_point_on_curve(point2)

        x1 = self.galois_field.element(point1.x)
        y1 = self.galois_field.element(point1.y)
        x2 = self.galois_field.element(point2.x)
        y2 = self.galois_field.element(point2.y)

        if point1.is_neutral:
            return Point(x2, y2)
        if point2.is_neutral:
            return Point(x1, y1)

        if x1 != x2:
            k = self.__calc_k_different_points(x1, y1, x2, y2)
        elif y1 != y2:
            return Point(is_neutral=True)
        else:
            k = self.__calc_k_equal_points(x1, y1)
        c = self.__calc_c(k, x1, y1)
        x3 = self.__calc_x3(k, x1, x2)
        y3 = self.__calc_y3(k, x3, c)
        y3_inverse = self.galois_field.add_inverse(y3)

        point_result = Point(x3, y3_inverse)

        assert self.is_point_on_curve(point_result)

        return point_result

    def mul(self, coefficient, point):
        """
        Applies the multiplication operation
        to two points of the elliptic curve

        :param coefficient: scalar factor
        :param point: point factor
        :return: product (elliptic curve point)
        """
        assert coefficient >= 0

        x = self.galois_field.element(point.x)
        y = self.galois_field.element(point.y)

        point_result = Point(is_neutral=True)
        current_term = Point(x, y)
        current_degree = 1
        while coefficient != 0:
            coefficient_potential = coefficient - current_degree
            next_degree = current_degree * 2
            if coefficient_potential % next_degree == 0:
                coefficient = coefficient_potential
                point_result = self.add(point_result, current_term)
            current_term = self.add(current_term, current_term)
            current_degree = next_degree
        return point_result

    def order(self, point):
        """
        Calculates order of the elliptic curve point

        :param point: point
        :return: order of the point
        """
        counter = 1
        point_res = Point(point.x, point.y)
        while not point_res.is_neutral:
            point_res = self.add(point_res, point)
            counter += 1
        return counter

    def point_generator(self):
        """
        Simple generator for points on the curve
        """
        for i in range(self.galois_field.p):
            for j in range(self.galois_field.p):
                if self.are_coordinates_on_curve(i, j):
                    yield Point(i, j)

    def __calc_discriminant(self):
        """
        4 * x1_coefficient ** 3 + 27 * x0_coefficient ** 2
        """
        return \
            self.galois_field.add(
                self.galois_field.mul(
                    4, self.galois_field.pow(self.x1_coefficient, 3)
                ),
                self.galois_field.mul(
                    27, self.galois_field.pow(self.x0_coefficient, 2)
                )
            )

    def __calc_elliptic_left(self, y):
        """
        y ** y
        """
        y = self.galois_field.element(y)
        return \
            self.galois_field.mul(y, y)

    def __calc_elliptic_right(self, x):
        """
        x ** 3 + x1_coefficient * x + x0_coefficient
        """
        x = self.galois_field.element(x)
        return \
            self.galois_field.add(
                self.galois_field.add(
                    self.galois_field.pow(x, 3),
                    self.galois_field.mul(self.x1_coefficient, x)
                ),
                self.x0_coefficient
            )

    def __calc_k_different_points(self, x1, y1, x2, y2):
        """
        (y2 - y1) / (x2 - x1)
        """
        return \
            self.galois_field.div(
                self.galois_field.sub(y2, y1),
                self.galois_field.sub(x2, x1)
            )

    def __calc_k_equal_points(self, x, y):
        """
        (3 * x ** 2 + x1_coefficient) / (2 * y)
        """
        return \
            self.galois_field.div(
                self.galois_field.add(
                    self.galois_field.mul(
                        3, self.galois_field.pow(x, 2)
                    ),
                    self.x1_coefficient
                ),
                self.galois_field.mul(2, y)
            )

    def __calc_c(self, k, x, y):
        """
        y - k * x
        """
        return \
            self.galois_field.sub(
                y, self.galois_field.mul(k, x)
            )

    def __calc_x3(self, k, x1, x2):
        """
        k ** k - (x1 + x2)
        """
        return \
            self.galois_field.sub(
                self.galois_field.pow(k, 2),
                self.galois_field.add(x1, x2)
            )

    def __calc_y3(self, k, x3, c):
        """
        k * x3 + c
        """
        return \
            self.galois_field.add(
                self.galois_field.mul(k, x3), c
            )
