import sympy as sp


class DerivativeCalculator:
    def __init__(self, expression, variable):
        self.expression = expression
        self.variable = variable

    def calculate_derivative(self, order=1):
        x = sp.symbols(self.variable)
        expr = sp.sympify(self.expression)
        derivative = sp.diff(expr, x, order)
        return derivative

    def find_roots(self):
        x = sp.symbols(self.variable)
        derivative = self.calculate_derivative()
        roots = sp.solve(derivative, x)
        numeric_roots = [float(root.evalf()) for root in roots if root.is_real]
        return numeric_roots


class FunctionSignDeterminer:
    def __init__(self, expression, variable, roots):
        self.expression = expression
        self.variable = variable
        self.roots = roots

    def evaluate_and_display_signs(self, x_values):
        x = sp.symbols(self.variable)
        expr = sp.sympify(self.expression)

        header = ["x", "Sign"]
        table = [header]

        for x_val in x_values:
            try:
                if x_val == sp.oo:
                    value_at_x = sp.limit(expr, x, x_val, dir='+')
                else:
                    value_at_x = expr.subs(x, x_val)
                sign_at_x = self.determine_sign(value_at_x)
            except TypeError:
                sign_at_x = "Невизначений (NaN)"
            table.append([x_val, sign_at_x])

        self.display_table(table)

    def determine_sign(self, value):
        if value > 0:
            return "(+)"
        elif value < 0:
            return "(-)"
        else:
            return "Нуль (0)"

    def display_table(self, table):
        col_width = max(len(str(word)) for row in table for word in row) + 2
        for row in table:
            print("".join(str(word).ljust(col_width) for word in row))


class NewtonMethod:
    def __init__(self, expression, variable, initial_guess, tolerance=1e-6):
        self.expression = expression
        self.variable = variable
        self.initial_guess = initial_guess
        self.tolerance = tolerance
        self.x = sp.symbols(variable)
        self.first_derivative = DerivativeCalculator(expression, variable)

    def run(self, max_iterations=100):
        expr = sp.sympify(self.expression)

        x_n = self.initial_guess
        x_n_1 = x_n - 2 * self.tolerance

        iterations = 0
        while iterations < max_iterations and abs(x_n - x_n_1) >= self.tolerance:
            f_x_n = expr.subs(self.x, x_n)
            f_prime_x_n = self.first_derivative.calculate_derivative().subs(self.x, x_n)

            x_n_1 = x_n
            x_n = x_n - f_x_n / f_prime_x_n

            iterations += 1

        return x_n


class BisectionMethod:
    def __init__(self, expression, variable, interval, tolerance=1e-6):
        self.expression = expression
        self.variable = variable
        self.interval = interval
        self.tolerance = tolerance
        self.x = sp.symbols(variable)

    def run(self, max_iterations=100):
        expr = sp.sympify(self.expression)
        a, b = self.interval

        if a >= b:
            raise ValueError("Invalid interval: a must be less than b.")

        f_a = expr.subs(self.x, a)
        f_b = expr.subs(self.x, b)

        if f_a * f_b > 0:
            raise ValueError("The function values at the interval endpoints must have opposite signs.")

        iterations = 0
        while iterations < max_iterations and abs(b - a) >= 2 * self.tolerance:
            c = (a + b) / 2
            f_c = expr.subs(self.x, c)

            if abs(f_c) < self.tolerance:
                return c

            if f_c * f_a < 0:
                b = c
            else:
                a = c

            iterations += 1

        return (a + b) / 2


# Введення рівняння та змінної для диференціювання
expression = input("Введіть рівняння: ")
variable = input("Введіть змінну: ")

# Визначення похідної та вивід її коренів
calculator = DerivativeCalculator(expression, variable)
calculator_result = calculator.calculate_derivative()
print(f"Похідна від {expression} по {variable}: {calculator_result}")

roots = calculator.find_roots()
print(f"Корені похідної: {roots}")

# Знаходження знаків функції
sign_determiner = FunctionSignDeterminer(expression, variable, roots)
x_values_to_evaluate = [-sp.oo] + sorted(sign_determiner.roots) + [sp.oo]
sign_determiner.evaluate_and_display_signs(x_values_to_evaluate)

# Зменшення проміжку
neg_inf_value = float(input("Введіть значення для -inf: "))
pos_inf_value = float(input("Введіть значення для +inf: "))

# Визначення інтервалів для коренів
roots_intervals = []
for root in sign_determiner.roots:
    if root <= 0:
        roots_intervals.append((neg_inf_value, 0))
    else:
        roots_intervals.append((9, pos_inf_value))

# Вивід інтервалів для коренів
for i, interval in enumerate(roots_intervals, start=1):
    print(f"Корінь x_{i} належить інтервалу: {interval}")

# Метод Ньютона для уточнення першого кореня
initial_guess_x1 = roots_intervals[0][0]
tolerance_newton = float(input("Введіть точність для методу Ньютона: "))
newton_method_x1 = NewtonMethod(expression, variable, initial_guess_x1, tolerance_newton)
result_x1 = newton_method_x1.run()
print(f"Уточнений x_1: {round(result_x1, 4)}")

# Метод бісекцій для уточнення другого кореня
interval_x2 = roots_intervals[1]
tolerance_bisection = float(input("Введіть точність для методу бісекції: "))
bisection_method_x2 = BisectionMethod(expression, variable, interval_x2, tolerance_bisection)
result_x2 = bisection_method_x2.run()
print(f"Уточнений x_2: {round(result_x2, 4)}")