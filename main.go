package main

import "fmt"

func main() {
	fmt.Println(SetF("x**3+4*x**2-10"))
	fmt.Println(SetDF("3*x**2+8*x"))
	fmt.Println(SetD2F("6*x+8"))
	fmt.Println(SetG("sqrt(10/(x+4))"))
	tol := float64(1e-15)
	fmt.Println(IncrementalSearch(0, 0.1, 50))
	fmt.Println(Bisection(1, 2, tol, 50))
	fmt.Println(bisTable)
	fmt.Println()
	fmt.Println(FalsePosition(1, 2, tol, 50))
	fmt.Println(falTable)
	fmt.Println()
	fmt.Println(FixedPoint(1, tol, 100))
	fmt.Println(fixTable)
	fmt.Println()
	fmt.Println(Newton(1, tol, 50))
	fmt.Println(newTable)
	fmt.Println()
	fmt.Println(Secant(1, 2, tol, 50))
	fmt.Println(secTable)
	fmt.Println()
	fmt.Println(MultipeRoot(1, tol, 50))
	fmt.Println(mulTable)
}
