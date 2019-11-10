package main

import "math"

type RootRange struct {
	a, b, root      float64
	isRoot, isRange bool
	iterations      uint
}

type RootAnswer struct {
	root                        float64
	isRoot, isAlmostRoot, badIn bool
	iterations                  uint
}

type Reg1 struct {
	xi, xs, xm, fxm float64
	it              int
}

type Reg2 struct {
	xm, fxm float64
	it      int
}

type Reg3 struct {
	xm, fxm, dfxm, d2fxm float64
	it                   int
}

type Reg4 struct {
	xm, fxm, dfxm float64
	it            int
}

var bisTable, falTable []Reg1
var fixTable, secTable []Reg2
var mulTable []Reg3
var newTable []Reg4

func IncrementalSearch(x0, dx float64, maxIt uint) (ans *RootRange) {
	ans = &RootRange{}
	xi := x0
	yi, _ := F(xi)
	if yi == 0 {
		ans.root = xi
		ans.isRoot = true
		return
	}
	xf := xi + dx
	yf, _ := F(xf)
	count := uint(1)

	for yi*yf > 0 && count <= maxIt {
		xi = xf
		yi = yf
		xf = xi + dx
		yf, _ = F(xf)
		count++
	}
	ans.iterations = count

	if yf == 0 {
		ans.root = xf
		ans.isRoot = true
		return
	}
	if yi*yf < 0 {
		if xi > xf {
			temp := xi
			xi = xf
			xf = temp
		}
		ans.a = xi
		ans.b = xf
		ans.isRange = true
		return
	}

	return
}

func Bisection(a, b, tolerance float64, maxIt uint) (ans *RootAnswer) {
	ans = &RootAnswer{}
	bisTable = make([]Reg1, int(maxIt))
	xi := a
	yi, _ := F(xi)
	if yi == 0 {
		ans.root = xi
		ans.isRoot = true
		return
	}

	xs := b
	ys, _ := F(xs)
	if ys == 0 {
		ans.root = xs
		ans.isRoot = true
		return
	}
	if yi*ys > 0 {
		ans.badIn = true
		return
	}
	xm := (xi + xs) / float64(2)
	ym, _ := F(xm)
	err := math.MaxFloat64
	count := uint(1)
	var temp float64

	bisTable = append(bisTable, Reg1{xi, xs, xm, ym, int(count)})
	for ym != 0 && err > tolerance && count < maxIt {
		if yi*ym < 0 {
			xs = xm
			ys = ym
		} else {
			xi = xm
			yi = ym
		}
		temp = xm
		xm = (xi + xs) / float64(2)
		ym, _ = F(xm)
		err = math.Abs(xm - temp)
		count++
		bisTable = append(bisTable, Reg1{xi, xs, xm, ym, int(count)})
	}
	ans.iterations = count
	if ym == 0 {
		ans.isRoot = true
		ans.root = xm
		return
	}
	if err <= tolerance {
		ans.isAlmostRoot = true
		ans.root = xm
		return
	}

	return
}

func falsePosition(a, b, tolerance float64, maxIt uint) (ans *RootAnswer) {
	ans = &RootAnswer{}
	falTable = make([]Reg1, int(maxIt))
	xi := a
	yi, _ := F(xi)
	if yi == 0 {
		ans.root = xi
		ans.isRoot = true
		return
	}

	xs := b
	ys, _ := F(xs)
	if ys == 0 {
		ans.root = xs
		ans.isRoot = true
		return
	}
	if yi*ys > 0 {
		ans.badIn = true
		return ans
	}
	xm := xi - (yi * (xs - xi) / (ys - yi))
	ym, _ := F(xm)
	err := math.MaxFloat64
	count := uint(1)
	var temp float64
	falTable = append(falTable, Reg1{xi, xs, xm, ym, int(count)})
	for ym != 0 && err > tolerance && count < maxIt {
		if yi*ym < 0 {
			xs = xm
			ys = ym
		} else {
			xi = xm
			yi = ym
		}
		temp = xm
		xm = xi - (yi * (xs - xi) / (ys - yi))
		ym, _ = F(xm)
		err = math.Abs(xm - temp)
		count++
		falTable = append(falTable, Reg1{xi, xs, xm, ym, int(count)})
	}
	ans.iterations = count
	if ym == 0 {
		ans.isRoot = true
		ans.root = xm
		return
	}
	if err <= tolerance {
		ans.isAlmostRoot = true
		ans.root = xm
		return
	}

	return
}

func fixedPoint(x0, tolerance float64, maxIt uint) (ans *RootAnswer) {
	ans = &RootAnswer{}
	fixTable = make([]Reg2, int(maxIt))
	xn := x0
	yn, _ := F(xn)
	err := math.MaxFloat64
	count := uint(0)
	var last float64
	fixTable = append(fixTable, Reg2{xn, yn, 0})

	for yn != 0 && err > tolerance && count < maxIt {
		last = xn
		xn, _ = G(xn)
		yn, _ = F(xn)
		err = math.Abs(xn - last)
		count++
		fixTable = append(fixTable, Reg2{xn, yn, int(count)})
	}
	ans.iterations = count
	if yn == 0 {
		ans.isRoot = true
		ans.root = xn
		return
	}
	if err <= tolerance {
		ans.isAlmostRoot = true
		ans.root = xn
		return
	}
	return
}

func newton(x0, tolerance float64, maxIt uint) (ans *RootAnswer) {
	ans = &RootAnswer{}
	newTable = make([]Reg4, int(maxIt))
	x := x0
	y, _ := F(x)
	dy, _ := DF(x)
	count := uint(0)
	err := math.MaxFloat64
	var x1 float64
	newTable = append(newTable, Reg4{x, y, dy, 0})
	for y != 0 && err > tolerance && dy != 0 && count < maxIt {
		x1 = x - y/dy
		y, _ = F(x1)
		dy, _ = DF(x1)
		err = math.Abs(x1 - x)
		x = x1
		count++
		newTable = append(newTable, Reg4{x, y, dy, int(count)})
	}
	ans.iterations = count
	if y == 0 {
		ans.isRoot = true
		ans.root = x
		return
	}
	if err <= tolerance {
		ans.isAlmostRoot = true
		ans.root = x
		return
	}
	if dy == 0 {
		ans.badIn = true
		return
	}
	return

}

func secant(a, b, tolerance float64, maxIt uint) (ans *RootAnswer) {
	ans = &RootAnswer{}
	secTable = make([]Reg2, int(maxIt))
	x0 := a
	y0, _ := F(x0)
	if y0 == 0 {
		ans.isRoot = true
		ans.root = x0
		return
	}
	x1 := b
	y1, _ := F(x1)
	count := uint(0)
	err := math.MaxFloat64
	var x2 float64
	for y1 != 0 && err > tolerance && y1 != y0 && count < maxIt {
		x2 = x1 - y1*(x1-x0)/(y1-y0)
		err = math.Abs(x2 - x1)
		x0 = x1
		y0 = y1
		x1 = x2
		y1, _ = F(x1)
		count++
		secTable = append(secTable, Reg2{x1, y1, int(count)})
	}
	ans.iterations = count
	if y1 == 0 {
		ans.isRoot = true
		ans.root = x1
		return
	}
	if err <= tolerance {
		ans.isAlmostRoot = true
		ans.root = x1
		return
	}
	if y1 == y0 {
		ans.badIn = true
	}
	return
}

func multipeRootAux(x0, tolerance float64, maxIt uint) (ans *RootAnswer) {
	ans = &RootAnswer{}
	mulTable = make([]Reg3, int(maxIt))
	x := x0
	y, _ := F(x)
	dy, _ := DF(x)
	dy2 := math.Pow(dy, 2)
	d2y, _ := D2F(x)
	count := uint(0)
	err := math.MaxFloat64
	var x1 float64
	mulTable = append(mulTable, Reg3{x, y, dy, d2y, 0})
	for y != 0 && err > tolerance && dy2 != y*d2y && count < maxIt {
		x1 = x - (y*dy)/(dy2-y*d2y)
		y, _ = F(x1)
		dy, _ = DF(x1)
		d2y, _ = D2F(x1)
		dy2 = math.Pow(dy, 2)
		err = math.Abs(x1 - x)
		x = x1
		count++
		mulTable = append(mulTable, Reg3{x, y, dy, d2y, int(count)})
	}
	ans.iterations = count
	if y == 0 {
		ans.isRoot = true
		ans.root = x
		return
	}
	if err <= tolerance {
		ans.isAlmostRoot = true
		ans.root = x
		return
	}
	if dy == 0 {
		ans.badIn = true
		return
	}
	return

}
