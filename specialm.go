package gof2

import (
	"fmt"
	"math/big"
)

// I is an immutable rectangular identity matrix defined such that I.At(r, c) is
// 1 if r == c and 0 otherwise.
type I struct {
	immutableM
	// r and c are the size of the matrix.
	r, c int
}

// Eye creates an identity matrix.
func Eye(rows, cols int) I {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("can't create %dx%d I: size must be positive", rows, cols))
	}
	return I{r: rows, c: cols}
}

// Size returns the size of the matrix.
func (i I) Size() (rows, cols int) {
	return i.r, i.c
}

// At returns a polynomial containing the element at the given one-based row
// and column. The returned value is a shared constant and must not be modified.
func (i I) At(r, c int) *big.Int {
	if r <= 0 || r > i.r || c <= 0 || c > i.c {
		panic(fmt.Sprintf("index (%d,%d) out of bounds (size %dx%d)", r, c, i.r, i.c))
	}
	return to01(r == c)
}

// Z is an immutable rectangular zero matrix.
type Z struct {
	immutableM
	// r and c are the size of the matrix.
	r, c int
}

// Zeros creates a zero matrix.
func Zeros(rows, cols int) Z {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("can't create %dx%d Z: size must be positive", rows, cols))
	}
	return Z{r: rows, c: cols}
}

// Size returns the size of the matrix.
func (z Z) Size() (rows, cols int) {
	return z.r, z.c
}

// At returns a polynomial containing the element at the given one-based row
// and column. The returned value is a shared constant and must not be modified.
func (z Z) At(r, c int) *big.Int {
	if r <= 0 || r > z.r || c <= 0 || c > z.c {
		panic(fmt.Sprintf("index (%d,%d) out of bounds (size %dx%d)", r, c, z.r, z.c))
	}
	return zeroP
}

// R is an immutable square matrix such that the multiplication X*R(n) gives Y
// where each row of X is rotated left n times to produce the corresponding row
// of Y. That is, R(1) is the matrix with the principle subdiagonal and the
// upper-right corner ones and all other elements zero, and R(n) is R(1)^n.
type R struct {
	immutableM
	// s is the size of the matrix.
	s int
	// n is the left-rotate amount.
	n int
}

// Rol creates a new left-rotation matrix of the given size and rotation amount.
func Rol(size, shift int) R {
	if size <= 0 {
		panic(fmt.Sprintf("can't create %dx%d R: size must be positive", size, size))
	}
	return R{s: size, n: shift % size}
}

// Size returns the size of the matrix.
func (rot R) Size() (rows, cols int) {
	return rot.s, rot.s
}

// At returns a polynomial containing the element at the given one-based row
// and column. The returned value is a shared constant and must not be modified.
func (rot R) At(r, c int) *big.Int {
	if r <= 0 || r > rot.s || c <= 0 || c > rot.s {
		panic(fmt.Sprintf("index (%d,%d) out of bounds (size %dx%d)", r, c, rot.s, rot.s))
	}
	r = (r-1+rot.n)%rot.s + 1
	if r < 0 {
		r += rot.s
	}
	return to01(r == c)
}

// S is an immutable square matrix such that the multiplication X*R(n), n > 0
// gives Y where each row of X is shifted left n times and filled in with zeros
// on the right to produce the corresponding row of Y; and for n < 0, the shift
// is done to the right instead. That is, S(1) is the matrix with the principle
// subdiagonal ones and all other elements zero, S(-1) is the same on the
// principle superdiagonal, S(n) with n > 0 is S(1)^n, and S(n) with n < 0 is
// S(-1)^(-n).
type S struct {
	immutableM
	// s is the size of the matrix.
	s int
	// n is the left-shift amount.
	n int
}

// Shl creates a new left-shift matrix of the given size and shift amount.
func Shl(size, shift int) S {
	if size <= 0 {
		panic(fmt.Sprintf("can't create %dx%d S: size must be positive", size, size))
	}
	return S{s: size, n: shift}
}

// Size returns the size of the matrix.
func (s S) Size() (rows, cols int) {
	return s.s, s.s // s :)
}

// At returns a polynomial containing the element at the given one-based row
// and column. The returned value is a shared constant and must not be modified.
func (s S) At(r, c int) *big.Int {
	if r <= 0 || r > s.s || c <= 0 || c > s.s {
		panic(fmt.Sprintf("index (%d,%d) out of bounds (size %dx%d)", r, c, s.s, s.s))
	}
	return to01(r+s.n == c)
}

type immutableM struct{}

// SetAt panics.
func (immutableM) SetAt(r, c int, p *big.Int) {
	panic("immutable matrix must be converted before modifying")
}

// AddAt panics.
func (immutableM) AddAt(r, c int, p *big.Int) *big.Int {
	panic("immutable matrix must be converted before modifying")
}

// MulAt panics.
func (immutableM) MulAt(r, c int, p *big.Int) *big.Int {
	panic("immutable matrix must be converted before modifying")
}
