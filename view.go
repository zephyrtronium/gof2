package gof2

import (
	"fmt"
	"math/big"
)

// V provides a submatrix view. It proxies methods to its matrix offset by the
// view area.
type V struct {
	m M
	// r, c are the view's row and column offsets into the parent matrix.
	r, c int
	// rr, cc are the view's size.
	rr, cc int
}

// View creates a rows x cols view of A with its top-left at (r, c). This does
// not panic if the view exceeds the bounds of A, but out-of-bounds accesses
// may panic.
func View(A M, r, c, rows, cols int) V {
	return V{A, r, c, rows, cols}
}

// Size returns the size of the view.
func (v V) Size() (rows, cols int) {
	return v.rr, v.cc
}

// At proxies to the viewed matrix's At method. This method will panic if the
// index is outside the view's bounds, and the matrix's method may panic if the
// adjusted index is outside its bounds.
func (v V) At(r, c int) *big.Int {
	return v.m.At(v.index(r, c))
}

// SetAt proxies to the viewed matrix's SetAt method. This method will panic if
// the index is outside the view's bounds, and the matrix's method may panic if
// the adjusted index is outside its bounds.
func (v V) SetAt(r, c int, p *big.Int) {
	r, c = v.index(r, c)
	v.m.SetAt(r, c, p)
}

// AddAt proxies to the viewed matrix's AddAt method. This method will panic if
// the index is outside the view's bounds, and the matrix's method may panic if
// the adjusted index is outside its bounds.
func (v V) AddAt(r, c int, p *big.Int) *big.Int {
	r, c = v.index(r, c)
	return v.m.AddAt(r, c, p)
}

// MulAt proxies to the viewed matrix's MulAt method. This method will panic if
// the index is outside the view's bounds, and the matrix's method may panic if
// the adjusted index is outside its bounds.
func (v V) MulAt(r, c int, p *big.Int) *big.Int {
	r, c = v.index(r, c)
	return v.m.MulAt(r, c, p)
}

// index panics if the given 1-based index is outside the view's bounds and
// returns the corresponding index into the viewed matrix.
func (v V) index(r, c int) (int, int) {
	if r <= 0 || r > v.rr {
		panic(fmt.Sprintf("row index %d out of bounds (size %dx%d)", r, v.rr, v.cc))
	}
	if c <= 0 || c > v.cc {
		panic(fmt.Sprintf("column index %d out of bounds (size %dx%d)", c, v.rr, v.cc))
	}
	return r - v.r + 1, c - v.c + 1
}
