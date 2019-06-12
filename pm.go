package gof2

import (
	"fmt"
	"math/big"
	"math/bits"
)

// PSM is a sparse polynomial matrix of size up to 65535x65535.
type PSM struct {
	// r and c are the size of the matrix.
	r, c uint16
	// v is the map of coordinates to polynomials. The column occupies the
	// upper 16 bits and the row the lower ones.
	v map[uint32]*big.Int
}

// NewPSparse creates a zero matrix of the given size. Panics if either size is
// non-positive or greater than 65535.
func NewPSparse(rows, cols int) *PSM {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: size must be positive", rows, cols))
	}
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	return &PSM{
		r: uint16(rows),
		c: uint16(cols),
		v: make(map[uint32]*big.Int),
	}
}

// PSparse converts any type of matrix to a sparse polynomial matrix. Panics if
// the argument is too large. Types SM, FM, PSM, PFM, I, Z, R, and S are
// special-cased. All other types are filled in O(mn) time.
func PSparse(m M) *PSM {
	rows, cols := m.Size()
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	B := PSM{uint16(rows), uint16(cols), make(map[uint32]*big.Int)}
	switch A := m.(type) {
	case *SM:
		for k, v := range A.v {
			if v != 0 {
				B.v[k] = big.NewInt(1)
			}
		}
	case *FM:
		k := 0
		for _, w := range A.v.Bits() {
			for w != 0 {
				b := k + bits.TrailingZeros(uint(w))
				r, c := b%int(B.r), b/int(B.r)
				// The bit vector in the full matrix has a bit set past the end
				// of the matrix data, so we need to make sure we don't include
				// that.
				if c < int(B.c) {
					B.v[uint32(c<<16)|uint32(r)] = big.NewInt(1)
				}
				w &= w - 1 // Mask off the low bit of w.
			}
			k += bits.UintSize
		}
	case *PSM:
		for k, v := range A.v {
			if v.Sign() != 0 {
				B.v[k] = new(big.Int).Set(v)
			}
		}
	case I:
		for r := 0; r < rows; r++ {
			B.v[uint32(r)*0x00010001] = big.NewInt(1)
		}
	case Z:
		// do nothing
	case R:
		for i := 0; i < rows; i++ {
			r := i + A.n%rows
			if r < 0 {
				r += rows
			}
			B.v[uint32(r<<16)|uint32(i)] = big.NewInt(1)
		}
	case S:
		if A.n >= 0 {
			for i := 0; i < rows-A.n; i++ {
				B.v[uint32(i+A.n)<<16|uint32(i)] = big.NewInt(1)
			}
		} else {
			for i := 0; i < rows+A.n; i++ {
				B.v[uint32(i)<<16|uint32(i-A.n)] = big.NewInt(1)
			}
		}
	}
	return &B
}

// Size returns the size of the matrix.
func (A *PSM) Size() (rows, cols int) {
	return int(A.r), int(A.c)
}

// At returns the polynomial at the given one-based index. The returned element
// is a reference if and only if it is nonzero.
func (A *PSM) At(r, c int) *big.Int {
	k := A.index(r, c)
	if p := A.v[k]; p == nil || p.Sign() == 0 {
		if p.Sign() == 0 {
			delete(A.v, k)
		}
		return p
	}
	return new(big.Int)
}

// SetAt sets the polynomial at the given one-based index. The polynomial is
// not copied.
func (A *PSM) SetAt(r, c int, p *big.Int) {
	k := A.index(r, c)
	if p.Sign() == 0 {
		delete(A.v, k)
	} else {
		A.v[k] = p
	}
}

// AddAt adds a polynomial to that in the given index. The returned value is
// always a reference, even if it is zero.
func (A *PSM) AddAt(r, c int, p *big.Int) *big.Int {
	k := A.index(r, c)
	q, ok := A.v[k]
	if !ok {
		q = new(big.Int)
		A.v[k] = q
	}
	return q.Xor(q, p)
}

// MulAt multiplies (i.e. convolves coefficients of) the polynomial in a given
// index by another. The returned value is always a reference, even if it is
// zero.
func (A *PSM) MulAt(r, c int, p *big.Int) *big.Int {
	k := A.index(r, c)
	q, ok := A.v[k]
	if !ok {
		q = new(big.Int)
		A.v[k] = q
		return q
	}
	return q.Mul(q, p)
}

// index panics if the given row or column indices are out of bounds and
// returns the corresponding sparse coordinate otherwise.
func (A *PSM) index(r, c int) uint32 {
	if r--; r < 0 || r >= int(A.r) {
		panic(fmt.Sprintf("row index %d out of bounds (size %dx%d)", r+1, A.r, A.c))
	}
	if c--; c < 0 || c >= int(A.c) {
		panic(fmt.Sprintf("column index %d out of bounds (size %dx%d)", c+1, A.r, A.c))
	}
	return uint32(uint16(c))<<16 | uint32(uint16(r))
}

// PFM is a full polynomial matrix of size up to 65535x65535.
type PFM struct {
	r, c uint16
	v    []*big.Int
}

// NewPFull creates a zero matrix of the given size. This allocates all
// elements. Panics if either size is non-positive or greater than 65535.
func NewPFull(rows, cols int) *PFM {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: size must be positive", rows, cols))
	}
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	v := make([]*big.Int, rows*cols)
	for i := range v {
		v[i] = new(big.Int)
	}
	return &PFM{
		r: uint16(rows),
		c: uint16(cols),
		v: v,
	}
}

// PFull converts any type of matrix to a full polynomial matrix. Panics if the
// argument is too large. There are no special cases; converting any matrix
// results in m*n calls to m.At().
func PFull(m M) *PFM {
	rows, cols := m.Size()
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	B := PFM{uint16(rows), uint16(cols), make([]*big.Int, rows*cols)}
	for c := 0; c < cols; c++ {
		for r := 0; r < rows; r++ {
			B.v[c*rows+r] = new(big.Int).Set(m.At(r, c))
		}
	}
	return &B
}

// Size returns the size of the matrix.
func (A *PFM) Size() (rows, cols int) {
	return int(A.r), int(A.c)
}

// At returns the element at the given one-based index. The returned element is
// always a reference.
func (A *PFM) At(r, c int) *big.Int {
	return A.v[A.index(r, c)]
}

// SetAt sets an element. The polynomial is not copied.
func (A *PFM) SetAt(r, c int, p *big.Int) {
	A.v[A.index(r, c)] = p
}

// AddAt adds a polynomial to that in the given index. The returned value is a
// a reference.
func (A *PFM) AddAt(r, c int, p *big.Int) *big.Int {
	k := A.index(r, c)
	return A.v[k].Xor(A.v[k], p)
}

// MulAt multiplies (i.e. convolves coefficients of) the polynomial in a given
// index by another. The returned value is a reference.
func (A *PFM) MulAt(r, c int, p *big.Int) *big.Int {
	k := A.index(r, c)
	return A.v[k].Mul(A.v[k], p)
}

// index panics if the given row or column indices are out of bounds and
// returns the corresponding bit vector coordinate otherwise. The vector is
// column-major.
func (A *PFM) index(r, c int) int {
	if r--; r < 0 || r >= int(A.r) {
		panic(fmt.Sprintf("row index %d out of bounds (size %dx%d)", r+1, A.r, A.c))
	}
	if c--; c < 0 || c >= int(A.c) {
		panic(fmt.Sprintf("column index %d out of bounds (size %dx%d)", c+1, A.r, A.c))
	}
	return c*int(A.r) + r
}
