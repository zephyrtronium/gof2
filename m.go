// Package gof2 implements some operations on large polynomials and matrices
// over GF(2).
//
// gof2 provides a variety of different matrix implementations with a single
// interface for working with them. Polynomials are standard big ints. The
// general focus is space efficiency rather than speed.
//
// The scope of gof2 is currently small, mostly serving as a way for me to
// calculate primitive polynomials of pseudo-random number generators.
package gof2

import (
	"fmt"
	"math/big"
	"math/bits"
)

// M represents a matrix.
type M interface {
	// Size returns the number of rows and columns in the matrix.
	Size() (rows, cols int)
	// At returns a polynomial containing the value at the given one-based row
	// and column. For matrices with binary elements (i.e., not polynomials),
	// the returned values are shared constants which must not be modified.
	// Out-of-bounds accesses panic.
	At(r, c int) *big.Int
	// SetAt sets the element at a one-based row and column index to the given
	// polynomial. If the matrix has only binary elements and the polynomial is
	// not 0 or 1, this will panic.
	SetAt(r, c int, p *big.Int)
	// AddAt adds a polynomial to the element at the given one-based row and
	// column. Panics if the polynomial is not 0 or 1 and the matrix has only
	// binary elements.
	AddAt(r, c int, p *big.Int) *big.Int
	// MulAt multiplies by a polynomial the element at the given one-based row
	// and column. Panics if the polynomial is not 0 or 1 and the matrix has
	// only binary elements.
	MulAt(r, c int, p *big.Int) *big.Int
}

// SM is a sparse matrix of size up to 65535x65535 of binary elements.
type SM struct {
	// r and c are the size of the matrix.
	r, c uint16
	// v is the map of zero-based coordinates to elements. The column occupies
	// the upper sixteen bits and the row the lower ones.
	v map[uint32]uint8
}

// NewSparse creates a zero matrix of the given size. Panics if either size is
// non-positive or greater than 65535.
func NewSparse(rows, cols int) *SM {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: size must be positive", rows, cols))
	}
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	return &SM{
		r: uint16(rows),
		c: uint16(cols),
		v: make(map[uint32]uint8),
	}
}

// Sparse converts any type of binary matrix to a new sparse matrix. Panics if
// the argument is a polynomial matrix with any element having degree higher
// than one, or if m is too large. Types SM, FM, I, Z, R, and S are
// special-cased. All other types are filled in O(mn) time.
func Sparse(m M) *SM {
	rows, cols := m.Size()
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	B := SM{uint16(rows), uint16(cols), make(map[uint32]uint8)}
	switch A := m.(type) {
	case *SM:
		for k, v := range A.v {
			if v != 0 {
				B.v[k] = v
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
					B.v[uint32(c<<16)|uint32(r)] = 1
				}
				w &= w - 1 // Mask off the low bit of w.
			}
			k += bits.UintSize
		}
	case I:
		for k := 0; k < rows; k++ {
			B.v[uint32(k)*0x00010001] = 1
		}
	case Z:
		// do nothing
	case R:
		for i := 0; i < rows; i++ {
			r := i + A.n%rows
			if r < 0 {
				r += rows
			}
			B.v[uint32(r<<16)|uint32(i)] = 1
		}
	case S:
		if A.n >= 0 {
			for i := 0; i < rows-A.n; i++ {
				B.v[uint32(i+A.n)<<16|uint32(i)] = 1
			}
		} else {
			for i := 0; i < rows+A.n; i++ {
				B.v[uint32(i)<<16|uint32(i-A.n)] = 1
			}
		}
	default:
		for r := 0; r < rows; r++ {
			for c := 0; c < cols; c++ {
				if check01(A.At(r+1, c+1)) != 0 {
					B.v[uint32(c<<16)|uint32(r)] = 1
				}
			}
		}
	}
	return &B
}

// Size returns the number of rows and columns in the matrix.
func (sm *SM) Size() (rows, cols int) {
	return int(sm.r), int(sm.c)
}

// At returns a polynomial containing the element at the given one-based row
// and column. The returned value is a shared constant and must not be
// modified.
func (sm *SM) At(r, c int) *big.Int {
	if sm.v[sm.index(r, c)] != 0 {
		return oneP
	}
	return zeroP
}

// SetAt sets the value at a one-based row and column index to the given
// polynomial. Panics if the index is out of bounds or if p is not 0 or 1.
func (sm *SM) SetAt(r, c int, p *big.Int) {
	sm.v[sm.index(r, c)] = check01(p)
}

// AddAt adds to the element at the given one-based row and column. Panics if
// the index is out of bounds or if p is not 0 or 1.
func (sm *SM) AddAt(r, c int, p *big.Int) *big.Int {
	k := sm.index(r, c)
	sm.v[k] ^= check01(p)
	return to01(sm.v[k] != 0)
}

// MulAt multiplies the element at the given one-based row and column. Panics
// if the index is out of bounds or if p is not 0 or 1.
func (sm *SM) MulAt(r, c int, p *big.Int) *big.Int {
	k := sm.index(r, c)
	sm.v[k] &= check01(p)
	return to01(sm.v[k] != 0)
}

// index panics if the given row or column indices are out of bounds and
// returns the corresponding sparse coordinate otherwise.
func (sm *SM) index(r, c int) uint32 {
	if r--; r < 0 || r >= int(sm.r) {
		panic(fmt.Sprintf("row index %d out of bounds (size %dx%d)", r+1, sm.r, sm.c))
	}
	if c--; c < 0 || c >= int(sm.c) {
		panic(fmt.Sprintf("column index %d out of bounds (size %dx%d)", c+1, sm.r, sm.c))
	}
	return uint32(uint16(c))<<16 | uint32(uint16(r))
}

// FM is a full matrix of size up to 65535x65535 of binary elements. In this
// implementation, FM is more space-efficient than SM by the time the SM
// reaches a load factor of 1/8.
type FM struct {
	r, c uint16
	v    *big.Int
}

// NewFull creates a zero matrix of the given size. Panics if either size is
// non-positive or greater than 65535.
func NewFull(rows, cols int) *FM {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: size must be positive", rows, cols))
	}
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	fm := FM{
		r: uint16(rows),
		c: uint16(cols),
		v: new(big.Int),
	}
	// Set the bit past the end of the matrix data so that we aren't
	// reallocating at random while doing work.
	fm.v.SetBit(fm.v, rows*cols, 1)
	return &fm
}

// Full converts any type of binary matrix to a new full matrix. Panics if
// the argument is a polynomial matrix with any element having degree higher
// than one, or if m is too large. Types SM, FM, I, R, and S are special-cased;
// all other types are filled in O(mn) time.
func Full(m M) *FM {
	rows, cols := m.Size()
	if rows > 65535 || cols > 65535 {
		panic(fmt.Sprintf("cannot make %dx%d matrix: maximum dimension is 65535", rows, cols))
	}
	B := FM{uint16(rows), uint16(cols), new(big.Int)}
	switch A := m.(type) {
	case *SM:
		B.v.SetBit(B.v, rows*cols, 1)
		for k, v := range A.v {
			if v != 0 {
				// Calculating indices manually lets us skip bounds checks.
				r, c := int(k&0xffff), int(k>>16)
				B.v.SetBit(B.v, c*rows+r, 1)
			}
		}
	case *FM:
		B.v.Set(A.v)
	case I:
		B.v.SetBit(B.v, rows*cols, 1)
		for k := 0; k < rows; k++ {
			B.v.SetBit(B.v, k*rows+k, 1)
		}
	default:
		B.v.SetBit(B.v, rows*cols, 1)
		for r := 0; r < rows; r++ {
			for c := 0; c < cols; c++ {
				p := check01(A.At(r+1, c+1))
				if p != 0 {
					B.v.SetBit(B.v, c*rows+r, 1)
				}
			}
		}
	}
	return &B
}

// Size returns the size of the matrix.
func (fm *FM) Size() (rows, cols int) {
	return int(fm.r), int(fm.c)
}

// At returns a polynomial containing the element at the given one-based row
// and column. The returned value is a shared constant and must not be
// modified.
func (fm *FM) At(r, c int) *big.Int {
	if fm.v.Bit(fm.index(r, c)) != 0 {
		return oneP
	}
	return zeroP
}

// SetAt sets the value at a one-based row and column index to the given
// polynomial. Panics if the index is out of bounds or if p is not 0 or 1.
func (fm *FM) SetAt(r, c int, p *big.Int) {
	fm.v.SetBit(fm.v, fm.index(r, c), uint(check01(p)))
}

// AddAt adds to the element at the given one-based row and column. Panics if
// the index is out of bounds or if p is not 0 or 1.
func (fm *FM) AddAt(r, c int, p *big.Int) *big.Int {
	k := fm.index(r, c)
	v := fm.v.Bit(k) ^ uint(check01(p))
	fm.v.SetBit(fm.v, k, v)
	return to01(v != 0)
}

// MulAt multiplies the element at the given one-based row and column. Panics
// if the index is out of bounds or if p is not 0 or 1.
func (fm *FM) MulAt(r, c int, p *big.Int) *big.Int {
	k := fm.index(r, c)
	v := fm.v.Bit(k) & uint(check01(p))
	fm.v.SetBit(fm.v, k, v)
	return to01(v != 0)
}

// index panics if the given row or column indices are out of bounds and
// returns the corresponding bit vector coordinate otherwise. The bit vector
// is column-major.
func (fm *FM) index(r, c int) int {
	if r--; r < 0 || r >= int(fm.r) {
		panic(fmt.Sprintf("row index %d out of bounds (size %dx%d)", r+1, fm.r, fm.c))
	}
	if c--; c < 0 || c >= int(fm.c) {
		panic(fmt.Sprintf("column index %d out of bounds (size %dx%d)", c+1, fm.r, fm.c))
	}
	return c*int(fm.r) + r
}

// check01 panics if the given polynomial is not zero or one and returns its
// value if it is.
func check01(p *big.Int) uint8 {
	if p.Sign() < 0 || p.BitLen() > 1 {
		s := "+"
		if p.Sign() < 0 {
			s = "-"
		}
		panic(fmt.Sprintf("cannot use polynomial %s%s in binary element matrix", s, p.Text(2)))
	}
	return uint8(p.Bit(0))
}

// to01 returns a polynomial representation of 1 if the argument is true or 0
// otherwise.
func to01(b bool) *big.Int {
	if b {
		return oneP
	}
	return zeroP
}

// Polynomial representations of F2's elements.
var zeroP, oneP = new(big.Int), big.NewInt(1)
