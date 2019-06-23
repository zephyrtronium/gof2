package gof2

import (
	"fmt"
)

// FMul multiplies two matrices in GF(2). If either argument is sparse, the
// result is SM. If either argument is Z, the result is Z. Panics if the inner
// dimensions of the matrices are not equal or if any element is not 0 or 1.
func FMul(A, B M) M {
	ar, ac := A.Size()
	br, bc := B.Size()
	if ac != br {
		panic(fmt.Sprintf("inner dimension mismatch: %dx%d * %dx%d", ar, ac, br, bc))
	}
	switch x := A.(type) {
	case Z:
		return Zeros(ar, bc)
	case *SM:
		return fMulSX(x, B)
	case *PSM:
		return fMulPSX(x, B)
	case I:
		switch B.(type) {
		case *SM:
			return Sparse(B)
		case *FM:
			return Full(B)
		case I:
			return Eye(ar, bc)
		}
	}
	switch x := B.(type) {
	case Z:
		return Zeros(ar, bc)
	case *SM:
		return fMulXS(A, x)
	case *PSM:
		return fMulXPS(A, x)
	case I:
		if _, ok := A.(*FM); ok {
			return Full(A)
		}
	}
	return fMulFull(A, B)
}

// fMulFull multiplies two matrices into a new FM.
func fMulFull(A, B M) *FM {
	ar, ac := A.Size()
	_, bc := B.Size()
	C := NewFull(ar, bc)
	for c := 1; c <= bc; c++ {
		for r := 1; r <= ar; r++ {
			var d uint8
			for i := 1; i <= ac; i++ {
				d ^= check01(A.At(r, i)) & check01(B.At(i, c))
			}
			C.SetAt(r, c, to01(d != 0))
		}
	}
	return C
}

// fMulSX multiplies a sparse matrix by another matrix into a new SM.
func fMulSX(A *SM, B M) *SM {
	ar, ac := A.Size()
	_, bc := B.Size()
	C := NewSparse(ar, bc)
	switch X := B.(type) {
	case *SM:
		for j, a := range A.v {
			if a == 0 {
				continue
			}
			for k, b := range X.v {
				if j>>16 == k&0xffff {
					// The column of the A element equals the row of the B
					// element. Their product is a term in the element of C at
					// the row of A and the column of B.
					C.v[k&0xffff0000|j&0x0000ffff] ^= a & b
				}
			}
		}
	case *PSM:
		for j, a := range A.v {
			if a == 0 {
				continue
			}
			for k, b := range X.v {
				if j>>16 == k&0xffff {
					C.v[k&0xffff0000|j&0x0000ffff] ^= a & check01(b)
				}
			}
		}
	case I:
		for j, a := range A.v {
			if a != 0 && j>>16 < uint32(bc) && j&0xffff < uint32(ar) {
				C.v[j] = 1
			}
		}
	case Z:
		// do nothing
	case R:
		for j, a := range A.v {
			if a != 0 {
				r, c := j&0xffff, j>>16
				cc := (int(c) - X.n) % bc
				if cc < 0 {
					cc += bc
				}
				C.v[uint32(cc)<<16|r] = 1
			}
		}
	case S:
		if X.n >= 0 {
			for j, a := range A.v {
				r, c := j&0xffff, j>>16
				cc := int(c) - X.n
				if a != 0 && cc >= 0 {
					C.v[uint32(cc)<<16|r] = 1
				}
			}
		} else {
			for j, a := range A.v {
				r, c := j&0xffff, j>>16
				cc := int(c) - X.n
				if a != 0 && cc < bc {
					C.v[uint32(cc)<<16|r] = 1
				}
			}
		}
	default:
		for j, a := range A.v {
			if a == 0 {
				continue
			}
			r, c := j&0xffff, int(j>>16)
			// This element multiplies with each element of the cth row of B
			// into the rth row and respective column of C.
			for i := 0; i < ac; i++ {
				b := check01(B.At(c+1, i+1))
				if b != 0 {
					C.v[uint32(i)<<16|r] ^= 1
				}
			}
		}
	}
	return C
}

// fMulPSX multiplies a sparse polynomial matrix by another matrix into a new
// SM.
func fMulPSX(A *PSM, B M) *SM {
	ar, ac := A.Size()
	_, bc := B.Size()
	C := NewSparse(ar, bc)
	switch X := B.(type) {
	case *SM:
		for j, a := range A.v {
			if a.Sign() == 0 {
				continue
			}
			for k, b := range X.v {
				if j>>16 == k&0xffff {
					// The column of the A element equals the row of the B
					// element. Their product is a term in the element of C at
					// the row of A and the column of B.
					C.v[k&0xffff0000|j&0x0000ffff] ^= check01(a) & b
				}
			}
		}
	case *PSM:
		for j, a := range A.v {
			if a.Sign() == 0 {
				continue
			}
			for k, b := range X.v {
				if j>>16 == k&0xffff {
					C.v[k&0xffff0000|j&0x0000ffff] ^= check01(a) & check01(b)
				}
			}
		}
	case I:
		for j, a := range A.v {
			if check01(a) != 0 && j>>16 < uint32(bc) && j&0xffff < uint32(ar) {
				C.v[j] = 1
			}
		}
	case Z:
		// do nothing
	case R:
		for j, a := range A.v {
			if check01(a) != 0 {
				r, c := j&0xffff, j>>16
				cc := (int(c) - X.n) % bc
				if cc < 0 {
					cc += bc
				}
				C.v[uint32(cc)<<16|r] = 1
			}
		}
	case S:
		if X.n >= 0 {
			for j, a := range A.v {
				r, c := j&0xffff, j>>16
				cc := int(c) - X.n
				if check01(a) != 0 && cc >= 0 {
					C.v[uint32(cc)<<16|r] = 1
				}
			}
		} else {
			for j, a := range A.v {
				r, c := j&0xffff, j>>16
				cc := int(c) - X.n
				if check01(a) != 0 && cc < bc {
					C.v[uint32(cc)<<16|r] = 1
				}
			}
		}
	default:
		for j, a := range A.v {
			if check01(a) == 0 {
				continue
			}
			r, c := j&0xffff, int(j>>16)
			// This element multiplies with each element of the cth row of B
			// into the rth row and respective column of C.
			for i := 0; i < ac; i++ {
				b := check01(B.At(c+1, i+1))
				if b != 0 {
					C.v[uint32(i)<<16|r] ^= 1
				}
			}
		}
	}
	return C
}

// fMulXS multiplies a matrix by an SM into a new SM.
func fMulXS(A M, B *SM) *SM {
	ar, ac := A.Size()
	_, bc := B.Size()
	C := NewSparse(ar, bc)
	switch X := A.(type) {
	case R:
		for j, a := range B.v {
			if a != 0 {
				r, c := j&0xffff, j>>16
				rr := (int(r) - X.n) % ar
				if rr < 0 {
					rr += ar
				}
				C.v[c<<16|uint32(rr)] = 1
			}
		}
	case S:
		if X.n >= 0 {
			for j, a := range B.v {
				r, c := j&0xffff, j>>16
				rr := int(r) + X.n
				if a != 0 && rr < ar {
					C.v[c<<16|uint32(rr)] = 1
				}
			}
		} else {
			for j, a := range B.v {
				r, c := j&0xfff, j>>16
				rr := int(r) + X.n
				if a != 0 && rr >= 0 {
					C.v[c<<16|uint32(rr)] = 1
				}
			}
		}
	default:
		for j, a := range B.v {
			if a == 0 {
				continue
			}
			r, c := int(j&0xffff), j>>16
			// This element multiplies with each element of the rth column of B
			// into the respective row and cth column of C.
			for i := 0; i < ac; i++ {
				b := check01(A.At(i+1, r+1))
				if b != 0 {
					C.v[c<<16|uint32(i)] ^= 1
				}
			}
		}
	}
	return C
}

// fMulXPS multiplies a matrix by a PSM into a new SM.
func fMulXPS(A M, B *PSM) *SM {
	ar, ac := A.Size()
	_, bc := B.Size()
	C := NewSparse(ar, bc)
	switch X := A.(type) {
	case R:
		for j, a := range B.v {
			if check01(a) != 0 {
				r, c := j&0xffff, j>>16
				rr := (int(r) - X.n) % ar
				if rr < 0 {
					rr += ar
				}
				C.v[c<<16|uint32(rr)] = 1
			}
		}
	case S:
		if X.n >= 0 {
			for j, a := range B.v {
				r, c := j&0xffff, j>>16
				rr := int(r) + X.n
				if check01(a) != 0 && rr < ar {
					C.v[c<<16|uint32(rr)] = 1
				}
			}
		} else {
			for j, a := range B.v {
				r, c := j&0xfff, j>>16
				rr := int(r) + X.n
				if check01(a) != 0 && rr >= 0 {
					C.v[c<<16|uint32(rr)] = 1
				}
			}
		}
	default:
		for j, a := range B.v {
			if check01(a) == 0 {
				continue
			}
			r, c := int(j&0xffff), j>>16
			// This element multiplies with each element of the rth column of B
			// into the respective row and cth column of C.
			for i := 0; i < ac; i++ {
				b := check01(A.At(i+1, r+1))
				if b != 0 {
					C.v[c<<16|uint32(i)] ^= 1
				}
			}
		}
	}
	return C
}
