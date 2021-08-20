package ntrulpr

// counter to trigger recompilation: 1

// #cgo CFLAGS: -I../c/common -I../c/ntrulpr -DPREFIX=ntrulpr_
// #include "../c/common/namespace.h"
// #include "../c/common/fips202.c"
// #include "../c/common/randombytes.c"
// #include "../c/common/verify.c"
// #include "../c/ntrulpr/ntrulpr.c"
// #include "../c/ntrulpr/mkem.c"
// #include "../c/ntrulpr/ntt_1907713_1536.c"
import "C"

import (
	"github.com/PQShield/chained-cmpke/mkem"

	"crypto/subtle"
	"unsafe"
)

type Group struct {
	packed [C.MKEM_GROUPBYTES]byte
}

type PrivateKey struct {
	group  *Group
	packed [C.MKEM_SECRETKEYBYTES]byte
}

type PublicKey struct {
	group  *Group
	packed [C.MKEM_PUBLICKEYBYTES]byte
}

type scheme struct {
}

func (*scheme) Name() string                   { return "ntrulpr" }
func (*scheme) GroupSize() int                 { return C.MKEM_GROUPBYTES }
func (*scheme) PrivateKeySize() int            { return C.MKEM_SECRETKEYBYTES }
func (*scheme) PublicKeySize() int             { return C.MKEM_PUBLICKEYBYTES }
func (*scheme) IndependentCiphertextSize() int { return C.MKEM_CTIBYTES }
func (*scheme) DependentCiphertextSize() int   { return C.MKEM_CTDBYTES }
func (*scheme) SharedSecretSize() int          { return C.MKEM_SSBYTES }

func (*scheme) GenerateGroup() mkem.Group {
	var ret Group
	C.mkem_group((*C.uchar)(&ret.packed[0]))
	return &ret
}

func (s *scheme) GroupFromBytes(group []byte) (mkem.Group, error) {
	if len(group) != s.GroupSize() {
		return nil, mkem.ErrWrongSize
	}
	var ret Group
	copy(ret.packed[:], group)
	return &ret, nil
}

func (g *Group) Encapsulate(ss []byte, cti []byte, ctd [][]byte,
	pks []mkem.PublicKey) error {
	if len(ctd) != len(pks) {
		return mkem.ErrWrongSize
	}
	if len(ss) != _scheme.SharedSecretSize() {
		return mkem.ErrWrongSize
	}
	if len(ctd) == 0 {
		return mkem.ErrNoRecips
	}
	if len(cti) != _scheme.IndependentCiphertextSize() {
		return mkem.ErrWrongSize
	}
	group, ok := pks[0].Group().(*Group)
	if group != g {
		return mkem.ErrGroupMismatch
	}
	if !ok {
		return mkem.ErrSchemeMismatch
	}

	pkPtrs := make([]unsafe.Pointer, len(pks))
	ctdPtrs := make([]unsafe.Pointer, len(pks))

	for i := 0; i < len(pks); i++ {
		if len(ctd[i]) != _scheme.DependentCiphertextSize() {
			return mkem.ErrWrongSize
		}
		packedPk, ok := pks[i].(*PublicKey)
		if !ok {
			return mkem.ErrSchemeMismatch
		}
		if pks[i].Group() != g {
			return mkem.ErrGroupMismatch
		}
		pkPtrs[i] = C.CBytes(packedPk.packed[:])
		ctdPtrs[i] = C.CBytes(ctd[i])
	}
	C.mkem_enc(
		(**C.uchar)((unsafe.Pointer)(&ctdPtrs[0])),
		(*C.uchar)(&cti[0]),
		(*C.uchar)(&ss[0]),
		(*C.uchar)(&g.packed[0]),
		(**C.uchar)((unsafe.Pointer)(&pkPtrs[0])),
		(C.uint)(len(pks)),
	)

	for i := 0; i < len(pks); i++ {
		copy(ctd[i], C.GoBytes(ctdPtrs[i], C.int(len(ctd[i]))))
		C.free(pkPtrs[i])
		C.free(ctdPtrs[i])
	}
	return nil
}

func (g *Group) GenerateKeypair() (mkem.PrivateKey, mkem.PublicKey) {
	pk := PublicKey{group: g}
	sk := PrivateKey{group: g}
	C.mkem_keypair(
		(*C.uchar)(&pk.packed[0]),
		(*C.uchar)(&sk.packed[0]),
		(*C.uchar)(&g.packed[0]),
	)
	return &sk, &pk
}

func (g *Group) Decapsulate(ss []byte, cti []byte, ctd []byte,
	sk mkem.PrivateKey, pk mkem.PublicKey) error {
	if len(ss) != _scheme.SharedSecretSize() ||
		len(cti) != _scheme.IndependentCiphertextSize() ||
		len(ctd) != _scheme.DependentCiphertextSize() {
		return mkem.ErrWrongSize
	}

	pk2, ok := pk.(*PublicKey)
	if !ok {
		return mkem.ErrSchemeMismatch
	}
	sk2, ok := sk.(*PrivateKey)
	if !ok {
		return mkem.ErrSchemeMismatch
	}
	if pk2.group != g || sk2.group != g {
		return mkem.ErrGroupMismatch
	}

	C.mkem_dec(
		(*C.uchar)(&ss[0]),
		(*C.uchar)(&cti[0]),
		(*C.uchar)(&ctd[0]),
		(*C.uchar)(&g.packed[0]),
		(*C.uchar)(&pk2.packed[0]),
		(*C.uchar)(&sk2.packed[0]),
	)
	return nil
}

func (*Group) Scheme() mkem.Scheme {
	return _scheme
}

func (pk *PublicKey) Group() mkem.Group {
	return pk.group
}
func (sk *PrivateKey) Group() mkem.Group {
	return sk.group
}

func (pk *PublicKey) Pack(buf []byte) {
	if len(buf) != _scheme.PublicKeySize() {
		panic(mkem.ErrWrongSize)
	}

	copy(buf, pk.packed[:])
}

func (sk *PrivateKey) Pack(buf []byte) {
	if len(buf) != _scheme.PrivateKeySize() {
		panic(mkem.ErrWrongSize)
	}

	copy(buf, sk.packed[:])
}

func (g *Group) Pack(buf []byte) {
	if len(buf) != _scheme.GroupSize() {
		panic(mkem.ErrWrongSize)
	}

	copy(buf, g.packed[:])
}

func (g *Group) UnpackPrivateKey(buf []byte) (mkem.PrivateKey, error) {
	if len(buf) != _scheme.PrivateKeySize() {
		return nil, mkem.ErrWrongSize
	}

	ret := PrivateKey{group: g}
	copy(ret.packed[:], buf)
	return &ret, nil
}

func (g *Group) UnpackPublicKey(buf []byte) (mkem.PublicKey, error) {
	if len(buf) != _scheme.PublicKeySize() {
		return nil, mkem.ErrWrongSize
	}

	ret := PublicKey{group: g}
	copy(ret.packed[:], buf)
	return &ret, nil
}

func (s *scheme) UnpackGroup(buf []byte) (mkem.Group, error) {
	if len(buf) != _scheme.GroupSize() {
		return nil, mkem.ErrWrongSize
	}
	var ret Group
	copy(ret.packed[:], buf)
	return &ret, nil
}

var _scheme mkem.Scheme = &scheme{}

func Scheme() mkem.Scheme {
	return _scheme
}

func (g *Group) Equals(other mkem.Group) bool {
	g2, ok := other.(*Group)
	if !ok {
		return false
	}
	return subtle.ConstantTimeCompare(g.packed[:], g2.packed[:]) == 1
}

func (sk *PrivateKey) Equals(other mkem.PrivateKey) bool {
	sk2, ok := other.(*PrivateKey)
	if !ok {
		return false
	}
	return subtle.ConstantTimeCompare(sk.packed[:], sk2.packed[:]) == 1
}

func (pk *PublicKey) Equals(other mkem.PublicKey) bool {
	pk2, ok := other.(*PublicKey)
	if !ok {
		return false
	}
	return subtle.ConstantTimeCompare(pk.packed[:], pk2.packed[:]) == 1
}
