package mkem

import (
	"errors"
)

type PrivateKey interface {
	Group() Group
	Pack(buf []byte)
	Equals(sk PrivateKey) bool
}

type PublicKey interface {
	Group() Group
	Pack(buf []byte)
	Equals(pk PublicKey) bool
}

type Group interface {
	Scheme() Scheme

	GenerateKeypair() (PrivateKey, PublicKey)

	Encapsulate(ss []byte, cti []byte, ctd [][]byte, pks []PublicKey) error
	Decapsulate(ss []byte, cti []byte, ctd []byte, sk PrivateKey,
		pk PublicKey) error

	UnpackPublicKey(buf []byte) (PublicKey, error)
	UnpackPrivateKey(buf []byte) (PrivateKey, error)
	Pack(buf []byte)
	Equals(g Group) bool
}

type Scheme interface {
	Name() string

	GroupSize() int
	PrivateKeySize() int
	PublicKeySize() int
	IndependentCiphertextSize() int
	DependentCiphertextSize() int
	SharedSecretSize() int

	GenerateGroup() Group
	UnpackGroup(buf []byte) (Group, error)
}

var (
	ErrWrongSize      = errors.New("an argument has the wrong size")
	ErrSchemeMismatch = errors.New("mixing objects of different schemes")
	ErrGroupMismatch  = errors.New("mixing objects of different groups")
	ErrNoRecips       = errors.New("empty list of recipients")
)
