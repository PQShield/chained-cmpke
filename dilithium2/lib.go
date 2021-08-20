package dilithium2

// counter to trigger recompilation: 1

// #include "api.h"
import "C"

import (
	"bytes"
	"crypto"
	"io"
)

const (
	PublicKeySize  = C.PQCLEAN_DILITHIUM2_CLEAN_CRYPTO_PUBLICKEYBYTES
	PrivateKeySize = C.PQCLEAN_DILITHIUM2_CLEAN_CRYPTO_SECRETKEYBYTES
	SignatureSize  = C.PQCLEAN_DILITHIUM2_CLEAN_CRYPTO_BYTES
)

type PublicKey [PublicKeySize]byte
type PrivateKey [PrivateKeySize]byte

func GenerateKey(io.Reader) (*PublicKey, *PrivateKey, error) {
	var pk PublicKey
	var sk PrivateKey

	C.PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_keypair(
		(*C.uchar)(&pk[0]),
		(*C.uchar)(&sk[0]),
	)

	return &pk, &sk, nil
}

func Verify(pk *PublicKey, msg, sig []byte) bool {
	return C.PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_verify(
		(*C.uchar)(&sig[0]),
		C.size_t(len(sig)),
		(*C.uchar)(&msg[0]),
		C.size_t(len(msg)),
		(*C.uchar)(&pk[0]),
	) == 0
}

func (sk *PrivateKey) Sign(r io.Reader, msg []byte, so crypto.SignerOpts) ([]byte, error) {
	var ret [SignatureSize]byte
	var slen C.size_t = C.size_t(len(ret))
	C.PQCLEAN_DILITHIUM2_CLEAN_crypto_sign_signature(
		(*C.uchar)(&ret[0]),
		&slen,
		(*C.uchar)(&msg[0]),
		C.size_t(len(msg)),
		(*C.uchar)(&sk[0]),
	)
	return ret[:], nil
}

func (pk *PublicKey) Equal(other *PublicKey) bool {
	return bytes.Equal(pk[:], other[:])
}

func (pk *PublicKey) Pack(out *[PublicKeySize]byte) {
	copy(out[:], pk[:])
}
