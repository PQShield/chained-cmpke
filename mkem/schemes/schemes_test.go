package schemes_test

import (
	"bytes"
	"testing"

	"github.com/PQShield/chained-cmpke/mkem"
	"github.com/PQShield/chained-cmpke/mkem/schemes"
)

func TestApi(t *testing.T) {
	allSchemes := schemes.All()

	for _, scheme := range allSchemes {
		scheme := scheme

		t.Run(scheme.Name(), func(t *testing.T) {
			if scheme == nil {
				t.Fatal()
			}

			g := scheme.GenerateGroup()
			if g == nil {
				t.Fatal()
			}

			sk, pk := g.GenerateKeypair()
			if sk == nil || pk == nil {
				t.Fatal()
			}

			ppk := make([]byte, scheme.PublicKeySize())
			psk := make([]byte, scheme.PrivateKeySize())
			pg := make([]byte, scheme.GroupSize())

			sk.Pack(psk)
			pk.Pack(ppk)
			g.Pack(pg)

			g2, err := scheme.UnpackGroup(pg)
			if err != nil || !g2.Equals(g) {
				t.Fatal(err)
			}
			pk2, err := g.UnpackPublicKey(ppk)
			if err != nil || !pk2.Equals(pk) {
				t.Fatal(err)
			}
			sk2, err := g.UnpackPrivateKey(psk)
			if err != nil || !sk2.Equals(sk) {
				t.Fatal(err)
			}

			sk2, pk2 = g.GenerateKeypair()

			ss := make([]byte, scheme.SharedSecretSize())
			ss2 := make([]byte, scheme.SharedSecretSize())
			cti := make([]byte, scheme.IndependentCiphertextSize())
			ctd := [][]byte{
				make([]byte, scheme.DependentCiphertextSize()),
				make([]byte, scheme.DependentCiphertextSize()),
			}
			pks := []mkem.PublicKey{pk, pk2}
			sks := []mkem.PrivateKey{sk, sk2}

			err = g.Encapsulate(
				ss,
				cti,
				ctd,
				pks,
			)
			if err != nil {
				t.Fatal(err)
			}

			for i := 0; i < len(pks); i++ {
				err = g.Decapsulate(ss2, cti, ctd[i], sks[i], pks[i])
				if err != nil {
					t.Fatal(err)
				}

				if !bytes.Equal(ss, ss2) {
					t.Fatal()
				}
			}
		})
	}
}
