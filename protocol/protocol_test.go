package protocol_test

import (
	"fmt"
	"strconv"
	"testing"

	"github.com/PQShield/chained-cmpke/mkem"
	"github.com/PQShield/chained-cmpke/mkem/schemes"
	"github.com/PQShield/chained-cmpke/protocol"
)

func BenchmarkKeypair(b *testing.B) {
	for _, scheme := range schemes.All() {
		b.Run(scheme.Name(), func(b *testing.B) {
			group := scheme.GenerateGroup()
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				group.GenerateKeypair()
			}
		})
	}
}

func BenchmarkCmDecrypt(b *testing.B) {
	for _, scheme := range schemes.All() {
		group := scheme.GenerateGroup()
		sk, pk := group.GenerateKeypair()
		M := make([]byte, 16)
		cti, ctd := protocol.CmEncrypt(M, []mkem.PublicKey{pk})

		b.Run(scheme.Name(), func(b *testing.B) {
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				protocol.CmDecrypt(cti, ctd[0], sk, pk)
			}
		})
	}
}

func BenchmarkCmEncrypt(b *testing.B) {
	for _, scheme := range schemes.All() {
		group := scheme.GenerateGroup()
		maxLogN := 10
		maxN := 1 << maxLogN
		pks := make([]mkem.PublicKey, maxN)
		for i := 0; i < maxN; i++ {
			_, pks[i] = group.GenerateKeypair()
		}
		M := make([]byte, 16)

		for logN := 0; logN <= maxLogN; logN++ {
			N := 1 << logN
			b.Run(fmt.Sprintf("%s-%d", scheme.Name(), N), func(b *testing.B) {
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					protocol.CmEncrypt(M, pks[:N])
				}
			})
		}
	}
}

func BenchmarkCommit(b *testing.B) {
	for _, scheme := range schemes.All() {
		for logN := 0; logN < 11; logN++ {
			N := 1 << logN
			b.Run(fmt.Sprintf("%s-%d", scheme.Name(), N), func(b *testing.B) {
				ctx := protocol.CreateContext(&protocol.ContextOpts{
					MKemScheme: scheme.Name(),
				})
				ctx.CreateMember("0")
				ctx.CreateGroup("0")
				ctx.CommitProcessAndJoin("0")

				for i := 1; i < N; i++ {
					ctx.CreateMember(strconv.Itoa(i))
					ctx.ProposeAdd("0", strconv.Itoa(i))
				}
				ctx.CommitProcessAndJoin("0")

				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					c0, cs, _, _ := ctx.GroupState("0").Commit(
						ctx, []*protocol.FramedProposal{})
					ctx.GroupState("0").Process(
						ctx, c0, cs["0"], []*protocol.FramedProposal{})
				}
			})
		}
	}
}

func BenchmarkProcess(b *testing.B) {
	for _, scheme := range schemes.All() {
		if scheme.Name() != "bilbo" {
			continue
		}
		for logN := 8; logN < 9; logN++ {
			N := 1 << logN
			b.Run(fmt.Sprintf("%s-%d", scheme.Name(), N), func(b *testing.B) {
				ctx := protocol.CreateContext(&protocol.ContextOpts{
					MKemScheme: scheme.Name(),
				})
				ctx.CreateMember("0")
				ctx.CreateGroup("0")
				ctx.CommitProcessAndJoin("0")

				for i := 1; i < N; i++ {
					ctx.CreateMember(strconv.Itoa(i))
					ctx.ProposeAdd("0", strconv.Itoa(i))
				}
				ctx.CommitProcessAndJoin("0")
				c0, cs, _, _ := ctx.GroupState("0").Commit(
					ctx, []*protocol.FramedProposal{})

				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					b.StopTimer()
					gs := ctx.GroupState("1").Clone()
					b.StartTimer()
					gs.Process(ctx, c0, cs["1"], []*protocol.FramedProposal{})
				}
			})
		}
	}

}

func TestProtocol(t *testing.T) {
	ctx := protocol.CreateContext(nil)
	ctx.CreateMember("creator")
	ctx.CreateGroup("creator")
	ctx.CommitProcessAndJoin("creator")

	ctx.CreateMember("joiner1")
	ctx.ProposeAdd("creator", "joiner1")
	ctx.CommitProcessAndJoin("creator")

	ctx.ProposeUpd("creator")
	ctx.CommitProcessAndJoin("joiner1")

	ctx.CreateMember("joiner2")
	ctx.CreateMember("joiner3")
	ctx.ProposeAdd("joiner1", "joiner2")
	ctx.ProposeAdd("creator", "joiner3")
	ctx.CommitProcessAndJoin("joiner1")

	ctx.CreateMember("joiner4")
	ctx.ProposeRem("joiner1", "joiner2")
	ctx.ProposeAdd("joiner1", "joiner4")

	ctx.CommitProcessAndJoin("joiner3")
	ctx.CommitProcessAndJoin("joiner4")

	ctx.ProposeUpd("joiner1")
	ctx.CommitProcessAndJoin("creator")
}
