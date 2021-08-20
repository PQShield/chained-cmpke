package protocol

// TODO cache intermediate computations in mKEM's Group, etc

import (
	"bytes"
	"crypto"
	"crypto/aes"
	cryptoRand "crypto/rand"
	"crypto/subtle"
	"encoding/binary"
	"fmt"
	"sort"

	"golang.org/x/crypto/sha3"
	sign "github.com/PQShield/chained-cmpke/dilithium2"

	"github.com/PQShield/chained-cmpke/mkem"
	mkemSchemes "github.com/PQShield/chained-cmpke/mkem/schemes"
)

const (
	Kappa = 16

	MemberHashId = iota
	ConfTransHashId
	PropIdHashId
	MacHashId
	CSKEHashId
	EpochKeysHashId
	JoinerSecretHashId
	InterimTransHashId
	FramedCommitKeyHashId
	KeyPackageHashId
	CommitSignId
	KeyPackageSignId
)

type GroupInfo struct {
	GroupId          []byte
	Epoch            uint
	MemberPublicInfo map[string]*KeyPackage
	MemberHash       []byte
	ConfTransHash    []byte
	InterimTransHash []byte
	ConfTag          []byte
	IdC              string
}

func (gi *GroupInfo) Pack(p *Pad) []byte {
	bufs := [][]byte{
		gi.GroupId,
		packUint(gi.Epoch),
		gi.MemberHash,
		gi.ConfTransHash,
		gi.InterimTransHash,
		gi.ConfTag,
		[]byte(gi.IdC),
	}
	ids := []string{}
	for id := range gi.MemberPublicInfo {
		ids = append(ids, id)
	}
	sort.Strings(ids)
	for _, id := range ids {
		bufs = append(bufs, []byte(id), gi.MemberPublicInfo[id].Hash(p))
	}

	return pack(bufs...)
}

// scratchpad to reduce allocations
type Pad struct {
	Svk   [sign.PublicKeySize]byte
	Ek    []byte
	Shake sha3.ShakeHash
}

func (ctx *Context) NewPad() *Pad {
	return &Pad{
		Shake: sha3.NewShake256(),
		Ek:    make([]byte, ctx.mkemScheme.PublicKeySize()),
	}
}

type WelcomeMessageI struct {
	GroupInfo *GroupInfo
	T         []byte
	Sig       []byte
}

type WelcomeMessageD struct {
	KpHash []byte
	Ctd    []byte
}

type MemberState struct {
	Id  string
	Ek  mkem.PublicKey
	Svk *sign.PublicKey
	Sig []byte
	Dk  mkem.PrivateKey
}

func (ms *MemberState) Kp() *KeyPackage {
	return &KeyPackage{
		Id:  ms.Id,
		Ek:  ms.Ek,
		Svk: ms.Svk,
		Sig: ms.Sig,
	}
}

type AS struct {
	Registered map[string]*sign.PublicKey
}

func NewAS() *AS {
	ret := AS{
		Registered: make(map[string]*sign.PublicKey),
	}
	return &ret
}

func (as *AS) RegisterSvk(id string, svk *sign.PublicKey) {
	if _, ok := as.Registered[id]; ok {
		panic("already registered")
	}
	as.Registered[id] = svk
}

func (as *AS) VerifyCert(id string, svk *sign.PublicKey) bool {
	return as.Registered[id].Equal(svk)
}

type KS struct {
	as *AS
	KP map[string]*KeyPackage
}

func NewKS(as *AS) *KS {
	return &KS{
		as: as,
		KP: make(map[string]*KeyPackage),
	}
}

func (ks *KS) RegisterKp(kp *KeyPackage) {
	if !ks.as.VerifyCert(kp.Id, kp.Svk) {
		panic("registerkp: svk mismatch")
	}
	if _, ok := ks.KP[kp.Id]; ok {
		panic("registerkp: kp already registered")
	}
	ks.KP[kp.Id] = kp
}

func (ks *KS) GetKp(id string) *KeyPackage {
	return ks.KP[id]
}

type GroupState struct {
	GroupId          []byte
	Epoch            uint
	ConfTransHash    []byte
	InterimTransHash []byte

	Member     map[string]*MemberState
	MemberHash []byte

	Id  string
	Ssk *sign.PrivateKey
	Svk *sign.PublicKey

	AppSecret  []byte
	MembKey    []byte
	InitSecret []byte

	MKemGroup mkem.Group
	CertSvks  map[string]*sign.PublicKey

	PendUpd map[string]*PendingUpdate
	PendCom map[string]*PendingCommit
}

type PendingUpdate struct {
	Ssk *sign.PrivateKey
	Dk  mkem.PrivateKey
}

type Commit struct {
	PropIds [][]byte
	Kp      *KeyPackage
	T       []byte
}

func (c *Commit) Pack(p *Pad) []byte {
	bufs := [][]byte{
		c.Kp.Hash(p),
		c.T,
	}
	bufs = append(bufs, c.PropIds...)
	return pack(bufs...)
}

type FramedCommit struct {
	GroupId          []byte
	Epoch            uint
	InterimTransHash []byte
	Id               string
	C0               *Commit
	Sig              []byte
	ConfTag          []byte
}

func (fc *FramedCommit) Key(pad *Pad) string {
	return string(hashPack(
		pad,
		FramedCommitKeyHashId,
		fc.GroupId,
		packUint(fc.Epoch),
		fc.InterimTransHash,
		[]byte(fc.Id),
		fc.C0.Pack(pad),
		fc.Sig,
		fc.ConfTag,
	))
}

func (gs *GroupState) Clone() *GroupState {
	ret := *gs
	ret.CertSvks = make(map[string]*sign.PublicKey)
	for id, svk := range gs.CertSvks {
		ret.CertSvks[id] = svk
	}
	ret.Member = make(map[string]*MemberState)
	for k, v := range gs.Member {
		v2 := *v
		ret.Member[k] = &v2
	}
	// TODO do we need to clone these as init-epoch clears them?
	// TODO are these cloned completely?
	// ret.PendUpd = make(map[string]*PendingUpdate)
	// for k, v := range gs.PendUpd {
	//     v2 := *v
	//     ret.PendUpd[k] = &v2
	// }
	// ret.PendCom = make(map[string]*PendingCommit)
	// for k, v := range gs.PendCom {
	//     v2 := *v
	//     ret.PendCom[k] = &v2
	// }
	return &ret
}

func (gs *GroupState) AssignKp(kp *KeyPackage) {
	member, ok := gs.Member[kp.Id]
	if !ok {
		member = &MemberState{Id: kp.Id}
		gs.Member[kp.Id] = member
	}
	member.Ek = kp.Ek
	member.Svk = kp.Svk
	member.Sig = kp.Sig
}

type KeyPackage struct {
	Id  string
	Ek  mkem.PublicKey
	Svk *sign.PublicKey
	Sig []byte
}

func (kp *KeyPackage) Hash(pad *Pad) []byte {
	kp.Svk.Pack(&pad.Svk)
	kp.Ek.Pack(pad.Ek)
	return hashPack(
		pad,
		KeyPackageHashId,
		[]byte(kp.Id),
		pad.Ek,
		pad.Svk[:],
		kp.Sig,
	)
}

func (kp *KeyPackage) VerifySig(pad *Pad) {
	kp.Svk.Pack(&pad.Svk)
	kp.Ek.Pack(pad.Ek)
	if !sign.Verify(
		kp.Svk,
		hashPack(pad, KeyPackageSignId, []byte(kp.Id), pad.Ek, pad.Svk[:]),
		kp.Sig,
	) {
		panic("keypackage signature invalid")
	}
}

func randomKey() []byte {
	ret := make([]byte, Kappa)
	_, err := cryptoRand.Read(ret)
	if err != nil {
		panic(err)
	}
	return ret
}

func hashPack(pad *Pad, sep byte, in ...[]byte) []byte {
	pad.Shake.Reset()
	h := pad.Shake
	var out [Kappa]byte
	h.Write([]byte{sep})
	var buf [4]byte

	for i := 0; i < len(in); i++ {
		binary.BigEndian.PutUint32(buf[:], uint32(len(in[i])))
		h.Write(buf[:])
	}

	for i := 0; i < len(in); i++ {
		h.Write(in[i][:])
	}

	h.Read(out[:])
	return out[:]
}

// Injection of 2** -> 2*
func pack(in ...[]byte) []byte {
	N := 0
	for i := 0; i < len(in); i++ {
		N += len(in[i]) + 4
		if len(in[i]) >= 2147483647 {
			panic("too big for pack")
		}
	}
	buf := make([]byte, N)
	ret := buf
	for i := 0; i < len(in); i++ {
		binary.BigEndian.PutUint32(buf[:], uint32(len(in[i])))
		buf = buf[4:]
	}
	for i := 0; i < len(in); i++ {
		copy(buf[:], in[i])
		buf = buf[len(in[i]):]
	}
	return ret
}

func genKp(id string, ssk *sign.PrivateKey, svk *sign.PublicKey,
	mg mkem.Group, pad *Pad) (*KeyPackage, mkem.PrivateKey) {
	dk, ek := mg.GenerateKeypair()
	ek.Pack(pad.Ek)
	svk.Pack(&pad.Svk)
	sig, err := ssk.Sign(
		nil,
		hashPack(pad, KeyPackageSignId, []byte(id), pad.Ek, pad.Svk[:]),
		crypto.Hash(0),
	)
	if err != nil {
		panic(err)
	}
	return &KeyPackage{
		Id:  id,
		Ek:  ek,
		Svk: svk,
		Sig: sig,
	}, dk
}

func (gs *GroupState) GenKp(pad *Pad) (*KeyPackage, mkem.PrivateKey) {
	return genKp(gs.Id, gs.Ssk, gs.Svk, gs.MKemGroup, pad)
}

func (gs *GroupState) ValidateKp(ctx *Context, kp *KeyPackage, id string, pad *Pad) {
	if id != kp.Id {
		panic(fmt.Sprintf("validate-pk: mismatch in ids: %s ≠ %s", id, kp.Id))
	}
	kp.VerifySig(pad) // panics if sig isn't correct
	svk, ok := gs.CertSvks[id]
	if ok {
		if !kp.Svk.Equal(svk) {
			panic("verified svk doesn't match kp svk")
		}
	} else {
		if !ctx.as.VerifyCert(id, kp.Svk) {
			panic("AS didn't verify certificate")
		}
		gs.CertSvks[id] = kp.Svk
	}
}

func (gs *GroupState) SetInterimTransHash(confTag []byte, pad *Pad) {
	gs.InterimTransHash = hashPack(
		pad,
		InterimTransHashId,
		gs.ConfTransHash,
		confTag,
	)
}

func (gs *GroupState) GroupCont() []byte {
	return pack(
		gs.GroupId,
		packUint(gs.Epoch),
		gs.MemberHash,
		gs.ConfTransHash,
	)
}

func (gs *GroupState) GroupContWInterim() []byte {
	return pack(
		gs.GroupId,
		packUint(gs.Epoch),
		gs.MemberHash,
		gs.ConfTransHash,
		gs.InterimTransHash,
	)
}

func (gs *GroupState) MemberIds() []string {
	ret := []string{}
	for id := range gs.Member {
		ret = append(ret, id)
	}
	return ret
}

func (gs *GroupState) DeriveMemberHash(pad *Pad) []byte {
	ids := gs.MemberIds()
	sort.Strings(ids)
	packedKps := [][]byte{}
	for _, id := range ids {
		packedKps = append(packedKps, gs.Member[id].Kp().Hash(pad))
	}
	return hashPack(pad, MemberHashId, packedKps...)
}

type Member struct {
	gs *GroupState

	// These are just the initial keys and might go out of date.
	ssk *sign.PrivateKey
	svk *sign.PublicKey
	dk  mkem.PrivateKey
}

type Context struct {
	as         *AS
	ks         *KS
	members    map[string]*Member
	mkemScheme mkem.Scheme
	mkemGroup  mkem.Group

	props []*FramedProposal
}

func (ctx *Context) GroupState(id string) *GroupState {
	return ctx.members[id].gs
}

type ContextOpts struct {
	MKemScheme string
}

func CreateContext(opts *ContextOpts) *Context {
	if opts == nil {
		opts = &ContextOpts{}
	}
	if opts.MKemScheme == "" {
		opts.MKemScheme = "ntrulpr"
	}
	mkemScheme := mkemSchemes.ByName(opts.MKemScheme)
	if mkemScheme == nil {
		panic("no such mkem scheme")
	}
	as := NewAS()
	return &Context{
		as:         as,
		ks:         NewKS(as),
		members:    make(map[string]*Member),
		mkemScheme: mkemScheme,
		props:      []*FramedProposal{},
	}
}

func (ctx *Context) ProposeUpd(id string) {
	ctx.props = append(ctx.props, ctx.members[id].gs.ProposeUpd(ctx))
}

func (ctx *Context) ProposeRem(idS, idT string) {
	ctx.props = append(ctx.props, ctx.members[idS].gs.ProposeRem(ctx, idT))
}

func (ctx *Context) ProposeAdd(idS, idT string) {
	ctx.props = append(ctx.props, ctx.members[idS].gs.ProposeAdd(ctx, idT))
}

func (gs *GroupState) ProposeRem(ctx *Context, id string) *FramedProposal {
	if _, ok := gs.Member[id]; !ok {
		panic("not member")
	}
	return gs.FrameProposal(&Proposal{
		Type: "rem",
		Id:   id,
	}, ctx.NewPad())
}

func (gs *GroupState) ProposeUpd(ctx *Context) *FramedProposal {
	// TODO update signing verification key
	pad := ctx.NewPad()
	kp, dk := genKp(gs.Id, gs.Ssk, gs.Svk, ctx.mkemGroup, pad)

	fp := gs.FrameProposal(&Proposal{
		Type: "upd",
		Kp:   kp,
	}, pad)
	gs.PendUpd[string(fp.PropId(pad))] = &PendingUpdate{
		Ssk: gs.Ssk,
		Dk:  dk,
	}
	return fp
}

func (gs *GroupState) ProposeAdd(ctx *Context, id string) *FramedProposal {
	if _, ok := gs.Member[id]; ok {
		panic("already member")
	}
	pad := ctx.NewPad()
	kp := ctx.ks.GetKp(id)
	gs.ValidateKp(ctx, kp, id, pad)
	return gs.FrameProposal(&Proposal{
		Type: "add",
		Kp:   kp,
	}, pad)
}

func (gs *GroupState) InitializeGroup(id string, info *GroupInfo) {
	gs.GroupId = info.GroupId
	gs.Epoch = info.Epoch
	gs.MemberHash = info.MemberHash
	gs.ConfTransHash = info.ConfTransHash
	gs.InterimTransHash = info.InterimTransHash

	gs.Member = make(map[string]*MemberState)

	for idM, kp := range info.MemberPublicInfo {
		gs.Member[idM] = &MemberState{
			Id:  idM,
			Sig: kp.Sig,
			Svk: kp.Svk,
			Ek:  kp.Ek,
		}
	}

	gs.CertSvks = make(map[string]*sign.PublicKey)
	gs.PendUpd = make(map[string]*PendingUpdate)
	gs.PendCom = make(map[string]*PendingCommit)
	gs.Id = id
}

func (ctx *Context) Join(w0 *WelcomeMessageI, w *WelcomeMessageD,
	id string) *GroupState {
	pad := ctx.NewPad()
	gs := &GroupState{}
	gs.InitializeGroup(id, w0.GroupInfo)
	gs.MKemGroup = ctx.mkemGroup
	if !sign.Verify(
		gs.Member[w0.GroupInfo.IdC].Svk,
		pack(w0.GroupInfo.Pack(pad), w0.T),
		w0.Sig,
	) {
		panic("group info sig invalid")
	}

	gs.VrfGroupState(ctx, pad)
	gs.Ssk = ctx.members[id].ssk
	gs.Svk = ctx.members[id].svk
	gs.Member[id].Dk = ctx.members[id].dk

	if !bytes.Equal(
		w.KpHash,
		gs.Member[id].Kp().Hash(pad),
	) {
		panic("kp mismatch")
	}
	joinerSecret := CmDecrypt(w0.T, w.Ctd, gs.Member[id].Dk, gs.Member[id].Ek)
	confKey := gs.DeriveEpochKeys(joinerSecret, pad)
	gs.VrfConfTag(confKey, w0.GroupInfo.ConfTag)

	return gs
}

func (gs *GroupState) Process(ctx *Context, c0 *FramedCommit,
	ctd []byte, fps []*FramedProposal) (idC string, changes *Changes) {
	pad := ctx.NewPad()
	idC, C0, sig, confTag := gs.UnframeCommit(c0, pad)
	if idC == gs.Id {
		pc := gs.PendCom[c0.Key(pad)]

		ok := len(fps) == len(pc.Fps)
		for i := 0; i < len(fps) && ok; i++ {
			ok = bytes.Equal(fps[i].PropId(pad), pc.Fps[i].PropId(pad))
		}
		if !ok {
			panic("proposals mismatch")
		}

		*gs = *pc.New
		return idC, pc.Changes
	}

	for i := 0; i < len(fps); i++ {
		if !bytes.Equal(C0.PropIds[i], fps[i].PropId(pad)) {
			panic("propids mismatch")
		}
	}

	old := gs.InitEpoch()
	changes = gs.ApplyProps(ctx, fps, old, pad)
	if changes.IsRemoved(idC) || changes.IsUpdated(idC) {
		panic("commiter can't remove or update themselves")
	}

	if changes.IsRemoved(gs.Id) {
		*gs = GroupState{}
		return idC, nil
	}

	gs.SetConfTransHash(old, idC, C0, sig, pad)
	comSecret := gs.ApplyRekey(ctx, idC, C0.Kp, C0.T, ctd, pad)
	gs.MemberHash = gs.DeriveMemberHash(pad)
	confKey, _ := gs.DeriveKeys(old, comSecret, pad)
	gs.VrfConfTag(confKey, confTag)
	gs.SetInterimTransHash(confTag, pad)
	return idC, changes
}

func (gs *GroupState) ApplyRekey(ctx *Context, idC string, kp *KeyPackage,
	T []byte, ctd []byte, pad *Pad) []byte {
	comSecret := CmDecrypt(T, ctd, gs.Member[gs.Id].Dk, gs.Member[gs.Id].Ek)
	gs.ValidateKp(ctx, kp, idC, pad)
	gs.AssignKp(kp)
	return comSecret
}

func (gs *GroupState) VrfGroupState(ctx *Context, pad *Pad) {
	if !bytes.Equal(gs.DeriveMemberHash(pad), gs.MemberHash) {
		panic(fmt.Sprintf("mismatched memberhash: %x %x",
			gs.DeriveMemberHash(pad), gs.MemberHash))
	}
	for id, member := range gs.Member {
		gs.ValidateKp(ctx, member.Kp(), id, pad)
	}
}

func (ctx *Context) CreateGroup(idCreator string) *GroupState {
	pad := ctx.NewPad()
	member, ok := ctx.members[idCreator]
	if !ok {
		panic("creator does not exist")
	}
	if member.gs != nil {
		panic("creator already in group")
	}
	groupState := &GroupState{
		GroupId:    randomKey(),
		InitSecret: randomKey(),
		Epoch:      0,

		MKemGroup: ctx.mkemScheme.GenerateGroup(),
		Id:        idCreator,
		Ssk:       member.ssk,

		Member:   make(map[string]*MemberState),
		CertSvks: make(map[string]*sign.PublicKey),
		PendUpd:  make(map[string]*PendingUpdate),
		PendCom:  make(map[string]*PendingCommit),
	}
	ctx.members[idCreator].gs = groupState
	groupState.Svk = member.svk
	kp, dk := groupState.GenKp(pad)
	groupState.Member[idCreator] = &MemberState{Id: idCreator}
	groupState.AssignKp(kp)
	groupState.Member[idCreator].Dk = dk

	ctx.mkemGroup = groupState.MKemGroup
	return groupState
}

func (gs *GroupState) InitEpoch() *GroupState {
	old := gs.Clone()

	gs.Epoch++

	gs.PendCom = make(map[string]*PendingCommit)
	gs.PendUpd = make(map[string]*PendingUpdate)

	return old
}

func (gs *GroupState) ApplyProps(ctx *Context, fps []*FramedProposal,
	old *GroupState, pad *Pad) *Changes {
	ret := Changes{
		Add:    []AddChange{},
		Update: []UpdateChange{},
		Remove: []RemoveChange{},
	}

	updateLut := make(map[string]struct{})

	// TODO first updates, then remove, then add

	for _, fp := range fps {
		idS, p := old.UnframeProposal(fp, pad)
		switch p.Type {
		case "upd":
			if _, ok := old.Member[idS]; !ok {
				panic("no such member")
			}

			if len(ret.Remove) != 0 || len(ret.Add) != 0 {
				panic("wrong order")
			}

			if _, ok := updateLut[idS]; ok {
				panic("already updated")
			}

			gs.ValidateKp(ctx, p.Kp, idS, pad)
			gs.AssignKp(p.Kp)

			if idS == gs.Id {
				pu := old.PendUpd[string(fp.PropId(pad))]
				gs.Ssk = pu.Ssk
				gs.Member[gs.Id].Dk = pu.Dk
				gs.Svk = p.Kp.Svk
			}

			updateLut[idS] = struct{}{}
			ret.Update = append(ret.Update, UpdateChange{
				Id:  idS,
				Svk: gs.Member[idS].Svk,
			})
		case "rem":
			idT := p.Id
			if idS == idT {
				panic("user can't remove themselves")
			}
			if _, ok := gs.Member[idT]; !ok {
				panic("no such member")
			}
			if _, ok := updateLut[idT]; ok {
				panic("can't remove updated member")
			}
			if len(ret.Add) != 0 {
				panic("wrong order")
			}
			delete(gs.Member, idT)
			ret.Remove = append(ret.Remove, RemoveChange{
				IdT: idT,
				IdS: idS,
			})
		case "add":
			idT := p.Kp.Id
			if _, ok := gs.Member[idT]; ok {
				panic("member already exists")
			}
			gs.ValidateKp(ctx, p.Kp, idT, pad)
			gs.AssignKp(p.Kp)
			ret.Add = append(ret.Add, AddChange{
				IdT: idT,
				IdS: idS,
				Svk: p.Kp.Svk,
			})
		default:
			panic("unknown type of proposal")
		}
	}

	return &ret
}

func (ctx *Context) CommitProcessAndJoin(idC string) {
	gs := ctx.members[idC].gs
	fps := ctx.props
	ctx.props = []*FramedProposal{}
	c0, cs, w0, ws := gs.Commit(ctx, fps)
	for id, w := range ws {
		ctx.members[id].gs = ctx.Join(w0, w, id)
	}
	for id, ctd := range cs {
		ctx.members[id].gs.Process(ctx, c0, ctd, fps)
	}
}

func (gs *GroupState) Commit(ctx *Context, fps []*FramedProposal) (
	c0 *FramedCommit,
	cs map[string][]byte,
	w0 *WelcomeMessageI,
	ws map[string]*WelcomeMessageD,
) {
	old := gs.InitEpoch()
	pad := ctx.NewPad()

	changes := gs.ApplyProps(ctx, fps, old, pad)

	if changes.IsRemoved(gs.Id) {
		panic("can't commit our own removal")
	}
	if changes.IsUpdated(gs.Id) {
		panic("can't commit our own update")
	}

	receivers := []string{}
	addedMem := []string{}
	for id := range gs.Member {
		if changes.IsAdded(id) {
			addedMem = append(addedMem, id)
		} else {
			receivers = append(receivers, id)
		}
	}

	comSecret, kp, T, ctd := gs.Rekey(receivers, pad)

	gs.MemberHash = gs.DeriveMemberHash(pad)

	C0 := Commit{
		Kp:      kp,
		T:       T,
		PropIds: [][]byte{},
	}

	for _, fp := range fps {
		C0.PropIds = append(C0.PropIds, fp.PropId(pad))
	}

	sig := old.SignCommit(&C0, pad)
	gs.SetConfTransHash(old, gs.Id, &C0, sig, pad)
	confKey, joinerSecret := gs.DeriveKeys(old, comSecret, pad)
	confTag := gs.GenConfTag(confKey)
	c0 = old.FrameCommit(&C0, sig, confTag)
	gs.SetInterimTransHash(confTag, pad)

	cs = make(map[string][]byte) // recipient ciphertexts
	for i, id := range receivers {
		cs[id] = ctd[i]
	}

	if len(changes.Add) != 0 {
		w0, ws = gs.WelcomeMsg(addedMem, joinerSecret, confTag, pad)
	}

	old.PendCom[c0.Key(pad)] = &PendingCommit{
		New:     gs.Clone(),
		Fps:     fps,
		Changes: changes,
	}

	*gs = *old

	return
}

func (gs *GroupState) WelcomeMsg(addedMem []string, joinerSecret,
	confTag []byte, pad *Pad) (w0 *WelcomeMessageI, ws map[string]*WelcomeMessageD) {
	var err error
	var ctd [][]byte

	eks := make([]mkem.PublicKey, len(addedMem))
	for i, id := range addedMem {
		eks[i] = gs.Member[id].Ek
	}

	w0 = new(WelcomeMessageI)
	w0.T, ctd = CmEncrypt(joinerSecret, eks)
	w0.GroupInfo = &GroupInfo{
		GroupId:          gs.GroupId,
		Epoch:            gs.Epoch,
		MemberPublicInfo: gs.MemberPublicInfo(),
		MemberHash:       gs.MemberHash,
		ConfTransHash:    gs.ConfTransHash,
		InterimTransHash: gs.InterimTransHash,
		ConfTag:          confTag,
		IdC:              gs.Id,
	}

	w0.Sig, err = gs.Ssk.Sign(
		nil,
		pack(w0.GroupInfo.Pack(pad), w0.T),
		crypto.Hash(0),
	)
	if err != nil {
		panic(err)
	}

	ws = make(map[string]*WelcomeMessageD)
	for i, id := range addedMem {
		ws[id] = &WelcomeMessageD{
			KpHash: gs.Member[id].Kp().Hash(pad),
			Ctd:    ctd[i],
		}
	}

	return
}

func (gs *GroupState) MemberPublicInfo() map[string]*KeyPackage {
	ret := make(map[string]*KeyPackage)
	for id, member := range gs.Member {
		ret[id] = member.Kp()
	}
	return ret
}

func (gs *GroupState) FrameCommit(commit *Commit, sig, confTag []byte) *FramedCommit {
	return &FramedCommit{
		GroupId:          gs.GroupId,
		Epoch:            gs.Epoch,
		InterimTransHash: gs.InterimTransHash,
		Id:               gs.Id,
		C0:               commit,
		Sig:              sig,
		ConfTag:          confTag,
	}
}

func (gs *GroupState) DeriveEpochKeys(joinerSecret []byte, pad *Pad) []byte {
	groupCont := gs.GroupCont()
	gs.AppSecret = hashPack(
		pad,
		EpochKeysHashId,
		joinerSecret,
		groupCont,
		[]byte("app"),
	)
	gs.MembKey = hashPack(
		pad,
		EpochKeysHashId,
		joinerSecret,
		groupCont,
		[]byte("memb"),
	)
	gs.InitSecret = hashPack(
		pad,
		EpochKeysHashId,
		joinerSecret,
		groupCont,
		[]byte("init"),
	)
	return hashPack(
		pad,
		EpochKeysHashId,
		joinerSecret,
		groupCont,
		[]byte("conf"),
	)
}

func (gs *GroupState) DeriveKeys(old *GroupState, comSecret []byte, pad *Pad) ([]byte, []byte) {
	joinerSecret := hashPack(
		pad,
		JoinerSecretHashId,
		old.InitSecret,
		comSecret,
	)
	confKey := gs.DeriveEpochKeys(joinerSecret, pad)
	return confKey, joinerSecret
}

type Proposal struct {
	Type string // add, rem, upd

	Kp *KeyPackage // in case of upd or add
	Id string      // in case of remove
}

type AddChange struct {
	IdS string
	IdT string
	Svk *sign.PublicKey
}

type UpdateChange struct {
	Id  string
	Svk *sign.PublicKey
}

type RemoveChange struct {
	IdS string
	IdT string
}

type Changes struct {
	Add    []AddChange
	Update []UpdateChange
	Remove []RemoveChange
}

func (cs *Changes) IsAdded(id string) bool {
	for _, c := range cs.Add {
		if c.IdT == id {
			return true
		}
	}
	return false
}

func (cs *Changes) IsRemoved(id string) bool {
	for _, c := range cs.Remove {
		if c.IdT == id {
			return true
		}
	}
	return false
}

func (cs *Changes) IsUpdated(id string) bool {
	for _, c := range cs.Update {
		if c.Id == id {
			return true
		}
	}
	return false
}

type PendingCommit struct {
	New     *GroupState
	Fps     []*FramedProposal
	Changes *Changes
}

func (p *Proposal) Pack(pad *Pad) []byte {
	var arg []byte
	if p.Type == "rem" {
		arg = []byte(p.Id)
	} else if p.Type == "add" || p.Type == "upd" {
		arg = p.Kp.Hash(pad)
	} else {
		panic("invalid type")
	}
	return pack(
		[]byte(p.Type),
		arg,
	)
}

func (gs *GroupState) SetConfTransHash(old *GroupState, idC string,
	C0 *Commit, sig []byte, pad *Pad) {
	gs.ConfTransHash = hashPack(
		pad,
		ConfTransHashId,
		old.InterimTransHash,
		old.GroupId,
		packUint(old.Epoch),
		[]byte(idC),
		[]byte("commit"),
		C0.Pack(pad),
		sig,
	)
}

func (gs *GroupState) SignCommit(c *Commit, pad *Pad) []byte {
	groupContWInterim := gs.GroupContWInterim()
	packedC := c.Pack(pad)
	comCont := hashPack(
		pad,
		CommitSignId,
		groupContWInterim,
		[]byte(gs.Id),
		[]byte("commit"),
		packedC,
	)
	sig, err := gs.Ssk.Sign(
		nil,
		comCont,
		crypto.Hash(0),
	)
	if err != nil {
		panic(err)
	}
	return sig
}

func (gs *GroupState) FrameProposal(p *Proposal, pad *Pad) *FramedProposal {
	propCont := pack(
		gs.GroupContWInterim(),
		[]byte(gs.Id),
		[]byte("proposal"),
		p.Pack(pad),
	)

	sig, err := gs.Ssk.Sign(
		nil,
		propCont,
		crypto.Hash(0),
	)
	if err != nil {
		panic(err)
	}

	membTag := macTag(gs.MembKey, pack(propCont, sig))

	return &FramedProposal{
		GroupId:          gs.GroupId,
		Epoch:            gs.Epoch,
		InterimTransHash: gs.InterimTransHash,
		Id:               gs.Id,
		P:                *p,
		Sig:              sig,
		MembTag:          membTag,
	}
}

func (gs *GroupState) UnframeCommit(fc *FramedCommit, pad *Pad) (idC string,
	C0 *Commit, sig []byte, confTag []byte) {
	idC = fc.Id
	C0 = fc.C0
	sig = fc.Sig
	confTag = fc.ConfTag
	if !bytes.Equal(fc.GroupId, gs.GroupId) {
		panic("groupid mismatch")
	}
	if fc.Epoch != gs.Epoch {
		panic(fmt.Sprintf("epoch mismatch: %d ≠ %d", fc.Epoch, gs.Epoch))
	}
	if !bytes.Equal(fc.InterimTransHash, gs.InterimTransHash) {
		panic("interimTransHash mismatch")
	}
	groupContWInterim := gs.GroupContWInterim()
	packedC0 := C0.Pack(pad)
	comCont := hashPack(
		pad,
		CommitSignId,
		groupContWInterim,
		[]byte(idC),
		[]byte("commit"),
		packedC0,
	)
	if !sign.Verify(
		gs.Member[idC].Svk,
		comCont,
		sig,
	) {
		panic("invalid sig")
	}
	return
}

func (gs *GroupState) UnframeProposal(f *FramedProposal, pad *Pad) (string, Proposal) {
	if !bytes.Equal(gs.GroupId, f.GroupId) {
		panic("proposal groupid does not match")
	}
	if gs.Epoch != f.Epoch {
		panic(fmt.Sprintf("frame epoch does not match %d ≠ %d", gs.Epoch,
			f.Epoch))
	}
	if !bytes.Equal(gs.InterimTransHash, f.InterimTransHash) {
		panic("frame interim hash does not match") // TODO this seems redundant
	}
	if _, ok := gs.Member[f.Id]; !ok {
		panic("don't know that member")
	}
	propCont := pack(
		gs.GroupContWInterim(),
		[]byte(f.Id),
		[]byte("proposal"),
		f.P.Pack(pad),
	)
	if !sign.Verify(
		gs.Member[f.Id].Svk,
		propCont,
		f.Sig,
	) {
		panic("invalid sig")
	}
	macVerify(gs.MembKey, pack(propCont, f.Sig), f.MembTag)
	return f.Id, f.P
}

func (gs *GroupState) GenConfTag(confKey []byte) []byte {
	return macTag(confKey, gs.ConfTransHash)
}

func (gs *GroupState) VrfConfTag(confKey, confTag []byte) {
	macVerify(confKey, gs.ConfTransHash, confTag)
}

func macVerify(K, M, T []byte) {
	if subtle.ConstantTimeCompare(macTag(K, M), T) != 1 {
		panic("mac does not verify")
	}
}

func macTag(K, M []byte) []byte {
	ret := make([]byte, Kappa)
	h := sha3.NewShake256()
	h.Write([]byte{MacHashId})
	h.Write(K)
	h.Write(M)
	h.Read(ret)
	return ret
}

type FramedProposal struct {
	GroupId          []byte
	Epoch            uint
	InterimTransHash []byte
	Id               string
	P                Proposal
	Sig              []byte
	MembTag          []byte
}

func packUint(x uint) []byte {
	var buf [8]byte
	binary.BigEndian.PutUint64(buf[:], uint64(x))
	return buf[:]
}

func (fp *FramedProposal) PropId(pad *Pad) []byte {
	return hashPack(
		pad,
		PropIdHashId,
		fp.GroupId,
		packUint(fp.Epoch),
		fp.InterimTransHash,
		[]byte(fp.Id),
		fp.P.Pack(pad),
		fp.Sig,
		fp.MembTag,
	)
}

func (gs *GroupState) Rekey(receivers []string, pad *Pad) (
	comSecret []byte,
	kp *KeyPackage,
	cti []byte,
	ctd [][]byte) {
	var dk mkem.PrivateKey
	kp, dk = gs.GenKp(pad)
	gs.AssignKp(kp)
	gs.Member[gs.Id].Dk = dk
	comSecret = randomKey()

	eks := make([]mkem.PublicKey, len(receivers))
	for i, id := range receivers {
		eks[i] = gs.Member[id].Ek
	}

	cti, ctd = CmEncrypt(comSecret, eks)

	return
}

func CmDecrypt(cti []byte, ctd []byte, dk mkem.PrivateKey,
	ek mkem.PublicKey) []byte {
	group := dk.Group()
	scheme := group.Scheme()

	ss := make([]byte, scheme.SharedSecretSize())

	kemCti := cti[2*Kappa:]
	symCti := cti[:2*Kappa]

	group.Decapsulate(ss, kemCti, ctd, dk, ek)
	return CSKEDec(symCti, ss)
}

func CmEncrypt(M []byte, eks []mkem.PublicKey) (cti []byte, ctd [][]byte) {
	mGroup := eks[0].Group()
	mScheme := mGroup.Scheme()

	ss := make([]byte, mScheme.SharedSecretSize())
	ctd = make([][]byte, len(eks))

	cti = make([]byte, mScheme.IndependentCiphertextSize()+2*Kappa)
	kemCti := cti[2*Kappa:]
	symCti := cti[:2*Kappa]

	for i := 0; i < len(eks); i++ {
		ctd[i] = make([]byte, mScheme.DependentCiphertextSize())
	}

	mGroup.Encapsulate(ss, kemCti, ctd, eks)

	CSKEEnc(symCti, ss, M)
	return
}

func CSKEDec(ct, K []byte) []byte {
	var K1, K2, pt [Kappa]byte

	h := sha3.NewShake256()
	h.Write([]byte{CSKEHashId})
	h.Write(K)
	h.Read(K1[:])
	h.Read(K2[:])
	cipher, err := aes.NewCipher(K1[:])
	if err != nil {
		panic(err)
	}
	cipher.Decrypt(pt[:], ct[:Kappa])
	if !bytes.Equal(K2[:], ct[Kappa:]) {
		panic("tag mismatch")
	}
	return pt[:]
}

func CSKEEnc(ct, K, pt []byte) {
	var K1 [Kappa]byte

	h := sha3.NewShake256()
	h.Write([]byte{CSKEHashId})
	h.Write(K)
	h.Read(K1[:])
	h.Read(ct[Kappa:])
	cipher, err := aes.NewCipher(K1[:])
	if err != nil {
		panic(err)
	}
	cipher.Encrypt(ct[:Kappa], pt)
}

func (ctx *Context) CreateMember(id string) {
	if _, ok := ctx.members[id]; ok {
		panic("already created")
	}
	pad := ctx.NewPad()
	svk, ssk, err := sign.GenerateKey(nil)
	if err != nil {
		panic(err)
	}
	m := &Member{
		ssk: ssk,
		svk: svk,
	}
	ctx.as.RegisterSvk(id, svk)
	ctx.members[id] = m
	if ctx.mkemGroup != nil {
		kp, dk := genKp(id, ssk, svk, ctx.mkemGroup, pad)
		ctx.ks.RegisterKp(kp)
		m.dk = dk
	}
}
