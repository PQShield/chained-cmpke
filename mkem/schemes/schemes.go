package schemes

import (
	"github.com/PQShield/chained-cmpke/mkem"
	"github.com/PQShield/chained-cmpke/mkem/bilbo"
	"github.com/PQShield/chained-cmpke/mkem/ilum"
	"github.com/PQShield/chained-cmpke/mkem/ntrulpr"
	"github.com/PQShield/chained-cmpke/mkem/sike"
	"strings"
)

var allSchemes = [...]mkem.Scheme{
	bilbo.Scheme(),
	ilum.Scheme(),
	sike.Scheme(),
	ntrulpr.Scheme(),
}

var allSchemesNames map[string]mkem.Scheme

func init() {
	allSchemesNames = make(map[string]mkem.Scheme)
	for _, scheme := range allSchemes {
		allSchemesNames[strings.ToLower(scheme.Name())] = scheme
	}
}

func ByName(name string) mkem.Scheme {
	return allSchemesNames[strings.ToLower(name)]
}

func All() []mkem.Scheme {
	a := allSchemes
	return a[:]
}
