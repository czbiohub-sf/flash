package main

import (
	"bufio"
	"log"
)

type twentymer uint64
type tenmer uint32

func basecode(c byte) twentymer {
	switch c {
	case 'A':
		return 0
	case 'T':
		return 3
	case 'C':
		return 1
	case 'G':
		return 2
	}
	panic("We support only ACGT.  No wildcards.")
	return 0xFFFFFFFFFF
}

// In this encoding, the LSBs of the code word correspoind to
// the most signifficant bases in the forward 20mer (i.e.,
// the bases nearest the NGG motif).

func encode(s string) twentymer {
	if len(s) != 20 {
		panic("received target query for string with length != 20")
	}
	var code twentymer = 0
	for i := 0; i < 20; i++ {
		code <<= 2
		code |= basecode(s[i])
	}
	return code
}

func split_in_halves(t twentymer) (tenmer, tenmer) {
	const lsb twentymer = 1
	const lower_half twentymer = (lsb << 20) - lsb
	return tenmer(t & lower_half), tenmer((t >> 20) & lower_half)
}

type OfftargetMatcher struct {
	index [][]tenmer
}

func NewOfftargetMatcher(scanner *bufio.Scanner) *OfftargetMatcher {
	m := OfftargetMatcher{}
	m.index = build_index(scanner)
	return &m
}

func build_index(scanner *bufio.Scanner) [][]tenmer {

	ind := make([][]tenmer, 1024*1024)

	for scanner.Scan() {
		s := scanner.Text()
		t := encode(s)
		ba, c := split_in_halves(t)
		ind[ba] = append(ind[ba], c)
	}

	log.Printf("computing stats\n")
	var s int
	var m int = -1
	var o int = 0
	for _, v := range ind {
		l := len(v)
		s += l
		if l > m {
			m = l
		}
		if l > 0 {
			o += 1
		}
	}
	log.Printf("index contains %d 20mers\n", s)
	log.Printf("max chain length is %d\n", m)
	log.Printf("occupied buckets %d (%d percent)\n", o, (100 * o / len(ind)))
	return ind
}

func count_bits(x tenmer) int {
	// There is an x86 instruction that implements this function.
	// It's available in GO as of 8/25/2017, but our version of
	// GO might predate that release.
	c := 0
	for x != 0 {
		x &= (x - 1)
		c++
	}
	return c
}

func find(haystack []tenmer, needle tenmer, max_differences int) bool {
	if max_differences == 0 {
		// fast path for exact match
		for _, hay := range haystack {
			if hay == needle {
				return true
			}
		}
		return false
	}
	const mask_even_bits tenmer = 0x55555
	const mask_odd_bits tenmer = mask_even_bits << 1
	for _, hay := range haystack {
		xor := hay ^ needle
		xor_even := xor & mask_even_bits
		xor_odd := xor & mask_odd_bits
		xor_any := (xor_odd >> 1) | xor_even
		// now the number of bits set to 1 in xor_any is equal to the number
		// of positional differences in the tenmers hay and needle
		if count_bits(xor_any) <= max_differences {
			return true
		}
	}
	return false
}

func (matcher *OfftargetMatcher) MatchForward(target string, lim_c5, lim_c10, lim_c20 int) bool {
	if lim_c5 != 5 || lim_c10 < 9 || lim_c10 > 10 || lim_c20 < lim_c10 || lim_c20 > 20 {
		panic("For now we only support 5_9_x and 5_10_x")
	}
	t := encode(target)
	ba, c := split_in_halves(t)
	max_diff := 20 - lim_c20
	if find(matcher.index[ba], c, max_diff) {
		return true
	}
	if lim_c10 == 9 {
		// generate every variant of b, total of 15
		for i := uint(5); i < 10; i++ {
			for j := uint(0); j < 4; j++ {
				mask := tenmer(3) << (2 * i)
				ba_prime := (ba &^ mask) | (tenmer(j) << (2 * i))
				if ba != ba_prime {
					if find(matcher.index[ba_prime], c, max_diff-1) {
						return true
					}
				}
			}
		}
	}
	return false
}
