package main

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestSplit20mer(t *testing.T) {
	a, b, c := split20mer("CCCCCCCCCCBBBBBAAAAA")
	assert.Equal(t, a, "AAAAA")
	assert.Equal(t, b, "BBBBB")
	assert.Equal(t, c, "CCCCCCCCCC")
}

func TestBuildIndex1(t *testing.T) {
	c := make(chan string)
	go func() {
		c <- "CCCCCCCCCCBBBBBAAAAA"
		close(c)
	}()
	i := build_index(c)
	fmt.Println(i)
	assert.Len(t, i, 1)
	assert.Contains(t, i, "AAAAA")
	assert.Contains(t, i["AAAAA"], "BBBBB")
	assert.Contains(t, i["AAAAA"]["BBBBB"], "CCCCCCCCCC")
}

func TestBuildIndex2(t *testing.T) {
	c := make(chan string)
	go func() {
		c <- "CCCCCCCCCCBBBBBAAAAA"
		c <- "CCCCCCCCCABBBBBAAAAA"
		close(c)
	}()
	i := build_index(c)
	fmt.Println(i)
	assert.Len(t, i, 1)
	assert.Contains(t, i, "AAAAA")
	assert.Len(t, i["AAAAA"], 1)
	assert.Contains(t, i["AAAAA"], "BBBBB")
	assert.Len(t, i["AAAAA"]["BBBBB"], 2)
	assert.Contains(t, i["AAAAA"]["BBBBB"], "CCCCCCCCCC")
	assert.Contains(t, i["AAAAA"]["BBBBB"], "CCCCCCCCCA")
}

func TestBuildIndex3(t *testing.T) {
	c := make(chan string)
	go func() {
		c <- "CCCCCCCCCCBBBBBAAAAA"
		c <- "CCCCCCCCCDBBBBDAAAAA"
		close(c)
	}()
	i := build_index(c)
	fmt.Println(i)
	assert.Len(t, i, 1)
	assert.Contains(t, i, "AAAAA")
	assert.Len(t, i["AAAAA"], 2)
	assert.Contains(t, i["AAAAA"], "BBBBB")
	assert.Contains(t, i["AAAAA"], "BBBBD")
	assert.Len(t, i["AAAAA"]["BBBBB"], 1)
	assert.Len(t, i["AAAAA"]["BBBBD"], 1)
	assert.Contains(t, i["AAAAA"]["BBBBB"], "CCCCCCCCCC")
	assert.Contains(t, i["AAAAA"]["BBBBD"], "CCCCCCCCCD")
}

func TestCountSame(t *testing.T) {
	var tests = []struct {
		l        string
		r        string
		expected int
	}{
		{"AAAAA", "AAAAA", 5},
		{"AAAAA", "AAAAB", 4},
		{"AAAAA", "AAACB", 3},
		{"AAAAA", "BBBBB", 0},
	}

	for _, tt := range tests {
		actual := countSame(tt.l, tt.r)
		assert.Equal(t, tt.expected, actual)
	}
}

func TestMatchForward(t *testing.T) {
	c := make(chan string)
	go func() {
		c <- "CCCCCCCCCCBBBBBAAAAA"
		c <- "CCCCCCCCCCBBBBBAAAAA"
		close(c)
	}()
	m := NewOfftargetMatcher(c)
	b := m.MatchForward("CCCCCCCCCCBBBBBAAAAA", 5, 10, 20)
	assert.True(t, b)
	b1 := m.MatchForward("DCCCCCCCCCBBBBBAAAAA", 5, 10, 19)
	assert.True(t, b1)
	b2 := m.MatchForward("DCCCCCCCCCBBBBBAAAAA", 5, 10, 20)
	assert.False(t, b2)
}
