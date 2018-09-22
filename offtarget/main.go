package main

import (
	"bufio"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"net/http"
	"net/url"
	"os"
	"strconv"
	"strings"
)

func main() {
	host := os.Getenv("HOST")
	u, e := url.Parse(host)
	if e != nil {
		panic("can't read HOST url")
	}

	var reader io.ReadCloser

	if u.Scheme == "http" {
		if e != nil {
			log.Fatalf("unable to parse host data url %s", host)
		}
		r, e2 := http.Get(u.String())

		if e2 != nil {
			log.Fatalf("unable to load url %s\n error: %#v\n", u.String(), e2)
		}

		log.Printf("building index from %s\n", u.String())
		reader = r.Body
	} else if u.Scheme == "file" {
		f, err := os.Open(u.Path)
		if err != nil {
			panic("can't read file")
		}
		reader = ioutil.NopCloser(bufio.NewReader(f))
	}

	log.Println("ingesting host 20mers")
	m := NewOfftargetMatcher(bufio.NewScanner(reader))
	reader.Close()

	http.HandleFunc("/search", func(w http.ResponseWriter, r *http.Request) {
		t := r.URL.Query().Get("targets")
		l := r.URL.Query().Get("limits")

		targets := strings.Split(t, ",")
		limits := strings.Split(l, ",")
		if len(targets) < 1 {
			http.Error(w, "please specify targets parameter as a comma separated list", 400)
			return
		}
		if len(limits) != 3 {
			http.Error(w, "please specify limits parmeter as a comma separated list of 3 items", 400)
			return
		}

		lim_c5, _ := strconv.Atoi(limits[0])
		lim_c10, _ := strconv.Atoi(limits[1])
		lim_c20, _ := strconv.Atoi(limits[2])

		for _, t := range targets {
			b := m.MatchForward(t, lim_c5, lim_c10, lim_c20)
			// logging here is too expensive
			// log.Printf("limits: %d %d %d targets %#v %t\n", lim_c5, lim_c10, lim_c20, t, b)
			fmt.Fprintf(w, "%s %t\n", t, b)
		}

	})
	log.Println("starting server")
	log.Fatal(http.ListenAndServe(":8080", nil))
}
