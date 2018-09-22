#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;


// This program scans its input for forward k-3 mers ending with GG,
// or reverse k-3 mers ending with CC.   It filters out guides that
// contain wildcards.  It supercedes crispr2_sites.py.
//
// Usage:
//
//    g++ -O3 --std=c++11 -o crispr_sites2 crispr_sites2.cpp
//    gzip -dc ../../HUGE_DOWNLOADS/hg38.fa.gz | ./crispr_sites2 >! human_guides.txt
//
// Takes about 2 minutes on 2017 MacBook Pro.

// Look for 20-mers at PAM sites.  Including NGG or CCN, k=23.
constexpr auto k = 23;

// Permit at most this many N characters per 23-mer.
//
// For example, 'ACATAGTGACGTTNAAACATNG' is permitted because the trailing TNG
// matches TGG which is a valid PAM for CRISPR.
//
// 'ACATAGTNACGTTNAAACATNG' contains 3 N characters and is not permitted.
//
// The choice max_N=2 is somewhat arbitrary.  Any value <5 will work pretty well
// for the human genome hg38.fa.   Values much greater than 5 should be avoided,
// as may result in signifficant explosion of output size.
//
constexpr int max_N = 2;

// Turn off expand_N_variants if you wish to see actual N characters in the output.
// For the default "true", each N will be expanded into all possible variants.
// For example, for 20-mer 'ACATAGTNACGTTNAAACA', there will be 16
// variants present in the output:
//
//      'ACATAGTAACGTTAAAACA'    'ACATAGTAACGTTCAAACA'
//      'ACATAGTAACGTTGAAACA'    'ACATAGTAACGTTTAAACA'
//      'ACATAGTCACGTTAAAACA'    'ACATAGTCACGTTCAAACA'
//      'ACATAGTCACGTTGAAACA'    'ACATAGTCACGTTTAAACA'
//      'ACATAGTGACGTTAAAACA'    'ACATAGTGACGTTCAAACA'
//      'ACATAGTGACGTTGAAACA'    'ACATAGTGACGTTTAAACA'
//      'ACATAGTTACGTTAAAACA'    'ACATAGTTACGTTCAAACA'
//      'ACATAGTTACGTTGAAACA'    'ACATAGTTACGTTTAAACA'
//
// This expansion permits the downstream pipeline not to need handling of N wildcards.
//
constexpr bool expand_N_variants = true;


// A 20-mer over the ACGTN alphabet is representable by a 60 bit integer,
// with each base letter represented by 3 consecutive bits.
constexpr int bits_per_base = 3;
constexpr int64_t lsb = 1;
constexpr int64_t base_mask = (lsb << bits_per_base) - lsb;


// Compute full-complement-mask (subtracting a code from that mask produces
// the code's complement).
constexpr int64_t fcm(int len) {
    return len == 1 ? 6 : (6 | (fcm(len - 1) << bits_per_base));
}

// This runs at compile time.
constexpr int64_t complement_mask = fcm(k - 3);


int64_t bitcode_for_base(const char c) {
    // The encodings are chosen such that you get the DNA complement
    // simply by subtracting from complement_mask.  Also lex order.
    switch (c) {
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'N':
            return 3;
        case 'G':
            return 4;
        case 'T':
            return 5;
    }
    assert(false);
}


char base_for_bitcode(const int64_t single_base_code) {
    switch (single_base_code) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 3:
            return 'N';
        case 4:
            return 'G';
        case 5:
            return 'T';
    }
    assert(false);
}


// This optimization reduces runtime for hg38.fa from 2m20sec to just 2m.
// But don't forget to run init_encoding.
int encoding[1 << (sizeof(char) * 8)];
void init_encoding() {
    encoding['A'] = bitcode_for_base('A');
    encoding['C'] = bitcode_for_base('C');
    encoding['N'] = bitcode_for_base('N');
    encoding['G'] = bitcode_for_base('G');
    encoding['T'] = bitcode_for_base('T');
}


template<int len>
int64_t encode(const char* buf) {
    int64_t code = 0;
    for (int i=0;  i < len;  ++i) {
        const int64_t base_code = encoding[buf[i]];
        code |= ((base_code & base_mask) << (bits_per_base * (len - i - 1)));
    }
    return code;
}


void decode(char* buf, const int len, const int64_t code) {
    for (int i=0;  i < len;  ++i) {
        const int64_t base_code = (code >> (bits_per_base * (len - i - 1))) & base_mask;
        buf[i] = base_for_bitcode(base_code);
    }
}


int64_t complement(const int64_t code) {
    return complement_mask - code;
}


constexpr bool forward_direction = false;
constexpr bool reverse_complement = true;


template <bool direction>
void emit(vector<int64_t>& results, const char* guide) {
    auto code = encode<k - 3>(guide);
    if (direction == reverse_complement) {
        code = complement(code);
    }
    results.push_back(code);
}


template <bool direction>
void emit_all_variants(vector<int64_t>& results, const char* guide, const int num_N_to_expand) {
    // assert count of 'N' characters in rc[0...k-3] equals num_N_to_expand
    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    char guide_variant[k - 3];  // not 0-terminated
    const int number_variants = (1 << (2 * num_N_to_expand));  // 4 power num_N
    for (int n_code = 0;  n_code < number_variants;  ++n_code) {
        int pos = k - 3;
        int code_pos = 0;
        do {
            --pos;
            char c = guide[pos];
            if (c == 'N') {
                // intepret n_code as a 4-digit number base 4, with each
                // digit representing a base
                c = bases[(n_code >> (code_pos * 2)) & 0x3];
                ++code_pos;
            }
            guide_variant[pos] = c;
        } while (pos);
        emit<direction>(results, guide_variant);
    }
}


template <bool direction>
int index(const int j) {
    return (direction == forward_direction) ? j : k - 1 - j;
}


template <bool direction, char cog>
void try_match(vector<int64_t>& results, const char* bufi) {
    char guide[k - 3];  // not 0 terminated
    int count[1 << (sizeof(char) * 8)];
    assert(cog == 'C' || cog == 'G');
    count[cog] = 0;
    count['N'] = 0;
    constexpr auto opposite_direction = !direction;
    count[bufi[index<opposite_direction>(0)]] += 1;
    count[bufi[index<opposite_direction>(1)]] += 1;
    if (count[cog] + count['N'] == 2) {
        const int pam_N = count['N'];
        for (int j = 0;  j < k - 3;  ++j) {
            const char c = bufi[index<direction>(j)];
            guide[j] = c;
            count[c] += 1;
        }
        // allow up to this many total N, including PAM site
        if (count['N'] <= max_N) {
            const int num_N_to_expand = count['N'] - pam_N;
            if (expand_N_variants && num_N_to_expand > 0) {
                emit_all_variants<direction>(results, guide, num_N_to_expand);
            } else {
                // common case, no Ns to expand in this 20-mer
                emit<direction>(results, guide);
            }
        }
    }
}


void scan_for_kmers(vector<int64_t>& results, const char* buf, size_t len) {
    assert(k <= 24);
    for (int i = 0;  i <= len - k;  ++i) {
        // match ...GG, or ...GN, or ...NG, or ...NN
        try_match<forward_direction, 'G'>(results, buf + i);
        // match CC..., or CN..., or NC..., or NN...
        try_match<reverse_complement, 'C'>(results, buf + i);
    }
}


// Return number of milliseconds elapsed since Jan 1, 1970 00:00 GMT.
long unixtime() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

void scan_stdin() {

    init_encoding();

    vector<int64_t> results;

    constexpr auto STRIDE_SIZE = 32 * 1024 * 1024;

    // to scan for k-mers, consecutive read windows must overlap by k-1 characters
    constexpr auto BUFFER_SIZE = STRIDE_SIZE + k - 1;

    // using c++ vector provides transparent memory management
    vector<char> buffer(BUFFER_SIZE);
    char* window = buffer.data();

    uintmax_t lines = 0;
    uintmax_t bases = 0;
    uintmax_t guides = 0;

    // chromosome comments are lines that begin with '>'
    // we filter those out
    bool chrm_comment = false;
    int overlap = 0;
    auto t = unixtime();
    auto t_last_print = t;
    auto t_start = t;

    while (true) {

        assert(0 <= overlap);
        assert(overlap < k);

        const ssize_t bytes_read = read(fileno(stdin), window + overlap, STRIDE_SIZE);

        // end of file
        if (bytes_read == 0) {
            break;
        }

        // Not sure when this happens.
        if (bytes_read == (ssize_t) -1) {
            throw runtime_error("ooops");
        }

        // Convert to uppercase and filter out chromosome comments.
        int len = overlap;
        for (int i = overlap;  i < bytes_read + overlap;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
                ++lines;
                chrm_comment = false;
            } else if (!(chrm_comment)) {
                if (c == '>') {
                    chrm_comment = true;
                } else {
                    // compact
                    window[len++] = c;
                }
            }
        }

        bases += (len - overlap);

        if (len < k) {
            // There are no k-mers in the current buffer.
            // This is likely the end of the file and the very last iteration.
            overlap = len;
        } else {
            // window now starts with the last k-1 bases from the previous read,
            // plus all bases from the current read
            scan_for_kmers(results, window, len);
            // overlap the last k-1 characters by moving them to the start of the window
            overlap = k - 1;
            for (int i=0;  i < overlap;  ++i) {
                window[i] = window[len - overlap + i];
            }
        }

        t = unixtime();
        if (t - t_last_print > 10000) {
            cerr << "Progress update "  << (t - t_start) / 1000 << " seconds after start." << endl;
            cerr << "Lines so far: "  << lines  << endl;
            cerr << "Bases so far: "  << bases  << endl;
            cerr << endl;
            t_last_print = t;
        }

    }

    cerr << "Finished reading input."  << endl;
    cerr << "Total lines: "  << lines  << endl;
    cerr << "Total bases: "  << bases  << endl;

    // If there are tons of duplicates, we may benefit from sorting each batch
    // and then merging incrementally with c++ algorithm set_union,
    // rather than doing a huge sort at the end.   Parallelizing, esp on GPU,
    // could yield phenomenal speedup if we ever need to run this program fast.
    cerr << "Sorting " << results.size() << " candidate guides." << endl;
    sort(results.begin(), results.end());

    // 0 is not a valid code
    int64_t last = 0;
    for (auto it = results.begin();  it != results.end();  ++it) {
        if (*it != last) {
            ++guides;
        }
        last = *it;
    }

    cerr << "Outputting " << guides << " unique guides." << endl;
    last = 0;
    char obuf[k-1];
    obuf[k-2] = 0;
    obuf[k-3] = '\n';
    for (auto it = results.begin();  it != results.end();  ++it) {
        if (*it != last) {
            decode(obuf, k-3, *it);
            cout << obuf;
        }
        last = *it;
    }

}


void silent_tests() {
    const char* kmer                = "ACGTGGTGGCAATGCACGGT";
    const char* kmer_complement     = "TGCACCACCGTTACGTGCCA";
    const char* kmer_2wc            = "ACGTGGTGGCAATNCACGGN";
    const char* kmer_2wc_complement = "TGCACCACCGTTANGTGCCN";
    const char* kmer_own_complement = "NNNNNNNNNNNNNNNNNNNN";
    const char* nonse               = "abcdefghijklmnopqrst";
    char buf[k-2];
    assert(strlen(kmer) == k - 3);
    assert(strlen(kmer_complement) == k - 3);
    assert(strlen(kmer_2wc) == k - 3);
    assert(strlen(kmer_2wc_complement) == k - 3);
    assert(strlen(nonse) == k - 3);
    assert(strlen(nonse) == sizeof(buf) - 1);
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    decode(buf, k-3, encode<k - 3>(kmer));
    assert(0 == strcmp(buf, kmer));
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    decode(buf, k-3, complement(encode<k - 3>(kmer)));
    assert(0 == strcmp(buf, kmer_complement));
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    decode(buf, k-3, complement(encode<k - 3>(kmer_complement)));
    assert(0 == strcmp(buf, kmer));
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    decode(buf, k-3, encode<k - 3>(kmer_2wc));
    assert(0 == strcmp(buf, kmer_2wc));
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    decode(buf, k-3, complement(encode<k - 3>(kmer_2wc)));
    assert(0 == strcmp(buf, kmer_2wc_complement));
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    decode(buf, k-3, complement(encode<k - 3>(kmer_2wc_complement)));
    assert(0 == strcmp(buf, kmer_2wc));
    // ----------------
    strcpy(buf, nonse);
    assert(0 == strcmp(buf, nonse));
    auto oc_code = encode<k - 3>(kmer_own_complement);
    assert(oc_code == complement(oc_code));
    // ----------------
    vector<int64_t> expansions;
    const char* expected[] = {
        "ACGTGGTGGCAATACACGGA",
        "ACGTGGTGGCAATACACGGC",
        "ACGTGGTGGCAATACACGGG",
        "ACGTGGTGGCAATACACGGT",
        "ACGTGGTGGCAATCCACGGA",
        "ACGTGGTGGCAATCCACGGC",
        "ACGTGGTGGCAATCCACGGG",
        "ACGTGGTGGCAATCCACGGT",
        "ACGTGGTGGCAATGCACGGA",
        "ACGTGGTGGCAATGCACGGC",
        "ACGTGGTGGCAATGCACGGG",
        "ACGTGGTGGCAATGCACGGT",
        "ACGTGGTGGCAATTCACGGA",
        "ACGTGGTGGCAATTCACGGC",
        "ACGTGGTGGCAATTCACGGG",
        "ACGTGGTGGCAATTCACGGT",
    };
    emit_all_variants<false>(expansions, kmer_2wc, 2);
    for (int i=0; i<16; ++i) {
        decode(buf, k - 3, expansions[i]);
        assert(0 == strcmp(expected[i], buf));
    }
    if (false) {
        // This is good but not silent.
        for (auto it = expansions.begin();   it != expansions.end();   ++it) {
            strcpy(buf, nonse);
            assert(0 == strcmp(buf, nonse));
            decode(buf, k - 3, *it);
            cout << buf << endl;
        }
    }
}


int main(int argc, char** argv){
    init_encoding();
    silent_tests();
    scan_stdin();
    return 0;
}
