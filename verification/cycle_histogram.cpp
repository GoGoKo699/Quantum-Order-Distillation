// Exact signed permutation-pair histogram for a multilinear order polynomial.
// Input: support_size degree, then coefficient followed by one-based word labels.
// Output: cycle_count signed_weight for every cycle count, including zero weights.
// Matrix products are read left to right. Uses C_word=(0,word_1,...,word_m).
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Term {
    std::int64_t coefficient;
    std::array<unsigned char, 32> cycle{};
    std::array<unsigned char, 32> inverse{};
};

int main() {
    int support, degree;
    if (!(std::cin >> support >> degree)) return 1;
    if (degree < 1 || degree > 30 || support < 1) return 2;
    std::vector<Term> terms(support);
    for (auto &term : terms) {
        std::cin >> term.coefficient;
        int previous=0;
        std::uint32_t seen=1;
        for (int k=0; k<degree; ++k) {
            int label; std::cin >> label;
            if (label < 1 || label > degree || (seen & (1u<<label))) return 3;
            term.cycle[previous]=static_cast<unsigned char>(label);
            seen |= 1u<<label;
            previous=label;
        }
        term.cycle[previous]=0;
        for (int k=0; k<=degree; ++k) term.inverse[term.cycle[k]]=k;
    }
    std::vector<std::int64_t> histogram(degree+2,0);
    for (int i=0; i<support; ++i) {
        for (int j=i; j<support; ++j) {
            std::uint32_t seen=0;
            int count=0;
            for (int k=0; k<=degree; ++k) {
                if (seen & (1u<<k)) continue;
                ++count;
                int next=k;
                do {
                    seen |= 1u<<next;
                    next=terms[i].cycle[terms[j].inverse[next]];
                } while (!(seen & (1u<<next)));
            }
            histogram[count] += (i==j ? 1 : 2)*terms[i].coefficient*terms[j].coefficient;
        }
    }
    for (int count=1; count<=degree+1; ++count)
        std::cout << count << " " << histogram[count] << "\n";
}
