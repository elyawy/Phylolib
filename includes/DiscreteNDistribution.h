#pragma once

#include <random>
#include <utility>
#include <vector>
#include <stack>
#include <iostream>

// Helper struct to split a single RNG call into uniform int and float
// K can only be a power of two, So we only need to mask the lower bits for the integer part, and use the upper bits for the float part.
// For example, if K=32, we can use the lower 5 bits for the integer (0-31) and the upper bits for the float (0.0-1.0).
// Since we get 64bits we can split up into two 32bits.
// each 32bits can be used for the integer part and the float part respectively, so we can use the lower 5 bits of the integer part for K=32, and the upper 27 bits for the float part, which gives us a very high precision for the float part.
// more specifically, we get at worst 1/(2^26) precision for the float part, which is more than enough for our use case.
template<uint8_t K>
struct SplitRNG {
    static_assert(K == 4 || K == 32 || K == 64, "K must be 4, 32, or 64");
    
    static constexpr uint64_t int_mask = K - 1;  // Mask to extract the integer part
    static constexpr uint8_t bits_for_int = __builtin_ctz(K);  // Number of bits needed for the integer part
    static constexpr double divisor = 1ULL << (64 - bits_for_int);  // Divisor to scale the float part to [0.0, 1.0)
    
    uint8_t integer;
    double uniform;

    
    explicit SplitRNG(uint64_t rng_val) {

        integer = rng_val & int_mask;                          // take low bits
        uniform= (rng_val >> bits_for_int) / divisor;  // take remaining high bits
    }
};

// Implementation of Vose's alias method, given a list of probabilities
// representing a finite distribution builds a datastructure in O(n) time
// that allows O(1) drawing from the distribution.
template<size_t N>
class DiscreteNDistribution
{
private:
    std::array<double,   N> probabilities_;
    std::array<uint8_t, N> alias_;

public:
    DiscreteNDistribution<N>(const std::vector<double> &probabilities, double normalizingFactor=1.0) {
        assert(probabilities.size() == N);  // safety check


        std::stack<std::pair<int, double>> small_;
        std::stack<std::pair<int, double>> large_;


        for(int i = 0; i < N; i++) {
            double scaled_prob = N*probabilities[i]*normalizingFactor;
            std::pair<int, double> current_prob(i, scaled_prob);
            if (current_prob.second < 1) {
                small_.push(current_prob);
            } else {
                large_.push(current_prob);
            }
        }
        while (!small_.empty() && !large_.empty()) {
            std::pair<int, double> s = small_.top();
            std::pair<int, double> l = large_.top();
            small_.pop();
            large_.pop();

            probabilities_[s.first] = s.second;
            alias_[s.first] = l.first;

            l.second = (l.second + s.second) - 1;
            if(l.second < 1) {
                small_.push(l);
            } else {
                large_.push(l);
            }
        }
        while (!large_.empty()) {
            std::pair<int, double> l = large_.top();
            large_.pop();
            probabilities_[l.first] = 1;
        }
        while (!small_.empty()) {
            std::pair<int, double> s = small_.top();
            small_.pop();
            probabilities_[s.first] = 1;
        }
    }

    template<typename RngType = std::mt19937_64>
    int drawSample(RngType &rng) {
        // split 64 bits into two 32 bits
        uint32_t intToSplit;

        uint64_t full64Bits = rng();

        SplitRNG<N> split(full64Bits);

        uint8_t die_roll = split.integer;  // Uses lower bits for N
        double coin_flip = split.uniform;  // Uses different bits
        
        int result = alias_[die_roll];
        int keep_original = (coin_flip < probabilities_[die_roll]); 

        int mask = -keep_original; // will be 0xFF..FF if keep_original is true, 0x00..00 if false
        result = (result & ~mask) | (die_roll & mask); // if keep_original is true, keep die_roll, else use alias result

        return result + 1;
    }


    void printTable() {
        for(auto &i: probabilities_){
            std::cout << i << " ";
        }
        std::cout << "\n";
        for(auto &i: alias_){
            std::cout << i << " ";
        }
        std::cout << "\n";
    }



    std::vector<std::pair<double, int>> getTable() {
        std::vector<std::pair<double, int>> prob_alias_table;
        for(int i=0; i < probabilities_.size(); i++){
            std::pair<double, int> current_item (probabilities_[i], alias_[i]);
            prob_alias_table.push_back(current_item);
        }
        return prob_alias_table;
    }

    double getProb(int i) {
        return probabilities_[i];
    }
};