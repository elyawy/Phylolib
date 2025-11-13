#pragma once

#include <random>
#include <utility>
#include <vector>
#include <stack>
#include <iostream>

// Helper struct to split a single RNG call into uniform int and float
template<uint32_t K>
struct SplitRNG {
    static_assert(K == 4 || K == 20 || K == 61, "K must be 4, 20, or 61");
    
    static constexpr int bits_for_int = 12;
    static constexpr uint64_t int_mask = (1ULL << bits_for_int) - 1;
    static constexpr int bits_for_float = 64 - bits_for_int;
    static constexpr double float_divisor = (double)(1ULL << bits_for_float);
    static constexpr uint64_t threshold = ((1ULL << bits_for_int) / K) * K;
    
    uint32_t integer;
    double uniform_float;
    
    template<typename RNG>
    explicit SplitRNG(RNG& rng) {
        uint64_t rng_val = rng();
        uint64_t candidate = rng_val & int_mask;
        
        while (candidate >= threshold) {
            rng_val = rng();
            candidate = rng_val & int_mask;
        }
        
        integer = (uint32_t)(candidate % K);
        uniform_float = (double)(rng_val >> bits_for_int) / float_divisor;
    }
};

// Implementation of Vose's alias method, given a list of probabilities
// representing a finite distribution builds a datastructure in O(n) time
// that allows O(1) drawing from the distribution.
template<size_t N>
class DiscreteNDistribution
{
private:
    std::array<double, N> probabilities_;
    std::array<int, N> alias_;

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
        auto split = SplitRNG<N>(rng);

        int die_roll = split.integer;  // Uses lower bits for N
        double coin_flip = split.uniform_float;  // Uses different bits
        
        int result = alias_[die_roll];
        if (coin_flip < probabilities_[die_roll]) result = die_roll;
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

    // static void setSeed(int seed) {
    //     rng_ = std::mt19937_64(seed);
    // }

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